! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2024 Fluid Numerics LLC
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
!    the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

program test

  implicit none
  integer :: exit_code

  exit_code = mappedtwopointvectordivergence_2d_linear()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function mappedtwopointvectordivergence_2d_linear() result(r)
    !! Verifies that the mapped 2-D split-form divergence is exact for
    !!   V = x * e_x + y * e_y
    !! on a uniform Cartesian block mesh, for which the exact divergence is 2.
    !!
    !! Physical two-point fluxes are set using arithmetic means of mesh node
    !! coordinates. On a Cartesian mesh the off-diagonal metric terms vanish,
    !! so:
    !!   interior(nn,i,j,iel,ivar,1) = (x_{i,j} + x_{nn,j}) / 2
    !!   interior(nn,i,j,iel,ivar,2) = (y_{i,j} + y_{i,nn}) / 2
    !!
    !! MappedDivergence projects these onto contravariant directions using
    !! averaged metric terms (Winters et al.) and divides by the Jacobian.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_MappedScalar_2D
    use SELF_MappedTwoPointVector_2D

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 1
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(MappedTwoPointVector2D) :: f
    type(MappedScalar2D) :: df
    character(LEN=255) :: WORKSPACE
    integer :: i,j,nn,iel,ivar

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS_LOBATTO, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D_mesh.h5")

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nelem)
    call df%Init(interp,nvar,mesh%nelem)
    call f%AssociateGeometry(geometry)

    ! Contravariant two-point fluxes for V = (x, y), divergence = 2.
    ! For V=(x,y), F_EC(sL,sR) = ((xL+xR)/2, (yL+yR)/2).
    ! Each direction r: Fc^r = avg(Ja^r) . F_EC with correct pair.
    do ivar = 1,nvar
      do iel = 1,mesh%nelem
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            do nn = 1,controlDegree+1
              ! xi^1: pair (i,j)-(nn,j)
              f%interior(nn,i,j,iel,ivar,1) = &
                0.5_prec*(geometry%dsdx%interior(i,j,iel,1,1,1)+ &
                          geometry%dsdx%interior(nn,j,iel,1,1,1))* &
                0.5_prec*(geometry%x%interior(i,j,iel,1,1)+ &
                          geometry%x%interior(nn,j,iel,1,1))+ &
                0.5_prec*(geometry%dsdx%interior(i,j,iel,1,2,1)+ &
                          geometry%dsdx%interior(nn,j,iel,1,2,1))* &
                0.5_prec*(geometry%x%interior(i,j,iel,1,2)+ &
                          geometry%x%interior(nn,j,iel,1,2))
              ! xi^2: pair (i,j)-(i,nn)
              f%interior(nn,i,j,iel,ivar,2) = &
                0.5_prec*(geometry%dsdx%interior(i,j,iel,1,1,2)+ &
                          geometry%dsdx%interior(i,nn,iel,1,1,2))* &
                0.5_prec*(geometry%x%interior(i,j,iel,1,1)+ &
                          geometry%x%interior(i,nn,iel,1,1))+ &
                0.5_prec*(geometry%dsdx%interior(i,j,iel,1,2,2)+ &
                          geometry%dsdx%interior(i,nn,iel,1,2,2))* &
                0.5_prec*(geometry%x%interior(i,j,iel,1,2)+ &
                          geometry%x%interior(i,nn,iel,1,2))
            enddo
          enddo
        enddo
      enddo
    enddo
    call f%UpdateDevice()

#ifdef ENABLE_GPU
    call f%MappedDivergence(df%interior_gpu)
#else
    call f%MappedDivergence(df%interior)
#endif
    call df%UpdateHost()

    ! Add surface term: (1/J)*M^{-1}*B^T * (F_phys . nhat * nScale)
    ! For V=(x,y): F_phys.nhat = x*nhat_x + y*nhat_y at each face node.
    block
      integer :: ii,jj,iiEl,iiVar
      real(prec) :: fbn_e,fbn_w,fbn_n,fbn_s,jac,xb,yb
      do concurrent(ii=1:controlDegree+1,jj=1:controlDegree+1, &
                    iiEl=1:mesh%nElem,iiVar=1:nvar)

        jac = geometry%J%interior(ii,jj,iiEl,1)

        ! East (side 2): boundary node j=jj at xi^1=+1
        xb = geometry%x%boundary(jj,2,iiEl,1,1)
        yb = geometry%x%boundary(jj,2,iiEl,1,2)
        fbn_e = (xb*geometry%nHat%boundary(jj,2,iiEl,1,1)+ &
                 yb*geometry%nHat%boundary(jj,2,iiEl,1,2))* &
                geometry%nScale%boundary(jj,2,iiEl,1)
        ! West (side 4)
        xb = geometry%x%boundary(jj,4,iiEl,1,1)
        yb = geometry%x%boundary(jj,4,iiEl,1,2)
        fbn_w = (xb*geometry%nHat%boundary(jj,4,iiEl,1,1)+ &
                 yb*geometry%nHat%boundary(jj,4,iiEl,1,2))* &
                geometry%nScale%boundary(jj,4,iiEl,1)
        ! North (side 3)
        xb = geometry%x%boundary(ii,3,iiEl,1,1)
        yb = geometry%x%boundary(ii,3,iiEl,1,2)
        fbn_n = (xb*geometry%nHat%boundary(ii,3,iiEl,1,1)+ &
                 yb*geometry%nHat%boundary(ii,3,iiEl,1,2))* &
                geometry%nScale%boundary(ii,3,iiEl,1)
        ! South (side 1)
        xb = geometry%x%boundary(ii,1,iiEl,1,1)
        yb = geometry%x%boundary(ii,1,iiEl,1,2)
        fbn_s = (xb*geometry%nHat%boundary(ii,1,iiEl,1,1)+ &
                 yb*geometry%nHat%boundary(ii,1,iiEl,1,2))* &
                geometry%nScale%boundary(ii,1,iiEl,1)

        df%interior(ii,jj,iiEl,iiVar) = df%interior(ii,jj,iiEl,iiVar)+ &
                                        (interp%bMatrix(ii,2)*fbn_e+interp%bMatrix(ii,1)*fbn_w)/(interp%qWeights(ii)*jac)+ &
                                        (interp%bMatrix(jj,2)*fbn_n+interp%bMatrix(jj,1)*fbn_s)/(interp%qWeights(jj)*jac)

      enddo
    endblock

    print*,"absmax (nabla.F) :",maxval(df%interior)
    df%interior = abs(df%interior-2.0_prec)
    print*,"absmax error     :",maxval(df%interior)

    if(maxval(df%interior) <= tolerance) then
      r = 0
    else
      r = 1
    endif

    call f%DissociateGeometry()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()
    call f%free()
    call df%free()

  endfunction mappedtwopointvectordivergence_2d_linear

endprogram test
