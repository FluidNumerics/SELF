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

  exit_code = mappedtwopointvectordivergence_3d_linear()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function mappedtwopointvectordivergence_3d_linear() result(r)
    !! Verifies that the mapped 3-D split-form divergence is exact for
    !!   V = x * e_x + y * e_y + z * e_z
    !! on a uniform Cartesian block mesh, for which the exact divergence is 3.
    !!
    !! Contravariant two-point fluxes are set using arithmetic means of mesh
    !! node coordinates projected onto the averaged metric terms (Winters et al.).
    !! Each direction r: Fc^r = sum_d avg(Ja^r_d) * F_EC_d with correct pair.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_MappedScalar_3D
    use SELF_MappedTwoPointVector_3D

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
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(MappedTwoPointVector3D) :: f
    type(MappedScalar3D) :: df
    character(LEN=255) :: WORKSPACE
    integer :: i,j,k,nn,iel,ivar

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS_LOBATTO, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block3D/Block3D_mesh.h5")

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nelem)
    call df%Init(interp,nvar,mesh%nelem)
    call f%AssociateGeometry(geometry)

    ! Contravariant two-point fluxes for V = (x, y, z), divergence = 3.
    ! For V=(x,y,z), F_EC_d(sL,sR) = avg(x_d) for each physical direction d.
    ! Each direction r: Fc^r = sum_d avg(Ja^r_d) * F_EC_d with correct pair.
    do ivar = 1,nvar
      do iel = 1,mesh%nelem
        do k = 1,controlDegree+1
          do j = 1,controlDegree+1
            do i = 1,controlDegree+1
              do nn = 1,controlDegree+1
                ! xi^1: pair (i,j,k)-(nn,j,k)
                f%interior(nn,i,j,k,iel,ivar,1) = &
                  0.5_prec*(geometry%dsdx%interior(i,j,k,iel,1,1,1)+ &
                            geometry%dsdx%interior(nn,j,k,iel,1,1,1))* &
                  0.5_prec*(geometry%x%interior(i,j,k,iel,1,1)+ &
                            geometry%x%interior(nn,j,k,iel,1,1))+ &
                  0.5_prec*(geometry%dsdx%interior(i,j,k,iel,1,2,1)+ &
                            geometry%dsdx%interior(nn,j,k,iel,1,2,1))* &
                  0.5_prec*(geometry%x%interior(i,j,k,iel,1,2)+ &
                            geometry%x%interior(nn,j,k,iel,1,2))+ &
                  0.5_prec*(geometry%dsdx%interior(i,j,k,iel,1,3,1)+ &
                            geometry%dsdx%interior(nn,j,k,iel,1,3,1))* &
                  0.5_prec*(geometry%x%interior(i,j,k,iel,1,3)+ &
                            geometry%x%interior(nn,j,k,iel,1,3))
                ! xi^2: pair (i,j,k)-(i,nn,k)
                f%interior(nn,i,j,k,iel,ivar,2) = &
                  0.5_prec*(geometry%dsdx%interior(i,j,k,iel,1,1,2)+ &
                            geometry%dsdx%interior(i,nn,k,iel,1,1,2))* &
                  0.5_prec*(geometry%x%interior(i,j,k,iel,1,1)+ &
                            geometry%x%interior(i,nn,k,iel,1,1))+ &
                  0.5_prec*(geometry%dsdx%interior(i,j,k,iel,1,2,2)+ &
                            geometry%dsdx%interior(i,nn,k,iel,1,2,2))* &
                  0.5_prec*(geometry%x%interior(i,j,k,iel,1,2)+ &
                            geometry%x%interior(i,nn,k,iel,1,2))+ &
                  0.5_prec*(geometry%dsdx%interior(i,j,k,iel,1,3,2)+ &
                            geometry%dsdx%interior(i,nn,k,iel,1,3,2))* &
                  0.5_prec*(geometry%x%interior(i,j,k,iel,1,3)+ &
                            geometry%x%interior(i,nn,k,iel,1,3))
                ! xi^3: pair (i,j,k)-(i,j,nn)
                f%interior(nn,i,j,k,iel,ivar,3) = &
                  0.5_prec*(geometry%dsdx%interior(i,j,k,iel,1,1,3)+ &
                            geometry%dsdx%interior(i,j,nn,iel,1,1,3))* &
                  0.5_prec*(geometry%x%interior(i,j,k,iel,1,1)+ &
                            geometry%x%interior(i,j,nn,iel,1,1))+ &
                  0.5_prec*(geometry%dsdx%interior(i,j,k,iel,1,2,3)+ &
                            geometry%dsdx%interior(i,j,nn,iel,1,2,3))* &
                  0.5_prec*(geometry%x%interior(i,j,k,iel,1,2)+ &
                            geometry%x%interior(i,j,nn,iel,1,2))+ &
                  0.5_prec*(geometry%dsdx%interior(i,j,k,iel,1,3,3)+ &
                            geometry%dsdx%interior(i,j,nn,iel,1,3,3))* &
                  0.5_prec*(geometry%x%interior(i,j,k,iel,1,3)+ &
                            geometry%x%interior(i,j,nn,iel,1,3))
              enddo
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
    ! For V=(x,y,z): F_phys.nhat = x*nhat_x + y*nhat_y + z*nhat_z at each face node.
    ! 3D sides: 1=Bottom, 2=South, 3=East, 4=North, 5=West, 6=Top
    block
      integer :: ii,jj,kk,iiEl,iiVar
      real(prec) :: fbn_e,fbn_w,fbn_n,fbn_s,fbn_t,fbn_b,jac,xb,yb,zb
      do concurrent(ii=1:controlDegree+1,jj=1:controlDegree+1, &
                    kk=1:controlDegree+1,iiEl=1:mesh%nElem,iiVar=1:nvar)

        jac = geometry%J%interior(ii,jj,kk,iiEl,1)

        ! East (side 3): boundary node (jj,kk) at xi^1=+1
        xb = geometry%x%boundary(jj,kk,3,iiEl,1,1)
        yb = geometry%x%boundary(jj,kk,3,iiEl,1,2)
        zb = geometry%x%boundary(jj,kk,3,iiEl,1,3)
        fbn_e = (xb*geometry%nHat%boundary(jj,kk,3,iiEl,1,1)+ &
                 yb*geometry%nHat%boundary(jj,kk,3,iiEl,1,2)+ &
                 zb*geometry%nHat%boundary(jj,kk,3,iiEl,1,3))* &
                geometry%nScale%boundary(jj,kk,3,iiEl,1)
        ! West (side 5): boundary node (jj,kk) at xi^1=-1
        xb = geometry%x%boundary(jj,kk,5,iiEl,1,1)
        yb = geometry%x%boundary(jj,kk,5,iiEl,1,2)
        zb = geometry%x%boundary(jj,kk,5,iiEl,1,3)
        fbn_w = (xb*geometry%nHat%boundary(jj,kk,5,iiEl,1,1)+ &
                 yb*geometry%nHat%boundary(jj,kk,5,iiEl,1,2)+ &
                 zb*geometry%nHat%boundary(jj,kk,5,iiEl,1,3))* &
                geometry%nScale%boundary(jj,kk,5,iiEl,1)
        ! North (side 4): boundary node (ii,kk) at xi^2=+1
        xb = geometry%x%boundary(ii,kk,4,iiEl,1,1)
        yb = geometry%x%boundary(ii,kk,4,iiEl,1,2)
        zb = geometry%x%boundary(ii,kk,4,iiEl,1,3)
        fbn_n = (xb*geometry%nHat%boundary(ii,kk,4,iiEl,1,1)+ &
                 yb*geometry%nHat%boundary(ii,kk,4,iiEl,1,2)+ &
                 zb*geometry%nHat%boundary(ii,kk,4,iiEl,1,3))* &
                geometry%nScale%boundary(ii,kk,4,iiEl,1)
        ! South (side 2): boundary node (ii,kk) at xi^2=-1
        xb = geometry%x%boundary(ii,kk,2,iiEl,1,1)
        yb = geometry%x%boundary(ii,kk,2,iiEl,1,2)
        zb = geometry%x%boundary(ii,kk,2,iiEl,1,3)
        fbn_s = (xb*geometry%nHat%boundary(ii,kk,2,iiEl,1,1)+ &
                 yb*geometry%nHat%boundary(ii,kk,2,iiEl,1,2)+ &
                 zb*geometry%nHat%boundary(ii,kk,2,iiEl,1,3))* &
                geometry%nScale%boundary(ii,kk,2,iiEl,1)
        ! Top (side 6): boundary node (ii,jj) at xi^3=+1
        xb = geometry%x%boundary(ii,jj,6,iiEl,1,1)
        yb = geometry%x%boundary(ii,jj,6,iiEl,1,2)
        zb = geometry%x%boundary(ii,jj,6,iiEl,1,3)
        fbn_t = (xb*geometry%nHat%boundary(ii,jj,6,iiEl,1,1)+ &
                 yb*geometry%nHat%boundary(ii,jj,6,iiEl,1,2)+ &
                 zb*geometry%nHat%boundary(ii,jj,6,iiEl,1,3))* &
                geometry%nScale%boundary(ii,jj,6,iiEl,1)
        ! Bottom (side 1): boundary node (ii,jj) at xi^3=-1
        xb = geometry%x%boundary(ii,jj,1,iiEl,1,1)
        yb = geometry%x%boundary(ii,jj,1,iiEl,1,2)
        zb = geometry%x%boundary(ii,jj,1,iiEl,1,3)
        fbn_b = (xb*geometry%nHat%boundary(ii,jj,1,iiEl,1,1)+ &
                 yb*geometry%nHat%boundary(ii,jj,1,iiEl,1,2)+ &
                 zb*geometry%nHat%boundary(ii,jj,1,iiEl,1,3))* &
                geometry%nScale%boundary(ii,jj,1,iiEl,1)

        df%interior(ii,jj,kk,iiEl,iiVar) = df%interior(ii,jj,kk,iiEl,iiVar)+ &
                                           (interp%bMatrix(ii,2)*fbn_e+interp%bMatrix(ii,1)*fbn_w)/(interp%qWeights(ii)*jac)+ &
                                           (interp%bMatrix(jj,2)*fbn_n+interp%bMatrix(jj,1)*fbn_s)/(interp%qWeights(jj)*jac)+ &
                                           (interp%bMatrix(kk,2)*fbn_t+interp%bMatrix(kk,1)*fbn_b)/(interp%qWeights(kk)*jac)

      enddo
    endblock

    print*,"absmax (nabla.F) :",maxval(df%interior)
    df%interior = abs(df%interior-3.0_prec)
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

  endfunction mappedtwopointvectordivergence_3d_linear

endprogram test
