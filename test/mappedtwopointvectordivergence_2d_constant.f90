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

  exit_code = mappedtwopointvectordivergence_2d_constant()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function mappedtwopointvectordivergence_2d_constant() result(r)
    !! Verifies that the mapped 2-D split-form divergence of a constant
    !! physical-space two-point vector field is identically zero on a
    !! uniform Cartesian mesh.
    !!
    !! A constant physical flux has a constant contravariant projection;
    !! each D-matrix column sum is zero, so the split-form sum must vanish.
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

    ! Set contravariant two-point fluxes for a constant physical flux F=(1,1).
    ! For constant F, F_EC(sL,sR) = F for all pairs.
    ! Contravariant: Fc^r(nn,i,j) = sum_d avg(Ja^r_d) * F_d
    ! On a uniform Cartesian mesh, avg(Ja^r_d) = Ja^r_d (constant per element).
    block
      integer :: nn,i,j,d,iEl,iVar
      real(prec) :: Fc
      do concurrent(nn=1:controlDegree+1,i=1:controlDegree+1, &
                    j=1:controlDegree+1,iEl=1:mesh%nElem,iVar=1:nvar)

        ! xi^1: pair (i,j)-(nn,j)
        Fc = 0.0_prec
        do d = 1,2
          Fc = Fc+0.5_prec*( &
               geometry%dsdx%interior(i,j,iEl,1,d,1)+ &
               geometry%dsdx%interior(nn,j,iEl,1,d,1))*1.0_prec
        enddo
        f%interior(nn,i,j,iEl,iVar,1) = Fc

        ! xi^2: pair (i,j)-(i,nn)
        Fc = 0.0_prec
        do d = 1,2
          Fc = Fc+0.5_prec*( &
               geometry%dsdx%interior(i,j,iEl,1,d,2)+ &
               geometry%dsdx%interior(i,nn,iEl,1,d,2))*1.0_prec
        enddo
        f%interior(nn,i,j,iEl,iVar,2) = Fc

      enddo
    endblock
    call f%UpdateDevice()

#ifdef ENABLE_GPU
    call f%MappedDivergence(df%interior_gpu)
#else
    call f%MappedDivergence(df%interior)
#endif
    call df%UpdateHost()

    ! Add surface term for constant flux on the mapped element.
    ! The boundary normal flux = F_phys . nhat * nScale at each face.
    ! For F=(1,1): fbn = (nhat_x + nhat_y)*nScale at each face node.
    block
      integer :: i,j,iEl,iVar
      real(prec) :: fbn_east,fbn_west,fbn_north,fbn_south,jac
      do concurrent(i=1:controlDegree+1,j=1:controlDegree+1, &
                    iEl=1:mesh%nElem,iVar=1:nvar)

        jac = geometry%J%interior(i,j,iEl,1)
        fbn_east = geometry%nScale%boundary(j,2,iEl,1)* &
                   (geometry%nHat%boundary(j,2,iEl,1,1)+geometry%nHat%boundary(j,2,iEl,1,2))
        fbn_west = geometry%nScale%boundary(j,4,iEl,1)* &
                   (geometry%nHat%boundary(j,4,iEl,1,1)+geometry%nHat%boundary(j,4,iEl,1,2))
        fbn_north = geometry%nScale%boundary(i,3,iEl,1)* &
                    (geometry%nHat%boundary(i,3,iEl,1,1)+geometry%nHat%boundary(i,3,iEl,1,2))
        fbn_south = geometry%nScale%boundary(i,1,iEl,1)* &
                    (geometry%nHat%boundary(i,1,iEl,1,1)+geometry%nHat%boundary(i,1,iEl,1,2))

        df%interior(i,j,iEl,iVar) = df%interior(i,j,iEl,iVar)+ &
                                    (interp%bMatrix(i,2)*fbn_east+interp%bMatrix(i,1)*fbn_west)/(interp%qWeights(i)*jac)+ &
                                    (interp%bMatrix(j,2)*fbn_north+interp%bMatrix(j,1)*fbn_south)/(interp%qWeights(j)*jac)

      enddo
    endblock

    df%interior = abs(df%interior-0.0_prec)

    print*,"absmax error :",maxval(df%interior)
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

  endfunction mappedtwopointvectordivergence_2d_constant

endprogram test
