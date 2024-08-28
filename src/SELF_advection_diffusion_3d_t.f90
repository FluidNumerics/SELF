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
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module self_advection_diffusion_3d_t

  use self_model
  use self_dgmodel3d
  use self_mesh

  implicit none

  type,extends(dgmodel3d) :: advection_diffusion_3d_t
    real(prec) :: nu ! diffusion coefficient
    real(prec) :: u ! constant x-component of velocity
    real(prec) :: v ! constant y-component of velocity
    real(prec) :: w ! constant z-component of velocity

  contains

    procedure :: setboundarycondition => setboundarycondition_advection_diffusion_3d_t
    procedure :: setgradientboundarycondition => setgradientboundarycondition_advection_diffusion_3d_t
    procedure :: riemannsolver => riemannsolver_advection_diffusion_3d_t
    procedure :: fluxmethod => fluxmethod_advection_diffusion_3d_t
    procedure :: CalculateEntropy => CalculateEntropy_advection_diffusion_3d_t

  endtype advection_diffusion_3d_t

contains

subroutine CalculateEntropy_advection_diffusion_3d_t(this)
    implicit none
    class(advection_diffusion_3d_t),intent(inout) :: this
    ! Local
    integer :: iel,i,j,k,ivar
    real(prec) :: e,s,jac

    e = 0.0_prec
    do ivar = 1,this%solution%nvar
      do iel = 1,this%geometry%nelem
        do k = 1,this%solution%interp%N+1
          do j = 1,this%solution%interp%N+1
            do i = 1,this%solution%interp%N+1
              jac = this%geometry%J%interior(i,j,k,iel,1)
              s = this%solution%interior(i,j,k,iel,ivar)
              e = e + 0.5_prec*s*s*jac
            enddo
          enddo
        enddo
      enddo
    enddo

    this%entropy = e

  endsubroutine CalculateEntropy_advection_diffusion_3d_t

  subroutine setboundarycondition_advection_diffusion_3d_t(this)
    ! Here, we use the pre-tendency method to calculate the
    ! derivative of the solution using a bassi-rebay method
    ! We then do a boundary interpolation and side exchange
    ! on the gradient field
    implicit none
    class(advection_diffusion_3d_t),intent(inout) :: this
    ! local
    integer :: i,j,ivar,iEl,k,e2

    !$omp target
    !$omp teams loop bind(teams) collapse(3)
    do ivar = 1,this%solution%nvar
      do iEl = 1,this%solution%nElem ! Loop over all elements
        do k = 1,6 ! Loop over all sides

          !bcid = this % mesh % sideInfo(5,k,iEl) ! Boundary Condition ID
          e2 = this%mesh%sideInfo(3,k,iEl) ! Neighboring Element ID

          if(e2 == 0) then
            !$omp loop bind(parallel) collapse(2)
            do j = 1,this%solution%interp%N+1 ! Loop over quadrature point
              do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
                this%solution%extBoundary(i,j,k,iEl,iVar) = 0.0_prec
              enddo
            enddo
          endif
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine setboundarycondition_advection_diffusion_3d_t

  subroutine setgradientboundarycondition_advection_diffusion_3d_t(this)
    ! Here, we set the boundary conditions for the
    ! solution and the solution gradient at the left
    ! and right most boundaries.
    !
    ! Here, we use periodic boundary conditions
    implicit none
    class(advection_diffusion_3d_t),intent(inout) :: this
    ! local
    integer :: i,j,ivar,iEl,k,e2

    !$omp target
    !$omp teams loop bind(teams) collapse(3)
    do ivar = 1,this%solution%nvar
      do iEl = 1,this%solution%nElem ! Loop over all elements
        do k = 1,6 ! Loop over all sides

          !bcid = this % mesh % sideInfo(5,k,iEl) ! Boundary Condition ID
          e2 = this%mesh%sideInfo(3,k,iEl) ! Neighboring Element ID

          if(e2 == 0) then
            !$omp loop bind(parallel) collapse(2)
            do j = 1,this%solution%interp%N+1 ! Loop over quadrature point
              do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
                this%solutionGradient%extBoundary(i,j,k,iEl,iVar,1:3) = &
                  this%solutionGradient%boundary(i,j,k,iEl,iVar,1:3)
              enddo
            enddo
          endif
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine setgradientboundarycondition_advection_diffusion_3d_t

  subroutine fluxmethod_advection_diffusion_3d_t(this)
    implicit none
    class(advection_diffusion_3d_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: ivar
    integer :: i
    integer :: j
    integer :: k
    real(prec) :: u,v,w,nu,f,dfdx,dfdy,dfdz

    u = this%u
    v = this%v
    w = this%w
    nu = this%nu
    !$omp target
    !$omp teams loop collapse(5)
    do ivar = 1,this%solution%nvar
      do iel = 1,this%mesh%nelem
        do k = 1,this%solution%interp%N+1
          do j = 1,this%solution%interp%N+1
            do i = 1,this%solution%interp%N+1

              f = this%solution%interior(i,j,k,iel,ivar)
              dfdx = this%solutionGradient%interior(i,j,k,iel,ivar,1)
              dfdy = this%solutionGradient%interior(i,j,k,iel,ivar,2)
              dfdz = this%solutionGradient%interior(i,j,k,iel,ivar,3)

              this%flux%interior(i,j,k,iel,ivar,1) = u*f-nu*dfdx ! advective flux + diffusive flux (x-component)
              this%flux%interior(i,j,k,iel,ivar,2) = v*f-nu*dfdy ! advective flux + diffusive flux (y-component)
              this%flux%interior(i,j,k,iel,ivar,3) = w*f-nu*dfdz ! advective flux + diffusive flux (z-component)

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine fluxmethod_advection_diffusion_3d_t

  subroutine riemannsolver_advection_diffusion_3d_t(this)
    ! this method uses an linear upwind solver for the
    ! advective flux and the bassi-rebay method for the
    ! diffusive fluxes
    implicit none
    class(advection_diffusion_3d_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: ivar
    integer :: k
    integer :: i,j
    real(prec) :: fin,fout,dfdn,un
    real(prec) :: nx,ny,nz,nmag

    !$omp target
    !$omp teams loop collapse(4)
    do ivar = 1,this%solution%nvar
      do iEl = 1,this%solution%nElem
        do k = 1,6
          do j = 1,this%solution%interp%N+1
            do i = 1,this%solution%interp%N+1

              ! Get the boundary normals on cell edges from the mesh geometry
              nx = this%geometry%nHat%boundary(i,j,k,iEl,1,1)
              ny = this%geometry%nHat%boundary(i,j,k,iEl,1,2)
              nz = this%geometry%nHat%boundary(i,j,k,iEl,1,3)

              un = this%u*nx+this%v*ny+this%w*nz
              dfdn = this%solutionGradient%boundary(i,j,k,iEl,iVar,1)*nx+ &
                     this%solutionGradient%boundary(i,j,k,iEl,iVar,2)*ny+ &
                     this%solutionGradient%boundary(i,j,k,iEl,iVar,3)*nz

              fin = this%solution%boundary(i,j,k,iEl,iVar) ! interior solution
              fout = this%solution%extboundary(i,j,k,iEl,iVar) ! exterior solution

              nmag = this%geometry%nScale%boundary(i,j,k,iEl,1)

              this%flux%boundaryNormal(i,j,k,iEl,1) = (0.5_prec*( &
                                                       un*(fin+fout)+abs(un)*(fin-fout))- & ! advective flux
                                                       this%nu*dfdn)*nmag
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine riemannsolver_advection_diffusion_3d_t

endmodule self_advection_diffusion_3d_t
