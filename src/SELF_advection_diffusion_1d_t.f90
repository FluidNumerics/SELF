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

module self_advection_diffusion_1d_t

  use self_model
  use self_dgmodel1d
  use self_mesh

  implicit none

  type,extends(dgmodel1d) :: advection_diffusion_1d_t
    real(prec) :: nu ! diffusion coefficient
    real(prec) :: u ! constant velocity

  contains
    procedure :: setboundarycondition => setboundarycondition_advection_diffusion_1d_t
    procedure :: setgradientboundarycondition => setgradientboundarycondition_advection_diffusion_1d_t
    procedure :: riemannsolver => riemannsolver_advection_diffusion_1d_t
    procedure :: fluxmethod => fluxmethod_advection_diffusion_1d_t
    procedure :: CalculateEntropy => CalculateEntropy_advection_diffusion_1d_t


  endtype advection_diffusion_1d_t

contains

  subroutine CalculateEntropy_advection_diffusion_1d_t(this)
    implicit none
    class(advection_diffusion_1d_t),intent(inout) :: this
    ! Local
    integer :: iel, i, ivar
    real(prec) :: e,s,J

    e = 0.0_prec
    do ivar = 1,this%solution%nvar
      do iel = 1,this%geometry%nelem
        do i = 1,this%solution%interp%N
          J = this%geometry%dxds%interior(i,iel,1)
          s = this%solution%interior(i,iel,ivar)
          e = e + 0.5_prec*s*s*J
        enddo
      enddo
    enddo

    this%entropy = e

  endsubroutine CalculateEntropy_advection_diffusion_1d_t

  subroutine setboundarycondition_advection_diffusion_1d_t(this)
    ! Here, we use the pre-tendency method to calculate the
    ! derivative of the solution using a bassi-rebay method
    ! We then do a boundary interpolation and side exchange
    ! on the gradient field
    implicit none
    class(advection_diffusion_1d_t),intent(inout) :: this
    ! local
    integer :: ivar
    integer :: N,nelem

    nelem = this%geometry%nelem ! number of elements in the mesh
    N = this%solution%interp%N ! polynomial degree

    do ivar = 1,this%solution%nvar

      ! left-most boundary
      this%solution%extBoundary(1,1,ivar) = &
        this%solution%boundary(2,nelem,ivar)

      ! right-most boundary
      this%solution%extBoundary(2,nelem,ivar) = &
        this%solution%boundary(1,1,ivar)

    enddo
 
  endsubroutine setboundarycondition_advection_diffusion_1d_t

  subroutine setgradientboundarycondition_advection_diffusion_1d_t(this)
    ! Here, we set the boundary conditions for the
    ! solution and the solution gradient at the left
    ! and right most boundaries.
    !
    ! Here, we use periodic boundary conditions
    implicit none
    class(advection_diffusion_1d_t),intent(inout) :: this
    ! local
    integer :: ivar
    integer :: nelem

    nelem = this%geometry%nelem ! number of elements in the mesh

    do ivar = 1,this%solution%nvar

      ! left-most boundary
      this%solutionGradient%extBoundary(1,1,ivar) = &
        this%solutionGradient%boundary(2,nelem,ivar)

      ! right-most boundary
      this%solutionGradient%extBoundary(2,nelem,ivar) = &
        this%solutionGradient%boundary(1,1,ivar)

    enddo


  endsubroutine setgradientboundarycondition_advection_diffusion_1d_t

  subroutine fluxmethod_advection_diffusion_1d_t(this)
    implicit none
    class(advection_diffusion_1d_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: ivar
    integer :: i
    real(prec) :: u,nu,f,dfdx

    u = this%u
    nu = this%nu
    !$omp target
    !$omp teams
    !$omp loop collapse(3)
    do ivar = 1,this%solution%nvar
      do iel = 1,this%mesh%nelem
        do i = 1,this%solution%interp%N+1

          f = this%solution%interior(i,iel,ivar)
          dfdx = this%solutionGradient%interior(i,iel,ivar)

          this%flux%interior(i,iel,ivar) = u*f-nu*dfdx ! advective flux + diffusive flux

        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target

  endsubroutine fluxmethod_advection_diffusion_1d_t

  subroutine riemannsolver_advection_diffusion_1d_t(this)
    ! this method uses an linear upwind solver for the
    ! advective flux and the bassi-rebay method for the
    ! diffusive fluxes
    implicit none
    class(advection_diffusion_1d_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: ivar
    integer :: iside
    real(prec) :: fin,fout,dfavg,u,nhat

    u = this%u
    !$omp target
    !$omp teams
    !$omp loop collapse(3)
    do ivar = 1,this%solution%nvar
      do iel = 1,this%mesh%nelem
        do iside = 1,2

          ! set the normal velocity
          if(iside == 1) then
            nhat = -1.0_prec
          else
            nhat = 1.0_prec
          endif

          fin = this%solution%boundary(iside,iel,ivar) ! interior solution
          fout = this%solution%extboundary(iside,iel,ivar) ! exterior solution
          dfavg = this%solutionGradient%avgboundary(iside,iel,ivar) ! average solution gradient (with direction taken into account)

          this%flux%boundarynormal(iside,iel,ivar) = 0.5_prec*(u*nhat*(fin+fout)+abs(u*nhat)*(fin-fout))- & ! advective flux
                                               this%nu*dfavg*nhat ! diffusive flux

        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target

  endsubroutine riemannsolver_advection_diffusion_1d_t

endmodule self_advection_diffusion_1d_t
