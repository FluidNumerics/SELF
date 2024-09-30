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

module self_advection_diffusion_2d_t

  use self_model
  use self_dgmodel2d
  use self_mesh

  implicit none

  type,extends(dgmodel2d) :: advection_diffusion_2d_t
    real(prec) :: nu ! diffusion coefficient
    real(prec) :: u ! constant x-component of velocity
    real(prec) :: v ! constant y-component of velocity

  contains
    procedure :: setboundarycondition => setboundarycondition_advection_diffusion_2d_t
    procedure :: setgradientboundarycondition => setgradientboundarycondition_advection_diffusion_2d_t
    procedure :: riemannsolver => riemannsolver_advection_diffusion_2d_t
    procedure :: fluxmethod => fluxmethod_advection_diffusion_2d_t
    procedure :: CalculateEntropy => CalculateEntropy_advection_diffusion_2d_t

  endtype advection_diffusion_2d_t

contains

  subroutine CalculateEntropy_advection_diffusion_2d_t(this)
    implicit none
    class(advection_diffusion_2d_t),intent(inout) :: this
    ! Local
    integer :: iel,i,j,ivar,ierror
    real(prec) :: e,s,jac

    e = 0.0_prec
    do ivar = 1,this%solution%nvar
      do iel = 1,this%geometry%nelem
        do j = 1,this%solution%interp%N+1
          do i = 1,this%solution%interp%N+1
            jac = this%geometry%J%interior(i,j,iel,1)
            s = this%solution%interior(i,j,iel,ivar)
            e = e+0.5_prec*s*s*jac
          enddo
        enddo
      enddo
    enddo

    if(this%mesh%decomp%mpiEnabled) then
      call mpi_allreduce(e, &
                         this%entropy, &
                         1, &
                         this%mesh%decomp%mpiPrec, &
                         MPI_SUM, &
                         this%mesh%decomp%mpiComm, &
                         iError)
    else
      this%entropy = e
    endif

  endsubroutine CalculateEntropy_advection_diffusion_2d_t

  subroutine setboundarycondition_advection_diffusion_2d_t(this)
    !! Boundary conditions for the solution are set to
    !! 0 for the external state to provide radiation type
    !! boundary conditions.
    implicit none
    class(advection_diffusion_2d_t),intent(inout) :: this
    ! local
    integer :: i,ivar,iEl,j,e2

    do ivar = 1,this%solution%nvar
      do iEl = 1,this%solution%nElem ! Loop over all elements
        do j = 1,4 ! Loop over all sides

          !bcid = this % mesh % sideInfo(5,j,iEl) ! Boundary Condition ID
          e2 = this%mesh%sideInfo(3,j,iEl) ! Neighboring Element ID

          if(e2 == 0) then
            do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
              this%solution%extBoundary(i,j,iEl,ivar) = 0.0_prec
            enddo
          endif

        enddo
      enddo
    enddo

  endsubroutine setboundarycondition_advection_diffusion_2d_t

  subroutine setgradientboundarycondition_advection_diffusion_2d_t(this)
    !! Boundary conditions on the solution gradient are set
    !! to prolong the solution gradient through the boundaries
    implicit none
    class(advection_diffusion_2d_t),intent(inout) :: this
    ! local
    integer :: i,ivar,iEl,j,e2

    do ivar = 1,this%solution%nvar
      do iEl = 1,this%solution%nElem ! Loop over all elements
        do j = 1,4 ! Loop over all sides

          !bcid = this % mesh % sideInfo(5,j,iEl) ! Boundary Condition ID
          e2 = this%mesh%sideInfo(3,j,iEl) ! Neighboring Element ID

          if(e2 == 0) then

            do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
              this%solutionGradient%extBoundary(i,j,iEl,ivar,1:2) = &
                this%solutionGradient%boundary(i,j,iEl,ivar,1:2)
            enddo

          endif

        enddo
      enddo
    enddo

  endsubroutine setgradientboundarycondition_advection_diffusion_2d_t

  ! function flux_func_advection_diffusion_2d( this, f, dfdx, dfdy ) result(flux)
  !   class(advection_diffusion_2d_t) :: this
  !   real(prec) :: f(1:this%nvar)
  !   real(prec) :: dfdx(1:this%nvar)
  !   real(prec) :: dfdy(1:this%nvar)
  !   real(prec) :: flux(1:this%nvar,1:2)

  !   flux(1:this%nvar,1) = this%u*f-this%nu*dfdx ! advective flux + diffusive flux (x-component)
  !   flux(1:this%nvar,1) = this%v*f-this%nu*dfdy ! advective flux + diffusive flux (x-component)

  ! end function

  subroutine fluxmethod_advection_diffusion_2d_t(this)
    implicit none
    class(advection_diffusion_2d_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: ivar
    integer :: i
    integer :: j
    real(prec) :: u,v,nu,f,dfdx,dfdy

    u = this%u
    v = this%v
    nu = this%nu
    do ivar = 1,this%solution%nvar
      do iel = 1,this%mesh%nelem
        do j = 1,this%solution%interp%N+1
          do i = 1,this%solution%interp%N+1

            f = this%solution%interior(i,j,iel,ivar)
            dfdx = this%solutionGradient%interior(i,j,iel,ivar,1)
            dfdy = this%solutionGradient%interior(i,j,iel,ivar,2)

            !this%flux%interior(i,j,iel,:,:) = this%flux_func( f, dfdx, dfdy )

            this%flux%interior(i,j,iel,ivar,1) = u*f-nu*dfdx ! advective flux + diffusive flux (x-component)
            this%flux%interior(i,j,iel,ivar,2) = v*f-nu*dfdy ! advective flux + diffusive flux (y-component)

          enddo
        enddo
      enddo
    enddo

  endsubroutine fluxmethod_advection_diffusion_2d_t

  subroutine riemannsolver_advection_diffusion_2d_t(this)
    ! this method uses an linear upwind solver for the
    ! advective flux and the bassi-rebay method for the
    ! diffusive fluxes
    implicit none
    class(advection_diffusion_2d_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: ivar
    integer :: j
    integer :: i
    real(prec) :: fin,fout,dfdn,un
    real(prec) :: nx,ny,nmag

    do ivar = 1,this%solution%nvar
      do iEl = 1,this%solution%nElem
        do j = 1,4
          do i = 1,this%solution%interp%N+1

            ! Get the boundary normals on cell edges from the mesh geometry
            nx = this%geometry%nHat%boundary(i,j,iEl,1,1)
            ny = this%geometry%nHat%boundary(i,j,iEl,1,2)

            un = this%u*nx+this%v*ny
            dfdn = this%solutionGradient%avgboundary(i,j,iEl,ivar,1)*nx+ &
                   this%solutionGradient%avgboundary(i,j,iEl,ivar,2)*ny

            fin = this%solution%boundary(i,j,iel,ivar) ! interior solution
            fout = this%solution%extboundary(i,j,iel,ivar) ! exterior solution

            nmag = this%geometry%nScale%boundary(i,j,iEl,1)

            this%flux%boundaryNormal(i,j,iEl,1) = (0.5_prec*( &
                                                   un*(fin+fout)+abs(un)*(fin-fout))- & ! advective flux
                                                   this%nu*dfdn)*nmag

          enddo
        enddo
      enddo
    enddo

  endsubroutine riemannsolver_advection_diffusion_2d_t

endmodule self_advection_diffusion_2d_t
