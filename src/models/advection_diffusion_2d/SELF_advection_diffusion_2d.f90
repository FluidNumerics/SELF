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

module self_advection_diffusion_2d

  use self_model
  use self_dgmodel2d
  use self_mesh

  implicit none

  type,extends(dgmodel2d) :: advection_diffusion_2d
    real(prec) :: nu ! diffusion coefficient
    real(prec) :: u ! constant x-component of velocity
    real(prec) :: v ! constant y-component of velocity

  contains
    procedure :: pretendency => pretendency_advection_diffusion_2d
    procedure :: setboundarycondition => setboundarycondition_advection_diffusion_2d
    procedure :: riemannsolver => riemannsolver_advection_diffusion_2d
    procedure :: fluxmethod => fluxmethod_advection_diffusion_2d

  endtype advection_diffusion_2d

  ! Remember, the order of operations for the tendency calculation is
  !
  !    solution % BoundaryInterp(this % gpuAccel)
  !    solution % SideExchange(this % mesh,this % decomp,this % gpuAccel)
  !    PreTendency()
  !    SetBoundaryCondition()
  !    SourceMethod()
  !    RiemannSolver()
  !    FluxMethod()
  !    CalculateFluxDivergence()

contains

  subroutine pretendency_advection_diffusion_2d(this)
    ! Here, we use the pre-tendency method to calculate the
    ! derivative of the solution using a bassi-rebay method
    ! We then do a boundary interpolation and side exchange
    ! on the gradient field
    implicit none
    class(advection_diffusion_2d),intent(inout) :: this
    ! local
    integer :: i,ivar,iEl,j,e2

    !$omp target map(from: this % mesh % sideInfo) map(tofrom: this % solution % extBoundary)
    !$omp teams loop bind(teams) collapse(3)
    do ivar = 1,this%solution%nvar
      do iEl = 1,this%solution%nElem ! Loop over all elements
        do j = 1,4 ! Loop over all sides

          !bcid = this % mesh % sideInfo(5,j,iEl) ! Boundary Condition ID
          e2 = this%mesh%sideInfo(3,j,iEl) ! Neighboring Element ID

          if(e2 == 0) then
            !$omp loop bind(parallel)
            do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
              this%solution%extBoundary(i,j,iEl,ivar) = 0.0_prec
            enddo
          endif

        enddo
      enddo
    enddo
    !$omp end target

    ! calculate the averages of the solutions on the element
    ! boundaries and store is this % solution % avgBoundary
    call this%solution%AverageSides()

    ! calculate the derivative using the weak form
    ! With the boundary attribute now storing the side average, this
    ! amounts to a bassi-rebay method for the gradient.
    this%solutionGradient%interior = this%solution%DGGradient(this%geometry)

    ! interpolate the solutiongradient to the element boundaries
    call this%solutionGradient%BoundaryInterp()

    ! perform the side exchange to populate the solutionGradient % extBoundary attribute
    call this%solutionGradient%SideExchange(this%mesh,this%decomp)

    ! Re-compute the solution%boundary attribute so that we don't use the avgboundary in the hyperbolic flux
    call this%solution%BoundaryInterp()

  endsubroutine pretendency_advection_diffusion_2d

  subroutine setboundarycondition_advection_diffusion_2d(this)
    ! Here, we set the boundary conditions for the
    ! solution and the solution gradient at the left
    ! and right most boundaries.
    !
    ! Here, we use periodic boundary conditions
    implicit none
    class(advection_diffusion_2d),intent(inout) :: this
    ! local
    integer :: i,ivar,iEl,j,e2

    !$omp target map(from: this % mesh % sideInfo) map(tofrom: this % solutionGradient % extBoundary)
    !$omp teams loop collapse(3)
    do ivar = 1,this%solution%nvar
      do iEl = 1,this%solution%nElem ! Loop over all elements
        do j = 1,4 ! Loop over all sides

          !bcid = this % mesh % sideInfo(5,j,iEl) ! Boundary Condition ID
          e2 = this%mesh%sideInfo(3,j,iEl) ! Neighboring Element ID

          if(e2 == 0) then

            do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
              this%solutionGradient%extBoundary(i,j,iEl,ivar,1:2) = 0.0_prec
            enddo

          endif

        enddo
      enddo
    enddo
    !$omp end target

    call this%solutionGradient%AverageSides()

  endsubroutine setboundarycondition_advection_diffusion_2d

  subroutine fluxmethod_advection_diffusion_2d(this)
    implicit none
    class(advection_diffusion_2d),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: ivar
    integer :: i
    integer :: j
    real(prec) :: u,v,nu,f,dfdx,dfdy

    u = this%u
    v = this%v
    nu = this%nu
    !$omp target map(to:this % solution % interior, this % solutionGradient % interior) map(from:this % flux % interior)
    !$omp teams loop collapse(4)
    do ivar = 1,this%solution%nvar
      do iel = 1,this%mesh%nelem
        do j = 1,this%solution%interp%N+1
          do i = 1,this%solution%interp%N+1

            f = this%solution%interior(i,j,iel,ivar)
            dfdx = this%solutionGradient%interior(i,j,iel,ivar,1)
            dfdy = this%solutionGradient%interior(i,j,iel,ivar,2)

            this%flux%interior(i,j,iel,ivar,1) = u*f-nu*dfdx ! advective flux + diffusive flux (x-component)
            this%flux%interior(i,j,iel,ivar,2) = v*f-nu*dfdy ! advective flux + diffusive flux (y-component)

          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine fluxmethod_advection_diffusion_2d

  subroutine riemannsolver_advection_diffusion_2d(this)
    ! this method uses an linear upwind solver for the
    ! advective flux and the bassi-rebay method for the
    ! diffusive fluxes
    implicit none
    class(advection_diffusion_2d),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: ivar
    integer :: j
    integer :: i
    real(prec) :: fin,fout,dfdn,un
    real(prec) :: nx,ny,nmag

    !$omp target map(to:this % geometry % nHat % boundary, this % solutionGradient % avgBoundary) &
    !$omp& map(to:this % solution % boundary, this % solution % extBoundary) &
    !$omp& map(to:this % geometry % nscale % boundary) map(from: this % flux % boundaryNormal)
    !$omp teams loop collapse(4)
    do ivar = 1,this%solution%nvar
      do iEl = 1,this%solution%nElem
        do j = 1,4
          do i = 1,this%solution%interp%N+1

            ! Get the boundary normals on cell edges from the mesh geometry
            nx = this%geometry%nHat%boundary(i,j,iEl,1,1)
            ny = this%geometry%nHat%boundary(i,j,iEl,1,2)

            un = this%u*nx+this%v*ny
            dfdn = this%solutionGradient%avgBoundary(i,j,iEl,ivar,1)*nx+ &
                   this%solutionGradient%avgBoundary(i,j,iEl,ivar,2)*ny

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
    !$omp end target

  endsubroutine riemannsolver_advection_diffusion_2d

endmodule self_advection_diffusion_2d
