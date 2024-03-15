module self_advection_diffusion_1d

use self_model
use self_dgmodel1d
use self_mesh

implicit none

  type, extends(dgmodel1d) :: advection_diffusion_1d
    real(prec) :: nu ! diffusion coefficient
    real(prec) :: u  ! constant velocity

    contains
    procedure :: pretendency => pretendency_advection_diffusion_1d
    procedure :: setboundarycondition => setboundarycondition_advection_diffusion_1d
    procedure :: riemannsolver => riemannsolver_advection_diffusion_1d
    procedure :: fluxmethod => fluxmethod_advection_diffusion_1d

  end type advection_diffusion_1d

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

  subroutine pretendency_advection_diffusion_1d(this)
  ! Here, we use the pre-tendency method to calculate the
  ! derivative of the solution using a bassi-rebay method
  ! We then do a boundary interpolation and side exchange
  ! on the gradient field
  implicit none
  class(advection_diffusion_1d), intent(inout) :: this
  ! local
  integer :: ivar
  integer :: N, nelem

    nelem = this % geometry % nelem ! number of elements in the mesh
    N = this % solution % interp % N ! polynomial degree

    do ivar = 1, this % solution % nvar

      ! left-most boundary
      this % solution % extBoundary(1,1,ivar) = &
        this % solution % boundary(2,nelem,ivar)

      ! right-most boundary
      this % solution % extBoundary(2,nelem,ivar) = &
        this % solution % boundary(1,1,ivar)

    enddo

    ! calculate the averages of the solutions on the element
    ! boundaries and store is this % solution % avgBoundary
    call this % solution % BassiRebaySides()

    ! calculate the derivative using the bassi-rebay form
    call this % solution % BRDerivative(this % geometry, &
            this % solutionGradient)

    ! interpolate the solutiongradient to the element boundaries
    call this % solutionGradient % BoundaryInterp()
  
    ! perform the side exchange to populate the 
    ! solutionGradient % extBoundary attribute
    call this % solutionGradient % SideExchange(this % mesh, &
           this % decomp) 

    call this % solutionGradient % BassiRebaySides()

  end subroutine pretendency_advection_diffusion_1d

  subroutine setboundarycondition_advection_diffusion_1d(this)
  ! Here, we set the boundary conditions for the 
  ! solution and the solution gradient at the left
  ! and right most boundaries.
  ! 
  ! Here, we use periodic boundary conditions
  implicit none
  class(advection_diffusion_1d), intent(inout) :: this
  ! local
  integer :: ivar
  integer :: N, nelem

    nelem = this % geometry % nelem ! number of elements in the mesh
    N = this % solution % interp % N ! polynomial degree

    do ivar = 1, this % solution % nvar

      ! left-most boundary
      this% solutionGradient % extBoundary(1,1,ivar) = &
        this % solutionGradient % boundary(2,nelem,ivar)

      this% solutionGradient % avgBoundary(1,1,ivar) = &
        -0.5_prec*(this% solutionGradient % extBoundary(1,1,ivar) + &
        this% solutionGradient % boundary(ivar,1,1))

      ! right-most boundary
      this % solutionGradient % extBoundary(2,nelem,ivar) = &
        this % solutionGradient % boundary(1,1,ivar)

      this% solutionGradient % avgBoundary(2,nelem,ivar) = &
        0.5_prec*(this% solutionGradient % extBoundary(2,nelem,ivar) + &
        this% solutionGradient % boundary(2,nelem,ivar))
    enddo

  end subroutine setboundarycondition_advection_diffusion_1d

  subroutine fluxmethod_advection_diffusion_1d(this)
  implicit none
  class(advection_diffusion_1d), intent(inout) :: this
  ! Local
  integer :: iel
  integer :: ivar
  integer :: i
  real(prec) :: u, nu, f, dfdx

    u = this % u
    nu = this % nu
    do ivar = 1, this % solution % nvar
      do iel = 1, this % mesh % nelem
        do i = 1, this % solution % interp % N+1

          f = this % solution % interior(i,iel,ivar)
          dfdx = this % solutionGradient % interior(i,iel,ivar)

          this % flux % interior(i,iel,ivar) = u*f - nu*dfdx  ! advective flux + diffusive flux

        enddo
      enddo
    enddo

  end subroutine fluxmethod_advection_diffusion_1d

  subroutine riemannsolver_advection_diffusion_1d(this)
  ! this method uses an linear upwind solver for the
  ! advective flux and the bassi-rebay method for the 
  ! diffusive fluxes
  implicit none
  class(advection_diffusion_1d), intent(inout) :: this
  ! Local
  integer :: iel
  integer :: ivar
  integer :: iside
  real(prec) :: fin, fout, dfavg, un

    call this % solutionGradient % BassiRebaySides()

    do ivar = 1, this % solution % nvar
      do iel = 1, this % mesh % nelem
        do iside = 1, 2

          ! set the normal velocity
          if( iside == 1 )then
            un = -this % u
          else
            un = this % u
          endif


            fin = this % solution % boundary(iside,iel,ivar) ! interior solution
            fout = this % solution % extboundary(iside,iel,ivar) ! exterior solution
            dfavg = this % solutionGradient % avgboundary(iside,iel,ivar) ! average solution gradient (with direction taken into account)

            this % flux % boundary(iside,iel,ivar) = 0.5_prec*( un*(fin + fout) + abs(un)*(fin - fout) ) -& ! advective flux
                                                                this % nu*dfavg ! diffusive flux

        enddo
      enddo
    enddo

  end subroutine riemannsolver_advection_diffusion_1d

  
end module self_advection_diffusion_1d