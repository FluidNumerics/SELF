module advection_diffusion_model_1d

use self_model
use self_model1d
use self_mesh

implicit none

  type, extends(model1d) :: advdiff1d
    real(prec) :: nu ! diffusion coefficient
    real(prec) :: u  ! constant velocity

    contains
    procedure :: pretendency => pretendency_advdiff1d
    procedure :: setboundarycondition => setboundarycondition_advdiff1d
!    procedure :: sourcemethod => sourcemethod_advdiff1d
    procedure :: riemannsolver => riemannsolver_advdiff1d
    procedure :: fluxmethod => fluxmethod_advdiff1d

  end type advdiff1d

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

  subroutine pretendency_advdiff1d(this)
  ! Here, we use the pre-tendency method to calculate the
  ! derivative of the solution using a bassi-rebay method
  ! We then do a boundary interpolation and side exchange
  ! on the gradient field
  implicit none
  class(advdiff1d), intent(inout) :: this
  ! local
  integer :: ivar
  integer :: N, nelem

    nelem = this % geometry % nelem ! number of elements in the mesh
    N = this % solution % interp % N ! polynomial degree

    do ivar = 1, this % solution % nvar

      ! left-most boundary
      this % solution % extBoundary % hostdata(ivar,1,1) = &
        this % solution % boundary % hostdata(ivar,2,nelem)

      ! right-most boundary
      this % solution % extBoundary % hostdata(ivar,2,nelem) = &
        this % solution % boundary % hostdata(ivar,1,1)

    enddo

    ! calculate the averages of the solutions on the element
    ! boundaries and store is this % solution % avgBoundary
    call this % solution % BassiRebaySides(this % gpuaccel)

    ! calculate the derivative using the bassi-rebay form
    call this % solution % Derivative(this % geometry, &
            this % solutionGradient, selfWeakBRForm, &
            this % gpuaccel)

    ! interpolate the solutiongradient to the element boundaries
    call this % solutionGradient % BoundaryInterp(this % gpuaccel)
  
    ! perform the side exchange to populate the 
    ! solutionGradient % extBoundary attribute
    call this % solutionGradient % SideExchange(this % mesh, &
           this % decomp, this % gpuaccel) 

    call this % solutionGradient % BassiRebaySides(this % gpuaccel)

  end subroutine pretendency_advdiff1d

  subroutine setboundarycondition_advdiff1d(this)
  ! Here, we set the boundary conditions for the 
  ! solution and the solution gradient at the left
  ! and right most boundaries.
  ! 
  ! Here, we use periodic boundary conditions
  implicit none
  class(advdiff1d), intent(inout) :: this
  ! local
  integer :: ivar
  integer :: N, nelem

    nelem = this % geometry % nelem ! number of elements in the mesh
    N = this % solution % interp % N ! polynomial degree

    do ivar = 1, this % solution % nvar

      ! left-most boundary
      this% solutionGradient % extBoundary % hostdata(ivar,1,1) = &
        this % solutionGradient % boundary % hostdata(ivar,2,nelem)

      this% solutionGradient % avgBoundary % hostdata(ivar,1,1) = &
        -0.5_prec*(this% solutionGradient % extBoundary % hostdata(ivar,1,1) + &
        this% solutionGradient % boundary % hostdata(ivar,1,1))

      ! right-most boundary
      this % solutionGradient % extBoundary % hostdata(ivar,2,nelem) = &
        this % solutionGradient % boundary % hostdata(ivar,1,1)

      this% solutionGradient % avgBoundary % hostdata(ivar,2,nelem) = &
        0.5_prec*(this% solutionGradient % extBoundary % hostdata(ivar,2,nelem) + &
        this% solutionGradient % boundary % hostdata(ivar,2,nelem))
    enddo

  end subroutine setboundarycondition_advdiff1d

  subroutine fluxmethod_advdiff1d(this)
  implicit none
  class(advdiff1d), intent(inout) :: this
  ! Local
  integer :: iel
  integer :: ivar
  integer :: i
  real(prec) :: u, nu, f, dfdx

    u = this % u
    nu = this % nu
    do iel = 1, this % mesh % nelem
      do ivar = 1, this % solution % nvar
        do i = 0, this % solution % interp % N

          f = this % solution % interior % hostdata(i,ivar,iel)
          dfdx = this % solutionGradient % interior % hostdata(i,ivar,iel)

          this % flux % interior % hostdata(i,ivar,iel) = u*f - nu*dfdx  ! advective flux + diffusive flux

        enddo
      enddo
    enddo

  end subroutine fluxmethod_advdiff1d

  subroutine riemannsolver_advdiff1d(this)
  ! this method uses an linear upwind solver for the
  ! advective flux and the bassi-rebay method for the 
  ! diffusive fluxes
  implicit none
  class(advdiff1d), intent(inout) :: this
  ! Local
  integer :: iel
  integer :: ivar
  integer :: iside
  real(prec) :: fin, fout, dfavg, un

    call this % solutionGradient % BassiRebaySides(this % gpuaccel)

    do iel = 1, this % mesh % nelem
      do iside = 1, 2

        ! set the normal velocity
        if( iside == 1 )then
          un = -this % u
        else
          un = this % u
        endif

        do ivar = 1, this % solution % nvar

          fin = this % solution % boundary % hostdata(ivar,iside,iel) ! interior solution
          fout = this % solution % extboundary % hostdata(ivar,iside,iel) ! exterior solution
          dfavg = this % solutionGradient % avgboundary % hostdata(ivar,iside,iel) ! average solution gradient (with direction taken into account)

          this % flux % boundary % hostdata(ivar,iside,iel) = 0.5_prec*( un*(fin + fout) + abs(un)*(fin - fout) ) -& ! advective flux
                                                              this % nu*dfavg ! diffusive flux

        enddo

      enddo
    enddo

  end subroutine riemannsolver_advdiff1d

  
end module advection_diffusion_model_1d

! ---- Main program below ---- !

program advection_diffusion_1d

  use self_data
  use advection_diffusion_model_1d

  implicit none
  character(SELF_INTEGRATOR_LENGTH), parameter :: integrator = 'rk3'
  integer, parameter :: nvar = 1
  integer, parameter :: nelem = 50
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  real(prec), parameter :: u = 1.0_prec ! velocity
  real(prec), parameter :: nu = 0.01_prec ! diffusivity
  real(prec), parameter :: dt = 1.0_prec*10.0_prec**(-4) ! time-step size
  real(prec), parameter :: endtime = 1.0_prec
  real(prec), parameter :: iointerval = 0.1_prec
  integer, parameter :: stepsperio = 1000
  type(advdiff1d) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh1D),target :: mesh
  type(Geometry1D),target :: geometry
  type(MPILayer),target :: decomp

  ! Create a mesh using the built-in
  ! uniform mesh generator.
  ! The domain is set to x in [0,1]
  ! We use `nelem` elements
  call mesh % UniformBlockMesh(nGeo=1,&
                               nElem=nelem,&
                               x=(/0.0_prec,1.0_prec/))

  ! We create a domain decomposition.
  call decomp % Init(enableMPI=.false.)
  call decomp % GenerateDecomposition(nelem,nelem+1)

  ! Create an interpolant
  call interp % Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

  ! Generate geometry (metric terms) from the mesh elements
  call geometry % Init(interp,mesh % nElem)
  call geometry % GenerateFromMesh(mesh)

  ! Initialize the model
  call modelobj % Init(nvar,mesh,geometry,decomp)

  ! Set the velocity
  modelobj % u = u
  !Set the diffusivity
  modelobj % nu = nu

  ! Set the initial condition
  call modelobj % solution % SetEquation( 1, 'f = 1.0')
  call modelobj % solution % SetInteriorFromEquation( geometry, 0.0_prec ) 

  print*, "min, max (interior)", &
    minval(modelobj % solution % interior % hostdata), &
    maxval(modelobj % solution % interior % hostdata)

  ! Set the model's time integration method
  call modelobj % SetTimeIntegrator( integrator )

  ! forward step the model to `endtime` using a time step
  ! of `dt` and outputing model data every `iointerval`
  call modelobj % ForwardStep(endtime,dt,iointerval)

  print*, "min, max (interior)", &
  minval(modelobj % solution % interior % hostdata), &
  maxval(modelobj % solution % interior % hostdata)
  

  ! Clean up
  call modelobj % free()
  call decomp % free()
  call mesh % free()
  call geometry % free()
  call interp % free()
  

end program advection_diffusion_1d
