PROGRAM ShallowWater_QuietFluid

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Mesh
USE SELF_Geometry
USE SELF_ShallowWater
USE SELF_CLI

  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 7 ! Polynomial degree of solution
  INTEGER, PARAMETER :: quadrature = GAUSS ! Quadrature
  INTEGER, PARAMETER :: M = 15 ! Number of points in the uniform plotting mesh
  INTEGER, PARAMETER :: nvar = 3 ! The number prognostic variables
  REAL(prec), PARAMETER :: dt = 0.001_prec ! Time step size
  REAL(prec), PARAMETER :: tn = 0.001_prec ! Total simulation time
  REAL(prec), PARAMETER :: ioInterval = 0.001_prec ! File IO interval
  REAL(prec), PARAMETER :: tolerance=1000.0_prec*epsilon(1.0_prec) ! Error tolerance

  REAL(prec) :: referenceEntropy
  REAL(prec) :: solutionMax(1:3)
  TYPE(Lagrange),TARGET :: interp
  TYPE(Mesh2D),TARGET :: mesh
  TYPE(SEMQuad),TARGET :: geometry
  TYPE(ShallowWater),TARGET :: semModel
  TYPE(MPILayer),TARGET :: decomp
  TYPE(CLI) :: args
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: initialCondition(1:nvar)
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: topography
  CHARACTER(LEN=255) :: SELF_PREFIX

    CALL get_environment_variable("SELF_PREFIX", SELF_PREFIX)
    CALL args % Init( TRIM(SELF_PREFIX)//"/etc/cli/default.json")
    CALL args % LoadFromCLI()

    ! Initialize a domain decomposition
    ! Here MPI is disabled, since scaling is currently
    ! atrocious with the uniform block mesh
    CALL decomp % Init(enableMPI=.FALSE.)

    ! Create an interpolant
    CALL interp % Init(N,quadrature,M,UNIFORM)

    ! Create a uniform block mesh
    CALL get_environment_variable("SELF_PREFIX", SELF_PREFIX)
    CALL mesh % Read_HOPr(TRIM(SELF_PREFIX)//"/etc/mesh/Block2D/Block2D_mesh.h5")

    ! Reset the boundary condition to reflecting
    !CALL mesh % ResetBoundaryConditionType(SELF_BC_NONORMALFLOW)
    CALL mesh % ResetBoundaryConditionType(SELF_BC_RADIATION)

    ! Generate a decomposition
     CALL decomp % GenerateDecomposition(mesh)

    ! Generate geometry (metric terms) from the mesh elements
    CALL geometry % Init(interp,mesh % nElem)
    CALL geometry % GenerateFromMesh(mesh)

    ! Initialize the semModel
    CALL semModel % Init(nvar,mesh,geometry,decomp)

    ! Enable GPU Acceleration (if a GPU is found) !
    CALL semModel % EnableGPUAccel()
 
    topography = "h = 1.0"
    CALL semModel % SetTopography(topography)

    ! Set the initial condition
    initialCondition = (/"u = 0.0", &
                         "v = 0.0", &
                         "H = 1.0"/)
    CALL semModel % SetSolution( initialCondition )
    referenceEntropy = semModel % entropy

    ! Write the initial condition to file
    CALL semModel % WriteModel()
    CALL semModel % WriteTecplot()

    ! Set the time integrator (euler, rk3)
    CALL semModel % SetTimeIntegrator("rk3")

    ! Set your time step
    semModel % dt = dt

    !! Forward step the semModel and do the file io
    CALL semModel % ForwardStep( tn = tn, ioInterval = ioInterval )

    !! Manually write the last semModel state
    CALL semModel % WriteModel('solution.pickup.h5')

    ! Error checking !
    IF( semModel % entropy /= semModel % entropy )THEN
      PRINT*, "Model entropy is not a number"
      STOP 2
    ENDIF

    IF( semModel % entropy >= HUGE(1.0_prec) )THEN
      PRINT*, "Model entropy is infinite."
      STOP 1
    ENDIF

    IF( semModel % entropy > referenceEntropy )THEN
      PRINT*, "Warning : final entropy greater than initial entropy"
      ! Currently do nothing in this situation, since
      ! conservative solvers in mapped geometries may
      ! not be entropy conservative.
      ! However, throwing this warning will bring some
      ! visibility
    ENDIF

    ! Check the solution !
    solutionMax = semModel % solution % AbsMaxInterior() 
    IF( solutionMax(1) > tolerance .OR. solutionMax(2) > tolerance)THEN
      PRINT*, "Non-zero velocity field detected for quiescent fluid."
      PRINT*, solutionMax
      STOP 1
    ENDIF

    ! Clean up
    CALL semModel % Free()
    CALL decomp % Free()
    CALL geometry % Free()
    CALL mesh % Free()
    CALL interp % Free()
    CALL args % Free()

END PROGRAM ShallowWater_QuietFluid
