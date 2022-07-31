PROGRAM LinearShallowWater_GravityWaveRelease

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Mesh
USE SELF_Geometry
USE SELF_CLI
USE SELF_LinearShallowWater

  IMPLICIT NONE

  INTEGER, PARAMETER :: nvar = 3 ! The number prognostic variables
  REAL(prec), PARAMETER :: g = 10.0 ! Acceleration of gravity (m/s^2)
  REAL(prec), PARAMETER :: H = 1000.0 ! Fluid depth (m)

  REAL(prec) :: dt
  REAL(prec) :: ioInterval
  REAL(prec) :: tn
  INTEGER :: N ! Control Degree
  INTEGER :: M ! Target degree
  INTEGER :: quadrature
  CHARACTER(LEN=self_QuadratureTypeCharLength) :: qChar
  LOGICAL :: mpiRequested
  LOGICAL :: gpuRequested

  REAL(prec) :: referenceEntropy
  TYPE(Lagrange),TARGET :: interp
  TYPE(Mesh2D),TARGET :: mesh
  TYPE(SEMQuad),TARGET :: geometry
  TYPE(LinearShallowWater),TARGET :: semModel
  TYPE(MPILayer),TARGET :: decomp
  TYPE(CLI) :: args
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: initialCondition(1:nvar)
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: topography
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: coriolis
  CHARACTER(LEN=255) :: SELF_PREFIX


    CALL get_environment_variable("SELF_PREFIX", SELF_PREFIX)
    CALL args % Init( TRIM(SELF_PREFIX)//"/etc/cli/default.json")
    CALL args % LoadFromCLI()

    CALL args % Get_CLI('--output-interval',ioInterval)
    CALL args % Get_CLI('--end-time',tn)
    CALL args % Get_CLI('--time-step',dt)
    CALL args % Get_CLI('--mpi',mpiRequested)
    CALL args % Get_CLI('--gpu',gpuRequested)
    CALL args % Get_CLI('--control-degree',N)
    CALL args % Get_CLI('--control-quadrature',qChar)
    quadrature = GetIntForChar(qChar)
    CALL args % Get_CLI('--target-degree',M)


    ! Initialize a domain decomposition
    ! Here MPI is disabled, since scaling is currently
    ! atrocious with the uniform block mesh
    CALL decomp % Init(enableMPI=mpiRequested)

    ! Create an interpolant
    CALL interp % Init(N,quadrature,M,UNIFORM)

    ! Create a uniform block mesh
    CALL mesh % Read_HOPr(TRIM(SELF_PREFIX)//"/etc/mesh/GeophysicalBlock2DMedium/Block2D_mesh.h5")

    ! Generate a decomposition
     CALL decomp % GenerateDecomposition(mesh)

    ! Generate geometry (metric terms) from the mesh elements
    CALL geometry % Init(interp,mesh % nElem)
    CALL geometry % GenerateFromMesh(mesh)
    
    ! Reset the boundary condition to reflecting
    CALL mesh % ResetBoundaryConditionType(SELF_BC_RADIATION)

    ! Initialize the semModel
    CALL semModel % Init(nvar,mesh,geometry,decomp)
    
    ! Set gravity acceleration and fluid depth
    semModel % g = g
    CALL semModel % SetBathymetry( H )

    ! Enable GPU Acceleration (if a GPU is found) !
    IF( gpuRequested )THEN
      CALL semModel % EnableGPUAccel()
    ENDIF

    ! Set the initial condition
    initialCondition = (/"u = 0.0                                             ", &
                         "v = 0.0                                             ", &
                         "n = 10^(-2)*exp( -( (x^2 + y^2 )/(2.0*10.0^(10)) ) )"/)
    CALL semModel % SetSolution( initialCondition )
    referenceEntropy = semModel % entropy

    ! Set the coriolis parameter
    coriolis = "f = 10^(-4) + 2.0*(10^(-11))*y"
    CALL semModel % SetCoriolis( coriolis )
    
    ! Get the geostrophic velocity from the free surface height field
    CALL semModel % DiagnoseGeostrophicVelocity()

    ! Write the initial condition to file
    CALL semModel % WriteModel()
    CALL semModel % WriteTecplot()

    ! Set the time integrator (euler, rk3, rk4)
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

    ! Clean up
    CALL semModel % Free()
    CALL decomp % Free()
    CALL geometry % Free()
    CALL mesh % Free()
    CALL interp % Free()
    CALL args % Free()

END PROGRAM LinearShallowWater_GravityWaveRelease
