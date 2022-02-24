PROGRAM LinearShallowWater_GravityWaveRelease

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Mesh
USE SELF_Geometry
USE SELF_LinearShallowWater

  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 7 ! Polynomial degree of solution
  INTEGER, PARAMETER :: quadrature = GAUSS ! Quadrature
  INTEGER, PARAMETER :: M = 15 ! Number of points in the uniform plotting mesh
  INTEGER, PARAMETER :: nvar = 3 ! The number prognostic variables
  REAL(prec), PARAMETER :: dt = 0.001_prec ! Time step size
  REAL(prec), PARAMETER :: tn = 2.0_prec ! Total simulation time
  REAL(prec), PARAMETER :: ioInterval = 0.025_prec ! File IO interval

  TYPE(Lagrange),TARGET :: interp
  TYPE(Mesh2D),TARGET :: mesh
  TYPE(SEMQuad),TARGET :: geometry
  TYPE(LinearShallowWater),TARGET :: semModel
  TYPE(MPILayer),TARGET :: decomp
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: initialCondition(1:nvar)
  CHARACTER(LEN=255) :: SELF_PREFIX

    ! Initialize a domain decomposition
    ! Here MPI is disabled, since scaling is currently
    ! atrocious with the uniform block mesh
    CALL decomp % Init(enableMPI=.FALSE.)

    ! Create an interpolant
    CALL interp % Init(N,quadrature,M,UNIFORM)

    ! Create a uniform block mesh
    CALL get_environment_variable("SELF_PREFIX", SELF_PREFIX)
    CALL mesh % Read_HOPr(TRIM(SELF_PREFIX)//"/etc/mesh/Circle/Circle_mesh.h5")

    ! Generate a decomposition
     CALL decomp % GenerateDecomposition(mesh)

    ! Generate geometry (metric terms) from the mesh elements
    CALL geometry % Init(interp,mesh % nElem)
    CALL geometry % GenerateFromMesh(mesh)

    ! Reset the boundary condition to reflecting
    CALL mesh % ResetBoundaryConditionType(SELF_BC_NONORMALFLOW)

    ! Initialize the semModel
    CALL semModel % Init(nvar,mesh,geometry,decomp)

    ! Enable GPU Acceleration (if a GPU is found) !
    CALL semModel % EnableGPUAccel()

    ! Set the initial condition
    initialCondition = (/"u = 0.0                                         ", &
                         "v = 0.0                                         ", &
                         "n = 0.01*exp( -( ((x+0.5)^2 + (y+0.5)^2 )/0.01 )"/)
    CALL semModel % SetSolution( initialCondition )

    ! Write the initial condition to file
    CALL semModel % WriteModel()
    CALL semModel % WriteTecplot()

    ! Set the time integrator (euler, rk3, rk4)
    CALL semModel % SetTimeIntegrator("Euler")

    ! Set your time step
    semModel % dt = dt

    !! Forward step the semModel and do the file io
    CALL semModel % ForwardStep( tn = tn, ioInterval = ioInterval )

    !! Manually write the last semModel state
    CALL semModel % WriteModel('solution.pickup.h5')

    ! Clean up
    CALL semModel % Free()
    CALL decomp % Free()
    CALL geometry % Free()
    CALL mesh % Free()
    CALL interp % Free()

END PROGRAM LinearShallowWater_GravityWaveRelease
