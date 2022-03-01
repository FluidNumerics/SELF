PROGRAM Advection_ConservativeForm_2D

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Mesh
USE SELF_Geometry
USE SELF_Advection2D

  INTEGER, PARAMETER :: N = 7 ! Polynomial degree of solution
  INTEGER, PARAMETER :: quadrature = GAUSS ! Quadrature
  INTEGER, PARAMETER :: M = 15 ! Number of points in the uniform plotting mesh
  INTEGER, PARAMETER :: nXe = 5 ! Number of elements in the x-direction
  INTEGER, PARAMETER :: nYe = 5 ! Number of elements in the y-direction
  INTEGER, PARAMETER :: nvar = 1 ! The number of tracer fields
  REAL(prec), PARAMETER :: Lx = 1.0_prec ! Length of the domain in the x-direction 
  REAL(prec), PARAMETER :: Ly = 1.0_prec ! Length of the domain in the y-direction 
  REAL(prec), PARAMETER :: dt = 0.001_prec ! Time step size
  REAL(prec), PARAMETER :: tn = 0.25_prec ! File time 
  REAL(prec), PARAMETER :: ioInterval = 0.25_prec ! File IO interval


  REAL(prec) :: referenceEntropy
  TYPE(Lagrange),TARGET :: interp
  TYPE(Mesh2D),TARGET :: mesh
  TYPE(SEMQuad),TARGET :: geometry
  TYPE(Advection2D),TARGET :: semModel
  TYPE(MPILayer),TARGET :: decomp
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: initialCondition(1:nvar)
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: velocityField(1:2)
  CHARACTER(LEN=255) :: SELF_PREFIX

    ! Initialize a domain decomposition
    ! Here MPI is disabled, since scaling is currently
    ! atrocious with the uniform block mesh
    CALL decomp % Init(enableMPI=.FALSE.)

    ! Create an interpolant
    CALL interp % Init(N,quadrature,M,UNIFORM)

    ! Create a uniform block mesh
    CALL get_environment_variable("SELF_PREFIX", SELF_PREFIX)
    CALL mesh % Read_HOPr(TRIM(SELF_PREFIX)//"/etc/mesh/Block2D/Block2D_mesh.h5")

    ! Generate a decomposition
    CALL decomp % GenerateDecomposition(mesh)

    ! Generate geometry (metric terms) from the mesh elements
    CALL geometry % Init(interp,mesh % nElem)
    CALL geometry % GenerateFromMesh(mesh)

    ! Initialize the semModel
    CALL semModel % Init(nvar,mesh,geometry,decomp)

    ! Set the solution name
    CALL semModel % solution % SetName(1,'Tracer')

    ! Enable GPU Acceleration (if a GPU is found) !
    CALL semModel % EnableGPUAccel()

    ! Set the velocity field
    velocityField = (/"vx=1.0","vy=1.0"/)
    CALL semModel % SetVelocityField( velocityField )

    ! Set the initial condition
    initialCondition = (/"s = exp( -( (x-0.5-t)^2 + (y-0.5-t)^2 )/0.01 )"/)
    CALL semModel % SetSolution( initialCondition )
    referenceEntropy = semModel % entropy

    ! Write the initial condition to file
    CALL semModel % WriteModel()
    CALL semModel % WriteTecplot()

    ! Set the time integrator (euler, rk3, rk4)
    CALL semModel % SetTimeIntegrator("Euler")

   ! ! Set your time step
    semModel % dt = dt

   ! ! Forward step the semModel and do the file io
    CALL semModel % ForwardStep( tn = tn, ioInterval = ioInterval )

   ! ! Manually write the last semModel state
    CALL semModel % WriteModel()
    CALL semModel % WriteTecplot()

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

END PROGRAM Advection_ConservativeForm_2D
