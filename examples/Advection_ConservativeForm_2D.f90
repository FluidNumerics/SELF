PROGRAM Advection_ConservativeForm_2D

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Mesh
USE SELF_Geometry
USE SELF_Advection2D

  INTEGER, PARAMETER :: N = 7 ! Polynomial degree of solution
  INTEGER, PARAMETER :: quadrature = GAUSS_LOBATTO ! Quadrature
  INTEGER, PARAMETER :: M = 15 ! Number of points in the uniform plotting mesh
  INTEGER, PARAMETER :: nXe = 10 ! Number of elements in the x-direction
  INTEGER, PARAMETER :: nYe = 10 ! Number of elements in the y-direction
  INTEGER, PARAMETER :: nvar = 1 ! The number of tracer fields
  REAL(prec), PARAMETER :: Lx = 1.0_prec ! Length of the domain in the x-direction 
  REAL(prec), PARAMETER :: Ly = 1.0_prec ! Length of the domain in the y-direction 
  REAL(prec), PARAMETER :: dt = 0.001_prec ! Time step size
  REAL(prec), PARAMETER :: tn = 0.5_prec ! File time 


  TYPE(Lagrange),TARGET :: interp
  TYPE(Mesh2D),TARGET :: mesh
  TYPE(SEMQuad),TARGET :: geometry
  TYPE(Advection2D),TARGET :: semModel
  TYPE(MPILayer),TARGET :: decomp
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: initialCondition(1:nvar)
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: velocityField(1:2)

    ! Initialize a domain decomposition
    ! Here MPI is disabled, since scaling is currently
    ! atrocious with the uniform block mesh
    CALL decomp % Init(enableMPI=.FALSE.)

    ! Create an interpolant
    CALL interp % Init(N,quadrature,M,UNIFORM)

    ! Create a uniform block mesh
    CALL mesh % UniformBlockMesh(N,(/nXe,nYe/),(/0.0_prec,Lx,0.0_prec,Ly/))

    ! Generate a decomposition
     CALL decomp % GenerateDecomposition(mesh)

    ! Generate geometry (metric terms) from the mesh elements
    CALL geometry % GenerateFromMesh(mesh,interp,meshQuadrature=GAUSS_LOBATTO)

    ! Initialize the semModel
    CALL semModel % Init(nvar,mesh,geometry,decomp)

    ! Enable GPU Acceleration (if a GPU is found) !
    CALL semModel % EnableGPUAccel()

    ! Set the velocity field
    velocityField = (/"vx=1.0","vy=1.0"/)
    ! CALL semModel % SetVelocityField( velocityField )

    ! Set the initial condition
    initialCondition = (/"s = exp( -( (x-0.5-t)^2 + (y-0.5-t)^2 )/0.1 )"/)
    CALL semModel % SetSolution( initialCondition )

    ! Set the boundary condition
    ! CALL semModel % SetPrescribedBoundaryCondition( initialCondition )

    ! Set the time integrator (euler, rk3, rk4)
    CALL semModel % SetTimeIntegrator("Euler")

    ! Set your time step
    semModel % dt = dt

    ! Set the file IO frequency
    ! semModel % 

    ! TO DO :  Set the formulation type (conservative or splitform)

    ! Forward step the semModel and do the file io
    CALL semModel % ForwardStep( tn = tn )


    ! Manually write the last semModel state
    CALL semModel % Write()

    ! Clean up
    CALL semModel % Free()
    CALL decomp % Free()
    CALL geometry % Free()
    CALL mesh % Free()
    CALL interp % Free()

END PROGRAM Advection_ConservativeForm_2D
