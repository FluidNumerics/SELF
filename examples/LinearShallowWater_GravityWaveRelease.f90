PROGRAM LinearShallowWater_GravityWaveRelease

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Mesh
USE SELF_Geometry
USE SELF_LinearShallowWater

  INTEGER, PARAMETER :: N = 7 ! Polynomial degree of solution
  INTEGER, PARAMETER :: quadrature = GAUSS ! Quadrature
  INTEGER, PARAMETER :: M = 15 ! Number of points in the uniform plotting mesh
  INTEGER, PARAMETER :: nvar = 3 ! The number prognostic variables
  REAL(prec), PARAMETER :: dt = 0.001_prec ! Time step size
  REAL(prec), PARAMETER :: tn = 0.5_prec ! File time 


  TYPE(Lagrange),TARGET :: interp
  TYPE(Mesh2D),TARGET :: mesh
  TYPE(SEMQuad),TARGET :: geometry
  TYPE(LinearShallowWater),TARGET :: semModel
  TYPE(MPILayer),TARGET :: decomp
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: initialCondition(1:nvar)

    ! Initialize a domain decomposition
    ! Here MPI is disabled, since scaling is currently
    ! atrocious with the uniform block mesh
    CALL decomp % Init(enableMPI=.FALSE.)

    ! Create an interpolant
    CALL interp % Init(N,quadrature,M,UNIFORM)

    ! Create a uniform block mesh
    !CALL mesh % UniformBlockMesh(N,(/nXe,nYe/),(/0.0_prec,Lx,0.0_prec,Ly/))
    CALL mesh % Read_ISMv2('./mesh/Circle.mesh')

    ! Generate a decomposition
     CALL decomp % GenerateDecomposition(mesh)

    ! Generate geometry (metric terms) from the mesh elements
    CALL geometry % GenerateFromMesh(mesh,interp,meshQuadrature=GAUSS_LOBATTO)

    ! Initialize the semModel
    CALL semModel % Init(nvar,mesh,geometry,decomp)

    ! Enable GPU Acceleration (if a GPU is found) !
    CALL semModel % EnableGPUAccel()

    ! Set the initial condition
    initialCondition = (/"u = 0.0                                      ", &
                         "v = 0.0                                      ", &
                         "n = exp( -( (x-0.5-t)^2 + (y-0.5-t)^2 )/0.1 )"/)
    CALL semModel % SetSolution( initialCondition )

    ! Write the initial condition to file
    CALL semModel % Write()
    CALL semModel % WriteTecplot()

    !! Set the time integrator (euler, rk3, rk4)
    !CALL semModel % SetTimeIntegrator("Euler")

    !! Set your time step
    !semModel % dt = dt

    !! Forward step the semModel and do the file io
    !CALL semModel % ForwardStep( tn = tn )

    !! Manually write the last semModel state
    !CALL semModel % Write()
    !CALL semModel % WriteTecplot()

    ! Clean up
    CALL semModel % Free()
    CALL decomp % Free()
    CALL geometry % Free()
    CALL mesh % Free()
    CALL interp % Free()

END PROGRAM LinearShallowWater_GravityWaveRelease
