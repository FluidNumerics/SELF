PROGRAM GradientBR_Test

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Mesh
USE SELF_Geometry
USE SELF_CLI
USE SELF_Model2D

  IMPLICIT NONE

  INTEGER, PARAMETER :: nvar = 1 ! The number prognostic variables

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
  TYPE(Model2D),TARGET :: semModel
  TYPE(Model2D),TARGET :: exact
  TYPE(MPILayer),TARGET :: decomp
  TYPE(CLI) :: args
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: initialCondition(1:nvar)
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: topography
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: coriolis
  CHARACTER(LEN=255) :: SELF_PREFIX
  CHARACTER(LEN=500) :: meshfile
  REAL(prec) :: absDiff, maxXerr, maxYerr
  INTEGER :: iEl, i, j 


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
    CALL args % Get_CLI('--mesh',meshfile)

    IF( TRIM(meshfile) == '')THEN
      meshfile = TRIM(SELF_PREFIX)//"/etc/mesh/GeophysicalBlock2DMedium/Block2D_mesh.h5"
    ENDIF

    ! Initialize a domain decomposition
    ! Here MPI is disabled, since scaling is currently
    ! atrocious with the uniform block mesh
    CALL decomp % Init(enableMPI=mpiRequested)

    ! Create an interpolant
    CALL interp % Init(N,quadrature,M,UNIFORM)

    ! Read the mesh file in
    CALL mesh % Read_HOPr(TRIM(meshfile),decomp)

    ! Generate geometry (metric terms) from the mesh elements
    CALL geometry % Init(interp,mesh % nElem)
    CALL geometry % GenerateFromMesh(mesh)
    
    ! Reset the boundary condition to reflecting
    CALL mesh % ResetBoundaryConditionType(SELF_BC_PRESCRIBED)

    ! Initialize the semModel
    CALL semModel % Init(nvar,mesh,geometry,decomp)
    CALL exact % Init(nvar,mesh,geometry,decomp)

    IF( gpuRequested )THEN
      CALL semModel % EnableGPUAccel()
      ! Update the device for the whole model
      ! This ensures that the mesh, geometry, and default state match on the GPU
      CALL semModel % UpdateDevice()
    ENDIF
    
    ! Set gravity acceleration and fluid depth
    ! Set the initial condition
    initialCondition = (/"s = sin(2.0*pi*x/1000000.0)*sin(2.0*pi*y/1000000.0)"/)
    CALL semModel % SetSolution( initialCondition )

    ! Set the exact solution
    CALL exact % solutionGradient % SetEquation( 1, 1, &
                     "s = (2.0*pi/1000000.0)*cos(2.0*pi*x/1000000.0)*sin(2.0*pi*y/1000000.0)" )

    CALL exact % solutionGradient % SetEquation( 2, 1, &
                     "s = (2.0*pi/1000000.0)*sin(2.0*pi*x/1000000.0)*cos(2.0*pi*y/1000000.0)" )

    CALL exact % solutionGradient % SetInteriorFromEquation( geometry, 0.0_prec )


    ! Estimate gradient using Bassi-Rebay method
    semModel % solution % extBoundary % hostData = 0.0_prec ! Set the boundary conditions
    CALL semModel % solution % SideExchange(mesh, decomp, semModel % gpuAccel) ! SideExchange (update interior edges)

    CALL semModel % solution % GradientBR(geometry, &
                                          semModel % solutionGradient, &
                                          semModel % gpuAccel)

    CALL semModel % WriteModel('solution.h5')
    CALL semModel % WriteTecplot('solution.tec')
    CALL exact % WriteModel('exact.h5')
    CALL exact % WriteTecplot('exact.tec')

    ! Clean up
    CALL semModel % Free()
    CALL decomp % Free()
    CALL geometry % Free()
    CALL mesh % Free()
    CALL interp % Free()
    CALL args % Free()
    CALL decomp % Finalize()

END PROGRAM GradientBR_Test
