PROGRAM Model1D_Write

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Mesh
USE SELF_Geometry
USE SELF_Model1D
USE SELF_CLI

  IMPLICIT NONE

  INTEGER, PARAMETER :: nvar = 1 ! The number prognostic variables
  INTEGER, PARAMETER :: nX = 2

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
  REAL(prec) :: solutionMax(1:1)
  TYPE(Lagrange),TARGET :: interp
  TYPE(Mesh1D),TARGET :: mesh
  TYPE(Geometry1D),TARGET :: geometry
  TYPE(Model1D),TARGET :: semModel
  TYPE(MPILayer),TARGET :: decomp
  TYPE(CLI) :: args
  CHARACTER(LEN=SELF_EQUATION_LENGTH) :: initialCondition(1:nvar)
  CHARACTER(LEN=255) :: SELF_PREFIX
  CHARACTER(LEN=3) :: Nchar

    CALL get_environment_variable("SELF_PREFIX", SELF_PREFIX)
    CALL args % Init( TRIM(SELF_PREFIX)//"/etc/cli/default.json")
    CALL args % LoadFromCLI()

    mpiRequested = .FALSE.
    gpuRequested = .FALSE.
    CALL args % Get_CLI('--control-degree',N)
    CALL args % Get_CLI('--control-quadrature',qChar)
    quadrature = GetIntForChar(qChar)
    M = N

    ! Initialize a domain decomposition
    ! Here MPI is disabled, since scaling is currently
    ! atrocious with the uniform block mesh
    CALL decomp % Init(enableMPI=.FALSE.)

    ! Create an interpolant
    CALL interp % Init(N,quadrature,M,UNIFORM)

    ! Create a uniform block mesh
    CALL mesh % UniformBlockMesh(1,nX,(/-1.0_prec,1.0_prec/))

    ! Generate a decomposition
    CALL decomp % GenerateDecomposition(mesh % nElem, 1)

    ! Generate geometry (metric terms) from the mesh elements
    CALL geometry % Init(interp,mesh % nElem)
    CALL geometry % GenerateFromMesh(mesh)

    ! Initialize the semModel
    CALL semModel % Init(nvar,mesh,geometry,decomp)

    WRITE(UNIT=Nchar,FMT='(I0.3)') N

    ! Write to file
    CALL semModel % WriteModel('sample1D_'//Nchar//'.h5')

    ! Clean up
    CALL semModel % Free()
    CALL decomp % Free()
    CALL geometry % Free()
    CALL mesh % Free()
    CALL interp % Free()
    CALL args % Free()
    CALL decomp % Finalize()

END PROGRAM Model1D_Write
