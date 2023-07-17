MODULE SELF_Main

  ! Core
  USE SELF_Constants
  USE SELF_SupportRoutines
  USE SELF_Mesh
  USE SELF_Geometry
  USE SELF_MappedData
  USE SELF_Config

  ! Models
  USE SELF_Model
  USE SELF_Model1D
  USE SELF_Model2D
  USE SELF_ECModel2D
  USE SELF_Burgers1D
  USE SELF_CompressibleIdealGas2D
  USE SELF_LinearShallowWater

  ! External
  USE ISO_FORTRAN_ENV
  USE HDF5
  USE FEQParse

  IMPLICIT NONE

#include "../src/SELF_Macros.h"

  TYPE(SELFConfig) :: config
  TYPE(Lagrange),TARGET :: interp
  TYPE(MPILayer),TARGET :: decomp

  CLASS(Model),POINTER :: selfModel
  CLASS(SEMMesh),POINTER :: selfMesh
  CLASS(SEMGeometry),POINTER :: selfGeometry

  TYPE(Burgers1D),TARGET,PRIVATE :: selfBurgers1D
  TYPE(CompressibleIdealGas2D),TARGET,PRIVATE :: selfCompressibleIdealGas2D
  TYPE(LinearShallowWater),TARGET,PRIVATE :: selfLinearShallowWater2D
  TYPE(Mesh1D),TARGET,PRIVATE :: selfMesh1D
  TYPE(Mesh2D),TARGET,PRIVATE :: selfMesh2D
  TYPE(Geometry1D),TARGET,PRIVATE :: selfGeometry1D
  TYPE(SEMQuad),TARGET,PRIVATE :: selfGeometry2D

  !REAL(prec),PUBLIC :: referenceEntropy

  INTEGER,PARAMETER :: MODEL_NAME_LENGTH = 50

CONTAINS

  SUBROUTINE InitializeSELF()
#undef __FUNC__
#define __FUNC__ "InitializeSELF"
    IMPLICIT NONE

    ! Local
    CHARACTER(LEN=MODEL_NAME_LENGTH) :: modelname

    ! TO DO : Check for --test flag
    ! IF( test ) testing=TRUE; return
    ! ELSE

    CALL config % Init()

    CALL config % Get("model_name",modelname)
    INFO("Setting model to "//TRIM(modelname))

    ! Select the model
    SELECT CASE (TRIM(modelname))
    CASE ("brg1d")

      CALL Init1DWorkspace()
      CALL InitBurgers1D()
      !CALL InitLinearShallowWater2D()

    CASE ("cns2d")

      CALL Init2DWorkspace()
      CALL InitCompressibleIdealGas2D()

    CASE ("lsw2d")

      CALL Init2DWorkspace()
      !CALL InitLinearShallowWater2D()

    CASE DEFAULT
    END SELECT

    CALL InitGPUAcceleration()

  END SUBROUTINE InitializeSELF

  SUBROUTINE InitGPUAcceleration()
#undef __FUNC__
#define __FUNC__ "InitGPUAcceleration"
    IMPLICIT NONE
    ! Local
    LOGICAL :: gpu

    CALL config % Get("deployment_options.gpu_accel",gpu)

    IF (gpu) THEN
      ! Enable GPU Acceleration (if a GPU is found)
      ! Update the device for the whole model
      ! This ensures that the mesh, geometry, and default state match on the GPU
      INFO("GPU acceleration enabled")

      CALL selfModel % EnableGPUAccel()

      IF (selfModel % gpuAccel) THEN

        SELECT TYPE (selfModel)

        TYPE IS (Burgers1D)
          !CALL selfModel % UpdateDevice()
          WARNING("GPU acceleration not implemented yet")

        TYPE IS (CompressibleIdealGas2D)
          CALL selfModel % UpdateDevice()

        TYPE IS (LinearShallowWater)
          CALL selfModel % UpdateDevice()

        END SELECT

        INFO("Model updated on device")
      END IF

    ELSE
      INFO("Running CPU-Only mode")
    END IF

  END SUBROUTINE InitGPUAcceleration

  SUBROUTINE Init1DWorkspace()
#undef __FUNC__
#define __FUNC__ "Init1DWorkspace"
    IMPLICIT NONE
    ! Local
    LOGICAL :: mpiRequested
    CHARACTER(LEN=self_QuadratureTypeCharLength) :: qChar
    CHARACTER(LEN=MODEL_NAME_LENGTH) :: meshfile
    INTEGER :: controlQuadrature
    INTEGER :: controlDegree
    INTEGER :: targetDegree
    INTEGER :: targetQuadrature
    REAL(prec) :: x0,xN
    INTEGER :: nX

    CALL config % Get("deployment_options.mpi",mpiRequested)
    IF (mpiRequested) THEN
      INFO("MPI domain decomposition not implemented. Disabling")
      mpiRequested = .FALSE.
    ELSE
      INFO("MPI domain decomposition disabled")
    END IF

    CALL config % Get("geometry.control_degree",controlDegree)
    CALL config % Get("geometry.target_degree",targetDegree)
    CALL config % Get("geometry.control_quadrature",qChar)
    controlQuadrature = GetIntForChar(TRIM(qChar))
    CALL config % Get("geometry.target_quadrature",qChar)
    targetQuadrature = GetIntForChar(TRIM(qChar))
    CALL config % Get("geometry.x0",x0)
    CALL config % Get("geometry.xN",xN)
    CALL config % Get("geometry.nX",nX)

    ! Initialize a domain decomposition
    CALL decomp % Init(enableMPI=mpiRequested)

    ! Create an interpolant
    CALL interp % Init(controlDegree, &
                       controlQuadrature, &
                       targetDegree, &
                       targetQuadrature)

    ! Read in mesh file and set the public mesh pointer to selfMesh2D
    CALL selfMesh1D % UniformBlockMesh(1,nX, (/x0,xN/))
    selfMesh => selfMesh1D

    ! Generate geometry (metric terms) from the mesh elements
    CALL selfGeometry1D % Init(interp,selfMesh1D % nElem)
    CALL selfGeometry1D % GenerateFromMesh(selfMesh1D)
    selfGeometry => selfGeometry1D

    ! Reset the boundary condition to prescribed

  END SUBROUTINE Init1DWorkspace

  SUBROUTINE Init2DWorkspace()
#undef __FUNC__
#define __FUNC__ "Init2DWorkspace"
    IMPLICIT NONE
    ! Local
    LOGICAL :: mpiRequested
    CHARACTER(LEN=self_QuadratureTypeCharLength) :: qChar
    CHARACTER(LEN=MODEL_NAME_LENGTH) :: meshfile
    INTEGER :: controlQuadrature
    INTEGER :: controlDegree
    INTEGER :: targetDegree
    INTEGER :: targetQuadrature

    CALL config % Get("deployment_options.mpi",mpiRequested)
    IF (mpiRequested) THEN
      INFO("MPI domain decomposition enabled")
    ELSE
      INFO("MPI domain decomposition disabled")
    END IF

    CALL config % Get("geometry.control_degree",controlDegree)
    CALL config % Get("geometry.target_degree",targetDegree)
    CALL config % Get("geometry.control_quadrature",qChar)
    controlQuadrature = GetIntForChar(TRIM(qChar))
    CALL config % Get("geometry.target_quadrature",qChar)
    targetQuadrature = GetIntForChar(TRIM(qChar))
    CALL config % Get("geometry.mesh_file",meshfile)

    INFO("Mesh file : "//TRIM(meshfile))

    ! Initialize a domain decomposition
    CALL decomp % Init(enableMPI=mpiRequested)

    ! Create an interpolant
    CALL interp % Init(controlDegree, &
                       controlQuadrature, &
                       targetDegree, &
                       targetQuadrature)

    ! Read in mesh file and set the public mesh pointer to selfMesh2D
    CALL selfMesh2D % Read_HOPr(TRIM(meshfile),decomp)
    selfMesh => selfMesh2D

    ! Generate geometry (metric terms) from the mesh elements
    CALL selfGeometry2D % Init(interp,selfMesh2D % nElem)
    CALL selfGeometry2D % GenerateFromMesh(selfMesh2D)
    selfGeometry => selfGeometry2D

    ! Reset the boundary condition to prescribed

  END SUBROUTINE Init2DWorkspace

  SUBROUTINE InitBurgers1D()
#undef __FUNC__
#define __FUNC__ "Init"
    IMPLICIT NONE

    INFO("Model set to Burgers1D")

    ! For now - default all boundaries to prescribed
    ! CALL selfMesh2D % ResetBoundaryConditionType(SELF_BC_PRESCRIBED)

    CALL selfBurgers1D % Init(1, &
                              selfMesh1D,selfGeometry1D,decomp)

    selfModel => selfBurgers1D

  END SUBROUTINE InitBurgers1D

  SUBROUTINE InitCompressibleIdealGas2D()
#undef __FUNC__
#define __FUNC__ "Init"
    IMPLICIT NONE

    INFO("Model set to CompressibleIdealGas2D")

    ! For now - default all boundaries to prescribed
    CALL selfMesh2D % ResetBoundaryConditionType(SELF_BC_PRESCRIBED)

    CALL selfCompressibleIdealGas2D % Init(5, &
                                           selfMesh2D,selfGeometry2D,decomp)

    selfModel => selfCompressibleIdealGas2D

  END SUBROUTINE InitCompressibleIdealGas2D

  SUBROUTINE FileIO()
#undef __FUNC__
#define __FUNC__ "FileIO"
    IMPLICIT NONE

    SELECT TYPE (selfModel)

    TYPE IS (Burgers1D)
    ! Write the initial condition to file
    CALL selfModel % WriteModel()
    CALL selfModel % WriteTecplot()

    TYPE IS (CompressibleIdealGas2D)
      ! Write the initial condition to file
      CALL selfModel % WriteModel()
      CALL selfModel % WriteTecplot()
    TYPE IS (LinearShallowWater)
      CALL selfModel % WriteModel()
      CALL selfModel % WriteTecplot()

    END SELECT
  END SUBROUTINE FileIO

  SUBROUTINE MainLoop()
#undef __FUNC__
#define __FUNC__ "MainLoop"
    IMPLICIT NONE
    REAL(prec) :: startTime
    REAL(prec) :: duration
    REAL(prec) :: endTime
    REAL(prec) :: ioInterval

    INFO("Starting main loop")
    CALL config % Get("time_options.start_time",startTime)
    CALL config % Get("time_options.duration",duration)
    endTime = startTime + duration
    CALL config % Get("time_options.io_interval",ioInterval)

    CALL selfModel % ForwardStep(tn=endTime,ioInterval=ioInterval)
!    referenceEntropy = selfModel % entropy

  !  !! Forward step the selfModel and do the file io
  !   SELECT TYPE (selfModel)

  !   TYPE IS (CompressibleIdealGas2D)
  !     CALL selfModel % ForwardStep(tn=endTime,ioInterval=ioInterval)
  !   TYPE IS (LinearShallowWater)
  !     CALL selfModel % ForwardStep(tn=endTime,ioInterval=ioInterval)
  !   END SELECT

    ! !! Manually write the last selfModel state
    !  CALL selfModel % WriteModel('solution.pickup.h5')

    !  ! Error checking !
    !  IF (selfModel % entropy /= selfModel % entropy) THEN
    !    PRINT *, "Model entropy is not a number"
    !    STOP 2
    !  END IF

    !  IF (selfModel % entropy >= HUGE(1.0_PREC)) THEN
    !    PRINT *, "Model entropy is infinite."
    !    STOP 1
    !  END IF

  END SUBROUTINE MainLoop

END MODULE SELF_Main

PROGRAM SELF

  ! Core
  USE SELF_SupportRoutines
  USE SELF_Mesh
  USE SELF_Geometry
  USE SELF_MappedData
  USE SELF_Config

  ! Models
  ! USE SELF_Model
  ! USE SELF_Model1D
  ! USE SELF_Model2D
  ! USE SELF_ECModel2D
  ! USE SELF_CompressibleIdealGas2D
  ! USE SELF_LinearShallowWater
  USE SELF_Main

  IMPLICIT NONE
  ! Public

  CALL InitializeSELF()

  ! Show us which model we're running
  CALL selfModel % PrintType()

  ! Set the model initial conditions
  CALL selfModel % SetInitialConditions(config)

  ! Do the initial file IO
  CALL FileIO()
  ! Main time integration loop
  CALL MainLoop()
  ! Do the last file IO
  CALL FileIO()

  ! Clean up [TO DO]
  ! CALL selfModel % Free()
  ! CALL selfGeometry % Free()
  ! CALL selfMesh % Free()
  ! CALL interp % Free()

END PROGRAM SELF
