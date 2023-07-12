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

  TYPE(CompressibleIdealGas2D),TARGET,PRIVATE :: selfCompressibleIdealGas2D
  TYPE(LinearShallowWater),TARGET,PRIVATE :: selfLinearShallowWater2D
  TYPE(Mesh1D),TARGET,PRIVATE :: selfMesh1D
  TYPE(Mesh2D),TARGET,PRIVATE :: selfMesh2D
  TYPE(Geometry1D),TARGET,PRIVATE :: selfGeometry1D
  TYPE(SEMQuad),TARGET,PRIVATE :: selfGeometry2D

  INTEGER,PARAMETER :: MODEL_NAME_LENGTH = 50

CONTAINS

  SUBROUTINE InitializeSimulation()
#undef __FUNC__
#define __FUNC__ "InitializeSimulation"
    IMPLICIT NONE

    ! Local
    CHARACTER(LEN=MODEL_NAME_LENGTH) :: modelname

    CALL config % Init()

    CALL config % Get("model_name",modelname)
    INFO("Setting model to "//TRIM(modelname))

    ! Select the model
    SELECT CASE (TRIM(modelname))
    CASE ("CompressibleIdealGas2D")

      CALL Init2DWorkspace()
      CALL InitCompressibleIdealGas2D()

    CASE ("LinearShallowWater2D")

      CALL Init2DWorkspace()
      !CALL InitLinearShallowWater2D()

    CASE DEFAULT
    END SELECT

    CALL InitGPUAcceleration()

  END SUBROUTINE InitializeSimulation

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

  END SUBROUTINE Init2DWorkspace

  SUBROUTINE InitCompressibleIdealGas2D()
#undef __FUNC__
#define __FUNC__ "Init"
    IMPLICIT NONE

    INFO("Model set to CompressibleIdealGas2D")
    CALL selfCompressibleIdealGas2D % Init(5, &
                                           selfMesh2D,selfGeometry2D,decomp)

    selfModel => selfCompressibleIdealGas2D

  END SUBROUTINE InitCompressibleIdealGas2D

END MODULE SELF_Main

PROGRAM SELF

  ! Core
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
  USE SELF_CompressibleIdealGas2D
  USE SELF_LinearShallowWater
  USE SELF_Main

  IMPLICIT NONE
  ! Public

  CALL InitializeSimulation()
  CALL selfModel % PrintType()
  CALL selfModel % SetInitialConditions(config)

END PROGRAM SELF
