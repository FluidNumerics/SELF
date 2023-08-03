module SELF_Main

  ! Core
  use SELF_Constants
  use SELF_SupportRoutines
  use SELF_Mesh
  use SELF_Geometry
  use SELF_MappedData
  use SELF_Config

  ! Models
  use SELF_Model
  use SELF_Model1D
  use SELF_Model2D
  use SELF_ECModel2D
  use SELF_Burgers1D
  use SELF_CompressibleIdealGas2D
  use SELF_LinearShallowWater

  ! External
  use iso_fortran_env
  use HDF5
  use FEQParse

  implicit none

#include "../src/SELF_Macros.h"

  type(SELFConfig) :: config
  type(Lagrange),target :: interp
  type(MPILayer),target :: decomp

  class(Model),pointer :: selfModel
  class(SEMMesh),pointer :: selfMesh
  class(SEMGeometry),pointer :: selfGeometry

  type(Burgers1D),target,private :: selfBurgers1D
  type(CompressibleIdealGas2D),target,private :: selfCompressibleIdealGas2D
  type(LinearShallowWater),target,private :: selfLinearShallowWater2D
  type(Mesh1D),target,private :: selfMesh1D
  type(Mesh2D),target,private :: selfMesh2D
  type(Geometry1D),target,private :: selfGeometry1D
  type(SEMQuad),target,private :: selfGeometry2D

  !REAL(prec),PUBLIC :: referenceEntropy

  integer,parameter :: MODEL_NAME_LENGTH = 50

contains

  subroutine InitializeSELF()
#undef __FUNC__
#define __FUNC__ "InitializeSELF"
    implicit none

    ! Local
    character(LEN=MODEL_NAME_LENGTH) :: modelname

    ! TO DO : Check for --test flag
    ! IF( test ) testing=TRUE; return
    ! ELSE

    call config % Init()

    call config % Get("model_name",modelname)
    INFO("Setting model to "//trim(modelname))

    ! Select the model
    select case (trim(modelname))
    case ("brg1d")

      call Init1DWorkspace()
      call InitBurgers1D()
      !CALL InitLinearShallowWater2D()

    case ("cns2d")

      call Init2DWorkspace()
      call InitCompressibleIdealGas2D()

    case ("lsw2d")

      call Init2DWorkspace()
      !CALL InitLinearShallowWater2D()

    case DEFAULT
    end select

    call InitGPUAcceleration()

  end subroutine InitializeSELF

  subroutine InitGPUAcceleration()
#undef __FUNC__
#define __FUNC__ "InitGPUAcceleration"
    implicit none
    ! Local
    logical :: gpu

    call config % Get("deployment_options.gpu_accel",gpu)

    if (gpu) then
      ! Enable GPU Acceleration (if a GPU is found)
      ! Update the device for the whole model
      ! This ensures that the mesh, geometry, and default state match on the GPU
      INFO("GPU acceleration enabled")

      call selfModel % EnableGPUAccel()

      if (selfModel % gpuAccel) then

        select type (selfModel)

        type is (Burgers1D)
          !CALL selfModel % UpdateDevice()
          WARNING("GPU acceleration not implemented yet")

        type is (CompressibleIdealGas2D)
          call selfModel % UpdateDevice()

        type is (LinearShallowWater)
          call selfModel % UpdateDevice()

        end select

        INFO("Model updated on device")
      end if

    else
      INFO("Running CPU-Only mode")
    end if

  end subroutine InitGPUAcceleration

  subroutine Init1DWorkspace()
#undef __FUNC__
#define __FUNC__ "Init1DWorkspace"
    implicit none
    ! Local
    logical :: mpiRequested
    character(LEN=self_QuadratureTypeCharLength) :: qChar
    character(LEN=MODEL_NAME_LENGTH) :: meshfile
    integer :: controlQuadrature
    integer :: controlDegree
    integer :: targetDegree
    integer :: targetQuadrature
    real(prec) :: x0,xN
    integer :: nX

    call config % Get("deployment_options.mpi",mpiRequested)
    if (mpiRequested) then
      INFO("MPI domain decomposition not implemented. Disabling")
      mpiRequested = .false.
    else
      INFO("MPI domain decomposition disabled")
    end if

    call config % Get("geometry.control_degree",controlDegree)
    call config % Get("geometry.target_degree",targetDegree)
    call config % Get("geometry.control_quadrature",qChar)
    controlQuadrature = GetIntForChar(trim(qChar))
    call config % Get("geometry.target_quadrature",qChar)
    targetQuadrature = GetIntForChar(trim(qChar))
    call config % Get("geometry.x0",x0)
    call config % Get("geometry.xN",xN)
    call config % Get("geometry.nX",nX)

    ! Initialize a domain decomposition
    call decomp % Init(enableMPI=mpiRequested)

    ! Create an interpolant
    call interp % Init(controlDegree, &
                       controlQuadrature, &
                       targetDegree, &
                       targetQuadrature)

    ! Read in mesh file and set the public mesh pointer to selfMesh2D
    call selfMesh1D % UniformBlockMesh(1,nX, (/x0,xN/))
    selfMesh => selfMesh1D

    ! Generate geometry (metric terms) from the mesh elements
    call selfGeometry1D % Init(interp,selfMesh1D % nElem)
    call selfGeometry1D % GenerateFromMesh(selfMesh1D)
    selfGeometry => selfGeometry1D

    ! Reset the boundary condition to prescribed

  end subroutine Init1DWorkspace

  subroutine Init2DWorkspace()
#undef __FUNC__
#define __FUNC__ "Init2DWorkspace"
    implicit none
    ! Local
    logical :: mpiRequested
    character(LEN=self_QuadratureTypeCharLength) :: qChar
    character(LEN=MODEL_NAME_LENGTH) :: meshfile
    character(LEN=MODEL_NAME_LENGTH) :: uniformBoundaryCondition
    integer :: controlQuadrature
    integer :: controlDegree
    integer :: targetDegree
    integer :: targetQuadrature
    integer :: bcFlag

    call config % Get("deployment_options.mpi",mpiRequested)
    if (mpiRequested) then
      INFO("MPI domain decomposition enabled")
    else
      INFO("MPI domain decomposition disabled")
    end if

    call config % Get("geometry.control_degree",controlDegree)
    call config % Get("geometry.target_degree",targetDegree)
    call config % Get("geometry.control_quadrature",qChar)
    controlQuadrature = GetIntForChar(trim(qChar))
    call config % Get("geometry.target_quadrature",qChar)
    targetQuadrature = GetIntForChar(trim(qChar))
    call config % Get("geometry.mesh_file",meshfile)
    call config % Get("geometry.uniform_boundary_condition",uniformBoundaryCondition)
    bcFlag = GetBCFlagForChar(uniformBoundaryCondition)

    INFO("Mesh file : "//trim(meshfile))

    ! Initialize a domain decomposition
    call decomp % Init(enableMPI=mpiRequested)

    ! Create an interpolant
    call interp % Init(controlDegree, &
                       controlQuadrature, &
                       targetDegree, &
                       targetQuadrature)

    ! Read in mesh file and set the public mesh pointer to selfMesh2D
    call selfMesh2D % Read_HOPr(trim(meshfile),decomp)
    selfMesh => selfMesh2D

    ! Generate geometry (metric terms) from the mesh elements
    call selfGeometry2D % Init(interp,selfMesh2D % nElem)
    call selfGeometry2D % GenerateFromMesh(selfMesh2D)

    ! Reset the boundary conditions, if requested
    if (bcFlag > 0) then
      call selfMesh2D % ResetBoundaryConditionType(bcFlag)
    end if
    selfGeometry => selfGeometry2D

  end subroutine Init2DWorkspace

  subroutine InitBurgers1D()
#undef __FUNC__
#define __FUNC__ "Init"
    implicit none

    INFO("Model set to Burgers1D")
    call selfBurgers1D % Init(1, &
                              selfMesh1D,selfGeometry1D,decomp)

    selfModel => selfBurgers1D

  end subroutine InitBurgers1D

  subroutine InitCompressibleIdealGas2D()
#undef __FUNC__
#define __FUNC__ "Init"
    implicit none

    INFO("Model set to CompressibleIdealGas2D")

    call selfCompressibleIdealGas2D % Init(5, &
                                           selfMesh2D,selfGeometry2D,decomp)

    selfModel => selfCompressibleIdealGas2D

  end subroutine InitCompressibleIdealGas2D

  subroutine FileIO()
#undef __FUNC__
#define __FUNC__ "FileIO"
    implicit none

    select type (selfModel)

    type is (Burgers1D)
      ! Write the initial condition to file
      call selfModel % WriteModel()
      call selfModel % WriteTecplot()

    type is (CompressibleIdealGas2D)
      ! Write the initial condition to file
      call selfModel % WriteModel()
      call selfModel % WriteTecplot()

    type is (LinearShallowWater)
      call selfModel % WriteModel()
      call selfModel % WriteTecplot()

    end select
  end subroutine FileIO

  subroutine MainLoop()
#undef __FUNC__
#define __FUNC__ "MainLoop"
    implicit none
    real(prec) :: startTime
    real(prec) :: duration
    real(prec) :: endTime
    real(prec) :: ioInterval

    INFO("Starting main loop")
    call config % Get("time_options.start_time",startTime)
    call config % Get("time_options.duration",duration)
    endTime = startTime + duration
    call config % Get("time_options.io_interval",ioInterval)

    call selfModel % ForwardStep(tn=endTime,ioInterval=ioInterval)
!    referenceEntropy = selfModel % entropy

    !  !! Forward step the selfModel and do the file io
    !   SELECT TYPE (selfModel)

    !   TYPE is (CompressibleIdealGas2D)
    !     CALL selfModel % ForwardStep(tn=endTime,ioInterval=ioInterval)
    !   TYPE is (LinearShallowWater)
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

  end subroutine MainLoop

end module SELF_Main

program SELF

  use SELF_Main

  implicit none
  ! Public

  call InitializeSELF()

  ! Show us which model we're running
  call selfModel % PrintType()

  ! Set the model initial conditions
  call selfModel % SetInitialConditions(config)

  ! Do the initial file IO
  call FileIO()
  ! Main time integration loop
  call MainLoop()
  ! Do the last file IO
  call FileIO()

  ! Clean up [TO DO]
  ! CALL selfModel % Free()
  ! CALL selfGeometry % Free()
  ! CALL selfMesh % Free()
  ! CALL interp % Free()

end program SELF
