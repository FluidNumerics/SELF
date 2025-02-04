
module SELF_Model_Interface

  ! Core
  use SELF_Constants
  use SELF_SupportRoutines
  use SELF_Mesh
  use SELF_Geometry_1D
  use SELF_Geometry_2D
  use SELF_Geometry_3D
  use SELF_JSON_Config

  ! Models
  use SELF_Model
  !use SELF_DGModel1D
  use SELF_DGModel2D
  ! use SELF_DGModel3D
  use self_LinearShallowWater2D

  ! External
  use iso_fortran_env
!   use HDF5
!   use FEQParse

  implicit none

  type(SELFConfig) :: config
  type(Lagrange),target,private :: interp
  !type(MPILayer),target :: decomp

  class(Model),pointer,private :: selfModel
  class(SEMMesh),pointer,private :: selfMesh
  !class(SEMGeometry),pointer,private :: selfGeometry

  ! Mesh
  ! type(Mesh1D),target,private :: selfMesh1D
  type(Mesh2D),target,private :: selfMesh2D

  ! Geometry
  ! type(Geometry1D),target,private :: selfGeometry1D
  type(SEMQuad),target,private :: selfGeometry2D

  ! Models
  ! type(Burgers1D),target,private :: selfBurgers1D
  ! type(CompressibleIdealGas2D),target,private :: selfCompressibleIdealGas2D
  type(LinearShallowWater2D),target,private :: selfLinearShallowWater2D

  integer,parameter,private :: MODEL_NAME_LENGTH = 50

  character(LEN=500),private :: model_configuration_file

  ! Interfaces
  public :: Initialize,ForwardStep,WritePickupFile!, !GetSolution, Finalize
  private :: GetBCFlagForChar,Init2DWorkspace,UpdateParameters,InitLinearShallowWater2D

contains

  function GetBCFlagForChar(charFlag) result(intFlag)
    !! This method is used to return the integer flag from a char for boundary conditions
    !!
    implicit none
    character(*),intent(in) :: charFlag
    integer :: intFlag

    select case(UpperCase(trim(charFlag)))

    case("PRESCRIBED")
      intFlag = SELF_BC_PRESCRIBED

    case("RADIATION")
      intFlag = SELF_BC_RADIATION

    case("NO_NORMAL_FLOW")
      intFlag = SELF_BC_NONORMALFLOW

    case DEFAULT
      intFlag = 0

    endselect

  endfunction GetBCFlagForChar
  function GetQFlagForChar(charFlag) result(intFlag)
    !! This method is used to return the integer flag from a char for boundary conditions
    !!
    implicit none
    character(*),intent(in) :: charFlag
    integer :: intFlag

    select case(UpperCase(trim(charFlag)))

    case("GAUSS")
      intFlag = GAUSS

    case("GAUSS-LOBATTO")
      intFlag = GAUSS_LOBATTO

    case DEFAULT
      intFlag = 0

    endselect

  endfunction GetQFlagForChar

  subroutine Initialize(config_file)
    implicit none
    character(LEN=*),intent(in) :: config_file
    ! local
    character(LEN=MODEL_NAME_LENGTH) :: modelname

    call config%Init(config_file)
    model_configuration_file = config_file

    ! Select the model
    select case(trim(modelname))

    case("linear-shallow-water-2d")

      call Init2DWorkspace()
      call InitLinearShallowWater2D()

    case default
    endselect

  endsubroutine Initialize

  subroutine Init2DWorkspace()
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

    call config%Get("geometry.control_degree",controlDegree)
    call config%Get("geometry.target_degree",targetDegree)
    call config%Get("geometry.control_quadrature",qChar)
    controlQuadrature = GetQFlagForChar(trim(qChar))
    call config%Get("geometry.target_quadrature",qChar)
    targetQuadrature = GetQFlagForChar(trim(qChar))
    call config%Get("geometry.mesh_file",meshfile)
    call config%Get("geometry.uniform_boundary_condition",uniformBoundaryCondition)
    bcFlag = GetBCFlagForChar(uniformBoundaryCondition)

    print*, "Using Mesh file : "//trim(meshfile)
    ! Read in mesh file and set the public mesh pointer to selfMesh2D
    call selfMesh2D%Read_HOPr(trim(meshfile))
    call selfMesh2D%ResetBoundaryConditionType(bcFlag)

    selfMesh => selfMesh2D

    ! Create an interpolant
    call interp%Init(controlDegree, &
                     controlQuadrature, &
                     targetDegree, &
                     targetQuadrature)

    ! Generate geometry (metric terms) from the mesh elements
    call selfGeometry2D%Init(interp,selfMesh2D%nElem)
    call selfGeometry2D%GenerateFromMesh(selfMesh2D)

!    selfGeometry => selfGeometry2D

  endsubroutine Init2DWorkspace

  subroutine InitLinearShallowWater2D()
    implicit none

    print*,"Model set to Linear Shallow Water (2D)"

    call selfLinearShallowWater2D%Init(selfMesh2D,selfGeometry2D)
    selfLinearShallowWater2D%prescribed_bcs_enabled = .false. ! Disables prescribed boundary condition block for gpu accelerated implementations
    selfLinearShallowWater2D%tecplot_enabled = .false. ! Disables tecplot output

    call UpdateParameters()

    selfModel => selfLinearShallowWater2D

  endsubroutine InitLinearShallowWater2D

  subroutine WritePickupFile()
    implicit none

    select type(selfModel)

    type is(LinearShallowWater2D)
      call selfModel%WriteModel()

    endselect
  endsubroutine WritePickupFile

  subroutine UpdateParameters()
    implicit none

    call config%Free()
    call config%Init(model_configuration_file)
    select type(selfModel)

    type is(LinearShallowWater2D)

      call config%Get("linear-shallow-water-2d.environment.g", &
                      selfLinearShallowWater2D%g)

      call config%Get("linear-shallow-water-2d.environment.H", &
                      selfLinearShallowWater2D%H)

      call config%Get("linear-shallow-water-2d.environment.Cd", &
                      selfLinearShallowWater2D%Cd)

      call config%Get("linear-shallow-water-2d.environment.f0", &
                      selfLinearShallowWater2D%f0)

      call config%Get("linear-shallow-water-2d.environment.beta", &
                      selfLinearShallowWater2D%beta)

    endselect

  endsubroutine UpdateParameters

  subroutine ForwardStep()
    implicit none

    real(prec) :: dt
    real(prec) :: targetTime
    integer :: updateInterval

    ! Reload config
    call UpdateParameters()

    call config%Get("time_options.update_interval",updateInterval)
    call config%Get("time_options.dt",dt)

    selfModel%dt = dt
    targetTime = selfModel%t+dt*real(updateInterval,prec)
    call selfModel%timeIntegrator(targetTime)
    selfModel%t = targetTime

  endsubroutine ForwardStep

endmodule SELF_Model_Interface

! program SELF

!   use SELF_Main

!   implicit none
!   ! Public

!   call InitializeSELF()

!   ! Show us which model we're running
!   call selfModel % PrintType()

!   ! Set the model initial conditions
!   call selfModel % SetInitialConditions(config)

!   ! Do the initial file IO
!   call FileIO()
!   ! Main time integration loop
!   call MainLoop()
!   ! Do the last file IO
!   call FileIO()

!   ! Clean up [TO DO]
!   ! CALL selfModel % Free()
!   ! CALL selfGeometry % Free()
!   ! CALL selfMesh % Free()
!   ! CALL interp % Free()

! endmodule SELF_Model_Interface
