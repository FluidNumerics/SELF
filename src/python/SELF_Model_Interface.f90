
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
  use SELF_DGModel1D
  use SELF_DGModel2D
  use SELF_DGModel3D
  use self_LinearShallowWater2D

  ! External
  use iso_fortran_env
  use iso_c_binding

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

  character(c_char,len=500),private :: model_configuration_file

  ! Interfaces
  public :: Initialize,ForwardStep,WritePickupFile,UpdateParameters,GetSolution,GetPrecision,GetVariableName,Finalize
  private :: GetBCFlagForChar,Init2DWorkspace,InitLinearShallowWater2D

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

  subroutine Initialize(config_file) bind(C,name="Initialize")
    implicit none
    character(kind=c_char,len=*),intent(in) :: config_file
    ! local
    character(len=MODEL_NAME_LENGTH) :: modelname

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

  subroutine Finalize() bind(C,name="Finalize")
    implicit none

    call config%Free()
    call selfModel%Free()
    selfModel => null()

    ! Free the interpolant
    call interp%Free()

    ! Free the mesh

    ! Free the geometry

  endsubroutine Finalize

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

    print*,"Using Mesh file : "//trim(meshfile)
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

  function WritePickupFile(case_directory) result(pickupFile) bind(C,name="WritePickupFile")
    implicit none
    character(kind=c_char,len=*) :: case_directory
    character(LEN=self_FileNameLength) :: pickupFile
    ! Local
    character(13) :: timeStampString

    write(timeStampString,'(I13.13)') this%ioIterate
    pickupFile = case_directory//'/solution.'//timeStampString//'.h5'
    call selfModel%WriteModel(pickupfile)

  endfunction WritePickupFile

  subroutine UpdateParameters() bind(c,name="UpdateParameters")
    implicit none
    character(len=self_IntegratorTypeCharLength) :: timeIntegrator

    call config%Free()
    call config%Init(model_configuration_file)

    ! Set the time integrator
    call config%Get("time_options.integrator",timeIntegrator)
    call selfModel%SetTimeIntegrator(trim(timeIntegrator))

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

  function ForwardStep(dt,updateInterval) result(err) bind(c,name="ForwardStep")
    implicit none
    real(c_prec) :: dt
    real(c_prec) :: updateInterval
    integer(c_int) :: err
    ! Local
    real(prec) :: targetTime

    selfModel%dt = real(dt,prec)
    targetTime = selfModel%t+selfModel%dt*real(updateInterval,prec)
    call selfModel%timeIntegrator(targetTime)
    selfModel%t = targetTime

    ! To do, check solution validity
    err = 0

  endfunction ForwardStep

  function GetPrecision() result(precision) bind(c,name="GetPrecision")
    integer(c_int) :: precision

    precision = prec

  endfunction GetPrecision

  subroutine GetSolution(solution,solshape,ndim) bind(C,name="GetSolution")
    type(c_ptr),intent(out) :: solution ! Pointer to data
    integer(c_int),intent(out) :: solshape(5) ! Shape array (max 4D)
    integer(c_int),intent(out) :: ndim ! Number of dimensions (3, 4, 5)

    select type(selfModel)

    class is(DGModel1D)
      solshape(1:3) = shape(selfModel%solution%interior)
      solshape(4:5) = 0
      ndim = 3
      call selfModel%solution%UpdateHost()
      solution = c_loc(selfModel%solution%interior)

    class is(DGModel2D)
      solshape(1:4) = shape(selfModel%solution%interior)
      solshape(5) = 0
      ndim = 4
      call selfModel%solution%UpdateHost()
      solution = c_loc(selfModel%solution%interior)

    class is(DGModel3D)
      solshape(1:5) = shape(selfModel%solution%interior)
      ndim = 5
      call selfModel%solution%UpdateHost()
      solution = c_loc(selfModel%solution%interior)

    endselect

  endsubroutine GetSolution

  function GetVariableName() result(name) bind(c,name="GetVariableName")
    character(kind=c_char,len=*) :: name

    select type(selfModel)

    class is(SELF_DGModel1D)
      name = selfModel%solution%meta(ivar)%name

    class is(SELF_DGModel2D)
      name = selfModel%solution%meta(ivar)%name

    class is(SELF_DGModel3D)
      name = selfModel%solution%meta(ivar)%name

    endselect

  endfunction GetVariableName

endmodule SELF_Model_Interface
