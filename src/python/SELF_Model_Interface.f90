
module SELF_Model_Interface

  ! Core
  use SELF_Constants
  use SELF_SupportRoutines
  use SELF_Metadata

  use SELF_Mesh
  use SELF_Mesh_1D
  use SELF_Mesh_2D
  use SELF_Mesh_3D

  use SELF_Geometry
  use SELF_Geometry_1D
  use SELF_Geometry_2D
  use SELF_Geometry_3D

  use SELF_JSON_Config

  ! Models
  use SELF_Model
  use SELF_DGModel1D
  use SELF_DGModel2D
  use SELF_DGModel3D

  use SELF_Burgers1D
  !use SELF_Burgers1D_Interface

  use SELF_LinearShallowWater2D
  use SELF_LinearShallowWater2D_Interface

  use SELF_LinearEuler2D
  !use SELF_LinearEuler2D_Interface

  use SELF_LinearEuler3D
  !use SELF_LinearEuler3D_Interface

  ! External
  use iso_fortran_env
  use iso_c_binding

  implicit none

  type(SELFConfig) :: config
  type(Lagrange),target,private :: interp
  !type(MPILayer),target :: decomp

  ! ========================================== !
  !             Top level pointers             !
  ! ========================================== !
  class(Model),pointer,private :: selfModel
  class(SEMMesh),pointer,private :: selfMesh
  class(SEMGeometry),pointer,private :: selfGeometry

  ! Mesh
  type(Mesh1D),target,private :: selfMesh1D
  type(Mesh2D),target,private :: selfMesh2D
  type(Mesh3D),target,private :: selfMesh3D

  ! Geometry
  type(Geometry1D),target,private :: selfGeometry1D
  type(SEMQuad),target,private :: selfGeometry2D
  type(SEMHex),target,private :: selfGeometry3D

  ! Models
  type(Burgers1D),target,private :: selfBurgers1D
  type(LinearShallowWater2D),target,private :: selfLinearShallowWater2D
  type(LinearEuler2D),target,private :: selfLinearEuler2D
  type(LinearEuler3D),target,private :: selfLinearEuler3D

  integer,parameter,private :: MODEL_NAME_LENGTH = 50

  character(kind=c_char,len=750),private :: model_configuration_file

  ! Interfaces
  public :: Initialize,ForwardStep,WritePickupFile,UpdateParameters,GetSolution,GetPrecision,GetVariableName,Finalize
  private :: GetBCFlagForChar,Init2DWorkspace
contains

  ! =================================================================
  !                      Public methods
  ! =================================================================

  function Initialize(config_file) result(error) bind(C,name="Initialize")
    implicit none
    character(kind=c_char,len=*),intent(in) :: config_file
    integer(c_int) :: error
    ! local
    character(len=MODEL_NAME_LENGTH) :: modelname

    call config%Init(config_file)
    model_configuration_file = config_file

    call config%Get("model_name",modelname)

    ! Select the model
    select case(trim(modelname))

    case("burgers-1d")

      print*,"Not implemented yet"
      error = -1
      ! call Init1DWorkspace()
      ! call Init_Burgers1D(selfBurgers1D,selfGeometry1D,selfMesh1D)
      ! selfModel => selfBurgers1D
      ! error = 0

    case("linear-shallow-water-2d")

      call Init2DWorkspace()
      call Init_LinearShallowWater2D(selfLinearShallowWater2D,selfGeometry2D,selfMesh2D)
      selfModel => selfLinearShallowWater2D
      error = 0

    case("linear-euler-2d")

      print*,"Not implemented yet"
      error = -1
      ! call Init2DWorkspace()
      ! call Init_LinearEuler2D(selfLinearEuler2D,selfGeometry2D,selfMesh2D)
      ! selfModel => selfLinearEuler2D
      ! error = 0

    case("linear-euler-3d")

      print*,"Not implemented yet"
      error = -1
      ! call Init2DWorkspace()
      ! call Init_LinearEuler2D(selfLinearEuler2D,selfGeometry2D,selfMesh2D)
      ! selfModel => selfLinearEuler2D
      ! error = 0

    case("gfdles-3d")

      print*,"Not implemented yet"
      error = -1
      ! call Init2DWorkspace()
      ! call Init_LinearEuler2D(selfLinearEuler2D,selfGeometry2D,selfMesh2D)
      ! selfModel => selfLinearEuler2D
      ! error = 0

    case default

    endselect

    ! Point the mesh and geometry top level pointers to the appropriate mesh and geometry objects
    select type(selfModel)

    class is(DGModel1D)
      selfMesh => selfMesh1D
      selfGeometry => selfGeometry1D

    class is(DGModel2D)
      selfMesh => selfMesh2D
      selfGeometry => selfGeometry2D

    class is(DGModel3D)
      selfMesh => selfMesh3D
      selfGeometry => selfGeometry3D

    endselect

    call UpdateParameters()

  endfunction Initialize

  subroutine Finalize() bind(C,name="Finalize")
    implicit none

    call config%Free()
    call selfModel%Free()
    call selfMesh%Free()
    call selfGeometry%Free()
    call interp%Free()

    ! Nullify the top level pointers
    selfModel => null()
    selfMesh => null()
    selfGeometry => null()

  endsubroutine Finalize

  subroutine WritePickupFile(case_directory,pickupFile) bind(C,name="WritePickupFile")
    implicit none
    character(kind=c_char,len=*),intent(in) :: case_directory
    character(kind=c_char,len=*),intent(out) :: pickupFile
    ! Local
    character(13) :: timeStampString

    write(timeStampString,'(I13.13)') selfModel%ioIterate
    pickupFile = trim(case_directory)//'/solution.'//timeStampString//'.h5'
    call selfModel%WriteModel(trim(pickupfile))

  endsubroutine WritePickupFile

  subroutine UpdateParameters() bind(c,name="UpdateParameters")
    implicit none
    character(len=self_IntegratorTypeCharLength) :: timeIntegrator

    call config%Free()
    call config%Init(model_configuration_file)

    ! Set the time integrator
    call config%Get("time_options.integrator",timeIntegrator)
    call selfModel%SetTimeIntegrator(trim(timeIntegrator))

    select type(selfModel)

    type is(Burgers1D)
      print*,"Not implemented yet"
      !call UpdateParameters_Burgers1D(selfModel,config)
    type is(LinearShallowWater2D)

      call UpdateParameters_LinearShallowWater2D(selfModel,config)

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

  subroutine GetVariableName(ivar,name) bind(c,name="GetVariableName")
    integer(c_int),intent(in) :: ivar
    character(kind=c_char,len=*),intent(out) :: name

    select type(selfModel)

    class is(DGModel1D)
      name = selfModel%solution%meta(ivar)%name

    class is(DGModel2D)
      name = selfModel%solution%meta(ivar)%name

    class is(DGModel3D)
      name = selfModel%solution%meta(ivar)%name

    endselect

  endsubroutine GetVariableName

  ! =================================================================
  !                      Private methods
  ! =================================================================

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

  subroutine Init2DWorkspace()
    implicit none
    ! Local
    logical :: mpiRequested
    character(len=self_QuadratureTypeCharLength) :: qChar
    character(len=MODEL_NAME_LENGTH) :: meshfile
    character(len=MODEL_NAME_LENGTH) :: uniformBoundaryCondition
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

endmodule SELF_Model_Interface
