module SELF_LinearShallowWater2D_Interface
! Core
  use SELF_Constants
  use SELF_SupportRoutines
  use SELF_Mesh
  use SELF_Geometry_2D
  use SELF_JSON_Config

! Models
  use self_LinearShallowWater2D

! External
  use iso_fortran_env
  use iso_c_binding

contains

  subroutine Init_LinearShallowWater2D(modelObj,geometry,mesh)
    implicit none
    type(LinearShallowWater2D),intent(inout) :: modelObj
    type(SEMQuad),intent(in) :: geometry
    type(Mesh2D),intent(in) :: mesh

    print*,"Model set to Linear Shallow Water (2D)"

    call modelObj%Init(mesh,geometry)
    modelObj%prescribed_bcs_enabled = .false. ! Disables prescribed boundary condition block for gpu accelerated implementations
    modelObj%tecplot_enabled = .false. ! Disables tecplot output

  endsubroutine Init_LinearShallowWater2D

  subroutine UpdateParameters_LinearShallowWater2D(modelObj,config)
    implicit none
    type(LinearShallowWater2D),intent(inout) :: modelObj
    type(SELFConfig),intent(inout) :: config

    call config%Get("linear-shallow-water-2d.environment.g", &
                    modelObj%g)

    call config%Get("linear-shallow-water-2d.environment.H", &
                    modelObj%H)

    call config%Get("linear-shallow-water-2d.environment.Cd", &
                    modelObj%Cd)

    call config%Get("linear-shallow-water-2d.environment.f0", &
                    modelObj%f0)

    call config%Get("linear-shallow-water-2d.environment.beta", &
                    modelObj%beta)

  endsubroutine UpdateParameters_LinearShallowWater2D

endmodule SELF_LinearShallowWater2D_Interface
