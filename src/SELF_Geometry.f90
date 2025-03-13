module SELF_Geometry

  use SELF_Constants
  use SELF_Lagrange

  implicit none

  type,abstract :: SEMGeometry
    integer :: nElem
  contains

    procedure(SELF_InitGeometry),deferred :: Init
    procedure(SELF_FreeGeometry),deferred:: Free

  endtype SEMGeometry

  interface
    subroutine SELF_InitGeometry(this,interp,nElem)
      import SEMGeometry
      import Lagrange
      implicit none
      class(SEMGeometry),intent(out) :: this
      type(Lagrange),pointer,intent(in) :: interp
      integer,intent(in) :: nElem
    endsubroutine SELF_InitGeometry
  endinterface

  interface
    subroutine SELF_FreeGeometry(this)
      import SEMGeometry
      implicit none
      class(SEMGeometry),intent(inout) :: this
    endsubroutine SELF_FreeGeometry
  endinterface

endmodule SELF_Geometry
