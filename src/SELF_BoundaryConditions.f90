! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2024 Fluid Numerics LLC
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
!    the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_BoundaryConditions

  use SELF_SupportRoutines
  use SELF_Metadata

  implicit none
  integer,parameter :: SELF_BCNAME_LENGTH = 32

  enum,bind(c)
    enumerator :: SELF_BC_STATE_CONTEXT = 0
    enumerator :: SELF_BC_GRADIENT_CONTEXT = 1
  endenum

  type SELF_BoundaryCondition
    procedure(SELF_BCFunction),pointer :: bcFunction => null() ! For state BCs
    procedure(SELF_BCgFunction),pointer :: bcgFunction => null() ! For gradient BCs
    integer :: bcid
    character(SELF_BCNAME_LENGTH) :: bcname
    integer :: context_enum
    type(SELF_BoundaryCondition),pointer :: next => null()
    type(SELF_BoundaryCondition),pointer :: prev => null()
  endtype SELF_BoundaryCondition

  type SELF_BoundaryConditionList
    type(SELF_BoundaryCondition),pointer :: current => null()
    type(SELF_BoundaryCondition),pointer :: head => null()
    type(SELF_BoundaryCondition),pointer :: tail => null()
    integer :: nbc

  contains
    procedure,public :: init => Init_BCList
    procedure,public :: free => Free_BCList
    procedure,private :: MoveNext
    procedure,private :: rewind
    procedure,public :: GetBCForID
    generic,public :: RegisterBoundaryCondition => RegisterBCFunction,RegisterBCGFunction
    procedure,private :: RegisterBCFunction
    procedure,private :: RegisterBCGFunction

  endtype SELF_BoundaryConditionList

  interface
    pure function SELF_BCFunction(this,s,dsdx,x,t,nhat,nvar,ndim) result(extstate)
      use SELF_Constants,only:prec
      import SELF_BoundaryCondition
      implicit none
      class(SELF_BoundaryCondition),intent(in) :: this
      integer,intent(in) :: nvar,ndim
      real(prec),intent(in) :: s(1:nvar)
      real(prec),intent(in) :: dsdx(1:nvar,1:ndim)
      real(prec),intent(in) :: x(1:ndim)
      real(prec),intent(in) :: nhat(1:ndim)
      real(prec),intent(in) :: t
      real(prec) :: extstate(1:nvar)
    endfunction SELF_BCFunction
  endinterface

  interface
    pure function SELF_BCGFunction(this,s,dsdx,x,t,nhat,nvar,ndim) result(extstate)
      use SELF_Constants,only:prec
      import SELF_BoundaryCondition
      implicit none
      class(SELF_BoundaryCondition),intent(in) :: this
      integer,intent(in) :: nvar,ndim
      real(prec),intent(in) :: s(1:nvar)
      real(prec),intent(in) :: dsdx(1:nvar,1:ndim)
      real(prec),intent(in) :: x(1:ndim)
      real(prec),intent(in) :: nhat(1:ndim)
      real(prec),intent(in) :: t
      real(prec) :: extstate(1:nvar,1:ndim)
    endfunction SELF_BCGFunction
  endinterface

contains

! //////////////////////////////////////////// !
!  Boundary Condition Methods
! ////////////////////////////////////////////// !

  subroutine Init_BCList(list)
    class(SELF_BoundaryConditionList),intent(inout) :: list
    list%head => null()
    list%tail => null()
    list%current => null()
    list%nbc = 0
  endsubroutine Init_BCList

  subroutine Free_BCList(list)
    class(SELF_BoundaryConditionList),intent(inout) :: list
    type(SELF_BoundaryCondition),pointer :: node,next_node

    node => list%head
    do while(associated(node))
      next_node => node%next
      nullify(node%bcFunction)
      deallocate(node)
      node => next_node
    enddo

    call Init_BCList(list)
  endsubroutine Free_BCList

  subroutine MoveNext(list)
    class(SELF_BoundaryConditionList),intent(inout) :: list
    if(associated(list%current%next)) then
      list%current => list%current%next
    else
      nullify(list%current)
    endif
  endsubroutine MoveNext

  subroutine rewind(list)
    class(SELF_BoundaryConditionList),intent(inout) :: list
    list%current => list%head
  endsubroutine rewind

  function GetBCForID(list,bcid) result(node)
    !! This function returns the node associated with the given bcid
    !! and context. If the bcid is not found, a null pointer is returned.
    class(SELF_BoundaryConditionList),intent(in) :: list
    integer,intent(in) :: bcid
    type(SELF_BoundaryCondition),pointer :: node

    node => list%head

    do while(associated(node))
      if(node%bcid == bcid) then
        return
      endif
      node => node%next
    enddo
    ! If we reach this point, the bcid was not found
    ! and we return a null pointer
    node => null()

  endfunction GetBCForID

  subroutine RegisterBCFunction(list,bcid,bcname,bcfunc)
    !! Register a boundary condition function
    !! with the given bcid and bcname. If the bcid
    !! is already registered, the function is updated.
    !! The function is expected to be a pointer to a
    !! SELF_BCFunction type.
    class(SELF_BoundaryConditionList),intent(inout) :: list
    integer,intent(in) :: bcid
    character(*),intent(in) :: bcname
    procedure(SELF_BCFunction),pointer,intent(in) :: bcfunc
    ! Local
    type(SELF_BoundaryCondition),pointer :: bc

    ! Check if bcid is registered
    bc => list%GetBCForID(bcid)
    if(associated(bc)) then
      ! If the bcid is already registered, we do not register it again
      print*,"Boundary condition with ID ",bcid," is already registered."
      print*,"Assigning new function to existing BC"
      bc%bcFunction => bcfunc
    else
      allocate(bc)
      bc%bcid = bcid
      bc%bcname = trim(bcname)
      bc%bcFunction => bcfunc
      nullify(bc%next)
      nullify(bc%prev)

      ! Insert at the tail
      if(.not. associated(list%head)) then
        ! First entry
        list%head => bc
        list%tail => bc
      else
        ! Append to tail
        bc%prev => list%tail
        list%tail%next => bc
        list%tail => bc
      endif

      list%nbc = list%nbc+1
      list%current => bc

    endif

  endsubroutine RegisterBCFunction

  subroutine RegisterBCGFunction(list,bcid,bcname,bcgfunc)
    !! Register a boundary condition function
    !! with the given bcid and bcname. If the bcid
    !! is already registered, the function is updated.
    !! The function is expected to be a pointer to a
    !! SELF_BCGFunction type.
    class(SELF_BoundaryConditionList),intent(inout) :: list
    integer,intent(in) :: bcid
    character(*),intent(in) :: bcname
    procedure(SELF_BCGFunction),pointer,intent(in) :: bcgfunc
    ! Local
    type(SELF_BoundaryCondition),pointer :: bc

    ! Check if bcid is registered
    bc => list%GetBCForID(bcid)
    if(associated(bc)) then
      ! If the bcid is already registered, we do not register it again
      print*,"Boundary condition with ID ",bcid," is already registered."
      print*,"Assigning new function to existing BC"
      bc%bcgFunction => bcgfunc
    else
      allocate(bc)
      bc%bcid = bcid
      bc%bcname = trim(bcname)
      bc%bcgFunction => bcgfunc
      nullify(bc%next)
      nullify(bc%prev)

      ! Insert at the tail
      if(.not. associated(list%head)) then
        ! First entry
        list%head => bc
        list%tail => bc
      else
        ! Append to tail
        bc%prev => list%tail
        list%tail%next => bc
        list%tail => bc
      endif

      list%nbc = list%nbc+1
      list%current => bc

    endif

  endsubroutine RegisterBCGFunction

endmodule SELF_BoundaryConditions
