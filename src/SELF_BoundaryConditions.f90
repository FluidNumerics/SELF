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
  use SELF_Model

  implicit none
  integer,parameter :: SELF_BCNAME_LENGTH = 32

  type BoundaryCondition
    procedure(SELF_bcMethod),pointer :: bcMethod => null() !
    integer :: bcid
    character(SELF_BCNAME_LENGTH) :: bcname
    integer :: nBoundaries ! Number of boundaries this BC applies to
    integer,allocatable :: elements(:) ! List of elements this BC applies to
    integer,allocatable :: sides(:) ! List of local sides this BC applies to
    type(BoundaryCondition),pointer :: next => null()
    type(BoundaryCondition),pointer :: prev => null()
  endtype BoundaryCondition

  type BoundaryConditionList
    type(BoundaryCondition),pointer :: current => null()
    type(BoundaryCondition),pointer :: head => null()
    type(BoundaryCondition),pointer :: tail => null()
    integer :: nbc

  contains
    procedure,public :: init => Init_BCList
    procedure,public :: free => Free_BCList
    procedure,private :: MoveNext
    procedure,private :: rewind
    procedure,public :: GetBCForID
    generic,public :: RegisterBoundaryCondition => RegisterbcMethod
    procedure,private :: RegisterbcMethod

  endtype BoundaryConditionList

  interface
    subroutine SELF_bcMethod(this,mymodel)
      use SELF_Constants,only:prec
      use SELF_Model,only:Model
      import BoundaryCondition
      implicit none
      class(BoundaryCondition),intent(in) :: this
      class(Model),intent(inout) :: mymodel
    endsubroutine SELF_bcMethod
  endinterface

contains

! //////////////////////////////////////////// !
!  Boundary Condition Methods
! ////////////////////////////////////////////// !

  subroutine Init_BCList(list)
    class(BoundaryConditionList),intent(inout) :: list
    list%head => null()
    list%tail => null()
    list%current => null()
    list%nbc = 0
  endsubroutine Init_BCList

  subroutine Free_BCList(list)
    class(BoundaryConditionList),intent(inout) :: list
    type(SELF_BoundaryCondition),pointer :: node,next_node

    node => list%head
    do while(associated(node))
      next_node => node%next
      nullify(node%bcMethod)
      if allocated(node%elements) deallocate(node%elements)
      if allocated(node%sides) deallocate(node%sides)
      deallocate(node)
      node => next_node
    enddo

    call Init_BCList(list)
  endsubroutine Free_BCList

  subroutine MoveNext(list)
    class(BoundaryConditionList),intent(inout) :: list
    if(associated(list%current%next)) then
      list%current => list%current%next
    else
      nullify(list%current)
    endif
  endsubroutine MoveNext

  subroutine rewind(list)
    class(BoundaryConditionList),intent(inout) :: list
    list%current => list%head
  endsubroutine rewind

  function GetBCForID(list,bcid) result(node)
    !! This function returns the node associated with the given bcid
    !! and context. If the bcid is not found, a null pointer is returned.
    class(BoundaryConditionList),intent(in) :: list
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

  subroutine RegisterbcMethod(list,bcid,bcname,bcfunc,nboundaries)
    !! Register a boundary condition function
    !! with the given bcid and bcname. If the bcid
    !! is already registered, the function is updated.
    !! The function is expected to be a pointer to a
    !! SELF_bcMethod type.
    class(BoundaryConditionList),intent(inout) :: list
    integer,intent(in) :: bcid
    character(*),intent(in) :: bcname
    procedure(SELF_bcMethod),pointer,intent(in) :: bcfunc
    integer,intent(in) :: nboundaries
    ! Local
    type(SELF_BoundaryCondition),pointer :: bc

    ! Check if bcid is registered
    bc => list%GetBCForID(bcid)
    if(associated(bc)) then
      ! If the bcid is already registered, we do not register it again
      print*,"Boundary condition with ID ",bcid," is already registered."
      print*,"Assigning new function to existing BC"
      bc%bcMethod => bcfunc
    else
      allocate(bc)
      bc%bcid = bcid
      bc%bcname = trim(bcname)
      bc%bcMethod => bcfunc
      allocate(bc%elements(1:nboundaries))
      allocate(bc%sides(1:nboundaries))
      bc%nBoundaries = nboundaries
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

  endsubroutine RegisterbcMethod

  subroutine RegisterBCGFunction(list,bcid,bcname,bcgfunc)
    !! Register a boundary condition function
    !! with the given bcid and bcname. If the bcid
    !! is already registered, the function is updated.
    !! The function is expected to be a pointer to a
    !! SELF_BCGFunction type.
    class(BoundaryConditionList),intent(inout) :: list
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
