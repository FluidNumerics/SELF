! SELF_LinkedList.f90
!
! Copyright 2017-2021 Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file SELF_LinkedList.f90
!! Contains the \ref SELF_LinkedList module, and <BR>
!! defines the \ref Record and \ref LinkedList data-structures.

!> \defgroup SELF_LinkedList SELF_LinkedList
!! This module defines the Record and \ref LinkedList data-structure and associated routines.
MODULE SELF_LinkedList

!  A linked-list data-structure with routines for construction, destruction, accession, and
!  modification are provided. This provides a template for a linked list of other data-structures.
!  Here, the data is taken as an integer and we also provide an integer key which can be useful
!  for accessing linked-list data for other data-structures.

  IMPLICIT NONE
!> \addtogroup SELF_LinkedList
!! @{

!> \struct Record
!! A template record entry for a Linked-List that consists of data, a key, and a pointer to the next
!! record
!!
!! This data structure is used in the LinkedList data structure. The data in this template is an
!! integer, though this template can be used to develop of Linked List of other data (e.g. see
!! Timing.f90, Node_Class.f90, QuadElement_Class.f90, HexElement_Class.f90, Edge_Class.f90)
!!
!! <H2> Record </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> listData <td> INTEGER  <td> Data for a single entry in a LinkedList
!!       <tr> <th>  key <td> INTEGER  <td> A key that can be used for conditional access, as in the
!!                                         HashTable_Class
!!       <tr> <th> next <td> Record  <td> Pointer to the next record in a LinkedList
!!    </table>
!!
!> @}
  TYPE Record
    INTEGER                 :: listData
    INTEGER                 :: key
    TYPE(Record),POINTER :: next

  END TYPE Record
!> \addtogroup SELF_LinkedList
!! @{

!> \struct LinkedList
!! A LinkedList of Records
!!
!!
!! <H2> LinkedList </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> head <td> Record  <td> A pointer to the first entry (the head) of the LinkedList
!!       <tr> <th> tail <td> Record  <td> A pointer to the last entry (the tail) of the LinkedList
!!       <tr> <th> current <td> Record  <td> A pointer to the current entry of the LinkedList. This
!!                                           pointer can be used to traverse the list.
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref SELF_LinkedList for more information. The first column lists the "call-name" and the second
!!    column lists the name of routine that is aliased onto the call-name. This is the list of
!!    PUBLIC type-bound procedures for the LinkedList data-structure.
!!
!!    <table>
!!       <tr> <th> Init <td> Init_LinkedList
!!       <tr> <th> Free <td> Free_LinkedList
!!       <tr> <th> GetData <td> GetData_LinkedList
!!       <tr> <th> SetData <td> SetData_LinkedList
!!       <tr> <th> GetKey <td> GetKey_LinkedList
!!       <tr> <th> SetKey <td> SetKey_LinkedList
!!       <tr> <th> ListIsEmpty <td> ListIsEmpty_LinkedList
!!       <tr> <th> AddToList <td> AddToList_LinkedList
!!       <tr> <th> RemoveCurrent <td> RemoveCurrent_LinkedList
!!       <tr> <th> MoveToNext <td> MoveToNext_LinkedList
!!       <tr> <th> MoveToHead <td> MoveToHead_LinkedList
!!       <tr> <th> MovetoTail <td> MoveToTail_LinkedList
!!       <tr> <th> PrintToScreen <td> PrintToScreen_LinkedList
!!    </table>
!> @}
  TYPE LinkedList
    TYPE(Record),POINTER :: head,tail,current

  CONTAINS

    PROCEDURE :: Init => Init_LinkedList
    PROCEDURE :: Free => Free_LinkedList

    PROCEDURE :: GetData => GetCurrentData_LinkedList
    PROCEDURE :: SetData => SetCurrentData_LinkedList
    PROCEDURE :: GetKey => GetCurrentKey_LinkedList
    PROCEDURE :: SetKey => SetCurrentKey_LinkedList

    PROCEDURE :: ListIsEmpty => ListIsEmpty_LinkedList
    PROCEDURE :: AddToList => AddToList_LinkedList
    PROCEDURE :: RemoveCurrent => RemoveCurrent_LinkedList

    PROCEDURE :: MoveToNext => MoveToNext_LinkedList
    PROCEDURE :: MoveToHead => MoveToHead_LinkedList
    PROCEDURE :: MoveToTail => MoveToTail_LinkedList
    PROCEDURE :: PrintToScreen => PrintToScreen_LinkedList

  END TYPE LinkedList

! Setting up some parameters pertaining to this module
  INTEGER,PARAMETER,PRIVATE :: keyInc = 1 ! The default increment in the Record Key
  INTEGER,PARAMETER,PRIVATE :: keyStart = 1 ! The default starting Record key

CONTAINS
!
!==================================================================================================!
!---------------------------- CONSTRUCTOR/DESTRUCTOR ROUTINES -------------------------------------!
!==================================================================================================!
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R Init
!
!> \fn Init_LinkedList
!! Initializes a LinkedList by nullifying the head, tail, and current pointers
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Init(  ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myList <td> LinkedList <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE Init_LinkedList(myList)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList

    myList % head => NULL()
    myList % tail => NULL()
    myList % current => NULL()

  END SUBROUTINE Init_LinkedList
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R Free
!
!> \fn Free_LinkedList
!!  Cycles through the Linked List, deallocates associated memory, and nullifies the LinkedList
!!  pointers.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Free( Inputs/Outputs ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myList <td> LinkedList <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE Free_LinkedList(myList)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList
    ! LOCAL
    TYPE(Record),POINTER :: pNext

    ! Set the current position of the list to the head
    myList % current => myList % head

    ! Scroll through the list until the current position is nullified
    DO WHILE (ASSOCIATED(myList % current))

      pNext => myList % current % next
      DEALLOCATE (myList % current)
      myList % current => pNext

    END DO

  END SUBROUTINE Free_LinkedList
!
!==================================================================================================!
!---------------------------------------- ACCESSORS -----------------------------------------------!
!==================================================================================================!
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R GetCurrentData
!
!> \fn GetCurrentData_LinkedList
!! Gets the current list data and returns the integer value.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!! <B>INTEGER</B>          :: outData <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetCurrentData( outData ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in <th> myList <td> LinkedList <td>
!!   <tr> <td> out <th> outData <td> INTEGER <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE GetCurrentData_LinkedList(myList,outData)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList
    INTEGER             :: outData

    outData = myList % current % listData

  END SUBROUTINE GetCurrentData_LinkedList
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R GetCurrentKey
!
!> \fn GetCurrentKey_LinkedList
!! Gets the current key and returns the integer value.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!! <B>INTEGER</B>          :: outKey <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetCurrentKey( outKey ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in <th> myList <td> LinkedList <td>
!!   <tr> <td> out <th> outKey <td> INTEGER <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE GetCurrentKey_LinkedList(myList,outKey)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList
    INTEGER             :: outKey

    outKey = myList % current % key

  END SUBROUTINE GetCurrentKey_LinkedList
  !
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R SetCurrentData
!
!> \fn SetCurrentData_LinkedList
!! Sets the current list data to the supplied integer value.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!! <B>INTEGER</B>          :: inData <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetCurrentData( inData ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myList <td> LinkedList <td>
!!   <tr> <td> in <th> inData <td> INTEGER <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE SetCurrentData_LinkedList(myList,inData)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList
    INTEGER             :: inData

    myList % current % listData = inData

  END SUBROUTINE SetCurrentData_LinkedList
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R SetCurrentKey
!
!> \fn SetCurrentKey_LinkedList
!! Sets the current list key to the supplied integer value.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!! <B>INTEGER</B>          :: inKey <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetCurrentKey( inKey ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myList <td> LinkedList <td>
!!   <tr> <td> in <th> inKey <td> INTEGER <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE SetCurrentKey_LinkedList(myList,inKey)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList
    INTEGER             :: inKey

    myList % current % key = inKey

  END SUBROUTINE SetCurrentKey_LinkedList
!
!
!==================================================================================================!
!-------------------------------- Linked-List Type Operations -------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R ListIsEmpty
!
!> \fn ListIsEmpty_LinkedList
!! Checks if the LinkedList is empty and returns a logical indicating the status of the list.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!! <B>LOGICAL</B>          :: listStatus <BR>
!!         .... <BR>
!!     listStatus =  this % ListIsEmpty( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in <th> myList <td> LinkedList <td>
!!   <tr> <td> out <th> isempty<td> LOGICAL <td>
!!                      If the list is empty, isempty=.TRUE., <BR>
!!                      If the list has entries, isempty=.FALSE., <BR>
!!  </table>
!!
! ================================================================================================ !
!>@}
  FUNCTION ListIsEmpty_LinkedList(myList) RESULT(isempty)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList
    LOGICAL             :: isempty

    isempty = .NOT. (ASSOCIATED(myList % head))

  END FUNCTION ListIsEmpty_LinkedList
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R AddToList
!
!> \fn AddToList_LinkedList
!! Adds a new entry to the LinkedList with the provided data and key
!!
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!! <B>INTEGER</B>          :: inData, inKey <BR>
!!         .... <BR>
!!     <B>CALL</B> this % AddToList( inData, inKey ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myList <td> LinkedList <td>
!!   <tr> <td> in <th> inData <td> INTEGER <td> Data to assign to the new entry in the list
!!   <tr> <td> in <th> inKey <td> INTEGER <td>  Key to assign to the new entry in the list
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE AddToList_LinkedList(myList,inData,inKey)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList
    INTEGER             :: inData
    INTEGER,OPTIONAL   :: inKey
    ! LOCAL
    TYPE(Record),POINTER :: previous
    INTEGER :: allocationStatus

    ! Check to see if this list is empty
    IF (myList % ListIsEmpty()) THEN

      ALLOCATE (myList % head,STAT=allocationStatus)
      IF (allocationStatus /= 0) THEN
        PRINT *, 'MODULE LinkedListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
        ! An exception handler should be built to handle these problems
      END IF

      ! Point the current position to the head
      myList % current => myList % head
      ! Set the data
      CALL myList % SetData(inData)

      IF (PRESENT(inKey)) THEN
        CALL myList % SetKey(inKey)
      ELSE
        CALL myList % SetKey(keyStart)
      END IF

      ! Point the next to null and the tail to current
      myList % current % next => NULL()
      myList % tail => myList % current

    ELSE ! the list is not empty

      ! Then we allocate space for the next item in the list
      ALLOCATE (myList % tail % next,STAT=allocationStatus)
      IF (allocationStatus /= 0) THEN
        PRINT *, 'MODULE LinkedListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
        ! An exception handler should be built to handle these problems
      END IF

      !  Temporarily hold onto the tail in case we need the key (if inKey is not present)
      previous => myList % tail
      ! Reassign the tail
      myList % tail => myList % tail % next

      ! Set the current to the tail
      myList % current => myList % tail

      ! Fill in the data
      CALL myList % SetData(inData)

      ! Fill in the key information
      IF (PRESENT(inKey)) THEN
        CALL myList % SetKey(inKey)
      ELSE
        CALL myList % SetKey(previous % key + keyInc)
      END IF

      ! Point the next to null and the tail to current
      myList % current % next => NULL()

    END IF

  END SUBROUTINE AddToList_LinkedList
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R RemoveCurrent
!
!> \fn RemoveCurrent_LinkedList
!! Removes the current entry from the list and patches the previous with the next.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % RemoveCurrent( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myList <td> LinkedList <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE RemoveCurrent_LinkedList(myList)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList
    ! LOCAL
    TYPE(Record),POINTER :: previous,pNext
    INTEGER                 :: currentKey,thisKey

    CALL myList % GetKey(currentKey)

    ! Check to see if this list is empty
    IF (myList % ListIsEmpty()) THEN

      PRINT *, 'Module LinkedListClass.f90 : S/R RemoveCurrent : List is empty. Nothing to remove'
      RETURN

    ELSE ! the list is not empty

      CALL myList % MoveToHead() ! Rewind the list

      ! Get the key for this list item
      CALL myList % GetKey(thisKey)

      ! Check if we are trying to remove the head of the list
      IF (thisKey == currentKey) THEN

        ! temporarily point to the next in the list
        pNext => myList % current % next

        ! Deallocate memory pointed to by current position
        DEALLOCATE (myList % current)

        ! Update current position
        myList % head => pNext ! Reset the head of the list

        RETURN
      END IF

      ! If the execution of the code has arrived here, then we are not removing the head of the
      ! list.
      ! Hang on to the head as the previous
      previous => myList % current
      CALL myList % MoveToNext()

      DO WHILE (ASSOCIATED(myList % current))

        ! Get the key for this list item
        CALL myList % GetKey(thisKey)

        ! Check if we are trying to remove the head of the list
        IF (thisKey == currentKey) THEN

          ! temporarily point to the next in the list
          pNext => myList % current % next

          ! Patch the previous item to the next item
          previous % next => pNext

          IF (.NOT. ASSOCIATED(pNext)) THEN
            myList % tail => previous
          END IF

          ! Deallocate memory pointed to by current position
          DEALLOCATE (myList % current)

          EXIT
        ELSE

          previous => myList % current
          CALL myList % moveToNext()

        END IF

      END DO

    END IF

  END SUBROUTINE RemoveCurrent_LinkedList
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R MoveToNext
!
!> \fn MoveToNext_LinkedList
!! Advances the current pointer in the LinkedList to the next position
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToNext( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myList <td> LinkedList <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE MoveToNext_LinkedList(myList)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList

    myList % current => myList % current % next

  END SUBROUTINE MoveToNext_LinkedList
!
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R MoveToHead
!
!> \fn MoveToHead_LinkedList
!! Advances the current pointer in the LinkedList to the start of the list  (list head)
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToHead( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myList <td> LinkedList <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE MoveToHead_LinkedList(myList)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList

    myList % current => myList % head

  END SUBROUTINE MoveToHead_LinkedList
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R MoveToTail
!
!> \fn MoveToTail_LinkedList
!! Advances the current pointer in the LinkedList to the end of the list (list tail)
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToTail( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myList <td> LinkedList <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE MoveToTail_LinkedList(myList)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList

    myList % current => myList % tail

  END SUBROUTINE MoveToTail_LinkedList
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R GetCount
!
!> \fn GetCount_LinkedList
!! Counts the number of associated entries in the list
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!! <B>INTEGER</B>          :: nItems <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetCount( nItems ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in <th> myList <td> LinkedList <td>
!!   <tr> <td> out <th> nItems <td> INTEGER <td> The number of entries in the list
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE GetCount_LinkedList(myList,numberOfListItems)
    IMPLICIT NONE
    CLASS(LinkedList)  :: myList
    INTEGER,INTENT(out) :: numberOfListItems

    numberOfListItems = 0 ! Initialize the number of list items
    ! Check to see if this list is empty
    IF (myList % ListIsEmpty()) THEN

      PRINT *, 'Module LinkedListClass.f90 : S/R ListCount : List is empty.'
      RETURN

    ELSE ! the list is not empty

      CALL myList % MoveToHead() ! Rewind the list

      DO WHILE (ASSOCIATED(myList % current))

        numberOfListItems = numberOfListItems + 1
        CALL myList % moveToNext()

      END DO

    END IF

  END SUBROUTINE GetCount_LinkedList
!
!> \addtogroup SELF_LinkedList
!! @{
! ================================================================================================ !
! S/R PrintToScreen
!
!> \fn PrintToScreen_LinkedList
!! Prints the list data and keys to the screen
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(LinkedList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % PrintToScreen(  ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in <th> myList <td> LinkedList <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE PrintToScreen_LinkedList(myList)
    IMPLICIT NONE
    CLASS(LinkedList) :: myList

    myList % current => myList % head

    PRINT *, '          Data        Key'
    DO WHILE (ASSOCIATED(myList % current))

      PRINT *, myList % current % listData,myList % current % key

      CALL myList % MoveToNext()

    END DO

  END SUBROUTINE PrintToScreen_LinkedList
!
END MODULE SELF_LinkedList
