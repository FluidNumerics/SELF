! QuadElement_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file QuadElement_Class.f90
!! Contains the \ref QuadElement_Class module, and <BR>
!! defines the QuadElement, QuadElementRecord, and QuadElementList data-structures.

!> \defgroup QuadElement_Class QuadElement_Class 
!! This module defines the QuadElement, QuadElementRecord, and QuadElementList data-structures
!! and its associated routines.

MODULE QuadElement_Class
 
! src/common/
USE ModelPrecision
USE ConstantsDictionary


IMPLICIT NONE

!> \addtogroup QuadElement_Class 
!! @{

!> \struct QuadElement
!!  The QuadElement data structure defines attributes needed to describe a quadrilateral element.
!!
!!  Nodes, elements, edges, and faces form the foundation of describing an unstructured mesh. 
!!  The relationship betweens nodes and elements and nodes and edges (or faces in 3-D) define the 
!!  connectivity in an unstructured mesh. In this data structure a QuadElement is defined through an
!!  integer ID, its four corner nodes, and its neighbors. 
!!  
!!
!! <H2> QuadElement </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> elementID <td> INTEGER  <td> ID that uniquely identifies the element in a mesh
!!       <tr> <th> nodeIDs(1:4) <td> INTEGER <td> ID's of the four corner nodes that comprise this element
!!       <tr> <th> neighbors(1:4) <td> INTEGER <td> ID's of the four elements that share an edge with 
!!                                                  this element
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref QuadElement_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_QuadElement
!!    </table>
!!

!>@}
   TYPE QuadElement
      INTEGER :: elementID
      INTEGER :: nodeIDs(1:4)   ! Corner QuadElement ID's
      INTEGER :: neighbors(1:4) ! Element IDs for the neighbors

      CONTAINS
      PROCEDURE :: Initialize => Initialize_QuadElement

   END TYPE QuadElement
!> \addtogroup QuadElement_Class 
!! @{

!> \struct QuadElementRecord
!!  An extension of the node class that includes a pointer to another QuadElementRecord so that a
!!  LinkedList style of data storage can be used.
!!
!!  This data structure inherits all of the attributes of the QuadElement, but also has a pointer to
!!  another QuadElementRecord ("next"). In the anticipation of an unstructured mesh generator being included
!!  with the SELF, it will be necessary to insert mesh primitives on-the-fly without knowing the memory
!!  requirements a'priori.
!!
!! <H2> QuadElement </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> next <td> QuadElementRecord, POINTER  <td> Pointer to the next QuadElementRecord in a Linked-List
!!                                                     of QuadElementRecords
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref QuadElement_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_QuadElementRecord
!!    </table>
!!

!>@}
   TYPE, EXTENDS( QuadElement ) :: QuadElementRecord 
      TYPE( QuadElementRecord ), POINTER :: next => NULL( )
      CONTAINS
      PROCEDURE :: Initialize => Initialize_QuadElementRecord
   END TYPE QuadElementRecord
!> \addtogroup QuadElement_Class 
!! @{

!> \struct QuadElementList
!!  A Linked-List of QuadElement Records.
!!
!!  This data structure inherits all of the attributes of the QuadElement, but also has a pointer to
!!  another QuadElementRecord ("next"). In the anticipation of an unstructured mesh generator being included
!!  with the SELF, it will be necessary to insert mesh primitives on-the-fly without the memory
!!  requirements a'priori.
!!
!! <H2> QuadElement </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> head <td> QuadElementRecord, POINTER  <td> Pointer to the first QuadElementRecord in the Linked-List
!!                                                     of QuadElementRecords
!!       <tr> <th> current <td> QuadElementRecord, POINTER  <td> Pointer to the current QuadElementRecord in the Linked-List
!!                                                     of QuadElementRecords
!!       <tr> <th> tail <td> QuadElementRecord, POINTER  <td> Pointer to the last QuadElementRecord in the Linked-List
!!                                                     of QuadElementRecords
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref QuadElement_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_QuadElementRecord
!!       <tr> <th> Trash <td> Trash_QuadElementRecord
!!       <tr> <th> GetElementID <td> GetElementID_QuadElementList
!!       <tr> <th> SetElementID <td> SetElementID_QuadElementList
!!       <tr> <th> GetNodeIDs <td> GetNodeIDs_QuadElementList
!!       <tr> <th> SetNodeIDs <td> SetNodeIDs_QuadElementList
!!       <tr> <th> GeNeighbors <td> GetNeighbors_QuadElementList
!!       <tr> <th> SetNeighbors <td> SetNeighbors_QuadElementList
!!       <tr> <th> ListIsEmpty <td> ListIsEmpty_QuadElementList
!!       <tr> <th> AddToList <td> AddToList_QuadElementList
!!       <tr> <th> RemoveCurrent <td> RemoveCurrent_QuadElementList
!!       <tr> <th> MoveToHead <td> MoveToHead_QuadElementList
!!       <tr> <th> MoveToNext <td> MoveToNext_QuadElementList
!!       <tr> <th> MoveToTail <td> MoveToTail_QuadElementList
!!       <tr> <th> GetCount<td> GetCount_QuadElementList
!!    </table>
!!

!>@}
   TYPE QuadElementList      
      TYPE( QuadElementRecord ), POINTER :: head, tail, current

      CONTAINS

      PROCEDURE :: Initialize => Initialize_QuadElementList
      PROCEDURE :: Trash      => Trash_QuadElementList

      PROCEDURE :: GetElementID => GetElementID_QuadElementList
      PROCEDURE :: SetElementID => SetElementID_QuadElementList
      PROCEDURE :: GetNodeIDs   => GetNodeIDs_QuadElementList
      PROCEDURE :: SetNodeIDs   => SetNodeIDs_QuadElementList
      PROCEDURE :: GetNeighbors => GetNeighbors_QuadElementList
      PROCEDURE :: SetNeighbors => SetNeighbors_QuadElementList
            
      PROCEDURE :: ListIsEmpty   => ListIsEmpty_QuadElementList
      PROCEDURE :: AddToList     => AddToList_QuadElementList
      PROCEDURE :: RemoveCurrent => RemoveCurrent_QuadElementList
      PROCEDURE :: MoveToHead    => MoveToHead_QuadElementList
      PROCEDURE :: MoveToNext    => MoveToNext_QuadElementList
      PROCEDURE :: MoveToTail    => MoveToTail_QuadElementList
      PROCEDURE :: GetCount      => GetCount_QuadElementList

   END TYPE QuadElementList


 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_QuadElement  
!! Initializes memory held by the attributes of the node to default values.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElement) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myQuadElement <td> QuadElement <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_QuadElement( myElement )
 
   IMPLICIT NONE
   CLASS( QuadElement ), INTENT(out) :: myElement
   
      myElement % nodeIDs   = 0
      myElement % neighbors = 0
      myElement % elementID = 0

 END SUBROUTINE Initialize_QuadElement
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_QuadElementRecord  
!! Initializes memory held by the attributes of the node to default values.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementRecord) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myQuadElement <td> QuadElementRecord <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_QuadElementRecord( myElement )
 
   IMPLICIT NONE
   CLASS( QuadElementRecord ), INTENT(out) :: myElement
   
      myElement % nodeIDs   = 0
      myElement % neighbors = 0
      myElement % elementID = 0
      myElement % next => NULL( )

 END SUBROUTINE Initialize_QuadElementRecord
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_QuadElementList
!! Nullifies the head, tail, and current pointers of the QuadElementList.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myQuadElement <td> QuadElementList <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_QuadElementList( myList )

   IMPLICIT NONE
   CLASS( QuadElementList ) :: myList
    
      myList % head => NULL( )
      myList % tail => NULL()
      myList % current => NULL( )
  
 END SUBROUTINE Initialize_QuadElementList
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_QuadElement  
!! Cycles through the QuadElementList, frees memory held by entries and the LinkedList, and nullifies 
!! the head, tail, and current pointers.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myQuadElement <td> QuadElementList <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_QuadElementList( myList )

   IMPLICIT NONE
   CLASS( QuadElementList ) :: myList
   ! LOCAL
   TYPE( QuadElementRecord ), POINTER :: pNext

      ! Set the current position of the list to the head
      myList % current => myList % head
     
      ! Scroll through the list until the current position is nullified
      DO WHILE ( ASSOCIATED( myList % current ) )
         pNext => myList % current % next 
         DEALLOCATE( myList % current ) 
         myList % current => pNext 
      ENDDO
  
 END SUBROUTINE Trash_QuadElementList
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R SetElementID
! 
!> \fn SetElementID_QuadElementList  
!! Sets the CURRENT elementID in the QuadElementList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!! <B>INTEGER</B>               :: elementID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetElementID( ElementID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> QuadElementList <td> On output, the current elementID is set to the 
!!                                                          incoming elementID. 
!!   <tr> <td> in <th> elementID <td> INTEGER <td> The unique identifier for the current element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetElementID_QuadElementList( myList, ElementID )
 
   IMPLICIT NONE
   CLASS( QuadElementList )   :: myList
   INTEGER, INTENT(in)        :: elementID

      myList % current % elementID = elementID

 END SUBROUTINE SetElementID_QuadElementList
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R GetElementID
! 
!> \fn GetElementID_QuadElementList  
!! Gets the CURRENT elementID in the QuadElementList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!! <B>INTEGER</B>               :: elementID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetElementID( ElementID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> QuadElementList <td> A LinkedList of QuadElementRecords
!!   <tr> <td> out <th> elementID <td> INTEGER <td> The unique identifier for the current element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetElementID_QuadElementList( myList, elementID )
 
   IMPLICIT NONE
   CLASS( QuadElementList )    :: myList
   INTEGER, INTENT(out)        :: elementID

      elementID = myList % current % elementID

 END SUBROUTINE GetElementID_QuadElementList
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R SetNodeIDs
! 
!> \fn SetNodeIDs_QuadElementList  
!! Sets the CURRENT nodeIDs in the QuadElementList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!! <B>INTEGER</B>               :: nodeIDs(1:4) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetNodeIDs( nodeIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> QuadElementList <td> On output, the current nodeIDs is set to the 
!!                                                          incoming nodeIDs. 
!!   <tr> <td> in <th> nodeIDs(1:4) <td> INTEGER <td> The unique identifiers for the four corner nodes
!!                                                    that define the current element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetNodeIDs_QuadElementList( myList, nodeIDs )
 
   IMPLICIT NONE
   CLASS( QuadElementList )   :: myList
   INTEGER, INTENT(in)        :: nodeIDs(1:4)

      myList % current % nodeIDs = nodeIDs

 END SUBROUTINE SetNodeIDs_QuadElementList
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R GetNodeIDs
! 
!> \fn GetNodeIDs_QuadElementList  
!! Gets the CURRENT nodeIDs in the QuadElementList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!! <B>INTEGER</B>               :: nodeIDs(1:4) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetNodeIDs( nodeIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> QuadElementList <td> A LinkedList of QuadElementRecords
!!   <tr> <td> out <th> nodeIDs(1:4) <td> INTEGER <td> The unique identifiers for the four corner nodes
!!                                                    that define the current element  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetNodeIDs_QuadElementList( myList, nodeIDs )
 
   IMPLICIT NONE
   CLASS( QuadElementList )    :: myList
   INTEGER, INTENT(out)        :: nodeIDs(1:4)

      nodeIDs = myList % current % nodeIDs

 END SUBROUTINE GetNodeIDs_QuadElementList
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R SetNeighbors
! 
!> \fn SetNeighbors_QuadElementList  
!! Sets the CURRENT neighbors in the QuadElementList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!! <B>INTEGER</B>               :: neighbors(1:4) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetNeighbors( neighbors ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> QuadElementList <td> On output, the current neighbors is set to the 
!!                                                          incoming neighbors. 
!!   <tr> <td> in <th> neighbors(1:4) <td> INTEGER <td> The unique identifiers for the four neighbors
!!                                                      of the current element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetNeighbors_QuadElementList( myList, neighbors )
 
   IMPLICIT NONE
   CLASS( QuadElementList )   :: myList
   INTEGER, INTENT(in)        :: neighbors(1:4)

      myList % current % neighbors = neighbors

 END SUBROUTINE SetNeighbors_QuadElementList
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R GetNeighbors
! 
!> \fn GetNeighbors_QuadElementList  
!! Gets the CURRENT neighbors in the QuadElementList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!! <B>INTEGER</B>               :: neighbors(1:4) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetNeighbors( Neighbors ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> QuadElementList <td> A LinkedList of QuadElementRecords
!!   <tr> <td> out <th> neighbors(1:4) <td> INTEGER <td> The unique identifiers for the four neighbors
!!                                                       of the current element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetNeighbors_QuadElementList( myList, neighbors )
 
   IMPLICIT NONE
   CLASS( QuadElementList )    :: myList
   INTEGER, INTENT(out)        :: neighbors(1:4)

      neighbors = myList % current % neighbors

 END SUBROUTINE GetNeighbors_QuadElementList
!
!
!==================================================================================================!
!----------------------------------- Linked List Routines  ----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! Function ListIsEmpty 
! 
!> \fn ListIsEmpty_QuadElementList  
!!  Tests if the QuadElementList has QuadElementRecords and returns a LOGICAL indicating the result.
!! 
!!  The List is deemed "empty" if the head of the list is not associated. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!! <B>LOGICAL</B>        :: TorF<BR>
!!         .... <BR>
!!     TorF = this % ListIsEmpty(  ) <BR>
!!     ! To use this function directly within a conditional, you can do the following <BR>
!!     IF( this % ListIsEmpty( ) )
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> QuadElementList <td> An initialized Linked-List of edges 
!!   <tr> <td> out <th> TorF <td> LOGICAL <td> 
!!                      <B>.TRUE.</B> if the head of the list is nullified (not associated),
!!                      <B>.FALSE.</B> if the list has at least one QuadElementRecord
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION ListIsEmpty_QuadElementList( myList ) RESULT( TorF )

   IMPLICIT NONE
   CLASS( QuadElementList ) :: myList
   LOGICAL           :: TorF

      TorF = .NOT.( ASSOCIATED( myList % head  ) )
     
 END FUNCTION ListIsEmpty_QuadElementList
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToNext 
! 
!> \fn MoveToNext_QuadElementList  
!! Shifts the current pointer in the LinkedList to the current % next position. 
!! 
!!  This routine provides a convenient means of traversing the QuadElementList one QuadElementRecord at a time.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToNext( ) <BR>
!! 
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE MoveToNext_QuadElementList( myList )

   IMPLICIT NONE
   CLASS( QuadElementList ) :: myList

     myList % current => myList % current % next

 END SUBROUTINE MoveToNext_QuadElementList
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToHead 
! 
!> \fn MoveToHead_QuadElementList  
!! Shifts the current pointer in the LinkedList to the head position. 
!! 
!!  This routine provides a convenient means of "rewinding" the QuadElementList
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToHead( ) <BR>
!! 
! ================================================================================================ ! 
!>@}
  SUBROUTINE MoveToHead_QuadElementList( myList )

   IMPLICIT NONE
   CLASS( QuadElementList ) :: myList

     myList % current => myList % head

 END SUBROUTINE MoveToHead_QuadElementList
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToTail 
! 
!> \fn MoveToTail_QuadElementList  
!! Shifts the current pointer in the LinkedList to the tail position. 
!! 
!!  This routine provides a convenient means of "fast-forwarding" the QuadElementList to the last 
!!  QuadElementRecord.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToTail( ) <BR>
!! 
! ================================================================================================ ! 
!>@}
 SUBROUTINE MoveToTail_QuadElementList( myList )

   IMPLICIT NONE
   CLASS( QuadElementList ) :: myList

      myList % current => myList % tail

 END SUBROUTINE MoveToTail_QuadElementList
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R AddToList 
! 
!> \fn AddToList_QuadElementList  
!! Adds an QuadElementRecord to the end of the QuadElementList 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!! <B> INTEGER</B>              :: elementID, nodeIDs(1:4), neighbors(1:4) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % AddToList( elementID, nodeIDs, neighbors ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> QuadElementList <td> 
!!       <tr> <td> elementID <td> INTEGER  <td> ID that uniquely identifies the element in a mesh
!!       <tr> <td> nodeIDs(1:4) <td> INTEGER <td> ID's of the four corner nodes that comprise this element
!!       <tr> <td> neighbors(1:4) <td> INTEGER <td> ID's of the four elements that share an edge with 
!!                                                  this element
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE AddToList_QuadElementList( myList, elementID, nodeIDs, neighbors )
   IMPLICIT NONE
   CLASS( QuadElementList ) :: myList
   INTEGER                  :: elementID, nodeIDs(1:4), neighbors(1:4)
   ! LOCAL
   TYPE( QuadElementRecord ), POINTER :: previous
   INTEGER                            :: allocationStatus

     ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
     
         ALLOCATE( myList % head, STAT = allocationStatus )
         IF( allocationStatus /=0 )THEN
            PRINT*, 'MODULE QuadElementListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
            ! An exception handler should be built to handle these problems
         ENDIF      
      
         ! Point the current position to the head
         myList % current => myList % head
         ! Set the data
         myList % current % elementID = elementID
         myList % current % nodeIDs   = nodeIDs
         myList % current % neighbors = neighbors
      
         ! Point the next to null and the tail to current
         myList % current % next => NULL( )
         myList % tail => myList % current
        
      ELSE ! the list is not empty
    
         ! Then we allocate space for the next item in the list    
         ALLOCATE( myList % tail % next, STAT = allocationStatus )
         IF( allocationStatus /=0 )THEN
            PRINT*, 'MODULE QuadElementListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
            ! An exception handler should be built to handle these problems
         ENDIF      
        
         !  Temporarily hold onto the tail in case we need the key (if inKey is not present)  
         previous => myList % tail
         ! Reassign the tail
         myList % tail => myList % tail % next
        
         ! Set the current to the tail
         myList % current => myList % tail
  
         ! Fill in the data
         myList % current % elementID = elementID
         myList % current % nodeIDs   = nodeIDs
         myList % current % neighbors = neighbors
        
         ! Point the next to null and the tail to current
         myList % current % next => NULL( )
        
      ENDIF

 END SUBROUTINE AddToList_QuadElementList
!
!> \addtogroup QuadElement_Class 
!! @{ 
! ================================================================================================ !
! S/R RemoveCurrent 
! 
!> \fn RemoveCurrent_QuadElementList  
!! Removes the current QuadElementRecord from the QuadElementList and "patches" the list by adjusting the "next" 
!! pointer of the previous QuadElementRecord. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % RemoveCurrent( Inputs/Outputs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> QuadElementList <td>  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE RemoveCurrent_QuadElementList( myList )
 
   IMPLICIT NONE
   CLASS( QuadElementList ) :: myList
   ! LOCAL
   TYPE( QuadElementRecord ), POINTER :: previous, pNext
   INTEGER                     :: currentKey, thisKey

      currentKey = myList % current % elementID
     
      ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
     
         PRINT*, 'Module QuadElementClass.f90 : S/R RemoveCurrent : List is empty. Nothing to remove'
         RETURN
        
      ELSE ! the list is not empty
    
         CALL myList % MoveToHead( ) ! Rewind the list
         
         ! Get the key for this list item
         thisKey = myList % current % elementID
        
         ! Check if we are trying to remove the head of the list
         IF( thisKey == currentKey )THEN 
         
            ! temporarily point to the next in the list
            pNext => myList % current % next 
           
            ! Deallocate memory pointed to by current position
            DEALLOCATE( myList % current ) 

            ! Update current position
            myList % head => pNext ! Reset the head of the list
          
            RETURN
         ENDIF
        
         ! If the execution of the code has arrived here, then we are not removing the head of the 
         ! list. 
         ! Hang on to the head as the previous
         previous => myList % current
         CALL myList % MoveToNext( )
        
         DO WHILE( ASSOCIATED( myList % current ) )
        
            ! Get the key for this list item
            thisKey = myList % current % elementID
           
            ! Check if we are trying to remove the head of the list
            IF( thisKey == currentKey )THEN 
           
               ! temporarily point to the next in the list
               pNext => myList % current % next 
            
               ! Patch the previous item to the next item
               previous % next => pNext
           
               IF( .NOT.ASSOCIATED(pNext)  )THEN
                  myList % tail => previous
               ENDIF
              
              ! Deallocate memory pointed to by current position
              DEALLOCATE( myList % current ) 
          
              EXIT
            ELSE
           
               previous => myList % current
               CALL myList % moveToNext( )

           ENDIF
        
        ENDDO

      ENDIF

 END SUBROUTINE RemoveCurrent_QuadElementList
!
!> \addtogroup QuadElement_Class
!! @{ 
! ================================================================================================ !
! Function GetCount
! 
!> \fn GetCount_QuadElementList  
!! Cycles through the QuadElementList and counts the number of QuadElementRecords 
!! 
!! Longer Description 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadElementList) :: this <BR>
!! <B>INTEGER</B>        :: nQuadElements <BR>
!!         .... <BR>
!!     nQuadElements = this % GetCount( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> QuadElementList <td> 
!!   <tr> <td> out <th> nQuadElements <td> The number of QuadElementRecords in the QuadElementList <td>  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION GetCount_QuadElementList( myList ) RESULT( numberOfQuadElements )

   IMPLICIT NONE
   CLASS( QuadElementList ) :: myList
   INTEGER           :: numberOfQuadElements

      numberOfQuadElements = 0 ! Initialize the number of list items
      ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
         PRINT*, 'Module QuadElementClass.f90 : S/R GetCount : List is empty.'
         RETURN
      ELSE ! the list is not empty
         CALL myList % MoveToHead( ) ! Rewind the list
         DO WHILE( ASSOCIATED( myList % current ) )
            numberOfQuadElements = numberOfQuadElements + 1
            CALL myList % moveToNext( )
         ENDDO
      ENDIF

 END FUNCTION GetCount_QuadElementList
!
END MODULE QuadElement_Class
