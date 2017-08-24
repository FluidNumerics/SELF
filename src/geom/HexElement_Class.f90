! HexElement_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file HexElement_Class.f90
!! Contains the \ref HexElement_Class module, and <BR>
!! defines the HexElement, HexElementRecord, and HexElementList data-structures.

!> \defgroup HexElement_Class HexElement_Class 
!! This module defines the HexElement, HexElementRecord, and HexElementList data-structures
!! and its associated routines.

MODULE HexElement_Class
 
! src/common/
USE ModelPrecision
USE ConstantsDictionary


IMPLICIT NONE

!> \addtogroup HexElement_Class 
!! @{

!> \struct HexElement
!!  The HexElement data structure defines attributes needed to describe a quadrilateral element.
!!
!!  HexElements, elements, edges, and faces form the foundation of describing an unstructured mesh. 
!!  The relationship betweens nodes and elements and nodes and edges (or faces in 3-D) define the 
!!  connectivity in an unstructured mesh. In this data structure a HexElement is defined through an
!!  integer ID, its four corner nodes, and its neighbors. 
!!  
!!
!! <H2> HexElement </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> elementID <td> INTEGER  <td> ID that uniquely identifies the element in a mesh
!!       <tr> <th> nodeIDs(1:8) <td> INTEGER <td> ID's of the eight corner nodes that comprise this element
!!       <tr> <th> neighbors(1:6) <td> INTEGER <td> ID's of the six elements that share an edge with 
!!                                                  this element
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref HexElement_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_HexElement
!!    </table>
!!

!>@}
   TYPE HexElement
      INTEGER :: elementID
      INTEGER :: nodeIDs(1:8)   ! Corner HexElement ID's
      INTEGER :: neighbors(1:6) ! Element IDs for the neighbors

      CONTAINS
      PROCEDURE :: Initialize => Initialize_HexElement

   END TYPE HexElement
!> \addtogroup HexElement_Class 
!! @{

!> \struct HexElementRecord
!!  An extension of the node class that includes a pointer to another HexElementRecord so that a
!!  LinkedList style of data storage can be used.
!!
!!  This data structure inherits all of the attributes of the HexElement, but also has a pointer to
!!  another HexElementRecord ("next"). In the anticipation of an unstructured mesh generator being included
!!  with the SELF, it will be necessary to insert mesh primitives on-the-fly without knowing the memory
!!  requirements a'priori.
!!
!! <H2> HexElement </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> next <td> HexElementRecord, POINTER  <td> Pointer to the next HexElementRecord in a Linked-List
!!                                                     of HexElementRecords
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref HexElement_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_HexElementRecord
!!    </table>
!!

!>@}
   TYPE, EXTENDS( HexElement ) :: HexElementRecord 
      TYPE( HexElementRecord ), POINTER :: next => NULL( )
      CONTAINS
      PROCEDURE :: Initialize => Initialize_HexElementRecord
   END TYPE HexElementRecord
!> \addtogroup HexElement_Class 
!! @{

!> \struct HexElementList
!!  A Linked-List of HexElement Records.
!!
!!  This data structure inherits all of the attributes of the HexElement, but also has a pointer to
!!  another HexElementRecord ("next"). In the anticipation of an unstructured mesh generator being included
!!  with the SELF, it will be necessary to insert mesh primitives on-the-fly without the memory
!!  requirements a'priori.
!!
!! <H2> HexElement </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> head <td> HexElementRecord, POINTER  <td> Pointer to the first HexElementRecord in the Linked-List
!!                                                     of HexElementRecords
!!       <tr> <th> current <td> HexElementRecord, POINTER  <td> Pointer to the current HexElementRecord in the Linked-List
!!                                                     of HexElementRecords
!!       <tr> <th> tail <td> HexElementRecord, POINTER  <td> Pointer to the last HexElementRecord in the Linked-List
!!                                                     of HexElementRecords
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref HexElement_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_HexElementRecord
!!       <tr> <th> Trash <td> Trash_HexElementRecord
!!       <tr> <th> GetElementID <td> GetElementID_HexElementList
!!       <tr> <th> SetElementID <td> SetElementID_HexElementList
!!       <tr> <th> GetNodeIDs <td> GetNodeIDs_HexElementList
!!       <tr> <th> SetNodeIDs <td> SetNodeIDs_HexElementList
!!       <tr> <th> GeNeighbors <td> GetNeighbors_HexElementList
!!       <tr> <th> SetNeighbors <td> SetNeighbors_HexElementList
!!       <tr> <th> ListIsEmpty <td> ListIsEmpty_HexElementList
!!       <tr> <th> AddToList <td> AddToList_HexElementList
!!       <tr> <th> RemoveCurrent <td> RemoveCurrent_HexElementList
!!       <tr> <th> MoveToHead <td> MoveToHead_HexElementList
!!       <tr> <th> MoveToNext <td> MoveToNext_HexElementList
!!       <tr> <th> MoveToTail <td> MoveToTail_HexElementList
!!       <tr> <th> GetCount<td> GetCount_HexElementList
!!    </table>
!!

!>@}
   TYPE HexElementList      
      TYPE( HexElementRecord ), POINTER :: head, tail, current

      CONTAINS

      PROCEDURE :: Initialize => Initialize_HexElementList
      PROCEDURE :: Trash      => Trash_HexElementList

      PROCEDURE :: GetElementID => GetElementID_HexElementList
      PROCEDURE :: SetElementID => SetElementID_HexElementList
      PROCEDURE :: GetNodeIDs   => GetNodeIDs_HexElementList
      PROCEDURE :: SetNodeIDs   => SetNodeIDs_HexElementList
      PROCEDURE :: GetNeighbors => GetNeighbors_HexElementList
      PROCEDURE :: SetNeighbors => SetNeighbors_HexElementList
            
      PROCEDURE :: ListIsEmpty   => ListIsEmpty_HexElementList
      PROCEDURE :: AddToList     => AddToList_HexElementList
      PROCEDURE :: RemoveCurrent => RemoveCurrent_HexElementList
      PROCEDURE :: MoveToHead    => MoveToHead_HexElementList
      PROCEDURE :: MoveToNext    => MoveToNext_HexElementList
      PROCEDURE :: MoveToTail    => MoveToTail_HexElementList
      PROCEDURE :: GetCount      => GetCount_HexElementList

   END TYPE HexElementList


 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_HexElement  
!! Initializes memory held by the attributes of the node to default values.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElement) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myHexElement <td> HexElement <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_HexElement( myElement )
 
   IMPLICIT NONE
   CLASS( HexElement ), INTENT(out) :: myElement
   
      myElement % nodeIDs   = 0
      myElement % neighbors = 0
      myElement % elementID = 0

 END SUBROUTINE Initialize_HexElement
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_HexElementRecord  
!! Initializes memory held by the attributes of the node to default values.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementRecord) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myHexElement <td> HexElementRecord <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_HexElementRecord( myElement )
 
   IMPLICIT NONE
   CLASS( HexElementRecord ), INTENT(out) :: myElement
   
      myElement % nodeIDs   = 0
      myElement % neighbors = 0
      myElement % elementID = 0
      myElement % next => NULL( )

 END SUBROUTINE Initialize_HexElementRecord
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_HexElementList
!! Nullifies the head, tail, and current pointers of the HexElementList.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myHexElement <td> HexElementList <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_HexElementList( myList )

   IMPLICIT NONE
   CLASS( HexElementList ) :: myList
    
      myList % head => NULL( )
      myList % tail => NULL()
      myList % current => NULL( )
  
 END SUBROUTINE Initialize_HexElementList
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_HexElement  
!! Cycles through the HexElementList, frees memory held by entries and the LinkedList, and nullifies 
!! the head, tail, and current pointers.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myHexElement <td> HexElementList <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_HexElementList( myList )

   IMPLICIT NONE
   CLASS( HexElementList ) :: myList
   ! LOCAL
   TYPE( HexElementRecord ), POINTER :: pNext

      ! Set the current position of the list to the head
      myList % current => myList % head
     
      ! Scroll through the list until the current position is nullified
      DO WHILE ( ASSOCIATED( myList % current ) )
         pNext => myList % current % next 
         DEALLOCATE( myList % current ) 
         myList % current => pNext 
      ENDDO
  
 END SUBROUTINE Trash_HexElementList
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R SetElementID
! 
!> \fn SetElementID_HexElementList  
!! Sets the CURRENT elementID in the HexElementList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!! <B>INTEGER</B>               :: elementID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetElementID( ElementID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> HexElementList <td> On output, the current elementID is set to the 
!!                                                          incoming elementID. 
!!   <tr> <td> in <th> elementID <td> INTEGER <td> The unique identifier for the current element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetElementID_HexElementList( myList, ElementID )
 
   IMPLICIT NONE
   CLASS( HexElementList )   :: myList
   INTEGER, INTENT(in)        :: elementID

      myList % current % elementID = elementID

 END SUBROUTINE SetElementID_HexElementList
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R GetElementID
! 
!> \fn GetElementID_HexElementList  
!! Gets the CURRENT elementID in the HexElementList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!! <B>INTEGER</B>               :: elementID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetElementID( ElementID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> HexElementList <td> A LinkedList of HexElementRecords
!!   <tr> <td> out <th> elementID <td> INTEGER <td> The unique identifier for the current element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetElementID_HexElementList( myList, elementID )
 
   IMPLICIT NONE
   CLASS( HexElementList )    :: myList
   INTEGER, INTENT(out)        :: elementID

      elementID = myList % current % elementID

 END SUBROUTINE GetElementID_HexElementList
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R SetNodeIDs
! 
!> \fn SetNodeIDs_HexElementList  
!! Sets the CURRENT nodeIDs in the HexElementList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!! <B>INTEGER</B>               :: nodeIDs(1:8) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetNodeIDs( nodeIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> HexElementList <td> On output, the current nodeIDs is set to the 
!!                                                          incoming nodeIDs. 
!!   <tr> <td> in <th> nodeIDs(1:8) <td> INTEGER <td> The unique identifiers for the four corner nodes
!!                                                    that define the current element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetNodeIDs_HexElementList( myList, nodeIDs )
 
   IMPLICIT NONE
   CLASS( HexElementList )   :: myList
   INTEGER, INTENT(in)       :: nodeIDs(1:8)

      myList % current % nodeIDs = nodeIDs

 END SUBROUTINE SetNodeIDs_HexElementList
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R GetNodeIDs
! 
!> \fn GetNodeIDs_HexElementList  
!! Gets the CURRENT nodeIDs in the HexElementList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!! <B>INTEGER</B>              :: nodeIDs(1:8) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetNodeIDs( nodeIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> HexElementList <td> A LinkedList of HexElementRecords
!!   <tr> <td> out <th> nodeIDs(1:8) <td> INTEGER <td> The unique identifiers for the four corner nodes
!!                                                    that define the current element  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetNodeIDs_HexElementList( myList, nodeIDs )
 
   IMPLICIT NONE
   CLASS( HexElementList )    :: myList
   INTEGER, INTENT(out)        :: nodeIDs(1:8)

      nodeIDs = myList % current % nodeIDs

 END SUBROUTINE GetNodeIDs_HexElementList
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R SetNeighbors
! 
!> \fn SetNeighbors_HexElementList  
!! Sets the CURRENT neighbors in the HexElementList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!! <B>INTEGER</B>               :: neighbors(1:6) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetNeighbors( neighbors ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> HexElementList <td> On output, the current neighbors is set to the 
!!                                                          incoming neighbors. 
!!   <tr> <td> in <th> neighbors(1:6) <td> INTEGER <td> The unique identifiers for the four neighbors
!!                                                      of the current element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetNeighbors_HexElementList( myList, neighbors )
 
   IMPLICIT NONE
   CLASS( HexElementList )   :: myList
   INTEGER, INTENT(in)        :: neighbors(1:6)

      myList % current % neighbors = neighbors

 END SUBROUTINE SetNeighbors_HexElementList
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R GetNeighbors
! 
!> \fn GetNeighbors_HexElementList  
!! Gets the CURRENT neighbors in the HexElementList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!! <B>INTEGER</B>               :: neighbors(1:6) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetNeighbors( Neighbors ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> HexElementList <td> A LinkedList of HexElementRecords
!!   <tr> <td> out <th> neighbors(1:6) <td> INTEGER <td> The unique identifiers for the four neighbors
!!                                                       of the current element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetNeighbors_HexElementList( myList, neighbors )
 
   IMPLICIT NONE
   CLASS( HexElementList )    :: myList
   INTEGER, INTENT(out)        :: neighbors(1:6)

      neighbors = myList % current % neighbors

 END SUBROUTINE GetNeighbors_HexElementList
!
!
!==================================================================================================!
!----------------------------------- Linked List Routines  ----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! Function ListIsEmpty 
! 
!> \fn ListIsEmpty_HexElementList  
!!  Tests if the HexElementList has HexElementRecords and returns a LOGICAL indicating the result.
!! 
!!  The List is deemed "empty" if the head of the list is not associated. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!! <B>LOGICAL</B>        :: TorF<BR>
!!         .... <BR>
!!     TorF = this % ListIsEmpty(  ) <BR>
!!     ! To use this function directly within a conditional, you can do the following <BR>
!!     IF( this % ListIsEmpty( ) )
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> HexElementList <td> An initialized Linked-List of edges 
!!   <tr> <td> out <th> TorF <td> LOGICAL <td> 
!!                      <B>.TRUE.</B> if the head of the list is nullified (not associated),
!!                      <B>.FALSE.</B> if the list has at least one HexElementRecord
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION ListIsEmpty_HexElementList( myList ) RESULT( TorF )

   IMPLICIT NONE
   CLASS( HexElementList ) :: myList
   LOGICAL           :: TorF

      TorF = .NOT.( ASSOCIATED( myList % head  ) )
     
 END FUNCTION ListIsEmpty_HexElementList
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToNext 
! 
!> \fn MoveToNext_HexElementList  
!! Shifts the current pointer in the LinkedList to the current % next position. 
!! 
!!  This routine provides a convenient means of traversing the HexElementList one HexElementRecord at a time.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToNext( ) <BR>
!! 
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE MoveToNext_HexElementList( myList )

   IMPLICIT NONE
   CLASS( HexElementList ) :: myList

     myList % current => myList % current % next

 END SUBROUTINE MoveToNext_HexElementList
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToHead 
! 
!> \fn MoveToHead_HexElementList  
!! Shifts the current pointer in the LinkedList to the head position. 
!! 
!!  This routine provides a convenient means of "rewinding" the HexElementList
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToHead( ) <BR>
!! 
! ================================================================================================ ! 
!>@}
  SUBROUTINE MoveToHead_HexElementList( myList )

   IMPLICIT NONE
   CLASS( HexElementList ) :: myList

     myList % current => myList % head

 END SUBROUTINE MoveToHead_HexElementList
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToTail 
! 
!> \fn MoveToTail_HexElementList  
!! Shifts the current pointer in the LinkedList to the tail position. 
!! 
!!  This routine provides a convenient means of "fast-forwarding" the HexElementList to the last 
!!  HexElementRecord.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToTail( ) <BR>
!! 
! ================================================================================================ ! 
!>@}
 SUBROUTINE MoveToTail_HexElementList( myList )

   IMPLICIT NONE
   CLASS( HexElementList ) :: myList

      myList % current => myList % tail

 END SUBROUTINE MoveToTail_HexElementList
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R AddToList 
! 
!> \fn AddToList_HexElementList  
!! Adds an HexElementRecord to the end of the HexElementList 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!! <B> INTEGER</B>              :: elementID, nodeIDs(1:8), neighbors(1:6) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % AddToList( elementID, nodeIDs, neighbors ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> HexElementList <td> 
!!       <tr> <td> elementID <td> INTEGER  <td> ID that uniquely identifies the element in a mesh
!!       <tr> <td> nodeIDs(1:8) <td> INTEGER <td> ID's of the four corner nodes that comprise this element
!!       <tr> <td> neighbors(1:6) <td> INTEGER <td> ID's of the four elements that share an edge with 
!!                                                  this element
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE AddToList_HexElementList( myList, elementID, nodeIDs, neighbors )
   IMPLICIT NONE
   CLASS( HexElementList ) :: myList
   INTEGER                  :: elementID, nodeIDs(1:8), neighbors(1:6)
   ! LOCAL
   TYPE( HexElementRecord ), POINTER :: previous
   INTEGER                            :: allocationStatus

     ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
     
         ALLOCATE( myList % head, STAT = allocationStatus )
         IF( allocationStatus /=0 )THEN
            PRINT*, 'MODULE HexElementListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
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
            PRINT*, 'MODULE HexElementListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
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

 END SUBROUTINE AddToList_HexElementList
!
!> \addtogroup HexElement_Class 
!! @{ 
! ================================================================================================ !
! S/R RemoveCurrent 
! 
!> \fn RemoveCurrent_HexElementList  
!! Removes the current HexElementRecord from the HexElementList and "patches" the list by adjusting the "next" 
!! pointer of the previous HexElementRecord. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % RemoveCurrent( Inputs/Outputs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> HexElementList <td>  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE RemoveCurrent_HexElementList( myList )
 
   IMPLICIT NONE
   CLASS( HexElementList ) :: myList
   ! LOCAL
   TYPE( HexElementRecord ), POINTER :: previous, pNext
   INTEGER                     :: currentKey, thisKey

      currentKey = myList % current % elementID
     
      ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
     
         PRINT*, 'Module HexElementClass.f90 : S/R RemoveCurrent : List is empty. Nothing to remove'
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

 END SUBROUTINE RemoveCurrent_HexElementList
!
!> \addtogroup HexElement_Class
!! @{ 
! ================================================================================================ !
! Function GetCount
! 
!> \fn GetCount_HexElementList  
!! Cycles through the HexElementList and counts the number of HexElementRecords 
!! 
!! Longer Description 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElementList) :: this <BR>
!! <B>INTEGER</B>        :: nHexElements <BR>
!!         .... <BR>
!!     nHexElements = this % GetCount( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> HexElementList <td> 
!!   <tr> <td> out <th> nHexElements <td> The number of HexElementRecords in the HexElementList <td>  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION GetCount_HexElementList( myList ) RESULT( numberOfHexElements )

   IMPLICIT NONE
   CLASS( HexElementList ) :: myList
   INTEGER           :: numberOfHexElements

      numberOfHexElements = 0 ! Initialize the number of list items
      ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
         PRINT*, 'Module HexElementClass.f90 : S/R GetCount : List is empty.'
         RETURN
      ELSE ! the list is not empty
         CALL myList % MoveToHead( ) ! Rewind the list
         DO WHILE( ASSOCIATED( myList % current ) )
            numberOfHexElements = numberOfHexElements + 1
            CALL myList % moveToNext( )
         ENDDO
      ENDIF

 END FUNCTION GetCount_HexElementList
!
END MODULE HexElement_Class
