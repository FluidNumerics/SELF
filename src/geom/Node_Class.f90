! Node_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>, The Florida State University
! Copyright 2016 Joseph Schoonover <jschoonover@lanl.gov>, Los Alamos National Laboratory
!
! The SELF and accompanying documentation were produced in part under the 
! support of Florida State University and the National Science Foundation 
! through Grant OCE-1049131 during 2015 and in part under the support of the 
! Center for Nonlinear Studies and the Department of Energy through the 
! LANL/LDRD program in 2016.
!
! Node_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
! 
! Licensed under the Apache License, Version 2.0 (the "License"); 
! You may obtain a copy of the License at 
!
! http://www.apache.org/licenses/LICENSE-2.0 
!
! Unless required by applicable law or agreed to in writing, software 
! distributed under the License is distributed on an "AS IS" BASIS, 
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and  
! limitations under the License.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file Node_Class.f90
!! Contains the \ref Node_Class module, and <BR>
!! defines the \ref Node, NodeRecord, and NodeList data-structures.

!> \defgroup Node_Class Node_Class 
!! This module defines the Node, NodeRecord, and NodeList data-structures and its associated routines.

MODULE Node_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE LinkedList_Class

IMPLICIT NONE

!> \addtogroup Node_Class 
!! @{

!> \struct Node
!!  The Node data structure defines attributes and type-bound procedures for working with the "node"
!!  mesh primitive in an unstructured mesh.
!!
!!  Nodes, elements, edges, and faces form the foundation of describing an unstructured mesh. The
!!  relationship betweens nodes and elements and nodes and edges (or faces in 3-D) define the 
!!  connectivity in an unstructured mesh. In this data structure a node is defined through an
!!  integer ID, its type (INTERIOR or BOUNDARY), and its position. 
!!  
!!
!! <H2> Node </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> nodeID <td> INTEGER  <td> ID that uniquely identifies the node in a mesh
!!       <tr> <th> nodeType <td> INTEGER <td> Flag that marks the node as being a node in the 
!!                                            interior of the mesh or on the boundary.
!!       <tr> <th> x <td> REAL(prec) <td> x-position of the node
!!       <tr> <th> y <td> REAL(prec) <td> x-position of the node
!!       <tr> <th> z <td> REAL(prec) <td> x-position of the node
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref Node_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_Node
!!       <tr> <th> ScalePosition <td> ScalePosition_Node
!!    </table>
!!

!>@}
   TYPE Node
      INTEGER    :: nodeID
      INTEGER    :: nodeType ! An INTEGER flag for INTERIOR or BOUNDARY
      REAL(prec) :: x, y, z

      CONTAINS

      PROCEDURE :: Initialize    => Initialize_Node
      PROCEDURE :: ScalePosition => ScalePosition_Node
      
   END TYPE Node
!> \addtogroup Node_Class 
!! @{

!> \struct NodeRecord
!!  An extension of the node class that includes a pointer to another NodeRecord so that a
!!  LinkedList style of data storage can be used.
!!
!!  This data structure inherits all of the attributes of the Node, but also has a pointer to
!!  another NodeRecord ("next"). In the anticipation of an unstructured mesh generator being included
!!  with the SELF, it will be necessary to insert mesh primitives on-the-fly without knowing the memory
!!  requirements a'priori.
!!
!! <H2> Node </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> next <td> NodeRecord, POINTER  <td> Pointer to the next NodeRecord in a Linked-List
!!                                                     of NodeRecords
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref Node_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_NodeRecord
!!    </table>
!!

!>@}
   TYPE, EXTENDS( Node ) :: NodeRecord 
      TYPE( NodeRecord ), POINTER :: next => NULL( )
      CONTAINS
      PROCEDURE :: Initialize => Initialize_NodeRecord
   END TYPE NodeRecord
!> \addtogroup Node_Class 
!! @{

!> \struct NodeList
!!  A Linked-List of Node Records.
!!
!!  This data structure inherits all of the attributes of the Node, but also has a pointer to
!!  another NodeRecord ("next"). In the anticipation of an unstructured mesh generator being included
!!  with the SELF, it will be necessary to insert mesh primitives on-the-fly without the memory
!!  requirements a'priori.
!!
!! <H2> Node </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> head <td> NodeRecord, POINTER  <td> Pointer to the first NodeRecord in the Linked-List
!!                                                     of NodeRecords
!!       <tr> <th> current <td> NodeRecord, POINTER  <td> Pointer to the current NodeRecord in the Linked-List
!!                                                     of NodeRecords
!!       <tr> <th> tail <td> NodeRecord, POINTER  <td> Pointer to the last NodeRecord in the Linked-List
!!                                                     of NodeRecords
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref Node_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_NodeRecord
!!       <tr> <th> Trash <td> Trash_NodeRecord
!!       <tr> <th> GetNodeID <td> GetNodeID_NodeList
!!       <tr> <th> SetNodeID <td> SetNodeID_NodeList
!!       <tr> <th> GetNodeType <td> GetNodeType_NodeList
!!       <tr> <th> SetNodeType <td> SetNodeType_NodeList
!!       <tr> <th> GetNodePosition <td> GetNodePosition_NodeList
!!       <tr> <th> SetNodePosition <td> SetNodePosition_NodeList
!!       <tr> <th> ListIsEmpty <td> ListIsEmpty_NodeList
!!       <tr> <th> AddToList <td> AddToList_NodeList
!!       <tr> <th> RemoveCurrent <td> RemoveCurrent_NodeList
!!       <tr> <th> MoveToHead <td> MoveToHead_NodeList
!!       <tr> <th> MoveToNext <td> MoveToNext_NodeList
!!       <tr> <th> MoveToTail <td> MoveToTail_NodeList
!!       <tr> <th> GetCount<td> GetCount_NodeList
!!    </table>
!!

!>@}
   TYPE NodeList      
      TYPE( NodeRecord ), POINTER :: head, tail, current

      CONTAINS

      PROCEDURE :: Initialize => Initialize_NodeList
      PROCEDURE :: Trash      => Trash_NodeList

      PROCEDURE :: GetNodeID       => GetNodeID_NodeList
      PROCEDURE :: SetNodeID       => SetNodeID_NodeList
      PROCEDURE :: GetNodeType     => GetNodeType_NodeList
      PROCEDURE :: SetNodeType     => SetNodeType_NodeList
      PROCEDURE :: GetNodePosition => GetNodePosition_NodeList
      PROCEDURE :: SetNodePosition => SetNodePosition_NodeList
            
      PROCEDURE :: ListIsEmpty   => ListIsEmpty_NodeList
      PROCEDURE :: AddToList     => AddToList_NodeList
      PROCEDURE :: RemoveCurrent => RemoveCurrent_NodeList
      PROCEDURE :: MoveToHead    => MoveToHead_NodeList
      PROCEDURE :: MoveToNext    => MoveToNext_NodeList
      PROCEDURE :: MoveToTail    => MoveToTail_NodeList
      PROCEDURE :: GetCount      => GetCount_NodeList

   END TYPE NodeList
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_Node  
!! Initializes memory held by the attributes of the node to default values.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Node) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myNode <td> Node <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
SUBROUTINE Initialize_Node( myNode )

   IMPLICIT NONE
   CLASS( Node ), INTENT(out) :: myNode

      myNode % x        = ZERO
      myNode % y        = ZERO
      myNode % z        = ZERO
      myNode % nodeType = INTERIOR
      myNode % nodeID   = 0
     
 END SUBROUTINE Initialize_Node
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_NodeRecord  
!! Initializes memory held by the attributes of the node to default values and nullifies the "next"
!! pointer. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeRecord) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myNode <td> NodeRecord <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_NodeRecord( myNode )

   IMPLICIT NONE
   CLASS( NodeRecord ), INTENT(out) :: myNode

      myNode % x        = ZERO
      myNode % y        = ZERO
      myNode % z        = ZERO
      myNode % nodeType = INTERIOR
      myNode % nodeID   = 0
      myNode % next => NULL()
     
 END SUBROUTINE Initialize_NodeRecord 
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_NodeList
!! Nullifies the head, tail, and current pointers of the NodeList.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myNode <td> NodeList <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_NodeList( myList )

   IMPLICIT NONE
   CLASS( NodeList ) :: myList
    
      myList % head => NULL( )
      myList % tail => NULL()
      myList % current => NULL( )
  
 END SUBROUTINE Initialize_NodeList
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_Node  
!! Cycles through the NodeList, frees memory held by entries and the LinkedList, and nullifies 
!! the head, tail, and current pointers.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myNode <td> NodeList <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_NodeList( myList )

   IMPLICIT NONE
   CLASS( NodeList ) :: myList
   ! LOCAL
   TYPE( NodeRecord ), POINTER :: pNext

      ! Set the current position of the list to the head
      myList % current => myList % head
     
      ! Scroll through the list until the current position is nullified
      DO WHILE ( ASSOCIATED( myList % current ) )
         pNext => myList % current % next 
         DEALLOCATE( myList % current ) 
         myList % current => pNext 
      ENDDO
  
 END SUBROUTINE Trash_NodeList
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R SetNodeID
! 
!> \fn SetNodeID_NodeList  
!! Sets the CURRENT nodeID in the NodeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!! <B>INTEGER</B>        :: nodeID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetNodeID( nodeID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> NodeList <td> On output, the current nodeID is set to the 
!!                                                   incoming nodeID. 
!!   <tr> <td> in <th> nodeID <td> INTEGER <td> The unique identifier for the current node 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetNodeID_NodeList( myList, nodeID )
 
   IMPLICIT NONE
   CLASS( NodeList )   :: myList
   INTEGER, INTENT(in) :: nodeID

      myList % current % nodeID = nodeID

 END SUBROUTINE SetNodeID_NodeList
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R GetNodeID
! 
!> \fn GetNodeID_NodeList  
!! Gets the CURRENT nodeID in the NodeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!! <B>INTEGER</B>        :: nodeID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetNodeID( nodeID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> NodeList <td> A LinkedList of NodeRecords
!!   <tr> <td> out <th> nodeID <td> INTEGER <td> The unique identifier for the current node 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetNodeID_NodeList( myList, nodeID )
 
   IMPLICIT NONE
   CLASS( NodeList )    :: myList
   INTEGER, INTENT(out) :: nodeID

      nodeID = myList % current % nodeID

 END SUBROUTINE GetNodeID_NodeList
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R SetNodeType
! 
!> \fn SetNodeType_NodeList  
!! Sets the CURRENT nodeType in the NodeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!! <B>INTEGER</B>        :: nodeType <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetNodeType( nodeType ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> NodeList <td> On output, the current nodeType is set to the 
!!                                                   incoming nodeType. 
!!   <tr> <td> in <th> nodeType <td> INTEGER <td> A flag designating whether the node is an interior
!!                                                or boundary node
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetNodeType_NodeList( myList, nodeType )
 
   IMPLICIT NONE
   CLASS( NodeList )   :: myList
   INTEGER, INTENT(in) :: nodeType

      myList % current % nodeType = nodeType

 END SUBROUTINE SetNodeType_NodeList
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R GetNodeType
! 
!> \fn GetNodeType_NodeList  
!! Gets the CURRENT nodeType in the NodeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!! <B>INTEGER</B>        :: nodeType <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetNodeType( nodeType ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> NodeList <td> A LinkedList of NodeRecords
!!   <tr> <td> out <th> nodeType <td> INTEGER <td> A flag designating whether the node is an interior
!!                                                 or boundary node
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetNodeType_NodeList( myList, nodeType )
 
   IMPLICIT NONE
   CLASS( NodeList )    :: myList
   INTEGER, INTENT(out) :: nodeType

      nodeType = myList % current % nodeType

 END SUBROUTINE GetNodeType_NodeList
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R SetNodePosition
! 
!> \fn SetNodePosition_NodeList  
!! Sets the CURRENT x, y, z in the NodeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!! <B>REAL(prec)</B>        :: x, y, z <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetNodePosition( x, y, z ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> NodeList <td> On output, the current x, y, z is set to the 
!!                                                   incoming x, y, z. 
!!   <tr> <td> in <th> x, y, z <td> REAL(prec) <td> Physical location of the node
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetNodePosition_NodeList( myList, x, y, z )
 
   IMPLICIT NONE
   CLASS( NodeList )   :: myList
   REAL(prec), INTENT(in) :: x, y, z

      myList % current % x = x
      myList % current % y = y
      myList % current % z = z

 END SUBROUTINE SetNodePosition_NodeList
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R GetNodePosition
! 
!> \fn GetNodePosition_NodeList  
!! Gets the CURRENT x, y, z in the NodeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!! <B>REAL(prec)</B>     :: x, y, z <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetNodePosition( x, y, z ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> NodeList <td> A LinkedList of NodeRecords
!!   <tr> <td> out <th> x, y, z <td> REAL(prec) <td> Physical location of the node
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetNodePosition_NodeList( myList, x, y, z )
 
   IMPLICIT NONE
   CLASS( NodeList )    :: myList
   REAL(prec), INTENT(out) :: x, y, z

      x = myList % current % x
      y = myList % current % y
      z = myList % current % z

 END SUBROUTINE GetNodePosition_NodeList
!
!
!==================================================================================================!
!--------------------------------- Type-Specific Routines  ----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Node_Class
!! @{ 
! ================================================================================================ !
! S/R ScalePosition
! 
!> \fn ScalePosition_Node
!!Brief Description 
!! 
!! Longer Description 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Node) :: this <BR>
!! <B>REAL</B>(prec) :: xScale(1:this % nDim) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ScalePosition( xScale ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> thisNode <td> Node <td> A previously constructed Node data structure
!!   <tr> <td> in <th> xScale(1:thisNode % nDim) <td> REAL(prec) <td> Array containing the factors
!!                     to scale the node position by. 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ScalePosition_Node( thisNode, xScale, yScale, zScale )

   IMPLICIT NONE
   CLASS( Node ), INTENT(inout) :: thisNode
   REAL(prec), INTENT(in)       :: xScale, yScale, zScale
   
      thisNode % x = xScale*thisNode % x
      thisNode % y = yScale*thisNode % y
      thisNode % z = zScale*thisNode % z

 END SUBROUTINE ScalePosition_Node
!
!
!==================================================================================================!
!----------------------------------- Linked List Routines  ----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! Function ListIsEmpty 
! 
!> \fn ListIsEmpty_NodeList  
!!  Tests if the NodeList has NodeRecords and returns a LOGICAL indicating the result.
!! 
!!  The List is deemed "empty" if the head of the list is not associated. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!! <B>LOGICAL</B>        :: TorF<BR>
!!         .... <BR>
!!     TorF = this % ListIsEmpty(  ) <BR>
!!     ! To use this function directly within a conditional, you can do the following <BR>
!!     IF( this % ListIsEmpty( ) )
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> NodeList <td> An initialized Linked-List of edges 
!!   <tr> <td> out <th> TorF <td> LOGICAL <td> 
!!                      <B>.TRUE.</B> if the head of the list is nullified (not associated),
!!                      <B>.FALSE.</B> if the list has at least one NodeRecord
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION ListIsEmpty_NodeList( myList ) RESULT( TorF )

   IMPLICIT NONE
   CLASS( NodeList ) :: myList
   LOGICAL           :: TorF

      TorF = .NOT.( ASSOCIATED( myList % head  ) )
     
 END FUNCTION ListIsEmpty_NodeList
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToNext 
! 
!> \fn MoveToNext_NodeList  
!! Shifts the current pointer in the LinkedList to the current % next position. 
!! 
!!  This routine provides a convenient means of traversing the NodeList one NodeRecord at a time.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToNext( ) <BR>
!! 
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE MoveToNext_NodeList( myList )

   IMPLICIT NONE
   CLASS( NodeList ) :: myList

     myList % current => myList % current % next

 END SUBROUTINE MoveToNext_NodeList
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToHead 
! 
!> \fn MoveToHead_NodeList  
!! Shifts the current pointer in the LinkedList to the head position. 
!! 
!!  This routine provides a convenient means of "rewinding" the NodeList
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToHead( ) <BR>
!! 
! ================================================================================================ ! 
!>@}
  SUBROUTINE MoveToHead_NodeList( myList )

   IMPLICIT NONE
   CLASS( NodeList ) :: myList

     myList % current => myList % head

 END SUBROUTINE MoveToHead_NodeList
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToTail 
! 
!> \fn MoveToTail_NodeList  
!! Shifts the current pointer in the LinkedList to the tail position. 
!! 
!!  This routine provides a convenient means of "fast-forwarding" the NodeList to the last 
!!  NodeRecord.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToTail( ) <BR>
!! 
! ================================================================================================ ! 
!>@}
 SUBROUTINE MoveToTail_NodeList( myList )

   IMPLICIT NONE
   CLASS( NodeList ) :: myList

      myList % current => myList % tail

 END SUBROUTINE MoveToTail_NodeList
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R AddToList 
! 
!> \fn AddToList_NodeList  
!! Adds an NodeRecord to the end of the NodeList 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!! <B>INTEGER</B>        :: nodeID, nodeType <BR>
!! <B>REAL</B>(prec)     :: x, y, z <BR>
!!         .... <BR>
!!     <B>CALL</B> this % AddToList( nodeID, nodeType, x, y, z ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> NodeList <td> 
!!       <tr> <td> nodeID <td> INTEGER  <td> ID that uniquely identifies the node in a mesh
!!       <tr> <td> nodeType <td> INTEGER <td> Flag that marks the node as being a node in the 
!!                                            interior of the mesh or on the boundary.
!!       <tr> <td> x <td> REAL(prec) <td> x-position of the node
!!       <tr> <td> y <td> REAL(prec) <td> x-position of the node
!!       <tr> <td> z <td> REAL(prec) <td> x-position of the node 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE AddToList_NodeList( myList, nodeID, nodeType, x, y, z )
   IMPLICIT NONE
   CLASS( NodeList ) :: myList
   INTEGER           :: nodeID, nodeType
   REAL(prec)        :: x, y, z
   ! LOCAL
   TYPE( NodeRecord ), POINTER :: previous
   INTEGER                     :: allocationStatus

     ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
     
         ALLOCATE( myList % head, STAT = allocationStatus )
         IF( allocationStatus /=0 )THEN
            PRINT*, 'MODULE NodeListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
            ! An exception handler should be built to handle these problems
         ENDIF      
      
         ! Point the current position to the head
         myList % current => myList % head
         ! Set the data
         myList % current % nodeID   = nodeID
         myList % current % nodeType = nodeType
         myList % current % x        = x
         myList % current % y        = y
         myList % current % z        = z
      
         ! Point the next to null and the tail to current
         myList % current % next => NULL( )
         myList % tail => myList % current
        
      ELSE ! the list is not empty
    
         ! Then we allocate space for the next item in the list    
         ALLOCATE( myList % tail % next, STAT = allocationStatus )
         IF( allocationStatus /=0 )THEN
            PRINT*, 'MODULE NodeListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
            ! An exception handler should be built to handle these problems
         ENDIF      
        
         !  Temporarily hold onto the tail in case we need the key (if inKey is not present)  
         previous => myList % tail
         ! Reassign the tail
         myList % tail => myList % tail % next
        
         ! Set the current to the tail
         myList % current => myList % tail
  
         ! Fill in the data
         myList % current % nodeID   = nodeID
         myList % current % nodeType = nodeType
         myList % current % x        = x
         myList % current % y        = y
         myList % current % z        = z
        
         ! Point the next to null and the tail to current
         myList % current % next => NULL( )
        
      ENDIF

 END SUBROUTINE AddToList_NodeList
!
!> \addtogroup Node_Class 
!! @{ 
! ================================================================================================ !
! S/R RemoveCurrent 
! 
!> \fn RemoveCurrent_NodeList  
!! Removes the current NodeRecord from the NodeList and "patches" the list by adjusting the "next" 
!! pointer of the previous NodeRecord. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % RemoveCurrent( Inputs/Outputs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> NodeList <td>  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE RemoveCurrent_NodeList( myList )
 
   IMPLICIT NONE
   CLASS( NodeList ) :: myList
   ! LOCAL
   TYPE( NodeRecord ), POINTER :: previous, pNext
   INTEGER                     :: currentKey, thisKey

      currentKey = myList % current % nodeID
     
      ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
     
         PRINT*, 'Module NodeClass.f90 : S/R RemoveCurrent : List is empty. Nothing to remove'
         RETURN
        
      ELSE ! the list is not empty
    
         CALL myList % MoveToHead( ) ! Rewind the list
         
         thisKey = myList % current % nodeID
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
            thisKey = myList % current % nodeID
           
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

 END SUBROUTINE RemoveCurrent_NodeList
!
!> \addtogroup Node_Class
!! @{ 
! ================================================================================================ !
! Function GetCount
! 
!> \fn GetCount_NodeList  
!! Cycles through the NodeList and counts the number of NodeRecords 
!! 
!! Longer Description 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodeList) :: this <BR>
!! <B>INTEGER</B>        :: nNodes <BR>
!!         .... <BR>
!!     nNodes = this % GetCount( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> NodeList <td> 
!!   <tr> <td> out <th> nNodes <td> The number of NodeRecords in the NodeList <td>  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION GetCount_NodeList( myList ) RESULT( numberOfNodes )

   IMPLICIT NONE
   CLASS( NodeList ) :: myList
   INTEGER           :: numberOfNodes

      numberOfNodes = 0 ! Initialize the number of list items
      ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
         PRINT*, 'Module NodeClass.f90 : S/R GetCount : List is empty.'
         RETURN
      ELSE ! the list is not empty
         CALL myList % MoveToHead( ) ! Rewind the list
         DO WHILE( ASSOCIATED( myList % current ) )
            numberOfNodes = numberOfNodes + 1
            CALL myList % moveToNext( )
         ENDDO
      ENDIF

 END FUNCTION GetCount_NodeList
!
END MODULE Node_Class
