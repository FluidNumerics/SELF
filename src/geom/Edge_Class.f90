! Edge_Class.f90
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
! Edge_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

!> \file Edge_Class.f90
!! Contains the \ref Edge_Class module, and <BR>
!! defines the \ref Edge, \ref EdgeRecord, and \ref EdgeList data-structures.

!> \defgroup Edge_Class Edge_Class 
!! This module defines the Edge, EdgeRecord, and EdgeList data-structures.
!!
!! The Edge data structure defines the attributes that make up an edge in an unstructured mesh.
!! An EdgeRecord appends a pointer to the edge data structure so that it can be used in a LinkedList,
!! and the EdgeList is the Linked-List of edges. A linked-list is included in the anticipation of
!! use in mesh generation algorithms. The Edge and EdgeRecords are kept separate so that a
!! spectral element solver on a fixed mesh can make use of arrays of Edges without having to carry
!! around a nullified pointer.
!!

MODULE Edge_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary

IMPLICIT NONE
!> \addtogroup Edge_Class 
!! @{

!> \struct Edge
!!  The Edge class defines the attributes necessary to describe an edge in an unstructured mesh.
!!
!!  <H3>edgeID</H3>
!!  An edge in an unstructured mesh should be given a unique ID. The edgeID is rather handy when 
!!  using a domain decomposition parallelization strategy; the edgeID can serve as a unique tag
!!  for an MPI_SEND/RECV pair when exchanging data across edges in the mesh. <BR>
!! 
!!  <H3>boundaryID</H3>
!!  When handling boundary conditions in the SELF, boundary edges are extracted from the mesh.
!!  The order in which the boundary edges are found provides them with a second ID (boundaryID) that
!!  enables quick cycling over the mesh boundary to enforce boundary conditions or MPI exchanges. <BR>
!!
!!  <H3>nodeIDs</H3>
!!  In the unstructured mesh, an edge joins two nodes; the node IDs of the terminating nodes is 
!!  stored in the edge data structure. <BR>
!!
!!  <H3>elementIDs</H3>
!!   For conformal unstructured meshes, two elements can share an edge. The "primary" and
!!  "secondary" element ID's are stored in the edge data structure. If the edge is a boundary edge,
!!  the secondary element ID is replaced with a boundary condition flag that indicates how to 
!!  implement the boundary condition or if and MPI_SEND/RECV pair is needed. <BR>
!!
!!  <H3>elementSides</H3>
!!  In order to exchange data across an edge, the local side ID's of the primary and secondary 
!!  elements is needed in order to pass the correct element-edge solution to its neighbor. <BR>
!!
!!  <H3>start and inc</H3>
!!  In spectral elements, the solution along an element-edge is a 1-D array of values. If two 
!!  neighboring elements do not have the same orientation, we need to be able to reverse the 
!!  ordering of one of the neighboring element's solutions. To handle this, the edge also contains 
!!  a "start" and "increment" value that indicates
!!  the edge solution ordering of the secondary element wrt to the primary element. <BR>
!!
!! <H2> Edge </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> edgeID <td> INTEGER  <td> Unique identifier for the edge
!!       <tr> <th> boundaryID <td> INTEGER  <td> Unique identifier for the edge (if it is a boundary edge)
!!       <tr> <th> nodeIDs(1:2) <td> INTEGER  <td> ID's for the two terminating nodes for this edge
!!       <tr> <th> elementIDs(1:2) <td> INTEGER  <td> ID's for the two abutting elements
!!       <tr> <th> elementSides(1:2) <td> INTEGER  <td> Local side ID's for the two abutting elements
!!       <tr> <th> start <td> INTEGER  <td> Loop start for solutions on the secondary element side
!!       <tr> <th> inc <td> INTEGER  <td> Loop increment for solutions on the secondary element side
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref Edge_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_Edge
!!    </table>
!!

!>@}
   TYPE Edge
      INTEGER      :: edgeID            ! The edge ID
      INTEGER      :: boundaryID        ! If the edge is part of the mesh boundary, the edge gets assigned a boundary edge ID
      INTEGER      :: nodeIDs(1:2)      ! Node IDs which start and terminate this edge
      INTEGER      :: elementIDs(1:2)   ! Neighboring elements IDs across the edge
      INTEGER      :: elementSides(1:2) ! Local side IDs for the neighboring elements
      INTEGER      :: start, inc        ! Loop start and increment for the secondary element side

      CONTAINS
      PROCEDURE :: Initialize => Initialize_Edge
   END TYPE Edge
!> \addtogroup Edge_Class 
!! @{

!> \struct EdgeRecord
!!  An extension of the edge class that includes a pointer to another EdgeRecord so that a
!!  LinkedList style of data storage can be used.
!!
!!  This data structure inherits all of the attributes of the Edge, but also has a pointer to
!!  another EdgeRecord ("next"). In the anticipation of an unstructured mesh generator being included
!!  with the SELF, it will be necessary to insert mesh primitives on-the-fly without knowing the memory
!!  requirements a'priori.
!!
!! <H2> Edge </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> next <td> EdgeRecord, POINTER  <td> Pointer to the next EdgeRecord in a Linked-List
!!                                                     of EdgeRecords
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref Edge_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_EdgeRecord
!!    </table>
!!

!>@}
   TYPE, EXTENDS( Edge ) :: EdgeRecord 
      TYPE( EdgeRecord ), POINTER :: next => NULL( )
      CONTAINS
      PROCEDURE :: Initialize => Initialize_EdgeRecord
   END TYPE EdgeRecord
!> \addtogroup Edge_Class 
!! @{

!> \struct EdgeList
!!  A Linked-List of Edge Records.
!!
!!  This data structure inherits all of the attributes of the Edge, but also has a pointer to
!!  another EdgeRecord ("next"). In the anticipation of an unstructured mesh generator being included
!!  with the SELF, it will be necessary to insert mesh primitives on-the-fly without the memory
!!  requirements a'priori.
!!
!! <H2> Edge </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> head <td> EdgeRecord, POINTER  <td> Pointer to the first EdgeRecord in the Linked-List
!!                                                     of EdgeRecords
!!       <tr> <th> current <td> EdgeRecord, POINTER  <td> Pointer to the current EdgeRecord in the Linked-List
!!                                                     of EdgeRecords
!!       <tr> <th> tail <td> EdgeRecord, POINTER  <td> Pointer to the last EdgeRecord in the Linked-List
!!                                                     of EdgeRecords
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref Edge_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_EdgeRecord
!!       <tr> <th> Trash <td> Trash_EdgeRecord
!!       <tr> <th> GetEdgeID <td> GetEdgeID_EdgeList
!!       <tr> <th> SetEdgeID <td> SetEdgeID_EdgeList
!!       <tr> <th> GetBoundaryID <td> GetBoundaryID_EdgeList
!!       <tr> <th> SetBoundaryID <td> SetBoundaryID_EdgeList
!!       <tr> <th> GetNodeIDs <td> GetNodeIDs_EdgeList
!!       <tr> <th> SetNodeIDs <td> SetNodeIDs_EdgeList
!!       <tr> <th> GetElementIDs <td> GetElementIDs_EdgeList
!!       <tr> <th> SetElementIDs <td> SetElementIDs_EdgeList
!!       <tr> <th> GetElementSides <td> GetElementSides_EdgeList
!!       <tr> <th> SetElementSides <td> SetElementSides_EdgeList
!!       <tr> <th> GetStartAndInc <td> GetStartAndInc_EdgeList
!!       <tr> <th> SetStartAndInc <td> SetStartAndInc_EdgeList
!!       <tr> <th> ListIsEmpty <td> ListIsEmpty_EdgeList
!!       <tr> <th> AddToList <td> AddToList_EdgeList
!!       <tr> <th> RemoveCurrent <td> RemoveCurrent_EdgeList
!!       <tr> <th> MoveToHead <td> MoveToHead_EdgeList
!!       <tr> <th> MoveToNext <td> MoveToNext_EdgeList
!!       <tr> <th> MoveToTail <td> MoveToTail_EdgeList
!!       <tr> <th> GetCount<td> GetCount_EdgeList
!!    </table>
!!

!>@}
   TYPE EdgeList      
      TYPE( EdgeRecord ), POINTER :: head, tail, current

      CONTAINS

      PROCEDURE :: Initialize => Initialize_EdgeList
      PROCEDURE :: Trash      => Trash_EdgeList

      PROCEDURE :: GetEdgeID       => GetEdgeID_EdgeList
      PROCEDURE :: SetEdgeID       => SetEdgeID_EdgeList
      PROCEDURE :: GetBoundaryID   => GetBoundaryID_EdgeList
      PROCEDURE :: SetBoundaryID   => SetBoundaryID_EdgeList
      PROCEDURE :: GetNodeIDs      => GetNodeIDs_EdgeList
      PROCEDURE :: SetNodeIDs      => SetNodeIDs_EdgeList
      PROCEDURE :: GetElementIDs   => GetElementIDs_EdgeList
      PROCEDURE :: SetElementIDs   => SetElementIDs_EdgeList
      PROCEDURE :: GetElementSides => GetElementSides_EdgeList
      PROCEDURE :: SetElementSides => SetElementSides_EdgeList
      PROCEDURE :: GetStartAndInc  => GetStartAndInc_EdgeList
      PROCEDURE :: SetStartAndInc  => SetStartAndInc_EdgeList
            
      PROCEDURE :: ListIsEmpty   => ListIsEmpty_EdgeList
      PROCEDURE :: AddToList     => AddToList_EdgeList
      PROCEDURE :: RemoveCurrent => RemoveCurrent_EdgeList
      PROCEDURE :: MoveToHead    => MoveToHead_EdgeList
      PROCEDURE :: MoveToNext    => MoveToNext_EdgeList
      PROCEDURE :: MoveToTail    => MoveToTail_EdgeList
      PROCEDURE :: GetCount      => GetCount_EdgeList

   END TYPE EdgeList
   

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_Edge  
!! Initializes memory held by the attributes of the edge to default values.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Edge) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myEdge <td> Edge <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_Edge( myEdge )

   IMPLICIT NONE
   CLASS( Edge ), INTENT(out) :: myEdge

      myEdge % edgeID       = 0 
      myEdge % nodeIDs      = NO_NORMAL_FLOW
      myEdge % elementIDs   = NO_NORMAL_FLOW
      myEdge % elementSides = 0
      myEdge % start        = 1 
      myEdge % inc          = 1
      myEdge % boundaryID   = 0
     
 END SUBROUTINE Initialize_Edge
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_EdgeRecord  
!! Initializes memory held by the attributes of the edge to default values and nullifies the "next"
!! pointer. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeRecord) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myEdge <td> EdgeRecord <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_EdgeRecord( myEdge )

   IMPLICIT NONE
   CLASS( EdgeRecord ), INTENT(out) :: myEdge

      myEdge % edgeID       = 0 
      myEdge % nodeIDs      = NO_NORMAL_FLOW
      myEdge % elementIDs   = NO_NORMAL_FLOW
      myEdge % elementSides = 0
      myEdge % start        = 1 
      myEdge % inc          = 1
      myEdge % boundaryID   = 0
      myEdge % next => NULL()
     
 END SUBROUTINE Initialize_EdgeRecord 
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_EdgeList
!! Nullifies the head, tail, and current pointers of the EdgeList.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myEdge <td> EdgeList <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_EdgeList( myList )

   IMPLICIT NONE
   CLASS( EdgeList ) :: myList
    
      myList % head => NULL( )
      myList % tail => NULL()
      myList % current => NULL( )
  
 END SUBROUTINE Initialize_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_Edge  
!! Cycles through the EdgeList, frees memory held by entries and the LinkedList, and nullifies 
!! the head, tail, and current pointers.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myEdge <td> EdgeList <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_EdgeList( myList )

   IMPLICIT NONE
   CLASS( EdgeList ) :: myList
   ! LOCAL
   TYPE( EdgeRecord ), POINTER :: pNext

      ! Set the current position of the list to the head
      myList % current => myList % head
     
      ! Scroll through the list until the current position is nullified
      DO WHILE ( ASSOCIATED( myList % current ) )
         pNext => myList % current % next 
         DEALLOCATE( myList % current ) 
         myList % current => pNext 
      ENDDO
  
 END SUBROUTINE Trash_EdgeList
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R SetEdgeID
! 
!> \fn SetEdgeID_EdgeList  
!! Sets the CURRENT edgeID in the EdgeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: edgeID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetEdgeID( edgeID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> EdgeList <td> On output, the current edgeID is set to the 
!!                                                   incoming edgeID. 
!!   <tr> <td> in <th> edgeID <td> INTEGER <td> The unique identifier for the current edge 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetEdgeID_EdgeList( myList, edgeID )
 
   IMPLICIT NONE
   CLASS( EdgeList )   :: myList
   INTEGER, INTENT(in) :: edgeID

      myList % current % edgeID = edgeID

 END SUBROUTINE SetEdgeID_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R GetEdgeID
! 
!> \fn GetEdgeID_EdgeList  
!! Gets the CURRENT edgeID in the EdgeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: edgeID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetEdgeID( edgeID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> EdgeList <td> A LinkedList of EdgeRecords
!!   <tr> <td> out <th> edgeID <td> INTEGER <td> The unique identifier for the current edge 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetEdgeID_EdgeList( myList, edgeID )
 
   IMPLICIT NONE
   CLASS( EdgeList )    :: myList
   INTEGER, INTENT(out) :: edgeID

      edgeID = myList % current % edgeID

 END SUBROUTINE GetEdgeID_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R SetBoundaryID
! 
!> \fn SetBoundaryID_EdgeList  
!! Sets the CURRENT boundaryID in the EdgeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: boundaryID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetBoundaryID( boundaryID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> EdgeList <td> On output, the current boundaryID is set to the 
!!                                                   incoming boundaryID. 
!!   <tr> <td> in <th> boundaryID <td> INTEGER <td> The unique boundary-edge identifier for the 
!!                                                  current edge 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetBoundaryID_EdgeList( myList, boundaryID )
 
   IMPLICIT NONE
   CLASS( EdgeList )   :: myList
   INTEGER, INTENT(in) :: boundaryID

      myList % current % boundaryID = boundaryID

 END SUBROUTINE SetBoundaryID_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R GetBoundaryID
! 
!> \fn GetBoundaryID_EdgeList  
!! Gets the CURRENT boundaryID in the EdgeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: boundaryID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetBoundaryID( boundaryID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> EdgeList <td> A LinkedList of EdgeRecords
!!   <tr> <td> out <th> boundaryID <td> INTEGER <td> The unique boundary identifier for the current 
!!                                                   edge 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetBoundaryID_EdgeList( myList, boundaryID )
 
   IMPLICIT NONE
   CLASS( EdgeList )    :: myList
   INTEGER, INTENT(out) :: boundaryID

      boundaryID = myList % current % boundaryID

 END SUBROUTINE GetBoundaryID_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R SetNodeIDs
! 
!> \fn SetNodeIDs_EdgeList  
!! Sets the CURRENT nodeIDs in the EdgeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: nodeIDs(1:2) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetNodeIDs( nodeIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> EdgeList <td> On output, the current nodeIDs are set to the 
!!                                                   incoming nodeIDs. 
!!   <tr> <td> in <th> nodeIDs(1:2) <td> INTEGER <td> Identifiers for the nodes terminating the edge
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetNodeIDs_EdgeList( myList, nodeIDs )
 
   IMPLICIT NONE
   CLASS( EdgeList )   :: myList
   INTEGER, INTENT(in) :: nodeIDs(1:2)

      myList % current % nodeIDs = nodeIDs

 END SUBROUTINE SetNodeIDs_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R GetNodeIDs
! 
!> \fn GetNodeIDs_EdgeList  
!! Gets the CURRENT nodeIDs in the EdgeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: nodeIDs(1:2) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetNodeIDs( nodeIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> EdgeList <td> A LinkedList of EdgeRecords
!!   <tr> <td> out <th> nodeIDs(1:2) <td> INTEGER <td> Identifiers for the nodes terminating the edge
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetNodeIDs_EdgeList( myList, nodeIDs )
 
   IMPLICIT NONE
   CLASS( EdgeList )    :: myList
   INTEGER, INTENT(out) :: nodeIDs(1:2)

      nodeIDs = myList % current % nodeIDs

 END SUBROUTINE GetNodeIDs_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R SetElementIDs
! 
!> \fn SetElementIDs_EdgeList  
!! Sets the CURRENT elementIDs in the EdgeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: elementIDs(1:2) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetElementIDs( elementIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> EdgeList <td> On output, the current elementIDs are set to the 
!!                                                   incoming elementIDs. 
!!   <tr> <td> in <th> elementIDs <td> INTEGER <td> Identifiers for the abutting elements
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetElementIDs_EdgeList( myList, elementIDs )
 
   IMPLICIT NONE
   CLASS( EdgeList )   :: myList
   INTEGER, INTENT(in) :: elementIDs(1:2)

      myList % current % elementIDs = elementIDs

 END SUBROUTINE SetElementIDs_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R GetElementIDs
! 
!> \fn GetElementIDs_EdgeList  
!! Gets the CURRENT elementIDs in the EdgeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: elementIDs(1:2) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetElementIDs( elementIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> EdgeList <td> A LinkedList of EdgeRecords
!!   <tr> <td> out <th> elementIDs <td> INTEGER <td> Identifiers for the abutting elements
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetElementIDs_EdgeList( myList, elementIDs )
 
   IMPLICIT NONE
   CLASS( EdgeList )    :: myList
   INTEGER, INTENT(out) :: elementIDs(1:2)

      elementIDs = myList % current % elementIDs

 END SUBROUTINE GetElementIDs_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R SetElementSides
! 
!> \fn SetElementSides_EdgeList  
!! Sets the CURRENT elementSides in the EdgeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: elementSides(1:2) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetElementSides( elementSides ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> EdgeList <td> On output, the current elementSides are set to the 
!!                                                   incoming elementSides. 
!!   <tr> <td> in <th> elementSides <td> INTEGER <td> Local side ID's for the abutting elements
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetElementSides_EdgeList( myList, elementSides )
 
   IMPLICIT NONE
   CLASS( EdgeList )   :: myList
   INTEGER, INTENT(in) :: elementSides(1:2)

      myList % current % elementSides = elementSides

 END SUBROUTINE SetElementSides_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R GetElementSides
! 
!> \fn GetElementSides_EdgeList  
!! Gets the CURRENT elementSides in the EdgeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: elementSides(1:2) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetElementSides( elementSides ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> EdgeList <td> A LinkedList of EdgeRecords
!!   <tr> <td> out <th> elementSides <td> INTEGER <td> Local side ID's for the abutting elements
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetElementSides_EdgeList( myList, elementSides )
 
   IMPLICIT NONE
   CLASS( EdgeList )    :: myList
   INTEGER, INTENT(out) :: elementSides(1:2)

      elementSides = myList % current % elementSides

 END SUBROUTINE GetElementSides_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R SetStartAndInc
! 
!> \fn SetStartAndInc_EdgeList  
!! Sets the CURRENT start and inc in the EdgeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: start, inc <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetStartAndInc( start, inc ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> EdgeList <td> On output, the current start and inc are set to
!!                                                   the incoming start and inc. 
!!   <tr> <td> in <th> start <td> INTEGER <td> Loop start index for the secondary element
!!   <tr> <td> in <th> inc <td> INTEGER <td> Loop increment for the secondary elemente
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetStartAndInc_EdgeList( myList, start, inc )
 
   IMPLICIT NONE
   CLASS( EdgeList )   :: myList
   INTEGER, INTENT(in) :: start, inc

      myList % current % start = start
      myList % current % inc   = inc

 END SUBROUTINE SetStartAndInc_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R GetStartAndInc
! 
!> \fn GetStartAndInc_EdgeList  
!! Gets the CURRENT start and inc in the EdgeList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: nodeIDs(1:2) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetStartAndInc( nodeIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> EdgeList <td> A LinkedList of EdgeRecords
!!   <tr> <td> out <th> start <td> INTEGER <td> Loop start index for the secondary element
!!   <tr> <td> out <th> inc <td> INTEGER <td> Loop increment for the secondary element
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetStartAndInc_EdgeList( myList, start, inc )
 
   IMPLICIT NONE
   CLASS( EdgeList )    :: myList
   INTEGER, INTENT(out) :: start, inc

      start = myList % current % start
      inc   = myList % current % inc

 END SUBROUTINE GetStartAndInc_EdgeList
!
!
!==================================================================================================!
!----------------------------------- Linked List Routines  ----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! Function ListIsEmpty 
! 
!> \fn ListIsEmpty_EdgeList  
!!  Tests if the EdgeList has EdgeRecords and returns a LOGICAL indicating the result.
!! 
!!  The List is deemed "empty" if the head of the list is not associated. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>LOGICAL</B>        :: TorF<BR>
!!         .... <BR>
!!     TorF = this % ListIsEmpty(  ) <BR>
!!     ! To use this function directly within a conditional, you can do the following <BR>
!!     IF( this % ListIsEmpty( ) )
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> EdgeList <td> An initialized Linked-List of edges 
!!   <tr> <td> out <th> TorF <td> LOGICAL <td> 
!!                      <B>.TRUE.</B> if the head of the list is nullified (not associated),
!!                      <B>.FALSE.</B> if the list has at least one EdgeRecord
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION ListIsEmpty_EdgeList( myList ) RESULT( TorF )

   IMPLICIT NONE
   CLASS( EdgeList ) :: myList
   LOGICAL           :: TorF

      TorF = .NOT.( ASSOCIATED( myList % head  ) )
     
 END FUNCTION ListIsEmpty_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToNext 
! 
!> \fn MoveToNext_EdgeList  
!! Shifts the current pointer in the LinkedList to the current % next position. 
!! 
!!  This routine provides a convenient means of traversing the EdgeList one EdgeRecord at a time.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToNext( ) <BR>
!! 
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE MoveToNext_EdgeList( myList )

   IMPLICIT NONE
   CLASS( EdgeList ) :: myList

     myList % current => myList % current % next

 END SUBROUTINE MoveToNext_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToHead 
! 
!> \fn MoveToHead_EdgeList  
!! Shifts the current pointer in the LinkedList to the head position. 
!! 
!!  This routine provides a convenient means of "rewinding" the EdgeList
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToHead( ) <BR>
!! 
! ================================================================================================ ! 
!>@}
  SUBROUTINE MoveToHead_EdgeList( myList )

   IMPLICIT NONE
   CLASS( EdgeList ) :: myList

     myList % current => myList % head

 END SUBROUTINE MoveToHead_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToTail 
! 
!> \fn MoveToTail_EdgeList  
!! Shifts the current pointer in the LinkedList to the tail position. 
!! 
!!  This routine provides a convenient means of "fast-forwarding" the EdgeList to the last 
!!  EdgeRecord.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToTail( ) <BR>
!! 
! ================================================================================================ ! 
!>@}
 SUBROUTINE MoveToTail_EdgeList( myList )

   IMPLICIT NONE
   CLASS( EdgeList ) :: myList

      myList % current => myList % tail

 END SUBROUTINE MoveToTail_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R AddToList 
! 
!> \fn AddToList_EdgeList  
!! Adds an EdgeRecord to the end of the EdgeList 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: edgeID, nodeIDs(1:2), elementIDs(1:2), elementSides(1:2), start, inc <BR>
!!         .... <BR>
!!     <B>CALL</B> this % AddToList( edgeID, nodeIDs, elementIDs, elementSides, start, inc ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> EdgeList <td> 
!!   <tr> <td> edgeID <td> INTEGER  <td> Unique identifier for the edge
!!   <tr> <td> boundaryID <td> INTEGER  <td> Unique identifier for the edge (if it is a boundary edge)
!!   <tr> <td> nodeIDs(1:2) <td> INTEGER  <td> ID's for the two terminating nodes for this edge
!!   <tr> <td> elementIDs(1:2) <td> INTEGER  <td> ID's for the two abutting elements
!!   <tr> <td> elementSides(1:2) <td> INTEGER  <td> Local side ID's for the two abutting elements
!!   <tr> <td> start <td> INTEGER  <td> Loop start for solutions on the secondary element side
!!   <tr> <td> inc <td> INTEGER  <td> Loop increment for solutions on the secondary element side 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE AddToList_EdgeList( myList, edgeID, nodeIDs, elementIDs, elementSides, start, inc )

   IMPLICIT NONE
   CLASS( EdgeList ) :: myList
   INTEGER           :: edgeID, nodeIDs(1:2), elementIDs(1:2), elementSides(1:2), start, inc
   ! LOCAL
   TYPE( EdgeRecord ), POINTER :: previous
   INTEGER                     :: allocationStatus

     ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
     
         ALLOCATE( myList % head, STAT = allocationStatus )
         IF( allocationStatus /=0 )THEN
            PRINT*, 'MODULE EdgeListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
            ! An exception handler should be built to handle these problems
         ENDIF      
      
         ! Point the current position to the head
         myList % current => myList % head
         ! Set the data
         myList % current % edgeID       = edgeID
         myList % current % nodeIDs      = nodeIDs
         myList % current % elementIDs   = elementIDs
         myList % current % elementSides = elementSides
         myList % current % start        = start
         myList % current % inc          = inc
      
         ! Point the next to null and the tail to current
         myList % current % next => NULL( )
         myList % tail => myList % current
        
      ELSE ! the list is not empty
    
         ! Then we allocate space for the next item in the list    
         ALLOCATE( myList % tail % next, STAT = allocationStatus )
         IF( allocationStatus /=0 )THEN
            PRINT*, 'MODULE EdgeListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
            ! An exception handler should be built to handle these problems
         ENDIF      
        
         !  Temporarily hold onto the tail in case we need the key (if inKey is not present)  
         previous => myList % tail
         ! Reassign the tail
         myList % tail => myList % tail % next
        
         ! Set the current to the tail
         myList % current => myList % tail
  
         ! Fill in the data
         myList % current % edgeID       = edgeID
         myList % current % nodeIDs      = nodeIDs
         myList % current % elementIDs   = elementIDs
         myList % current % elementSides = elementSides
         myList % current % start        = start
         myList % current % inc          = inc
        
         ! Point the next to null and the tail to current
         myList % current % next => NULL( )
        
      ENDIF

 END SUBROUTINE AddToList_EdgeList
!
!> \addtogroup Edge_Class 
!! @{ 
! ================================================================================================ !
! S/R RemoveCurrent 
! 
!> \fn RemoveCurrent_EdgeList  
!! Removes the current EdgeRecord from the EdgeList and "patches" the list by adjusting the "next" 
!! pointer of the previous EdgeRecord. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % RemoveCurrent( Inputs/Outputs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> EdgeList <td>  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE RemoveCurrent_EdgeList( myList )
 
   IMPLICIT NONE
   CLASS( EdgeList ) :: myList
   ! LOCAL
   TYPE( EdgeRecord ), POINTER :: previous, pNext
   INTEGER               :: currentKey, thisKey

      currentKey = myList % current % edgeID
     
      ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
     
         PRINT*, 'Module EdgeClass.f90 : S/R RemoveCurrent : List is empty. Nothing to remove'
         RETURN
        
      ELSE ! the list is not empty
    
         CALL myList % MoveToHead( ) ! Rewind the list
         
         ! Get the key for this list item
         thisKey = myList % current % edgeID
        
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
            thisKey = myList % current % edgeID
           
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

 END SUBROUTINE RemoveCurrent_EdgeList
!
!> \addtogroup Edge_Class
!! @{ 
! ================================================================================================ !
! Function GetCount
! 
!> \fn GetCount_EdgeList  
!! Cycles through the EdgeList and counts the number of EdgeRecords 
!! 
!! Longer Description 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(EdgeList) :: this <BR>
!! <B>INTEGER</B>        :: nEdges <BR>
!!         .... <BR>
!!     nEdges = this % GetCount( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> EdgeList <td> 
!!   <tr> <td> out <th> nEdges <td> The number of EdgeRecords in the EdgeList <td>  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION GetCount_EdgeList( myList ) RESULT( numberOfEdges )

   IMPLICIT NONE
   CLASS( EdgeList ) :: myList
   INTEGER           :: numberOfEdges

      numberOfEdges = 0 ! Initialize the number of list items
      ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
         PRINT*, 'Module EdgeClass.f90 : S/R GetCount : List is empty.'
         RETURN
      ELSE ! the list is not empty
         CALL myList % MoveToHead( ) ! Rewind the list
         DO WHILE( ASSOCIATED( myList % current ) )
            numberOfEdges = numberOfEdges + 1
            CALL myList % moveToNext( )
         ENDDO
      ENDIF

 END FUNCTION GetCount_EdgeList

END MODULE Edge_Class
