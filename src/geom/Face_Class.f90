! Face_Class.f90
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
! Face_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

!> \file Face_Class.f90
!! Contains the \ref Face_Class module, and <BR>
!! defines the \ref Face, FaceRecord, and FaceList data-structures.

!> \defgroup Face_Class Face_Class 
!! This module defines the Face, FaceRecord, and FaceList data-structures.
!!
!! The Face data structure defines the attributes that make up an face in an unstructured mesh.
!! An FaceRecord appends a pointer to the face data structure so that it can be used in a LinkedList,
!! and the FaceList is the Linked-List of faces. A linked-list is included in the anticipation of
!! use in mesh generation algorithms. The Face and FaceRecords are kept separate so that a
!! spectral element solver on a fixed mesh can make use of arrays of Faces without having to carry
!! around a nullified pointer.
!!

MODULE Face_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary

IMPLICIT NONE
!> \addtogroup Face_Class 
!! @{

!> \struct Face
!!  The Face class defines the attributes necessary to describe an face in an unstructured mesh.
!!
!!  <H3>faceID</H3>
!!  An face in an unstructured mesh should be given a unique ID. The faceID is rather handy when 
!!  using a domain decomposition parallelization strategy; the faceID can serve as a unique tag
!!  for an MPI_SEND/RECV pair when exchanging data across faces in the mesh. <BR>
!! 
!!  <H3>boundaryID</H3>
!!  When handling boundary conditions in the SELF, boundary faces are extracted from the mesh.
!!  The order in which the boundary faces are found provides them with a second ID (boundaryID) that
!!  enables quick cycling over the mesh boundary to enforce boundary conditions or MPI exchanges. <BR>
!!
!!  <H3>nodeIDs</H3>
!!  In the unstructured mesh, an face joins two nodes; the node IDs of the terminating nodes is 
!!  stored in the face data structure. <BR>
!!
!!  <H3>elementIDs</H3>
!!   For conformal unstructured meshes, two elements can share an face. The "primary" and
!!  "secondary" element ID's are stored in the face data structure. If the face is a boundary face,
!!  the secondary element ID is replaced with a boundary condition flag that indicates how to 
!!  implement the boundary condition or if and MPI_SEND/RECV pair is needed. <BR>
!!
!!  <H3>elementSides</H3>
!!  In order to exchange data across an face, the local side ID's of the primary and secondary 
!!  elements is needed in order to pass the correct element-face solution to its neighbor. <BR>
!!
!!  <H3>start and inc</H3>
!!  In spectral elements, the solution along an element-face is a 1-D array of values. If two 
!!  neighboring elements do not have the same orientation, we need to be able to reverse the 
!!  ordering of one of the neighboring element's solutions. To handle this, the face also contains 
!!  a "start" and "increment" value that indicates
!!  the face solution ordering of the secondary element wrt to the primary element. <BR>
!!
!! <H2> Face </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> faceID <td> INTEGER  <td> Unique identifier for the face
!!       <tr> <th> boundaryID <td> INTEGER  <td> Unique identifier for the face (if it is a boundary face)
!!       <tr> <th> nodeIDs(1:4) <td> INTEGER  <td> ID's for the two terminating nodes for this face
!!       <tr> <th> elementIDs(1:2) <td> INTEGER  <td> ID's for the two abutting elements
!!       <tr> <th> elementSides(1:2) <td> INTEGER  <td> Local side ID's for the two abutting elements
!!       <tr> <th> iStart <td> INTEGER  <td> Loop start for solutions on the secondary element side,
!!                                           in the first computational direction
!!       <tr> <th> iInc <td> INTEGER  <td> Loop increment for solutions on the secondary element side,
!!                                           in the first computational direction
!!       <tr> <th> jStart <td> INTEGER  <td> Loop start for solutions on the secondary element side,
!!                                           in the second computational direction
!!       <tr> <th> jInc <td> INTEGER  <td> Loop increment for solutions on the secondary element side,
!!                                           in the second computational direction
!!       <tr> <th> swapDimensions <td> INTEGER <td> A flag used to swap the computational directions
!!                                                  between the primary and secondary elements
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref Face_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_Face
!!    </table>
!!

!>@}
   TYPE Face
      INTEGER      :: faceID            ! The face ID
      INTEGER      :: boundaryID        ! If the face is part of the mesh boundary, the face gets assigned a boundary face ID
      INTEGER      :: nodeIDs(1:4)      ! Node IDs which start and terminate this face
      INTEGER      :: elementIDs(1:2)   ! Neighboring elements IDs across the face
      INTEGER      :: elementSides(1:2) ! Local side IDs for the neighboring elements
      INTEGER      :: iStart, iInc      ! Loop start and increment for the secondary element side (1st computational direction)
      INTEGER      :: jStart, jInc      ! Loop start and increment for the secondary element side (2nd computational direction)
      INTEGER      :: swapDimensions    ! A flag used to swap the computational directions
!                                         between the primary and secondary elements

      CONTAINS
      PROCEDURE :: Initialize => Initialize_Face
   END TYPE Face
!> \addtogroup Face_Class 
!! @{

!> \struct FaceRecord
!!  An extension of the face class that includes a pointer to another FaceRecord so that a
!!  LinkedList style of data storage can be used.
!!
!!  This data structure inherits all of the attributes of the Face, but also has a pointer to
!!  another FaceRecord ("next"). In the anticipation of an unstructured mesh generator being included
!!  with the SELF, it will be necessary to insert mesh primitives on-the-fly without knowing the memory
!!  requirements a'priori.
!!
!! <H2> Face </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> next <td> FaceRecord, POINTER  <td> Pointer to the next FaceRecord in a Linked-List
!!                                                     of FaceRecords
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref Face_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_FaceRecord
!!    </table>
!!

!>@}
   TYPE, EXTENDS( Face ) :: FaceRecord 
      TYPE( FaceRecord ), POINTER :: next => NULL( )
      CONTAINS
      PROCEDURE :: Initialize => Initialize_FaceRecord
   END TYPE FaceRecord
!> \addtogroup Face_Class 
!! @{

!> \struct FaceList
!!  A Linked-List of Face Records.
!!
!!  This data structure inherits all of the attributes of the Face, but also has a pointer to
!!  another FaceRecord ("next"). In the anticipation of an unstructured mesh generator being included
!!  with the SELF, it will be necessary to insert mesh primitives on-the-fly without the memory
!!  requirements a'priori.
!!
!! <H2> Face </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> head <td> FaceRecord, POINTER  <td> Pointer to the first FaceRecord in the Linked-List
!!                                                     of FaceRecords
!!       <tr> <th> current <td> FaceRecord, POINTER  <td> Pointer to the current FaceRecord in the Linked-List
!!                                                     of FaceRecords
!!       <tr> <th> tail <td> FaceRecord, POINTER  <td> Pointer to the last FaceRecord in the Linked-List
!!                                                     of FaceRecords
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref Face_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_FaceRecord
!!       <tr> <th> Trash <td> Trash_FaceRecord
!!       <tr> <th> GetFaceID <td> GetFaceID_FaceList
!!       <tr> <th> SetFaceID <td> SetFaceID_FaceList
!!       <tr> <th> GetBoundaryID <td> GetBoundaryID_FaceList
!!       <tr> <th> SetBoundaryID <td> SetBoundaryID_FaceList
!!       <tr> <th> GetNodeIDs <td> GetNodeIDs_FaceList
!!       <tr> <th> SetNodeIDs <td> SetNodeIDs_FaceList
!!       <tr> <th> GetElementIDs <td> GetElementIDs_FaceList
!!       <tr> <th> SetElementIDs <td> SetElementIDs_FaceList
!!       <tr> <th> GetElementSides <td> GetElementSides_FaceList
!!       <tr> <th> SetElementSides <td> SetElementSides_FaceList
!!       <tr> <th> GetStartAndInc <td> GetStartAndInc_FaceList
!!       <tr> <th> SetStartAndInc <td> SetStartAndInc_FaceList
!!       <tr> <th> GetSwapFlag <td> GetSwapFlag_FaceList
!!       <tr> <th> SetSwapFlag <td> SetSwapFlag_FaceList
!!       <tr> <th> ListIsEmpty <td> ListIsEmpty_FaceList
!!       <tr> <th> AddToList <td> AddToList_FaceList
!!       <tr> <th> RemoveCurrent <td> RemoveCurrent_FaceList
!!       <tr> <th> MoveToHead <td> MoveToHead_FaceList
!!       <tr> <th> MoveToNext <td> MoveToNext_FaceList
!!       <tr> <th> MoveToTail <td> MoveToTail_FaceList
!!       <tr> <th> GetCount<td> GetCount_FaceList
!!    </table>
!!

!>@}
   TYPE FaceList      
      TYPE( FaceRecord ), POINTER :: head, tail, current

      CONTAINS

      PROCEDURE :: Initialize => Initialize_FaceList
      PROCEDURE :: Trash      => Trash_FaceList

      PROCEDURE :: GetFaceID       => GetFaceID_FaceList
      PROCEDURE :: SetFaceID       => SetFaceID_FaceList
      PROCEDURE :: GetBoundaryID   => GetBoundaryID_FaceList
      PROCEDURE :: SetBoundaryID   => SetBoundaryID_FaceList
      PROCEDURE :: GetNodeIDs      => GetNodeIDs_FaceList
      PROCEDURE :: SetNodeIDs      => SetNodeIDs_FaceList
      PROCEDURE :: GetElementIDs   => GetElementIDs_FaceList
      PROCEDURE :: SetElementIDs   => SetElementIDs_FaceList
      PROCEDURE :: GetElementSides => GetElementSides_FaceList
      PROCEDURE :: SetElementSides => SetElementSides_FaceList
      PROCEDURE :: GetStartAndInc  => GetStartAndInc_FaceList
      PROCEDURE :: SetStartAndInc  => SetStartAndInc_FaceList
            
      PROCEDURE :: ListIsEmpty   => ListIsEmpty_FaceList
      PROCEDURE :: AddToList     => AddToList_FaceList
      PROCEDURE :: RemoveCurrent => RemoveCurrent_FaceList
      PROCEDURE :: MoveToHead    => MoveToHead_FaceList
      PROCEDURE :: MoveToNext    => MoveToNext_FaceList
      PROCEDURE :: MoveToTail    => MoveToTail_FaceList
      PROCEDURE :: GetCount      => GetCount_FaceList

   END TYPE FaceList
   

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_Face  
!! Initializes memory held by the attributes of the face to default values.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Face) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myFace <td> Face <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_Face( myFace )

   IMPLICIT NONE
   CLASS( Face ), INTENT(out) :: myFace

      myFace % faceID         = 0 
      myFace % nodeIDs        = NO_NORMAL_FLOW
      myFace % elementIDs     = NO_NORMAL_FLOW
      myFace % elementSides   = 0
      myFace % iStart         = 1 
      myFace % iInc           = 1
      myFace % jStart         = 1 
      myFace % jInc           = 1
      myFace % swapDimensions = 1 
      myFace % boundaryID     = 0
     
 END SUBROUTINE Initialize_Face
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_FaceRecord  
!! Initializes memory held by the attributes of the face to default values and nullifies the "next"
!! pointer. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceRecord) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myFace <td> FaceRecord <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_FaceRecord( myFace )

   IMPLICIT NONE
   CLASS( FaceRecord ), INTENT(out) :: myFace

      myFace % faceID         = 0 
      myFace % nodeIDs        = NO_NORMAL_FLOW
      myFace % elementIDs     = NO_NORMAL_FLOW
      myFace % elementSides   = 0
      myFace % iStart         = 1 
      myFace % iInc           = 1
      myFace % jStart         = 1 
      myFace % jInc           = 1
      myFace % swapDimensions = 1 
      myFace % boundaryID     = 0
      myFace % next => NULL()
     
 END SUBROUTINE Initialize_FaceRecord 
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_FaceList
!! Nullifies the head, tail, and current pointers of the FaceList.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myFace <td> FaceList <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_FaceList( myList )

   IMPLICIT NONE
   CLASS( FaceList ) :: myList
    
      myList % head => NULL( )
      myList % tail => NULL()
      myList % current => NULL( )
  
 END SUBROUTINE Initialize_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_Face  
!! Cycles through the FaceList, frees memory held by entries and the LinkedList, and nullifies 
!! the head, tail, and current pointers.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myFace <td> FaceList <td>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_FaceList( myList )

   IMPLICIT NONE
   CLASS( FaceList ) :: myList
   ! LOCAL
   TYPE( FaceRecord ), POINTER :: pNext

      ! Set the current position of the list to the head
      myList % current => myList % head
     
      ! Scroll through the list until the current position is nullified
      DO WHILE ( ASSOCIATED( myList % current ) )
         pNext => myList % current % next 
         DEALLOCATE( myList % current ) 
         myList % current => pNext 
      ENDDO
  
 END SUBROUTINE Trash_FaceList
!
!
!==================================================================================================!
!--------------------------------------- Accessors ------------------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R SetFaceID
! 
!> \fn SetFaceID_FaceList  
!! Sets the CURRENT faceID in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: faceID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetFaceID( faceID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> FaceList <td> On output, the current faceID is set to the 
!!                                                   incoming faceID. 
!!   <tr> <td> in <th> faceID <td> INTEGER <td> The unique identifier for the current face 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetFaceID_FaceList( myList, faceID )
 
   IMPLICIT NONE
   CLASS( FaceList )   :: myList
   INTEGER, INTENT(in) :: faceID

      myList % current % faceID = faceID

 END SUBROUTINE SetFaceID_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R GetFaceID
! 
!> \fn GetFaceID_FaceList  
!! Gets the CURRENT faceID in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: faceID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetFaceID( faceID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> FaceList <td> A LinkedList of FaceRecords
!!   <tr> <td> out <th> faceID <td> INTEGER <td> The unique identifier for the current face 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetFaceID_FaceList( myList, faceID )
 
   IMPLICIT NONE
   CLASS( FaceList )    :: myList
   INTEGER, INTENT(out) :: faceID

      faceID = myList % current % faceID

 END SUBROUTINE GetFaceID_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R SetBoundaryID
! 
!> \fn SetBoundaryID_FaceList  
!! Sets the CURRENT boundaryID in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: boundaryID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetBoundaryID( boundaryID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> FaceList <td> On output, the current boundaryID is set to the 
!!                                                   incoming boundaryID. 
!!   <tr> <td> in <th> boundaryID <td> INTEGER <td> The unique boundary-face identifier for the 
!!                                                  current face 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetBoundaryID_FaceList( myList, boundaryID )
 
   IMPLICIT NONE
   CLASS( FaceList )   :: myList
   INTEGER, INTENT(in) :: boundaryID

      myList % current % boundaryID = boundaryID

 END SUBROUTINE SetBoundaryID_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R GetBoundaryID
! 
!> \fn GetBoundaryID_FaceList  
!! Gets the CURRENT boundaryID in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: boundaryID <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetBoundaryID( boundaryID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> FaceList <td> A LinkedList of FaceRecords
!!   <tr> <td> out <th> boundaryID <td> INTEGER <td> The unique boundary identifier for the current 
!!                                                   face 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetBoundaryID_FaceList( myList, boundaryID )
 
   IMPLICIT NONE
   CLASS( FaceList )    :: myList
   INTEGER, INTENT(out) :: boundaryID

      boundaryID = myList % current % boundaryID

 END SUBROUTINE GetBoundaryID_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R SetNodeIDs
! 
!> \fn SetNodeIDs_FaceList  
!! Sets the CURRENT nodeIDs in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: nodeIDs(1:4) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetNodeIDs( nodeIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> FaceList <td> On output, the current nodeIDs are set to the 
!!                                                   incoming nodeIDs. 
!!   <tr> <td> in <th> nodeIDs(1:4) <td> INTEGER <td> Identifiers for the nodes terminating the face
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetNodeIDs_FaceList( myList, nodeIDs )
 
   IMPLICIT NONE
   CLASS( FaceList )   :: myList
   INTEGER, INTENT(in) :: nodeIDs(1:4)

      myList % current % nodeIDs = nodeIDs

 END SUBROUTINE SetNodeIDs_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R GetNodeIDs
! 
!> \fn GetNodeIDs_FaceList  
!! Gets the CURRENT nodeIDs in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: nodeIDs(1:4) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetNodeIDs( nodeIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> FaceList <td> A LinkedList of FaceRecords
!!   <tr> <td> out <th> nodeIDs(1:4) <td> INTEGER <td> Identifiers for the nodes terminating the face
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetNodeIDs_FaceList( myList, nodeIDs )
 
   IMPLICIT NONE
   CLASS( FaceList )    :: myList
   INTEGER, INTENT(out) :: nodeIDs(1:4)

      nodeIDs = myList % current % nodeIDs

 END SUBROUTINE GetNodeIDs_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R SetElementIDs
! 
!> \fn SetElementIDs_FaceList  
!! Sets the CURRENT elementIDs in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: elementIDs(1:2) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetElementIDs( elementIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> FaceList <td> On output, the current elementIDs are set to the 
!!                                                   incoming elementIDs. 
!!   <tr> <td> in <th> elementIDs <td> INTEGER <td> Identifiers for the abutting elements
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetElementIDs_FaceList( myList, elementIDs )
 
   IMPLICIT NONE
   CLASS( FaceList )   :: myList
   INTEGER, INTENT(in) :: elementIDs(1:2)

      myList % current % elementIDs = elementIDs

 END SUBROUTINE SetElementIDs_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R GetElementIDs
! 
!> \fn GetElementIDs_FaceList  
!! Gets the CURRENT elementIDs in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: elementIDs(1:2) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetElementIDs( elementIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> FaceList <td> A LinkedList of FaceRecords
!!   <tr> <td> out <th> elementIDs <td> INTEGER <td> Identifiers for the abutting elements
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetElementIDs_FaceList( myList, elementIDs )
 
   IMPLICIT NONE
   CLASS( FaceList )    :: myList
   INTEGER, INTENT(out) :: elementIDs(1:2)

      elementIDs = myList % current % elementIDs

 END SUBROUTINE GetElementIDs_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R SetElementSides
! 
!> \fn SetElementSides_FaceList  
!! Sets the CURRENT elementSides in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: elementSides(1:2) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetElementSides( elementSides ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> FaceList <td> On output, the current elementSides are set to the 
!!                                                   incoming elementSides. 
!!   <tr> <td> in <th> elementSides <td> INTEGER <td> Local side ID's for the abutting elements
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetElementSides_FaceList( myList, elementSides )
 
   IMPLICIT NONE
   CLASS( FaceList )   :: myList
   INTEGER, INTENT(in) :: elementSides(1:2)

      myList % current % elementSides = elementSides

 END SUBROUTINE SetElementSides_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R GetElementSides
! 
!> \fn GetElementSides_FaceList  
!! Gets the CURRENT elementSides in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: elementSides(1:2) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetElementSides( elementSides ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> FaceList <td> A LinkedList of FaceRecords
!!   <tr> <td> out <th> elementSides <td> INTEGER <td> Local side ID's for the abutting elements
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetElementSides_FaceList( myList, elementSides )
 
   IMPLICIT NONE
   CLASS( FaceList )    :: myList
   INTEGER, INTENT(out) :: elementSides(1:2)

      elementSides = myList % current % elementSides

 END SUBROUTINE GetElementSides_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R SetStartAndInc
! 
!> \fn SetStartAndInc_FaceList  
!! Sets the CURRENT start and inc in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: start, inc <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetStartAndInc( start, inc ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> FaceList <td> On output, the current start and inc are set to
!!                                                   the incoming start and inc. 
!!   <tr> <td> in <th> iStart <td> INTEGER <td> Loop start index for the secondary element (1st Computational Direction)
!!   <tr> <td> in <th> iInc <td> INTEGER <td> Loop increment for the secondary element (1st Computational Direction)
!!   <tr> <td> in <th> jStart <td> INTEGER <td> Loop start index for the secondary element (2nd Computational Direction)
!!   <tr> <td> in <th> jInc <td> INTEGER <td> Loop increment for the secondary element (2nd Computational Direction)
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetStartAndInc_FaceList( myList, iStart, iInc, jStart, jInc )
 
   IMPLICIT NONE
   CLASS( FaceList )   :: myList
   INTEGER, INTENT(in) :: iStart, iInc, jStart, jInc

      myList % current % iStart = iStart
      myList % current % iInc = iInc 
      myList % current % jStart = jStart
      myList % current % jInc = jInc 

 END SUBROUTINE SetStartAndInc_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R GetStartAndInc
! 
!> \fn GetStartAndInc_FaceList  
!! Gets the CURRENT start and inc in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: nodeIDs(1:2) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetStartAndInc( nodeIDs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> FaceList <td> A LinkedList of FaceRecords
!!   <tr> <td> out <th> iStart <td> INTEGER <td> Loop start index for the secondary element (1st Computational Direction)
!!   <tr> <td> out <th> iInc <td> INTEGER <td> Loop increment for the secondary element (1st Computational Direction)
!!   <tr> <td> out <th> jStart <td> INTEGER <td> Loop start index for the secondary element (2nd Computational Direction)
!!   <tr> <td> out <th> jInc <td> INTEGER <td> Loop increment for the secondary element (2nd Computational Direction)
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetStartAndInc_FaceList( myList, iStart, iInc, jStart, jInc )
 
   IMPLICIT NONE
   CLASS( FaceList )    :: myList
   INTEGER, INTENT(out) :: iStart, iInc, jStart, jInc

      iStart = myList % current % iStart
      iInc   = myList % current % iInc
      jStart = myList % current % jStart
      jInc   = myList % current % jInc

 END SUBROUTINE GetStartAndInc_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R SetSwapFlag
! 
!> \fn SetSwapFlag_FaceList  
!! Sets the CURRENT swapDimensions in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: swapDimensions <BR>
!!         .... <BR>
!!     <B>CALL</B> this % SetSwapFlag( swapDimensions ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> FaceList <td> On output, the current swapDimensions is set to the 
!!                                                   incoming swapDimensions. 
!!   <tr> <td> in <th> swapDimensions <td> INTEGER <td> The unique boundary-face identifier for the 
!!                                                  current face 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetSwapFlag_FaceList( myList, swapDimensions )
 
   IMPLICIT NONE
   CLASS( FaceList )   :: myList
   INTEGER, INTENT(in) :: swapDimensions

      myList % current % swapDimensions = swapDimensions

 END SUBROUTINE SetSwapFlag_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R GetSwapFlag
! 
!> \fn GetSwapFlag_FaceList  
!! Gets the CURRENT swapDimensions in the FaceList. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: swapDimensions <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetSwapFlag( swapDimensions ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> FaceList <td> A LinkedList of FaceRecords
!!   <tr> <td> out <th> swapDimensions <td> INTEGER <td> The unique boundary identifier for the current 
!!                                                   face 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetSwapFlag_FaceList( myList, swapDimensions )
 
   IMPLICIT NONE
   CLASS( FaceList )    :: myList
   INTEGER, INTENT(out) :: swapDimensions

      swapDimensions = myList % current % swapDimensions

 END SUBROUTINE GetSwapFlag_FaceList
!
!
!==================================================================================================!
!----------------------------------- Linked List Routines  ----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! Function ListIsEmpty 
! 
!> \fn ListIsEmpty_FaceList  
!!  Tests if the FaceList has FaceRecords and returns a LOGICAL indicating the result.
!! 
!!  The List is deemed "empty" if the head of the list is not associated. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>LOGICAL</B>        :: TorF<BR>
!!         .... <BR>
!!     TorF = this % ListIsEmpty(  ) <BR>
!!     ! To use this function directly within a conditional, you can do the following <BR>
!!     IF( this % ListIsEmpty( ) )
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myList <td> FaceList <td> An initialized Linked-List of faces 
!!   <tr> <td> out <th> TorF <td> LOGICAL <td> 
!!                      <B>.TRUE.</B> if the head of the list is nullified (not associated),
!!                      <B>.FALSE.</B> if the list has at least one FaceRecord
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION ListIsEmpty_FaceList( myList ) RESULT( TorF )

   IMPLICIT NONE
   CLASS( FaceList ) :: myList
   LOGICAL           :: TorF

      TorF = .NOT.( ASSOCIATED( myList % head  ) )
     
 END FUNCTION ListIsEmpty_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToNext 
! 
!> \fn MoveToNext_FaceList  
!! Shifts the current pointer in the LinkedList to the current % next position. 
!! 
!!  This routine provides a convenient means of traversing the FaceList one FaceRecord at a time.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToNext( ) <BR>
!! 
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE MoveToNext_FaceList( myList )

   IMPLICIT NONE
   CLASS( FaceList ) :: myList

     myList % current => myList % current % next

 END SUBROUTINE MoveToNext_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToHead 
! 
!> \fn MoveToHead_FaceList  
!! Shifts the current pointer in the LinkedList to the head position. 
!! 
!!  This routine provides a convenient means of "rewinding" the FaceList
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToHead( ) <BR>
!! 
! ================================================================================================ ! 
!>@}
  SUBROUTINE MoveToHead_FaceList( myList )

   IMPLICIT NONE
   CLASS( FaceList ) :: myList

     myList % current => myList % head

 END SUBROUTINE MoveToHead_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R MoveToTail 
! 
!> \fn MoveToTail_FaceList  
!! Shifts the current pointer in the LinkedList to the tail position. 
!! 
!!  This routine provides a convenient means of "fast-forwarding" the FaceList to the last 
!!  FaceRecord.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % MoveToTail( ) <BR>
!! 
! ================================================================================================ ! 
!>@}
 SUBROUTINE MoveToTail_FaceList( myList )

   IMPLICIT NONE
   CLASS( FaceList ) :: myList

      myList % current => myList % tail

 END SUBROUTINE MoveToTail_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R AddToList 
! 
!> \fn AddToList_FaceList  
!! Adds an FaceRecord to the end of the FaceList 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: faceID, nodeIDs(1:2), elementIDs(1:2), elementSides(1:2) <BR>
!! <B>INTEGER</B>        :: iStart, iInc, jStart, jInc <BR>
!!         .... <BR>
!!     <B>CALL</B> this % AddToList( faceID, nodeIDs, elementIDs, elementSides, iStart, iInc, jStart, jInc ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> FaceList <td> 
!!       <tr> <th> faceID <td> INTEGER  <td> Unique identifier for the face
!!       <tr> <th> boundaryID <td> INTEGER  <td> Unique identifier for the face (if it is a boundary face)
!!       <tr> <th> nodeIDs(1:4) <td> INTEGER  <td> ID's for the two terminating nodes for this face
!!       <tr> <th> elementIDs(1:2) <td> INTEGER  <td> ID's for the two abutting elements
!!       <tr> <th> elementSides(1:2) <td> INTEGER  <td> Local side ID's for the two abutting elements
!!       <tr> <th> iStart <td> INTEGER  <td> Loop start for solutions on the secondary element side,
!!                                           in the first computational direction
!!       <tr> <th> iInc <td> INTEGER  <td> Loop increment for solutions on the secondary element side,
!!                                           in the first computational direction
!!       <tr> <th> jStart <td> INTEGER  <td> Loop start for solutions on the secondary element side,
!!                                           in the second computational direction
!!       <tr> <th> jInc <td> INTEGER  <td> Loop increment for solutions on the secondary element side,
!!                                           in the second computational direction
!!       <tr> <th> swapDimensions <td> INTEGER <td> A flag used to swap the computational directions
!!                                                  between the primary and secondary elements
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE AddToList_FaceList( myList, faceID, nodeIDs, elementIDs, elementSides,&
                                        iStart, iInc, jStart, jInc, swapDimensions )

   IMPLICIT NONE
   CLASS( FaceList ) :: myList
   INTEGER           :: faceID, nodeIDs(1:4), elementIDs(1:2), elementSides(1:2)
   INTEGER           :: iStart, iInc, jStart, jInc, swapDimensions
   ! LOCAL
   TYPE( FaceRecord ), POINTER :: previous
   INTEGER                     :: allocationStatus

     ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
     
         ALLOCATE( myList % head, STAT = allocationStatus )
         IF( allocationStatus /=0 )THEN
            PRINT*, 'MODULE FaceListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
            ! An exception handler should be built to handle these problems
         ENDIF      
      
         ! Point the current position to the head
         myList % current => myList % head
         ! Set the data
         myList % current % faceID         = faceID
         myList % current % nodeIDs        = nodeIDs
         myList % current % elementIDs     = elementIDs
         myList % current % elementSides   = elementSides
         myList % current % iStart         = iStart
         myList % current % iInc           = iInc
         myList % current % jStart         = jStart
         myList % current % jInc           = jInc
         myList % current % swapDimensions = swapDimensions
         
         ! Point the next to null and the tail to current
         myList % current % next => NULL( )
         myList % tail => myList % current
        
      ELSE ! the list is not empty
    
         ! Then we allocate space for the next item in the list    
         ALLOCATE( myList % tail % next, STAT = allocationStatus )
         IF( allocationStatus /=0 )THEN
            PRINT*, 'MODULE FaceListClass.f90 : S/R AddToList : Memory not allocated for next entry in list.'
            ! An exception handler should be built to handle these problems
         ENDIF      
        
         !  Temporarily hold onto the tail in case we need the key (if inKey is not present)  
         previous => myList % tail
         ! Reassign the tail
         myList % tail => myList % tail % next
        
         ! Set the current to the tail
         myList % current => myList % tail
  
         ! Fill in the data
         myList % current % faceID         = faceID
         myList % current % nodeIDs        = nodeIDs
         myList % current % elementIDs     = elementIDs
         myList % current % elementSides   = elementSides
         myList % current % iStart         = iStart
         myList % current % iInc           = iInc
         myList % current % jStart         = jStart
         myList % current % jInc           = jInc
         myList % current % swapDimensions = swapDimensions
        
         ! Point the next to null and the tail to current
         myList % current % next => NULL( )
        
      ENDIF

 END SUBROUTINE AddToList_FaceList
!
!> \addtogroup Face_Class 
!! @{ 
! ================================================================================================ !
! S/R RemoveCurrent 
! 
!> \fn RemoveCurrent_FaceList  
!! Removes the current FaceRecord from the FaceList and "patches" the list by adjusting the "next" 
!! pointer of the previous FaceRecord. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % RemoveCurrent( Inputs/Outputs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> FaceList <td>  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE RemoveCurrent_FaceList( myList )
 
   IMPLICIT NONE
   CLASS( FaceList ) :: myList
   ! LOCAL
   TYPE( FaceRecord ), POINTER :: previous, pNext
   INTEGER               :: currentKey, thisKey

      currentKey = myList % current % FaceID
     
      ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
     
         PRINT*, 'Module FaceClass.f90 : S/R RemoveCurrent : List is empty. Nothing to remove'
         RETURN
        
      ELSE ! the list is not empty
    
         CALL myList % MoveToHead( ) ! Rewind the list
         
         ! Get the key for this list item
         thisKey = myList % current % FaceID
        
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
            thisKey = myList % current % FaceID
           
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

 END SUBROUTINE RemoveCurrent_FaceList
!
!> \addtogroup Face_Class
!! @{ 
! ================================================================================================ !
! Function GetCount
! 
!> \fn GetCount_FaceList  
!! Cycles through the FaceList and counts the number of FaceRecords 
!! 
!! Longer Description 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(FaceList) :: this <BR>
!! <B>INTEGER</B>        :: nFaces <BR>
!!         .... <BR>
!!     nFaces = this % GetCount( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myList <td> FaceList <td> 
!!   <tr> <td> out <th> nFaces <td> The number of FaceRecords in the FaceList <td>  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION GetCount_FaceList( myList ) RESULT( numberOfFaces )

   IMPLICIT NONE
   CLASS( FaceList ) :: myList
   INTEGER           :: numberOfFaces

      numberOfFaces = 0 ! Initialize the number of list items
      ! Check to see if this list is empty
      IF( myList % ListIsEmpty() )THEN
         PRINT*, 'Module FaceClass.f90 : S/R GetCount : List is empty.'
         RETURN
      ELSE ! the list is not empty
         CALL myList % MoveToHead( ) ! Rewind the list
         DO WHILE( ASSOCIATED( myList % current ) )
            numberOfFaces = numberOfFaces + 1
            CALL myList % moveToNext( )
         ENDDO
      ENDIF

 END FUNCTION GetCount_FaceList

END MODULE Face_Class
