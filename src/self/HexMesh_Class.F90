! HexMesh_Class.f90
!
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE HexMesh_Class

  USE ModelPrecision
  USE ConstantsDictionary
  USE LinkedList_Class
  USE KeyRing_Class
  USE Quadrature
  USE Lagrange_Class
  USE NodalDG_Class
  USE Surfaces_Class
  USE HexElements_Class
  USE Edges_Class
  USE Faces_Class
  USE Nodes_Class
  USE ModelParameters_Class
  USE BoundaryCommunicator_Class
  USE Geom_EquationParser_Class


  IMPLICIT NONE

! HexMesh
!  The HexMesh DATA structure defines attributes needed to describe a conformal unstructured
!  spectral element mesh.
!
!  The HexMesh DATA-structure brings together nodes, elements, faces, and mapped-geometry
!  DATA-structures to completely describe a 3-D spectral element mesh. TYPE-bound procedures are
!  provided for filling in mesh connectivity information including element-to-element neighbors,
!  edge-to-element IDs, element-to-node IDs, and edge-to-node IDs.
!

  TYPE HexMesh
    TYPE( HexElements ) :: elements
    TYPE( Nodes  )      :: nodes
    TYPE( Edges )       :: edges 
    TYPE( Faces )       :: faces
    INTEGER             :: cornerMap(1:3,1:8)
    INTEGER             :: sideMap(1:6)
    INTEGER             :: faceMap(1:4,1:6)
    INTEGER             :: edgeMap(1:2,1:12)
    INTEGER             :: edgeFaceMap(1:2,1:4)

#ifdef HAVE_CUDA
    INTEGER, DEVICE, ALLOCATABLE :: cornerMap_dev(:,:)
    INTEGER, DEVICE, ALLOCATABLE :: sideMap_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: faceMap_dev(:,:)
    INTEGER, DEVICE, ALLOCATABLE :: edgeFaceMap_dev(:,:)
    INTEGER, DEVICE, ALLOCATABLE :: edgeMap_dev(:,:)
#endif

  CONTAINS

    PROCEDURE :: Build => Build_HexMesh
    PROCEDURE :: Trash => Trash_HexMesh

#ifdef HAVE_CUDA
    PROCEDURE :: UpdateDevice => UpdateDevice_HexMesh
    PROCEDURE :: UpdateHost   => UpdateHost_HexMesh
#endif

    ! Built in "strucured" mesh construction
    PROCEDURE :: ConstructStructuredMesh => ConstructStructuredMesh_HexMesh

    ! Mesh Format readers/writers
    PROCEDURE :: ReadSELFMeshFile    => ReadSELFMeshFile_HexMesh
    PROCEDURE :: WriteSELFMeshFile   => WriteSELFMeshFile_HexMesh
    PROCEDURE :: ReadUCDMeshFile     => ReadUCDMeshFile_HexMesh
    PROCEDURE :: ReadTrellisUCDMeshFile => ReadTrellisUCDMeshFile_HexMesh

    ! Connectivity Routines
    PROCEDURE :: ConstructFaces               => ConstructFaces_HexMesh
    PROCEDURE :: ConstructStructuredFaces     => ConstructStructuredFaces_HexMesh
    PROCEDURE :: ConstructDoublyPeriodicFaces => ConstructDoublyPeriodicFaces_HexMesh
    PROCEDURE :: ConstructElementNeighbors    => ConstructElementNeighbors_HexMesh
    PROCEDURE :: DetermineOrientation         => DetermineOrientation_HexMesh
    PROCEDURE :: ScaleTheMesh                 => ScaleTheMesh_HexMesh
    PROCEDURE :: PartitionStructuredElementsAndNodes

    ! Visualization I/O Routines
    PROCEDURE :: WriteTecplot         => WriteTecplot_Hexmesh
    PROCEDURE :: WriteMaterialTecplot => WriteMaterialTecplot_Hexmesh

    PROCEDURE :: NumberOfBoundaryFaces

  END TYPE HexMesh

  INTEGER, PRIVATE, PARAMETER    :: BoundaryFlagStructured = NO_NORMAL_FLOW



CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup HexMesh_Class
!! @{
! ================================================================================================ !
! S/R Build
!
!> \fn Build_HexMesh
!! Builds the attributes and "convenience" arrays of the HexMesh DATA-structure.
!!
!!
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(HexMesh) :: this <BR>
!! <B>INTEGER</B>        :: nNodes, nElements, nFaces, N <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( nNodes, nElements, nFaces, N ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> out <th> myHexMesh <td> HexMesh <td>
!!   <tr> <td> in <th> nNodes <td> INTEGER <td> The number of nodes in the mesh
!!   <tr> <td> in <th> nElements <td> INTEGER <td> The number of elements in the mesh
!!   <tr> <td> in <th> nFaces <td> INTEGER <td> The number of faces in the mesh
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree for the geometry within each element
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE Build_HexMesh( myHexMesh, nNodes, nElements, nFaces, N )

    IMPLICIT NONE
    CLASS(HexMesh), INTENT(out) :: myHexMesh
    INTEGER, INTENT(in)         :: nNodes, nElements, nFaces, N

    ! A hexahedron element (hex-element for short) has six faces. Each face has geometry that
    ! requires the USE of two computational coordinates. The third computational coordinate is
    ! fixed. The sideMap gives the value of the remaining computational coordinate for each face
    ! of the hex-element.
    !
    ! The face ordering is  1 - South, 2 - East, 3 - North, 4 - West, 5 - Bottom, 6 - Top
    myHexmesh % sideMap(1:6) = (/ 0, N, N, 0, 0, N /)

    ! The eight corner nodes in an approximation that USEs Gauss-Lobatto points, typically CG-TYPE
    ! methods, have fixed computational coordinates. The corner node numbering  starts in the
    ! southwest corner (in the computational grid) of the bottom face, proceeds counter clockwise
    ! around the face, and THEN repeats for the top face. This gives the local ID for the corner
    ! nodes as
    !
    ! Bottom, SouthWest = 1
    ! Bottom, SouthEast = 2
    ! Bottom, NorthEast = 3
    ! Bottom, NorthWest = 4
    ! Top, SouthWest = 5
    ! Top, SouthEast = 6
    ! Top, NorthEast = 7
    ! Top, NorthWest = 8
    !
    ! The computational coordinates for the corner nodes is given assuming a Gauss-Lobatto
    ! computational mesh is USEd. Note that for a Gauss mesh, the corner nodes are not included.

    myHexmesh % cornerMap(1, 1:8) = (/ 0, N, N, 0, 0, N, N, 0 /)
    myHexmesh % cornerMap(2, 1:8) = (/ 0, 0, N, N, 0, 0, N, N /)
    myHexmesh % cornerMap(3, 1:8) = (/ 0, 0, 0, 0, N, N, N, N /)

    ! Mesh construction usually begins with the specIFication of elements and the corner nodes, in
    ! addition to the element geometry. From the element-to-node connectivity, we need to construct
    ! the unique faces in the mesh and specIFy the abutting elements and their relative orientation.
    ! This procedure is aided by a "convenience array" that lists the local corner node IDs in the
    ! local counter-clockwise direction beginning in the local southwest corner of the face. When
    ! two elements share a face, the global node IDs for each element can be found using this
    ! convenience array (called "faceMap") and the relative orientation of the neighboring elements
    ! can be determined. The first index cycles over the nodes which make up the face in a
    ! counterclockwise direction. The second index cycles over the faces in the element.

    myHexMesh % faceMap(1:4, south)  = (/ 1, 2, 6, 5 /)
    myHexMesh % faceMap(1:4, east)   = (/ 2, 3, 7, 6 /)
    myHexMesh % faceMap(1:4, north)  = (/ 4, 3, 7, 8 /)
    myHexMesh % faceMap(1:4, west)   = (/ 1, 4, 8, 5 /)
    myHexMesh % faceMap(1:4, bottom) = (/ 1, 2, 3, 4 /)
    myHexMesh % faceMap(1:4, top)    = (/ 5, 6, 7, 8 /)

    ! Each of the faces can be identIFied by their four corner nodes. The geometry of the faces
    ! is described using two computational coordinates between [-1,1]X[-1,1]. This 2-D computational
    ! grid has its own "southwest", "southeast", "northeast", and "northwest" identIFications. The
    ! The corner nodes of the faces are labeled in the order mentioned in the previous sentence
    ! ( counterclockwise starting from southwest ). For quick referencing when producing tri-linear
    ! elements, a book-keeping array is USEful for ensuring that we reference each edge of the
    ! face in the order of increasing computational coordinate. This is identical to what is DOne
    ! in the HexMeshCLASS.f90 for the "edgeMap". Here, it is called the "edgeFaceMap".
    ! The first index references the starting(1) or ending(2) node. The second index references
    ! the edge of the face with the first being the southern edge and increasing the second index
    ! proceeds counter-clockwise.

    myHexmesh % edgeFaceMap(1, 1:4) = (/ 1, 2, 4, 1 /)
    myHexmesh % edgeFaceMap(2, 1:4) = (/ 2, 3, 3, 4 /)

    ! the edge map is used to order the 12 edges that make up the hex cell by using the local
    ! corner node ID's
    myHexMesh % edgeMap(1:2, 1)  = (/ 1, 2 /)
    myHexMesh % edgeMap(1:2, 2)  = (/ 2, 3 /)
    myHexMesh % edgeMap(1:2, 3)  = (/ 3, 4 /)
    myHexMesh % edgeMap(1:2, 4)  = (/ 4, 1 /)
    myHexMesh % edgeMap(1:2, 5)  = (/ 5, 6 /)
    myHexMesh % edgeMap(1:2, 6)  = (/ 6, 7 /)
    myHexMesh % edgeMap(1:2, 7)  = (/ 7, 8 /)
    myHexMesh % edgeMap(1:2, 8)  = (/ 8, 5 /)
    myHexMesh % edgeMap(1:2, 9)  = (/ 5, 1 /)
    myHexMesh % edgeMap(1:2, 10) = (/ 2, 6 /)
    myHexMesh % edgeMap(1:2, 11) = (/ 3, 7 /)
    myHexMesh % edgeMap(1:2, 12) = (/ 8, 4 /)

    ! The number of nodes, the number of elements, and the number of faces are stored in this DATA
    ! structure for convenience. In another implementation (planned for the next version), the
    ! number of elements, nodes, and faces is dynamic; that implementation
    ! requires the USE of dynamic storage, e.g. a linked-list like structure for elements, edges,
    ! and nodes.

    CALL myHexMesh % elements % Build( N, nElements )
    CALL myHexmesh % nodes % Build( nNodes )
    CALL myHexMesh % faces % Build( nFaces, N )

#ifdef HAVE_CUDA

    ALLOCATE( myHexMesh % cornerMap_dev(1:3,1:8), &
              myHexMesh % sideMap_dev(1:6), &
              myHexMesh % faceMap_dev(1:4,1:6), &
              myHexMesh % edgeFaceMap_dev(1:2,1:4), &
              myHexMesh % edgeMap_dev(1:2,1:12) )

    myHexMesh % cornerMap_dev   = myHexMesh % cornerMap
    myHexMesh % sideMap_dev     = myHexMesh % sideMap
    myHexMesh % faceMap_dev     = myHexMesh % faceMap
    myHexMesh % edgeFaceMap_dev = myHexMesh % edgeFaceMap
    myHexMesh % edgeMap_dev     = myHexMesh % edgeMap

#endif


  END SUBROUTINE Build_HexMesh
!
!> \addtogroup HexMesh_Class
!! @{
! ================================================================================================ !
! S/R Trash
!
!> \fn Trash_HexMesh
!! Frees memory associated with each of the attributes of the HexMesh DATA structure
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(HexMesh) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myHexMesh <td> HexMesh <td> On output, the memory associated with
!!                         attributes of the HexMesh structure have been freed
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE Trash_HexMesh( myHexMesh )

    IMPLICIT NONE
    CLASS(HexMesh), INTENT(inout) :: myHexMesh

    CALL myHexMesh % elements % Trash( )
    CALL myHexmesh % nodes % Trash( )
    CALL myHexMesh % faces % Trash( )

#ifdef HAVE_CUDA
    DEALLOCATE( myHexMesh % cornerMap_dev, &
                myHexMesh % sideMap_dev, &
                myHexMesh % faceMap_dev, &
                myHexMesh % edgeFaceMap_dev, &
                myHexMesh % edgeMap_dev )
#endif

  END SUBROUTINE Trash_HexMesh
!
#ifdef HAVE_CUDA
  SUBROUTINE UpdateDevice_HexMesh( myHexMesh )
    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh

    CALL myHexMesh % faces % UpdateDevice( )
    CALL myHexMesh % elements % UpdateDevice( )
    CALL myHexMesh % nodes % UpdateDevice( )

  END SUBROUTINE UpdateDevice_HexMesh

  SUBROUTINE UpdateHost_HexMesh( myHexMesh )
    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh

    CALL myHexMesh % faces % UpdateHost( )
    CALL myHexMesh % elements % UpdateHost( )
    CALL myHexMesh % nodes % UpdateHost( )

  END SUBROUTINE UpdateHost_HexMesh
#endif
!
!
!==================================================================================================!
!--------------------------------- TYPE SpecIFic Routines -----------------------------------------!
!==================================================================================================!
!
  SUBROUTINE ConstructNodeToElementConnectivity( myHexMesh )
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh
    ! Local
    INTEGER :: iEl, j, localNodeID, globalNodeID
    INTEGER :: nodeLogic(1:myHexMesh % nodes % nNodes)

   
      nodeLogic = 0
      DO iEl = 1, myHexMesh % elements % nElements 
        DO localNodeID = 1, 8

          globalNodeID = myHexMesh % elements % nodeIDs(localNodeID,iEl)

          myHexMesh % nodes % nElements(globalNodeID) = myHexMesh % nodes % nElements(globalNodeID) + 1
          j = myHexMesh % nodes % nElements(globalNodeID)

          IF( j <= maxNodeValence )THEN
            myHexMesh % nodes % elementIDs(j,globalNodeID)   = iEl
            myHexMesh % nodes % elementNodes(j,globalNodeID) = localNodeID
          ENDIF

        ENDDO
      ENDDO
    
  END SUBROUTINE ConstructNodeToElementConnectivity


!
!> \addtogroup HexMesh_Class
!! @{
! ================================================================================================ !
! S/R ConstructFaces
!
!> \fn ConstructFaces_HexMesh
!! USEs the element-to-node connectivity to construct all of the unique faces in a mesh
!!
!! Similar to the edge construction in a QuadMesh, we loop over the elements and construct a face
!! using the element-to-node connectivity and a convenience array that relates the local nodes
!! to each face. A face is identIFied by its four corner nodes. The four corner node IDs for each
!! elements face are gathered and compared against a list of already identIFied faces. IF the
!! face is not in the list, THEN a new face is added to the list and the element that added
!! the face to the list is designated the "primary" element. IF the face is on the list,
!! the orientation of the secondary element relative to the primary element is determined.
!!
!! This routine depends on <BR>
!!     Module HexMesh_Class : \ref DetermineOrientation_HexMesh
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(DATATYPE) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ConstructFaces( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myHexMesh <td> HexMesh <td>
!!                         On <B>input</B>, the element and node information has been filled in, <BR>
!!                         On <B>output</B>, the unique faces have been identIFied and the face
!!                         information has been filled in.
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE ConstructFaces_HexMesh( myHexMesh )

    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh
    ! LOCAL
    TYPE( KeyRing ) :: KeyCabinet(1:myHexMesh % nodes % nNodes)
    INTEGER :: nEls, nNodes, iEl, nFaces, i, j, k, ii, jj
    INTEGER :: locNodeIDs(1:4), globNodeIDs(1:4)
    INTEGER :: keyRingID, globFaceID
    INTEGER :: N
    LOGICAL :: keyExists

    nNodes = myHexMesh % nodes % nNodes
    nEls   = myHexMesh % elements % nElements
    N      = myHexMesh % elements % N
    nFaces = 0

    DO k = 1, nNodes
      CALL keyCabinet(k) % Build( )
    ENDDO

    DO iEl = 1, nEls ! Loop over the elements in the mesh

      DO k = 1, 6 ! Loop over the faces of each element

        ! In this first step, we want to identIFy the unique global node ID's for this face
        ! To DO this, we start by gathering the local (to the element) node ID's.
        DO j = 1, 4
          locNodeIDs(j) = myHexMesh % faceMap(j,k) ! starting local node for this Face
        ENDDO

        ! Now, we extract, for this element, the global node ID's for each local node ID
        DO j = 1, 4
          globNodeIDs(j) = myHexMesh % elements % nodeIDs( locNodeIDs(j), iEl )
        ENDDO
        ! Our key cabinet has many key-rings with each key-ring containing a set of notched keys.
        ! Each notched key corresponds to a unique face in the mesh. To enable fast searching
        ! for a unique face, we address our key-rings according to the minimum global node ID
        ! that resides on the current face. IF another face shares this minimum global node ID,
        ! THEN it is possible that the face has already been generated. IF not, our search is
        ! limited only to the key-ring with the same key-ring ID. Here, we grab the key-ring ID.
        keyRingID = MINVAL( globNodeIDs )

        ! Now, we check to see IF a notched key already exists with the same set of global ID's
        ! This is DOne by making a call to "FindDATAForNotches". This routine searches through
        ! the current key ring for the key which has the same global ID's (though not
        ! necessarily in the same order). IF the face has been built already, the face ID is
        ! returned in globFaceID and the LOGICAL, "keyExists", is set to TRUE. Otherwise,
        ! globFaceID is set to zero and "keyExsists" is set to FALSE.
        CALL KeyCabinet(keyRingID) % FindDATAForNotches( globNodeIDS, 4, &
          globFaceID, keyExists )
        ! Here is where the conditional processing begins.
        !
        ! IF this face has already been found THEN we DO nothing
        !
        ! IF this is a new face, we increment the number of faces and add to the key-ring.

        ! Add element to corner-node connectivity list
        IF( .NOT.(keyExists) )THEN ! this is a new face

          nFaces = nFaces + 1
          CALL KeyCabinet(keyRingID) % AddToList( nFaces, globNodeIDs, 4 )

        ENDIF


      ENDDO ! k, Loop over the faces of each element

    ENDDO! iEl, Loop over the elements in the mesh


    DO k = 1, nNodes
      CALL keyCabinet(k) % Trash( ) ! Trash the Facetable
      CALL keyCabinet(k) % Build( ) ! and rebuild a blank cabinet
    ENDDO

    ! Re-allocate space for the mesh Faces

    CALL myHexMesh % Faces % Trash( )
    CALL myHexMesh % Faces % Build( nFaces, N )
    nFaces = 0

    DO iEl = 1, nEls ! Loop over the elements in the mesh

      DO k = 1, 6 ! Loop over the faces of each element

        ! In this first step, we want to identIFy the unique global node ID's for this face
        ! To DO this, we start by gathering the local (to the element) node ID's.
        DO j = 1, 4
          locNodeIDs(j) = myHexMesh % faceMap(j,k) ! starting local node for this Face
        ENDDO

        ! Now, we extract, for this element, the global node ID's for each local node ID
        DO j = 1, 4
          globNodeIDs(j) = myHexMesh % elements % nodeIDs( locNodeIDs(j), iEl )
        ENDDO

        ! Our key cabinet has many key-rings with each key-ring containing a set of notched keys.
        ! Each notched key corresponds to a unique face in the mesh. To enable fast searching
        ! for a unique face, we address our key-rings according to the minimum global node ID
        ! that resides on the current face. IF another face shares this minimum global node ID,
        ! THEN it is possible that the face has already been generated. IF not, our search is
        ! limited only to the key-ring with the same key-ring ID. Here, we grab the key-ring ID.
        keyRingID = MINVAL( globNodeIDs )

        ! Now, we check to see IF a notched key already exists with the same set of global ID's
        ! This is DOne by making a call to "FindDATAForNotches". This routine searches through
        ! the current key ring for the key which has the same global ID's (though not
        ! necessarily in the same order). IF the face has been built already, the face ID is
        ! returned in globFaceID and the LOGICAL, "keyExists", is set to TRUE. Otherwise,
        ! globFaceID is set to zero and "keyExsists" is set to FALSE.
        CALL KeyCabinet(keyRingID) % FindDATAForNotches( globNodeIDS, 4, &
          globFaceID, keyExists )

        ! Here is where the conditional processing begins.
        !
        ! IF this face has already been found THEN we need to determine the face orientation
        ! of the secondary element relative to the primary element.
        !
        ! IF this is a new face, we set the primary element information, and default the
        ! secondary element information.

        ! Add element to corner-node connectivity list
        IF( keyExists )THEN ! this face has already been found

          ! Since this face exists, we need to compare the relative face orientation of the
          ! secondary element to the primary element. .

          CALL myHexMesh % DetermineOrientation( globFaceID, globNodeIDs )

          myHexMesh % faces % elementIDs(2, globFaceID)   = iEl
          myHexMesh % faces % elementSides(2, globFaceID) = k

        ELSE ! This is a new face

          ! First, we store the key-ring information
          nFaces = nFaces + 1
          CALL KeyCabinet(keyRingID) % AddToList( nFaces, globNodeIDs, 4 )

          ! Now, we set the primary element information
          myHexMesh % faces % faceID(nFaces)         = nFaces
          myHexMesh % faces % nodeIDs(1:4,nFaces)    = globNodeIDs
          myHexMesh % faces % elementIDs(1,nFaces)   = iEl
          myHexMesh % faces % elementSides(1,nFaces) = k
          ! Now we default the secondary element information and the swap flag
          myHexMesh % faces % elementIDs(2,nFaces)   = BoundaryFlagStructured
          myHexMesh % faces % elementSides(2,nFaces) = k
          myHexMesh % faces % iStart(nFaces)         = 0
          myHexMesh % faces % iInc(nFaces)           = 1
          myHexMesh % faces % jStart(nFaces)         = 0
          myHexMesh % faces % jInc(nFaces)           = 1
          myHexMesh % faces % swapDimensions(nFaces) = 0

        ENDIF

      ENDDO ! k, Loop over the faces of each element

    ENDDO! iEl, Loop over the elements in the mesh
    DO k = 1, nNodes
      CALL keyCabinet(k) % Trash( ) ! Trash the Facetable
    ENDDO

    DO k = 1, nFaces
      DO j = 0, N
        DO i = 0, N

          IF( i == 0 )THEN
            IF( j == 0 )THEN
              ii = (1-myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % iStart(k)) + &
                (myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % jStart(k))
              jj = (1-myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % jStart(k)) + &
                (myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % iStart(k))
            ELSE
              ii = myHexMesh % faces % swapDimensions(k)*&
                (ii+myHexMesh % faces % jInc(k)) + &
                (1-myHexMesh % faces % swapDimensions(k))*&
                myHexMesh % faces % iStart(k)
              jj = (1-myHexMesh % faces % swapDimensions(k))*&
                (jj+myHexMesh % faces % jInc(k)) +&
                myHexMesh % faces % swapDimensions(k)*&
                myHexMesh % faces % jStart(k)
            ENDIF
          ELSE
            ii = (1-myHexMesh % faces % swapDimensions(k))*&
              (ii + myHexMesh % faces % iInc(k)) +&
              myHexMesh % faces % swapDimensions(k)*ii
            jj = myHexMesh % faces % swapDimensions(k)*&
              (jj+myHexMesh % faces % iInc(k)) + &
              (1-myHexMesh % faces % swapDimensions(k))*jj
          ENDIF

          myHexMesh % faces % iMap(i,j,k) = ii
          myHexMesh % faces % jMap(i,j,k) = jj

        ENDDO
      ENDDO
    ENDDO

#ifdef HAVE_CUDA

    CALL myHexMesh % faces % UpdateDevice( )

#endif

  END SUBROUTINE ConstructFaces_HexMesh
!
  SUBROUTINE ConstructStructuredFaces_HexMesh( myHexMesh, nXElem, nYElem, nZElem )

! Assumes structured mesh
    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh
    INTEGER, INTENT(in) :: nXElem, nYElem, nZElem
    ! LOCAL
    INTEGER ::  e1, e2, s1, s2, nFaces, i, j, k, l, IFace, ii, jj

    nFaces = (nZElem+1)*nXElem*nYElem + (nXElem+1)*nYElem*nZElem + (nYElem+1)*nXElem*nZElem

    ! Re-allocate space for the mesh Faces
    CALL myHexMesh % faces % Trash( )
    CALL myHexMesh % faces % Build( nFaces, myHexMesh % elements % N )

    IFace = 0

    DO k = 1, nZElem
      DO j = 1, nYElem
        DO i = 1, nXElem

          e1 = i + nXElem*( j-1 + nYElem*(k-1) ) ! Primary element ID
          ! Element e1's southern boundary
          s1 = SOUTH
          IF( j==1 )THEN

            IFace = IFace + 1
            e2 = NO_NORMAL_FLOW ! Enforce boundary condition on south boundary of the domain
            s2 = NORTH
            myHexMesh % faces % faceID(IFace) = IFace
            DO l = 1, 4
              myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,SOUTH), e1 )
            ENDDO
            myHexMesh % faces % elementIDs(1,IFace)   = e1
            myHexMesh % faces % elementIDs(2,IFace)   = e2
            myHexMesh % faces % elementSides(1,IFace) = s1
            myHexMesh % faces % elementSides(2,IFace) = s2
            myHexMesh % faces % iStart(IFace)         = 0
            myHexMesh % faces % iInc(IFace)           = 1
            myHexMesh % faces % jStart(IFace)         = 0
            myHexMesh % faces % jInc(IFace)           = 1
            myHexMesh % faces % swapDimensions(IFace) = 0

          ELSE

            IFace = IFace + 1
            e2 = i + nXElem*( j-2 + (nYElem)*(k-1) )
            s2 = NORTH
            myHexMesh % faces % faceID(IFace) = IFace
            DO l = 1, 4
              myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,SOUTH), e1 )
            ENDDO
            myHexMesh % faces % elementIDs(1,IFace)   = e1
            myHexMesh % faces % elementIDs(2,IFace)   = e2
            myHexMesh % faces % elementSides(1,IFace) = s1
            myHexMesh % faces % elementSides(2,IFace) = s2
            myHexMesh % faces % iStart(IFace)         = 0
            myHexMesh % faces % iInc(IFace)           = 1
            myHexMesh % faces % jStart(IFace)         = 0
            myHexMesh % faces % jInc(IFace)           = 1
            myHexMesh % faces % swapDimensions(IFace) = 0

          ENDIF

          ! Element e1's western boundary
          s1 = WEST
          IF( i==1 )THEN

            IFace = IFace + 1
            e2 = NO_NORMAL_FLOW ! Enforce boundary condition on west domain boundary
            s2 = EAST
            myHexMesh % faces % faceID(IFace) = IFace
            DO l = 1, 4
              myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,WEST), e1 )
            ENDDO
            myHexMesh % faces % elementIDs(1,IFace)   = e1
            myHexMesh % faces % elementIDs(2,IFace)   = e2
            myHexMesh % faces % elementSides(1,IFace) = s1
            myHexMesh % faces % elementSides(2,IFace) = s2
            myHexMesh % faces % iStart(IFace)         = 0
            myHexMesh % faces % iInc(IFace)           = 1
            myHexMesh % faces % jStart(IFace)         = 0
            myHexMesh % faces % jInc(IFace)           = 1
            myHexMesh % faces % swapDimensions(IFace) = 0

          ELSE
            IFace = IFace + 1
            e2 = i-1 + nXElem*( j-1 + (nYElem)*(k-1) )
            s2 = EAST
            myHexMesh % faces % faceID(IFace) = IFace
            DO l = 1, 4
              myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,WEST), e1 )
            ENDDO
            myHexMesh % faces % elementIDs(1,IFace)   = e1
            myHexMesh % faces % elementIDs(2,IFace)   = e2
            myHexMesh % faces % elementSides(1,IFace) = s1
            myHexMesh % faces % elementSides(2,IFace) = s2
            myHexMesh % faces % iStart(IFace)         = 0
            myHexMesh % faces % iInc(IFace)           = 1
            myHexMesh % faces % jStart(IFace)         = 0
            myHexMesh % faces % jInc(IFace)           = 1
            myHexMesh % faces % swapDimensions(IFace) = 0

          ENDIF

          ! Element e1's bottom boundary
          s1 = BOTTOM
          IF( k==1 )THEN

            IFace = IFace + 1
            e2 = NO_NORMAL_FLOW ! Enforce boundary condition on bottom domain boundary
            s2 = TOP
            myHexMesh % faces % faceID(IFace) = IFace
            DO l = 1, 4
              myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,BOTTOM), e1 )
            ENDDO
            myHexMesh % faces % elementIDs(1,IFace)   = e1
            myHexMesh % faces % elementIDs(2,IFace)   = e2
            myHexMesh % faces % elementSides(1,IFace) = s1
            myHexMesh % faces % elementSides(2,IFace) = s2
            myHexMesh % faces % iStart(IFace)         = 0
            myHexMesh % faces % iInc(IFace)           = 1
            myHexMesh % faces % jStart(IFace)         = 0
            myHexMesh % faces % jInc(IFace)           = 1
            myHexMesh % faces % swapDimensions(IFace) = 0

          ELSE

            IFace = IFace + 1
            e2 = i + nXElem*( j-1 + (nYElem)*(k-2) )
            s2 = TOP
            myHexMesh % faces % faceID(IFace) = IFace
            DO l = 1, 4
              myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,BOTTOM), e1 )
            ENDDO
            myHexMesh % faces % elementIDs(1,IFace)   = e1
            myHexMesh % faces % elementIDs(2,IFace)   = e2
            myHexMesh % faces % elementSides(1,IFace) = s1
            myHexMesh % faces % elementSides(2,IFace) = s2
            myHexMesh % faces % iStart(IFace)         = 0
            myHexMesh % faces % iInc(IFace)           = 1
            myHexMesh % faces % jStart(IFace)         = 0
            myHexMesh % faces % jInc(IFace)           = 1
            myHexMesh % faces % swapDimensions(IFace) = 0

          ENDIF

        ENDDO ! i

        e1 = nXElem + nXElem*( j-1 + (nYElem)*(k-1) )
        s1 = EAST
        IFace = IFace + 1
        e2 = NO_NORMAL_FLOW ! Enforce boundary condition on east domain boundary
        s2 = WEST
        myHexMesh % faces % faceID(IFace) = IFace
        DO l = 1, 4
          myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,EAST), e1 )
        ENDDO
        myHexMesh % faces % elementIDs(1,IFace)   = e1
        myHexMesh % faces % elementIDs(2,IFace)   = e2
        myHexMesh % faces % elementSides(1,IFace) = s1
        myHexMesh % faces % elementSides(2,IFace) = s2
        myHexMesh % faces % iStart(IFace)         = 0
        myHexMesh % faces % iInc(IFace)           = 1
        myHexMesh % faces % jStart(IFace)         = 0
        myHexMesh % faces % jInc(IFace)           = 1
        myHexMesh % faces % swapDimensions(IFace) = 0

      ENDDO ! j

      DO i = 1, nXElem

        IFace = IFace + 1
        e1 = i + nXElem*( nYElem-1 + (nYElem)*(k-1) )
        s1 = NORTH
        e2 = NO_NORMAL_FLOW ! Enforce boundary condition on domain north boundary
        s2 = SOUTH
        myHexMesh % faces % faceID(IFace) = IFace
        DO l = 1, 4
          myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,NORTH), e1 )
        ENDDO
        myHexMesh % faces % elementIDs(1,IFace)    = e1
        myHexMesh % faces % elementIDs(2,IFace)    = e2
        myHexMesh % faces % elementSides(1,IFace)  = s1
        myHexMesh % faces % elementSides(2,IFace)  = s2
        myHexMesh % faces % iStart(IFace)          = 0
        myHexMesh % faces % iInc(IFace)            = 1
        myHexMesh % faces % jStart(IFace)          = 0
        myHexMesh % faces % jInc(IFace)            = 1
        myHexMesh % faces % swapDimensions(IFace)  = 0

      ENDDO

    ENDDO ! k

    DO j = 1, nYElem
      DO i = 1, nXElem

        e1 = i + nXElem*( j-1 + (nYElem)*(nZElem-1) ) ! Primary element ID
        IFace = IFace + 1
        e2 = NO_NORMAL_FLOW ! Enforce boundary condition on domain top
        s1 = TOP
        s2 = s1
        myHexMesh % faces % faceID(IFace) = IFace
        DO l = 1, 4
          myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,TOP), e1 )
        ENDDO
        myHexMesh % faces % elementIDs(1,IFace)   = e1
        myHexMesh % faces % elementIDs(2,IFace)   = e2
        myHexMesh % faces % elementSides(1,IFace) = s1
        myHexMesh % faces % elementSides(2,IFace) = s2
        myHexMesh % faces % iStart(IFace)         = 0
        myHexMesh % faces % iInc(IFace)           = 1
        myHexMesh % faces % jStart(IFace)         = 0
        myHexMesh % faces % jInc(IFace)           = 1
        myHexMesh % faces % swapDimensions(IFace) = 0

      ENDDO
    ENDDO

    DO k = 1, myHexMesh % faces % nFaces
      DO j = 0, myHexMesh % elements % N
        DO i = 0, myHexMesh % elements % N

          IF( i == 0 )THEN
            IF( j == 0 )THEN
              ii = (1-myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % iStart(k)) + &
                (myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % jStart(k))
              jj = (1-myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % jStart(k)) + &
                (myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % iStart(k))
            ELSE
              ii = myHexMesh % faces % swapDimensions(k)*&
                (ii+myHexMesh % faces % jInc(k)) + &
                (1-myHexMesh % faces % swapDimensions(k))*&
                myHexMesh % faces % iStart(k)
              jj = (1-myHexMesh % faces % swapDimensions(k))*&
                (jj+myHexMesh % faces % jInc(k)) +&
                myHexMesh % faces % swapDimensions(k)*&
                myHexMesh % faces % jStart(k)
            ENDIF
          ELSE
            ii = (1-myHexMesh % faces % swapDimensions(k))*&
              (ii + myHexMesh % faces % iInc(k)) +&
              myHexMesh % faces % swapDimensions(k)*ii
            jj = myHexMesh % faces % swapDimensions(k)*&
              (jj+myHexMesh % faces % iInc(k)) + &
              (1-myHexMesh % faces % swapDimensions(k))*jj
          ENDIF

          myHexMesh % faces % iMap(i,j,k) = ii
          myHexMesh % faces % jMap(i,j,k) = jj

        ENDDO
      ENDDO
    ENDDO

#ifdef HAVE_CUDA

    CALL myHexMesh % faces % UpdateDevice( )

#endif

  END SUBROUTINE ConstructStructuredFaces_HexMesh
!  
  SUBROUTINE ConstructDoublyPeriodicFaces_HexMesh( myHexMesh, nXElem, nYElem, nZElem )

! Assumes structured mesh
    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh
    INTEGER, INTENT(in) :: nXElem, nYElem, nZElem
    ! LOCAL
    INTEGER ::  e1, e2, s1, s2, nFaces, i, j, k, l, IFace, ii, jj

    nFaces = (nZElem+1)*nXElem*nYElem + (nXElem+1)*nYElem*nZElem + (nYElem+1)*nXElem*nZElem

    ! Re-allocate space for the mesh Faces
    CALL myHexMesh % faces % Trash( )
    CALL myHexMesh % faces % Build( nFaces, myHexMesh % elements % N )

    IFace = 0

    DO k = 1, nZElem
      DO j = 1, nYElem
        DO i = 1, nXElem

          e1 = i + nXElem*( j-1 + nYElem*(k-1) ) ! Primary element ID
          ! Element e1's southern boundary
          s1 = SOUTH
          IF( j==1 )THEN

            IFace = IFace + 1
            e2 = i + nXElem*( nYElem-1 + (nYElem)*(k-1) )! Enforce periodicity with the "northern" most element
            s2 = NORTH
            myHexMesh % faces % faceID(IFace) = IFace
            DO l = 1, 4
              myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,SOUTH), e1 )
            ENDDO
            myHexMesh % faces % elementIDs(1,IFace)   = e1
            myHexMesh % faces % elementIDs(2,IFace)   = e2
            myHexMesh % faces % elementSides(1,IFace) = s1
            myHexMesh % faces % elementSides(2,IFace) = s2
            myHexMesh % faces % iStart(IFace)         = 0
            myHexMesh % faces % iInc(IFace)           = 1
            myHexMesh % faces % jStart(IFace)         = 0
            myHexMesh % faces % jInc(IFace)           = 1
            myHexMesh % faces % swapDimensions(IFace) = 0

          ELSE

            IFace = IFace + 1
            e2 = i + nXElem*( j-2 + (nYElem)*(k-1) )
            s2 = NORTH
            myHexMesh % faces % faceID(IFace) = IFace
            DO l = 1, 4
              myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,SOUTH), e1 )
            ENDDO
            myHexMesh % faces % elementIDs(1,IFace)   = e1
            myHexMesh % faces % elementIDs(2,IFace)   = e2
            myHexMesh % faces % elementSides(1,IFace) = s1
            myHexMesh % faces % elementSides(2,IFace) = s2
            myHexMesh % faces % iStart(IFace)         = 0
            myHexMesh % faces % iInc(IFace)           = 1
            myHexMesh % faces % jStart(IFace)         = 0
            myHexMesh % faces % jInc(IFace)           = 1
            myHexMesh % faces % swapDimensions(IFace) = 0

          ENDIF

          ! Element e1's western boundary
          s1 = WEST
          IF( i==1 )THEN

            IFace = IFace + 1
            e2 = nXElem + nXElem*( j-1 + (nYElem)*(k-1) ) ! Enforce periodicity with the "eastern" most element
            s2 = EAST
            myHexMesh % faces % faceID(IFace) = IFace
            DO l = 1, 4
              myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,WEST), e1 )
            ENDDO
            myHexMesh % faces % elementIDs(1,IFace)   = e1
            myHexMesh % faces % elementIDs(2,IFace)   = e2
            myHexMesh % faces % elementSides(1,IFace) = s1
            myHexMesh % faces % elementSides(2,IFace) = s2
            myHexMesh % faces % iStart(IFace)         = 0
            myHexMesh % faces % iInc(IFace)           = 1
            myHexMesh % faces % jStart(IFace)         = 0
            myHexMesh % faces % jInc(IFace)           = 1
            myHexMesh % faces % swapDimensions(IFace) = 0

          ELSE
            IFace = IFace + 1
            e2 = i-1 + nXElem*( j-1 + (nYElem)*(k-1) )
            s2 = EAST
            myHexMesh % faces % faceID(IFace) = IFace
            DO l = 1, 4
              myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,WEST), e1 )
            ENDDO
            myHexMesh % faces % elementIDs(1,IFace)   = e1
            myHexMesh % faces % elementIDs(2,IFace)   = e2
            myHexMesh % faces % elementSides(1,IFace) = s1
            myHexMesh % faces % elementSides(2,IFace) = s2
            myHexMesh % faces % iStart(IFace)         = 0
            myHexMesh % faces % iInc(IFace)           = 1
            myHexMesh % faces % jStart(IFace)         = 0
            myHexMesh % faces % jInc(IFace)           = 1
            myHexMesh % faces % swapDimensions(IFace) = 0

          ENDIF

          ! Element e1's bottom boundary
          s1 = BOTTOM
          IF( k==1 )THEN

            IFace = IFace + 1
            e2 = NO_NORMAL_FLOW
            s2 = TOP
            myHexMesh % faces % faceID(IFace) = IFace
            DO l = 1, 4
              myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,BOTTOM), e1 )
            ENDDO
            myHexMesh % faces % elementIDs(1,IFace)   = e1
            myHexMesh % faces % elementIDs(2,IFace)   = e2
            myHexMesh % faces % elementSides(1,IFace) = s1
            myHexMesh % faces % elementSides(2,IFace) = s2
            myHexMesh % faces % iStart(IFace)         = 0
            myHexMesh % faces % iInc(IFace)           = 1
            myHexMesh % faces % jStart(IFace)         = 0
            myHexMesh % faces % jInc(IFace)           = 1
            myHexMesh % faces % swapDimensions(IFace) = 0

          ELSE

            IFace = IFace + 1
            e2 = i + nXElem*( j-1 + (nYElem)*(k-2) )
            s2 = TOP
            myHexMesh % faces % faceID(IFace) = IFace
            DO l = 1, 4
              myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,BOTTOM), e1 )
            ENDDO
            myHexMesh % faces % elementIDs(1,IFace)   = e1
            myHexMesh % faces % elementIDs(2,IFace)   = e2
            myHexMesh % faces % elementSides(1,IFace) = s1
            myHexMesh % faces % elementSides(2,IFace) = s2
            myHexMesh % faces % iStart(IFace)         = 0
            myHexMesh % faces % iInc(IFace)           = 1
            myHexMesh % faces % jStart(IFace)         = 0
            myHexMesh % faces % jInc(IFace)           = 1
            myHexMesh % faces % swapDimensions(IFace) = 0

          ENDIF

        ENDDO ! i

        e1 = nXElem + nXElem*( j-1 + (nYElem)*(k-1) )
        s1 = EAST
        IFace = IFace + 1
        e2 = 1 + nXElem*( j-1 + (nYElem)*(k-1) ) ! Enforce periodicity with the "western" most element
        s2 = WEST
        myHexMesh % faces % faceID(IFace) = IFace
        DO l = 1, 4
          myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,EAST), e1 )
        ENDDO
        myHexMesh % faces % elementIDs(1,IFace)   = e1
        myHexMesh % faces % elementIDs(2,IFace)   = e2
        myHexMesh % faces % elementSides(1,IFace) = s1
        myHexMesh % faces % elementSides(2,IFace) = s2
        myHexMesh % faces % iStart(IFace)         = 0
        myHexMesh % faces % iInc(IFace)           = 1
        myHexMesh % faces % jStart(IFace)         = 0
        myHexMesh % faces % jInc(IFace)           = 1
        myHexMesh % faces % swapDimensions(IFace) = 0

      ENDDO ! j

      DO i = 1, nXElem

        IFace = IFace + 1
        e1 = i + nXElem*( nYElem-1 + (nYElem)*(k-1) )
        s1 = NORTH
        e2 = i + nXElem*( nYElem*(k-1) )
        s2 = SOUTH
        myHexMesh % faces % faceID(IFace) = IFace
        DO l = 1, 4
          myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,NORTH), e1 )
        ENDDO
        myHexMesh % faces % elementIDs(1,IFace)    = e1
        myHexMesh % faces % elementIDs(2,IFace)    = e2
        myHexMesh % faces % elementSides(1,IFace)  = s1
        myHexMesh % faces % elementSides(2,IFace)  = s2
        myHexMesh % faces % iStart(IFace)          = 0
        myHexMesh % faces % iInc(IFace)            = 1
        myHexMesh % faces % jStart(IFace)          = 0
        myHexMesh % faces % jInc(IFace)            = 1
        myHexMesh % faces % swapDimensions(IFace)  = 0

      ENDDO

    ENDDO ! k

    DO j = 1, nYElem
      DO i = 1, nXElem

        e1 = i + nXElem*( j-1 + (nYElem)*(nZElem-1) ) ! Primary element ID
        IFace = IFace + 1
        e2 = PRESCRIBED
        s1 = TOP
        s2 = s1
        myHexMesh % faces % faceID(IFace) = IFace
        DO l = 1, 4
          myHexMesh % faces % nodeIDs(l,IFace) = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(l,TOP), e1 )
        ENDDO
        myHexMesh % faces % elementIDs(1,IFace)   = e1
        myHexMesh % faces % elementIDs(2,IFace)   = e2
        myHexMesh % faces % elementSides(1,IFace) = s1
        myHexMesh % faces % elementSides(2,IFace) = s2
        myHexMesh % faces % iStart(IFace)         = 0
        myHexMesh % faces % iInc(IFace)           = 1
        myHexMesh % faces % jStart(IFace)         = 0
        myHexMesh % faces % jInc(IFace)           = 1
        myHexMesh % faces % swapDimensions(IFace) = 0

      ENDDO
    ENDDO

    DO k = 1, myHexMesh % faces % nFaces
      DO j = 0, myHexMesh % elements % N
        DO i = 0, myHexMesh % elements % N

          IF( i == 0 )THEN
            IF( j == 0 )THEN
              ii = (1-myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % iStart(k)) + &
                (myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % jStart(k))
              jj = (1-myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % jStart(k)) + &
                (myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % iStart(k))
            ELSE
              ii = myHexMesh % faces % swapDimensions(k)*&
                (ii+myHexMesh % faces % jInc(k)) + &
                (1-myHexMesh % faces % swapDimensions(k))*&
                myHexMesh % faces % iStart(k)
              jj = (1-myHexMesh % faces % swapDimensions(k))*&
                (jj+myHexMesh % faces % jInc(k)) +&
                myHexMesh % faces % swapDimensions(k)*&
                myHexMesh % faces % jStart(k)
            ENDIF
          ELSE
            ii = (1-myHexMesh % faces % swapDimensions(k))*&
              (ii + myHexMesh % faces % iInc(k)) +&
              myHexMesh % faces % swapDimensions(k)*ii
            jj = myHexMesh % faces % swapDimensions(k)*&
              (jj+myHexMesh % faces % iInc(k)) + &
              (1-myHexMesh % faces % swapDimensions(k))*jj
          ENDIF

          myHexMesh % faces % iMap(i,j,k) = ii
          myHexMesh % faces % jMap(i,j,k) = jj

        ENDDO
      ENDDO
    ENDDO

#ifdef HAVE_CUDA

    CALL myHexMesh % faces % UpdateDevice( )

#endif

  END SUBROUTINE ConstructDoublyPeriodicFaces_HexMesh
!
!> \addtogroup HexMesh_Class
!! @{
! ================================================================================================ !
! S/R DetermineOrientation
!
!> \fn DetermineOrientation_HexMesh
!! Given a face-ID and an array of global node ID's, this routine determines the relative orientation
!! of the face to a face represented by the incoming global node ID's.
!!
!!  This support routine takes as input the mesh, a face ID number, and a list of node IDs.
!!  This routine assumes that the primary element information for the given face ID has been
!!  filled upon issuing a call to this routine. IF the primary element shares this face with
!!  another element, this routine is called and the secondary element node ID's are passed in.
!!  The order of the secondary node ID's relative to the order of the primary node ID's determines
!!  the orientation of the secondary element relative to the primary element. We need to know
!!  IF the roles of the computational coordinates are flipped, and (simultaneously) IF the
!!  incrementing of the computational coordinates are reversed (or not).
!!  This routine determines the orientation by determining how many times the ordered secondary
!!  node ID's need to be shIFted in order to match the ordered primary node ID's.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(HexMesh) :: this <BR>
!! <B>INTEGER</B>       :: faceID <BR>
!! <B>INTEGER</B>       :: secondaryNodes(1:4) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % DetermineOrientation( faceID, secondaryNodes ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myHexMesh <td> HexMesh <td>
!!                         On <B>input</B>, HexMesh with myHexMesh % faces(faceID) primary element
!!                         information filled in, <BR>
!!                         On <B>output</B>, secondary element information is filled in for this face
!!   <tr> <td> in <th> faceID <td> INTEGER <td> ID for the face in question
!!   <tr> <td> in <th> secondaryNodes(1:4) <td> INTEGER <td> Global node ID's associated with the
!!                                                           secondary element.
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE DetermineOrientation_HexMesh( myHexMesh, faceID, secondaryNodes )

    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh
    INTEGER, INTENT(in)             :: faceID
    INTEGER, INTENT(in)             :: secondaryNodes(1:4)
    ! Local
    INTEGER :: primaryNodes(1:4)
    INTEGER :: nShIFts, i, N
    LOGICAL :: theyMatch

    primaryNodes = myHexMesh % faces % nodeIDs(1:4,faceID)
    N = myHexMesh % elements % N
    nShIFts = 0
    theyMatch = .FALSE.

    DO i = 1, 4

      ! First, we compare the primary and secondary nodes. This routine returns a zero
      ! IF the arrays match, and a one IF they DO not match.
      theyMatch = CompareArray( primaryNodes, secondaryNodes, 4 )

      IF( theyMatch )THEN
        EXIT
      ELSE
        nShIFts = nShIFts + 1
        CALL ForwardShIFt( primaryNodes, 4 )
      ENDIF

    ENDDO

    IF( theyMatch )THEN

      SELECT CASE ( nShIFts )

       CASE (0)
        myHexMesh % faces % iStart(faceID)          = 0
        myHexMesh % faces % iInc(faceID)            = 1
        myHexMesh % faces % jStart(faceID)          = 0
        myHexMesh % faces % jInc(faceID)            = 1
        myHexMesh % faces % swapDimensions(faceID)  = 0
       CASE (1)
        myHexMesh % faces % iStart(faceID)          = 0
        myHexMesh % faces % iInc(faceID)            = 1
        myHexMesh % faces % jStart(faceID)          = N
        myHexMesh % faces % jInc(faceID)            = -1
        myHexMesh % faces % swapDimensions(faceID)  = 1
       CASE (2)
        myHexMesh % faces % iStart(faceID)          = N
        myHexMesh % faces % iInc(faceID)            = -1
        myHexMesh % faces % jStart(faceID)          = N
        myHexMesh % faces % jInc(faceID)            = -1
        myHexMesh % faces % swapDimensions(faceID)  = 0
       CASE (3)
        myHexMesh % faces % iStart(faceID)          = N
        myHexMesh % faces % iInc(faceID)            = -1
        myHexMesh % faces % jStart(faceID)          = 0
        myHexMesh % faces % jInc(faceID)            = 1
        myHexMesh % faces % swapDimensions(faceID)  = 1
       CASE DEFAULT

      END SELECT

    ENDIF


  END SUBROUTINE DetermineOrientation_HexMesh
!
!> \addtogroup HexMesh_Class
!! @{
! ================================================================================================ !
! S/R ConstructElementNeighbors
!
!> \fn sConstructElementNeighbors_HexMesh
!! USEs the edge-to-element connectivity to construct the element neighbors attribute.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(myHexMesh) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ConstructElementNeighbors( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myHexMesh <td> HexMesh <td>
!!                         On <B>input</B>, the edge-to-element information has been constructed, <BR>
!!                         On <B>output</B>, the element neighbors information has been filled in.
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE ConstructElementNeighbors_HexMesh( myHexMesh )

    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh
    ! LOCAL
    INTEGER :: e(1:2), s(1:2), IFace

    DO IFace = 1, myHexMesh % faces % nFaces
      e = myHexMesh % faces % elementIDs(1:2,IFace)
      s = myHexMesh % faces % elementSides(1:2,IFace)
      IF( e(2) > 0 )THEN
        myHexMesh % elements % neighbors(s(1), e(1))      = e(2)
        myHexMesh % elements % neighbors(ABS(s(2)), e(2)) = e(1)
      ENDIF
    ENDDO


  END SUBROUTINE ConstructElementNeighbors_HexMesh
!
!> \addtogroup HexMesh_Class
!! @{
! ================================================================================================ !
! S/R ScaleTheMesh_HexMesh
!
!> \fn ScaleTheMesh
!! Scales the element geometry and corner node positions by the provided x-scale, y-scale, and z-scale.
!!
!! This routine depend on <BR>
!!   Module \ref MappedGeometry_2D_Class, S/R ScaleGeometry_MappedGeometry_2D <BR>
!!   Module \ref Node_Class, S/R ScalePosition_Node
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(HexMesh)    :: this <BR>
!! <B>TYPE</B>(Lagrange) :: interp <BR>
!! <B>REAL</B>(prec)        :: xScale, yScale, zScale <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ScaleGeometry( interp, xScale, yScale, zScale ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myHexMesh <td> HexMesh <td>
!!   <tr> <td> in <th> interp <td> Lagrange <td>
!!   <tr> <td> in <th> xScale <td> REAL(prec) <td>
!!   <tr> <td> in <th> yScale <td> REAL(prec) <td>
!!   <tr> <td> in <th> zScale <td> REAL(prec) <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE ScaleTheMesh_HexMesh( myHexMesh, interp, xScale, yScale, zScale  )

    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh
    TYPE( Lagrange ), INTENT(in)    :: interp
    REAL(prec), INTENT(in)          :: xScale, yScale, zScale
#ifdef HAVE_CUDA
    INTEGER :: istat
#endif

    CALL myHexMesh % elements % ScaleGeometry( interp, xScale, yScale, zScale )
    CALL myHexMesh % nodes % ScalePosition( xScale, yScale, zScale )

#ifdef HAVE_CUDA

    CALL myHexMesh % elements % UpdateDevice( )
    CALL myHexMesh % nodes % UpdateDevice( )
    istat = cudaDeviceSynchronize( )

#endif

  END SUBROUTINE ScaleTheMesh_HexMesh
!
!> \addtogroup HexMesh_Class
!! @{
! ================================================================================================ !
! S/R ConstructStructuredMesh
!
!> \fn ConstructStructuredMesh_HexMesh
!! Constructs a "structured" spectral element mesh with nXelem-by-nYelem elements.
!!
!! This routine builds a mesh with nXelem elements in the x-direction,  nYelem elements in the
!! y-direction, and nZelem elements in the z-direction. The mesh nodes are between [0,1]x[0,1]x[0,1].
!!  After construction, the USEr can call "ScaleTheMesh" to change the physical extents of the DOmain.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(HexMesh)    :: this <BR>
!! <B>TYPE</B>(Lagrange) :: interp <BR>
!! <B>INTEGER</B>           :: nXElem, nYElem, nZElem <BR>
!!         .... <BR>
!!     <B>CALL</B> this % RoutineName( interp, nXElem, nYElem, nZElem ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myHexMesh <td> HexMesh <td> On output, CONTAINS the "structured" mesh
!!                                                       in the unstructured spectral element
!!                                                       mesh format.
!!   <tr> <td> in <th> interp <td> Lagrange <td> 2-D interpolant that stores the computational
!!                                                  quadrature mesh for each element.
!!   <tr> <td> in <th> nXElem <td> INTEGER <td> The number of desired elements in the x-direction
!!   <tr> <td> in <th> nYElem <td> INTEGER <td> The number of desired elements in the y-direction
!!   <tr> <td> in <th> nZElem <td> INTEGER <td> The number of desired elements in the z-direction
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE ConstructStructuredMesh_HexMesh( myHexMesh, interp, nXelem, nYelem, nZelem, geomparser, DoublyPeriodic  )

    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout)  :: myHexMesh
    TYPE( Lagrange ), INTENT(in)     :: interp
    INTEGER, INTENT(in)              :: nXelem, nYelem, nZelem
    TYPE( Geom_EquationParser ), INTENT(in) :: geomparser
    LOGICAL, INTENT(in)              :: DoublyPeriodic
    ! LOGICAL
    TYPE( Surfaces ) :: boundSurfs
    REAL(prec) :: x, y, z, zb, zi, zu, zip1, dxElem, dyElem, dzElem
    REAL(prec) :: x1(1:3), x2(1:3), x3(1:3), x4(1:3)
    REAL(prec) :: c1(1:3), c2(1:3), xE(1:3)
    REAL(prec), ALLOCATABLE :: xc(:,:,:,:), s(:), weights(:)

    INTEGER :: nNodes, nElements, nFaces, gPolyDeg, nSurf
    INTEGER :: nodes(1:8)
    INTEGER :: n1, n2, n3, n4
    INTEGER :: iNode, iEl, iSide, iX, iY, iZ, i, j, iSurf

    dxElem = 1.0_prec/nXElem
    dyElem = 1.0_prec/nYElem
    dzElem = 1.0_prec/nZElem

    ! ** "Hard-wired" values for a structured mesh with no holes ** !
    nNodes    = (nXElem+1)*(nYElem+1)*(nZElem+1)
    nElements = (nXElem)*(nYElem)*(nZElem)
    nFaces    = (nXElem)*(nYElem)*(nZElem+1) + (nXElem)*(nZElem)*(nYElem+1) + (nYElem)*(nZElem)*(nXElem+1)
    gPolyDeg  = interp % N
    nSurf     = 6*nElements
    ! ************************************************************************* !


    ! Generate the Legendre-Gauss Lobatto points of order gPolyDeg
    ! These are the points USEd to define the parametric
    ! curves for the element boundaries

    ALLOCATE( s(0:gPolyDeg), xc(0:gPolyDeg,0:gPolyDeg,1:3,1:nSurf), weights(0:gpolyDeg) )
    CALL LegendreQuadrature( gPolyDeg, s, weights, GAUSS_LOBATTO )

    ! ---- Build the mesh (empty) ---- !
    CALL myHexMesh % Build( nNodes, nElements, nFaces, interp % N )


    ! ---- Read in the corner nodes ---- !
    DO iZ = 1, nZElem + 1

      zi = dZElem*(REAL(iZ-1,prec))

      DO iY = 1, nYElem+1

        y = dYElem*(REAL(iY-1,prec))

        DO iX = 1, nXElem+1

          iNode = iX + (iY-1)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)
          x = dXElem*(REAL(iX-1,prec))

          xE(1:3) = (/ x, y, 0.0_prec /)
          zb = geomparser % topography % Evaluate( xE )

          z = zb*(1.0_prec-zi) + 1.0_prec*zi

          myHexMesh % nodes % x(1:3,iNode)  = (/ x, y, z /)
          myHexMesh % nodes % nodeID(iNode) = iNode

        ENDDO
      ENDDO
    ENDDO

    ! DO the element information
    xc = 0.0_prec
    CALL boundSurfs % Build( s, gPolyDeg, nSurf )

    DO iZ = 1, nZElem

      zi   = dZElem*(REAL(iZ-1,prec))
      zip1 = dZElem*(REAL(iZ,prec))

      DO iY = 1, nYElem
        DO iX = 1, nXElem

          iEl = iX + (iY-1)*(nXElem) + (iZ-1)*(nXElem)*(nYElem)
          ! Calculate the global node IDs for this element.
          nodes(1) = iX + (iY-1)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)    ! Southwest
          nodes(2) = iX + 1 + (iY-1)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)! SouthEast
          nodes(3) = iX + 1 + (iY)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)  ! NorthEast
          nodes(4) = iX + (iY)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)      ! NorthWest

          nodes(5) = iX + (iY-1)*(nXElem+1) + (iZ)*(nXElem+1)*(nYElem+1)      ! Southwest
          nodes(6) = iX + 1 + (iY-1)*(nXElem+1) + (iZ)*(nXElem+1)*(nYElem+1)  ! SouthEast
          nodes(7) = iX + 1 + (iY)*(nXElem+1) + (iZ)*(nXElem+1)*(nYElem+1)    ! NorthEast
          nodes(8) = iX + (iY)*(nXElem+1) + (iZ)*(nXElem+1)*(nYElem+1)        ! NorthWest

          myHexMesh % elements % nodeIDs(1:8,iEl) = nodes
          myHexMesh % elements % elementID(iEl) = iEl

            DO iSide = 1, 6 ! Loop over the sides of the quads

              iSurf = iSide + 6*(iEl-1)
              ! To build the current face, we construct a plane that passes through
              ! the four corner nodes. Here, we grab the global node ID's for the four
              ! corner nodes.
              n1 = nodes( myHexMesh % faceMap(1,iSide) )
              n2 = nodes( myHexMesh % faceMap(2,iSide) )
              n3 = nodes( myHexMesh % faceMap(3,iSide) )
              n4 = nodes( myHexMesh % faceMap(4,iSide) )

              x1(1) = myHexMesh % nodes % x(1,n1)
              x1(2) = myHexMesh % nodes % x(2,n1)
              x1(3) = myHexMesh % nodes % x(3,n1)

              x2(1) = myHexMesh % nodes % x(1,n2)
              x2(2) = myHexMesh % nodes % x(2,n2)
              x2(3) = myHexMesh % nodes % x(3,n2)

              x3(1) = myHexMesh % nodes % x(1,n3)
              x3(2) = myHexMesh % nodes % x(2,n3)
              x3(3) = myHexMesh % nodes % x(3,n3)

              x4(1) = myHexMesh % nodes % x(1,n4)
              x4(2) = myHexMesh % nodes % x(2,n4)
              x4(3) = myHexMesh % nodes % x(3,n4)

              DO j = 0, gPolyDeg
                DO i = 0, gPolyDeg
                  ! Transfinite inerpolation with linear blending is USEd to construct the face
                  c1 = ( 0.5_prec*(x2-x1)*(1.0_prec+s(i)) + x1 )
                  c2 = ( 0.5_prec*(x3-x4)*(1.0_prec+s(i)) + x4 )
                  xc(i,j,1:3,iSurf) = 0.5_prec*(c2-c1)*(1.0_prec+s(j)) + c1

                ENDDO
              ENDDO

                IF( iSide == BOTTOM )THEN
  
                  zu = zi
                  DO j = 0, gPolyDeg
                    DO i = 0, gPolyDeg
  
                      xE(1:3) = (/ xc(i,j,1,iSurf), xc(i,j,2,iSurf), 0.0_prec /)
                      zb = geomparser % topography % Evaluate( xE )
                      xc(i,j,3,iSurf) = zb*(1.0_prec-zu) + 1.0_prec*zu
  
                    ENDDO
                  ENDDO
  
                ELSEIF( iSide == TOP )THEN
  
                  zu = zip1
                  DO j = 0, gPolyDeg
                    DO i = 0, gPolyDeg
  
                      xE(1:3) = (/ xc(i,j,1,iSurf), xc(i,j,2,iSurf), 0.0_prec /)
                      zb = geomparser % topography % Evaluate( xE )
                      xc(i,j,3,iSurf) = zb*(1.0_prec-zu) + 1.0_prec*zu
  
                    ENDDO
                  ENDDO
  
  
  
                ELSE
  
                  DO j = 0, gPolyDeg
                    DO i = 0, gPolyDeg
  
                      zu = 0.5_prec*( zi*( s(j) - 1.0_prec ) + zip1*( s(j) + 1.0_prec ) ) 
                      xE(1:3) = (/ xc(i,j,1,iSurf), xc(i,j,2,iSurf), 0.0_prec /)
                      zb = geomparser % topography % Evaluate( xE )
                      xc(i,j,3,iSurf) = zb*(1.0_prec-zu) + 1.0_prec*zu
  
                    ENDDO
                  ENDDO
  
                ENDIF
  

            ENDDO

        ENDDO
      ENDDO
    ENDDO ! iEl, cycle over the elements

    CALL boundSurfs % Set_Surfaces( xc )
    CALL myHexMesh % elements % GenerateMesh( interp, boundSurfs )
    CALL myHexMesh % elements % GenerateMetrics( interp )

    IF( DoublyPeriodic )THEN
      CALL myHexMesh % ConstructDoublyPeriodicFaces( nXElem, nYElem, nZElem )
    ELSE
      CALL myHexMesh % ConstructStructuredFaces( nXElem, nYElem, nZElem )
    ENDIF


    CALL myHexMesh % ConstructElementNeighbors( )

    ! Clear up memory
    DEALLOCATE( s, xc, weights )

    CALL boundSurfs % Trash( )

#ifdef HAVE_CUDA
    CALL myHexMesh % UpdateDevice( )
#endif

  END SUBROUTINE ConstructStructuredMesh_HexMesh

  SUBROUTINE PartitionStructuredElementsAndNodes( myMesh, params, partitions, nElPerProc, globalToLocal, nodeLogic, nNodePerProc, globalToLocalNode, nProc )
    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(in)        :: myMesh
    TYPE( ModelParameters ), INTENT(in) :: params
    INTEGER, INTENT(in)                 :: nProc
    INTEGER, INTENT(out)                :: partitions(1:myMesh % elements % nElements)
    INTEGER, INTENT(out)                :: nElPerProc(0:nProc-1)
    INTEGER, INTENT(out)                :: globalToLocal(1:myMesh % elements % nElements,1:2)
    INTEGER, INTENT(out)                :: nodeLogic(1:myMesh % nodes % nNodes,0:nProc-1)
    INTEGER, INTENT(out)                :: nNodePerProc(0:nProc-1)
    INTEGER, INTENT(out)                :: globalToLocalNode(1:myMesh % nodes % nNodes,0:nProc-1)
   
    ! Local
!    INTEGER :: nxp, nyp, nzp
    INTEGER :: iPz, iPy, iPx
    INTEGER :: iZp, iYp, iXp
    INTEGER :: iZ, iY, iX
    INTEGER :: iNode, nID
    INTEGER :: iEl, procID
    INTEGER :: nxp(1:params % nProcX)
    INTEGER :: nyp(1:params % nProcY)
    INTEGER :: nzp(1:params % nProcZ)
    

      partitions        = 0
      nElPerProc        = 0
      globalToLocal     = 0
      nodeLogic         = 0
      nNodePerProc      = 0
      globalToLocalNode = 0

      IF( nProc > 1 )THEN

                  
         nxp = params % nXelem/params % nProcX
         IF( SUM(nxp) < params % nXElem )THEN

           DO iPx = 1, params % nProcX
             nxp(iPx) = nxp(iPx) + 1
             IF( SUM(nxp) == params % nXElem )THEN
               EXIT
             ENDIF
           ENDDO

         ENDIF

         nyp = params % nYelem/params % nProcY
         IF( SUM(nyp) < params % nYElem )THEN

           DO iPy = 1, params % nProcY
             nyp(iPy) = nyp(iPy) + 1
             IF( SUM(nyp) == params % nYElem )THEN
               EXIT
             ENDIF
           ENDDO

         ENDIF

         nzp = params % nZelem/params % nProcZ
         IF( SUM(nzp) < params % nZElem )THEN

           DO iPz = 1, params % nProcZ
             nzp(iPz) = nzp(iPz) + 1
             IF( SUM(nzp) == params % nZElem )THEN
               EXIT
             ENDIF
           ENDDO

         ENDIF

         iZ = 0
         DO iPz = 1, params % nProcZ
           DO iZp = 1, nzp(iPz)

             iZ = iZ+1
             iY = 0

             DO iPy = 1, params % nProcY
               DO iYp = 1, nyp(iPy)

                 iY = iY+1
                 iX = 0

                 DO iPx = 1, params % nProcX
                   DO iXp = 1, nxp(iPx)

                     iX = iX +1               
                        
                     iEl = iX + params % nXelem*( iY-1 + params % nYelem*( iZ-1 ) )
                     IF( iEl <= myMesh % elements % nElements )THEN
                       partitions(iEl) = iPx-1 + params % nProcX*( iPy-1 + params % nProcY*(iPz-1) )
                     ENDIF

                   ENDDO
                 ENDDO

               ENDDO
             ENDDO

           ENDDO
         ENDDO
      
         DO iEl = 1, myMesh % elements % nElements
      
            procID               = partitions(iEl)
            nElPerProc(procID)   = nElPerProc(procID) + 1 ! Increment the local element ID for this process
            globalToLocal(iEl,1) = nElPerProc(procID)     ! Store the local element ID on this process
            globalToLocal(iEl,2) = procID                 ! Store the process ID
         
            DO iNode = 1, 8
               nID = myMesh % elements % nodeIDs(iNode,iEl) 
               nodeLogic(nID,procID) = 1
            ENDDO
         
         ENDDO
         DO nID = 1, myMesh % nodes % nNodes
            nNodePerProc = nNodePerProc + nodeLogic(nID,:)
            globalToLocalNode(nID,:) = nNodePerProc
         ENDDO

      ELSE

         DO iEl = 1, myMesh % elements % nElements
            nElPerProc(0) = myMesh % elements % nElements
            globalToLocal(iEl,1) = iEl ! Local ID is global ID
            globalToLocal(iEl,2) = 0   ! process ID
            
            DO iNode = 1, 8
               nID = myMesh % elements % nodeIDs(iNode,iEl) 
               nodeLogic(nID,0) = 1
            ENDDO
         ENDDO
         
         DO nID = 1, myMesh % nodes % nNodes
            nNodePerProc = nNodePerProc + nodeLogic(nID,:)
            globalToLocalNode(nID,:) = nNodePerProc
         ENDDO

      ENDIF

  END SUBROUTINE PartitionStructuredElementsAndNodes
!
!==================================================================================================!
!--------------------------------- Mesh File I/O Routines -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup HexMesh_Class
!! @{
! ================================================================================================ !
! S/R ReadSELFMeshFile
!
!> \fn ReadSELFMeshFile_HexMesh
!! Reads a SELFMesh file and constructs the HexMesh DATA structure.
!!
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(HexMesh)    :: this <BR>
!! <B>CHARACTER</B>        :: filename <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ReadSELFMeshFile( filename ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> out <th> myHexMesh <td> HexMesh <td>
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the SELFMesh file
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE ReadSELFMeshFile_HexMesh( myHexMesh, filename )

    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(out)   :: myHexMesh
    CHARACTER(*), INTENT(in)        :: filename
    ! LOCAL
    INTEGER :: nNodes, nElements, nFaces, N
    INTEGER :: IFace, iNode, iEl
    INTEGER :: fUnit, k, i, j, l, row, col, ii, jj
#ifdef HAVE_CUDA
    INTEGER :: istat
#endif


    PRINT(MsgFMT), 'Reading '//TRIM( filename )//'.mesh'
    ! Get a new file unit
    OPEN( UNIT    = NEWUNIT(fUnit), &
      FILE    = TRIM( filename )//'.mesh', &
      FORM    = 'UNFORMATTED',&
      STATUS  = 'OLD', &
      ACCESS  = 'DIRECT', &
      CONVERT = 'BIG_ENDIAN', &
      RECL    = SIZEOF(nNodes) ) ! How to DO variable record length

    ! ---- Gather the number of nodes, number of elements, and number of edges ---- !
    k = 1
    READ( fUnit, rec=k )nNodes
    k = k+1
    READ( fUnit, rec=k )nElements
    k = k+1
    READ( fUnit, rec=k )nFaces
    k = k+1
    READ( fUnit, rec=k )N
    k = k+1

    ! ---- Build the quadrature mesh (empty) ---- !
    CALL myHexMesh % Build( nNodes, nElements, nFaces, N )

    ! ---- Read in the element connectivity ---- !
    DO iEl = 1, nElements
      READ( fUnit, rec=k ) myHexMesh % elements % elementID(iEl)
      k = k+1
      DO i = 1, 8
        READ( fUnit, rec=k ) myHexMesh % elements % nodeIDs(i,iEl)
        k = k+1
      ENDDO
    ENDDO

    ! ---- Read in the face information ---- !

    DO IFace = 1, nFaces
      READ( fUnit, rec=k ) myHexMesh % faces % faceID(IFace)
      k = k+1
      READ( fUnit, rec=k ) myHexMesh % faces % boundaryID(IFace)
      k = k+1
      DO i = 1, 4
        READ( fUnit, rec=k ) myHexMesh % faces % nodeIDs(i,IFace)
        k = k+1
      ENDDO
      DO i = 1, 2
        READ( fUnit, rec=k ) myHexMesh % faces % elementIDs(i,IFace)
        k = k+1
        READ( fUnit, rec=k ) myHexMesh % faces % elementSides(i,IFace)
        k = k+1
      ENDDO
      READ( fUnit, rec=k ) myHexMesh % faces % iStart(IFace)
      k = k+1
      READ( fUnit, rec=k ) myHexMesh % faces % iInc(IFace)
      k = k+1
      READ( fUnit, rec=k ) myHexMesh % faces % jStart(IFace)
      k = k+1
      READ( fUnit, rec=k ) myHexMesh % faces % jInc(IFace)
      k = k+1
      READ( fUnit, rec=k ) myHexMesh % faces % swapDimensions(IFace)
      k = k+1
    ENDDO

    CLOSE( fUnit )
    ! Get a new file unit
    OPEN( UNIT    = NEWUNIT(fUnit), &
      FILE    = TRIM( filename )//'.geom', &
      FORM    = 'UNFORMATTED',&
      STATUS  = 'OLD', &
      ACCESS  = 'DIRECT', &
      CONVERT = 'BIG_ENDIAN', &
      RECL    = prec ) ! How to DO variable record length

    ! ---- Read in the corner nodes ---- !
    k = 1
    DO iNode = 1, nNodes  ! Loop over the nodes in the file
      READ( fUnit, rec=k ) myHexMesh % nodes % x(1,iNode)
      k = k+1
      READ( fUnit, rec=k ) myHexMesh % nodes % x(2,iNode)
      k = k+1
      READ( fUnit, rec=k ) myHexMesh % nodes % x(3,iNode)
      k = k+1
    ENDDO

    ! ---- Read in the element information ---- !
    DO iEl = 1, nElements
      DO l = 0, N
        DO j = 0, N
          DO i = 0, N
            READ( fUnit, rec=k ) myHexMesh % elements % x(i,j,l,1,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % x(i,j,l,2,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % x(i,j,l,3,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % dxds(i,j,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % dxdp(i,j,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % dxdq(i,j,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % dyds(i,j,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % dydp(i,j,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % dydq(i,j,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % dzds(i,j,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % dzdp(i,j,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % dzdq(i,j,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % J(i,j,l,iEl)
            k = k+1
            DO col = 1, 3
              DO row = 1, 3
                READ( fUnit, rec=k ) myHexMesh % elements % Ja(i,j,l,row,col,iEl)
                k = k+1
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO l = 1, 6
        DO j = 0, N
          DO i = 0, N
            READ( fUnit, rec=k ) myHexMesh % elements % boundaryLengthScale(i,j,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % xBound(i,j,1,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % xBound(i,j,2,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % xBound(i,j,3,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % nHat(1,i,j,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % nHat(2,i,j,l,iEl)
            k = k+1
            READ( fUnit, rec=k ) myHexMesh % elements % nHat(3,i,j,l,iEl)
            k = k+1
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CLOSE( fUnit )

    CALL myHexMesh % ConstructElementNeighbors( )

    DO k = 1, myHexMesh % faces % nFaces
      DO j = 0, N
        DO i = 0, N

          IF( i == 0 )THEN
            IF( j == 0 )THEN
              ii = (1-myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % iStart(k)) + &
                (myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % jStart(k))
              jj = (1-myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % jStart(k)) + &
                (myHexMesh % faces % swapDimensions(k))*&
                (myHexMesh % faces % iStart(k))
            ELSE
              ii = myHexMesh % faces % swapDimensions(k)*&
                (ii+myHexMesh % faces % jInc(k)) + &
                (1-myHexMesh % faces % swapDimensions(k))*&
                myHexMesh % faces % iStart(k)
              jj = (1-myHexMesh % faces % swapDimensions(k))*&
                (jj+myHexMesh % faces % jInc(k)) +&
                myHexMesh % faces % swapDimensions(k)*&
                myHexMesh % faces % jStart(k)
            ENDIF
          ELSE
            ii = (1-myHexMesh % faces % swapDimensions(k))*&
              (ii + myHexMesh % faces % iInc(k)) +&
              myHexMesh % faces % swapDimensions(k)*ii
            jj = myHexMesh % faces % swapDimensions(k)*&
              (jj+myHexMesh % faces % iInc(k)) + &
              (1-myHexMesh % faces % swapDimensions(k))*jj
          ENDIF

          myHexMesh % faces % iMap(i,j,k) = ii
          myHexMesh % faces % jMap(i,j,k) = jj

        ENDDO
      ENDDO
    ENDDO

#ifdef HAVE_CUDA
    CALL myHexMesh % UpdateDevice( )
    istat = cudaDeviceSynchronize( )
#endif

  END SUBROUTINE ReadSELFMeshFile_HexMesh
!
!> \addtogroup HexMesh_Class
!! @{
! ================================================================================================ !
! S/R WriteSELFMeshFile
!
!> \fn WriteSELFMeshFile_HexMesh
!! Writes a SELFMesh file using the HexMesh DATA structure.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(HexMesh)    :: this <BR>
!! <B>CHARACTER</B>        :: filename <BR>
!!         .... <BR>
!!     <B>CALL</B> this % WriteSELFMeshFile( filename ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in <th> myHexMesh <td> HexMesh <td>
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the SELFMesh file
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE WriteSELFMeshFile_HexMesh( myHexMesh, filename )

    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(in)   :: myHexMesh
    CHARACTER(*), INTENT(in)       :: filename
    ! LOCAL
    INTEGER :: nNodes, nElements, nFaces, N
    INTEGER :: IFace, iNode, iEl
    INTEGER :: fUnit, k, i, j, l, row, col


    nNodes = 1
    ! Get a new file unit
    OPEN( UNIT    = NEWUNIT(fUnit), &
      FILE    = TRIM( filename )//'.mesh', &
      FORM    = 'UNFORMATTED',&
      STATUS  = 'REPLACE', &
      ACCESS  = 'DIRECT', &
      CONVERT = 'BIG_ENDIAN', &
      RECL    = SIZEOF(nNodes) ) ! How to DO variable record length


    ! ---- Gather the number of nodes, number of elements, and number of edges ---- !
    nNodes = myHexMesh % nodes % nNodes
    nElements = myHexMesh % elements % nElements
    nFaces = myHexMesh % faces % nFaces
    N      = myHexMesh % elements % N
    k = 1
    WRITE( fUnit, rec=k )nNodes
    k = k+1
    WRITE( fUnit, rec=k )nElements
    k = k+1
    WRITE( fUnit, rec=k )nFaces
    k = k+1
    WRITE( fUnit, rec=k )N
    k = k+1


    ! ---- Read in the element connectivity ---- !
    DO iEl = 1, nElements
      WRITE( fUnit, rec=k ) myHexMesh % elements % elementID(iEl)
      k = k+1
      DO i = 1, 8
        WRITE( fUnit, rec=k ) myHexMesh % elements % nodeIDs(i,iEl)
        k = k+1
      ENDDO
    ENDDO

    ! ---- Read in the face information ---- !

    DO IFace = 1, nFaces
      WRITE( fUnit, rec=k ) myHexMesh % faces % faceID(IFace)
      k = k+1
      WRITE( fUnit, rec=k ) myHexMesh % faces % boundaryID(IFace)
      k = k+1
      DO i = 1, 4
        WRITE( fUnit, rec=k ) myHexMesh % faces % nodeIDs(i,IFace)
        k = k+1
      ENDDO
      DO i = 1, 2
        WRITE( fUnit, rec=k ) myHexMesh % faces % elementIDs(i,IFace)
        k = k+1
        WRITE( fUnit, rec=k ) myHexMesh % faces % elementSides(i,IFace)
        k = k+1
      ENDDO
      WRITE( fUnit, rec=k ) myHexMesh % faces % iStart(IFace)
      k = k+1
      WRITE( fUnit, rec=k ) myHexMesh % faces % iInc(IFace)
      k = k+1
      WRITE( fUnit, rec=k ) myHexMesh % faces % jStart(IFace)
      k = k+1
      WRITE( fUnit, rec=k ) myHexMesh % faces % jInc(IFace)
      k = k+1
      WRITE( fUnit, rec=k ) myHexMesh % faces % swapDimensions(IFace)
      k = k+1
    ENDDO

    CLOSE( fUnit )
    ! Get a new file unit
    OPEN( UNIT    = NEWUNIT(fUnit), &
      FILE    = TRIM( filename )//'.geom', &
      FORM    = 'UNFORMATTED',&
      STATUS  = 'REPLACE', &
      ACCESS  = 'DIRECT', &
      CONVERT = 'BIG_ENDIAN', &
      RECL    = prec ) ! How to DO variable record length

    ! ---- Read in the corner nodes ---- !
    k = 1
    DO iNode = 1, nNodes  ! Loop over the nodes in the file
      WRITE( fUnit, rec=k ) myHexMesh % nodes % x(1,iNode)
      k = k+1
      WRITE( fUnit, rec=k ) myHexMesh % nodes % x(2,iNode)
      k = k+1
      WRITE( fUnit, rec=k ) myHexMesh % nodes % x(3,iNode)
      k = k+1
    ENDDO

    ! ---- Read in the element information ---- !
    DO iEl = 1, nElements
      DO l = 0, N
        DO j = 0, N
          DO i = 0, N
            WRITE( fUnit, rec=k ) myHexMesh % elements % x(i,j,l,1,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % x(i,j,l,2,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % x(i,j,l,3,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % dxds(i,j,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % dxdp(i,j,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % dxdq(i,j,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % dyds(i,j,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % dydp(i,j,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % dydq(i,j,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % dzds(i,j,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % dzdp(i,j,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % dzdq(i,j,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % J(i,j,l,iEl)
            k = k+1
            DO col = 1, 3
              DO row = 1, 3
                WRITE( fUnit, rec=k ) myHexMesh % elements % Ja(i,j,l,row,col,iEl)
                k = k+1
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO l = 1, 6
        DO j = 0, N
          DO i = 0, N
            WRITE( fUnit, rec=k ) myHexMesh % elements % boundaryLengthScale(i,j,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % xBound(i,j,1,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % xBound(i,j,2,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % xBound(i,j,3,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % nHat(1,i,j,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % nHat(2,i,j,l,iEl)
            k = k+1
            WRITE( fUnit, rec=k ) myHexMesh % elements % nHat(3,i,j,l,iEl)
            k = k+1
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CLOSE( fUnit )

  END SUBROUTINE WriteSELFMeshFile_HexMesh

  SUBROUTINE ReadTrellisUCDMeshFile_HexMesh( myHexMesh, interp, filename )

    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(out)   :: myHexMesh
    TYPE( Lagrange ), INTENT(in)    :: interp
    CHARACTER(*), INTENT(in)        :: filename
    ! Local 
    INTEGER :: fUnit, readStatus, nNodes, nElements, id, gPolyDeg
    INTEGER :: n1, n2, n3, n4, iEl, iSurf, iSide, i, j
    LOGICAL :: workingOnNodes, workingOnElements, withinFile
    CHARACTER(100) :: lineInTheFile
    TYPE( Surfaces ) :: boundSurfs
    REAL(prec) :: x, y, z, zb, zi, zu, zip1, dxElem, dyElem, dzElem
    REAL(prec) :: x1(1:3), x2(1:3), x3(1:3), x4(1:3)
    REAL(prec) :: c1(1:3), c2(1:3)
    REAL(prec), ALLOCATABLE :: xc(:,:,:,:), s(:), weights(:)



    OPEN( UNIT    = NEWUNIT(fUnit), &
      FILE    = TRIM( filename ), &
      FORM    = 'FORMATTED',&
      STATUS  = 'OLD', &
      ACCESS  = 'SEQUENTIAL' )

    withinFile        = .TRUE.
    workingOnNodes    = .FALSE.
    workingOnElements = .FALSE.

    nNodes    = 0
    nElements = 0
    DO WHILE ( withinFile )

      READ( fUnit, '(A100)', IOSTAT=readStatus ) lineInTheFile

      IF( readStatus == 0 )THEN
        
        IF( lineInTheFile(1:5) == '*NODE' ) THEN

          workingOnNodes    = .TRUE.
          workingOnElements = .FALSE.

        ELSEIF( lineInTheFile(1:5) == '*ELEM' ) THEN

          workingOnNodes    = .FALSE.
          workingOnElements = .TRUE.

        ELSEIF( lineInTheFile(1:2) == '**' ) THEN

          ! Skip comment lines
          workingOnNodes    = .FALSE.
          workingOnElements = .FALSE.

        ELSE

          IF( workingOnNodes )THEN

            nNodes = nNodes + 1       

          ELSEIF( workingOnElements )THEN

            nElements = nElements + 1

          ENDIF


        ENDIF

      ELSE

        withinFile = .FALSE.

      ENDIF

    ENDDO

    REWIND( fUnit )

    CALL myHexMesh % Build( nNodes, nElements, 1, interp % N ) 

    withinFile        = .TRUE.
    workingOnNodes    = .FALSE.
    workingOnElements = .FALSE.
   
    nNodes    = 0
    nElements = 0

    DO WHILE ( withinFile )

      READ( fUnit, '(A100)', IOSTAT=readStatus ) lineInTheFile

      IF( readStatus == 0 )THEN
        
        IF( lineInTheFile(1:5) == '*NODE' ) THEN

          workingOnNodes    = .TRUE.
          workingOnElements = .FALSE.

        ELSEIF( lineInTheFile(1:5) == '*ELEM' ) THEN

          workingOnNodes    = .FALSE.
          workingOnElements = .TRUE.

        ELSEIF( lineInTheFile(1:2) == '**' ) THEN

          ! Skip comment lines
          workingOnNodes    = .FALSE.
          workingOnElements = .FALSE.

        ELSE

          IF( workingOnNodes )THEN

            nNodes = nNodes + 1       
            READ( lineInTheFile, '(I6,3(",",2x,E14.6))' ) id, myHexMesh % nodes % x(1:3,nNodes)

          ELSEIF( workingOnElements )THEN

            nElements = nElements + 1
            READ( lineInTheFile, '(9(2x,I6))' ) id, myHexMesh % elements % nodeIDs(1:8,nElements)

          ENDIF


        ENDIF

      ELSE

        withinFile = .FALSE.

      ENDIF

    ENDDO

    CLOSE( fUnit )


    gPolyDeg = interp % N
    ALLOCATE( s(0:gPolyDeg), xc(0:gPolyDeg,0:gPolyDeg,1:3,1:6*nElements), weights(0:gpolyDeg) )
    CALL LegendreQuadrature( gPolyDeg, s, weights, GAUSS_LOBATTO )

    xc = 0.0_prec
    CALL boundSurfs % Build( s, gPolyDeg, 6*nElements )

    DO iEl = 1, nElements

      myHexMesh % elements % elementID(iEl) = iEl

      DO iSide = 1, 6 ! Loop over the sides of the quads

        iSurf = iSide + 6*(iEl-1)
        ! To build the current face, we construct a plane that passes through
        ! the four corner nodes. Here, we grab the global node ID's for the four
        ! corner nodes.
        n1 = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(1,iSide),iEl )
        n2 = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(2,iSide),iEl )
        n3 = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(3,iSide),iEl )
        n4 = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(4,iSide),iEl )

        x1(1) = myHexMesh % nodes % x(1,n1)
        x1(2) = myHexMesh % nodes % x(2,n1)
        x1(3) = myHexMesh % nodes % x(3,n1)

        x2(1) = myHexMesh % nodes % x(1,n2)
        x2(2) = myHexMesh % nodes % x(2,n2)
        x2(3) = myHexMesh % nodes % x(3,n2)

        x3(1) = myHexMesh % nodes % x(1,n3)
        x3(2) = myHexMesh % nodes % x(2,n3)
        x3(3) = myHexMesh % nodes % x(3,n3)

        x4(1) = myHexMesh % nodes % x(1,n4)
        x4(2) = myHexMesh % nodes % x(2,n4)
        x4(3) = myHexMesh % nodes % x(3,n4)

        DO j = 0, gPolyDeg
          DO i = 0, gPolyDeg
            ! Transfinite inerpolation with linear blending is USEd to construct the face
            c1 = ( 0.5_prec*(x2-x1)*(1.0_prec+s(i)) + x1 )
            c2 = ( 0.5_prec*(x3-x4)*(1.0_prec+s(i)) + x4 )
            xc(i,j,1:3,iSurf) = 0.5_prec*(c2-c1)*(1.0_prec+s(j)) + c1

          ENDDO
        ENDDO

      ENDDO

    ENDDO

    CALL boundSurfs % Set_Surfaces( xc )
    CALL myHexMesh % elements % GenerateMesh( interp, boundSurfs )
    CALL myHexMesh % elements % GenerateMetrics( interp )

    CALL myHexMesh % ConstructFaces( )


    CALL myHexMesh % ConstructElementNeighbors( )

    DEALLOCATE( s, xc, weights )
    CALL boundSurfs % Trash( )

  END SUBROUTINE ReadTrellisUCDMeshFile_HexMesh

!
!> \addtogroup HexMesh_Class
!! @{
! ================================================================================================ !
! S/R ReadUCDMeshFile
!
!> \fn ReadUCDMeshFile_HexMesh
!! Reads in a ucd unstructured mesh file. It is assumed that the nodes and elements are numbered
!! between 1 to nNodes and 1 to nElements, respectively. To ensure this in Trellis or CuBit, TYPE
!! "compress all" before exporting ucd file.
!!
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(HexMesh)    :: this <BR>
!! <B>CHARACTER</B>        :: filename <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ReadUCDMeshFile( filename ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> out <th> myHexMesh <td> HexMesh <td>
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the SELFMesh file
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE ReadUCDMeshFile_HexMesh( myHexMesh, interp, filename )

    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(out)   :: myHexMesh
    TYPE( Lagrange ), INTENT(in)    :: interp
    CHARACTER(*), INTENT(in)        :: filename
    ! LOCAL
    CHARACTER(60) :: longDummy
    CHARACTER(3)  :: shortDummy
    INTEGER :: nNodes, nElements, nFaces
    INTEGER :: IFace, iNode, iEl, iSide
    INTEGER :: fUnit, k, i, j, l, row, col, n1, n2, n3, n4, iSurf
    REAL(prec) :: s(0:1)
    REAL(prec), ALLOCATABLE :: xc(:,:,:,:)
    REAL(prec) :: x1(1:3), x2(1:3), x3(1:3), x4(1:3), c1(1:3), c2(1:3)
    TYPE( Surfaces ) :: boundSurfs


    PRINT(MsgFMT), 'Mesh File : '//TRIM( filename )

    ! Get a new file unit
    OPEN( UNIT    = NEWUNIT(fUnit), &
      FILE    = TRIM( filename ), &
      FORM    = 'FORMATTED',&
      STATUS  = 'OLD', &
      ACCESS  = 'SEQUENTIAL' )

    READ( fUnit, * ) longDummy
    WRITE( *, * ) longDummy

    READ( fUnit, * ) longDummy
    WRITE( *, * ) longDummy

    READ( fUnit, * ) longDummy
    WRITE( *, * ) longDummy

    READ( fUnit, * ) nNodes, nElements, i, j, l ! i,j, and l are just fillers for now.


    ! ---- Gather the number of nodes, number of elements, and number of edges ---- !
    ! Structured nFaces = 1 for initial build
    nFaces = 1

    ! ---- Build the quadrature mesh (empty) ---- !
    CALL myHexMesh % Build( nNodes, nElements, nFaces, interp % N )

    ALLOCATE( xc(0:1,0:1,1:3,1:6*nElements) )

    ! Read in the corner node positions
    DO iNode = 1, nNodes
      READ( fUnit, * ) myHexMesh % nodes % nodeID(iNode), &
        myHexMesh % nodes % x(1,iNode), &
        myHexMesh % nodes % x(2,iNode), &
        myHexMesh % nodes % x(3,iNode)
    ENDDO

    ! ---- Read in the element connectivity ---- !
    DO iEl = 1, nElements
      myHexMesh % elements % elementID(iEl) = iEl
      READ( fUnit, * ) i, j, shortDummy, &
        myHexMesh % elements % nodeIDs(1:8,iEl)
    ENDDO
    CLOSE( fUnit )

    CALL myHexMesh % ConstructFaces( )
    CALL myHexMesh % ConstructElementNeighbors( )

    ! DO the element information
    s(0) = -1.0_prec
    s(1) = 1.0_prec
    xc   = 0.0_prec

    CALL boundSurfs % Build( s, 1, 6*nElements )


    ! For now, we assume tri-linear mappings so that the geometry can be constructed
    ! using only the corner node locations. These corner nodes will be USEd to generate
    ! the bounding surfaces for each element, and transfinite interpolation with linear
    ! blending is USEd to construct the internal mesh geometry.
    DO iEl = 1, nElements

      DO iSide = 1, 6 ! Loop over the sides of the quads

        iSurf = iSide + (iEl-1)*6
        ! To build the current face, we construct a plane that passes through
        ! the four corner nodes. Here, we grab the global node ID's for the four
        ! corner nodes.
        n1 = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(1,iSide), iEl )
        n2 = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(2,iSide), iEl )
        n3 = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(3,iSide), iEl )
        n4 = myHexMesh % elements % nodeIDs( myHexMesh % faceMap(4,iSide), iEl )

        x1(1) = myHexMesh % nodes % x(1,n1)
        x1(2) = myHexMesh % nodes % x(2,n1)
        x1(3) = myHexMesh % nodes % x(3,n1)

        x2(1) = myHexMesh % nodes % x(1,n2)
        x2(2) = myHexMesh % nodes % x(2,n2)
        x2(3) = myHexMesh % nodes % x(3,n2)

        x3(1) = myHexMesh % nodes % x(1,n3)
        x3(2) = myHexMesh % nodes % x(2,n3)
        x3(3) = myHexMesh % nodes % x(3,n3)

        x4(1) = myHexMesh % nodes % x(1,n4)
        x4(2) = myHexMesh % nodes % x(2,n4)
        x4(3) = myHexMesh % nodes % x(3,n4)

        DO j = 0, 1
          DO i = 0, 1
            ! Transfinite inerpolation with linear blending is USEd to construct the face
            c1 = ( 0.5_prec*(x2-x1)*(1.0_prec+s(i)) + x1 )
            c2 = ( 0.5_prec*(x3-x4)*(1.0_prec+s(i)) + x4 )
            xc(i,j,1:3,iSurf) = 0.5_prec*(c2-c1)*(1.0_prec+s(j)) + c1
          ENDDO
        ENDDO

      ENDDO
    ENDDO

    CALL boundSurfs % Set_Surfaces( xc )
    CALL myHexMesh % elements % GenerateMesh( interp, boundSurfs )
    CALL myHexMesh % elements % GenerateMetrics( interp )


    CALL boundSurfs % Trash( )

    DEALLOCATE( xc )


#ifdef HAVE_CUDA
    CALL myHexMesh % UpdateDevice( )
#endif


  END SUBROUTINE ReadUCDMeshFile_HexMesh
!
!> \addtogroup HexMesh_Class
!! @{
! ================================================================================================ !
! S/R WriteTecplot
!
!> \fn WriteTecplot_HexMesh
!! Writes a tecplot file of the mesh geometry.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(HexMesh) :: this <BR>
!! <B>CHARACTER</B>      :: filename <BR>
!!         .... <BR>
!!     ! To write a file with a specIFied name <BR>
!!     <B>CALL</B> this % WriteTecplot( filename ) <BR>
!!     ! Or, the file "mesh.tec" will be USEd with the following calling sequence <BR>
!!     <B>CALL</B> this % WriteTecplot( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in <th> myHexMesh <td> HexMesh <td>
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the tecplot file
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE WriteTecplot_Hexmesh( myHexMesh, filename )

    IMPLICIT NONE
    CLASS(HexMesh), INTENT(inout)     :: myHexMesh
    CHARACTER(*), INTENT(in), OPTIONAL :: filename
    ! Local
    INTEGER :: i,j,k, N, iEl, fUnit, eID
    CHARACTER(7) :: zoneID

    N = myHexMesh % elements % N

    IF( PRESENT(filename) )THEN

      OPEN( UNIT=NEWUNIT(fUnit), &
        FILE= TRIM(filename)//'.tec', &
        FORM='formatted')
    ELSE

      OPEN( UNIT=NEWUNIT(fUnit), &
        FILE= 'mesh.tec', &
        FORM='formatted')

    ENDIF

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "Jacobian", "dxds", "dxdp", "dxdq", "dyds", "dydp", "dydq", "dzds", "dzdp", "dzdq" '


    DO iEl = 1, myHexMesh % elements % nElements

      eID = myHexMesh % elements % elementID(iEl)
      WRITE(zoneID,'(I7.7)') eID
      WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',N+1,', J=',N+1,', K=',N+1,',F=POINT'

      DO k = 0, N
        DO j = 0, N
          DO i = 0,N

            WRITE(fUnit,'(13(E15.7,1x))')  myHexMesh % elements % x(i,j,k,1,iEl), &
              myHexMesh % elements % x(i,j,k,2,iEl), &
              myHexMesh % elements % x(i,j,k,3,iEl), &
              myHexMesh % elements % J(i,j,k,iEl), &
              myHexMesh % elements % dxds(i,j,k,iEl), &
              myHexMesh % elements % dxdp(i,j,k,iEl), &
              myHexMesh % elements % dxdq(i,j,k,iEl), &
              myHexMesh % elements % dyds(i,j,k,iEl), &
              myHexMesh % elements % dydp(i,j,k,iEl), &
              myHexMesh % elements % dydq(i,j,k,iEl), &
              myHexMesh % elements % dzds(i,j,k,iEl), &
              myHexMesh % elements % dzdp(i,j,k,iEl), &
              myHexMesh % elements % dzdq(i,j,k,iEl)
          ENDDO
        ENDDO
      ENDDO

    ENDDO

    CLOSE(UNIT=fUnit)

    RETURN

  END SUBROUTINE WriteTecplot_Hexmesh
!
!> \addtogroup HexMesh_Class
!! @{
! ================================================================================================ !
! S/R WriteMaterialTecplot
!
!> \fn WriteMaterialTecplot_HexMesh
!! Writes a tecplot file of the mesh geometry and a given "material" field.
!!
!! When performing DOmain decomposition (e.g. with the DecomposeHexMesh.f90 program), the
!! decomposition can be visualized by assigning each element a number that corresponds to the
!! process ID it has been assigned to. This routine takes in an array of "material ID's" that
!! identIFy each element.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(HexMesh) :: this <BR>
!! <B>REAL</B>(prec)     :: pID(1:this % nElements) <BR>
!! <B>CHARACTER</B>      :: filename <BR>
!!         .... <BR>
!!     ! To write a file with a specIFied name <BR>
!!     <B>CALL</B> this % WriteMaterialTecplot( pID, filename ) <BR>
!!     ! Or, the file "mesh.tec" will be USEd with the following calling sequence <BR>
!!     <B>CALL</B> this % WriteMaterialTecplot( pID ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in <th> myHexMesh <td> HexMesh <td>
!!   <tr> <td> in <th> materialIDs(1:myHexMesh % nElements) <td> REAL(prec) <td>
!!                     Array of values that are assigned to each element.
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the MaterialTecplot file
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE WriteMaterialTecplot_Hexmesh( myHexMesh, materialIDs, filename )

    IMPLICIT NONE
    CLASS(HexMesh), INTENT(inout)      :: myHexMesh
    REAL(prec), INTENT(in)             :: materialIDs(1:myHexMesh % elements % nElements)
    CHARACTER(*), INTENT(in), OPTIONAL :: filename
    ! Local
    INTEGER :: i, j, k, N, iEl,fUnit
    CHARACTER(7) :: zoneID

    N = myHexMesh % elements % N

    IF( PRESENT(filename) )THEN
      OPEN( UNIT=NEWUNIT(fUnit), &
        FILE= TRIM(filename)//'.tec', &
        FORM='formatted')
    ELSE
      OPEN( UNIT=NEWUNIT(fUnit), &
        FILE= 'mesh.tec', &
        FORM='formatted')
    ENDIF

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "materialIDs" '

    DO iEl = 1, myHexMesh % elements % nElements

      WRITE(zoneID,'(I7.7)') myHexMesh % elements % elementID(iEl)
      WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',N+1,', J=',N+1,', K=',N+1,',F=POINT'

      DO k = 0, N
        DO j = 0, N
          DO i = 0, N
            WRITE(fUnit,*) myHexMesh % elements % x(i,j,k,1,iEl), &
              myHexMesh % elements % x(i,j,k,2,iEl), &
              myHexMesh % elements % x(i,j,k,3,iEl), &
              materialIDs(iEl)
          ENDDO
        ENDDO
      ENDDO

    ENDDO

    CLOSE(UNIT=fUnit)
    RETURN

  END SUBROUTINE WriteMaterialTecplot_Hexmesh
!
  FUNCTION NumberOfBoundaryFaces( mesh ) RESULT( nBFaces )
  
    IMPLICIT NONE
    CLASS( HexMesh ) :: mesh
    INTEGER          :: nBFaces
    ! Local
    INTEGER :: iFace

    nBFaces = 0

    DO iFace = 1, mesh % faces % nFaces

      IF( mesh % faces % elementIDs(2,iFace) < 0 )THEN

        nBFaces = nBFaces + 1 

      ENDIF

    ENDDO   
    
  END FUNCTION NumberOfBoundaryFaces
  
 SUBROUTINE StructuredMeshGenerator_3D( paramFile, setupSuccess )
 IMPLICIT NONE
 LOGICAL, INTENT(out)                      :: setupSuccess
 CHARACTER(*), INTENT(in)                  :: paramFile
  ! Local
 TYPE( NodalDG )                           :: nodal
 TYPE( HexMesh )                           :: mesh
 TYPE( HexMesh )                           :: procMesh
 TYPE( ModelParameters )                   :: params
 TYPE( Geom_EquationParser )               :: geomParser
 TYPE( BoundaryCommunicator ), ALLOCATABLE :: bcom(:)
 INTEGER, ALLOCATABLE                      :: faceProcCount(:), faceProcTable(:,:), faceProcOwners(:,:), faceBoundaryIDs(:,:)
 INTEGER                      :: procID, mpiErr
 INTEGER                      :: nElems, nProc, pID, i, j
 INTEGER                      :: iEl, jEl, iSide, iNode, nAdj, nID, elID, nMPI
 INTEGER                      :: e1, e2, p1, p2, iFace, iFaceLocal, jFace, localID
 INTEGER                      :: nBe, npFaces
 INTEGER                      :: globalFaceID, localFaceID, bID, extProc, extBID
 CHARACTER(4)                 :: pIDChar
 INTEGER, ALLOCATABLE         :: globalToLocal(:,:), nElPerProc(:), nodeLogic(:,:), nNodePerProc(:), partitions(:)
 INTEGER, ALLOCATABLE         :: globalToLocalNode(:,:), nLocMPI(:)
 REAL(prec), ALLOCATABLE      :: materials(:)
 

      ! Read in the parameters
      CALL params % Build( TRIM(paramFile), setupSuccess )

#ifdef HAVE_MPI
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, params % nProc, mpiErr )

      IF( params % nProc > params % nXElem*params % nYElem*params % nZElem )THEN

        PRINT*, ' Number of processes exceeds number of elements.'
        setupSuccess = .FALSE.

      ELSEIF( params % nProc == params % nXElem*params % nYElem*params % nZElem )THEN

        params % nProcX = params % nXElem
        params % nProcY = params % nYElem
        params % nProcZ = params % nZElem

      ELSE

        IF( floorSQRT(params % nProc)**2 == params % nProc )THEN

          params % nProcZ = 1
          params % nProcY = floorSQRT(params % nProc)
          params % nProcX = floorSQRT(params % nProc)

        ELSEIF(  floorCURT( params % nProc )**3 ==  params % nProc )THEN

          params % nProcZ = floorCURT( params % nProc )
          params % nProcY = floorCURT( params % nProc )
          params % nProcX = floorCURT( params % nProc )

        ELSE

          params % nProcZ = 1 
          params % nProcY = 1 
          params % nProcX = params % nProc 

        ENDIF

          

        DO WHILE( params % nProcX > params % nXElem )

          DO i = 2, params % nProcX
            IF( MOD( params % nProcX, i ) == 0 )THEN
              ! nProcX is divisible by i
              j = i
              EXIT
            ENDIF
          ENDDO

          params % nProcX = params % nProcX/j

        ENDDO

        IF( params % nProcX*params % nProcY*params % nProcZ < params % nProc )THEN

          params % nProcY = params % nProc/params % nProcX
          IF( params % nProcY > 1 )THEN
            DO WHILE( params % nProcY > params % nYElem )
  
              DO i = 2, params % nProcY
                IF( MOD( params % nProcY, i ) == 0 )THEN
                  ! nProcY is divisible by i
                  j = i
                  EXIT
                ENDIF
              ENDDO
  
              params % nProcY = params % nProcY/j
  
            ENDDO
  
          ENDIF

        ENDIF

        IF( params % nProcX*params % nProcY*params % nProcZ < params % nProc )THEN

          params % nProcZ = params % nProc/(params % nProcX*params % nProcY)
          IF( params % nProcZ > 1 )THEN
            DO WHILE( params % nProcZ > params % nZElem )
  
              DO i = 2, params % nProcZ
                IF( MOD( params % nProcZ, i ) == 0 )THEN
                  ! nProcZ is divisible by i
                  j = i
                  EXIT
                ENDIF
              ENDDO
  
              params % nProcZ = params % nProcZ/j
  
            ENDDO
  
          ENDIF

        ENDIF

      ENDIF

      PRINT*, ' nProcX :', params % nProcX
      PRINT*, ' nProcY :', params % nProcY
      PRINT*, ' nProcZ :', params % nProcZ
      
#endif
      CALL geomParser % Build( 'self.equations' )

      IF( setupSuccess )THEN
         ! Build an interpolant
         CALL nodal % Build( targetPoints = UniformPoints(-1.0_prec,1.0_prec,params % nPlot), &
                             N = params % polyDeg, &
                             nTargetPoints = params % nPlot, &
                             quadrature = GAUSS_LOBATTO )
   
         ! Build the Geometry
         IF( TRIM( params % UCDMeshFile ) == '' )THEN

           IF( params % MeshType == DoublyPeriodic )THEN
              CALL mesh % ConstructStructuredMesh( nodal % interp, &
                params % nXelem, &
                params % nYelem, &
                params % nZelem, &
                geomParser, &
                .TRUE. )
           ELSE
              CALL mesh % ConstructStructuredMesh( nodal % interp, &
                params % nXelem, &
                params % nYelem, &
                params % nZelem, &
                geomParser, &
                .FALSE. )
           ENDIF

         ELSE   

           CALL mesh % ReadTrellisUCDMeshFile( nodal % interp, TRIM( params % UCDMeshFile ) )           
#ifdef TECPLOT
           CALL mesh % WriteTecplot( 'mesh' )
#endif

         ENDIF

         nElems = mesh % elements % nElements
         nProc  = params % nProcX*params % nProcY*params % nProcZ
         IF( nProc == 0 )THEN
            nProc = 1
         ENDIF
         ALLOCATE( materials(1:nElems), &
                   globalToLocal(1:nElems,1:2), &
                   nElPerProc(0:nProc-1), &
                   nodeLogic(1:mesh % nodes % nNodes, 0:nProc-1), &
                   nNodePerProc(0:nProc-1), &
                   globalToLocalNode(1:mesh % nodes % nNodes, 0:nProc-1), &
                   nLocMPI(0:nProc-1) )
                   
         materials     = 0.0_prec
         globalToLocal = 0
         nElPerProc    = 0
         nodeLogic     = 0
         nNodePerProc  = 0
         
         ALLOCATE( partitions(1:nElems) )
         ALLOCATE( bcom(0:nProc-1) ) 
         ALLOCATE( faceProcCount(1:mesh % faces % nFaces), &
                   faceProcOwners(1:mesh % faces % nFaces,1:2), &
                   faceBoundaryIDs(1:mesh % faces % nFaces,1:2), &
                   faceProcTable(1:mesh % faces % nFaces, 0:nProc-1) )
         faceProcCount   = 0
         faceProcOwners  = -1
         faceBoundaryIDs = 0
         faceProcTable   = -5000

         CALL mesh % PartitionStructuredElementsAndNodes( params, partitions, nElPerProc, &
                                                          globalToLocal, nodeLogic, nNodePerProc, &
                                                          globalToLocalNode, nProc )

   
            
         ! Now we generate the local mesh for each process

         DO procID = 0, nProc-1
            CALL procMesh % Build( nNodePerProc(procID), &
                                   nElPerProc(procID), &
                                   1, params % polyDeg )
! ----------------------------------------------------------------------------- !

            DO nID = 1, mesh % nodes % nNodes
            
               ! We only assign node information if it is "owned" by this process
               IF( nodeLogic(nID,procID) == 1 )THEN 

                  localID = globalToLocalNode( nID, procID )
                  procMesh % nodes % nodeID(localID)   = nID
                  procMesh % nodes % nodeType(localID) = mesh % nodes % nodeTYPE(nID)
                  procMesh % nodes % x(1:3,localID)    = mesh % nodes % x(1:3,nID)

               ENDIF
               
            ENDDO
   
! ----------------------------------------------------------------------------- !

            ! Set the element-node connectivity
            DO iEl = 1, mesh % elements % nElements
               ! We only assign element information if it is "owned" by this process
               IF( globalToLocal(iEl,2) == procID )THEN
                  elID = globalToLocal(iEl,1) ! Obtain the local element ID from the global-to-local array
                  
                  DO iNode = 1, 8
       
                    ! First we grab the global node ID for this element's corner node
                    nID = mesh % elements % nodeIDs(iNode,iEl)
                    ! THEN we convert this to a local node ID for this process
                    localID = globalToLocalNode( nID, procID )
                    ! And now we assign the local Node ID to the local element's corner node ID's
                    procMesh % elements % nodeIDs(iNode,elID)  = localID
       
                  ENDDO
                  procMesh % elements % elementID(elID) = iEl
                     
                  ! Set the geometry
                  procMesh % elements % nHat(:,:,:,:,elID)   = mesh % elements % nHat(:,:,:,:,iEl)
                  procMesh % elements % xBound(:,:,:,:,elID) = mesh % elements % xBound(:,:,:,:,iEl)
                  procMesh % elements % x(:,:,:,:,elID)      = mesh % elements % x(:,:,:,:,iEl)
                  procMesh % elements % J(:,:,:,elID)        = mesh % elements % J(:,:,:,iEl)
                  procMesh % elements % dxds(:,:,:,elID)     = mesh % elements % dxds(:,:,:,iEl)
                  procMesh % elements % dxdp(:,:,:,elID)     = mesh % elements % dxdp(:,:,:,iEl)
                  procMesh % elements % dxdq(:,:,:,elID)     = mesh % elements % dxdq(:,:,:,iEl)
                  procMesh % elements % dyds(:,:,:,elID)     = mesh % elements % dyds(:,:,:,iEl)
                  procMesh % elements % dydp(:,:,:,elID)     = mesh % elements % dydp(:,:,:,iEl)
                  procMesh % elements % dydq(:,:,:,elID)     = mesh % elements % dydq(:,:,:,iEl)
                  procMesh % elements % dzds(:,:,:,elID)     = mesh % elements % dzds(:,:,:,iEl)
                  procMesh % elements % dzdp(:,:,:,elID)     = mesh % elements % dzdp(:,:,:,iEl)
                  procMesh % elements % dzdq(:,:,:,elID)     = mesh % elements % dzdq(:,:,:,iEl)
                  procMesh % elements % Ja(:,:,:,:,:,elID)   = mesh % elements % Ja(:,:,:,:,:,iEl)
       
               ENDIF
               
            ENDDO
   
! ----------------------------------------------------------------------------- !

            npFaces = 0 
            nBe     = 0
            iFaceLocal = 0
            nMPI       = 0
            DO iFace = 1, mesh % faces % nFaces

               e1 = mesh % faces % elementIDs(1,iFace)
               e2 = mesh % faces % elementIDs(2,iFace)

               p1 = globalToLocal(e1,2)
               IF( e2 > 0 )THEN
                  p2 = globalToLocal(e2,2)
               ELSE
                  p2  = p1
               ENDIF
               IF( p1 == procID .OR. p2 == procID )THEN

                  faceProcTable(iFace,p1) = 1
                  faceProcTable(iFace,p2) = 1
                  iFaceLocal = iFaceLocal+1

                  faceProcCount(iFace) = faceProcCount(iFace) + 1
                  faceProcOwners(iFace,faceProcCount(iFace)) = procID 

                  npFaces = npFaces + 1

                  IF( e2 < 0 .AND. p2 == p1 )THEN

                     nBe = nBe + 1
                     faceBoundaryIDs(iFace, faceProcCount(iFace)) = nBe

                  ELSEIF( p2 /= p1 .AND. p1 == procID )THEN

                     nMPI = nMPI + 1 
                     nBe = nBe + 1
                     faceBoundaryIDs(iFace, faceProcCount(iFace)) = nBe

                  ELSEIF( p2 /= p1 .AND. p2 == procID )THEN

                     nMPI = nMPI + 1 
                     nBe = nBe + 1
                     faceBoundaryIDs(iFace, faceProcCount(iFace)) = nBe

                  ENDIF

               ENDIF

            ENDDO

            PRINT*, '      Process ID        nMPI       nFaces       nBFaces       nElements' 
            PRINT*,  procID, nMPI, npFaces, nBe, nElPerProc(procID)

            CALL procMesh % faces % Trash( )
            CALL procMesh % faces % Build( npFaces, params % polyDeg ) 
            CALL bCom(procID) % Build( nBe )

            iFaceLocal = 0
            nBe        = 0
            nMPI       = 0
            DO iFace = 1, mesh % faces % nFaces

               e1 = mesh % faces % elementIDs(1,iFace)
               e2 = mesh % faces % elementIDs(2,iFace)

               p1 = globalToLocal(e1,2)
               IF( e2 > 0 )THEN
                  p2 = globalToLocal(e2,2)
               ELSE
                  p2  = p1
               ENDIF

               IF( p1 == procID .OR. p2 == procID )THEN

                  iFaceLocal = iFaceLocal + 1


                  IF( p2 == p1 .AND. e2 > 0 )THEN ! Internal Face
                    
                     procMesh % faces % faceID(iFaceLocal)         = iFace
                     procMesh % faces % elementIDs(1,iFaceLocal)   = globalToLocal( e1, 1 )
                     procMesh % faces % elementIDs(2,iFaceLocal)   = globalToLocal( e2, 1 )
                     procMesh % faces % elementSides(1,iFaceLocal) = mesh % faces % elementSides(1,iFace)
                     procMesh % faces % elementSides(2,iFaceLocal) = mesh % faces % elementSides(2,iFace)
                     procMesh % faces % iStart(iFaceLocal)         = mesh % faces % iStart(iFace)
                     procMesh % faces % iInc(iFaceLocal)           = mesh % faces % iInc(iFace)
                     procMesh % faces % jStart(iFaceLocal)         = mesh % faces % jStart(iFace)
                     procMesh % faces % jInc(iFaceLocal)           = mesh % faces % jInc(iFace)
                     procMesh % faces % swapDimensions(iFaceLocal) = mesh % faces % swapDimensions(iFace)
                     procMesh % faces % iMap(:,:,iFaceLocal)       = mesh % faces % iMap(:,:,iFace)
                     procMesh % faces % jMap(:,:,iFaceLocal)       = mesh % faces % jMap(:,:,iFace)

                  ELSEIF( p2 == p1 .AND. e2 < 0 )THEN ! Physical boundary

                     procMesh % faces % faceID(iFaceLocal)         = iFace
                     procMesh % faces % elementIDs(1,iFaceLocal)   = globalToLocal( e1, 1 )
                     procMesh % faces % elementIDs(2,iFaceLocal)   = e2
                     procMesh % faces % elementSides(1,iFaceLocal) = mesh % faces % elementSides(1,iFace)
                     procMesh % faces % elementSides(2,iFaceLocal) = mesh % faces % elementSides(2,iFace)
                     procMesh % faces % iStart(iFaceLocal)         = mesh % faces % iStart(iFace)
                     procMesh % faces % iInc(iFaceLocal)           = mesh % faces % iInc(iFace)
                     procMesh % faces % jStart(iFaceLocal)         = mesh % faces % jStart(iFace)
                     procMesh % faces % jInc(iFaceLocal)           = mesh % faces % jInc(iFace)
                     procMesh % faces % swapDimensions(iFaceLocal) = mesh % faces % swapDimensions(iFace)
                     procMesh % faces % iMap(:,:,iFaceLocal)       = mesh % faces % iMap(:,:,iFace)
                     procMesh % faces % jMap(:,:,iFaceLocal)       = mesh % faces % jMap(:,:,iFace)

                     nBe = nBe + 1
                     procMesh % faces % boundaryID(iFaceLocal) = -nBe

                     ! Physical boundary ID's will be set to negative for physical boundaries
                     ! so that boundaries requiring communication can be separated from physical
                     ! boundaries. This will allow more operations to be run concurrently with
                     ! communication.
                     bCom(procID) % boundaryIDs(nBe) = iFaceLocal
                     bCom(procID) % boundaryGlobalIDs(nBe) = iFace
                     bCom(procID) % extProcIDs(nBe)  = procID

                  ELSEIF( p2 /= p1 .AND. procID == p1 )THEN ! MPI Boundary
                
                     procMesh % faces % faceID(iFaceLocal)         = iFace
                     procMesh % faces % elementIDs(1,iFaceLocal)   = globalToLocal( e1, 1 )
                     procMesh % faces % elementIDs(2,iFaceLocal)   = -e2
                     procMesh % faces % elementSides(1,iFaceLocal) = mesh % faces % elementSides(1,iFace)
                     procMesh % faces % elementSides(2,iFaceLocal) = mesh % faces % elementSides(2,iFace)
                     procMesh % faces % iStart(iFaceLocal)         = mesh % faces % iStart(iFace)
                     procMesh % faces % iInc(iFaceLocal)           = mesh % faces % iInc(iFace)
                     procMesh % faces % jStart(iFaceLocal)         = mesh % faces % jStart(iFace)
                     procMesh % faces % jInc(iFaceLocal)           = mesh % faces % jInc(iFace)
                     procMesh % faces % swapDimensions(iFaceLocal) = mesh % faces % swapDimensions(iFace)
                     procMesh % faces % iMap(:,:,iFaceLocal)       = mesh % faces % iMap(:,:,iFace)
                     procMesh % faces % jMap(:,:,iFaceLocal)       = mesh % faces % jMap(:,:,iFace)

                     nBe = nBe + 1
                     procMesh % faces % boundaryID(iFaceLocal) = nBe

                     ! Boundary ID's associated with MPI boundaries are reported as positive integers
                     ! so that they can be distinguished from physical boundaries that do not require
                     ! communication with neighboring processes. 
                     bCom(procID) % boundaryIDs(nBe) = iFaceLocal
                     bCom(procID) % boundaryGlobalIDs(nBe) = iFace
                     bCom(procID) % extProcIDs(nBe)  = p2


                  ELSEIF( p2 /= p1 .AND. procID == p2 )THEN ! MPI Boundary
                  
                  
                     procMesh % faces % faceID(iFaceLocal)         = iFace
                     procMesh % faces % elementIDs(1,iFaceLocal)   = globalToLocal( e2, 1 )
                     procMesh % faces % elementIDs(2,iFaceLocal)   = -e1
                     procMesh % faces % elementSides(1,iFaceLocal) = mesh % faces % elementSides(2,iFace)
                     procMesh % faces % elementSides(2,iFaceLocal) = mesh % faces % elementSides(1,iFace)
                     procMesh % faces % iStart(iFaceLocal)         = mesh % faces % iStart(iFace)
                     procMesh % faces % iInc(iFaceLocal)           = mesh % faces % iInc(iFace)
                     procMesh % faces % jStart(iFaceLocal)         = mesh % faces % jStart(iFace)
                     procMesh % faces % jInc(iFaceLocal)           = mesh % faces % jInc(iFace)
                     procMesh % faces % swapDimensions(iFaceLocal) = mesh % faces % swapDimensions(iFace)
                     procMesh % faces % iMap(:,:,iFaceLocal)       = mesh % faces % iMap(:,:,iFace)
                     procMesh % faces % jMap(:,:,iFaceLocal)       = mesh % faces % jMap(:,:,iFace)
                  
                     nBe = nBe + 1
                     procMesh % faces % boundaryID(iFaceLocal) = nBe

                     ! Boundary ID's associated with MPI boundaries are reported as positive integers
                     ! so that they can be distinguished from physical boundaries that do not require
                     ! communication with neighboring processes. 
                     bCom(procID) % boundaryIDs(nBe) = iFaceLocal
                     bCom(procID) % boundaryGlobalIDs(nBe) = iFace
                     bCom(procID) % extProcIDs(nBe)  = p1


                  ENDIF
               ENDIF

            ENDDO
            
! ----------------------------------------------------------------------------- !

            WRITE( pIDChar, '(I4.4)' ) procID
#ifdef TECPLOT
            CALL procMesh % WriteTecplot( 'mesh.'//pIDChar )
#endif
            CALL procMesh % WriteSELFMeshFile( TRIM(params % SELFMeshFile)//'.'//pIDChar )
            
            CALL procMesh % Trash( )

         ENDDO
   


      
      DO procID = 0, nProc-1

         nlocMPI = 0
         DO bID = 1, bCom(procID) % nBoundaries

            extProc = bCom(procID) % extProcIDs(bID)
            IF( extProc /= procID )THEN

               ! nMPI is an array that keeps track of the order in which MPI messages will be sent
               ! to procID's neighbors
               nlocMPI(extProc) = nlocMPI(extProc) + 1
               globalFaceID  = bCom(procID) % boundaryGlobalIDs(bID)

               IF( faceProcOwners(globalFaceID, 1) == procID )THEN

                  extBID = faceBoundaryIDs(globalFaceID,2)

               ELSEIF( faceProcOwners(globalFaceID,2) == procID )THEN

                  extBID = faceBoundaryIDs(globalFaceID,1)

               ELSE
                  PRINT*, '  SetupCommTables : Something catastrophic happened !'
                  setupSuccess = .FALSE.
                  RETURN
               ENDIF


               bCom(extProc) % unPackMap(extBID) = nlocMPI(extProc)
                  

            ENDIF

         ENDDO

      ENDDO   

      DO procID = 0, nProc-1
         WRITE( pIDChar, '(I4.4)' ) procID
         CALL bCom(procID) % WritePickup( 'ExtComm.'//pIDChar )
      ENDDO
         
! ----------------------------------------------------------------------------- !

      DO procID = 0, nProc-1
         CALL bCom(procID) % Trash( )
      ENDDO
      DEALLOCATE( bCom )
      

      ! For visualizing the decomposition
      materials = REAL( partitions, prec )
#ifdef TECPLOT
      CALL mesh % WriteMaterialTecplot( materials )
#endif
      

      ! Clean up memory !
      DEALLOCATE( materials, &
                  faceProcCount, &
                  faceProcOwners, &
                  faceBoundaryIDs, &
                  nLocMPI )
      
      CALL mesh % Trash( )
      CALL nodal % Trash( )

      ENDIF

 END SUBROUTINE StructuredMeshGenerator_3D

!
 FUNCTION NodesAreTheSame( nlist, mlist, N ) RESULT( NATS )
   IMPLICIT NONE
   INTEGER :: N
   INTEGER :: nlist(1:N), mlist(1:N)
   LOGICAL :: NATS
   ! LOCAL
   INTEGER :: nlistSorted(1:N), mlistSorted(1:N), i

      CALL InsertionSort( nlist, nlistSorted, N )
      CALL InsertionSort( mlist, mlistSorted, N )

      NATS = .TRUE.

      DO i = 1, N
         IF( .NOT.( nlistSorted(i) == mlistSorted(i) ) )THEN
            NATS = .FALSE.
            EXIT
         ENDIF
      ENDDO

 END FUNCTION NodesAreTheSame

 END MODULE HexMesh_Class
