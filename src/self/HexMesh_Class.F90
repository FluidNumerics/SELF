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
  USE Boundaries_Class
  USE Geom_EquationParser_Class

  USE HDF5

  IMPLICIT NONE


#include "self_macros.h"

! HexMesh
!  The HexMesh DATA structure defines attributes needed to describe a conformal unstructured
!  spectral element mesh.
!
!  The HexMesh DATA-structure brings together nodes, elements, faces, and mapped-geometry
!  DATA-structures to completely describe a 3-D spectral element mesh. TYPE-bound procedures are
!  provided for filling in mesh connectivity information including element-to-element neighbors,
!  edge-to-element IDs, element-to-node IDs, and edge-to-node IDs.
!

  TYPE MeshObjectList
    INTEGER              :: nNodes, nFaces, nElements
    INTEGER, ALLOCATABLE :: nodeids(:)
    INTEGER, ALLOCATABLE :: faceids(:)
    INTEGER, ALLOCATABLE :: elementids(:)
  END TYPE MeshObjectList

  TYPE DomainDecomposition
    INTEGER                           :: nBlocks, nGlobalElements
    TYPE(MeshObjectList), ALLOCATABLE :: mesh_obj(:) 
    INTEGER, ALLOCATABLE              :: element_to_blockID(:) 
    INTEGER, ALLOCATABLE              :: global_to_localID(:) 

    CONTAINS
    PROCEDURE, PRIVATE :: Build => Build_DomainDecomposition

  END TYPE DomainDecomposition

  TYPE HexMesh
    TYPE( HexElements )         :: elements
    TYPE( Nodes  )              :: nodes
    TYPE( Edges )               :: edges 
    TYPE( Faces )               :: faces
    TYPE( Boundaries )          :: boundaryMap 
    TYPE( DomainDecomposition ) :: decomp 
    INTEGER                     :: cornerMap(1:3,1:8)
    INTEGER                     :: sideMap(1:6)
    INTEGER                     :: faceMap(1:4,1:6)
    INTEGER                     :: edgeMap(1:2,1:12)
    INTEGER                     :: edgeFaceMap(1:2,1:4)

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

    PROCEDURE :: Load_SELFMesh

    ! Mesh Format readers/writers
    PROCEDURE          :: Write_TecplotMesh
    PROCEDURE, PRIVATE :: Read_SELFMesh
    PROCEDURE, PRIVATE :: Write_SELFMesh
    PROCEDURE, PRIVATE :: Read_UCDMesh
    PROCEDURE, PRIVATE :: Read_TrellisUCDMesh
  
    ! Support Routines
    PROCEDURE, PRIVATE :: ConstructStructuredMesh
    PROCEDURE, PRIVATE :: Write_MeshElements
    PROCEDURE, PRIVATE :: Write_MeshFaces
    PROCEDURE, PRIVATE :: Write_MeshNodes
    PROCEDURE, PRIVATE :: Write_MeshGeometry
    PROCEDURE, PRIVATE :: Write_MeshDecomp
    PROCEDURE, PRIVATE :: Read_MeshDecomp
    PROCEDURE, PRIVATE :: Read_MeshElements
    PROCEDURE, PRIVATE :: ConstructFaces               
    PROCEDURE, PRIVATE :: ConstructStructuredFaces     
    PROCEDURE, PRIVATE :: ConstructDoublyPeriodicFaces 
    PROCEDURE, PRIVATE :: ConstructElementNeighbors    
    PROCEDURE, PRIVATE :: DetermineOrientation         
    PROCEDURE, PRIVATE :: ScaleTheMesh                 
    PROCEDURE, PRIVATE :: SetupFaceMaps

    ! Visualization I/O Routines

  END TYPE HexMesh

  INTEGER, PRIVATE, PARAMETER    :: BoundaryFlagStructured = NO_NORMAL_FLOW
  PRIVATE :: MeshDecompose
  PRIVATE :: SetupProcessBoundaryMap
  PRIVATE :: CreateMeshObjLists
  PRIVATE :: StructuredMeshToBlocks
  PRIVATE :: StructuredDecompose
  PRIVATE :: UnstructuredDecompose
  PRIVATE :: Add_FloatMeshObj_to_HDF5
  PRIVATE :: Add_IntMeshObj_to_HDF5
  PRIVATE :: SetGlobalToLocalMapping
  PRIVATE :: Get_HDF5_Obj_Dimensions

CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
  SUBROUTINE Build_HexMesh( myHexMesh, nNodes, nElements, nFaces, N )
#undef __FUNC__
#define __FUNC__ "ConstructFaces"
    IMPLICIT NONE
    CLASS(HexMesh), INTENT(inout) :: myHexMesh
    INTEGER, INTENT(in)           :: nNodes, nElements, nFaces, N

    INFO('Start')
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


    INFO('End')

  END SUBROUTINE Build_HexMesh
!
  SUBROUTINE Trash_HexMesh( myHexMesh )

    IMPLICIT NONE
    CLASS(HexMesh), INTENT(inout) :: myHexMesh
    ! Local
    INTEGER :: i 

    CALL myHexMesh % elements % Trash( )
    CALL myHexmesh % nodes % Trash( )
    CALL myHexMesh % faces % Trash( )

    CALL myHexMesh % boundaryMap % Trash( )

#ifdef HAVE_CUDA
    DEALLOCATE( myHexMesh % cornerMap_dev, &
                myHexMesh % sideMap_dev, &
                myHexMesh % faceMap_dev, &
                myHexMesh % edgeFaceMap_dev, &
                myHexMesh % edgeMap_dev )
#endif

    IF( ALLOCATED(myHexMesh % decomp % mesh_obj) )THEN
      DO i = 0, myHexMesh % decomp % nBlocks-1
        IF( ALLOCATED(myHexMesh % decomp % mesh_obj(i) % nodeids) )DEALLOCATE(myHexMesh % decomp % mesh_obj(i) % nodeids)
        IF( ALLOCATED(myHexMesh % decomp % mesh_obj(i) % faceids) )DEALLOCATE(myHexMesh % decomp % mesh_obj(i) % faceids)
        IF( ALLOCATED(myHexMesh % decomp % mesh_obj(i) % elementids) )DEALLOCATE(myHexMesh % decomp % mesh_obj(i) % elementids)
      ENDDO
      DEALLOCATE(myHexMesh % decomp % mesh_obj)
    ENDIF
    IF( ALLOCATED(myHexMesh % decomp % element_to_blockID) )DEALLOCATE(myHexMesh % decomp % element_to_blockID)
    IF( ALLOCATED(myHexMesh % decomp % global_to_localID) )DEALLOCATE(myHexMesh % decomp % global_to_localID)

  END SUBROUTINE Trash_HexMesh
!
  SUBROUTINE  Build_DomainDecomposition( decomp, nGlobalElements, nBlocks, nLocalElements, nLocalFaces, nLocalNodes )

    IMPLICIT NONE
    CLASS( DomainDecomposition ), INTENT(out) :: decomp
    INTEGER, INTENT(in)                       :: nGlobalElements, nBlocks
    INTEGER, INTENT(in)                       :: nLocalElements(0:nBlocks-1), nLocalFaces(0:nBlocks-1), nLocalNodes(0:nBlocks-1)
    ! Local
    INTEGER :: i 

     decomp % nGlobalElements = nGlobalElements
     decomp % nBlocks         = nBlocks
     ALLOCATE( decomp % element_to_blockID(1:nGlobalElements) )
     ALLOCATE( decomp % global_to_localID(1:nGlobalElements) )
     decomp % element_to_blockID = 0
     decomp % global_to_localID = 0

     ALLOCATE( decomp % mesh_obj(0:nBlocks-1) )
     DO i = 0, nBlocks-1
       decomp % mesh_obj(i) % nNodes    = nLocalNodes(i)
       decomp % mesh_obj(i) % nElements = nLocalElements(i)
       decomp % mesh_obj(i) % nFaces    = nLocalFaces(i)

       ALLOCATE( decomp % mesh_obj(i) % nodeids(1:nLocalNodes(i)) )
       ALLOCATE( decomp % mesh_obj(i) % elementids(1:nLocalElements(i)) )
       ALLOCATE( decomp % mesh_obj(i) % faceids(1:nLocalFaces(i)) )

       decomp % mesh_obj(i) % nodeids(1:nLocalNodes(i)) = 0
       decomp % mesh_obj(i) % elementids(1:nLocalElements(i)) = 0
       decomp % mesh_obj(i) % faceids(1:nLocalFaces(i)) = 0
     ENDDO
 
  END SUBROUTINE  Build_DomainDecomposition
!
#ifdef HAVE_CUDA
  SUBROUTINE UpdateDevice_HexMesh( myHexMesh )
    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh

    CALL myHexMesh % faces % UpdateDevice( )
    CALL myHexMesh % elements % UpdateDevice( )
    CALL myHexMesh % nodes % UpdateDevice( )
    CALL myHexMesh % boundaryMap % UpdateDevice()

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
    IMPLICIT NONE
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
  SUBROUTINE SetupFaceMaps( mesh )
#undef __FUNC__
#define __FUNC__ "ConstructFaces"
  CLASS( HexMesh ), INTENT(inout) :: mesh
  ! Local
  INTEGER :: i, j, k, ii, jj

  INFO('Start')
    DO k = 1, mesh % faces % nFaces
      DO j = 0, mesh % faces % N
        DO i = 0, mesh % faces % N

          IF( i == 0 )THEN
            IF( j == 0 )THEN
              ii = (1-mesh % faces % swapDimensions(k))*&
                (mesh % faces % iStart(k)) + &
                (mesh % faces % swapDimensions(k))*&
                (mesh % faces % jStart(k))
              jj = (1-mesh % faces % swapDimensions(k))*&
                (mesh % faces % jStart(k)) + &
                (mesh % faces % swapDimensions(k))*&
                (mesh % faces % iStart(k))
            ELSE
              ii = mesh % faces % swapDimensions(k)*&
                (ii+mesh % faces % jInc(k)) + &
                (1-mesh % faces % swapDimensions(k))*&
                mesh % faces % iStart(k)
              jj = (1-mesh % faces % swapDimensions(k))*&
                (jj+mesh % faces % jInc(k)) +&
                mesh % faces % swapDimensions(k)*&
                mesh % faces % jStart(k)
            ENDIF
          ELSE
            ii = (1-mesh % faces % swapDimensions(k))*&
              (ii + mesh % faces % iInc(k)) +&
              mesh % faces % swapDimensions(k)*ii
            jj = mesh % faces % swapDimensions(k)*&
              (jj+mesh % faces % iInc(k)) + &
              (1-mesh % faces % swapDimensions(k))*jj
          ENDIF

          mesh % faces % iMap(i,j,k) = ii
          mesh % faces % jMap(i,j,k) = jj

        ENDDO
      ENDDO
    ENDDO

  INFO('End')

  END SUBROUTINE SetupFaceMaps
!
  SUBROUTINE ConstructFaces( myHexMesh )
#undef __FUNC__
#define __FUNC__ "ConstructFaces"
    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh
    ! LOCAL
    TYPE( KeyRing ) :: KeyCabinet(1:myHexMesh % nodes % nNodes)
    INTEGER :: nEls, nNodes, iEl, nFaces, i, j, k, ii, jj
    INTEGER :: locNodeIDs(1:4), globNodeIDs(1:4)
    INTEGER :: keyRingID, globFaceID
    INTEGER :: N
    LOGICAL :: keyExists

    INFO('Start')

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

    CALL myHexMesh % SetupFaceMaps( )

#ifdef HAVE_CUDA

    CALL myHexMesh % faces % UpdateDevice( )

#endif
    INFO('End')

  END SUBROUTINE ConstructFaces
!
  SUBROUTINE ConstructStructuredFaces( myHexMesh, nXElem, nYElem, nZElem, boundaryconditionflags )
#undef __FUNC__
#define __FUNC__ "ConstructStructuredFaces"
    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh
    INTEGER, INTENT(in) :: nXElem, nYElem, nZElem
    INTEGER, INTENT(in), OPTIONAL :: boundaryconditionflags(1:6)
    ! LOCAL
    INTEGER ::  e1, e2, s1, s2, nFaces, i, j, k, l, IFace, ii, jj
    INTEGER :: bcflags(1:6)

    INFO('Start')

    IF( PRESENT(boundaryconditionflags) )THEN
      bcflags = boundaryconditionflags
    ELSE
      bcflags = NO_NORMAL_FLOW
    ENDIF

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
            e2 = bcFlags(1) ! Enforce boundary condition on south boundary of the domain
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
            e2 = bcFlags(4) ! Enforce boundary condition on west domain boundary
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
            e2 = bcFlags(5) ! Enforce boundary condition on bottom domain boundary
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
        e2 = bcFlags(2) ! Enforce boundary condition on east domain boundary
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
        e2 = bcFlags(3) ! Enforce boundary condition on domain north boundary
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
        e2 = bcFlags(6) ! Enforce boundary condition on domain top
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

    CALL myHexMesh % SetupFaceMaps( )

#ifdef HAVE_CUDA

    CALL myHexMesh % faces % UpdateDevice( )

#endif
    INFO('End')

  END SUBROUTINE ConstructStructuredFaces
!  
  SUBROUTINE ConstructDoublyPeriodicFaces( myHexMesh, nXElem, nYElem, nZElem, boundaryconditionflags )
#undef __FUNC__
#define __FUNC__ "ConstructDoublyPeriodicFaces"
    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout) :: myHexMesh
    INTEGER, INTENT(in) :: nXElem, nYElem, nZElem
    INTEGER, INTENT(in), OPTIONAL :: boundaryconditionflags(1:6)
    ! LOCAL
    INTEGER ::  e1, e2, s1, s2, nFaces, i, j, k, l, IFace, ii, jj
    INTEGER :: bcflags(1:6)

    INFO('Start')

    IF( PRESENT(boundaryconditionflags) )THEN
      bcflags = boundaryconditionflags
    ELSE
      bcflags = NO_NORMAL_FLOW
      bcflags(6) = PRESCRIBED
    ENDIF

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
            e2 = bcFlags(5)
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
        e2 = bcFlags(6)
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

    CALL myHexMesh % SetupFaceMaps( )

#ifdef HAVE_CUDA

    CALL myHexMesh % faces % UpdateDevice( )

#endif
    INFO('End')

  END SUBROUTINE ConstructDoublyPeriodicFaces
!
  SUBROUTINE DetermineOrientation( myHexMesh, faceID, secondaryNodes )

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


  END SUBROUTINE DetermineOrientation
!
  SUBROUTINE ConstructElementNeighbors( myHexMesh )

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


  END SUBROUTINE ConstructElementNeighbors
!
  SUBROUTINE ScaleTheMesh( myHexMesh, interp, xScale, yScale, zScale  )

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

  END SUBROUTINE ScaleTheMesh

  SUBROUTINE ConstructStructuredMesh( myHexMesh, interp, nXelem, nYelem, nZelem, geomparser  )
#undef __FUNC__
#define __FUNC__ "ConstructStructuredMesh"
    IMPLICIT NONE
    CLASS( HexMesh ), INTENT(inout)  :: myHexMesh
    TYPE( Lagrange ), INTENT(in)     :: interp
    INTEGER, INTENT(in)              :: nXelem, nYelem, nZelem
    TYPE( Geom_EquationParser ), INTENT(in) :: geomparser
    ! LOGICAL
    TYPE( Surfaces ) :: boundSurfs
    REAL(prec) :: x, y, z, zb, zi, zu, zip1, dxElem, dyElem, dzElem
    REAL(prec) :: x1(1:3), x2(1:3), x3(1:3), x4(1:3)
    REAL(prec) :: c1(1:3), c2(1:3), xE(1:3)
    REAL(prec), ALLOCATABLE :: xc(:,:,:,:), s(:), weights(:)
    LOGICAL :: DoublyPeriodic
    INTEGER :: nNodes, nElements, nFaces, gPolyDeg, nSurf
    INTEGER :: nodes(1:8)
    INTEGER :: n1, n2, n3, n4
    INTEGER :: iNode, iEl, iSide, iX, iY, iZ, i, j, iSurf

    INFO('Start')
    IF( geomparser % boundaryconditionflags(1) == PERIODIC ) DoublyPeriodic = .TRUE.
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
    ! curves for the element boundaryMap

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
      CALL myHexMesh % ConstructDoublyPeriodicFaces( nXElem, nYElem, nZElem, geomparser % boundaryconditionflags )
    ELSE
      CALL myHexMesh % ConstructStructuredFaces( nXElem, nYElem, nZElem, geomparser % boundaryconditionflags )
    ENDIF


    CALL myHexMesh % ConstructElementNeighbors( )

    ! Clear up memory
    DEALLOCATE( s, xc, weights )

    CALL boundSurfs % Trash( )

#ifdef HAVE_CUDA
    CALL myHexMesh % UpdateDevice( )
#endif
    INFO('End')

  END SUBROUTINE ConstructStructuredMesh
!
!==================================================================================================!
!--------------------------------- Mesh File I/O Routines -----------------------------------------!
!==================================================================================================!
!
!
!  SUBROUTINE ReadSELFMeshFile_HexMesh( myHexMesh, filename )
!
!#undef __FUNC__
!#define __FUNC__ "ReadSELFMeshFile"
!
!    CLASS( HexMesh ), INTENT(out)   :: myHexMesh
!    CHARACTER(*), INTENT(in)        :: filename
!    ! LOCAL
!    INTEGER :: nNodes, nElements, nFaces, N, nBoundaries
!    INTEGER :: IFace, iNode, iEl
!    INTEGER :: fUnit, k, i, j, l, row, col, ii, jj
!#ifdef HAVE_CUDA
!    INTEGER :: istat
!#endif
!
!
!    INFO('Start')
!    INFO('Reading '//TRIM( filename )//'.mesh')
!    ! Get a new file unit
!    OPEN( UNIT    = NEWUNIT(fUnit), &
!      FILE    = TRIM( filename )//'.mesh', &
!      FORM    = 'UNFORMATTED',&
!      STATUS  = 'OLD', &
!      ACCESS  = 'DIRECT', &
!      CONVERT = 'BIG_ENDIAN', &
!      RECL    = SIZEOF(nNodes) ) ! How to DO variable record length
!
!    ! ---- Gather the number of nodes, number of elements, and number of edges ---- !
!    k = 1
!    READ( fUnit, rec=k )nNodes
!    k = k+1
!    READ( fUnit, rec=k )nElements
!    k = k+1
!    READ( fUnit, rec=k )nFaces
!    k = k+1
!    READ( fUnit, rec=k )N
!    k = k+1
!    READ( fUnit, rec=k )nBoundaries
!    k = k+1
!
!    ! ---- Build the quadrature mesh (empty) ---- !
!    CALL myHexMesh % Build( nNodes, nElements, nFaces, N )
!    CALL myHexMesh % boundaryMap % Build( nBoundaries, 0 )
!
!    ! ---- Read in the element connectivity ---- !
!    DO iEl = 1, nElements
!      READ( fUnit, rec=k ) myHexMesh % elements % elementID(iEl)
!      k = k+1
!      DO i = 1, 8
!        READ( fUnit, rec=k ) myHexMesh % elements % nodeIDs(i,iEl)
!        k = k+1
!      ENDDO
!    ENDDO
!
!    ! ---- Read in the face information ---- !
!
!    DO IFace = 1, nFaces
!      READ( fUnit, rec=k ) myHexMesh % faces % faceID(IFace)
!      k = k+1
!      READ( fUnit, rec=k ) myHexMesh % faces % boundaryID(IFace)
!      k = k+1
!      DO i = 1, 4
!        READ( fUnit, rec=k ) myHexMesh % faces % nodeIDs(i,IFace)
!        k = k+1
!      ENDDO
!      DO i = 1, 2
!        READ( fUnit, rec=k ) myHexMesh % faces % elementIDs(i,IFace)
!        k = k+1
!        READ( fUnit, rec=k ) myHexMesh % faces % elementSides(i,IFace)
!        k = k+1
!      ENDDO
!      READ( fUnit, rec=k ) myHexMesh % faces % iStart(IFace)
!      k = k+1
!      READ( fUnit, rec=k ) myHexMesh % faces % iInc(IFace)
!      k = k+1
!      READ( fUnit, rec=k ) myHexMesh % faces % jStart(IFace)
!      k = k+1
!      READ( fUnit, rec=k ) myHexMesh % faces % jInc(IFace)
!      k = k+1
!      READ( fUnit, rec=k ) myHexMesh % faces % swapDimensions(IFace)
!      k = k+1
!    ENDDO
!
!    DO iface = 1, myHexMesh % boundaryMap % nBoundaries
!      READ( fUnit, rec=k ) myHexMesh % boundaryMap % extProcIDs(iface)
!      k = k+1
!      READ( fUnit, rec=k ) myHexMesh % boundaryMap % boundaryIDs(iface)
!      k = k+1
!      READ( fUnit, rec=k ) myHexMesh % boundaryMap % boundaryGlobalIDs(iface)
!      k = k+1
!      READ( fUnit, rec=k ) myHexMesh % boundaryMap % boundaryCondition(iface)
!      k = k+1
!    ENDDO
!
!
!    CLOSE( fUnit )
!    ! Get a new file unit
!    OPEN( UNIT    = NEWUNIT(fUnit), &
!      FILE    = TRIM( filename )//'.geom', &
!      FORM    = 'UNFORMATTED',&
!      STATUS  = 'OLD', &
!      ACCESS  = 'DIRECT', &
!      CONVERT = 'BIG_ENDIAN', &
!      RECL    = prec ) ! How to DO variable record length
!
!    ! ---- Read in the corner nodes ---- !
!    k = 1
!    DO iNode = 1, nNodes  ! Loop over the nodes in the file
!      READ( fUnit, rec=k ) myHexMesh % nodes % x(1,iNode)
!      k = k+1
!      READ( fUnit, rec=k ) myHexMesh % nodes % x(2,iNode)
!      k = k+1
!      READ( fUnit, rec=k ) myHexMesh % nodes % x(3,iNode)
!      k = k+1
!    ENDDO
!
!    ! ---- Read in the element information ---- !
!    DO iEl = 1, nElements
!      DO l = 0, N
!        DO j = 0, N
!          DO i = 0, N
!            READ( fUnit, rec=k ) myHexMesh % elements % x(i,j,l,1,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % x(i,j,l,2,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % x(i,j,l,3,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % dxds(i,j,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % dxdp(i,j,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % dxdq(i,j,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % dyds(i,j,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % dydp(i,j,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % dydq(i,j,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % dzds(i,j,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % dzdp(i,j,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % dzdq(i,j,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % J(i,j,l,iEl)
!            k = k+1
!            DO col = 1, 3
!              DO row = 1, 3
!                READ( fUnit, rec=k ) myHexMesh % elements % Ja(i,j,l,row,col,iEl)
!                k = k+1
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO
!      DO l = 1, 6
!        DO j = 0, N
!          DO i = 0, N
!            READ( fUnit, rec=k ) myHexMesh % elements % boundaryLengthScale(i,j,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % xBound(i,j,1,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % xBound(i,j,2,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % xBound(i,j,3,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % nHat(1,i,j,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % nHat(2,i,j,l,iEl)
!            k = k+1
!            READ( fUnit, rec=k ) myHexMesh % elements % nHat(3,i,j,l,iEl)
!            k = k+1
!          ENDDO
!        ENDDO
!      ENDDO
!    ENDDO
!
!    CLOSE( fUnit )
!
!    CALL myHexMesh % ConstructElementNeighbors( )
!
!    DO k = 1, myHexMesh % faces % nFaces
!      DO j = 0, N
!        DO i = 0, N
!
!          IF( i == 0 )THEN
!            IF( j == 0 )THEN
!              ii = (1-myHexMesh % faces % swapDimensions(k))*&
!                (myHexMesh % faces % iStart(k)) + &
!                (myHexMesh % faces % swapDimensions(k))*&
!                (myHexMesh % faces % jStart(k))
!              jj = (1-myHexMesh % faces % swapDimensions(k))*&
!                (myHexMesh % faces % jStart(k)) + &
!                (myHexMesh % faces % swapDimensions(k))*&
!                (myHexMesh % faces % iStart(k))
!            ELSE
!              ii = myHexMesh % faces % swapDimensions(k)*&
!                (ii+myHexMesh % faces % jInc(k)) + &
!                (1-myHexMesh % faces % swapDimensions(k))*&
!                myHexMesh % faces % iStart(k)
!              jj = (1-myHexMesh % faces % swapDimensions(k))*&
!                (jj+myHexMesh % faces % jInc(k)) +&
!                myHexMesh % faces % swapDimensions(k)*&
!                myHexMesh % faces % jStart(k)
!            ENDIF
!          ELSE
!            ii = (1-myHexMesh % faces % swapDimensions(k))*&
!              (ii + myHexMesh % faces % iInc(k)) +&
!              myHexMesh % faces % swapDimensions(k)*ii
!            jj = myHexMesh % faces % swapDimensions(k)*&
!              (jj+myHexMesh % faces % iInc(k)) + &
!              (1-myHexMesh % faces % swapDimensions(k))*jj
!          ENDIF
!
!          myHexMesh % faces % iMap(i,j,k) = ii
!          myHexMesh % faces % jMap(i,j,k) = jj
!
!        ENDDO
!      ENDDO
!    ENDDO
!
!#ifdef HAVE_CUDA
!    CALL myHexMesh % UpdateDevice( )
!    istat = cudaDeviceSynchronize( )
!#endif
!    INFO('End')
!
!  END SUBROUTINE ReadSELFMeshFile_HexMesh
!!
!  SUBROUTINE WriteSELFMeshFile_HexMesh( myHexMesh, filename )
!
!    IMPLICIT NONE
!    CLASS( HexMesh ), INTENT(in)   :: myHexMesh
!    CHARACTER(*), INTENT(in)       :: filename
!    ! LOCAL
!    INTEGER :: nNodes, nElements, nFaces, N
!    INTEGER :: IFace, iNode, iEl
!    INTEGER :: fUnit, k, i, j, l, row, col
!
!
!    nNodes = 1
!    ! Get a new file unit
!    OPEN( UNIT    = NEWUNIT(fUnit), &
!      FILE    = TRIM( filename )//'.mesh', &
!      FORM    = 'UNFORMATTED',&
!      STATUS  = 'REPLACE', &
!      ACCESS  = 'DIRECT', &
!      CONVERT = 'BIG_ENDIAN', &
!      RECL    = SIZEOF(nNodes) ) ! How to DO variable record length
!
!
!    ! ---- Gather the number of nodes, number of elements, and number of edges ---- !
!    nNodes = myHexMesh % nodes % nNodes
!    nElements = myHexMesh % elements % nElements
!    nFaces = myHexMesh % faces % nFaces
!    N      = myHexMesh % elements % N
!    k = 1
!    WRITE( fUnit, rec=k )nNodes
!    k = k+1
!    WRITE( fUnit, rec=k )nElements
!    k = k+1
!    WRITE( fUnit, rec=k )nFaces
!    k = k+1
!    WRITE( fUnit, rec=k )N
!    k = k+1
!    WRITE( fUnit, rec=k )myHexMesh % boundaryMap % nBoundaries
!    k = k+1
!
!
!    ! ---- Read in the element connectivity ---- !
!    DO iEl = 1, nElements
!      WRITE( fUnit, rec=k ) myHexMesh % elements % elementID(iEl)
!      k = k+1
!      DO i = 1, 8
!        WRITE( fUnit, rec=k ) myHexMesh % elements % nodeIDs(i,iEl)
!        k = k+1
!      ENDDO
!    ENDDO
!
!    ! ---- Read in the face information ---- !
!
!    DO IFace = 1, nFaces
!      WRITE( fUnit, rec=k ) myHexMesh % faces % faceID(IFace)
!      k = k+1
!      WRITE( fUnit, rec=k ) myHexMesh % faces % boundaryID(IFace)
!      k = k+1
!      DO i = 1, 4
!        WRITE( fUnit, rec=k ) myHexMesh % faces % nodeIDs(i,IFace)
!        k = k+1
!      ENDDO
!      DO i = 1, 2
!        WRITE( fUnit, rec=k ) myHexMesh % faces % elementIDs(i,IFace)
!        k = k+1
!        WRITE( fUnit, rec=k ) myHexMesh % faces % elementSides(i,IFace)
!        k = k+1
!      ENDDO
!      WRITE( fUnit, rec=k ) myHexMesh % faces % iStart(IFace)
!      k = k+1
!      WRITE( fUnit, rec=k ) myHexMesh % faces % iInc(IFace)
!      k = k+1
!      WRITE( fUnit, rec=k ) myHexMesh % faces % jStart(IFace)
!      k = k+1
!      WRITE( fUnit, rec=k ) myHexMesh % faces % jInc(IFace)
!      k = k+1
!      WRITE( fUnit, rec=k ) myHexMesh % faces % swapDimensions(IFace)
!      k = k+1
!    ENDDO
!
!    DO iface = 1, myHexMesh % boundaryMap % nBoundaries
!      WRITE( fUnit, rec=k ) myHexMesh % boundaryMap % extProcIDs(iface)
!      k = k+1
!      WRITE( fUnit, rec=k ) myHexMesh % boundaryMap % boundaryIDs(iface)
!      k = k+1
!      WRITE( fUnit, rec=k ) myHexMesh % boundaryMap % boundaryGlobalIDs(iface)
!      k = k+1
!      WRITE( fUnit, rec=k ) myHexMesh % boundaryMap % boundaryCondition(iface)
!      k = k+1
!    ENDDO
!
!    CLOSE( fUnit )
!    ! Get a new file unit
!    OPEN( UNIT    = NEWUNIT(fUnit), &
!      FILE    = TRIM( filename )//'.geom', &
!      FORM    = 'UNFORMATTED',&
!      STATUS  = 'REPLACE', &
!      ACCESS  = 'DIRECT', &
!      CONVERT = 'BIG_ENDIAN', &
!      RECL    = prec ) ! How to DO variable record length
!
!    ! ---- Read in the corner nodes ---- !
!    k = 1
!    DO iNode = 1, nNodes  ! Loop over the nodes in the file
!      WRITE( fUnit, rec=k ) myHexMesh % nodes % x(1,iNode)
!      k = k+1
!      WRITE( fUnit, rec=k ) myHexMesh % nodes % x(2,iNode)
!      k = k+1
!      WRITE( fUnit, rec=k ) myHexMesh % nodes % x(3,iNode)
!      k = k+1
!    ENDDO
!
!    ! ---- Read in the element information ---- !
!    DO iEl = 1, nElements
!      DO l = 0, N
!        DO j = 0, N
!          DO i = 0, N
!            WRITE( fUnit, rec=k ) myHexMesh % elements % x(i,j,l,1,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % x(i,j,l,2,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % x(i,j,l,3,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % dxds(i,j,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % dxdp(i,j,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % dxdq(i,j,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % dyds(i,j,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % dydp(i,j,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % dydq(i,j,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % dzds(i,j,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % dzdp(i,j,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % dzdq(i,j,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % J(i,j,l,iEl)
!            k = k+1
!            DO col = 1, 3
!              DO row = 1, 3
!                WRITE( fUnit, rec=k ) myHexMesh % elements % Ja(i,j,l,row,col,iEl)
!                k = k+1
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO
!      DO l = 1, 6
!        DO j = 0, N
!          DO i = 0, N
!            WRITE( fUnit, rec=k ) myHexMesh % elements % boundaryLengthScale(i,j,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % xBound(i,j,1,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % xBound(i,j,2,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % xBound(i,j,3,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % nHat(1,i,j,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % nHat(2,i,j,l,iEl)
!            k = k+1
!            WRITE( fUnit, rec=k ) myHexMesh % elements % nHat(3,i,j,l,iEl)
!            k = k+1
!          ENDDO
!        ENDDO
!      ENDDO
!    ENDDO
!
!    CLOSE( fUnit )
!
!  END SUBROUTINE WriteSELFMeshFile_HexMesh

  SUBROUTINE Read_TrellisUCDMesh( myHexMesh, interp, filename )

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

  END SUBROUTINE Read_TrellisUCDMesh
!
  SUBROUTINE Read_UCDMesh( myHexMesh, interp, filename )
#undef __FUNC__
#define __FUNC__ "Read_UCDMesh"
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


    INFO('Start')
    INFO('Mesh File : '//TRIM( filename ))

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

    INFO('End')

  END SUBROUTINE Read_UCDMesh
!
  SUBROUTINE Write_TecplotMesh( myHexMesh, filename )

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

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "Jacobian", "dxds", "dxdp", "dxdq", "dyds", "dydp", "dydq", "dzds", "dzdp", "dzdq", "process-id" '


    DO iEl = 1, myHexMesh % elements % nElements

      eID = myHexMesh % elements % elementID(iEl)
      WRITE(zoneID,'(I7.7)') eID
      WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',N+1,', J=',N+1,', K=',N+1,',F=POINT'

      DO k = 0, N
        DO j = 0, N
          DO i = 0,N

            WRITE(fUnit,'(13(E15.7,1x), F10.4)')  myHexMesh % elements % x(i,j,k,1,iEl), &
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
              myHexMesh % elements % dzdq(i,j,k,iEl), &
              FLOAT(myHexMesh % decomp % element_to_blockID(iEl))
          ENDDO
        ENDDO
      ENDDO

    ENDDO

    CLOSE(UNIT=fUnit)

    RETURN

 END SUBROUTINE Write_TecplotMesh
!
 SUBROUTINE Read_SELFMesh( mesh, meshfile, polydeg, my_RankID, nMPI_Ranks, mpiCommunicator )
#undef __FUNC__
#define __FUNC__ "Read_SELFMesh"
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: mesh 
   CHARACTER(*), INTENT(in)        :: meshfile
   INTEGER, INTENT(in)             :: polydeg
   INTEGER, INTENT(in)             :: my_RankID
   INTEGER, INTENT(in)             :: nMPI_Ranks
   INTEGER, INTENT(in)             :: mpiCommunicator
   ! Local
   INTEGER(HID_T) :: file_id, memspace, filespace, group_id
   INTEGER(HID_T) :: plist_id
   INTEGER        :: error
  
   INFO('Start')

   CALL h5open_f(error)  

#ifdef HAVE_MPI
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, mpiCommunicator, MPI_INFO_NULL, error)
    CALL h5fopen_f(TRIM(meshfile), H5F_ACC_RDWR_F, file_id, error, access_prp=plist_id)
    CALL h5pclose_f(plist_id,error)
#else
    CALL h5fopen_f(TRIM(meshfile), H5F_ACC_RDWR_F, file_id, error)
#endif

    IF( error /= 0 )THEN
      INFO('HDF5 file open failed')
      STOP
    ENDIF

    CALL mesh % Read_MeshDecomp( file_id, my_RankID )

    CALL mesh % Build( mesh % decomp % mesh_obj(0) % nNodes, &
                       mesh % decomp % mesh_obj(0) % nElements, &
                       mesh % decomp % mesh_obj(0) % nFaces, polydeg )

    CALL mesh % Read_MeshElements( file_id )

STOP



   INFO('End')


 END SUBROUTINE Read_SELFMesh

 SUBROUTINE Write_SELFMesh( global_mesh, meshfile, nMPI_Ranks )
#undef __FUNC__
#define __FUNC__ "Write_SELFMesh"
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: global_mesh 
   CHARACTER(*), INTENT(in)        :: meshfile
   INTEGER, INTENT(in)             :: nMPI_Ranks
   ! Local 
   INTEGER(HID_T)   :: file_id, memspace, filespace, group_id
   INTEGER          :: error

     INFO('Start')

       CALL h5open_f(error)  
  
       CALL h5fcreate_f(TRIM(meshfile), H5F_ACC_TRUNC_F, file_id, error)

       CALL h5gcreate_f( file_id, "/mesh", group_id, error )
       CALL h5gclose_f( group_id, error )

       CALL h5gcreate_f( file_id, "/mesh/global", group_id, error )
       CALL h5gclose_f( group_id, error )

       CALL global_mesh % Write_MeshElements( file_id )
       CALL global_mesh % Write_MeshFaces( file_id )
       CALL global_mesh % Write_MeshNodes( file_id )
       CALL global_mesh % Write_MeshGeometry( file_id )
       CALL global_mesh % Write_MeshDecomp( file_id )

       CALL h5fclose_f( file_id, error )
       CALL h5close_f( error )

     INFO('End')


 END SUBROUTINE Write_SELFMesh

 SUBROUTINE Write_MeshElements( global_mesh, file_id )
#undef __FUNC__
#define __FUNC__ "Write_MeshElements"
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: global_mesh 
   INTEGER(HID_T), INTENT(in)      :: file_id
   !
   INTEGER(HSIZE_T) :: dimensions(1:2)
   INTEGER(HID_T)   :: group_id
   INTEGER          :: error


   INFO('Start')
       CALL h5gcreate_f( file_id, "/mesh/global/elements", group_id, error )
       CALL h5gclose_f( group_id, error )

       dimensions(1:2) = (/8, global_mesh % elements % nElements /)  
       CALL Add_IntMeshObj_to_HDF5( rank=2,&
                                    dimensions=dimensions(1:2),&
                                    variable_name='/mesh/global/elements/node-ids', &
                                    int_variable=global_mesh % elements % nodeIDs,&
                                    file_id=file_id )

   INFO('End')

 END SUBROUTINE Write_MeshElements

 SUBROUTINE Read_MeshElements( local_mesh, file_id )
#undef __FUNC__
#define __FUNC__ "Read_MeshElements"
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: local_mesh 
   INTEGER(HID_T), INTENT(in)      :: file_id
   !
   INTEGER(HSIZE_T) :: dimensions(1:2)
   INTEGER(HID_T)   :: dataset_id, filespace
   INTEGER          :: error
   INTEGER          :: node_ids(1:8,1:local_mesh % decomp % nGlobalElements)
   INTEGER          :: iEl, iEl_global 
   CHARACTER(50)    :: msg


   INFO('Start')

       dimensions(1:2) = (/8, local_mesh % decomp % nGlobalElements /)  
       CALL h5dopen_f(file_id, '/mesh/global/elements/node-ids', dataset_id, error)
       CALL h5dget_space_f( dataset_id, filespace, error )
       CALL h5dread_f( dataset_id, H5T_STD_I32LE, node_ids, dimensions, error, H5S_ALL_F, filespace )

       IF( local_mesh % decomp % nBlocks >  1 )THEN

         DO iEl = 1, local_mesh % decomp % mesh_obj(0) % nElements
           iEl_global = local_mesh % decomp % mesh_obj(0) % elementIDs(iEl) 
           local_mesh % elements % nodeIDs(1:8,iEl) = node_ids(1:8,iEl_global)
         ENDDO
  
       ELSE

         local_mesh % elements % nodeIDs = node_ids

       ENDIF
       CALL h5sclose_f( filespace, error)
       CALL h5dclose_f( dataset_id, error )

   INFO('End')

 END SUBROUTINE Read_MeshElements

 SUBROUTINE Write_MeshFaces( global_mesh, file_id )
#undef __FUNC__
#define __FUNC__ "Write_MeshFaces"
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: global_mesh 
   INTEGER(HID_T), INTENT(in)      :: file_id
   !
   INTEGER(HSIZE_T) :: dimensions(1:3)
   INTEGER(HID_T)   :: group_id
   INTEGER          :: error

   INFO('Start')

       CALL h5gcreate_f( file_id, "/mesh/global/faces", group_id, error )
       CALL h5gclose_f( group_id, error )

       dimensions(1) = global_mesh % faces % nFaces
       CALL Add_IntMeshObj_to_HDF5( rank=1,&
                                    dimensions=dimensions(1:1),&
                                    variable_name='/mesh/global/faces/boundary-ids', &
                                    int_variable=global_mesh % faces % boundaryID,&
                                    file_id=file_id )

       CALL Add_IntMeshObj_to_HDF5( rank=1,&
                                    dimensions=dimensions(1:1),&
                                    variable_name='/mesh/global/faces/i-start', &
                                    int_variable=global_mesh % faces % iStart,&
                                    file_id=file_id )

       CALL Add_IntMeshObj_to_HDF5( rank=1,&
                                    dimensions=dimensions(1:1),&
                                    variable_name='/mesh/global/faces/j-start', &
                                    int_variable=global_mesh % faces % jStart,&
                                    file_id=file_id )

       CALL Add_IntMeshObj_to_HDF5( rank=1,&
                                    dimensions=dimensions(1:1),&
                                    variable_name='/mesh/global/faces/i-inc', &
                                    int_variable=global_mesh % faces % iInc,&
                                    file_id=file_id )

       CALL Add_IntMeshObj_to_HDF5( rank=1,&
                                    dimensions=dimensions(1:1),&
                                    variable_name='/mesh/global/faces/j-inc', &
                                    int_variable=global_mesh % faces % jInc,&
                                    file_id=file_id )

       CALL Add_IntMeshObj_to_HDF5( rank=1,&
                                    dimensions=dimensions(1:1),&
                                    variable_name='/mesh/global/faces/swap-dimension', &
                                    int_variable=global_mesh % faces % swapDimensions,&
                                    file_id=file_id )

       dimensions(1:2) = (/4,global_mesh % faces % nFaces /)  
       CALL Add_IntMeshObj_to_HDF5( rank=2,&
                                    dimensions=dimensions(1:2),&
                                    variable_name='/mesh/global/faces/node-ids', &
                                    int_variable=global_mesh % faces % nodeIds,&
                                    file_id=file_id )

       dimensions(1:2) = (/2, global_mesh % faces % nFaces /)  
       CALL Add_IntMeshObj_to_HDF5( rank=2,&
                                    dimensions=dimensions(1:2),&
                                    variable_name='/mesh/global/faces/element-ids', &
                                    int_variable=global_mesh % faces % elementIds,&
                                    file_id=file_id )

       CALL Add_IntMeshObj_to_HDF5( rank=2,&
                                    dimensions=dimensions(1:2),&
                                    variable_name='/mesh/global/faces/element-sides', &
                                    int_variable=global_mesh % faces % elementSides,&
                                    file_id=file_id )

   INFO('End')

 END SUBROUTINE Write_MeshFaces

 SUBROUTINE Write_MeshNodes( global_mesh, file_id )
#undef __FUNC__
#define __FUNC__ "Write_MeshNodes"
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: global_mesh 
   INTEGER(HID_T), INTENT(in)      :: file_id
   !
   INTEGER(HSIZE_T) :: dimensions(1:2)
   INTEGER(HID_T)   :: group_id
   INTEGER          :: error


   INFO('Start')
       CALL h5gcreate_f( file_id, "/mesh/global/nodes", group_id, error )
       CALL h5gclose_f( group_id, error )

       dimensions(1:2) = (/3, global_mesh % nodes % nNodes /)  
       CALL Add_FloatMeshObj_to_HDF5( rank=2,&
                                      dimensions=dimensions,&
                                      variable_name='/mesh/global/nodes/positions', &
                                      float_variable=global_mesh % nodes % x,&
                                      file_id=file_id )
   INFO('End')

 END SUBROUTINE Write_MeshNodes

 SUBROUTINE Write_MeshGeometry( global_mesh, file_id )
#undef __FUNC__
#define __FUNC__ "Write_MeshGeometry"
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: global_mesh 
   INTEGER(HID_T), INTENT(in)      :: file_id
   !
   INTEGER(HSIZE_T) :: dimensions(1:5)
   INTEGER(HID_T)   :: group_id
   INTEGER          :: error


   INFO('Start')

       CALL h5gcreate_f( file_id, "/mesh/global/geometry", group_id, error )
       CALL h5gclose_f( group_id, error )

       dimensions(1:5) = (/global_mesh % elements % N+1, global_mesh % elements % N+1, global_mesh % elements % N+1, 3, global_mesh % elements % nElements /)  
       CALL Add_FloatMeshObj_to_HDF5( rank=5,&
                                      dimensions=dimensions,&
                                      variable_name='/mesh/global/geometry/positions', &
                                      float_variable=global_mesh % elements % x,&
                                      file_id=file_id )

       dimensions=(/global_mesh % elements % N+1, global_mesh % elements % N+1, 3, 6, global_mesh % elements % nElements /)
       CALL Add_FloatMeshObj_to_HDF5( rank=5,&
                                      dimensions=dimensions,&
                                      variable_name='/mesh/global/geometry/boundary-positions', &
                                      float_variable=global_mesh % elements % xBound,&
                                      file_id=file_id )

   INFO('End')

 END SUBROUTINE Write_MeshGeometry

 SUBROUTINE Write_MeshDecomp( global_mesh, file_id )
#undef __FUNC__
#define __FUNC__ "Write_MeshDecomp"
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: global_mesh 
   INTEGER(HID_T), INTENT(in)      :: file_id
   !
   INTEGER(HSIZE_T) :: dimensions(1)
   INTEGER(HID_T)   :: group_id
   INTEGER          :: error, blockid
   CHARACTER(6)     :: blockchar


   INFO('Start')

       CALL h5gcreate_f( file_id, "/mesh/decomp", group_id, error )
       CALL h5gclose_f( group_id, error )

       dimensions(1:1) = (/global_mesh % elements % nElements /)  
       CALL Add_IntMeshObj_to_HDF5( rank=1,&
                                    dimensions=dimensions(1:1),&
                                    variable_name='/mesh/decomp/element-to-blockid', &
                                    int_variable=global_mesh % decomp % element_to_blockID,&
                                    file_id=file_id )

       CALL Add_IntMeshObj_to_HDF5( rank=1,&
                                    dimensions=dimensions(1:1),&
                                    variable_name='/mesh/decomp/global-to-localid', &
                                    int_variable=global_mesh % decomp % global_to_localID,&
                                    file_id=file_id )

       DO blockid = 0, global_mesh % decomp % nBlocks-1
         WRITE(blockchar,'(I6.6)')blockid

         CALL h5gcreate_f( file_id, '/mesh/decomp/proc'//TRIM(blockchar), group_id, error )
         CALL h5gclose_f( group_id, error )

         dimensions(1:1) = (/global_mesh % decomp % mesh_obj(blockid) % nElements /)  
         CALL Add_IntMeshObj_to_HDF5( rank=1,&
                                      dimensions=dimensions(1:1),&
                                      variable_name='/mesh/decomp/proc'//TRIM(blockchar)//'/element-ids', &
                                      int_variable=global_mesh % decomp % mesh_obj(blockid) % elementids,&
                                      file_id=file_id )

         dimensions(1:1) = (/global_mesh % decomp % mesh_obj(blockid) % nNodes /)  
         CALL Add_IntMeshObj_to_HDF5( rank=1,&
                                      dimensions=dimensions(1:1),&
                                      variable_name='/mesh/decomp/proc'//TRIM(blockchar)//'/node-ids', &
                                      int_variable=global_mesh % decomp % mesh_obj(blockid) % nodeids,&
                                      file_id=file_id )
        
         dimensions(1:1) = (/global_mesh % decomp % mesh_obj(blockid) % nFaces /)  
         CALL Add_IntMeshObj_to_HDF5( rank=1,&
                                      dimensions=dimensions(1:1),&
                                      variable_name='/mesh/decomp/proc'//TRIM(blockchar)//'/face-ids', &
                                      int_variable=global_mesh % decomp % mesh_obj(blockid) % faceids,&
                                      file_id=file_id )
        
        
       ENDDO

   INFO('End')

 END SUBROUTINE Write_MeshDecomp

 SUBROUTINE Read_MeshDecomp( local_mesh, file_id, my_RankID )
#undef __FUNC__
#define __FUNC__ "Read_MeshDecomp"
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout) :: local_mesh 
   INTEGER(HID_T), INTENT(in)      :: file_id
   INTEGER, INTENT(in)             :: my_RankID
   !
   INTEGER(HSIZE_T) :: dimensions(1:1)
   INTEGER(HSIZE_T) :: maxdims(1:1)
   INTEGER(HID_T)   :: dataset_id, filespace
   INTEGER          :: error, blockid, ndims
   INTEGER          :: nGlobalElements, nLocalNodes(0:0), nLocalFaces(0:0), nLocalElements(0:0)
   CHARACTER(6)     :: blockchar
   CHARACTER(50)    :: msg


   INFO('Start')


       CALL Get_HDF5_Obj_Dimensions( file_id, '/mesh/decomp/element-to-blockid', 1, dimensions )
       nGlobalElements = dimensions(1)

       WRITE(blockchar,'(I6.6)')my_RankID
       CALL Get_HDF5_Obj_Dimensions( file_id, '/mesh/decomp/proc'//TRIM(blockchar)//'/element-ids', 1, dimensions )
       nLocalElements(0) = dimensions(1)

       CALL Get_HDF5_Obj_Dimensions( file_id, '/mesh/decomp/proc'//TRIM(blockchar)//'/node-ids', 1, dimensions )
       nLocalNodes(0) = dimensions(1)

       CALL Get_HDF5_Obj_Dimensions( file_id, '/mesh/decomp/proc'//TRIM(blockchar)//'/face-ids', 1, dimensions )
       nLocalFaces(0) = dimensions(1)

       CALL local_mesh % decomp % Build( nGlobalElements, 1, nLocalElements, nLocalFaces, nLocalNodes )

       WRITE(msg,'(I5)')nGlobalElements
       msg = 'Number of Global Elements : '//TRIM(msg)
       INFO(msg)

       WRITE(msg,'(I5)')nLocalElements(0)
       msg = 'Number of Local Elements : '//TRIM(msg)
       INFO(msg)

       WRITE(msg,'(I5)')nLocalFaces(0)
       msg = 'Number of Faces : '//TRIM(msg)
       INFO(msg)

       WRITE(msg,'(I5)')nLocalNodes(0)
       msg = 'Number of Nodes : '//TRIM(msg)
       INFO(msg)
       
       dimensions(1:1) = nGlobalElements 
       CALL h5dopen_f(file_id, '/mesh/decomp/element-to-blockid', dataset_id, error)
       CALL h5dget_space_f( dataset_id, filespace, error )
       CALL h5dread_f( dataset_id, H5T_STD_I32LE, &
                       local_mesh % decomp % element_to_blockID, &
                       dimensions, error)

       CALL h5dopen_f(file_id, '/mesh/decomp/global-to-localid', dataset_id, error)
       CALL h5dget_space_f( dataset_id, filespace, error )
       CALL h5dread_f( dataset_id, H5T_STD_I32LE, &
                       local_mesh % decomp % global_to_localID, &
                       dimensions, error)

       dimensions(1:1) = nLocalElements(0:0)  
       CALL h5dopen_f(file_id, '/mesh/decomp/proc'//TRIM(blockchar)//'/element-ids', dataset_id, error)
       CALL h5dget_space_f( dataset_id, filespace, error )
       CALL h5dread_f( dataset_id, H5T_STD_I32LE, &
                       local_mesh % decomp % mesh_obj(0) % elementids, &
                       dimensions, error)

       dimensions(1:1) = nLocalNodes(0:0)
       CALL h5dopen_f(file_id, '/mesh/decomp/proc'//TRIM(blockchar)//'/node-ids', dataset_id, error)
       CALL h5dget_space_f( dataset_id, filespace, error )
       CALL h5dread_f( dataset_id, H5T_STD_I32LE, &
                       local_mesh % decomp % mesh_obj(0) % nodeids, &
                       dimensions, error)
       
       dimensions(1:1) = nLocalFaces(0:0)
       CALL h5dopen_f(file_id, '/mesh/decomp/proc'//TRIM(blockchar)//'/face-ids', dataset_id, error)
       CALL h5dget_space_f( dataset_id, filespace, error )
       CALL h5dread_f( dataset_id, H5T_STD_I32LE, &
                       local_mesh % decomp % mesh_obj(0) % faceids, &
                       dimensions, error)
       
        
   INFO('End')

 END SUBROUTINE Read_MeshDecomp

 SUBROUTINE Get_HDF5_Obj_Dimensions( file_id, variable_name, rank, dimensions )
#undef __FUNC__
#define __FUNC__ "Get_HDF5_Obj_Dimensions"
   IMPLICIT NONE
   INTEGER(HID_T), INTENT(in)    :: file_id
   CHARACTER(*)                  :: variable_name
   INTEGER, INTENT(in)           :: rank
   INTEGER(HSIZE_T), INTENT(out) :: dimensions(1:rank)
   ! Local
   INTEGER          :: error
   INTEGER(HID_T)   :: dataset_id
   INTEGER(HSIZE_T) :: maxdims(1:rank)
   INTEGER(HID_T)   :: filespace

        CALL h5dopen_f(file_id, TRIM(variable_name), dataset_id, error)
        CALL h5dget_space_f(dataset_id, filespace, error)
        CALL h5sget_simple_extent_dims_f(filespace, dimensions, maxdims, error) 
        
 END SUBROUTINE Get_HDF5_Obj_Dimensions

 SUBROUTINE Add_FloatMeshObj_to_HDF5( rank, dimensions, variable_name, float_variable, file_id )
#undef __FUNC__
#define __FUNC__ "Add_FloatMeshObj_to_HDF5"
   IMPLICIT NONE
   INTEGER, INTENT(in)          :: rank
   INTEGER(HSIZE_T), INTENT(in) :: dimensions(1:rank)
   CHARACTER(*)                 :: variable_name
   REAL(prec), INTENT(in)       :: float_variable(*)
   INTEGER(HID_T), INTENT(in)   :: file_id
   ! local
   INTEGER(HID_T)   :: memspace, dataset_id
   INTEGER          :: error

       CALL h5screate_simple_f(rank, dimensions, memspace, error)
       CALL h5dcreate_f( file_id, TRIM(variable_name), H5T_IEEE_F32LE, memspace, dataset_id, error)
       CALL h5dwrite_f( dataset_id, H5T_IEEE_F32LE, float_variable, dimensions, error)
       CALL h5dclose_f( dataset_id, error )
       CALL h5sclose_f( memspace, error )

 END SUBROUTINE Add_FloatMeshObj_to_HDF5

 SUBROUTINE Add_IntMeshObj_to_HDF5( rank, dimensions, variable_name, int_variable, file_id )
#undef __FUNC__
#define __FUNC__ "Add_IntMeshObj_to_HDF5"
   IMPLICIT NONE
   INTEGER, INTENT(in)          :: rank
   INTEGER(HSIZE_T), INTENT(in) :: dimensions(1:rank)
   CHARACTER(*)                 :: variable_name
   INTEGER, INTENT(in)          :: int_variable(*)
   INTEGER(HID_T), INTENT(in)   :: file_id
   ! local
   INTEGER(HID_T)   :: memspace, dataset_id
   INTEGER          :: error

       CALL h5screate_simple_f(rank, dimensions, memspace, error)
       CALL h5dcreate_f( file_id, TRIM(variable_name), H5T_STD_I32LE, memspace, dataset_id, error)
       CALL h5dwrite_f( dataset_id, H5T_STD_I32LE , int_variable, dimensions, error)
       CALL h5dclose_f( dataset_id, error )
       CALL h5sclose_f( memspace, error )

 END SUBROUTINE Add_IntMeshObj_to_HDF5
  
 SUBROUTINE Load_SELFMesh( mesh, params, my_RankID, nMPI_Ranks, mpiCommunicator )
#undef __FUNC__
#define __FUNC__ "Load_SELFMesh"
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout)     :: mesh 
   TYPE( ModelParameters ), INTENT(in) :: params
   INTEGER, INTENT(in)                 :: my_RankID
   INTEGER, INTENT(in)                 :: nMPI_Ranks
   INTEGER, INTENT(in)                 :: mpiCommunicator

   ! Read From HDF5
   !  - If serial, load in the global mesh. Then create the boundary map
   !  - If parallel, load in the local process mesh information, the element_to_blockID array and the decomp objects into decomp % mesh_obj(0). 
   !       Then adjust the face information and create the boundary map 

      CALL mesh % Read_SELFMesh( params % SELFMeshFile, params % polydeg, my_RankID, nMPI_Ranks, mpiCommunicator )
      CALL mesh % SetupFaceMaps( )
      CALL SetupProcessBoundaryMap( mesh, my_RankID )


 END SUBROUTINE Load_SELFMesh

 SUBROUTINE Generate_SELFMesh(paramFile, equationFile, nMPI_Ranks)
#undef __FUNC__
#define __FUNC__ "Generate_SELFMesh"
   ! Generates a global SELFMeshFile with a boundary map and
   ! decomposition for a structured mesh or a UCD mesh file.
   IMPLICIT NONE
   CHARACTER(*), INTENT(in) :: paramFile
   CHARACTER(*), INTENT(in) :: equationFile
   INTEGER, INTENT(in)      :: nMPI_Ranks
   ! Local
   LOGICAL                     :: setupSuccess
   TYPE( NodalDG )             :: nodal
   TYPE( HexMesh )             :: mesh
   TYPE( ModelParameters )     :: params
   TYPE( Geom_EquationParser ) :: geomParser
   INTEGER :: nbf, i

   INFO('Start')
      
      CALL params % Build( TRIM(paramFile), setupSuccess )
      CALL geomParser % Build( equationFile )

      ! Build an interpolant
      CALL nodal % Build( targetPoints = UniformPoints(-1.0_prec,1.0_prec,params % nPlot), &
                          N = params % polyDeg, &
                          nTargetPoints = params % nPlot, &
                          quadrature = GAUSS )
   
      ! Build the Geometry
      IF( TRIM( params % UCDMeshFile ) == '' )THEN
        CALL mesh % ConstructStructuredMesh( nodal % interp, params % nXelem, params % nYelem, params % nZelem,geomParser )
        ! TO DO : push the mesh scaling inside the structured mesh construction
        CALL mesh % ScaleTheMesh( nodal % interp, params % xScale, params % yScale, params % zScale  )
      ELSE   
        CALL mesh % Read_TrellisUCDMesh( nodal % interp, TRIM( params % UCDMeshFile ) )           
      ENDIF


      CALL mesh % boundaryMap % Build(mesh % faces % nBoundaryFaces(),0)

      ! Set the boundary map for the global mesh
      nbf = 0
      DO i = 1, mesh % faces % nFaces
        IF( mesh % faces % elementIDs(2,i) < 0 )THEN
          nbf = nbf+1
          mesh % faces % boundaryID(i)                = nbf  
          mesh % boundaryMap % boundaryIDs(nbf)       = i
          mesh % boundaryMap % boundaryGlobalIDs(nbf) = i
          mesh % boundaryMap % extProcIDs(nbf)        = 0
          mesh % boundaryMap % boundarycondition(nbf) = mesh % faces % elementIDs(2,i)
        ENDIF
      ENDDO

      CALL MeshDecompose(mesh, params, nMPI_Ranks )

      CALL mesh % Write_TecplotMesh(TRIM(params % SELFMeshFile))
      CALL mesh % Write_SELFMesh(TRIM(params % SELFMeshFile), nMPI_Ranks)
      CALL nodal % Trash( )
      CALL mesh % Trash( )

    INFO('End')

 END SUBROUTINE Generate_SELFMesh
!
 SUBROUTINE MeshDecompose( mesh, params, nMPI_Ranks )
#undef __FUNC__
#define __FUNC__ "MeshDecompose"
   IMPLICIT NONE
   TYPE( HexMesh ), INTENT(inout)        :: mesh
   TYPE( ModelParameters ),INTENT(inout) :: params  
   INTEGER, INTENT(in)                   :: nMPI_Ranks
   ! Local 
   LOGICAL                 :: setupSuccess

   INFO('Start')

     ALLOCATE( mesh % decomp % element_to_blockID(1:mesh % elements % nElements) )
     ALLOCATE( mesh % decomp % global_to_localID(1:mesh % elements % nElements) )
     mesh % decomp % element_to_blockID = 0
     mesh % decomp % global_to_localID = 0
     mesh % decomp % nGlobalElements = mesh % elements % nElements
     mesh % decomp % nBlocks = nMPI_Ranks

     IF( TRIM( params % UCDMeshFile ) /= '' )THEN
       IF( nMPI_Ranks > 1 )THEN
         CALL UnstructuredDecompose( mesh, nMPI_Ranks )
       ENDIF
     ELSE
       CALL StructuredDecompose( params, nMPI_Ranks )
       CALL StructuredMeshToBlocks( mesh, params )
     ENDIF
     CALL SetGlobalToLocalMapping( mesh, nMPI_Ranks )
     CALL CreateMeshObjLists( mesh, nMPI_Ranks )

   INFO('End')

 END SUBROUTINE MeshDecompose

 SUBROUTINE UnstructuredDecompose( global_mesh, nMPI_Ranks )
#undef __FUNC__
#define __FUNC__ "UnstructuredDecompose"
   IMPLICIT NONE
   CLASS( HexMesh ), INTENT(inout)        :: global_mesh
   INTEGER, INTENT(in)                    :: nMPI_Ranks
   ! Local
   INTEGER :: i, iel, n
   CHARACTER(60)  :: msg
   INTEGER :: n_partitions, n_nodes, n_elements
   INTEGER :: eptr(1:global_mesh % elements % nElements+1)
   INTEGER :: eind(1:8*global_mesh % elements % nElements)
   INTEGER :: epart(1:global_mesh % elements % nElements)
   INTEGER :: npart(1:global_mesh % nodes % nNodes)
   INTEGER, POINTER :: vwgt =>null(), vsize => null(), options => null()
   REAL(8), POINTER :: tpwgts=>null()
   EXTERNAL :: METIS_PartMeshNodal

   INFO('Start')

   ! Count the number of graph edges (number of internal faces)
   n_nodes      = global_mesh % nodes % nNodes
   n_elements   = global_mesh % elements % nElements
   n_partitions = nMPI_Ranks

   INFO('Passing off to METIS_PartMeshNodal')

   DO iel = 1, global_mesh % elements % nElements
     eptr(iel)   = 8*(iel-1)
     eptr(iel+1) = 8*(iel)
     DO i = 1, 8
       eind(i+8*(iel-1)) = global_mesh % elements % nodeIDs(i,iel)-1
     ENDDO
   ENDDO

   CALL METIS_PartMeshNodal(n_elements, n_nodes, eptr, eind, vwgt, vsize, n_partitions, tpwgts, options, n, epart, npart)

   INFO('METIS_PartMeshNodal complete')

   DO iel = 1, global_mesh % elements % nElements
     global_mesh % decomp % element_to_blockID(iel) = epart(iel)
   ENDDO

   INFO('End')

 END SUBROUTINE UnstructuredDecompose

 SUBROUTINE SetupProcessBoundaryMap( mesh, my_RankID )
#undef __FUNC__
#define __FUNC__ "SetupProcessBoundaryMap"
 ! This routine takes in a mesh block assigned to a given MPI rank, along with
 ! the filled in decomp object. When the decomp object is created on an MPI rank
 ! only the 0th array element is allocated
 ! This routine also assumes that the decomp element_to_blockID and global_to_localID arrays
 ! have been filled in
 !
 ! On input, the mesh has the node, element, and face information filled in. However, the
 ! face element IDs refer to global element IDs and must be converted to local element IDs
 !
   IMPLICIT NONE
   TYPE( HexMesh ), INTENT(inout) :: mesh
   INTEGER,  INTENT(in)           :: my_RankID
   ! Local
   INTEGER :: iFace, e1, e2, p1, p2, e1Local, e2Local, s1, s2
   INTEGER :: nbf, nmpi
   CHARACTER(30) :: msg
   
   INFO('Start')

     ! Count the number of boundary faces first
     nbf  = 0
     nmpi = 0
     DO iFace = 1, mesh % faces % nFaces

       e1 = mesh % faces % elementIDs(1,iFace)
       e2 = mesh % faces % elementIDs(2,iFace)
       p1 = mesh % decomp % element_to_blockID(e1)

       IF( e2 > 0 )THEN ! Global internal face
         p2 = mesh % decomp % element_to_blockID(e2)

         IF( p1 == my_RankID .AND. p2 == my_RankID )THEN ! Local internal face
          ! mesh % faces % elementIDs(1,iFace) = mesh % decomp % global_to_localID(e1)
          ! mesh % faces % elementIDs(2,iFace) = mesh % decomp % global_to_localID(e2)
         ELSE ! MPI Boundary
            nbf = nbf + 1
            nmpi = nmpi + 1
         ENDIF

       ELSE ! Physical Boundary

         nbf = nbf + 1

       ENDIF

     ENDDO

     WRITE(msg,'(A,I5)') 'Number of boundary faces : ',nbf
     INFO( TRIM(msg) )
     WRITE(msg,'(A,I5)') 'Number of MPI exchanges : ',nmpi
     INFO( TRIM(msg) )

     CALL mesh % boundaryMap % build( nbf, nmpi )

     nbf  = 0
     nmpi = 0
     DO iFace = 1, mesh % faces % nFaces

       e1 = mesh % faces % elementIDs(1,iFace)
       e2 = mesh % faces % elementIDs(2,iFace)
       p1 = mesh % decomp % element_to_blockID(e1)

       IF( e2 > 0 )THEN ! Global internal face
         p2 = mesh % decomp % element_to_blockID(e2)

         IF( p1 == my_RankID .AND. p2 == my_RankID )THEN ! Local internal face
           mesh % faces % elementIDs(1,iFace) = mesh % decomp % global_to_localID(e1)
           mesh % faces % elementIDs(2,iFace) = mesh % decomp % global_to_localID(e2)
         ELSE ! MPI Boundary
            nbf = nbf + 1
            nmpi = nmpi + 1

            e1local = mesh % decomp % global_to_localID(e1)
            e2local = mesh % decomp % global_to_localID(e2)
            s1 = mesh % faces % elementSides(1,iFace)
            s2 = mesh % faces % elementSides(2,iFace)

            IF( p1 == my_RankID )THEN

              mesh % faces % elementIDs(1,iFace)   = e1local
              mesh % faces % elementIDs(2,iFace)   = -e2 ! Store the negative of the global element ID as the secondary element
              mesh % faces % elementSides(1,iFace) = s1  
              mesh % faces % elementSides(2,iFace) = s2  
              mesh % faces % boundaryID(iFace)     = nbf  

              mesh % boundaryMap % extProcIDs(nbf)        = p2
              mesh % boundaryMap % boundaryIDs(nbf)       = iFace
              mesh % boundaryMap % boundaryGlobalIDs(nbf) = mesh % faces % faceID(iFace) 
              mesh % boundaryMap % boundaryCondition(nbf) = MPI_BOUNDARY 
              mesh % boundaryMap % boundary_to_mpi(nmpi)  = nbf

            ELSE ! this process owns the secondary element - need to swap the ids and sides

              mesh % faces % elementIDs(1,iFace)   = e2local
              mesh % faces % elementIDs(2,iFace)   = -e1 ! Store the negative of the global element ID as the secondary element
              mesh % faces % elementSides(1,iFace) = s2  
              mesh % faces % elementSides(2,iFace) = s1  
              mesh % faces % boundaryID(iFace)     = nbf  

              mesh % boundaryMap % extProcIDs(nbf)        = p1
              mesh % boundaryMap % boundaryIDs(nbf)       = iFace
              mesh % boundaryMap % boundaryGlobalIDs(nbf) = mesh % faces % faceID(iFace) 
              mesh % boundaryMap % boundaryCondition(nbf) = MPI_BOUNDARY 
              mesh % boundaryMap % boundary_to_mpi(nmpi)  = nbf

            ENDIF
         ENDIF

       ELSE ! Physical Boundary

         nbf = nbf + 1
         mesh % faces % elementIDs(1,iFace)   = e1local
         mesh % faces % elementSides(1,iFace) = s1  
         mesh % faces % boundaryID(iFace)     = nbf  

         mesh % boundaryMap % extProcIDs(nbf)        = p1
         mesh % boundaryMap % boundaryIDs(nbf)       = iFace
         mesh % boundaryMap % boundaryGlobalIDs(nbf) = mesh % faces % faceID(iFace) 
         mesh % boundaryMap % boundaryCondition(nbf) = e2 ! Retain the boundary condition from the face information 

       ENDIF

     ENDDO

   INFO('End')

 END SUBROUTINE SetupProcessBoundaryMap

 SUBROUTINE CreateMeshObjLists( global_mesh, nMPI_Ranks )
#undef __FUNC__
#define __FUNC__ "CreateMeshObjLists"
   IMPLICIT NONE
   TYPE( HexMesh ), INTENT(inout) :: global_mesh
   INTEGER, INTENT(in)            :: nMPI_Ranks
   ! Local
   INTEGER :: iEl, iNode, iFace, e1, e2
   INTEGER :: procID, nodeID, faceID
   INTEGER :: nNodes, nFaces
   INTEGER :: nEl(0:nMPI_Ranks-1)
   INTEGER :: nodeLogic(1:global_mesh % nodes % nNodes,0:nMPI_Ranks-1)
   INTEGER :: faceLogic(1:global_mesh % faces % nFaces,0:nMPI_Ranks-1)
   CHARACTER(30) :: msg

   INFO('Start')

     nEl   = 0
     nodeLogic = 0
     faceLogic = 0
     ! Set the logic mask for the nodes and count the number of elements
     DO iEl = 1, global_mesh % elements % nElements
       procID = global_mesh % decomp % element_to_blockID(iEl)
       nEl(procID) = nEl(procID) + 1
       DO iNode = 1, 8
         nodeID = global_mesh % elements % nodeIDs(iNode,iEl) 
         nodeLogic(nodeID,procID) = 1
       ENDDO
     ENDDO
    
     ! Set the logic mask for the faces
     DO procID = 0, nMPI_Ranks -1
       DO iFace = 1, global_mesh % faces % nFaces
         e1 = global_mesh % faces % elementIDs(1,iFace)
         e2 = global_mesh % faces % elementIDs(2,iFace)
         IF( global_mesh % decomp % element_to_blockID(e1) == procID .OR. &
             global_mesh % decomp % element_to_blockID(e2) == procID )THEN

           faceLogic(iFace,procID) = 1
 
         ENDIF
       ENDDO
     ENDDO

     ! Allocate space for the decomposition mesh object lists
     ! TO DO : Constructor for decomp
     ALLOCATE( global_mesh % decomp % mesh_obj(0:nMPI_Ranks-1) )
     global_mesh % decomp % nBlocks = nMPI_Ranks
     global_mesh % decomp % nGlobalElements = global_mesh % elements % nElements 
     DO procID = 0, nMPI_Ranks-1
       ALLOCATE( global_mesh % decomp % mesh_obj(procID) % elementids(1:nEl(procID)) ) 
       nNodes = SUM( nodeLogic(:,procID) )
       ALLOCATE( global_mesh % decomp % mesh_obj(procID) % nodeids(1:nNodes) ) 
       nFaces = SUM( faceLogic(:,procID) )
       ALLOCATE( global_mesh % decomp % mesh_obj(procID) % faceids(1:nFaces) ) 
       global_mesh % decomp % mesh_obj(procID) % nElements = nEl(procID)
       global_mesh % decomp % mesh_obj(procID) % nNodes    = nNodes
       global_mesh % decomp % mesh_obj(procID) % nFaces    = nFaces

       WRITE(msg,'(A,I5)') 'Process ID : ',procID
       INFO( TRIM(msg) )
       WRITE(msg,'(A,I5)') 'Number of elements : ',nEl(procID)
       INFO( TRIM(msg) )
       WRITE(msg,'(A,I5)') 'Number of nodes : ',nNodes
       INFO( TRIM(msg) )
       WRITE(msg,'(A,I5)') 'Number of faces : ',nFaces
       INFO( TRIM(msg) )

     ENDDO


     ! Set the element lists
     nEl   = 0
     DO iEl = 1, global_mesh % elements % nElements
       procID = global_mesh % decomp % element_to_blockID(iEl)
       nEl(procID) = nEl(procID) + 1
       global_mesh % decomp % mesh_obj(procID) % elementids(nEl(procID)) = iEl
     ENDDO

     ! Set the node lists
     DO procID = 0, nMPI_Ranks -1
       nNodes = 0
       DO iNode = 1, global_mesh % nodes % nNodes
         IF( nodeLogic(iNode,procID) == 1 )THEN
           nNodes = nNodes + 1
           global_mesh % decomp % mesh_obj(procID) % nodeids(nNodes) = iNode
         ENDIF
       ENDDO
     ENDDO

     ! Set the face lists
     DO procID = 0, nMPI_Ranks -1
       nFaces = 0
       DO iFace = 1, global_mesh % faces % nFaces
         IF( faceLogic(iFace,procID) == 1 )THEN
           nFaces = nFaces + 1
           global_mesh % decomp % mesh_obj(procID) % faceids(nFaces) = iFace
         ENDIF
       ENDDO
     ENDDO

   INFO('End')

 END SUBROUTINE CreateMeshObjLists

 SUBROUTINE SetGlobalToLocalMapping( global_mesh, nMPI_Ranks )
#undef __FUNC__
#define __FUNC__ "SetGlobalToLocalMapping"
   IMPLICIT NONE
   TYPE( HexMesh ), INTENT(inout) :: global_mesh
   INTEGER, INTENT(in)            :: nMPI_Ranks
   ! Local
   INTEGER :: nEl(0:nMPI_Ranks), procID, iEl

   INFO('Start')

     nEl = 0
     DO iEl = 1, global_mesh % elements % nElements
       procID = global_mesh % decomp % element_to_blockID(iEl)
       nEl(procID) = nEl(procID) + 1
       global_mesh % decomp % global_to_localID(iEl) = nEl(procID)
     ENDDO

   INFO('End')

 END SUBROUTINE SetGlobalToLocalMapping

 SUBROUTINE StructuredMeshToBlocks( global_mesh, params )
#undef __FUNC__
#define __FUNC__ "StructuredMeshToBlocks"
   IMPLICIT NONE
   TYPE( HexMesh ), INTENT(inout)        :: global_mesh
   TYPE( ModelParameters ),INTENT(inout) :: params  
   ! Local
   INTEGER :: iPz, iPy, iPx
   INTEGER :: iZp, iYp, iXp
   INTEGER :: iZ, iY, iX
   INTEGER :: iNode, nID
   INTEGER :: iEl, procID
   INTEGER :: nxp(1:params % nProcX)
   INTEGER :: nyp(1:params % nProcY)
   INTEGER :: nzp(1:params % nProcZ)
   INTEGER :: nEl(0:params % nProc-1)
    
   INFO('Start')


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
                    IF( iEl <= global_mesh % elements % nElements )THEN
                      global_mesh % decomp % element_to_blockID(iEl) = iPx-1 + params % nProcX*( iPy-1 + params % nProcY*(iPz-1) )
                    ENDIF

                  ENDDO
                ENDDO

              ENDDO
            ENDDO

          ENDDO
        ENDDO


   INFO('End')

 END SUBROUTINE StructuredMeshToBlocks

 SUBROUTINE StructuredDecompose( params, nMPI_Ranks )
#undef __FUNC__
#define __FUNC__ "StructuredDecompose"
   IMPLICIT NONE
   TYPE( ModelParameters ), INTENT(inout) :: params
   INTEGER, INTENT(in)                    :: nMPI_Ranks
   ! Local
   INTEGER :: procDim(1:3), meshDim(1:3), dimIndex(1:3), temp
   INTEGER :: i, j 
   CHARACTER(50) :: msg

   INFO('Start')
      params % nProc = nMPI_Ranks

      IF( nMPI_Ranks > 1 )THEN

         IF( params % nProc > params % nXElem*params % nYElem*params % nZElem )THEN

           INFO('Number of processes exceeds number of elements.')
           INFO('Cannot decompose.')
           INFO('Stopping.')
           STOP

         ELSEIF( params % nProc == params % nXElem*params % nYElem*params % nZElem )THEN

           params % nProcX = params % nXElem
           params % nProcY = params % nYElem
           params % nProcZ = params % nZElem

         ELSE

           meshDim(1:3)  = (/ params % nXElem, params % nYElem, params % nZElem /) 
           dimIndex(1:3) = (/ 1, 2, 3 /) 

           ! Sort from largest mesh dimension to smallest with insertion sort
           DO i = 2,  3
             j = i
             DO WHILE( j > 1 )
               IF( meshDim(j-1) < meshDim(j) )THEN

                 temp         = meshDim(j)
                 meshDim(j)   = meshDim(j-1)
                 meshDim(j-1) = temp 

                 temp          = dimIndex(j)
                 dimIndex(j)   = dimIndex(j-1)
                 dimIndex(j-1) = temp 
                 j = j-1

               ELSE
                 EXIT
               ENDIF
             ENDDO
           ENDDO
           
           procDim(1:3)  = (/ params % nProc, 1, 1 /)

           DO WHILE( procDim(1) > meshDim(1) )

             DO i = 2, procDim(1)
               IF( MOD( procDim(1), i ) == 0 )THEN
                 j = i
                 EXIT
               ENDIF
             ENDDO

             procDim(1) = procDim(1)/j

           ENDDO

           IF( procDim(1)*procDim(2)*procDim(3) < params % nProc )THEN

             procDim(2) = params % nProc/procDim(1)

             IF( procDim(2) > 1 )THEN
               DO WHILE( procDim(2) > meshDim(2) )
  
                 DO i = 2, procDim(2)
                   IF( MOD( procDim(2), i ) == 0 )THEN
                     j = i
                     EXIT
                   ENDIF
                 ENDDO
  
                 procDim(2) = procDim(2)/j
  
               ENDDO
  
             ENDIF

           ENDIF

           IF( procDim(1)*procDim(2)*procDim(3) < params % nProc )THEN

             procDim(3) = params % nProc/(procDim(1)*procDim(2))
             IF( procDim(3) > 1 )THEN
               DO WHILE( procDim(3) > meshDim(3) )
  
                 DO i = 2, procDim(3)
                   IF( MOD( procDim(3), i ) == 0 )THEN
                     j = i
                     EXIT
                   ENDIF
                 ENDDO
  
                 procDim(3) = procDim(3)/j
  
               ENDDO
  
             ENDIF

           ENDIF

           DO i = 1, 3

             IF( dimIndex(i) == 1 )THEN
               params % nProcX = procDim(i)
             ELSEIF( dimIndex(i) == 2 )THEN
               params % nProcY = procDim(i)
             ELSE
               params % nProcZ = procDim(i)
             ENDIF

           ENDDO

         ENDIF

       ELSE

         params % nProcX = 1
         params % nProcY = 1
         params % nProcZ = 1
 
       ENDIF

       WRITE(msg,'(A,I5)')'nProcX :', params % nProcX
       INFO( TRIM(msg) )
       WRITE(msg,'(A,I5)')'nProcY :', params % nProcY
       INFO( TRIM(msg) )
       WRITE(msg,'(A,I5)')'nProcZ :', params % nProcZ
       INFO( TRIM(msg) )

   INFO('End')
      
 END SUBROUTINE StructuredDecompose     

! SUBROUTINE StructuredMeshGenerator_3D( paramFile, equationFile, setupSuccess, nMPI_Ranks )
!#undef __FUNC__
!#define __FUNC__ "StructuredMeshGenerator_3D"
! LOGICAL, INTENT(out)                      :: setupSuccess
! CHARACTER(*), INTENT(in)                  :: paramFile
! CHARACTER(*), INTENT(in)                  :: equationFile
! INTEGER, INTENT(in)                       :: nMPI_Ranks
!  ! Local
! TYPE( NodalDG )                           :: nodal
! TYPE( HexMesh )                           :: mesh
! TYPE( HexMesh )                           :: procMesh
! TYPE( ModelParameters )                   :: params
! TYPE( Geom_EquationParser )               :: geomParser
! INTEGER, ALLOCATABLE                      :: faceProcCount(:), faceProcTable(:,:), faceProcOwners(:,:), faceBoundaryIDs(:,:)
! INTEGER                      :: procID, mpiErr
! INTEGER                      :: nElems, nProc, pID, i, j
! INTEGER                      :: iEl, jEl, iSide, iNode, nAdj, nID, elID, nMPI
! INTEGER                      :: e1, e2, p1, p2, iFace, iFaceLocal, jFace, localID
! INTEGER                      :: nBe, npFaces
! INTEGER                      :: globalFaceID, localFaceID, bID, extProc, extBID
! CHARACTER(4)                 :: pIDChar
! INTEGER, ALLOCATABLE         :: globalToLocal(:,:), nElPerProc(:), nodeLogic(:,:), nNodePerProc(:), partitions(:)
! INTEGER, ALLOCATABLE         :: globalToLocalNode(:,:), nLocMPI(:)
! REAL(prec), ALLOCATABLE      :: materials(:)
! INTEGER :: procDim(1:3), meshDim(1:3), dimIndex(1:3), temp
! CHARACTER(50) :: msg
! 
!
!      INFO('Start')
!      ! Read in the PARAMETERs
!      CALL params % Build( TRIM(paramFile), setupSuccess )
!
!#ifdef HAVE_MPI
!      params % nProc = nMPI_Ranks
!
!      IF( params % nProc > params % nXElem*params % nYElem*params % nZElem )THEN
!
!        INFO('Number of processes exceeds number of elements.')
!        setupSuccess = .FALSE.
!
!      ELSEIF( params % nProc == params % nXElem*params % nYElem*params % nZElem )THEN
!
!        params % nProcX = params % nXElem
!        params % nProcY = params % nYElem
!        params % nProcZ = params % nZElem
!
!      ELSE
!
!
!        meshDim(1:3)  = (/ params % nXElem, params % nYElem, params % nZElem /) 
!        dimIndex(1:3) = (/ 1, 2, 3 /) 
!
!        ! Sort from largest mesh dimension to smallest with insertion sort
!        DO i = 2,  3
!          j = i
!          DO WHILE( j > 1 )
!            IF( meshDim(j-1) < meshDim(j) )THEN
!
!              temp         = meshDim(j)
!              meshDim(j)   = meshDim(j-1)
!              meshDim(j-1) = temp 
!
!              temp          = dimIndex(j)
!              dimIndex(j)   = dimIndex(j-1)
!              dimIndex(j-1) = temp 
!              j = j-1
!
!            ELSE
!              EXIT
!            ENDIF
!          ENDDO
!        ENDDO
!        
!        procDim(1:3)  = (/ params % nProc, 1, 1 /)
!
!        DO WHILE( procDim(1) > meshDim(1) )
!
!          DO i = 2, procDim(1)
!            IF( MOD( procDim(1), i ) == 0 )THEN
!              j = i
!              EXIT
!            ENDIF
!          ENDDO
!
!          procDim(1) = procDim(1)/j
!
!        ENDDO
!
!        IF( procDim(1)*procDim(2)*procDim(3) < params % nProc )THEN
!
!          procDim(2) = params % nProc/procDim(1)
!
!          IF( procDim(2) > 1 )THEN
!            DO WHILE( procDim(2) > meshDim(2) )
!  
!              DO i = 2, procDim(2)
!                IF( MOD( procDim(2), i ) == 0 )THEN
!                  j = i
!                  EXIT
!                ENDIF
!              ENDDO
!  
!              procDim(2) = procDim(2)/j
!  
!            ENDDO
!  
!          ENDIF
!
!        ENDIF
!
!        IF( procDim(1)*procDim(2)*procDim(3) < params % nProc )THEN
!
!          procDim(3) = params % nProc/(procDim(1)*procDim(2))
!          IF( procDim(3) > 1 )THEN
!            DO WHILE( procDim(3) > meshDim(3) )
!  
!              DO i = 2, procDim(3)
!                IF( MOD( procDim(3), i ) == 0 )THEN
!                  j = i
!                  EXIT
!                ENDIF
!              ENDDO
!  
!              procDim(3) = procDim(3)/j
!  
!            ENDDO
!  
!          ENDIF
!
!        ENDIF
!
!        DO i = 1, 3
!
!          IF( dimIndex(i) == 1 )THEN
!            params % nProcX = procDim(i)
!          ELSEIF( dimIndex(i) == 2 )THEN
!            params % nProcY = procDim(i)
!          ELSE
!            params % nProcZ = procDim(i)
!          ENDIF
!
!        ENDDO
!
!      ENDIF
!
!      WRITE(msg,'(A,I5)')'nProcX :', params % nProcX
!      INFO( TRIM(msg) )
!      WRITE(msg,'(A,I5)')'nProcY :', params % nProcY
!      INFO( TRIM(msg) )
!      WRITE(msg,'(A,I5)')'nProcZ :', params % nProcZ
!      INFO( TRIM(msg) )
!      
!#endif
!      CALL geomParser % Build( equationFile )
!
!      IF( setupSuccess )THEN
!         ! Build an interpolant
!         CALL nodal % Build( targetPoints = UniformPoints(-1.0_prec,1.0_prec,params % nPlot), &
!                             N = params % polyDeg, &
!                             nTargetPoints = params % nPlot, &
!                             quadrature = GAUSS_LOBATTO )
!   
!         ! Build the Geometry
!         IF( TRIM( params % UCDMeshFile ) == '' )THEN
!              CALL mesh % ConstructStructuredMesh( nodal % interp, &
!                params % nXelem, &
!                params % nYelem, &
!                params % nZelem, &
!                geomParser )
!         ELSE   
!
!           CALL mesh % ReadTrellisUCDMeshFile( nodal % interp, TRIM( params % UCDMeshFile ) )           
!#ifdef TECPLOT
!           CALL mesh % WriteTecplot( 'mesh' )
!#endif
!
!         ENDIF
!
!         nElems = mesh % elements % nElements
!         nProc  = params % nProcX*params % nProcY*params % nProcZ
!         IF( nProc == 0 )THEN
!            nProc = 1
!         ENDIF
!         ALLOCATE( materials(1:nElems), &
!                   globalToLocal(1:nElems,1:2), &
!                   nElPerProc(0:nProc-1), &
!                   nodeLogic(1:mesh % nodes % nNodes, 0:nProc-1), &
!                   nNodePerProc(0:nProc-1), &
!                   globalToLocalNode(1:mesh % nodes % nNodes, 0:nProc-1), &
!                   nLocMPI(0:nProc-1) )
!                   
!         materials     = 0.0_prec
!         globalToLocal = 0
!         nElPerProc    = 0
!         nodeLogic     = 0
!         nNodePerProc  = 0
!         
!         ALLOCATE( partitions(1:nElems) )
!         ALLOCATE( faceProcCount(1:mesh % faces % nFaces), &
!                   faceProcOwners(1:mesh % faces % nFaces,1:2), &
!                   faceBoundaryIDs(1:mesh % faces % nFaces,1:2), &
!                   faceProcTable(1:mesh % faces % nFaces, 0:nProc-1) )
!         faceProcCount   = 0
!         faceProcOwners  = -1
!         faceBoundaryIDs = 0
!         faceProcTable   = -5000
!
!         CALL mesh % PartitionStructuredElementsAndNodes( params, partitions, nElPerProc, &
!                                                          globalToLocal, nodeLogic, nNodePerProc, &
!                                                          globalToLocalNode, nProc )
!
!   
!            
!         ! Now we generate the local mesh for each process
!
!         DO procID = 0, nProc-1
!            CALL procMesh % Build( nNodePerProc(procID), &
!                                   nElPerProc(procID), &
!                                   1, params % polyDeg )
!! ----------------------------------------------------------------------------- !
!
!            DO nID = 1, mesh % nodes % nNodes
!            
!               ! We only assign node information if it is "owned" by this process
!               IF( nodeLogic(nID,procID) == 1 )THEN 
!
!                  localID = globalToLocalNode( nID, procID )
!                  procMesh % nodes % nodeID(localID)   = nID
!                  procMesh % nodes % nodeType(localID) = mesh % nodes % nodeTYPE(nID)
!                  procMesh % nodes % x(1:3,localID)    = mesh % nodes % x(1:3,nID)
!
!               ENDIF
!               
!            ENDDO
!   
!! ----------------------------------------------------------------------------- !
!
!            ! Set the element-node connectivity
!            DO iEl = 1, mesh % elements % nElements
!               ! We only assign element information if it is "owned" by this process
!               IF( globalToLocal(iEl,2) == procID )THEN
!                  elID = globalToLocal(iEl,1) ! Obtain the local element ID from the global-to-local array
!                  
!                  DO iNode = 1, 8
!       
!                    ! First we grab the global node ID for this element's corner node
!                    nID = mesh % elements % nodeIDs(iNode,iEl)
!                    ! THEN we convert this to a local node ID for this process
!                    localID = globalToLocalNode( nID, procID )
!                    ! And now we assign the local Node ID to the local element's corner node ID's
!                    procMesh % elements % nodeIDs(iNode,elID)  = localID
!       
!                  ENDDO
!                  procMesh % elements % elementID(elID) = iEl
!                     
!                  ! Set the geometry
!                  procMesh % elements % nHat(:,:,:,:,elID)   = mesh % elements % nHat(:,:,:,:,iEl)
!                  procMesh % elements % xBound(:,:,:,:,elID) = mesh % elements % xBound(:,:,:,:,iEl)
!                  procMesh % elements % x(:,:,:,:,elID)      = mesh % elements % x(:,:,:,:,iEl)
!                  procMesh % elements % J(:,:,:,elID)        = mesh % elements % J(:,:,:,iEl)
!                  procMesh % elements % dxds(:,:,:,elID)     = mesh % elements % dxds(:,:,:,iEl)
!                  procMesh % elements % dxdp(:,:,:,elID)     = mesh % elements % dxdp(:,:,:,iEl)
!                  procMesh % elements % dxdq(:,:,:,elID)     = mesh % elements % dxdq(:,:,:,iEl)
!                  procMesh % elements % dyds(:,:,:,elID)     = mesh % elements % dyds(:,:,:,iEl)
!                  procMesh % elements % dydp(:,:,:,elID)     = mesh % elements % dydp(:,:,:,iEl)
!                  procMesh % elements % dydq(:,:,:,elID)     = mesh % elements % dydq(:,:,:,iEl)
!                  procMesh % elements % dzds(:,:,:,elID)     = mesh % elements % dzds(:,:,:,iEl)
!                  procMesh % elements % dzdp(:,:,:,elID)     = mesh % elements % dzdp(:,:,:,iEl)
!                  procMesh % elements % dzdq(:,:,:,elID)     = mesh % elements % dzdq(:,:,:,iEl)
!                  procMesh % elements % Ja(:,:,:,:,:,elID)   = mesh % elements % Ja(:,:,:,:,:,iEl)
!       
!               ENDIF
!               
!            ENDDO
!   
!! ----------------------------------------------------------------------------- !
!
!            npFaces = 0 
!            nBe     = 0
!            iFaceLocal = 0
!            nMPI       = 0
!            DO iFace = 1, mesh % faces % nFaces
!
!               e1 = mesh % faces % elementIDs(1,iFace)
!               e2 = mesh % faces % elementIDs(2,iFace)
!
!               p1 = globalToLocal(e1,2)
!               IF( e2 > 0 )THEN
!                  p2 = globalToLocal(e2,2)
!               ELSE
!                  p2  = p1
!               ENDIF
!               IF( p1 == procID .OR. p2 == procID )THEN
!
!                  faceProcTable(iFace,p1) = 1
!                  faceProcTable(iFace,p2) = 1
!                  iFaceLocal = iFaceLocal+1
!
!                  faceProcCount(iFace) = faceProcCount(iFace) + 1
!                  faceProcOwners(iFace,faceProcCount(iFace)) = procID 
!
!                  npFaces = npFaces + 1
!
!                  IF( e2 < 0 .AND. p2 == p1 )THEN
!
!                     nBe = nBe + 1
!                     faceBoundaryIDs(iFace, faceProcCount(iFace)) = nBe
!
!                  ELSEIF( p2 /= p1 .AND. p1 == procID )THEN
!
!                     nMPI = nMPI + 1 
!                     nBe = nBe + 1
!                     faceBoundaryIDs(iFace, faceProcCount(iFace)) = nBe
!
!                  ELSEIF( p2 /= p1 .AND. p2 == procID )THEN
!
!                     nMPI = nMPI + 1 
!                     nBe = nBe + 1
!                     faceBoundaryIDs(iFace, faceProcCount(iFace)) = nBe
!
!                  ENDIF
!
!               ENDIF
!
!            ENDDO
!
!            WRITE(msg,'(A,I5)') 'Process ID : ',procID
!            INFO( TRIM(msg) )
!            WRITE(msg,'(A,I5)') 'Number of MPI Communications : ',nMPI
!            INFO( TRIM(msg) )
!            WRITE(msg,'(A,I5)') 'Number of faces : ',npFaces
!            INFO( TRIM(msg) )
!            WRITE(msg,'(A,I5)') 'Number of boundary faces : ',nBe
!            INFO( TRIM(msg) )
!            WRITE(msg,'(A,I5)') 'Number of elements : ',nElPerProc(procID)
!            INFO( TRIM(msg) )
!
!            CALL procMesh % faces % Trash( )
!            CALL procMesh % faces % Build( npFaces, params % polyDeg ) 
!            CALL procMesh % boundaryMap % Build( nBe, nMPI )
!
!            iFaceLocal = 0
!            nBe        = 0
!            nMPI       = 0
!            DO iFace = 1, mesh % faces % nFaces
!
!               e1 = mesh % faces % elementIDs(1,iFace)
!               e2 = mesh % faces % elementIDs(2,iFace)
!
!               p1 = globalToLocal(e1,2)
!               IF( e2 > 0 )THEN
!                  p2 = globalToLocal(e2,2)
!               ELSE
!                  p2  = p1
!               ENDIF
!
!               IF( p1 == procID .OR. p2 == procID )THEN
!
!                  iFaceLocal = iFaceLocal + 1
!
!
!                  IF( p2 == p1 .AND. e2 > 0 )THEN ! Internal Face
!                    
!                     procMesh % faces % faceID(iFaceLocal)         = iFace
!                     procMesh % faces % elementIDs(1,iFaceLocal)   = globalToLocal( e1, 1 )
!                     procMesh % faces % elementIDs(2,iFaceLocal)   = globalToLocal( e2, 1 )
!                     procMesh % faces % elementSides(1,iFaceLocal) = mesh % faces % elementSides(1,iFace)
!                     procMesh % faces % elementSides(2,iFaceLocal) = mesh % faces % elementSides(2,iFace)
!                     procMesh % faces % iStart(iFaceLocal)         = mesh % faces % iStart(iFace)
!                     procMesh % faces % iInc(iFaceLocal)           = mesh % faces % iInc(iFace)
!                     procMesh % faces % jStart(iFaceLocal)         = mesh % faces % jStart(iFace)
!                     procMesh % faces % jInc(iFaceLocal)           = mesh % faces % jInc(iFace)
!                     procMesh % faces % swapDimensions(iFaceLocal) = mesh % faces % swapDimensions(iFace)
!                     procMesh % faces % iMap(:,:,iFaceLocal)       = mesh % faces % iMap(:,:,iFace)
!                     procMesh % faces % jMap(:,:,iFaceLocal)       = mesh % faces % jMap(:,:,iFace)
!
!                  ELSEIF( p2 == p1 .AND. e2 < 0 )THEN ! Physical boundary
!
!                     procMesh % faces % faceID(iFaceLocal)         = iFace
!                     procMesh % faces % elementIDs(1,iFaceLocal)   = globalToLocal( e1, 1 )
!                     procMesh % faces % elementIDs(2,iFaceLocal)   = e2
!                     procMesh % faces % elementSides(1,iFaceLocal) = mesh % faces % elementSides(1,iFace)
!                     procMesh % faces % elementSides(2,iFaceLocal) = mesh % faces % elementSides(2,iFace)
!                     procMesh % faces % iStart(iFaceLocal)         = mesh % faces % iStart(iFace)
!                     procMesh % faces % iInc(iFaceLocal)           = mesh % faces % iInc(iFace)
!                     procMesh % faces % jStart(iFaceLocal)         = mesh % faces % jStart(iFace)
!                     procMesh % faces % jInc(iFaceLocal)           = mesh % faces % jInc(iFace)
!                     procMesh % faces % swapDimensions(iFaceLocal) = mesh % faces % swapDimensions(iFace)
!                     procMesh % faces % iMap(:,:,iFaceLocal)       = mesh % faces % iMap(:,:,iFace)
!                     procMesh % faces % jMap(:,:,iFaceLocal)       = mesh % faces % jMap(:,:,iFace)
!
!                     nBe = nBe + 1
!                     procMesh % faces % boundaryID(iFaceLocal) = -nBe
!
!                     ! Physical boundary ID's will be set to negative for physical boundaryMap
!                     ! so that boundaryMap requiring communication can be separated from physical
!                     ! boundaryMap. This will allow more operations to be run concurrently with
!                     ! communication.
!                     procMesh % boundaryMap % boundaryIDs(nBe)       = iFaceLocal
!                     procMesh % boundaryMap % boundaryGlobalIDs(nBe) = iFace
!                     procMesh % boundaryMap % extProcIDs(nBe)        = procID
!
!                  ELSEIF( p2 /= p1 .AND. procID == p1 )THEN ! MPI Boundary
!                
!                     procMesh % faces % faceID(iFaceLocal)         = iFace
!                     procMesh % faces % elementIDs(1,iFaceLocal)   = globalToLocal( e1, 1 )
!                     procMesh % faces % elementIDs(2,iFaceLocal)   = -e2
!                     procMesh % faces % elementSides(1,iFaceLocal) = mesh % faces % elementSides(1,iFace)
!                     procMesh % faces % elementSides(2,iFaceLocal) = mesh % faces % elementSides(2,iFace)
!                     procMesh % faces % iStart(iFaceLocal)         = mesh % faces % iStart(iFace)
!                     procMesh % faces % iInc(iFaceLocal)           = mesh % faces % iInc(iFace)
!                     procMesh % faces % jStart(iFaceLocal)         = mesh % faces % jStart(iFace)
!                     procMesh % faces % jInc(iFaceLocal)           = mesh % faces % jInc(iFace)
!                     procMesh % faces % swapDimensions(iFaceLocal) = mesh % faces % swapDimensions(iFace)
!                     procMesh % faces % iMap(:,:,iFaceLocal)       = mesh % faces % iMap(:,:,iFace)
!                     procMesh % faces % jMap(:,:,iFaceLocal)       = mesh % faces % jMap(:,:,iFace)
!
!                     nBe = nBe + 1
!                     procMesh % faces % boundaryID(iFaceLocal) = nBe
!
!                     ! Boundary ID's associated with MPI boundaryMap are reported as positive INTEGERs
!                     ! so that they can be distinguished from physical boundaryMap that do not require
!                     ! communication with neighboring processes. 
!                     procMesh % boundaryMap % boundaryIDs(nBe) = iFaceLocal
!                     procMesh % boundaryMap % boundaryGlobalIDs(nBe) = iFace
!                     procMesh % boundaryMap % extProcIDs(nBe)  = p2
!
!
!                  ELSEIF( p2 /= p1 .AND. procID == p2 )THEN ! MPI Boundary
!                  
!                  
!                     procMesh % faces % faceID(iFaceLocal)         = iFace
!                     procMesh % faces % elementIDs(1,iFaceLocal)   = globalToLocal( e2, 1 )
!                     procMesh % faces % elementIDs(2,iFaceLocal)   = -e1
!                     procMesh % faces % elementSides(1,iFaceLocal) = mesh % faces % elementSides(2,iFace)
!                     procMesh % faces % elementSides(2,iFaceLocal) = mesh % faces % elementSides(1,iFace)
!                     procMesh % faces % iStart(iFaceLocal)         = mesh % faces % iStart(iFace)
!                     procMesh % faces % iInc(iFaceLocal)           = mesh % faces % iInc(iFace)
!                     procMesh % faces % jStart(iFaceLocal)         = mesh % faces % jStart(iFace)
!                     procMesh % faces % jInc(iFaceLocal)           = mesh % faces % jInc(iFace)
!                     procMesh % faces % swapDimensions(iFaceLocal) = mesh % faces % swapDimensions(iFace)
!                     procMesh % faces % iMap(:,:,iFaceLocal)       = mesh % faces % iMap(:,:,iFace)
!                     procMesh % faces % jMap(:,:,iFaceLocal)       = mesh % faces % jMap(:,:,iFace)
!                  
!                     nBe = nBe + 1
!                     procMesh % faces % boundaryID(iFaceLocal) = nBe
!
!                     ! Boundary ID's associated with MPI boundaryMap are reported as positive INTEGERs
!                     ! so that they can be distinguished from physical boundaryMap that do not require
!                     ! communication with neighboring processes. 
!                     procMesh % boundaryMap % boundaryIDs(nBe) = iFaceLocal
!                     procMesh % boundaryMap % boundaryGlobalIDs(nBe) = iFace
!                     procMesh % boundaryMap % extProcIDs(nBe)  = p1
!
!
!                  ENDIF
!               ENDIF
!
!            ENDDO
!            
!! ----------------------------------------------------------------------------- !
!
!            WRITE( pIDChar, '(I4.4)' ) procID
!#ifdef TECPLOT
!            CALL procMesh % WriteTecplot( 'mesh.'//pIDChar )
!#endif
!            CALL procMesh % WriteSELFMeshFile( TRIM(params % SELFMeshFile)//'.'//pIDChar )
!            
!            CALL procMesh % Trash( )
!
!         ENDDO
!   
!
!      ! For visualizing the decomposition
!      materials = REAL( partitions, prec )
!#ifdef TECPLOT
!      CALL mesh % WriteMaterialTecplot( materials )
!#endif
!      
!
!      ! Clean up memory !
!      DEALLOCATE( materials, &
!                  faceProcCount, &
!                  faceProcOwners, &
!                  faceBoundaryIDs, &
!                  nLocMPI )
!      
!      CALL mesh % Trash( )
!      CALL nodal % Trash( )
!
!      ENDIF
!
!      INFO('End')
!
! END SUBROUTINE StructuredMeshGenerator_3D

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
