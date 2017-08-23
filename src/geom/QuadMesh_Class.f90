! QuadMesh_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file QuadMesh_Class.f90
!! Contains the \ref QuadMesh_Class module, and <BR>
!! defines the QuadMesh data-structure.

!> \defgroup QuadMesh_Class QuadMesh_Class 
!! This module defines the QuadMesh data-structure and its associated routines.
!!
!! There are some routines derived from D.A. Kopriva, 2009, "Implementing Spectral Methods for 
!! Partial Differential Equations: Algorithms for Scientists and Engineers"
!!
MODULE QuadMesh_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE LinkedList_Class
USE HashTable_Class
! src/interp/
USE Quadrature
USE Lagrange_Class
! src/geom/
USE Curve_Class
USE MappedGeometry_2D_Class
USE QuadElement_Class
USE Edge_Class
USE Node_Class

IMPLICIT NONE

!> \addtogroup QuadMesh_Class 
!! @{

!> \struct QuadMesh
!!  The QuadMesh data structure defines attributes needed to describe a conformal unstructured
!!  spectral element mesh. 
!!
!!  The QuadMesh data-structure brings together nodes, elements, edges, and mapped-geometry 
!!  data-structures to completely describe a spectral element mesh. Type-bound procedures are 
!!  provided for filling in mesh connectivity information including element-to-element neighbors, 
!!  edge-to-element IDs, element-to-node IDs, and edge-to-node IDs.
!!
!! <H2> QuadMesh </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> nElems <td> INTEGER  <td> Number of elements in the mesh
!!       <tr> <th> nNodes <td> INTEGER <td>  Number of nodes in the mesh
!!       <tr> <th> nEdges <td> INTEGER <td>  Number of edges in the mesh
!!       <tr> <th> elements(1:nElems) <td> QuadElement <td> Element-to-node and element-to-neighbor
!!                                                           information
!!       <tr> <th> geom(1:nElems) <td> MappedGeometry_2D <td> Computational-to-physical space
!!                                                            mapping information.
!!       <tr> <th> nodes(1:nNodes) <td> Node <td> Node ID's with node position
!!       <tr> <th> edges(1:nEdges) <td> Edge <td> Edge-to-node and edge-to-element information
!!       <tr> <th> cornerMap(1:2,1:4) <td> INTEGER <td> "Convenience array" for relating the 
!!                                          computational coordinates to local corner-node ID's
!!                                          (primarily for use use GAUSS_LOBATTO quadratures)
!!       <tr> <th> sideMap(1:4) <td> INTEGER <td> "Convenience array" for relating the 
!!                                          computational coordinates to local side ID's
!!                                          (primarily for use use GAUSS_LOBATTO quadratures)
!!       <tr> <th> edgeMap(1:2,1:4) <td> INTEGER <td>  "Convenience array" for relating the 
!!                                          local edge ID's to local corner-node ID's
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref QuadMesh_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_QuadMesh
!!       <tr> <th> Trash <td> Trash_QuadMesh
!!       <tr> <th> ConstructElementNeighbors <td> ConstructElementNeighbors_QuadMesh
!!       <tr> <th> ConstructEdges <td> ConstructEdges_QuadMesh
!!       <tr> <th> ScaleTheMesh <td> ScaleTheMesh_QuadMesh
!!       <tr> <th> LoadDefaultMesh <td> LoadDefaultMesh_QuadMesh
!!       <tr> <th> ConstructLinearMeshFromOther <td> ConstructLinearMeshFromOther_QuadMesh
!!       <tr> <th> ReadSpecMeshFile <td> ReadSpecMeshFile_QuadMesh
!!       <tr> <th> ReadPeaceMeshFile <td> ReadPeaceMeshFile_QuadMesh
!!       <tr> <th> WritePeaceMeshFile <td> WritePeaceMeshFile_QuadMesh
!!       <tr> <th> WriteTecplot <td> WriteTecplot_QuadMesh
!!       <tr> <th> WriteMaterialTecplot <td> WriteMaterialTecplot_QuadMesh
!!    </table>
!!

!>@}
   TYPE QuadMesh 
      INTEGER                                :: nElems, nNodes, nEdges
      TYPE( QuadElement ), ALLOCATABLE       :: elements(:)
      TYPE( MappedGeometry_2D ), ALLOCATABLE :: geom(:)
      TYPE( Node ), ALLOCATABLE              :: nodes(:)  
      TYPE( Edge ), ALLOCATABLE              :: edges(:)
      INTEGER, ALLOCATABLE                   :: nodeToElement(:,:,:) 
      INTEGER, ALLOCATABLE                   :: nodalValence(:)
      INTEGER                                :: cornerMap(1:2,1:4) 
      INTEGER                                :: sideMap(1:4) 
      INTEGER                                :: edgeMap(1:2,1:4) 

      CONTAINS

      PROCEDURE :: Initialize => Initialize_QuadMesh
      PROCEDURE :: Trash => Trash_QuadMesh
       
      PROCEDURE :: ConstructElementNeighbors    => ConstructElementNeighbors_QuadMesh
      PROCEDURE :: ConstructEdges               => ConstructEdges_QuadMesh
      PROCEDURE :: GetNodeToElementConnectivity => GetNodeToElementConnectivity_QuadMesh
      PROCEDURE :: ScaleTheMesh                 => ScaleTheMesh_QuadMesh
      PROCEDURE :: LoadDefaultMesh              => LoadDefaultMesh_QuadMesh
      PROCEDURE :: ConstructLinearMeshFromOther => ConstructLinearMeshFromOther_QuadMesh       

      PROCEDURE :: ReadSpecMeshFile     => ReadSpecMeshFile_QuadMesh
      PROCEDURE :: ReadPeaceMeshFile    => ReadPeaceMeshFile_QuadMesh
      PROCEDURE :: WritePeaceMeshFile   => WritePeaceMeshFile_QuadMesh
      PROCEDURE :: WriteTecplot         => WriteTecplot_Quadmesh
      PROCEDURE :: WriteMaterialTecplot => WriteMaterialTecplot_Quadmesh
       
   END TYPE QuadMesh

 INTEGER, PRIVATE, PARAMETER    :: nXElemDefault = 5
 INTEGER, PRIVATE, PARAMETER    :: nYElemDefault = 5
 INTEGER, PRIVATE, PARAMETER    :: BoundaryFlagDefault = NO_NORMAL_FLOW
 REAL(prec), PRIVATE, PARAMETER :: dXDefault = ONE/(nXElemDefault)
 REAL(prec), PRIVATE, PARAMETER :: dYDefault = ONE/(nYElemDefault)


 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup QuadMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_QuadMesh 
!! Initializes the attributes and "convenience" arrays of the QuadMesh data-structure.
!! 
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadMesh) :: this <BR>
!! <B>INTEGER</B>        :: nNodes, nElems, nEdges, N <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize( nNodes, nElems, nEdges, N ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myQuadMesh <td> QuadMesh <td> 
!!   <tr> <td> in <th> nNodes <td> INTEGER <td> The number of nodes in the mesh
!!   <tr> <td> in <th> nElems <td> INTEGER <td> The number of elements in the mesh
!!   <tr> <td> in <th> nEdges <td> INTEGER <td> The number of edges in the mesh 
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree for the geometry within each element
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_QuadMesh( myQuadMesh, nNodes, nElems, nEdges, N )

   IMPLICIT NONE
   CLASS(QuadMesh), INTENT(out) :: myQuadMesh
   INTEGER, INTENT(in)          :: nNodes, nElems, nEdges, N
   !LOCAL
   INTEGER :: i

      myQuadmesh % sideMap(1:4)      = (/ 0, N, N, 0 /)
      myQuadmesh % cornerMap(1, 1:4) = (/ 0, N, N, 0 /)
      myQuadmesh % cornerMap(2, 1:4) = (/ 0, 0, N, N /)
      myQuadmesh % edgeMap(1, 1:4)   = (/ 1, 2, 4, 1 /)
      myQuadmesh % edgeMap(2, 1:4)   = (/ 2, 3, 3, 4 /)

      myQuadmesh % nNodes = nNodes
      myQuadmesh % nElems = nElems
      myQuadMesh % nEdges = nEdges
      
      ALLOCATE( myQuadmesh % nodes(1:nNodes), &
                myQuadmesh % edges(1:nEdges), &
                myQuadMesh % geom(1:nElems), &
                myQuadMesh % elements(1:nElems), &
                myQuadMesh % nodeToElement(1:2,0:maxNodalValence,1:nNodes), &
                myQuadMesh % nodalValence(1:nNodes) )
                
      myQuadMesh % nodeToElement = 0
      myQuadMesh % nodalValence  = 0
      
      DO i = 1, myQuadmesh % nNodes
         CALL myQuadmesh % nodes(i) % Initialize( )
      ENDDO 
      DO i = 1, myQuadMesh % nElems
         CALL myQuadMesh % elements(i) % Initialize( )
         CALL myQuadMesh % geom(i) % Initialize( N )
      ENDDO
      DO i = 1, myQuadMesh % nEdges
         CALL myQuadMesh % edges(i) % Initialize( )
      ENDDO

 END SUBROUTINE Initialize_QuadMesh
!
!> \addtogroup QuadMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_QuadMesh 
!! Frees memory associated with each of the attributes of the QuadMesh data structure
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadMesh) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myQuadMesh <td> QuadMesh <td> On output, the memory associated with 
!!                         attributes of the QuadMesh structure have been freed 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_QuadMesh( myQuadMesh )

   IMPLICIT NONE
   CLASS(QuadMesh), INTENT(inout) :: myQuadMesh
  ! LOCAL
   INTEGER :: i, nEl
      
      nEl = myQuadMesh % nElems
      DO i = 1, nEl
         CALL myQuadMesh % geom(i) % Trash( )
      ENDDO

      DEALLOCATE( myQuadMesh % nodes, &
                  myQuadMesh % edges, &
                  myQuadMesh % elements, &
                  myQuadMesh % nodeToElement, &
                  myQuadMesh % nodalValence, &
                  myQuadMesh % geom )
     
 END SUBROUTINE Trash_QuadMesh
!
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup QuadMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R ConstructEdges
! 
!> \fn ConstructEdges_QuadMesh  
!! Uses the element-to-node connectivity to construct all of the unique edges in a mesh
!! 
!! To construct the unique edges in the mesh, we cycle over all of the elements and over each 
!! element-side. An edge is identified by the two nodes that mark its endpoints - these we obtain
!! from the element to node connectivity and the "edgeMap" convenience array. The two nodes are stored
!! as "keys" in a hash-table with the smaller node ID being the primary key. If the keys do not
!! exist in the hash-table, a new edge is created and the element that caused the edge creation
!! is designated the primary element. If the key already exists in the hash-table, then the secondary
!! element information is filled in. The node ordering for the secondary element edge determines
!! the "start" and "inc" attributes. By default, the secondary element information is set to a
!! physical boundary condition flag. This is overwritten in the "ReadSpecMesh" and "ReadPeaceMesh"
!! routines.
!! 
!! This routine is derived from Alg. 148 on pg. 383 of D.A. Kopriva, 2009.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadMesh) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ConstructEdges( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myQuadMesh <td> QuadMesh <td> 
!!                         On <B>input</B>, the element and node information has been filled in, <BR>
!!                         On <B>output</B>, the unique edges have been identified and the edge 
!!                         information has been filled in. 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
SUBROUTINE ConstructEdges_QuadMesh( myQuadMesh )

   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   ! LOCAL
   TYPE( HashTable ) :: edgeTable
   INTEGER :: nEls, nNodes, iEl, nEdges, k  
   INTEGER :: l1, l2, startID, endID, key1, key2
   INTEGER :: e1, s1, edgeID, n1

      nNodes = myQuadMesh % nNodes  
      nEls   = myQuadMesh % nElems
      nEdges = 0

      ! First, just count the number of edges
      CALL edgeTable % Build( nNodes )

      DO iEl = 1, nEls ! Loop over the elements in the mesh
         DO k = 1, 4 ! Loop over the sides of each element

            l1 = myQuadMesh % edgeMap(1,k) ! starting local node for this edge
            l2 = myQuadMesh % edgeMap(2,k) ! ending local node for this edge
            
            startID = myQuadMesh % elements(iEl) % nodeIDs(l1)
            endID   = myQuadMesh % elements(iEl) % nodeIDs(l2)
            
            key1 = min( startID, endID )
            key2 = max( startID, endID )

            ! Add element to corner-node connectivity list
            IF( edgeTable % ContainsKeys( key1, key2 ) .EQV. .FALSE. )then ! this is a new edge
               
               ! Add the edge to the list
               nEdges = nEdges + 1
               CALL edgeTable % AddDataForKeys( nEdges, key1, key2 )
               
            ENDIF

         ENDDO ! k, Loop over the sides of each element
      ENDDO ! iEl, Loop over the elements in the mesh
 
      CALL edgeTable % Trash( ) ! TRASH the edgetable

      ! And rebuild it
      CALL edgeTable % Build( nNodes )
      
      ! Re-allocate space for the mesh edges

      DEALLOCATE( myQuadMesh % edges )

      ALLOCATE( myQuadMesh % edges( 1:nEdges ) )

      nEdges = 0 ! restart the edge counting

      DO iEl = 1, nEls ! Loop over the elements in the mesh
         DO k = 1, 4 ! Loop over the sides of each element

            l1 = myQuadMesh % edgeMap(1,k) ! starting local node for this edge
            l2 = myQuadMesh % edgeMap(2,k) ! ending local node for this edge

            startID = myQuadMesh % elements(iEl) % nodeIDs(l1)
            endID   = myQuadMesh % elements(iEl) % nodeIDs(l2)

            key1 = min( startID, endID )
            key2 = max( startID, endID )

            IF( edgeTable % ContainsKeys( key1, key2 )  )then ! this edge already exists
               
               !Get the edgeID
               CALL edgeTable % GetDataForKeys( edgeID, key1, key2 )
               ! Find the primary element and the starting node for this element's edge
               ! This is compared with the secondary element's starting node to infer
               ! the relative orientation of the two elements.
               e1 = myQuadMesh % edges(edgeID) % elementIDs(1)
               s1 = myQuadMesh % edges(edgeID) % elementSides(1)
               
               l1 = myQuadMesh % edgeMap(1,s1)
               n1 = myQuadMesh % elements(e1) % nodeIDs(l1)

               ! Set the secondary element information
               myQuadMesh % edges(edgeID) % elementIDs(2) = iEl
               
               !PRINT*, 'edgeID, primary, secondary:', edgeID, e1, iEl
               IF( startID == n1 ) then ! the elements are oriented the same direction
               
                  myQuadMesh % edges(edgeID) % elementSides(2) = k

               ELSE ! the elements are oriented in the opposite direction

                  ! For these edges, we mark the side ID as negative
                  myQuadMesh % edges(edgeID) % elementSides(2) = -k
                 
               ENDIF


            ELSE ! this is a new edge

               ! Add the edge to the list
               nEdges = nEdges + 1
               
               edgeID = nEdges
               CALL myQuadMesh % edges(edgeID) % Initialize()
               myQuadMesh % edges(edgeID) % elementIDs(1) = iEl
               myQuadMesh % edges(edgeID) % elementSides(1) = k
               myQuadMesh % edges(edgeID) % nodeIDs = (/startID, endID/)

               ! Default the secondary information
               myQuadMesh % edges(edgeID) % elementIDs(2) = BoundaryFlagDefault
               CALL edgeTable % AddDataForKeys( edgeID, key1, key2 )
               
            ENDIF
         enddo ! k, Loop over the sides of each element
        
      enddo ! iEl, Loop over the elements in the mesh

      CALL edgeTable % Trash( )
      myQuadMesh % nEdges = nEdges

 END SUBROUTINE ConstructEdges_QuadMesh
!
!> \addtogroup QuadMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R GetNodeToElementConnectivity 
! 
!> \fn GetNodeToElementConnectivity_QuadMesh 
!! Cycles over the elements and local nodes to generate the node to element connectivity list. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadMesh) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GetNodeToElementConnectivity( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myQuadMesh <td> QuadMesh <td>
!!                         On <B>input</B>, the mesh elements have been constructed, <BR>
!!                         On <B>output</B>, the node-to-element connectivity list is constructed.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetNodeToElementConnectivity_QuadMesh( myQuadMesh )

   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   ! LOCAL
   INTEGER :: nEls, nNodes, iEl, k, i 
   INTEGER :: nID

      nNodes = myQuadMesh % nNodes 
      nEls   = myQuadMesh % nElems

      myQuadMesh % nodeToElement = 0

      DO iEl = 1, nEls
         DO k = 1, 4 
            nID = myQuadMesh % elements(iEl) % nodeIDs(k)
            ! The 0th entry of the second dimension of the "nodeToElement" array
            ! provides the nodal valence for this node.
            myQuadMesh % nodalValence(nID) = myQuadMesh % nodalValence(nID) + 1
            i = myQuadMesh % nodalValence(nID)
            ! And here we add the element ID for the node-to-element connectivity list
            myQuadMesh % nodeToElement(1,i,nID) = iEl
            ! And here we add the local node ID for the node-to-element connectivity list
            myQuadMesh % nodeToElement(2,i,nID) = k
         ENDDO 
      ENDDO 

 END SUBROUTINE GetNodeToElementConnectivity_QuadMesh
!
!> \addtogroup QuadMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R ConstructElementNeighbors 
! 
!> \fn sConstructElementNeighbors_QuadMesh 
!! Uses the edge-to-element connectivity to construct the element neighbors attribute.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(myQuadMesh) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ConstructElementNeighbors( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myQuadMesh <td> QuadMesh <td>
!!                         On <B>input</B>, the edge-to-element information has been constructed, <BR>
!!                         On <B>output</B>, the element neighbors information has been filled in. 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ConstructElementNeighbors_QuadMesh( myQuadMesh )

   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   ! LOCAL
   INTEGER :: e(1:2), s(1:2), nEdges, iEdge
  
      nEdges = myQuadMesh % nEdges
      DO iEdge = 1, nEdges ! Loop over the edges

         e = myQuadMesh % edges(iEdge) % elementIDs
         s = myQuadMesh % edges(iEdge) % elementSides

         IF( e(2) > 0 )THEN
            myQuadMesh % elements(e(1)) % neighbors( s(1) )      = e(2)
            myQuadMesh % elements(e(2)) % neighbors( ABS(s(2)) ) = e(1)
         ENDIF

      ENDDO ! iEdge, Loop over the edges

     
 END SUBROUTINE ConstructElementNeighbors_QuadMesh
!
!> \addtogroup QuadMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R ScaleTheMesh_QuadMesh 
! 
!> \fn ScaleTheMesh  
!! Scales the element geometry and corner node positions by the provided x-scale and y-scale. 
!! 
!! This routine depend on <BR>
!!   Module \ref MappedGeometry_2D_Class, S/R ScaleGeometry_MappedGeometry_2D <BR>
!!   Module \ref Node_Class, S/R ScalePosition_Node 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadMesh)    :: this <BR>
!! <B>TYPE</B>(Lagrange) :: interp <BR>
!! <B>REAL</B>(prec)        :: xScale, yScale <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ScaleGeometry( interp, xScale, yScale ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myQuadMesh <td> QuadMesh <td> 
!!   <tr> <td> in <th> interp <td> Lagrange <td> 
!!   <tr> <td> in <th> xScale <td> REAL(prec) <td>
!!   <tr> <td> in <th> yScale <td> REAL(prec) <td>  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ScaleTheMesh_QuadMesh( myQuadMesh, interp, xScale, yScale  )
 
   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout) :: myQuadMesh
   TYPE( Lagrange ), INTENT(in)  :: interp
   REAL(prec), INTENT(in)           :: xScale, yScale
   ! LOCAL
   INTEGER :: nElems, nNodes, iEl, iNode
   
      nElems = myQuadMesh % nElems
      nNodes = myQuadMesh % nNodes 

      DO iEl = 1, nElems
         CALL myQuadMesh % geom(iEl) % ScaleGeometry( interp, xScale, yScale )
      ENDDO
      
      DO iNode = 1, nNodes
         CALL myQuadMesh % nodes(iNode) % ScalePosition( xScale, yScale, ZERO )
      ENDDO

 END SUBROUTINE ScaleTheMesh_QuadMesh
!
!> \addtogroup QuadMesh_Class
!! @{ 
! ================================================================================================ !
! S/R LoadDefaultMesh
! 
!> \fn LoadDefaultMesh_QuadMesh  
!! Constructs a "structured" spectral element mesh with nXelem-by-nYelem elements. 
!! 
!! This routine builds a mesh with nXelem elements in the x-direction and nYelem elements in the
!! y-direction. The mesh nodes are between [0,1]x[0,1]. After construction, the user can call
!! "ScaleTheMesh" to change the physical extents of the domain.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadMesh)    :: this <BR>
!! <B>TYPE</B>(Lagrange) :: interp <BR>
!! <B>INTEGER</B>           :: nXElem, nYElem <BR>
!!         .... <BR>
!!     <B>CALL</B> this % RoutineName( interp, nXElem, nYElem ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myQuadMesh <td> QuadMesh <td> On output, contains the "structured" mesh
!!                                                       in the unstructured spectral element
!!                                                       mesh format.
!!   <tr> <td> in <th> interp <td> Lagrange <td> 2-D interpolant that stores the computational
!!                                                  quadrature mesh for each element.
!!   <tr> <td> in <th> nXElem <td> INTEGER <td> The number of desired elements in the x-direction
!!   <tr> <td> in <th> nYElem <td> INTEGER <td> The number of desired elements in the y-direction
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE LoadDefaultMesh_QuadMesh( myQuadMesh, interp, nXelem, nYelem  )

   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout)  :: myQuadMesh
   TYPE( Lagrange ), INTENT(in)   :: interp
   INTEGER, INTENT(in)               :: nXelem, nYelem
   ! LOCAL
   TYPE( Curve ) :: elBoundCurves(1:4)

   REAL(prec) :: x, y, dxElem, dyElem
   REAL(prec) :: x1, x2, y1, y2
   REAL(prec), ALLOCATABLE :: xc(:,:), s(:)

   INTEGER :: nNodes, nElems, nEdges, gPolyDeg
   INTEGER :: nodes(1:4), nids(1:2), eID
   INTEGER :: s2, n1, n2
   INTEGER :: iEdge, iNode, iEl, iSide, iX, iY
      
      dxElem = ONE/nXElem
      dyElem = ONE/nYElem
      
      ! ** "Hard-wired" values for a structured mesh with no holes ** !
      nNodes   = (nXElem+1)*(nYElem+1)
      nElems   = (nXElem)*(nYElem)
      nEdges   = (nXElem)*(nYElem+1) + (nXElem+1)*(nYElem)
      gPolyDeg = 1
      ! ************************************************************************* !

      PRINT*, 'nNodes    : ', nNodes
      PRINT*, 'nElems    : ', nElems
      PRINT*, 'gPolyDeg  : ', gPolyDeg

      ! Generate the chebyshev points of order gPolyDeg
      ! These are the points used to define the parametric
      ! curves for the element boundaries

      ALLOCATE( s(0:gPolyDeg), xc(0:gPolyDeg,1:2) )
      CALL ChebyshevQuadrature( gPolyDeg, Gauss_Lobatto, s, xc(:,1) )! ** xc is a placeholder here only.

      ! ---- Initialize the quadrature mesh (empty) ---- !
      CALL myQuadMesh % Initialize( nNodes, nElems, nEdges, interp % N ) 
      
      ! ---- Read in the corner nodes ---- !
      DO iY = 1, nYElem+1
      
         y = dYElem*(REAL(iY-1,prec))
         DO iX = 1, nXElem+1
            
            iNode = iX + (iY-1)*(nXElem+1)
            x = dXElem*(REAL(iX-1,prec))
            myQuadMesh % nodes(iNode) % x = x
            myQuadMesh % nodes(iNode) % y = y
            
         ENDDO
      ENDDO
  
      ! Do the element information
      xc = ZERO
      ! Do the initial build for the parametric curves
      DO iSide = 1, 4
         CALL elBoundCurves(iSide) % Build( xc, s, gPolyDeg, 2 ) 
      ENDDO
     
      DO iY = 1, nYElem
         DO iX = 1, nXElem

            iEl = iX + (iY-1)*(nXElem)
            ! Calculate the global node IDs for this element.
            nodes(1) = iX + (iY-1)*(nXElem+1)    ! Southwest
            nodes(2) = iX + 1 + (iY-1)*(nXElem+1)! SouthEast
            nodes(3) = iX + 1 + (iY)*(nXElem+1)  ! NorthEast
            nodes(4) = iX + (iY)*(nXElem+1)      ! NorthWest
         
            DO iSide = 1, 4 ! Loop over the sides of the quads

               ! Initialize a straight curve for this side
               n1 = nodes( myQuadMesh % edgeMap(1,iSide) )
               n2 = nodes( myQuadMesh % edgeMap(2,iSide) )

               x1 = myQuadMesh % nodes(n1) % x
               y1 = myQuadMesh % nodes(n1) % y
               x2 = myQuadMesh % nodes(n2) % x
               y2 = myQuadMesh % nodes(n2) % y
               
               DO iNode = 0, gPolyDeg
                  xc(iNode,1) = x1 + (x2-x1)*HALF*( s(iNode) + ONE )
                  xc(iNode,2) = y1 + (y2-y1)*HALF*( s(iNode) + ONE )
               ENDDO
               elBoundCurves(iSide) % x = xc

            ENDDO
            myQuadMesh % elements(iEl) % nodeIDs(1:4) = nodes
            CALL myQuadMesh % geom(iEl) % GenerateMesh( interp, elBoundCurves )
            CALL myQuadMesh % geom(iEl) % GenerateMetrics( interp )
         ENDDO
      ENDDO ! iEl, cycle over the elements

      CALL myQuadMesh % ConstructEdges( )
      nEdges = myQuadMesh % nEdges
      PRINT*, 'nEdges    : ', nEdges
      
      ! Set the start and increment for the secondary element 
      DO iEdge = 1, nEdges
            s2 = myQuadMesh % edges(iEdge) % elementSides(2)
            IF(s2 < 0)THEN
               myQuadMesh % edges(iEdge) % start = interp % N-1
               myQuadMesh % edges(iEdge) % inc   = -1
            ELSE
               myQuadMesh % edges(iEdge) % start = 1
               myQuadMesh % edges(iEdge) % inc   = 1
            ENDIF
            eID = myQuadMesh % edges(iEdge) % elementIDs(2)
            IF( eID < 0 )THEN
               myQuadMesh % edges(iEdge) % elementIDs(2) = eID
               nids = myQuadMesh % edges(iEdge) % nodeIDs
               myQuadMesh % nodes(nids(1)) % nodeType = eID
               myQuadMesh % nodes(nids(2)) % nodeType = eID
            ENDIF
      ENDDO

      CALL myQuadMesh % GetNodeToElementConnectivity( )
      CALL myQuadMesh % ConstructElementNeighbors( )

      ! Clear up memory
      DEALLOCATE( s, xc )

      DO iSide = 1, 4
         CALL elBoundCurves(iSide) % Trash( )
      ENDDO
  
 END SUBROUTINE LoadDefaultMesh_QuadMesh
!
!> \addtogroup QuadMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R ConstructLinearMeshFromOther 
! 
!> \fn ConstructLinearMeshFromOther_QuadMesh  
!! Creates a new mesh by forming bilinear elements from the quadrature points of all of the elements
!! in a mesh.
!!
!! ** Sensible output is produced only for Gauss-Lobatto quadrature; this is the only SELF-supported
!!    quadrature that contains integration nodes that include element boundaries.
!! Two integer arrays are also provided for mapping solution storage structures between the two
!! meshes.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadMesh) :: newMesh, oldMesh <BR>
!! <B>INTEGER</B>        :: newToOld(:,:,:,:) <BR>
!! <B>INTEGER</B>        :: oldToNew(:,:,:,:) <BR>
!!         .... <BR>
!!     <B>CALL</B> newMesh % ConstructLinearMeshFromOther( oldMesh, newToOld, oldToNew ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> linMesh <td> QuadMesh <td> Spectral element mesh comprised of bilinear 
!!                                                 elements that are generated from a higher order mesh
!!   <tr> <td> in <th> otherMesh <td> QuadMesh <td> Higher order mesh from which to construct
!!                                                  the bilinear spectral element mesh
!!   <tr> <td> out <th> linToOther(:,:,:,:) <td> INTEGER <td> 4-D array that provides the mapping
!!                                               of a solution storage structure from the bilinear
!!                                               spectral element mesh to the higher order incoming
!!                                               mesh
!!   <tr> <td> out <th> otherToLin(:,:,:,:) <td> INTEGER <td> 4-D array that provides the mapping
!!                                               of a solution storage structure from the higher
!!                                               order incoming mesh to the bilinear spectral 
!!                                               element mesh
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ConstructLinearMeshFromOther_QuadMesh( linMesh, otherMesh, linToOther, otherToLin )

   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(inout)    :: linMesh
   TYPE( QuadMesh ), INTENT(in)        :: otherMesh
   INTEGER, INTENT(inout), ALLOCATABLE :: linToOther(:,:,:,:), otherToLin(:,:,:,:)
   ! LOCAL
   TYPE( Lagrange )     :: interp
   TYPE( Curve )           :: elBoundCurves(1:4)
   REAL(prec)              :: x, y, xc, yc
   REAL(prec)              :: x1, x2, y1, y2
   REAL(prec), ALLOCATABLE :: cnodelocs(:,:)
   REAL(prec), ALLOCATABLE :: xx(:,:),s(:)
   INTEGER, ALLOCATABLE    :: uniqueNids(:)
   INTEGER :: nodes(1:4), nids(1:2)
   INTEGER :: nEl, nP, n1, n2, nEdges
   INTEGER :: nLinEl, N, nLinNodes, nUnique, inc, nID, eID
   INTEGER :: iNode, i, j, iS, iP, iEl, iSide, iEdge

      N  = otherMesh % geom(1) % N
      nP  = N
      nEl = otherMesh % nElems
      
      N      = 1 ! Polynomial order for the linear mesh
      nLinEl = nEl*N*nP

      ALLOCATE( linToOther(0:N,0:N,1:nLinEl,1:3), otherToLin(0:N,0:nP,1:nEl,1:3) )
      linToOther = 0
      otherToLin = 0

      nLinNodes = nEl*(N+1)*(nP+1)
      ALLOCATE( cnodelocs(1:nLinNodes,1:2), uniqueNids(1:nLinNodes) )

      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, N

               iNode = iS+1 + (N+1)*( iP + (nP+1)*(iEl-1) )
               x = otherMesh % geom(iEl) % x(iS,iP)
               y = otherMesh % geom(iEl) % y(iS,iP)
               cnodelocs(iNode,1) = x
               cnodelocs(iNode,2) = y

            ENDDO
         ENDDO
      ENDDO 

      nUnique = 0
      ! Now we determine the unique nodes by comparing their location. The list is traversed in 
      ! increasing order of node ID. The node is compared to all of the previous nodes. If this node
      ! is different from all of the previously identified node, the number of unique nodes is
      ! incremented and the node is assigned that number as it's ID. An integer array is used
      ! to map the global node list ID to the unique node ID, so that we can establish a mapping
      ! between the linear mesh and the original mesh.
      DO i = 1, nLinNodes

         x = cnodelocs(i,1)
         y = cnodelocs(i,2)

         inc = 1 ! reset the unique node ID increment
         DO j = 1,i-1

            ! Grab the nodes to compare with
            xc = cnodelocs(j,1)
            yc = cnodelocs(j,2)

            ! If the current node (node "i") does not match a previous node (node "j"), then we record
            ! the unique ID counter in the integer array "uniqueNids". 
            IF( AlmostEqual(x,xc) .AND. AlmostEqual(y,yc) )THEN
               uniqueNids(i) = uniqueNids(j)
               inc = 0
               EXIT
            ENDIF

         ENDDO

         nUnique = nUnique + inc
         IF( inc == 1 )THEN
            uniqueNids(i) = nUnique 
         ENDIF

      ENDDO
 
      PRINT*, 'S/R ConstructLinearMeshFromOther : Found ', nUnique, ' unique nodes in the linear mesh.'
      CALL linMesh % Initialize( nUnique, nLinEl, 1, N ) 

      ALLOCATE( s(0:N), xx(0:N,1:2) )
      s(0) = -ONE
      s(1) = ONE

      CALL interp % Build( N, N, s, s )
      ! Now the linear mesh construction and mapping can begin.
     
      ! First we start with the nodes
      DO iEl = 1, nEl
         DO iP = 0, nP
            DO iS = 0, N

               iNode = iS+1 + (N+1)*( iP + (nP+1)*(iEl-1) )
               nID   = uniqueNids(iNode)  ! Convert iNode to a unique node ID             
               x     = cnodelocs(iNode,1)
               y     = cnodelocs(iNode,2)
               linmesh % nodes(nID) % x = x
               linmesh % nodes(nID) % y = y 

            ENDDO
         ENDDO
      ENDDO

      ! Do the element information
 
      xx = ZERO
      ! Do the initial build for the parametric curves
      DO iSide = 1, 4
         CALL elBoundCurves(iSide) % Build( xx, s, N, 2 ) 
      ENDDO
   
     
      DO iEl = 1, nEl
         DO iP = 1, nP
            DO iS = 1, N

               eID = iS + N*( (iP-1) + (iEl-1)*nP )
               ! Linear to original mesh mapping
               ! Southwest
               linToOther(0,0,eID,1) = iS-1 
               linToOther(0,0,eID,2) = iP-1
               linToOther(0,0,eID,3) = iEl
               ! Southeast
               linToOther(1,0,eID,1) = iS 
               linToOther(1,0,eID,2) = iP-1
               linToOther(1,0,eID,3) = iEl
               ! Northeast
               linToOther(1,1,eID,1) = iS 
               linToOther(1,1,eID,2) = iP
               linToOther(1,1,eID,3) = iEl
               ! NorthWest
               linToOther(0,1,eID,1) = iS-1 
               linToOther(0,1,eID,2) = iP
               linToOther(0,1,eID,3) = iEl

               ! Original to linear mesh mapping
               ! Southwest
               otherToLin(iS-1,iP-1,iEl,1) = 0 
               otherToLin(iS-1,iP-1,iEl,2) = 0
               otherToLin(iS-1,iP-1,iEl,3) = eID
               ! Southeast
               otherToLin(iS,iP-1,iEl,1) = 1 
               otherToLin(iS,iP-1,iEl,2) = 0
               otherToLin(iS,iP-1,iEl,3) = eID
               ! Northeast
               otherToLin(iS,iP,iEl,1) = 1 
               otherToLin(iS,iP,iEl,2) = 1
               otherToLin(iS,iP,iEl,3) = eID
               ! NorthWest
               otherToLin(iS-1,iP,iEl,1) = 0 
               otherToLin(iS-1,iP,iEl,2) = 1
               otherToLin(iS-1,iP,iEl,3) = eID
   
               ! Calculate the global node IDs for this element.
               nodes(1) = uniqueNids( iS + (N+1)*( iP - 1 + (nP+1)*(iEl-1) ) )   ! Southwest
               nodes(2) = uniqueNids( iS+1 + (N+1)*( iP - 1 + (nP+1)*(iEl-1) ) ) ! SouthEast
               nodes(3) = uniqueNids( iS+1 + (N+1)*( iP + (nP+1)*(iEl-1) ) )     ! NorthEast
               nodes(4) = uniqueNids( iS + (N+1)*( iP + (nP+1)*(iEl-1) ) )       ! NorthWest
         
               DO iSide = 1, 4 ! Loop over the sides of the quads

                  ! Initialize a straight curve for this side
                  n1 = nodes( linmesh % edgeMap(1,iSide) )
                  n2 = nodes( linmesh % edgeMap(2,iSide) )
                  x1 = linmesh % nodes(n1) % x
                  y1 = linmesh % nodes(n1) % y
                  x2 = linmesh % nodes(n2) % x
                  y2 = linmesh % nodes(n2) % y
               
                  DO iNode = 0, N
                     xx(iNode,1) = x1 + (x2-x1)*HALF*( s(iNode) + ONE )
                     xx(iNode,2) = y1 + (y2-y1)*HALF*( s(iNode) + ONE )
                  ENDDO
                  elBoundCurves(iSide) % x = xx
               ENDDO
               linmesh % elements(eID) % nodeIDs(1:4) = nodes
               linmesh % elements(eID) % elementID     = eID
               
               CALL linmesh % geom(eID) % GenerateMesh( interp, elBoundCurves )
               CALL linmesh % geom(eID) % GenerateMetrics( interp )
               
            ENDDO
         ENDDO
      ENDDO

      CALL linmesh % ConstructEdges( )
      nEdges = linmesh % nEdges
      PRINT*, 'nEdges    : ', nEdges
      
      ! Set the start and increment for the secondary element 
      DO iEdge = 1, nEdges
            iSide = linmesh % edges(iEdge) % elementSides(2)
            IF(iSide < 0)THEN
               linmesh % edges(iEdge) % start = N-1
               linmesh % edges(iEdge) % inc   = -1
            ELSE
               linmesh % edges(iEdge) % start = 1
               linmesh % edges(iEdge) % inc   = 1
            ENDIF
            
            eID = linmesh % edges(iEdge) % elementIDs(2)
            IF( eID < 0 )THEN
               linmesh % edges(iEdge) % elementIDs(2) = eID
               nids = linmesh % edges(iEdge) % nodeIDs
               linmesh % nodes(nids(1)) % nodeType = eID
               linmesh % nodes(nids(2)) % nodeType = eID
            ENDIF
      ENDDO

      CALL linmesh % GetNodeToElementConnectivity( )

      ! Clear up memory
      DEALLOCATE( s, xx )

      DO iSide = 1, 4
         CALL elBoundCurves(iSide) % Trash( )
      ENDDO
  
      
 END SUBROUTINE ConstructLinearMeshFromOther_QuadMesh
!
!
!==================================================================================================!
!--------------------------------- Mesh File I/O Routines -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup QuadMesh_Class
!! @{ 
! ================================================================================================ !
! S/R ReadSpecMeshFile 
! 
!> \fn ReadSpecMeshFile_QuadMesh  
!! Reads in a SpecMesh file (format ISM-v2). 
!! 
!! For the formatting requirements of a SpecMesh file, see the SpecMesh manual included with the
!! SELF. SpecMesh is software provided by David Kopriva and is included with the SELF. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadMesh)    :: this <BR>
!! <B>TYPE</B>(Lagrange) :: interp <BR>
!! <B>CHARACTER</B>         :: filename <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ReadSpecMeshFile( interp, filename ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myQuadMesh <td> QuadMesh <td>
!!   <tr> <td> in <th> interp <td> Lagrange <td> 2-D Lagrange interpolation structure that 
!!                                                  contains the computational quadrature mesh
!!                                                  for each element
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Name of the SpecMesh file
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ReadSpecMeshFile_QuadMesh( myQuadMesh, interp, filename )

   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(out)  :: myQuadMesh
   TYPE( Lagrange ), INTENT(in)    :: interp
   CHARACTER(*), INTENT(in)        :: filename
   ! LOCAL
   TYPE( Curve ) :: elBoundCurves(1:4)

   REAL(prec) :: x, y, z
   REAL(prec) :: x1, x2, y1, y2
   REAL(prec), ALLOCATABLE :: xc(:,:), s(:)

   INTEGER :: nNodes, nElems, nEdges, gPolyDeg
   INTEGER :: bFlags(1:4), nodes(1:4)
   INTEGER :: e1, s1, s2, n1, n2, n(1:2), e(1:2), si(1:2)
   INTEGER :: iEdge, iNode, iEl, iSide
   INTEGER :: fUnit, iC, jC
   INTEGER, ALLOCATABLE :: sideFlags(:,:)

   CHARACTER(20) :: ISMversion
   CHARACTER(100) :: edgeNames
   CHARACTER(100) :: thisEdge

      PRINT*, 'Mesh File : '//TRIM( filename )
      
      ! Get a new file unit
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE= TRIM( filename ), &
            FORM='FORMATTED',&
            STATUS='OLD' )

      
      ! ---- Read in the file version ----- !

      READ( fUnit, '(A20)' ) ISMversion

      PRINT*, 'Version   : '//TRIM( ISMversion )


      ! ---- Gather the number of nodes, number of elements, and number of edges ---- !
      
      READ( fUnit, * ) nNodes, nEdges, nElems, gPolyDeg

      PRINT*, 'nNodes    : ', nNodes
      PRINT*, 'nElems    : ', nElems
      PRINT*, 'nEdges    : ', nEdges
      PRINT*, 'gPolyDeg  : ', gPolyDeg

      ! Generate the chebyshev points of order gPolyDeg
      ! These are the points used to define the parametric
      ! curves for the element boundaries

      ALLOCATE( s(0:gPolyDeg), xc(0:gPolyDeg,1:2) )
      ALLOCATE( sideFlags(1:nElems,1:4) )
      CALL ChebyshevQuadrature( gPolyDeg, Gauss_Lobatto, s, xc(:,1) )! ** xc is a placeholder here only.

      ! ---- Initialize the quadrature mesh (empty) ---- !
  
      CALL myQuadMesh % Initialize( nNodes, nElems, nEdges, interp % N ) 
      
      ! ---- Read in the corner nodes ---- !

      DO iNode = 1, nNodes  ! Loop over the nodes in the file
         READ( fUnit, * ) x, y, z
         myQuadMesh % nodes(iNode) % x = x
         myQuadMesh % nodes(iNode) % y = y
      ENDDO


      ! ---- Read in the edge information ---- !
      DO iEdge = 1, nEdges
         READ( fUnit, * ) n, e, si 
      ENDDO
      ! ---- Read in the element information ---- !
 
      xc = ZERO
      ! Do the initial build for the parametric curves
      DO iSide = 1, 4
         CALL elBoundCurves(iSide) % Build( xc, s, gPolyDeg, 2 ) 
      ENDDO
   
      sideFlags = INTERIOR
      DO iEl = 1, nElems

         READ( fUnit, * ) nodes(1:4)
         READ( fUnit, * ) bFlags(1:4)
         
         DO iSide = 1, 4 ! Loop over the sides of the quads

            IF( bFlags(iSide) == 0 )then ! this is an interior edge

               ! Initialize a straight curve for this side
               n1 = nodes( myQuadMesh % edgeMap(1,iSide) )
               n2 = nodes( myQuadMesh % edgeMap(2,iSide) )
               
               x1 = myQuadMesh % nodes(n1) % x
               y1 = myQuadMesh % nodes(n1) % y
               x2 = myQuadMesh % nodes(n2) % x
               y2 = myQuadMesh % nodes(n2) % y
               
               DO iNode = 0, gPolyDeg

                  xc(iNode,1) = x1 + (x2-x1)*HALF*( s(iNode) + ONE )
                  xc(iNode,2) = y1 + (y2-y1)*HALF*( s(iNode) + ONE )
  
               ENDDO
               elBoundCurves(iSide) % x = xc

            ELSEIF( bFlags(iSide) == 1 )then ! this is a boundary edge

                ! Read in the parametric curve
                DO iNode = 0, gPolyDeg
                   READ( fUnit, * ) x, y, z
                   xc(iNode,1) = x
                   xc(iNode,2) = y
                ENDDO
                elBoundCurves(iSide) % x = xc

             ELSE
 
                PRINT*,' S/R ReadSpecMeshFile : Impossible element boundary flag '

             ENDIF

         ENDDO

         ! Initialize this element's geometry
         myQuadMesh % elements(iEl) % nodeIDs(1:4) = nodes
         myQuadMesh % elements(iEl) % elementID    = iEl
         CALL myQuadMesh % geom(iEl) % GenerateMesh( interp, elBoundCurves )
         CALL myQuadMesh % geom(iEl) % GenerateMetrics( interp )

         ! Read in and parse the edge names
         READ( fUnit, '(1x, A100)' )  edgeNames

         ! Parse the edge names into four edge names
         iSide = 1
         iC = 1
         jC = 1
         thisEdge = ' ' 
         DO while( iSide <= 4 )

            IF( edgeNames(iC:iC) == ' ' )then ! we have reached a blank space
     
               n1 = nodes( myQuadMesh % edgeMap(1,iSide) )
               n2 = nodes( myQuadMesh % edgeMap(2,iSide) )
               
               IF( TRIM(thisEdge) == 'DIRICHLET' )then
                  
                  sideFlags(iEl,iSide) = DIRICHLET
                  
               ELSEIF( TRIM(thisEdge) == 'DIRICHLET_INFLOW' )then
                  
                  sideFlags(iEl,iSide) = DIRICHLET_INFLOW

               ELSEIF( TRIM(thisEdge) == 'DIRICHLET_OUTFLOW' )then
                  
                  sideFlags(iEl,iSide) = DIRICHLET_OUTFLOW

               ELSEIF( TRIM(thisEdge) == 'ROBIN' )then

                  sideFlags(iEl,iSide) = ROBIN

               ELSEIF( TRIM(thisEdge) == 'ROBIN_FORCED' )then

                  sideFlags(iEl,iSide) = ROBIN_FORCED

               ELSEIF( TRIM(thisEdge) == 'HOMOGENEOUS_NEUMANN' )then

                  sideFlags(iEl,iSide) = HOMOGENEOUS_NEUMANN

               ELSEIF( TRIM(thisEdge) == 'NEUMANN_WALL' )then

                  sideFlags(iEl,iSide) = NEUMANN_WALL

               ELSEIF( TRIM(thisEdge) == 'NEUMANN' )then

                  sideFlags(iEl,iSide) = NEUMANN
 
               ELSEIF( TRIM(thisEdge) == 'NO_NORMAL_FLOW' )then
                  
                  sideFlags(iEl,iSide) = NO_NORMAL_FLOW

               ELSEIF( TRIM(thisEdge) == 'PRESCRIBED' )then
               
                  sideFlags(iEl,iSide) = PRESCRIBED

               ELSEIF( TRIM(thisEdge) == 'RADIATION' )then
                  
                  sideFlags(iEl,iSide) = RADIATION

               ENDIF

               ! Reset thisEdge
               thisEdge = ' '
               jC = 1 
               iC = iC + 1
               ! Increment the side
               iSide = iSide + 1

            ELSE
            
               
               thisEdge(jC:jC) = edgeNames(iC:iC)
               iC = iC + 1
               jC = jC + 1

            ENDIF

         ENDDO ! while (iSide <= 4 )


      ENDDO ! iEl, cycle over the elements

      CLOSE( fUnit )
      CALL myQuadMesh % ConstructEdges( )
      
      ! Set the secondary element on boundary edges to the boundary flag
      DO iEdge = 1, nEdges

         e1 = myQuadMesh % edges(iEdge) % elementIDs(1)
         s1 = myQuadMesh % edges(iEdge) % elementSides(1)
      
         IF( sideFlags(e1,s1) == DIRICHLET )then

            myQuadMesh % edges(iEdge) % elementIDs(2) = DIRICHLET
            n = myQuadMesh % edges(iEdge) % nodeIDs
            myQuadMesh % nodes( n(1) ) % nodeType = DIRICHLET
            myQuadMesh % nodes( n(2) ) % nodeType = DIRICHLET
   
         ELSEIF( sideFlags(e1,s1) == DIRICHLET_INFLOW )then

            myQuadMesh % edges(iEdge) % elementIDs(2) = DIRICHLET_INFLOW
            n = myQuadMesh % edges(iEdge) % nodeIDs
            myQuadMesh % nodes( n(1) ) % nodeType = DIRICHLET_INFLOW
            myQuadMesh % nodes( n(2) ) % nodeType = DIRICHLET_INFLOW
            
         ELSEIF( sideFlags(e1,s1) == DIRICHLET_OUTFLOW )then

            myQuadMesh % edges(iEdge) % elementIDs(2) = DIRICHLET_OUTFLOW
            n = myQuadMesh % edges(iEdge) % nodeIDs
            myQuadMesh % nodes( n(1) ) % nodeType = DIRICHLET_OUTFLOW
            myQuadMesh % nodes( n(2) ) % nodeType = DIRICHLET_OUTFLOW
            
         ELSEIF( sideFlags(e1,s1) == ROBIN )then

            myQuadMesh % edges(iEdge) % elementIDs(2) = ROBIN
            n = myQuadMesh % edges(iEdge) % nodeIDs
            myQuadMesh % nodes( n(1) ) % nodeType = ROBIN
            myQuadMesh % nodes( n(2) ) % nodeType = ROBIN

         ELSEIF( sideFlags(e1,s1) == ROBIN_FORCED )then

            myQuadMesh % edges(iEdge) % elementIDs(2) = ROBIN_FORCED
            n = myQuadMesh % edges(iEdge) % nodeIDs
            myQuadMesh % nodes( n(1) ) % nodeType = ROBIN_FORCED
            myQuadMesh % nodes( n(2) ) % nodeType = ROBIN_FORCED

         ELSEIF( sideFlags(e1,s1) == HOMOGENEOUS_NEUMANN )then

            myQuadMesh % edges(iEdge) % elementIDs(2) = HOMOGENEOUS_NEUMANN
            n = myQuadMesh % edges(iEdge) % nodeIDs
            myQuadMesh % nodes( n(1) ) % nodeType = HOMOGENEOUS_NEUMANN
            myQuadMesh % nodes( n(2) ) % nodeType = HOMOGENEOUS_NEUMANN

         ELSEIF( sideFlags(e1,s1) == NEUMANN_WALL )then

            myQuadMesh % edges(iEdge) % elementIDs(2) = NEUMANN_WALL
            n = myQuadMesh % edges(iEdge) % nodeIDs
            myQuadMesh % nodes( n(1) ) % nodeType = NEUMANN_WALL
            myQuadMesh % nodes( n(2) ) % nodeType = NEUMANN_WALL

         ELSEIF( sideFlags(e1,s1) == NEUMANN )then

            myQuadMesh % edges(iEdge) % elementIDs(2) = NEUMANN
            n = myQuadMesh % edges(iEdge) % nodeIDs
            myQuadMesh % nodes( n(1) ) % nodeType = NEUMANN
            myQuadMesh % nodes( n(2) ) % nodeType = NEUMANN

         ELSEIF( sideFlags(e1,s1) == NO_NORMAL_FLOW )then

            myQuadMesh % edges(iEdge) % elementIDs(2) = NO_NORMAL_FLOW
            n = myQuadMesh % edges(iEdge) % nodeIDs
            myQuadMesh % nodes( n(1) ) % nodeType = NO_NORMAL_FLOW
            myQuadMesh % nodes( n(2) ) % nodeType = NO_NORMAL_FLOW

         ELSEIF( sideFlags(e1,s1) == PRESCRIBED )then

            myQuadMesh % edges(iEdge) % elementIDs(2) = PRESCRIBED
            n = myQuadMesh % edges(iEdge) % nodeIDs
            myQuadMesh % nodes( n(1) ) % nodeType = PRESCRIBED
            myQuadMesh % nodes( n(2) ) % nodeType = PRESCRIBED
            
         ELSEIF( sideFlags(e1,s1) == RADIATION )then

            myQuadMesh % edges(iEdge) % elementIDs(2) = RADIATION
            n = myQuadMesh % edges(iEdge) % nodeIDs
            myQuadMesh % nodes( n(1) ) % nodeType = RADIATION
            myQuadMesh % nodes( n(2) ) % nodeType = RADIATION
         ENDIF
 
      ENDDO

      DO iEdge = 1, nEdges
            s2 = myQuadMesh % edges(iEdge) % elementSides(2)
            IF(s2 < 0)THEN
               myQuadMesh % edges(iEdge) % start = interp % N -1
               myQuadMesh % edges(iEdge) % inc   = -1
            ELSE
               myQuadMesh % edges(iEdge) % start = 1
               myQuadMesh % edges(iEdge) % inc   = 1
            ENDIF
            
      ENDDO

      ! Get the node to element connectivity
      CALL myQuadMesh % GetNodeToElementConnectivity( )
      CALL myQuadMesh % ConstructElementNeighbors( )

      ! Clear up memory
      DEALLOCATE( s, xc, sideFlags )

      DO iSide = 1, 4
         CALL elBoundCurves(iSide) % Trash( )
      ENDDO
  
 END SUBROUTINE ReadSpecMeshFile_QuadMesh
!
!> \addtogroup QuadMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R ReadPeaceMeshFile 
! 
!> \fn ReadPeaceMeshFile_QuadMesh 
!! Reads a PeaceMesh file and constructs the QuadMesh data structure.
!! 
!! In constrast to the ReadSpecMesh routine, an interpolant is not required to build the mesh.
!! The geometry information is written into the 
!! The "PeaceMesh" file format is described in the SELF Technical documentation. It was initially
!! written in order to handle domain decomposition for MPI parallelism. For those wondering 
!! "PeaceMesh" was initially meant to be "PieceMesh" but due to initially developing this file
!! format and the associated I/O routines late at night, this naming convention went unnoticed 
!! initially. For no other reasons than "I like the convention PeaceMesh" and "I appreciate the reminder
!! that we all make mistakes", this naming has been left alone.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadMesh)    :: this <BR>
!! <B>CHARACTER</B>         :: filename <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ReadPeaceMeshFile( filename ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myQuadMesh <td> QuadMesh <td>
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the PeacMesh file
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ReadPeaceMeshFile_QuadMesh( myQuadMesh, filename )

   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(out)  :: myQuadMesh
   CHARACTER(*), INTENT(in)        :: filename
   ! LOCAL
   INTEGER :: nNodes, nElems, nEdges, N
   INTEGER :: iEdge, iNode, iEl
   INTEGER :: fUnit, k, i, j


      !PRINT*, 'Mesh File : '//TRIM( filename )//'.pc.mesh'
      
      ! Get a new file unit
      OPEN( UNIT    = NEWUNIT(fUnit), &
            FILE    = TRIM( filename )//'.pc.mesh', &
            FORM    = 'UNFORMATTED',&
            STATUS  = 'OLD', &
            ACCESS  = 'DIRECT', &
            CONVERT = 'BIG_ENDIAN', &
            RECL    = SIZEOF(nNodes) ) ! How to do variable record length

      ! ---- Gather the number of nodes, number of elements, and number of edges ---- !
      k = 1
      READ( fUnit, rec=k )nNodes
      k = k+1
      READ( fUnit, rec=k )nElems
      k = k+1
      READ( fUnit, rec=k )nEdges
      k = k+1
      READ( fUnit, rec=k )N
      k = k+1

      PRINT*, 'nNodes    : ', nNodes
      PRINT*, 'nElems    : ', nElems
      PRINT*, 'nEdges    : ', nEdges
      PRINT*, 'N         : ', N
      ! ---- Initialize the quadrature mesh (empty) ---- !
      CALL myQuadMesh % Initialize( nNodes, nElems, nEdges, N ) 
      
      ! ---- Read in the element connectivity ---- !
      DO iEl = 1, nElems
         READ( fUnit, rec=k ) myQuadMesh % elements(iEl) % elementID 
         k = k+1
         DO i = 1, 4
            READ( fUnit, rec=k ) myQuadMesh % elements(iEl) % nodeIDs(i)
            k = k+1
         ENDDO
      ENDDO 
      
      ! ---- Read in the edge information ---- !

      DO iEdge = 1, nEdges
         READ( fUnit, rec=k ) myQuadMesh % edges(iEdge) % edgeID
         k = k+1
         READ( fUnit, rec=k ) myQuadMesh % edges(iEdge) % boundaryID
         k = k+1
         DO i = 1, 2
            READ( fUnit, rec=k ) myQuadMesh % edges(iEdge) % nodeIDs(i)
            k = k+1
            READ( fUnit, rec=k ) myQuadMesh % edges(iEdge) % elementIDs(i)
            k = k+1
            READ( fUnit, rec=k ) myQuadMesh % edges(iEdge) % elementSides(i)
            k = k+1
         ENDDO
         READ( fUnit, rec=k ) myQuadMesh % edges(iEdge) % start
         k = k+1
         READ( fUnit, rec=k ) myQuadMesh % edges(iEdge) % inc
         k = k+1
      ENDDO
      
      CLOSE( fUnit )
      ! Get a new file unit
      OPEN( UNIT    = NEWUNIT(fUnit), &
            FILE    = TRIM( filename )//'.pc.geom', &
            FORM    = 'UNFORMATTED',&
            STATUS  = 'OLD', &
            ACCESS  = 'DIRECT', &
            CONVERT = 'BIG_ENDIAN', &
            RECL    = prec ) ! How to do variable record length
      
      ! ---- Read in the corner nodes ---- !
      k = 1
      DO iNode = 1, nNodes  ! Loop over the nodes in the file
         READ( fUnit, rec=k ) myQuadMesh % nodes(iNode) % x
         k = k+1
         READ( fUnit, rec=k ) myQuadMesh % nodes(iNode) % y
         k = k+1
      ENDDO

      ! ---- Read in the element information ---- !
      DO iEl = 1, nElems
         DO j = 0, N
            DO i = 0, N
               READ( fUnit, rec=k ) myQuadMesh % geom(iEl) % x(i,j)
               k = k+1
               READ( fUnit, rec=k ) myQuadMesh % geom(iEl) % y(i,j)
               k = k+1
               READ( fUnit, rec=k ) myQuadMesh % geom(iEl) % dxds(i,j)
               k = k+1
               READ( fUnit, rec=k ) myQuadMesh % geom(iEl) % dxdp(i,j)
               k = k+1
               READ( fUnit, rec=k ) myQuadMesh % geom(iEl) % dyds(i,j)
               k = k+1
               READ( fUnit, rec=k ) myQuadMesh % geom(iEl) % dydp(i,j)
               k = k+1
               READ( fUnit, rec=k ) myQuadMesh % geom(iEl) % J(i,j)
               k = k+1
            ENDDO
         ENDDO
         DO j = 1, nQuadEdges
            DO i = 0, N
               READ( fUnit, rec=k ) myQuadMesh % geom(iEl) % xBound(i,j)
               k = k+1
               READ( fUnit, rec=k ) myQuadMesh % geom(iEl) % yBound(i,j)
               k = k+1
               READ( fUnit, rec=k ) myQuadMesh % geom(iEl) % nHat(1,i,j)
               k = k+1
               READ( fUnit, rec=k ) myQuadMesh % geom(iEl) % nHat(2,i,j)
               k = k+1
            ENDDO
         ENDDO
      ENDDO 

      CLOSE( fUnit )
     
      ! Get the node to element connectivity
      CALL myQuadMesh % GetNodeToElementConnectivity( )
      CALL myQuadMesh % ConstructElementNeighbors( )
  
 END SUBROUTINE ReadPeaceMeshFile_QuadMesh
!
!> \addtogroup QuadMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R WritePeaceMeshFile 
! 
!> \fn WritePeaceMeshFile_QuadMesh 
!! Writes a PeaceMesh file using the QuadMesh data structure.
!! 
!! The "PeaceMesh" file format is described in the SELF Technical documentation. It was initially
!! written in order to handle domain decomposition for MPI parallelism. For those wondering 
!! "PeaceMesh" was initially meant to be "PieceMesh" but due to initially developing this file
!! format and the associated I/O routines late at night, this naming convention went unnoticed 
!! initially. For no other reasons than "I like the convention PeaceMesh" and "I appreciate the reminder
!! that we all make mistakes", this naming has been left alone.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadMesh)    :: this <BR>
!! <B>CHARACTER</B>         :: filename <BR>
!!         .... <BR>
!!     <B>CALL</B> this % WritePeaceMeshFile( filename ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myQuadMesh <td> QuadMesh <td>
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the PeaceMesh file
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
SUBROUTINE WritePeaceMeshFile_QuadMesh( myQuadMesh, filename )

   IMPLICIT NONE
   CLASS( QuadMesh ), INTENT(in)  :: myQuadMesh
   CHARACTER(*), INTENT(in)        :: filename
   ! LOCAL
   INTEGER :: nNodes, nElems, nEdges, N
   INTEGER :: iEdge, iNode, iEl
   INTEGER :: fUnit, k, i, j

      !PRINT*, 'Mesh File : '//TRIM( filename )//'.pc.mesh'
      
      ! Get a new file unit
      OPEN( UNIT    = NEWUNIT(fUnit), &
            FILE    = TRIM( filename )//'.pc.mesh', &
            FORM    = 'UNFORMATTED',&
            STATUS  = 'REPLACE', &
            ACCESS  = 'DIRECT', &
            CONVERT = 'BIG_ENDIAN', &
            RECL    = SIZEOF(nNodes) ) ! How to do variable record length
      
      nNodes = myQuadMesh % nNodes
      nElems = myQuadMesh % nElems
      nEdges = myQuadMesh % nEdges
      N      = myQuadMesh % geom(1) % N
      ! ---- Gather the number of nodes, number of elements, and number of edges ---- !
      k = 1
      WRITE( fUnit, rec=k ) nNodes
      k = k+1
      WRITE( fUnit, rec=k ) nElems
      k = k+1
      WRITE( fUnit, rec=k ) nEdges
      k = k+1
      WRITE( fUnit, rec=k ) N
      k = k+1

      ! ---- WRITE in the element connectivity ---- !
      DO iEl = 1, nElems
         WRITE( fUnit, rec=k ) myQuadMesh % elements(iEl) % elementID 
         k = k+1
         DO i = 1, 4
            WRITE( fUnit, rec=k ) myQuadMesh % elements(iEl) % nodeIDs(i)
            k = k+1
         ENDDO
      ENDDO 
      
      ! ---- WRITE in the edge information ---- !

      DO iEdge = 1, nEdges
         WRITE( fUnit, rec=k ) myQuadMesh % edges(iEdge) % edgeID
         k = k+1
         WRITE( fUnit, rec=k ) myQuadMesh % edges(iEdge) % boundaryID
         k = k+1
         DO i = 1, 2
            WRITE( fUnit, rec=k ) myQuadMesh % edges(iEdge) % nodeIDs(i)
            k = k+1
            WRITE( fUnit, rec=k ) myQuadMesh % edges(iEdge) % elementIDs(i)
            k = k+1
            WRITE( fUnit, rec=k ) myQuadMesh % edges(iEdge) % elementSides(i)
            k = k+1
         ENDDO
         WRITE( fUnit, rec=k ) myQuadMesh % edges(iEdge) % start
         k = k+1
         WRITE( fUnit, rec=k ) myQuadMesh % edges(iEdge) % inc
         k = k+1
      ENDDO
      
      CLOSE( fUnit )
      ! Get a new file unit
      OPEN( UNIT    = NEWUNIT(fUnit), &
            FILE    = TRIM( filename )//'.pc.geom', &
            FORM    = 'UNFORMATTED',&
            STATUS  = 'REPLACE', &
            ACCESS  = 'DIRECT', &
            CONVERT = 'BIG_ENDIAN', &
            RECL    = prec ) ! How to do variable record length
      
      ! ---- WRITE in the corner nodes ---- !
      k = 1
      DO iNode = 1, nNodes  ! Loop over the nodes in the file
         WRITE( fUnit, rec=k ) myQuadMesh % nodes(iNode) % x
         k = k+1
         WRITE( fUnit, rec=k ) myQuadMesh % nodes(iNode) % y
         k = k+1
      ENDDO

      ! ---- WRITE in the element information ---- !
      DO iEl = 1, nElems
         DO j = 0, N
            DO i = 0, N
               WRITE( fUnit, rec=k ) myQuadMesh % geom(iEl) % x(i,j)
               k = k+1
               WRITE( fUnit, rec=k ) myQuadMesh % geom(iEl) % y(i,j)
               k = k+1
               WRITE( fUnit, rec=k ) myQuadMesh % geom(iEl) % dxds(i,j)
               k = k+1
               WRITE( fUnit, rec=k ) myQuadMesh % geom(iEl) % dxdp(i,j)
               k = k+1
               WRITE( fUnit, rec=k ) myQuadMesh % geom(iEl) % dyds(i,j)
               k = k+1
               WRITE( fUnit, rec=k ) myQuadMesh % geom(iEl) % dydp(i,j)
               k = k+1
               WRITE( fUnit, rec=k ) myQuadMesh % geom(iEl) % J(i,j)
               k = k+1
            ENDDO
         ENDDO
         DO j = 1, nQuadEdges
            DO i = 0, N
               WRITE( fUnit, rec=k ) myQuadMesh % geom(iEl) % xBound(i,j)
               k = k+1
               WRITE( fUnit, rec=k ) myQuadMesh % geom(iEl) % yBound(i,j)
               k = k+1
               WRITE( fUnit, rec=k ) myQuadMesh % geom(iEl) % nHat(1,i,j)
               k = k+1
               WRITE( fUnit, rec=k ) myQuadMesh % geom(iEl) % nHat(2,i,j)
               k = k+1
            ENDDO
         ENDDO
      ENDDO 

      CLOSE( fUnit )
     
 END SUBROUTINE WritePeaceMeshFile_QuadMesh
!
!> \addtogroup QuadMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R WriteTecplot
! 
!> \fn WriteTecplot_QuadMesh  
!! Writes a tecplot file of the mesh geometry. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadMesh) :: this <BR>
!! <B>CHARACTER</B>      :: filename <BR>
!!         .... <BR>
!!     ! To write a file with a specified name <BR>
!!     <B>CALL</B> this % WriteTecplot( filename ) <BR>
!!     ! Or, the file "mesh.tec" will be used with the following calling sequence <BR>
!!     <B>CALL</B> this % WriteTecplot( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myQuadMesh <td> QuadMesh <td>
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the tecplot file
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE WriteTecplot_Quadmesh( myQuadMesh, filename )

   IMPLICIT NONE
   CLASS(QuadMesh), INTENT(inout)     :: myQuadMesh
   CHARACTER(*), INTENT(in), OPTIONAL :: filename  
   ! Local
   INTEGER :: iS, iP, N, iEl,fUnit
   CHARACTER(7) :: zoneID

      N = myQuadMesh % geom(1) % N

      IF( PRESENT(filename) )THEN
         OPEN( UNIT=NEWUNIT(fUnit), &
               FILE= TRIM(filename)//'.tec', &
               FORM='formatted')
      ELSE
         OPEN( UNIT=NEWUNIT(fUnit), &
               FILE= 'mesh.tec', &
               FORM='formatted')
      ENDIF
    
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Jacobian", "dxds", "dxdp", "dyds", "dydp" '

      DO iEl = 1, myQuadMesh % nElems

         WRITE(zoneID,'(I7.7)') myQuadMesh % elements(iEl) % elementID
         WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',N+1,', J=', N+1,',F=POINT'

         DO iP = 0, N
            DO iS = 0, N
               WRITE(fUnit,*)  myQuadMesh % geom(iEl) % x(iS,iP), &
                               myQuadMesh % geom(iEl) % y(iS,iP), &
                               myQuadMesh % geom(iEl) % J(iS,iP), &
                               myQuadMesh % geom(iEl) % dxds(iS,iP), &
                               myQuadMesh % geom(iEl) % dxdp(iS,iP), &
                               myQuadMesh % geom(iEl) % dyds(iS,iP), &
                               myQuadMesh % geom(iEl) % dydp(iS,iP)
            ENDDO
         ENDDO
      ENDDO
    
      CLOSE(UNIT=fUnit)

    RETURN

 END SUBROUTINE WriteTecplot_Quadmesh
!
!> \addtogroup QuadMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R WriteMaterialTecplot
! 
!> \fn WriteMaterialTecplot_QuadMesh  
!! Writes a tecplot file of the mesh geometry and a given "material" field. 
!! 
!! When performing domain decomposition (e.g. with the DecomposeQuadMesh.f90 program), the 
!! decomposition can be visualized by assigning each element a number that corresponds to the 
!! process ID it has been assigned to. This routine takes in an array of "material ID's" that
!! identify each element.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(QuadMesh) :: this <BR>
!! <B>REAL</B>(prec)     :: pID(1:this % nElems) <BR>
!! <B>CHARACTER</B>      :: filename <BR>
!!         .... <BR>
!!     ! To write a file with a specified name <BR>
!!     <B>CALL</B> this % WriteMaterialTecplot( pID, filename ) <BR>
!!     ! Or, the file "mesh.tec" will be used with the following calling sequence <BR>
!!     <B>CALL</B> this % WriteMaterialTecplot( pID ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myQuadMesh <td> QuadMesh <td>
!!   <tr> <td> in <th> materialIDs(1:myQuadMesh % nElems) <td> REAL(prec) <td> 
!!                     Array of values that are assigned to each element.
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the MaterialTecplot file
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE WriteMaterialTecplot_Quadmesh( myQuadMesh, materialIDs, filename )

  IMPLICIT NONE
  CLASS(QuadMesh), INTENT(inout)     :: myQuadMesh
  REAL(prec), INTENT(in)             :: materialIDs(1:myQuadMesh % nElems)
  CHARACTER(*), INTENT(in), OPTIONAL :: filename  
  ! Local
  INTEGER :: iS, iP, N, iEl,fUnit
  CHARACTER(7) :: zoneID

     N = myQuadMesh % geom(1) % N

    IF( PRESENT(filename) )THEN
    
       OPEN( UNIT=NEWUNIT(fUnit), &
             FILE= TRIM(filename)//'.tec', &
             FORM='formatted')
    ELSE
    
       OPEN( UNIT=NEWUNIT(fUnit), &
             FILE= 'mesh.tec', &
             FORM='formatted')

    ENDIF
    
    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "materialIDs" '


    DO iEl = 1, myQuadMesh % nElems

       WRITE(zoneID,'(I7.7)') myQuadMesh % elements(iEl) % elementID
       WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',N+1,', J=', N+1,',F=POINT'

       DO iP = 0, N
          DO iS = 0, N
             WRITE(fUnit,*)  myQuadMesh % geom(iEl) % x(iS,iP), &
                             myQuadMesh % geom(iEl) % y(iS,iP), &
                             materialIDs(iEl)
          ENDDO
      ENDDO

    ENDDO
    
    CLOSE(UNIT=fUnit)

    RETURN

 END SUBROUTINE WriteMaterialTecplot_Quadmesh
!
END MODULE QuadMesh_Class
