
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.

PROGRAM DecomposeHexMesh

 ! src/common/
 USE ModelPrecision
 USE ConstantsDictionary
 USE CommonRoutines
 USE Timing
 ! src/nodal/
 USE NodalStorage_Class
 ! src/geom/
 USE HexMesh_Class
 USE TopographicShapes
 USE Params_Class
 ! src/boundary/
 USE BoundaryCommunicator_Class
 
 ! To make use of the METIS data types and intrinsics
  
 IMPLICIT NONE

 ! ===== METIS Variables ===== !
 INTEGER(KIND=8)              :: ncon, objval, nparts, nvtxs
 INTEGER(KIND=8), POINTER     :: vwgt=>null(), vsize=>null(), adjwgt=>null()
 INTEGER(KIND=8)              :: metisOptions(0:40)
 REAL(prec), POINTER          :: tpwgts=>null(), ubvec=>null()
 INTEGER(KIND=8), ALLOCATABLE :: xadj(:), adjacency(:) , partitions(:)
 ! =========================== !

 TYPE( NodalStorage )         :: nodal
 TYPE( HexMesh )              :: mesh, procMesh
 TYPE( RunParams )            :: params
 TYPE( MultiTimers )          :: timers
 TYPE( BoundaryCommunicator ) :: bcom
 INTEGER(KIND=8)              :: procID
 INTEGER                      :: pfnodes(1:nQuadNodes), gfNodes(1:nQuadNodes), nIDs(1:nQuadNodes)
 INTEGER                      :: nElems, nProc, pID, i
 INTEGER                      :: iEl, jEl, iSide, iNode, nAdj, nID, elID, nBe, nMPI
 INTEGER                      :: e2, iFace, jFace, localID
 INTEGER, ALLOCATABLE         :: globalToLocal(:,:), nElPerProc(:), nodeLogic(:,:), nNodePerProc(:)
 INTEGER, ALLOCATABLE         :: globalToLocalNode(:,:)
 REAL(prec), ALLOCATABLE      :: materials(:)
 CHARACTER(4)                 :: pIDChar
 
      ! Read in the parameters
      CALL params % Build( )
      CALL timers % Build( )
      
      CALL timers % AddTimer( 'Spectral Partition Time', 1 )
      ! Build an interpolant
      CALL nodal % Build( N = params % polyDeg, &
                          nPlot = params % nPlot, &
                          quadrature = GAUSS,  &
                          approxForm = DG )
                          
      IF( params % topographicShape == Gaussian )THEN
         TopographicShape => GaussianHill
      ELSE
         TopographicShape => DefaultTopography
      ENDIF
      
      ! Build the Geometry
      IF( TRIM( params % UCDMeshFile ) == nada )THEN
         IF( params % MeshType == DoublyPeriodic )THEN
            PRINT*,' Loading doubly periodic mesh.'
            CALL mesh % LoadDoublyPeriodicMesh( nodal % interp, &
                                                params % nXelem, &
                                                params % nYelem, &
                                                params % nZelem )
         ELSE
            PRINT*,' Loading default mesh.'
            CALL mesh % LoadDefaultMesh( nodal % interp, &
                                         params % nXelem, &
                                         params % nYelem, &
                                         params % nZelem )
         ENDIF
      ELSE
      ! Builds the lateral mesh
         PRINT*, 'Reading mesh from '//trim(params % UCDMeshFile)//'.'
         CALL  mesh % ReadUCDMeshFile( nodal % interp, &
                                       params % UCDMeshFile )
      ENDIF

      nElems = mesh % nElems
      nProc  = params % nProc
      ALLOCATE( materials(1:nElems), &
                globalToLocal(1:nElems,1:2), &
                nElPerProc(0:nProc-1), &
                nodeLogic(1:mesh % nNodes, 0:nProc-1), &
                nNodePerProc(0:nProc-1), &
                globalToLocalNode(1:mesh % nNodes, 0:nProc-1) )
                
      materials     = ZERO
      globalToLocal = 0
      nElPerProc    = 0
      nodeLogic     = 0
      nNodePerProc  = 0
      
      !========= Assign METIS parameters ===========!
      ncon   = 1
      nvtxs  = nElems
      nparts = nProc
      ! ============================================ !
      ! To allocate space for the adjacency, excess space is allocated
      ALLOCATE( xadj(0:nElems), partitions(0:nElems-1) )
      xadj = 0
      partitions = 0
      
      ! Determine the number of faces
      nAdj = 0
      xadj(0) = 0
      DO iEl = 1, nElems
         DO iSide = 1, 6
            jEl = mesh % elements(iEl) % neighbors(iSide)
            IF( jEl > 0 )THEN
               ! Increment the number of adjacent elements
               nAdj = nAdj + 1
            ENDIF
         ENDDO
         xadj(iEl) = nAdj
      ENDDO
      ! To allocate space for the adjacency, excess space is allocated
      ALLOCATE( adjacency(0:nAdj-1) )
      adjacency = 0
      
      CALL timers % StartTimer( 1 )
      
      nAdj = 0
      ! Now we set up the adjacency graph
      DO iEl = 1, nElems
         DO iSide = 1, 6
            jEl = mesh % elements(iEl) % neighbors(iSide)
            IF( jEl > 0 )THEN
               ! Add the element ID to the adjacency list
               adjacency(nAdj) = jEl-1
               ! Increment the number of adjacent elements
               nAdj = nAdj + 1
            ENDIF
         ENDDO
      ENDDO
      
     
      ! Now call METIS
      CALL METIS_SetDefaultOptions(metisOptions)
      metisoptions(17) = 0 ! METIS_OPTION_NUMBERING
      
      ! The the number of processes is greater than one, then we need to call then
      ! K-way graph partitioning.
      IF ( nProc > 1 )THEN
         PRINT*,' DecomposeHexMesh : Calling METIS_PartGraphKWay.'
     
         CALL METIS_PartGraphKWay( nvtxs, ncon, xadj, adjacency, &
                                   vwgt, vsize, adjwgt, nparts, tpwgts, &
                                   ubvec, metisoptions, objval, partitions )
                            
      
         PRINT*,' DecomposeHexMesh :  Building Connectivity.'
         ! Now, we need to build each mesh and and the associated connecivity
         DO iEl = 1, nElems
      
            procID               = partitions(iEl-1)
            nElPerProc(procID)   = nElPerProc(procID) + 1 ! Increment the local element ID for this process
            globalToLocal(iEl,1) = nElPerProc(procID)     ! Store the local element ID on this process
            globalToLocal(iEl,2) = procID                 ! Store the process ID
         
            DO iNode = 1, 8
               nID = mesh % elements(iEl) % nodeIDs(iNode) 
               nodeLogic(nID,procID) = 1
            ENDDO
         
         ENDDO
         DO nID = 1, mesh % nNodes
            nNodePerProc = nNodePerProc + nodeLogic(nID,:)
            globalToLocalNode(nID,:) = nNodePerProc
         ENDDO
      
         DO procID = 0, nProc-1
            PRINT*, 'Process ID :',procID, ', nElems :', nElPerProc(procID)
            PRINT*, 'Process ID :',procID, ', nNodes :', nNodePerProc(procID)
         ENDDO
      
         CALL timers % StopTimer( 1 )
      
      ELSE  ! This is a single rank job
         DO iEl = 1, nElems
            nElPerProc(0) = mesh % nElems
            globalToLocal(iEl,1) = iEl ! Local ID is global ID
            globalToLocal(iEl,2) = 0   ! process ID
            
            DO iNode = 1, 8
               nID = mesh % elements(iEl) % nodeIDs(iNode) 
               nodeLogic(nID,0) = 1
            ENDDO
         ENDDO
         
         DO nID = 1, mesh % nNodes
            nNodePerProc = nNodePerProc + nodeLogic(nID,:)
            globalToLocalNode(nID,:) = nNodePerProc
         ENDDO
         
      ENDIF
         
         ! Add in the node logic!!!
         
      ! Now we generate the local mesh for each process
      DO procID = 0, nProc-1
      
         CALL procMesh % Initialize(  nNodePerProc(procID), &
                                      nElPerProc(procID), &
                                      1, params % polyDeg )
                                      
         !PRINT*, 'Setting up local nodes'
         ! The nodes on the local node mesh are assigned using the global mesh and the "nodeLogic"
         ! array
         DO nID = 1, mesh % nNodes
         
            ! We only assign node information if it is "owned" by this process
            IF( nodeLogic(nID,procID) == 1 )THEN 
               ! Obtain the local node ID from the "global-to-local" array
               localID = globalToLocalNode( nID, procID )
               ! Assign the local node attributes
               procMesh % nodes(localID) % nodeID   = nID
               procMesh % nodes(localID) % nodeType = mesh % nodes(nID) % nodeType
               procMesh % nodes(localID) % x        = mesh % nodes(nID) % x
               procMesh % nodes(localID) % y        = mesh % nodes(nID) % y
               procMesh % nodes(localID) % z        = mesh % nodes(nID) % z
            ENDIF
            
         ENDDO

         ! Set the element-node connectivity
         DO iEl = 1, mesh % nElems
            ! We only assign element information if it is "owned" by this process
            IF( globalToLocal(iEl,2) == procID )THEN
               elID = globalToLocal(iEl,1) ! Obtain the local element ID from the global-to-local array
               
               DO iNode = 1, 8
                  ! First we grab the global node ID for this element's corner node
                  nID = mesh % elements(iEl) % nodeIDs(iNode)
                  ! Then we convert this to a local node ID for this process
                  localID = globalToLocalNode( nID, procID )
                  ! And now we assign the local Node ID to the local element's corner node ID's
                  procMesh % elements(elID) % nodeIDs(iNode)  = localID
               ENDDO
               procMesh % elements(elID) % elementID = iEl
                  
               ! Set the geometry
               procMesh % geom(elID) % nHat   = mesh % geom(iEl) % nHat
               procMesh % geom(elID) % xBound = mesh % geom(iEl) % xBound
               procMesh % geom(elID) % yBound = mesh % geom(iEl) % yBound
               procMesh % geom(elID) % zBound = mesh % geom(iEl) % zBound
               procMesh % geom(elID) % x      = mesh % geom(iEl) % x
               procMesh % geom(elID) % y      = mesh % geom(iEl) % y
               procMesh % geom(elID) % z      = mesh % geom(iEl) % z
               procMesh % geom(elID) % J      = mesh % geom(iEl) % J
               procMesh % geom(elID) % dxds   = mesh % geom(iEl) % dxds
               procMesh % geom(elID) % dxdp   = mesh % geom(iEl) % dxdp
               procMesh % geom(elID) % dxdq   = mesh % geom(iEl) % dxdq
               procMesh % geom(elID) % dyds   = mesh % geom(iEl) % dyds
               procMesh % geom(elID) % dydp   = mesh % geom(iEl) % dydp
               procMesh % geom(elID) % dydq   = mesh % geom(iEl) % dydq
               procMesh % geom(elID) % dzds   = mesh % geom(iEl) % dzds
               procMesh % geom(elID) % dzdp   = mesh % geom(iEl) % dzdp
               procMesh % geom(elID) % dzdq   = mesh % geom(iEl) % dzdq
               procMesh % geom(elID) % Ja     = mesh % geom(iEl) % Ja

            ENDIF
            
         ENDDO
         
         ! Set the faces
         ! First we count the local faces and fill in the face connectivity
         
         CALL procMesh % ConstructFaces( )
         PRINT*, 'Process ID :',procID, ', nFaces :', procMesh % nFaces
         
         
         ! Now we need to cycle through the faces and find the boundary faces
         ! If we find a boundary face (unassigned secondary element), then we need to search
         ! through the global mesh to find the matching face to determine the appropriate 
         ! boundary condition. The boundary condition could either be a physical boundary condition,
         ! numerical boundary condition, or may need to indicate an MPI_SEND/RECV pair.
         
         ! Count the number of boundary faces
         nBe = 0
         DO iFace = 1, procMesh % nFaces
            e2 = procMesh % faces(iFace) % ElementIDs(2)
            IF( e2 <= 0 )THEN
               nBe = nBe + 1
            ENDIF
         ENDDO
         PRINT*, 'Process ID :',procID, ', nBFace :', nBe
         !PRINT*, 'Setting up boundary communications'
         
         ! And initialize the boundary communicator
         CALL bCom % Initialize( nBe )
         ! Assign the boundary face ID to the
         nBe = 0
         nMPI = 0
         DO iFace = 1, procMesh % nFaces
         
            e2 = procMesh % faces(iFace) % ElementIDs(2)
            IF( e2 <= 0 )THEN ! This face is a boundary of this process's face
            
               nBe = nBe + 1
               bCom % boundaryIDs(nBe) = iFace
               ! Initially fill in the boundary condition as provided, and assign the external
               ! process ID to the current process ID (implying no MPI communication)
               bCom % extElemIDs(nBe) = e2
               bCom % extProcIDs(nBe) = procID 
               procMesh % faces(iFace) % boundaryID = nBe
               
               ! Now we search through the global mesh for an identical face
               ! Recall that the face is identified by its terminating nodes.
               ! Here, we must use the global node ID's
               nIDs = procMesh % faces(iFace) % nodeIDs ! local node ID
               DO i = 1, nQuadNodes
                  pfnodes(i) = procMesh % nodes(nIDs(i)) % nodeID
               ENDDO
               DO jFace = 1, mesh % nFaces ! Cycle through the global mesh and search for an identical face
                  ! Obtain the terminating global node IDs from the global mesh
                  gfnodes = mesh % faces(jFace) % nodeIDs
              
                  IF( NodesAreTheSame( pfnodes, gfnodes, nQuadNodes ) )THEN ! we have found an identical face
                     
                     ! In this case, we have found an identical face 
                     ! We need to copy over the secondary element information
                     iEl = mesh % faces(jFace) % elementIDs(2)
                     IF( iEl > 0 )THEN ! If iEl > 0, then this face is internal to the global mesh
                                       ! implying that an MPI communication is needed
                        nMPI = nMPI + 1
                        ! For this case, an MPI_SEND/RECV pair is needed and all secondary element
                        ! information must be copied from the global mesh
                        elID = globalToLocal(iEl, 1) ! obtain the process local element ID
                        pID  = globalToLocal(iEl, 2) ! obtain the process ID
                        
                        ! If the process ID of the global-secondary element is owned by the current 
                        ! process, then we need to assign the global primary element as the local
                        ! secondary element
                        IF( pID == procID )THEN
                           iEl = mesh % faces(jFace) % elementIDs(1)
                           elID = globalToLocal(iEl, 1) ! obtain the process local element ID
                           pID  = globalToLocal(iEl, 2) ! obtain the process ID
                        ENDIF
                        
                     !   print*,iEl, elID, pID
                        procMesh % faces(iFace) % faceID = jFace
                        procMesh % faces(iFace) % elementIDs(2) = -elID
                        procMesh % faces(iFace) % iStart = mesh % faces(jFace) % iStart
                        procMesh % faces(iFace) % iInc = mesh % faces(jFace) % iInc
                        procMesh % faces(iFace) % jStart = mesh % faces(jFace) % jStart
                        procMesh % faces(iFace) % jInc = mesh % faces(jFace) % jInc
                        procMesh % faces(iFace) % swapDimensions = mesh % faces(jFace) % swapDimensions
                        
                        ! Update the boundary communicator information
                        bCom % extElemIDs(nBe) = elID
                        bCom % extProcIDs(nBe) = pID
                        
                     ELSE
                        ! In this case, a physical or numerical boundary is found and we only 
                        ! need to assign the secondary element ID
                        procMesh % faces(iFace) % elementIDs(2) = iEl
                     ENDIF
                     
                     ! In either case, we have identified the face and we can exit this inner loop
                     EXIT
                     
                  ENDIF
               
               ENDDO 
               
            ENDIF ! If this is a boundary face
            
         ENDDO ! loop over the local faces
         PRINT*, 'Process ID :',procID, ', nMPI   :', nMPI
         
         ! Now we need to write a peace-mesh file and and communicator file
         WRITE( pIDChar, '(I4.4)' ) procID
         CALL procmesh % WriteTecplot( 'mesh.'//pIDChar )
         CALL procmesh % WritePeaceMeshFile( TRIM(params % PeaceMeshFile)//'.'//pIDChar )
         CALL bCom % WritePickup( 'ExtComm.'//pIDChar )
         
         CALL procMesh % Trash( )
         CALL bCom % Trash( )
         
      ENDDO
      
      ! For visualizing the decomposition
      materials = REAL( partitions, prec )
      CALL mesh % WriteMaterialTecplot( materials )
      
      CALL timers % Write_MultiTimers( )

      ! Clean up memory !
      DEALLOCATE( materials )
      
   !   CALL mesh % WriteTecplot( )
      CALL mesh % Trash( )
      CALL nodal % Trash( )
      CALL timers % Trash( )

 CONTAINS

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

END PROGRAM DecomposeHexMesh
