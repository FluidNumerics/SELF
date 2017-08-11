PROGRAM DecomposeStructuredHexMesh

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
 
  
 IMPLICIT NONE

 TYPE( NodalStorage )                      :: nodal
 TYPE( HexMesh )                           :: mesh,
 TYPE( HexMesh ), ALLOCATABLE              :: procMesh(:)
 TYPE( RunParams )                         :: params
 TYPE( MultiTimers )                       :: timers
 TYPE( BoundaryCommunicator ), ALLOCATABLE :: bcom(:)
 INTEGER, ALLOCATABLE                      :: faceProcCount(:), faceProcOwners(:,:), faceBoundaryIDs(:,:)
 INTEGER(KIND=8)               :: procID
 INTEGER                      :: pfnodes(1:nQuadNodes), gfNodes(1:nQuadNodes), nIDs(1:nQuadNodes)
 INTEGER                      :: nElems, nProc, pID, i
 INTEGER                      :: iEl, jEl, iSide, iNode, nAdj, nID, elID, nBe, nMPI
 INTEGER                      :: e2, iFace, jFace, localID, eIDs(1:2), pIDs(1:2)
 INTEGER                      :: nxp, nyp, nzp, iPx, iPy, iPz, iXp, iYp, iZp, iX, iY, iZ
 INTEGER, ALLOCATABLE         :: globalToLocal(:,:), nElPerProc(:), nodeLogic(:,:), nNodePerProc(:), partitions(:)
 INTEGER, ALLOCATABLE         :: globalToLocalNode(:,:)
 REAL(prec), ALLOCATABLE      :: materials(:)
 

      CALL Setup( )

      CALL PartitionElementsAndNodes_Structured( )
      
         
      ! Now we generate the local mesh for each process
      DO procID = 0, nProc-1
      
         CALL procMesh(procID) % Initialize(  nNodePerProc(procID), &
                                      nElPerProc(procID), &
                                      1, params % polyDeg )
          
         CALL DecomposeNodes( procID )

         CALL DecomposeElements( procID )                            

         CALL DecomposeFaces( procID ) 
         
         CALL FileIO( procID ) 
         
      ENDDO

      CALL SetupCommTables( )
      
      CALL Cleanup( )

 CONTAINS

 SUBROUTINE Setup( )
    IMPLICIT NONE

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

      nElems = mesh % nElems
      nProc  = params % nProcX*params % nProcY*params % nProcZ
      PRINT*, nProc
      IF( nProc == 0 )THEN
         nProc = 1
      ENDIF
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
      
      ALLOCATE( partitions(1:nElems) )
      partitions = 0
      ALLOCATE( bcom(0:nProcs-1), procMesh(0:nProcs-1) ) 
      ALLOCATE( faceProcCount(1:mesh % nFaces), &
                faceProcOwners(1:mesh % nFaces,1:2), &
                faceBoundaryIDs(1:mesh % nFaces,1:2) )

      faceProcCount   = 0
      faceProcOwners  = -1
      faceBoundaryIDs = 0

 END SUBROUTINE Setup
!
 SUBROUTINE Cleanup( )
    IMPLICIT NONE
    
      DO procID = 0, nProc-1
         CALL bCom(procID) % Trash( )
         CALL procMesh(procID) % Trash( )
      ENDDO
      DEALLOCATE( bCom, procMesh )
      

      ! For visualizing the decomposition
      materials = REAL( partitions, prec )
      CALL mesh % WriteMaterialTecplot( materials )
      
      CALL timers % Write_MultiTimers( )

      ! Clean up memory !
      DEALLOCATE( materials, &
                  faceProcCount, &
                  faceProcOwners, &
                  faceBoundaryIDs )
      
      CALL mesh % Trash( )
      CALL nodal % Trash( )
      CALL timers % Trash( )

 END SUBROUTINE Cleanup
!
 SUBROUTINE FileIO( procID )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: procID
    ! Local 
    CHARACTER(4) :: pIDChar


         ! Now we need to write a peace-mesh file and and communicator file
         WRITE( pIDChar, '(I4.4)' ) procID
         CALL procmesh(procID) % WriteTecplot( 'mesh.'//pIDChar )
         CALL procmesh(procID) % WritePeaceMeshFile( TRIM(params % PeaceMeshFile)//'.'//pIDChar )
         !CALL bCom(procID) % WritePickup( 'ExtComm.'//pIDChar )
         
         CALL procMesh(procID) % Trash( )
 
 END SUBROUTINE FileIO
!
 SUBROUTINE PartionElementsAndNodes_Structured( )
    IMPLICIT NONE
    

      IF( nProc > 1 )THEN
         nxp = params % nXelem/params % nProcX
         nyp = params % nYelem/params % nProcY
         nzp = params % nZelem/params % nProcZ
         DO iPz = 1, params % nProcZ
            DO iPy = 1, params % nProcY
               DO iPx = 1, params % nProcX
               
                  DO iZp = 1, nzp
                     DO iYp = 1, nyp
                        DO iXp = 1, nxp
                        
                           iX = iXp + (iPx-1)*nxp
                           iY = iYp + (iPy-1)*nyp
                           iZ = iZp + (iPz-1)*nzp
                           
                           iEl = iX + params % nXelem*( iY-1 + params % nYelem*( iZ-1 ) )
                           partitions(iEl) = iPx-1 + params % nProcX*( iPy-1 + params % nProcY*(iPz-1) )
                           
                        ENDDO
                     ENDDO
                  ENDDO
                  
               ENDDO
            ENDDO
         ENDDO
      
         PRINT*,' DecomposeHexMesh :  Building Connectivity.'
         ! Now, we need to build each mesh and and the associated connecivity
         DO iEl = 1, nElems
      
            procID               = partitions(iEl)
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

      ELSE

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
 
      
 END SUBROUTINE PartitionElementsAndNodes_Structured
!
 SUBROUTINE DecomposeNodes( procID )
   IMPLICIT NONE
   INTEGER, INTENT(in) :: procID
   ! Local 
   INTEGER :: nID, localID

     DO nID = 1, mesh % nNodes
     
        ! We only assign node information if it is "owned" by this process
        IF( nodeLogic(nID,procID) == 1 )THEN 
           ! Obtain the local node ID from the "global-to-local" array
           localID = globalToLocalNode( nID, procID )
           ! Assign the local node attributes
           procMesh(procID) % nodes(localID) % nodeID   = nID
           procMesh(procID) % nodes(localID) % nodeType = mesh % nodes(nID) % nodeType
           procMesh(procID) % nodes(localID) % x        = mesh % nodes(nID) % x
           procMesh(procID) % nodes(localID) % y        = mesh % nodes(nID) % y
           procMesh(procID) % nodes(localID) % z        = mesh % nodes(nID) % z
        ENDIF
        
     ENDDO

 END SUBROUTINE DecomposeNodes
!
 SUBROUTINE DecomposeElements( procID )
   IMPLICIT NONE
   INTEGER, INTENT(in) :: procID
   ! Local
   INTEGER :: iEl, elID, nID, localID

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
              procMesh(procID) % elements(elID) % nodeIDs(iNode)  = localID
           ENDDO
           procMesh(procID) % elements(elID) % elementID = iEl
              
           ! Set the geometry
           procMesh(procID) % geom(elID) % nHat   = mesh % geom(iEl) % nHat
           procMesh(procID) % geom(elID) % xBound = mesh % geom(iEl) % xBound
           procMesh(procID) % geom(elID) % yBound = mesh % geom(iEl) % yBound
           procMesh(procID) % geom(elID) % zBound = mesh % geom(iEl) % zBound
           procMesh(procID) % geom(elID) % x      = mesh % geom(iEl) % x
           procMesh(procID) % geom(elID) % y      = mesh % geom(iEl) % y
           procMesh(procID) % geom(elID) % z      = mesh % geom(iEl) % z
           procMesh(procID) % geom(elID) % J      = mesh % geom(iEl) % J
           procMesh(procID) % geom(elID) % dxds   = mesh % geom(iEl) % dxds
           procMesh(procID) % geom(elID) % dxdp   = mesh % geom(iEl) % dxdp
           procMesh(procID) % geom(elID) % dxdq   = mesh % geom(iEl) % dxdq
           procMesh(procID) % geom(elID) % dyds   = mesh % geom(iEl) % dyds
           procMesh(procID) % geom(elID) % dydp   = mesh % geom(iEl) % dydp
           procMesh(procID) % geom(elID) % dydq   = mesh % geom(iEl) % dydq
           procMesh(procID) % geom(elID) % dzds   = mesh % geom(iEl) % dzds
           procMesh(procID) % geom(elID) % dzdp   = mesh % geom(iEl) % dzdp
           procMesh(procID) % geom(elID) % dzdq   = mesh % geom(iEl) % dzdq
           procMesh(procID) % geom(elID) % Ja     = mesh % geom(iEl) % Ja

        ENDIF
        
     ENDDO

 END SUBROUTINE DecomposeElements
!
 SUBROUTINE DecomposeFaces( procID )
   IMPLICIT NONE
   INTEGER, INTENT(in) :: procID 
   INTEGER             :: e1, e2, p1, p2, iFace, iFaceLocal
   INTEGER             :: nBe, nMPI

      procMesh(procID) % nFaces = 0 
      nBe               = 0
      DO iFace = 1, mesh % nFaces

         e1 = mesh % faces(iFace) % elementIDs(1)
         e2 = mesh % faces(iFace) % elementIDs(2)

         p1 = globalToLocal(e1,2)
         IF( e2 > 0 )THEN
            p2 = globalToLocal(e2,2)
         ELSE
            p2  = p1
         ENDIF

         IF( p1 == procID .OR. p2 == procID )THEN
            procMesh(procID) % nFaces = procMesh(procID) % nFaces + 1
            IF( e2 < 0 .OR. p2 /= p1 )THEN
               nBe = nBe + 1
            ENDIF
         ENDIF

      ENDDO

      PRINT*, 'Process ID :',procID, ', nFaces :', procMesh(procID) % nFaces
      PRINT*, 'Process ID :',procID, ', nBFace :', nBe
      DEALLOCATE( procMesh(procID) % faces )
      ALLOCATE( procMesh(procID) % faces(1:proceMesh % nFaces) ) 
      CALL bCom(procID) % Initialize( nBe )

      iFaceLocal = 0
      nBe        = 0
      nMPI       = 0
      DO iFace = 1, mesh % nFaces

         e1 = mesh % faces(iFace) % elementIDs(1)
         e2 = mesh % faces(iFace) % elementIDs(2)

         p1 = globalToLocal(e1,2)
         IF( e2 > 0 )THEN
            p2 = globalToLocal(e2,2)
         ELSE
            p2  = p1
         ENDIF

         IF( p1 == procID .OR. p2 == procID )THEN
            iFaceLocal = iFaceLocal + 1
            faceProcCount(iFace) = faceProcCount(iFace) + 1
            IF( faceProcCount(iFace) > 2 )THEN
               PRINT*, ' DecomposeFaces : Something catastrophic happened !' 
               CALL Cleanup( )
               STOP 
            ENDIF
            faceProcOwners(iFace, faceProcCount(iFace)) = procID

            IF( p2 == p1 .AND. e2 > 0 )THEN ! Internal Face
              
               procMesh(procID) % faces(iFaceLocal) % faceID          = iFace
               procMesh(procID) % faces(iFaceLocal) % elementIDs(1)   = globalToLocal( e1, 1 ) 
               procMesh(procID) % faces(iFaceLocal) % elementIDs(2)   = globalToLocal( e2, 1 ) 
               procMesh(procID) % faces(iFaceLocal) % elementSides(1) = mesh % faces(iFace) % elementSides(1) 
               procMesh(procID) % faces(iFaceLocal) % elementSides(2) = mesh % faces(iFace) % elementSides(2)
               procMesh(procID) % faces(iFaceLocal) % iStart          = mesh % faces(iFace) % iStart
               procMesh(procID) % faces(iFaceLocal) % iInc            = mesh % faces(iFace) % iInc
               procMesh(procID) % faces(iFaceLocal) % jStart          = mesh % faces(iFace) % jStart
               procMesh(procID) % faces(iFaceLocal) % jInc            = mesh % faces(iFace) % jInc
               procMesh(procID) % faces(iFaceLocal) % swapDimensions  = mesh % faces(iFace) % swapDimensions

            ELSEIF( p2 == p1 .AND. e2 < 0 )THEN ! Physical boundary

               procMesh(procID) % faces(iFaceLocal) % faceID          = iFace
               procMesh(procID) % faces(iFaceLocal) % elementIDs(1)   = globalToLocal( e1, 1 ) 
               procMesh(procID) % faces(iFaceLocal) % elementIDs(2)   = e2 
               procMesh(procID) % faces(iFaceLocal) % elementSides(1) = mesh % faces(iFace) % elementSides(1) 
               procMesh(procID) % faces(iFaceLocal) % elementSides(2) = mesh % faces(iFace) % elementSides(2)
               procMesh(procID) % faces(iFaceLocal) % iStart          = mesh % faces(iFace) % iStart
               procMesh(procID) % faces(iFaceLocal) % iInc            = mesh % faces(iFace) % iInc
               procMesh(procID) % faces(iFaceLocal) % jStart          = mesh % faces(iFace) % jStart
               procMesh(procID) % faces(iFaceLocal) % jInc            = mesh % faces(iFace) % jInc
               procMesh(procID) % faces(iFaceLocal) % swapDimensions  = mesh % faces(iFace) % swapDimensions

               nBe = nBe + 1
               bCom(procID) % boundaryIDs(nBe) = iFaceLocal
               bCom(procID) % extProcIDs(nBe)  = procID 

               ! Physical boundary ID's will be set to negative for physical boundaries
               ! so that boundaries requiring communication can be separated from physical
               ! boundaries. This will allow more operations to be run concurrently with
               ! communication.
               procMesh(procID) % faces(iFaceLocal) % boundaryID = -nBe

            ELSEIF( p2 /= p1 .AND. p1 == procID )THEN ! MPI Boundary
 
               nMPI = nMPI + 1 
          
               procMesh(procID) % faces(iFaceLocal) % faceID          = iFace
               procMesh(procID) % faces(iFaceLocal) % elementIDs(1)   = globalToLocal( e1, 1 ) 
               procMesh(procID) % faces(iFaceLocal) % elementIDs(2)   = -e2 
               procMesh(procID) % faces(iFaceLocal) % elementSides(1) = mesh % faces(iFace) % elementSides(1) 
               procMesh(procID) % faces(iFaceLocal) % elementSides(2) = mesh % faces(iFace) % elementSides(2)
               procMesh(procID) % faces(iFaceLocal) % iStart          = mesh % faces(iFace) % iStart
               procMesh(procID) % faces(iFaceLocal) % iInc            = mesh % faces(iFace) % iInc
               procMesh(procID) % faces(iFaceLocal) % jStart          = mesh % faces(iFace) % jStart
               procMesh(procID) % faces(iFaceLocal) % jInc            = mesh % faces(iFace) % jInc
               procMesh(procID) % faces(iFaceLocal) % swapDimensions  = mesh % faces(iFace) % swapDimensions

               nBe = nBe + 1
               bCom(procID) % boundaryIDs(nBe) = iFaceLocal
               bCom(procID) % extProcIDs(nBe)  = p2 

               ! Boundary ID's associated with MPI boundaries are reported as positive integers
               ! so that they can be distinguished from physical boundaries that do not require
               ! communication with neighboring processes. 
               procMesh(procID) % faces(iFaceLocal) % boundaryID = nBe

               faceBoundaryIDs(iFace, faceProcCount(iFace)) = nBe

            ELSEIF( p2 /= p1 .AND. p2 == procID )THEN ! MPI Boundary
            
               nMPI = nMPI + 1 
            
               procMesh(procID) % faces(iFaceLocal) % faceID          = iFace
               procMesh(procID) % faces(iFaceLocal) % elementIDs(1)   = globalToLocal( e2, 1 ) 
               procMesh(procID) % faces(iFaceLocal) % elementIDs(2)   = -e1 
               procMesh(procID) % faces(iFaceLocal) % elementSides(1) = mesh % faces(iFace) % elementSides(1) 
               procMesh(procID) % faces(iFaceLocal) % elementSides(2) = mesh % faces(iFace) % elementSides(2)
               procMesh(procID) % faces(iFaceLocal) % iStart          = mesh % faces(iFace) % iStart
               procMesh(procID) % faces(iFaceLocal) % iInc            = mesh % faces(iFace) % iInc
               procMesh(procID) % faces(iFaceLocal) % jStart          = mesh % faces(iFace) % jStart
               procMesh(procID) % faces(iFaceLocal) % jInc            = mesh % faces(iFace) % jInc
               procMesh(procID) % faces(iFaceLocal) % swapDimensions  = mesh % faces(iFace) % swapDimensions
            
               nBe = nBe + 1
               bCom(procID) % boundaryIDs(nBe) = iFaceLocal
               bCom(procID) % extProcIDs(nBe)  = p1 

               ! Boundary ID's associated with MPI boundaries are reported as positive integers
               ! so that they can be distinguished from physical boundaries that do not require
               ! communication with neighboring processes. 
               procMesh(procID) % faces(iFaceLocal) % boundaryID = nBe

               faceBoundaryIDs(iFace, faceProcCount(iFace)) = nBe

         ENDIF

      ENDDO
      PRINT*, 'Process ID :',procID, ', nMPI   :', nMPI

 END SUBROUTINE DecomposeFaces
!
 SUBROUTINE SetupCommTables( )
   IMPLICIT NONE
   ! Local
   INTEGER :: procID, bID, extProc
   INTEGER :: nMPI(0:nProc-1)


      
      DO procID = 0, nProc-1

         nMPI = 0
         DO bID = 1, bCom(jProc) % nBoundaries

            extProc = bCom(procID) % extProcIDs(bID)
            IF( extProc /= procID )THEN

               ! nMPI is an array that keeps track of the order in which MPI messages will be sent
               ! to procID's neighbors
               nMPI(procID) = nMPI(procID) + 1

               localFaceID  = bCom(procID) % boundaryIDs(bID)
               globalFaceID = procMesh(procID) % faces(localFaceID) % faceID

               IF( faceProcOwners(globalFaceID, 1) == procID )THEN
                  extBID = faceBoundaryIDs(globalFaceID,2)
               ELSEIF( faceProcOwners(globalFaceID,2) == procID )THEN
                  extBID = faceBoundaryIDs(globalFaceID,1)
               ELSE
                  PRINT*, '  SetupCommTables : Something catastrophic happened !'
                  CALL Cleanup( )
                  STOP
               ENDIF


               bCom(extProc) % unPackMap(extBID) = nMPI(procID)
                  

            ENDIF

         ENDDO

      ENDDO   


 END SUBROUTINE SetupCommTables
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

END PROGRAM DecomposeStructuredHexMesh
