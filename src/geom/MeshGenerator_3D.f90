
! Copyright 2018 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.

PROGRAM MeshGenerator_3D

 USE ModelPrecision
 USE ConstantsDictionary
 USE CommonRoutines
 USE Timing
 USE NodalDG_Class
 USE HexMesh_Class
 USE TopographicShapes
 USE ModelParameters_Class
 USE BoundaryCommunicator_Class
 
  
 IMPLICIT NONE

 TYPE( NodalDG )                           :: nodal
 TYPE( HexMesh )                           :: mesh
 TYPE( HexMesh )                           :: procMesh
 TYPE( ModelParameters )                   :: params
 TYPE( MultiTimers )                       :: timers
 TYPE( BoundaryCommunicator ), ALLOCATABLE :: bcom(:)
 INTEGER, ALLOCATABLE                      :: faceProcCount(:), faceProcTable(:,:), faceProcOwners(:,:), faceBoundaryIDs(:,:)
 INTEGER                      :: procID
 INTEGER                      :: nElems, nProc, pID, i
 INTEGER                      :: iEl, jEl, iSide, iNode, nAdj, nID, elID, nMPI
 INTEGER                      :: e1, e2, p1, p2, iFace, iFaceLocal, jFace, localID
 INTEGER                      :: nBe, npFaces
 INTEGER                      :: globalFaceID, localFaceID, bID, extProc, extBID
 CHARACTER(4)                 :: pIDChar
 INTEGER, ALLOCATABLE         :: globalToLocal(:,:), nElPerProc(:), nodeLogic(:,:), nNodePerProc(:), partitions(:)
 INTEGER, ALLOCATABLE         :: globalToLocalNode(:,:), nLocMPI(:)
 REAL(prec), ALLOCATABLE      :: materials(:)
 LOGICAL                      :: setupSuccess
 

      ! Read in the parameters
      CALL params % Build( setupSuccess )

      IF( setupSuccess )THEN

         CALL timers % Build( )
         
         CALL timers % AddTimer( 'Total Time', 1 )
         ! Build an interpolant
         CALL nodal % Build( targetPoints = UniformPoints(-1.0_prec,1.0_prec,params % nPlot), &
                             N = params % polyDeg, &
                             nTargetPoints = params % nPlot, &
                             quadrature = GAUSS_LOBATTO )
   
         IF( params % topographicShape == Gaussian )THEN
            TopographicShape => GaussianHill
         ELSE
            TopographicShape => DefaultTopography
         ENDIF
         ! Build the Geometry
         IF( TRIM( params % UCDMeshFile ) == '' )THEN

           IF( params % MeshType == DoublyPeriodic )THEN
              PRINT*,' Loading doubly periodic mesh.'
              CALL mesh % ConstructStructuredMesh( nodal % interp, &
                params % nXelem, &
                params % nYelem, &
                params % nZelem, &
                .TRUE. )
           ELSE
              PRINT*,' Loading default mesh.'
              CALL mesh % ConstructStructuredMesh( nodal % interp, &
                params % nXelem, &
                params % nYelem, &
                params % nZelem, &
                .FALSE. )
           ENDIF

         ELSE   

           CALL mesh % ReadTrellisUCDMeshFile( nodal % interp, TRIM( params % UCDMeshFile ) )           
           CALL mesh % WriteTecplot( 'mesh' )

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
         

         CALL timers % StartTimer( 1 )
   
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

 !        IF( TRIM( params % UCDMeshFile ) == '' )THEN

           CALL mesh % PartitionStructuredElementsAndNodes( params, partitions, nElPerProc, &
                                                            globalToLocal, nodeLogic, nNodePerProc, &
                                                            globalToLocalNode, nProc )

 !        ELSE


 !        ENDIF        
   
            
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

            PRINT*, '========================================================================='
            PRINT*, 'Process ID :',procID, ', nMPI   :', nMPI
            PRINT*, 'Process ID :',procID, ', nFaces :', npFaces
            PRINT*, 'Process ID :',procID, ', nBFace :', nBe

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
            CALL procMesh % WriteTecplot( 'mesh.'//pIDChar )
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
                  STOP
               ENDIF


               bCom(extProc) % unPackMap(extBID) = nlocMPI(extProc)
                  

            ENDIF

         ENDDO

      ENDDO   

      DO procID = 0, nProc-1
         WRITE( pIDChar, '(I4.4)' ) procID
         CALL bCom(procID) % WritePickup( 'ExtComm.'//pIDChar )
      ENDDO
   
         CALL timers % StopTimer( 1 )
         
! ----------------------------------------------------------------------------- !

      DO procID = 0, nProc-1
         CALL bCom(procID) % Trash( )
      ENDDO
      DEALLOCATE( bCom )
      

      ! For visualizing the decomposition
      materials = REAL( partitions, prec )
      CALL mesh % WriteMaterialTecplot( materials )
      
      CALL timers % Write_MultiTimers( )

      ! Clean up memory !
      DEALLOCATE( materials, &
                  faceProcCount, &
                  faceProcOwners, &
                  faceBoundaryIDs, &
                  nLocMPI )
      
      CALL mesh % Trash( )
      CALL nodal % Trash( )
      CALL timers % Trash( )


      ENDIF

 CONTAINS

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

END PROGRAM MeshGenerator_3D
