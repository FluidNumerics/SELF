! MeshGenerator_3D
!
! Copyright 2018 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.

PROGRAM MeshGenerator_3D

  USE ModelPrecision
  USE ConstantsDictionary
  USE commonroutines
  USE Timing
  USE NodalDG_CLASS
  USE HexMesh_CLASS
  USE TopographicShapes
  USE FluidParams_CLASS
  USE BoundaryCommunicator_CLASS


  IMPLICIT NONE

  TYPE( NodalDG )                           :: nodal
  TYPE( HexMesh )                           :: mesh
  TYPE( HexMesh ), ALLOCATABLE              :: procMesh(:)
  TYPE( FluidParams )                       :: params
  TYPE( MultiTimers )                       :: timers
  TYPE( BoundaryCommunicator ), ALLOCATABLE :: bcom(:)
  INTEGER, ALLOCATABLE                      :: faceProcCount(:), faceProcOwners(:,:), faceBoundaryIDs(:,:)
  INTEGER                      :: procID
  INTEGER                      :: pfnodes(1:nQuadNodes), gfNodes(1:nQuadNodes), nIDs(1:nQuadNodes)
  INTEGER                      :: nElems, nProc, pID, i
  INTEGER                      :: iEl, jEl, iSide, iNode, nAdj, nID, elID, nBe, nMPI
  INTEGER                      :: e2, IFace, jFace, localID, eIDs(1:2), pIDs(1:2)
  INTEGER                      :: nxp, nyp, nzp, iPx, iPy, iPz, iXp, iYp, iZp, iX, iY, iZ
  INTEGER, ALLOCATABLE         :: globalToLocal(:,:), nElPerProc(:), nodeLogic(:,:), nNodePerProc(:), partitions(:)
  INTEGER, ALLOCATABLE         :: globalToLocalNode(:,:)
  REAL(prec), ALLOCATABLE      :: materials(:)
  LOGICAL                      :: setupSuccess


  CALL Setup( )

  IF( setupSuccess )THEN

    CALL timers % StartTimer( 1 )

    CALL PartitionElementsAndNodes_Structured( )


    ! Now we generate the local mesh for each process
    DO procID = 0, nProc-1

      CALL procMesh(procID) % Build(  nNodePerProc(procID), &
        nElPerProc(procID), &
        1, params % polyDeg )

      CALL DecomposeNodes( procID )

      CALL DecomposeElements( procID )

      CALL DecomposeFaces( procID )

      CALL FileIO( procID )

    ENDDO

    CALL SetupCommTables( )

    CALL timers % STOPTimer( 1 )

    CALL Cleanup( )
  ENDIF

CONTAINS

  SUBROUTINE Setup( )
    IMPLICIT NONE

    ! Read in the parameters
    CALL params % Build( setupSuccess )
    IF( .NOT. setupSuccess )THEN
      RETURN
    ENDIF
    CALL timers % Build( )

    CALL timers % AddTimer( 'Total Time', 1 )
    ! Build an interpolant
    CALL nodal % Build( targetPoints = UnIFormPoints(-1.0_prec,1.0_prec,params % nPlot), &
      N = params % polyDeg, &
      nTargetPoints = params % nPlot, &
      quadrature = GAUSS )

    IF( params % topographicShape == Gaussian )THEN
      TopographicShape => GaussianHill
    ELSE
      TopographicShape => DefaultTopography
    ENDIF
    ! Build the Geometry
    IF( params % MeshTYPE == DOublyPeriodic )THEN
      PRINT*,' Constructing DOubly periodic Structured mesh.'
      CALL mesh % ConstructStructuredMesh( nodal % interp, &
        params % nXelem, &
        params % nYelem, &
        params % nZelem, &
        .TRUE. )
    ELSE
      PRINT*,' Constructing Structured mesh.'
      CALL mesh % ConstructStructuredMesh( nodal % interp, &
        params % nXelem, &
        params % nYelem, &
        params % nZelem, &
        .FALSE. )
    ENDIF

    nElems = mesh % elements % nElements
    nProc  = params % nProcX*params % nProcY*params % nProcZ
    PRINT*, nProc
    IF( nProc == 0 )THEN
      nProc = 1
    ENDIF
    ALLOCATE( materials(1:nElems), &
      globalToLocal(1:nElems,1:2), &
      nElPerProc(0:nProc-1), &
      nodeLogic(1:mesh % nodes % nNodes, 0:nProc-1), &
      nNodePerProc(0:nProc-1), &
      globalToLocalNode(1:mesh % nodes % nNodes, 0:nProc-1) )

    materials     = 0.0_prec
    globalToLocal = 0
    nElPerProc    = 0
    nodeLogic     = 0
    nNodePerProc  = 0

    ALLOCATE( partitions(1:nElems) )
    partitions = 0
    ALLOCATE( bcom(0:nProc-1), procMesh(0:nProc-1) )
    ALLOCATE( faceProcCount(1:mesh % faces % nFaces), &
      faceProcOwners(1:mesh % faces % nFaces,1:2), &
      faceBoundaryIDs(1:mesh % faces % nFaces,1:2) )

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
    CALL procmesh(procID) % WriteSELFMeshFile( TRIM(params % PeaceMeshFile)//'.'//pIDChar )


  END SUBROUTINE FileIO
!
  SUBROUTINE PartitionElementsAndNodes_Structured( )
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
        nElPerProc(procID)   = nElPerProc(procID) + 1
        globalToLocal(iEl,1) = nElPerProc(procID)
        globalToLocal(iEl,2) = procID

        DO iNode = 1, 8
          nID = mesh % elements % nodeIDs(iNode,iEl)
          nodeLogic(nID,procID) = 1
        ENDDO

      ENDDO
      DO nID = 1, mesh % nodes % nNodes
        nNodePerProc = nNodePerProc + nodeLogic(nID,:)
        globalToLocalNode(nID,:) = nNodePerProc
      ENDDO

      DO procID = 0, nProc-1
        PRINT*, 'Process ID :',procID, ', nElems :', nElPerProc(procID)
        PRINT*, 'Process ID :',procID, ', nNodes :', nNodePerProc(procID)
      ENDDO

    ELSE

      DO iEl = 1, nElems
        nElPerProc(0) = mesh % elements % nElements
        globalToLocal(iEl,1) = iEl ! Local ID is global ID
        globalToLocal(iEl,2) = 0   ! process ID

        DO iNode = 1, 8
          nID = mesh % elements % nodeIDs(iNode,iEl)
          nodeLogic(nID,0) = 1
        ENDDO
      ENDDO

      DO nID = 1, mesh % nodes % nNodes
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

    DO nID = 1, mesh % nodes % nNodes

      ! We only assign node information IF it is "owned" by this process
      IF( nodeLogic(nID,procID) == 1 )THEN
        ! Obtain the local node ID from the "global-to-local" array
        localID = globalToLocalNode( nID, procID )
        ! Assign the local node attributes
        procMesh(procID) % nodes % nodeID(localID)   = nID
        procMesh(procID) % nodes % nodeTYPE(localID) = mesh % nodes % nodeTYPE(nID)
        procMesh(procID) % nodes % x(1:3,localID)    = mesh % nodes % x(1:3,nID)
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
    DO iEl = 1, mesh % elements % nElements

      ! We only assign element information IF it is "owned" by this process

      IF( globalToLocal(iEl,2) == procID )THEN
        elID = globalToLocal(iEl,1) ! Obtain the local element ID from the global-to-local array

        DO iNode = 1, 8

          ! First we grab the global node ID for this element's corner node
          nID = mesh % elements % nodeIDs(iNode,iEl)
          ! THEN we convert this to a local node ID for this process
          localID = globalToLocalNode( nID, procID )
          ! And now we assign the local Node ID to the local element's corner node ID's
          procMesh(procID) % elements % nodeIDs(iNode,elID)  = localID

        ENDDO
        procMesh(procID) % elements % elementID(elID) = iEl

        ! Set the geometry
        procMesh(procID) % elements % nHat(:,:,:,:,elID)   = mesh % elements % nHat(:,:,:,:,iEl)
        procMesh(procID) % elements % xBound(:,:,:,:,elID) = mesh % elements % xBound(:,:,:,:,iEl)
        procMesh(procID) % elements % x(:,:,:,:,elID)      = mesh % elements % x(:,:,:,:,iEl)
        procMesh(procID) % elements % J(:,:,:,elID)        = mesh % elements % J(:,:,:,iEl)
        procMesh(procID) % elements % dxds(:,:,:,elID)     = mesh % elements % dxds(:,:,:,iEl)
        procMesh(procID) % elements % dxdp(:,:,:,elID)     = mesh % elements % dxdp(:,:,:,iEl)
        procMesh(procID) % elements % dxdq(:,:,:,elID)     = mesh % elements % dxdq(:,:,:,iEl)
        procMesh(procID) % elements % dyds(:,:,:,elID)     = mesh % elements % dyds(:,:,:,iEl)
        procMesh(procID) % elements % dydp(:,:,:,elID)     = mesh % elements % dydp(:,:,:,iEl)
        procMesh(procID) % elements % dydq(:,:,:,elID)     = mesh % elements % dydq(:,:,:,iEl)
        procMesh(procID) % elements % dzds(:,:,:,elID)     = mesh % elements % dzds(:,:,:,iEl)
        procMesh(procID) % elements % dzdp(:,:,:,elID)     = mesh % elements % dzdp(:,:,:,iEl)
        procMesh(procID) % elements % dzdq(:,:,:,elID)     = mesh % elements % dzdq(:,:,:,iEl)
        procMesh(procID) % elements % Ja(:,:,:,:,:,elID)   = mesh % elements % Ja(:,:,:,:,:,iEl)

      ENDIF

    ENDDO

  END SUBROUTINE DecomposeElements
!
  SUBROUTINE DecomposeFaces( procID )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: procID
    INTEGER             :: e1, e2, p1, p2, IFace, IFaceLocal
    INTEGER             :: nBe, nMPI, nLocalFaces

    nLocalFaces = 0
    nBe         = 0
    DO IFace = 1, mesh % faces % nFaces

      e1 = mesh % faces % elementIDs(1,IFace)
      e2 = mesh % faces % elementIDs(2,IFace)

      p1 = globalToLocal(e1,2)
      IF( e2 > 0 )THEN
        p2 = globalToLocal(e2,2)
      ELSE
        p2  = p1
      ENDIF

      IF( p1 == procID .OR. p2 == procID )THEN
        nLocalFaces = nLocalFaces + 1
        IF( e2 < 0 .OR. p2 /= p1 )THEN
          nBe = nBe + 1
        ENDIF
      ENDIF

    ENDDO

    PRINT*, '========================================================================='
    PRINT*, 'Process ID :',procID, ', nFaces :', nLocalFaces
    PRINT*, 'Process ID :',procID, ', nBFace :', nBe

    CALL procMesh(procID) % faces % Trash( )
    CALL procMesh(procID) % faces % Build( nLocalFaces, mesh % elements % N )
    CALL bCom(procID) % Initialize( nBe )

    IFaceLocal = 0
    nBe        = 0
    nMPI       = 0
    DO IFace = 1, mesh % faces % nFaces

      e1 = mesh % faces % elementIDs(1,IFace)
      e2 = mesh % faces % elementIDs(2,IFace)

      p1 = globalToLocal(e1,2)
      IF( e2 > 0 )THEN
        p2 = globalToLocal(e2,2)
      ELSE
        p2  = p1
      ENDIF

      IF( p1 == procID .OR. p2 == procID )THEN
        IFaceLocal = IFaceLocal + 1
        faceProcCount(IFace) = faceProcCount(IFace) + 1
        IF( faceProcCount(IFace) > 2 )THEN
          PRINT*, ' DecomposeFaces : Something catastrophic happened !'
          CALL Cleanup( )
          STOP
        ENDIF
        faceProcOwners(IFace, faceProcCount(IFace)) = procID

        IF( p2 == p1 .AND. e2 > 0 )THEN ! Internal Face

          procMesh(procID) % faces % faceID(IFaceLocal)         = IFace
          procMesh(procID) % faces % elementIDs(1,IFaceLocal)   = globalToLocal( e1, 1 )
          procMesh(procID) % faces % elementIDs(2,IFaceLocal)   = globalToLocal( e2, 1 )
          procMesh(procID) % faces % elementSides(1,IFaceLocal) = mesh % faces % elementSides(1,IFace)
          procMesh(procID) % faces % elementSides(2,IFaceLocal) = mesh % faces % elementSides(2,IFace)
          procMesh(procID) % faces % iStart(IFaceLocal)         = mesh % faces % iStart(IFace)
          procMesh(procID) % faces % iInc(IFaceLocal)           = mesh % faces % iInc(IFace)
          procMesh(procID) % faces % jStart(IFaceLocal)         = mesh % faces % jStart(IFace)
          procMesh(procID) % faces % jInc(IFaceLocal)           = mesh % faces % jInc(IFace)
          procMesh(procID) % faces % swapDimensions(IFaceLocal) = mesh % faces % swapDimensions(IFace)
          procMesh(procID) % faces % iMap(:,:,IFaceLocal)       = mesh % faces % iMap(:,:,IFace)
          procMesh(procID) % faces % jMap(:,:,IFaceLocal)       = mesh % faces % jMap(:,:,IFace)

        ELSEIF( p2 == p1 .AND. e2 < 0 )THEN ! Physical boundary

          procMesh(procID) % faces % faceID(IFaceLocal)         = IFace
          procMesh(procID) % faces % elementIDs(1,IFaceLocal)   = globalToLocal( e1, 1 )
          procMesh(procID) % faces % elementIDs(2,IFaceLocal)   = e2
          procMesh(procID) % faces % elementSides(1,IFaceLocal) = mesh % faces % elementSides(1,IFace)
          procMesh(procID) % faces % elementSides(2,IFaceLocal) = mesh % faces % elementSides(2,IFace)
          procMesh(procID) % faces % iStart(IFaceLocal)         = mesh % faces % iStart(IFace)
          procMesh(procID) % faces % iInc(IFaceLocal)           = mesh % faces % iInc(IFace)
          procMesh(procID) % faces % jStart(IFaceLocal)         = mesh % faces % jStart(IFace)
          procMesh(procID) % faces % jInc(IFaceLocal)           = mesh % faces % jInc(IFace)
          procMesh(procID) % faces % swapDimensions(IFaceLocal) = mesh % faces % swapDimensions(IFace)
          procMesh(procID) % faces % iMap(:,:,IFaceLocal)       = mesh % faces % iMap(:,:,IFace)
          procMesh(procID) % faces % jMap(:,:,IFaceLocal)       = mesh % faces % jMap(:,:,IFace)

          nBe = nBe + 1
          bCom(procID) % boundaryIDs(nBe) = IFaceLocal
          bCom(procID) % extProcIDs(nBe)  = procID

          ! Physical boundary ID's will be set to negative for physical boundaries
          ! so that boundaries requiring communication can be separated from physical
          ! boundaries. This will allow more operations to be run concurrently with
          ! communication.
          procMesh(procID) % faces % boundaryID(IFaceLocal) = -nBe

        ELSEIF( p2 /= p1 .AND. procID == p1 )THEN ! MPI Boundary

          nMPI = nMPI + 1

          procMesh(procID) % faces % faceID(IFaceLocal)         = IFace
          procMesh(procID) % faces % elementIDs(1,IFaceLocal)   = globalToLocal( e1, 1 )
          procMesh(procID) % faces % elementIDs(2,IFaceLocal)   = -e2
          procMesh(procID) % faces % elementSides(1,IFaceLocal) = mesh % faces % elementSides(1,IFace)
          procMesh(procID) % faces % elementSides(2,IFaceLocal) = mesh % faces % elementSides(2,IFace)
          procMesh(procID) % faces % iStart(IFaceLocal)         = mesh % faces % iStart(IFace)
          procMesh(procID) % faces % iInc(IFaceLocal)           = mesh % faces % iInc(IFace)
          procMesh(procID) % faces % jStart(IFaceLocal)         = mesh % faces % jStart(IFace)
          procMesh(procID) % faces % jInc(IFaceLocal)           = mesh % faces % jInc(IFace)
          procMesh(procID) % faces % swapDimensions(IFaceLocal) = mesh % faces % swapDimensions(IFace)
          procMesh(procID) % faces % iMap(:,:,IFaceLocal)       = mesh % faces % iMap(:,:,IFace)
          procMesh(procID) % faces % jMap(:,:,IFaceLocal)       = mesh % faces % jMap(:,:,IFace)


          nBe = nBe + 1
          bCom(procID) % boundaryIDs(nBe) = IFaceLocal
          bCom(procID) % extProcIDs(nBe)  = p2

          ! Boundary ID's associated with MPI boundaries are reported as positive INTEGERs
          ! so that they can be distinguished from physical boundaries that DO not require
          ! communication with neighboring processes.
          procMesh(procID) % faces % boundaryID(IFaceLocal) = nBe

          faceBoundaryIDs(IFace, faceProcCount(IFace)) = nBe

        ELSEIF( p2 /= p1 .AND. procID == p2 )THEN ! MPI Boundary

          nMPI = nMPI + 1

          procMesh(procID) % faces % faceID(IFaceLocal)         = IFace
          procMesh(procID) % faces % elementIDs(1,IFaceLocal)   = globalToLocal( e2, 1 )
          procMesh(procID) % faces % elementIDs(2,IFaceLocal)   = -e1
          procMesh(procID) % faces % elementSides(1,IFaceLocal) = mesh % faces % elementSides(2,IFace)
          procMesh(procID) % faces % elementSides(2,IFaceLocal) = mesh % faces % elementSides(1,IFace)
          procMesh(procID) % faces % iStart(IFaceLocal)         = mesh % faces % iStart(IFace)
          procMesh(procID) % faces % iInc(IFaceLocal)           = mesh % faces % iInc(IFace)
          procMesh(procID) % faces % jStart(IFaceLocal)         = mesh % faces % jStart(IFace)
          procMesh(procID) % faces % jInc(IFaceLocal)           = mesh % faces % jInc(IFace)
          procMesh(procID) % faces % swapDimensions(IFaceLocal) = mesh % faces % swapDimensions(IFace)
          procMesh(procID) % faces % iMap(:,:,IFaceLocal)       = mesh % faces % iMap(:,:,IFace)
          procMesh(procID) % faces % jMap(:,:,IFaceLocal)       = mesh % faces % jMap(:,:,IFace)

          nBe = nBe + 1
          bCom(procID) % boundaryIDs(nBe) = IFaceLocal
          bCom(procID) % extProcIDs(nBe)  = p1

          ! Boundary ID's associated with MPI boundaries are reported as positive INTEGERs
          ! so that they can be distinguished from physical boundaries that DO not require
          ! communication with neighboring processes.
          procMesh(procID) % faces % boundaryID(IFaceLocal) = nBe

          faceBoundaryIDs(IFace, faceProcCount(IFace)) = nBe

        ENDIF
      ENDIF

    ENDDO
    PRINT*, 'Process ID :',procID, ', nMPI   :', nMPI

  END SUBROUTINE DecomposeFaces
!
  SUBROUTINE SetupCommTables( )
    IMPLICIT NONE
    ! Local
    INTEGER :: procID, bID, extProc, extBID
    INTEGER :: globalFaceID, localFaceID
    INTEGER :: nMPI(0:nProc-1)
    CHARACTER(4) :: pIDChar



    DO procID = 0, nProc-1

      nMPI = 0
      DO bID = 1, bCom(procID) % nBoundaries

        extProc = bCom(procID) % extProcIDs(bID)
        IF( extProc /= procID )THEN

          ! nMPI is an array that keeps track of the order in which MPI messages will be sent
          ! to procID's neighbors
          nMPI(extProc) = nMPI(extProc) + 1
          localFaceID  = bCom(procID) % boundaryIDs(bID)
          globalFaceID = procMesh(procID) % faces % faceID(localFaceID)

          IF( faceProcOwners(globalFaceID, 1) == procID )THEN
            extBID = faceBoundaryIDs(globalFaceID,2)
          ELSEIF( faceProcOwners(globalFaceID,2) == procID )THEN
            extBID = faceBoundaryIDs(globalFaceID,1)
          ELSE
            PRINT*, '  SetupCommTables : Something catastrophic happened !'
            CALL Cleanup( )
            STOP
          ENDIF


          bCom(extProc) % unPackMap(extBID) = nMPI(extProc)


        ENDIF

      ENDDO

    ENDDO

    DO procID = 0, nProc-1
      WRITE( pIDChar, '(I4.4)' ) procID
      CALL bCom(procID) % WritePickup( 'ExtComm.'//pIDChar )
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

ENDPROGRAM MeshGenerator_3D
