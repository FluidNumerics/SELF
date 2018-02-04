! CommunicationTable_CLASS.f90
!
! Copyright 2018 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE CommunicationTable_CLASS

  USE ModelPrecision
  USE BoundaryCommunicator_CLASS


  IMPLICIT NONE

  INCLUDE 'mpIF.h'

  TYPE CommunicationTable

    INTEGER :: myRank, nProc, nNeighbors
    INTEGER :: maxBufferSize
    INTEGER :: MPI_COMM, MPI_PREC, mpiErr

    INTEGER, ALLOCATABLE :: bufferMap(:)
    INTEGER, ALLOCATABLE :: neighborRank(:)
    INTEGER, ALLOCATABLE :: bufferSize(:)
    INTEGER, ALLOCATABLE :: rankTable(:)

#ifdef HAVE_CUDA
    INTEGER, DEVICE, ALLOCATABLE :: bufferMap_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: neighborRank_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: bufferSize_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: rankTable_dev(:)
#endif

  CONTAINS

    PROCEDURE :: Build => Build_CommunicationTable
    PROCEDURE :: Trash => Trash_CommunicationTable

    PROCEDURE :: ConstructCommTables

  END TYPE CommunicationTable

CONTAINS

  SUBROUTINE Build_CommunicationTable( myComm, extComm, setupSuccess )

    IMPLICIT NONE
    CLASS( CommunicationTable ), INTENT(out)    :: myComm
    TYPE( BoundaryCommunicator ), INTENT(inout) :: extComm
    LOGICAL, INTENT(inout)                      :: setupSuccess


    myComm % MPI_COMM = MPI_COMM_WORLD

    IF( prec == sp )THEN
      myComm % MPI_PREC = MPI_FLOAT
    ELSE
      myComm % MPI_PREC = MPI_DOUBLE
    ENDIF

    CALL MPI_INIT( myComm % mpiErr )
    CALL MPI_COMM_RANK( myComm % MPI_COMM, myComm % myRank, myComm % mpiErr )
    CALL MPI_COMM_SIZE( myComm % MPI_COMM, myComm % nProc, myComm % mpiErr )

    PRINT*, '    S/R Build_CommunicationTable : Greetings from Process ', myComm % myRank+1, ' of ', myComm % nProc

    ALLOCATE( myComm % bufferMap(1:extComm % nBoundaries), &
      myComm % rankTable(0:myComm % nProc-1) )

#ifdef HAVE_CUDA

    ALLOCATE( myComm % bufferMap_dev(1:extComm % nBoundaries), &
      myComm % rankTable_dev(0:myComm % nProc-1) )

#endif

    CALL myComm % ConstructCommTables( extComm )

#ifdef HAVE_CUDA

    ALLOCATE( myComm % myRank_dev, &
      myComm % nProc_dev, &
      myComm % nNeighbors_dev )

    myComm % myRank_dev      = myComm % myRank
    myComm % nProc_dev       = myComm % nProc
    myComm % nNeighbors_dev  = myComm % nNeighbors

#endif


  END SUBROUTINE Build_CommunicationTable

!

  SUBROUTINE Trash_CommunicationTable( myComm )

    IMPLICIT NONE
    CLASS( CommunicationTable ), INTENT(inout) :: myComm


    PRINT*, '    S/R Trash_CommunicationTable : Clearing memory.'

    DEALLOCATE( myComm % neighborRank, &
      myComm % bufferSize, &
      myComm % bufferMap, &
      myComm % rankTable )

#ifdef HAVE_CUDA

    DEALLOCATE( myComm % neighborRank_dev, &
      myComm % bufferSize_dev, &
      myComm % bufferMap_dev, &
      myComm % rankTable_dev )


#endif

    CALL MPI_FINALIZE( myComm % mpiErr )


  END SUBROUTINE Trash_CommunicationTable

!

  SUBROUTINE ConstructCommTables( myComm, extComm )
    IMPLICIT NONE
    CLASS( CommunicationTable ), INTENT(inout)  :: myComm
    TYPE( BoundaryCommunicator ), INTENT(inout) :: extComm
    ! Local
    INTEGER, ALLOCATABLE :: bufferCounter(:)
    INTEGER :: sharedFaceCount(0:myComm % nProc-1)
    INTEGER :: IFace, bID, iNeighbor
    INTEGER :: tag, ierror
    INTEGER :: e1, e2, s1, p2, nmsg, maxFaceCount
    INTEGER :: fUnit


    ! Count up the number of neighboring ranks
    myComm % rankTable = 0
    sharedFaceCount     = 0
    DO bID = 1, extComm % nBoundaries

      p2 = extComm % extProcIDS(bID)

      IF( p2 /= myComm % myRank )THEN

        myComm % rankTable(p2) = 1
        sharedFaceCount(p2) = sharedFaceCount(p2)+1

      ENDIF

    ENDDO


    myComm % nNeighbors = SUM( myComm % rankTable )
    PRINT*, '  S/R ConstructCommTables : Found', myComm % nNeighbors, 'neighbors for Rank', myComm % myRank+1

    ALLOCATE( myComm % neighborRank(1:myComm % nNeighbors), &
      myComm % bufferSize(1:myComm % nNeighbors), &
      bufferCounter(1:myComm % nNeighbors) )


    ! For each neighbor, set the neighbor's rank
    iNeighbor = 0
    DO p2 = 0, myComm % nProc-1

      IF( myComm % rankTable(p2) == 1 )THEN

        iNeighbor = iNeighbor + 1
        myComm % neighborRank(iNeighbor) = p2
        myComm % rankTable(p2) = iNeighbor

      ENDIF

    ENDDO


    maxFaceCount = MAXVAL( sharedFaceCount )
    DO iNeighbor = 1, myComm % nNeighbors

      p2 = myComm % neighborRank(iNeighbor)
      myComm % bufferSize(iNeighbor) = sharedFaceCount(p2)

    ENDDO


    myComm % maxBufferSize = maxFaceCount
    bufferCounter = 0



    myComm % bufferMap = 0

    DO bID = 1, extComm % nBoundaries

      p2 = extComm % extProcIDs(bID)

      ! In the event that the external process ID (p2) is identical to the current rank (p1),
      ! THEN this boundary edge involves a physical boundary condition and DOes not require a
      ! message exchange

      IF( p2 /= myComm % myRank )THEN

        iNeighbor = myComm % rankTable(p2)

        bufferCounter(iNeighbor) = bufferCounter(iNeighbor) + 1
        myComm % bufferMap(bID)   = bufferCounter(iNeighbor)

      ENDIF

    ENDDO

    DEALLOCATE( bufferCounter )

#ifdef HAVE_CUDA
    ALLOCATE( myComm % neighborRank_dev(1:myComm % nNeighbors), &
      myComm % bufferSize_dev(1:myComm % nNeighbors) )

    myDGSEM % mpiPackets % rankTable_dev    = myDGSEM % mpiPackets % rankTable
    myDGSEM % mpiPackets % neighborRank_dev = myDGSEM % mpiPackets % neighborRank
    myDGSEM % mpiPackets % bufferSize_dev   = myDGSEM % mpiPackets % bufferSize
    myDGSEM % mpiPackets % bufferMap_dev    = myDGSEM % mpiPackets % bufferMap
#endif

  END SUBROUTINE ConstructCommTables

END MODULE CommunicationTable_CLASS


