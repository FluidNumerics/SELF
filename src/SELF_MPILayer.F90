MODULE SELF_MPILayer

USE SELF_Data
USE SELF_MappedData
USE SELF_Mesh


IMPLICIT NONE 

  TYPE MPILayer

    INTEGER :: nMessages
    INTEGER :: mpiComm
    INTEGER :: mpiPrec
    INTEGER :: myRank
    INTEGER :: nProc
    INTEGER, ALLOCATABLE :: requestHandle(:), requestStats(:,:) 

    CONTAINS

      PROCEDURE :: Build => Build_MPILayer
      PROCEDURE :: Trash => Trash_MPILayer

      GENERIC, PUBLIC :: MPI_Exchange => MPIExchange_MappedScalar2D
      PROCEDURE, PUBLIC :: MPIExchange_MappedScalar2D

      PROCEDURE :: Finalize_MPI_Exchange

  END TYPE MPILayer

  PUBLIC :: SELF_CalculateRankDistributions

CONTAINS


  FUNCTION Int2Str( aNumber ) RESULT( aString )
    IMPLICIT NONE
    INTEGER :: aNumber
    CHARACTER(12) :: aString

      WRITE(aString, '(I)') aNumber

  END FUNCTION Int2Str

  SUBROUTINE SELF_CalculateRankDistributions( nObj, nRanks, objBounds )
#undef __FUNC__
#define __FUNC__ SELF_CalculateRankDistributions
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nObj
    INTEGER, INTENT(in) :: nRanks
    INTEGER, INTENT(out) :: objBounds(1:2,0:nRanks-1)
    ! Local
    INTEGER :: rank, nObjPerRank, remainder, lastUpperBound

      nObjPerRank = nObj/nRanks
      remainder = nObj - nObjPerRank*nRanks

      lastUpperBound = 0
      DO rank = 0, nRanks - 1

        objBounds(1,rank) = 1 + lastUpperBound
        IF( rank < remainder )THEN
          objBounds(2,rank) = objBounds(1,rank) + nObjPerRank + 1
        ELSE
          objBounds(2,rank) = objBounds(1,rank) + nObjPerRank
        ENDIF
        lastUpperBound = objBounds(2,rank)

        INFO( 'Rank '//TRIM(Int2Str(rank))//&
              ' [lower,upper/total] = '//&
              TRIM(Int2Str(objBounds(1,rank)))//','//&
              TRIM(Int2Str(objBounds(2,rank)))//'/'//&
              TRIM(Int2Str(nObj)) )

      ENDDO

  END SUBROUTINE SELF_CalculateRankDistributions

  SUBROUTINE MPIExchange_MappedScalar2D( mpiHandler, scalar, mesh, workScalar )
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(MappedScalar2D), INTENT(in) :: scalar
    TYPE(Mesh2D), INTENT(in) :: mesh
    TYPE(MappedScalar2D), INTENT(inout) :: scalar
    ! Local
    INTEGER :: msgCnt, e1, s1, e2, s2
    INTEGER :: globalSideId, externalProcId

    msgCnt = 0
    DO e1 = 1, mesh % nElem
      s1 = 1
      DO sideId = mesh % elemInfo % hostData(3,e1)+1, mesh % elemInfo % hostData(4,e1)
        ! Secondary element ID for this face
        e2 = mesh % sideInfo % hostData(3,sideId) 

        ! In SELF, we read in HOPR pre-processed mesh information. Upon reading in and
        ! performing data decomposition, we set e2 = -e2 if the neighboring element is
        ! owned by another rank
        IF( e2 < 0 )THEN 
          s2 = mesh % sideInfo % hostData(4,sideId)/10
          globalSideId = ABS(mesh % sideInfo % hostdata(2,sideId))

          ! Assume that mesh has been pre-processed with HOPR and
          ! elements are distributed to ranks by simply dividing the list [1:nElements] 
          ! evenly and in order to all ranks
          ! TODO: externalProcId = GetRank( e2, mesh % nElem, mpiHandler % nRanks ) 

          msgCnt = msgCnt + 1
          ! Receive data on this rank's workScalar
          CALL MPI_IRECV( workScalar % boundary % hostData(:,:,s1,e1), &
                          (scalar % N+1)*scalar % nVar, &
                          mpiHandler % mpiPrec, &
                          externalProcId, globalSideId,  &
                          mpiHandler % mpiComm, &
                          mpiHandler % requestHandle(msgCnt*2-1), iError )
  

          ! Send data from this rank's scalar
          CALL MPI_ISEND( workScalar % boundary % hostData(:,:,s1,e1), &
                          (scalar % N+1)*scalar % nVar, &
                          mpiHandler % mpiPrec, &
                          externalProcID, globalSideId, &
                          mpiHandler % mpiComm, &
                          mpiHandler % requestHandle(msgCnt*2), iError)

        ENDIF
      
      ENDDO
    ENDDO

  END FUNCTION MPIExchange_MappedScalar2D

  SUBROUTINE Finalize_MPIExchange( mpiHandler )
    CLASS( MPILayer ), INTENT(inout) :: mpiHandler
    ! Local
    INTEGER    :: ierror

    CALL MPI_WaitAll( mpiHandler % nMessages*2, &
                      mpiHandler % requestHandle, &
                      mpiHandler % requestStats, &
                      iError)

  END SUBROUTINE Finalize_MPI_Exchange

END MODULE SELF_MPILayer
