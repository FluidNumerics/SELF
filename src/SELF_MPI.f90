!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_MPI

  USE SELF_Constants
  USE SELF_Memory

  IMPLICIT NONE

#include "SELF_Macros.h"

  TYPE MPILayer
    LOGICAL :: mpiEnabled
    INTEGER :: mpiComm
    INTEGER :: mpiPrec
    INTEGER :: rankId
    INTEGER :: nRanks
    INTEGER :: nElem
    TYPE(hfInt32_r1) :: elemToRank
    TYPE(hfInt32_r1) :: offSetElem

    CONTAINS

      PROCEDURE :: Init => Init_MPILayer
!      PROCEDURE :: Free => Free_MPILayer

      PROCEDURE :: SetElemToRank
!
!      GENERIC, PUBLIC :: MPI_Exchange => MPIExchange_MappedScalar2D
!      PROCEDURE, PUBLIC :: MPIExchange_MappedScalar2D
!
!      PROCEDURE :: Finalize_MPI_Exchange
!
  END TYPE MPILayer


CONTAINS

  SUBROUTINE Init_MPILayer(this)
#undef __FUNC__
#define __FUNC__ "Init_MPILayer"
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(out) :: this
    ! Local
    INTEGER       :: ierror
    CHARACTER(30) :: msg

    this % mpiComm = 0
    this % mpiPrec = prec
    this % rankId = 0
    this % nRanks = 1
    this % nElem = 0
    this % mpiEnabled = .FALSE.
#ifdef MPI
    this % mpiEnabled = .TRUE.
    this % mpiComm = MPI_COMM_WORLD
    CALL MPI_INIT( ierror )
    CALL MPI_COMM_RANK( this % mpiComm, this % rankId, ierror )
    CALL MPI_COMM_SIZE( this % mpiComm, this % nRanks,  ierror )

    IF( prec == sp )THEN
      this % mpiPrec=MPI_FLOAT
    ELSE
      this % mpiPrec=MPI_DOUBLE
    ENDIF
#endif

    CALL this % offSetElem % Alloc(0,this % nRanks)

    WRITE(msg,'(I5)')this % rankId
    msg="Greetings from rank "//TRIM(msg)//"."
    INFO(TRIM(msg))
      
  END SUBROUTINE Init_MPILayer

  SUBROUTINE SetElemToRank(this, nElem)
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: this
    INTEGER, INTENT(in) :: nElem
    ! Local
    INTEGER :: iel
    
      this % nElem = nElem
      CALL this % elemToRank % Alloc(1,nElem)
      CALL DomainDecomp(nElem,&
                        this % nRanks,&
                        this % offSetElem % hostData)

      DO iel = 1, nElem
        CALL ElemToRank(this % nRanks, &
                        this % offSetElem % hostData, &
                        iel, &
                        this % elemToRank % hostData(iel))
      ENDDO

#ifdef GPU
      CALL this % offSetElem % UpdateDevice()
      CALL this % elemToRank % UpdateDevice()
#endif

  END SUBROUTINE SetElemToRank

  SUBROUTINE DomainDecomp(nElems,nDomains,offSetElem)
  ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 4
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nElems
    INTEGER, INTENT(in) :: nDomains
    INTEGER, INTENT(out) :: offsetElem(0:nDomains)
    ! Local
    INTEGER :: nLocalElems
    INTEGER :: remainElems
    INTEGER :: iDom

      nLocalElems = nElems/nDomains
      remainElems = nElems - nLocalElems*nDomains
      DO iDom = 0, nDomains-1
        offSetElem(iDom) = iDom*nLocalElems + MIN(iDom,remainElems)
      ENDDO
      offSetElem(nDomains) = nElems

  END SUBROUTINE DomainDecomp

  SUBROUTINE ElemToRank(nDomains,offsetElem,elemID,domain)
  ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 7
  !   "Find domain containing element index"
  !
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nDomains
    INTEGER, INTENT(in) :: offsetElem(0:nDomains)
    INTEGER, INTENT(in) :: elemID
    INTEGER, INTENT(out) :: domain
    ! Local
    INTEGER :: maxSteps
    INTEGER :: low,up,mid
    INTEGER :: i

      domain = 0
      maxSteps = INT(LOG10(REAL(nDomains))/LOG10(2.0))+1
      low = 0
      up = nDomains-1

      IF(offsetElem(low) < elemID .AND. elemID <= offsetElem(low+1))THEN
        domain = low
      ELSEIF(offsetElem(up) < elemID .AND. elemID <= offsetElem(up+1))THEN
        domain = up
      ELSE
        DO i = 1, maxSteps
          mid = (up-low)/2+low
          IF(offsetElem(mid) < elemID .AND. elemID <= offsetElem(mid+1))THEN
            domain = mid
            RETURN
          ELSEIF(elemID > offsetElem(mid+1))THEN
            low = mid+1
          ELSE
            up = mid
          ENDIF
        ENDDO
      ENDIF

  END SUBROUTINE ElemToRank
!  SUBROUTINE MPIExchange_MappedScalar2D( mpiHandler, scalar, mesh, workScalar )
!    IMPLICIT NONE
!    CLASS(MPILayer), INTENT(inout) :: mpiHandler
!    TYPE(MappedScalar2D), INTENT(in) :: scalar
!    TYPE(Mesh2D), INTENT(in) :: mesh
!    TYPE(MappedScalar2D), INTENT(inout) :: scalar
!    ! Local
!    INTEGER :: msgCnt, e1, s1, e2, s2
!    INTEGER :: globalSideId, externalProcId
!
!    msgCnt = 0
!    DO e1 = 1, mesh % nElem
!      s1 = 1
!      DO sideId = mesh % hopr_elemInfo % hostData(3,e1)+1, mesh % hopr_elemInfo % hostData(4,e1)
!        ! Secondary element ID for this face
!        e2 = mesh % hopr_sideInfo % hostData(3,sideId)
!
!        ! In SELF, we read in HOPR pre-processed mesh information. Upon reading in and
!        ! performing data decomposition, we set e2 = -e2 if the neighboring element is
!        ! owned by another rank
!        IF( e2 < 0 )THEN
!          s2 = mesh % hopr_sideInfo % hostData(4,sideId)/10
!          globalSideId = ABS(mesh % hopr_sideInfo % hostdata(2,sideId))
!
!          ! Assume that mesh has been pre-processed with HOPR and
!          ! elements are distributed to ranks by simply dividing the list [1:nElements]
!          ! evenly and in order to all ranks
!          ! TODO: externalProcId = GetRank( e2, mesh % nElem, mpiHandler % nRanks )
!
!          msgCnt = msgCnt + 1
!          ! Receive data on this rank's workScalar
!          CALL MPI_IRECV( workScalar % boundary % hostData(:,:,s1,e1), &
!                          (scalar % N+1)*scalar % nVar, &
!                          mpiHandler % mpiPrec, &
!                          externalProcId, globalSideId,  &
!                          mpiHandler % mpiComm, &
!                          mpiHandler % requestHandle(msgCnt*2-1), iError )
!
!
!          ! Send data from this rank's scalar
!          CALL MPI_ISEND( workScalar % boundary % hostData(:,:,s1,e1), &
!                          (scalar % N+1)*scalar % nVar, &
!                          mpiHandler % mpiPrec, &
!                          externalProcID, globalSideId, &
!                          mpiHandler % mpiComm, &
!                          mpiHandler % requestHandle(msgCnt*2), iError)
!
!        ENDIF
!
!      ENDDO
!    ENDDO
!
!  END FUNCTION MPIExchange_MappedScalar2D
!
!  SUBROUTINE Finalize_MPIExchange( mpiHandler )
!    CLASS( MPILayer ), INTENT(inout) :: mpiHandler
!    ! Local
!    INTEGER    :: ierror
!
!    CALL MPI_WaitAll( mpiHandler % nMessages*2, &
!                      mpiHandler % requestHandle, &
!                      mpiHandler % requestStats, &
!                      iError)
!
!  END SUBROUTINE Finalize_MPI_Exchange

END MODULE SELF_MPI
