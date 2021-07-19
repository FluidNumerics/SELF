!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_MPI

  USE SELF_Constants
  USE SELF_Memory
  USE ISO_C_BINDING

  IMPLICIT NONE

#include "SELF_Macros.h"

  TYPE MPILayer
    LOGICAL :: mpiEnabled
    INTEGER :: mpiComm
    INTEGER :: mpiPrec
    INTEGER :: rankId
    INTEGER :: nRanks
    INTEGER :: nElem
    INTEGER :: maxMsg
    INTEGER :: msgCount
    TYPE(hfInt32_r1) :: elemToRank
    TYPE(hfInt32_r1) :: offSetElem
    TYPE(hfInt32_r2) :: requests


    CONTAINS

      PROCEDURE :: Init => Init_MPILayer
      PROCEDURE :: Free => Free_MPILayer

      PROCEDURE :: SetElemToRank

      PROCEDURE,PUBLIC :: FinalizeMPIExchangeAsync

  END TYPE MPILayer



CONTAINS

  SUBROUTINE Init_MPILayer(this,maxMsg)
#undef __FUNC__
#define __FUNC__ "Init_MPILayer"
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(out) :: this
    INTEGER, INTENT(in) :: maxMsg
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
    CALL MPI_INIT(ierror)
    CALL MPI_COMM_RANK(this % mpiComm, this % rankId, ierror)
    CALL MPI_COMM_SIZE(this % mpiComm, this % nRanks,  ierror)

    IF(prec == sp)THEN
      this % mpiPrec=MPI_FLOAT
    ELSE
      this % mpiPrec=MPI_DOUBLE
    ENDIF
#endif

    CALL this % offSetElem % Alloc(0,this % nRanks)
    CALL this % requests % Alloc((/1,1/),&
                                 (/maxMsg,2/))
    this % maxMsg = maxMsg

    WRITE(msg,'(I5)')this % rankId
    msg="Greetings from rank "//TRIM(msg)//"."
    INFO(TRIM(msg))

  END SUBROUTINE Init_MPILayer

  SUBROUTINE Free_MPILayer(this)
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: this

      CALL this % offSetElem % Free()
      CALL this % requests % Free()
      CALL this % elemToRank % Free()

  END SUBROUTINE Free_MPILayer

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
  
  SUBROUTINE FinalizeMPIExchangeAsync(mpiHandler)
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    ! Local
    INTEGER :: ierror

#ifdef MPI
    CALL MPI_WaitAll(mpiHandler % msgCount, &
                      mpiHandler % requests(1:mpiHandler % msgCount,1), &
                      mpiHandler % requests(1:mpiHandler % msgCount,2), &
                      iError)
#endif

  END SUBROUTINE FinalizeMPIExchangeAsync

END MODULE SELF_MPI
