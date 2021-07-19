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
      GENERIC,PUBLIC :: MPIExchangeAsync => MPIExchangeAsync_MappedScalar2D,&
                                            MPIExchangeAsync_MappedVector2D,&
                                            MPIExchangeAsync_MappedTensor2D,&
                                            MPIExchangeAsync_MappedScalar3D,&
                                            MPIExchangeAsync_MappedVector3D,&
                                            MPIExchangeAsync_MappedTensor3D
      PROCEDURE,PRIVATE :: MPIExchangeAsync_MappedScalar2D
      PROCEDURE,PRIVATE :: MPIExchangeAsync_MappedVector2D
      PROCEDURE,PRIVATE :: MPIExchangeAsync_MappedTensor2D
      PROCEDURE,PRIVATE :: MPIExchangeAsync_MappedScalar3D
      PROCEDURE,PRIVATE :: MPIExchangeAsync_MappedVector3D
      PROCEDURE,PRIVATE :: MPIExchangeAsync_MappedTensor3D
      PROCEDURE,PUBLIC :: FinalizeMPIExchangeAsync

  END TYPE MPILayer

#ifdef GPU

  INTERFACE
    SUBROUTINE ApplyFlip_MappedScalar2D_gpu_wrapper(selfSideInfo,elemToRank,extBoundary,rankId,N,nVar,nEl) &
      bind(c,name="ApplyFlip_MappedScalar2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: selfSideInfo,elemToRank,extBoundary
      INTEGER,VALUE :: rankId,N,nVar,nEl
    END SUBROUTINE ApplyFlip_MappedScalar2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ApplyFlip_MappedVector2D_gpu_wrapper(selfSideInfo,elemToRank,extBoundary,rankId,N,nVar,nEl) &
      bind(c,name="ApplyFlip_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: selfSideInfo,elemToRank,extBoundary
      INTEGER,VALUE :: rankId,N,nVar,nEl
    END SUBROUTINE ApplyFlip_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ApplyFlip_MappedTensor2D_gpu_wrapper(selfSideInfo,elemToRank,extBoundary,rankId,N,nVar,nEl) &
      bind(c,name="ApplyFlip_MappedTensor2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: selfSideInfo,elemToRank,extBoundary
      INTEGER,VALUE :: rankId,N,nVar,nEl
    END SUBROUTINE ApplyFlip_MappedTensor2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ApplyFlip_MappedScalar3D_gpu_wrapper(selfSideInfo,elemToRank,extBoundary,rankId,N,nVar,nEl) &
      bind(c,name="ApplyFlip_MappedScalar3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: selfSideInfo,elemToRank,extBoundary
      INTEGER,VALUE :: rankId,N,nVar,nEl
    END SUBROUTINE ApplyFlip_MappedScalar3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ApplyFlip_MappedVector3D_gpu_wrapper(selfSideInfo,elemToRank,extBoundary,rankId,N,nVar,nEl) &
      bind(c,name="ApplyFlip_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: selfSideInfo,elemToRank,extBoundary
      INTEGER,VALUE :: rankId,N,nVar,nEl
    END SUBROUTINE ApplyFlip_MappedVector3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ApplyFlip_MappedTensor3D_gpu_wrapper(selfSideInfo,elemToRank,extBoundary,rankId,N,nVar,nEl) &
      bind(c,name="ApplyFlip_MappedTensor3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: selfSideInfo,elemToRank,extBoundary
      INTEGER,VALUE :: rankId,N,nVar,nEl
    END SUBROUTINE ApplyFlip_MappedTensor3D_gpu_wrapper
  END INTERFACE

#endif


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
  
  SUBROUTINE MPIExchangeAsync_MappedScalar2D(mpiHandler,selfSideInfo,scalar,resetCount,useDevicePtr)
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(hfInt32_r3), INTENT(in) :: selfSideInfo
    TYPE(MappedScalar2D), INTENT(inout) :: scalar
    LOGICAL, INTENT(in) :: resetCount
    LOGICAL, INTENT(in) :: useDevicePtr
    ! Local
    INTEGER :: e1, s1, e2, s2
    INTEGER :: globalSideId, externalProcId
    INTEGER :: msgCount

#ifdef MPI
    IF(resetCount)THEN
      msgCount = 0
    ELSE
      msgCount = mpiHandler % msgCount 
    ENDIF

    DO e1 = 1, scalar % nElem
      DO s1 = 1, 4

        e2 = selfSideInfo % hostData(3,s1,e1) ! Neighbor Element
        r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

        IF(r2 /= mpiHandler % rankId)THEN

          s2 = selfSideInfo % hostData(4,s1,e1)/10
          globalSideId = selfSideInfo % hostdata(2,s1,e1)

          IF(useDevicePtr)THEN
            msgCount = msgCount + 1
            CALL MPI_IRECV(scalar % extBoundary % deviceData(:,:,s1,e1), &
                            (scalar % N+1)*scalar % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcId, globalSideId,  &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

            msgCount = msgCount +1
            CALL MPI_ISEND(scalar % boundary % deviceData(:,:,s1,e1), &
                            (scalar % N+1)*scalar % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcID, globalSideId, &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)
          ELSE
            msgCount = msgCount + 1
            CALL MPI_IRECV(scalar % extBoundary % hostData(:,:,s1,e1), &
                            (scalar % N+1)*scalar % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcId, globalSideId,  &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

            msgCount = msgCount +1
            CALL MPI_ISEND(scalar % boundary % hostData(:,:,s1,e1), &
                            (scalar % N+1)*scalar % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcID, globalSideId, &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)
          ENDIF
        ENDIF

      ENDDO
    ENDDO

    mpiHandler % msgCount = msgCount
#endif

  END SUBROUTINE MPIExchangeAsync_MappedScalar2D
!
  SUBROUTINE ApplyFlip_MappedScalar2D(mpiHandler,selfSideInfo,scalar,gpuAccel)
  ! Apply side flips to sides where MPI exchanges took place.
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(hfInt32_r3), INTENT(in) :: selfSideInfo
    TYPE(MappedScalar2D), INTENT(inout) :: scalar
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1, s1, e2, s2
    INTEGER :: i, i2
    INTEGER :: globalSideId, externalProcId
    REAL(prec) :: extBuff(0:scalar % N)


    IF(gpuAccel)THEN

      CALL ApplyFlip_MappedScalar2D_gpu_wrapper(selfSideInfo % deviceData, &
                                                mpiHandler % elemToRank % deviceData, &
                                                scalar % extBoundary % deviceData, &
                                                mpiHandler % rankId, &
                                                scalar % N, &
                                                scalar % nVar, &
                                                scalar % nElem)


    ELSE
      DO e1 = 1, scalar % nElem
        DO s1 = 1, 4

          e2 = selfSideInfo % hostData(3,s1,e1) ! Neighbor Element
          r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

          IF(r2 /= mpiHandler % rankId)THEN

            s2 = selfSideInfo % hostData(4,s1,e1)/10
            flip = selfSideInfo % hostData(4,s1,e1) - s2*10
            globalSideId = selfSideInfo % hostdata(2,s1,e1)

            ! Need to update extBoundary with flip applied
            IF(flip == 1)THEN

              DO ivar = 1, scalar % nvar
                DO i = 0, scalar % N
                  i2 = scalar % N - i
                  extBuff(i) = scalar % extBoundary % hostData(i2,ivar,s1,e1)
                ENDDO
                DO i = 0, scalar % N
                  scalar % extBoundary % hostData(i,ivar,s1,e1) = extBuff(i)
                ENDDO
              ENDDO
  
            ENDIF
          ENDIF

        ENDDO
      ENDDO

  END SUBROUTINE ApplyFlip_MappedScalar2D

  SUBROUTINE MPIExchangeAsync_MappedVector2D(mpiHandler,selfSideInfo,vector,resetCount,useDevicePtr)
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(hfInt32_r3), INTENT(in) :: selfSideInfo
    TYPE(MappedVector2D), INTENT(inout) :: vector
    LOGICAL, INTENT(in) :: resetCount
    LOGICAL, INTENT(in) :: useDevicePtr
    ! Local
    INTEGER :: e1, s1, e2, s2
    INTEGER :: globalSideId, externalProcId
    INTEGER :: msgCount

#ifdef MPI
    IF(resetCount)THEN
      msgCount = 0
    ELSE
      msgCount = mpiHandler % msgCount 
    ENDIF

    DO e1 = 1, scalar % nElem
      DO s1 = 1, 4

        e2 = selfSideInfo % hostData(3,s1,e1) ! Neighbor Element
        r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

        IF(r2 /= mpiHandler % rankId)THEN

          s2 = selfSideInfo % hostData(4,s1,e1)/10
          globalSideId = selfSideInfo % hostdata(2,s1,e1)

          IF(useDevicePtr)THEN
            msgCount = msgCount + 1
            CALL MPI_IRECV(vector % extBoundary % deviceData(:,:,:,s1,e1), &
                           2*(vector % N+1)*vector % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcId, globalSideId,  &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

            msgCount = msgCount +1
            CALL MPI_ISEND(vector % boundary % deviceData(:,:,:,s1,e1), &
                            2*(vector % N+1)*vector % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcID, globalSideId, &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)
          ELSE
            msgCount = msgCount + 1
            CALL MPI_IRECV(vector % extBoundary % hostData(:,:,:,s1,e1), &
                           2*(vector % N+1)*vector % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcId, globalSideId,  &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

            msgCount = msgCount +1
            CALL MPI_ISEND(vector % boundary % hostData(:,:,:,s1,e1), &
                            2*(vector % N+1)*vector % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcID, globalSideId, &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)
          ENDIF

        ENDIF

      ENDDO
    ENDDO

    mpiHandler % msgCount = msgCount
#endif

  END SUBROUTINE MPIExchangeAsync_MappedVector2D

  SUBROUTINE ApplyFlip_MappedVector2D(mpiHandler,selfSideInfo,vector,gpuAccel)
  ! Apply side flips to sides where MPI exchanges took place.
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(hfInt32_r3), INTENT(in) :: selfSideInfo
    TYPE(MappedVector2D), INTENT(inout) :: vector
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1, s1, e2, s2
    INTEGER :: i, i2
    INTEGER :: globalSideId, externalProcId
    REAL(prec) :: extBuff(1:2,0:vector % N)


    IF(gpuAccel)THEN

      CALL ApplyFlip_MappedVector2D_gpu_wrapper(selfSideInfo % deviceData, &
                                                mpiHandler % elemToRank % deviceData, &
                                                vector % extBoundary % deviceData, &
                                                mpiHandler % rankId, &
                                                vector % N, &
                                                vector % nVar, &
                                                vector % nElem)


    ELSE
      DO e1 = 1, vector % nElem
        DO s1 = 1, 4

          e2 = selfSideInfo % hostData(3,s1,e1) ! Neighbor Element
          r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

          IF(r2 /= mpiHandler % rankId)THEN

            s2 = selfSideInfo % hostData(4,s1,e1)/10
            flip = selfSideInfo % hostData(4,s1,e1) - s2*10
            globalSideId = selfSideInfo % hostdata(2,s1,e1)

            ! Need to update extBoundary with flip applied
            IF(flip == 1)THEN

              DO ivar = 1, vector % nvar
                DO i = 0, vector % N
                  i2 = vector % N - i
                  extBuff(1:2,i) = vector % extBoundary % hostData(1:2,i2,ivar,s1,e1)
                ENDDO
                DO i = 0, vector % N
                  vector % extBoundary % hostData(1:2,i,ivar,s1,e1) = extBuff(1:2,i)
                ENDDO
              ENDDO
  
            ENDIF
          ENDIF

        ENDDO
      ENDDO

  END SUBROUTINE ApplyFlip_MappedVector2D

  SUBROUTINE MPIExchangeAsync_MappedTensor2D(mpiHandler,selfSideInfo,tensor,resetCount,useDevicePtr)
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(hfInt32_r3), INTENT(in) :: selfSideInfo
    TYPE(MappedTensor2D), INTENT(inout) :: tensor
    LOGICAL, INTENT(in) :: resetCount
    LOGICAL, INTENT(in) :: useDevicePtr
    ! Local
    INTEGER :: e1, s1, e2, s2
    INTEGER :: globalSideId, externalProcId
    INTEGER :: msgCount

#ifdef MPI
    IF(resetCount)THEN
      msgCount = 0
    ELSE
      msgCount = mpiHandler % msgCount 
    ENDIF

    DO e1 = 1, scalar % nElem
      DO s1 = 1, 4

        e2 = selfSideInfo % hostData(3,s1,e1) ! Neighbor Element
        r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

        IF(r2 /= mpiHandler % rankId)THEN

          s2 = selfSideInfo % hostData(4,s1,e1)/10
          globalSideId = selfSideInfo % hostdata(2,s1,e1)

          IF(useDevicePtr)THEN
            msgCount = msgCount + 1
            CALL MPI_IRECV(tensor % extBoundary % deviceData(:,:,:,:,s1,e1), &
                           4*(tensor % N+1)*tensor % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcId, globalSideId,  &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

            msgCount = msgCount +1
            CALL MPI_ISEND(tensor % boundary % deviceData(:,:,:,:,s1,e1), &
                            4*(tensor % N+1)*tensor % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcID, globalSideId, &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)
          ELSE
            msgCount = msgCount + 1
            CALL MPI_IRECV(tensor % extBoundary % hostData(:,:,:,:,s1,e1), &
                           4*(tensor % N+1)*tensor % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcId, globalSideId,  &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

            msgCount = msgCount +1
            CALL MPI_ISEND(tensor % boundary % hostData(:,:,:,:,s1,e1), &
                            4*(tensor % N+1)*tensor % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcID, globalSideId, &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)
          ENDIF

        ENDIF

      ENDDO
    ENDDO

    mpiHandler % msgCount = msgCount
#endif

  END SUBROUTINE MPIExchangeAsync_MappedTensor2D

  SUBROUTINE ApplyFlip_MappedTensor2D(mpiHandler,selfSideInfo,tensor,gpuAccel)
  ! Apply side flips to sides where MPI exchanges took place.
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(hfInt32_r3), INTENT(in) :: selfSideInfo
    TYPE(MappedTensor2D), INTENT(inout) :: tensor
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1, s1, e2, s2
    INTEGER :: i, i2
    INTEGER :: globalSideId, externalProcId
    REAL(prec) :: extBuff(1:2,1:2,0:tensor % N)


    IF(gpuAccel)THEN

      CALL ApplyFlip_MappedTensor2D_gpu_wrapper(selfSideInfo % deviceData, &
                                                mpiHandler % elemToRank % deviceData, &
                                                tensor % extBoundary % deviceData, &
                                                mpiHandler % rankId, &
                                                tensor % N, &
                                                tensor % nVar, &
                                                tensor % nElem)


    ELSE
      DO e1 = 1, tensor % nElem
        DO s1 = 1, 4

          e2 = selfSideInfo % hostData(3,s1,e1) ! Neighbor Element
          r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

          IF(r2 /= mpiHandler % rankId)THEN

            s2 = selfSideInfo % hostData(4,s1,e1)/10
            flip = selfSideInfo % hostData(4,s1,e1) - s2*10
            globalSideId = selfSideInfo % hostdata(2,s1,e1)

            ! Need to update extBoundary with flip applied
            IF(flip == 1)THEN

              DO ivar = 1, tensor % nvar
                DO i = 0, tensor % N
                  i2 = tensor % N - i
                  extBuff(1:2,1:2,i) = tensor % extBoundary % hostData(1:2,1:2,i2,ivar,s1,e1)
                ENDDO
                DO i = 0, tensor % N
                  tensor % extBoundary % hostData(1:2,1:2,i,ivar,s1,e1) = extBuff(1:2,1:2,i)
                ENDDO
              ENDDO
  
            ENDIF
          ENDIF

        ENDDO
      ENDDO

  END SUBROUTINE ApplyFlip_MappedTensor2D

  SUBROUTINE MPIExchangeAsync_MappedScalar3D(mpiHandler,selfSideInfo,scalar,resetCount,useDevicePtr)
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(hfInt32_r3), INTENT(in) :: selfSideInfo
    TYPE(MappedScalar3D), INTENT(inout) :: scalar
    LOGICAL, INTENT(in) :: resetCount
    LOGICAL, INTENT(in) :: useDevicePtr
    ! Local
    INTEGER :: e1, s1, e2, s2
    INTEGER :: globalSideId, externalProcId
    INTEGER :: msgCount

#ifdef MPI
    IF(resetCount)THEN
      msgCount = 0
    ELSE
      msgCount = mpiHandler % msgCount 
    ENDIF

    DO e1 = 1, scalar % nElem
      DO s1 = 1, 6

        e2 = selfSideInfo % hostData(3,s1,e1) ! Neighbor Element
        r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

        IF(r2 /= mpiHandler % rankId)THEN

          s2 = selfSideInfo % hostData(4,s1,e1)/10
          globalSideId = selfSideInfo % hostdata(2,s1,e1)

          IF(useDevicePtr)THEN
            msgCount = msgCount + 1
            CALL MPI_IRECV(scalar % extBoundary % deviceData(:,:,:,s1,e1), &
                            (scalar % N+1)*(scalar % N+1)*scalar % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcId, globalSideId,  &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

            msgCount = msgCount +1
            CALL MPI_ISEND(scalar % boundary % deviceData(:,:,:,s1,e1), &
                            (scalar % N+1)*(scalar % N+1)*scalar % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcID, globalSideId, &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)
          ELSE
            msgCount = msgCount + 1
            CALL MPI_IRECV(scalar % extBoundary % hostData(:,:,:,s1,e1), &
                            (scalar % N+1)*(scalar % N+1)*scalar % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcId, globalSideId,  &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

            msgCount = msgCount +1
            CALL MPI_ISEND(scalar % boundary % hostData(:,:,:,s1,e1), &
                            (scalar % N+1)*(scalar % N+1)*scalar % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcID, globalSideId, &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)
          ENDIF

        ENDIF

      ENDDO
    ENDDO

    mpiHandler % msgCount = msgCount
#endif

  END SUBROUTINE MPIExchangeAsync_MappedScalar3D
!
  SUBROUTINE ApplyFlip_MappedScalar3D(mpiHandler,selfSideInfo,scalar,gpuAccel)
  ! Apply side flips to sides where MPI exchanges took place.
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(hfInt32_r3), INTENT(in) :: selfSideInfo
    TYPE(MappedScalar3D), INTENT(inout) :: scalar
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1, s1, e2, s2
    INTEGER :: i, i2, j, j2
    INTEGER :: globalSideId, externalProcId
    REAL(prec) :: extBuff(0:scalar % N,0:scalar % N)


    IF(gpuAccel)THEN

      CALL ApplyFlip_MappedScalar3D_gpu_wrapper(selfSideInfo % deviceData, &
                                                mpiHandler % elemToRank % deviceData, &
                                                scalar % extBoundary % deviceData, &
                                                mpiHandler % rankId, &
                                                scalar % N, &
                                                scalar % nVar, &
                                                scalar % nElem)


    ELSE
      DO e1 = 1, scalar % nElem
        DO s1 = 1, 6

          e2 = selfSideInfo % hostData(3,s1,e1) ! Neighbor Element
          r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

          IF(r2 /= mpiHandler % rankId)THEN

            s2 = selfSideInfo % hostData(4,s1,e1)/10
            flip = selfSideInfo % hostData(4,s1,e1) - s2*10
            globalSideId = selfSideInfo % hostdata(2,s1,e1)

            ! Need to update extBoundary with flip applied
            IF(flip == 2)THEN

              DO ivar = 1, scalar % nvar
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    i2 = scalar % N - j
                    j2 = i
                    extBuff(i,j) = scalar % extBoundary % hostData(i2,j2,ivar,s1,e1)
                  ENDDO
                ENDDO
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    scalar % extBoundary % hostData(i,j,ivar,s1,e1) = extBuff(i,j)
                  ENDDO
                ENDDO
              ENDDO

            ELSEIF(flip == 3)THEN

              DO ivar = 1, scalar % nvar
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    i2 = scalar % N - i
                    j2 = scalar % N - j
                    extBuff(i,j) = scalar % extBoundary % hostData(i2,j2,ivar,s1,e1)
                  ENDDO
                ENDDO
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    scalar % extBoundary % hostData(i,j,ivar,s1,e1) = extBuff(i,j)
                  ENDDO
                ENDDO
              ENDDO
                    
            ELSEIF(flip == 4)THEN
                    
              DO ivar = 1, scalar % nvar
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    i2 = j
                    j2 = scalar % N - i
                    extBuff(i,j) = scalar % extBoundary % hostData(i2,j2,ivar,s1,e1)
                  ENDDO
                ENDDO
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    scalar % extBoundary % hostData(i,j,ivar,s1,e1) = extBuff(i,j)
                  ENDDO
                ENDDO
              ENDDO
  
            ENDIF
          ENDIF

        ENDDO
      ENDDO

  END SUBROUTINE ApplyFlip_MappedScalar3D

  SUBROUTINE MPIExchangeAsync_MappedVector3D(mpiHandler,selfSideInfo,vector,resetCount,useDevicePtr)
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(hfInt32_r3), INTENT(in) :: selfSideInfo
    TYPE(MappedVector3D), INTENT(inout) :: vector
    LOGICAL, INTENT(in) :: resetCount
    LOGICAL, INTENT(in) :: useDevicePtr
    ! Local
    INTEGER :: e1, s1, e2, s2
    INTEGER :: globalSideId, externalProcId
    INTEGER :: msgCount

#ifdef MPI
    IF(resetCount)THEN
      msgCount = 0
    ELSE
      msgCount = mpiHandler % msgCount 
    ENDIF

    DO e1 = 1, scalar % nElem
      DO s1 = 1, 6

        e2 = selfSideInfo % hostData(3,s1,e1) ! Neighbor Element
        r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

        IF(r2 /= mpiHandler % rankId)THEN

          s2 = selfSideInfo % hostData(4,s1,e1)/10
          globalSideId = selfSideInfo % hostdata(2,s1,e1)

          IF(useDevicePtr)THEN
            msgCount = msgCount + 1
            CALL MPI_IRECV(vector % extBoundary % deviceData(:,:,:,:,s1,e1), &
                            3*(vector % N+1)*(vector % N+1)*vector % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcId, globalSideId,  &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

            msgCount = msgCount +1
            CALL MPI_ISEND(vector % boundary % deviceData(:,:,:,:,s1,e1), &
                            3*(vector % N+1)*(vector % N+1)*vector % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcID, globalSideId, &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)
          ELSE
            msgCount = msgCount + 1
            CALL MPI_IRECV(vector % extBoundary % hostData(:,:,:,:,s1,e1), &
                            3*(vector % N+1)*(vector % N+1)*vector % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcId, globalSideId,  &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

            msgCount = msgCount +1
            CALL MPI_ISEND(vector % boundary % hostData(:,:,:,:,s1,e1), &
                            3*(vector % N+1)*(vector % N+1)*vector % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcID, globalSideId, &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)
          ENDIF
        ENDIF

      ENDDO
    ENDDO

    mpiHandler % msgCount = msgCount
#endif

  END SUBROUTINE MPIExchangeAsync_MappedVector3D

  SUBROUTINE ApplyFlip_MappedVector3D(mpiHandler,selfSideInfo,scalar,gpuAccel)
  ! Apply side flips to sides where MPI exchanges took place.
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(hfInt32_r3), INTENT(in) :: selfSideInfo
    TYPE(MappedVector3D), INTENT(inout) :: scalar
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1, s1, e2, s2
    INTEGER :: i, i2, j, j2
    INTEGER :: globalSideId, externalProcId
    REAL(prec) :: extBuff(1:3,0:scalar % N,0:scalar % N)


    IF(gpuAccel)THEN

      CALL ApplyFlip_MappedVector3D_gpu_wrapper(selfSideInfo % deviceData, &
                                                mpiHandler % elemToRank % deviceData, &
                                                scalar % extBoundary % deviceData, &
                                                mpiHandler % rankId, &
                                                scalar % N, &
                                                scalar % nVar, &
                                                scalar % nElem)


    ELSE
      DO e1 = 1, scalar % nElem
        DO s1 = 1, 6

          e2 = selfSideInfo % hostData(3,s1,e1) ! Neighbor Element
          r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

          IF(r2 /= mpiHandler % rankId)THEN

            s2 = selfSideInfo % hostData(4,s1,e1)/10
            flip = selfSideInfo % hostData(4,s1,e1) - s2*10
            globalSideId = selfSideInfo % hostdata(2,s1,e1)

            IF(flip == 2)THEN

              DO ivar = 1, scalar % nvar
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    i2 = scalar % N - j
                    j2 = i
                    extBuff(1:3,i,j) = scalar % extBoundary % hostData(1:3,i2,j2,ivar,s1,e1)
                  ENDDO
                ENDDO
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    scalar % extBoundary % hostData(1:3,i,j,ivar,s1,e1) = extBuff(1:3,i,j)
                  ENDDO
                ENDDO
              ENDDO

            ELSEIF(flip == 3)THEN

              DO ivar = 1, scalar % nvar
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    i2 = scalar % N - i
                    j2 = scalar % N - j
                    extBuff(1:3,i,j) = scalar % extBoundary % hostData(1:3,i2,j2,ivar,s1,e1)
                  ENDDO
                ENDDO
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    scalar % extBoundary % hostData(1:3,i,j,ivar,s1,e1) = extBuff(1:3,i,j)
                  ENDDO
                ENDDO
              ENDDO
                    
            ELSEIF(flip == 4)THEN
                    
              DO ivar = 1, scalar % nvar
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    i2 = j
                    j2 = scalar % N - i
                    extBuff(1:3,i,j) = scalar % extBoundary % hostData(1:3,i2,j2,ivar,s1,e1)
                  ENDDO
                ENDDO
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    scalar % extBoundary % hostData(1:3,i,j,ivar,s1,e1) = extBuff(1:3,i,j)
                  ENDDO
                ENDDO
              ENDDO
  
            ENDIF
          ENDIF
        ENDDO
      ENDDO

  END SUBROUTINE ApplyFlip_MappedVector3D

  SUBROUTINE MPIExchangeAsync_MappedTensor3D(mpiHandler,selfSideInfo,tensor,resetCount,useDevicePtr)
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(hfInt32_r3), INTENT(in) :: selfSideInfo
    TYPE(MappedTensor3D), INTENT(inout) :: tensor
    LOGICAL, INTENT(in) :: resetCount
    LOGICAL, INTENT(in) :: useDevicePtr
    ! Local
    INTEGER :: e1, s1, e2, s2
    INTEGER :: globalSideId, externalProcId
    INTEGER :: msgCount

#ifdef MPI
    IF(resetCount)THEN
      msgCount = 0
    ELSE
      msgCount = mpiHandler % msgCount 
    ENDIF

    DO e1 = 1, scalar % nElem
      DO s1 = 1, 6

        e2 = selfSideInfo % hostData(3,s1,e1) ! Neighbor Element
        r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

        IF(r2 /= mpiHandler % rankId)THEN

          s2 = selfSideInfo % hostData(4,s1,e1)/10
          globalSideId = selfSideInfo % hostdata(2,s1,e1)

          IF(useDevicePtr)THEN
            msgCount = msgCount + 1
            CALL MPI_IRECV(tensor % extBoundary % deviceData(:,:,:,:,:,s1,e1), &
                            9*(tensor % N +1)*(tensor % N+1)*tensor % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcId, globalSideId,  &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

            msgCount = msgCount +1
            CALL MPI_ISEND(tensor % boundary % deviceData(:,:,:,:,:,s1,e1), &
                            9*(tensor % N+1)*(tensor % N+1)*tensor % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcID, globalSideId, &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)
          ELSE
            msgCount = msgCount + 1
            CALL MPI_IRECV(tensor % extBoundary % hostData(:,:,:,:,:,s1,e1), &
                            9*(tensor % N +1)*(tensor % N+1)*tensor % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcId, globalSideId,  &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

            msgCount = msgCount +1
            CALL MPI_ISEND(tensor % boundary % hostData(:,:,:,:,:,s1,e1), &
                            9*(tensor % N+1)*(tensor % N+1)*tensor % nVar, &
                            mpiHandler % mpiPrec, &
                            externalProcID, globalSideId, &
                            mpiHandler % mpiComm, &
                            mpiHandler % requests(msgCount,1), iError)

          ENDIF

        ENDIF

      ENDDO
    ENDDO

    mpiHandler % msgCount = msgCount
#endif

  END SUBROUTINE MPIExchangeAsync_MappedTensor3D

  SUBROUTINE ApplyFlip_MappedTensor3D(mpiHandler,selfSideInfo,scalar,gpuAccel)
  ! Apply side flips to sides where MPI exchanges took place.
    IMPLICIT NONE
    CLASS(MPILayer), INTENT(inout) :: mpiHandler
    TYPE(hfInt32_r3), INTENT(in) :: selfSideInfo
    TYPE(MappedTensor3D), INTENT(inout) :: scalar
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1, s1, e2, s2
    INTEGER :: i, i2, j, j2
    INTEGER :: globalSideId, externalProcId
    REAL(prec) :: extBuff(1:3,1:3,0:scalar % N,0:scalar % N)


    IF(gpuAccel)THEN

      CALL ApplyFlip_MappedTensor3D_gpu_wrapper(selfSideInfo % deviceData, &
                                                mpiHandler % elemToRank % deviceData, &
                                                scalar % extBoundary % deviceData, &
                                                mpiHandler % rankId, &
                                                scalar % N, &
                                                scalar % nVar, &
                                                scalar % nElem)


    ELSE
      DO e1 = 1, scalar % nElem
        DO s1 = 1, 6

          e2 = selfSideInfo % hostData(3,s1,e1) ! Neighbor Element
          r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

          IF(r2 /= mpiHandler % rankId)THEN

            s2 = selfSideInfo % hostData(4,s1,e1)/10
            flip = selfSideInfo % hostData(4,s1,e1) - s2*10
            globalSideId = selfSideInfo % hostdata(2,s1,e1)

            IF(flip == 2)THEN

              DO ivar = 1, scalar % nvar
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    i2 = scalar % N - j
                    j2 = i
                    extBuff(1:3,1:3,i,j) = scalar % extBoundary % hostData(1:3,1:3,i2,j2,ivar,s1,e1)
                  ENDDO
                ENDDO
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    scalar % extBoundary % hostData(1:3,1:3,i,j,ivar,s1,e1) = extBuff(1:3,1:3,i,j)
                  ENDDO
                ENDDO
              ENDDO

            ELSEIF(flip == 3)THEN

              DO ivar = 1, scalar % nvar
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    i2 = scalar % N - i
                    j2 = scalar % N - j
                    extBuff(1:3,1:3,i,j) = scalar % extBoundary % hostData(1:3,1:3,i2,j2,ivar,s1,e1)
                  ENDDO
                ENDDO
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    scalar % extBoundary % hostData(1:3,1:3,i,j,ivar,s1,e1) = extBuff(1:3,1:3,i,j)
                  ENDDO
                ENDDO
              ENDDO
                    
            ELSEIF(flip == 4)THEN
                    
              DO ivar = 1, scalar % nvar
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    i2 = j
                    j2 = scalar % N - i
                    extBuff(1:3,1:3,i,j) = scalar % extBoundary % hostData(1:3,1:3,i2,j2,ivar,s1,e1)
                  ENDDO
                ENDDO
                DO j = 0, scalar % N
                  DO i = 0, scalar % N
                    scalar % extBoundary % hostData(1:3,1:3,i,j,ivar,s1,e1) = extBuff(1:3,1:3,i,j)
                  ENDDO
                ENDDO
              ENDDO
  
            ENDIF
          ENDIF

        ENDDO
      ENDDO

  END SUBROUTINE ApplyFlip_MappedTensor3D

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
