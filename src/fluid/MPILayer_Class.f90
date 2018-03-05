! MPILayer_Class.f90
!
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE MPILayer_Class

  USE ModelPrecision
  USE NodalDGSolution_3D_Class
  USE Faces_Class
  USE BoundaryCommunicator_Class

#ifdef HAVE_CUDA
  USE cudafor
#endif


  IMPLICIT NONE


  TYPE MPILayer

    ! For the unstructured mesh, I opt for building my own DATA structure that bundles messages between neighboring
    ! ranks, rather than attempting to build an MPI DATA structure. IF a structured mesh is USEd, one optimization
    ! would be to USE MPI-DATA TYPEs to handle the message passing.

    INTEGER :: nVars, N

    INTEGER, ALLOCATABLE :: requestHandle(:)
    INTEGER, ALLOCATABLE :: requestStats(:,:)

#ifdef HAVE_CUDA
    REAL(prec), DEVICE, ALLOCATABLE :: sendBuffer_dev(:,:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: recvBuffer_dev(:,:,:,:,:)

    REAL(prec), PINNED, ALLOCATABLE :: sendBuffer(:,:,:,:,:)
    REAL(prec), PINNED, ALLOCATABLE :: recvBuffer(:,:,:,:,:)

#else
    REAL(prec), ALLOCATABLE :: sendBuffer(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: recvBuffer(:,:,:,:,:)
#endif

    CONTAINS

      PROCEDURE :: Build => Build_MPILayer
      PROCEDURE :: Trash => Trash_MPILayer
      PROCEDURE :: MPI_Exchange
      PROCEDURE :: Finalize_MPI_Exchange

  END TYPE MPILayer

CONTAINS

  SUBROUTINE Build_MPILayer( myMPI, extComm, N, nVars )

    IMPLICIT NONE
    CLASS( MPILayer ), INTENT(out)           :: myMPI
    TYPE( BoundaryCommunicator ), INTENT(in) :: extComm
    INTEGER, INTENT(in)                      :: N, nVars


    myMPI % N     = N
    myMPI % nVars = nVars

    ALLOCATE( myMPI % requestHandle(1:extComm % nNeighbors*2), &
      myMPI % requestStats(MPI_STATUS_SIZE,1:extComm % nNeighbors*2), &
      myMPI % recvBuffer(0:N, 0:N, 1:nVars, 1:extComm % maxBufferSize, 1:extComm % nNeighbors), &
      myMPI % sendBuffer(0:N, 0:N, 1:nVars, 1:extComm % maxBufferSize, 1:extComm % nNeighbors) )

    myMPI % requestHandle = 0
    myMPI % requestStats  = 0
    myMPI % recvBuffer    = 0.0_prec
    myMPI % sendBuffer    = 0.0_prec


#ifdef HAVE_CUDA

    ALLOCATE( myMPI % recvBuffer_dev(0:N, 0:N, 1:nVars, 1:extComm % maxBufferSize, 1:extComm % nNeighbors), &
              myMPI % sendBuffer_dev(0:N, 0:N, 1:nVars, 1:extComm % maxBufferSize, 1:extComm % nNeighbors) )
    myMPI % recvBuffer_dev = 0.0_prec
    myMPI % sendBuffer_dev = 0.0_prec

#endif


  END SUBROUTINE Build_MPILayer
!
  SUBROUTINE Trash_MPILayer( myMPI )

    IMPLICIT NONE
    CLASS( MPILayer ), INTENT(inout) :: myMPI

    PRINT*, '    S/R Trash_MPILayer : Clearing memory.'

    DEALLOCATE( myMPI % requestHandle, &
                myMPI % requestStats, &
                myMPI % recvBuffer, &
                myMPI % sendBuffer )

#ifdef HAVE_CUDA


    DEALLOCATE( myMPI % recvBuffer_dev, myMPI % sendBuffer_dev )

#endif


  END SUBROUTINE Trash_MPILayer

  SUBROUTINE MPI_Exchange( myMPI, state, meshFaces, extComm )

    IMPLICIT NONE
    CLASS( MPILayer ), INTENT(inout)          :: myMPI
    TYPE( NodalDGSolution_3D ), INTENT(inout) :: state
    TYPE( Faces ), INTENT(in)                 :: meshFaces
    TYPE( BoundaryCommunicator ), INTENT(in)  :: extComm
    ! Local
    INTEGER    :: iNeighbor, iError
#ifdef HAVE_CUDA
    TYPE(dim3) :: grid, tBlock


    tBlock = dim3( state % N+1, state % N+1, 1 )
    grid = dim3( state % nEquations, extComm % nBoundaries, 1 )

    CALL BoundaryToBuffer_CUDAKernel<<<grid, tBlock>>>( myMPI % sendBuffer_dev, &
                                                        state % boundarySolution_dev, &
                                                        meshFaces % elementIDs_dev, &
                                                        meshFaces % elementSides_dev, &
                                                        extComm % boundaryIDs_dev, &
                                                        extComm % extProcIDs_dev, &
                                                        extComm % rankTable_dev,&
                                                        extComm % bufferMap_dev,&
                                                        meshFaces % nFaces, extComm % nBoundaries, &
                                                        extComm % nProc, extComm % myRank, &
                                                        state % N, state % nEquations,&
                                                        extComm % nNeighbors, extComm % maxBufferSize, &
                                                        state % nElements )

#ifdef GPU_DIRECT
    iError = cudaDeviceSynchronize( )
    DO iNeighbor = 1, extComm % nNeighbors


      CALL MPI_IRECV( myMPI % recvBuffer_dev(:,:,:,:,iNeighbor), &
        (myMPI % N+1)*(myMPI % N+1)*myMPI % nVars*extComm % bufferSize(iNeighbor), &
        extComm % MPI_PREC,   &
        extComm % neighborRank(iNeighbor), 0,  &
        extComm % MPI_COMM,   &
        myMPI % requestHandle((iNeighbor-1)*2+1), iError )

      CALL MPI_ISEND( myMPI % sendBuffer_dev(:,:,:,:,iNeighbor), &
        (myMPI % N+1)*(myMPI % N+1)*myMPI % nVars*extComm % bufferSize(iNeighbor), &
        extComm % MPI_PREC, &
        extComm % neighborRank(iNeighbor), 0, &
        extComm % MPI_COMM, &
        myMPI % requestHandle(iNeighbor*2), iError)

    ENDDO

#else
    myMPI % sendBuffer= myMPI % sendBuffer_dev
    iError = cudaDeviceSynchronize( )

    DO iNeighbor = 1, extComm % nNeighbors

      CALL MPI_IRECV( myMPI % recvBuffer(:,:,:,:,iNeighbor), &
        (myMPI % N+1)*(myMPI % N+1)*myMPI % nVars*extComm % bufferSize(iNeighbor), &
        extComm % MPI_PREC,   &
        extComm % neighborRank(iNeighbor), 0,  &
        extComm % MPI_COMM,   &
        myMPI % requestHandle(iNeighbor*2-1), iError )

      CALL MPI_ISEND( myMPI % sendBuffer(:,:,:,:,iNeighbor), &
        (myMPI % N+1)*(myMPI % N+1)*myMPI % nVars*extComm % bufferSize(iNeighbor), &
        extComm % MPI_PREC, &
        extComm % neighborRank(iNeighbor), 0, &
        extComm % MPI_COMM, &
        myMPI % requestHandle(iNeighbor*2), iError)

    ENDDO


#endif

#else
    INTEGER    :: IFace, bID
    INTEGER    :: tag
    INTEGER    :: e1, e2, s1, p2

    DO bID = 1, extComm % nBoundaries

      IFace     = extComm % boundaryIDs( bID )
      p2        = extComm % extProcIDs(bID)

      IF( p2 /= extComm % myRank )THEN

        e1        = meshFaces % elementIDs(1,IFace)
        s1        = meshFaces % elementSides(1,IFace)
        iNeighbor = extComm % rankTable(p2)

        myMPI % sendBuffer(:,:,:,extComm % bufferMap(bID), iNeighbor ) = state % boundarySolution(:,:,:,s1,e1)

      ENDIF
    ENDDO

    DO iNeighbor = 1, extComm % nNeighbors

      CALL MPI_IRECV( myMPI % recvBuffer(:,:,:,:,iNeighbor), &
        (myMPI % N+1)*(myMPI % N+1)*myMPI % nVars*extComm % bufferSize(iNeighbor), &
        extComm % MPI_PREC,   &
        extComm % neighborRank(iNeighbor), 0,  &
        extComm % MPI_COMM,   &
        myMPI % requestHandle(iNeighbor*2-1), iError )

      CALL MPI_ISEND( myMPI % sendBuffer(:,:,:,:,iNeighbor), &
        (myMPI % N+1)*(myMPI % N+1)*myMPI % nVars*extComm % bufferSize(iNeighbor), &
        extComm % MPI_PREC, &
        extComm % neighborRank(iNeighbor), 0, &
        extComm % MPI_COMM, &
        myMPI % requestHandle(iNeighbor*2), iError)

    ENDDO
#endif


  END SUBROUTINE MPI_Exchange
!
  SUBROUTINE Finalize_MPI_Exchange( myMPI, state, meshFaces, extComm )

    IMPLICIT NONE
    CLASS( MPILayer ), INTENT(inout)          :: myMPI
    TYPE( NodalDGSolution_3D ), INTENT(inout) :: state
    TYPE( Faces ), INTENT(in)                 :: meshFaces
    TYPE( BoundaryCommunicator ), INTENT(in)  :: extComm
    ! Local
    INTEGER    :: IFace, bID
    INTEGER    :: tag, ierror
    INTEGER    :: e1, e2, s1, p2, iNeighbor, jUnpack
    INTEGER    :: fUnit
#ifdef HAVE_CUDA
    TYPE(dim3) :: grid, tBlock
#endif

    CALL MPI_WaitAll( extComm % nNeighbors*2, &
                      myMPI % requestHandle, &
                      myMPI % requestStats, &
                      iError)

#ifdef HAVE_CUDA

#ifndef GPU_DIRECT
    myMPI % recvBuffer_dev = myMPI % recvBuffer
#endif
    tBlock = dim3(myMPI % N+1, &
                  myMPI % N+1, &
                  1 )
    grid = dim3( myMPI % nVars, extComm % nBoundaries,1)

    CALL BufferToBoundary_CUDAKernel<<<grid, tBlock>>>( myMPI % recvBuffer_dev, &
      state % externalState_dev, &
      extComm % boundaryIDs_dev, &
      extComm % extProcIDs_dev, &
      extComm % rankTable_dev, &
      extComm % unPackMap_dev, &
      meshFaces % nFaces, &
      extComm % nBoundaries, &
      extComm % nProc, extComm % myRank, &
      myMPI % N, myMPI % nVars, extComm % nNeighbors, &
      extComm % maxBufferSize )

#else

    DO bID = 1, extComm % nBoundaries

      p2        = extComm % extProcIDs(bID)

      IF( p2 /= extComm % myRank )THEN

        iNeighbor = extComm % rankTable(p2)
        jUnpack   = extComm % unpackMap(bID)

        state % externalState(:,:,:,bID) = myMPI % recvBuffer(:,:,:,jUnpack,iNeighbor)

      ENDIF

    ENDDO

#endif


  END SUBROUTINE Finalize_MPI_Exchange
!
#ifdef HAVE_CUDA
  ATTRIBUTES(Global) SUBROUTINE BoundaryToBuffer_CUDAKernel( sendBuffer, boundarySolution, faceToElement, faceToSide, boundaryToFaceID, &
    boundaryToProcID, rankToNeighbor, boundaryToBuffer, nFaces, nBoundaries, nRanks, myRank, N, numEq, nNeighbors, &
    bufferSize, nElements )

    IMPLICIT NONE
    INTEGER, VALUE, INTENT(in)        :: nFaces, nBoundaries, nRanks, myRank
    INTEGER, VALUE, INTENT(in)        :: N, numEq, nNeighbors, nElements, bufferSize
    REAL(prec), DEVICE, INTENT(out)   :: sendBuffer(0:N,0:N,1:numEq,1:bufferSize,1:nNeighbors)
    REAL(prec), DEVICE, INTENT(in)    :: boundarySolution(0:N,0:N,1:numEq,1:6,1:nElements)
    INTEGER, DEVICE, INTENT(in)       :: faceToElement(1:2,1:nFaces)
    INTEGER, DEVICE, INTENT(in)       :: faceToSide(1:2,1:nFaces)
    INTEGER, DEVICE, INTENT(in)       :: boundaryToFaceID(1:nBoundaries)
    INTEGER, DEVICE, INTENT(in)       :: boundaryToProcID(1:nBoundaries)
    INTEGER, DEVICE, INTENT(in)       :: rankToNeighbor(0:nRanks-1)
    INTEGER, DEVICE, INTENT(in)       :: boundaryToBuffer(1:nBoundaries)
    ! Local
    INTEGER :: bID, i, j, iEq
    INTEGER :: IFace, p2, neighborID, bufferID, elementID, sideID

    iEq = blockIdx % x
    bID = blockIdx % y

    i   = threadIdx % x - 1
    j   = threadIdx % y - 1

    p2 = boundaryToProcID(bID)

  IF( p2 /= myRank )THEN

    neighborID = rankToNeighbor(p2)
    bufferID   = boundaryToBuffer(bID)
    iFace      = boundaryToFaceID(bID)
    elementID  = faceToElement(1,iFace)
    sideID     = faceToSide(1,iFace)

    sendBuffer(i,j,iEq,bufferID,neighborID) = boundarySolution(i,j,iEq,sideID,elementID)
  ENDIF

  END SUBROUTINE BoundaryToBuffer_CUDAKernel
!
ATTRIBUTES(Global) SUBROUTINE BufferToBoundary_CUDAKernel( recvBuffer, externalSolution, boundaryToFaceID, &
  boundaryToProcID, rankTable, unPackMap, nFaces, nBoundaries, nRanks, myRank, N, numEq, nNeighbors, &
  bufferSize )
IMPLICIT NONE
INTEGER, VALUE, INTENT(in)        :: nFaces, nBoundaries, nRanks, myRank
INTEGER, VALUE, INTENT(in)        :: N, numEq, nNeighbors, bufferSize
INTEGER, DEVICE, INTENT(in)       :: boundaryToFaceID(1:nBoundaries)
INTEGER, DEVICE, INTENT(in)       :: boundaryToProcID(1:nBoundaries)
INTEGER, DEVICE, INTENT(in)       :: rankTable(0:nRanks-1)
INTEGER, DEVICE, INTENT(in)       :: unPackMap(1:nBoundaries)
REAL(prec), DEVICE, INTENT(in)    :: recvBuffer(0:N,0:N,1:numEq,1:bufferSize,1:nNeighbors)
REAL(prec), DEVICE, INTENT(inout) :: externalSolution(0:N,0:N,1:numEq,1:nBoundaries)
 ! Local
INTEGER :: bID, i, j, iEq
INTEGER :: p2

iEq = blockIdx % x
bID = blockIdx % y

i   = threadIdx % x - 1
j   = threadIdx % y - 1

p2        = boundaryToProcID(bID)


IF( p2 /= myRank )THEN
  externalSolution(i,j,iEq,bID) = recvBuffer(i,j,iEq,unpackMap(bID),rankTable(p2))
ENDIF


END SUBROUTINE BufferToBoundary_CUDAKernel

#endif

END MODULE MPILayer_Class


