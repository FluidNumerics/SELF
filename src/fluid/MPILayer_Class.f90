! MPILayer_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE MPI_Layer_Class

USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
USE FluidParams_Class
USE Lagrange_Class
USE dgStorage_Class
USE RollOffFilter_Class
USE NodalDGSolution_3D_Class
USE HexMesh_Class
USE BoundaryCommunicator_Class

#ifdef HAVE_CUDA
USE cudafor
#else


IMPLICIT NONE

INCLUDE 'mpif.h'

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

    END TYPE CommunicationTable

    TYPE MPILayer

       ! For the unstructured mesh, I opt for building my own data structure that bundles messages between neighboring
       ! ranks, rather than attempting to build an MPI data structure. If a structured mesh is used, one optimization 
       ! would be to use MPI-data types to handle the message passing.

       INTEGER :: nEq

       INTEGER, ALLOCATABLE :: requestHandle(:)
       INTEGER, ALLOCATABLE :: requestStats(:,:)

#ifdef HAVE_CUDA
       INTEGER, DEVICE, ALLOCATABLE :: nState_dev, nStress_dev, nSGS_dev
       INTEGER, DEVICE, ALLOCATABLE :: myRank_dev, nProc_dev, nNeighbors_dev
       REAL(prec), DEVICE, ALLOCATABLE :: sendBuffer_dev(:,:,:,:,:) ! (0:N,0:N,1:nEq,1:nSharedFaces)
       REAL(prec), DEVICE, ALLOCATABLE :: recvBuffer_dev(:,:,:,:,:)

       REAL(prec), PINNED, ALLOCATABLE :: sendBuffer(:,:,:,:,:) ! (0:N,0:N,1:nEq,1:nSharedFaces)
       REAL(prec), PINNED, ALLOCATABLE :: recvBuffer(:,:,:,:,:)

#else
       REAL(prec), ALLOCATABLE :: sendBuffer(:,:,:,:,:) ! (0:N,0:N,1:nEq,1:nSharedFaces)
       REAL(prec), ALLOCATABLE :: recvBuffer(:,:,:,:,:)
#endif

    END TYPE MPILayer

CONTAINS

SUBROUTINE Build_MPILayer( myMPI, params, nState, nStress, nSGS, setupSuccess )

   IMPLICIT NONE
   CLASS( MPILayer ), INTENT(out)  :: myMPI
   TYPE( ModelParams ), INTENT(in) :: params
   LOGICAL, INTENT(inout)          :: setupSuccess


      myMPI % myRank = 0
      myMPI % nProc  = 1
      myMPI % MPI_COMM = MPI_COMM_WORLD

      IF( prec == sp )THEN
         myMPI % MPI_PREC = MPI_FLOAT
      ELSE
         myMPI % MPI_PREC = MPI_DOUBLE
      ENDIF

      CALL MPI_INIT( myMPI % mpiErr )
      CALL MPI_COMM_RANK( myMPI % MPI_COMM, myMPI % myRank, mpiErr )
      CALL MPI_COMM_SIZE( myMPI % MPI_COMM, myMPI % nProc, mpiErr )

      PRINT*, '    S/R Build_MPILayer : Greetings from Process ', myMPI % myRank+1, ' of ', myMPI % nProc

      CALL myMPI % ConstructCommTables(  )

      ALLOCATE( stateReqHandle(1:myMPI % nNeighbors*2), &
                stressReqHandle(1:myMPI % nNeighbors*2), &
                SGSReqHandle(1:myMPI % nNeighbors*2), &
                stateStats(MPI_STATUS_SIZE,1:myMPI % nNeighbors*2), &
                stressStats(MPI_STATUS_SIZE,1:myMPI % nNeighbors*2), &
                SGSStats(MPI_STATUS_SIZE,1:myMPI % nNeighbors*2) )

      myMPI % externalSGS = params % viscosity

#ifdef HAVE_CUDA

      ALLOCATE( myMPI % nState_dev, &
                myMPI % nStress_dev, &
                myMPI % nSGS_dev, &
                myMPI % myRank_dev, &
                myMPI % nProc_dev, &
                myMPI % nNeighbors_dev )

       myMPI % nState_dev      = nState
       myMPI % nStress_dev     = nStress
       myMPI % nSGS_dev        = nSGS
       myMPI % myRank_dev      = myMPI % myRank
       myMPI % nProc_dev       = myMPI % nProc
       myMPI % nNeighbors_dev  = myMPI % nNeighbors
       myMPI % externalSGS_dev = params % viscosity

#endif


 END SUBROUTINE Build_MPILayer
!
  SUBROUTINE Trash_MPILayer( myMPI )

   IMPLICIT NONE
   CLASS( MPILayer ), INTENT(inout) :: myMPI
   
      PRINT*, 'S/R Trash_MPILayer : Clearing memory.'

      DEALLOCATE( myMPI % neighborRank, & 
                  myMPI % bufferSize, & 
                  myMPI % bufferMap, & 
                  myMPI % sendStateBuffer, & 
                  myMPI % recvStateBuffer, &  
                  myMPI % sendStressBuffer, &  
                  myMPI % recvStressBuffer, &  
                  myMPI % sendSGSBuffer, &  
                  myMPI % recvSGSBuffer, &  
                  myMPI % rankTable, &  
                  myMPI % stateReqHandle, &
                  myMPI % stressReqHandle, &
                  myMPI % SGSReqHandle, &
                  myMPI % stateStats, &
                  myMPI % stressStats, &
                  myMPI % SGSStats )
#ifdef HAVE_CUDA

      DEALLOCATE( myMPI % sendStateBuffer_dev, & 
                  myMPI % recvStateBuffer_dev, & 
                  myMPI % sendStressBuffer_dev, & 
                  myMPI % recvStressBuffer_dev, & 
                  myMPI % sendSGSBuffer_dev, & 
                  myMPI % recvSGSBuffer_dev, & 
                  myMPI % neighborRank_dev, &
                  myMPI % bufferSize_dev, &
                  myMPI % bufferMap_dev, &
                  myMPI % rankTable_dev )

#endif

     CALL MPI_FINALIZE( myMPI % mpiErr )
      
 END SUBROUTINE Trash_MPILayer
!
 SUBROUTINE ConstructCommTables( myMPI, extComm, N, nEq )
   IMPLICIT NONE
   CLASS( myMPI ), INTENT(inout)               :: myMPI
   TYPE( BoundaryCommunicator ), INTENT(inout) :: extComm
   INTEGER, INTENT(in)                         :: N, nEq
   ! Local
   INTEGER, ALLOCATABLE :: bufferCounter(:)
   INTEGER :: sharedFaceCount(0:myMPI % nProc-1)
   INTEGER :: iFace, bID, iNeighbor
   INTEGER :: tag, ierror
   INTEGER :: e1, e2, s1, p2, nmsg, maxFaceCount
   INTEGER :: fUnit


      ALLOCATE( myMPI % rankTable(0:myMPI % nProc-1) )

      ! Count up the number of neighboring ranks
      myMPI % rankTable = 0
      sharedFaceCount     = 0
      DO bID = 1, extComm % nBoundaries

         p2 = extComm % extProcIDS(bID)

         IF( p2 /= myMPI % myRank )THEN

            mpiLayer % rankTable(p2) = 1
            sharedFaceCount(p2) = sharedFaceCount(p2)+1

         ENDIF

      ENDDO


      myMPI % nNeighbors = SUM( myMPI % rankTable )
      PRINT*, '  S/R ConstructCommTables : Found', myMPI % nNeighbors, 'neighbors for Rank', myMPI % myRank+1
      
      ALLOCATE( myMPI % neighborRank(1:myMPI % nNeighbors), &
                myMPI % bufferSize(1:myMPI % nNeighbors), &
                bufferCounter(1:myMPI % nNeighbors) )


      ! For each neighbor, set the neighbor's rank
      iNeighbor = 0
      DO p2 = 0, myMPI % nProc-1

         IF( myMPI % rankTable(p2) == 1 )THEN

            iNeighbor = iNeighbor + 1
            myMPI % neighborRank(iNeighbor) = p2
            myMPI % rankTable(p2) = iNeighbor

         ENDIF

      ENDDO

      
      maxFaceCount = MAXVAL( sharedFaceCount )
      DO iNeighbor = 1, myMPI % nNeighbors

         p2 = myMPI % neighborRank(iNeighbor)
         myMPI % bufferSize(iNeighbor) = sharedFaceCount(p2)

      ENDDO


      myMPI % maxBufferSize = maxFaceCount
         
      ALLOCATE( myMPI % recvStateBuffer(0:N, 0:N, 1:nEq, 1:maxFaceCount, 1:myMPI % nNeighbors), &
                myMPI % sendStateBuffer(0:N, 0:N, 1:nEq, 1:maxFaceCount, 1:myMPI % nNeighbors), &
                myMPI % recvStressBuffer(0:N, 0:N, 1:(nEq-1)*3, 1:maxFaceCount, 1:myMPI % nNeighbors), &
                myMPI % sendStressBuffer(0:N, 0:N, 1:(nEq-1)*3, 1:maxFaceCount, 1:myMPI % nNeighbors), &
                myMPI % sendSGSBuffer(0:N, 0:N, 1:(nEq-1), 1:maxFaceCount, 1:myMPI % nNeighbors), &
                myMPI % recvSGSBuffer(0:N, 0:N, 1:(nEq-1), 1:maxFaceCount, 1:myMPI % nNeighbors) )
                
      myMPI % recvStateBuffer  = 0.0_prec
      myMPI % sendStateBuffer  = 0.0_prec
      myMPI % recvStressBuffer = 0.0_prec
      myMPI % sendStressBuffer = 0.0_prec
      myMPI % recvSGSBuffer    = 0.0_prec
      myMPI % sendSGSBuffer    = 0.0_prec

      bufferCounter = 0
       
      ALLOCATE( myMPI % bufferMap(1:extComm % nBoundaries) )


      myMPI % bufferMap = 0

      DO bID = 1, extComm % nBoundaries
      
         p2 = extComm % extProcIDs(bID)

         ! In the event that the external process ID (p2) is identical to the current rank (p1),
         ! then this boundary edge involves a physical boundary condition and does not require a 
         ! message exchange

         IF( p2 /= myMPI % myRank )THEN 

            iNeighbor = myMPI % rankTable(p2)

            bufferCounter(iNeighbor) = bufferCounter(iNeighbor) + 1
            myMPI % bufferMap(bID)   = bufferCounter(iNeighbor)   

         ENDIF

      ENDDO

      DEALLOCATE( bufferCounter )

#ifdef HAVE_CUDA
      ALLOCATE( myMPI % rankTable_dev(0:myMPI % nProc-1), &
                myMPI % neighborRank_dev(1:myMPI % nNeighbors), &
                myMPI % bufferSize_dev(1:myMPI % nNeighbors), &  
                myMPI % recvStateBuffer_dev(0:N, 0:N, 1:nEq, 1:maxFaceCount, 1:myMPI % nNeighbors), &
                myMPI % sendStateBuffer_dev(0:N, 0:N, 1:nEq, 1:maxFaceCount, 1:myMPI % nNeighbors), &
                myMPI % recvStressBuffer_dev(0:N, 0:N, 1:(nEq-1)*3, 1:maxFaceCount, 1:myMPI % nNeighbors), &
                myMPI % sendStressBuffer_dev(0:N, 0:N, 1:(nEq-1)*3, 1:maxFaceCount, 1:myMPI % nNeighbors), &
                myMPI % sendSGSBuffer_dev(0:N, 0:N, 1:(nEq-1), 1:maxFaceCount, 1:myMPI % nNeighbors), &
                myMPI % recvSGSBuffer_dev(0:N, 0:N, 1:(nEq-1), 1:maxFaceCount, 1:myMPI % nNeighbors), &  
                myMPI % bufferMap_dev(1:extComm % nBoundaries) )

      myDGSEM % mpiPackets % rankTable_dev    = myDGSEM % mpiPackets % rankTable
      myDGSEM % mpiPackets % neighborRank_dev = myDGSEM % mpiPackets % neighborRank
      myDGSEM % mpiPackets % bufferSize_dev   = myDGSEM % mpiPackets % bufferSize
      myDGSEM % mpiPackets % bufferMap_dev    = myDGSEM % mpiPackets % bufferMap
#endif

 END SUBROUTINE ConstructCommTables

 SUBROUTINE MPI_State_Exchange( myMPI, fluidState, meshFaces, extComm ) 

   IMPLICIT NONE
   CLASS( MPILayer ), INTENT(inout)       :: myMPI
   TYPE( NodalDGSolution_3D ), INTENT(in) :: fluidState
   TYPE( Faces ), INTENT(in)              :: meshFaces
   ! Local
   INTEGER    :: iNeighbor
#ifdef HAVE_CUDA
   TYPE(dim3) :: grid, tBlock

  
      tBlock = dim3(4*(ceiling( REAL(N+1)/4 ) ), &
                    4*(ceiling( REAL(N+1)/4 ) ) , &
                    nEq )
      grid = dim3(extComm % nBoundaries,1,1) 

      CALL BoundaryToBuffer_CUDAKernel<<<grid, tBlock>>>( myMPI % sendStateBuffer_dev, &
                                                          fluidState % boundarySolution_dev, &
                                                          meshFaces % elementIDs_dev, &
                                                          meshFaces % elementSides_dev, &
                                                          extComm % boundaryIDs_dev, &
                                                          extComm % extProcIDs_dev, &
                                                          mpiLayer % rankTable_dev,&
                                                          mpiLayer % bufferMap_dev,&
                                                          meshFaces % nFaces_dev, extComm % nBoundaries_dev, &
                                                          myMPI % nProc_dev, myMPI % myRank_dev, &
                                                          fluidState % N_dev, fluidState % nEquations_dev,&
                                                          myMPI % nNeighbors_dev, myMPI % maxBufferSize_dev, &
                                                          fluidState % nElements_dev )


#ifdef CUDA_DIRECT
      iError = cudaDeviceSynchronize( )
      DO iNeighbor = 1, myMPI % nNeighbors 

            
            CALL MPI_IRECV( myMPI % recvStateBuffer_dev(:,:,:,:,iNeighbor), & 
                           (fluidState % N+1)*(fluidState % N+1)*fluidState % nEq*myMPI % bufferSize(iNeighbor), &                
                           myMPI % MPI_PREC,   &                      
                           myMPI % mpiPackets % neighborRank(iNeighbor), 0,  &                        
                           myMPI % MPI_COMM,   &                
                           stateReqHandle((iNeighbor-1)*2+1), iError )           

            CALL MPI_ISEND( myMPI % sendStateBuffer_dev(:,:,:,:,iNeighbor), & 
                           (fluidState % N+1)*(fluidState % N+1)*fluidState % nEq*myMPI % bufferSize(iNeighbor), &                
                           myMPI % MPI_PREC, &      
                           myMPI % neighborRank(iNeighbor), 0, &       
                           myMPI % MPI_COMM, &
                           stateReqHandle(iNeighbor*2), iError)  
                           
      ENDDO

#else

      myMPI % sendStateBuffer= myMPI % sendStateBuffer_dev
      iError = cudaDeviceSynchronize( )

      DO iNeighbor = 1, myMPI % nNeighbors 
            
            CALL MPI_IRECV( myMPI % recvStateBuffer(:,:,:,:,iNeighbor), & 
                           (fluidState % N+1)*(fluidState % N+1)*fluidState % nEq*myMPI % bufferSize(iNeighbor), &                
                           myMPI % MPI_PREC,   &                      
                           myMPI % mpiPackets % neighborRank(iNeighbor), 0,  &                        
                           myMPI % MPI_COMM,   &                
                           myMPI % stateReqHandle((iNeighbor-1)*2+1), iError )           

            CALL MPI_ISEND( myMPI % sendStateBuffer(:,:,:,:,iNeighbor), & 
                           (fluidState % N+1)*(fluidState % N+1)*fluidState % nEq*myMPI % bufferSize(iNeighbor), &                
                           myMPI % MPI_PREC, &      
                           myMPI % neighborRank(iNeighbor), 0, &       
                           myMPI % MPI_COMM, &
                           myMPI % stateReqHandle(iNeighbor*2), iError)  
                         
      ENDDO
      
      
#endif

#else
   INTEGER    :: iFace, bID
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, iNeighbor

      DO bID = 1, extComm % nBoundaries
      
         iFace     = extComm % boundaryIDs( bID )
         p2        = extComm % extProcIDs(bID)
         
         IF( p2 /= myMPI % myRank )THEN 

            e1        = meshFaces % elementIDs(1,iFace)
            s1        = meshFaces % elementSides(1,iFace)
            iNeighbor = myMPI % rankTable(p2)
         
            myMPI % sendStateBuffer(:,:,:,myMPI % bufferMap(bID), iNeighbor ) = boundaryState(:,:,:,s1,e1) 

         ENDIF
      ENDDO

      DO iNeighbor = 1, myMPI % nNeighbors 
            
            CALL MPI_IRECV( myMPI % recvStateBuffer(:,:,:,:,iNeighbor), & 
                           (fluidState % N+1)*(fluidState % N+1)*fluidState % nEq*myMPI % bufferSize(iNeighbor), &                
                           myMPI % MPI_PREC,   &                      
                           myMPI % mpiPackets % neighborRank(iNeighbor), 0,  &                        
                           myMPI % MPI_COMM,   &                
                           myMPI % stateReqHandle((iNeighbor-1)*2+1), iError )           

            CALL MPI_ISEND( myMPI % sendStateBuffer(:,:,:,:,iNeighbor), & 
                           (fluidState % N+1)*(fluidState % N+1)*fluidState % nEq*myMPI % bufferSize(iNeighbor), &                
                           myMPI % MPI_PREC, &      
                           myMPI % neighborRank(iNeighbor), 0, &       
                           myMPI % MPI_COMM, &
                           myMPI % stateReqHandle(iNeighbor*2), iError)  
                           
      ENDDO
#endif


 END SUBROUTINE MPIState_Exchange
!
 SUBROUTINE Finalize_MPIState_Exchange( myMPI, boundaryConditions, extComm ) 

   IMPLICIT NONE
   CLASS( MPILayer ), INTENT(inout)               :: myMPI
   TYPE( FluidBoundaryConditions ), INTENT(inout) :: boundaryConditions
   REAL
   ! Local
   INTEGER    :: iFace, bID
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, iNeighbor, jUnpack
   INTEGER    :: fUnit
#ifdef HAVE_CUDA
   TYPE(dim3) :: grid, tBlock
#endif

      CALL MPI_WaitAll(myMPI % nNeighbors*2,stateReqHandle,stateStats,iError)

#ifdef HAVE_CUDA
  
#ifndef CUDA_DIRECT
      myMPI % recvStateBuffer_dev = myMPI % recvStateBuffer
#endif
      tBlock = dim3(4*(ceiling( REAL(myMPI % N+1)/4 ) ), &
                    4*(ceiling( REAL(myMPI % N+1)/4 ) ) , &
                    myMPI % nState )
      grid = dim3(extComm % nBoundaries,1,1) 

      CALL BufferToBoundary_CUDAKernel<<<grid, tBlock>>>( myMPI % recvStateBuffer_dev, &
                                                          boundaryConditions % externalState_dev, &
                                                          extComm % boundaryIDs_dev, &
                                                          extComm % extProcIDs_dev, &
                                                          myMPI % rankTable_dev, &
                                                          extComm % unPackMap_dev, &
                                                          meshFaces % nFaces_dec, &
                                                          extComm % nBoundaries_dev, &
                                                          myMPI % nProc_dev, myMPI % myRank_dev, 7
                                                          my % N, nEq, myDGSEM % nNeighbors, &
                                                          myDGSEM % mpiPackets % maxBufferSize, myDGSEM % mesh % nElems )

#else

      DO bID = 1, myDGSEM % extComm % nBoundaries

         p2        = myDGSEM % extComm % extProcIDs(bID)
         
         IF( p2 /= myDGSEM % myRank )THEN 
      
            iNeighbor = myDGSEM % mpiPackets % rankTable(p2)
            jUnpack   = myDGSEM % extComm % unpackMap(bID)

            myDGSEM % externalState(:,:,:,bID) = myDGSEM % mpiPackets % recvStateBuffer(:,:,:,jUnpack,iNeighbor)
         ENDIF

      ENDDO

#endif


 END SUBROUTINE FinalizeMPI_StateExchange_Fluid
!
 SUBROUTINE MPI_StressExchange_Fluid( myDGSEM ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iFace, bID
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, iNeighbor, jUnpack
#ifdef HAVE_CUDA
   TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                    (nEq-1)*3 )
      grid = dim3(myDGSEM % extComm % nBoundaries,1,1) 

      CALL BoundaryToBuffer_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mpiPackets % sendStressBuffer_dev, &
                                                          myDGSEM % stressTensor % boundarySolution_dev, &
                                                          myDGSEM % mesh % faces_dev % elementIDs, &
                                                          myDGSEM % mesh % faces_dev % elementSides, &
                                                          myDGSEM % extComm % boundaryIDs_dev, &
                                                          myDGSEM % extComm % extProcIDs_dev, &
                                                          myDGSEM % mpiPackets % rankTable_dev,&
                                                          myDGSEM % mpiPackets % bufferMap_dev,&
                                                          myDGSEM % mesh % nFaces, myDGSEM % extComm % nBoundaries, &
                                                          myDGSEM % nProc, myDGSEM % myRank, myDGSEM % N, (nEq-1)*3,&
                                                          myDGSEM % nNeighbors, myDGSEM % mpiPackets % maxBufferSize, &
                                                          myDGSEM % mesh % nElems )

#ifdef CUDA_DIRECT
      iError = cudaDeviceSynchronize( )
      DO iNeighbor = 1, myDGSEM % nNeighbors 
            
            CALL MPI_IRECV( myDGSEM % mpiPackets % recvStressBuffer_dev(:,:,:,:,iNeighbor), & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*3*myDGSEM % mpiPackets % bufferSize(iNeighbor), &                
                           MPI_PREC,   &                      
                           myDGSEM % mpiPackets % neighborRank(iNeighbor), 0,  &                        
                           MPI_COMM_WORLD,   &                
                           stressReqHandle((iNeighbor-1)*2+1), iError )           

            CALL MPI_ISEND( myDGSEM % mpiPackets % sendStressBuffer_dev(:,:,:,:,iNeighbor), & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*3*myDGSEM % mpiPackets % bufferSize(iNeighbor), &       
                           MPI_PREC, &      
                           myDGSEM % mpiPackets % neighborRank(iNeighbor), 0, &       
                           MPI_COMM_WORLD, &
                           stressReqHandle(iNeighbor*2), iError)  
                           
      ENDDO

#else

      myDGSEM % mpiPackets % sendStressBuffer= myDGSEM % mpiPackets % sendStressBuffer_dev
      iError = cudaDeviceSynchronize( )
      
      DO iNeighbor = 1, myDGSEM % nNeighbors 
            
            CALL MPI_IRECV( myDGSEM % mpiPackets % recvStressBuffer(:,:,:,:,iNeighbor), & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*3*myDGSEM % mpiPackets % bufferSize(iNeighbor), &                
                           MPI_PREC,   &                      
                           myDGSEM % mpiPackets % neighborRank(iNeighbor), 0,  &                        
                           MPI_COMM_WORLD,   &                
                           stressReqHandle((iNeighbor-1)*2+1), iError )           

            CALL MPI_ISEND( myDGSEM % mpiPackets % sendStressBuffer(:,:,:,:,iNeighbor), & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*3*myDGSEM % mpiPackets % bufferSize(iNeighbor), &       
                           MPI_PREC, &      
                           myDGSEM % mpiPackets % neighborRank(iNeighbor), 0, &       
                           MPI_COMM_WORLD, &
                           stressReqHandle(iNeighbor*2), iError)  
                           
      ENDDO
#endif

#else
      DO bID = 1, myDGSEM % extComm % nBoundaries
      
         iFace     = myDGSEM % extComm % boundaryIDs( bID )
         p2        = myDGSEM % extComm % extProcIDs(bID)
         
         IF( p2 /= myDGSEM % myRank )THEN 

            e1        = myDGSEM % mesh % Faces(iFace) % elementIDs(1)
            s1        = myDGSEM % mesh % Faces(iFace) % elementSides(1)
            iNeighbor = myDGSEM % mpiPackets % rankTable(p2)
         
            myDGSEM % mpiPackets % sendStressBuffer(:,:,:,myDGSEM % mpiPackets % bufferMap(bID), iNeighbor ) =&
               myDGSEM % stressTensor % boundarySolution(:,:,:,s1,e1) 

         ENDIF
      ENDDO
      DO iNeighbor = 1, myDGSEM % nNeighbors 
            
            CALL MPI_IRECV( myDGSEM % mpiPackets % recvStressBuffer(:,:,:,:,iNeighbor), & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*3*myDGSEM % mpiPackets % bufferSize(iNeighbor), &                
                           MPI_PREC,   &                      
                           myDGSEM % mpiPackets % neighborRank(iNeighbor), 0,  &                        
                           MPI_COMM_WORLD,   &                
                           stressReqHandle((iNeighbor-1)*2+1), iError )           

            CALL MPI_ISEND( myDGSEM % mpiPackets % sendStressBuffer(:,:,:,:,iNeighbor), & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*3*myDGSEM % mpiPackets % bufferSize(iNeighbor), &       
                           MPI_PREC, &      
                           myDGSEM % mpiPackets % neighborRank(iNeighbor), 0, &       
                           MPI_COMM_WORLD, &
                           stressReqHandle(iNeighbor*2), iError)  
                           
      ENDDO
#endif
              


 END SUBROUTINE MPI_StressExchange_Fluid
!
 SUBROUTINE FinalizeMPI_StressExchange_Fluid( myDGSEM ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iFace, bID
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, iNeighbor, jUnpack
#ifdef HAVE_CUDA
   TYPE(dim3) :: grid, tBlock
   REAL(prec) :: t2, t1
#endif

      CALL MPI_WaitAll(myDGSEM % nNeighbors*2,stressReqHandle,stressStats,iError)

#ifdef HAVE_CUDA
  
#ifndef CUDA_DIRECT
      myDGSEM % mpiPackets % recvStressBuffer_dev = myDGSEM % mpiPackets % recvStressBuffer
#endif
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                   4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                   (nEq-1)*3 )
      grid = dim3(myDGSEM % extComm % nBoundaries,1,1) 

      CALL BufferToBoundary_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mpiPackets % recvStressBuffer_dev, &
                                                          myDGSEM % externalStress_dev, &
                                                          myDGSEM % extComm % boundaryIDs_dev, &
                                                          myDGSEM % extComm % extProcIDs_dev, &
                                                          myDGSEM % mpiPackets % rankTable_dev, &
                                                          myDGSEM % extComm % unPackMap_dev, &
                                                          myDGSEM % mesh % nFaces, myDGSEM % extComm % nBoundaries, &
                                                          myDGSEM % nProc, myDGSEM % myRank, myDGSEM % N, (nEq-1)*3, myDGSEM % nNeighbors, &
                                                          myDGSEM % mpiPackets % maxBufferSize, myDGSEM % mesh % nElems )

#else

      DO bID = 1, myDGSEM % extComm % nBoundaries

         p2        = myDGSEM % extComm % extProcIDs(bID)
         
         IF( p2 /= myDGSEM % myRank )THEN 
      
            iNeighbor = myDGSEM % mpiPackets % rankTable(p2)
            jUnpack   = myDGSEM % extComm % unpackMap(bID)

            myDGSEM % externalStress(:,:,:,bID) = myDGSEM % mpiPackets % recvStressBuffer(:,:,:,jUnpack,iNeighbor)
         ENDIF

      ENDDO

#endif

 END SUBROUTINE FinalizeMPI_StressExchange_Fluid
!
 SUBROUTINE MPI_SGSExchange_Fluid( myDGSEM ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iFace, bID
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, iNeighbor, jUnpack
#ifdef HAVE_CUDA
   TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                   4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                   nEq-1 )
      grid = dim3(myDGSEM % extComm % nBoundaries,1,1) 

      CALL BoundaryToBuffer_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mpiPackets % sendSGSBuffer_dev, &
                                                          myDGSEM % sgsCoeffs % boundarySolution_dev, &
                                                          myDGSEM % mesh % faces_dev % elementIDs, &
                                                          myDGSEM % mesh % faces_dev % elementSides, &
                                                          myDGSEM % extComm % boundaryIDs_dev, &
                                                          myDGSEM % extComm % extProcIDs_dev, &
                                                          myDGSEM % mpiPackets % rankTable_dev,&
                                                          myDGSEM % mpiPackets % bufferMap_dev,&
                                                          myDGSEM % mesh % nFaces, myDGSEM % extComm % nBoundaries, &
                                                          myDGSEM % nProc, myDGSEM % myRank, myDGSEM % N, nEq-1,&
                                                          myDGSEM % nNeighbors, myDGSEM % mpiPackets % maxBufferSize, &
                                                          myDGSEM % mesh % nElems )

#ifdef CUDA_DIRECT
      iError = cudaDeviceSynchronize( )
      DO iNeighbor = 1, myDGSEM % nNeighbors 
            
            CALL MPI_IRECV( myDGSEM % mpiPackets % recvSGSBuffer_dev(:,:,:,:,iNeighbor), & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*myDGSEM % mpiPackets % bufferSize(iNeighbor), &                
                           MPI_PREC,   &                      
                           myDGSEM % mpiPackets % neighborRank(iNeighbor), 0,  &                        
                           MPI_COMM_WORLD,   &                
                           sgsReqHandle((iNeighbor-1)*2+1), iError )           

            CALL MPI_ISEND( myDGSEM % mpiPackets % sendSGSBuffer_dev(:,:,:,:,iNeighbor), & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*myDGSEM % mpiPackets % bufferSize(iNeighbor), &       
                           MPI_PREC, &      
                           myDGSEM % mpiPackets % neighborRank(iNeighbor), 0, &       
                           MPI_COMM_WORLD, &
                           sgsReqHandle(iNeighbor*2), iError)  
                           
      ENDDO

#else

      myDGSEM % mpiPackets % sendSGSBuffer= myDGSEM % mpiPackets % sendSGSBuffer_dev
      iError = cudaDeviceSynchronize( )
      DO iNeighbor = 1, myDGSEM % nNeighbors 
            
            CALL MPI_IRECV( myDGSEM % mpiPackets % recvSGSBuffer(:,:,:,:,iNeighbor), & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*myDGSEM % mpiPackets % bufferSize(iNeighbor), &                
                           MPI_PREC,   &                      
                           myDGSEM % mpiPackets % neighborRank(iNeighbor), 0,  &                        
                           MPI_COMM_WORLD,   &                
                           SGSReqHandle((iNeighbor-1)*2+1), iError )           

            CALL MPI_ISEND( myDGSEM % mpiPackets % sendSGSBuffer(:,:,:,:,iNeighbor), & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*myDGSEM % mpiPackets % bufferSize(iNeighbor), &       
                           MPI_PREC, &      
                           myDGSEM % mpiPackets % neighborRank(iNeighbor), 0, &       
                           MPI_COMM_WORLD, &
                           SGSReqHandle(iNeighbor*2), iError)  
                           
      ENDDO
#endif

#else
      DO bID = 1, myDGSEM % extComm % nBoundaries
      
         iFace     = myDGSEM % extComm % boundaryIDs( bID )
         p2        = myDGSEM % extComm % extProcIDs(bID)
         
         IF( p2 /= myDGSEM % myRank )THEN 

            e1        = myDGSEM % mesh % Faces(iFace) % elementIDs(1)
            s1        = myDGSEM % mesh % Faces(iFace) % elementSides(1)
            iNeighbor = myDGSEM % mpiPackets % rankTable(p2)
         
            myDGSEM % mpiPackets % sendSGSBuffer(:,:,:,myDGSEM % mpiPackets % bufferMap(bID),iNeighbor ) =&
               myDGSEM % sgsCoeffs % boundarySolution(:,:,:,s1,e1) 

         ENDIF
      ENDDO

      DO iNeighbor = 1, myDGSEM % nNeighbors 
            
            CALL MPI_IRECV( myDGSEM % mpiPackets % recvSGSBuffer(:,:,:,:,iNeighbor), & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*myDGSEM % mpiPackets % bufferSize(iNeighbor), &                
                           MPI_PREC,   &                      
                           myDGSEM % mpiPackets % neighborRank(iNeighbor), 0,  &                        
                           MPI_COMM_WORLD,   &                
                           SGSReqHandle((iNeighbor-1)*2+1), iError )           

            CALL MPI_ISEND( myDGSEM % mpiPackets % sendSGSBuffer(:,:,:,:,iNeighbor), & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*myDGSEM % mpiPackets % bufferSize(iNeighbor), &       
                           MPI_PREC, &      
                           myDGSEM % mpiPackets % neighborRank(iNeighbor), 0, &       
                           MPI_COMM_WORLD, &
                           SGSReqHandle(iNeighbor*2), iError)  
                           
      ENDDO
#endif

   
 END SUBROUTINE MPI_SGSExchange_Fluid
!
 SUBROUTINE FinalizeMPI_SGSExchange_Fluid( myDGSEM ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: iFace, bID
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, iNeighbor, jUnpack
#ifdef HAVE_CUDA
   TYPE(dim3) :: grid, tBlock
   REAL(prec) :: t2, t1
#endif

      CALL MPI_WaitAll(myDGSEM % nNeighbors*2,SGSReqHandle,sgsStats,iError)

#ifdef HAVE_CUDA
  
#ifndef CUDA_DIRECT
      myDGSEM % mpiPackets % recvSGSBuffer_dev = myDGSEM % mpiPackets % recvSGSBuffer
#endif
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                   4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                   nEq-1 )
      grid = dim3(myDGSEM % extComm % nBoundaries,1,1) 

      CALL BufferToBoundary_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mpiPackets % recvSGSBuffer_dev, &
                                                          myDGSEM % externalSGS_dev, &
                                                          myDGSEM % extComm % boundaryIDs_dev, &
                                                          myDGSEM % extComm % extProcIDs_dev, &
                                                          myDGSEM % mpiPackets % rankTable_dev, &
                                                          myDGSEM % extComm % unPackMap_dev, &
                                                          myDGSEM % mesh % nFaces, myDGSEM % extComm % nBoundaries, &
                                                          myDGSEM % nProc, myDGSEM % myRank, myDGSEM % N, nEq-1, myDGSEM % nNeighbors, &
                                                          myDGSEM % mpiPackets % maxBufferSize, myDGSEM % mesh % nElems )

#else

      DO bID = 1, myDGSEM % extComm % nBoundaries

         p2        = myDGSEM % extComm % extProcIDs(bID)
         
         IF( p2 /= myDGSEM % myRank )THEN 
      
            iNeighbor = myDGSEM % mpiPackets % rankTable(p2)
            jUnpack   = myDGSEM % extComm % unpackMap(bID)

            myDGSEM % externalSGS(:,:,:,bID) = myDGSEM % mpiPackets % recvSGSBuffer(:,:,:,jUnpack,iNeighbor)

         ENDIF

      ENDDO

#endif


   
 END SUBROUTINE FinalizeMPI_SGSExchange_Fluid

#endif
!
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!                                                                                                 !
!  This section of code contains routines for computing the fluxes through the element faces      !
!                                                                                                 !
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!
 SUBROUTINE InternalFaceFlux_Fluid( myDGSEM )

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                    1 )
      grid = dim3(myDGSEM % mesh % nFaces,1,1)  
      
      CALL InternalFaceFlux_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mesh % faces_dev % elementIDs, &
                                                  myDGSEM % mesh % faces_dev % elementSides, &
                                                  myDGSEM % mesh % faces_dev % boundaryID, &
                                                  myDGSEM % mesh % faces_dev % iMap, &
                                                  myDGSEM % mesh % faces_dev % jMap, &
                                                  myDGSEM % mesh % geom_dev % nHat_dev, &
                                                  myDGSEM % state % boundarySolution_dev, &
                                                  myDGSEM % static % boundarySolution_dev, &
                                                  myDGSEM % externalState_dev, &
                                                  myDGSEM % state % boundaryFlux_dev, &
                                                  myDGSEM % stressTensor % boundaryFlux_dev )
#else
   ! Local
   INTEGER :: iEl, iFace
   INTEGER    :: i, j, k, m, iEq, jEq
   INTEGER    :: ii, jj, bID
   INTEGER    :: e1, s1, e2, s2
   REAL(prec) :: nHat(1:3), norm
   REAL(prec) :: uOut, uIn, cIn, cOut, T
   REAL(prec) :: jump(1:nEq-1), aS(1:nEq-1)
   REAL(prec) :: fac, hCapRatio, rC

      hCapRatio = ( myDGSEM % params % R + myDGSEM % params % Cv ) / myDGSEM % params % Cv
      rC        =   myDGSEM % params % R / ( myDGSEM % params % R + myDGSEM % params % Cv )

      !$OMP DO PRIVATE( jump, aS )
      DO iFace = 1, myDGSEM % mesh % nFaces

         
         e1 = myDGSEM % mesh % faces(iFace) % elementIDs(1)
         s1 = myDGSEM % mesh % faces(iFace) % elementSides(1)
         e2 = myDGSEM % mesh % faces(iFace) % elementIDs(2)
         s2 = ABS(myDGSEM % mesh % faces(iFace) % elementSides(2))

         IF( e2 > 0 )THEN

            DO j = 0, myDGSEM % N  
               DO i = 0, myDGSEM % N
   
                  IF( i == 0 )THEN
                     IF( j == 0 )THEN
                        ii = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % iStart) + &
                             (myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % jStart)
                        jj = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % jStart) + &
                             (myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % iStart)
                     ELSE
                        ii = myDGSEM % mesh % faces(iFace) % swapDimensions*&
                             (ii+myDGSEM % mesh % faces(iFace) % jInc) + &
                             (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             myDGSEM % mesh % faces(iFace) % iStart
                        jj = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (jj+myDGSEM % mesh % faces(iFace) % jInc) +&
                             myDGSEM % mesh % faces(iFace) % swapDimensions*&
                             myDGSEM % mesh % faces(iFace) % jStart
                     ENDIF
                  ELSE
                     ii = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                          (ii + myDGSEM % mesh % faces(iFace) % iInc) +&
                          myDGSEM % mesh % faces(iFace) % swapDimensions*ii
                     jj = myDGSEM % mesh % faces(iFace) % swapDimensions*&
                          (jj+myDGSEM % mesh % faces(iFace) % iInc) + &
                          (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*jj
                  ENDIF
   
                  norm = sqrt( myDGSEM % mesh % geom(e1) % nHat(1,i,j,s1)*myDGSEM % mesh % geom(e1) % nHat(1,i,j,s1) + &
                               myDGSEM % mesh % geom(e1) % nHat(2,i,j,s1)*myDGSEM % mesh % geom(e1) % nHat(2,i,j,s1) + &
                               myDGSEM % mesh % geom(e1) % nHat(3,i,j,s1)*myDGSEM % mesh % geom(e1) % nHat(3,i,j,s1) )
   
                  DO k = 1, 3
                     nHat(k) = myDGSEM % mesh % geom(e1) % nHat(k,i,j,s1)/norm
                  ENDDO
               
                  DO iEq = 1, nEq-1
                  jump(iEq)  = myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2) - &
                               myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) !outState - inState
                  ENDDO
                            
                  T =   (myDGSEM % static % boundarySolution(ii,jj,5,s2,e2) + myDGSEM % state % boundarySolution(ii,jj,5,s2,e2))/&
                        (myDGSEM % static % boundarySolution(ii,jj,4,s2,e2) + myDGSEM % state % boundarySolution(ii,jj,4,s2,e2))
                          
                  ! Sound speed estimate for the external and internal states
                  cOut = sqrt( myDGSEM % params % R *T* &
                              ( (myDGSEM % state % boundarySolution(ii,jj,6,s2,e2)+&
                                 myDGSEM % static % boundarySolution(ii,jj,6,s2,e2))/&
                                 myDGSEM % params % P0 )**rC   )
                        
                  T =   (myDGSEM % static % boundarySolution(i,j,5,s1,e1) + &
                         myDGSEM % state % boundarySolution(i,j,5,s1,e1))/&
                        (myDGSEM % static % boundarySolution(i,j,4,s1,e1) + &
                         myDGSEM % state % boundarySolution(i,j,4,s1,e1) )        
                             
                  cIn  = sqrt( myDGSEM % params % R*T* &
                              ( (myDGSEM % state % boundarySolution(i,j,6,s1,e1)+&
                                 myDGSEM % static % boundarySolution(i,j,6,s1,e1))/&
                                 myDGSEM % params % P0 )**rC  )
                               
                  ! External normal velocity component

                  uOut = ( myDGSEM % state % boundarySolution(ii,jj,1,s2,e2)*nHat(1) + &
                           myDGSEM % state % boundarySolution(ii,jj,2,s2,e2)*nHat(2) + &
                           myDGSEM % state % boundarySolution(ii,jj,3,s2,e2)*nHat(3) )/& 
                         ( myDGSEM % state % boundarySolution(ii,jj,4,s2,e2) + &
                           myDGSEM % static % boundarySolution(ii,jj,4,s2,e2) )
                           
                  ! Internal normal velocity component
                  uIn  = ( myDGSEM % state % boundarySolution(i,j,1,s1,e1)*nHat(1) + &
                           myDGSEM % state % boundarySolution(i,j,2,s1,e1)*nHat(2) + &
                           myDGSEM % state % boundarySolution(i,j,3,s1,e1)*nHat(3) )/& 
                         ( myDGSEM % state % boundarySolution(i,j,4,s1,e1) + &
                           myDGSEM % static % boundarySolution(i,j,4,s1,e1) ) 
                         
                  ! Lax-Friedrich's estimate of the magnitude of the flux jacobian matrix
                  fac = max( abs(uIn+cIn), abs(uIn-cIn), abs(uOut+cOut), abs(uOut-cOut) )
                  !fac = max( abs(uIn), abs(uOut) )

                  ! Advective flux
                  DO iEq = 1, nEq-1
                        aS(iEq) = uIn*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) + &
                                        myDGSEM % static % boundarySolution(i,j,iEq,s1,e1) ) +&
                                 uOut*( myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2) + &
                                        myDGSEM % static % boundarySolution(ii,jj,iEq,s2,e2) )
                  ENDDO
                  
                  DO k = 1, 3
                  ! Momentum flux due to pressure
                  aS(k) = aS(k) + (myDGSEM % state % boundarySolution(i,j,6,s1,e1) + &
                                   myDGSEM % state % boundarySolution(ii,jj,6,s2,e2))*nHat(k)
                  ENDDO    
      
                 
                  DO iEq = 1, nEq-1
                     myDGSEM % state % boundaryFlux(i,j,iEq,s1,e1) =  0.5_prec*( aS(iEq) - fac*jump(iEq) )*norm
                     myDGSEM % state % boundaryFlux(ii,jj,iEq,s2,e2) = -myDGSEM % state % boundaryFlux(i,j,iEq,s1,e1)
                     IF( iEq == 4 )THEN
                        DO k = 1, 3
                           jEq = k+(iEq-1)*3
                           ! Calculate the LDG flux for the stress tensor.
                           myDGSEM % stressTensor % boundaryFlux(i,j,jEq,s1,e1) = 0.5_prec*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) +&
                                                                                             myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2) )*&
                                                                                             myDGSEM % mesh % geom(e1) % nHat(k,i,j,s1)
                                                                                        
                           myDGSEM % stressTensor % boundaryFlux(ii,jj,jEq,s2,e2) = -myDGSEM % stressTensor % boundaryFlux(i,j,jEq,s1,e1)
                        ENDDO
                     ELSE
                        DO k = 1, 3
                           jEq = k+(iEq-1)*3
                           ! Calculate the LDG flux for the stress tensor.
                           myDGSEM % stressTensor % boundaryFlux(i,j,jEq,s1,e1) = 0.5_prec*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1)/&
                                                                                             (myDGSEM % state % boundarySolution(i,j,4,s1,e1) +&
                                                                                              myDGSEM % static % boundarySolution(i,j,4,s1,e1))+&
                                                                                             myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2)/&
                                                                                             (myDGSEM % state % boundarySolution(ii,jj,4,s2,e2) +&
                                                                                              myDGSEM % static % boundarySolution(ii,jj,4,s2,e2)) )*&
                                                                                             myDGSEM % mesh % geom(e1) % nHat(k,i,j,s1)
                                                                                        
                           myDGSEM % stressTensor % boundaryFlux(ii,jj,jEq,s2,e2) = -myDGSEM % stressTensor % boundaryFlux(i,j,jEq,s1,e1)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
             
         ENDIF 
         
      ENDDO 
      !$OMP ENDDO
      

#endif

 END SUBROUTINE InternalFaceFlux_Fluid
!
 SUBROUTINE BoundaryFaceFlux_Fluid( myDGSEM )

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                    1 )
      grid = dim3(myDGSEM % mesh % nFaces,1,1)  
      
      CALL BoundaryFaceFlux_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mesh % faces_dev % elementIDs, &
                                                  myDGSEM % mesh % faces_dev % elementSides, &
                                                  myDGSEM % mesh % faces_dev % boundaryID, &
                                                  myDGSEM % mesh % faces_dev % iMap, &
                                                  myDGSEM % mesh % faces_dev % jMap, &
                                                  myDGSEM % mesh % geom_dev % nHat_dev, &
                                                  myDGSEM % state % boundarySolution_dev, &
                                                  myDGSEM % static % boundarySolution_dev, &
                                                  myDGSEM % externalState_dev, &
                                                  myDGSEM % state % boundaryFlux_dev, &
                                                  myDGSEM % stressTensor % boundaryFlux_dev )
#else
   ! Local
   INTEGER :: iEl, iFace
   INTEGER    :: i, j, k, m, iEq, jEq
   INTEGER    :: ii, jj, bID
   INTEGER    :: e1, s1, e2, s2
   REAL(prec) :: nHat(1:3), norm
   REAL(prec) :: uOut, uIn, cIn, cOut, T
   REAL(prec) :: jump(1:nEq-1), aS(1:nEq-1)
   REAL(prec) :: fac, hCapRatio, rC

      hCapRatio = ( myDGSEM % params % R + myDGSEM % params % Cv ) / myDGSEM % params % Cv
      rC        =   myDGSEM % params % R / ( myDGSEM % params % R + myDGSEM % params % Cv )

      !$OMP DO PRIVATE( jump, aS )
      DO iFace = 1, myDGSEM % mesh % nFaces

         
         e1 = myDGSEM % mesh % faces(iFace) % elementIDs(1)
         s1 = myDGSEM % mesh % faces(iFace) % elementSides(1)
         e2 = myDGSEM % mesh % faces(iFace) % elementIDs(2)
         s2 = ABS(myDGSEM % mesh % faces(iFace) % elementSides(2))

         IF( e2 < 0 )THEN

            DO j = 0, myDGSEM % N  
               DO i = 0, myDGSEM % N
   
                  IF( i == 0 )THEN
                     IF( j == 0 )THEN
                        ii = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % iStart) + &
                             (myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % jStart)
                        jj = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % jStart) + &
                             (myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % iStart)
                     ELSE
                        ii = myDGSEM % mesh % faces(iFace) % swapDimensions*&
                             (ii+myDGSEM % mesh % faces(iFace) % jInc) + &
                             (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             myDGSEM % mesh % faces(iFace) % iStart
                        jj = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (jj+myDGSEM % mesh % faces(iFace) % jInc) +&
                             myDGSEM % mesh % faces(iFace) % swapDimensions*&
                             myDGSEM % mesh % faces(iFace) % jStart
                     ENDIF
                  ELSE
                     ii = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                          (ii + myDGSEM % mesh % faces(iFace) % iInc) +&
                          myDGSEM % mesh % faces(iFace) % swapDimensions*ii
                     jj = myDGSEM % mesh % faces(iFace) % swapDimensions*&
                          (jj+myDGSEM % mesh % faces(iFace) % iInc) + &
                          (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*jj
                  ENDIF
   
                  norm = sqrt( myDGSEM % mesh % geom(e1) % nHat(1,i,j,s1)*myDGSEM % mesh % geom(e1) % nHat(1,i,j,s1) + &
                               myDGSEM % mesh % geom(e1) % nHat(2,i,j,s1)*myDGSEM % mesh % geom(e1) % nHat(2,i,j,s1) + &
                               myDGSEM % mesh % geom(e1) % nHat(3,i,j,s1)*myDGSEM % mesh % geom(e1) % nHat(3,i,j,s1) )
   
                  DO k = 1, 3
                     nHat(k) = myDGSEM % mesh % geom(e1) % nHat(k,i,j,s1)/norm
                  ENDDO
 
                  bID  = ABS(myDGSEM % mesh % faces(iFace) % boundaryID)
                  DO iEq = 1, nEq-1                 
                  jump(iEq)  = myDGSEM % externalState(ii,jj,iEq,bID) - &
                               myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) !outState - inState
                  ENDDO

                  ! Sound speed estimate for the external and internal states
                  
                  T = (myDGSEM % static % boundarySolution(i,j,5,s1,e1)+myDGSEM % externalState(ii,jj,5,bID))/&
                      (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % externalState(ii,jj,4,bID))
                           
                  cOut = sqrt( myDGSEM % params % R*T* &
                              ( (myDGSEM % externalState(ii,jj,6,bID)+&
                                 myDGSEM % static % boundarySolution(i,j,6,s1,e1) )/&
                                 myDGSEM % params % P0 )**rC   )
                   
                  T = (myDGSEM % static % boundarySolution(i,j,5,s1,e1)+myDGSEM % state % boundarySolution(i,j,5,s1,e1))/&
                      (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % boundarySolution(i,j,4,s1,e1))  
                                      
                  cIn  = sqrt( myDGSEM % params % R*T* &
                              ( (myDGSEM % state % boundarySolution(i,j,6,s1,e1)+&
                                 myDGSEM % static % boundarySolution(i,j,6,s1,e1) )/&
                                 myDGSEM % params % P0 )**rC  )
                               
                  ! External normal velocity component
                  uOut = ( myDGSEM % externalState(ii,jj,1,bID)*nHat(1) + &
                           myDGSEM % externalState(ii,jj,2,bID)*nHat(2) + &
                           myDGSEM % externalState(ii,jj,3,bID)*nHat(3) )/&
                           (myDGSEM % externalState(ii,jj,4,bID)+&
                            myDGSEM % static % boundarySolution(i,j,4,s1,e1))
                            
                  ! Internal normal velocity component
                  uIn  = ( myDGSEM % state % boundarySolution(i,j,1,s1,e1)*nHat(1) + &
                           myDGSEM % state % boundarySolution(i,j,2,s1,e1)*nHat(2) + &
                           myDGSEM % state % boundarySolution(i,j,3,s1,e1)*nHat(3) )/& 
                           (myDGSEM % state % boundarySolution(i,j,4,s1,e1)+&
                            myDGSEM % static % boundarySolution(i,j,4,s1,e1) )
                            
                  fac = max( abs(uIn+cIn), abs(uIn-cIn), abs(uOut+cOut), abs(uOut-cOut) )
                 ! fac = max( abs(uIn), abs(uOut) )
                  
                  ! Advective flux
                  DO iEq = 1, nEq-1
                        aS(iEq) = uIn*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) + myDGSEM % static % boundarySolution(i,j,iEq,s1,e1) ) +&
                                 uOut*( myDGSEM % externalState(ii,jj,iEq,bID) + myDGSEM % static % boundarySolution(i,j,iEq,s1,e1) )
                  ENDDO
                  
                  DO k = 1, 3
                  ! Momentum flux due to pressure
                  aS(k) = aS(k) + (myDGSEM % state % boundarySolution(i,j,6,s1,e1) + &
                                   myDGSEM % externalState(ii,jj,6,bID))*nHat(k)
                  ENDDO  
                  
                  
                  DO iEq = 1, nEq-1
                     myDGSEM % state % boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*( aS(iEq) - fac*jump(iEq) )*norm
                     IF( iEq == 4 )THEN
                        DO k = 1, 3
                           jEq = k+(iEq-1)*3
                           ! Calculate the Bassi-Rebay flux for the stress tensor.
                           myDGSEM % stressTensor % boundaryFlux(i,j,jEq,s1,e1) = 0.5_prec*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) +&
                                                                                             myDGSEM % externalState(ii,jj,iEq,bID) )*&
                                                                                             myDGSEM % mesh % geom(e1) % nHat(k,i,j,s1)
                        ENDDO
                     ELSE
                        DO k = 1, 3
                           jEq = k+(iEq-1)*3
                           ! Calculate the Bassi-Rebay flux for the stress tensor.
                           myDGSEM % stressTensor % boundaryFlux(i,j,jEq,s1,e1) = 0.5_prec*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1)/&
                                                                                             (myDGSEM % state % boundarySolution(i,j,4,s1,e1) +&
                                                                                              myDGSEM % static % boundarySolution(i,j,4,s1,e1)) +&
                                                                                             myDGSEM % externalState(ii,jj,iEq,bID)/&
                                                                                             (myDGSEM % externalState(ii,jj,4,bID)+&
                                                                                              myDGSEM % static % boundarySolution(i,j,4,s1,e1)) )*&
                                                                                             myDGSEM % mesh % geom(e1) % nHat(k,i,j,s1)
                        ENDDO
                     ENDIF
                  ENDDO
                  
               ENDDO
            ENDDO 

         ENDIF 
         
      ENDDO 
      !$OMP ENDDO
      

#endif

 END SUBROUTINE BoundaryFaceFlux_Fluid

!
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!                                                                                                 !
!  This section of code contains routines for computing tendency from the internal and Riemann    !
!  fluxes and the buoyancy source term.                                                           !
!                                                                                                 !
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!
 SUBROUTINE MappedTimeDerivative_Fluid( myDGSEM )

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) )
      grid = dim3(myDGSEM % mesh % nElems,nEq-1,1)
      
      CALL MappedTimeDerivative_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                             myDGSEM % static % solution_dev, &
                                                             myDGSEM % state % boundaryFlux_dev, &
                                                             myDGSEM % dragProfile_dev, &
                                                             myDGSEM % mesh % geom_dev % Ja_dev, &
                                                             myDGSEM % mesh % geom_dev % J_dev, &
                                                             myDGSEM % dgStorage % bMat_dev, &
                                                             myDGSEM % dgStorage % quadratureWeights_dev, &
                                                             myDGSEM % dgStorage % dMatP_dev, &
                                                             myDGSEM % state % tendency_dev )
#else
   ! Local
   INTEGER    :: iEl, i, j, k, m, iEq, row, col
   REAL(prec) :: F
   REAL(prec) :: pContFlux(0:myDGSEM % N,0:myDGSEM % N,0:myDGSEM % N)
   REAL(prec) :: sContFlux(0:myDGSEM % N,0:myDGSEM % N,0:myDGSEM % N)
   REAL(prec) :: qContFlux(0:myDGSEM % N,0:myDGSEM % N,0:myDGSEM % N)


      !$OMP DO PRIVATE( pContFlux, sContFlux, qContFlux )
      DO iEl = 1, myDGSEM % mesh % nElems
         DO iEq = 1, nEq-1
         
            DO k = 0, myDGSEM % N 
               DO j = 0, myDGSEM % N 
                  DO i = 0, myDGSEM % N
                     sContFlux(i,j,k) = 0.0_prec
                     pContFlux(i,j,k) = 0.0_prec
                     qContFlux(i,j,k) = 0.0_prec
                  ENDDO
               ENDDO 
            ENDDO
         
         ! Here the flux tensor in physical space is calculated and rotated to give the 
         ! contravariant flux tensor in the reference computational domain.
         
         !//////////////////////////////// Advection ///////////////////////////////////////!
            DO col = 1, 3
               DO k = 0, myDGSEM % N 
                  DO j = 0, myDGSEM % N 
                     DO i = 0, myDGSEM % N

                        F = myDGSEM % state % solution(i,j,k,col,iEl)*&
                           (myDGSEM % state % solution(i,j,k,iEq,iEl)+&
                            myDGSEM % static % solution(i,j,k,iEq,iEl))/& 
                           (myDGSEM % state % solution(i,j,k,4,iEl)+&
                            myDGSEM % static % solution(i,j,k,4,iEl)) 
                            
               
                        sContFlux(i,j,k) = sContFlux(i,j,k) + &
                                               myDGSEM % mesh % geom(iEl) % Ja(i,j,k,col,1)*F
                                       
                        pContFlux(i,j,k) = pContFlux(i,j,k) + &
                                               myDGSEM % mesh % geom(iEl) % Ja(i,j,k,col,2)*F 
                                      
                        qContFlux(i,j,k) = qContFlux(i,j,k) + &
                                              myDGSEM % mesh % geom(iEl) % Ja(i,j,k,col,3)*F
                                              
                     ENDDO
                  ENDDO
               ENDDO 
            ENDDO

        ! //////////////////// Pressure (Momentum only) /////////////////////////// !
            IF( iEq <= 3 )THEN
               DO k = 0, myDGSEM % N 
                  DO j = 0, myDGSEM % N 
                     DO i = 0, myDGSEM % N

                        sContFlux(i,j,k) = sContFlux(i,j,k) + &
                                               myDGSEM % mesh % geom(iEl) % Ja(i,j,k,iEq,1)*&
                                               myDGSEM % state % solution(i,j,k,6,iEl)
                                       
                        pContFlux(i,j,k) = pContFlux(i,j,k) + &
                                               myDGSEM % mesh % geom(iEl) % Ja(i,j,k,iEq,2)*&
                                               myDGSEM % state % solution(i,j,k,6,iEl) 
                                      
                        qContFlux(i,j,k) = qContFlux(i,j,k) + &
                                               myDGSEM % mesh % geom(iEl) % Ja(i,j,k,iEq,3)*&
                                               myDGSEM % state % solution(i,j,k,6,iEl)
                                              
                     ENDDO
                  ENDDO
               ENDDO 
               
            ENDIF
      ! Now, the flux divergence is computed by multiplying the internally calculated fluxes by the
      ! DG-Derivative matrix and adding the boundary weighted fluxes.
            DO k = 0, myDGSEM % N
               DO j = 0, myDGSEM % N
                  DO i = 0, myDGSEM % N
      
                     myDGSEM % state % tendency(i,j,k,iEq,iEl) = 0.0_prec
               
                     DO m = 0, myDGSEM % N
                        myDGSEM % state % tendency(i,j,k,iEq,iEl) = myDGSEM % state % tendency(i,j,k,iEq,iEl) + &
                                          myDGSEM % dgStorage % dMatP(m,i)*sContFlux(m,j,k) + &
                                          myDGSEM % dgStorage % dMatP(m,j)*pContFlux(i,m,k) + &
                                          myDGSEM % dgStorage % dMatP(m,k)*qContFlux(i,j,m)
                     ENDDO
                                    
                     myDGSEM % state % tendency(i,j,k,iEq,iEl) = -( myDGSEM % state % tendency(i,j,k,iEq,iEl) + &
                                    ( myDGSEM % state % boundaryFlux(i,k,iEq,SOUTH,iEl)*&
                                      myDGSEM % dgStorage % bmat(j,left-1) + &
                                      myDGSEM % state % boundaryFlux(i,k,iEq,NORTH,iEl)*&
                                      myDGSEM % dgStorage % bMat(j,right-1) )/&
                                      myDGSEM % dgStorage % quadratureWeights(j) + &
                                    ( myDGSEM % state % boundaryFlux(j,k,iEq,WEST,iEl)*&
                                      myDGSEM % dgStorage % bMat(i,left-1) + &
                                      myDGSEM % state % boundaryFlux(j,k,iEq,EAST,iEl)*&
                                      myDGSEM % dgStorage % bMat(i,right-1) )/&
                                      myDGSEM % dgStorage % quadratureWeights(i) + &
                                    ( myDGSEM % state % boundaryFlux(i,j,iEq,BOTTOM,iEl)*&
                                      myDGSEM % dgStorage % bMat(k,left-1) + &
                                      myDGSEM % state % boundaryFlux(i,j,iEq,TOP,iEl)*&
                                      myDGSEM % dgStorage % bMat(k,right-1) )/&
                                      myDGSEM % dgStorage % quadratureWeights(k) )/myDGSEM % mesh % geom(iEl) % J(i,j,k)
                  ENDDO
               ENDDO
            ENDDO
            
            IF( iEq == 1 )THEN
               DO k = 0, myDGSEM % N  
                  DO j = 0, myDGSEM % N 
                     DO i = 0, myDGSEM % N
                        F = sqrt( myDGSEM % state % solution(i,j,k,1,iEl)**2 + &
                                  myDGSEM % state % solution(i,j,k,2,iEl)**2 + &
                                  myDGSEM % state % solution(i,j,k,3,iEl)**2 )/&
                                  ( myDGSEM % static % solution(i,j,k,4,iEl) + myDGSEM % state % solution(i,j,k,4,iEl) )
                        myDGSEM % state % tendency(i,j,k,1,iEl) = myDGSEM % state % tendency(i,j,k,1,iEl) -&
                                                                  myDGSEM % dragProfile(i,j,k,iEl)*myDGSEM % state % solution(i,j,k,1,iEl)*F -&
                                                                  myDGSEM % state % solution(i,j,k,3,iEl)*myDGSEM % params % fRotY +&
                                                                  myDGSEM % state % solution(i,j,k,2,iEl)*myDGSEM % params % fRotZ
                     ENDDO
                  ENDDO
               ENDDO
            
            ELSEIF( iEq == 2 )THEN
               DO k = 0, myDGSEM % N  
                  DO j = 0, myDGSEM % N 
                     DO i = 0, myDGSEM % N
                        F = sqrt( myDGSEM % state % solution(i,j,k,1,iEl)**2 + &
                                  myDGSEM % state % solution(i,j,k,2,iEl)**2 + &
                                  myDGSEM % state % solution(i,j,k,3,iEl)**2 )/&
                                  ( myDGSEM % static % solution(i,j,k,4,iEl) + myDGSEM % state % solution(i,j,k,4,iEl) )
                        myDGSEM % state % tendency(i,j,k,1,iEl) = myDGSEM % state % tendency(i,j,k,1,iEl) -&
                                                                  myDGSEM % dragProfile(i,j,k,iEl)*myDGSEM % state % solution(i,j,k,2,iEl)*F - &
                                                                  myDGSEM % state % solution(i,j,k,1,iEl)*myDGSEM % params % fRotZ +&
                                                                  myDGSEM % state % solution(i,j,k,3,iEl)*myDGSEM % params % fRotX
                     ENDDO
                  ENDDO
               ENDDO
               
            ELSEIF( iEq == 3 )THEN
               DO k = 0, myDGSEM % N  
                  DO j = 0, myDGSEM % N 
                     DO i = 0, myDGSEM % N
                        F = sqrt( myDGSEM % state % solution(i,j,k,1,iEl)**2 + &
                                  myDGSEM % state % solution(i,j,k,2,iEl)**2 + &
                                  myDGSEM % state % solution(i,j,k,3,iEl)**2 )/&
                                  ( myDGSEM % static % solution(i,j,k,4,iEl) + myDGSEM % state % solution(i,j,k,4,iEl) )
                        myDGSEM % state % tendency(i,j,k,3,iEl) = myDGSEM % state % tendency(i,j,k,3,iEl) -&
                                                                  myDGSEM % dragProfile(i,j,k,iEl)*myDGSEM % state % solution(i,j,k,3,iEl)*F - &
                                                                  myDGSEM % state % solution(i,j,k,2,iEl)*myDGSEM % params % fRotX +&
                                                                  myDGSEM % state % solution(i,j,k,1,iEl)*myDGSEM % params % fRotY-&
                                                                 ( myDGSEM % state % solution(i,j,k,4,iEl) )*myDGSEM % params % g  !& 
                     ENDDO
                  ENDDO 
               ENDDO
            ENDIF
         
         ENDDO
      ENDDO
      !$OMP ENDDO
      
#endif

 END SUBROUTINE MappedTimeDerivative_Fluid
!
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!                                                                                                 !
!  This section of code contains routines for computing the gradients of the prognostic variables !
!                                                                                                 !
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!
 SUBROUTINE CalculateStressTensor_Fluid( myDGSEM )
 
   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) )
      grid = dim3(myDGSEM % mesh % nElems,nEq-1,3)
   
      CALL CalculateStressTensor_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                              myDGSEM % static % solution_dev, &
                                                              myDGSEM % dgStorage % dMatP_dev, &
                                                              myDGSEM % dgStorage % bMat_dev, &
                                                              myDGSEM % dgStorage % quadratureWeights_dev, &
                                                              myDGSEM % mesh % geom_dev % Ja_dev, &
                                                              myDGSEM % mesh % geom_dev % J_dev, &
                                                              myDGSEM % stressTensor % boundaryFlux_dev, &
                                                              myDGSEM % stressTensor % solution_dev )
#else
   ! Local
   INTEGER :: iEl, iEq, idir, i, j, k, m, jEq
   REAL(prec) :: pContFlux(0:myDGSEM % N,0:myDGSEM % N,0:myDGSEM % N)
   REAL(prec) :: sContFlux(0:myDGSEM % N,0:myDGSEM % N,0:myDGSEM % N)
   REAL(prec) :: qContFlux(0:myDGSEM % N,0:myDGSEM % N,0:myDGSEM % N)


      !$OMP DO PRIVATE( pContFlux, sContFlux, qContFlux )
      DO iEl = 1, myDGSEM % mesh % nElems
         DO iEq = 1, nEq-1
            DO idir = 1, 3
         
               DO k = 0, myDGSEM % N 
                  DO j = 0, myDGSEM % N 
                     DO i = 0, myDGSEM % N
                        sContFlux(i,j,k) = 0.0_prec
                        pContFlux(i,j,k) = 0.0_prec
                        qContFlux(i,j,k) = 0.0_prec
                     ENDDO
                  ENDDO 
               ENDDO
         
         ! Here the flux tensor in physical space is calculated and rotated to give the 
         ! contravariant flux tensor in the reference computational domain.
               IF( iEq == 4 )THEN
	               DO k = 0, myDGSEM % N 
	                  DO j = 0, myDGSEM % N 
	                     DO i = 0, myDGSEM % N
	
	                        sContFlux(i,j,k) = myDGSEM % mesh % geom(iEl) % Ja(i,j,k,idir,1)*&
	                                           myDGSEM % state % solution(i,j,k,iEq,iEl)
	                                       
	                        pContFlux(i,j,k) = myDGSEM % mesh % geom(iEl) % Ja(i,j,k,idir,2)*&
	                                           myDGSEM % state % solution(i,j,k,iEq,iEl) 
	                                      
	                        qContFlux(i,j,k) = myDGSEM % mesh % geom(iEl) % Ja(i,j,k,idir,3)*&
	                                           myDGSEM % state % solution(i,j,k,iEq,iEl)
	                                              
	                     ENDDO
	                  ENDDO
	               ENDDO 
               ELSE
	               DO k = 0, myDGSEM % N 
	                  DO j = 0, myDGSEM % N 
	                     DO i = 0, myDGSEM % N
	
	                        sContFlux(i,j,k) = myDGSEM % mesh % geom(iEl) % Ja(i,j,k,idir,1)*&
	                                           myDGSEM % state % solution(i,j,k,iEq,iEl)/&
	                                          (myDGSEM % state % solution(i,j,k,4,iEl)+&
	                                           myDGSEM % static % solution(i,j,k,4,iEl) )
	                                       
	                        pContFlux(i,j,k) = myDGSEM % mesh % geom(iEl) % Ja(i,j,k,idir,2)*&
	                                           myDGSEM % state % solution(i,j,k,iEq,iEl)/&
	                                           (myDGSEM % state % solution(i,j,k,4,iEl)+&
	                                           myDGSEM % static % solution(i,j,k,4,iEl) )
	                                      
	                        qContFlux(i,j,k) = myDGSEM % mesh % geom(iEl) % Ja(i,j,k,idir,3)*&
	                                           myDGSEM % state % solution(i,j,k,iEq,iEl)/&
	                                           (myDGSEM % state % solution(i,j,k,4,iEl)+&
	                                           myDGSEM % static % solution(i,j,k,4,iEl) )
	                                              
	                     ENDDO
	                  ENDDO
	               ENDDO 
               ENDIF
	               
      ! Now, the flux divergence is computed by multiplying the internally calculated fluxes by the
      ! DG-Derivative matrix and adding the boundary weighted fluxes.
               jEq = idir + (iEq-1)*3
               DO k = 0, myDGSEM % N
                  DO j = 0, myDGSEM % N
                     DO i = 0, myDGSEM % N
      
                        myDGSEM % stressTensor % solution(i,j,k,jEq,iEl) = 0.0_prec
                
                        DO m = 0, myDGSEM % N
                           myDGSEM % stressTensor % solution(i,j,k,jEq,iEl) = myDGSEM % stressTensor % solution(i,j,k,jEq,iEl) + &
                                          myDGSEM % dgStorage % dMatP(m,i)*sContFlux(m,j,k) + &
                                          myDGSEM % dgStorage % dMatP(m,j)*pContFlux(i,m,k) + &
                                          myDGSEM % dgStorage % dMatP(m,k)*qContFlux(i,j,m)
                        ENDDO
         
                        myDGSEM % stressTensor % solution(i,j,k,jEq,iEl) = ( myDGSEM % stressTensor % solution(i,j,k,jEq,iEl) + &
                                    ( myDGSEM % stressTensor % boundaryFlux(i,k,jEq,SOUTH,iEl)*&
                                      myDGSEM % dgStorage % bmat(j,left-1) + &
                                      myDGSEM % stressTensor % boundaryFlux(i,k,jEq,NORTH,iEl)*&
                                      myDGSEM % dgStorage % bMat(j,right-1) )/&
                                      myDGSEM % dgStorage % quadratureWeights(j) + &
                                    ( myDGSEM % stressTensor % boundaryFlux(j,k,jEq,WEST,iEl)*&
                                      myDGSEM % dgStorage % bMat(i,left-1) + &
                                      myDGSEM % stressTensor % boundaryFlux(j,k,jEq,EAST,iEl)*&
                                      myDGSEM % dgStorage % bMat(i,right-1) )/&
                                      myDGSEM % dgStorage % quadratureWeights(i) + &
                                    ( myDGSEM % stressTensor % boundaryFlux(i,j,jEq,BOTTOM,iEl)*&
                                      myDGSEM % dgStorage % bMat(k,left-1) + &
                                      myDGSEM % stressTensor % boundaryFlux(i,j,jEq,TOP,iEl)*&
                                      myDGSEM % dgStorage % bMat(k,right-1) )/&
                                      myDGSEM % dgStorage % quadratureWeights(k) )/myDGSEM % mesh % geom(iEl) % J(i,j,k)
                     ENDDO
                  ENDDO
               ENDDO
            
            ENDDO
         ENDDO    
      ENDDO
      !$OMP ENDDO
      
#endif

 END SUBROUTINE CalculateStressTensor_Fluid
!
 SUBROUTINE CalculateBoundaryStress_Fluid( myDGSEM ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
     tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                   4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                   15 )
     grid = dim3(myDGSEM % mesh % nElems, 1, 1)  
     
     CALL CalculateBoundaryStress_CUDAKernel<<<grid, tBlock>>>( myDGSEM % stressTensor % solution_dev, &
                                                                myDGSEM % dgStorage % bMat_dev, &
                                                                myDGSEM % stressTensor % boundarySolution_dev )
#else
   ! Local
   INTEGER :: iEq, iEl, i, j, k

      !$OMP DO
      DO iEl = 1, myDGSEM % mesh % nElems
         DO iEq = 1, (myDGSEM % nEq-1)*3
            DO k = 0, myDGSEM % N
               DO j = 0, myDGSEM % N
               
                  myDGSEM % stressTensor % boundarySolution(j,k,iEq,WEST,iEl)   = 0.0_prec
                  myDGSEM % stressTensor % boundarySolution(j,k,iEq,EAST,iEl)   = 0.0_prec
                  myDGSEM % stressTensor % boundarySolution(j,k,iEq,SOUTH,iEl)  = 0.0_prec
                  myDGSEM % stressTensor % boundarySolution(j,k,iEq,NORTH,iEl)  = 0.0_prec
                  myDGSEM % stressTensor % boundarySolution(j,k,iEq,BOTTOM,iEl) = 0.0_prec
                  myDGSEM % stressTensor % boundarySolution(j,k,iEq,TOP,iEL)    = 0.0_prec
    
                  DO i = 0, myDGSEM % N

                     myDGSEM % stressTensor % boundarySolution(j,k,iEq,WEST,iEl)  = &
                        myDGSEM % stressTensor % boundarySolution(j,k,iEq,WEST,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,0)*myDGSEM % stressTensor % solution(i,j,k,iEq,iEl)

                     myDGSEM % stressTensor % boundarySolution(j,k,iEq,EAST,iEl)  = &
                        myDGSEM % stressTensor % boundarySolution(j,k,iEq,EAST,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,1)*myDGSEM % stressTensor % solution(i,j,k,iEq,iEl)

                     myDGSEM % stressTensor % boundarySolution(j,k,iEq,SOUTH,iEl)  = &
                        myDGSEM % stressTensor % boundarySolution(j,k,iEq,SOUTH,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,0)*myDGSEM % stressTensor % solution(j,i,k,iEq,iEl)

                     myDGSEM % stressTensor % boundarySolution(j,k,iEq,NORTH,iEl)   = &
                        myDGSEM % stressTensor % boundarySolution(j,k,iEq,NORTH,iEl)  +  &
                        myDGSEM % dgStorage % bMat(i,1)*myDGSEM % stressTensor % solution(j,i,k,iEq,iEl)

                     myDGSEM % stressTensor % boundarySolution(j,k,iEq,BOTTOM,iEl)  = &
                        myDGSEM % stressTensor % boundarySolution(j,k,iEq,BOTTOM,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,0)*myDGSEM % stressTensor % solution(j,k,i,iEq,iEl)

                     myDGSEM % stressTensor % boundarySolution(j,k,iEq,TOP,iEl)  = &
                        myDGSEM % stressTensor % boundarySolution(j,k,iEq,TOP,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,1)*myDGSEM % stressTensor % solution(j,k,i,iEq,iEl)

                  ENDDO
               
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !$OMP ENDDO
      
#endif
      
 END SUBROUTINE CalculateBoundaryStress_Fluid
!
 SUBROUTINE UpdateExternalStress_Fluid( myDGSEM, tn ) ! ////////// !

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)      :: tn
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
     tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                   4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                   1 )
     grid = dim3(myDGSEM % nBoundaryFaces,(nEq-1)*3,1)  
     
     CALL UpdateExternalStress_CUDAKernel<<<grid, tBlock>>>( myDGSEM % extComm % boundaryIDs_dev, &            ! I
                                                            myDGSEM % mesh % faces_dev % elementIDs, &         ! I
                                                            myDGSEM % mesh % faces_dev % elementSides, &       ! I
                                                            myDGSEM % extComm % extProcIDs_dev, &              ! I
                                                            myDGSEM % externalStress_dev, &                    ! O
                                                            myDGSEM % stressTensor % boundarySolution_dev, &   ! I  
                                                            myDGSEM % prescribedStress_dev, &                  ! I
                                                            myDGSEM % mesh % geom_dev % nHat_dev )             ! I

#else
   ! Local
   INTEGER    :: iEl, iFace, bFaceID, i, j, k, iEq
   INTEGER    :: bID, p2
   INTEGER    :: e1, e2, s1, s2

      !$OMP DO
      DO bID = 1, myDGSEM % extComm % nBoundaries

         iFace = myDGSEM % extComm % boundaryIDs( bID ) ! Obtain the process-local face id for this boundary-face id
         e1    = myDGSEM % mesh % Faces(iFace) % elementIDs(1)
         s1    = myDGSEM % mesh % Faces(iFace) % elementSides(1)
         e2    = myDGSEM % mesh % Faces(iFace) % elementIDs(2)
         p2    = myDGSEM % extComm % extProcIDs( bID )
         
         IF( p2 == myDGSEM % myRank )THEN ! Enforce no boundary flux due to the fluid stress
            DO j = 0, myDGSEM % N 
               DO i = 0, myDGSEM % N
                  DO iEq = 1, (myDGSEM % nEq-1)*3
                        myDGSEM % externalStress(i,j,iEq,bID) = -myDGSEM % stressTensor % boundarySolution(i,j,iEq,s1,e1)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         
      ENDDO
      !$OMP ENDDO

#endif     

 END SUBROUTINE UpdateExternalStress_Fluid
!
 SUBROUTINE InternalStressFlux_Fluid( myDGSEM )

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                    1 )
      grid = dim3(myDGSEM % mesh % nFaces,nEq-1,1)  
      
      CALL InternalStressFlux_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mesh % faces_dev % elementIDs, &
                                                  myDGSEM % mesh % faces_dev % elementSides, &
                                                  myDGSEM % mesh % faces_dev % boundaryID, &
                                                  myDGSEM % mesh % faces_dev % iMap, &
                                                  myDGSEM % mesh % faces_dev % jMap, &
                                                  myDGSEM % mesh % geom_dev % nHat_dev, &
                                                  myDGSEM % state % boundarySolution_dev, &
                                                  myDGSEM % static % boundarySolution_dev, &
                                                  myDGSEM % stressTensor % boundarySolution_dev, &
                                                  myDGSEM % sgsCoeffs % boundarySolution_dev, &
                                                  myDGSEM % externalState_dev, &
                                                  myDGSEM % externalStress_dev, &
                                                  myDGSEM % stressTensor % boundaryFlux_dev )
#else
   ! Local
   REAL(prec) :: norm, rhoIn, rhoOut
   INTEGER :: iEl, iFace
   INTEGER    :: i, j, k, m, iEq, jEq
   INTEGER    :: ii, jj, bID
   INTEGER    :: e1, s1, e2, s2

      !$OMP DO
      DO iFace = 1, myDGSEM % mesh % nFaces

         
         e1 = myDGSEM % mesh % faces(iFace) % elementIDs(1)
         s1 = myDGSEM % mesh % faces(iFace) % elementSides(1)
         e2 = myDGSEM % mesh % faces(iFace) % elementIDs(2)
         s2 = ABS(myDGSEM % mesh % faces(iFace) % elementSides(2))

         IF( e2 > 0 )THEN

            DO j = 0, myDGSEM % N  
               DO i = 0, myDGSEM % N

                  IF( i == 0 )THEN
                     IF( j == 0 )THEN
                        ii = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % iStart) + &
                             (myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % jStart)
                        jj = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % jStart) + &
                             (myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % iStart)
                     ELSE
                        ii = myDGSEM % mesh % faces(iFace) % swapDimensions*&
                             (ii+myDGSEM % mesh % faces(iFace) % jInc) + &
                             (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             myDGSEM % mesh % faces(iFace) % iStart
                        jj = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (jj+myDGSEM % mesh % faces(iFace) % jInc) +&
                             myDGSEM % mesh % faces(iFace) % swapDimensions*&
                             myDGSEM % mesh % faces(iFace) % jStart
                     ENDIF
                  ELSE
                     ii = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                          (ii + myDGSEM % mesh % faces(iFace) % iInc) +&
                          myDGSEM % mesh % faces(iFace) % swapDimensions*ii
                     jj = myDGSEM % mesh % faces(iFace) % swapDimensions*&
                          (jj+myDGSEM % mesh % faces(iFace) % iInc) + &
                          (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*jj
                  ENDIF

                   norm = sqrt( myDGSEM % mesh % geom(e1) % nHat(1,i,j,s1)**2 + &
                               myDGSEM % mesh % geom(e1) % nHat(2,i,j,s1)**2 + &
                               myDGSEM % mesh % geom(e1) % nHat(3,i,j,s1)**2 )
               
                   DO iEq = 1, nEq-1
                   
                      myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = &
	                          0.5_prec*( myDGSEM % sgsCoeffs % boundarysolution(ii,jj,iEq,s2,e2)*myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2)-&
                                        myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) )/&
	                            myDGSEM % params % viscLengthScale*norm
                               
                      IF( iEq == 4 )THEN
	                      DO m = 1, 3    
	                         jEq = m + (iEq-1)*3  
	                         myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) + &
	                           0.5_prec*(myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*myDGSEM % stressTensor % boundarysolution(i,j,jEq,s1,e1)+&
	                                  myDGSEM % sgsCoeffs % boundarysolution(ii,jj,iEq,s2,e2)*myDGSEM % stressTensor % boundarysolution(ii,jj,jEq,s2,e2))*&
	                          myDGSEM % mesh % geom(e1) % nHat(m,i,j,s1)
	                      ENDDO
	                      
	                   ELSE
                      
	                      rhoOut = (myDGSEM % static % boundarySolution(ii,jj,4,s2,e2)+myDGSEM % state % boundarySolution(ii,jj,4,s2,e2) )
	                      rhoIn  = (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % boundarySolution(i,j,4,s1,e1) )
	                   
	                      DO m = 1, 3    
	                         jEq = m + (iEq-1)*3  
	                         myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) + &
                              0.5_prec*( rhoIn*myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*myDGSEM % stressTensor % boundarysolution(i,j,jEq,s1,e1)+&
	                                       rhoOut*myDGSEM % sgsCoeffs % boundarysolution(ii,jj,iEq,s2,e2)*myDGSEM % stressTensor % boundarysolution(ii,jj,jEq,s2,e2))*&
	                          myDGSEM % mesh % geom(e1) % nHat(m,i,j,s1)
	                      ENDDO
	                   ENDIF
	                      
                      myDGSEM % stressTensor % boundaryFlux(ii,jj,iEq,s2,e2) = -myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1)
                   ENDDO
                ENDDO
             ENDDO     

         ENDIF 
         
      ENDDO 
      !$OMP ENDDO
      
#endif

 END SUBROUTINE InternalStressFlux_Fluid
!
 SUBROUTINE BoundaryStressFlux_Fluid( myDGSEM )

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                    1 )
      grid = dim3(myDGSEM % mesh % nFaces,nEq-1,1)  
      
      CALL BoundaryStressFlux_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mesh % faces_dev % elementIDs, &
                                                  myDGSEM % mesh % faces_dev % elementSides, &
                                                  myDGSEM % mesh % faces_dev % boundaryID, &
                                                  myDGSEM % mesh % faces_dev % iMap, &
                                                  myDGSEM % mesh % faces_dev % jMap, &
                                                  myDGSEM % mesh % geom_dev % nHat_dev, &
                                                  myDGSEM % state % boundarySolution_dev, &
                                                  myDGSEM % static % boundarySolution_dev, &
                                                  myDGSEM % stressTensor % boundarySolution_dev, &
                                                  myDGSEM % sgsCoeffs % boundarySolution_dev, &
                                                  myDGSEM % externalSGS_dev, &
                                                  myDGSEM % externalState_dev, &
                                                  myDGSEM % externalStress_dev, &
                                                  myDGSEM % stressTensor % boundaryFlux_dev )
#else
   ! Local
   REAL(prec) :: norm, rhoIn, rhoOut
   INTEGER :: iEl, iFace
   INTEGER    :: i, j, k, m, iEq, jEq
   INTEGER    :: ii, jj, bID
   INTEGER    :: e1, s1, e2, s2

      !$OMP DO
      DO iFace = 1, myDGSEM % mesh % nFaces

         
         e1 = myDGSEM % mesh % faces(iFace) % elementIDs(1)
         s1 = myDGSEM % mesh % faces(iFace) % elementSides(1)
         e2 = myDGSEM % mesh % faces(iFace) % elementIDs(2)
         s2 = ABS(myDGSEM % mesh % faces(iFace) % elementSides(2))

         IF( e2 < 0 )THEN

            DO j = 0, myDGSEM % N  
               DO i = 0, myDGSEM % N

                  IF( i == 0 )THEN
                     IF( j == 0 )THEN
                        ii = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % iStart) + &
                             (myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % jStart)
                        jj = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % jStart) + &
                             (myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (myDGSEM % mesh % faces(iFace) % iStart)
                     ELSE
                        ii = myDGSEM % mesh % faces(iFace) % swapDimensions*&
                             (ii+myDGSEM % mesh % faces(iFace) % jInc) + &
                             (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             myDGSEM % mesh % faces(iFace) % iStart
                        jj = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                             (jj+myDGSEM % mesh % faces(iFace) % jInc) +&
                             myDGSEM % mesh % faces(iFace) % swapDimensions*&
                             myDGSEM % mesh % faces(iFace) % jStart
                     ENDIF
                  ELSE
                     ii = (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*&
                          (ii + myDGSEM % mesh % faces(iFace) % iInc) +&
                          myDGSEM % mesh % faces(iFace) % swapDimensions*ii
                     jj = myDGSEM % mesh % faces(iFace) % swapDimensions*&
                          (jj+myDGSEM % mesh % faces(iFace) % iInc) + &
                          (1-myDGSEM % mesh % faces(iFace) % swapDimensions)*jj
                  ENDIF

                   norm = sqrt( myDGSEM % mesh % geom(e1) % nHat(1,i,j,s1)**2 + &
                               myDGSEM % mesh % geom(e1) % nHat(2,i,j,s1)**2 + &
                               myDGSEM % mesh % geom(e1) % nHat(3,i,j,s1)**2 )
               
                  bID  = myDGSEM % mesh % faces(iFace) % boundaryID
                  IF( bID < 0 )THEN ! Physical Boundary
                     bID = ABS(bID)
                     DO iEq = 1, nEq-1
                        myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = &
                                0.5_prec*myDGSEM % sgsCoeffs % boundarySolution(i,j,iEq,s1,e1)*&
                                ( myDGSEM % externalState(ii,jj,iEq,bID) - myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) )/&
                                myDGSEM % params % viscLengthScale*norm
                                 
                        IF( iEq == 4 )THEN
                           DO m = 1, 3    
                              jEq = m + (iEq-1)*3  
                              myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = &
                                      myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) + &
                                      0.5_prec*myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*&
                                      ( myDGSEM % stressTensor % boundarysolution(i,j,jEq,s1,e1)+ myDGSEM % externalStress(ii,jj,jEq,bID) )*&
                                      myDGSEM % mesh % geom(e1) % nHat(m,i,j,s1)
                           ENDDO
                        
                        ELSE
                        
                           rhoOut = (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % externalState(ii,jj,4,bID) )
                           rhoIn  = (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % boundarySolution(i,j,4,s1,e1) )
                     
                           DO m = 1, 3    
                              jEq = m + (iEq-1)*3  
                              myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = &
                                      myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) + &
                                      0.5_prec*myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*&
                                      ( rhoIn*myDGSEM % stressTensor % boundarysolution(i,j,jEq,s1,e1)+&
                                        rhoOut*myDGSEM % externalStress(ii,jj,jEq,bID) )*&
                                      myDGSEM % mesh % geom(e1) % nHat(m,i,j,s1)
                           ENDDO
                        ENDIF
                     ENDDO
                  ELSE ! Neighboring process
                     DO iEq = 1, nEq-1
                        myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = &
                                0.5_prec*(myDGSEM % externalSGS(ii,jj,iEq,bID)*myDGSEM % externalState(ii,jj,iEq,bID) -&
                                          myDGSEM % sgsCoeffs % boundarySolution(i,j,iEq,s1,e1)*myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) )/&
                                myDGSEM % params % viscLengthScale*norm
                                 
                        IF( iEq == 4 )THEN
                           DO m = 1, 3    
                              jEq = m + (iEq-1)*3  
                              myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = &
                                      myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) + &
                                      0.5_prec*( myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*&
                                                 myDGSEM % stressTensor % boundarysolution(i,j,jEq,s1,e1)+ &
                                                 myDGSEM % externalSGS(ii,jj,iEq,bID )*myDGSEM % externalStress(ii,jj,jEq,bID) )*&
                                      myDGSEM % mesh % geom(e1) % nHat(m,i,j,s1)
                           ENDDO
                        
                        ELSE
                        
                           rhoOut = (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % externalState(ii,jj,4,bID) )
                           rhoIn  = (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % boundarySolution(i,j,4,s1,e1) )
                     
                           DO m = 1, 3    
                              jEq = m + (iEq-1)*3  
                              myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = &
                                      myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) + &
                                      0.5_prec*( myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*&
                                                 rhoIn*myDGSEM % stressTensor % boundarysolution(i,j,jEq,s1,e1)+&
                                                 myDGSEM % externalSGS(ii,jj,iEq,bID)*&
                                                 rhoOut*myDGSEM % externalStress(ii,jj,jEq,bID) )*&
                                      myDGSEM % mesh % geom(e1) % nHat(m,i,j,s1)
                           ENDDO
                        ENDIF
                     ENDDO

                  ENDIF
               ENDDO
            ENDDO 

         ENDIF 
         
      ENDDO 
      !$OMP ENDDO
      
#endif

 END SUBROUTINE BoundaryStressFlux_Fluid
!
 SUBROUTINE StressDivergence_Fluid( myDGSEM )

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) )
      grid = dim3(myDGSEM % mesh % nElems,nEq-1,1)
      
      CALL StressDivergence_CUDAKernel<<<grid,tBlock>>>( myDGSEM % stressTensor % solution_dev, &
                                                             myDGSEM % stressTensor % boundaryFlux_dev, &
                                                             myDGSEM % state % solution_dev, &
                                                             myDGSEM % static % solution_dev, &
                                                             myDGSEM % sgsCoeffs % solution_dev, &
                                                             myDGSEM % mesh % geom_dev % Ja_dev, &
                                                             myDGSEM % mesh % geom_dev % J_dev, &
                                                             myDGSEM % dgStorage % bMat_dev, &
                                                             myDGSEM % dgStorage % quadratureWeights_dev, &
                                                             myDGSEM % dgStorage % dMatP_dev, &
                                                             myDGSEM % stressTensor % tendency_dev )
#else
   ! Local
   INTEGER    :: iEl, i, j, k, m, iEq, row, col, jEq
   REAL(prec) :: F
   REAL(prec) :: pContFlux(0:myDGSEM % N,0:myDGSEM % N,0:myDGSEM % N)
   REAL(prec) :: sContFlux(0:myDGSEM % N,0:myDGSEM % N,0:myDGSEM % N)
   REAL(prec) :: qContFlux(0:myDGSEM % N,0:myDGSEM % N,0:myDGSEM % N)

      !$OMP DO PRIVATE( pContFlux, sContFlux, qContFlux )
      DO iEl = 1, myDGSEM % mesh % nElems
         DO iEq = 1, nEq-1
         
            DO k = 0, myDGSEM % N 
               DO j = 0, myDGSEM % N 
                  DO i = 0, myDGSEM % N
                     sContFlux(i,j,k) = 0.0_prec
                     pContFlux(i,j,k) = 0.0_prec
                     qContFlux(i,j,k) = 0.0_prec
                  ENDDO
               ENDDO 
            ENDDO
         
         ! Here the flux tensor in physical space is calculated and rotated to give the 
         ! contravariant flux tensor in the reference computational domain.
         
            DO col = 1, 3
               jEq = col + (iEq-1)*3
               DO k = 0, myDGSEM % N 
                  DO j = 0, myDGSEM % N 
                     DO i = 0, myDGSEM % N

                        IF( iEq == 4 )THEN
	                        F = myDGSEM % stressTensor % solution(i,j,k,jEq,iEl)*&
	                            myDGSEM % sgsCoeffs % solution(i,j,k,iEq,iEl)
                            
                        ELSE
	                        F = myDGSEM % stressTensor % solution(i,j,k,jEq,iEl)*&
	                            (myDGSEM % state % solution(i,j,k,4,iEl)+&
	                             myDGSEM % static % solution(i,j,k,4,iEl))*&
	                            myDGSEM % sgsCoeffs % solution(i,j,k,iEq,iEl)
                        ENDIF
                        sContFlux(i,j,k) = sContFlux(i,j,k) + &
                                               myDGSEM % mesh % geom(iEl) % Ja(i,j,k,col,1)*F
                                       
                        pContFlux(i,j,k) = pContFlux(i,j,k) + &
                                               myDGSEM % mesh % geom(iEl) % Ja(i,j,k,col,2)*F 
                                      
                        qContFlux(i,j,k) = qContFlux(i,j,k) + &
                                              myDGSEM % mesh % geom(iEl) % Ja(i,j,k,col,3)*F
                                              
                     ENDDO
                  ENDDO
               ENDDO 
            ENDDO

      ! Now, the flux divergence is computed by multiplying the internally calculated fluxes by the
      ! DG-Derivative matrix and adding the boundary weighted fluxes.
            DO k = 0, myDGSEM % N
               DO j = 0, myDGSEM % N
                  DO i = 0, myDGSEM % N
      
                     myDGSEM % stressTensor % tendency(i,j,k,iEq,iEl) = 0.0_prec
               
                     DO m = 0, myDGSEM % N
                        myDGSEM % stressTensor % tendency(i,j,k,iEq,iEl) = myDGSEM % stressTensor % tendency(i,j,k,iEq,iEl) + &
                                          myDGSEM % dgStorage % dMatP(m,i)*sContFlux(m,j,k) + &
                                          myDGSEM % dgStorage % dMatP(m,j)*pContFlux(i,m,k) + &
                                          myDGSEM % dgStorage % dMatP(m,k)*qContFlux(i,j,m)
                     ENDDO
         
                     myDGSEM % stressTensor % tendency(i,j,k,iEq,iEl) = (myDGSEM % stressTensor % tendency(i,j,k,iEq,iEl) + &
                                    ( myDGSEM % stressTensor % boundaryFlux(i,k,iEq,SOUTH,iEl)*&
                                      myDGSEM % dgStorage % bmat(j,left-1) + &
                                      myDGSEM % stressTensor % boundaryFlux(i,k,iEq,NORTH,iEl)*&
                                      myDGSEM % dgStorage % bMat(j,right-1) )/&
                                      myDGSEM % dgStorage % quadratureWeights(j) + &
                                    ( myDGSEM % stressTensor % boundaryFlux(j,k,iEq,WEST,iEl)*&
                                      myDGSEM % dgStorage % bMat(i,left-1) + &
                                      myDGSEM % stressTensor % boundaryFlux(j,k,iEq,EAST,iEl)*&
                                      myDGSEM % dgStorage % bMat(i,right-1) )/&
                                      myDGSEM % dgStorage % quadratureWeights(i) + &
                                    ( myDGSEM % stressTensor % boundaryFlux(i,j,iEq,BOTTOM,iEl)*&
                                      myDGSEM % dgStorage % bMat(k,left-1) + &
                                      myDGSEM % stressTensor % boundaryFlux(i,j,iEq,TOP,iEl)*&
                                      myDGSEM % dgStorage % bMat(k,right-1) )/&
                                      myDGSEM % dgStorage % quadratureWeights(k) )/myDGSEM % mesh % geom(iEl) % J(i,j,k)
                  ENDDO
               ENDDO
            ENDDO
            
         ENDDO
      ENDDO
      !$OMP ENDDO
      
#endif
       
 END SUBROUTINE StressDivergence_Fluid
!
 SUBROUTINE EquationOfState_Fluid( myDGSEM )
   ! This routine calculates the anomalous pressure referenced to the static state.
   ! The pressure is calculated using the ideal gas law.
   IMPLICIT NONE
   CLASS( Fluid ), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
   ! Local
    TYPE(dim3) :: grid, tBlock
  
      ! How should we pick the thread and block size
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) )
      grid = dim3(myDGSEM % mesh % nElems,1,1) 
      
      CALL EquationOfState_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                        myDGSEM % static % solution_dev ) 
                                                        
#else
   ! Local
   INTEGER :: iEl, i, j, k
   REAL(prec) :: hCapRatio, rC, rhoT


      hCapRatio = ( myDGSEM % params % R + myDGSEM % params % Cv ) / myDGSEM % params % Cv
      rC        = myDGSEM % params % R / ( myDGSEM % params % R + myDGSEM % params % Cv )

      !$OMP DO
      DO iEl = 1, myDGSEM % mesh % nElems
         DO k = 0, myDGSEM % N
            DO j = 0, myDGSEM % N
               DO i = 0, myDGSEM % N

                  ! Pressure = rho*e*R/Cv (Ideal Gas Law, P = rho*R*T, thermo ~> T = theta*(P/P0)^r)
                  ! Then P = (rho*theta*R/P0^rC)^(Cp/Cv)
                  ! And P' = P - P_static
                  rhoT = (myDGSEM % static % solution(i,j,k,5,iEl) + myDGSEM % state % solution(i,j,k,5,iEl) )
                  myDGSEM % state % solution(i,j,k,6,iEl) = myDGSEM % params % P0*( rhoT*myDGSEM % params % R/myDGSEM % params % P0 )**hCapRatio -&
                                                              myDGSEM % static % solution(i,j,k,6,iEl)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !$OMP ENDDO
      
#endif  
    
 END SUBROUTINE EquationOfState_Fluid
!
 SUBROUTINE CalculateStaticState_Fluid( myDGSEM )
   IMPLICIT NONE
   CLASS( Fluid ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: i, j, k, iEl,iEq
   REAL(prec) :: z, H, P0, Cp, T, T0, dTdz, P, rC, g, R

      R    = myDGSEM % params % R
      Cp   = (R + myDGSEM % params % Cv)
      rC   = R/Cp
      g    = myDGSEM % params % g
      H    = myDGSEM % params % zScale
      T0   = myDGSEM % params % T0
      P0   = myDGSEM % params % P0
      dTdz = myDGSEM % params % dTdz
      
      ! /////////////////////  Build the Static/Background State ///////////////////////// !

      !$OMP DO
      DO iEl = 1, myDGSEM % mesh % nElems
         DO iEq = 1, nEq
            DO k = 0, myDGSEM % N
               DO j = 0, myDGSEM % N
                  DO i = 0, myDGSEM % N
                     myDGSEM % state % solution(i,j,k,iEq,iEl)  = 0.0_prec
                     myDGSEM % static % solution(i,j,k,iEq,iEl) = 0.0_prec
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !$OMP ENDDO
      
      !$OMP DO
      DO iEl = 1, myDGSEM % mesh % nElems
         DO k = 0, myDGSEM % N
            DO j = 0, myDGSEM % N
               DO i = 0, myDGSEM % N
            
                  z = myDGSEM % mesh % geom(iEl) % z(i,j,k)
                  
                  ! The static profile is determined from hydrostatic balance, the equation of state,
                  ! and a prescribed potential temperature profile. 
                  ! ** The potential temperature is assumed to vary linearly with z ** 
                  
                  T = T0 + dTdz*z ! Potential temperature
                  IF( dTdz == 0.0_prec )THEN
                     P = P0*( 1.0_prec - g*z/(T0*Cp) )**(Cp/R)
                  ELSE
                     P = P0*( 1.0_prec - g*rC/R*log( (T/T0)**(1.0_prec/dTdz) ) )**(Cp/R)
                  ENDIF
                  
                  ! Density
                  myDGSEM % static % solution(i,j,k,4,iEl) = (P/( T*R*(P/P0)**rC) )
                  
                  ! Potential Temperature (weighted with density)
                  myDGSEM % static % solution(i,j,k,5,iEl) = myDGSEM % static % solution(i,j,k,4,iEl)*T
                  
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !$OMP ENDDO
      
#ifdef HAVE_CUDA
      myDGSEM % static % solution_dev = myDGSEM % static % solution
#endif
      ! This routine Calculates the pressure
      CALL myDGSEM % EquationOfState( )
#ifdef HAVE_CUDA
      myDGSEM % state % solution = myDGSEM % state % solution_dev
#endif

      !$OMP MASTER
      DO iEl = 1, myDGSEM % mesh % nElems
         myDGSEM % static % solution(:,:,:,6,iEl) = myDGSEM % state % solution(:,:,:,6,iEl)
      ENDDO
      !$OMP END MASTER
      
#ifdef HAVE_CUDA
      myDGSEM % static % solution_dev = myDGSEM % static % solution
      myDGSEM % state % solution_dev  = 0.0_prec
#endif
      
      !$OMP MASTER
      myDGSEM % state % solution = 0.0_prec
      !$OMP END MASTER
      ! Add a CUDA Kernel !
      
 END SUBROUTINE CalculateStaticState_Fluid
!                      
!
!
!==================================================================================================!
!-------------------------------------- FILE I/O ROUTINES -----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE WriteTecplot_Fluid( myDGSEM )

  IMPLICIT NONE
 
  CLASS( Fluid ), INTENT(inout) :: myDGsem
  !LOCAL
  REAL(prec)  :: x(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot)
  REAL(prec)  :: y(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot)
  REAL(prec)  :: z(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot)
  REAL(prec)  ::  sol(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:nEq)
  REAL(prec)  :: bsol(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:nEq)

  INTEGER       :: i, j, k, iEl, iEq, fUnit
  CHARACTER(5)  :: zoneID
  CHARACTER(4)  :: rankChar
  REAL(prec)    :: hCapRatio, c, T
  CHARACTER(13) :: timeStampString
   
      timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )
  
      hCapRatio = ( myDGSEM % params % R + myDGSEM % params % Cv ) / myDGSEM % params % Cv

      WRITE(rankChar,'(I4.4)') myDGSEM % myRank
      
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE= 'State.'//rankChar//'.'//timeStampString//'.tec', &
            FORM='formatted', &
            STATUS='replace')  
    
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "u", "v", "w", "rho", "Pot. Temp.", "Pressure",'//&
                                ' "u_b", "v_b", "w_b", "rho_b", "Pot. Temp._b", "Pressure_b", "Drag", "c" '
 
      DO iEl = 1, myDGsem % mesh % nElems

         x = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % geom(iEl) % x )
         y = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % geom(iEl) % y )
         z = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % geom(iEl) % z )
           
         DO iEq = 1, nEq
            sol(:,:,:,iEq) = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % state % solution(:,:,:,iEq,iEl) )
         ENDDO
      
         DO iEq = 1, nEq
            bsol(:,:,:,iEq) = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % static % solution(:,:,:,iEq,iEl) )
         ENDDO
         
         WRITE(zoneID,'(I5.5)') myDGSEM % mesh % elements(iEl) % elementID
         WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',myDGSEM % params % nPlot+1,&
                                                     ', J=',myDGSEM % params % nPlot+1,&
                                                     ', K=',myDGSEM % params % nPlot+1,',F=POINT'

         DO k = 0, myDGSEM % params % nPlot
            DO j = 0, myDGSEM % params % nPlot
               DO i = 0, myDGSEM % params % nPlot
                  T =   (bsol(i,j,k,5) + sol(i,j,k,5))/(bsol(i,j,k,4)+sol(i,j,k,4) )
                          
                  ! Sound speed estimate for the external and internal states
                  c = sqrt( myDGSEM % params % R*T*( ( sol(i,j,k,6) + bsol(i,j,k,6) )/myDGSEM % params % P0 )**hCapRatio   )
                  WRITE(fUnit,'(17(E15.7,1x))') x(i,j,k), y(i,j,k), z(i,j,k),&
                                  sol(i,j,k,1)/( sol(i,j,k,4) + bsol(i,j,k,4) ), &
                                  sol(i,j,k,2)/( sol(i,j,k,4) + bsol(i,j,k,4) ), &
                                  sol(i,j,k,3)/( sol(i,j,k,4) + bsol(i,j,k,4) ), &
                                  sol(i,j,k,4), &
                                  (sol(i,j,k,5) + bsol(i,j,k,5))/( sol(i,j,k,4) + bsol(i,j,k,4) ), &
                                  sol(i,j,k,6), &
                                  bsol(i,j,k,1)/( bsol(i,j,k,4) ), &
                                  bsol(i,j,k,2)/( bsol(i,j,k,4) ), &
                                  bsol(i,j,k,3)/( bsol(i,j,k,4) ), &
                                  bsol(i,j,k,4), &
                                  bsol(i,j,k,5)/( bsol(i,j,k,4) ),&
                                  bsol(i,j,k,6),&
                                  myDGSEM % dragProfile(i,j,k,iEl), c
               ENDDO
            ENDDO
         ENDDO
        
      ENDDO

      CLOSE(UNIT=fUnit)

 END SUBROUTINE WriteTecplot_Fluid
!
#ifdef DIAGNOSTICS
 SUBROUTINE OpenDiagnosticsFiles_Fluid( myDGSEM, fileUnits )
   IMPLICIT NONE
   CLASS( Fluid  ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(out)           :: fileUnits(1:nDiagnostics) 
   ! Local
   CHARACTER(13) :: timeStampString
   

      myDGSEM % volume = 0.0_prec
      myDGSEM % mass   = 0.0_prec
      myDGSEM % KE     = 0.0_prec
      myDGSEM % PE     = 0.0_prec
      myDGSEM % heat   = 0.0_prec

      IF( myDGSEM % myRank == 0 )THEN
        timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )


        OPEN( UNIT=NewUnit(fileUnits(1)), &
              FILE='Mass.'//timeStampString//'.curve', &
              FORM='FORMATTED', &
              STATUS='REPLACE' )
        WRITE(fileUnits(1),*) '#TotalMass'

        OPEN( UNIT=NewUnit(fileUnits(2)), &
              FILE='KineticEnergy.'//timeStampString//'.curve', &
              FORM='FORMATTED', &
              STATUS='REPLACE' )
        WRITE(fileUnits(2),*) '#TotalKineticEnergy'

        OPEN( UNIT=NewUnit(fileUnits(3)), &
              FILE='PotentialEnergy.'//timeStampString//'.curve', &
              FORM='FORMATTED', &
              STATUS='REPLACE' )
        WRITE(fileUnits(3),*) '#TotalPotentialEnergy'

        OPEN( UNIT=NewUnit(fileUnits(4)), &
              FILE='Heat.'//timeStampString//'.curve', &
              FORM='FORMATTED', &
              STATUS='REPLACE' )
        WRITE(fileUnits(4),*) '#TotalHeat'

        OPEN( UNIT=NewUnit(fileUnits(5)), &
              FILE='Volume.'//timeStampString//'.curve', &
              FORM='FORMATTED', &
              STATUS='REPLACE' )
        WRITE(fileUnits(5),*) '#TotalVolume'
      ENDIF

 END SUBROUTINE OpenDiagnosticsFiles_Fluid
!
 SUBROUTINE WriteDiagnostics_Fluid( myDGSEM, fileUnits )
   IMPLICIT NONE
   CLASS( Fluid ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)        :: fileUnits(1:nDiagnostics)

      IF( myDGSEM % myRank == 0 )THEN
         WRITE(fileUnits(1),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % mass
         WRITE(fileUnits(2),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % KE
         WRITE(fileUnits(3),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % PE
         WRITE(fileUnits(4),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % heat
         WRITE(fileUnits(5),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % volume
      ENDIF

 END SUBROUTINE WriteDiagnostics_Fluid
!
 SUBROUTINE CloseDiagnosticsFiles_Fluid( myDGSEM, fileUnits )
   IMPLICIT NONE
   CLASS( Fluid  ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)            :: fileUnits(1:nDiagnostics) 
   

      IF( myDGSEM % myRank == 0 ) THEN
         CLOSE( UNIT=fileUnits(1) )
         CLOSE( UNIT=fileUnits(2) )
         CLOSE( UNIT=fileUnits(3) )
         CLOSE( UNIT=fileUnits(4) )
         CLOSE( UNIT=fileUnits(5) )
      ENDIF

 END SUBROUTINE CloseDiagnosticsFiles_Fluid
!
 SUBROUTINE Diagnostics_Fluid( myDGSEM )
  IMPLICIT NONE
  CLASS( Fluid ), INTENT(inout) :: myDGSEM 
  ! Local
  INTEGER    :: iEl, i, j, k 
  REAL(prec) :: volume, mass, KE, PE, heat
#ifdef HAVE_MPI
  INTEGER    :: mpiErr
#endif


     volume = 0.0_prec
     mass   = 0.0_prec
     KE     = 0.0_prec
     PE     = 0.0_prec
     heat   = 0.0_prec

     DO iEl = 1, myDGSEM % mesh % nElems
        DO k = 0, myDGSEM % N
           DO j = 0, myDGSEM % N
              DO i = 0, myDGSEM % N

                 volume = volume + myDGSEM % mesh % geom(iEl) % J(i,j,k)*&
                                             myDGSEM % dgStorage % quadratureWeights(i)*&
                                             myDGSEM % dgStorage % quadratureWeights(j)*&
                                             myDGSEM % dgStorage % quadratureWeights(k)

                 mass = mass + ( myDGSEM % state % solution(i,j,k,4,iEl)+&
                                           myDGSEM % static % solution(i,j,k,4,iEl) )*&
                                         myDGSEM % mesh % geom(iEl) % J(i,j,k)*&
                                             myDGSEM % dgStorage % quadratureWeights(i)*&
                                             myDGSEM % dgStorage % quadratureWeights(j)*&
                                             myDGSEM % dgStorage % quadratureWeights(k)

                 KE   = KE + ( myDGSEM % state % solution(i,j,k,1,iEl)**2 +&
                                         myDGSEM % state % solution(i,j,k,2,iEl)**2 +&
                                         myDGSEM % state % solution(i,j,k,3,iEl)**2 )/&
                                       ( myDGSEM % state % solution(i,j,k,4,iEl)+&
                                         myDGSEM % static % solution(i,j,k,4,iEl) )*&
                                       myDGSEM % mesh % geom(iEl) % J(i,j,k)*&
                                             myDGSEM % dgStorage % quadratureWeights(i)*&
                                             myDGSEM % dgStorage % quadratureWeights(j)*&
                                             myDGSEM % dgStorage % quadratureWeights(k)

                 PE   = PE - myDGSEM % state % solution(i,j,k,4,iEl)*&
                                       myDGSEM % params % g*&
                                       myDGSEM % mesh % geom(iEl) % z(i,j,k)*&
                                       myDGSEM % mesh % geom(iEl) % J(i,j,k)*&
                                             myDGSEM % dgStorage % quadratureWeights(i)*&
                                             myDGSEM % dgStorage % quadratureWeights(j)*&
                                             myDGSEM % dgStorage % quadratureWeights(k)

                 heat = heat + ( myDGSEM % static % solution(i,j,k,5,iEl) + &
                                 myDGSEM % state % solution(i,j,k,5,iEl) )/&
                               ( myDGSEM % static % solution(i,j,k,4,iEl) +&
                                 myDGSEM % state % solution(i,j,k,4,iEl) )*&
                                         myDGSEM % mesh % geom(iEl) % J(i,j,k)*& 
                                             myDGSEM % dgStorage % quadratureWeights(i)*&
                                             myDGSEM % dgStorage % quadratureWeights(j)*&
                                             myDGSEM % dgStorage % quadratureWeights(k)

              ENDDO
           ENDDO
        ENDDO
     ENDDO

     heat = heat*myDGSEM % params % Cv

     myDGSEM % volume = volume
     myDGSEM % mass   = mass
     myDGSEM % KE     = KE
     myDGSEM % PE     = PE
     myDGSEM % heat   = heat

#ifdef HAVE_MPI
     CALL MPI_ALLREDUCE( volume, myDGSEM % volume, 1, MPI_PREC, MPI_SUM, MPI_COMM_WORLD, mpiErr ) 
     CALL MPI_ALLREDUCE( mass, myDGSEM % mass, 1, MPI_PREC, MPI_SUM, MPI_COMM_WORLD, mpiErr ) 
     CALL MPI_ALLREDUCE( KE, myDGSEM % KE, 1, MPI_PREC, MPI_SUM, MPI_COMM_WORLD, mpiErr ) 
     CALL MPI_ALLREDUCE( PE, myDGSEM % PE, 1, MPI_PREC, MPI_SUM, MPI_COMM_WORLD, mpiErr ) 
     CALL MPI_ALLREDUCE( heat, myDGSEM % heat, 1, MPI_PREC, MPI_SUM, MPI_COMM_WORLD, mpiErr ) 
#endif

 END SUBROUTINE Diagnostics_Fluid
#endif
 SUBROUTINE ObtainPlottingMesh_Fluid( myDGSEM, x, y, z )

  IMPLICIT NONE
 
  CLASS( Fluid ), INTENT(inout) :: myDGsem
  REAL(prec), INTENT(out)       :: x(0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     1:myDGSEM % mesh % nElems)
  REAL(prec), INTENT(out)       :: y(0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     1:myDGSEM % mesh % nElems)
  REAL(prec), INTENT(out)       :: z(0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     1:myDGSEM % mesh % nElems)

  INTEGER       :: i, j, k, iEl
  

      
      DO iEl = 1, myDGsem % mesh % nElems

         x(:,:,:,iEl) = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % geom(iEl) % x )
         y(:,:,:,iEl) = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % geom(iEl) % y )
         z(:,:,:,iEl) = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % geom(iEl) % z )
           
      ENDDO

 END SUBROUTINE ObtainPlottingMesh_Fluid
!
 SUBROUTINE FluidStateAtPlottingPoints_Fluid( myDGSEM, u, v, w, density, potentialTemp, pressure )

  IMPLICIT NONE
 
  CLASS( Fluid ), INTENT(inout) :: myDGsem
  REAL(prec), INTENT(out)       :: u(0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     1:myDGSEM % mesh % nElems)
  REAL(prec), INTENT(out)       :: v(0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     1:myDGSEM % mesh % nElems)
  REAL(prec), INTENT(out)       :: w(0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     1:myDGSEM % mesh % nElems)
  REAL(prec), INTENT(out)       :: density(0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     1:myDGSEM % mesh % nElems)
  REAL(prec), INTENT(out)       :: potentialTemp(0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     1:myDGSEM % mesh % nElems)
  REAL(prec), INTENT(out)       :: pressure(0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     0:myDGSEM % params % nPlot, &
                                     1:myDGSEM % mesh % nElems)
  ! Local
  REAL(prec) :: sol(0:myDGSEM % params % nPlot, &
                    0:myDGSEM % params % nPlot, &
                    0:myDGSEM % params % nPlot, &
                    1:nEq )
  REAL(prec) :: bsol(0:myDGSEM % params % nPlot, &
                     0:myDGSEM % params % nPlot, &
                     0:myDGSEM % params % nPlot, &
                     1:nEq )
  INTEGER       :: i, j, k, iEl, iEq
  

 
      DO iEl = 1, myDGsem % mesh % nElems

           
         DO iEq = 1, nEq
            sol(:,:,:,iEq) = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % state % solution(:,:,:,iEq,iEl) )
         ENDDO
      
         DO iEq = 1, nEq
            bsol(:,:,:,iEq) = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % static % solution(:,:,:,iEq,iEl) )
         ENDDO
         

         DO k = 0, myDGSEM % params % nPlot
            DO j = 0, myDGSEM % params % nPlot
               DO i = 0, myDGSEM % params % nPlot
                          
                  u(i,j,k,iEl)             = sol(i,j,k,1)/( sol(i,j,k,4) + bsol(i,j,k,4) )
                  v(i,j,k,iEl)             = sol(i,j,k,2)/( sol(i,j,k,4) + bsol(i,j,k,4) )
                  w(i,j,k,iEl)             = sol(i,j,k,3)/( sol(i,j,k,4) + bsol(i,j,k,4) )
                  density(i,j,k,iEl)       = sol(i,j,k,4)
                  potentialTemp(i,j,k,iEl) = (sol(i,j,k,5) + bsol(i,j,k,5))/( sol(i,j,k,4) + bsol(i,j,k,4) )
                  pressure(i,j,k,iEl)      = sol(i,j,k,6)

               ENDDO
            ENDDO
         ENDDO
        
      ENDDO


 END SUBROUTINE FluidStateAtPlottingPoints_Fluid
!
 SUBROUTINE WriteSmoothedTecplot_Fluid( myDGSEM, iter, nPlot )

  IMPLICIT NONE
 
  CLASS( Fluid ), INTENT(inout) :: myDGsem
  INTEGER, INTENT(in)           :: iter, nPlot
  !LOCAL
  REAL(prec)  :: x(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: y(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: z(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: sol(0:nPlot,0:nPlot,0:nPlot,1:nEq)
  REAL(prec)  :: bsol(0:nPlot,0:nPlot,0:nPlot,1:nEq)

  INTEGER       :: i, j, k, iEl, iEq, fUnit
  CHARACTER(5)  :: zoneID
  CHARACTER(4)  :: rankChar
  CHARACTER(13) :: timeStampString
   
      timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )

      WRITE(rankChar,'(I4.4)') myDGSEM % myRank
      
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE= 'State-smoothed.'//rankChar//'.'//timeStampString//'.tec', &
            FORM='formatted', &
            STATUS='replace')  
    
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "u", "v", "w", "rho", "Pot. Temp.", "SGS KE"'
 
      DO iEl = 1, myDGsem % mesh % nElems

         x = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % geom(iEl) % x )
         y = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % geom(iEl) % y )
         z = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % geom(iEl) % z )
           
         DO iEq = 1, nEq
            sol(:,:,:,iEq) = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % smoothstate % solution(:,:,:,iEq,iEl) )
         ENDDO
         
         DO iEq = 1, nEq
            bsol(:,:,:,iEq) = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % static % solution(:,:,:,iEq,iEl) )
         ENDDO
         WRITE(zoneID,'(I5.5)') myDGSEM % mesh % elements(iEl) % elementID
         WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,&
                                                     ', J=',nPlot+1,&
                                                     ', K=',nPlot+1,',F=POINT'

         DO k = 0, nPlot
            DO j = 0, nPlot
               DO i = 0, nPlot
                  WRITE(fUnit,'(9(F15.7,1x))') x(i,j,k), y(i,j,k), z(i,j,k),&
                                  sol(i,j,k,1)/( sol(i,j,k,4) + bsol(i,j,k,4) ), &
                                  sol(i,j,k,2)/( sol(i,j,k,4) + bsol(i,j,k,4) ), &
                                  sol(i,j,k,3)/( sol(i,j,k,4) + bsol(i,j,k,4) ), &
                                  sol(i,j,k,4), &
                                  sol(i,j,k,5)/( sol(i,j,k,4) + bsol(i,j,k,4) ), &
                                  sol(i,j,k,6)
               ENDDO
            ENDDO
         ENDDO
        
      ENDDO

      CLOSE(UNIT=fUnit)

 END SUBROUTINE WriteSmoothedTecplot_Fluid
!
 SUBROUTINE WriteSGSTecplot_Fluid( myDGSEM, iter )

  IMPLICIT NONE
 
  CLASS( Fluid ), INTENT(inout) :: myDGsem
  INTEGER, INTENT(in)           :: iter
  !LOCAL
  INTEGER       :: i, j, k, iEl, iEq, fUnit
  CHARACTER(5)  :: zoneID
  CHARACTER(4)  :: rankChar
  CHARACTER(13) :: timeStampString
   
      timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )

      WRITE(rankChar,'(I4.4)') myDGSEM % myRank
      
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE= 'State-SGS.'//rankChar//'.'//timeStampString//'.tec', &
            FORM='formatted', &
            STATUS='replace')  
    
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "kappa_u", "kappa_v", "kappa_w", "kappa_rho", "kappa_T"'
 
      DO iEl = 1, myDGsem % mesh % nElems

         WRITE(zoneID,'(I5.5)') myDGSEM % mesh % elements(iEl) % elementID
         WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',myDGSEM % N+1,&
                                                     ', J=',myDGSEM % N+1,&
                                                     ', K=',myDGSEM % N+1,',F=POINT'

         DO k = 0, myDGSEM % N
            DO j = 0, myDGSEM % N
               DO i = 0, myDGSEM % N
                  WRITE(fUnit,'(8(F15.7,1x))') myDGSEM % mesh % geom(iEl) % x(i,j,k),&
                                               myDGSEM % mesh % geom(iEl) % y(i,j,k),&
                                               myDGSEM % mesh % geom(iEl) % z(i,j,k),&
                                               myDGSEM % sgsCoeffs % solution(i,j,k,1:5,iEl)
               ENDDO
            ENDDO
         ENDDO
        
      ENDDO

      CLOSE(UNIT=fUnit)

 END SUBROUTINE WriteSGSTecplot_Fluid
!
 SUBROUTINE WriteStressTensorTecplot_Fluid( myDGSEM, iter, nPlot )

  IMPLICIT NONE
 
  CLASS( Fluid ), INTENT(inout) :: myDGsem
  INTEGER, INTENT(in)           :: iter, nPlot
  !LOCAL
  REAL(prec)  :: x(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: y(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: z(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: sol(0:nPlot,0:nPlot,0:nPlot,1:15)

  INTEGER       :: i, j, k, iEl, iEq, fUnit
  CHARACTER(5)  :: zoneID
  CHARACTER(4)  :: rankChar
  CHARACTER(13) :: timeStampString
   
      timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )

      WRITE(rankChar,'(I4.4)') myDGSEM % myRank
      
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE= 'Stress.'//rankChar//'.'//timeStampString//'.tec', &
            FORM='formatted', &
            STATUS='replace')  
    
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "dmx/dx", "dmx/dy", "dmx/dz",'//&
                                                '"dmy/dx", "dmy/dy", "dmy/dz",'//&
                                                '"dmz/dx", "dmz/dy", "dmz/dz",'//&
                                                '"drho/dx", "drho/dy", "drho/dz",'//&
                                                '"drT/dx", "drT/dy", "drT/dz",'
      DO iEl = 1, myDGsem % mesh % nElems

         x = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % geom(iEl) % x )
         y = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % geom(iEl) % y )
         z = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % geom(iEl) % z )
           
         DO iEq = 1, (nEq-1)*3
            sol(:,:,:,iEq) = myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % stressTensor % solution(:,:,:,iEq,iEl) )
         ENDDO
      
         
         WRITE(zoneID,'(I5.5)') myDGSEM % mesh % elements(iEl) % elementID
         WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,&
                                                     ', J=',nPlot+1,&
                                                     ', K=',nPlot+1,',F=POINT'

         DO k = 0, nPlot
            DO j = 0, nPlot
               DO i = 0, nPlot
                  WRITE(fUnit,'(18(F15.7,1x))') x(i,j,k), y(i,j,k), z(i,j,k),sol(i,j,k,1:15)
               ENDDO
            ENDDO
         ENDDO
        
      ENDDO

      CLOSE(UNIT=fUnit)

 END SUBROUTINE WriteStressTensorTecplot_Fluid
!
 SUBROUTINE WritePickup_Fluid( myDGSEM )

   IMPLICIT NONE
   CLASS( Fluid ), INTENT(in) :: myDGSEM
  ! LOCAL
   CHARACTER(4)  :: rankChar
   INTEGER       :: iEl
   INTEGER       :: thisRec, fUnit
   INTEGER       :: iEq, N
   CHARACTER(13) :: timeStampString
   
      timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )

      N = myDGSEM % N
     
      WRITE(rankChar,'(I4.4)') myDGSEM % myRank
      PRINT(MsgFMT), ' S/R WritePickup_Fluid : Writing output file :  State.'//&
                       rankChar//'.'//timeStampString//'.pickup' 

      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE='State.'//rankChar//'.'//timeStampString//'.pickup', &
            FORM='UNFORMATTED',&
            ACCESS='DIRECT',&
            STATUS='REPLACE',&
            ACTION='WRITE',&
            CONVERT='BIG_ENDIAN',&
            RECL=prec*(N+1)*(N+1)*(N+1) )

      thisRec = 1 
      DO iEl = 1, myDGSEM % mesh % nElems
        
         DO iEq = 1, nEq
            WRITE( fUnit, REC=thisRec )myDGSEM % state % solution(:,:,:,iEq,iEl)
            thisRec = thisRec+1
         ENDDO
         DO iEq = 1, nEq
            WRITE( fUnit, REC=thisRec )myDGSEM % static % solution(:,:,:,iEq,iEl) 
            thisRec = thisRec+1
         ENDDO
        
         WRITE( fUnit, REC=thisRec )myDGSEM % dragProfile(:,:,:,iEl)
         thisRec = thisRec+1
         
      ENDDO

      CLOSE(UNIT=fUnit)
     
      OPEN( UNIT   = NEWUNIT(fUnit), &
            FILE   = 'State.'//rankChar//'.'//timeStampString//'.exs', &
            FORM   ='UNFORMATTED',&
            ACCESS ='DIRECT',&
            STATUS ='REPLACE',&
            ACTION ='WRITE', &
            CONVERT='BIG_ENDIAN',&
            RECL   = prec*(myDGSEM % N+1)*(myDGSEM % N+1)*(nEq)*(myDGSEM % nBoundaryFaces) )
      WRITE( fUnit, rec = 1 ) myDGSEM % externalState
      WRITE( fUnit, rec = 2 ) myDGSEM % prescribedState
      CLOSE(fUnit) 


 END SUBROUTINE WritePickup_Fluid
!
 SUBROUTINE ReadPickup_Fluid( myDGSEM )

   IMPLICIT NONE
   CLASS( Fluid ), INTENT(inout) :: myDGSEM
  ! LOCAL
   CHARACTER(4)  :: rankChar
   INTEGER       :: iEl, istat
   INTEGER       :: thisRec, fUnit
   INTEGER       :: iEq, N
   LOGICAL       :: itExists
   CHARACTER(13) :: timeStampString
   
      timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )
   
      N = myDGSEM % N
      
      WRITE(rankChar,'(I4.4)') myDGSEM % myRank
      INQUIRE( FILE='State.'//rankChar//'.'//timeStampString//'.pickup', EXIST = itExists )
     
      IF( itExists )THEN
      
         PRINT*, '  Opening State.'//rankChar//'.'//timeStampString//'.pickup'
     
         OPEN( UNIT=NEWUNIT(fUnit), &
               FILE='State.'//rankChar//'.'//timeStampString//'.pickup', &
               FORM='unformatted',&
               ACCESS='direct',&
               STATUS='old',&
               ACTION='READ',&
               CONVERT='big_endian',&
               RECL=prec*(N+1)*(N+1)*(N+1) )
        
         thisRec = 1
         DO iEl = 1, myDGSEM % mesh % nElems
           
            DO iEq = 1, nEq
               READ( fUnit, REC=thisRec )myDGSEM % state % solution(:,:,:,iEq,iEl) 
               thisRec = thisRec+1
            ENDDO
            DO iEq = 1, nEq
               READ( fUnit, REC=thisRec )myDGSEM % static % solution(:,:,:,iEq,iEl)
               thisRec = thisRec+1
            ENDDO
            
            READ( fUnit, REC=thisRec )myDGSEM % dragProfile(:,:,:,iEl)
            thisRec = thisRec+1
         
         ENDDO
         
         CLOSE(UNIT=fUnit)
         
      ENDIF

#ifdef HAVE_CUDA
      myDGSEM % dragProfile_dev       = myDGSEM % dragProfile
      myDGSEM % state % solution_dev  = myDGSEM % state % solution
      myDGSEM % static % solution_dev = myDGSEM % static % solution
#endif
      
	  CALL myDGSEM % extComm % ReadPickup( 'ExtComm.'//rankChar )
	  myDGSEM % nBoundaryFaces = myDGSEM % extComm % nBoundaries
	  
	  ALLOCATE( myDGSEM % externalState(0:N,0:N,1:nEq,1:myDGSEM % nBoundaryFaces) )
	  ALLOCATE( myDGSEM % externalStress(0:N,0:N,1:(nEq-1)*3,1:myDGSEM % nBoundaryFaces) )
	  ALLOCATE( myDGSEM % externalSGS(0:N,0:N,1:nEq-1,1:myDGSEM % nBoundaryFaces) )
	  ALLOCATE( myDGSEM % prescribedState(0:N,0:N,1:nEq,1:myDGSEM % nBoundaryFaces) )
	  ALLOCATE( myDGSEM % prescribedStress(0:N,0:N,1:(nEq-1)*3,1:myDGSEM % nBoundaryFaces) )
	  myDGSEM % externalState    = 0.0_prec
	  myDGSEM % externalStress   = 0.0_prec
	  myDGSEM % externalSGS      = 0.0_prec
	  myDGSEM % prescribedState  = 0.0_prec
	  myDGSEM % prescribedStress = 0.0_prec

#ifdef HAVE_CUDA	 
	  ALLOCATE( myDGSEM % externalState_dev(0:N,0:N,1:nEq,1:myDGSEM % nBoundaryFaces) )
	  ALLOCATE( myDGSEM % externalStress_dev(0:N,0:N,1:(nEq-1)*3,1:myDGSEM % nBoundaryFaces) )
	  ALLOCATE( myDGSEM % externalSGS_dev(0:N,0:N,1:nEq-1,1:myDGSEM % nBoundaryFaces) )
	  ALLOCATE( myDGSEM % prescribedState_dev(0:N,0:N,1:nEq,1:myDGSEM % nBoundaryFaces) )
	  ALLOCATE( myDGSEM % prescribedStress_dev(0:N,0:N,1:(nEq-1)*3,1:myDGSEM % nBoundaryFaces) )
	  myDGSEM % externalState_dev    = 0.0_prec
	  myDGSEM % externalStress_dev   = 0.0_prec
	  myDGSEM % externalSGS_dev      = 0.0_prec
	  myDGSEM % prescribedState_dev  = 0.0_prec
	  myDGSEM % prescribedStress_dev = 0.0_prec
#endif	  
     
      INQUIRE( FILE='State.'//rankChar//'.'//timeStampString//'.exs', EXIST = itExists )
     
      IF( itExists )THEN
      
         OPEN( UNIT   = NEWUNIT(fUnit), &
                FILE   = 'State.'//rankChar//'.'//timeStampString//'.exs', &
                FORM   ='UNFORMATTED',&
                ACCESS ='DIRECT',&
                STATUS ='OLD',&
                ACTION ='READ', &
                CONVERT='BIG_ENDIAN',&
                RECL   = prec*(myDGSEM % N+1)*(myDGSEM % N+1)*(nEq)*(myDGSEM % nBoundaryFaces) )
         READ( fUnit, rec = 1 ) myDGSEM % externalState
         READ( fUnit, rec = 2 ) myDGSEM % prescribedState
         CLOSE(fUnit) 
         
#ifdef HAVE_CUDA
         ! Copy the external state and the prescribed-state to the device
         myDGSEM % externalState_dev   = myDGSEM % externalState
         myDGSEM % prescribedState_dev = myDGSEM % prescribedState
#endif
         
      ENDIF
      
      PRINT*, 'S/R ReadPickup : Done.'
      
      ! Interpolate the static state to the element boundaries
      CALL myDGSEM % CalculateStaticBoundarySolution( )  
      
 END SUBROUTINE ReadPickup_Fluid
!
#ifdef HAVE_CUDA
! ============================================================================================================================ !
!------------------------------------------- CUDA Kernels Below -------------------------------------------------------------- !
! ============================================================================================================================ !
 ATTRIBUTES(Global) SUBROUTINE UpdateG3D_CUDAKernel( G3D, a, g, solution, tendency, diffusiveTendency )
    
   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(inout) :: G3D(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)    :: a, g
   REAL(prec), DEVICE, INTENT(inout) :: solution(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)    :: tendency(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)    :: diffusivetendency(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:15,1:nEl_dev)
   ! Local
   INTEGER :: i, j, k, iEq, iEl
   
      iEl = blockIDx % x
      iEq = blockIDx % y
      
      i = threadIdx % x - 1
      j = threadIdx % y - 1
      k = threadIdx % z - 1
      
      G3D(i,j,k,iEq,iEl)      = a*G3D(i,j,k,iEq,iEl) + tendency(i,j,k,iEq,iEl) + diffusivetendency(i,j,k,iEq,iEl) 
      solution(i,j,k,iEq,iEl) = solution(i,j,k,iEq,iEl) + dt_dev*g*G3D(i,j,k,iEq,iEl)

 END SUBROUTINE UpdateG3D_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE CalculateSmoothedState_CUDAKernel( solution, filterMat, smoothSolution )
 
   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(in)    :: solution(0:polydeg_dev, 0:polydeg_dev, 0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)    :: filterMat(0:polydeg_dev, 0:polydeg_dev)
   REAL(prec), DEVICE, INTENT(inout) :: smoothsolution(0:polydeg_dev, 0:polydeg_dev, 0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   ! Local
   INTEGER :: iEl, iEq, i, j, k, ii, jj, kk
   REAL(prec) :: uijk, uij, ui
   
      iEl = blockIdx % x
      iEq = blockIdx % y
      
      i = threadIdx % x-1
      j = threadIdx % y-1
      k = threadIdx % z-1
   
	  uijk = 0.0_prec
	  DO kk = 0, polydeg_dev
	 
		 uij = 0.0_prec
		 DO jj = 0, polydeg_dev
		   
		    ui = 0.0_prec
		    DO ii = 0, polydeg_dev
			   ui = ui + filterMat(ii,i)*solution(ii,jj,kk,iEq,iEl)
		    ENDDO
		   
		    uij = uij + filterMat(jj,j)*ui
		 ENDDO
		
 		 uijk = uijk + filterMat(kk,k)*uij
		
      ENDDO
	 
	  smoothSolution(i,j,k,iEq,iEl) = uijk
                     
 END SUBROUTINE CalculateSmoothedState_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE CalculateSGSCoefficients_CUDAKernel( solution, static, smoothState, filterMat, sgsCoeffs )
 
   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(in)    :: solution(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)    :: static(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(inout) :: smoothState(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)    :: filterMat(0:polydeg_dev,0:polydeg_dev)
   REAL(prec), DEVICE, INTENT(inout) :: sgsCoeffs(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev-1,1:nEl_dev)
   ! Local
   INTEGER :: iEl, i, j, k, m, ii, jj, kk
   REAL(prec) :: sgsKE, uijk, uij, ui
#ifdef VIZ
   REAL(prec), SHARED :: KE(0:7,0:7,0:7) 
#endif

      iEl = blockIDx % x
      
      i = threadIdx % x-1
      j = threadIdx % y-1
      k = threadIdx % z-1
      
      ! Here, the SGS Kinetic energy is calculated using the 
      ! "high wavenumber" component of the velocity field.
      ! This component is defined (here) as the difference
      ! between the full solution and the smoothed solution.
     
      sgsKE = 0.0_prec
      DO m = 1, 3  
         sgsKE = sgsKE + &
        ( solution(i,j,k,m,iEl)/( solution(i,j,k,4,iEl) + static(i,j,k,4,iEl))- &
        smoothState(i,j,k,m,iEl)/(smoothState(i,j,k,4,iEl)+static(i,j,k,4,iEl)) )**2
      ENDDO

#ifdef VIZ
      KE(i,j,k) = 0.5_prec*sgsKE
                  
      CALL syncthreads( )
      
      ! Smooth the subgrid scale Kinetic energy
      uijk = 0.0_prec
      DO kk = 0, polydeg_dev
 
         uij = 0.0_prec
         DO jj = 0, polydeg_dev

            ui = 0.0_prec
            DO ii = 0, polydeg_dev
               ui = ui + filterMat(ii,i)*KE(ii,jj,kk)
            ENDDO

            uij = uij + filterMat(jj,j)*ui
         ENDDO
 
         uijk = uijk + filterMat(kk,k)*uij
 
      ENDDO
                  
      ! Here, we store the smoothed SGS kinetic energy, in
      ! case we would like to visualize the data later
      smoothState(i,j,k,6,iEl) = ABS(uijk)
#endif
 
         ! Now we calculate the viscosity and diffusivities (currently assumes isotropic and low mach number)
      DO m = 1, nEq_dev-1
         !IF( m == 4 )THEN
         !   sgsCoeffs(i,j,k,m,iEl) = 0.0_prec ! No density diffusion
         !ELSE
            ! This is the parameterization used in Jeremy Sauer's dissertation ... citation ?!
            !** Note that the filtering process may not preserve positivity of the EKE.. hence 
            !   we need to take the absolute value of uijk
            sgsCoeffs(i,j,k,m,iEl) = 0.09_prec*viscLengthScale_dev*sqrt( sgsKE )
         ! ENDIF
      ENDDO
 
 END SUBROUTINE CalculateSGSCoefficients_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE CalculateBoundarySGS_CUDAKernel( sgsCoeffs, bMat, boundarySGSCoeffs ) 

   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(in)  :: sgsCoeffs(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev-1,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: bMat(0:polydeg_dev,0:1)
   REAL(prec), DEVICE, INTENT(out) :: boundarysgsCoeffs(0:polydeg_dev,0:polydeg_dev,1:nEq_dev-1,1:6,1:nEl_dev)
   ! Local
   INTEGER :: iEq, iEl, i, j, k
   REAL(prec) :: bSol(1:6)

      iEl = blockIdx % x
      
      iEq = threadIdx % z
      k   = threadIdx % y-1
      j   = threadIdx % x-1
      
	  bSol(1:6) = 0.0_prec

	  DO i = 0, polydeg_dev
		 bSol(1) = bSol(1) + bMat(i,0)*sgsCoeffs(j,i,k,iEq,iEl) ! south
			
		 bSol(2) = bSol(2) + bMat(i,1)*sgsCoeffs(i,j,k,iEq,iEl) ! east
			
		 bSol(3) = bSol(3) + bMat(i,1)*sgsCoeffs(j,i,k,iEq,iEl) ! north
															
		 bSol(4) = bSol(4) + bMat(i,0)*sgsCoeffs(i,j,k,iEq,iEl) ! west

		 bSol(5) = bSol(5) + bMat(i,0)*sgsCoeffs(j,k,i,iEq,iEl) ! botom

		 bSol(6) = bSol(6) + bMat(i,1)*sgsCoeffs(j,k,i,iEq,iEl) ! top
	  ENDDO
               
      DO i = 1, 6
         boundarysgsCoeffs(j,k,iEq,i,iEl) = ABS(bSol(i)) !Ensure positivity of the viscosity/diffusivity coefficients
      ENDDO
      
 END SUBROUTINE CalculateBoundarySGS_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE UpdateExternalSGSCoeffs_CUDAKernel( boundaryIDs, elementIDs, &
                                                               elementSides, procIDs, &
                                                               externalsgsCoeffs, &
                                                               sgsCoeffsBsols,nHat ) 
   IMPLICIT NONE
   INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nBoundaryFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: procIDs(1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(out) :: externalsgsCoeffs(0:polydeg_dev,0:polydeg_dev,1:nEq_dev-1,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: sgsCoeffsBsols(0:polydeg_dev,0:polydeg_dev,1:nEq_dev-1,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: nhat(1:3,0:polydeg_dev,0:polydeg_dev,1:6,1:nEl_dev)
   ! Local
   INTEGER    :: iEq, iFace, i, j, k
   INTEGER    :: iFace2, p2
   INTEGER    :: e1, e2, s1, s2, m
   
      iFace = blockIdx % x
      iEq   = blockIDx % y
      ! ////////////////////////////////////////////////////////////////////////// !
      i   = threadIdx % x-1
      j   = threadIdx % y-1

	  iFace2 = boundaryIDs( iFace ) ! Obtain the process-local face id for this boundary-face id
	  e1     = elementIDs(1,iFace2)
	  s1     = elementSides(1,iFace2)
	  e2     = elementIDs(2,iFace2)
	  p2     = procIDs( iFace )
	 
	  IF( i <= polydeg_dev .AND. j <= polydeg_dev )THEN
	 
		 IF( p2 == myRank_dev )THEN
		    externalsgsCoeffs(i,j,iEq,iFace) = sgsCoeffsBsols(i,j,iEq,s1,e1)
		 ENDIF
		
	 ENDIF
                  
 END SUBROUTINE UpdateExternalSGSCoeffs_CUDAKernel  
!
 ATTRIBUTES(Global) SUBROUTINE CalculateBoundarySolution_CUDAKernel( solution, bMat, boundarySolution ) 

   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(in)  :: solution(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: bMat(0:polydeg_dev,0:1)
   REAL(prec), DEVICE, INTENT(out) :: boundarySolution(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   ! Local
   INTEGER :: iEq, iEl, i, j, k
   REAL(prec) :: bSol(1:6)

      iEl = blockIdx % x
      
      iEq = threadIdx % z
      k   = threadIdx % y-1
      j   = threadIdx % x-1
      
	  bSol(1:6) = 0.0_prec

	  DO i = 0, polydeg_dev
		 bSol(1) = bSol(1) + bMat(i,0)*solution(j,i,k,iEq,iEl) ! south
			
		 bSol(2) = bSol(2) + bMat(i,1)*solution(i,j,k,iEq,iEl) ! east
			
		 bSol(3) = bSol(3) + bMat(i,1)*solution(j,i,k,iEq,iEl) ! north
															
		 bSol(4) = bSol(4) + bMat(i,0)*solution(i,j,k,iEq,iEl) ! west

		 bSol(5) = bSol(5) + bMat(i,0)*solution(j,k,i,iEq,iEl) ! botom

		 bSol(6) = bSol(6) + bMat(i,1)*solution(j,k,i,iEq,iEl) ! top
	  ENDDO
               
      DO i = 1, 6
         boundarySolution(j,k,iEq,i,iEl) = bSol(i)
      ENDDO
      
 END SUBROUTINE CalculateBoundarySolution_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE UpdateExternalState_CUDAKernel( boundaryIDs, elementIDs, &
                                                               elementSides, procIDs, &
                                                               externalState, &
                                                               stateBsols, &
                                                               prescribedState, nHat ) 
   IMPLICIT NONE
   INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nBoundaryFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: procIDs(1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(out) :: externalState(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: stateBsols(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: prescribedState(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: nhat(1:3,0:polydeg_dev,0:polydeg_dev,1:6,1:nEl_dev)
   ! Local
   INTEGER    :: iEl, iFace, bFaceID, i, j, k, iEq
   INTEGER    :: iFace2, p2
   INTEGER    :: e1, e2, s1, s2
   REAL(prec) :: norm, un, ut, us, speed
   REAL(prec) :: nx, ny, nz
   REAL(prec) :: sx, sy, sz
   REAL(prec) :: tx, ty, tz
   
      iFace = blockIdx % x
      i     = threadIdx % x-1
      j     = threadIdx % y-1

      IF( iFace <= nBoundaryFaces_dev )THEN
      
         iFace2 = boundaryIDs( iFace ) ! Obtain the process-local face id for this boundary-face id
         e1     = elementIDs(1,iFace2)
         s1     = elementSides(1,iFace2)
         e2     = elementIDs(2,iFace2)
         p2     = procIDs( iFace )
         
         IF( i <= polydeg_dev .AND. j <= polydeg_dev .AND. p2 == myRank_dev)THEN
         
            IF( e2 == PRESCRIBED )THEN
               
               DO iEq = 1, nEq_dev
                  externalState(i,j,iEq,iFace) = prescribedState(i,j,iEq,iFace)
               ENDDO
                  
            ELSEIF( e2 == RADIATION )THEN
                        
               DO iEq = 1, nEq_dev
                  externalState(i,j,iEq,iFace) = 0.0_prec
               ENDDO
                 
            ELSEIF( e2 == NO_NORMAL_FLOW )THEN
                             
               ! normal
               nx = nHat(1,i,j,s1,e1) !**
               ny = nHat(2,i,j,s1,e1)
               nz = nHat(3,i,j,s1,e1)
               norm = sqrt( nx*nx + ny*ny + nz*nz )
               nx = nx/norm
               ny = ny/norm
               nz = nz/norm
      
               ! tangent (built by performing 90 deg rotation in y - if zero, performs rotation in x)
               IF( nz == 0.0_prec .AND. ny == 0.0_prec )THEN ! rotate about y-axis
                  sx = -nz
                  sy = 0.0_prec
                  sz = nx
               ELSE
                  sx = 0.0_prec
                  sy = nz
                  sz = -ny
               ENDIF
  
               norm = sqrt( sx*sx + sy*sy + sz*sz )
               sx = sx/norm
               sy = sy/norm
               sz = sz/norm
                       
               !binormal
               tx = sy*nz - sz*ny
               ty = nx*sz - nz*sx
               tz = sx*ny - nx*sy
               norm = sqrt( tx*tx + ty*ty + tz*tz )
               tx = tx/norm
               ty = ty/norm
               tz = tz/norm
                        
               un = stateBsols(i,j,1,s1,e1)*nx + &
                    stateBsols(i,j,2,s1,e1)*ny + &
                    stateBsols(i,j,3,s1,e1)*nz
               us = stateBsols(i,j,1,s1,e1)*sx    + &
                    stateBsols(i,j,2,s1,e1)*sy    + &
                    stateBsols(i,j,3,s1,e1)*sz
               ut = stateBsols(i,j,1,s1,e1)*tx  + &
                    stateBsols(i,j,2,s1,e1)*ty  + &
                    stateBsols(i,j,3,s1,e1)*tz
                        
               externalState(i,j,1,iFace) = -nx*un + us*sx + ut*tx ! u
               externalState(i,j,2,iFace) = -ny*un + us*sy + ut*ty ! v
               externalState(i,j,3,iFace) = -nz*un + us*sz + ut*tz ! w
               externalState(i,j,4,iFace) =  stateBsols(i,j,4,s1,e1) ! rho
               externalState(i,j,5,iFace) =  stateBsols(i,j,5,s1,e1) ! potential temperature
               externalState(i,j,6,iFace) =  stateBsols(i,j,6,s1,e1) ! P
                        
            ELSEIF( e2 == DRAG_SLIP )THEN
                             
               ! normal
               nx = nHat(1,i,j,s1,e1) !**
               ny = nHat(2,i,j,s1,e1)
               nz = nHat(3,i,j,s1,e1)
               norm = sqrt( nx*nx + ny*ny + nz*nz )
               nx = nx/norm
               ny = ny/norm
               nz = nz/norm
      
               ! tangent (built by performing 90 deg rotation in y - if zero, performs rotation in x)
               IF( nz == 0.0_prec .AND. ny == 0.0_prec )THEN ! rotate about y-axis
                  sx = -nz
                  sy = 0.0_prec
                  sz = nx
               ELSE
                  sx = 0.0_prec
                  sy = nz
                  sz = -ny
               ENDIF
  
               norm = sqrt( sx*sx + sy*sy + sz*sz )
               sx = sx/norm
               sy = sy/norm
               sz = sz/norm
                       
               !binormal
               tx = sy*nz - sz*ny
               ty = nx*sz - nz*sx
               tz = sx*ny - nx*sy
               norm = sqrt( tx*tx + ty*ty + tz*tz )
               tx = tx/norm
               ty = ty/norm
               tz = tz/norm
                       
               speed = ( stateBsols(i,j,1,s1,e1)**2 +&
                         stateBsols(i,j,2,s1,e1)**2 +&
                         stateBsols(i,j,3,s1,e1)**2 )/&
                         stateBsols(i,j,4,s1,e1)
               
               un = stateBsols(i,j,1,s1,e1)*nx + &
                    stateBsols(i,j,2,s1,e1)*ny + &
                    stateBsols(i,j,3,s1,e1)*nz
                    
               us = ( stateBsols(i,j,1,s1,e1)*sx + &
                      stateBsols(i,j,2,s1,e1)*sy + &
                      stateBsols(i,j,3,s1,e1)*sz )*&
                      (1.0_prec-Cd_dev*dScale_dev*speed)
                      
               ut = ( stateBsols(i,j,1,s1,e1)*tx + &
                      stateBsols(i,j,2,s1,e1)*ty + &
                      stateBsols(i,j,3,s1,e1)*tz )*&
                      (1.0_prec-Cd_dev*dScale_dev*speed)
                        
               externalState(i,j,1,iFace) = -nx*un + us*sx + ut*tx ! u
               externalState(i,j,2,iFace) = -ny*un + us*sy + ut*ty ! v
               externalState(i,j,3,iFace) = -nz*un + us*sz + ut*tz ! w
               externalState(i,j,4,iFace) =  stateBsols(i,j,4,s1,e1) ! rho
               externalState(i,j,5,iFace) =  stateBsols(i,j,5,s1,e1) ! potential temperature
               externalState(i,j,6,iFace) =  stateBsols(i,j,6,s1,e1) ! P
                        
                     
            ENDIF
            
         ENDIF
                  
      ENDIF 
     

 END SUBROUTINE UpdateExternalState_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE InternalFaceFlux_CUDAKernel( elementIDs, elementSides, boundaryIDs, iMap, jMap, &
                                                 nHat, boundarySolution, boundarySolution_static, &
                                                 externalState, boundaryFlux, stressFlux )

   IMPLICIT NONE
   INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: iMap(0:polydeg_dev,0:polydeg_dev,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: jMap(0:polydeg_dev,0:polydeg_dev,1:nFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:polydeg_dev,0:polydeg_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundarySolution(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundarySolution_static(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: externalState(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(out) :: boundaryFlux(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(out) :: stressFlux(0:polydeg_dev,0:polydeg_dev,1:15,1:6,1:nEl_dev)
   ! Local
   INTEGER    :: iEl, iFace, jEq
   INTEGER    :: i, j, k, iEq
   INTEGER    :: ii, jj, bID
   INTEGER    :: e1, s1, e2, s2
   REAL(prec) :: uOut, uIn, cIn, cOut, norm, T
   REAL(prec) :: jump(1:5), aS(1:5)
   REAL(prec) :: fac


      iFace = blockIdx % x
      j     = threadIdx % y - 1
      i     = threadIdx % x -1
     
         e1 = elementIDs(1,iFace)
         s1 = elementSides(1,iFace)
         e2 = elementIDs(2,iFace)
         s2 = ABS(elementSides(2,iFace))
         bID  = ABS(boundaryIDs(iFace))

               ii = iMap(i,j,iFace)
               jj = jMap(i,j,iFace)
               
               norm = sqrt( nHat(1,i,j,s1,e1)*nHat(1,i,j,s1,e1) + &
                            nHat(2,i,j,s1,e1)*nHat(2,i,j,s1,e1) + &
                            nHat(3,i,j,s1,e1)*nHat(3,i,j,s1,e1) )

               
               IF( e2 > 0 )THEN
               
                  DO iEq = 1, nEq_dev-1
                     jump(iEq)  = boundarySolution(ii,jj,iEq,s2,e2) - &
                                  boundarySolution(i,j,iEq,s1,e1) !outState - inState
                  ENDDO

                  
                  T =   (boundarySolution_static(ii,jj,5,s2,e2) + boundarySolution(ii,jj,5,s2,e2))/&
                          (boundarySolution(ii,jj,4,s2,e2)+boundarySolution_static(ii,jj,4,s2,e2) )
                          
                  ! Sound speed estimate for the external and internal states
                  cOut = sqrt( R_dev*T* &
                              ( (boundarySolution(ii,jj,6,s2,e2)+boundarySolution_static(ii,jj,6,s2,e2))/ P0_dev )**rC_dev   )
                        
                  T =   (boundarySolution_static(i,j,5,s1,e1) + boundarySolution(i,j,5,s1,e1))/&
                          (boundarySolution(i,j,4,s1,e1)+boundarySolution_static(i,j,4,s1,e1) )        
                             
                  cIn  = sqrt( R_dev*T* &
                              ( (boundarySolution(i,j,6,s1,e1)+boundarySolution_static(i,j,6,s1,e1))/P0_dev )**rC_dev  )
                               
                  ! External normal velocity component
                  uOut = ( boundarySolution(ii,jj,1,s2,e2)*nHat(1,i,j,s1,e1)/norm + &
                           boundarySolution(ii,jj,2,s2,e2)*nHat(2,i,j,s1,e1)/norm + &
                           boundarySolution(ii,jj,3,s2,e2)*nHat(3,i,j,s1,e1)/norm )/& 
                         ( boundarySolution(ii,jj,4,s2,e2) + boundarySolution_static(ii,jj,4,s2,e2) )
                           
                  ! Internal normal velocity component
                  uIn  = ( boundarySolution(i,j,1,s1,e1)*nHat(1,i,j,s1,e1)/norm + &
                           boundarySolution(i,j,2,s1,e1)*nHat(2,i,j,s1,e1)/norm + &
                           boundarySolution(i,j,3,s1,e1)*nHat(3,i,j,s1,e1)/norm )/& 
                         ( boundarySolution(i,j,4,s1,e1) + boundarySolution_static(i,j,4,s1,e1) ) 
                           
                  ! Lax-Friedrich's estimate of the magnitude of the flux jacobian matrix
                  fac = max( abs(uIn+cIn), abs(uIn-cIn), abs(uOut+cOut), abs(uOut-cOut) )
                  !fac = max( abs(uIn),  abs(uOut) )

                  ! Advective flux
                  DO iEq = 1, nEq_dev-1
                        aS(iEq) = uIn*( boundarySolution(i,j,iEq,s1,e1) + boundarySolution_static(i,j,iEq,s1,e1) ) +&
                                 uOut*( boundarySolution(ii,jj,iEq,s2,e2) + boundarySolution_static(ii,jj,iEq,s2,e2) )
                  ENDDO
                  
                  DO k = 1, 3
                  ! Momentum flux due to pressure
                  aS(k) = aS(k) + (boundarySolution(i,j,6,s1,e1) + &
                                   boundarySolution(ii,jj,6,s2,e2))*nHat(k,i,j,s1,e1)/norm
                  ENDDO    
      
                         
                  DO iEq = 1, nEq_dev-1
                     boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*( aS(iEq) - fac*jump(iEq) )*norm
                     boundaryFlux(ii,jj,iEq,s2,e2) = -boundaryFlux(i,j,iEq,s1,e1)
                     IF( iEq == 4 )THEN
                        DO k = 1, 3
                           jEq = k+(iEq-1)*3
                           ! Calculate the LDG flux for the stress tensor.
                           stressFlux(i,j,jEq,s1,e1) = 0.5_prec*( boundarySolution(i,j,iEq,s1,e1) +&
                                                                  boundarySolution(ii,jj,iEq,s2,e2))*& 
                                                                  nHat(k,i,j,s1,e1)
                                                                                        
                           stressFlux(ii,jj,jEq,s2,e2) = -stressFlux(i,j,jEq,s1,e1)
                        ENDDO
                     ELSE
                        DO k = 1, 3
                           jEq = k+(iEq-1)*3
                           ! Calculate the LDG flux for the stress tensor.
                           stressFlux(i,j,jEq,s1,e1) = 0.5_prec*( boundarySolution(i,j,iEq,s1,e1)/&
                                                                  (boundarySolution(i,j,4,s1,e1)+&
                                                                   boundarySolution_static(i,j,4,s1,e1)) +&
                                                                  boundarySolution(ii,jj,iEq,s2,e2)/&
                                                                  (boundarySolution(ii,jj,4,s2,e2)+&
                                                                   boundarySolution_static(ii,jj,4,s2,e2)) )*& 
                                                                  nHat(k,i,j,s1,e1)
                                                                                        
                           stressFlux(ii,jj,jEq,s2,e2) = -stressFlux(i,j,jEq,s1,e1)
                        ENDDO
                     ENDIF
                     
                  ENDDO
             
               ENDIF 


 END SUBROUTINE InternalFaceFlux_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE BoundaryFaceFlux_CUDAKernel( elementIDs, elementSides, boundaryIDs, iMap, jMap, &
                                                 nHat, boundarySolution, boundarySolution_static, &
                                                 externalState, boundaryFlux, stressFlux )

   IMPLICIT NONE
   INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: iMap(0:polydeg_dev,0:polydeg_dev,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: jMap(0:polydeg_dev,0:polydeg_dev,1:nFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:polydeg_dev,0:polydeg_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundarySolution(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundarySolution_static(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: externalState(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(out) :: boundaryFlux(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(out) :: stressFlux(0:polydeg_dev,0:polydeg_dev,1:15,1:6,1:nEl_dev)
   ! Local
   INTEGER    :: iEl, iFace, jEq
   INTEGER    :: i, j, k, iEq
   INTEGER    :: ii, jj, bID
   INTEGER    :: e1, s1, e2, s2
   REAL(prec) :: uOut, uIn, cIn, cOut, norm, T
   REAL(prec) :: jump(1:5), aS(1:5)
   REAL(prec) :: fac


      iFace = blockIdx % x
      j     = threadIdx % y - 1
      i     = threadIdx % x -1
     
         e1 = elementIDs(1,iFace)
         s1 = elementSides(1,iFace)
         e2 = elementIDs(2,iFace)
         s2 = ABS(elementSides(2,iFace))
         bID  = ABS(boundaryIDs(iFace))

               ii = iMap(i,j,iFace)
               jj = jMap(i,j,iFace)
               
               norm = sqrt( nHat(1,i,j,s1,e1)*nHat(1,i,j,s1,e1) + &
                            nHat(2,i,j,s1,e1)*nHat(2,i,j,s1,e1) + &
                            nHat(3,i,j,s1,e1)*nHat(3,i,j,s1,e1) )

               
               IF( e2 < 0 )THEN
               
                  
                  DO iEq = 1, nEq_dev-1              
                  jump(iEq)  = externalState(ii,jj,iEq,bID) - &
                               boundarySolution(i,j,iEq,s1,e1) !outState - inState
                  ENDDO
                 
                  T =   (boundarySolution_static(i,j,5,s1,e1) + externalState(ii,jj,5,bID))/&
                          (externalState(ii,jj,4,bID)+boundarySolution_static(i,j,4,s1,e1) )
                 ! Sound speed estimate for the external and internal states
                  cOut = sqrt( R_dev*T* &
                              ( (externalState(ii,jj,6,bID)+boundarySolution_static(i,j,6,s1,e1))/ P0_dev )**rC_dev   )
                  
                  T =   (boundarySolution_static(i,j,5,s1,e1) + boundarySolution(i,j,5,s1,e1))/&
                          (boundarySolution(i,j,4,s1,e1)+boundarySolution_static(i,j,4,s1,e1) )  
                                   
                  cIn  = sqrt( R_dev*T* &
                              ( (boundarySolution(i,j,6,s1,e1)+boundarySolution_static(i,j,6,s1,e1))/P0_dev )**rC_dev  )
                               
                  ! External normal velocity component
                  uOut = ( externalState(ii,jj,1,bID)*nHat(1,i,j,s1,e1)/norm + &
                           externalState(ii,jj,2,bID)*nHat(2,i,j,s1,e1)/norm + &
                           externalState(ii,jj,3,bID)*nHat(3,i,j,s1,e1)/norm )/& 
                         ( externalState(ii,jj,4,bID) + boundarySolution_static(i,j,4,s1,e1) )
                  ! Internal normal velocity component
                  uIn  = ( boundarySolution(i,j,1,s1,e1)*nHat(1,i,j,s1,e1)/norm + &
                           boundarySolution(i,j,2,s1,e1)*nHat(2,i,j,s1,e1)/norm + &
                           boundarySolution(i,j,3,s1,e1)*nHat(3,i,j,s1,e1)/norm )/& 
                         ( boundarySolution(i,j,4,s1,e1) + boundarySolution_static(i,j,4,s1,e1) )


                  fac = max( abs(uIn+cIn), abs(uIn-cIn), abs(uOut+cOut), abs(uOut-cOut) )

                  DO iEq = 1, nEq_dev-1
                        aS(iEq) = uIn*( boundarySolution(i,j,iEq,s1,e1) + boundarySolution_static(i,j,iEq,s1,e1) ) +&
                                 uOut*( externalState(ii,jj,iEq,bID) + boundarySolution_static(i,j,iEq,s1,e1) )
                  ENDDO
                  
                  ! Pressure !
                  DO k = 1, 3         
                  aS(k) = aS(k) + (boundarySolution(i,j,6,s1,e1)+externalState(ii,jj,6,bID))*nHat(k,i,j,s1,e1)/norm
                  ENDDO
                  
                          
                  DO iEq = 1, nEq_dev-1
                     boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*( aS(iEq) - fac*jump(iEq) )*norm
                     IF( iEq == 4 )THEN
                        DO k = 1, 3
                           jEq = k+(iEq-1)*3
                           ! Calculate the Bassi-Rebay flux for the stress tensor.
                           stressFlux(i,j,jEq,s1,e1) = 0.5_prec*( boundarySolution(i,j,iEq,s1,e1) +&
                                                                  externalState(ii,jj,iEq,bID)  )*& 
                                                                  nHat(k,i,j,s1,e1)
                        ENDDO
                     ELSE
                        DO k = 1, 3
                           jEq = k+(iEq-1)*3
                           ! Calculate the Bassi-Rebay flux for the stress tensor.
                           stressFlux(i,j,jEq,s1,e1) = 0.5_prec*( boundarySolution(i,j,iEq,s1,e1)/&
                                                                  (boundarySolution(i,j,4,s1,e1)+&
                                                                   boundarySolution_static(i,j,4,s1,e1)) +&
                                                                  externalState(ii,jj,iEq,bID)/&
                                                                  (externalState(ii,jj,4,bID)+&
                                                                   boundarySolution_static(i,j,4,s1,e1))  )*& 
                                                                  nHat(k,i,j,s1,e1)
                        ENDDO
                     ENDIF
                     
                  ENDDO
                  
               ENDIF 


 END SUBROUTINE BoundaryFaceFlux_CUDAKernel

!
 ATTRIBUTES(Global) SUBROUTINE MappedTimeDerivative_CUDAKernel( solution, static, boundaryFlux, drag, &
                                                                Ja, Jac, bMat, quadratureWeights, dMatP, tendency )

   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(in)  :: solution(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: static(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundaryFlux(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: drag(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: Ja(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:3,1:3,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: Jac(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: bMat(0:polydeg_dev,0:1)
   REAL(prec), DEVICE, INTENT(in)  :: quadratureWeights(0:polydeg_dev)
   REAL(prec), DEVICE, INTENT(in)  :: dMatP(0:polydeg_dev,0:polydeg_dev)
   REAL(prec), DEVICE, INTENT(out) :: tendency(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   ! Local
   INTEGER            :: i, j, k, row, col
   INTEGER, SHARED    :: iEl, iEq
   REAL(prec), SHARED :: contFlux(0:7,0:7,0:7,1:3)
   REAL(prec)         :: tend, F

      iEl = blockIDx % x
      iEq = blockIDx % y
      
      i = threadIdx % x - 1
      j = threadIdx % y - 1
      k = threadIdx % z - 1

      ! Here the flux tensor in physical space is calculated and rotated to give the 
      ! contravariant flux tensor in the reference computational domain.
      DO col = 1, 3
         contFlux(i,j,k,col) = 0.0_prec
         DO row = 1, 3
         !//////////////////////////////// Advection ///////////////////////////////////////!
               contFlux(i,j,k,col) = contFlux(i,j,k,col) +&
                                     Ja(i,j,k,row,col,iEl)*&
                                     solution(i,j,k,row,iEl)*&
                                      (solution(i,j,k,iEq,iEl) + static(i,j,k,iEq,iEl))/&    ! Density weighted variable being advected
                                  (solution(i,j,k,4,iEl) + static(i,j,k,4,iEl) )      
         ENDDO
     ! //////////////////// Pressure (Momentum only) /////////////////////////// !
         IF( iEq <= 3 )THEN
            contFlux(i,j,k,col) = contFlux(i,j,k,col) + Ja(i,j,k,iEq,col,iEl)*solution(i,j,k,6,iEl)
         ENDIF
         
      ENDDO                                  
            
      CALL syncthreads( )
      ! Now, the flux divergence is computed by multiplying the internally calculated fluxes by the
      ! DG-Derivative matrix and adding the boundary weighted fluxes.
      
      tend = 0.0_prec
      DO row = 0, polydeg_dev
         tend = tend + dMatP(row,i)*contFlux(row,j,k,1) + &
                       dMatP(row,j)*contFlux(i,row,k,2) + &
                       dMatP(row,k)*contFlux(i,j,row,3)
      ENDDO
       
      tend = -( tend + &
                ( boundaryFlux(i,k,iEq,1,iEl)*bmat(j,0) + &
                  boundaryFlux(i,k,iEq,3,iEl)*bMat(j,1) )/&
                quadratureWeights(j) + &
                ( boundaryFlux(j,k,iEq,4,iEl)*bMat(i,0) + &
                  boundaryFlux(j,k,iEq,2,iEl)*bMat(i,1) )/&
                quadratureWeights(i) + &
                ( boundaryFlux(i,j,iEq,5,iEl)*bMat(k,0) + &
                  boundaryFlux(i,j,iEq,6,iEl)*bMat(k,1) )/&
                quadratureWeights(k) )/Jac(i,j,k,iEl)
                      
             
      tendency(i,j,k,iEq,iEl) = tend
      F = sqrt( solution(i,j,k,1,iEl)**2 + &
                solution(i,j,k,2,iEl)**2 + &
                solution(i,j,k,3,iEl)**2 ) /&
                  (solution(i,j,k,4,iEl) + static(i,j,k,4,iEl))
                  
      IF( iEq == 1 )THEN
         tendency(i,j,k,1,iEl) = tendency(i,j,k,1,iEl) -&
                                 drag(i,j,k,iEl)*solution(i,j,k,1,iEl)*F-&
                                 solution(i,j,k,3,iEl)*fRotY_dev +&
                                 solution(i,j,k,2,iEl)*fRotz_dev 
      
      ELSEIF( iEq == 2 )THEN
         tendency(i,j,k,2,iEl) = tendency(i,j,k,2,iEl) -&
                                 drag(i,j,k,iEl)*solution(i,j,k,2,iEl)*F -&
                                 solution(i,j,k,1,iEl)*fRotZ_dev +&
                                 solution(i,j,k,3,iEl)*fRotX_dev 
      ELSEIF( iEq == 3 )THEN ! Add in the buoyancy acceleration
         tendency(i,j,k,3,iEl) = tendency(i,j,k,3,iEl) -&
                                 drag(i,j,k,iEl)*solution(i,j,k,3,iEl)*F -&
                                 solution(i,j,k,2,iEl)*fRotX_dev +&
                                 solution(i,j,k,1,iEl)*fRotY_dev -&
                                 solution(i,j,k,4,iEl)*g_dev
      ENDIF
         

      
       
 END SUBROUTINE MappedTimeDerivative_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE CalculateStressTensor_CUDAKernel( solution, static, dMatP, bmat, quadratureWeights, Ja, Jac, stressFlux, stressTensor ) 
 
   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(in)  :: solution(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: static(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: dMatP(0:polydeg_dev,0:polydeg_dev)
   REAL(prec), DEVICE, INTENT(in)  :: bmat(0:polydeg_dev,0:1)
   REAL(prec), DEVICE, INTENT(in)  :: quadratureWeights(0:polydeg_dev)
   REAL(prec), DEVICE, INTENT(in)  :: Ja(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:3,1:3,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: Jac(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: stressFlux(0:polydeg_dev,0:polydeg_dev,1:15,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(out) :: stressTensor(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:15,1:nEl_dev)
   ! Local
   INTEGER :: iEl, iEq, idir, i, j, k, m
   REAL(prec), SHARED :: contFlux(0:7,0:7,0:7,1:3)
   REAL(prec) :: strTens


      iEl = blockIDx % x
      iEq = blockIDx % y
      idir = blockIDx % z
      
      i = threadIDx % x-1
      j = threadIDx % y-1
      k = threadIDx % z-1
      
         
         ! Here the flux tensor in physical space is calculated and rotated to give the 
         ! contravariant flux tensor in the reference computational domain.
         IF( iEq == 4 )THEN
			 DO m = 1, 3
				contFlux(i,j,k,m) = Ja(i,j,k,idir,m,iEl)*solution(i,j,k,iEq,iEl)
			 ENDDO
	     ELSE
	         DO m = 1, 3
				contFlux(i,j,k,m) = Ja(i,j,k,idir,m,iEl)*solution(i,j,k,iEq,iEl)/&
				                      ( solution(i,j,k,4,iEl)+static(i,j,k,4,iEl) )
			 ENDDO
	     ENDIF
	     CALL syncthreads( ) 
	            
      ! Now, the flux divergence is computed by multiplying the internally calculated fluxes by the
      ! DG-Derivative matrix and adding the boundary weighted fluxes.
      
         ! Reduction for the stress tensor
         strTens = 0.0_prec
         DO m = 0, polydeg_dev
            strTens = strTens + dMatP(m,i)*contFlux(m,j,k,1) + &
                                dMatP(m,j)*contFlux(i,m,k,2) + &
                                dMatP(m,k)*contFlux(i,j,m,3)
         ENDDO


         stressTensor(i,j,k,idir + (iEq-1)*3,iEl) = ( strTens + &
                ( stressFlux(i,k,idir + (iEq-1)*3,1,iEl)*bmat(j,0) + &
                  stressFlux(i,k,idir + (iEq-1)*3,3,iEl)*bMat(j,1) )/&
                quadratureWeights(j) + &
                ( stressFlux(j,k,idir + (iEq-1)*3,4,iEl)*bMat(i,0) + &
                  stressFlux(j,k,idir + (iEq-1)*3,2,iEl)*bMat(i,1) )/&
                quadratureWeights(i) + &
                ( stressFlux(i,j,idir + (iEq-1)*3,5,iEl)*bMat(k,0) + &
                  stressFlux(i,j,idir + (iEq-1)*3,6,iEl)*bMat(k,1) )/&
                quadratureWeights(k) )/Jac(i,j,k,iEl)

                     
 END SUBROUTINE CalculateStressTensor_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE CalculateBoundaryStress_CUDAKernel( solution, bMat, boundarySolution ) 

   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(in)  :: solution(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:15,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: bMat(0:polydeg_dev,0:1)
   REAL(prec), DEVICE, INTENT(out) :: boundarySolution(0:polydeg_dev,0:polydeg_dev,1:15,1:6,1:nEl_dev)
   ! Local
   INTEGER :: iEq, iEl, i, j, k
   REAL(prec) :: bSol(1:6)

      iEl = blockIdx % x
      
      iEq = threadIdx % z
      k   = threadIdx % y-1
      j   = threadIdx % x-1
      
	  bSol(1:6) = 0.0_prec

	  DO i = 0, polydeg_dev
		 bSol(1) = bSol(1) + bMat(i,0)*solution(j,i,k,iEq,iEl) ! south
			
		 bSol(2) = bSol(2) + bMat(i,1)*solution(i,j,k,iEq,iEl) ! east
			
		 bSol(3) = bSol(3) + bMat(i,1)*solution(j,i,k,iEq,iEl) ! north
															
		 bSol(4) = bSol(4) + bMat(i,0)*solution(i,j,k,iEq,iEl) ! west

		 bSol(5) = bSol(5) + bMat(i,0)*solution(j,k,i,iEq,iEl) ! botom

		 bSol(6) = bSol(6) + bMat(i,1)*solution(j,k,i,iEq,iEl) ! top
	  ENDDO
               
      DO i = 1, 6
         boundarySolution(j,k,iEq,i,iEl) = bSol(i)
      ENDDO
      
 END SUBROUTINE CalculateBoundaryStress_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE UpdateExternalStress_CUDAKernel( boundaryIDs, elementIDs, &
                                                               elementSides, procIDs, &
                                                               externalStress, &
                                                               stressBsols, &
                                                               prescribedStress, nHat ) 
   IMPLICIT NONE
   INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nBoundaryFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: procIDs(1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(out) :: externalStress(0:polydeg_dev,0:polydeg_dev,1:15,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: stressBsols(0:polydeg_dev,0:polydeg_dev,1:15,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: prescribedStress(0:polydeg_dev,0:polydeg_dev,1:15,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: nhat(1:3,0:polydeg_dev,0:polydeg_dev,1:6,1:nEl_dev)
   ! Local
   INTEGER    :: iEq, iFace, i, j, k
   INTEGER    :: iFace2, p2
   INTEGER    :: e1, e2, s1, s2, m
   
      iFace = blockIdx % x
      iEq   = blockIDx % y
      ! ////////////////////////////////////////////////////////////////////////// !
      i   = threadIdx % x-1
      j   = threadIdx % y-1

	  iFace2 = boundaryIDs( iFace ) ! Obtain the process-local face id for this boundary-face id
	  e1     = elementIDs(1,iFace2)
	  s1     = elementSides(1,iFace2)
	  e2     = elementIDs(2,iFace2)
	  p2     = procIDs( iFace )
	 
	  IF( i <= polydeg_dev .AND. j <= polydeg_dev )THEN
	 
		 IF( p2 == myRank_dev )THEN
		    externalStress(i,j,iEq,iFace) = -stressBsols(i,j,iEq,s1,e1)
		 ENDIF
		
	 ENDIF
                  
 END SUBROUTINE UpdateExternalStress_CUDAKernel 
!
 ATTRIBUTES(Global) SUBROUTINE InternalStressFlux_CUDAKernel( elementIDs, elementSides, boundaryIDs, iMap, jMap, &
                                                      nHat, boundaryState, static, boundaryStress, sgsCoeffs, &
                                                      externalState, externalStress, boundaryFlux )

   IMPLICIT NONE
   INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: iMap(0:polydeg_dev,0:polydeg_dev,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: jMap(0:polydeg_dev,0:polydeg_dev,1:nFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:polydeg_dev,0:polydeg_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundaryState(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundaryStress(0:polydeg_dev,0:polydeg_dev,1:15,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: sgsCoeffs(0:polydeg_dev,0:polydeg_dev,1:5,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: externalState(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: static(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: externalStress(0:polydeg_dev,0:polydeg_dev,1:15,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(out) :: boundaryFlux(0:polydeg_dev,0:polydeg_dev,1:15,1:6,1:nEl_dev)
   ! Local
   INTEGER    :: iEl, iFace
   INTEGER    :: i, j, iEq
   INTEGER    :: ii, jj, bID
   INTEGER    :: e1, s1, e2, s2
   INTEGER    :: m, jEq
   REAL(prec) :: norm, rhoOut, rhoIn

      iFace = blockIdx % x
      iEq   = blockIdx % y
      j     = threadIdx % y-1
      i     = threadIdx % x-1
     
      e1 = elementIDs(1,iFace)
      s1 = elementSides(1,iFace)
      e2 = elementIDs(2,iFace)
      s2 = ABS(elementSides(2,iFace))
      bID  = ABS(boundaryIDs(iFace))

      ii = iMap(i,j,iFace)
      jj = jMap(i,j,iFace)
               
      norm = sqrt( nHat(1,i,j,s1,e1)**2 + nHat(2,i,j,s1,e1)**2 + nHat(3,i,j,s1,e1)**2 )
      
      IF( e2 > 0 )THEN
      
         boundaryFlux(i,j,iEq,s1,e1) =  0.5_prec*( sgsCoeffs(ii,jj,iEq,s2,e2)*boundaryState(ii,jj,iEq,s2,e2)-&
						      sgsCoeffs(i,j,iEq,s1,e1)*boundaryState(i,j,iEq,s1,e1))/viscLengthScale_dev*norm
         
		   IF( iEq == 4 )THEN
       
            
			   DO m = 1, 3    
				   jEq = m + (iEq-1)*3  
				   boundaryFlux(i,j,iEq,s1,e1) = boundaryFlux(i,j,iEq,s1,e1) + &
				     0.5_prec*(sgsCoeffs(i,j,iEq,s1,e1)*boundaryStress(i,j,jEq,s1,e1)+&
					    	     sgsCoeffs(ii,jj,iEq,s2,e2)*boundaryStress(ii,jj,jEq,s2,e2) )*nHat(m,i,j,s1,e1)
			 ENDDO
          
		   ELSE
                        
          rhoIn = static(i,j,4,s1,e1) + boundaryState(i,j,4,s1,e1)
          rhoOut = static(ii,jj,4,s2,e2) + boundaryState(ii,jj,4,s2,e2)
          
            DO m = 1, 3    
				  jEq = m + (iEq-1)*3  
				  boundaryFlux(i,j,iEq,s1,e1) = boundaryFlux(i,j,iEq,s1,e1) + &
				    0.5_prec*(rhoIn*sgsCoeffs(i,j,iEq,s1,e1)*boundaryStress(i,j,jEq,s1,e1)+&
                          rhoOut*sgsCoeffs(ii,jj,iEq,s2,e2)*boundaryStress(ii,jj,jEq,s2,e2))*nHat(m,i,j,s1,e1)
			   ENDDO
          
		   ENDIF
			 
         boundaryFlux(ii,jj,iEq,s2,e2) = -boundaryFlux(i,j,iEq,s1,e1)
             
      ENDIF

 END SUBROUTINE InternalStressFlux_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE BoundaryStressFlux_CUDAKernel( elementIDs, elementSides, boundaryIDs, iMap, jMap, &
                                                      nHat, boundaryState, static, boundaryStress, sgsCoeffs, externalSGS, &
                                                      externalState, externalStress, boundaryFlux )

   IMPLICIT NONE
   INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: iMap(0:polydeg_dev,0:polydeg_dev,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: jMap(0:polydeg_dev,0:polydeg_dev,1:nFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:polydeg_dev,0:polydeg_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundaryState(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundaryStress(0:polydeg_dev,0:polydeg_dev,1:15,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: sgsCoeffs(0:polydeg_dev,0:polydeg_dev,1:5,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: externalSGS(0:polydeg_dev,0:polydeg_dev,1:5,nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: externalState(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: static(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: externalStress(0:polydeg_dev,0:polydeg_dev,1:15,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(out) :: boundaryFlux(0:polydeg_dev,0:polydeg_dev,1:15,1:6,1:nEl_dev)
   ! Local
   INTEGER    :: iEl, iFace
   INTEGER    :: i, j, iEq
   INTEGER    :: ii, jj, bID
   INTEGER    :: e1, s1, e2, s2
   INTEGER    :: m, jEq
   REAL(prec) :: norm, rhoOut, rhoIn

      iFace = blockIdx % x
      iEq   = blockIdx % y
      j     = threadIdx % y-1
      i     = threadIdx % x-1
     
      e1 = elementIDs(1,iFace)
      s1 = elementSides(1,iFace)
      e2 = elementIDs(2,iFace)
      s2 = ABS(elementSides(2,iFace))
      bID  = boundaryIDs(iFace)

      ii = iMap(i,j,iFace)
      jj = jMap(i,j,iFace)
               
      norm = sqrt( nHat(1,i,j,s1,e1)**2 + nHat(2,i,j,s1,e1)**2 + nHat(3,i,j,s1,e1)**2 )
      IF( e2 < 0 )THEN !Physical boundary
         IF( bID < 0 )THEN
         
            bID = ABS(bID)
            boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*sgsCoeffs(i,j,iEq,s1,e1)*&
                                        (externalState(ii,jj,iEq,bID)-boundaryState(i,j,iEq,s1,e1))/viscLengthScale_dev*norm

           IF( iEq == 4 )THEN
             
             DO m = 1, 3    
               jEq = m + (iEq-1)*3  
               boundaryFlux(i,j,iEq,s1,e1) = boundaryFlux(i,j,iEq,s1,e1) + &
                                             0.5_prec*sgsCoeffs(i,j,iEq,s1,e1)*(boundaryStress(i,j,jEq,s1,e1)+&
                                             externalStress(ii,jj,jEq,bID) )*nHat(m,i,j,s1,e1)
             ENDDO
             
           ELSE
                           
             rhoIn  = static(i,j,4,s1,e1) + boundaryState(i,j,4,s1,e1)
             rhoOut = static(i,j,4,s1,e1) + externalState(ii,jj,4,bID)
             
             DO m = 1, 3    
               jEq = m + (iEq-1)*3  
               boundaryFlux(i,j,iEq,s1,e1) = boundaryFlux(i,j,iEq,s1,e1) + &
                                             0.5_prec*sgsCoeffs(i,j,iEq,s1,e1)*( rhoIn*boundaryStress(i,j,jEq,s1,e1)+&
                                                                           rhoOut*externalStress(ii,jj,jEq,bID) )*nHat(m,i,j,s1,e1)
             ENDDO
             
           ENDIF

         ELSE
         
            boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*( externalSGS(ii,jj,iEq,bID)*externalState(ii,jj,iEq,bID)-&
                                                     sgsCoeffs(i,j,iEq,s1,e1)*boundaryState(i,j,iEq,s1,e1))/viscLengthScale_dev*norm

           IF( iEq == 4 )THEN
             
             DO m = 1, 3    
               jEq = m + (iEq-1)*3  
               boundaryFlux(i,j,iEq,s1,e1) = boundaryFlux(i,j,iEq,s1,e1) + &
                                             0.5_prec*( sgsCoeffs(i,j,iEq,s1,e1)*boundaryStress(i,j,jEq,s1,e1)+&
                                                        externalSGS(ii,jj,iEq,bID)*externalStress(ii,jj,jEq,bID) )*nHat(m,i,j,s1,e1)
             ENDDO
             
           ELSE
                           
             rhoIn  = static(i,j,4,s1,e1) + boundaryState(i,j,4,s1,e1)
             rhoOut = static(i,j,4,s1,e1) + externalState(ii,jj,4,bID)
             
             DO m = 1, 3    
               jEq = m + (iEq-1)*3  
               boundaryFlux(i,j,iEq,s1,e1) = boundaryFlux(i,j,iEq,s1,e1) + &
                                             0.5_prec*( sgsCoeffs(i,j,iEq,s1,e1)*rhoIn*boundaryStress(i,j,jEq,s1,e1)+&
                                                        externalSGS(ii,jj,iEq,bID)*rhoOut*externalStress(ii,jj,jEq,bID) )*nHat(m,i,j,s1,e1)
             ENDDO
             
           ENDIF
         ENDIF

      ENDIF

 END SUBROUTINE BoundaryStressFlux_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE StressDivergence_CUDAKernel( stress, stressFlux, state, static, sgsCoeffs, &
                                                             Ja, Jac, bMat, quadratureWeights, dMatP, tendency ) ! ///////////////////// !

   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(in)  :: stress(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:15,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: stressFlux(0:polydeg_dev,0:polydeg_dev,1:15,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: state(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: static(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: sgsCoeffs(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev-1,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: Ja(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:3,1:3,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: Jac(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: bMat(0:polydeg_dev,0:1)
   REAL(prec), DEVICE, INTENT(in)  :: quadratureWeights(0:polydeg_dev)
   REAL(prec), DEVICE, INTENT(in)  :: dMatP(0:polydeg_dev,0:polydeg_dev)
   REAL(prec), DEVICE, INTENT(out) :: tendency(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:15,1:nEl_dev)
   ! Local
   INTEGER            :: i, j, k, row, col
   INTEGER, SHARED    :: iEl, iEq, jEq
   REAL(prec), SHARED :: contFlux(0:7,0:7,0:7,1:3)
   REAL(prec)         :: tend

      iEl = blockIDx % x
      iEq = blockIDx % y
      
      i = threadIdx % x - 1
      j = threadIdx % y - 1
      k = threadIdx % z - 1

      ! Here the flux tensor in physical space is calculated and rotated to give the 
      ! contravariant flux tensor in the reference computational domain.
      IF( iEq == 4 )THEN
         DO col = 1, 3
            contFlux(i,j,k,col) = 0.0_prec
            DO row = 1, 3
               jEq = row + (iEq-1)*3
               contFlux(i,j,k,col) = contFlux(i,j,k,col) +&
                                    Ja(i,j,k,row,col,iEl)*&
                                    stress(i,j,k,jEq,iEl)*&
                                    sgsCoeffs(i,j,k,iEq,iEl)
                        
            ENDDO
         ENDDO                                  
      ELSE
        DO col = 1, 3
            contFlux(i,j,k,col) = 0.0_prec
            DO row = 1, 3
               jEq = row + (iEq-1)*3
               contFlux(i,j,k,col) = contFlux(i,j,k,col) +&
                                    Ja(i,j,k,row,col,iEl)*&
                                    stress(i,j,k,jEq,iEl)*&
                                    (state(i,j,k,4,iEl)+static(i,j,k,4,iEl))*&
                                    sgsCoeffs(i,j,k,iEq,iEl)
                        
           ENDDO
        ENDDO
      
      ENDIF    
      CALL syncthreads( )
      ! Now, the flux divergence is computed by multiplying the internally calculated fluxes by the
      ! DG-Derivative matrix and adding the boundary weighted fluxes.
      
      tend = 0.0_prec
      DO row = 0, polydeg_dev
         tend = tend + dMatP(row,i)*contFlux(row,j,k,1) + &
                       dMatP(row,j)*contFlux(i,row,k,2) + &
                       dMatP(row,k)*contFlux(i,j,row,3)
      ENDDO
       
      tend = ( tend + &
                ( stressFlux(i,k,iEq,1,iEl)*bmat(j,0) + &
                  stressFlux(i,k,iEq,3,iEl)*bMat(j,1) )/&
                quadratureWeights(j) + &
                ( stressFlux(j,k,iEq,4,iEl)*bMat(i,0) + &
                  stressFlux(j,k,iEq,2,iEl)*bMat(i,1) )/&
                quadratureWeights(i) + &
                ( stressFlux(i,j,iEq,5,iEl)*bMat(k,0) + &
                  stressFlux(i,j,iEq,6,iEl)*bMat(k,1) )/&
                quadratureWeights(k) )/Jac(i,j,k,iEl)
                      
             
      tendency(i,j,k,iEq,iEl) = tend
       
 END SUBROUTINE StressDivergence_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE EquationOfState_CUDAKernel( solution, static )
   ! This routine calculates the anomalous pressure referenced to the static state.
   ! The pressure is calculated using the ideal gas law.
   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(inout) :: solution(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)    :: static(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   ! Local
   INTEGER :: i, j, k, iEl
   REAL(prec) :: rhoT

      iEl = blockIdx % x
      i   = threadIdx % x - 1
      j   = threadIdx % y - 1
      k   = threadIdx % z - 1

      ! Pressure = rho*e*R/Cv (Ideal Gas Law, P = rho*R*T, thermo ~> T = theta*(P/P0)^r)
      ! Then P = P0*(rho*theta*R/P0)^(Cp/Cv)
      ! And P' = P - P_static
      rhoT = static(i,j,k,5,iEl) + solution(i,j,k,5,iEl)
      solution(i,j,k,6,iEl) = P0_dev*( rhoT*R_dev/P0_dev )**hCapRatio_dev - static(i,j,k,6,iEl)

 END SUBROUTINE EquationOfState_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE BoundaryToBuffer_CUDAKernel( sendBuffer, boundarySolution, faceToElement, faceToSide, boundaryToFaceID, &
                                                            boundaryToProcID, rankToNeighbor, boundaryToBuffer, nFaces, nBoundaries, nRanks, myRank, N, numEq, nNeighbors, &
                                                            bufferSize, nElements )
   IMPLICIT NONE
   INTEGER, VALUE, INTENT(in)        :: nFaces, nBoundaries, nRanks, myRank
   INTEGER, VALUE, INTENT(in)        :: N, numEq, nNeighbors, nElements, bufferSize
   INTEGER, DEVICE, INTENT(in)       :: faceToElement(1:2,1:nFaces)
   INTEGER, DEVICE, INTENT(in)       :: faceToSide(1:2,1:nFaces)
   INTEGER, DEVICE, INTENT(in)       :: boundaryToFaceID(1:nBoundaries)
   INTEGER, DEVICE, INTENT(in)       :: boundaryToProcID(1:nBoundaries)
   INTEGER, DEVICE, INTENT(in)       :: rankToNeighbor(0:nRanks-1)
   INTEGER, DEVICE, INTENT(in)       :: boundaryToBuffer(1:nBoundaries)
   REAL(prec), DEVICE, INTENT(inout) :: sendBuffer(0:N,0:N,1:numEq,1:bufferSize,1:nNeighbors)
   REAL(prec), DEVICE, INTENT(in)    :: boundarySolution(0:N,0:N,1:numEq,1:6,1:nElements)
   ! Local
   INTEGER :: bID, i, j, iEq
   INTEGER :: iFace, p2, neighborID, bufferID, elementID, sideID

      bID = blockIdx % x
      i   = threadIdx % x - 1
      j   = threadIdx % y - 1
      iEq = threadIdx % z 
      
      p2 = boundaryToProcID(bID)
      
      IF( p2 /= myRank )THEN 

         neighborID = rankToNeighbor(p2)
         bufferID   = boundaryToBuffer(bID)
         iFace      = boundaryToFaceID(bID)
         elementID  = faceToElement(1,iFace)
         sideID     = faceToSide(1,iFace)

         !PRINT*, "ACCESS", i, j, iEq, bufferID, neighborID
         sendBuffer(i,j,iEq,bufferID,neighborID) = boundarySolution(i,j,iEq,sideID,elementID) 
      ENDIF

 END SUBROUTINE BoundaryToBuffer_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE BufferToBoundary_CUDAKernel( recvBuffer, externalSolution, boundaryToFaceID, &
                                                            boundaryToProcID, rankTable, unPackMap, nFaces, nBoundaries, nRanks, myRank, N, numEq, nNeighbors, &
                                                            bufferSize, nElements )
   IMPLICIT NONE
   INTEGER, VALUE, INTENT(in)        :: nFaces, nBoundaries, nRanks, myRank
   INTEGER, VALUE, INTENT(in)        :: N, numEq, nNeighbors, nElements, bufferSize
   INTEGER, DEVICE, INTENT(in)       :: boundaryToFaceID(1:nBoundaries)
   INTEGER, DEVICE, INTENT(in)       :: boundaryToProcID(1:nBoundaries)
   INTEGER, DEVICE, INTENT(in)       :: rankTable(0:nRanks-1)
   INTEGER, DEVICE, INTENT(in)       :: unPackMap(1:nBoundaries)
   REAL(prec), DEVICE, INTENT(in)    :: recvBuffer(0:N,0:N,1:numEq,1:bufferSize,1:nNeighbors)
   REAL(prec), DEVICE, INTENT(inout) :: externalSolution(0:N,0:N,1:numEq,1:nBoundaries)
   ! Local
   INTEGER :: bID, i, j, iEq
   INTEGER :: p2

      bID = blockIdx % x
      i   = threadIdx % x - 1
      j   = threadIdx % y - 1
      iEq = threadIdx % z 
      
      p2        = boundaryToProcID(bID)
      

      IF( p2 /= myRank )THEN 
         externalSolution(i,j,iEq,bID) = recvBuffer(i,j,iEq,unpackMap(bID),rankTable(p2))
      ENDIF


 END SUBROUTINE BufferToBoundary_CUDAKernel

#endif
 END MODULE Fluid_Class


