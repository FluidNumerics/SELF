! Fluid_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE Fluid_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
#ifdef TESTING
USE ModelDataInstances_Class
#endif
! src/highend/Fluid/
USE FluidParams_Class

#ifdef HAVE_CUDA
! CUDA libraries
! src/spectralops/cuda/
USE Lagrange_Cuda_Class
USE NodalStorage_Cuda_Class
! src/filters/cuda/
USE RollOffFilter_Cuda_Class
! src/solutionstorage/cuda/
USE DGSEMSolutionStorage_3D_Cuda_Class
! src/geometry/cuda/
USE Face_Cuda_Class
USE Element_Cuda_Class  
USE HexMesh_Cuda_Class 
USE BoundaryCommunicator_Cuda_Class

USE cudafor
#else

! src/spectralops/
USE Lagrange_Class
USE NodalStorage_Class
! src/filters/
USE RollOffFilter_Class
! src/solutionstorage/
USE DGSEMSolutionStorage_3D_Class
! src/geometry/
USE HexMesh_Class
USE BoundaryCommunicator_Class

#endif
IMPLICIT NONE

#ifdef HAVE_MPI
INCLUDE 'mpif.h'

    TYPE PairWiseMPIPacket
       ! For the unstructured mesh, I opt for building my own data structure that bundles messages between neighboring
       ! ranks, rather than attempting to build an MPI data structure. If a structured mesh is used, one optimization 
       ! would be to use MPI-data types to handle the message passing.
       !
       ! :: Attributes ::  
       ! unPackMap(1:nSharedFaces)   >> Maps the local shared boundary face ID in the message to the correct boundaryID
        
       INTEGER :: neighborRank
       INTEGER :: bufferSize
       INTEGER :: nSharedFaces
       REAL(prec), ALLOCATABLE :: sendStateBuffer(:,:,:,:) ! (0:N,0:N,1:nEq,1:nSharedFaces)
       REAL(prec), ALLOCATABLE :: recvStateBuffer(:,:,:,:)
       REAL(prec), ALLOCATABLE :: sendStressBuffer(:,:,:,:)
       REAL(prec), ALLOCATABLE :: recvStressBuffer(:,:,:,:)
       REAL(prec), ALLOCATABLE :: sendSGSBuffer(:,:,:,:)
       REAL(prec), ALLOCATABLE :: recvSGSBuffer(:,:,:,:)
       INTEGER                 :: bufferCounter
    END TYPE PairWiseMPIPacket
#endif

    TYPE Fluid
      INTEGER                                :: nEq, N, nBoundaryFaces, nNeighbors
      REAL(prec)                             :: simulationTime
      TYPE( FluidParams )                    :: params
      REAL(prec), ALLOCATABLE                :: dragProfile(:,:,:,:)
#ifdef HAVE_MPI
      TYPE( PairWiseMPIPacket ), ALLOCATABLE :: mpiPackets(:)
#endif

#ifdef HAVE_CUDA
      REAL(prec), DEVICE, ALLOCATABLE        :: dragProfile_dev(:,:,:,:)
      TYPE( HexMesh_Cuda )                   :: mesh
      TYPE( NodalStorage_Cuda )              :: dGStorage
      TYPE( RollOffFilter_Cuda )             :: filter
      TYPE( DGSEMSolution_3D_Cuda )          :: state
      TYPE( LightDGSEMSolution_3D_Cuda )     :: smoothState
      TYPE( DGSEMSolution_3D_Cuda )          :: static
      TYPE( DGSEMSolution_3D_Cuda )          :: stressTensor
      TYPE( LightDGSEMSolution_3D_Cuda )     :: sgsCoeffs
#else
      TYPE( HexMesh )                        :: mesh
      TYPE( NodalStorage )                   :: dGStorage
      TYPE( RollOffFilter )                  :: filter
      TYPE( DGSEMSolution_3D )               :: state
      TYPE( LightDGSEMSolution_3D )          :: smoothState
      TYPE( DGSEMSolution_3D )               :: static
      TYPE( DGSEMSolution_3D )               :: stressTensor
      TYPE( LightDGSEMSolution_3D )          :: sgsCoeffs
#endif
      ! ////////////  Boundary communication information  //////////// !
#ifdef HAVE_CUDA
      TYPE( BoundaryCommunicator_Cuda )       :: extComm
#else
      TYPE( BoundaryCommunicator )            :: extComm
#endif

      INTEGER, ALLOCATABLE                    :: rankTable(:)
      
      REAL(prec), ALLOCATABLE                 :: prescribedState(:,:,:,:)
      REAL(prec), ALLOCATABLE                 :: externalState(:,:,:,:)
      REAL(prec), ALLOCATABLE                 :: prescribedStress(:,:,:,:)
      REAL(prec), ALLOCATABLE                 :: externalStress(:,:,:,:)
      REAL(prec), ALLOCATABLE                 :: externalSGS(:,:,:,:)
#ifdef HAVE_CUDA
      REAL(prec), DEVICE, ALLOCATABLE         :: prescribedState_dev(:,:,:,:)
      REAL(prec), DEVICE, ALLOCATABLE         :: externalState_dev(:,:,:,:)
      REAL(prec), DEVICE, ALLOCATABLE         :: prescribedStress_dev(:,:,:,:)
      REAL(prec), DEVICE, ALLOCATABLE         :: externalStress_dev(:,:,:,:)
      REAL(prec), DEVICE, ALLOCATABLE         :: externalSGS_dev(:,:,:,:)
#endif
      ! ////////////////////////////////////////////////////////////// !

      CONTAINS

      PROCEDURE :: Build => Build_Fluid
      PROCEDURE :: Trash => Trash_Fluid
      PROCEDURE, PRIVATE :: BuildHexMesh => BuildHexMesh_Fluid
      
      PROCEDURE :: CalculateStaticState => CalculateStaticState_Fluid
      ! Time integrators
      PROCEDURE :: ForwardStepRK3        => ForwardStepRK3_Fluid

      ! ////////////////////////// Equation Of State ////////////////////////////////////////////////////// !
      PROCEDURE :: EquationOfState        => EquationOfState_Fluid
      
      ! /////////////////////////////////////////////////////////////////////////////////////////////////// !
      PROCEDURE :: GlobalTimeDerivative         => GlobalTimeDerivative_Fluid
      
      
      !  /////////// Routines for interpolating state and static attributes to the element faces ////////// !
      PROCEDURE :: CalculateStaticBoundarySolution => CalculateStaticBoundarySolution_Fluid
      PROCEDURE :: CalculateBoundarySolution       => CalculateBoundarySolution_Fluid
      PROCEDURE :: UpdateExternalState             => UpdateExternalState_Fluid
      PROCEDURE :: InternalFaceFlux                => InternalFaceFlux_Fluid
      PROCEDURE :: BoundaryFaceFlux                => BoundaryFaceFlux_Fluid
      PROCEDURE :: MappedTimeDerivative            => MappedTimeDerivative_Fluid
      
#ifdef HAVE_MPI
      PROCEDURE :: ConstructCommTables             => ConstructCommTables_Fluid
      PROCEDURE :: MPI_StateExchange               => MPI_StateExchange_Fluid
      PROCEDURE :: FinalizeMPI_StateExchange       => FinalizeMPI_StateExchange_Fluid
      PROCEDURE :: MPI_StressExchange              => MPI_StressExchange_Fluid
      PROCEDURE :: FinalizeMPI_StressExchange      => MPI_StressExchange_Fluid
      PROCEDURE :: MPI_SGSExchange                 => MPI_SGSExchange_Fluid
      PROCEDURE :: FinalizeMPI_SGSExchange         => MPI_SGSExchange_Fluid
#endif
      ! /////////////////////////////////////////////////////////////////////////////////////////////////// !
      ! Routines for the fluid-stress model and
      ! Routines for spectral filtering and calculating the sgs coefficients
      PROCEDURE :: CalculateSmoothedState   => CalculateSmoothedState_Fluid
      PROCEDURE :: CalculateSGSCoefficients => CalculateSGSCoefficients_Fluid
      PROCEDURE :: CalculateBoundarySGS     => CalculateBoundarySGS_Fluid
      PROCEDURE :: UpdateExternalSGS        => UpdateExternalSGS_Fluid
      
      PROCEDURE :: CalculateStressTensor   => CalculateStressTensor_Fluid
      PROCEDURE :: CalculateBoundaryStress => CalculateBoundaryStress_Fluid
      PROCEDURE :: UpdateExternalStress    => UpdateExternalStress_Fluid
      PROCEDURE :: InternalStressFlux      => InternalStressFlux_Fluid
      PROCEDURE :: BoundaryStressFlux      => BoundaryStressFlux_Fluid
      PROCEDURE :: StressDivergence        => StressDivergence_Fluid
      
      !
      
      ! /////////////////////////////////////////////////////////////////////////////////////////////////// !
      PROCEDURE :: WriteTecplot => WriteTecplot_Fluid
      PROCEDURE :: WriteSmoothedTecplot => WriteSmoothedTecplot_Fluid
      PROCEDURE :: WriteSGSTecplot      => WriteSGSTecplot_Fluid
      PROCEDURE :: WriteStressTensorTecplot => WriteStressTensorTecplot_Fluid
      PROCEDURE :: WritePickup => WritePickup_Fluid
      PROCEDURE :: ReadPickup => ReadPickup_Fluid
  !    PROCEDURE :: QuickDiagnostics => QuickDiagnostics_Fluid
      PROCEDURE :: FluidStateAtPlottingPoints => FluidStateAtPlottingPoints_Fluid
      PROCEDURE :: ObtainPlottingMesh         => ObtainPlottingMesh_Fluid

    END TYPE Fluid

 INTEGER, PARAMETER   :: nEq = 6

#ifdef HAVE_CUDA
 INTEGER, CONSTANT    :: nEq_dev
 INTEGER, CONSTANT    :: polyDeg_dev
 INTEGER, CONSTANT    :: nEl_dev
 INTEGER, CONSTANT    :: myRank_dev
 REAL(prec), CONSTANT :: R_dev, Cv_dev, P0_dev, hCapRatio_dev, rC_dev, g_dev
 REAL(prec), CONSTANT :: rk3_a_dev, rk3_g_dev, dt_dev
 REAL(prec), CONSTANT :: viscLengthScale_dev, dScale_dev, Cd_dev
 REAL(prec), CONSTANT :: fRotX_dev, fRotY_dev, fRotZ_dev
#endif
 
#ifdef HAVE_MPI
 INTEGER :: MPI_PREC
 INTEGER, ALLOCATABLE :: stateReqHandle(:), stressReqHandle(:), SGSReqHandle(:)
 INTEGER, ALLOCATABLE :: stateStats(:,:), stressStats(:,:), SGSStats(:,:)
#endif

 INTEGER :: callID

#ifdef TESTING
 TYPE( ModelDataInstances ) :: mdi
#endif

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
  SUBROUTINE Build_Fluid( myDGSEM, myRank, nProc )

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)         :: myRank, nProc
   !
   integer :: istat
#ifdef HAVE_CUDA
   integer(kind=cuda_count_kind) :: freebytes, totalbytes
#endif   

#ifdef HAVE_MPI
      IF( prec == sp )THEN
         MPI_PREC = MPI_FLOAT
      ELSE
         MPI_PREC = MPI_DOUBLE
      ENDIF
#endif
      callid = 0
      CALL myDGSEM % params % Build( )
      myDGSEM % N   = myDGSEM % params % polyDeg
      myDGSEM % nEq = nEq
      myDGSEM % simulationTime = myDGSEM % params % iterInit*myDGSEM % params % dt

#ifdef TESTING
      PRINT*, '  Module Fluid_Class.f90 : S/R Build_Fluid :'
      PRINT*, '    Testing is enabled. Restricting number of time steps to 1.'
      myDGSEM % params % nTimeSteps = 1
      myDGSEM % params % dumpFreq   = 1
#endif
      
      ! Construct the data structure that holds the derivative and interpolation matrices
      ! and the quadrature weights. This call will also perform the device copies.
      CALL myDGSEM % dGStorage % Build( myDGSEM % N, myDGSEM % params % nPlot, GAUSS, DG )
      
      ! Construct the roll-off filter matrix
      CALL myDGSEM % filter % Build( myDGSEM % dgStorage % interp % s,&
                                     myDGSEM % dgStorage % qWeight, &
                                     myDGSEM % N, myDGSEM % params % nCutoff )
                                     
      ! Load the mesh from the pc-mesh file and copy the mesh to the GPU
      CALL myDGSEM % BuildHexMesh( myRank )

      ALLOCATE( myDGSEM % dragProfile(0:myDGSEM % N, 0:myDGSEM % N, 0:myDGSEM % N, 1:myDGSEM % mesh % nElems) )
      myDGSEM % dragProfile = 0.0_prec

      ! Initialize the state and static "solutionstorage" data structures
      CALL myDGSEM % state % Build( myDGSEM % N, nEq, &
                                    myDGSEM % mesh % nElems )
                                    
      CALL myDGSEM % smoothState % Build( myDGSEM % N, nEq, &
                                          myDGSEM % mesh % nElems ) ! > Contains the smoothed SGS KE
                                          
      ! The "sgsCoeffs" attribute contains coefficients for the 
      ! subgrid scale parameterization. Currently, these coefficients
      ! are the eddy viscosities for the momentum equations (assuming
      ! isotropic turbulence), and the eddy diffusivities for the 
      ! potential temperature and density equations.
      CALL myDGSEM % sgsCoeffs % Build( myDGSEM % N, nEq-1, &
                                          myDGSEM % mesh % nElems )  
      ! Initially set all of the SGS coefficients to the "viscosity". In the event
      ! the Laplacian model is used, this will be the laplacian coefficient that is
      ! used for the momentum, potential temperature, and density equations.
      ! If another SGS model is used (e.g. SpectralEKE ), then these values will be
      ! overwritten in the "CalculateSGSCoefficients" routine.
      
      myDGSEM % sgsCoeffs % solution         = myDGSEM % params % viscosity
      myDGSEM % sgsCoeffs % boundarySolution = myDGSEM % params % viscosity

      CALL myDGSEM % static % Build( myDGSEM % N, nEq, &
                                     myDGSEM % mesh % nElems )

      ! Place conditionals on the construction of this attribute
      CALL myDGSEM % stressTensor % Build( myDGSEM % N, (nEq-1)*3, &
                                           myDGSEM % mesh % nElems )
                                           
#ifdef HAVE_CUDA

      ALLOCATE( myDGSEM % dragProfile_dev(0:myDGSEM % N, 0:myDGSEM % N, 0:myDGSEM % N, 1:myDGSEM % mesh % nElems) )
      myDGSEM % dragProfile_dev = 0.0_prec
      
      ! Copy values to the GPU variables
      myDGSEM % sgsCoeffs % solution_dev         = myDGSEM % sgsCoeffs % solution
      myDGSEM % sgsCoeffs % boundarySolution_dev = myDGSEM % sgsCoeffs % boundarySolution

      ! Set Device Constants
      R_dev         = myDGSEM % params % R
      Cv_dev        = myDGSEM % params % Cv
      P0_dev        = myDGSEM % params % P0
      hCapRatio_dev = ( myDGSEM % params % R + myDGSEM % params % Cv ) / myDGSEM % params % Cv
      rC_dev        = myDGSEM % params % R / ( myDGSEM % params % R + myDGSEM % params % Cv )
      g_dev         = myDGSEM % params % g
      dt_dev        = myDGSEM % params % dt
      viscLengthScale_dev = myDGSEM % params % viscLengthScale
      fRotX_dev     = myDGSEM % params % fRotX
      fRotY_dev     = myDGSEM % params % fRotY
      fRotZ_dev     = myDGSEM % params % fRotZ
      dScale_dev    = myDGSEM % params % dragScale
      Cd_dev        = myDGSEM % params % Cd

      nEq_dev     = nEq
      polydeg_dev = myDGSEM % N
      nEl_dev     = myDGSEM % mesh % nElems
      myRank_dev  = myRank

#endif

      ! Read the initial conditions, static state, and the boundary communicator
      CALL myDGSEM % ReadPickup( myDGSEM % params % iterInit, myRank )
      
      
#ifdef HAVE_CUDA
      istat = cudaDeviceSynchronize( )
      istat = cudaMemGetInfo( freebytes, totalbytes )
      PRINT*, '                                                         '
      PRINT*, '========================================================='
      PRINT*, '=============== GPU Device memory status ================'
      PRINT*, '                                                         '
      PRINT*, 'Total: ', REAL(totalbytes)/10.0**9, ' GB'
      PRINT*, 'Free : ', REAL(freebytes)/10.0**9, ' GB'
      PRINT*, 'Used : ', REAL(totalbytes-freebytes)/10.0**9, ' GB'
      PRINT*, '========================================================='
      PRINT*, '                                                         '
#endif

#ifdef HAVE_MPI
      CALL myDGSEM % ConstructCommTables( myRank, nProc )
      ALLOCATE( stateReqHandle(1:myDGSEM % nNeighbors*2), &
                stressReqHandle(1:myDGSEM % nNeighbors*2), &
                SGSReqHandle(1:myDGSEM % nNeighbors*2), &
                stateStats(MPI_STATUS_SIZE,1:myDGSEM % nNeighbors*2), &
                stressStats(MPI_STATUS_SIZE,1:myDGSEM % nNeighbors*2), &
                SGSStats(MPI_STATUS_SIZE,1:myDGSEM % nNeighbors*2) )
#endif

 END SUBROUTINE Build_Fluid
!
  SUBROUTINE Trash_Fluid( myDGSEM )

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   ! LOCAL
   INTEGER :: i
   
      PRINT*, 'S/R Trash_Fluid : Clearing memory.'
      CALL myDGSEM % state % Trash( )
      CALL myDGSEM % smoothState % Trash( )
      CALL myDGSEM % sgsCoeffs % Trash( )
      CALL myDGSEM % static % Trash( )
      CALL myDGSEM % stressTensor % Trash( )
#ifdef HAVE_MPI      
      DO i = 1, myDGSEM % nNeighbors
         DEALLOCATE( myDGSEM % mpiPackets(i) % sendStateBuffer )
         DEALLOCATE( myDGSEM % mpiPackets(i) % recvStateBuffer )
         DEALLOCATE( myDGSEM % mpiPackets(i) % sendStressBuffer )
         DEALLOCATE( myDGSEM % mpiPackets(i) % recvStressBuffer )
         DEALLOCATE( myDGSEM % mpiPackets(i) % sendSGSBuffer )
         DEALLOCATE( myDGSEM % mpiPackets(i) % recvSGSBuffer )
      ENDDO
      DEALLOCATE( myDGSEM % mpiPackets )
      DEALLOCATE( stateReqHandle, &
                  stressReqHandle, &
                  SGSReqHandle, &
                  stateStats, &
                  stressStats, &
                  SGSStats )
#endif
      CALL myDGSEM % dGStorage % Trash( )
      CALL myDGSEM % mesh % Trash( )
      CALL myDGSEM % extComm % Trash( )

      DEALLOCATE( myDGSEM % dragProfile )
      
      DEALLOCATE( myDGSEM % prescribedState )
      DEALLOCATE( myDGSEM % externalState )
      DEALLOCATE( myDGSEM % prescribedStress )
      DEALLOCATE( myDGSEM % externalStress )
      DEALLOCATE( myDGSEM % externalSGS )
      
#ifdef HAVE_CUDA
      DEALLOCATE( myDGSEM % dragProfile_dev )
      DEALLOCATE( myDGSEM % prescribedState_dev )
      DEALLOCATE( myDGSEM % externalState_dev )
      DEALLOCATE( myDGSEM % prescribedStress_dev )
      DEALLOCATE( myDGSEM % externalStress_dev )
      DEALLOCATE( myDGSEM % externalSGS_dev )
#endif

#ifdef TESTING
      CALL mdi % Trash( )
#endif

 END SUBROUTINE Trash_Fluid
!
 SUBROUTINE BuildHexMesh_Fluid( myDGSEM, myRank )

   IMPLICIT NONE
   CLASS( Fluid ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)           :: myRank
   ! Local
   CHARACTER(4) :: rankChar

      WRITE( rankChar, '(I4.4)' )myRank
      PRINT*,'Module FluidClass.f90 : S/R BuildHexMesh :'

      PRINT*, 'Reading mesh from '//TRIM(myDGSEM % params % PeaceMeshFile)//'.'//rankChar//'.pc.mesh '
      
      ! This loads in the mesh from the "pc-mesh file" and sets up the device arrays for the mesh
      CALL myDGSEM % mesh % ReadPeaceMeshFile( TRIM(myDGSEM % params % PeaceMeshFile)//'.'//rankChar )
      
      ! Multiply the mesh positions to scale the size of the mesh
      CALL myDGSEM % mesh % ScaleTheMesh( myDGSEM % dgStorage % interp, &
                                          myDGSEM % params % xScale, &
                                          myDGSEM % params % yScale, &
                                          myDGSEM % params % zScale )

#ifdef HAVE_CUDA
      ! Here we update the device arrays on the mesh
      CALL myDGSEM % mesh % UpdateDevice( myDGSEM % N )
#endif

 END SUBROUTINE BuildHexMesh_Fluid
!
#ifdef HAVE_MPI
 SUBROUTINE ConstructCommTables_Fluid( myDGSEM, myRank, nProc )
   IMPLICIT NONE
   CLASS( Fluid ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)           :: myRank, nProc
   ! Local
   INTEGER    :: sharedFaceCount(0:nProc-1)
   INTEGER    :: iFace, bID, iNeighbor
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, nmsg
   

      ! Count up the number of neighboring ranks
      ALLOCATE( myDGSEM % rankTable(0:nProc-1) )
      myDGSEM % rankTable = 0
      sharedFaceCount     = 0
      DO bID = 1, myDGSEM % extComm % nBoundaries
         p2 = myDGSEM % extComm % extProcIDS(bID)
         IF( p2 /= myRank )THEN
            myDGSEM % rankTable(p2) = 1
            sharedFaceCount(p2) = sharedFaceCount(p2)+1
         ENDIF
      ENDDO
      myDGSEM % nNeighbors = SUM( myDGSEM % rankTable )
      PRINT*, '  S/R ConstructCommTables : Found', myDGSEM % nNeighbors, 'neighbors for Rank', myRank
      
      ALLOCATE( myDGSEM % mpiPackets(1:myDGSEM % nNeighbors) )
      
      ! For each neighbor, set the neighbor's rank
      iNeighbor = 0
      DO p2 = 0, nProc-1
         IF( myDGSEM % rankTable(p2) == 1 )THEN
            iNeighbor = iNeighbor + 1
            myDGSEM % mpiPackets(iNeighbor) % neighborRank = p2
            myDGSEM % rankTable(p2) = iNeighbor
         ENDIF
      ENDDO
      
      DO iNeighbor = 1, myDGSEM % nNeighbors
      
         p2 = myDGSEM % mpiPackets(iNeighbor) % neighborRank

         myDGSEM % mpiPackets(iNeighbor) % bufferSize = sharedFaceCount(p2)
         
         ALLOCATE( myDGSEM % mpiPackets(iNeighbor) % recvStateBuffer(0:myDGSEM % N, 0:myDGSEM % N, 1:nEq,1:sharedFaceCount(p2)), &
                   myDGSEM % mpiPackets(iNeighbor) % sendStateBuffer(0:myDGSEM % N, 0:myDGSEM % N, 1:nEq,1:sharedFaceCount(p2)), &
                   myDGSEM % mpiPackets(iNeighbor) % recvStressBuffer(0:myDGSEM % N, 0:myDGSEM % N, 1:(nEq-1)*3,1:sharedFaceCount(p2)), &
                   myDGSEM % mpiPackets(iNeighbor) % sendStressBuffer(0:myDGSEM % N, 0:myDGSEM % N, 1:(nEq-1)*3,1:sharedFaceCount(p2)), &
                   myDGSEM % mpiPackets(iNeighbor) % sendSGSBuffer(0:myDGSEM % N, 0:myDGSEM % N, 1:(nEq-1),1:sharedFaceCount(p2)), &
                   myDGSEM % mpiPackets(iNeighbor) % recvSGSBuffer(0:myDGSEM % N, 0:myDGSEM % N, 1:(nEq-1),1:sharedFaceCount(p2)) )
                   
         myDGSEM % mpiPackets(iNeighbor) % bufferCounter = 0
      ENDDO
      
   
 END SUBROUTINE ConstructCommTables_Fluid
#endif
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
  SUBROUTINE ForwardStepRK3_Fluid( myDGSEM, nT, myRank )

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)         :: nT
   INTEGER, INTENT(in)         :: myRank
#ifdef HAVE_CUDA
   ! LOCAL
   REAL(prec) :: t, dt
   REAL(prec),DEVICE :: G3D(0:myDGSEM % N,&
                            0:myDGSEM % N,&
                            0:myDGSEM % N,&
                            1:nEq,&
                            1:myDGSEM % mesh % nElems) 
   REAL(prec), DEVICE :: rk3_a_dev(1:3), rk3_g_dev(1:3)
   INTEGER    :: iT, m, iStat
   TYPE(dim3) :: grid, tBlock
  
     ! How should we pick the thread and block size

      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) )
      grid = dim3(myDGSEM % mesh % nElems, nEq-1, 1) ! Because nElems can reach a rather large size, it must remain the first dimension of the grid (see pgaccelinfo for MAximum Grid Size)
     

      dt = myDGSEM % params % dt
      
      rk3_a_dev = rk3_a
      rk3_g_dev = rk3_g
      DO iT = 1, nT
     
         G3D  = 0.0_prec
         DO m = 1,3 ! Loop over RK3 steps
            
            t = myDGSEM % simulationTime + rk3_b(m)*dt
            CALL myDGSEM % GlobalTimeDerivative( t, myRank )
            
            !
            CALL UpdateG3D_CUDAKernel<<<grid,tBlock>>>( G3D, rk3_a_dev(m), rk3_g_dev(m), &
                                                        myDGSEM % state % solution_dev, &
                                                        myDGSEM % state % tendency_dev, &
                                                        myDGSEM % stressTensor % tendency_dev )
            
            ! Calculate the pressure
            CALL myDGSEM % EquationOfState( )
            
         ENDDO ! m, loop over the RK3 steps
      
      ENDDO
      myDGSEM % simulationTime = myDGSEM % simulationTime + nT*dt
#else
   REAL(prec) :: t, dt, rk3_a_local, rk3_g_local
   REAL(prec) :: G3D(0:myDGSEM % N,&
                     0:myDGSEM % N,&
                     0:myDGSEM % N,&
                     1:myDGSEM % nEq-1,&
                     1:myDGSEM % mesh % nElems) 
   INTEGER    :: m, iEl, iT, i, j, k, iEq

      dt = myDGSEM % params % dt
      

      DO iT = 1, nT
     
         !$OMP DO
         DO iEl = 1, myDGSEM % mesh % nElems
            DO iEq = 1, nEq-1
              DO k = 0, myDGSEM % N
                 DO j = 0, myDGSEM % N
                    DO i = 0, myDGSEM % N
                       G3D(i,j,k,iEq,iEl) = ZERO
                    ENDDO
                 ENDDO
              ENDDO
            ENDDO
         ENDDO
         !$OMP ENDDO
         
         DO m = 1,3 ! Loop over RK3 steps
            
            t = myDGSEM % simulationTime + rk3_b(m)*dt
            CALL myDGSEM % GlobalTimeDerivative( t, myRank )
            
            rk3_a_local = rk3_a(m)
            rk3_g_local = rk3_g(m)

            !$OMP DO
            DO iEl = 1, myDGSEM % mesh % nElems
               DO iEq = 1, nEq-1
                  DO k = 0, myDGSEM % N
                     DO j = 0, myDGSEM % N
                        DO i = 0, myDGSEM % N
                           G3D(i,j,k,iEq,iEl) = rk3_a_local*G3D(i,j,k,iEq,iEl) + myDGSEM % state % tendency(i,j,k,iEq,iEl) +&
                                                                                 myDGSEM % stressTensor % tendency(i,j,k,iEq,iEl)
                           myDGSEM % state % solution(i,j,k,iEq,iEl) = myDGSEM % state % solution(i,j,k,iEq,iEl) + &
                                                    rk3_g_local*dt*G3D(i,j,k,iEq,iEl)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO ! iEl, loop over all of the elements
            !$OMP ENDDO
            
            ! Calculate the internal energy and the pressure
            CALL myDGSEM % EquationOfState( )
            
         ENDDO ! m, loop over the RK3 steps
            
      ENDDO
      myDGSEM % simulationTime = myDGSEM % simulationTime + nT*myDGSEM % params % dt 
#endif          


 END SUBROUTINE ForwardStepRK3_Fluid
!
 SUBROUTINE GlobalTimeDerivative_Fluid( myDGSEM, tn, myRank ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)      :: tn
   INTEGER, INTENT(in)         :: myRank


! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
!  If SpectralFiltering is used as the subgridscale model, then the spectral
!  filtering matrix (specified in src/filtering/RollOffFilter_Class.f90) is used
!  to smooth the solution variables before proceeding.

      IF( myDGSEM % params % SubGridModel == SpectralFiltering )THEN
         CALL myDGSEM % CalculateSmoothedState( .TRUE. )

#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % state % solution = myDGSEM % state % solution_dev
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'CalculateSmoothedState', &
                            'Smooth State for Spectral Filtering', &
                             SIZE(myDGSEM % state % solution), &
                             PACK(myDGSEM % state % solution,.TRUE.) )
      ENDIF
#endif
      ENDIF

! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
!  Here, the solution within each element is interpolated to the faces of each
!  element in order to prepare for computing the external state for enforcing
!  boundary conditions, Riemann Fluxes, and MPI data exchanges that need to
!  occur.

      CALL myDGSEM % CalculateBoundarySolution( ) 

#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % state % boundarySolution = myDGSEM % state % boundarySolution_dev
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'CalculateBoundarySolution', &
                            'Interpolation to element boundaries', &
                             SIZE(myDGSEM % state % boundarySolution), &
                             PACK(myDGSEM % state % boundarySolution,.TRUE.) )
      ENDIF
#endif

! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! When MPI is used, the boundary solutions that are stored on faces shared with
! a neighboring rank are passed to that neighboring rank. Additionally, the
! perceived "external state" is received from the neighboring rank. Calling this
! routine is dependent on the result of CalculateBoundarySolutio, but can be
! done at the same time as UpdateExternalState; doing so should hide some
! communication costs.

#ifdef HAVE_MPI
      CALL myDGSEM % MPI_StateExchange( myRank ) 
#endif
! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! The boundary solutions are used to calculate the external states that, when
! accompanied with a Riemann Solver, enforce boundary conditions. Calling this
! routine is dependent on the result of CalculateBoundarySolution

      CALL myDGSEM % UpdateExternalState( tn, myRank ) 
#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % state % externalState = myDGSEM % state % externalState_dev
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'UpdateExternalState', &
                            'Update of Boundary Conditions', &
                             SIZE(myDGSEM % externalState), &
                             PACK(myDGSEM % externalState,.TRUE.) )
      ENDIF
#endif


! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! The inviscid fluxes (advection + pressure) through the faces of each element
! are estimated here using a (linear) Lax-Friedrich's upwind solver. In order to
! call this routine, CalculateBoundarySolution, UpdateExternalState, and
! MPI_StateExchange must have been completed.

      CALL myDGSEM % InternalFaceFlux( )

#ifdef HAVE_MPI
      CALL myDGSEM % FinalizeMPI_StateExchange( myRank ) 
#endif

      CALL myDGSEM % BoundaryFaceFlux( )

#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % state % boundaryFlux = myDGSEM % state % boundaryFlux_dev
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'FaceFlux', &
                            'Update of boundary fluxes', &
                             SIZE(myDGSEM % state % boundaryFlux), &
                             PACK(myDGSEM % state % boundaryFlux,.TRUE.) )
      ENDIF
#endif
      
! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! If the SpectralEKE or Laplacian subgridscale models are used, a Laplacian-like
! operator is used to diffuse momentum and heat. When the Laplacian model is
! used a fixed viscosity is specified in runtime.params and a Rayleigh number of
! 1 is assumed (viscosity = diffusivity). When the SpectralEKE model is used,
! the viscosity coefficient is diagnosed using a Smagorinksy closure, where the
! unresolved kinetic energy is diagnosed from a highpass spectral filter.
!
      IF( myDGSEM % params % SubGridModel == SpectralEKE .OR. &
          myDGSEM % params % SubGridModel == Laplacian ) THEN
          
         IF( myDGSEM % params % SubGridModel == SpectralEKE )THEN ! 

! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! In the Spectral-EKE model, the under-resolved state is diagnosed from a high
! pass spectral filter (the difference between the state and smoothed state).
! Here, we first calculate the smoothed state and store it in the smoothedState
! attribute. This subroutine call has no dependence to any other within this
! subroutine.

            CALL myDGSEM % CalculateSmoothedState( .FALSE. )              

#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % smoothState % solution = myDGSEM % smoothState % solution_dev
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'CalculateSmoothedState', &
                            'Smooth State for Spectral EKE', &
                             SIZE(myDGSEM % smoothState % solution), &
                             PACK(myDGSEM % smoothState % solution,.TRUE.) )
      ENDIF
#endif

! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! The high passed solution is used to diagnose an isotropic viscosity
! coefficient, similar to a Smagorinksy closure and similar to the closure in 
!
!  J. Sauer (2013), "Towards Improved Capability and Confidence in Coupled
!  Atmospheric and Wildland Fire Modeling"
!
! The main difference in this work, is in the diagnosis of the SGS Kinetic
! Energy from the high pass filtered solution.
! This routine depends on the results from CalculateSmoothedState.

            CALL myDGSEM % CalculateSGSCoefficients( ) 

#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % sgsCoeffs % solution = myDGSEM % sgsCoeffs % solution
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'CalculateSGSCoefficients', &
                            'Estimate viscosity and diffusivity', &
                             SIZE(myDGSEM % sgsCoeffs % solution), &
                             PACK(myDGSEM % sgsCoeffs % solution,.TRUE.) )
      ENDIF
#endif

! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
! The viscosity coefficient that is calculated is now interpolated to the faces
! of each element so that the viscous flux can later be computed. This routine
! depends on the result of CalculateSGSCoefficients.

            CALL myDGSEM % CalculateBoundarySGS( ) 

#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % sgsCoeffs % boundarySolution = myDGSEM % sgsCoeffs % boundarySolution_dev
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'CalculateBoundarySGS', &
                            'Interpolate viscosity to element faces', &
                             SIZE(myDGSEM % sgsCoeffs % boundarySolution), &
                             PACK(myDGSEM % sgsCoeffs % boundarySolution,.TRUE.) )
      ENDIF
#endif

! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! The viscosity coefficients are exchanged with neighboring ranks that share
! common faces. MPI_SGSExchange can be run simulataneously with
! CalculateStressTensor, CalculateBoundaryStress, UpdateExternalStress, and the
! MPI_StressExchange. The viscosity coefficients that are exchanged are not
! needed until StressFlux

#ifdef HAVE_MPI
            CALL myDGSEM % MPI_SGSExchange( myRank )  
#endif
         ENDIF
! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! Now, the internal solution and the stress tensor boundary flux can be pieced
! together to calculate gradients in the velocity and potential temperature.
! This routine depends on the result of FaceFlux (stressTensor % boundaryFlux)   

         CALL myDGSEM % CalculateStressTensor( )                  

#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % stressTensor % solution = myDGSEM % stressTensor % solution_dev
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'CalculateStressTensor', &
                            'Gradients of velocity and temperature', &
                             SIZE(myDGSEM % stressTensor % solution), &
                             PACK(myDGSEM % stressTensor % solution,.TRUE.) )
      ENDIF
#endif


#ifdef HAVE_MPI
      IF( myDGSEM % params % SubGridModel == SpectralEKE )THEN ! 
         CALL myDGSEM % FinalizeMPI_SGSExchange( myRank )  
      ENDIF
#endif

! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! The stress tensor values are interpolated to the faces of each element to
! prepare for the calculation of the divergence of the viscous fluxes. This
! routine depends on the result of CalculateStressTensor.

         CALL myDGSEM % CalculateBoundaryStress( )

#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % stressTensor % boundarySolution = myDGSEM % stressTensor % boundarySolution_dev
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'CalculateBoundaryStress', &
                            'Interpolate stress tensor to element faces', &
                             SIZE(myDGSEM % stressTensor % boundarySolution), &
                             PACK(myDGSEM % stressTensor % boundarySolution,.TRUE.) )
      ENDIF
#endif

! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! Stress tensor values are exchanged with neighboring ranks along shared faces.
! This routine depends on the result of CalculateBoundaryStress, but can be run
! at the same time as UpdateExternalStress.

#ifdef HAVE_MPI
         CALL myDGSEM % MPI_StressExchange( myRank )
#endif

! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! Now that the stress tensor is available on element faces, boundary conditions
! can be applied by setting the external stress tensor state. This routine
! depends on the result of CalculateBoundaryStress. Note that this routine can
! be run simultaneously with the MPI_StressExchange

         CALL myDGSEM % UpdateExternalStress( tn, myRank ) 

#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % externalStress = myDGSEM % externalStress_dev
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'UpdateExternalStress', &
                            'Apply Stress Boundary Conditions', &
                             SIZE(myDGSEM % externalStress), &
                             PACK(myDGSEM % externalStress,.TRUE.) )
      ENDIF
#endif


! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! Using the boundary and the external state for the stress tensor, the viscous
! fluxes are estimated using a Bassi-Rebay flux that averages neighboring values
! of the stress tensor plus the jump in the solution weighted by a spatial
! wave-number. This routine depends on the result of the UpdateExternalStress
! and the MPI_StressExchange.

         CALL myDGSEM % InternalStressFlux( )
#ifdef HAVE_MPI
         CALL myDGSEM % FinalizeMPI_StressExchange( myRank )
#endif
         CALL myDGSEM % BoundaryStressFlux( )
 
#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % stressTensor % boundaryFlux = myDGSEM % stressTensor % boundaryFlux_dev
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'StressFlux', &
                            'Estimate Viscous Stress Flux', &
                             SIZE(myDGSEM % stressTensor % boundaryFlux), &
                             PACK(myDGSEM % stressTensor % boundaryFlux,.TRUE.) )
      ENDIF
#endif

! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! With the boundary stress flux and the internal stress tensor values, the
! divergence of the stress tensor can be calculated, giving the viscous tendency
! for the momentum and the potential temperature. This routine depends on the
! result of StressFlux (and the dependencies of StressFlux), but can be done
! simultaneously with the MappedTimeDerivative.

         CALL myDGSEM % StressDivergence( ) 

#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % stressTensor % tendency = myDGSEM % stressTensor % tendency_dev
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'StressDivergence', &
                            'Tendency due to viscous terms', &
                             SIZE(myDGSEM % stressTensor % tendency), &
                             PACK(myDGSEM % stressTensor % tendency,.TRUE.) )
      ENDIF
#endif

! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 

      ENDIF
      
! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 
!
! Once the inviscid fluxes through the faces are calculated, and the internal
! state is known, the tendency due to the inviscid flux terms and
! nonconservative source terms is calculated here. This routine depends on the
! result of FaceFlux, but can be done at the same time as StressDivergence
      
      CALL myDGSEM % MappedTimeDerivative( )

#ifdef TESTING
      IF( myRank == 0 )THEN
#ifdef CUDA
         myDGSEM % state % tendency = myDGSEM % state % tendency_dev
#endif
         CALL mdi % Update( 'Fluid_Class.f90', &
                            'MappedTimeDerivative', &
                            'Tendency due to inviscid and source terms', &
                             SIZE(myDGSEM % state % tendency), &
                             PACK(myDGSEM % state % tendency,.TRUE.) )
      ENDIF
#endif

! ----------------------------------------------------------------------------- ! 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- ! 

#ifdef TESTING
      IF( myRank == 0 )THEN
         CALL mdi % Write_ModelDataInstances( 'SELF-Fluid' ) 
      ENDIF
#endif
      
      
 END SUBROUTINE GlobalTimeDerivative_Fluid
!
  SUBROUTINE CalculateSmoothedState_Fluid( myDGSEM, overwriteState )
 
   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM 
   LOGICAL, INTENT(in)         :: overWriteState
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) )
      grid = dim3(myDGSEM % mesh % nElems, nEq-1, 1)
  
      IF( overWriteState )THEN 
         CALL CalculateSmoothedState_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                                  myDGSEM % filter % filterMat_dev, &
                                                                  myDGSEM % state % solution_dev )
      ELSE
         CALL CalculateSmoothedState_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                                  myDGSEM % filter % filterMat_dev, &
                                                                  myDGSEM % smoothState % solution_dev )
      ENDIF

#else
   ! Local
   INTEGER :: iEl, iEq, i, j, k, ii, jj, kk
   REAL(prec) :: uijk, uij, ui
   
      IF( overWriteState )THEN

         !$OMP DO
         DO iEl = 1, myDGSEM % mesh % nElems
            DO iEq = 1, nEq-1
               DO k = 0, myDGSEM % N
                  DO j = 0, myDGSEM % N
                     DO i = 0, myDGSEM % N
                     
                        uijk = 0.0_prec
                        DO kk = 0, myDGSEM % N
                        
                           uij = 0.0_prec
                           DO jj = 0, myDGSEM % N
                              
                              ui = 0.0_prec
                              DO ii = 0, myDGSEM % N
                                 ui = ui + myDGSEM % filter % filterMat(ii,i)*&
                                           myDGSEM % state % solution(ii,jj,kk,iEq,iEl)
                              ENDDO
                              
                              uij = uij + myDGSEM % filter % filterMat(jj,j)*ui
                           ENDDO
                           
                           uijk = uijk + myDGSEM % filter % filterMat(kk,k)*uij
                           
                        ENDDO
                        
                        myDGSEM % state % solution(i,j,k,iEq,iEl) = uijk
                        
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
        !$OMP ENDDO
 
     ELSE

         !$OMP DO
         DO iEl = 1, myDGSEM % mesh % nElems
            DO iEq = 1, nEq-1
               DO k = 0, myDGSEM % N
                  DO j = 0, myDGSEM % N
                     DO i = 0, myDGSEM % N
                     
                        uijk = 0.0_prec
                        DO kk = 0, myDGSEM % N
                        
                           uij = 0.0_prec
                           DO jj = 0, myDGSEM % N
                              
                              ui = 0.0_prec
                              DO ii = 0, myDGSEM % N
                                 ui = ui + myDGSEM % filter % filterMat(ii,i)*&
                                           myDGSEM % state % solution(ii,jj,kk,iEq,iEl)
                              ENDDO
                              
                              uij = uij + myDGSEM % filter % filterMat(jj,j)*ui
                           ENDDO
                           
                           uijk = uijk + myDGSEM % filter % filterMat(kk,k)*uij
                           
                        ENDDO
                        
                        myDGSEM % smoothState % solution(i,j,k,iEq,iEl) = uijk
                        
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
        !$OMP ENDDO

     ENDIF
      
#endif

 END SUBROUTINE CalculateSmoothedState_Fluid
! 
 SUBROUTINE CalculateSGSCoefficients_Fluid( myDGSEM )
 
   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                    4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) )
      grid = dim3(myDGSEM % mesh % nElems, 1, 1)
   
      CALL CalculateSGSCoefficients_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                                 myDGSEM % static % solution_dev, &
                                                                 myDGSEM % smoothState % solution_dev, &
                                                                 myDGSEM % filter % filterMat_dev, &
                                                                 myDGSEM % sgsCoeffs % solution_dev )
#else
   ! Local
   INTEGER :: iEl, i, j, k, m, ii, jj, kk
   REAL(prec) :: sgsKE, uijk, uij, ui
   REAL(prec) :: KE(0:myDGSEM % N,0:myDGSEM % N,0:myDGSEM % N)
   

      !$OMP DO PRIVATE( KE )
      DO iEl = 1, myDGSEM % mesh % nElems
      
         ! Here, the SGS Kinetic energy is calculated using the 
         ! "high wavenumber" component of the velocity field.
         ! This component is defined (here) as the difference
         ! between the full solution and the smoothed solution.
         DO k = 0, myDGSEM % N
            DO j = 0, myDGSEM % N
               DO i = 0, myDGSEM % N
               
                  sgsKE = 0.0_prec
                  DO m = 1, 3  
                     sgsKE = sgsKE + &
                             ( myDGSEM % state % solution(i,j,k,m,iEl)/&
                               (myDGSEM % state % solution(i,j,k,4,iEl)+myDGSEM % static % solution(i,j,k,4,iEl)) - &
                               myDGSEM % smoothState % solution(i,j,k,m,iEl)/&
                               (myDGSEM % smoothState % solution(i,j,k,4,iEl)+myDGSEM % static % solution(i,j,k,4,iEl)) )**2
                  ENDDO
                  KE(i,j,k) = 0.5_prec*sgsKE
                  
               ENDDO
            ENDDO
         ENDDO

#ifdef VIZ         
         ! Smooth the subgrid scale Kinetic energy
         DO k = 0, myDGSEM % N
            DO j = 0, myDGSEM % N
               DO i = 0, myDGSEM % N
                  
                  uijk = 0.0_prec
                  DO kk = 0, myDGSEM % N
                  
                     uij = 0.0_prec
                     DO jj = 0, myDGSEM % N
                        
                        ui = 0.0_prec
                        DO ii = 0, myDGSEM % N
                           ui = ui + myDGSEM % filter % filterMat(ii,i)*&
                                     KE(ii,jj,kk)
                        ENDDO
                        
                        uij = uij + myDGSEM % filter % filterMat(jj,j)*ui
                     ENDDO
                     
                     uijk = uijk + myDGSEM % filter % filterMat(kk,k)*uij
                     
                  ENDDO
                  
                  ! Here, we store the smoothed SGS kinetic energy, in
                  ! case we would like to visualize the data later
                  myDGSEM % smoothState % solution(i,j,k,6,iEl) = ABS(uijk)
                  
               ENDDO
            ENDDO
         ENDDO
#endif         
         
         ! Now we calculate the viscosity and diffusivities (currently assumes isotropic and low mach number)
         DO k = 0, myDGSEM % N
            DO j = 0, myDGSEM % N
               DO i = 0, myDGSEM % N
                  DO m = 1, nEq-1
                  
                     myDGSEM % sgsCoeffs % solution(i,j,k,m,iEl) = 0.09_prec*&
                                    myDGSEM % params % viscLengthScale*sqrt( KE(i,j,k) )
                                      
                  ENDDO
               ENDDO
            ENDDO
         ENDDO 
      ENDDO
      !$OMP ENDDO

#endif

 END SUBROUTINE CalculateSGSCoefficients_Fluid
!
 SUBROUTINE CalculateBoundarySGS_Fluid( myDGSEM ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
     tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                   4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                   nEq-1 )
     grid = dim3(myDGSEM % mesh % nElems, 1, 1)  
     CALL CalculateBoundarySGS_CUDAKernel<<<grid, tBlock>>>( myDGSEM % sgsCoeffs % solution_dev, &
                                                              myDGSEM % dgStorage % bMat_dev, &
                                                              myDGSEM % sgsCoeffs % boundarySolution_dev )
#else
   ! Local
   INTEGER :: iEq, iEl, i, j, k


      !$OMP DO
      DO iEl = 1, myDGSEM % mesh % nElems
         DO iEq = 1, myDGSEM % nEq-1
            DO k = 0, myDGSEM % N
               DO j = 0, myDGSEM % N
               
                  myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,WEST,iEl)   = 0.0_prec
                  myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,EAST,iEl)   = 0.0_prec
                  myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,SOUTH,iEl)  = 0.0_prec
                  myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,NORTH,iEl)  = 0.0_prec
                  myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,BOTTOM,iEl) = 0.0_prec
                  myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,TOP,iEL)    = 0.0_prec
    
                  DO i = 0, myDGSEM % N

                     myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,WEST,iEl)  = &
                        myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,WEST,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,0)*myDGSEM % sgsCoeffs % solution(i,j,k,iEq,iEl)

                     myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,EAST,iEl)  = &
                        myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,EAST,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,1)*myDGSEM % sgsCoeffs % solution(i,j,k,iEq,iEl)

                     myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,SOUTH,iEl)  = &
                        myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,SOUTH,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,0)*myDGSEM % sgsCoeffs % solution(j,i,k,iEq,iEl)

                     myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,NORTH,iEl)   = &
                        myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,NORTH,iEl)  +  &
                        myDGSEM % dgStorage % bMat(i,1)*myDGSEM % sgsCoeffs % solution(j,i,k,iEq,iEl)

                     myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,BOTTOM,iEl)  = &
                        myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,BOTTOM,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,0)*myDGSEM % sgsCoeffs % solution(j,k,i,iEq,iEl)

                     myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,TOP,iEl)  = &
                        myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,TOP,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,1)*myDGSEM % sgsCoeffs % solution(j,k,i,iEq,iEl)

                  ENDDO
                  
                  ! Ensure positivity of the subgridscale viscosity/diffusivity coefficients
                  myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,1:6,iEl) = ABS(myDGSEM % sgsCoeffs % boundarySolution(j,k,iEq,1:6,iEl))
               
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !$OMP ENDDO
      

#endif
      
 END SUBROUTINE CalculateBoundarySGS_Fluid
!
 SUBROUTINE UpdateExternalSGS_Fluid( myDGSEM, myRank ) ! ////////// !

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)         :: myRank
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
     tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                   4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                   1 )
     grid = dim3(myDGSEM % nBoundaryFaces,nEq-1,1)  
     
     CALL UpdateExternalSGSCoeffs_CUDAKernel<<<grid, tBlock>>>( myDGSEM % extComm % boundaryIDs_dev, &       ! I
                                                            myDGSEM % mesh % faces_dev % elementIDs, &   ! I
                                                            myDGSEM % mesh % faces_dev % elementSides, & ! I
                                                            myDGSEM % extComm % extProcIDs_dev, &           ! I
                                                            myDGSEM % externalSGS_dev, &                    ! O
                                                            myDGSEM % sgsCoeffs % boundarySolution_dev, &   ! I  
                                                            myDGSEM % mesh % geom_dev % nHat_dev )           ! I
#else
   ! Local
   INTEGER    :: iEl, iFace, bFaceID, i, j, k, iEq
   INTEGER    :: iFace2, p2
   INTEGER    :: e1, e2, s1, s2
   
      !$OMP DO
      DO iFace = 1, myDGSEM % nBoundaryFaces

         iFace2 = myDGSEM % extComm % boundaryIDs( iFace ) ! Obtain the process-local face id for this boundary-face id
         e1     = myDGSEM % mesh % Faces(iFace2) % elementIDs(1)
         s1     = myDGSEM % mesh % Faces(iFace2) % elementSides(1)
         e2     = myDGSEM % mesh % Faces(iFace2) % elementIDs(2)
         p2     = myDGSEM % extComm % extProcIDs( iFace )
         
         IF( p2 == myRank )THEN ! Enforce no boundary flux due to the fluid stress
            DO j = 0, myDGSEM % N 
               DO i = 0, myDGSEM % N
                  DO iEq = 1, myDGSEM % nEq-1
                        myDGSEM % externalSGS(i,j,iEq,iFace) = myDGSEM % sgsCoeffs % boundarySolution(i,j,iEq,s1,e1)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         
      ENDDO
      !$OMP ENDDO 

#endif

 END SUBROUTINE UpdateExternalSGS_Fluid
!
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!                                                                                                 !
!  This section of code contains routines for interpolating the the solutions within each element !
!  to the element faces. Included are routines for the static and state attributes of the Fluid   !
!  data structure.                                                                                !
!                                                                                                 !
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!
  SUBROUTINE CalculateStaticBoundarySolution_Fluid( myDGSEM )

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
     tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                   4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                   nEq )
     grid = dim3(myDGSEM % mesh % nElems, 1, 1) 
      
     CALL CalculateBoundarySolution_CUDAKernel<<<grid, tBlock>>>( myDGSEM % static % solution_dev, &
                                                                  myDGSEM % dgStorage % bMat_dev, &
                                                                  myDGSEM % static % boundarySolution_dev )
#else
   ! Local
   INTEGER :: iEq, iEl, i, j, k

      !$OMP DO
      DO iEl = 1, myDGSEM % mesh % nElems
         DO iEq = 1, myDGSEM % nEq
            DO k = 0, myDGSEM % N
               DO j = 0, myDGSEM % N
               
                  myDGSEM % static % boundarySolution(j,k,iEq,WEST,iEl)   = 0.0_prec
                  myDGSEM % static % boundarySolution(j,k,iEq,EAST,iEl)   = 0.0_prec
                  myDGSEM % static % boundarySolution(j,k,iEq,SOUTH,iEl)  = 0.0_prec
                  myDGSEM % static % boundarySolution(j,k,iEq,NORTH,iEl)  = 0.0_prec
                  myDGSEM % static % boundarySolution(j,k,iEq,BOTTOM,iEl) = 0.0_prec
                  myDGSEM % static % boundarySolution(j,k,iEq,TOP,iEL)    = 0.0_prec
    
                  DO i = 0, myDGSEM % N

                     myDGSEM % static % boundarySolution(j,k,iEq,WEST,iEl)  = &
                        myDGSEM % static % boundarySolution(j,k,iEq,WEST,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,0)*myDGSEM % static % solution(i,j,k,iEq,iEl)

                     myDGSEM % static % boundarySolution(j,k,iEq,EAST,iEl)  = &
                        myDGSEM % static % boundarySolution(j,k,iEq,EAST,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,1)*myDGSEM % static % solution(i,j,k,iEq,iEl)

                     myDGSEM % static % boundarySolution(j,k,iEq,SOUTH,iEl)  = &
                        myDGSEM % static % boundarySolution(j,k,iEq,SOUTH,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,0)*myDGSEM % static % solution(j,i,k,iEq,iEl)

                     myDGSEM % static % boundarySolution(j,k,iEq,NORTH,iEl)   = &
                        myDGSEM % static % boundarySolution(j,k,iEq,NORTH,iEl)  +  &
                        myDGSEM % dgStorage % bMat(i,1)*myDGSEM % static % solution(j,i,k,iEq,iEl)

                     myDGSEM % static % boundarySolution(j,k,iEq,BOTTOM,iEl)  = &
                        myDGSEM % static % boundarySolution(j,k,iEq,BOTTOM,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,0)*myDGSEM % static % solution(j,k,i,iEq,iEl)

                     myDGSEM % static % boundarySolution(j,k,iEq,TOP,iEl)  = &
                        myDGSEM % static % boundarySolution(j,k,iEq,TOP,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,1)*myDGSEM % static % solution(j,k,i,iEq,iEl)

                  ENDDO
               
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !$OMP ENDDO

#endif

 END SUBROUTINE CalculateStaticBoundarySolution_Fluid
!
 SUBROUTINE CalculateBoundarySolution_Fluid( myDGSEM ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
     tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                   4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                   nEq )
     grid = dim3(myDGSEM % mesh % nElems, 1, 1)  
     CALL CalculateBoundarySolution_CUDAKernel<<<grid, tBlock>>>( myDGSEM % state % solution_dev, &
                                                                  myDGSEM % dgStorage % bMat_dev, &
                                                                  myDGSEM % state % boundarySolution_dev )
#else
   ! Local
   INTEGER :: iEq, iEl, i, j, k

      !$OMP DO
      DO iEl = 1, myDGSEM % mesh % nElems
         DO iEq = 1, myDGSEM % nEq
            DO k = 0, myDGSEM % N
               DO j = 0, myDGSEM % N
               
                  myDGSEM % state % boundarySolution(j,k,iEq,WEST,iEl)   = 0.0_prec
                  myDGSEM % state % boundarySolution(j,k,iEq,EAST,iEl)   = 0.0_prec
                  myDGSEM % state % boundarySolution(j,k,iEq,SOUTH,iEl)  = 0.0_prec
                  myDGSEM % state % boundarySolution(j,k,iEq,NORTH,iEl)  = 0.0_prec
                  myDGSEM % state % boundarySolution(j,k,iEq,BOTTOM,iEl) = 0.0_prec
                  myDGSEM % state % boundarySolution(j,k,iEq,TOP,iEL)    = 0.0_prec
    
                  DO i = 0, myDGSEM % N

                     myDGSEM % state % boundarySolution(j,k,iEq,WEST,iEl)  = &
                        myDGSEM % state % boundarySolution(j,k,iEq,WEST,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,0)*myDGSEM % state % solution(i,j,k,iEq,iEl)

                     myDGSEM % state % boundarySolution(j,k,iEq,EAST,iEl)  = &
                        myDGSEM % state % boundarySolution(j,k,iEq,EAST,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,1)*myDGSEM % state % solution(i,j,k,iEq,iEl)

                     myDGSEM % state % boundarySolution(j,k,iEq,SOUTH,iEl)  = &
                        myDGSEM % state % boundarySolution(j,k,iEq,SOUTH,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,0)*myDGSEM % state % solution(j,i,k,iEq,iEl)

                     myDGSEM % state % boundarySolution(j,k,iEq,NORTH,iEl)   = &
                        myDGSEM % state % boundarySolution(j,k,iEq,NORTH,iEl)  +  &
                        myDGSEM % dgStorage % bMat(i,1)*myDGSEM % state % solution(j,i,k,iEq,iEl)

                     myDGSEM % state % boundarySolution(j,k,iEq,BOTTOM,iEl)  = &
                        myDGSEM % state % boundarySolution(j,k,iEq,BOTTOM,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,0)*myDGSEM % state % solution(j,k,i,iEq,iEl)

                     myDGSEM % state % boundarySolution(j,k,iEq,TOP,iEl)  = &
                        myDGSEM % state % boundarySolution(j,k,iEq,TOP,iEl)  + &
                        myDGSEM % dgStorage % bMat(i,1)*myDGSEM % state % solution(j,k,i,iEq,iEl)

                  ENDDO
               
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !$OMP ENDDO

#endif

 END SUBROUTINE CalculateBoundarySolution_Fluid
!
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!                                                                                                 !
!  This section of code contains routines for applying boundary conditions along physical         !
!  boundaries.                                                                                    !
!                                                                                                 !
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!
 SUBROUTINE UpdateExternalState_Fluid( myDGSEM, tn, myRank ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)      :: tn
   INTEGER, INTENT(in)         :: myRank
#ifdef HAVE_CUDA
   ! Local
   TYPE(dim3) :: grid, tBlock
  
     tBlock = dim3(4*(ceiling( REAL(myDGSEM % N+1)/4 ) ), &
                   4*(ceiling( REAL(myDGSEM % N+1)/4 ) ) , &
                   1 )
     grid = dim3(myDGSEM % nBoundaryFaces,1,1) 
     
     CALL UpdateExternalState_CUDAKernel<<<grid, tBlock>>>( myDGSEM % extComm % boundaryIDs_dev, &       ! I
                                                            myDGSEM % mesh % faces_dev % elementIDs, &   ! I
                                                            myDGSEM % mesh % faces_dev % elementSides, & ! I
                                                            myDGSEM % extComm % extProcIDs_dev, &           ! I
                                                            myDGSEM % externalState_dev, &               ! O
                                                            myDGSEM % state % boundarySolution_dev, &    ! I
                                                            myDGSEM % prescribedState_dev, &             ! I
                                                            myDGSEM % mesh % geom_dev % nHat_dev )           ! I
#else
   ! Local
   INTEGER    :: iEl, iFace, bFaceID, i, j, k, iEq
   INTEGER    :: iFace2, p2
   INTEGER    :: e1, e2, s1, s2
   REAL(prec) :: norm, un, ut, us, speed
   REAL(prec) :: nx, ny, nz
   REAL(prec) :: sx, sy, sz
   REAL(prec) :: tx, ty, tz
   
      !$OMP DO
      DO iFace = 1, myDGSEM % nBoundaryFaces

         iFace2 = myDGSEM % extComm % boundaryIDs( iFace ) ! Obtain the process-local face id for this boundary-face id
         e1     = myDGSEM % mesh % Faces(iFace2) % elementIDs(1)
         s1     = myDGSEM % mesh % Faces(iFace2) % elementSides(1)
         e2     = myDGSEM % mesh % Faces(iFace2) % elementIDs(2)
         p2     = myDGSEM % extComm % extProcIDs( iFace )
            DO j = 0, myDGSEM % N 
               DO i = 0, myDGSEM % N
               
                     IF( e2 == PRESCRIBED .AND. p2 == myRank )THEN
                     !   PRINT*,'PRESCRIBED',myDGSEM % prescribedState(i,j,1,iFace)
                        DO iEq = 1, myDGSEM % nEq
                           myDGSEM % externalState(i,j,iEq,iFace) = myDGSEM % prescribedState(i,j,iEq,iFace)
                        ENDDO
                     ELSEIF( e2 == RADIATION .AND. p2 == myRank )THEN
                        
                        ! momentum
                        ! rho*u
                        myDGSEM % externalState(i,j,1,iFace) = 0.0_prec
                        ! rho*v
                        myDGSEM % externalState(i,j,2,iFace) = 0.0_prec
                        ! rho*w
                        myDGSEM % externalState(i,j,3,iFace) = 0.0_prec
                        ! Density is set to the static density field
                        myDGSEM % externalState(i,j,4,iFace) = 0.0_prec
                        ! Potential Temperature anomaly (multiplied by density) is set to its static state
                        myDGSEM % externalState(i,j,5,iFace) = 0.0_prec
                        ! Pressure anomaly is set to zero                               
                        myDGSEM % externalState(i,j,6,iFace) = 0.0_prec
                        
                     ELSEIF( e2 == NO_NORMAL_FLOW .AND. p2 == myRank )THEN
                              
                        ! normal
                        nx = myDGSEM % mesh % geom(e1) % nHat(1,i,j,s1)
                        ny = myDGSEM % mesh % geom(e1) % nHat(2,i,j,s1)
                        nz = myDGSEM % mesh % geom(e1) % nHat(3,i,j,s1)
                        norm = sqrt( nx*nx + ny*ny + nz*nz )
                        nx = nx/norm
                        ny = ny/norm
                        nz = nz/norm
      
                        ! tangent
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
                        
                        un = myDGSEM % state % boundarySolution(i,j,1,s1,e1)*nx + &
                             myDGSEM % state % boundarySolution(i,j,2,s1,e1)*ny + &
                             myDGSEM % state % boundarySolution(i,j,3,s1,e1)*nz
                        us = myDGSEM % state % boundarySolution(i,j,1,s1,e1)*sx    + &
                             myDGSEM % state % boundarySolution(i,j,2,s1,e1)*sy    + &
                             myDGSEM % state % boundarySolution(i,j,3,s1,e1)*sz
                        ut = myDGSEM % state % boundarySolution(i,j,1,s1,e1)*tx  + &
                             myDGSEM % state % boundarySolution(i,j,2,s1,e1)*ty  + &
                             myDGSEM % state % boundarySolution(i,j,3,s1,e1)*tz
                             
                        
                        
                        myDGSEM % externalState(i,j,1,iFace) = -nx*un + us*sx + ut*tx ! u
                        myDGSEM % externalState(i,j,2,iFace) = -ny*un + us*sy + ut*ty ! v
                        myDGSEM % externalState(i,j,3,iFace) = -nz*un + us*sz + ut*tz ! w
                        myDGSEM % externalState(i,j,4,iFace) =  myDGSEM % state % boundarySolution(i,j,4,s1,e1) ! rho
                        myDGSEM % externalState(i,j,5,iFace) =  myDGSEM % state % boundarySolution(i,j,5,s1,e1) ! potential temperature
                        myDGSEM % externalState(i,j,6,iFace) =  myDGSEM % state % boundarySolution(i,j,6,s1,e1) ! P
                        
                     ELSEIF( e2 == DRAG_SLIP.AND. p2 == myRank )THEN
                              
                        ! normal
                        nx = myDGSEM % mesh % geom(e1) % nHat(1,i,j,s1)
                        ny = myDGSEM % mesh % geom(e1) % nHat(2,i,j,s1)
                        nz = myDGSEM % mesh % geom(e1) % nHat(3,i,j,s1)
                        norm = sqrt( nx*nx + ny*ny + nz*nz )
                        nx = nx/norm
                        ny = ny/norm
                        nz = nz/norm
      
                        ! tangent
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
                        
                        speed = sqrt( myDGSEM % state % boundarySolution(i,j,1,s1,e1)**2 + &
                                      myDGSEM % state % boundarySolution(i,j,2,s1,e1)**2 + &
                                      myDGSEM % state % boundarySolution(i,j,3,s1,e1)**2 )/&
                                      myDGSEM % state % boundarySolution(i,j,4,s1,e1)
                                      
                        un = myDGSEM % state % boundarySolution(i,j,1,s1,e1)*nx + &
                             myDGSEM % state % boundarySolution(i,j,2,s1,e1)*ny + &
                             myDGSEM % state % boundarySolution(i,j,3,s1,e1)*nz
                             
                        us = ( myDGSEM % state % boundarySolution(i,j,1,s1,e1)*sx + &
                               myDGSEM % state % boundarySolution(i,j,2,s1,e1)*sy + &
                               myDGSEM % state % boundarySolution(i,j,3,s1,e1)*sz )*&
                               (1.0_prec-myDGSEM % params % Cd*myDGSEM % params % dragScale*speed)
                             
                        ut = ( myDGSEM % state % boundarySolution(i,j,1,s1,e1)*tx + &
                               myDGSEM % state % boundarySolution(i,j,2,s1,e1)*ty + &
                               myDGSEM % state % boundarySolution(i,j,3,s1,e1)*tz )*&
                               (1.0_prec-myDGSEM % params % Cd*myDGSEM % params % dragScale*speed)
                        
                        myDGSEM % externalState(i,j,1,iFace) = -nx*un + us*sx + ut*tx ! u
                        myDGSEM % externalState(i,j,2,iFace) = -ny*un + us*sy + ut*ty ! v
                        myDGSEM % externalState(i,j,3,iFace) = -nz*un + us*sz + ut*tz ! w
                        myDGSEM % externalState(i,j,4,iFace) =  myDGSEM % state % boundarySolution(i,j,4,s1,e1) ! rho
                        myDGSEM % externalState(i,j,5,iFace) =  myDGSEM % state % boundarySolution(i,j,5,s1,e1) ! potential temperature
                        myDGSEM % externalState(i,j,6,iFace) =  myDGSEM % state % boundarySolution(i,j,6,s1,e1) ! P
                        
                     
                     ENDIF
               ENDDO
            ENDDO
      ENDDO
      !$OMP ENDDO 
      
#endif

 END SUBROUTINE UpdateExternalState_Fluid
! 
#ifdef HAVE_MPI
 SUBROUTINE MPI_StateExchange_Fluid( myDGSEM, myRank ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)         :: myRank
   ! Local
   INTEGER    :: iFace, bID
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, iNeighbor, jUnpack

#ifdef HAVE_CUDA
      ! For now, we update the CPU with the device boundary solution data before message passing
      myDGSEM % state % boundarySolution = myDGSEM % state % boundarySolution_dev
#endif
      DO iNeighbor = 1, myDGSEM % nNeighbors
         myDGSEM % mpiPackets(iNeighbor) % bufferCounter = 0
      ENDDO
 
      DO bID = 1, myDGSEM % extComm % nBoundaries
      
         iFace     = myDGSEM % extComm % boundaryIDs( bID )
         p2        = myDGSEM % extComm % extProcIDs(bID)
         
         ! In the event that the external process ID (p2) is identical to the current rank (p1),
         ! then this boundary edge involves a physical boundary condition and does not require a 
         ! message exchange
         IF( p2 /= myRank )THEN 

            e1        = myDGSEM % mesh % Faces(iFace) % elementIDs(1)
            s1        = myDGSEM % mesh % Faces(iFace) % elementSides(1)
            iNeighbor = myDGSEM % rankTable(p2)
         
            myDGSEM % mpiPackets(iNeighbor) % bufferCounter = myDGSEM % mpiPackets(iNeighbor) % bufferCounter + 1
            myDGSEM % mpiPackets(iNeighbor) % sendStateBuffer(:,:,:,myDGSEM % mpiPackets(iNeighbor) % bufferCounter ) =&
               myDGSEM % state % boundarySolution(:,:,:,s1,e1) 

         ENDIF
      ENDDO

      DO iNeighbor = 1, myDGSEM % nNeighbors 
            
            ! In the event that the external process ID (p2) is not equal to the current rank (p1),
            ! then this boundary requires a message exchange between ranks p1 and p2.
         
            ! We need to send the internal state along the shared edge to process p2.
            ! A unique "tag" for the message is the global edge ID
            CALL MPI_IRECV( myDGSEM % mpiPackets(iNeighbor) % recvStateBuffer, & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*nEq*myDGSEM % mpiPackets(iNeighbor) % bufferSize, &                
                           MPI_PREC,   &                      
                           myDGSEM % mpiPackets(iNeighbor) % neighborRank, 0,  &                        
                           MPI_COMM_WORLD,   &                
                           stateReqHandle((iNeighbor-1)*2+1), iError )           

            CALL MPI_ISEND( myDGSEM % mpiPackets(iNeighbor) % sendStateBuffer, & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*nEq*myDGSEM % mpiPackets(iNeighbor) % bufferSize, &       
                           MPI_PREC, &      
                           myDGSEM % mpiPackets(iNeighbor) % neighborRank, 0, &       
                           MPI_COMM_WORLD, &
                           stateReqHandle(iNeighbor*2), iError)  
                           
      ENDDO


 END SUBROUTINE MPI_StateExchange_Fluid
!
 SUBROUTINE FinalizeMPI_StateExchange_Fluid( myDGSEM, myRank ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)         :: myRank
   ! Local
   INTEGER    :: iFace, bID
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, iNeighbor, jUnpack


      CALL MPI_WaitAll(myDGSEM % nNeighbors*2,stateReqHandle,stateStats,iError)


      DO bID = 1, myDGSEM % extComm % nBoundaries

         iFace     = myDGSEM % extComm % boundaryIDs( bID )
         p2        = myDGSEM % extComm % extProcIDs(bID)
         
         ! In the event that the external process ID (p2) is identical to the current rank (p1),
         ! then this boundary edge involves a physical boundary condition and does not require a 
         ! message exchange
         IF( p2 /= myRank )THEN 
      
            iNeighbor = myDGSEM % rankTable(p2)
            jUnpack   = myDGSEM % extComm % unpackMap(bID)
            IF( jUnpack == 0 )THEN
              PRINT*, 'Something catastrophic happenend!'
              CALL myDGSEM % Trash( )
              STOP
            ENDIF
            myDGSEM % externalState(:,:,:,bID) = myDGSEM % mpiPackets(iNeighbor) % recvStateBuffer(:,:,:,jUnpack)

         ENDIF

      ENDDO
#ifdef HAVE_CUDA
      ! For now, we update the device with data received from message passing
      myDGSEM % externalState_dev = myDGSEM % externalState
#endif


 END SUBROUTINE FinalizeMPI_StateExchange_Fluid
!
 SUBROUTINE MPI_StressExchange_Fluid( myDGSEM, myRank ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)         :: myRank
   ! Local
   INTEGER    :: iFace, bID
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, iNeighbor, jUnpack

#ifdef HAVE_CUDA
      ! For now, we update the CPU with the device boundary solution data before message passing
      myDGSEM % stressTensor % boundarySolution = myDGSEM % stressTensor % boundarySolution_dev
#endif
      DO iNeighbor = 1, myDGSEM % nNeighbors
         myDGSEM % mpiPackets(iNeighbor) % bufferCounter = 0
      ENDDO
 
      DO bID = 1, myDGSEM % extComm % nBoundaries
      
         iFace     = myDGSEM % extComm % boundaryIDs( bID )
         p2        = myDGSEM % extComm % extProcIDs(bID)
         
         ! In the event that the external process ID (p2) is identical to the current rank (p1),
         ! then this boundary edge involves a physical boundary condition and does not require a 
         ! message exchange
         IF( p2 /= myRank )THEN 

            e1        = myDGSEM % mesh % Faces(iFace) % elementIDs(1)
            s1        = myDGSEM % mesh % Faces(iFace) % elementSides(1)
            iNeighbor = myDGSEM % rankTable(p2)
         
            myDGSEM % mpiPackets(iNeighbor) % bufferCounter = myDGSEM % mpiPackets(iNeighbor) % bufferCounter + 1
            myDGSEM % mpiPackets(iNeighbor) % sendStressBuffer(:,:,:,myDGSEM % mpiPackets(iNeighbor) % bufferCounter ) =&
               myDGSEM % stressTensor % boundarySolution(:,:,:,s1,e1) 

         ENDIF
      ENDDO

      DO iNeighbor = 1, myDGSEM % nNeighbors 
            
            ! In the event that the external process ID (p2) is not equal to the current rank (p1),
            ! then this boundary requires a message exchange between ranks p1 and p2.
         
            ! We need to send the internal state along the shared edge to process p2.
            ! A unique "tag" for the message is the global edge ID
            CALL MPI_IRECV( myDGSEM % mpiPackets(iNeighbor) % recvStressBuffer, & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*3*myDGSEM % mpiPackets(iNeighbor) % bufferSize, &                
                           MPI_PREC,   &                      
                           myDGSEM % mpiPackets(iNeighbor) % neighborRank, 0,  &                        
                           MPI_COMM_WORLD,   &                
                           stressReqHandle((iNeighbor-1)*2+1), iError )           

            CALL MPI_ISEND( myDGSEM % mpiPackets(iNeighbor) % sendStressBuffer, & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*3*myDGSEM % mpiPackets(iNeighbor) % bufferSize, &       
                           MPI_PREC, &      
                           myDGSEM % mpiPackets(iNeighbor) % neighborRank, 0, &       
                           MPI_COMM_WORLD, &
                           stressReqHandle(iNeighbor*2), iError)  
                           
      ENDDO

 END SUBROUTINE MPI_StressExchange_Fluid
!
 SUBROUTINE FinalizeMPI_StressExchange_Fluid( myDGSEM, myRank ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)         :: myRank
   ! Local
   INTEGER    :: iFace, bID
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, iNeighbor, jUnpack


      CALL MPI_WaitAll(myDGSEM % nNeighbors*2,stressReqHandle,stressStats,iError)


      DO bID = 1, myDGSEM % extComm % nBoundaries

         iFace     = myDGSEM % extComm % boundaryIDs( bID )
         p2        = myDGSEM % extComm % extProcIDs(bID)
         
         ! In the event that the external process ID (p2) is identical to the current rank (p1),
         ! then this boundary edge involves a physical boundary condition and does not require a 
         ! message exchange
         IF( p2 /= myRank )THEN 
      
            iNeighbor = myDGSEM % rankTable(p2)
            jUnpack   = myDGSEM % extComm % unpackMap(bID)
            IF( jUnpack == 0 )THEN
              PRINT*, 'Something catastrophic happenend!'
              CALL myDGSEM % Trash( )
              STOP
            ENDIF
            myDGSEM % externalStress(:,:,:,bID) = myDGSEM % mpiPackets(iNeighbor) % recvStressBuffer(:,:,:,jUnpack)

         ENDIF

      ENDDO
#ifdef HAVE_CUDA
      ! For now, we update the device with data received from message passing
      myDGSEM % externalStress_dev = myDGSEM % externalStress
#endif

 END SUBROUTINE FinalizeMPI_StressExchange_Fluid
!
 SUBROUTINE MPI_SGSExchange_Fluid( myDGSEM, myRank ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)         :: myRank
   ! Local
   INTEGER    :: iFace, bID
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, iNeighbor, jUnpack
   INTEGER    :: reqHandle(1:myDGSEM % nNeighbors*2)
   INTEGER    :: theStats(MPI_STATUS_SIZE,1:myDGSEM % nNeighbors*2)

#ifdef HAVE_CUDA
      ! For now, we update the CPU with the device boundary solution data before message passing
      myDGSEM % sgsCoeffs % boundarySolution = myDGSEM % sgsCoeffs % boundarySolution_dev
#endif
      DO iNeighbor = 1, myDGSEM % nNeighbors
         myDGSEM % mpiPackets(iNeighbor) % bufferCounter = 0
      ENDDO
 
      DO bID = 1, myDGSEM % extComm % nBoundaries
      
         iFace     = myDGSEM % extComm % boundaryIDs( bID )
         p2        = myDGSEM % extComm % extProcIDs(bID)
         
         ! In the event that the external process ID (p2) is identical to the current rank (p1),
         ! then this boundary edge involves a physical boundary condition and does not require a 
         ! message exchange
         IF( p2 /= myRank )THEN 

            e1        = myDGSEM % mesh % Faces(iFace) % elementIDs(1)
            s1        = myDGSEM % mesh % Faces(iFace) % elementSides(1)
            iNeighbor = myDGSEM % rankTable(p2)
         
            myDGSEM % mpiPackets(iNeighbor) % bufferCounter = myDGSEM % mpiPackets(iNeighbor) % bufferCounter + 1
            myDGSEM % mpiPackets(iNeighbor) % sendSGSBuffer(:,:,:,myDGSEM % mpiPackets(iNeighbor) % bufferCounter ) =&
               myDGSEM % sgsCoeffs % boundarySolution(:,:,:,s1,e1) 

         ENDIF
      ENDDO

      DO iNeighbor = 1, myDGSEM % nNeighbors 
            
            ! In the event that the external process ID (p2) is not equal to the current rank (p1),
            ! then this boundary requires a message exchange between ranks p1 and p2.
         
            ! We need to send the internal state along the shared edge to process p2.
            ! A unique "tag" for the message is the global edge ID
            CALL MPI_IRECV( myDGSEM % mpiPackets(iNeighbor) % recvSGSBuffer, & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*myDGSEM % mpiPackets(iNeighbor) % bufferSize, &                
                           MPI_PREC,   &                      
                           myDGSEM % mpiPackets(iNeighbor) % neighborRank, 0,  &                        
                           MPI_COMM_WORLD,   &                
                           SGSReqHandle((iNeighbor-1)*2+1), iError )           

            CALL MPI_ISEND( myDGSEM % mpiPackets(iNeighbor) % sendSGSBuffer, & 
                           (myDGSEM % N+1)*(myDGSEM % N+1)*(nEq-1)*myDGSEM % mpiPackets(iNeighbor) % bufferSize, &       
                           MPI_PREC, &      
                           myDGSEM % mpiPackets(iNeighbor) % neighborRank, 0, &       
                           MPI_COMM_WORLD, &
                           SGSReqHandle(iNeighbor*2), iError)  
                           
      ENDDO

   
 END SUBROUTINE MPI_SGSExchange_Fluid
!
 SUBROUTINE FinalizeMPI_SGSExchange_Fluid( myDGSEM, myRank ) 

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)         :: myRank
   ! Local
   INTEGER    :: iFace, bID
   INTEGER    :: tag, ierror
   INTEGER    :: e1, e2, s1, p2, iNeighbor, jUnpack


      CALL MPI_WaitAll(myDGSEM % nNeighbors*2,SGSReqHandle,SGSStats,iError)


      DO bID = 1, myDGSEM % extComm % nBoundaries

         iFace     = myDGSEM % extComm % boundaryIDs( bID )
         p2        = myDGSEM % extComm % extProcIDs(bID)
         
         ! In the event that the external process ID (p2) is identical to the current rank (p1),
         ! then this boundary edge involves a physical boundary condition and does not require a 
         ! message exchange
         IF( p2 /= myRank )THEN 
      
            iNeighbor = myDGSEM % rankTable(p2)
            jUnpack   = myDGSEM % extComm % unpackMap(bID)
            IF( jUnpack == 0 )THEN
              PRINT*, 'Something catastrophic happenend!'
              CALL myDGSEM % Trash( )
              STOP
            ENDIF
            myDGSEM % externalSGS(:,:,:,bID) = myDGSEM % mpiPackets(iNeighbor) % recvSGSBuffer(:,:,:,jUnpack)

         ENDIF

      ENDDO
#ifdef HAVE_CUDA
      ! For now, we update the device with data received from message passing
      myDGSEM % externalSGS_dev = myDGSEM % externalSGS
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
                                                             myDGSEM % dgStorage % qWeight_dev, &
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
                                      myDGSEM % dgStorage % qWeight(j) + &
                                    ( myDGSEM % state % boundaryFlux(j,k,iEq,WEST,iEl)*&
                                      myDGSEM % dgStorage % bMat(i,left-1) + &
                                      myDGSEM % state % boundaryFlux(j,k,iEq,EAST,iEl)*&
                                      myDGSEM % dgStorage % bMat(i,right-1) )/&
                                      myDGSEM % dgStorage % qWeight(i) + &
                                    ( myDGSEM % state % boundaryFlux(i,j,iEq,BOTTOM,iEl)*&
                                      myDGSEM % dgStorage % bMat(k,left-1) + &
                                      myDGSEM % state % boundaryFlux(i,j,iEq,TOP,iEl)*&
                                      myDGSEM % dgStorage % bMat(k,right-1) )/&
                                      myDGSEM % dgStorage % qWeight(k) )/myDGSEM % mesh % geom(iEl) % J(i,j,k)
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
                                                              myDGSEM % dgStorage % qWeight_dev, &
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
                                      myDGSEM % dgStorage % qWeight(j) + &
                                    ( myDGSEM % stressTensor % boundaryFlux(j,k,jEq,WEST,iEl)*&
                                      myDGSEM % dgStorage % bMat(i,left-1) + &
                                      myDGSEM % stressTensor % boundaryFlux(j,k,jEq,EAST,iEl)*&
                                      myDGSEM % dgStorage % bMat(i,right-1) )/&
                                      myDGSEM % dgStorage % qWeight(i) + &
                                    ( myDGSEM % stressTensor % boundaryFlux(i,j,jEq,BOTTOM,iEl)*&
                                      myDGSEM % dgStorage % bMat(k,left-1) + &
                                      myDGSEM % stressTensor % boundaryFlux(i,j,jEq,TOP,iEl)*&
                                      myDGSEM % dgStorage % bMat(k,right-1) )/&
                                      myDGSEM % dgStorage % qWeight(k) )/myDGSEM % mesh % geom(iEl) % J(i,j,k)
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
 SUBROUTINE UpdateExternalStress_Fluid( myDGSEM, tn, myRank ) ! ////////// !

   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)      :: tn
   INTEGER, INTENT(in)         :: myRank
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
   INTEGER    :: iFace2, p2
   INTEGER    :: e1, e2, s1, s2

      !$OMP DO
      DO iFace = 1, myDGSEM % nBoundaryFaces

         iFace2 = myDGSEM % extComm % boundaryIDs( iFace ) ! Obtain the process-local face id for this boundary-face id
         e1     = myDGSEM % mesh % Faces(iFace2) % elementIDs(1)
         s1     = myDGSEM % mesh % Faces(iFace2) % elementSides(1)
         e2     = myDGSEM % mesh % Faces(iFace2) % elementIDs(2)
         p2     = myDGSEM % extComm % extProcIDs( iFace )
         
         IF( p2 == myRank )THEN ! Enforce no boundary flux due to the fluid stress
            DO j = 0, myDGSEM % N 
               DO i = 0, myDGSEM % N
                  DO iEq = 1, (myDGSEM % nEq-1)*3
                        myDGSEM % externalStress(i,j,iEq,iFace) = -myDGSEM % stressTensor % boundarySolution(i,j,iEq,s1,e1)
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
               
                  bID  = ABS(myDGSEM % mesh % faces(iFace) % boundaryID)
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
                                                             myDGSEM % dgStorage % qWeight_dev, &
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
                                      myDGSEM % dgStorage % qWeight(j) + &
                                    ( myDGSEM % stressTensor % boundaryFlux(j,k,iEq,WEST,iEl)*&
                                      myDGSEM % dgStorage % bMat(i,left-1) + &
                                      myDGSEM % stressTensor % boundaryFlux(j,k,iEq,EAST,iEl)*&
                                      myDGSEM % dgStorage % bMat(i,right-1) )/&
                                      myDGSEM % dgStorage % qWeight(i) + &
                                    ( myDGSEM % stressTensor % boundaryFlux(i,j,iEq,BOTTOM,iEl)*&
                                      myDGSEM % dgStorage % bMat(k,left-1) + &
                                      myDGSEM % stressTensor % boundaryFlux(i,j,iEq,TOP,iEl)*&
                                      myDGSEM % dgStorage % bMat(k,right-1) )/&
                                      myDGSEM % dgStorage % qWeight(k) )/myDGSEM % mesh % geom(iEl) % J(i,j,k)
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
 SUBROUTINE WriteTecplot_Fluid( myDGSEM, iter, nPlot, myRank )

  IMPLICIT NONE
 
  CLASS( Fluid ), INTENT(inout) :: myDGsem
  INTEGER, INTENT(in)           :: iter, nPlot
  INTEGER, INTENT(in)           :: myRank
  !LOCAL
  REAL(prec)  :: x(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: y(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: z(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: sol(0:nPlot,0:nPlot,0:nPlot,1:nEq)
  REAL(prec)  :: bsol(0:nPlot,0:nPlot,0:nPlot,1:nEq)

  INTEGER       :: i, j, k, iEl, iEq, fUnit
  CHARACTER(5)  :: zoneID
  CHARACTER(10) :: iterchar
  CHARACTER(4)  :: rankChar
  REAL(prec) :: hCapRatio, c, T
  
  hCapRatio = ( myDGSEM % params % R + myDGSEM % params % Cv ) / myDGSEM % params % Cv

      WRITE(iterChar,'(I10.10)') iter
      WRITE(rankChar,'(I4.4)') myRank
      
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE= 'State.'//rankChar//'.'//iterChar//'.tec', &
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
         WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',nPlot+1,&
                                                     ', J=',nPlot+1,&
                                                     ', K=',nPlot+1,',F=POINT'

         DO k = 0, nPlot
            DO j = 0, nPlot
               DO i = 0, nPlot
                  T =   (bsol(i,j,k,5) + sol(i,j,k,5))/(bsol(i,j,k,4)+sol(i,j,k,4) )
                          
                  ! Sound speed estimate for the external and internal states
                  c = sqrt( myDGSEM % params % R*T*( ( sol(i,j,k,6) + bsol(i,j,k,6) )/myDGSEM % params % P0 )**hCapRatio   )
                  WRITE(fUnit,'(17(F15.7,1x))') x(i,j,k), y(i,j,k), z(i,j,k),&
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
 SUBROUTINE WriteSmoothedTecplot_Fluid( myDGSEM, iter, nPlot, myRank )

  IMPLICIT NONE
 
  CLASS( Fluid ), INTENT(inout) :: myDGsem
  INTEGER, INTENT(in)           :: iter, nPlot
  INTEGER, INTENT(in)           :: myRank
  !LOCAL
  REAL(prec)  :: x(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: y(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: z(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: sol(0:nPlot,0:nPlot,0:nPlot,1:nEq)
  REAL(prec)  :: bsol(0:nPlot,0:nPlot,0:nPlot,1:nEq)

  INTEGER       :: i, j, k, iEl, iEq, fUnit
  CHARACTER(5)  :: zoneID
  CHARACTER(10) :: iterchar
  CHARACTER(4)  :: rankChar

      WRITE(iterChar,'(I10.10)') iter
      WRITE(rankChar,'(I4.4)') myRank
      
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE= 'State-smoothed.'//rankChar//'.'//iterChar//'.tec', &
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
 SUBROUTINE WriteSGSTecplot_Fluid( myDGSEM, iter, myRank )

  IMPLICIT NONE
 
  CLASS( Fluid ), INTENT(inout) :: myDGsem
  INTEGER, INTENT(in)           :: iter
  INTEGER, INTENT(in)           :: myRank
  !LOCAL
  INTEGER       :: i, j, k, iEl, iEq, fUnit
  CHARACTER(5)  :: zoneID
  CHARACTER(10) :: iterchar
  CHARACTER(4)  :: rankChar

      WRITE(iterChar,'(I10.10)') iter
      WRITE(rankChar,'(I4.4)') myRank
      
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE= 'State-SGS.'//rankChar//'.'//iterChar//'.tec', &
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
 SUBROUTINE WriteStressTensorTecplot_Fluid( myDGSEM, iter, nPlot, myRank )

  IMPLICIT NONE
 
  CLASS( Fluid ), INTENT(inout) :: myDGsem
  INTEGER, INTENT(in)           :: iter, nPlot
  INTEGER, INTENT(in)           :: myRank
  !LOCAL
  REAL(prec)  :: x(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: y(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: z(0:nPlot,0:nPlot,0:nPlot)
  REAL(prec)  :: sol(0:nPlot,0:nPlot,0:nPlot,1:15)

  INTEGER       :: i, j, k, iEl, iEq, fUnit
  CHARACTER(5)  :: zoneID
  CHARACTER(10) :: iterchar
  CHARACTER(4)  :: rankChar

      WRITE(iterChar,'(I10.10)') iter
      WRITE(rankChar,'(I4.4)') myRank
      
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE= 'Stress.'//rankChar//'.'//iterChar//'.tec', &
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
 SUBROUTINE WritePickup_Fluid( myDGSEM, iter, myRank )

   IMPLICIT NONE
   CLASS( Fluid ), INTENT(in) :: myDGSEM
   INTEGER, INTENT(in)        :: iter
   INTEGER, INTENT(in)        :: myRank
  ! LOCAL
   CHARACTER(10) :: iterChar
   CHARACTER(4)  :: rankChar
   INTEGER       :: iEl
   INTEGER       :: thisRec, fUnit
   INTEGER       :: iEq, N

      N = myDGSEM % N
     
      WRITE(iterChar,'(I10.10)') iter
      WRITE(rankChar,'(I4.4)') myRank

      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE='State.'//rankChar//'.'//iterChar//'.pickup', &
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
            FILE   = 'State.'//rankChar//'.'//iterChar//'.exs', &
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
 SUBROUTINE ReadPickup_Fluid( myDGSEM, iter, myRank )

   IMPLICIT NONE
   CLASS( Fluid ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)           :: iter
   INTEGER, INTENT(in)           :: myRank
  ! LOCAL
   CHARACTER(10) :: iterChar
   CHARACTER(4)  :: rankChar
   INTEGER       :: iEl, istat
   INTEGER       :: thisRec, fUnit
   INTEGER       :: iEq, N
   LOGICAL       :: itExists
   
      N = myDGSEM % N
      
      WRITE(iterChar,'(I10.10)') iter
      WRITE(rankChar,'(I4.4)') myRank
      INQUIRE( FILE='State.'//rankChar//'.'//iterChar//'.pickup', EXIST = itExists )
     
      IF( itExists )THEN
     
         OPEN( UNIT=NEWUNIT(fUnit), &
               FILE='State.'//rankChar//'.'//iterChar//'.pickup', &
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
     
      INQUIRE( FILE='State.'//rankChar//'.'//iterChar//'.exs', EXIST = itExists )
     
      IF( itExists )THEN
      
         OPEN( UNIT   = NEWUNIT(fUnit), &
                FILE   = 'State.'//rankChar//'.'//iterChar//'.exs', &
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
! SUBROUTINE QuickDiagnostics_Fluid( myDGSEM, tn )
!   IMPLICIT NONE
!   CLASS( Fluid ), INTENT(in) :: myDGSEM 
!   REAL(prec), INTENT(in)     :: tn
!   ! Local
!   INTEGER :: fUnit, i, j, k, iEl
!   LOGICAL :: itExists
!   REAL(prec) :: U(0:myDGSEM % N, 0:myDGSEM % N, 0:myDGSEM % N, 1:myDGSEM % mesh % nElems)
!   REAL(prec) :: V(0:myDGSEM % N, 0:myDGSEM % N, 0:myDGSEM % N, 1:myDGSEM % mesh % nElems)
!   REAL(prec) :: W(0:myDGSEM % N, 0:myDGSEM % N, 0:myDGSEM % N, 1:myDGSEM % mesh % nElems)
!   REAL(prec) :: P(0:myDGSEM % N, 0:myDGSEM % N, 0:myDGSEM % N, 1:myDGSEM % mesh % nElems)
!   REAL(prec) :: c(0:myDGSEM % N, 0:myDGSEM % N, 0:myDGSEM % N, 1:myDGSEM % mesh % nElems)
!   REAL(prec) :: theta(0:myDGSEM % N, 0:myDGSEM % N, 0:myDGSEM % N, 1:myDGSEM % mesh % nElems)
!   REAL(prec) :: rho(0:myDGSEM % N, 0:myDGSEM % N, 0:myDGSEM % N, 1:myDGSEM % mesh % nElems)
!   REAL(prec) :: hCapRatio, rC
   
!      INQUIRE( FILE = "QuickDiagnostics.txt", EXIST = itExists )
!      IF( itExists )THEN
!          OPEN( UNIT = NewUnit(fUnit), &
!                FILE = "QuickDiagnostics.txt", &
!                STATUS ="OLD", &
!                POSITION = "APPEND", &
!                ACTION = "WRITE" )
!      ELSE
!         OPEN( UNIT = NewUnit(fUnit), &
!               FILE = "QuickDiagnostics.txt", &
!               STATUS = "NEW", &
!               ACTION = "WRITE" )
               

!      ENDIF
!      ! Sound speed estimate for the external and internal states
!      hCapRatio = ( myDGSEM % params % R + myDGSEM % params % Cv ) / myDGSEM % params % Cv
!      rC        =   myDGSEM % params % R / ( myDGSEM % params % R + myDGSEM % params % Cv )
                  
      
!      DO iEl = 1, myDGSEM % mesh % nElems
!         DO k = 0, myDGSEM % N
!            DO j = 0, myDGSEM % N
!               DO i = 0, myDGSEM % N
               
!                  U(i,j,k,iEl) = myDGSEM % state(iEl) % solution(i,j,k,1)/&
!                               ( myDGSEM % state(iEl) % solution(i,j,k,4) )
                                 
!                  V(i,j,k,iEl) = myDGSEM % state(iEl) % solution(i,j,k,2)/&
!                               ( myDGSEM % state(iEl) % solution(i,j,k,4) )
                                 
!                  W(i,j,k,iEl) = myDGSEM % state(iEl) % solution(i,j,k,3)/&
!                               ( myDGSEM % state(iEl) % solution(i,j,k,4) )
                                 
!                  P(i,j,k,iEl) = myDGSEM % state(iEl) % solution(i,j,k,6)
                  
!                  rho(i,j,k,iEl) = myDGSEM % state(iEl) % solution(i,j,k,4)
                  
!                  theta(i,j,k,iEl) = myDGSEM % state(iEl) % solution(i,j,k,5)/ &
!                                    (myDGSEM % state(iEl) % solution(i,j,k,4))
                       
!                  c(i,j,k,iEl)  = sqrt(  myDGSEM % params % R*theta(i,j,k,iEl)*( &
!                                        ( P(i,j,k,iEl) + myDGSEM % static(iEl) % solution(i,j,k,3) )/ &
!                                          myDGSEM % params % P0 )**rC  ) 

!               ENDDO
!            ENDDO
!         ENDDO
!      ENDDO
      
!      WRITE(fUnit, *) "==============================================================="                
!      WRITE(fUnit,*) "time : ", tn
!      WRITE(fUnit,*) "Max/Min (U)          : ", MAXVAL(U), MINVAL(U)
!      WRITE(fUnit,*) "Max/Min (V)          : ", MAXVAL(V), MINVAL(V)
!      WRITE(fUnit,*) "Max/Min (W)          : ", MAXVAL(W), MINVAL(W)
!      WRITE(fUnit,*) "Max/Min (rho)        : ", MAXVAL(rho), MINVAL(rho)
!      WRITE(fUnit,*) "Max/Min (P)          : ", MAXVAL(P), MINVAL(P)
!      WRITE(fUnit,*) "Max/Min (c)          : ", MAXVAL(c), MINVAL(c)
!      WRITE(fUnit,*) "Max/Min (Pot. Temp.) : ", MAXVAL(theta), MINVAL(theta)
      
!      CLOSE(fUnit)
   
! END SUBROUTINE QuickDiagnostics_Fluid
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
      ! ////////////////////////////////////////////////////////////////////////// !
      i   = threadIdx % x-1
      j   = threadIdx % y-1

      IF( iFace <= nBoundaryFaces_dev )THEN
      
         iFace2 = boundaryIDs( iFace ) ! Obtain the process-local face id for this boundary-face id
         e1     = elementIDs(1,iFace2)
         s1     = elementSides(1,iFace2)
         e2     = elementIDs(2,iFace2)
         p2     = procIDs( iFace )
         
         IF( i <= polydeg_dev .AND. j <= polydeg_dev )THEN
         
            IF( e2 == PRESCRIBED .AND. p2 == myRank_dev )THEN
               
               DO iEq = 1, nEq_dev
                  externalState(i,j,iEq,iFace) = prescribedState(i,j,iEq,iFace)
               ENDDO
                  
            ELSEIF( e2 == RADIATION .AND. p2 == myRank_dev )THEN
                        
               DO iEq = 1, nEq_dev
                  externalState(i,j,iEq,iFace) = 0.0_prec
               ENDDO
                 
            ELSEIF( e2 == NO_NORMAL_FLOW .AND. p2 == myRank_dev )THEN
                             
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
                        
            ELSEIF( e2 == DRAG_SLIP .AND. p2 == myRank_dev )THEN
                             
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
                          (boundarySolution(i,j,4,s2,e2)+boundarySolution_static(i,j,4,s1,e1) )  
                                   
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
                                                                Ja, Jac, bMat, qWeight, dMatP, tendency )

   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(in)  :: solution(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: static(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundaryFlux(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: drag(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: Ja(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:3,1:3,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: Jac(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: bMat(0:polydeg_dev,0:1)
   REAL(prec), DEVICE, INTENT(in)  :: qWeight(0:polydeg_dev)
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
                qWeight(j) + &
                ( boundaryFlux(j,k,iEq,4,iEl)*bMat(i,0) + &
                  boundaryFlux(j,k,iEq,2,iEl)*bMat(i,1) )/&
                qWeight(i) + &
                ( boundaryFlux(i,j,iEq,5,iEl)*bMat(k,0) + &
                  boundaryFlux(i,j,iEq,6,iEl)*bMat(k,1) )/&
                qWeight(k) )/Jac(i,j,k,iEl)
                      
             
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
 ATTRIBUTES(Global) SUBROUTINE CalculateStressTensor_CUDAKernel( solution, static, dMatP, bmat, qWeight, Ja, Jac, stressFlux, stressTensor ) 
 
   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(in)  :: solution(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: static(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: dMatP(0:polydeg_dev,0:polydeg_dev)
   REAL(prec), DEVICE, INTENT(in)  :: bmat(0:polydeg_dev,0:1)
   REAL(prec), DEVICE, INTENT(in)  :: qWeight(0:polydeg_dev)
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
                qWeight(j) + &
                ( stressFlux(j,k,idir + (iEq-1)*3,4,iEl)*bMat(i,0) + &
                  stressFlux(j,k,idir + (iEq-1)*3,2,iEl)*bMat(i,1) )/&
                qWeight(i) + &
                ( stressFlux(i,j,idir + (iEq-1)*3,5,iEl)*bMat(k,0) + &
                  stressFlux(i,j,idir + (iEq-1)*3,6,iEl)*bMat(k,1) )/&
                qWeight(k) )/Jac(i,j,k,iEl)

                     
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
      
      IF( e2 < 0 )THEN
      
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
      ENDIF

 END SUBROUTINE BoundaryStressFlux_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE StressDivergence_CUDAKernel( stress, stressFlux, state, static, sgsCoeffs, &
                                                             Ja, Jac, bMat, qWeight, dMatP, tendency ) ! ///////////////////// !

   IMPLICIT NONE
   REAL(prec), DEVICE, INTENT(in)  :: stress(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:15,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: stressFlux(0:polydeg_dev,0:polydeg_dev,1:15,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: state(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: static(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: sgsCoeffs(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev-1,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: Ja(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:3,1:3,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: Jac(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: bMat(0:polydeg_dev,0:1)
   REAL(prec), DEVICE, INTENT(in)  :: qWeight(0:polydeg_dev)
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
                qWeight(j) + &
                ( stressFlux(j,k,iEq,4,iEl)*bMat(i,0) + &
                  stressFlux(j,k,iEq,2,iEl)*bMat(i,1) )/&
                qWeight(i) + &
                ( stressFlux(i,j,iEq,5,iEl)*bMat(k,0) + &
                  stressFlux(i,j,iEq,6,iEl)*bMat(k,1) )/&
                qWeight(k) )/Jac(i,j,k,iEl)
                      
             
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
#endif

 END MODULE Fluid_Class


