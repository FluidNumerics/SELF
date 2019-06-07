! Fluid_Class.f90
!
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE Fluid_Class

  USE ModelPrecision
  USE ConstantsDictionary
  USE CommonRoutines
  USE Timing
  USE ModelParameters_Class
  USE Lagrange_Class
  USE NodalDG_Class
  USE NodalDGSolution_3D_Class
  USE HexMesh_Class
  USE MPILayer_Class

  USE BodyForces_Class
  USE Fluid_EquationParser_Class

#ifdef HAVE_CUDA
  USE cudafor
#endif

  USE HDF5

  IMPLICIT NONE

#include "self_macros.h"

  TYPE Fluid

    REAL(prec) :: simulationTime
    REAL(prec) :: volume, mass, KE, PE, heat
  
    TYPE( Fluid_EquationParser ) :: fluidEquations
    TYPE( MultiTimers )          :: timers
    TYPE( ModelParameters )      :: params
    TYPE( HexMesh )              :: mesh
    TYPE( NodalDG )              :: dGStorage
    TYPE( BodyForces )           :: sourceTerms
    TYPE( NodalDGSolution_3D )   :: static
    TYPE( NodalDGSolution_3D )   :: state
    TYPE( NodalDGSolution_3D )   :: stressTensor

    TYPE( MPILayer )             :: mpiStateHandler
    TYPE( MPILayer )             :: mpiStressHandler

  CONTAINS

    PROCEDURE :: Build => Build_Fluid
    PROCEDURE :: Trash => Trash_Fluid

    PROCEDURE :: SetInitialConditions => SetInitialConditions_Fluid
    PROCEDURE :: SetPrescribedState   => SetPrescribedState_Fluid
    PROCEDURE :: CalculateStaticState => CalculateStaticState_Fluid

    PROCEDURE :: Diagnostics           => Diagnostics_Fluid
    PROCEDURE :: WriteDiagnostics      => WriteDiagnostics_Fluid

    PROCEDURE :: IO => IO_Fluid
    PROCEDURE :: Write_to_HDF5
    PROCEDURE :: Read_from_HDF5
    PROCEDURE :: WriteTecplot => WriteTecplot_Fluid

    ! (Soon to be) PRIVATE Routines
    PROCEDURE :: ForwardStepRK3        => ForwardStepRK3_Fluid
!    PROCEDURE :: CrankNicholsonRHS     => CrankNicholsonRHS_Fluid

    PROCEDURE :: GlobalTimeDerivative   => GlobalTimeDerivative_Fluid
    PROCEDURE :: EquationOfState        => EquationOfState_Fluid
    PROCEDURE :: Update_FluidState_BCs
    PROCEDURE :: InternalFace_StateFlux => InternalFace_StateFlux_Fluid
    PROCEDURE :: BoundaryFace_StateFlux => BoundaryFace_StateFlux_Fluid
    PROCEDURE :: MappedTimeDerivative   => MappedTimeDerivative_Fluid

    PROCEDURE :: StressFluxDivergence              => StressFluxDivergence_Fluid
    PROCEDURE :: CalculateStateAtBoundaries        => CalculateStateAtBoundaries_Fluid
    PROCEDURE :: CalculateSolutionGradient         => CalculateSolutionGradient_Fluid
    PROCEDURE :: CalculateNormalStressAtBoundaries => CalculateNormalStressAtBoundaries_Fluid
    PROCEDURE :: CalculateStressFlux               => CalculateStressFlux_Fluid
    PROCEDURE :: Update_FluidStress_BCs
    PROCEDURE :: InternalFace_StressFlux           => InternalFace_StressFlux_Fluid
    PROCEDURE :: BoundaryFace_StressFlux           => BoundaryFace_StressFlux_Fluid

    PROCEDURE ::  Update_FluidStatics_BCs 

  END TYPE Fluid


  INTEGER, PARAMETER, PRIVATE :: nDiagnostics = 5
  INTEGER, PRIVATE            :: diagUnits

  INTEGER, PARAMETER :: nEquations   = 7

#ifdef HAVE_CUDA
 INTEGER, CONSTANT    :: nEq_dev
 INTEGER, CONSTANT    :: nStress_dev
 INTEGER, CONSTANT    :: nEl_dev
 INTEGER, CONSTANT    :: myRank_dev
 INTEGER, CONSTANT    :: nProc_dev
 INTEGER, CONSTANT    :: nNeighbors_dev
 INTEGER, CONSTANT    :: nFaces_dev
 INTEGER, CONSTANT    :: nBoundaryFaces_dev
#endif

CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
  SUBROUTINE Build_Fluid( myDGSEM, equationFile, paramFile, setupSuccess )
  IMPLICIT NONE
#undef __FUNC__
#define __FUNC__ "Build_Fluid"
    CLASS(Fluid), INTENT(inout) :: myDGSEM
    CHARACTER(*), INTENT(in)    :: equationFile
    CHARACTER(*), INTENT(in)    :: paramFile
    LOGICAL, INTENT(inout)      :: setupSuccess
    ! Local
#ifdef HAVE_CUDA
    INTEGER(KIND=cuda_count_KIND) :: freebytes, totalbytes
    INTEGER                       :: iStat, cudaDeviceNumber, nDevices
#endif
    INTEGER :: nB, nEl, nFace, threadID, remainder

    INFO('Start')
    CALL myDGSEM % params % Build( TRIM( paramFile), setupSuccess )
    myDGSEM % simulationTime = myDGSEM % params % startTime

    IF( .NOT. SetupSuccess ) THEN
      INFO( 'Model parameter setup not successful' )
      INFO( 'Halting before building' )
      RETURN
    ENDIF

    CALL myDGSEM % fluidEquations % Build( TRIM( equationFile ) )

#ifdef HAVE_CUDA
    CALL UpdateDeviceDictionary( )
#endif

    CALL myDGSEM % dGStorage % Build( UniformPoints(-1.0_prec, 1.0_prec,myDGSEM % params % nPlot), &
      myDGSEM % params % polyDeg, myDGSEM % params % nPlot, GAUSS )

    CALL myDGSEM % mesh % Load_SELFMesh( myDGSEM % params, myRank, nProc, mpiComm )

    CALL myDGSEM % sourceTerms % Build( myDGSEM % params % polyDeg, nEquations, &
      myDGSEM % mesh % elements % nElements )

    CALL myDGSEM % state % Build( myDGSEM % params % polyDeg, nEquations, &
      myDGSEM % mesh % elements % nElements, myDGSEM % mesh % decomp % nBoundaryFaces )

    CALL myDGSEM % static % Build( myDGSEM % params % polyDeg, nEquations, &
      myDGSEM % mesh % elements % nElements, myDGSEM % mesh % decomp % nBoundaryFaces )


    CALL myDGSEM % stressTensor % Build( myDGSEM % params % polyDeg, &
                                         myDGSEM % state % nEquations-1, &
                                         myDGSEM % mesh % elements % nElements, &
                                         myDGSEM % mesh % decomp % nBoundaryFaces )

#ifdef HAVE_MPI
    CALL myDGSEM % mpiStateHandler % Build( myDGSEM % params % polyDeg, &
                                            myDGSEM % state % nEquations, &
                                            myDGSEM % mesh % decomp % nMPI_Messages )

    CALL myDGSEM % mpiStressHandler % Build( myDGSEM % params % polyDeg, &
                                             myDGSEM % stressTensor % nEquations, &
                                             myDGSEM % mesh % decomp % nMPI_Messages )
#endif

#ifdef HAVE_CUDA
      nEq_dev            = myDGSEM % state % nEquations
      nStress_dev        = myDGSEM % stressTensor % nEquations
      nEl_dev            = myDGSEM % mesh % elements % nElements
      nFaces_dev         = myDGSEM % mesh % faces % nFaces
      nBoundaryFaces_dev = myDGSEM % mesh % decomp % nBoundaryFaces
      myRank_dev         = myRank
      nProc_dev          = nProc
#endif


    INFO('End')

  END SUBROUTINE Build_Fluid
!
  SUBROUTINE Trash_Fluid( myDGSEM )
  IMPLICIT NONE
#undef __FUNC__
#define __FUNC__ "Trash_Fluid"
    CLASS(Fluid), INTENT(inout) :: myDGSEM
    INFO('Start')

    CALL myDGSEM % mesh % Trash(  )
    CALL myDGSEM % dGStorage % Trash( )
    CALL myDGSEM % sourceTerms % Trash( )
    CALL myDGSEM % static % Trash( )
    CALL myDGSEM % state % Trash( )
    CALL myDGSEM % stressTensor % Trash( )

#ifdef HAVE_MPI
    CALL myDGSEM % mpiStateHandler % Trash( )
    CALL myDGSEM % mpiStressHandler % Trash( )
#endif
    INFO('End')

  END SUBROUTINE Trash_Fluid
!
  SUBROUTINE SetInitialConditions_Fluid( myDGSEM )
#undef __FUNC__
#define __FUNC__ "SetInitialConditions"
    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM
    ! Local
    INTEGER    :: i, j, k, iEl, istat
    INTEGER    :: iFace, bID, e1, s1, e2
    REAL(prec) :: x(1:3)
    REAL(prec) :: T, Tbar, u, v, w, rho, rhobar, s, s0

    INFO('Start')

    CALL myDGSEM % CalculateStaticState( ) !! CPU Kernel

    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO k = 0, myDGSEM % params % polyDeg
        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg

            x(1:3) = myDGSEM % mesh % elements % x(i,j,k,1:3,iEl)

            IF( myDGSEM % fluidEquations % calculate_density_from_T )THEN
               
              u  = myDGSEM % fluidEquations % u % evaluate( x ) 
              v  = myDGSEM % fluidEquations % v % evaluate( x ) 
              w  = myDGSEM % fluidEquations % w % evaluate( x ) 
              T  = myDGSEM % fluidEquations % t % evaluate( x ) ! Potential temperature anomaly
              s  = myDGSEM % fluidEquations % tracer % evaluate( x )
              s0 = myDGSEM % fluidEquations % staticTracer % evaluate( x )

              Tbar = myDGSEM % static % solution(i,j,k,5,iEl)/myDGSEM % static % solution(i,j,k,4,iEl)

              myDGSEM % state % solution(i,j,k,4,iEl) = -myDGSEM % static % solution(i,j,k,4,iEl)*T/(Tbar + T)
              ! In the in-situ temperature formulation, the potential temperature is calculated from the
              ! equations file, and we convert to the in-situ temperature here
              ! Since (rho*T)' = 0 in this setting, the pressure anomaly is also zero.
              T  = T*( (myDGSEM % static % solution(i,j,k,7,iEl))/myDGSEM % params % P0 )**( myDGSEM % params % R/( myDGSEM % params % R + myDGSEM % params % Cv ) ) 

              myDGSEM % state % solution(i,j,k,1,iEl) = ( myDGSEM % state % solution(i,j,k,4,iEl) + myDGSEM % static % solution(i,j,k,4,iEl) )*u
              myDGSEM % state % solution(i,j,k,2,iEl) = ( myDGSEM % state % solution(i,j,k,4,iEl) + myDGSEM % static % solution(i,j,k,4,iEl) )*v
              myDGSEM % state % solution(i,j,k,3,iEl) = ( myDGSEM % state % solution(i,j,k,4,iEl) + myDGSEM % static % solution(i,j,k,4,iEl) )*w
              myDGSEM % state % solution(i,j,k,6,iEl) = ( myDGSEM % state % solution(i,j,k,4,iEl) + myDGSEM % static % solution(i,j,k,4,iEl) )*s
              myDGSEM % static % solution(i,j,k,6,iEl) = ( myDGSEM % state % solution(i,j,k,4,iEl) + myDGSEM % static % solution(i,j,k,4,iEl) )*s0

            ELSE

              ! In this case, it is assumed that the temperature passed in is consistent with the CPP flag for the selection of the temperature variable
              u   = myDGSEM % fluidEquations % u % evaluate( x ) 
              v   = myDGSEM % fluidEquations % v % evaluate( x ) 
              w   = myDGSEM % fluidEquations % w % evaluate( x ) 
              rho = myDGSEM % fluidEquations % rho % evaluate( x ) 
              T   = myDGSEM % fluidEquations % t % evaluate( x ) ! Potential temperature anomaly
              s   = myDGSEM % fluidEquations % tracer % evaluate( x )
              s0  = myDGSEM % fluidEquations % staticTracer % evaluate( x )


              myDGSEM % state % solution(i,j,k,4,iEl) = rho
              myDGSEM % state % solution(i,j,k,1,iEl) = ( myDGSEM % state % solution(i,j,k,4,iEl) + myDGSEM % static % solution(i,j,k,4,iEl) )*u
              myDGSEM % state % solution(i,j,k,2,iEl) = ( myDGSEM % state % solution(i,j,k,4,iEl) + myDGSEM % static % solution(i,j,k,4,iEl) )*v
              myDGSEM % state % solution(i,j,k,3,iEl) = ( myDGSEM % state % solution(i,j,k,4,iEl) + myDGSEM % static % solution(i,j,k,4,iEl) )*w
              myDGSEM % state % solution(i,j,k,5,iEl) = ( myDGSEM % state % solution(i,j,k,4,iEl) + myDGSEM % static % solution(i,j,k,4,iEl) )*T
              myDGSEM % state % solution(i,j,k,6,iEl) = ( myDGSEM % state % solution(i,j,k,4,iEl) + myDGSEM % static % solution(i,j,k,4,iEl) )*s
              myDGSEM % static % solution(i,j,k,6,iEl) = ( myDGSEM % state % solution(i,j,k,4,iEl) + myDGSEM % static % solution(i,j,k,4,iEl) )*s0

            ENDIF

            myDGSEM % sourceTerms % drag(i,j,k,iEl) = myDGSEM % fluidEquations % drag % evaluate( x )

          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CALL myDGSEM % SetPrescribedState( ) ! CPU Kernel

#ifdef HAVE_CUDA
    CALL myDGSEM % state % UpdateDevice( )
    CALL myDGSEM % static % UpdateDevice( )
    CALL myDGSEM % sourceTerms % UpdateDevice( )
    istat = cudaDeviceSynchronize( )
#endif

    CALL myDGSEM % EquationOfState( ) ! GPU Kernel (if CUDA)

    CALL myDGSEM % Update_FluidStatics_BCs( ) ! GPU Kernel (if CUDA)

#ifdef HAVE_CUDA
    CALL myDGSEM % state % UpdateHost( )
    istat = cudaDeviceSynchronize( )
#endif
    INFO('End')

  END SUBROUTINE SetInitialConditions_Fluid

  SUBROUTINE SetPrescribedState_Fluid( myDGSEM )
#undef __FUNC__
#define __FUNC__ "SetPrescribedState"
    IMPLICIT NONE

    CLASS( Fluid ), INTENT(inout) :: myDGSEM 
    ! Local
    INTEGER    :: i, j, k, iEl
    INTEGER    :: iFace, bID, e1, s1, e2, p2
    REAL(prec) :: x(1:3)
    REAL(prec) :: T, Tbar, u, v, w, rho, rhobar, s, s0

    DO iFace = 1, myDGSEM % mesh % faces % nFaces

       e1    = myDGSEM % mesh % faces % elementIDs(1,iFace)
       s1    = myDGSEM % mesh % faces % elementSides(1,iFace)
       e2    = myDGSEM % mesh % faces % elementIDs(2,iFace)
       bID   = myDGSEM % mesh % faces % boundaryID(iFace)

       IF( e2 == PRESCRIBED )THEN


         IF( myDGSEM % fluidEquations % calculate_density_from_T )THEN

           DO j = 0, myDGSEM % params % polyDeg
             DO i = 0, myDGSEM % params % polyDeg

               x(1:3) = myDGSEM % mesh % elements % xBound(i,j,1:3,s1,e1)

               u = myDGSEM % fluidEquations % u % evaluate( x ) 
               v = myDGSEM % fluidEquations % v % evaluate( x ) 
               w = myDGSEM % fluidEquations % w % evaluate( x ) 
               T = myDGSEM % fluidEquations % t % evaluate( x ) ! Potential temperature anomaly
               s = myDGSEM % fluidEquations % tracer % evaluate( x )
      
               Tbar = myDGSEM % static % boundarySolution(i,j,5,s1,e1)/myDGSEM % static % boundarySolution(i,j,4,s1,e1)
               ! In the in-situ temperature formulation, the potential temperature is calculated from the
               ! equations file, and we convert to the in-situ temperature here
               ! Since (rho*T)' = 0 in this setting, the pressure anomaly is also zero.
               T  = T*( (myDGSEM % static % boundarySolution(i,j,7,s1,e1))/myDGSEM % params % P0 )**( myDGSEM % params % R/( myDGSEM % params % R + myDGSEM % params % Cv ) ) 
      
               myDGSEM % state % prescribedState(i,j,4,bID) = -myDGSEM % static % boundarySolution(i,j,4,s1,e1)*T/(Tbar + T)
               myDGSEM % state % prescribedState(i,j,1,bID) = ( myDGSEM % state % prescribedState(i,j,4,bID) + myDGSEM % static % boundarySolution(i,j,4,s1,e1) )*u
               myDGSEM % state % prescribedState(i,j,2,bID) = ( myDGSEM % state % prescribedState(i,j,4,bID) + myDGSEM % static % boundarySolution(i,j,4,s1,e1) )*v
               myDGSEM % state % prescribedState(i,j,3,bID) = ( myDGSEM % state % prescribedState(i,j,4,bID) + myDGSEM % static % boundarySolution(i,j,4,s1,e1) )*w
               myDGSEM % state % prescribedState(i,j,6,bID) = ( myDGSEM % state % prescribedState(i,j,4,bID) + myDGSEM % static % boundarySolution(i,j,4,s1,e1) )*s

             ENDDO
           ENDDO
          
         ELSE

           DO j = 0, myDGSEM % params % polyDeg
             DO i = 0, myDGSEM % params % polyDeg

               x(1:3) = myDGSEM % mesh % elements % xBound(i,j,1:3,s1,e1)

               u   = myDGSEM % fluidEquations % u % evaluate( x ) 
               v   = myDGSEM % fluidEquations % v % evaluate( x ) 
               w   = myDGSEM % fluidEquations % w % evaluate( x ) 
               T   = myDGSEM % fluidEquations % t % evaluate( x ) ! Potential temperature anomaly
               rho = myDGSEM % fluidEquations % rho % evaluate( x ) 
               s   = myDGSEM % fluidEquations % tracer % evaluate( x )
  
               Tbar = myDGSEM % static % boundarySolution(i,j,5,s1,e1)/myDGSEM % static % boundarySolution(i,j,4,e1,s1)
  
               myDGSEM % state % prescribedState(i,j,4,bID) = rho
               myDGSEM % state % prescribedState(i,j,1,bID) = ( rho + myDGSEM % static % boundarySolution(i,j,4,s1,e1) )*u
               myDGSEM % state % prescribedState(i,j,2,bID) = ( rho + myDGSEM % static % boundarySolution(i,j,4,s1,e1) )*v
               myDGSEM % state % prescribedState(i,j,3,bID) = ( rho + myDGSEM % static % boundarySolution(i,j,4,s1,e1) )*w
               myDGSEM % state % prescribedState(i,j,5,bID) = ( rho + myDGSEM % static % boundarySolution(i,j,4,s1,e1) )*T
               myDGSEM % state % prescribedState(i,j,6,bID) = ( rho + myDGSEM % static % boundarySolution(i,j,4,s1,e1) )*s

             ENDDO
           ENDDO

         ENDIF

       ENDIF

    ENDDO

#ifdef HAVE_CUDA
    CALL myDGSEM % state % UpdateDevice( )
#endif

  END SUBROUTINE SetPrescribedState_Fluid
!
  SUBROUTINE ForwardStepRK3_Fluid( myDGSEM, nT )
#undef __FUNC__
#define __FUNC__ "ForwardStepRK3"
    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
    INTEGER, INTENT(in)         :: nT
    ! LOCAL
    REAL(prec)                  :: t0
#ifdef HAVE_CUDA
    REAL(prec)           :: t, dt
    REAL(prec), DEVICE   :: t_dev
    REAL(prec), DEVICE :: G3D(0:myDGSEM % params % polyDeg,&
                             0:myDGSEM % params % polyDeg,&
                             0:myDGSEM % params % polyDeg,&
                             1:myDGSEM % state % nEquations,&
                             1:myDGSEM % mesh % elements % nElements)
    INTEGER            :: iT, m, iStat
    TYPE(dim3)         :: fgrid, grid, tBlock
#else
    REAL(prec) :: t, dt, rk3_a_local, rk3_g_local
    REAL(prec), ALLOCATABLE :: G3D(:,:,:,:,:)
    INTEGER    :: m, iEl, iT, i, j, k, iEq

#endif

    INFO('Start')


#ifdef HAVE_CUDA

    tBlock = dim3( myDGSEM % params % polyDeg+1, &
                   myDGSEM % params % polyDeg+1, &
                   myDGSEM % params % polyDeg+1  )

    grid  = dim3(myDGSEM % state % nElements, myDGSEM % state % nEquations-1, 1)
    fgrid = dim3(myDGSEM % state % nElements, myDGSEM % state % nEquations, 1)

    t0 = myDGSEM % simulationTime
    dt = myDGSEM % params % dt

    DO iT = 1, nT

      G3D  = 0.0_prec

      DO m = 1,3

        t = myDGSEM % simulationTime + rk3_b(m)*dt

        CALL myDGSEM % EquationOfState( )
        CALL myDGSEM % GlobalTimeDerivative( t )
        CALL UpdateG3D_CUDAKernel<<<grid,tBlock>>>( G3D, rk3_a_dev(m), rk3_g_dev(m), &
                                                    myDGSEM % state % solution_dev, &
                                                    myDGSEM % state % fluxDivergence_dev, &
                                                    myDGSEM % state % source_dev, &
                                                    myDGSEM % stressTensor % fluxDivergence_dev )


      ENDDO

      myDGSEM % simulationTime = t0 + REAL(iT,prec)*dt

    ENDDO

    ! Determine if we need to take another step with reduced time step to get the solution
    ! at exactly t0+outputFrequency
    IF( .NOT. AlmostEqual( myDGSEM % simulationTime, t0+myDGSEM % params % outputFrequency ) )THEN

      dt = t0+myDGSEM % params % outputFrequency - myDGSEM % simulationTime
      G3D  = 0.0_prec

      DO m = 1,3


        t = myDGSEM % simulationTime  + rk3_b(m)*dt

        CALL myDGSEM % EquationOfState( )


        CALL myDGSEM % GlobalTimeDerivative( t )
        CALL UpdateG3D_CUDAKernel<<<grid,tBlock>>>( G3D, rk3_a_dev(m), rk3_g_dev(m), &
                                                    myDGSEM % state % solution_dev, &
                                                    myDGSEM % state % fluxDivergence_dev, &
                                                    myDGSEM % state % source_dev, &
                                                    myDGSEM % stressTensor % fluxDivergence_dev )


      ENDDO 


      myDGSEM % simulationTime = myDGSEM % simulationTime + dt

    ENDIF

    ! Update the host copy of the solution
    myDGSEM % state % solution = myDGSEM % state % solution_dev

#else

    ALLOCATE( G3D(0:myDGSEM % params % polyDeg,&
                  0:myDGSEM % params % polyDeg,&
                  0:myDGSEM % params % polyDeg,&
                  1:myDGSEM % state % nEquations,&
                  1:myDGSEM % mesh % elements % nElements) )

    t0 = myDGSEM % simulationTime
    dt = myDGSEM % params % dt


    DO iT = 1, nT


      DO iEl = 1, myDGSEM % mesh % elements % nElements
        DO iEq = 1, myDGSEM % state % nEquations-1
          DO k = 0, myDGSEM % params % polyDeg
            DO j = 0, myDGSEM % params % polyDeg
              DO i = 0, myDGSEM % params % polyDeg
                G3D(i,j,k,iEq,iEl) = 0.0_prec
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,3 ! Loop over RK3 steps


        t = myDGSEM % simulationTime + rk3_b(m)*dt

        CALL myDGSEM % EquationOfState( )
        CALL myDGSEM % GlobalTimeDerivative( t )
       !$OMP BARRIER

        DO iEl = 1, myDGSEM % mesh % elements % nElements
          DO iEq = 1, myDGSEM % state % nEquations-1
            DO k = 0, myDGSEM % params % polyDeg
              DO j = 0, myDGSEM % params % polyDeg
                DO i = 0, myDGSEM % params % polyDeg

                  G3D(i,j,k,iEq,iEl) = rk3_a(m)*G3D(i,j,k,iEq,iEl) - ( myDGSEM % state % fluxDivergence(i,j,k,iEq,iEl) -&
                    myDGSEM % stressTensor % fluxDivergence(i,j,k,iEq,iEl) ) + &
                    myDGSEM % state % source(i,j,k,iEq,iEl)


                  myDGSEM % state % solution(i,j,k,iEq,iEl) = myDGSEM % state % solution(i,j,k,iEq,iEl) + &
                    rk3_g(m)*dt*G3D(i,j,k,iEq,iEl)

                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO 

      ENDDO

      myDGSEM % simulationTime = myDGSEM % simulationTime + dt

    ENDDO
    !$OMP BARRIER

    ! Determine IF we need to take another step with reduced time step to get the solution
    ! at exactly t0+outputFrequency
    IF( .NOT. AlmostEqual( myDGSEM % simulationTime, t0+myDGSEM % params % outputFrequency ) )THEN

      dt = t0+myDGSEM % params % outputFrequency - myDGSEM % simulationTime
      DO iEl = 1, myDGSEM % mesh % elements % nElements
        DO iEq = 1, myDGSEM % state % nEquations-1
          DO k = 0, myDGSEM % params % polyDeg
            DO j = 0, myDGSEM % params % polyDeg
              DO i = 0, myDGSEM % params % polyDeg
                G3D(i,j,k,iEq,iEl) = 0.0_prec
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,3

        t = myDGSEM % simulationTime + rk3_b(m)*dt

        CALL myDGSEM % EquationOfState( )
        CALL myDGSEM % GlobalTimeDerivative( t )

        DO iEl = 1, myDGSEM % mesh % elements % nElements
          DO iEq = 1, myDGSEM % state % nEquations-1
            DO k = 0, myDGSEM % params % polyDeg
              DO j = 0, myDGSEM % params % polyDeg
                DO i = 0, myDGSEM % params % polyDeg

                  G3D(i,j,k,iEq,iEl) = rk3_a(m)*G3D(i,j,k,iEq,iEl) - ( myDGSEM % state % fluxDivergence(i,j,k,iEq,iEl) -&
                    myDGSEM % stressTensor % fluxDivergence(i,j,k,iEq,iEl) ) + &
                    myDGSEM % state % source(i,j,k,iEq,iEl)

                  myDGSEM % state % solution(i,j,k,iEq,iEl) = myDGSEM % state % solution(i,j,k,iEq,iEl) + &
                    rk3_g(m)*dt*G3D(i,j,k,iEq,iEl)

                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDDO

      myDGSEM % simulationTime = myDGSEM % simulationTime +  dt

    ENDIF

    CALL myDGSEM % EquationOfState( )


    DEALLOCATE( G3D )

#endif

    INFO('End')

  END SUBROUTINE ForwardStepRK3_Fluid
!
! ============================================================================= !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ============================================================================= !
!
! Crank Nicholson time integrator routines

 SUBROUTINE CrankNicholsonBiCGStab_Fluid( myDGSEM, snk, explicitTendency )
   IMPLICIT NONE
   CLASS(Fluid), INTENT(inout) :: myDGSEM
   REAL(prec), INTENT(in)      :: snk(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)
   REAL(prec), INTENT(in)      :: explicitTendency(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)
   ! Local
   REAL(prec)   :: r(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)
   REAL(prec)   :: ds(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)
   REAL(prec)   :: v(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)
   REAL(prec)   :: p(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)
   REAL(prec)   :: t(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)
   REAL(prec)   :: rho, alpha, omega, beta

      ! Calculate the initial residual
      ! Assumes an initial guess of ds=0
      !r = explicitTendency + myDGSEM % CrankNicholsonRHS( snk )



 END SUBROUTINE CrankNicholsonBiCGStab_Fluid
!
 FUNCTION CrankNicholsonRHS_Fluid( myDGSEM, snk ) RESULT( b )
   ! Given
   IMPLICIT NONE
   CLASS(Fluid) :: myDGSEM
   REAL(prec)   :: snk(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)
   REAL(prec)   :: b(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)


      CALL myDGSEM % GlobalTimeDerivative( myDGSEM % simulationTime )
      b = -( snk - 0.5_prec*myDGSEM % params % dt*myDGSEM % state % tendency )


 END FUNCTION CrankNicholsonRHS_Fluid
!
 FUNCTION CrankNicholsonJacobianAction_Fluid( myDGSEM, s, ds, Fs ) RESULT( Jds )
   IMPLICIT NONE
   CLASS(Fluid) :: myDGSEM
   REAL(prec)   :: s(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)
   REAL(prec)   :: ds(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)
   REAL(prec)   :: Fs(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)
   REAL(prec)   :: Jds(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEquations, 1:myDGSEM % mesh % elements % nElements)


      myDGSEM % state % solution = s + myDGSEM % params % jacobianStepSize*ds

      CALL myDGSEM % GlobalTimeDerivative( myDGSEM % simulationTime )

      ! J*ds = (I - (dt/2)* dF/ds )*ds
      Jds = ds - 0.5_prec*myDGSEM % params % dt*( myDGSEM % state % tendency - Fs )/myDGSEM % params % jacobianStepSize

 END FUNCTION CrankNicholsonJacobianAction_Fluid
!
! ============================================================================= !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ============================================================================= !
!
  SUBROUTINE GlobalTimeDerivative_Fluid( myDGSEM, tn )

#undef __FUNC__
#define __FUNC__ "GlobalTimeDerivative"
!  Here, the solution within each element is interpolated to the faces of each
!  element in order to prepare for computing the external state for enforcing
!  boundary conditions, Riemann Fluxes, and MPI DATA exchanges that need to
!
! When MPI is USEd, the boundary solutions that are stored on faces shared with
! a neighboring rank are passed to that neighboring rank. Additionally, the
! perceived "external state" is received from the neighboring rank. Calling this
! routine is dependent on the result of CalculateBoundarySolutio, but can be
! DOne at the same time as Update_FluidState_BCs; DOing so should hide some
! communication costs.
!  occur.
!
! The boundary solutions are USEd to calculate the external states that, when
! accompanied with a Riemann Solver, enforce boundary conditions. Calling this
! routine is dependent on the result of CalculateBoundarySolution
!
! The inviscid fluxes (advection + pressure) through the faces of each element
! are estimated here using a (linear) Lax-Friedrich's upwind solver. In order to
! call this routine, CalculateBoundarySolution, Update_FluidState_BCs, and
! MPI_StateExchange must have been completed.
!
! IF the SpectralEKE or Laplacian subgridscale models are USEd, a Laplacian-like
! operator is USEd to dIFfUSE momentum and heat. When the Laplacian model is
! USEd a fixed viscosity is specIFied in runtime.params and a Rayleigh number of
! 1 is assumed (viscosity = dIFfusivity). When the SpectralEKE model is USEd,
! the viscosity coefficient is diagnosed using a Smagorinksy closure, where the
! unresolved kinetic energy is diagnosed from a highpass spectral filter.
!
! In the Spectral-EKE model, the under-resolved state is diagnosed from a high
! pass spectral filter (the dIFference between the state and smoothed state).
! Here, we first calculate the smoothed state and store it in the smoothedState
! attribute. This SUBROUTINE call has no dependence to any other within this
! SUBROUTINE.
!
! The high passed solution is USEd to diagnose an isotropic viscosity
! coefficient, similar to a Smagorinksy closure and similar to the closure in
!
!  J. Sauer (2013), "Towards Improved Capability and Confidence in Coupled
!  Atmospheric and Wildland Fire Modeling"
!
! The main dIFference in this work, is in the diagnosis of the SGS Kinetic
! Energy from the high pass filtered solution.
! This routine depends on the results from CalculateSmoothedState.
!
! The viscosity coefficient that is calculated is now interpolated to the faces
! of each element so that the viscous flux can later be computed. This routine
! depends on the result of CalculateSGSCoefficients.
!
! The viscosity coefficients are exchanged with neighboring ranks that share
! COMMON faces. MPI_SGSExchange can be run simulataneously with
! CalculateStressTensor, CalculateBoundaryStress, Update_FluidStress_BCs, and the
! MPI_StressExchange. The viscosity coefficients that are exchanged are not
! needed until StressFlux
!
! Now, the internal solution and the boundaryGradientFlux can be pieced
! together to calculate gradients in the velocity and potential temperature.
! This routine depends on the result of FaceFlux ( state % boundaryGradientFlux )
!
! The solution gradient values are interpolated to the faces of each element and
! projected onto the face normal direction. This
! routine depends on the result of CalculateStressTensor.
!
! Stress tensor values are exchanged with neighboring ranks along shared faces.
! This routine depends on the result of CalculateBoundaryStress, but can be run
! at the same time as Update_FluidStress_BCs.
!
!  Using the solution gradient and the eddy-viscosities/diffusivities, the
!  stress flux is calculated
!
! Now that the stress tensor is available on element faces, boundary conditions
! can be applied by setting the external stress tensor state. This routine
! depends on the result of CalculateBoundaryStress. Note that this routine can
! be run simultaneously with the MPI_StressExchange
!
! Using the boundary and the external state for the stress tensor, the viscous
! fluxes are estimated using a Bassi-Rebay flux that averages neighboring values
! of the stress tensor plus the jump in the solution weighted by a spatial
! wave-number. This routine depends on the result of the Update_FluidStress_BCs
! and the MPI_StressExchange.
!
! With the boundary stress flux and the internal stress tensor values, the
! divergence of the stress tensor can be calculated, giving the viscous tendency
! for the momentum and the potential temperature. This routine depends on the
! result of StressFlux (and the dependencies of StressFlux), but can be DOne
! simultaneously with the MappedTimeDerivative.
!
! Once the inviscid fluxes through the faces are calculated, and the internal
! state is known, the tendency due to the inviscid flux terms and
! nonconservative source terms is calculated here. This routine depends on the
! result of FaceFlux, but can be DOne at the same time as StressDivergence

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
    REAL(prec), INTENT(in)      :: tn
#ifdef HAVE_CUDA
    TYPE(dim3) :: tBlock, grid
    INTEGER    :: istat
#endif

    CALL myDGSEM % CalculateStateAtBoundaries( )
#ifdef HAVE_MPI
    CALL myDGSEM % mpiStateHandler % MPI_Exchange( myDGSEM % state, myDGSEM % mesh )
#endif
    CALL myDGSEM % Update_FluidState_BCs( tn )
    CALL myDGSEM % InternalFace_StateFlux( )

#ifdef HAVE_MPI
    CALL myDGSEM % mpiStateHandler % Finalize_MPI_Exchange( )
#endif

    CALL myDGSEM % BoundaryFace_StateFlux( )

    CALL myDGSEM % CalculateSolutionGradient( )

    CALL myDGSEM % CalculateNormalStressAtBoundaries( )

#ifdef HAVE_MPI
      CALL myDGSEM % mpiStressHandler % MPI_Exchange( myDGSEM % stressTensor, myDGSEM % mesh )
#endif
      CALL myDGSEM % CalculateStressFlux( )

      CALL myDGSEM % Update_FluidStress_BCs( tn )
#ifdef HAVE_MPI
      CALL myDGSEM % mpiStressHandler % Finalize_MPI_Exchange( )
#endif
      CALL myDGSEM % BoundaryFace_StressFlux( )

      CALL myDGSEM % StressFluxDivergence( )


      CALL myDGSEM % MappedTimeDerivative( )

  END SUBROUTINE GlobalTimeDerivative_Fluid
!
  SUBROUTINE StressFluxDivergence_Fluid( myDGSEM )
    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM 
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

      tBlock = dim3(myDGSEM % params % polyDeg+1, &
                    myDGSEM % params % polyDeg+1, &
                    myDGSEM % params % polyDeg+1 )

      grid = dim3(myDGSEM % stressTensor % nEquations, myDGSEM % state % nElements, 1)  

      CALL Stress_Mapped_DG_Divergence_3D_CUDAKernel<<<grid, tBlock>>>( myDGSEM % stressTensor % flux_dev, &
                                                          myDGSEM % stressTensor % boundaryFlux_dev, &
                                                          myDGSEM % stressTensor % fluxDivergence_dev, &
                                                          myDGSEM % dgStorage % boundaryInterpolationMatrix_dev, &
                                                          myDGSEM % dgStorage % dgDerivativeMatrixTranspose_dev, &
                                                          myDGSEM % dgStorage % quadratureWeights_dev, &
                                                          myDGSEM % mesh % elements % J_dev )

#else
    ! Local
    INTEGER    :: ii, i, j, k, iVar, iEl
    REAL(prec) :: df

      DO iEl = 1, myDGSEM % stressTensor % nElements
        DO iVar = 1, myDGSEM % stressTensor % nEquations
          DO k = 0, myDGSEM % stressTensor % N
            DO j = 0, myDGSEM % stresstensor % N
              DO i = 0, myDGSEM % stressTensor % N   
 
                df = 0.0_prec
                DO ii = 0, myDGSEM % dgStorage % N
                  df = df + myDGSEM % dgStorage % dgDerivativeMatrixTranspose(ii,i)*myDGSEM % stressTensor % flux(1,ii,j,k,iVar,iEl) + &
                            myDGSEM % dgStorage % dgDerivativeMatrixTranspose(ii,j)*myDGSEM % stressTensor % flux(2,i,ii,k,iVar,iEl) + &
                            myDGSEM % dgStorage % dgDerivativeMatrixTranspose(ii,k)*myDGSEM % stressTensor % flux(3,i,j,ii,iVar,iEl)
                ENDDO
                 
                myDGSEM % stressTensor % fluxDivergence(i,j,k,iVar,iEl) =  ( df+ ( myDGSEM % stressTensor % boundaryFlux(i,k,iVar,1,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(j,0) + &
                                                myDGSEM % stressTensor % boundaryFlux(i,k,iVar,3,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(j,1) )/&
                                              myDGSEM % dgStorage % quadratureWeights(j) + &
                                              ( myDGSEM % stressTensor % boundaryFlux(j,k,iVar,4,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(i,0) + &
                                                myDGSEM % stressTensor % boundaryFlux(j,k,iVar,2,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(i,1) )/&
                                              myDGSEM % dgStorage % quadratureWeights(i) + &
                                              ( myDGSEM % stressTensor % boundaryFlux(i,j,iVar,5,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(k,0) + &
                                                myDGSEM % stressTensor % boundaryFlux(i,j,iVar,6,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(k,1) )/&
                                              myDGSEM % dgStorage % quadratureWeights(k) )/myDGSEM % mesh % elements % J(i,j,k,iEl) 

              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

#endif


  END SUBROUTINE StressFluxDivergence_Fluid

  SUBROUTINE CalculateStateAtBoundaries_Fluid( myDGSEM ) 
    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM 
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1 , &
                  1 )
    grid = dim3(myDGSEM % state % nEquations, myDGSEM % state % nElements, 1)  

    CALL CalculateStateAtBoundaries_CUDAKernel<<<grid, tBlock>>>( myDGSEM % state % solution_dev, &
                                                                  myDGSEM % state % boundarySolution_dev,  &
                                                                  myDGSEM % dgStorage % boundaryInterpolationMatrix_dev )
#else
    ! Local
    INTEGER :: i, j, k, iVar, iEl
    REAL(prec) :: fb(1:6)


      DO iEl = 1, myDGSEM % state % nElements
        DO iVar = 1, myDGSEM % state % nEquations
          DO j = 0, myDGSEM % state % N
            DO i = 0, myDGSEM % state % N
            
              fb(1:6) = 0.0_prec
              
              DO k = 0, myDGSEM % state % N
               
                fb(1) = fb(1) + myDGSEM % dgStorage % boundaryInterpolationMatrix(k,0)*myDGSEM % state % solution(i,k,j,iVar,iEl) ! South
                fb(2) = fb(2) + myDGSEM % dgStorage % boundaryInterpolationMatrix(k,1)*myDGSEM % state % solution(k,i,j,iVar,iEl) ! East
                fb(3) = fb(3) + myDGSEM % dgStorage % boundaryInterpolationMatrix(k,1)*myDGSEM % state % solution(i,k,j,iVar,iEl) ! North
                fb(4) = fb(4) + myDGSEM % dgStorage % boundaryInterpolationMatrix(k,0)*myDGSEM % state % solution(k,i,j,iVar,iEl) ! West
                fb(5) = fb(5) + myDGSEM % dgStorage % boundaryInterpolationMatrix(k,0)*myDGSEM % state % solution(i,j,k,iVar,iEl) ! Bottom
                fb(6) = fb(6) + myDGSEM % dgStorage % boundaryInterpolationMatrix(k,1)*myDGSEM % state % solution(i,j,k,iVar,iEl) ! Top
               
              ENDDO

              myDGSEM % state % boundarySolution(i,j,iVar,1:6,iEl) = fb(1:6)
              
            ENDDO
          ENDDO
        ENDDO
      ENDDO

#endif

  END SUBROUTINE CalculateStateAtBoundaries_Fluid

! /////////////////////////////////////////////////////////////////////////////////////////////// !
!                                                                                                 !
!  This section of code CONTAINS routines for applying boundary conditions along physical         !
!  boundaries.                                                                                    !
!                                                                                                 !
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!
  SUBROUTINE Update_FluidState_BCs( myDGSEM, tn )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
    REAL(prec), INTENT(in)      :: tn
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1, &
                  1 )
    grid = dim3(myDGSEM % mesh % faces % nFaces,1,1)

    CALL Update_FluidState_BCs_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mesh % decomp % element_to_blockID_dev, &       ! I
                                                             myDGSEM % mesh % faces % boundaryID_dev, &
                                                             myDGSEM % mesh % faces % elementIDs_dev, &   ! I
                                                             myDGSEM % mesh % faces % elementSides_dev, & ! I
                                                             myDGSEM % state % externalState_dev, &               ! O
                                                             myDGSEM % state % boundarySolution_dev, &    ! I
                                                             myDGSEM % state % prescribedState_dev, &             ! I
                                                             myDGSEM % mesh % elements % nHat_dev )
#else
    ! Local
    INTEGER    :: bID, i, j, k, iEq
    INTEGER    :: iFace, p2
    INTEGER    :: e1, e2, s1, s2
    REAL(prec) :: norm, un, ut, us, speed
    REAL(prec) :: nx, ny, nz
    REAL(prec) :: sx, sy, sz
    REAL(prec) :: tx, ty, tz


    DO iFace = 1, myDGSEM % mesh % faces % nFaces

      e1  = myDGSEM % mesh % faces % elementIDs(1,iFace)
      s1  = myDGSEM % mesh % faces % elementSides(1,iFace)
      e2  = myDGSEM % mesh % faces % elementIDs(2,iFace)
      bID = myDGSEM % mesh % faces % boundaryID(iFace)

      IF( e2 == PRESCRIBED )THEN

        DO iEq = 1, myDGSEM % state % nEquations
          DO j = 0, myDGSEM % params % polyDeg
            DO i = 0, myDGSEM % params % polyDeg
              myDGSEM % state % externalState(i,j,iEq,bID) = myDGSEM % state % prescribedState(i,j,iEq,bID)
            ENDDO
          ENDDO
        ENDDO

      ELSEIF( e2 == RADIATION )THEN

        DO iEq = 1, myDGSEM % state % nEquations
          DO j = 0, myDGSEM % params % polyDeg
            DO i = 0, myDGSEM % params % polyDeg
              myDGSEM % state % externalState(i,j,iEq,bID) = 0.0_prec
            ENDDO
          ENDDO
        ENDDO

      ELSEIF( e2 == NO_NORMAL_FLOW )THEN

        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg
            ! normal
            nx = myDGSEM % mesh % elements % nHat(1,i,j,s1,e1)
            ny = myDGSEM % mesh % elements % nHat(2,i,j,s1,e1)
            nz = myDGSEM % mesh % elements % nHat(3,i,j,s1,e1)
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

            myDGSEM % state % externalState(i,j,1,bID) = -nx*un + us*sx + ut*tx ! u
            myDGSEM % state % externalState(i,j,2,bID) = -ny*un + us*sy + ut*ty ! v
            myDGSEM % state % externalState(i,j,3,bID) = -nz*un + us*sz + ut*tz ! w
            myDGSEM % state % externalState(i,j,4,bID) =  myDGSEM % state % boundarySolution(i,j,4,s1,e1) ! rho
            myDGSEM % state % externalState(i,j,5,bID) =  myDGSEM % state % boundarySolution(i,j,5,s1,e1) ! potential temperature
            myDGSEM % state % externalState(i,j,6,bID) =  myDGSEM % state % boundarySolution(i,j,6,s1,e1) 

            myDGSEM % state % externalState(i,j,nEquations,bID) =  myDGSEM % state % boundarySolution(i,j,nEquations,s1,e1) ! P

          ENDDO
        ENDDO

      ELSEIF( e2 == INFLOW )THEN

        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg
            ! normal
            nx = myDGSEM % mesh % elements % nHat(1,i,j,s1,e1)
            ny = myDGSEM % mesh % elements % nHat(2,i,j,s1,e1)
            nz = myDGSEM % mesh % elements % nHat(3,i,j,s1,e1)
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
                 myDGSEM % state % boundarySolution(i,j,3,s1,e1)*nz + &
                 myDGSEM % state % prescribedState(i,j,1,bID)
           
            us = myDGSEM % state % boundarySolution(i,j,1,s1,e1)*sx    + &
              myDGSEM % state % boundarySolution(i,j,2,s1,e1)*sy    + &
              myDGSEM % state % boundarySolution(i,j,3,s1,e1)*sz

            ut = myDGSEM % state % boundarySolution(i,j,1,s1,e1)*tx  + &
              myDGSEM % state % boundarySolution(i,j,2,s1,e1)*ty  + &
              myDGSEM % state % boundarySolution(i,j,3,s1,e1)*tz


            myDGSEM % state % externalState(i,j,1,bID) = -nx*un + us*sx + ut*tx ! u
            myDGSEM % state % externalState(i,j,2,bID) = -ny*un + us*sy + ut*ty ! v
            myDGSEM % state % externalState(i,j,3,bID) = -nz*un + us*sz + ut*tz ! w
            myDGSEM % state % externalState(i,j,4,bID) =  myDGSEM % state % prescribedState(i,j,4,bID) ! rho
            myDGSEM % state % externalState(i,j,5,bID) =  myDGSEM % state % prescribedState(i,j,5,bID) ! potential temperature
            myDGSEM % state % externalState(i,j,6,bID) =  myDGSEM % state % prescribedState(i,j,6,bID)
            myDGSEM % state % externalState(i,j,nEquations,bID) =  myDGSEM % state % prescribedState(i,j,nEquations,bID) ! P
          ENDDO
        ENDDO

      ENDIF

    ENDDO

#endif

  END SUBROUTINE Update_FluidState_BCs
!
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!                                                                                                 !
!  This section of code CONTAINS routines for computing the fluxes through the element faces      !
!                                                                                                 !
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!
  SUBROUTINE InternalFace_StateFlux_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(myDGSEM % params % polyDeg + 1, &
                  myDGSEM % params % polyDeg + 1, &
                  1 )
    grid = dim3(myDGSEM % mesh % faces % nFaces,1,1)

    CALL InternalFace_StateFlux_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mesh % faces % elementIDs_dev, &
                                                              myDGSEM % mesh % faces % elementSides_dev, &
                                                              myDGSEM % mesh % faces % boundaryID_dev, &
                                                              myDGSEM % mesh % faces % iMap_dev, &
                                                              myDGSEM % mesh % faces % jMap_dev, &
                                                              myDGSEM % mesh % elements % nHat_dev, &
                                                              myDGSEM % state % boundarySolution_dev, &
                                                              myDGSEM % static % boundarySolution_dev, &
                                                              myDGSEM % state % boundaryFlux_dev, &
                                                              myDGSEM % state % boundaryGradientFlux_dev )

#else
    ! Local
    INTEGER    :: iFace
    INTEGER    :: i, j, k, m, iEq
    INTEGER    :: ii, jj
    INTEGER    :: e1, s1, e2, s2, bID
    REAL(prec) :: nHat(1:3), norm
    REAL(prec) :: uOut, uIn, cIn, cOut, T
    REAL(prec) :: jump(1:myDGSEM % state % nEquations-1), aS(1:myDGSEM % state % nEquations-1)
    REAL(prec) :: fac, rC


    DO iFace = 1, myDGSEM % mesh % faces % nFaces

      e1 = myDGSEM % mesh % faces % elementIDs(1,iFace)
      s1 = myDGSEM % mesh % faces % elementSides(1,iFace)
      e2 = myDGSEM % mesh % faces % elementIDs(2,iFace)
      s2 = ABS(myDGSEM % mesh % faces % elementSides(2,iFace))
      bID = myDGSEM % mesh % faces % boundaryID(iFace)

      IF( bID == 0 )THEN

        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg

            ii = myDGSEM % mesh % faces % iMap(i,j,iFace)
            jj = myDGSEM % mesh % faces % jMap(i,j,iFace)

            norm = sqrt( myDGSEM % mesh % elements % nHat(1,i,j,s1,e1)*myDGSEM % mesh % elements % nHat(1,i,j,s1,e1) + &
              myDGSEM % mesh % elements % nHat(2,i,j,s1,e1)*myDGSEM % mesh % elements % nHat(2,i,j,s1,e1) + &
              myDGSEM % mesh % elements % nHat(3,i,j,s1,e1)*myDGSEM % mesh % elements % nHat(3,i,j,s1,e1) )

            DO k = 1, 3
              nHat(k) = myDGSEM % mesh % elements % nHat(k,i,j,s1,e1)/norm
            ENDDO

            DO iEq = 1, myDGSEM % state % nEquations-1
              jump(iEq)  = myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2) - &
                           myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) !outState - inState
            ENDDO

            T =   (myDGSEM % static % boundarySolution(ii,jj,5,s2,e2) + myDGSEM % state % boundarySolution(ii,jj,5,s2,e2))/&
                  (myDGSEM % static % boundarySolution(ii,jj,4,s2,e2) + myDGSEM % state % boundarySolution(ii,jj,4,s2,e2))

            ! Sound speed estimate for the external and internal states
            cOut = sqrt( myDGSEM % params % R*T )

            T =   (myDGSEM % static % boundarySolution(i,j,5,s1,e1) + myDGSEM % state % boundarySolution(i,j,5,s1,e1))/&
                  (myDGSEM % static % boundarySolution(i,j,4,s1,e1) + myDGSEM % state % boundarySolution(i,j,4,s1,e1))

            cIn = sqrt( myDGSEM % params % R*T )
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

            fac = max( abs(uIn+cIn), abs(uIn-cIn), abs(uOut+cOut), abs(uOut-cOut) )

            ! Advective flux
            DO iEq = 1, myDGSEM % state % nEquations-1
              aS(iEq) = uIn*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) + &
                              myDGSEM % static % boundarySolution(i,j,iEq,s1,e1) ) +&
                       uOut*( myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2) + &
                              myDGSEM % static % boundarySolution(ii,jj,iEq,s2,e2) )
            ENDDO

            DO k = 1, 3
              ! Momentum flux due to pressure
              aS(k) = aS(k) + (myDGSEM % state % boundarySolution(i,j,nEquations,s1,e1) + &
                               myDGSEM % state % boundarySolution(ii,jj,nEquations,s2,e2))*nHat(k)
            ENDDO


            DO iEq = 1, myDGSEM % state % nEquations-1

              myDGSEM % state % boundaryFlux(i,j,iEq,s1,e1) =  0.5_prec*( aS(iEq) - fac*jump(iEq) )*norm
              myDGSEM % state % boundaryFlux(ii,jj,iEq,s2,e2) = -myDGSEM % state % boundaryFlux(i,j,iEq,s1,e1)

              IF( iEq == 4 )THEN

                DO k = 1, 3
                  ! Calculate the LDG flux for the stress tensor.
                  myDGSEM % state % boundaryGradientFlux(k,i,j,iEq,s1,e1) = 0.5_prec*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) +&
                                                                                       myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2) )*&
                                                                                       myDGSEM % mesh % elements % nHat(k,i,j,s1,e1)

                  myDGSEM % state % boundaryGradientFlux(k,ii,jj,iEq,s2,e2) = -myDGSEM % state % boundaryGradientFlux(k,i,j,iEq,s1,e1)
                ENDDO

              ELSE
                DO k = 1, 3
                  ! Calculate the LDG flux for the stress tensor.
                  myDGSEM % state % boundaryGradientFlux(k,i,j,iEq,s1,e1) = 0.5_prec*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1)/&
                                                                                       (myDGSEM % state % boundarySolution(i,j,4,s1,e1) +&
                                                                                        myDGSEM % static % boundarySolution(i,j,4,s1,e1))+&
                                                                                       myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2)/&
                                                                                       (myDGSEM % state % boundarySolution(ii,jj,4,s2,e2) +&
                                                                                        myDGSEM % static % boundarySolution(ii,jj,4,s2,e2)) )*&
                                                                                       myDGSEM % mesh % elements % nHat(k,i,j,s1,e1)

                  myDGSEM % state % boundaryGradientFlux(k,ii,jj,iEq,s2,e2) = -myDGSEM % state % boundaryGradientFlux(k,i,j,iEq,s1,e1)
                ENDDO
              ENDIF

            ENDDO
          ENDDO
        ENDDO

      ENDIF

    ENDDO

#endif

  END SUBROUTINE InternalFace_StateFlux_Fluid
!
  SUBROUTINE BoundaryFace_StateFlux_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1, &
                  1 )
    grid = dim3(myDGSEM % mesh % faces % nFaces,1,1)

    CALL BoundaryFace_StateFlux_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mesh % faces % elementIDs_dev, &
                                                              myDGSEM % mesh % faces % elementSides_dev, &
                                                              myDGSEM % mesh % faces % boundaryID_dev, &
                                                              myDGSEM % mesh % faces % iMap_dev, &
                                                              myDGSEM % mesh % faces % jMap_dev, &
                                                              myDGSEM % mesh % elements % nHat_dev, &
                                                              myDGSEM % state % boundarySolution_dev, &
                                                              myDGSEM % static % boundarySolution_dev, &
                                                              myDGSEM % state % externalState_dev, &
                                                              myDGSEM % static % externalState_dev, &
                                                              myDGSEM % state % boundaryFlux_dev, &
                                                              myDGSEM % state % boundaryGradientFlux_dev )


#else
    ! Local
    INTEGER :: iEl, iFace
    INTEGER    :: i, j, k, m, iEq
    INTEGER    :: ii, jj, bID
    INTEGER    :: e1, s1, e2, s2
    REAL(prec) :: nHat(1:3), norm
    REAL(prec) :: uOut, uIn, cIn, cOut, T
    REAL(prec) :: jump(1:myDGSEM % state % nEquations-1), aS(1:myDGSEM % state % nEquations-1)
    REAL(prec) :: fac, hCapRatio, rC

    DO iFace = 1, myDGSEM % mesh % faces % nFaces

      e1 = myDGSEM % mesh % faces % elementIDs(1,iFace)
      s1 = myDGSEM % mesh % faces % elementSides(1,iFace)
      e2 = myDGSEM % mesh % faces % elementIDs(2,iFace)
      s2 = ABS(myDGSEM % mesh % faces % elementSides(2,iFace))
      bID  = myDGSEM % mesh % faces % boundaryID(iFace)

      IF( e2 < 0 )THEN

        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg

            ii = myDGSEM % mesh % faces % iMap(i,j,iFace)
            jj = myDGSEM % mesh % faces % jMap(i,j,iFace)

            norm = sqrt( myDGSEM % mesh % elements % nHat(1,i,j,s1,e1)*myDGSEM % mesh % elements % nHat(1,i,j,s1,e1) + &
              myDGSEM % mesh % elements % nHat(2,i,j,s1,e1)*myDGSEM % mesh % elements % nHat(2,i,j,s1,e1) + &
              myDGSEM % mesh % elements % nHat(3,i,j,s1,e1)*myDGSEM % mesh % elements % nHat(3,i,j,s1,e1) )

            DO k = 1, 3
              nHat(k) = myDGSEM % mesh % elements % nHat(k,i,j,s1,e1)/norm
            ENDDO

            DO iEq = 1, myDGSEM % state % nEquations-1
              jump(iEq)  = myDGSEM % state % externalState(ii,jj,iEq,bID) - &
                           myDGSEM % state % boundarySolution(i,j,iEq,s1,e1)
            ENDDO

            ! Sound speed estimate for the external and internal states

            T = (myDGSEM % static % externalState(ii,jj,5,bID)+myDGSEM % state % externalState(ii,jj,5,bID))/&
                (myDGSEM % static % externalState(ii,jj,4,bID)+myDGSEM % state % externalState(ii,jj,4,bID))

            cOut = sqrt( myDGSEM % params % R*T )

            T = (myDGSEM % static % boundarySolution(i,j,5,s1,e1)+myDGSEM % state % boundarySolution(i,j,5,s1,e1))/&
                (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % boundarySolution(i,j,4,s1,e1))

            cIn = sqrt( myDGSEM % params % R*T )

            uOut = ( myDGSEM % state % externalState(ii,jj,1,bID)*nHat(1) + &
                     myDGSEM % state % externalState(ii,jj,2,bID)*nHat(2) + &
                     myDGSEM % state % externalState(ii,jj,3,bID)*nHat(3) )/&
                   ( myDGSEM % state % externalState(ii,jj,4,bID)+&
                     myDGSEM % static % externalState(ii,jj,4,bID) )

            ! Internal normal velocity component
            uIn  = ( myDGSEM % state % boundarySolution(i,j,1,s1,e1)*nHat(1) + &
                     myDGSEM % state % boundarySolution(i,j,2,s1,e1)*nHat(2) + &
                     myDGSEM % state % boundarySolution(i,j,3,s1,e1)*nHat(3) )/&
                   ( myDGSEM % state % boundarySolution(i,j,4,s1,e1)+&
                     myDGSEM % static % boundarySolution(i,j,4,s1,e1) )

            fac = max( abs(uIn+cIn), abs(uIn-cIn), abs(uOut+cOut), abs(uOut-cOut) )

            DO iEq = 1, myDGSEM % state % nEquations-1
              aS(iEq) = uIn*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) + myDGSEM % static % boundarySolution(i,j,iEq,s1,e1) ) +&
                        uOut*( myDGSEM % state % externalState(ii,jj,iEq,bID) + myDGSEM % static % externalState(ii,jj,iEq,bID) )
            ENDDO

            DO k = 1, 3
              ! Momentum flux due to pressure
              aS(k) = aS(k) + (myDGSEM % state % boundarySolution(i,j,nEquations,s1,e1) + &
                myDGSEM % state % externalState(ii,jj,nEquations,bID))*nHat(k)
            ENDDO


            DO iEq = 1, myDGSEM % state % nEquations-1

              myDGSEM % state % boundaryFlux(i,j,iEq,s1,e1) =  0.5_prec*( aS(iEq) - fac*jump(iEq) )*norm

              DO k = 1, 3

                IF( iEq == 4 )THEN
                  ! Calculate the Bassi-Rebay (LDG) flux for the stress tensor.
                  myDGSEM % state % boundaryGradientFlux(k,i,j,iEq,s1,e1) = 0.5_prec*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) +&
                                                                                    myDGSEM % state % externalState(ii,jj,iEq,bID) )*&
                                                                                    myDGSEM % mesh % elements % nHat(k,i,j,s1,e1)
                                                                                   
                ELSE
                  ! Calculate the Bassi-Rebay (LDG) flux for the stress tensor.
                  myDGSEM % state % boundaryGradientFlux(k,i,j,iEq,s1,e1) = 0.5_prec*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1)/&
                                                                                    (myDGSEM % state % boundarySolution(i,j,4,s1,e1) +&
                                                                                     myDGSEM % static % boundarySolution(i,j,4,s1,e1))+&
                                                                                    myDGSEM % state % externalState(ii,jj,iEq,bID)/&
                                                                                    (myDGSEM % state % externalState(ii,jj,4,bID) +&
                                                                                     myDGSEM % static % externalState(ii,jj,4,bID)) )*&
                                                                                    myDGSEM % mesh % elements % nHat(k,i,j,s1,e1)
                                                                               
                ENDIF
                     
              ENDDO
            ENDDO

          ENDDO
        ENDDO

      ENDIF

    ENDDO


#endif

  END SUBROUTINE BoundaryFace_StateFlux_Fluid

!
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!                                                                                                 !
!  This section of code CONTAINS routines for computing tendency from the internal and Riemann    !
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
#else
    INTEGER    :: iEl, i, j, k, m, iEq, row, col, ii
    REAL(prec) :: df, F
#endif



#ifdef HAVE_CUDA

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1 )
    grid = dim3(myDGSEM % state % nEquations-1,myDGSEM % state % nElements,1)

    CALL CalculateFlux_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                    myDGSEM % static % solution_dev, &
                                                    myDGSEM % mesh % elements % Ja_dev, &
                                                    myDGSEM % mesh % elements % J_dev, &
                                                    myDGSEM % state % flux_dev )

    CALL State_Mapped_DG_Divergence_3D_CUDAKernel<<<grid, tBlock>>>( myDGSEM % state % flux_dev, &
                                                               myDGSEM % state % boundaryFlux_dev, &
                                                               myDGSEM % state % fluxDivergence_dev, &
                                                               myDGSEM % dgStorage % boundaryInterpolationMatrix_dev, &
                                                               myDGSEM % dgStorage % dgDerivativeMatrixTranspose_dev, &
                                                               myDGSEM % dgStorage % quadratureWeights_dev, &
                                                               myDGSEM % mesh % elements % J_dev ) 

    CALL CalculateSourceTerms_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                           myDGSEM % state % solutionGradient_dev, &
                                                           myDGSEM % static % solution_dev, &
                                                           myDGSEM % static % source_dev, &
                                                           myDGSEM % state % source_dev, &
                                                           myDGSEM % sourceTerms % drag_dev )
  
#else


    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO iEq = 1, myDGSEM % state % nEquations-1

        DO k = 0, myDGSEM % params % polyDeg
          DO j = 0, myDGSEM % params % polyDeg
            DO i = 0, myDGSEM % params % polyDeg
              myDGSEM % state % flux(1,i,j,k,iEq,iEl) = 0.0_prec
              myDGSEM % state % flux(2,i,j,k,iEq,iEl) = 0.0_prec
              myDGSEM % state % flux(3,i,j,k,iEq,iEl) = 0.0_prec
            ENDDO
          ENDDO
        ENDDO

        ! Here the flux tensor in physical space is calculated and rotated to give the
        ! contravariant flux tensor in the reference computational DOmain.

        !//////////////////////////////// Advection ///////////////////////////////////////!
        DO col = 1, 3
          DO k = 0, myDGSEM % params % polyDeg
            DO j = 0, myDGSEM % params % polyDeg
              DO i = 0, myDGSEM % params % polyDeg

                F = myDGSEM % state % solution(i,j,k,col,iEl)*&
                  (myDGSEM % state % solution(i,j,k,iEq,iEl)+&
                  myDGSEM % static % solution(i,j,k,iEq,iEl))/&
                  (myDGSEM % state % solution(i,j,k,4,iEl)+&
                  myDGSEM % static % solution(i,j,k,4,iEl))


                myDGSEM % state % flux(1,i,j,k,iEq,iEl) = myDGSEM % state % flux(1,i,j,k,iEq,iEl)+ &
                  myDGSEM % mesh % elements % Ja(i,j,k,col,1,iEl)*F

                myDGSEM % state % flux(2,i,j,k,iEq,iEl) = myDGSEM % state % flux(2,i,j,k,iEq,iEl) + &
                  myDGSEM % mesh % elements % Ja(i,j,k,col,2,iEl)*F

                myDGSEM % state % flux(3,i,j,k,iEq,iEl) = myDGSEM % state % flux(3,i,j,k,iEq,iEl)  + &
                  myDGSEM % mesh % elements % Ja(i,j,k,col,3,iEl)*F

              ENDDO
            ENDDO
          ENDDO
        ENDDO

        ! //////////////////// Pressure (Momentum only) /////////////////////////// !
        IF( iEq <= 3 )THEN

          DO k = 0, myDGSEM % params % polyDeg
            DO j = 0, myDGSEM % params % polyDeg
              DO i = 0, myDGSEM % params % polyDeg

                myDGSEM % state % flux(1,i,j,k,iEq,iEl) = myDGSEM % state % flux(1,i,j,k,iEq,iEl) + &
                  myDGSEM % mesh % elements % Ja(i,j,k,iEq,1,iEl)*&
                  myDGSEM % state % solution(i,j,k,nEquations,iEl)

                myDGSEM % state % flux(2,i,j,k,iEq,iEl) = myDGSEM % state % flux(2,i,j,k,iEq,iEl) + &
                  myDGSEM % mesh % elements % Ja(i,j,k,iEq,2,iEl)*&
                  myDGSEM % state % solution(i,j,k,nEquations,iEl)

                myDGSEM % state % flux(3,i,j,k,iEq,iEl) = myDGSEM % state % flux(3,i,j,k,iEq,iEl) + &
                  myDGSEM % mesh % elements % Ja(i,j,k,iEq,3,iEl)*&
                  myDGSEM % state % solution(i,j,k,nEquations,iEl)

              ENDDO
            ENDDO
          ENDDO

        ENDIF

      ENDDO
    ENDDO

    DO iEl = 1, myDGSEM % state % nElements
      DO iEq = 1, myDGSEM % state % nEquations
        DO k = 0, myDGSEM % state % N
          DO j = 0, myDGSEM % state % N
            DO i = 0, myDGSEM % state % N   
 
              df = 0.0_prec
              DO ii = 0, myDGSEM % dgStorage % N
                df = df + myDGSEM % dgStorage % dgDerivativeMatrixTranspose(ii,i)*myDGSEM % state % flux(1,ii,j,k,iEq,iEl) + &
                          myDGSEM % dgStorage % dgDerivativeMatrixTranspose(ii,j)*myDGSEM % state % flux(2,i,ii,k,iEq,iEl) + &
                          myDGSEM % dgStorage % dgDerivativeMatrixTranspose(ii,k)*myDGSEM % state % flux(3,i,j,ii,iEq,iEl)
              ENDDO
               
              myDGSEM % state % fluxDivergence(i,j,k,iEq,iEl) =  ( df+ ( myDGSEM % state % boundaryFlux(i,k,iEq,1,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(j,0) + &
                                              myDGSEM % state % boundaryFlux(i,k,iEq,3,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(j,1) )/&
                                            myDGSEM % dgStorage % quadratureWeights(j) + &
                                            ( myDGSEM % state % boundaryFlux(j,k,iEq,4,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(i,0) + &
                                              myDGSEM % state % boundaryFlux(j,k,iEq,2,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(i,1) )/&
                                            myDGSEM % dgStorage % quadratureWeights(i) + &
                                            ( myDGSEM % state % boundaryFlux(i,j,iEq,5,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(k,0) + &
                                              myDGSEM % state % boundaryFlux(i,j,iEq,6,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(k,1) )/&
                                            myDGSEM % dgStorage % quadratureWeights(k) )/myDGSEM % mesh % elements % J(i,j,k,iEl) 

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO iEq = 1, myDGSEM % state % nEquations-1

        IF( iEq == 1 )THEN

          DO k = 0, myDGSEM % params % polyDeg
            DO j = 0, myDGSEM % params % polyDeg
              DO i = 0, myDGSEM % params % polyDeg

                F = sqrt( myDGSEM % state % solution(i,j,k,1,iEl)**2 + &
                  myDGSEM % state % solution(i,j,k,2,iEl)**2 + &
                  myDGSEM % state % solution(i,j,k,3,iEl)**2 )/&
                  ( myDGSEM % static % solution(i,j,k,4,iEl) + myDGSEM % state % solution(i,j,k,4,iEl) )

                myDGSEM % state % source(i,j,k,1,iEl) = -myDGSEM % sourceTerms % drag(i,j,k,iEl)*myDGSEM % state % solution(i,j,k,1,iEl)*F -&
                                                        myDGSEM % state % solution(i,j,k,3,iEl)*myDGSEM % params % fRotY +&
                                                        myDGSEM % state % solution(i,j,k,2,iEl)*myDGSEM % params % fRotZ

              ENDDO
            ENDDO
          ENDDO

        ELSEIF( iEq == 2 )THEN

          DO k = 0, myDGSEM % params % polyDeg
            DO j = 0, myDGSEM % params % polyDeg
              DO i = 0, myDGSEM % params % polyDeg

                F = sqrt( myDGSEM % state % solution(i,j,k,1,iEl)**2 + &
                  myDGSEM % state % solution(i,j,k,2,iEl)**2 + &
                  myDGSEM % state % solution(i,j,k,3,iEl)**2 )/&
                  ( myDGSEM % static % solution(i,j,k,4,iEl) + myDGSEM % state % solution(i,j,k,4,iEl) )

                myDGSEM % state % source(i,j,k,2,iEl) = -myDGSEM % sourceTerms % drag(i,j,k,iEl)*myDGSEM % state % solution(i,j,k,2,iEl)*F - &
                                                        myDGSEM % state % solution(i,j,k,1,iEl)*myDGSEM % params % fRotZ +&
                                                        myDGSEM % state % solution(i,j,k,3,iEl)*myDGSEM % params % fRotX

              ENDDO
            ENDDO
          ENDDO

        ELSEIF( iEq == 3 )THEN

          DO k = 0, myDGSEM % params % polyDeg
            DO j = 0, myDGSEM % params % polyDeg
              DO i = 0, myDGSEM % params % polyDeg

                F = sqrt( myDGSEM % state % solution(i,j,k,1,iEl)**2 + &
                  myDGSEM % state % solution(i,j,k,2,iEl)**2 + &
                  myDGSEM % state % solution(i,j,k,3,iEl)**2 )/&
                  ( myDGSEM % static % solution(i,j,k,4,iEl) + myDGSEM % state % solution(i,j,k,4,iEl) )

                myDGSEM % state % source(i,j,k,3,iEl) = -myDGSEM % sourceTerms % drag(i,j,k,iEl)*myDGSEM % state % solution(i,j,k,3,iEl)*F - &
                                                        myDGSEM % state % solution(i,j,k,2,iEl)*myDGSEM % params % fRotX +&
                                                        myDGSEM % state % solution(i,j,k,1,iEl)*myDGSEM % params % fRotY-&
                                                        ( myDGSEM % state % solution(i,j,k,4,iEl) )*myDGSEM % params % g  !&

              ENDDO
            ENDDO
          ENDDO

        ! When the in-situ temperature formulation is used, we must add in the adiabatic heating term
        ! due to fluid convergences. This term is -( P_{total}/C_v )*div( u )
        ELSEIF( iEq == 5 )THEN

          DO k = 0, myDGSEM % params % polyDeg
            DO j = 0, myDGSEM % params % polyDeg
              DO i = 0, myDGSEM % params % polyDeg

                myDGSEM % state % source(i,j,k,5,iEl) = -( myDGSEM % static % solution(i,j,k,nEquations,iEl) + myDGSEM % state % solution(i,j,k,nEquations,iEl) )*&
                                                         ( myDGSEM % state % solutionGradient(1,i,j,k,1,iEl) + &
                                                           myDGSEM % state % solutionGradient(2,i,j,k,2,iEl) + &
                                                           myDGSEM % state % solutionGradient(3,i,j,k,3,iEl) )/myDGSEM % params % Cv

              ENDDO
            ENDDO
          ENDDO
   
        ENDIF

      ENDDO
    ENDDO
#endif


  END SUBROUTINE MappedTimeDerivative_Fluid
!
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!                                                                                                 !
!  This section of code CONTAINS routines for computing the gradients of the prognostic variables !
!                                                                                                 !
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!
  SUBROUTINE CalculateSolutionGradient_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock
#else
    INTEGER :: iEl, iEq, idir, i, j, k, m, ii
    REAL(prec) :: f(1:3,0:myDGSEM % params % polydeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg)
    REAL(prec) :: df
#endif



#ifdef HAVE_CUDA
    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1 )
    grid = dim3(myDGSEM % state % nEquations,myDGSEM % mesh % elements % nElements,3)

    CALL CalculateSolutionGradient_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solutionGradient_dev, &
                                                                myDGSEM % state % solution_dev, &
                                                                myDGSEM % static % solution_dev, &
                                                                myDGSEM % state % boundaryGradientFlux_dev, &
                                                                myDGSEM % mesh % elements % Ja_dev, &
                                                                myDGSEM % mesh % elements % J_dev, &
                                                                myDGSEM % dgStorage % dgDerivativeMatrixTranspose_dev, &
                                                                myDGSEM % dgStorage % boundaryInterpolationMatrix_dev, &
                                                                myDGSEM % dgStorage % quadratureWeights_dev )
#else

    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO iEq = 1, myDGSEM % state % nEquations

        DO idir = 1, 3

          IF( iEq == 4 )THEN

            DO k = 0, myDGSEM % params % polyDeg
              DO j = 0, myDGSEM % params % polyDeg
                DO i = 0, myDGSEM % params % polyDeg

                  f(1,i,j,k) = myDGSEM % state % solution(i,j,k,iEq,iEl)*myDGSEM % mesh % elements % Ja(i,j,k,idir,1,iEl)
                  f(2,i,j,k) = myDGSEM % state % solution(i,j,k,iEq,iEl)*myDGSEM % mesh % elements % Ja(i,j,k,idir,2,iEl)
                  f(3,i,j,k) = myDGSEM % state % solution(i,j,k,iEq,iEl)*myDGSEM % mesh % elements % Ja(i,j,k,idir,3,iEl)

                ENDDO
              ENDDO
            ENDDO

          ELSE

            DO k = 0, myDGSEM % params % polyDeg
              DO j = 0, myDGSEM % params % polyDeg
                DO i = 0, myDGSEM % params % polyDeg

                  f(1,i,j,k) = myDGSEM % mesh % elements % Ja(i,j,k,idir,1,iEl)*&
                               myDGSEM % state % solution(i,j,k,iEq,iEl)/&
                               (myDGSEM % state % solution(i,j,k,4,iEl)+&
                               myDGSEM % static % solution(i,j,k,4,iEl) )


                  f(2,i,j,k) = myDGSEM % mesh % elements % Ja(i,j,k,idir,2,iEl)*&
                               myDGSEM % state % solution(i,j,k,iEq,iEl)/&
                               (myDGSEM % state % solution(i,j,k,4,iEl)+&
                               myDGSEM % static % solution(i,j,k,4,iEl) )

                  f(3,i,j,k) = myDGSEM % mesh % elements % Ja(i,j,k,idir,3,iEl)*&
                               myDGSEM % state % solution(i,j,k,iEq,iEl)/&
                               (myDGSEM % state % solution(i,j,k,4,iEl)+&
                               myDGSEM % static % solution(i,j,k,4,iEl) )

                ENDDO
              ENDDO
            ENDDO

          ENDIF

          DO k = 0, myDGSEM % params % polyDeg
            DO j = 0, myDGSEM % params % polyDeg
              DO i = 0, myDGSEM % params % polyDeg

                df = 0.0_prec
                DO ii = 0, myDGSEM % params % polyDeg
                  df = df + myDGSEM % dgStorage % dgDerivativeMatrixTranspose(ii,i)*f(1,ii,j,k) + &
                            myDGSEM % dgStorage % dgDerivativeMatrixTranspose(ii,j)*f(2,i,ii,k) + &
                            myDGSEM % dgStorage % dgDerivativeMatrixTranspose(ii,k)*f(3,i,j,ii)
                ENDDO

                myDGSEM % state % solutionGradient(idir,i,j,k,iEq,iEl) =  ( df+ ( myDGSEM % state % boundaryGradientFlux(idir,i,k,iEq,1,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(j,0) + &
                                                                              myDGSEM % state % boundaryGradientFlux(idir,i,k,iEq,3,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(j,1) )/&
                                                                            myDGSEM % dgStorage % quadratureWeights(j) + &
                                                                            ( myDGSEM % state % boundaryGradientFlux(idir,j,k,iEq,4,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(i,0) + &
                                                                              myDGSEM % state % boundaryGradientFlux(idir,j,k,iEq,2,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(i,1) )/&
                                                                            myDGSEM % dgStorage % quadratureWeights(i) + &
                                                                            ( myDGSEM % state % boundaryGradientFlux(idir,i,j,iEq,5,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(k,0) + &
                                                                              myDGSEM % state % boundaryGradientFlux(idir,i,j,iEq,6,iEl)*myDGSEM % dgStorage % boundaryInterpolationMatrix(k,1) )/&
                                                                            myDGSEM % dgStorage % quadratureWeights(k) )/myDGSEM % mesh % elements % J(i,j,k,iEl)

              ENDDO
            ENDDO
          ENDDO

        ENDDO
      ENDDO
    ENDDO

#endif

  END SUBROUTINE CalculateSolutionGradient_Fluid
!
  SUBROUTINE CalculateNormalStressAtBoundaries_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock
#else
    INTEGER :: iEl, iEq, idir, i, j, k, m
    REAL(prec) :: fAtBoundaries(1:3,1:6)
#endif
  
#ifdef HAVE_CUDA

      tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                    1 )
      grid = dim3(myDGSEM % stressTensor % nEquations, myDGSEM % stressTensor % nElements, 1)  

      CALL CalculateNormalStressAtBoundaries_CUDAKernel<<<grid, tBlock>>>( myDGSEM % state % solutionGradient_dev, &
                                                                           myDGSEM % mesh % elements % nHat_dev, &
                                                                           myDGSEM % stressTensor % boundarySolution_dev,  &
                                                                           myDGSEM % dgStorage % boundaryInterpolationMatrix_dev )

#else


    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO iEq = 1, myDGSEM % stressTensor % nEquations

          DO j = 0, myDGSEM % params % polyDeg
            DO i = 0, myDGSEM % params % polyDeg
            
              fAtBoundaries(1:3,1:6) = 0.0_prec
              
              DO k = 0, myDGSEM % params % polyDeg
               
                fAtBoundaries(1:3,1) = fAtBoundaries(1:3,1) + myDGSEM % dgStorage % boundaryInterpolationMatrix(k,0)*myDGSEM % state % solutionGradient(1:3,i,k,j,iEq,iEl) ! South
                fAtBoundaries(1:3,2) = fAtBoundaries(1:3,2) + myDGSEM % dgStorage % boundaryInterpolationMatrix(k,1)*myDGSEM % state % solutionGradient(1:3,k,i,j,iEq,iEl) ! East
                fAtBoundaries(1:3,3) = fAtBoundaries(1:3,3) + myDGSEM % dgStorage % boundaryInterpolationMatrix(k,1)*myDGSEM % state % solutionGradient(1:3,i,k,j,iEq,iEl) ! North
                fAtBoundaries(1:3,4) = fAtBoundaries(1:3,4) + myDGSEM % dgStorage % boundaryInterpolationMatrix(k,0)*myDGSEM % state % solutionGradient(1:3,k,i,j,iEq,iEl) ! West
                fAtBoundaries(1:3,5) = fAtBoundaries(1:3,5) + myDGSEM % dgStorage % boundaryInterpolationMatrix(k,0)*myDGSEM % state % solutionGradient(1:3,i,j,k,iEq,iEl) ! Bottom
                fAtBoundaries(1:3,6) = fAtBoundaries(1:3,6) + myDGSEM % dgStorage % boundaryInterpolationMatrix(k,1)*myDGSEM % state % solutionGradient(1:3,i,j,k,iEq,iEl) ! Top
               
              ENDDO

              DO k = 1, 6

                myDGSEM % stressTensor % boundarySolution(i,j,iEq,k,iEl) = myDGSEM % params % viscosity*&
                                                                           ( fAtBoundaries(1,k)*myDGSEM % mesh % elements % nHat(1,i,j,k,iEL) + &
                                                                             fAtBoundaries(2,k)*myDGSEM % mesh % elements % nHat(2,i,j,k,iEL) + &
                                                                             fAtBoundaries(3,k)*myDGSEM % mesh % elements % nHat(3,i,j,k,iEL) )
          
              
              ENDDO

            ENDDO
          ENDDO

      ENDDO
    ENDDO
  

#endif
  END SUBROUTINE CalculateNormalStressAtBoundaries_Fluid
!
  SUBROUTINE CalculateStressFlux_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock
#else
    INTEGER :: iEl, iEq, i, j, k, idir, jdir, ii 
    REAL(prec) :: flux(1:3)
#endif


#ifdef HAVE_CUDA
    tBlock = dim3( myDGSEM % params % polyDeg+1, &
                   myDGSEM % params % polyDeg+1, &
                   myDGSEM % params % polyDeg+1 )

    grid = dim3( myDGSEM % stressTensor % nEquations, myDGSEM % stressTensor % nElements, 1)  

    CALL CalculateStressFlux_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solutionGradient_dev, &
                                                          myDGSEM % state % solution_dev, &
                                                          myDGSEM % static % solution_dev, &
                                                          myDGSEM % mesh % elements % Ja_dev, &
                                                          myDGSEM % stressTensor % flux_dev )

#else

    ! FLOPS = 3*[ (polyDeg+1)^3 ]*[3]*[5]*[nElements]
    ! # of array reads = 12*[ (polyDeg+1)^3 ]*[3]*[5]*[nElements] 
    !
    ! > For the thermal bubble demo @ single precisions, default settings,
    !     FLOPS = 3,456,000
    !     reads = 13,824,000
    !     writes = 1,152,000
    !
    !     Oryx, Intel I7, TPP ~ 12.84 GFLOPs/s, (Peak compute bandwidth of 51.31 GB/s @ single precision)
    !
    !     Current serial runtime is 0.776 s, 
    !         >> Performance is          4.45 FLOPs/s
    !         >> read/write bandwidth of 19.30 MB/s
    !
    !      Performance ~ Arithmetic Intensity * Bandwidth
    !      FLOPS/s ~ FLOPS/Bytes * ( Bytes/second )
    !      Theoretical arithmetic intensity for this routine is  ~ 0.231 (FLOPS/Byte)
    !
    !
    !

    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO iEq = 1, myDGSEM % stressTensor % nEquations
        DO k = 0, myDGSEM % params % polyDeg
          DO j = 0, myDGSEM % params % polyDeg
            DO i = 0, myDGSEM % params % polyDeg

              flux(1:3) = 0.0_prec

              DO idir = 1, 3


                flux(1) = flux(1) + myDGSEM % mesh % elements % Ja(i,j,k,idir,1,iEl)*&
                                          myDGSEM % state % solutionGradient(idir,i,j,k,iEq,iEl)*&
                                          myDGSEM % params % viscosity

                flux(2) = flux(2) + myDGSEM % mesh % elements % Ja(i,j,k,idir,2,iEl)*&
                                          myDGSEM % state % solutionGradient(idir,i,j,k,iEq,iEl)*&
                                          myDGSEM % params % viscosity

                flux(3) = flux(3) + myDGSEM % mesh % elements % Ja(i,j,k,idir,3,iEl)*&
                                          myDGSEM % state % solutionGradient(idir,i,j,k,iEq,iEl)*&
                                          myDGSEM % params % viscosity
              ENDDO

              myDGSEM % stressTensor % flux(1:3,i,j,k,iEq,iEl) = flux(1:3)

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

#endif

  END SUBROUTINE CalculateStressFlux_Fluid
!
!
  SUBROUTINE Update_FluidStress_BCs( myDGSEM, tn ) ! ////////// !

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
    REAL(prec), INTENT(in)      :: tn
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1, &
                  1 )
    grid = dim3(myDGSEM % mesh % faces % nFaces,myDGSEM % stressTensor % nEquations,1)

    CALL Update_FluidStress_BCs_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mesh % decomp % element_to_blockID_dev, &            ! I
                                                              myDGSEM % mesh % faces % boundaryID_dev, &
                                                              myDGSEM % mesh % faces % elementIDs_dev, &         ! I
                                                              myDGSEM % mesh % faces % elementSides_dev, &       ! I
                                                              myDGSEM % stressTensor % externalState_dev, &                    ! O
                                                              myDGSEM % stressTensor % boundarySolution_dev, &   ! I
                                                              myDGSEM % stressTensor % prescribedState_dev )
#else
    ! Local
    INTEGER    :: iEl, iFace, i, j, k, iEq
    INTEGER    :: bID, p2
    INTEGER    :: e1, e2, s1, s2

    DO iFace = 1, myDGSEM % mesh % faces % nFaces

      e1  = myDGSEM % mesh % faces % elementIDs(1,iFace)
      s1  = myDGSEM % mesh % faces % elementSides(1,iFace)
      e2  = myDGSEM % mesh % faces % elementIDs(2,iFace)
      bID = myDGSEM % mesh % faces % boundaryID(iFace)

      IF( e2 < 0 )THEN

        DO iEq = 1, myDGSEM % stressTensor % nEquations
          DO j = 0, myDGSEM % params % polyDeg
            DO i = 0, myDGSEM % params % polyDeg
              myDGSEM % stressTensor % externalState(i,j,iEq,bID) = myDGSEM % stressTensor % boundarySolution(i,j,iEq,s1,e1)
            ENDDO
          ENDDO
        ENDDO

      ENDIF

    ENDDO

#endif

  END SUBROUTINE Update_FluidStress_BCs
!
  SUBROUTINE InternalFace_StressFlux_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1, &
                  1 )
    grid = dim3(myDGSEM % mesh % faces % nFaces,myDGSEM % stressTensor % nEquations,1)

    CALL InternalFace_StressFlux_CUDAKernel<<<grid, tBlock>>>( myDGSEM % state % boundarySolution_dev, &
                                                               myDGSEM % stressTensor % boundarySolution_dev, &
                                                               myDGSEM % static % boundarySolution_dev, &
                                                               myDGSEM % mesh % elements % boundaryLengthScale_dev, &
                                                               myDGSEM % mesh % elements % nHat_dev, &
                                                               myDGSEM % stressTensor % boundaryFlux_dev, &
                                                               myDGSEM % mesh % faces % elementIDs_dev, &
                                                               myDGSEM % mesh % faces % elementSides_dev, &
                                                               myDGSEM % mesh % faces % boundaryID_dev, &
                                                               myDGSEM % mesh % faces % iMap_dev, &
                                                               myDGSEM % mesh % faces % jMap_dev )
#else
    ! Local
    REAL(prec) :: norm, rhoIn, rhoOut
    INTEGER :: iEl, iFace
    INTEGER    :: i, j, k, m, iEq
    INTEGER    :: ii, jj, bID
    INTEGER    :: e1, s1, e2, s2, p2

    DO iFace = 1, myDGSEM % mesh % faces % nFaces


      e1 = myDGSEM % mesh % faces % elementIDs(1,iFace)
      s1 = myDGSEM % mesh % faces % elementSides(1,iFace)
      e2 = myDGSEM % mesh % faces % elementIDs(2,iFace)
      s2 = ABS(myDGSEM % mesh % faces % elementSides(2,iFace))
      IF( e2 > 0 )THEN
        p2 = myDGSEM % mesh % decomp % element_to_blockID(e2)
      ELSE
        p2 = myRank 
      ENDIF

      IF( e2 > 0 .AND. p2 == myRank )THEN

        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg

            ii = myDGSEM % mesh % faces % iMap(i,j,iFace)
            jj = myDGSEM % mesh % faces % jMap(i,j,iFace)

            norm = sqrt( myDGSEM % mesh % elements % nHat(1,i,j,s1,e1)**2 + &
                         myDGSEM % mesh % elements % nHat(2,i,j,s1,e1)**2 + &
                         myDGSEM % mesh % elements % nHat(3,i,j,s1,e1)**2 )

            DO iEq = 1, myDGSEM % state % nEquations-1

              IF( iEq == 4 )THEN

                myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*( myDGSEM % stressTensor % boundarySolution(i,j,iEq,s1,e1)- &
                                                                                    myDGSEM % stressTensor % boundarySolution(ii,jj,iEq,s2,e2) )

              ELSE

                rhoOut = (myDGSEM % static % boundarySolution(ii,jj,4,s2,e2)+myDGSEM % state % boundarySolution(ii,jj,4,s2,e2) )
                rhoIn  = (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % boundarySolution(i,j,4,s1,e1) )

                myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*( rhoIn*myDGSEM % stressTensor % boundarySolution(i,j,iEq,s1,e1)-&
                                                                                  rhoOut*myDGSEM % stressTensor % boundarySolution(ii,jj,iEq,s2,e2) )


              ENDIF

              myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) +&
                0.5_prec*( myDGSEM % params % viscosity*myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2)-&
                           myDGSEM % params % viscosity*myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) )/&
                ( 0.5_prec*( myDGSEM % mesh % elements % boundaryLengthScale(i,j,s1,e1) + myDGSEM % mesh % elements % boundaryLengthScale(ii,jj,s2,e2) ) )*norm


              myDGSEM % stressTensor % boundaryFlux(ii,jj,iEq,s2,e2) = -myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1)

            ENDDO
          ENDDO
        ENDDO

      ENDIF

    ENDDO

#endif

  END SUBROUTINE InternalFace_StressFlux_Fluid
!
  SUBROUTINE BoundaryFace_StressFlux_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1, &
                  1 )
    grid = dim3(myDGSEM % mesh % faces % nFaces,myDGSEM % stressTensor % nEquations,1)

    CALL BoundaryFace_StressFlux_CUDAKernel<<<grid, tBlock>>>( myDGSEM % state % boundarySolution_dev, &
                                                               myDGSEM % state % externalState_dev, &
                                                               myDGSEM % stressTensor % boundarySolution_dev, &
                                                               myDGSEM % stressTensor % externalState_dev, &
                                                               myDGSEM % static % boundarySolution_dev, &
                                                               myDGSEM % static % externalState_dev, &
                                                               myDGSEM % mesh % elements % boundaryLengthScale_dev, &
                                                               myDGSEM % mesh % elements % nHat_dev, &
                                                               myDGSEM % stressTensor % boundaryFlux_dev, &
                                                               myDGSEM % mesh % faces % elementIDs_dev, &
                                                               myDGSEM % mesh % faces % elementSides_dev, &
                                                               myDGSEM % mesh % faces % boundaryID_dev, &
                                                               myDGSEM % mesh % faces % iMap_dev, &
                                                               myDGSEM % mesh % faces % jMap_dev )

#else
    ! Local
    REAL(prec) :: norm, rhoIn, rhoOut
    INTEGER :: iEl, iFace
    INTEGER    :: i, j, k, m, iEq
    INTEGER    :: ii, jj, bID
    INTEGER    :: e1, s1, e2, s2

    DO iFace = 1, myDGSEM % mesh % faces % nFaces


      e1  = myDGSEM % mesh % faces % elementIDs(1,iFace)
      s1  = myDGSEM % mesh % faces % elementSides(1,iFace)
      e2  = myDGSEM % mesh % faces % elementIDs(2,iFace)
      s2  = ABS(myDGSEM % mesh % faces % elementSides(2,iFace))
      bID = myDGSEM % mesh % faces % boundaryID(iFace)

      IF( e2 < 0 )THEN

        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg


            ii = myDGSEM % mesh % faces % iMap(i,j,iFace)
            jj = myDGSEM % mesh % faces % jMap(i,j,iFace)

            norm = sqrt( myDGSEM % mesh % elements % nHat(1,i,j,s1,e1)**2 + &
                         myDGSEM % mesh % elements % nHat(2,i,j,s1,e1)**2 + &
                         myDGSEM % mesh % elements % nHat(3,i,j,s1,e1)**2 )


            DO iEq = 1, myDGSEM % state % nEquations-1

              IF( iEq == 4 )THEN

               
                  myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*( myDGSEM % stressTensor % boundarySolution(i,j,iEq,s1,e1)- &
                                                                                      myDGSEM % stressTensor % externalState(ii,jj,iEq,bID) )
                

              ELSE

                rhoOut = (myDGSEM % static % externalState(ii,jj,4,bID)+myDGSEM % state % externalState(ii,jj,4,bID) )
                rhoIn  = (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % boundarySolution(i,j,4,s1,e1) )

                

                  myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*( rhoIn*myDGSEM % stressTensor % boundarySolution(i,j,iEq,s1,e1)-&
                                                                                    rhoOut*myDGSEM % stressTensor % externalState(ii,jj,iEq,bID) )

               

              ENDIF

              myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) +&
                ( myDGSEM % params % viscosity*myDGSEM % state % externalState(ii,jj,iEq,bID)-&
                  myDGSEM % params % viscosity*myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) )/&
                ( myDGSEM % mesh % elements % boundaryLengthScale(i,j,s1,e1) )*norm


            ENDDO

          ENDDO
        ENDDO

      ENDIF

    ENDDO

#endif

  END SUBROUTINE BoundaryFace_StressFlux_Fluid
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
    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1 )
    grid = dim3(myDGSEM % mesh % elements % nElements,1,1)

    CALL EquationOfState_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                      myDGSEM % static % solution_dev )

#else
    ! Local
    INTEGER :: iEl, i, j, k

    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO k = 0, myDGSEM % params % polyDeg
        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg

            ! Pressure = rho*R*T
            ! And P' = P - P_static = R*(rho*T)'
            myDGSEM % state % solution(i,j,k,nEquations,iEl) = myDGSEM % state % solution(i,j,k,5,iEl)*myDGSEM % params % R

          ENDDO
        ENDDO
      ENDDO
    ENDDO
#endif

  END SUBROUTINE EquationOfState_Fluid
!
  SUBROUTINE CalculateStaticState_Fluid( myDGSEM )
#undef __FUNC__
#define __FUNC__ "CalculateStaticState"
    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM
    ! Local
    INTEGER    :: i, j, k, iEl,iEq
    REAL(prec) :: z, H, P0, Cp, T, T0, dTdz, P, rC, g, R, hCapRatio, fs
#ifdef HAVE_CUDA
    INTEGER :: istat
#endif

    INFO('Start')
    R    = myDGSEM % params % R
    Cp   = (R + myDGSEM % params % Cv)
    rC   = R/Cp
    hcapratio = myDGSEM % params % Cv/Cp
    g    = myDGSEM % params % g
    T0   = myDGSEM % params % T0

    ! /////////////////////  Build the Static/Background State ///////////////////////// !
    !
    ! The static state is constructed using the assumptions of an isentropic, neutrally
    ! stable (constant potential temperature) atmosphere.
    ! 

    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO iEq = 1, myDGSEM % state % nEquations
        DO k = 0, myDGSEM % params % polyDeg
          DO j = 0, myDGSEM % params % polyDeg
            DO i = 0, myDGSEM % params % polyDeg
              myDGSEM % state % solution(i,j,k,iEq,iEl)  = 0.0_prec
              myDGSEM % static % solution(i,j,k,iEq,iEl) = 0.0_prec
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO k = 0, myDGSEM % params % polyDeg
        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg

            z = myDGSEM % mesh % elements % x(i,j,k,3,iEl)


            fs = (1.0_prec - (1.0_prec - hCapRatio)*g*z/( R*T0 ) )**( 1.0_prec/(1.0_prec - hCapRatio) ) 
              
            ! Density
            myDGSEM % static % solution(i,j,k,4,iEl) = myDGSEM % params % rho0*( fs )**hCapRatio

            myDGSEM % static % solution(i,j,k,5,iEl) = myDGSEM % static % solution(i,j,k,4,iEl)*T0*( fs )**( R/Cp )
            myDGSEM % static % solution(i,j,k,nEquations,iEl) = myDGSEM % static % solution(i,j,k,5,iEl)*myDGSEM % params % R

          ENDDO
        ENDDO
      ENDDO
    ENDDO

#ifdef HAVE_CUDA
    myDGSEM % static % solution_dev = myDGSEM % static % solution
    istat = cudaDeviceSynchronize( )
#endif

    CALL myDGSEM % Update_FluidStatics_BCs( )
    INFO('End')

  END SUBROUTINE CalculateStaticState_Fluid
!
  SUBROUTINE IO_Fluid( myDGSEM )
    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM 
    ! Local
    CHARACTER(50) :: filename
    CHARACTER(13) :: timeStampString
    CHARACTER(4)  :: rankChar

      timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )
      WRITE(rankChar,'(I4.4)') myRank

      filename = "State."//timeStampString//".h5"
      CALL myDGSEM % Write_to_HDF5( filename )

#ifdef TECPLOT
      filename = "State."//rankChar//"."//timeStampString//".tec"
      CALL myDGSEM % WriteTecplot( filename )
#endif

#ifdef DIAGNOSTICS
      CALL myDGSEM % Diagnostics( )
      CALL myDGSEM % WriteDiagnostics( )
#endif

      
  END SUBROUTINE IO_Fluid
!
  SUBROUTINE WriteTecplot_Fluid( myDGSEM, filename )

    IMPLICIT NONE

    CLASS( Fluid ), INTENT(inout) :: myDGsem
    CHARACTER(*), INTENT(in)      :: filename
    !LOCAL
    REAL(prec)  :: x(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:3,1:myDGSEM % mesh % elements % nElements)
    REAL(prec)  :: sol(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:myDGSEM % state % nEquations,1:myDGSEM % mesh % elements % nElements)
    REAL(prec)  :: bsol(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:myDGSEM % state % nEquations, 1:myDGSEM % mesh % elements % nElements)
#ifdef HAVE_CUDA
    INTEGER :: istat
#endif
    INTEGER       :: i, j, k, iEl, iEq, fUnit
    CHARACTER(5)  :: zoneID
    CHARACTER(4)  :: rankChar
    REAL(prec)    :: hCapRatio, c, T, insitu, pottemp


#ifdef HAVE_CUDA
    CALL myDGSEM % state % UpdateHost( )
    CALL myDGSEM % static % UpdateHost( )
    istat = cudaDeviceSynchronize( )
#endif


    IF( myDGSEM % params % nPlot == myDGSEM % params % polyDeg )THEN
      sol = myDGSEM % state % solution
      bsol = myDGSEM % static % solution
      x    = myDGSEM % mesh % elements % x
    ELSE
      sol = ApplyInterpolationMatrix_3D_Lagrange( myDGSEM % dgStorage % interp, &
                                                  myDGSEM % state % solution, &
                                                  myDGSEM % state % nEquations, &
                                                  myDGSEM % mesh % elements % nElements )

      bsol = ApplyInterpolationMatrix_3D_Lagrange( myDGSEM % dgStorage % interp, myDGSEM % static % solution, &
                                                   myDGSEM % static % nEquations, &
                                                   myDGSEM % mesh % elements % nElements )

      x = ApplyInterpolationMatrix_3D_Lagrange( myDGSEM % dgStorage % interp, myDGSEM % mesh % elements % x, &
                                                3, &
                                                myDGSEM % mesh % elements % nElements )
    ENDIF


    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(filename), &
      FORM='formatted', &
      STATUS='replace')


    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "u", "v", "w", "rho", "rho-bar", "In-situ Temp.", "Pot. Temp.", "Tracer", "Pressure"'


    DO iEl = 1, myDGsem % mesh % elements % nElements

      WRITE(zoneID,'(I5.5)') myDGSEM % mesh % elements % elementID(iEl)
      WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',myDGSEM % params % nPlot+1,&
                                                 ', J=',myDGSEM % params % nPlot+1,&
                                                 ', K=',myDGSEM % params % nPlot+1,',F=POINT'

      DO k = 0, myDGSEM % params % nPlot
        DO j = 0, myDGSEM % params % nPlot
          DO i = 0, myDGSEM % params % nPlot

            insitu  = ( sol(i,j,k,5,iEl) + bsol(i,j,k,5,iEl) )/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) )
            pottemp = insitu*( (bsol(i,j,k,7,iEl) + sol(i,j,k,7,iEl))/myDGSEM % params % P0 )**(-myDGSEM % params % R/( myDGSEM % params % R + myDGSEM % params % Cv ) ) 


            WRITE(fUnit,'(17(E15.7,1x))') x(i,j,k,1,iEl), &
                                          x(i,j,k,2,iEl), &
                                          x(i,j,k,3,iEl), &
                                          sol(i,j,k,1,iEl)/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ), &
                                          sol(i,j,k,2,iEl)/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ), &
                                          sol(i,j,k,3,iEl)/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ), &
                                          sol(i,j,k,4,iEl), &
                                          bsol(i,j,k,4,iEl), &
                                          insitu, pottemp, &
                                          ( sol(i,j,k,6,iEl) + bsol(i,j,k,6,iEl) )/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ),  &
                                          sol(i,j,k,nEquations,iEl)

          ENDDO
        ENDDO
      ENDDO

    ENDDO

    CLOSE(UNIT=fUnit)

  END SUBROUTINE WriteTecplot_Fluid
!
  SUBROUTINE WriteDiagnostics_Fluid( myDGSEM )
    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM
    ! Local
    LOGICAL :: fileExists


    IF( myRank == 0 )THEN

      INQUIRE( file='Diagnostics.curve', EXIST=fileExists )

      OPEN( UNIT=NewUnit(diagUnits), &
        FILE='Diagnostics.curve', &
        FORM='FORMATTED', &
        ACCESS='APPEND' )

      IF( .NOT. fileExists )THEN
        WRITE(diagUnits,*) 'Time, Total Mass, Total Kinetic Energy, Total Potential Energy, Total heat content, Total volume'
      ENDIF

    ENDIF

    IF( myRank == 0 )THEN
      WRITE(diagUnits,'(5(E15.5,",",2x),E15.5)') myDGSEM % simulationTime, myDGSEM % mass, myDGSEM % KE, myDGSEM % PE, myDGSEM % heat, myDGSEM % volume
    ENDIF

    IF( myRank == 0 ) THEN
      CLOSE( UNIT=diagUnits )
    ENDIF

  END SUBROUTINE WriteDiagnostics_Fluid
!
  SUBROUTINE Diagnostics_Fluid( myDGSEM )

#undef __FUNC__
#define __FUNC__ "Diagnostics"

    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM
    ! Local
    INTEGER    :: iEl, i, j, k
    REAL(prec) :: volume, mass, KE, PE, heat
#ifdef HAVE_MPI
    INTEGER    :: mpiErr
#endif
    CHARACTER(40) :: diagString

    INFO('Start')

    volume = 0.0_prec
    mass   = 0.0_prec
    KE     = 0.0_prec
    PE     = 0.0_prec
    heat   = 0.0_prec

    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO k = 0, myDGSEM % params % polyDeg
        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg

            volume = volume + myDGSEM % mesh % elements % J(i,j,k,iEl)*&
              myDGSEM % dgStorage % quadratureWeights(i)*&
              myDGSEM % dgStorage % quadratureWeights(j)*&
              myDGSEM % dgStorage % quadratureWeights(k)

            mass = mass + ( myDGSEM % state % solution(i,j,k,4,iEl)+&
              myDGSEM % static % solution(i,j,k,4,iEl) )*&
              myDGSEM % mesh % elements % J(i,j,k,iEl)*&
              myDGSEM % dgStorage % quadratureWeights(i)*&
              myDGSEM % dgStorage % quadratureWeights(j)*&
              myDGSEM % dgStorage % quadratureWeights(k)

            KE   = KE + ( myDGSEM % state % solution(i,j,k,1,iEl)**2 +&
              myDGSEM % state % solution(i,j,k,2,iEl)**2 +&
              myDGSEM % state % solution(i,j,k,3,iEl)**2 )/&
              ( myDGSEM % state % solution(i,j,k,4,iEl)+&
              myDGSEM % static % solution(i,j,k,4,iEl) )*&
              myDGSEM % mesh % elements % J(i,j,k,iEl)*&
              myDGSEM % dgStorage % quadratureWeights(i)*&
              myDGSEM % dgStorage % quadratureWeights(j)*&
              myDGSEM % dgStorage % quadratureWeights(k)

            PE   = PE - myDGSEM % state % solution(i,j,k,4,iEl)*&
              myDGSEM % params % g*&
              myDGSEM % mesh % elements % x(i,j,k,3,iEl)*&
              myDGSEM % mesh % elements % J(i,j,k,iEl)*&
              myDGSEM % dgStorage % quadratureWeights(i)*&
              myDGSEM % dgStorage % quadratureWeights(j)*&
              myDGSEM % dgStorage % quadratureWeights(k)

            heat = heat + ( myDGSEM % static % solution(i,j,k,5,iEl) + &
              myDGSEM % state % solution(i,j,k,5,iEl) )/&
              ( myDGSEM % static % solution(i,j,k,4,iEl) +&
              myDGSEM % state % solution(i,j,k,4,iEl) )*&
              myDGSEM % mesh % elements % J(i,j,k,iEl)*&
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
    CALL MPI_ALLREDUCE( volume, myDGSEM % volume, 1, mpiPrec, MPI_SUM, mpiComm, mpiErr )
    CALL MPI_ALLREDUCE( mass, myDGSEM % mass, 1, mpiPrec, MPI_SUM, mpiComm, mpiErr )
    CALL MPI_ALLREDUCE( KE, myDGSEM % KE, 1, mpiPrec, MPI_SUM, mpiComm, mpiErr )
    CALL MPI_ALLREDUCE( PE, myDGSEM % PE, 1, mpiPrec, MPI_SUM, mpiComm, mpiErr )
    CALL MPI_ALLREDUCE( heat, myDGSEM % heat, 1, mpiPrec, MPI_SUM, mpiComm, mpiErr )
#endif

    IF( myRank == 0 )THEN
  
      WRITE(diagString, '(E15.5,", ",E15.5)') myDGSEM % simulationTime, myDGSEM % volume
      INFO('Time (s), Volume (m^3) : '//diagString )

      WRITE(diagString, '(E15.5,", ",E15.5)') myDGSEM % simulationTime, myDGSEM % mass
      INFO('Time (s), Mass (kg) : '//diagString )

      WRITE(diagString, '(E15.5,", ",E15.5)') myDGSEM % simulationTime, myDGSEM % KE
      INFO('Time (s), Kinetic Energy (N) : '//diagString )

      WRITE(diagString, '(E15.5,", ",E15.5)') myDGSEM % simulationTime, myDGSEM % PE
      INFO('Time (s), Potential Energy (N) : '//diagString )

      WRITE(diagString, '(E15.5,", ",E15.5)') myDGSEM % simulationTime, myDGSEM % heat
      INFO('Time (s), Heat (Km^3) : '//diagString )

    ENDIF
    IF( IsNaN( myDGSEM % KE ) .OR. IsInf( myDGSEM % KE/ myDGSEM % volume) )THEN

      INFO('Model bust!')
      INFO('Consider reducing the time step or increasing dissipation.')

      CALL myDGSEM % WriteDiagnostics( )
      CALL myDGSEM % Trash( )

      STOP
    ENDIF
    INFO('End')

  END SUBROUTINE Diagnostics_Fluid
!
  SUBROUTINE Add_Variable_to_HDF5( file_id, variable_name, variable, filespace, memspace, plist_id, dimensions, elementIDs, N, nElements, error )
    IMPLICIT NONE
    INTEGER(HID_T), INTENT(in)   :: file_id
    CHARACTER(*), INTENT(in)     :: variable_name
    INTEGER(HID_T), INTENT(in)   :: filespace
    INTEGER(HID_T), INTENT(in)   :: memspace
    INTEGER(HID_T), INTENT(in)   :: plist_id
    INTEGER, INTENT(in)          :: N, nElements
    REAL(prec), INTENT(in)       :: variable(0:N,0:N,0:N,1:nElements)
    INTEGER(HSIZE_T), INTENT(in) :: dimensions(1:4)
    INTEGER, INTENT(in)          :: elementIDs(1:nElements)
    INTEGER, INTENT(out)         :: error
    ! Local
    INTEGER(HID_T) :: dataset_id 
    INTEGER        :: iEl, elID
    INTEGER(HSIZE_T)        :: starts(1:4), counts(1:4), strides(1:4)
    



#ifdef DOUBLE_PRECISION
    CALL h5dcreate_f( file_id, TRIM(variable_name), H5T_IEEE_F64LE, filespace, dataset_id, error, plist_id)
#else
    CALL h5dcreate_f( file_id, TRIM(variable_name), H5T_IEEE_F32LE, filespace, dataset_id, error, plist_id)
#endif

    DO iEl = 1, nElements
      elID = elementIDs(iEl)
      starts = (/ 0, 0, 0, elID-1 /)
      counts = (/ 1, 1, 1, 1 /)
      strides = (/ 1, 1, 1, 1 /)

      CALL h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, starts, counts, error, strides, dimensions )

#ifdef DOUBLE_PRECISION
      CALL h5dwrite_f( dataset_id, H5T_IEEE_F64LE, &
                       variable(0:N,0:N,0:N,iEl:iEl), &
                       dimensions, error, memspace, filespace )
#else
      CALL h5dwrite_f( dataset_id, H5T_IEEE_F32LE, &
                       variable(0:N,0:N,0:N,iEl:iEl), &
                       dimensions, error, memspace, filespace )
#endif

    ENDDO

#else

#ifdef DOUBLE_PRECISION
    CALL h5dcreate_f( file_id, TRIM(variable_name), &
                       H5T_IEEE_F64LE, memspace, dataset_id, error)
    CALL h5dwrite_f( dataset_id, H5T_IEEE_F64LE, &
                     variable, dimensions, error)
#else
    CALL h5dcreate_f( file_id, TRIM(variable_name), &
                       H5T_IEEE_F32LE, memspace, dataset_id, error)
    CALL h5dwrite_f( dataset_id, H5T_IEEE_F32LE, &
                     variable, dimensions, error)
#endif

#endif

!    CALL h5sclose_f( filespace, error)
    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Add_Variable_to_HDF5
!
  SUBROUTINE Get_Variable_from_HDF5( file_id, variable_name, variable, filespace, memspace, plist_id, dimensions, elementIDs, N, nElements, error )
#undef __FUNC__
#define __FUNC__ "Get_Variable_from_HDF5"
    IMPLICIT NONE
    INTEGER(HID_T), INTENT(in)   :: file_id
    CHARACTER(*), INTENT(in)     :: variable_name
    INTEGER(HID_T), INTENT(inout):: filespace
    INTEGER(HID_T), INTENT(in)   :: memspace
    INTEGER(HID_T), INTENT(in)   :: plist_id
    INTEGER, INTENT(in)          :: N, nElements
    REAL(prec), INTENT(inout)    :: variable(0:N,0:N,0:N,1:nElements)
    INTEGER(HSIZE_T), INTENT(in) :: dimensions(1:4)
    INTEGER, INTENT(in)          :: elementIDs(1:nElements)
    INTEGER, INTENT(out)         :: error
    ! Local
    INTEGER(HID_T) :: dataset_id 
    INTEGER        :: iEl, elID
    INTEGER(HSIZE_T) :: starts(1:4), counts(1:4), strides(1:4)

    CALL h5dopen_f(file_id, TRIM(variable_name), dataset_id, error)
    CALL h5dget_space_f( dataset_id, filespace, error )

#ifdef HAVE_MPI

    DO iEl = 1, nElements
      elID = elementIDs(iEl)
      starts = (/ 0, 0, 0, elID-1 /)
      counts = (/ 1, 1, 1, 1 /)
      strides = (/ 1, 1, 1, 1 /)

      CALL h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, starts, counts, error, strides, dimensions )
#ifdef DOUBLE_PRECISION
      CALL h5dread_f( dataset_id, H5T_IEEE_F64LE, &
                       variable(0:N,0:N,0:N,iEl:iEl), &
                       dimensions, error, memspace, filespace)
#else
      CALL h5dread_f( dataset_id, H5T_IEEE_F32LE, &
                       variable(0:N,0:N,0:N,iEl:iEl), &
                       dimensions, error, memspace, filespace)
#endif
    IF( error /= 0 )THEN
      INFO('h5dread_f error :'//TRIM(variable_name) )
    ENDIF

    ENDDO

#else

#ifdef DOUBLE_PRECISION
    CALL h5dread_f( dataset_id, H5T_IEEE_F64LE, &
                    variable(0:N,0:N,0:N,1:nElements), dimensions, error)
#else
    CALL h5dread_f( dataset_id, H5T_IEEE_F32LE, &
                    variable(0:N,0:N,0:N,1:nElements), dimensions, error)
    IF( error /= 0 )THEN
      INFO('h5dread_f error :'//TRIM(variable_name) )
    ENDIF
#endif

#endif

    CALL h5sclose_f( filespace, error)
    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Get_Variable_from_HDF5
!
  SUBROUTINE Write_to_HDF5( myDGSEM, filename )
#undef __FUNC__
#define __FUNC__ "Write_to_HDF5"
    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM
    CHARACTER(*), INTENT(in)      :: filename
    ! Local
    CHARACTER(13)    :: timeStampString
    CHARACTER(10)    :: zoneID
    INTEGER          :: iEl, N, rank, m_rank, error, istat, nEl, elID
    INTEGER(HSIZE_T) :: dimensions(1:4), global_dimensions(1:4)
    INTEGER(HSIZE_T) :: starts(1:4), counts(1:4), strides(1:4)
    INTEGER(HID_T)   :: file_id, memspace, dataset_id, filespace
    INTEGER(HID_T)   :: group_id
    INTEGER(HID_T)   :: conditions_group_id, conditions_element_group_id, plist_id
    INTEGER(HID_T)   :: create_plist_id, access_plist_id, transfer_plist_id

    INFO('Start')
#ifdef HAVE_CUDA
    CALL myDGSEM % state % UpdateHost( )
    CALL myDGSEM % static % UpdateHost( )
    istat = cudaDeviceSynchronize( )
#endif

    timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )
    INFO('Writing output file : '//TRIM(filename) )


    N = myDGSEM % params % polyDeg
#ifdef HAVE_MPI
    CALL MPI_BARRIER( mpiComm, error )
    CALL MPI_ALLREDUCE( myDGSEM % mesh % elements % nElements, nEl, 1, MPI_INTEGER, MPI_SUM, mpiComm, error )
    rank = 4
    ! Local Dimensions
    dimensions = (/ N+1, N+1, N+1, 1 /)
    global_dimensions = (/ N+1, N+1, N+1, nEl /)
#else
    nEl = myDGSEM % mesh % elements % nElements
    rank = 4
    dimensions = (/ N+1, N+1, N+1, nEl /)
    global_dimensions = (/ N+1, N+1, N+1, nEl /)
#endif


    CALL h5open_f(error)  
  
#ifdef HAVE_MPI

    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, mpiComm, MPI_INFO_NULL, error)

    ! Create a new file using default properties.
    CALL h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
    CALL h5pclose_f(plist_id, error)

    CALL h5screate_simple_f(rank, global_dimensions, filespace, error)
    CALL h5screate_simple_f(rank, dimensions, memspace, error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    CALL h5pset_chunk_f(plist_id, rank, dimensions, error)

#else

    CALL h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, file_id, error)
    CALL h5screate_simple_f(rank, dimensions, memspace, error)

#endif

    ! Create groups 
    CALL h5gcreate_f( file_id, "/model_output", group_id, error )
    CALL h5gclose_f( group_id, error )

    CALL h5gcreate_f( file_id, "/static", group_id, error )
    CALL h5gclose_f( group_id, error )

    CALL h5gcreate_f( file_id, "/model_conditions", group_id, error )
    CALL h5gclose_f( group_id, error )

    ! Create variables
    CALL Add_Variable_to_HDF5( file_id, "/model_output/x_momentum", myDGSEM % state % solution(:,:,:,1,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/model_output/y_momentum", myDGSEM % state % solution(:,:,:,2,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/model_output/z_momentum", myDGSEM % state % solution(:,:,:,3,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/model_output/density", myDGSEM % state % solution(:,:,:,4,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/model_output/density_weighted_temperature", myDGSEM % state % solution(:,:,:,5,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/model_output/density_weighted_tracer", myDGSEM % state % solution(:,:,:,6,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/model_output/pressure", myDGSEM % state % solution(:,:,:,7,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/static/x_momentum", myDGSEM % static % solution(:,:,:,1,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/static/y_momentum", myDGSEM % static % solution(:,:,:,2,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/static/z_momentum", myDGSEM % static % solution(:,:,:,3,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/static/density", myDGSEM % static % solution(:,:,:,4,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/static/density_weighted_temperature", myDGSEM % static % solution(:,:,:,5,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/static/density_weighted_tracer", myDGSEM % static % solution(:,:,:,6,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/static/pressure", myDGSEM % static % solution(:,:,:,7,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Add_Variable_to_HDF5( file_id, "/model_conditions/drag", myDGSEM % sourceTerms % drag, &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )


    ! Clean up property list and dataspace handles
#ifdef HAVE_MPI
    CALL h5pclose_f( plist_id, error )
    CALL h5sclose_f( filespace, error )
#endif
    CALL h5sclose_f( memspace, error )

    ! Close the file
    CALL h5fclose_f( file_id, error )
    ! Close access to HDF5
    CALL h5close_f( error )
    INFO('Finished writing output file : '//TRIM(filename))
    INFO('End')

  END SUBROUTINE Write_to_HDF5

  SUBROUTINE Read_from_HDF5( myDGSEM, itExists, filename )
#undef __FUNC__
#define __FUNC__ "Read_from_HDF5"
    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout)      :: myDGSEM
    CHARACTER(*), OPTIONAL, INTENT(in) :: filename
    LOGICAL, INTENT(out)               :: itExists
    ! Local
    CHARACTER(100)   :: fname, groupname
    CHARACTER(13)    :: timeStampString
    CHARACTER(10)    :: zoneID
    INTEGER          :: iEl, N, rank, m_rank, error, nEl, elID, istat
    INTEGER(HSIZE_T) :: dimensions(1:4), global_dimensions(1:4)
    INTEGER(HSIZE_T) :: starts(1:4), counts(1:4), strides(1:4)
    INTEGER(HID_T)   :: file_id, memspace, filespace, dataset_id
    INTEGER(HID_T)   :: model_group_id, static_group_id, static_element_group_id, element_group_id
    INTEGER(HID_T)   :: conditions_group_id, conditions_element_group_id
    INTEGER(HID_T)   :: plist_id

    INFO('Start')
    IF( PRESENT( filename ) )THEN
      fname = TRIM(filename)
    ELSE
      timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )
      fname = "State."//timeStampString//".h5"
    ENDIF

    INQUIRE( FILE=TRIM(fname), EXIST = itExists )
     
    IF( .NOT. itExists ) THEN
      RETURN
    ENDIF

    INFO('Reading '//TRIM(fname))

    N = myDGSEM % params % polyDeg
#ifdef HAVE_MPI
    CALL MPI_ALLREDUCE( myDGSEM % mesh % elements % nElements, nEl, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, error )
    rank = 4
    dimensions = (/ N+1, N+1, N+1, 1 /)
    global_dimensions = (/ N+1, N+1, N+1, nEl /)
#else
    nEl = myDGSEM % mesh % elements % nElements
    rank = 4
    dimensions = (/ N+1, N+1, N+1, nEl /)
    global_dimensions = (/ N+1, N+1, N+1, nEl /)

#endif

    CALL h5open_f(error)  

#ifdef HAVE_MPI
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, mpiComm, MPI_INFO_NULL, error)

    CALL h5fopen_f(TRIM(fname), H5F_ACC_RDWR_F, file_id, error, access_prp=plist_id)
    CALL h5pclose_f(plist_id,error)


#else

    CALL h5fopen_f(TRIM(fname), H5F_ACC_RDWR_F, file_id, error)
    IF( error /= 0 ) STOP

#endif

    CALL h5screate_simple_f(rank, dimensions, memspace, error)

    groupname = "/model_output"
    CALL h5gopen_f( file_id, TRIM(groupname), model_group_id, error )
    IF( error /= 0 ) STOP

    groupname = "/static"
    CALL h5gopen_f( file_id, TRIM(groupname), static_group_id, error )
    IF( error /= 0 ) STOP

    groupname = "/model_conditions"
    CALL h5gopen_f( file_id, TRIM(groupname), conditions_group_id, error )
    IF( error /= 0 ) STOP


    CALL Get_Variable_from_HDF5( file_id, "/model_output/x_momentum", myDGSEM % state % solution(:,:,:,1,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/model_output/y_momentum", myDGSEM % state % solution(:,:,:,2,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/model_output/z_momentum", myDGSEM % state % solution(:,:,:,3,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/model_output/density", myDGSEM % state % solution(:,:,:,4,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/model_output/density_weighted_temperature", myDGSEM % state % solution(:,:,:,5,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/model_output/density_weighted_tracer", myDGSEM % state % solution(:,:,:,6,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/model_output/pressure", myDGSEM % state % solution(:,:,:,7,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/static/x_momentum", myDGSEM % static % solution(:,:,:,1,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/static/y_momentum", myDGSEM % static % solution(:,:,:,2,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/static/z_momentum", myDGSEM % static % solution(:,:,:,3,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/static/density", myDGSEM % static % solution(:,:,:,4,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/static/density_weighted_temperature", myDGSEM % static % solution(:,:,:,5,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/static/density_weighted_tracer", myDGSEM % static % solution(:,:,:,6,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/static/pressure", myDGSEM % static % solution(:,:,:,7,:), &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )

    CALL Get_Variable_from_HDF5( file_id, "/model_conditions/drag", myDGSEM % sourceTerms % drag, &
                               filespace, memspace, plist_id, dimensions, &
                               myDGSEM % mesh % elements % elementID, N, &
                               myDGSEM % mesh % elements % nElements, error )


    IF( error /= 0 ) STOP

#ifdef HAVE_MPI
    CALL MPI_BARRIER( mpiComm, error )
#endif
    CALL h5gclose_f( model_group_id, error )
    CALL h5gclose_f( static_group_id, error )
    CALL h5gclose_f( conditions_group_id, error )
    CALL h5sclose_f( memspace, error )
    CALL h5fclose_f( file_id, error )
    CALL h5close_f( error )


    CALL myDGSEM % SetPrescribedState( )
    CALL myDGSEM % Update_FluidStatics_BCs( )
#ifdef HAVE_CUDA
    CALL myDGSEM % state % UpdateDevice( )
    CALL myDGSEM % static % UpdateDevice( )
    CALL myDGSEM % sourceTerms % UpdateDevice( )
#endif

    INFO('End')

  END SUBROUTINE Read_from_HDF5

  SUBROUTINE Update_FluidStatics_BCs( myDGSEM )
#undef __FUNC__
#define __FUNC__ "Update_FluidStatics_BCs"
    IMPLICIT NONE
    CLASS( Fluid ), INTENT (inout) :: myDGSEM 
    ! Local
    INTEGER :: bID, iFace, e1, s1, e2, i, j, iEq, error
#ifdef HAVE_CUDA
    INTEGER :: iStat
#endif

     INFO('Start')
     myDGSEM % static % boundarySolution = CalculateFunctionsAtBoundaries_3D_NodalDG( myDGSEM % dgStorage, &
                                                                                      myDGSEM % static % solution, &
                                                                                      myDGSEM % static % nEquations, &
                                                                                      myDGSEM % static % nElements )
    ! Update the static external states
    DO iFace = 1, myDGSEM % mesh % faces % nFaces
  
      e1  = myDGSEM % mesh % faces % elementIDs(1,iFace)
      s1  = myDGSEM % mesh % faces % elementSides(1,iFace)
      e2  = myDGSEM % mesh % faces % elementIDs(2,iFace)
      bID = myDGSEM % mesh % faces % boundaryID(iFace)
  
      IF( e2 < 0 )THEN

        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg
            DO iEq = 1, myDGSEM % static % nEquations
              myDGSEM % static % externalState(i,j,iEq,bID) = myDGSEM % static % boundarySolution(i,j,iEq,s1,e1)
            ENDDO
          ENDDO
        ENDDO

      ENDIF
  
    ENDDO

#ifdef HAVE_MPI
    CALL myDGSEM % mpiStateHandler % MPI_Exchange( myDGSEM % static, myDGSEM % mesh )
    CALL myDGSEM % mpiStateHandler % Finalize_MPI_Exchange( )

#endif

#ifdef HAVE_CUDA
    myDGSEM % static % externalState_dev = myDGSEM % static % externalState
    myDGSEM % static % boundarySolution_dev = myDGSEM % static % boundarySolution
    istat = cudaDeviceSynchronize( )
#endif

     
     INFO('End')

  END SUBROUTINE Update_FluidStatics_BCs

!
#ifdef HAVE_CUDA
! ============================================================================================================================ !
!------------------------------------------- CUDA Kernels Below -------------------------------------------------------------- !
! ============================================================================================================================ !
  ATTRIBUTES(Global) SUBROUTINE UpdateG3D_CUDAKernel( G3D, a, g, solution, fluxDivergence, source, diffusiveFluxDivergence )
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(inout) :: G3D(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: a, g
    REAL(prec), DEVICE, INTENT(inout) :: solution(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: source(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: fluxDivergence(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: diffusiveFluxDivergence(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:nEl_dev)
    ! Local
    INTEGER :: i, j, k, iEq, iEl
  
    iEl = blockIDx % x
    iEq = blockIDx % y
  
    i = threadIdx % x - 1
    j = threadIdx % y - 1
    k = threadIdx % z - 1
  
    G3D(i,j,k,iEq,iEl)      = a*G3D(i,j,k,iEq,iEl) - fluxDivergence(i,j,k,iEq,iEl) + diffusiveFluxDivergence(i,j,k,iEq,iEl) + source(i,j,k,iEq,iEl)
  
    solution(i,j,k,iEq,iEl) = solution(i,j,k,iEq,iEl) + dt_dev*g*G3D(i,j,k,iEq,iEl)

  
  END SUBROUTINE UpdateG3D_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE Update_FluidState_BCs_CUDAKernel( element_to_blockID, boundaryID, elementIDs, elementSides, externalState, stateBsols, prescribedState, nHat )
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)       :: element_to_blockID(1:nGlobalElements_dev)
    INTEGER, DEVICE, INTENT(in)       :: boundaryID(1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)       :: elementIDs(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)       :: elementSides(1:2,1:nFaces_dev)
    REAL(prec), DEVICE, INTENT(inout) :: externalState(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)    :: stateBsols(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: prescribedState(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)    :: nhat(1:3,0:polyDeg_dev,0:polyDeg_dev,1:6,1:nEl_dev)
     ! Local
    INTEGER    :: iEl, iFace, bID, i, j, k, iEq
    INTEGER    :: p2
    INTEGER    :: e1, e2, s1, s2
    REAL(prec) :: norm, un, ut, us, speed
    REAL(prec) :: nx, ny, nz
    REAL(prec) :: sx, sy, sz
    REAL(prec) :: tx, ty, tz
    
    iFace = blockIdx % x
    i     = threadIdx % x-1
    j     = threadIdx % y-1
    bID   = boundaryID(iFace)
    
    e1    = elementIDs(1,iFace)
    s1    = elementSides(1,iFace)
    e2    = elementIDs(2,iFace)
   
    IF( e2 == PRESCRIBED )THEN

      DO iEq = 1, nEq_dev
        externalState(i,j,iEq,bID) = prescribedState(i,j,iEq,bID)
      ENDDO
    
    ELSEIF( e2 == RADIATION )THEN

      DO iEq = 1, nEq_dev
        externalState(i,j,iEq,bID) = 0.0_prec
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
    
      ! tangent (built by performing 90 deg rotation in y - IF zero, performs rotation in x)
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
    
      externalState(i,j,1,bID) = -nx*un + us*sx + ut*tx ! u
      externalState(i,j,2,bID) = -ny*un + us*sy + ut*ty ! v
      externalState(i,j,3,bID) = -nz*un + us*sz + ut*tz ! w
      externalState(i,j,4,bID) =  stateBsols(i,j,4,s1,e1) ! rho
      externalState(i,j,5,bID) =  stateBsols(i,j,5,s1,e1) ! potential temperature
      externalState(i,j,6,bID) =  stateBsols(i,j,6,s1,e1) ! tracer
      externalState(i,j,nEq_dev,bID) =  stateBsols(i,j,nEq_dev,s1,e1) ! P
    
    ELSEIF( e2 == INFLOW )THEN
    
      ! normal
      nx = nHat(1,i,j,s1,e1) !**
      ny = nHat(2,i,j,s1,e1)
      nz = nHat(3,i,j,s1,e1)
      norm = sqrt( nx*nx + ny*ny + nz*nz )
      nx = nx/norm
      ny = ny/norm
      nz = nz/norm
    
      ! tangent (built by performing 90 deg rotation in y - IF zero, performs rotation in x)
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
           stateBsols(i,j,3,s1,e1)*nz + &
           prescribedState(i,j,1,bID)

      us = stateBsols(i,j,1,s1,e1)*sx    + &
        stateBsols(i,j,2,s1,e1)*sy    + &
        stateBsols(i,j,3,s1,e1)*sz
      ut = stateBsols(i,j,1,s1,e1)*tx  + &
        stateBsols(i,j,2,s1,e1)*ty  + &
        stateBsols(i,j,3,s1,e1)*tz
    
      externalState(i,j,1,bID) = -nx*un + us*sx + ut*tx ! u
      externalState(i,j,2,bID) = -ny*un + us*sy + ut*ty ! v
      externalState(i,j,3,bID) = -nz*un + us*sz + ut*tz ! w
      externalState(i,j,4,bID) =  prescribedState(i,j,4,bID) ! rho
      externalState(i,j,5,bID) =  prescribedState(i,j,5,bID) ! potential temperature
      externalState(i,j,6,bID) =  prescribedState(i,j,6,bID) ! tracer
      externalState(i,j,nEq_dev,bID) =  prescribedState(i,j,nEq_dev,bID) ! P
    
    ENDIF
    
  END SUBROUTINE Update_FluidState_BCs_CUDAKernel
!
ATTRIBUTES(Global) SUBROUTINE InternalFace_StateFlux_CUDAKernel( elementIDs, elementSides, boundaryID, iMap, jMap, &
                                                                 nHat, boundarySolution, boundarySolution_static, &
                                                                 boundaryFlux, stressFlux )
   IMPLICIT NONE
   INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: boundaryID(1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: iMap(0:polydeg_dev,0:polydeg_dev,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: jMap(0:polydeg_dev,0:polydeg_dev,1:nFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:polydeg_dev,0:polydeg_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundarySolution(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundarySolution_static(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(inout) :: boundaryFlux(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(inout) :: stressFlux(1:3,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   ! Local

   INTEGER    :: iEl, iFace
   INTEGER    :: i, j, k, iEq
   INTEGER    :: ii, jj
   INTEGER    :: e1, s1, e2, s2, bID
   REAL(prec) :: uOut, uIn, cIn, cOut, norm
   REAL(prec) :: aS(1:6)
   REAL(prec) :: fac, jump, T
   

      iFace = blockIdx % x
      j     = threadIdx % y - 1
      i     = threadIdx % x - 1

      e1 = elementIDs(1,iFace)
      e2 = elementIDs(2,iFace)
      s1 = elementSides(1,iFace)
      s2 = ABS(elementSides(2,iFace))
      bID = boundaryID(iFace)

      IF( bID == 0 )THEN

        ii = iMap(i,j,iFace)
        jj = jMap(i,j,iFace)

        norm = sqrt( nHat(1,i,j,s1,e1)*nHat(1,i,j,s1,e1) + &
                     nHat(2,i,j,s1,e1)*nHat(2,i,j,s1,e1) + &
                     nHat(3,i,j,s1,e1)*nHat(3,i,j,s1,e1) )
        
        ! Sound speed estimate for the external and internal states
        T =   (boundarySolution_static(ii,jj,5,s2,e2) + boundarySolution(ii,jj,5,s2,e2))/&
                 (boundarySolution(ii,jj,4,s2,e2)+boundarySolution_static(ii,jj,4,s2,e2) )
                 
        ! Sound speed estimate for the external and internal states
        cOut = sqrt( R_dev*T )

        T =   (boundarySolution_static(i,j,5,s1,e1) + boundarySolution(i,j,5,s1,e1))/&
                (boundarySolution(i,j,4,s1,e1)+boundarySolution_static(i,j,4,s1,e1) )        

        cIn = sqrt( R_dev*T )
                     
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
          aS(k) = aS(k) + (boundarySolution(i,j,7,s1,e1) + boundarySolution(ii,jj,7,s2,e2))*nHat(k,i,j,s1,e1)/norm
        ENDDO    

        DO iEq = 1, nEq_dev-1

           jump = boundarySolution(ii,jj,iEq,s2,e2)-boundarySolution(i,j,iEq,s1,e1)
           boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*( aS(iEq) - fac*jump )*norm
           boundaryFlux(ii,jj,iEq,s2,e2) = -boundaryFlux(i,j,iEq,s1,e1)

           IF( iEq == 4 )THEN
              DO k = 1, 3
                 ! Calculate the LDG flux for the stress tensor.
                 stressFlux(k,i,j,iEq,s1,e1) = 0.5_prec*( boundarySolution(i,j,iEq,s1,e1) +&
                                                          boundarySolution(ii,jj,iEq,s2,e2))*nHat(k,i,j,s1,e1)
                                                                              
                 stressFlux(k,ii,jj,iEq,s2,e2) = -stressFlux(k,i,j,iEq,s1,e1)
              ENDDO

           ELSE
              DO k = 1, 3
                 ! Calculate the LDG flux for the stress tensor.
                 stressFlux(k,i,j,iEq,s1,e1) = 0.5_prec*( boundarySolution(i,j,iEq,s1,e1)/&
                                                        (boundarySolution(i,j,4,s1,e1)+&
                                                         boundarySolution_static(i,j,4,s1,e1)) +&
                                                        boundarySolution(ii,jj,iEq,s2,e2)/&
                                                        (boundarySolution(ii,jj,4,s2,e2)+&
                                                         boundarySolution_static(ii,jj,4,s2,e2)) )*nHat(k,i,j,s1,e1)
                                                                              
                 stressFlux(k,ii,jj,iEq,s2,e2) = -stressFlux(k,i,j,iEq,s1,e1)

              ENDDO
           ENDIF
           
        ENDDO
                         
      ENDIF 


 END SUBROUTINE InternalFace_StateFlux_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE BoundaryFace_StateFlux_CUDAKernel( elementIDs, elementSides, boundaryID, iMap, jMap, &
                                                                 nHat, boundarySolution, boundarySolution_static, &
                                                                 externalState, externalStatic, boundaryFlux, stressFlux )

   IMPLICIT NONE
   INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: boundaryID(1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: iMap(0:polydeg_dev,0:polydeg_dev,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: jMap(0:polydeg_dev,0:polydeg_dev,1:nFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:polydeg_dev,0:polydeg_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundarySolution(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: boundarySolution_static(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(in)  :: externalState(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(in)  :: externalStatic(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
   REAL(prec), DEVICE, INTENT(inout) :: boundaryFlux(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(inout) :: stressFlux(1:3,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   ! Local
   INTEGER    :: iEl, iFace
   INTEGER    :: i, j, k, iEq
   INTEGER    :: ii, jj, bID
   INTEGER    :: e1, s1, e2, s2
   REAL(prec) :: uOut, uIn, cIn, cOut, norm, T
   REAL(prec) :: jump(1:6), aS(1:6)
   REAL(prec) :: fac


      iFace = blockIdx % x
      j     = threadIdx % y - 1
      i     = threadIdx % x -1
     
      e1 = elementIDs(1,iFace)
      s1 = elementSides(1,iFace)
      e2 = elementIDs(2,iFace)
      s2 = ABS(elementSides(2,iFace))
      bID  = boundaryID(iFace)

      ii = iMap(i,j,iFace)
      jj = jMap(i,j,iFace)
      
      norm = sqrt( nHat(1,i,j,s1,e1)*nHat(1,i,j,s1,e1) + &
                   nHat(2,i,j,s1,e1)*nHat(2,i,j,s1,e1) + &
                   nHat(3,i,j,s1,e1)*nHat(3,i,j,s1,e1) )

      
      IF( e2 < 0 )THEN
      
         
         DO iEq = 1, nEq_dev-1              
         jump(iEq)  = externalState(ii,jj,iEq,bID)-boundarySolution(i,j,iEq,s1,e1)
         ENDDO
        
         T =   (externalStatic(ii,jj,5,bID) + externalState(ii,jj,5,bID))/&
                 (externalState(ii,jj,4,bID)+externalStatic(ii,jj,4,bID) )

         cOut = sqrt( R_dev*T )
         
         T =   (boundarySolution_static(i,j,5,s1,e1) + boundarySolution(i,j,5,s1,e1))/&
                 (boundarySolution(i,j,4,s1,e1)+boundarySolution_static(i,j,4,s1,e1) )  
                          
         cIn = sqrt( R_dev*T )
                      
         ! External normal velocity component
         uOut = ( externalState(ii,jj,1,bID)*nHat(1,i,j,s1,e1)/norm + &
                  externalState(ii,jj,2,bID)*nHat(2,i,j,s1,e1)/norm + &
                  externalState(ii,jj,3,bID)*nHat(3,i,j,s1,e1)/norm )/& 
                ( externalState(ii,jj,4,bID) + externalStatic(ii,jj,4,bID) )
         ! Internal normal velocity component
         uIn  = ( boundarySolution(i,j,1,s1,e1)*nHat(1,i,j,s1,e1)/norm + &
                  boundarySolution(i,j,2,s1,e1)*nHat(2,i,j,s1,e1)/norm + &
                  boundarySolution(i,j,3,s1,e1)*nHat(3,i,j,s1,e1)/norm )/& 
                ( boundarySolution(i,j,4,s1,e1) + boundarySolution_static(i,j,4,s1,e1) )


         fac = max( abs(uIn+cIn), abs(uIn-cIn), abs(uOut+cOut), abs(uOut-cOut) )

         DO iEq = 1, nEq_dev-1
               aS(iEq) = uIn*( boundarySolution(i,j,iEq,s1,e1) + boundarySolution_static(i,j,iEq,s1,e1) ) +&
                        uOut*( externalState(ii,jj,iEq,bID) + externalStatic(i,j,iEq,bID) )
         ENDDO
         
         ! Pressure !
         DO k = 1, 3         
         aS(k) = aS(k) + (boundarySolution(i,j,7,s1,e1)+externalState(ii,jj,7,bID))*nHat(k,i,j,s1,e1)/norm
         ENDDO
         
                 
         DO iEq = 1, nEq_dev-1

            boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*( aS(iEq) - fac*jump(iEq) )*norm

            IF( iEq == 4 )THEN
               DO k = 1, 3
                  ! Calculate the Bassi-Rebay flux for the stress tensor.
                  stressFlux(k,i,j,iEq,s1,e1) = 0.5_prec*( boundarySolution(i,j,iEq,s1,e1) +&
                                                           externalState(ii,jj,iEq,bID)  )*& 
                                                         nHat(k,i,j,s1,e1)
               ENDDO
            ELSE
               DO k = 1, 3
                  ! Calculate the Bassi-Rebay flux for the stress tensor.
                  stressFlux(k,i,j,iEq,s1,e1) = 0.5_prec*( boundarySolution(i,j,iEq,s1,e1)/&
                                                         (boundarySolution(i,j,4,s1,e1)+&
                                                          boundarySolution_static(i,j,4,s1,e1)) +&
                                                         externalState(ii,jj,iEq,bID)/&
                                                         (externalState(ii,jj,4,bID)+&
                                                          externalStatic(ii,jj,4,bID))  )*& 
                                                         nHat(k,i,j,s1,e1)
               ENDDO
            ENDIF
            
         ENDDO
         
      ENDIF 


 END SUBROUTINE BoundaryFace_StateFlux_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE CalculateFlux_CUDAKernel( solution, static, Ja, Jac, flux )
  
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)   :: solution(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)   :: static(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)   :: Ja(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:3,1:3,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)   :: Jac(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(out)  :: flux(1:3,0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
     ! Local
    INTEGER            :: i, j, k, row
    INTEGER            :: iEl, iEq
    REAL(prec)         :: F
    REAL(prec)         :: flux_local(1:3)
    
    iEq = blockIDx % x
    iEl = blockIDx % y
    
    i = threadIdx % x - 1
    j = threadIdx % y - 1
    k = threadIdx % z - 1

    flux_local(1) = 0.0_prec
    flux_local(2) = 0.0_prec
    flux_local(3) = 0.0_prec

    ! Here the flux tensor in physical space is calculated and rotated to give the
    ! contravariant flux tensor in the reference computational DOmain.

    !//////////////////////////////// Advection ///////////////////////////////////////!
    DO row = 1, 3

      F = solution(i,j,k,row,iEl)*( solution(i,j,k,iEq,iEl) + static(i,j,k,iEq,iEl) )/&
                                  ( solution(i,j,k,4,iEl) + static(i,j,k,4,iEl) )


      flux_local(1) = flux_local(1) + Ja(i,j,k,row,1,iEl)*F

      flux_local(2) = flux_local(2) + Ja(i,j,k,row,2,iEl)*F

      flux_local(3) = flux_local(3) + Ja(i,j,k,row,3,iEl)*F

    ENDDO

    ! //////////////////// Pressure (Momentum only) /////////////////////////// !

    IF( iEq <= 3 )THEN

      flux(1,i,j,k,iEq,iEl) = flux_local(1) + Ja(i,j,k,iEq,1,iEl)*solution(i,j,k,nEq_dev,iEl)

      flux(2,i,j,k,iEq,iEl) = flux_local(2) + Ja(i,j,k,iEq,2,iEl)*solution(i,j,k,nEq_dev,iEl)

      flux(3,i,j,k,iEq,iEl) = flux_local(3) + Ja(i,j,k,iEq,3,iEl)*solution(i,j,k,nEq_dev,iEl)

    ELSE

      flux(1,i,j,k,iEq,iEl) = flux_local(1)
      flux(2,i,j,k,iEq,iEl) = flux_local(2)
      flux(3,i,j,k,iEq,iEl) = flux_local(3)

    ENDIF
  
  END SUBROUTINE CalculateFlux_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE CalculateSourceTerms_CUDAKernel( solution, solutionGradient, static, staticSource, source, drag )
  
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)    :: solution(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: solutionGradient(1:3,0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: static(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: drag(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: staticSource(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(inout) :: source(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
     ! Local
    INTEGER    :: i, j, k, row, col
    INTEGER    :: iEl, iEq
    REAL(prec) :: F
    
    iEq = blockIDx % x
    iEl = blockIDx % y
    
    i = threadIdx % x - 1
    j = threadIdx % y - 1
    k = threadIdx % z - 1
    
    
    F = sqrt( solution(i,j,k,1,iEl)**2 + &
      solution(i,j,k,2,iEl)**2 + &
      solution(i,j,k,3,iEl)**2 ) /&
      (solution(i,j,k,4,iEl) + static(i,j,k,4,iEl))
    
    IF( iEq == 1 )THEN

      source(i,j,k,1,iEl) = -drag(i,j,k,iEl)*solution(i,j,k,1,iEl)*F-&
                            solution(i,j,k,3,iEl)*fRotY_dev +&
                            solution(i,j,k,2,iEl)*fRotz_dev
    
    ELSEIF( iEq == 2 )THEN

      source(i,j,k,2,iEl) = -drag(i,j,k,iEl)*solution(i,j,k,2,iEl)*F -&
                            solution(i,j,k,1,iEl)*fRotZ_dev +&
                            solution(i,j,k,3,iEl)*fRotX_dev

    ELSEIF( iEq == 3 )THEN ! Add in the buoyancy acceleration

      source(i,j,k,3,iEl) = -drag(i,j,k,iEl)*solution(i,j,k,3,iEl)*F -&
                            solution(i,j,k,2,iEl)*fRotX_dev +&
                            solution(i,j,k,1,iEl)*fRotY_dev -&
                            solution(i,j,k,4,iEl)*g_dev

    ! When the in-situ temperature formulation is used, we must add in the adiabatic heating term
    ! due to fluid convergences. This term is -( P_{total}/C_v )*div( u )
    ELSEIF( iEq == 5 )THEN 

      source(i,j,k,5,iEl) = -( static(i,j,k,nEq_dev,iEl) + solution(i,j,k,nEq_dev,iEl) )*&
                             ( solutionGradient(1,i,j,k,1,iEl) + &
                               solutionGradient(2,i,j,k,2,iEl) + &
                               solutionGradient(3,i,j,k,3,iEl) )/Cv_dev


    ENDIF
  
  END SUBROUTINE CalculateSourceTerms_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE Update_FluidStress_BCs_CUDAKernel( element_to_blockID, boundaryID, elementIDs, elementSides, &
                                                                 externalStress, stressBsols, prescribedStress )
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: element_to_blockID(1:nGlobalElements_dev)
    INTEGER, DEVICE, INTENT(in)     :: boundaryID(1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
    REAL(prec), DEVICE, INTENT(inout) :: externalStress(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)  :: stressBsols(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: prescribedStress(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:nBoundaryFaces_dev)
     ! Local
    INTEGER    :: iEq, bID, i, j, k
    INTEGER    :: iFace, p2
    INTEGER    :: e1, e2, s1, s2, m
    
    iFace = blockIdx % x
    iEq   = blockIDx % y
    i     = threadIdx % x-1
    j     = threadIdx % y-1
    
    bID   = boundaryID(iFace) 
    e1    = elementIDs(1,iFace)
    e2    = elementIDs(2,iFace)
    s1    = elementSides(1,iFace)

    IF( e2 < 0 )THEN
      externalStress(i,j,iEq,bID) = stressBsols(i,j,iEq,s1,e1)
    ENDIF
  
  END SUBROUTINE Update_FluidStress_BCs_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE InternalFace_StressFlux_CUDAKernel( boundarySolution, boundaryStress, &
                                                                    staticBoundarySolution, lengthScale, nHat, boundaryStressFlux, &
                                                                    elementIDs, elementSides, boundaryID, iMap, jMap )
  
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: boundarySolution(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryStress(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: staticBoundarySolution(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:polyDeg_dev,0:polyDeg_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: lengthScale(0:polyDeg_dev,0:polyDeg_dev,1:6,1:nEl_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: boundaryID(1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: iMap(0:polyDeg_dev,0:polyDeg_dev,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: jMap(0:polyDeg_dev,0:polyDeg_dev,1:nFaces_dev)
    REAL(prec), DEVICE, INTENT(inout) :: boundaryStressFlux(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:6,1:nEl_dev)
     ! Local
    INTEGER    :: iEl, iFace
    INTEGER    :: i, j, iEq
    INTEGER    :: ii, jj
    INTEGER    :: e1, s1, e2, s2, bID
    INTEGER    :: m
    REAL(prec) :: norm, rhoOut, rhoIn
    
    iFace = blockIdx % x
    iEq   = blockIdx % y
    j     = threadIdx % y-1
    i     = threadIdx % x-1
    
    e1 = elementIDs(1,iFace)
    s1 = elementSides(1,iFace)
    e2 = elementIDs(2,iFace)
    s2 = ABS(elementSides(2,iFace))
    bID = boundaryID(iFace)
    
    ii = iMap(i,j,iFace)
    jj = jMap(i,j,iFace)
    
    norm = sqrt( nHat(1,i,j,s1,e1)**2 + nHat(2,i,j,s1,e1)**2 + nHat(3,i,j,s1,e1)**2 )
    IF( bID == 0 )THEN

      IF( iEq == 4 )THEN

        boundaryStressFlux(i,j,iEq,s1,e1) = 0.5_prec*( boundaryStress(i,j,iEq,s1,e1)- &
                                                       boundaryStress(ii,jj,iEq,s2,e2) )! +&



      ELSE

        rhoOut = (staticBoundarySolution(ii,jj,4,s2,e2)+boundarySolution(ii,jj,4,s2,e2) )
        rhoIn  = (staticBoundarySolution(i,j,4,s1,e1)+boundarySolution(i,j,4,s1,e1) )

        boundaryStressFlux(i,j,iEq,s1,e1) = 0.5_prec*( rhoIn*boundaryStress(i,j,iEq,s1,e1)-&
                                                       rhoOut*boundaryStress(ii,jj,iEq,s2,e2) )


      ENDIF

      boundaryStressFlux(i,j,iEq,s1,e1) = boundaryStressFlux(i,j,iEq,s1,e1) +&
                                          ( viscosity_dev*boundarySolution(ii,jj,iEq,s2,e2)-&
                                            viscosity_dev*boundarySolution(i,j,iEq,s1,e1) )/&
                                          ( 0.5_prec*(lengthScale(i,j,s1,e1)+lengthScale(ii,jj,s2,e2)) )*norm


      boundaryStressFlux(ii,jj,iEq,s2,e2) = -boundaryStressFlux(i,j,iEq,s1,e1)

    ENDIF

  END SUBROUTINE InternalFace_StressFlux_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE BoundaryFace_StressFlux_CUDAKernel( boundarySolution, externalSolution, boundaryStress, externalStress, &
                                                                    staticBoundarySolution, externalStaticSolution, lengthScale, nHat, boundaryStressFlux, &
                                                                    elementIDs, elementSides, boundaryID, iMap, jMap )
  
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: boundarySolution(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: externalSolution(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryStress(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: externalStress(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)  :: staticBoundarySolution(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: externalStaticSolution(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)  :: lengthScale(0:polyDeg_dev,0:polyDeg_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:polyDeg_dev,0:polyDeg_dev,1:6,1:nEl_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: boundaryID(1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: iMap(0:polyDeg_dev,0:polyDeg_dev,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: jMap(0:polyDeg_dev,0:polyDeg_dev,1:nFaces_dev)
    REAL(prec), DEVICE, INTENT(inout) :: boundaryStressFlux(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:6,1:nEl_dev)
     ! Local
    INTEGER    :: iEl, iFace
    INTEGER    :: i, j, iEq
    INTEGER    :: ii, jj, bID
    INTEGER    :: e1, s1, e2, s2
    INTEGER    :: m
    REAL(prec) :: norm, rhoOut, rhoIn
    
    iFace = blockIdx % x
    iEq   = blockIdx % y
    j     = threadIdx % y-1
    i     = threadIdx % x-1
    
    e1 = elementIDs(1,iFace)
    s1 = elementSides(1,iFace)
    e2 = elementIDs(2,iFace)
    s2 = ABS(elementSides(2,iFace))
    bID  = ABS(boundaryID(iFace))
    
    ii = iMap(i,j,iFace)
    jj = jMap(i,j,iFace)
    
    norm = sqrt( nHat(1,i,j,s1,e1)**2 + nHat(2,i,j,s1,e1)**2 + nHat(3,i,j,s1,e1)**2 )
    IF( e2 < 0 )THEN

      IF( iEq == 4 )THEN

        boundaryStressFlux(i,j,iEq,s1,e1) = 0.5_prec*( boundaryStress(i,j,iEq,s1,e1) - &
                                                       externalStress(ii,jj,iEq,bID) )



      ELSE

        rhoOut = (externalStaticSolution(ii,jj,4,bID)+externalSolution(ii,jj,4,bID) )
        rhoIn  = (staticBoundarySolution(i,j,4,s1,e1)+boundarySolution(i,j,4,s1,e1) )

        boundaryStressFlux(i,j,iEq,s1,e1) = 0.5_prec*( rhoIn*boundaryStress(i,j,iEq,s1,e1)-&
                                                       rhoOut*externalStress(ii,jj,iEq,bID) )


      ENDIF

      boundaryStressFlux(i,j,iEq,s1,e1) = boundaryStressFlux(i,j,iEq,s1,e1) +&
                                          ( viscosity_dev*externalSolution(ii,jj,iEq,bID)-&
                                            viscosity_dev*boundarySolution(i,j,iEq,s1,e1) )/(lengthScale(i,j,s1,e1))*norm



    ENDIF
  
  END SUBROUTINE BoundaryFace_StressFlux_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE EquationOfState_CUDAKernel( solution, static )
    ! This routine calculates the anomalous pressure referenced to the static state.
    ! The pressure is calculated using the ideal gas law.
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(inout) :: solution(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: static(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
     ! Local
    INTEGER :: i, j, k, iEl
    
    iEl = blockIdx % x
    i   = threadIdx % x - 1
    j   = threadIdx % y - 1
    k   = threadIdx % z - 1
    
     ! Pressure = rho*R*T
     solution(i,j,k,nEq_dev,iEl) = solution(i,j,k,5,iEl)*R_dev

  END SUBROUTINE EquationOfState_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE State_Mapped_DG_Divergence_3D_CUDAKernel( f, fnAtBoundaries, divF, boundaryMatrix, dgDerivativeMatrixTranspose, quadratureWeights, Jac )
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: f(1:3,0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: fnAtBoundaries(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: Jac(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryMatrix(0:polydeg_dev,0:1)
    REAL(prec), DEVICE, INTENT(in)  :: dgDerivativeMatrixTranspose(0:polydeg_dev,0:polydeg_dev)
    REAL(prec), DEVICE, INTENT(in)  :: quadratureWeights(0:polydeg_dev)
    REAL(prec), DEVICE, INTENT(out) :: divF(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
    ! Local
    INTEGER            :: i, j, k, iVar, iEl, ii
    REAL(prec)         :: df
    REAL(prec), SHARED :: fLocal(1:3,0:7,0:7,0:7)
    
    
      iVar = blockIDx % x
      iEl  = blockIDx % y
      
      i = threadIdx % x - 1
      j = threadIdx % y - 1
      k = threadIdx % z - 1
      
      fLocal(1,i,j,k) = f(1,i,j,k,iVar,iEl)
      fLocal(2,i,j,k) = f(2,i,j,k,iVar,iEl)
      fLocal(3,i,j,k) = f(3,i,j,k,iVar,iEl)
      
      CALL syncthreads( )
      
      
        df = 0.0_prec
        DO ii = 0, polydeg_dev
          df = df + dgDerivativeMatrixTranspose(ii,i)*fLocal(1,ii,j,k) + &
                    dgDerivativeMatrixTranspose(ii,j)*fLocal(2,i,ii,k) + &
                    dgDerivativeMatrixTranspose(ii,k)*fLocal(3,i,j,ii)
        ENDDO
       
        divF(i,j,k,iVar,iEl) = ( df+ ( fnAtBoundaries(i,k,iVar,1,iEl)*boundaryMatrix(j,0) + &
                                        fnAtBoundaries(i,k,iVar,3,iEl)*boundaryMatrix(j,1) )/&
                                      quadratureWeights(j) + &
                                      ( fnAtBoundaries(j,k,iVar,4,iEl)*boundaryMatrix(i,0) + &
                                        fnAtBoundaries(j,k,iVar,2,iEl)*boundaryMatrix(i,1) )/&
                                      quadratureWeights(i) + &
                                      ( fnAtBoundaries(i,j,iVar,5,iEl)*boundaryMatrix(k,0) + &
                                        fnAtBoundaries(i,j,iVar,6,iEl)*boundaryMatrix(k,1) )/&
                                      quadratureWeights(k) )/Jac(i,j,k,iEl)
                                      

  END SUBROUTINE State_Mapped_DG_Divergence_3D_CUDAKernel

  ATTRIBUTES(Global) SUBROUTINE Stress_Mapped_DG_Divergence_3D_CUDAKernel( f, fnAtBoundaries, divF, boundaryMatrix, dgDerivativeMatrixTranspose, quadratureWeights, Jac )
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: f(1:3,0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nStress_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: fnAtBoundaries(0:polydeg_dev,0:polydeg_dev,1:nStress_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: Jac(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryMatrix(0:polydeg_dev,0:1)
    REAL(prec), DEVICE, INTENT(in)  :: dgDerivativeMatrixTranspose(0:polydeg_dev,0:polydeg_dev)
    REAL(prec), DEVICE, INTENT(in)  :: quadratureWeights(0:polydeg_dev)
    REAL(prec), DEVICE, INTENT(out) :: divF(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nStress_dev,1:nEl_dev)
    ! Local
    INTEGER            :: i, j, k, iVar, iEl, ii
    REAL(prec)         :: df
    REAL(prec), SHARED :: fLocal(1:3,0:7,0:7,0:7)
    
    
      iVar = blockIDx % x
      iEl  = blockIDx % y
      
      i = threadIdx % x - 1
      j = threadIdx % y - 1
      k = threadIdx % z - 1
    
      fLocal(1,i,j,k) = f(1,i,j,k,iVar,iEl)
      fLocal(2,i,j,k) = f(2,i,j,k,iVar,iEl)
      fLocal(3,i,j,k) = f(3,i,j,k,iVar,iEl)
      
      CALL syncthreads( )
      
      df = 0.0_prec
      DO ii = 0, polydeg_dev
        df = df + dgDerivativeMatrixTranspose(ii,i)*fLocal(1,ii,j,k) + &
                  dgDerivativeMatrixTranspose(ii,j)*fLocal(2,i,ii,k) + &
                  dgDerivativeMatrixTranspose(ii,k)*fLocal(3,i,j,ii)
      ENDDO
      
      divF(i,j,k,iVar,iEl) = ( df+ ( fnAtBoundaries(i,k,iVar,1,iEl)*boundaryMatrix(j,0) + &
                                      fnAtBoundaries(i,k,iVar,3,iEl)*boundaryMatrix(j,1) )/&
                                    quadratureWeights(j) + &
                                    ( fnAtBoundaries(j,k,iVar,4,iEl)*boundaryMatrix(i,0) + &
                                      fnAtBoundaries(j,k,iVar,2,iEl)*boundaryMatrix(i,1) )/&
                                    quadratureWeights(i) + &
                                    ( fnAtBoundaries(i,j,iVar,5,iEl)*boundaryMatrix(k,0) + &
                                      fnAtBoundaries(i,j,iVar,6,iEl)*boundaryMatrix(k,1) )/&
                                    quadratureWeights(k) )/Jac(i,j,k,iEl)
                                      

  END SUBROUTINE Stress_Mapped_DG_Divergence_3D_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE CalculateStateAtBoundaries_CUDAKernel( f, fAtBoundaries, boundaryMatrix ) 
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: f(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryMatrix(0:polydeg_dev,0:1)
    REAL(prec), DEVICE, INTENT(out) :: fAtBoundaries(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
    ! Local
    INTEGER    :: iEq, iEl, i, j, k
    REAL(prec) :: bSol(1:6)

      iEq = blockIdx % x
      iEl = blockIdx % y

      k   = threadIdx % y-1
      j   = threadIdx % x-1
      
      bSol(1:6) = 0.0_prec

      DO i = 0, polydeg_dev

        bSol(1) = bSol(1) + boundaryMatrix(i,0)*f(j,i,k,iEq,iEl) ! south
        bSol(2) = bSol(2) + boundaryMatrix(i,1)*f(i,j,k,iEq,iEl) ! east
        bSol(3) = bSol(3) + boundaryMatrix(i,1)*f(j,i,k,iEq,iEl) ! north
        bSol(4) = bSol(4) + boundaryMatrix(i,0)*f(i,j,k,iEq,iEl) ! west
        bSol(5) = bSol(5) + boundaryMatrix(i,0)*f(j,k,i,iEq,iEl) ! botom
        bSol(6) = bSol(6) + boundaryMatrix(i,1)*f(j,k,i,iEq,iEl) ! top

      ENDDO
               
      DO i = 1, 6
        fAtBoundaries(j,k,iEq,i,iEl) = bSol(i)
      ENDDO
      
      
  END SUBROUTINE CalculateStateAtBoundaries_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE CalculateNormalStressAtBoundaries_CUDAKernel( solutionGradient, nHat, boundaryStress, boundaryMatrix ) 
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: solutionGradient(1:3,0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:polyDeg_dev,0:polyDeg_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryMatrix(0:polydeg_dev,0:1)
    REAL(prec), DEVICE, INTENT(out) :: boundaryStress(0:polydeg_dev,0:polydeg_dev,1:nStress_dev,1:6,1:nEl_dev)
    ! Local
    INTEGER    :: iEq, iEl, i, j, k
    REAL(prec) :: fAtBoundaries(1:3,1:6)

      iEq = blockIdx % x
      iEl = blockIdx % y

      i   = threadIdx % x-1
      j   = threadIdx % y-1
            
      fAtBoundaries(1:3,1:6) = 0.0_prec
      
      DO k = 0, polyDeg_dev
       
        fAtBoundaries(1:3,1) = fAtBoundaries(1:3,1) + boundaryMatrix(k,0)*solutionGradient(1:3,i,k,j,iEq,iEl) ! South
        fAtBoundaries(1:3,2) = fAtBoundaries(1:3,2) + boundaryMatrix(k,1)*solutionGradient(1:3,k,i,j,iEq,iEl) ! East
        fAtBoundaries(1:3,3) = fAtBoundaries(1:3,3) + boundaryMatrix(k,1)*solutionGradient(1:3,i,k,j,iEq,iEl) ! North
        fAtBoundaries(1:3,4) = fAtBoundaries(1:3,4) + boundaryMatrix(k,0)*solutionGradient(1:3,k,i,j,iEq,iEl) ! West
        fAtBoundaries(1:3,5) = fAtBoundaries(1:3,5) + boundaryMatrix(k,0)*solutionGradient(1:3,i,j,k,iEq,iEl) ! Bottom
        fAtBoundaries(1:3,6) = fAtBoundaries(1:3,6) + boundaryMatrix(k,1)*solutionGradient(1:3,i,j,k,iEq,iEl) ! Top
       
      ENDDO

      DO k = 1, 6

        boundaryStress(i,j,iEq,k,iEl) = viscosity_dev*( fAtBoundaries(1,k)*nHat(1,i,j,k,iEl) + &
                                                        fAtBoundaries(2,k)*nHat(2,i,j,k,iEl) + &
                                                        fAtBoundaries(3,k)*nHat(3,i,j,k,iEl) )

      ENDDO

      
  END SUBROUTINE CalculateNormalStressAtBoundaries_CUDAKernel

  ATTRIBUTES(Global) SUBROUTINE CalculateSolutionGradient_CUDAKernel( solutionGradient, solution, static, boundaryGradientFlux, Ja, Jac, dgDerivativeMatrixTranspose, boundaryInterpolationMatrix, quadratureWeights )

    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(out) :: solutionGradient(1:3,0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev) 
    REAL(prec), DEVICE, INTENT(in)  :: solution(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev) 
    REAL(prec), DEVICE, INTENT(in)  :: static(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev) 
    REAL(prec), DEVICE, INTENT(in)  :: boundaryGradientFlux(1:3,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:6,1:nEl_dev) 
    REAL(prec), DEVICE, INTENT(in)  :: Ja(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:3,1:3,1:nEl_dev) 
    REAL(prec), DEVICE, INTENT(in)  :: Jac(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEl_dev) 
    REAL(prec), DEVICE, INTENT(in)  :: dgDerivativeMatrixTranspose(0:polyDeg_dev,0:polyDeg_dev) 
    REAL(prec), DEVICE, INTENT(in)  :: boundaryInterpolationMatrix(0:polyDeg_dev,0:1) 
    REAL(prec), DEVICE, INTENT(in)  :: quadratureWeights(0:polyDeg_dev) 
    ! Local
    INTEGER            :: iEq, iEl, idir, i, j, k, ii
    REAL(prec), SHARED :: f(1:3,0:7,0:7,0:7)
    REAL(prec)         :: df!, bgf(1:6)
    
 
    iEq  = blockIdx % x
    iEl  = blockIdx % y
    idir = blockIdx % z

    i   = threadIdx % x-1
    j   = threadIdx % y-1
    k   = threadIdx % z-1

      IF( iEq == 4 )THEN

        f(1,i,j,k) = solution(i,j,k,iEq,iEl)*Ja(i,j,k,idir,1,iEl)
        f(2,i,j,k) = solution(i,j,k,iEq,iEl)*Ja(i,j,k,idir,2,iEl)
        f(3,i,j,k) = solution(i,j,k,iEq,iEl)*Ja(i,j,k,idir,3,iEl)

      ELSE

        f(1,i,j,k) = Ja(i,j,k,idir,1,iEl)*&
                     solution(i,j,k,iEq,iEl)/&
                     (solution(i,j,k,4,iEl)+&
                     static(i,j,k,4,iEl) )


        f(2,i,j,k) = Ja(i,j,k,idir,2,iEl)*&
                     solution(i,j,k,iEq,iEl)/&
                     (solution(i,j,k,4,iEl)+&
                     static(i,j,k,4,iEl) )

        f(3,i,j,k) = Ja(i,j,k,idir,3,iEl)*&
                     solution(i,j,k,iEq,iEl)/&
                     (solution(i,j,k,4,iEl)+&
                     static(i,j,k,4,iEl) )

      ENDIF

      CALL syncthreads( )
 
      df = 0.0_prec
      DO ii = 0, polyDeg_dev
        df = df + dgDerivativeMatrixTranspose(ii,i)*f(1,ii,j,k) + &
                  dgDerivativeMatrixTranspose(ii,j)*f(2,i,ii,k) + &
                  dgDerivativeMatrixTranspose(ii,k)*f(3,i,j,ii)
      ENDDO
 
      solutionGradient(idir,i,j,k,iEq,iEl) =  ( df+ ( boundaryGradientFlux(idir,i,k,ieq,1,iel)*boundaryInterpolationMatrix(j,0) + &
                                                      boundaryGradientFlux(idir,i,k,ieq,3,iel)*boundaryInterpolationMatrix(j,1) )/&
                                                    quadratureWeights(j) + &
                                                    ( boundaryGradientFlux(idir,j,k,ieq,4,iel)*boundaryInterpolationMatrix(i,0) + &
                                                      boundaryGradientFlux(idir,j,k,ieq,2,iel)*boundaryInterpolationMatrix(i,1) )/&
                                                    quadratureWeights(i) + &
                                                    ( boundaryGradientFlux(idir,i,j,ieq,5,iel)*boundaryInterpolationMatrix(k,0) + &
                                                      boundaryGradientFlux(idir,i,j,ieq,6,iel)*boundaryInterpolationMatrix(k,1) )/&
                                                    quadratureWeights(k) )/Jac(i,j,k,iEl)


  END SUBROUTINE CalculateSolutionGradient_CUDAKernel

  ATTRIBUTES(Global) SUBROUTINE CalculateStressFlux_CUDAKernel( solutionGradient, state, static, Ja, stressFlux )

    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: solutionGradient(1:3,0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: state(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: static(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: Ja(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:3,1:3,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(out) :: stressFlux(1:3,0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:nEl_dev)
    ! Local
    INTEGER    :: iEl, iEq, i, j, k, idir
    REAL(prec) :: sf(1:3), rho

    iEq  = blockIdx % x
    iEl  = blockIdx % y

    i   = threadIdx % x-1
    j   = threadIdx % y-1
    k   = threadIdx % z-1

    sf(1) = 0.0_prec
    sf(2) = 0.0_prec
    sf(3) = 0.0_prec


    IF( iEq == 4 )THEN

      DO idir = 1, 3
  
        sf(1) = sf(1) + Ja(i,j,k,idir,1,iEl)*solutionGradient(idir,i,j,k,iEq,iEl)*viscosity_dev
  
        sf(2) = sf(2) + Ja(i,j,k,idir,2,iEl)*solutionGradient(idir,i,j,k,iEq,iEl)*viscosity_dev
  
        sf(3) = sf(3) + Ja(i,j,k,idir,3,iEl)*solutionGradient(idir,i,j,k,iEq,iEl)*viscosity_dev
  
      ENDDO

    ELSE

      rho = state(i,j,k,4,iEl) + static(i,j,k,4,iEl)
      DO idir = 1, 3
  
        sf(1) = sf(1) + Ja(i,j,k,idir,1,iEl)*rho*solutionGradient(idir,i,j,k,iEq,iEl)*viscosity_dev
  
        sf(2) = sf(2) + Ja(i,j,k,idir,2,iEl)*rho*solutionGradient(idir,i,j,k,iEq,iEl)*viscosity_dev
  
        sf(3) = sf(3) + Ja(i,j,k,idir,3,iEl)*rho*solutionGradient(idir,i,j,k,iEq,iEl)*viscosity_dev
  
      ENDDO

    ENDIF

    stressflux(1,i,j,k,iEq,iEl) = sf(1)
    stressflux(2,i,j,k,iEq,iEl) = sf(2)
    stressflux(3,i,j,k,iEq,iEl) = sf(3)

  END SUBROUTINE CalculateStressFlux_CUDAKernel

#endif


END MODULE Fluid_Class


