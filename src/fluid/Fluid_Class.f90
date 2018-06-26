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
  USE SpectralFilter_Class
  USE NodalDGSolution_3D_Class
  USE HexMesh_Class
  USE BoundaryCommunicator_Class
  USE BodyForces_Class
#ifdef HAVE_MPI
  USE MPILayer_Class
#endif

#ifdef HAVE_CUDA
  USE cudafor
#endif


  IMPLICIT NONE


  TYPE Fluid

    REAL(prec) :: simulationTime
    REAL(prec) :: volume, mass, KE, PE, heat
  
    TYPE( MultiTimers )          :: timers
    TYPE( ModelParameters )      :: params
    TYPE( HexMesh )              :: mesh
    TYPE( BoundaryCommunicator ) :: extComm
    TYPE( NodalDG )              :: dGStorage
    TYPE( SpectralFilter )       :: filter
    TYPE( BodyForces )           :: sourceTerms
    TYPE( NodalDGSolution_3D )   :: static
    TYPE( NodalDGSolution_3D )   :: state
    TYPE( NodalDGSolution_3D )   :: smoothState
    TYPE( NodalDGSolution_3D )   :: stressTensor
    TYPE( NodalDGSolution_3D )   :: sgsCoeffs

#ifdef HAVE_MPI
    TYPE( MPILayer )             :: mpiStateHandler
    TYPE( MPILayer )             :: mpiStressHandler
    TYPE( MPILayer )             :: mpiSGSHandler
#endif

  CONTAINS

    PROCEDURE :: Build => Build_Fluid
    PROCEDURE :: Trash => Trash_Fluid

    PROCEDURE :: BuildHexMesh => BuildHexMesh_Fluid
    PROCEDURE :: CalculateStaticState => CalculateStaticState_Fluid

    ! Time integrators
    PROCEDURE :: ForwardStepRK3        => ForwardStepRK3_Fluid
!    PROCEDURE :: CrankNicholsonRHS     => CrankNicholsonRHS_Fluid

    PROCEDURE :: GlobalTimeDerivative   => GlobalTimeDerivative_Fluid
    PROCEDURE :: EquationOfState        => EquationOfState_Fluid
    PROCEDURE :: UpdateExternalState    => UpdateExternalState_Fluid
    PROCEDURE :: InternalFace_StateFlux => InternalFace_StateFlux_Fluid
    PROCEDURE :: BoundaryFace_StateFlux => BoundaryFace_StateFlux_Fluid
    PROCEDURE :: MappedTimeDerivative   => MappedTimeDerivative_Fluid

    PROCEDURE :: CalculateSGSCoefficients => CalculateSGSCoefficients_Fluid
    PROCEDURE :: UpdateExternalSGS        => UpdateExternalSGS_Fluid

    PROCEDURE :: CalculateSolutionGradient         => CalculateSolutionGradient_Fluid
    PROCEDURE :: CalculateNormalStressAtBoundaries => CalculateNormalStressAtBoundaries_Fluid
    PROCEDURE :: CalculateStressFlux               => CalculateStressFlux_Fluid
    PROCEDURE :: UpdateExternalStress              => UpdateExternalStress_Fluid
    PROCEDURE :: InternalFace_StressFlux           => InternalFace_StressFlux_Fluid
    PROCEDURE :: BoundaryFace_StressFlux           => BoundaryFace_StressFlux_Fluid

    PROCEDURE :: Diagnostics           => Diagnostics_Fluid
    PROCEDURE :: OpenDiagnosticsFiles  => OpenDiagnosticsFiles_Fluid
    PROCEDURE :: WriteDiagnostics      => WriteDiagnostics_Fluid
    PROCEDURE :: CloseDiagnosticsFiles => CloseDiagnosticsFiles_Fluid

    PROCEDURE :: WriteTecplot => WriteTecplot_Fluid
    PROCEDURE :: WritePickup  => WritePickup_Fluid
    PROCEDURE :: ReadPickup   => ReadPickup_Fluid

  END TYPE Fluid


  INTEGER, PARAMETER, PRIVATE :: nDiagnostics = 5
  INTEGER, PRIVATE            :: diagUnits(1:5)

#ifdef PASSIVE_TRACERS
  INTEGER, PARAMETER :: nEquations   = 7
#else
  INTEGER, PARAMETER :: nEquations   = 6
#endif

#ifdef HAVE_CUDA
 INTEGER, CONSTANT    :: nEq_dev
 INTEGER, CONSTANT    :: nStress_dev
 INTEGER, CONSTANT    :: nSGS_dev
 INTEGER, CONSTANT    :: polyDeg_dev
 INTEGER, CONSTANT    :: nEl_dev
 INTEGER, CONSTANT    :: myRank_dev
 INTEGER, CONSTANT    :: nProc_dev
 INTEGER, CONSTANT    :: nNeighbors_dev
 INTEGER, CONSTANT    :: nFaces_dev
 INTEGER, CONSTANT    :: nBoundaryFaces_dev
 INTEGER, CONSTANT    :: bufferSize_dev
 REAL(prec), CONSTANT :: R_dev, Cv_dev, P0_dev, hCapRatio_dev, rC_dev, g_dev, dt_dev
 REAL(prec), CONSTANT :: viscLengthScale_dev, dScale_dev, Cd_dev
 REAL(prec), CONSTANT :: fRotX_dev, fRotY_dev, fRotZ_dev
#endif

CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
  SUBROUTINE Build_Fluid( myDGSEM, setupSuccess )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
    LOGICAL, INTENT(inout)      :: setupSuccess
    ! Local
#ifdef HAVE_CUDA
    INTEGER(KIND=cuda_count_KIND) :: freebytes, totalbytes
    INTEGER                       :: iStat, cudaDeviceNumber, nDevices
#endif


    CALL myDGSEM % params % Build( setupSuccess )
    myDGSEM % simulationTime = myDGSEM % params % startTime

    IF( .NOT. SetupSuccess ) THEN
      PRINT(MsgFMT), 'S/R Build_Fluid : Halting before building,'
      RETURN
    ENDIF

    ! This call to the extComm % ReadPickup reads in the external communicator
    ! data. If MPI is enabled, MPI is initialized. If CUDA and MPI are enabled
    ! the device for each rank is also set here.
    CALL myDGSEM % extComm % ReadPickup(  )

#ifdef HAVE_CUDA
    CALL UpdateDeviceDictionary( )
#endif

    ! Construct the DATA structure that holds the derivative and interpolation matrices
    ! and the quadrature weights. This call will also perform the device copies.
    CALL myDGSEM % dGStorage % Build( UniformPoints(-1.0_prec, 1.0_prec,myDGSEM % params % nPlot), &
      myDGSEM % params % polyDeg, myDGSEM % params % nPlot, GAUSS )

    CALL myDGSEM % filter % Build( myDGSEM % dgStorage % interp % interpolationPoints,&
      myDGSEM % dgStorage % quadratureWeights, &
      myDGSEM % params % polyDeg, myDGSEM % params % nCutoff, &
      myDGSEM % params % filter_a, &
      myDGSEM % params % filter_b, &
      myDGSEM % params % filterType )

    CALL myDGSEM % BuildHexMesh(  )


    CALL myDGSEM % sourceTerms % Build( myDGSEM % params % polyDeg, nEquations, &
      myDGSEM % mesh % elements % nElements )

    CALL myDGSEM % state % Build( myDGSEM % params % polyDeg, nEquations, &
      myDGSEM % mesh % elements % nElements, myDGSEM % extComm % nBoundaries )

    CALL myDGSEM % static % Build( myDGSEM % params % polyDeg, nEquations, &
      myDGSEM % mesh % elements % nElements, myDGSEM % extComm % nBoundaries )

    CALL myDGSEM % smoothState % Build( myDGSEM % params % polyDeg, nEquations, &
      myDGSEM % mesh % elements % nElements, myDGSEM % extComm % nBoundaries )


    ! The "sgsCoeffs" attribute CONTAINS coefficients for the
    ! subgrid scale parameterization. Currently, these coefficients
    ! are the eddy viscosities for the momentum equations (assuming
    ! isotropic turbulence), and the eddy dIFfusivities for the
    ! potential temperature and density equations.

    CALL myDGSEM % sgsCoeffs % Build( myDGSEM % params % polyDeg, &
                                      myDGSEM % state % nEquations-1, &
                                      myDGSEM % mesh % elements % nElements, &
                                      myDGSEM % extComm % nBoundaries )

    ! Initially set all of the SGS coefficients to the "viscosity". In the event
    ! the Laplacian model is USEd, this will be the laplacian coefficient that is
    ! USEd for the momentum, potential temperature, and density equations.
    ! IF another SGS model is USEd (e.g. SpectralEKE ), THEN these values will be
    ! overwritten in the "CalculateSGSCoefficients" routine.

    CALL myDGSEM % stressTensor % Build( myDGSEM % params % polyDeg, &
                                         myDGSEM % state % nEquations-1, &
                                         myDGSEM % mesh % elements % nElements, &
                                         myDGSEM % extComm % nBoundaries )

    myDGSEM % sgsCoeffs % solution         = myDGSEM % params % viscosity
    myDGSEM % sgsCoeffs % boundarySolution = myDGSEM % params % viscosity
    myDGSEM % sgsCoeffs % externalState    = myDGSEM % params % viscosity

#ifdef HAVE_MPI
    CALL myDGSEM % mpiStateHandler % Build( myDGSEM % extComm, myDGSEM % params % polyDeg, myDGSEM % state % nEquations )
    CALL myDGSEM % mpiStressHandler % Build( myDGSEM % extComm, myDGSEM % params % polyDeg, myDGSEM % stressTensor % nEquations )
    CALL myDGSEM % mpiSGSHandler % Build( myDGSEM % extComm, myDGSEM % params % polyDeg, myDGSEM % sgsCoeffs % nEquations )
#endif

#ifdef HAVE_CUDA
    CALL myDGSEM % sgsCoeffs % UpdateDevice( )

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

      nEq_dev     = myDGSEM % state % nEquations
      nStress_dev = myDGSEM % stressTensor % nEquations
      nSGS_dev = myDGSEM % sgsCoeffs % nEquations
      polydeg_dev = myDGSEM % params % polyDeg
      nEl_dev     = myDGSEM % mesh % elements % nElements
      nFaces_dev     = myDGSEM % mesh % faces % nFaces
      nBoundaryFaces_dev     = myDGSEM % extComm % nBoundaries
      myRank_dev  = myDGSEM % extComm % myRank
      nProc_dev   = myDGSEM % extComm % nProc
#ifdef HAVE_MPI
      nNeighbors_dev = myDGSEM % extComm % nNeighbors
      bufferSize_dev = myDGSEM % extComm % maxbufferSize
#endif      

#endif

    ! Read the initial conditions, static state, and the boundary communicator
    CALL myDGSEM % ReadPickup( )

#ifdef TIMING
    ! Setup timers
    CALL myDGSEM % timers % Build( )
    CALL myDGSEM % timers % AddTimer( 'Forward_Step_RK3', 1 )
    CALL myDGSEM % timers % AddTimer( 'Calculate_State_At_Boundaries', 2 )
    CALL myDGSEM % timers % AddTimer( 'MPI_State_Exchange', 3 )
    CALL myDGSEM % timers % AddTimer( 'Update_External_State', 4 )
    CALL myDGSEM % timers % AddTimer( 'Internal_Face_State_Flux', 5 )
    CALL myDGSEM % timers % AddTimer( 'Finalize_MPI_State_Exchange', 6 )
    CALL myDGSEM % timers % AddTimer( 'Boundary_Face_State_Flux', 7 )
    CALL myDGSEM % timers % AddTimer( 'Filter3D_for_SmoothedState', 8 )
    CALL myDGSEM % timers % AddTimer( 'Calculate_SGS_Coefficients', 9 )
    CALL myDGSEM % timers % AddTimer( 'Calculate_SGS_At_Boundaries', 10 )
    CALL myDGSEM % timers % AddTimer( 'MPI_SGS_Exchange', 11 )
    CALL myDGSEM % timers % AddTimer( 'Calculate_Solution_Gradient', 12 )
    CALL myDGSEM % timers % AddTimer( 'Finalize_SGS_Exchange', 13 )
    CALL myDGSEM % timers % AddTimer( 'Calculate_Normal_Stress_At_Boundaries', 14 )
    CALL myDGSEM % timers % AddTimer( 'MPI_Stress_Exchange', 15 )
    CALL myDGSEM % timers % AddTimer( 'Calculate_Stress_Flux', 16 )
    CALL myDGSEM % timers % AddTimer( 'Update_External_Stress', 17 )
    CALL myDGSEM % timers % AddTimer( 'Finalize_Stress_Exchange', 18 )
    CALL myDGSEM % timers % AddTimer( 'Boundary_Face_Stress_Flux', 19 )
    CALL myDGSEM % timers % AddTimer( 'Calculate_Stress_Flux_Divergence', 20 )
    CALL myDGSEM % timers % AddTimer( 'Mapped_Time_Derivative', 21 )
#endif

  END SUBROUTINE Build_Fluid
!
  SUBROUTINE Trash_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM

#ifdef HAVE_MPI
    CALL MPI_BARRIER( myDGSEM % extComm % mpi_comm )
#endif

    PRINT*, '    S/R Trash_Fluid : Clearing memory.'

    CALL myDGSEM % params % Trash( )
    CALL myDGSEM % mesh % Trash(  )
    CALL myDGSEM % extComm % Trash( )
    CALL myDGSEM % dGStorage % Trash( )
    CALL myDGSEM % filter % Trash( )
    CALL myDGSEM % sourceTerms % Trash( )
    CALL myDGSEM % static % Trash( )
    CALL myDGSEM % state % Trash( )
    CALL myDGSEM % smoothState % Trash( )
    CALL myDGSEM % stressTensor % Trash( )
    CALL myDGSEM % sgsCoeffs % Trash( )

#ifdef HAVE_MPI
    CALL myDGSEM % mpiStateHandler % Trash( )
    CALL myDGSEM % mpiStressHandler % Trash( )
    CALL myDGSEM % mpiSGSHandler % Trash( )
#endif

#ifdef TIMING
    CALL myDGSEM % timers % Write_MultiTimers( )
    CALL myDGSEM % timers % Trash( )
#endif

  END SUBROUTINE Trash_Fluid
!
  SUBROUTINE BuildHexMesh_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM
    ! Local
    CHARACTER(4) :: rankChar

    WRITE( rankChar, '(I4.4)' )myDGSEM % extComm % myRank
    PRINT*,'    Module Fluid_Class.f90 : S/R BuildHexMesh :'

    PRINT*,'      Reading mesh from '//TRIM(myDGSEM % params % SELFMeshFile)//'.'//rankChar//'.smesh '

    ! This loads in the mesh from the "pc-mesh file" and sets up the device arrays for the mesh
    CALL myDGSEM % mesh % ReadSELFMeshFile( TRIM(myDGSEM % params % SELFMeshFile)//'.'//rankChar )

    ! Multiply the mesh positions to scale the size of the mesh
    CALL myDGSEM % mesh % ScaleTheMesh( myDGSEM % dgStorage % interp, &
                                        myDGSEM % params % xScale, &
                                        myDGSEM % params % yScale, &
                                        myDGSEM % params % zScale )


  END SUBROUTINE BuildHexMesh_Fluid
!
  SUBROUTINE ForwardStepRK3_Fluid( myDGSEM, nT )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
    INTEGER, INTENT(in)         :: nT
    ! LOCAL
    REAL(prec)                  :: t0
#ifdef HAVE_CUDA
    REAL(prec)         :: t, dt
    REAL(prec), DEVICE :: t_dev
    REAL(prec), DEVICE :: G3D(0:myDGSEM % params % polyDeg,&
                             0:myDGSEM % params % polyDeg,&
                             0:myDGSEM % params % polyDeg,&
                             1:myDGSEM % state % nEquations,&
                             1:myDGSEM % mesh % elements % nElements)
    INTEGER            :: iT, m, iStat
    TYPE(dim3)         :: fgrid, grid, tBlock
#else
    REAL(prec) :: t, dt, rk3_a_local, rk3_g_local
    REAL(prec) :: G3D(0:myDGSEM % params % polyDeg,&
                      0:myDGSEM % params % polyDeg,&
                      0:myDGSEM % params % polyDeg,&
                      1:myDGSEM % state % nEquations,&
                      1:myDGSEM % mesh % elements % nElements)
    INTEGER    :: m, iEl, iT, i, j, k, iEq

#endif


#ifdef TIMING
    CALL myDGSEM % timers % StartTimer( 1 )
#endif


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
      dt_dev = dt

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


      myDGSEM % simulationTime = myDGSEM % simulationTime + dt

    ENDIF

#else

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

        DO iEl = 1, myDGSEM % mesh % elements % nElements
          DO iEq = 1, myDGSEM % state % nEquations-1
            DO k = 0, myDGSEM % params % polyDeg
              DO j = 0, myDGSEM % params % polyDeg
                DO i = 0, myDGSEM % params % polyDeg

                  G3D(i,j,k,iEq,iEl) = rk3_a(m)*G3D(i,j,k,iEq,iEl) - ( myDGSEM % state % fluxDivergence(i,j,k,iEq,iEl) -&
                    myDGSEM % stressTensor % fluxDivergence(i,j,k,iEq,iEl) )/myDGSEM % mesh % elements % J(i,j,k,iEl) + &
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
                    myDGSEM % stressTensor % fluxDivergence(i,j,k,iEq,iEl) )/myDGSEM % mesh % elements % J(i,j,k,iEl) + &
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
#endif

#ifdef TIMING
    CALL myDGSEM % timers % StopTimer( 1 )
#endif

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

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
    REAL(prec), INTENT(in)      :: tn
#ifdef HAVE_CUDA
    TYPE(dim3) :: tBlock, grid
    INTEGER    :: istat
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
!  Here, the solution within each element is interpolated to the faces of each
!  element in order to prepare for computing the external state for enforcing
!  boundary conditions, Riemann Fluxes, and MPI DATA exchanges that need to
!  occur.

#ifdef TIMING
    CALL myDGSEM % timers % StartTimer( 2 )
#endif

#ifdef HAVE_CUDA

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1 , &
                  1 )
    grid = dim3(myDGSEM % state % nEquations, myDGSEM % state % nElements, 1)  

    CALL CalculateStateAtBoundaries_CUDAKernel<<<grid, tBlock>>>( myDGSEM % state % solution_dev, &
                                                                  myDGSEM % state % boundarySolution_dev,  &
                                                                  myDGSEM % dgStorage % boundaryInterpolationMatrix_dev )
#else

    CALL myDGSEM % state % Calculate_Solution_At_Boundaries( myDGSEM % dgStorage )

#endif

#ifdef TIMING
    CALL myDGSEM % timers % StopTimer( 2 )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! When MPI is USEd, the boundary solutions that are stored on faces shared with
! a neighboring rank are passed to that neighboring rank. Additionally, the
! perceived "external state" is received from the neighboring rank. Calling this
! routine is dependent on the result of CalculateBoundarySolutio, but can be
! DOne at the same time as UpdateExternalState; DOing so should hide some
! communication costs.

#ifdef TIMING
    CALL myDGSEM % timers % StartTimer( 3 )
#endif

#ifdef HAVE_MPI
    CALL myDGSEM % mpiStateHandler % MPI_Exchange( myDGSEM % state, &
                                                   myDGSEM % mesh % faces, &
                                                   myDGSEM % extComm )
#endif

#ifdef TIMING
    CALL myDGSEM % timers % StopTimer( 3 )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! The boundary solutions are USEd to calculate the external states that, when
! accompanied with a Riemann Solver, enforce boundary conditions. Calling this
! routine is dependent on the result of CalculateBoundarySolution

#ifdef TIMING
    CALL myDGSEM % timers % StartTimer( 4 )
#endif

    CALL myDGSEM % UpdateExternalState( tn )

#ifdef TIMING
    CALL myDGSEM % timers % StopTimer( 4 )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! The inviscid fluxes (advection + pressure) through the faces of each element
! are estimated here using a (linear) Lax-Friedrich's upwind solver. In order to
! call this routine, CalculateBoundarySolution, UpdateExternalState, and
! MPI_StateExchange must have been completed.

#ifdef TIMING
    CALL myDGSEM % timers % StartTimer( 5 )
#endif

    CALL myDGSEM % InternalFace_StateFlux( )

#ifdef TIMING
    CALL myDGSEM % timers % StopTimer( 5 )
#endif

#ifdef TIMING
    CALL myDGSEM % timers % StartTimer( 6 )
#endif

#ifdef HAVE_MPI
    CALL myDGSEM % mpiStateHandler % Finalize_MPI_Exchange( myDGSEM % state, &
                                                            myDGSEM % mesh % faces, &
                                                            myDGSEM % extComm )
#endif

#ifdef TIMING
    CALL myDGSEM % timers % StopTimer( 6 )
#endif


#ifdef TIMING
    CALL myDGSEM % timers % StartTimer( 7 )
#endif

    CALL myDGSEM % BoundaryFace_StateFlux( )

#ifdef TIMING
    CALL myDGSEM % timers % StopTimer( 7 )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! IF the SpectralEKE or Laplacian subgridscale models are USEd, a Laplacian-like
! operator is USEd to dIFfUSE momentum and heat. When the Laplacian model is
! USEd a fixed viscosity is specIFied in runtime.params and a Rayleigh number of
! 1 is assumed (viscosity = dIFfusivity). When the SpectralEKE model is USEd,
! the viscosity coefficient is diagnosed using a Smagorinksy closure, where the
! unresolved kinetic energy is diagnosed from a highpass spectral filter.
!
    IF( myDGSEM % params % SubGridModel == SpectralEKE .OR. &
        myDGSEM % params % SubGridModel == Laplacian ) THEN

      IF( myDGSEM % params % SubGridModel == SpectralEKE )THEN !
!
! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! In the Spectral-EKE model, the under-resolved state is diagnosed from a high
! pass spectral filter (the dIFference between the state and smoothed state).
! Here, we first calculate the smoothed state and store it in the smoothedState
! attribute. This SUBROUTINE call has no dependence to any other within this
! SUBROUTINE.

#ifdef TIMING
        CALL myDGSEM % timers % StartTimer( 8 )
#endif

#ifdef HAVE_CUDA
        CALL myDGSEM % filter % Filter3D( myDGSEM % state % solution_dev, &
                                          myDGSEM % smoothState % solution_dev, &
                                          myDGSEM % state % nEquations_dev, &
                                          myDGSEM % mesh % elements % nElements_dev )
#else
        CALL myDGSEM % filter % Filter3D( myDGSEM % state % solution, &
                                          myDGSEM % smoothState % solution, &
                                          myDGSEM % state % nEquations, &
                                          myDGSEM % mesh % elements % nElements )
#endif

#ifdef TIMING
        CALL myDGSEM % timers % StopTimer( 8 )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
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

#ifdef TIMING
        CALL myDGSEM % timers % StartTimer( 9 )
#endif

        CALL myDGSEM % CalculateSGSCoefficients( )

#ifdef TIMING
        CALL myDGSEM % timers % StopTimer( 9 )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
! The viscosity coefficient that is calculated is now interpolated to the faces
! of each element so that the viscous flux can later be computed. This routine
! depends on the result of CalculateSGSCoefficients.

#ifdef TIMING
        CALL myDGSEM % timers % StartTimer( 10 )
#endif

#ifdef HAVE_CUDA

        tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                      4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                      1 )
        grid = dim3(myDGSEM % sgsCoeffs % nEquations, myDGSEM % sgsCoeffs % nElements, 1) 
         
        CALL CalculateSGSAtBoundaries_3D_CUDAKernel<<<grid, tBlock>>>( myDGSEM % sgsCoeffs % solution_dev, &
                                                                       myDGSEM % sgsCoeffs % boundarySolution_dev,  &
                                                                       myDGSEM % dgStorage % boundaryInterpolationMatrix_dev )

#else

        CALL myDGSEM % sgsCoeffs % Calculate_Solution_At_Boundaries( myDGSEM % dgStorage )

#endif

#ifdef TIMING
        CALL myDGSEM % timers % StopTimer( 10 )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! The viscosity coefficients are exchanged with neighboring ranks that share
! COMMON faces. MPI_SGSExchange can be run simulataneously with
! CalculateStressTensor, CalculateBoundaryStress, UpdateExternalStress, and the
! MPI_StressExchange. The viscosity coefficients that are exchanged are not
! needed until StressFlux
!

#ifdef TIMING
        CALL myDGSEM % timers % StartTimer( 11 )
#endif

#ifdef HAVE_MPI
        CALL myDGSEM % mpiSGSHandler % MPI_Exchange( myDGSEM % sgsCoeffs, &
                                                     myDGSEM % mesh % faces, &
                                                     myDGSEM % extComm )
#endif

        CALL myDGSEM % timers % StopTimer( 11 )
      ENDIF
! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! Now, the internal solution and the boundaryGradientFlux can be pieced
! together to calculate gradients in the velocity and potential temperature.
! This routine depends on the result of FaceFlux ( state % boundaryGradientFlux )

#ifdef TIMING
      CALL myDGSEM % timers % StartTimer( 12 )
#endif

      CALL myDGSEM % CalculateSolutionGradient( )

#ifdef TIMING
      CALL myDGSEM % timers % StopTimer( 12 )
#endif

#ifdef TIMING
      CALL myDGSEM % timers % StartTimer( 13 )
#endif

#ifdef HAVE_MPI
      IF( myDGSEM % params % SubGridModel == SpectralEKE )THEN !
        CALL myDGSEM % mpiSGSHandler % Finalize_MPI_Exchange( myDGSEM % sgsCoeffs, &
                                                              myDGSEM % mesh % faces, &
                                                              myDGSEM % extComm )
      ENDIF
#endif

#ifdef TIMING
      CALL myDGSEM % timers % StopTimer( 13 )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! The solution gradient values are interpolated to the faces of each element and
! projected onto the face normal direction. This
! routine depends on the result of CalculateStressTensor.

#ifdef TIMING
      CALL myDGSEM % timers % StartTimer( 14 )
#endif

      CALL myDGSEM % CalculateNormalStressAtBoundaries( )

#ifdef TIMING
      CALL myDGSEM % timers % StopTimer( 14 )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! Stress tensor values are exchanged with neighboring ranks along shared faces.
! This routine depends on the result of CalculateBoundaryStress, but can be run
! at the same time as UpdateExternalStress.

#ifdef TIMING
      CALL myDGSEM % timers % StartTimer( 15 )
#endif

#ifdef HAVE_MPI
      CALL myDGSEM % mpiStressHandler % MPI_Exchange( myDGSEM % stressTensor, &
                                                      myDGSEM % mesh % faces, &
                                                      myDGSEM % extComm )
#endif

#ifdef TIMING
      CALL myDGSEM % timers % StopTimer( 15 )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
!  Using the solution gradient and the eddy-viscosities/diffusivities, the
!  stress flux is calculated

#ifdef TIMING
      CALL myDGSEM % timers % StartTimer( 16 )
#endif

      CALL myDGSEM % CalculateStressFlux( )

#ifdef TIMING
      CALL myDGSEM % timers % StopTimer( 16 )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! Now that the stress tensor is available on element faces, boundary conditions
! can be applied by setting the external stress tensor state. This routine
! depends on the result of CalculateBoundaryStress. Note that this routine can
! be run simultaneously with the MPI_StressExchange

#ifdef TIMING
      CALL myDGSEM % timers % StartTimer( 17 )
#endif

      CALL myDGSEM % UpdateExternalStress( tn )

#ifdef TIMING
      CALL myDGSEM % timers % StopTimer( 17 )
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

#ifdef TIMING
      CALL myDGSEM % timers % StartTimer( 18 )
#endif

#ifdef HAVE_MPI
      CALL myDGSEM % mpiStressHandler % Finalize_MPI_Exchange( myDGSEM % stressTensor,&
                                                               myDGSEM % mesh % faces, &
                                                               myDGSEM % extComm )
#endif

#ifdef TIMING
      CALL myDGSEM % timers % StopTimer( 18 )
#endif

#ifdef TIMING
      CALL myDGSEM % timers % StartTimer( 19 )
#endif

      CALL myDGSEM % BoundaryFace_StressFlux( )
      
#ifdef TIMING
      CALL myDGSEM % timers % StopTimer( 19 )
#endif      

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><>><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! With the boundary stress flux and the internal stress tensor values, the
! divergence of the stress tensor can be calculated, giving the viscous tendency
! for the momentum and the potential temperature. This routine depends on the
! result of StressFlux (and the dependencies of StressFlux), but can be DOne
! simultaneously with the MappedTimeDerivative.
!

#ifdef TIMING
      CALL myDGSEM % timers % StartTimer( 20 )
#endif

#ifdef HAVE_CUDA
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
      CALL myDGSEM % stressTensor % Calculate_Weak_Flux_Divergence( myDGSEM % dgStorage )
#endif

#ifdef TIMING
      CALL myDGSEM % timers % StopTimer( 20 )
#endif


    ENDIF


! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! Once the inviscid fluxes through the faces are calculated, and the internal
! state is known, the tendency due to the inviscid flux terms and
! nonconservative source terms is calculated here. This routine depends on the
! result of FaceFlux, but can be DOne at the same time as StressDivergence

#ifdef TIMING
      CALL myDGSEM % timers % StartTimer( 21 )
#endif

    CALL myDGSEM % MappedTimeDerivative( )

#ifdef TIMING
      CALL myDGSEM % timers % StopTimer( 21 )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !


  END SUBROUTINE GlobalTimeDerivative_Fluid
!
  SUBROUTINE CalculateSGSCoefficients_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1 , &
                  myDGSEM % params % polyDeg+1 )
    grid = dim3(myDGSEM % mesh % elements % nElements, 1, 1)

    CALL CalculateSGSCoefficients_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                               myDGSEM % static % solution_dev, &
                                                               myDGSEM % smoothState % solution_dev, &
                                                               myDGSEM % filter % filterMat_dev, &
                                                               myDGSEM % sgsCoeffs % solution_dev )
#else
    ! Local
    INTEGER :: iEl, i, j, k, m, ii, jj, kk
    REAL(prec) :: sgsKE, uijk, uij, ui
    REAL(prec) :: KE(0:myDGSEM % params % polyDeg,0:myDGSEM % params % polyDeg,0:myDGSEM % params % polyDeg)


    DO iEl = 1, myDGSEM % mesh % elements % nElements

      ! Here, the SGS Kinetic energy is calculated using the
      ! "high wavenumber" component of the velocity field.
      ! This component is defined (here) as the dIFference
      ! between the full solution and the smoothed solution.
      DO k = 0, myDGSEM % params % polyDeg
        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg

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

      ! Now we calculate the viscosity and dIFfusivities (currently assumes isotropic and low mach number)
      DO k = 0, myDGSEM % params % polyDeg
        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg
            DO m = 1, myDGSEM % state % nEquations-1

              myDGSEM % sgsCoeffs % solution(i,j,k,m,iEl) = 0.09_prec*&
                myDGSEM % params % viscLengthScale*sqrt( KE(i,j,k) )

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

#endif

  END SUBROUTINE CalculateSGSCoefficients_Fluid
!
  SUBROUTINE UpdateExternalSGS_Fluid( myDGSEM ) ! ////////// !

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1, &
                  1 )
    grid = dim3(myDGSEM % extComm % nBoundaries,myDGSEM % sgsCoeffs % nEquations,1)

    CALL UpdateExternalSGSCoeffs_CUDAKernel<<<grid, tBlock>>>( myDGSEM % extComm % boundaryIDs_dev, &       ! I
                                                               myDGSEM % mesh % faces % elementIDs_dev, &   ! I
                                                               myDGSEM % mesh % faces % elementSides_dev, & ! I
                                                               myDGSEM % extComm % extProcIDs_dev, &           ! I
                                                               myDGSEM % sgsCoeffs % externalState_dev, &                    ! O
                                                               myDGSEM % sgsCoeffs % boundarySolution_dev, &   ! I
                                                               myDGSEM % mesh % elements % nHat_dev )
#else
    ! Local
    INTEGER    :: iEl, bID, bFaceID, i, j, k, iEq
    INTEGER    :: iFace2, p2
    INTEGER    :: e1, e2, s1, s2

    DO bID = 1, myDGSEM % extComm % nBoundaries

      iFace2 = myDGSEM % extComm % boundaryIDs( bID ) ! Obtain the process-local face id for this boundary-face id
      e1     = myDGSEM % mesh % faces % elementIDs(1,iFace2)
      s1     = myDGSEM % mesh % faces % elementSides(1,iFace2)
      e2     = myDGSEM % mesh % faces % elementIDs(2,iFace2)
      p2     = myDGSEM % extComm % extProcIDs( bID )

      IF( p2 == myDGSEM % extComm % myRank )THEN ! Enforce no boundary flux due to the fluid stress

        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg
            DO iEq = 1, myDGSEM % state % nEquations-1
              myDGSEM % sgsCoeffs % externalState(i,j,iEq,bID) = myDGSEM % sgsCoeffs % boundarySolution(i,j,iEq,s1,e1)
            ENDDO
          ENDDO
        ENDDO

      ENDIF

    ENDDO

#endif

  END SUBROUTINE UpdateExternalSGS_Fluid
!
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!                                                                                                 !
!  This section of code CONTAINS routines for applying boundary conditions along physical         !
!  boundaries.                                                                                    !
!                                                                                                 !
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!
  SUBROUTINE UpdateExternalState_Fluid( myDGSEM, tn )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
    REAL(prec), INTENT(in)      :: tn
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1, &
                  1 )
    grid = dim3(myDGSEM % extComm % nBoundaries,1,1)

    CALL UpdateExternalState_CUDAKernel<<<grid, tBlock>>>( myDGSEM % extComm % boundaryIDs_dev, &       ! I
                                                           myDGSEM % mesh % faces % elementIDs_dev, &   ! I
                                                           myDGSEM % mesh % faces % elementSides_dev, & ! I
                                                           myDGSEM % extComm % extProcIDs_dev, &           ! I
                                                           myDGSEM % state % externalState_dev, &               ! O
                                                           myDGSEM % state % boundarySolution_dev, &    ! I
                                                           myDGSEM % state % prescribedState_dev, &             ! I
                                                           myDGSEM % mesh % elements % nHat_dev )
#else
    ! Local
    INTEGER    :: iEl, bID, bFaceID, i, j, k, iEq
    INTEGER    :: iFace2, p2
    INTEGER    :: e1, e2, s1, s2
    REAL(prec) :: norm, un, ut, us, speed
    REAL(prec) :: nx, ny, nz
    REAL(prec) :: sx, sy, sz
    REAL(prec) :: tx, ty, tz

    DO bID = 1, myDGSEM % extComm % nBoundaries

      iFace2 = myDGSEM % extComm % boundaryIDs( bID ) ! Obtain the process-local face id for this boundary-face id
      p2     = myDGSEM % extComm % extProcIDs( bID )
      e1     = myDGSEM % mesh % Faces % elementIDs(1,iFace2)
      s1     = myDGSEM % mesh % Faces % elementSides(1,iFace2)
      e2     = myDGSEM % mesh % Faces % elementIDs(2,iFace2)

      DO j = 0, myDGSEM % params % polyDeg
        DO i = 0, myDGSEM % params % polyDeg

          IF( e2 == PRESCRIBED .AND. p2 == myDGSEM % extComm % myRank )THEN
            DO iEq = 1, myDGSEM % state % nEquations
              myDGSEM % state % externalState(i,j,iEq,bID) = myDGSEM % state % prescribedState(i,j,iEq,bID)
            ENDDO
          ELSEIF( e2 == RADIATION .AND. p2 == myDGSEM % extComm % myRank )THEN

            DO iEq = 1, myDGSEM % state % nEquations
              myDGSEM % state % externalState(i,j,iEq,bID) = 0.0_prec
            ENDDO

          ELSEIF( e2 == NO_NORMAL_FLOW .AND. p2 == myDGSEM % extComm % myRank )THEN

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

#ifdef PASSIVE_TRACERS
            myDGSEM % state % externalState(i,j,6,bID) =  myDGSEM % state % boundarySolution(i,j,6,s1,e1) 
#endif

            myDGSEM % state % externalState(i,j,nEquations,bID) =  myDGSEM % state % boundarySolution(i,j,nEquations,s1,e1) ! P

          ELSEIF( e2 == INFLOW .AND. p2 == myDGSEM % extComm % myRank )THEN

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
#ifdef PASSIVE_TRACERS
            myDGSEM % state % externalState(i,j,6,bID) =  myDGSEM % state % prescribedState(i,j,6,bID)
#endif
            myDGSEM % state % externalState(i,j,nEquations,bID) =  myDGSEM % state % prescribedState(i,j,nEquations,bID) ! P

          ENDIF

        ENDDO
      ENDDO
    ENDDO

#endif

  END SUBROUTINE UpdateExternalState_Fluid
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
    INTEGER :: iEl, iFace
    INTEGER    :: i, j, k, m, iEq
    INTEGER    :: ii, jj, bID
    INTEGER    :: e1, s1, e2, s2
    REAL(prec) :: nHat(1:3), norm
    REAL(prec) :: uOut, uIn, cIn, cOut, T
    REAL(prec) :: jump(1:myDGSEM % state % nEquations-1), aS(1:myDGSEM % state % nEquations-1)
    REAL(prec) :: fac, rC

    DO iFace = 1, myDGSEM % mesh % faces % nFaces


      e1 = myDGSEM % mesh % faces % elementIDs(1,iFace)
      s1 = myDGSEM % mesh % faces % elementSides(1,iFace)
      e2 = myDGSEM % mesh % faces % elementIDs(2,iFace)
      s2 = ABS(myDGSEM % mesh % faces % elementSides(2,iFace))

      IF( e2 > 0 )THEN

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
            cOut = sqrt( myDGSEM % params % R *T* &
              ( (myDGSEM % state % boundarySolution(ii,jj,nEquations,s2,e2)+&
              myDGSEM % static % boundarySolution(ii,jj,nEquations,s2,e2))/&
              myDGSEM % params % P0 )**myDGSEM % params % rC   )

            T =   (myDGSEM % static % boundarySolution(i,j,5,s1,e1) + &
              myDGSEM % state % boundarySolution(i,j,5,s1,e1))/&
              (myDGSEM % static % boundarySolution(i,j,4,s1,e1) + &
              myDGSEM % state % boundarySolution(i,j,4,s1,e1) )

            cIn  = sqrt( myDGSEM % params % R*T* &
              ( (myDGSEM % state % boundarySolution(i,j,nEquations,s1,e1)+&
              myDGSEM % static % boundarySolution(i,j,nEquations,s1,e1))/&
              myDGSEM % params % P0 )**myDGSEM % params % rC  )

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

            bID  = ABS(myDGSEM % mesh % faces % boundaryID(iFace))
            DO iEq = 1, myDGSEM % state % nEquations-1
              jump(iEq)  = myDGSEM % state % externalState(ii,jj,iEq,bID) - &
                           myDGSEM % state % boundarySolution(i,j,iEq,s1,e1)
            ENDDO

            ! Sound speed estimate for the external and internal states

            T = (myDGSEM % static % boundarySolution(i,j,5,s1,e1)+myDGSEM % state % externalState(ii,jj,5,bID))/&
              (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % externalState(ii,jj,4,bID))

            cOut = sqrt( myDGSEM % params % R*T* &
              ( (myDGSEM % state % externalState(ii,jj,nEquations,bID)+&
              myDGSEM % static % boundarySolution(i,j,nEquations,s1,e1) )/&
              myDGSEM % params % P0 )**myDGSEM % params % rC   )

            T = (myDGSEM % static % boundarySolution(i,j,5,s1,e1)+myDGSEM % state % boundarySolution(i,j,5,s1,e1))/&
                (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % boundarySolution(i,j,4,s1,e1))

            cIn  = sqrt( myDGSEM % params % R*T* &
              ( (myDGSEM % state % boundarySolution(i,j,nEquations,s1,e1)+&
              myDGSEM % static % boundarySolution(i,j,nEquations,s1,e1) )/&
              myDGSEM % params % P0 )**myDGSEM % params % rC  )

            ! External normal velocity component
            uOut = ( myDGSEM % state % externalState(ii,jj,1,bID)*nHat(1) + &
                     myDGSEM % state % externalState(ii,jj,2,bID)*nHat(2) + &
                     myDGSEM % state % externalState(ii,jj,3,bID)*nHat(3) )/&
                   ( myDGSEM % state % externalState(ii,jj,4,bID)+&
                     myDGSEM % static % boundarySolution(i,j,4,s1,e1) )

            ! Internal normal velocity component
            uIn  = ( myDGSEM % state % boundarySolution(i,j,1,s1,e1)*nHat(1) + &
                     myDGSEM % state % boundarySolution(i,j,2,s1,e1)*nHat(2) + &
                     myDGSEM % state % boundarySolution(i,j,3,s1,e1)*nHat(3) )/&
                   ( myDGSEM % state % boundarySolution(i,j,4,s1,e1)+&
                     myDGSEM % static % boundarySolution(i,j,4,s1,e1) )

            fac = max( abs(uIn+cIn), abs(uIn-cIn), abs(uOut+cOut), abs(uOut-cOut) )

            ! Advective flux
            DO iEq = 1, myDGSEM % state % nEquations-1
              aS(iEq) = uIn*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) + myDGSEM % static % boundarySolution(i,j,iEq,s1,e1) ) +&
                uOut*( myDGSEM % state % externalState(ii,jj,iEq,bID) + myDGSEM % static % boundarySolution(i,j,iEq,s1,e1) )
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
                                                                                     myDGSEM % static % boundarySolution(i,j,4,s1,e1)) )*&
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
    INTEGER    :: iEl, i, j, k, m, iEq, row, col
    REAL(prec) :: F
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

     CALL myDGSEM % state % Calculate_Weak_Flux_Divergence( myDGSEM % dgStorage )  

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
    grid = dim3(myDGSEM % sgsCoeffs % nEquations,myDGSEM % mesh % elements % nElements,3)

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

    !$OMP PARALLEL
    !$OMP DO PRIVATE( f, df )
    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO iEq = 1, myDGSEM % sgsCoeffs % nEquations

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
    !$OMP ENDDO
    !$OMP END PARALLEL

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
                                                                           myDGSEM % sgsCoeffs % boundarySolution_dev,  &
                                                                           myDGSEM % mesh % elements % nHat_dev, &
                                                                           myDGSEM % stressTensor % boundarySolution_dev,  &
                                                                           myDGSEM % dgStorage % boundaryInterpolationMatrix_dev )

#else


    !$OMP PARALLEL
    !$OMP DO PRIVATE( fAtBoundaries )
    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO iEq = 1, myDGSEM % sgsCoeffs % nEquations

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

                myDGSEM % stressTensor % boundarySolution(i,j,iEq,k,iEl) = myDGSEM % sgsCoeffs % boundarySolution(i,j,iEq,k,iEl)*&
                                                                           ( fAtBoundaries(1,k)*myDGSEM % mesh % elements % nHat(1,i,j,k,iEL) + &
                                                                             fAtBoundaries(2,k)*myDGSEM % mesh % elements % nHat(2,i,j,k,iEL) + &
                                                                             fAtBoundaries(3,k)*myDGSEM % mesh % elements % nHat(3,i,j,k,iEL) )
          
              
              ENDDO

            ENDDO
          ENDDO

      ENDDO
    ENDDO
    !$OMP ENDDO
    !$OMP END PARALLEL
  

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
                                                          myDGSEM % sgsCoeffs % solution_dev, &
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

    !$OMP PARALLEL
    !$OMP DO PRIVATE( flux )
    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO iEq = 1, myDGSEM % sgsCoeffs % nEquations
        DO k = 0, myDGSEM % params % polyDeg
          DO j = 0, myDGSEM % params % polyDeg
            DO i = 0, myDGSEM % params % polyDeg

              flux(1:3) = 0.0_prec

              DO idir = 1, 3


                flux(1) = flux(1) + myDGSEM % mesh % elements % Ja(i,j,k,idir,1,iEl)*&
                                          myDGSEM % state % solutionGradient(idir,i,j,k,iEq,iEl)*&
                                          myDGSEM % sgsCoeffs % solution(i,j,k,iEq,iEl)

                flux(2) = flux(2) + myDGSEM % mesh % elements % Ja(i,j,k,idir,2,iEl)*&
                                          myDGSEM % state % solutionGradient(idir,i,j,k,iEq,iEl)*&
                                          myDGSEM % sgsCoeffs % solution(i,j,k,iEq,iEl) 

                flux(3) = flux(3) + myDGSEM % mesh % elements % Ja(i,j,k,idir,3,iEl)*&
                                          myDGSEM % state % solutionGradient(idir,i,j,k,iEq,iEl)*&
                                          myDGSEM % sgsCoeffs % solution(i,j,k,iEq,iEl) 
              ENDDO

              myDGSEM % stressTensor % flux(1:3,i,j,k,iEq,iEl) = flux(1:3)

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !$OMP ENDDO
    !$OMP END PARALLEL

#endif

  END SUBROUTINE CalculateStressFlux_Fluid
!
!
  SUBROUTINE UpdateExternalStress_Fluid( myDGSEM, tn ) ! ////////// !

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
    REAL(prec), INTENT(in)      :: tn
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1, &
                  1 )
    grid = dim3(myDGSEM % extComm % nBoundaries,myDGSEM % stressTensor % nEquations,1)

    CALL UpdateExternalStress_CUDAKernel<<<grid, tBlock>>>( myDGSEM % extComm % boundaryIDs_dev, &            ! I
                                                            myDGSEM % mesh % faces % elementIDs_dev, &         ! I
                                                            myDGSEM % mesh % faces % elementSides_dev, &       ! I
                                                            myDGSEM % extComm % extProcIDs_dev, &              ! I
                                                            myDGSEM % stressTensor % externalState_dev, &                    ! O
                                                            myDGSEM % stressTensor % boundarySolution_dev, &   ! I
                                                            myDGSEM % stressTensor % prescribedState_dev )
#else
    ! Local
    INTEGER    :: iEl, iFace, bFaceID, i, j, k, iEq
    INTEGER    :: bID, p2
    INTEGER    :: e1, e2, s1, s2

    DO bID = 1, myDGSEM % extComm % nBoundaries

      iFace = myDGSEM % extComm % boundaryIDs( bID ) ! Obtain the process-local face id for this boundary-face id
      e1    = myDGSEM % mesh % faces % elementIDs(1,iFace)
      s1    = myDGSEM % mesh % faces % elementSides(1,iFace)
      e2    = myDGSEM % mesh % faces % elementIDs(2,iFace)
      p2    = myDGSEM % extComm % extProcIDs( bID )

      IF( p2 == myDGSEM % extComm % myRank )THEN ! Enforce no boundary flux due to the fluid stress
        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg
            DO iEq = 1, myDGSEM % stressTensor % nEquations
              myDGSEM % stressTensor % externalState(i,j,iEq,bID) = myDGSEM % stressTensor % boundarySolution(i,j,iEq,s1,e1)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

    ENDDO

#endif

  END SUBROUTINE UpdateExternalStress_Fluid
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
                                                               myDGSEM % sgsCoeffs % boundarySolution_dev, &
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
    INTEGER    :: e1, s1, e2, s2

    DO iFace = 1, myDGSEM % mesh % faces % nFaces


      e1 = myDGSEM % mesh % faces % elementIDs(1,iFace)
      s1 = myDGSEM % mesh % faces % elementSides(1,iFace)
      e2 = myDGSEM % mesh % faces % elementIDs(2,iFace)
      s2 = ABS(myDGSEM % mesh % faces % elementSides(2,iFace))

      IF( e2 > 0 )THEN

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
                0.5_prec*( myDGSEM % sgsCoeffs % boundarysolution(ii,jj,iEq,s2,e2)*myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2)-&
                           myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) )/&
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
                                                               myDGSEM % sgsCoeffs % boundarySolution_dev, &
                                                               myDGSEM % sgsCoeffs % externalState_dev, &
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


      e1 = myDGSEM % mesh % faces % elementIDs(1,iFace)
      s1 = myDGSEM % mesh % faces % elementSides(1,iFace)
      e2 = myDGSEM % mesh % faces % elementIDs(2,iFace)
      s2 = ABS(myDGSEM % mesh % faces % elementSides(2,iFace))

      IF( e2 < 0 )THEN

        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg


            ii = myDGSEM % mesh % faces % iMap(i,j,iFace)
            jj = myDGSEM % mesh % faces % jMap(i,j,iFace)

            norm = sqrt( myDGSEM % mesh % elements % nHat(1,i,j,s1,e1)**2 + &
                         myDGSEM % mesh % elements % nHat(2,i,j,s1,e1)**2 + &
                         myDGSEM % mesh % elements % nHat(3,i,j,s1,e1)**2 )

            bID  = ABS(myDGSEM % mesh % faces % boundaryID(iFace))

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
                ( myDGSEM % sgsCoeffs % externalState(ii,jj,iEq,bID)*myDGSEM % state % externalState(ii,jj,iEq,bID)-&
                  myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) )/&
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
    REAL(prec) :: hCapRatio, rC, rhoT


    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO k = 0, myDGSEM % params % polyDeg
        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg

            ! Pressure = rho*e*R/Cv (Ideal Gas Law, P = rho*R*T, thermo ~> T = theta*(P/P0)^r)
            ! THEN P = (rho*theta*R/P0^rC)^(Cp/Cv)
            ! And P' = P - P_static
            rhoT = (myDGSEM % static % solution(i,j,k,5,iEl) + myDGSEM % state % solution(i,j,k,5,iEl) )
            myDGSEM % state % solution(i,j,k,nEquations,iEl) = myDGSEM % params % P0*( rhoT*myDGSEM % params % R/myDGSEM % params % P0 )**myDGSEM % params % hCapRatio -&
                                                      myDGSEM % static % solution(i,j,k,nEquations,iEl)

          ENDDO
        ENDDO
      ENDDO
    ENDDO

#endif

  END SUBROUTINE EquationOfState_Fluid
!
  SUBROUTINE CalculateStaticState_Fluid( myDGSEM )
    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM
    ! Local
    INTEGER    :: i, j, k, iEl,iEq
    REAL(prec) :: z, H, P0, Cp, T, T0, dTdz, P, rC, g, R
#ifdef HAVE_CUDA
    INTEGER :: istat
#endif

    R    = myDGSEM % params % R
    Cp   = (R + myDGSEM % params % Cv)
    rC   = R/Cp
    g    = myDGSEM % params % g
    H    = myDGSEM % params % zScale
    T0   = myDGSEM % params % T0
    P0   = myDGSEM % params % P0
    dTdz = myDGSEM % params % dTdz

    ! /////////////////////  Build the Static/Background State ///////////////////////// !

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

            ! The static profile is determined from hydrostatic balance, the equation of state,
            ! and a prescribed potential temperature profile.
            ! ** The potential temperature is assumed to vary linearly with z **

            T = T0 + dTdz*z ! Potential temperature
            IF( AlmostEqual( dTdz, 0.0_prec ) )THEN
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

#ifdef HAVE_CUDA
    myDGSEM % static % solution_dev = myDGSEM % static % solution
#endif

    ! This routine Calculates the pressure
    CALL myDGSEM % EquationOfState( )

#ifdef HAVE_CUDA
    myDGSEM % state % solution = myDGSEM % state % solution_dev
#endif

    DO iEl = 1, myDGSEM % mesh % elements % nElements
      myDGSEM % static % solution(:,:,:,nEquations,iEl) = myDGSEM % state % solution(:,:,:,nEquations,iEl)
    ENDDO

#ifdef HAVE_CUDA
    myDGSEM % static % solution_dev = myDGSEM % static % solution
    myDGSEM % state % solution_dev  = 0.0_prec
    istat = cudaDeviceSynchronize( )
#endif

    myDGSEM % state % solution = 0.0_prec

      PRINT*,   MINVAL( myDGSEM % static % solution(:,:,:,4,:) ), MAXVAL( myDGSEM % static % solution(:,:,:,4,:) )

  END SUBROUTINE CalculateStaticState_Fluid
!
  SUBROUTINE WriteTecplot_Fluid( myDGSEM )

    IMPLICIT NONE

    CLASS( Fluid ), INTENT(inout) :: myDGsem
    !LOCAL
    REAL(prec)  :: x(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:3,1:myDGSEM % mesh % elements % nElements)
    REAL(prec)  :: sol(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:myDGSEM % state % nEquations,1:myDGSEM % mesh % elements % nElements)
    REAL(prec)  :: bsol(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:myDGSEM % state % nEquations, 1:myDGSEM % mesh % elements % nElements)
    REAL(prec)  :: drag(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1, 1:myDGSEM % mesh % elements % nElements)
    REAL(prec)  :: drag_init(0:myDGSEM % params % polydeg,0:myDGSEM % params % polydeg,0:myDGSEM % params % polydeg,1, 1:myDGSEM % mesh % elements % nElements)
#ifdef HAVE_CUDA
    INTEGER :: istat
#endif
    INTEGER       :: i, j, k, iEl, iEq, fUnit
    CHARACTER(5)  :: zoneID
    CHARACTER(4)  :: rankChar
    REAL(prec)    :: hCapRatio, c, T
    CHARACTER(13) :: timeStampString

    timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )

#ifdef HAVE_CUDA
    CALL myDGSEM % state % UpdateHost( )
    CALL myDGSEM % static % UpdateHost( )
    istat = cudaDeviceSynchronize( )

#endif

    drag_init(:,:,:,1,:) = myDGSEM % sourceTerms % drag
    sol = ApplyInterpolationMatrix_3D_Lagrange( myDGSEM % dgStorage % interp, &
                                                myDGSEM % state % solution, &
                                                myDGSEM % state % nEquations, &
                                                myDGSEM % mesh % elements % nElements )

    bsol = ApplyInterpolationMatrix_3D_Lagrange( myDGSEM % dgStorage % interp, myDGSEM % static % solution, &
                                                 myDGSEM % static % nEquations, &
                                                 myDGSEM % mesh % elements % nElements )

    drag = ApplyInterpolationMatrix_3D_Lagrange( myDGSEM % dgStorage % interp, drag_init, &
                                                 1, &
                                                 myDGSEM % mesh % elements % nElements )

    x = ApplyInterpolationMatrix_3D_Lagrange( myDGSEM % dgStorage % interp, myDGSEM % mesh % elements % x, &
                                              3, &
                                              myDGSEM % mesh % elements % nElements )

    WRITE(rankChar,'(I4.4)') myDGSEM % extComm % myRank

    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= 'State.'//rankChar//'.'//timeStampString//'.tec', &
      FORM='formatted', &
      STATUS='replace')

#ifdef PASSIVE_TRACERS

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "u", "v", "w", "rho", "Pot. Temp.", "Tracer", "Pressure",'//&
      ' "rho_b", "Pot. Temp._b", "Pressure_b", "Drag", "c" '

#else

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "u", "v", "w", "rho", "Pot. Temp.", "Pressure",'//&
      ' "rho_b", "Pot. Temp._b", "Pressure_b", "Drag", "c" '

#endif

    DO iEl = 1, myDGsem % mesh % elements % nElements

      WRITE(zoneID,'(I5.5)') myDGSEM % mesh % elements % elementID(iEl)
      WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',myDGSEM % params % nPlot+1,&
                                                 ', J=',myDGSEM % params % nPlot+1,&
                                                 ', K=',myDGSEM % params % nPlot+1,',F=POINT'

      DO k = 0, myDGSEM % params % nPlot
        DO j = 0, myDGSEM % params % nPlot
          DO i = 0, myDGSEM % params % nPlot
            T =   (bsol(i,j,k,5,iEl) + sol(i,j,k,5,iEl))/(bsol(i,j,k,4,iEl)+sol(i,j,k,4,iEl) )

            ! Sound speed estimate for the external and internal states
            c = sqrt( myDGSEM % params % R*T*( ( sol(i,j,k,nEquations,iEl) + bsol(i,j,k,nEquations,iEl) )/myDGSEM % params % P0 )**myDGSEM % params % hCapRatio   )
            WRITE(fUnit,'(17(E15.7,1x))') x(i,j,k,1,iEl), x(i,j,k,2,iEl), x(i,j,k,3,iEl),&
              sol(i,j,k,1,iEl)/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ), &
              sol(i,j,k,2,iEl)/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ), &
              sol(i,j,k,3,iEl)/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ), &
              sol(i,j,k,4,iEl), &
              (sol(i,j,k,5,iEl) + bsol(i,j,k,5,iEl))/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ), &
#ifdef PASSIVE_TRACERS
              ( bsol(i,j,k,6,iEl) + sol(i,j,k,6,iEl) )/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ), &
#endif
              sol(i,j,k,nEquations,iEl), &
              bsol(i,j,k,4,iEl), &
              bsol(i,j,k,5,iEl)/( bsol(i,j,k,4,iEl) ),&
              bsol(i,j,k,nEquations,iEl),&
              drag(i,j,k,1,iEl), c


          ENDDO
        ENDDO
      ENDDO

    ENDDO

    CLOSE(UNIT=fUnit)

  END SUBROUTINE WriteTecplot_Fluid
!
  SUBROUTINE OpenDiagnosticsFiles_Fluid( myDGSEM )
    IMPLICIT NONE
    CLASS( Fluid  ), INTENT(inout) :: myDGSEM
    ! Local
    CHARACTER(13) :: timeStampString


    myDGSEM % volume = 0.0_prec
    myDGSEM % mass   = 0.0_prec
    myDGSEM % KE     = 0.0_prec
    myDGSEM % PE     = 0.0_prec
    myDGSEM % heat   = 0.0_prec

    IF( myDGSEM % extComm % myRank == 0 )THEN
      timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )


      OPEN( UNIT=NewUnit(diagUnits(1)), &
        FILE='Mass.'//timeStampString//'.curve', &
        FORM='FORMATTED', &
        STATUS='REPLACE', &
        ACCESS='APPEND' )
      WRITE(diagUnits(1),*) '#TotalMass'

      OPEN( UNIT=NewUnit(diagUnits(2)), &
        FILE='KineticEnergy.'//timeStampString//'.curve', &
        FORM='FORMATTED', &
        STATUS='REPLACE', &
        ACCESS='APPEND' )
      WRITE(diagUnits(2),*) '#TotalKineticEnergy'

      OPEN( UNIT=NewUnit(diagUnits(3)), &
        FILE='PotentialEnergy.'//timeStampString//'.curve', &
        FORM='FORMATTED', &
        STATUS='REPLACE', &
        ACCESS='APPEND' )
      WRITE(diagUnits(3),*) '#TotalPotentialEnergy'

      OPEN( UNIT=NewUnit(diagUnits(4)), &
        FILE='Heat.'//timeStampString//'.curve', &
        FORM='FORMATTED', &
        STATUS='REPLACE', &
        ACCESS='APPEND' )
      WRITE(diagUnits(4),*) '#TotalHeat'

      OPEN( UNIT=NewUnit(diagUnits(5)), &
        FILE='Volume.'//timeStampString//'.curve', &
        FORM='FORMATTED', &
        STATUS='REPLACE', &
        ACCESS='APPEND' )
      WRITE(diagUnits(5),*) '#TotalVolume'
    ENDIF

  END SUBROUTINE OpenDiagnosticsFiles_Fluid
!
  SUBROUTINE WriteDiagnostics_Fluid( myDGSEM )
    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM

    CALL myDGSEM % OpenDiagnosticsFiles( )

    IF( myDGSEM % extComm % myRank == 0 )THEN
      WRITE(diagUnits(1),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % mass
      WRITE(diagUnits(2),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % KE
      WRITE(diagUnits(3),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % PE
      WRITE(diagUnits(4),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % heat
      WRITE(diagUnits(5),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % volume
    ENDIF

    CALL myDGSEM % CloseDiagnosticsFiles( )

  END SUBROUTINE WriteDiagnostics_Fluid
!
  SUBROUTINE CloseDiagnosticsFiles_Fluid( myDGSEM )
    IMPLICIT NONE
    CLASS( Fluid  ), INTENT(inout) :: myDGSEM


    IF( myDGSEM % extComm % myRank == 0 ) THEN
      CLOSE( UNIT=diagUnits(1) )
      CLOSE( UNIT=diagUnits(2) )
      CLOSE( UNIT=diagUnits(3) )
      CLOSE( UNIT=diagUnits(4) )
      CLOSE( UNIT=diagUnits(5) )
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
    CALL MPI_ALLREDUCE( volume, myDGSEM % volume, 1, myDGSEM % extComm % MPI_PREC, MPI_SUM, MPI_COMM_WORLD, mpiErr )
    CALL MPI_ALLREDUCE( mass, myDGSEM % mass, 1, myDGSEM % extComm % MPI_PREC, MPI_SUM, MPI_COMM_WORLD, mpiErr )
    CALL MPI_ALLREDUCE( KE, myDGSEM % KE, 1, myDGSEM % extComm % MPI_PREC, MPI_SUM, MPI_COMM_WORLD, mpiErr )
    CALL MPI_ALLREDUCE( PE, myDGSEM % PE, 1, myDGSEM % extComm % MPI_PREC, MPI_SUM, MPI_COMM_WORLD, mpiErr )
    CALL MPI_ALLREDUCE( heat, myDGSEM % heat, 1, myDGSEM % extComm % MPI_PREC, MPI_SUM, MPI_COMM_WORLD, mpiErr )
#endif

    IF( IsNaN( myDGSEM % KE ) .OR. IsInf( myDGSEM % KE/ myDGSEM % volume) )THEN

      PRINT*, '  Model bust at simulation time : ', myDGSEM % simulationTime
      PRINT*, '  Total Kinetic Energy :', myDGSEM % KE
      PRINT*, '  Consider reducing the time step or increasing dissipation.'
      PRINT*,   MINVAL( myDGSEM % static % solution(:,:,:,1,:) ), MAXVAL( myDGSEM % static % solution(:,:,:,1,:) )
      PRINT*,   MINVAL( myDGSEM % static % solution(:,:,:,2,:) ), MAXVAL( myDGSEM % static % solution(:,:,:,2,:) )
      PRINT*,   MINVAL( myDGSEM % static % solution(:,:,:,3,:) ), MAXVAL( myDGSEM % static % solution(:,:,:,3,:) )
      PRINT*,   MINVAL( myDGSEM % static % solution(:,:,:,4,:) ), MAXVAL( myDGSEM % static % solution(:,:,:,4,:) )
      PRINT*,   MINVAL( myDGSEM % static % solution(:,:,:,5,:) ), MAXVAL( myDGSEM % static % solution(:,:,:,5,:) )
      PRINT*,   MINVAL( myDGSEM % static % solution(:,:,:,6,:) ), MAXVAL( myDGSEM % static % solution(:,:,:,6,:) )

      CALL myDGSEM % WriteDiagnostics( )
      CALL myDGSEM % Trash( )

      STOP
    ENDIF

  END SUBROUTINE Diagnostics_Fluid
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

    N = myDGSEM % params % polyDeg

    WRITE(rankChar,'(I4.4)') myDGSEM % extComm % myRank
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
    DO iEl = 1, myDGSEM % mesh % elements % nElements

      DO iEq = 1, myDGSEM % state % nEquations
        WRITE( fUnit, REC=thisRec )myDGSEM % state % solution(:,:,:,iEq,iEl)
        thisRec = thisRec+1
      ENDDO
      DO iEq = 1, myDGSEM % state % nEquations
        WRITE( fUnit, REC=thisRec )myDGSEM % static % solution(:,:,:,iEq,iEl)
        thisRec = thisRec+1
      ENDDO
      DO iEq = 1, myDGSEM % state % nEquations
        WRITE( fUnit, REC=thisRec )myDGSEM % static % source(:,:,:,iEq,iEl)
        thisRec = thisRec+1
      ENDDO

      WRITE( fUnit, REC=thisRec )myDGSEM % sourceTerms % drag(:,:,:,iEl)
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
      RECL   = prec*(myDGSEM % params % polyDeg+1)*(myDGSEM % params % polyDeg+1)*(myDGSEM % state % nEquations)*(myDGSEM % extComm % nBoundaries) )
    WRITE( fUnit, rec = 1 ) myDGSEM % state % externalState
    WRITE( fUnit, rec = 2 ) myDGSEM % state % prescribedState
    CLOSE(fUnit)


  END SUBROUTINE WritePickup_Fluid
!
  SUBROUTINE ReadPickup_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM
    ! LOCAL
    CHARACTER(4)  :: rankChar
    INTEGER       :: iEl, istat
    INTEGER       :: thisRec, fUnit, i, j 
    INTEGER       :: iEq, N, iFace, e1, s1, bID
    LOGICAL       :: itExists
    CHARACTER(13) :: timeStampString
#ifdef HAVE_CUDA
    TYPE(dim3) :: tBlock, grid
#endif
   
    timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )

    N = myDGSEM % params % polyDeg

    WRITE(rankChar,'(I4.4)') myDGSEM % extComm % myRank
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
      DO iEl = 1, myDGSEM % mesh % elements % nElements

        DO iEq = 1, myDGSEM  % state % nEquations
          READ( fUnit, REC=thisRec )myDGSEM % state % solution(:,:,:,iEq,iEl)
          thisRec = thisRec+1
        ENDDO
        DO iEq = 1, myDGSEM % state % nEquations
          READ( fUnit, REC=thisRec )myDGSEM % static % solution(:,:,:,iEq,iEl)
          thisRec = thisRec+1
        ENDDO
        DO iEq = 1, myDGSEM % state % nEquations
          READ( fUnit, REC=thisRec )myDGSEM % static % source(:,:,:,iEq,iEl)
          thisRec = thisRec+1
        ENDDO

        READ( fUnit, REC=thisRec )myDGSEM % sourceTerms % drag(:,:,:,iEl)
        thisRec = thisRec+1

      ENDDO

      CLOSE(UNIT=fUnit)

    ENDIF

    INQUIRE( FILE='State.'//rankChar//'.'//timeStampString//'.exs', EXIST = itExists )

    IF( itExists )THEN

      OPEN( UNIT   = NEWUNIT(fUnit), &
        FILE   = 'State.'//rankChar//'.'//timeStampString//'.exs', &
        FORM   ='UNFORMATTED',&
        ACCESS ='DIRECT',&
        STATUS ='OLD',&
        ACTION ='READ', &
        CONVERT='BIG_ENDIAN',&
        RECL   = prec*(myDGSEM % params % polyDeg+1)*(myDGSEM % params % polyDeg+1)*(myDGSEM % state % nEquations)*(myDGSEM % extComm % nBoundaries) )
      READ( fUnit, rec = 1 ) myDGSEM % state % externalState
      READ( fUnit, rec = 2 ) myDGSEM % state % prescribedState
      CLOSE(fUnit)

    ENDIF


#ifdef HAVE_CUDA
    CALL myDGSEM % static % UpdateDevice( )
    CALL myDGSEM % state % UpdateDevice( )
    CALL myDGSEM % sourceTerms % UpdateDevice( )
    istat = cudaDeviceSynchronize( )

    tBlock = dim3(myDGSEM % params % polyDeg+1, &
                  myDGSEM % params % polyDeg+1 , &
                  1 )
    grid = dim3(myDGSEM % state % nEquations, myDGSEM % state % nElements, 1)  

    CALL CalculateStateAtBoundaries_CUDAKernel<<<grid, tBlock>>>( myDGSEM % static % solution_dev, &
                                                                         myDGSEM % static % boundarySolution_dev,  &
                                                                         myDGSEM % dgStorage % boundaryInterpolationMatrix_dev )

   CALL myDGSEM % static % UpdateHost( )
   istat = cudaDeviceSynchronize( )
#else

    ! Interpolate the static state to the element boundaries
    CALL myDGSEM % static % Calculate_Solution_At_Boundaries( myDGSEM % dgStorage )

#endif

  ! Update the static external states
  DO bID = 1, myDGSEM % extComm % nBoundaries

    iFace = myDGSEM % extComm % boundaryIDs( bID ) ! Obtain the process-local face id for this boundary-face id
    e1    = myDGSEM % mesh % faces % elementIDs(1,iFace)
    s1    = myDGSEM % mesh % faces % elementSides(1,iFace)

    DO j = 0, myDGSEM % params % polyDeg
      DO i = 0, myDGSEM % params % polyDeg
        DO iEq = 1, myDGSEM % static % nEquations
          myDGSEM % static % externalState(i,j,iEq,bID) = myDGSEM % static % boundarySolution(i,j,iEq,s1,e1)
        ENDDO
      ENDDO
    ENDDO

  ENDDO

#ifdef HAVE_CUDA
  myDGSEM % static % externalState_dev = myDGSEM % static % externalState
  istat = cudaDeviceSynchronize( )
#endif

    PRINT*, 'S/R ReadPickup : Done.'


  END SUBROUTINE ReadPickup_Fluid
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
  ATTRIBUTES(Global) SUBROUTINE CalculateSGSCoefficients_CUDAKernel( solution, static, smoothState, filterMat, sgsCoeffs  )
  
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)    :: solution(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: static(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(inout) :: smoothState(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: filterMat(0:polyDeg_dev,0:polyDeg_dev)
    REAL(prec), DEVICE, INTENT(inout) :: sgsCoeffs(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nSGS_dev,1:nEl_dev)
     ! Local
    INTEGER :: iEl, i, j, k, m, ii, jj, kk
    REAL(prec) :: sgsKE, uijk, uij, ui
    
    iEl = blockIDx % x
    
    i = threadIdx % x-1
    j = threadIdx % y-1
    k = threadIdx % z-1
    
     ! Here, the SGS Kinetic energy is calculated using the
     ! "high wavenumber" component of the velocity field.
     ! This component is defined (here) as the dIFference
     ! between the full solution and the smoothed solution.
    
    sgsKE = 0.0_prec
    DO m = 1, 3
      sgsKE = sgsKE + &
        ( solution(i,j,k,m,iEl)/( solution(i,j,k,4,iEl) + static(i,j,k,4,iEl))- &
        smoothState(i,j,k,m,iEl)/(smoothState(i,j,k,4,iEl)+static(i,j,k,4,iEl)) )**2
    ENDDO
    
     ! Now we calculate the viscosity and dIFfusivities (currently assumes isotropic and low mach number)
    DO m = 1, nSGS_dev
      sgsCoeffs(i,j,k,m,iEl) = 0.09_prec*viscLengthScale_dev*sqrt( sgsKE )
    ENDDO
  
  END SUBROUTINE CalculateSGSCoefficients_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE UpdateExternalSGSCoeffs_CUDAKernel( boundaryIDs, elementIDs, elementSides, procIDs, &
                                                                    externalsgsCoeffs, sgsCoeffsBsols, nHat )
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nBoundaryFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: procIDs(1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(inout) :: externalsgsCoeffs(0:polydeg_dev,0:polydeg_dev,1:nSGS_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)  :: sgsCoeffsBsols(0:polydeg_dev,0:polydeg_dev,1:nSGS_dev,1:6,1:nEl_dev)
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
        externalsgsCoeffs(i,j,iEq,iFace2) = sgsCoeffsBsols(i,j,iEq,s1,e1)
      ENDIF
    
    ENDIF
  
  END SUBROUTINE UpdateExternalSGSCoeffs_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE UpdateExternalState_CUDAKernel( boundaryIDs, elementIDs, elementSides, procIDs, &
                                                                externalState, stateBsols, prescribedState, nHat )
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)       :: boundaryIDs(1:nBoundaryFaces_dev)
    INTEGER, DEVICE, INTENT(in)       :: elementIDs(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)       :: elementSides(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)       :: procIDs(1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(inout) :: externalState(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)    :: stateBsols(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: prescribedState(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)    :: nhat(1:3,0:polyDeg_dev,0:polyDeg_dev,1:6,1:nEl_dev)
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
    
      iFace2 = boundaryIDs( iFace ) ! Obtain the process-local face id for this boundary-face id
      e1     = elementIDs(1,iFace2)
      s1     = elementSides(1,iFace2)
      e2     = elementIDs(2,iFace2)
      p2     = procIDs( iFace )
    
      IF( p2 == myRank_dev )THEN
    
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
    
          externalState(i,j,1,iFace) = -nx*un + us*sx + ut*tx ! u
          externalState(i,j,2,iFace) = -ny*un + us*sy + ut*ty ! v
          externalState(i,j,3,iFace) = -nz*un + us*sz + ut*tz ! w
          externalState(i,j,4,iFace) =  stateBsols(i,j,4,s1,e1) ! rho
          externalState(i,j,5,iFace) =  stateBsols(i,j,5,s1,e1) ! potential temperature
#ifdef PASSIVE_TRACERS
          externalState(i,j,6,iFace) =  stateBsols(i,j,6,s1,e1) ! tracer
#endif
          externalState(i,j,nEq_dev,iFace) =  stateBsols(i,j,nEq_dev,s1,e1) ! P
    
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
               prescribedState(i,j,1,iFace)

          us = stateBsols(i,j,1,s1,e1)*sx    + &
            stateBsols(i,j,2,s1,e1)*sy    + &
            stateBsols(i,j,3,s1,e1)*sz
          ut = stateBsols(i,j,1,s1,e1)*tx  + &
            stateBsols(i,j,2,s1,e1)*ty  + &
            stateBsols(i,j,3,s1,e1)*tz
    
          externalState(i,j,1,iFace) = -nx*un + us*sx + ut*tx ! u
          externalState(i,j,2,iFace) = -ny*un + us*sy + ut*ty ! v
          externalState(i,j,3,iFace) = -nz*un + us*sz + ut*tz ! w
          externalState(i,j,4,iFace) =  prescribedState(i,j,4,iFace) ! rho
          externalState(i,j,5,iFace) =  prescribedState(i,j,5,iFace) ! potential temperature
#ifdef PASSIVE_TRACERS
          externalState(i,j,6,iFace) =  prescribedState(i,j,6,iFace) ! tracer
#endif
          externalState(i,j,nEq_dev,iFace) =  prescribedState(i,j,nEq_dev,iFace) ! P
    
        ENDIF

      ENDIF
    
  END SUBROUTINE UpdateExternalState_CUDAKernel
!
ATTRIBUTES(Global) SUBROUTINE InternalFace_StateFlux_CUDAKernel( elementIDs, elementSides, boundaryIDs, iMap, jMap, &
                                                                 nHat, boundarySolution, boundarySolution_static, &
                                                                 boundaryFlux, stressFlux )

   IMPLICIT NONE
   INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
   INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nFaces_dev)
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
   INTEGER    :: ii, jj, bID
   INTEGER    :: e1, s1, e2, s2
   REAL(prec) :: uOut, uIn, cIn, cOut, norm, T
#ifdef PASSIVE_TRACERS
   REAL(prec) :: jump(1:6), aS(1:6)
#else
   REAL(prec) :: jump(1:5), aS(1:5)
#endif
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
                              ( (boundarySolution(ii,jj,nEq_dev,s2,e2)+boundarySolution_static(ii,jj,nEq_dev,s2,e2))/ P0_dev )**rC_dev   )
                        
                  T =   (boundarySolution_static(i,j,5,s1,e1) + boundarySolution(i,j,5,s1,e1))/&
                          (boundarySolution(i,j,4,s1,e1)+boundarySolution_static(i,j,4,s1,e1) )        
                             
                  cIn  = sqrt( R_dev*T* &
                              ( (boundarySolution(i,j,nEq_dev,s1,e1)+boundarySolution_static(i,j,nEq_dev,s1,e1))/P0_dev )**rC_dev  )
                               
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
                  aS(k) = aS(k) + (boundarySolution(i,j,nEq_dev,s1,e1) + &
                                   boundarySolution(ii,jj,nEq_dev,s2,e2))*nHat(k,i,j,s1,e1)/norm
                  ENDDO    
      
                         
                  DO iEq = 1, nEq_dev-1

                     boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*( aS(iEq) - fac*jump(iEq) )*norm
                     boundaryFlux(ii,jj,iEq,s2,e2) = -boundaryFlux(i,j,iEq,s1,e1)

                     IF( iEq == 4 )THEN
                        DO k = 1, 3
                           ! Calculate the LDG flux for the stress tensor.
                           stressFlux(k,i,j,iEq,s1,e1) = 0.5_prec*( boundarySolution(i,j,iEq,s1,e1) +&
                                                                    boundarySolution(ii,jj,iEq,s2,e2))*& 
                                                                    nHat(k,i,j,s1,e1)
                                                                                        
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
                                                                   boundarySolution_static(ii,jj,4,s2,e2)) )*& 
                                                                  nHat(k,i,j,s1,e1)
                                                                                        
                           stressFlux(k,ii,jj,iEq,s2,e2) = -stressFlux(k,i,j,iEq,s1,e1)

                        ENDDO
                     ENDIF
                     
                  ENDDO
             
               ENDIF 


 END SUBROUTINE InternalFace_StateFlux_CUDAKernel
!
 ATTRIBUTES(Global) SUBROUTINE BoundaryFace_StateFlux_CUDAKernel( elementIDs, elementSides, boundaryIDs, iMap, jMap, &
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
   REAL(prec), DEVICE, INTENT(inout) :: boundaryFlux(0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   REAL(prec), DEVICE, INTENT(inout) :: stressFlux(1:3,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:6,1:nEl_dev)
   ! Local
   INTEGER    :: iEl, iFace
   INTEGER    :: i, j, k, iEq
   INTEGER    :: ii, jj, bID
   INTEGER    :: e1, s1, e2, s2
   REAL(prec) :: uOut, uIn, cIn, cOut, norm, T
#ifdef PASSIVE_TRACERS
   REAL(prec) :: jump(1:6), aS(1:6)
#else
   REAL(prec) :: jump(1:5), aS(1:5)
#endif
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
                  jump(iEq)  = externalState(ii,jj,iEq,bID)-boundarySolution(i,j,iEq,s1,e1)
                  ENDDO
                 
                  T =   (boundarySolution_static(i,j,5,s1,e1) + externalState(ii,jj,5,bID))/&
                          (externalState(ii,jj,4,bID)+boundarySolution_static(i,j,4,s1,e1) )
                 ! Sound speed estimate for the external and internal states
                  cOut = sqrt( R_dev*T* &
                              ( (externalState(ii,jj,nEq_dev,bID)+boundarySolution_static(i,j,nEq_dev,s1,e1))/ P0_dev )**rC_dev   )
                  
                  T =   (boundarySolution_static(i,j,5,s1,e1) + boundarySolution(i,j,5,s1,e1))/&
                          (boundarySolution(i,j,4,s1,e1)+boundarySolution_static(i,j,4,s1,e1) )  
                                   
                  cIn  = sqrt( R_dev*T* &
                              ( (boundarySolution(i,j,nEq_dev,s1,e1)+boundarySolution_static(i,j,nEq_dev,s1,e1))/P0_dev )**rC_dev  )
                               
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
                  aS(k) = aS(k) + (boundarySolution(i,j,nEq_dev,s1,e1)+externalState(ii,jj,nEq_dev,bID))*nHat(k,i,j,s1,e1)/norm
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
                                                                   boundarySolution_static(i,j,4,s1,e1))  )*& 
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
  ATTRIBUTES(Global) SUBROUTINE CalculateSourceTerms_CUDAKernel( solution, static, staticSource, source, drag )
  
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)    :: solution(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
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

    ENDIF
  
  END SUBROUTINE CalculateSourceTerms_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE UpdateExternalStress_CUDAKernel( boundaryIDs, elementIDs, elementSides, procIDs, &
                                                                 externalStress, stressBsols, prescribedStress )
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nBoundaryFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: procIDs(1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(inout) :: externalStress(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)  :: stressBsols(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: prescribedStress(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:nBoundaryFaces_dev)
     ! Local
    INTEGER    :: iEq, bID, i, j, k
    INTEGER    :: iFace, p2
    INTEGER    :: e1, s1, s2, m
    
    bID = blockIdx % x
    iEq = blockIDx % y
     ! ////////////////////////////////////////////////////////////////////////// !
    i   = threadIdx % x-1
    j   = threadIdx % y-1
    
    iFace = boundaryIDs( bID ) ! Obtain the process-local face id for this boundary-face id
    e1    = elementIDs(1,iFace )
    s1    = elementSides(1,iFace)
    p2    = procIDs( bID )
    
      IF( p2 == myRank_dev )THEN
        externalStress(i,j,iEq,bID) = stressBsols(i,j,iEq,s1,e1)
      ENDIF
    
  
  END SUBROUTINE UpdateExternalStress_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE InternalFace_StressFlux_CUDAKernel( boundarySolution, boundaryStress, boundaryViscosity, &
                                                                    staticBoundarySolution, lengthScale, nHat, boundaryStressFlux, &
                                                                    elementIDs, elementSides, boundaryIDs, iMap, jMap )
  
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: boundarySolution(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryStress(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryViscosity(0:polyDeg_dev,0:polyDeg_dev,1:nSGS_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: staticBoundarySolution(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:polyDeg_dev,0:polyDeg_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: lengthScale(0:polyDeg_dev,0:polyDeg_dev,1:6,1:nEl_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: iMap(0:polyDeg_dev,0:polyDeg_dev,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: jMap(0:polyDeg_dev,0:polyDeg_dev,1:nFaces_dev)
    REAL(prec), DEVICE, INTENT(inout) :: boundaryStressFlux(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:6,1:nEl_dev)
     ! Local
    INTEGER    :: iEl, iFace
    INTEGER    :: i, j, iEq
    INTEGER    :: ii, jj
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
    
    ii = iMap(i,j,iFace)
    jj = jMap(i,j,iFace)
    
    norm = sqrt( nHat(1,i,j,s1,e1)**2 + nHat(2,i,j,s1,e1)**2 + nHat(3,i,j,s1,e1)**2 )
    IF( e2 > 0 )THEN

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
                                          ( boundaryViscosity(ii,jj,iEq,s2,e2)*boundarySolution(ii,jj,iEq,s2,e2)-&
                                            boundaryViscosity(i,j,iEq,s1,e1)*boundarySolution(i,j,iEq,s1,e1) )/&
                                          ( 0.5_prec*(lengthScale(i,j,s1,e1)+lengthScale(ii,jj,s2,e2)) )*norm


      boundaryStressFlux(ii,jj,iEq,s2,e2) = -boundaryStressFlux(i,j,iEq,s1,e1)

    ENDIF

  END SUBROUTINE InternalFace_StressFlux_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE BoundaryFace_StressFlux_CUDAKernel( boundarySolution, externalSolution, boundaryStress, externalStress, boundaryViscosity, externalViscosity, &
                                                                    staticBoundarySolution, externalStaticSolution, lengthScale, nHat, boundaryStressFlux, &
                                                                    elementIDs, elementSides, boundaryIDs, iMap, jMap )
  
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: boundarySolution(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: externalSolution(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryStress(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: externalStress(0:polyDeg_dev,0:polyDeg_dev,1:nStress_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryViscosity(0:polyDeg_dev,0:polyDeg_dev,1:nSGS_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: externalViscosity(0:polyDeg_dev,0:polyDeg_dev,1:nSGS_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)  :: staticBoundarySolution(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: externalStaticSolution(0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nBoundaryFaces_dev)
    REAL(prec), DEVICE, INTENT(in)  :: lengthScale(0:polyDeg_dev,0:polyDeg_dev,1:6,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:polyDeg_dev,0:polyDeg_dev,1:6,1:nEl_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces_dev)
    INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nFaces_dev)
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
    bID  = ABS(boundaryIDs(iFace))
    
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
                                          ( externalViscosity(ii,jj,iEq,bID)*externalSolution(ii,jj,iEq,bID)-&
                                            boundaryViscosity(i,j,iEq,s1,e1)*boundarySolution(i,j,iEq,s1,e1) )/(lengthScale(i,j,s1,e1))*norm



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
    REAL(prec) :: rhoT
    
    iEl = blockIdx % x
    i   = threadIdx % x - 1
    j   = threadIdx % y - 1
    k   = threadIdx % z - 1
    
     ! Pressure = rho*e*R/Cv (Ideal Gas Law, P = rho*R*T, thermo ~> T = theta*(P/P0)^r)
     ! THEN P = P0*(rho*theta*R/P0)^(Cp/Cv)
     ! And P' = P - P_static
    rhoT = static(i,j,k,5,iEl) + solution(i,j,k,5,iEl)
    solution(i,j,k,nEq_dev,iEl) = P0_dev*( rhoT*R_dev/P0_dev )**hCapRatio_dev - static(i,j,k,nEq_dev,iEl)
  
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
  ATTRIBUTES(Global) SUBROUTINE CalculateSGSAtBoundaries_3D_CUDAKernel( f, fAtBoundaries, boundaryMatrix ) 
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: f(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nSGS_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryMatrix(0:polydeg_dev,0:1)
    REAL(prec), DEVICE, INTENT(out) :: fAtBoundaries(0:polydeg_dev,0:polydeg_dev,1:nSGS_dev,1:6,1:nEl_dev)
    ! Local
    INTEGER    :: iEq, iEl, i, j, k
    REAL(prec) :: bSol(1:6)

      iEq = blockIdx % x
      iEl = blockIdx % y

      k   = threadIdx % y-1
      j   = threadIdx % x-1
      
      IF( j <= polyDeg_dev .AND. k <= polyDeg_dev )THEN
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
      
      ENDIF
      
  END SUBROUTINE CalculateSGSAtBoundaries_3D_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE CalculateNormalStressAtBoundaries_CUDAKernel( solutionGradient, boundaryViscosity, nHat, boundaryStress, boundaryMatrix ) 
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: solutionGradient(1:3,0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryViscosity(0:polydeg_dev,0:polydeg_dev,1:nSGS_dev,1:6,1:nEl_dev)
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

        boundaryStress(i,j,iEq,k,iEl) = boundaryViscosity(i,j,iEq,k,iEl)*( fAtBoundaries(1,k)*nHat(1,i,j,k,iEl) + &
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

  ATTRIBUTES(Global) SUBROUTINE CalculateStressFlux_CUDAKernel( solutionGradient, viscosity, state, static, Ja, stressFlux )

    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: solutionGradient(1:3,0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: viscosity(0:polyDeg_dev,0:polyDeg_dev,0:polyDeg_dev,1:nSGS_dev,1:nEl_dev)
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
  
        sf(1) = sf(1) + Ja(i,j,k,idir,1,iEl)*solutionGradient(idir,i,j,k,iEq,iEl)*viscosity(i,j,k,iEq,iEl)
  
        sf(2) = sf(2) + Ja(i,j,k,idir,2,iEl)*solutionGradient(idir,i,j,k,iEq,iEl)*viscosity(i,j,k,iEq,iEl) 
  
        sf(3) = sf(3) + Ja(i,j,k,idir,3,iEl)*solutionGradient(idir,i,j,k,iEq,iEl)*viscosity(i,j,k,iEq,iEl) 
  
      ENDDO

    ELSE

      rho = state(i,j,k,4,iEl) + static(i,j,k,4,iEl)
      DO idir = 1, 3
  
        sf(1) = sf(1) + Ja(i,j,k,idir,1,iEl)*rho*solutionGradient(idir,i,j,k,iEq,iEl)*viscosity(i,j,k,iEq,iEl)
  
        sf(2) = sf(2) + Ja(i,j,k,idir,2,iEl)*rho*solutionGradient(idir,i,j,k,iEq,iEl)*viscosity(i,j,k,iEq,iEl) 
  
        sf(3) = sf(3) + Ja(i,j,k,idir,3,iEl)*rho*solutionGradient(idir,i,j,k,iEq,iEl)*viscosity(i,j,k,iEq,iEl) 
  
      ENDDO

    ENDIF

    stressflux(1,i,j,k,iEq,iEl) = sf(1)
    stressflux(2,i,j,k,iEq,iEl) = sf(2)
    stressflux(3,i,j,k,iEq,iEl) = sf(3)

  END SUBROUTINE CalculateStressFlux_CUDAKernel

  ATTRIBUTES( Global ) SUBROUTINE FilterState_CUDAKernel( f, filterMatrix )
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(inout) :: f(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nEq_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)    :: filterMatrix(0:polydeg_dev,0:polydeg_dev)
    ! Local
    INTEGER            :: i, j, k, iVar, iEl
    INTEGER            :: ii, jj, kk
    REAL(prec)         :: uijk, uij, ui
    REAL(prec), SHARED :: fLocal(0:7,0:7,0:7)


      iVar = blockIdx % x
      iEl  = blockIdx % y

      i = threadIdx % x-1
      j = threadIdx % y-1
      k = threadIdx % z-1

      IF( i <= polydeg_dev .AND. j <= polydeg_dev .AND. k <= polydeg_dev )THEN

        fLocal(i,j,k) = f(i,j,k,iVar,iEl)

        CALL syncthreads( )

        uijk = 0.0_prec

        DO kk = 0, polydeg_dev

          uij = 0.0_prec

          DO jj = 0, polydeg_dev

            ui = 0.0_prec

            DO ii = 0, polydeg_dev

              ui = ui + filterMatrix(ii,i)*fLocal(ii,jj,kk)

            ENDDO

            uij = uij + filterMatrix(jj,j)*ui

          ENDDO

          uijk = uijk + filterMatrix(kk,k)*uij

        ENDDO

        f(i,j,k,iVar,iEl) = uijk

      ENDIF


 END SUBROUTINE FilterState_CUDAKernel

#endif


END MODULE Fluid_Class


