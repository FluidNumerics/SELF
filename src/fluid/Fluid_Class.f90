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

    PROCEDURE :: GlobalTimeDerivative            => GlobalTimeDerivative_Fluid
    PROCEDURE :: EquationOfState                 => EquationOfState_Fluid
    PROCEDURE :: UpdateExternalState             => UpdateExternalState_Fluid
    PROCEDURE :: InternalFaceFlux                => InternalFaceFlux_Fluid
    PROCEDURE :: BoundaryFaceFlux                => BoundaryFaceFlux_Fluid
    PROCEDURE :: MappedTimeDerivative            => MappedTimeDerivative_Fluid

    PROCEDURE :: CalculateSGSCoefficients => CalculateSGSCoefficients_Fluid
    PROCEDURE :: UpdateExternalSGS        => UpdateExternalSGS_Fluid

    PROCEDURE :: CalculateStressTensor   => CalculateStressTensor_Fluid
    PROCEDURE :: UpdateExternalStress    => UpdateExternalStress_Fluid
    PROCEDURE :: InternalStressFlux      => InternalStressFlux_Fluid
    PROCEDURE :: BoundaryStressFlux      => BoundaryStressFlux_Fluid

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
  INTEGER, PARAMETER, PRIVATE :: nEquations   = 6

#ifdef HAVE_CUDA
 INTEGER, CONSTANT    :: nEq_dev
 INTEGER, CONSTANT    :: nStress_dev
 INTEGER, CONSTANT    :: nSGS_dev
 INTEGER, CONSTANT    :: polyDeg_dev
 INTEGER, CONSTANT    :: nEl_dev
 INTEGER, CONSTANT    :: myRank_dev
 INTEGER, CONSTANT    :: nFaces_dev
 INTEGER, CONSTANT    :: nBoundaryFaces_dev
 REAL(prec), CONSTANT :: R_dev, Cv_dev, P0_dev, hCapRatio_dev, rC_dev, g_dev
! REAL(prec), CONSTANT :: rk3_a_dev, rk3_g_dev, dt_dev
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
    !
#ifdef HAVE_CUDA
    INTEGER(KIND=cuda_count_KIND) :: freebytes, totalbytes
    INTEGER                       :: iStat, cudaDeviceNumber, nDevices

    CALL UpdateDeviceDictionary( )
#endif

    CALL myDGSEM % params % Build( setupSuccess )
    myDGSEM % simulationTime = myDGSEM % params % startTime

    IF( .NOT. SetupSuccess ) THEN
      PRINT(MsgFMT), 'S/R Build_Fluid : Halting before building,'
      RETURN
    ENDIF

    CALL myDGSEM % extComm % ReadPickup(  )

#ifdef HAVE_MPI
    CALL myDGSEM % extComm % ConstructCommTables(  )
#endif

#ifdef HAVE_CUDA

    CALL myDGSEM % extComm % UpdateDevice( )


    ! Assuming the number of GPU's and the number of ranks per node is unIForm,
    ! each rank is assigned to it's own GPU.
    iStat = cudaGetDeviceCount( nDevices )
    cudaDeviceNumber = MOD( myDGSEM % extComm % myRank, nDevices )

    PRINT*, '    S/R Build_Fluid : Rank :', &
      myDGSEM % extComm % myRank, ': Getting Device # ', cudaDeviceNumber

    iStat = cudaSetDevice( cudaDeviceNumber )

#endif

    ! Construct the DATA structure that holds the derivative and interpolation matrices
    ! and the quadrature weights. This call will also perform the device copies.
    CALL myDGSEM % dGStorage % Build( UniformPoints(-1.0_prec,1.0_prec,myDGSEM % params % nPlot), &
      myDGSEM % params % polyDeg, myDGSEM % params % nPlot, GAUSS )

    CALL myDGSEM % filter % Build( myDGSEM % dgStorage % interp % interpolationPoints,&
      myDGSEM % dgStorage % quadratureWeights, &
      myDGSEM % params % polyDeg, myDGSEM % params % nCutoff, &
      TanhRollOff )

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
                                         (myDGSEM % state % nEquations-1)*3, &
                                         myDGSEM % mesh % elements % nElements, &
                                         myDGSEM % extComm % nBoundaries )

    myDGSEM % sgsCoeffs % solution         = myDGSEM % params % viscosity
    myDGSEM % sgsCoeffs % boundarySolution = myDGSEM % params % viscosity
    myDGSEM % sgsCoeffs % externalState = myDGSEM % params % viscosity

#ifdef HAVE_MPI
    CALL myDGSEM % mpiStateHandler % Build( myDGSEM % extComm, myDGSEM % params % polyDeg, myDGSEM % state % nEquations )
    CALL myDGSEM % mpiStressHandler % Build( myDGSEM % extComm, myDGSEM % params % polyDeg, myDGSEM % stressTensor % nEquations )
    CALL myDGSEM % mpiSGSHandler % Build( myDGSEM % extComm, myDGSEM % params % polyDeg, myDGSEM % sgsCoeffs % nEquations )
#endif

#ifdef HAVE_CUDA
    CALL myDGSEM % sgsCoeffs % UpdateDevice( )
    CALL myDGSEM % sgsCoeffs % UpdateDevice( )
#endif

    ! Read the initial conditions, static state, and the boundary communicator
    CALL myDGSEM % ReadPickup( )

#ifdef HAVE_CUDA

      ! Set Device Constants
      R_dev         = myDGSEM % params % R
      Cv_dev        = myDGSEM % params % Cv
      P0_dev        = myDGSEM % params % P0
      hCapRatio_dev = ( myDGSEM % params % R + myDGSEM % params % Cv ) / myDGSEM % params % Cv
      rC_dev        = myDGSEM % params % R / ( myDGSEM % params % R + myDGSEM % params % Cv )
      g_dev         = myDGSEM % params % g
      !dt_dev        = myDGSEM % params % dt
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

#endif

  END SUBROUTINE Build_Fluid
!
  SUBROUTINE Trash_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM

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
    REAL(prec), DEVICE :: t_dev, dt_dev
    REAL(prec), DEVICE :: G3D(0:myDGSEM % params % polyDeg,&
                             0:myDGSEM % params % polyDeg,&
                             0:myDGSEM % params % polyDeg,&
                             1:myDGSEM % state % nEquations,&
                             1:myDGSEM % mesh % elements % nElements)
    INTEGER            :: iT, m, iStat
    TYPE(dim3)         :: grid, tBlock


    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) )
    grid = dim3(myDGSEM % state % nElements, myDGSEM % state % nEquations-1, 1)

    t0 = myDGSEM % simulationTime
    dt = myDGSEM % params % dt
    dt_dev = dt

    DO iT = 1, nT

      G3D  = 0.0_prec

      DO m = 1,3
        t = myDGSEM % simulationTime + rk3_b(m)*dt
        CALL myDGSEM % GlobalTimeDerivative( t )

        CALL UpdateG3D_CUDAKernel<<<grid,tBlock>>>( G3D, rk3_a_dev(m), rk3_g_dev(m), dt_dev, &
                                                    myDGSEM % state % solution_dev, &
                                                    myDGSEM % state % fluxDivergence_dev, &
                                                    myDGSEM % state % source_dev, &
                                                    myDGSEM % stressTensor % fluxDivergence_dev, &
                                                    myDGSEM % mesh % elements % J_dev, &
                                                    myDGSEM % params % polyDeg_dev, &
                                                    myDGSEM % state % nEquations_dev, &
                                                    myDGSEM % stressTensor % nEquations_dev, &
                                                    myDGSEM % mesh % elements % nElements_dev )

        CALL myDGSEM % EquationOfState( )

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
        CALL myDGSEM % GlobalTimeDerivative( t )

        CALL UpdateG3D_CUDAKernel<<<grid,tBlock>>>( G3D, rk3_a_dev(m), rk3_g_dev(m), dt_dev, &
                                                    myDGSEM % state % solution_dev, &
                                                    myDGSEM % state % fluxDivergence_dev, &
                                                    myDGSEM % state % source_dev, &
                                                    myDGSEM % stressTensor % fluxDivergence_dev, &
                                                    myDGSEM % mesh % elements % J_dev, &
                                                    myDGSEM % params % polyDeg_dev, &
                                                    myDGSEM % state % nEquations_dev, &
                                                    myDGSEM % stressTensor % nEquations_dev, &
                                                    myDGSEM % mesh % elements % nElements_dev )

        CALL myDGSEM % EquationOfState( )

      ENDDO 

      myDGSEM % simulationTime = myDGSEM % simulationTime + dt

    ENDIF

#else
    REAL(prec) :: t, dt, rk3_a_local, rk3_g_local
    REAL(prec) :: G3D(0:myDGSEM % params % polyDeg,&
                      0:myDGSEM % params % polyDeg,&
                      0:myDGSEM % params % polyDeg,&
                      1:myDGSEM % state % nEquations,&
                      1:myDGSEM % mesh % elements % nElements)
    INTEGER    :: m, iEl, iT, i, j, k, iEq

    t0 = myDGSEM % simulationTime
    dt = myDGSEM % params % dt


    DO iT = 1, nT

      !$OMP DO
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
      !$OMP ENDDO

      DO m = 1,3 ! Loop over RK3 steps

        t = myDGSEM % simulationTime + rk3_b(m)*dt
        CALL myDGSEM % GlobalTimeDerivative( t )

        !$OMP DO
        DO iEl = 1, myDGSEM % mesh % elements % nElements
          DO iEq = 1, myDGSEM % state % nEquations-1
            DO k = 0, myDGSEM % params % polyDeg
              DO j = 0, myDGSEM % params % polyDeg
                DO i = 0, myDGSEM % params % polyDeg

                  G3D(i,j,k,iEq,iEl) = rk3_a(m)*G3D(i,j,k,iEq,iEl) + ( myDGSEM % state % fluxDivergence(i,j,k,iEq,iEl) +&
                    myDGSEM % stressTensor % fluxDivergence(i,j,k,iEq,iEl) )/myDGSEM % mesh % elements % J(i,j,k,iEl) + &
                    myDGSEM % state % source(i,j,k,iEq,iEl)

                  myDGSEM % state % solution(i,j,k,iEq,iEl) = myDGSEM % state % solution(i,j,k,iEq,iEl) + &
                    rk3_g(m)*dt*G3D(i,j,k,iEq,iEl)

                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO 
        !$OMP ENDDO

        CALL myDGSEM % EquationOfState( )

      ENDDO

      myDGSEM % simulationTime = myDGSEM % simulationTime + dt

    ENDDO

    ! Determine IF we need to take another step with reduced time step to get the solution
    ! at exactly t0+outputFrequency
    IF( .NOT. AlmostEqual( myDGSEM % simulationTime, t0+myDGSEM % params % outputFrequency ) )THEN

      dt = t0+myDGSEM % params % outputFrequency - myDGSEM % simulationTime

      !$OMP DO
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
      !$OMP ENDDO

      DO m = 1,3

        t = myDGSEM % simulationTime + rk3_b(m)*dt
        CALL myDGSEM % GlobalTimeDerivative( t )


        !$OMP DO
        DO iEl = 1, myDGSEM % mesh % elements % nElements
          DO iEq = 1, myDGSEM % state % nEquations-1
            DO k = 0, myDGSEM % params % polyDeg
              DO j = 0, myDGSEM % params % polyDeg
                DO i = 0, myDGSEM % params % polyDeg

                  G3D(i,j,k,iEq,iEl) = rk3_a(m)*G3D(i,j,k,iEq,iEl) + ( myDGSEM % state % fluxDivergence(i,j,k,iEq,iEl) +&
                    myDGSEM % stressTensor % fluxDivergence(i,j,k,iEq,iEl) )/myDGSEM % mesh % elements % J(i,j,k,iEl) + &
                    myDGSEM % state % source(i,j,k,iEq,iEl)

                  myDGSEM % state % solution(i,j,k,iEq,iEl) = myDGSEM % state % solution(i,j,k,iEq,iEl) + &
                    rk3_g(m)*dt*G3D(i,j,k,iEq,iEl)

                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        !$OMP ENDDO

        CALL myDGSEM % EquationOfState( )

      ENDDO

      myDGSEM % simulationTime = myDGSEM % simulationTime +  dt

    ENDIF
#endif



  END SUBROUTINE ForwardStepRK3_Fluid
!
! ============================================================================= !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ============================================================================= !
!
! Crank Nicholson time integrator routines

! SUBROUTINE CrankNicholsonBiCGStab_Fluid( myDGSEM, snk, explicitTendency )
!   IMPLICIT NONE
!   CLASS(Fluid), INTENT(inout) :: myDGSEM
!   REAL(prec), INTENT(in)      :: snk(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!   REAL(prec), INTENT(in)      :: explicitTendency(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!   ! Local
!   REAL(prec)   :: r(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!   REAL(prec)   :: ds(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!   REAL(prec)   :: v(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!   REAL(prec)   :: p(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!   REAL(prec)   :: t(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!   REAL(prec)   :: rho, alpha, omega, beta
!
!      ! Calculate the initial residual
!      ! Assumes an initial guess of ds=0
!      r = explicitTendency + myDGSEM % CrankNicholsonRHS( snk )
!
!
!
! END SUBROUTINE CrankNicholsonBiCGStab_Fluid
!!
! FUNCTION CrankNicholsonRHS_Fluid( myDGSEM, snk ) RESULT( b )
!   ! Given
!   IMPLICIT NONE
!   CLASS(Fluid) :: myDGSEM
!   REAL(prec)   :: snk(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!   REAL(prec)   :: b(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!
!
!      CALL myDGSEM % GlobalTimeDerivative( myDGSEM % simulationTime )
!      b = -( snk - 0.5_prec*myDGSEM % params % dt*myDGSEM % state % tendency )
!
!
! END FUNCTION CrankNicholsonRHS_Fluid
!!
! FUNCTION CrankNicholsonJacobianAction_Fluid( myDGSEM, s, ds, Fs ) RESULT( Jds )
!   IMPLICIT NONE
!   CLASS(Fluid) :: myDGSEM
!   REAL(prec)   :: s(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!   REAL(prec)   :: ds(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!   REAL(prec)   :: Fs(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!   REAL(prec)   :: Jds(0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 0:myDGSEM % params % polyDeg, 1:nEq, 1:myDGSEM % mesh % elements % nElements)
!
!
!      myDGSEM % state % solution = s + myDGSEM % params % jacobianStepSize*ds
!
!      CALL myDGSEM % GlobalTimeDerivative( myDGSEM % simulationTime )
!
!      ! J*ds = (I - (dt/2)* dF/ds )*ds
!      Jds = ds - 0.5_prec*myDGSEM % params % dt*( myDGSEM % state % tendency - Fs )/myDGSEM % params % jacobianStepSize
!
! END FUNCTION CrankNicholsonJacobianAction_Fluid
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
!  IF SpectralFiltering is USEd as the subgridscale model, THEN the spectral
!  filtering matrix (specIFied in src/filtering/SpectralFilter_Class.f90) is USEd
!  to smooth the solution variables before proceeding.

    IF( myDGSEM % params % SubGridModel == SpectralFiltering )THEN

#ifdef HAVE_CUDA
      CALL myDGSEM % filter % Filter3D( myDGSEM % state % solution_dev, &
                                        myDGSEM % state % solution_dev, &
                                        myDGSEM % state % nEquations_dev, &
                                        myDGSEM % state % nElements_dev )

#else
      CALL myDGSEM % filter % Filter3D( myDGSEM % state % solution, &
                                        myDGSEM % state % solution, &
                                        myDGSEM % state % nEquations, &
                                        myDGSEM % state % nElements )
#endif

    ENDIF

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
!  Here, the solution within each element is interpolated to the faces of each
!  element in order to prepare for computing the external state for enforcing
!  boundary conditions, Riemann Fluxes, and MPI DATA exchanges that need to
!  occur.
#ifdef HAVE_CUDA

    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                  1 )
    grid = dim3(myDGSEM % state % nEquations, myDGSEM % state % nElements, 1)  

    CALL CalculateStateAtBoundaries_3D_CUDAKernel<<<grid, tBlock>>>( myDGSEM % state % solution_dev, &
                                                                         myDGSEM % state % boundarySolution_dev,  &
                                                                         myDGSEM % dgStorage % boundaryInterpolationMatrix_dev )
#else

    CALL myDGSEM % state % Calculate_Solution_At_Boundaries( myDGSEM % dgStorage )

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

#ifdef HAVE_MPI
    CALL myDGSEM % mpiStateHandler % MPI_Exchange( myDGSEM % state, &
                                                   myDGSEM % mesh % faces, &
                                                   myDGSEM % extComm )
#endif
! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! The boundary solutions are USEd to calculate the external states that, when
! accompanied with a Riemann Solver, enforce boundary conditions. Calling this
! routine is dependent on the result of CalculateBoundarySolution

    CALL myDGSEM % UpdateExternalState( tn )

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
    CALL myDGSEM % mpiStateHandler % Finalize_MPI_Exchange( myDGSEM % state, &
                                                            myDGSEM % mesh % faces, &
                                                            myDGSEM % extComm )
#endif

    CALL myDGSEM % BoundaryFaceFlux( )

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

        CALL myDGSEM % CalculateSGSCoefficients( )

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
! The viscosity coefficient that is calculated is now interpolated to the faces
! of each element so that the viscous flux can later be computed. This routine
! depends on the result of CalculateSGSCoefficients.

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

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! The viscosity coefficients are exchanged with neighboring ranks that share
! COMMON faces. MPI_SGSExchange can be run simulataneously with
! CalculateStressTensor, CalculateBoundaryStress, UpdateExternalStress, and the
! MPI_StressExchange. The viscosity coefficients that are exchanged are not
! needed until StressFlux

#ifdef HAVE_MPI
        CALL myDGSEM % mpiSGSHandler % MPI_Exchange( myDGSEM % sgsCoeffs, &
                                                     myDGSEM % mesh % faces, &
                                                     myDGSEM % extComm )
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

#ifdef HAVE_MPI
      IF( myDGSEM % params % SubGridModel == SpectralEKE )THEN !
        CALL myDGSEM % mpiSGSHandler % Finalize_MPI_Exchange( myDGSEM % sgsCoeffs, &
                                                              myDGSEM % mesh % faces, &
                                                              myDGSEM % extComm )
      ENDIF
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! The stress tensor values are interpolated to the faces of each element to
! prepare for the calculation of the divergence of the viscous fluxes. This
! routine depends on the result of CalculateStressTensor.

#ifdef HAVE_CUDA

      tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                    1 )
      grid = dim3(myDGSEM % stressTensor % nEquations, myDGSEM % stressTensor % nElements, 1)  
      CALL CalculateStressAtBoundaries_3D_CUDAKernel<<<grid, tBlock>>>( myDGSEM % stressTensor % solution_dev, &
                                                                           myDGSEM % stressTensor % boundarySolution_dev,  &
                                                                           myDGSEM % dgStorage % boundaryInterpolationMatrix_dev )

#else
  
      CALL myDGSEM % stressTensor % Calculate_Solution_At_Boundaries( myDGSEM % dgStorage )

#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! Stress tensor values are exchanged with neighboring ranks along shared faces.
! This routine depends on the result of CalculateBoundaryStress, but can be run
! at the same time as UpdateExternalStress.

#ifdef HAVE_MPI
      CALL myDGSEM % mpiStressHandler % MPI_Exchange( myDGSEM % stressTensor, &
                                                      myDGSEM % mesh % faces, &
                                                      myDGSEM % extComm )
#endif

! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! Now that the stress tensor is available on element faces, boundary conditions
! can be applied by setting the external stress tensor state. This routine
! depends on the result of CalculateBoundaryStress. Note that this routine can
! be run simultaneously with the MPI_StressExchange

      CALL myDGSEM % UpdateExternalStress( tn )

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
      CALL myDGSEM % mpiStressHandler % Finalize_MPI_Exchange( myDGSEM % stressTensor,&
                                                               myDGSEM % mesh % faces, &
                                                               myDGSEM % extComm )
#endif
      CALL myDGSEM % BoundaryStressFlux( )

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
#ifdef HAVE_CUDA
      tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                    4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                    4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) )
      grid = dim3(myDGSEM % state % nEquations, myDGSEM % stressTensor % nElements, 1)  

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


    ENDIF


! ----------------------------------------------------------------------------- !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< !
! ----------------------------------------------------------------------------- !
!
! Once the inviscid fluxes through the faces are calculated, and the internal
! state is known, the tendency due to the inviscid flux terms and
! nonconservative source terms is calculated here. This routine depends on the
! result of FaceFlux, but can be DOne at the same time as StressDivergence

    CALL myDGSEM % MappedTimeDerivative( )

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

    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
      4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
      4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) )
    grid = dim3(myDGSEM % mesh % elements % nElements, 1, 1)

    CALL CalculateSGSCoefficients_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                               myDGSEM % static % solution_dev, &
                                                               myDGSEM % smoothState % solution_dev, &
                                                               myDGSEM % filter % filterMat_dev, &
                                                               myDGSEM % sgsCoeffs % solution_dev, &
                                                               myDGSEM % params % viscLengthScale_dev, &
                                                               myDGSEM % params % polyDeg_dev, &
                                                               myDGSEM % state % nEquations_dev, &
                                                               myDGSEM % mesh % elements % nElements_dev )
#else
    ! Local
    INTEGER :: iEl, i, j, k, m, ii, jj, kk
    REAL(prec) :: sgsKE, uijk, uij, ui
    REAL(prec) :: KE(0:myDGSEM % params % polyDeg,0:myDGSEM % params % polyDeg,0:myDGSEM % params % polyDeg)


    !$OMP DO PRIVATE( KE )
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
    !$OMP ENDDO

#endif

  END SUBROUTINE CalculateSGSCoefficients_Fluid
!
  SUBROUTINE UpdateExternalSGS_Fluid( myDGSEM ) ! ////////// !

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                  1 )
    grid = dim3(myDGSEM % extComm % nBoundaries,myDGSEM % state % nEquations-1,1)

    CALL UpdateExternalSGSCoeffs_CUDAKernel<<<grid, tBlock>>>( myDGSEM % extComm % boundaryIDs_dev, &       ! I
                                                               myDGSEM % mesh % faces % elementIDs_dev, &   ! I
                                                               myDGSEM % mesh % faces % elementSides_dev, & ! I
                                                               myDGSEM % extComm % extProcIDs_dev, &           ! I
                                                               myDGSEM % sgsCoeffs % externalState_dev, &                    ! O
                                                               myDGSEM % sgsCoeffs % boundarySolution_dev, &   ! I
                                                               myDGSEM % mesh % elements % nHat_dev, &
                                                               myDGSEM % extComm % myRank_dev, &
                                                               myDGSEM % params % polyDeg_dev, &
                                                               myDGSEM % sgsCoeffs % nEquations_dev, &   ! I
                                                               myDGSEM % extComm % nBoundaries_dev, &
                                                               myDGSEM % mesh % faces % nFaces_dev, &
                                                               myDGSEM % mesh % elements % nElements_dev )           ! I
#else
    ! Local
    INTEGER    :: iEl, bID, bFaceID, i, j, k, iEq
    INTEGER    :: iFace2, p2
    INTEGER    :: e1, e2, s1, s2

    !$OMP DO
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
    !$OMP ENDDO

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

    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                  1 )
    grid = dim3(myDGSEM % extComm % nBoundaries,1,1)

    CALL UpdateExternalState_CUDAKernel<<<grid, tBlock>>>( myDGSEM % extComm % boundaryIDs_dev, &       ! I
                                                           myDGSEM % mesh % faces % elementIDs_dev, &   ! I
                                                           myDGSEM % mesh % faces % elementSides_dev, & ! I
                                                           myDGSEM % extComm % extProcIDs_dev, &           ! I
                                                           myDGSEM % state % externalState_dev, &               ! O
                                                           myDGSEM % state % boundarySolution_dev, &    ! I
                                                           myDGSEM % state % prescribedState_dev, &             ! I
                                                           myDGSEM % mesh % elements % nHat_dev, &
                                                           myDGSEM % extComm % myRank_dev, &
                                                           myDGSEM % params % Cd_dev, &
                                                           myDGSEM % params % dragScale_dev, &
                                                           myDGSEM % params % polydeg_dev, &
                                                           myDGSEM % state % nEquations_dev, &    ! I
                                                           myDGSEM % extComm % nBoundaries_dev, &
                                                           myDGSEM % mesh % faces % nFaces_dev, &
                                                           myDGSEM % mesh % elements % nElements_dev )           ! I
#else
    ! Local
    INTEGER    :: iEl, bID, bFaceID, i, j, k, iEq
    INTEGER    :: iFace2, p2
    INTEGER    :: e1, e2, s1, s2
    REAL(prec) :: norm, un, ut, us, speed
    REAL(prec) :: nx, ny, nz
    REAL(prec) :: sx, sy, sz
    REAL(prec) :: tx, ty, tz

    !$OMP DO
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

            ! momentum
            ! rho*u
            myDGSEM % state % externalState(i,j,1,bID) = 0.0_prec
            ! rho*v
            myDGSEM % state % externalState(i,j,2,bID) = 0.0_prec
            ! rho*w
            myDGSEM % state % externalState(i,j,3,bID) = 0.0_prec
            ! Density is set to the static density field
            myDGSEM % state % externalState(i,j,4,bID) = 0.0_prec
            ! Potential Temperature anomaly (multiplied by density) is set to its static state
            myDGSEM % state % externalState(i,j,5,bID) = 0.0_prec
            ! Pressure anomaly is set to zero
            myDGSEM % state % externalState(i,j,6,bID) = 0.0_prec

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
            myDGSEM % state % externalState(i,j,6,bID) =  myDGSEM % state % boundarySolution(i,j,6,s1,e1) ! P

          ELSEIF( e2 == DRAG_SLIP.AND. p2 == myDGSEM % extComm % myRank )THEN

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

            myDGSEM % state % externalState(i,j,1,bID) = -nx*un + us*sx + ut*tx ! u
            myDGSEM % state % externalState(i,j,2,bID) = -ny*un + us*sy + ut*ty ! v
            myDGSEM % state % externalState(i,j,3,bID) = -nz*un + us*sz + ut*tz ! w
            myDGSEM % state % externalState(i,j,4,bID) =  myDGSEM % state % boundarySolution(i,j,4,s1,e1) ! rho
            myDGSEM % state % externalState(i,j,5,bID) =  myDGSEM % state % boundarySolution(i,j,5,s1,e1) ! potential temperature
            myDGSEM % state % externalState(i,j,6,bID) =  myDGSEM % state % boundarySolution(i,j,6,s1,e1) ! P


          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !$OMP ENDDO

#endif

  END SUBROUTINE UpdateExternalState_Fluid
!
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!                                                                                                 !
!  This section of code CONTAINS routines for computing the fluxes through the element faces      !
!                                                                                                 !
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!
  SUBROUTINE InternalFaceFlux_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg + 1)/4 ) ), &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg + 1)/4 ) ) , &
                  1 )
    grid = dim3(myDGSEM % mesh % faces % nFaces,1,1)

    CALL InternalFaceFlux_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mesh % faces % elementIDs_dev, &
                                                        myDGSEM % mesh % faces % elementSides_dev, &
                                                        myDGSEM % mesh % faces % boundaryID_dev, &
                                                        myDGSEM % mesh % faces % iMap_dev, &
                                                        myDGSEM % mesh % faces % jMap_dev, &
                                                        myDGSEM % mesh % elements % nHat_dev, &
                                                        myDGSEM % state % boundarySolution_dev, &
                                                        myDGSEM % static % boundarySolution_dev, &
                                                        myDGSEM % state % externalState_dev, &
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
    REAL(prec) :: jump(1:myDGSEM % state % nEquations-1), aS(1:myDGSEM % state % nEquations-1)
    REAL(prec) :: fac, rC

    !$OMP DO PRIVATE( jump, aS )
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
              ( (myDGSEM % state % boundarySolution(ii,jj,6,s2,e2)+&
              myDGSEM % static % boundarySolution(ii,jj,6,s2,e2))/&
              myDGSEM % params % P0 )**myDGSEM % params % rC   )

            T =   (myDGSEM % static % boundarySolution(i,j,5,s1,e1) + &
              myDGSEM % state % boundarySolution(i,j,5,s1,e1))/&
              (myDGSEM % static % boundarySolution(i,j,4,s1,e1) + &
              myDGSEM % state % boundarySolution(i,j,4,s1,e1) )

            cIn  = sqrt( myDGSEM % params % R*T* &
              ( (myDGSEM % state % boundarySolution(i,j,6,s1,e1)+&
              myDGSEM % static % boundarySolution(i,j,6,s1,e1))/&
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
              aS(k) = aS(k) + (myDGSEM % state % boundarySolution(i,j,6,s1,e1) + &
                myDGSEM % state % boundarySolution(ii,jj,6,s2,e2))*nHat(k)
            ENDDO


            DO iEq = 1, myDGSEM % state % nEquations-1
              myDGSEM % state % boundaryFlux(i,j,iEq,s1,e1) =  0.5_prec*( aS(iEq) - fac*jump(iEq) )*norm
              myDGSEM % state % boundaryFlux(ii,jj,iEq,s2,e2) = -myDGSEM % state % boundaryFlux(i,j,iEq,s1,e1)
!              IF( iEq == 4 )THEN
!
!                DO k = 1, 3
!                  ! Calculate the LDG flux for the stress tensor.
!                  myDGSEM % state % boundaryGradientFlux(k,i,j,iEq,s1,e1) = 0.5_prec*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) +&
!                    myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2) )*&
!                    myDGSEM % mesh % elements % nHat(k,i,j,s1,e1)
!
!                  myDGSEM % state % boundaryGradientFlux(k,ii,jj,iEq,s2,e2) = -myDGSEM % state % boundaryGradientFlux(k,i,j,jEq,s1,e1)
!                ENDDO
!
!              ELSE
!                DO k = 1, 3
!                  jEq = k+(iEq-1)*3
!                  ! Calculate the LDG flux for the stress tensor.
!                  myDGSEM % stressTensor % boundaryFlux(i,j,jEq,s1,e1) = 0.5_prec*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1)/&
!                    (myDGSEM % state % boundarySolution(i,j,4,s1,e1) +&
!                    myDGSEM % static % boundarySolution(i,j,4,s1,e1))+&
!                    myDGSEM % state % boundarySolution(ii,jj,iEq,s2,e2)/&
!                    (myDGSEM % state % boundarySolution(ii,jj,4,s2,e2) +&
!                    myDGSEM % static % boundarySolution(ii,jj,4,s2,e2)) )*&
!                    myDGSEM % mesh % elements % nHat(k,i,j,s1,e1)
!
!                  myDGSEM % stressTensor % boundaryFlux(ii,jj,jEq,s2,e2) = -myDGSEM % stressTensor % boundaryFlux(i,j,jEq,s1,e1)
!                ENDDO
!              ENDIF
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

    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
      4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
      1 )
    grid = dim3(myDGSEM % mesh % faces % nFaces,1,1)

    CALL BoundaryFaceFlux_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mesh % faces % elementIDs_dev, &
                                                        myDGSEM % mesh % faces % elementSides_dev, &
                                                        myDGSEM % mesh % faces % boundaryID_dev, &
                                                        myDGSEM % mesh % faces % iMap_dev, &
                                                        myDGSEM % mesh % faces % jMap_dev, &
                                                        myDGSEM % mesh % elements % nHat_dev, &
                                                        myDGSEM % state % boundarySolution_dev, &
                                                        myDGSEM % static % boundarySolution_dev, &
                                                        myDGSEM % state % externalState_dev, &
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
    REAL(prec) :: jump(1:myDGSEM % state % nEquations-1), aS(1:myDGSEM % state % nEquations-1)
    REAL(prec) :: fac, hCapRatio, rC

    !$OMP DO PRIVATE( jump, aS )
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
              ( (myDGSEM % state % externalState(ii,jj,6,bID)+&
              myDGSEM % static % boundarySolution(i,j,6,s1,e1) )/&
              myDGSEM % params % P0 )**myDGSEM % params % rC   )

            T = (myDGSEM % static % boundarySolution(i,j,5,s1,e1)+myDGSEM % state % boundarySolution(i,j,5,s1,e1))/&
                (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % boundarySolution(i,j,4,s1,e1))

            cIn  = sqrt( myDGSEM % params % R*T* &
              ( (myDGSEM % state % boundarySolution(i,j,6,s1,e1)+&
              myDGSEM % static % boundarySolution(i,j,6,s1,e1) )/&
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
              aS(k) = aS(k) + (myDGSEM % state % boundarySolution(i,j,6,s1,e1) + &
                myDGSEM % state % externalState(ii,jj,6,bID))*nHat(k)
            ENDDO


            DO iEq = 1, myDGSEM % state % nEquations-1
              myDGSEM % state % boundaryFlux(i,j,iEq,s1,e1) = 0.5_prec*( aS(iEq) - fac*jump(iEq) )*norm
              
!              DO k = 1, 3
!                jEq = k + (iEq-1)*3
!                     IF( iEq == 4 )THEN
!                           ! Calculate the Bassi-Rebay (LDG) flux for the stress tensor.
!                           myDGSEM % stressTensor % boundaryFlux(i,j,jEq,s1,e1) = 0.5_prec*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) +&
!                                                                                             myDGSEM % state % externalState(ii,jj,iEq,bID) )*&
!                                                                                             myDGSEM % mesh % elements % nHat(k,i,j,s1,e1)
!                                                                                        
!                     ELSE
!                           ! Calculate the Bassi-Rebay (LDG) flux for the stress tensor.
!                           myDGSEM % stressTensor % boundaryFlux(i,j,jEq,s1,e1) = 0.5_prec*( myDGSEM % state % boundarySolution(i,j,iEq,s1,e1)/&
!                                                                                             (myDGSEM % state % boundarySolution(i,j,4,s1,e1) +&
!                                                                                              myDGSEM % static % boundarySolution(i,j,4,s1,e1))+&
!                                                                                             myDGSEM % state % externalState(ii,jj,iEq,bID)/&
!                                                                                             (myDGSEM % state % externalState(ii,jj,4,bID) +&
!                                                                                              myDGSEM % static % boundarySolution(i,j,4,s1,e1)) )*&
!                                                                                             myDGSEM % mesh % elements % nHat(k,i,j,s1,e1)
!                                                                                        
!                     ENDIF
!                     
!              ENDDO
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

    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) )
    grid = dim3(myDGSEM % state % nEquations-1,myDGSEM % mesh % elements % nElements,1)

    CALL CalculateFlux_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                    myDGSEM % static % solution_dev, &
                                                    myDGSEM % mesh % elements % Ja_dev, &
                                                    myDGSEM % mesh % elements % J_dev, &
                                                    myDGSEM % state % flux_dev, &
                                                    myDGSEM % params % polyDeg_dev, &
                                                    myDGSEM % state % nEquations_dev, &
                                                    myDGSEM % mesh % elements % nElements_dev )

    CALL State_Mapped_DG_Divergence_3D_CUDAKernel<<<grid, tBlock>>>( myDGSEM % state % flux_dev, &
                                                               myDGSEM % state % boundaryFlux_dev, &
                                                               myDGSEM % state % fluxDivergence_dev, &
                                                               myDGSEM % dgStorage % boundaryInterpolationMatrix_dev, &
                                                               myDGSEM % dgStorage % dgDerivativeMatrixTranspose_dev, &
                                                               myDGSEM % dgStorage % quadratureWeights_dev, &
                                                               myDGSEM % mesh % elements % J_dev ) 

!    CALL DG_Divergence_3D_CUDAKernel<<<grid, tBlock>>>( myDGSEM % state % flux_dev, &
!                                                        myDGSEM % state % boundaryFlux_dev, &
!                                                        myDGSEM % state % fluxDivergence_dev, &
!                                                        myDGSEM % dgStorage % boundaryInterpolationMatrix_dev, &
!                                                        myDGSEM % dgStorage % dgDerivativeMatrixTranspose_dev, &
!                                                        myDGSEM % dgStorage % quadratureWeights_dev, &
!                                                        myDGSEM % params % polydeg_dev, &
!                                                        myDGSEM % state % nEquations_dev, &
!                                                        myDGSEM % state % nElements_dev )
!
    CALL CalculateSourceTerms_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                           myDGSEM % static % solution_dev, &
                                                           myDGSEM % state % source_dev, &
                                                           myDGSEM % sourceTerms % drag_dev, &
                                                           myDGSEM % params % fRotX_dev, &
                                                           myDGSEM % params % fRotY_dev, &
                                                           myDGSEM % params % fRotZ_dev, &
                                                           myDGSEM % params % g_dev, &
                                                           myDGSEM % params % polydeg_dev, &
                                                           myDGSEM % state % nEquations_dev, &
                                                           myDGSEM % mesh % elements % nElements_dev )
  
#else


    !$OMP DO
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
                  myDGSEM % state % solution(i,j,k,6,iEl)

                myDGSEM % state % flux(2,i,j,k,iEq,iEl) = myDGSEM % state % flux(2,i,j,k,iEq,iEl) + &
                  myDGSEM % mesh % elements % Ja(i,j,k,iEq,2,iEl)*&
                  myDGSEM % state % solution(i,j,k,6,iEl)

                myDGSEM % state % flux(3,i,j,k,iEq,iEl) = myDGSEM % state % flux(3,i,j,k,iEq,iEl) + &
                  myDGSEM % mesh % elements % Ja(i,j,k,iEq,3,iEl)*&
                  myDGSEM % state % solution(i,j,k,6,iEl)

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

                myDGSEM % state % source(i,j,k,1,iEl) = myDGSEM % sourceTerms % drag(i,j,k,iEl)*myDGSEM % state % solution(i,j,k,1,iEl)*F -&
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

                myDGSEM % state % source(i,j,k,2,iEl) = myDGSEM % sourceTerms % drag(i,j,k,iEl)*myDGSEM % state % solution(i,j,k,2,iEl)*F - &
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

                myDGSEM % state % source(i,j,k,3,iEl) = myDGSEM % sourceTerms % drag(i,j,k,iEl)*myDGSEM % state % solution(i,j,k,3,iEl)*F - &
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
!  This section of code CONTAINS routines for computing the gradients of the prognostic variables !
!                                                                                                 !
! /////////////////////////////////////////////////////////////////////////////////////////////// !
!
  SUBROUTINE CalculateStressTensor_Fluid( myDGSEM )

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock
#else
    INTEGER :: iEl, iEq, idir, i, j, k, m, jEq
    REAL(prec) :: F
#endif



#ifdef HAVE_CUDA
    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) )
    grid = dim3(myDGSEM % sgsCoeffs % nEquations,myDGSEM % mesh % elements % nElements,3)

    CALL CalculateStressTensorFlux_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                                myDGSEM % static % solution_dev, &
                                                                myDGSEM % mesh % elements % Ja_dev, &
                                                                myDGSEM % stressTensor % flux_dev, &
                                                                myDGSEM % params % polyDeg_dev, &
                                                                myDGSEM % state % nEquations_dev, &
                                                                myDGSEM % stressTensor % nEquations_dev, &
                                                                myDGSEM % mesh % elements % nElements_dev )

    grid = dim3(myDGSEM % state % nEquations,myDGSEM % mesh % elements % nElements,3)
    
    CALL Stress_Mapped_DG_Divergence_3D_CUDAKernel<<<grid, tBlock>>>( myDGSEM % stressTensor % flux_dev, &
                                                        myDGSEM % stressTensor % boundaryFlux_dev, &
                                                        myDGSEM % stressTensor % fluxDivergence_dev, &
                                                        myDGSEM % dgStorage % boundaryInterpolationMatrix_dev, &
                                                        myDGSEM % dgStorage % dgDerivativeMatrixTranspose_dev, &
                                                        myDGSEM % dgStorage % quadratureWeights_dev, &
                                                        myDGSEM % mesh % elements % J_dev )

    grid = dim3(myDGSEM % sgsCoeffs % nEquations,myDGSEM % mesh % elements % nElements,1)

    CALL CalculateStressTensor_CUDAKernel<<<grid,tBlock>>>( myDGSEM % stressTensor % solution_dev, &
                                                            myDGSEM % stressTensor % fluxDivergence_dev, &
                                                            myDGSEM % stressTensor % flux_dev, &
                                                            myDGSEM % sgsCoeffs % solution_dev, &
                                                            myDGSEM % state % solution_dev, &
                                                            myDGSEM % static % solution_dev, &
                                                            myDGSEM % mesh % elements % J_dev, &
                                                            myDGSEM % mesh % elements % Ja_dev, &
                                                            myDGSEM % params % polyDeg_dev, &
                                                            myDGSEM % state % nEquations_dev, &
                                                            myDGSEM % stressTensor % nEquations_dev, &
                                                            myDGSEM % sgsCoeffs % nEquations_dev, &
                                                            myDGSEM % mesh % elements % nElements_dev )

#else

    !$OMP DO
    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO iEq = 1, myDGSEM % state % nEquations-1
        DO idir = 1, 3

          jEq = idir + (iEq-1)*3

          DO k = 0, myDGSEM % params % polyDeg
            DO j = 0, myDGSEM % params % polyDeg
              DO i = 0, myDGSEM % params % polyDeg
                myDGSEM % stressTensor % flux(1,i,j,k,jEq,iEl) = 0.0_prec
                myDGSEM % stressTensor % flux(2,i,j,k,jEq,iEl) = 0.0_prec
                myDGSEM % stressTensor % flux(3,i,j,k,jEq,iEl) = 0.0_prec
              ENDDO
            ENDDO
          ENDDO

          ! Here the flux tensor in physical space is calculated and rotated to give the
          ! contravariant flux tensor in the reference computational DOmain.
          IF( iEq == 4 )THEN

            DO k = 0, myDGSEM % params % polyDeg
              DO j = 0, myDGSEM % params % polyDeg
                DO i = 0, myDGSEM % params % polyDeg

                  myDGSEM % stressTensor % flux(1,i,j,k,jEq,iEl) = myDGSEM % mesh % elements % Ja(i,j,k,idir,1,iEl)*&
                    myDGSEM % state % solution(i,j,k,iEq,iEl)

                  myDGSEM % stressTensor % flux(2,i,j,k,jEq,iEl) = myDGSEM % mesh % elements % Ja(i,j,k,idir,2,iEl)*&
                    myDGSEM % state % solution(i,j,k,iEq,iEl)

                  myDGSEM % stressTensor % flux(3,i,j,k,jEq,iEl) = myDGSEM % mesh % elements % Ja(i,j,k,idir,3,iEl)*&
                    myDGSEM % state % solution(i,j,k,iEq,iEl)

                ENDDO
              ENDDO
            ENDDO

          ELSE

            DO k = 0, myDGSEM % params % polyDeg
              DO j = 0, myDGSEM % params % polyDeg
                DO i = 0, myDGSEM % params % polyDeg

                  myDGSEM % stressTensor % flux(1,i,j,k,jEq,iEl) = myDGSEM % mesh % elements % Ja(i,j,k,idir,1,iEl)*&
                    myDGSEM % state % solution(i,j,k,iEq,iEl)/&
                    (myDGSEM % state % solution(i,j,k,4,iEl)+&
                    myDGSEM % static % solution(i,j,k,4,iEl) )

                  myDGSEM % stressTensor % flux(2,i,j,k,jEq,iEl) = myDGSEM % mesh % elements % Ja(i,j,k,idir,2,iEl)*&
                    myDGSEM % state % solution(i,j,k,iEq,iEl)/&
                    (myDGSEM % state % solution(i,j,k,4,iEl)+&
                    myDGSEM % static % solution(i,j,k,4,iEl) )

                  myDGSEM % stressTensor % flux(3,i,j,k,jEq,iEl) = myDGSEM % mesh % elements % Ja(i,j,k,idir,3,iEl)*&
                    myDGSEM % state % solution(i,j,k,iEq,iEl)/&
                    (myDGSEM % state % solution(i,j,k,4,iEl)+&
                    myDGSEM % static % solution(i,j,k,4,iEl) )

                ENDDO
              ENDDO
            ENDDO

          ENDIF

        ENDDO
      ENDDO
    ENDDO
    !$OMP ENDDO

    CALL myDGSEM % stressTensor % Calculate_Weak_Flux_Divergence( myDGSEM % dgStorage ) 

    !$OMP DO
    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO iEq = 1, myDGSEM % state % nEquations-1
        DO idir = 1, 3

          jEq = idir + (iEq-1)*3

          DO k = 0, myDGSEM % params % polyDeg
            DO j = 0, myDGSEM % params % polyDeg
              DO i = 0, myDGSEM % params % polyDeg

                myDGSEM % stressTensor % solution(i,j,k,jEq,iEl) = myDGSEM % stressTensor % fluxDivergence(i,j,k,jEq,iEl)/myDGSEM % mesh % elements % J(i,j,k,iEl)

                IF( iEq == 4 )THEN

                  F = myDGSEM % stressTensor % solution(i,j,k,jEq,iEl)*&
                      myDGSEM % sgsCoeffs % solution(i,j,k,iEq,iEl)

                ELSE

                  F = myDGSEM % stressTensor % solution(i,j,k,jEq,iEl)*&
                     (myDGSEM % state % solution(i,j,k,4,iEl)+&
                      myDGSEM % static % solution(i,j,k,4,iEl))*&
                     myDGSEM % sgsCoeffs % solution(i,j,k,iEq,iEl)

                ENDIF

                myDGSEM % stressTensor % flux(1,i,j,k,iEq,iEl) = myDGSEM % stressTensor % flux(1,i,j,k,iEq,iEl) +&
                                                                 myDGSEM % mesh % elements % Ja(i,j,k,idir,1,iEl)*F

                myDGSEM % stressTensor % flux(2,i,j,k,iEq,iEl) = myDGSEM % stressTensor % flux(2,i,j,k,iEq,iEl) +&
                                                                 myDGSEM % mesh % elements % Ja(i,j,k,idir,2,iEl)*F

                myDGSEM % stressTensor % flux(3,i,j,k,iEq,iEl) = myDGSEM % stressTensor % flux(3,i,j,k,iEq,iEl) +&
                                                                myDGSEM % mesh % elements % Ja(i,j,k,idir,3,iEl)*F

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
  SUBROUTINE UpdateExternalStress_Fluid( myDGSEM, tn ) ! ////////// !

    IMPLICIT NONE
    CLASS(Fluid), INTENT(inout) :: myDGSEM
    REAL(prec), INTENT(in)      :: tn
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                  1 )
    grid = dim3(myDGSEM % extComm % nBoundaries,myDGSEM % stressTensor % nEquations,1)

    CALL UpdateExternalStress_CUDAKernel<<<grid, tBlock>>>( myDGSEM % extComm % boundaryIDs_dev, &            ! I
                                                            myDGSEM % mesh % faces % elementIDs_dev, &         ! I
                                                            myDGSEM % mesh % faces % elementSides_dev, &       ! I
                                                            myDGSEM % extComm % extProcIDs_dev, &              ! I
                                                            myDGSEM % stressTensor % externalState_dev, &                    ! O
                                                            myDGSEM % stressTensor % boundarySolution_dev, &   ! I
                                                            myDGSEM % stressTensor % prescribedState_dev, &                  ! I
                                                            myDGSEM % mesh % elements % nHat_dev, &
                                                            myDGSEM % extComm % myRank_dev, &
                                                            myDGSEM % params % polyDeg_dev, &
                                                            myDGSEM % stressTensor % nEquations_dev, &
                                                            myDGSEM % extComm % nBoundaries_dev, &
                                                            myDGSEM % mesh % faces % nFaces_dev, &
                                                            myDGSEM % mesh % elements % nElements_dev )             ! I
#else
    ! Local
    INTEGER    :: iEl, iFace, bFaceID, i, j, k, iEq
    INTEGER    :: bID, p2
    INTEGER    :: e1, e2, s1, s2

    !$OMP DO
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
              myDGSEM % stressTensor % externalState(i,j,iEq,bID) = -myDGSEM % stressTensor % boundarySolution(i,j,iEq,s1,e1)
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

    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                  1 )
    grid = dim3(myDGSEM % mesh % faces % nFaces,myDGSEM % state % nEquations-1,1)

    CALL InternalStressFlux_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mesh % faces % elementIDs_dev, &
                                                          myDGSEM % mesh % faces % elementSides_dev, &
                                                          myDGSEM % mesh % faces % boundaryID_dev, &
                                                          myDGSEM % mesh % faces % iMap_dev, &
                                                          myDGSEM % mesh % faces % jMap_dev, &
                                                          myDGSEM % mesh % elements % nHat_dev, &
                                                          myDGSEM % state % boundarySolution_dev, &
                                                          myDGSEM % static % boundarySolution_dev, &
                                                          myDGSEM % stressTensor % boundarySolution_dev, &
                                                          myDGSEM % sgsCoeffs % boundarySolution_dev, &
                                                          myDGSEM % state % externalState_dev, &
                                                          myDGSEM % stressTensor % externalState_dev, &
                                                          myDGSEM % stressTensor % boundaryFlux_dev, &
                                                          myDGSEM % params % viscLengthScale_dev, &
                                                          myDGSEM % params % polyDeg_dev, &
                                                          myDGSEM % state % nEquations_dev, &
                                                          myDGSEM % stressTensor % nEquations_dev, &
                                                          myDGSEM % sgsCoeffs % nEquations_dev, & 
                                                          myDGSEM % extComm % nBoundaries_dev, &
                                                          myDGSEM % mesh % faces % nFaces_dev, &
                                                          myDGSEM % mesh % elements % nElements_dev )
#else
    ! Local
    REAL(prec) :: norm, rhoIn, rhoOut
    INTEGER :: iEl, iFace
    INTEGER    :: i, j, k, m, iEq, jEq
    INTEGER    :: ii, jj, bID
    INTEGER    :: e1, s1, e2, s2

    !$OMP DO
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
                    myDGSEM % mesh % elements % nHat(m,i,j,s1,e1)
                ENDDO

              ELSE

                rhoOut = (myDGSEM % static % boundarySolution(ii,jj,4,s2,e2)+myDGSEM % state % boundarySolution(ii,jj,4,s2,e2) )
                rhoIn  = (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % boundarySolution(i,j,4,s1,e1) )

                DO m = 1, 3
                  jEq = m + (iEq-1)*3
                  myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) + &
                    0.5_prec*( rhoIn*myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*myDGSEM % stressTensor % boundarysolution(i,j,jEq,s1,e1)+&
                    rhoOut*myDGSEM % sgsCoeffs % boundarysolution(ii,jj,iEq,s2,e2)*myDGSEM % stressTensor % boundarysolution(ii,jj,jEq,s2,e2))*&
                    myDGSEM % mesh % elements % nHat(m,i,j,s1,e1)
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

    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                  1 )
    grid = dim3(myDGSEM % mesh % faces % nFaces,myDGSEM % state % nEquations-1,1)

    CALL BoundaryStressFlux_CUDAKernel<<<grid, tBlock>>>( myDGSEM % mesh % faces % elementIDs_dev, &
                                                          myDGSEM % mesh % faces % elementSides_dev, &
                                                          myDGSEM % mesh % faces % boundaryID_dev, &
                                                          myDGSEM % mesh % faces % iMap_dev, &
                                                          myDGSEM % mesh % faces % jMap_dev, &
                                                          myDGSEM % mesh % elements % nHat_dev, &
                                                          myDGSEM % state % boundarySolution_dev, &
                                                          myDGSEM % static % boundarySolution_dev, &
                                                          myDGSEM % stressTensor % boundarySolution_dev, &
                                                          myDGSEM % sgsCoeffs % boundarySolution_dev, &
                                                          myDGSEM % sgsCoeffs % externalState_dev, &
                                                          myDGSEM % state % externalState_dev, &
                                                          myDGSEM % stressTensor % externalState_dev, &
                                                          myDGSEM % stressTensor % boundaryFlux_dev, &
                                                          myDGSEM % params % viscLengthScale_dev, &
                                                          myDGSEM % params % polyDeg_dev, &
                                                          myDGSEM % state % nEquations_dev, &
                                                          myDGSEM % stressTensor % nEquations_dev, &
                                                          myDGSEM % sgsCoeffs % nEquations_dev, & 
                                                          myDGSEM % extComm % nBoundaries_dev, &
                                                          myDGSEM % mesh % faces % nFaces_dev, &
                                                          myDGSEM % mesh % elements % nElements_dev )
#else
    ! Local
    REAL(prec) :: norm, rhoIn, rhoOut
    INTEGER :: iEl, iFace
    INTEGER    :: i, j, k, m, iEq, jEq
    INTEGER    :: ii, jj, bID
    INTEGER    :: e1, s1, e2, s2

    !$OMP DO
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

            bID  = myDGSEM % mesh % faces % boundaryID(iFace)
            IF( bID < 0 )THEN ! Physical Boundary

              bID = ABS(bID)
              DO iEq = 1, myDGSEM % state % nEquations-1
                myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = &
                  0.5_prec*myDGSEM % sgsCoeffs % boundarySolution(i,j,iEq,s1,e1)*&
                  ( myDGSEM % state % externalState(ii,jj,iEq,bID) - myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) )/&
                  myDGSEM % params % viscLengthScale*norm

                IF( iEq == 4 )THEN

                  DO m = 1, 3
                    jEq = m + (iEq-1)*3
                    myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) + &
                                                                           0.5_prec*myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*&
                                                                           ( myDGSEM % stressTensor % boundarysolution(i,j,jEq,s1,e1) + &
                                                                             myDGSEM % stressTensor % externalState(ii,jj,jEq,bID) )*&
                                                                           myDGSEM % mesh % elements % nHat(m,i,j,s1,e1)
                  ENDDO

                ELSE

                  rhoOut = (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % externalState(ii,jj,4,bID) )
                  rhoIn  = (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % boundarySolution(i,j,4,s1,e1) )

                  DO m = 1, 3
                    jEq = m + (iEq-1)*3
                    myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) + &
                                                                           0.5_prec*myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*&
                                                                          ( rhoIn*myDGSEM % stressTensor % boundarysolution(i,j,jEq,s1,e1)+&
                                                                            rhoOut*myDGSEM % stressTensor % externalState(ii,jj,jEq,bID) )*&
                                                                          myDGSEM % mesh % elements % nHat(m,i,j,s1,e1)
                  ENDDO

                ENDIF
              ENDDO

            ELSE ! Neighboring process

              DO iEq = 1, myDGSEM % state % nEquations-1
                myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = &
                  0.5_prec*(myDGSEM % sgsCoeffs % externalState(ii,jj,iEq,bID)*myDGSEM % state % externalState(ii,jj,iEq,bID) -&
                  myDGSEM % sgsCoeffs % boundarySolution(i,j,iEq,s1,e1)*myDGSEM % state % boundarySolution(i,j,iEq,s1,e1) )/&
                  myDGSEM % params % viscLengthScale*norm

                IF( iEq == 4 )THEN

                  DO m = 1, 3
                    jEq = m + (iEq-1)*3
                    myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) + &
                                                                           0.5_prec*( myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*&
                                                                           myDGSEM % stressTensor % boundarysolution(i,j,jEq,s1,e1)+ &
                                                                           myDGSEM % sgsCoeffs % externalState(ii,jj,iEq,bID )*&
                                                                           myDGSEM % stressTensor % externalState(ii,jj,jEq,bID) )*&
                                                                           myDGSEM % mesh % elements % nHat(m,i,j,s1,e1)
                  ENDDO

                ELSE

                  rhoOut = (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % externalState(ii,jj,4,bID) )
                  rhoIn  = (myDGSEM % static % boundarySolution(i,j,4,s1,e1)+myDGSEM % state % boundarySolution(i,j,4,s1,e1) )

                  DO m = 1, 3

                    jEq = m + (iEq-1)*3

                    myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) = myDGSEM % stressTensor % boundaryFlux(i,j,iEq,s1,e1) + &
                                                                           0.5_prec*( myDGSEM % sgsCoeffs % boundarysolution(i,j,iEq,s1,e1)*&
                                                                           rhoIn*myDGSEM % stressTensor % boundarysolution(i,j,jEq,s1,e1)+&
                                                                           myDGSEM % sgsCoeffs % externalState(ii,jj,iEq,bID)*&
                                                                           rhoOut*myDGSEM % stressTensor % externalState(ii,jj,jEq,bID) )*&
                                                                           myDGSEM % mesh % elements % nHat(m,i,j,s1,e1)

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
  SUBROUTINE EquationOfState_Fluid( myDGSEM )
    ! This routine calculates the anomalous pressure referenced to the static state.
    ! The pressure is calculated using the ideal gas law.
    IMPLICIT NONE
    CLASS( Fluid ), INTENT(inout) :: myDGSEM
#ifdef HAVE_CUDA
    ! Local
    TYPE(dim3) :: grid, tBlock

    ! How should we pick the thread and block size
    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) )
    grid = dim3(myDGSEM % mesh % elements % nElements,1,1)

    CALL EquationOfState_CUDAKernel<<<grid,tBlock>>>( myDGSEM % state % solution_dev, &
                                                      myDGSEM % static % solution_dev, &
                                                      myDGSEM % params % P0_dev, &
                                                      myDGSEM % params % R_dev, &
                                                      myDGSEM % params % hCapRatio_dev, &
                                                      myDGSEM % params % polyDeg_dev, &
                                                      myDGSEM % state % nEquations_dev, &
                                                      myDGSEM % mesh % elements % nElements_dev )

#else
    ! Local
    INTEGER :: iEl, i, j, k
    REAL(prec) :: hCapRatio, rC, rhoT


    !$OMP DO
    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO k = 0, myDGSEM % params % polyDeg
        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg

            ! Pressure = rho*e*R/Cv (Ideal Gas Law, P = rho*R*T, thermo ~> T = theta*(P/P0)^r)
            ! THEN P = (rho*theta*R/P0^rC)^(Cp/Cv)
            ! And P' = P - P_static
            rhoT = (myDGSEM % static % solution(i,j,k,5,iEl) + myDGSEM % state % solution(i,j,k,5,iEl) )
            myDGSEM % state % solution(i,j,k,6,iEl) = myDGSEM % params % P0*( rhoT*myDGSEM % params % R/myDGSEM % params % P0 )**myDGSEM % params % hCapRatio -&
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

    PRINT*, dTdz, AlmostEqual( dTdz,0.0_prec )
 
    ! /////////////////////  Build the Static/Background State ///////////////////////// !

    !$OMP DO
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
    !$OMP ENDDO

    !$OMP DO
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
    DO iEl = 1, myDGSEM % mesh % elements % nElements
      myDGSEM % static % solution(:,:,:,6,iEl) = myDGSEM % state % solution(:,:,:,6,iEl)
    ENDDO
    !$OMP END MASTER

#ifdef HAVE_CUDA
    myDGSEM % static % solution_dev = myDGSEM % static % solution
    myDGSEM % state % solution_dev  = 0.0_prec
    istat = cudaDeviceSynchronize( )
#endif

    !$OMP MASTER
    myDGSEM % state % solution = 0.0_prec
    !$OMP END MASTER


  END SUBROUTINE CalculateStaticState_Fluid
!
  SUBROUTINE WriteTecplot_Fluid( myDGSEM )

    IMPLICIT NONE

    CLASS( Fluid ), INTENT(inout) :: myDGsem
    !LOCAL
    REAL(prec)  :: x(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:3,1:myDGSEM % mesh % elements % nElements)
    REAL(prec)  :: sol(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:myDGSEM % state % nEquations,1:myDGSEM % mesh % elements % nElements)
    REAL(prec)  :: bsol(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:myDGSEM % state % nEquations, 1:myDGSEM % mesh % elements % nElements)
#ifdef HAVE_CUDA
!    REAL(prec), DEVICE  :: x_dev(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:3,1:myDGSEM % mesh % elements % nElements)
!    REAL(prec), DEVICE  :: sol_dev(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:myDGSEM % state % nEquations,1:myDGSEM % mesh % elements % nElements)
!    REAL(prec), DEVICE  :: bsol_dev(0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,0:myDGSEM % params % nPlot,1:myDGSEM % state % nEquations, 1:myDGSEM % mesh % elements % nElements)
!    INTEGER,  DEVICE     :: nDim_dev
    INTEGER :: istat
#endif
    INTEGER       :: i, j, k, iEl, iEq, fUnit
    CHARACTER(5)  :: zoneID
    CHARACTER(4)  :: rankChar
    REAL(prec)    :: hCapRatio, c, T
    CHARACTER(13) :: timeStampString

    timeStampString = TimeStamp( myDGSEM % simulationTime, 's' )

!#ifdef HAVE_CUDA
!
!
!    nDim_dev = 3
!
!    CALL myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % state % solution_dev, &
!                                                                     sol_dev, myDGSEM % state % nEquations_dev, &
!                                                                     myDGSEM % mesh % elements % nElements_dev )  
!
!    CALL myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % static % solution_dev, &
!                                                                     bsol_dev, myDGSEM % static % nEquations_dev, &
!                                                                     myDGSEM % mesh % elements % nElements_dev )  
!
!    CALL myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % elements % x_dev, &
!                                                                     x_dev, nDim_dev, &
!                                                                     myDGSEM % mesh % elements % nElements_dev )  
!
!    x = x_dev
!    sol = sol_dev
!    bsol = bsol_dev
!    istat = cudaDeviceSynchronize( )
!
!#else
!
!    CALL myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % state % solution, &
!                                                                     sol, myDGSEM % state % nEquations, &
!                                                                     myDGSEM % mesh % elements % nElements )  
!
!    CALL myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % static % solution, &
!                                                                     bsol, myDGSEM % static % nEquations, &
!                                                                     myDGSEM % mesh % elements % nElements )  
!
!    CALL myDGSEM % dgStorage % interp % ApplyInterpolationMatrix_3D( myDGSEM % mesh % elements % x, &
!                                                                     x, 3, &
!                                                                     myDGSEM % mesh % elements % nElements )  
!#endif

#ifdef HAVE_CUDA
    CALL myDGSEM % state % UpdateHost( )
    CALL myDGSEM % static % UpdateHost( )
    istat = cudaDeviceSynchronize( )

#endif
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

    WRITE(rankChar,'(I4.4)') myDGSEM % extComm % myRank

    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= 'State.'//rankChar//'.'//timeStampString//'.tec', &
      FORM='formatted', &
      STATUS='replace')

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "u", "v", "w", "rho", "Pot. Temp.", "Pressure",'//&
      ' "u_b", "v_b", "w_b", "rho_b", "Pot. Temp._b", "Pressure_b", "Drag", "c" '

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
            c = sqrt( myDGSEM % params % R*T*( ( sol(i,j,k,6,iEl) + bsol(i,j,k,6,iEl) )/myDGSEM % params % P0 )**myDGSEM % params % hCapRatio   )
            WRITE(fUnit,'(17(E15.7,1x))') x(i,j,k,1,iEl), x(i,j,k,2,iEl), x(i,j,k,3,iEl),&
              sol(i,j,k,1,iEl)/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ), &
              sol(i,j,k,2,iEl)/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ), &
              sol(i,j,k,3,iEl)/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ), &
              sol(i,j,k,4,iEl), &
              (sol(i,j,k,5,iEl) + bsol(i,j,k,5,iEl))/( sol(i,j,k,4,iEl) + bsol(i,j,k,4,iEl) ), &
              sol(i,j,k,6,iEl), &
              bsol(i,j,k,1,iEl)/( bsol(i,j,k,4,iEl) ), &
              bsol(i,j,k,2,iEl)/( bsol(i,j,k,4,iEl) ), &
              bsol(i,j,k,3,iEl)/( bsol(i,j,k,4,iEl) ), &
              bsol(i,j,k,4,iEl), &
              bsol(i,j,k,5,iEl)/( bsol(i,j,k,4,iEl) ),&
              bsol(i,j,k,6,iEl),&
              myDGSEM % sourceTerms % drag(i,j,k,iEl), c
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
        STATUS='REPLACE' )
      WRITE(diagUnits(1),*) '#TotalMass'

      OPEN( UNIT=NewUnit(diagUnits(2)), &
        FILE='KineticEnergy.'//timeStampString//'.curve', &
        FORM='FORMATTED', &
        STATUS='REPLACE' )
      WRITE(diagUnits(2),*) '#TotalKineticEnergy'

      OPEN( UNIT=NewUnit(diagUnits(3)), &
        FILE='PotentialEnergy.'//timeStampString//'.curve', &
        FORM='FORMATTED', &
        STATUS='REPLACE' )
      WRITE(diagUnits(3),*) '#TotalPotentialEnergy'

      OPEN( UNIT=NewUnit(diagUnits(4)), &
        FILE='Heat.'//timeStampString//'.curve', &
        FORM='FORMATTED', &
        STATUS='REPLACE' )
      WRITE(diagUnits(4),*) '#TotalHeat'

      OPEN( UNIT=NewUnit(diagUnits(5)), &
        FILE='Volume.'//timeStampString//'.curve', &
        FORM='FORMATTED', &
        STATUS='REPLACE' )
      WRITE(diagUnits(5),*) '#TotalVolume'
    ENDIF

  END SUBROUTINE OpenDiagnosticsFiles_Fluid
!
  SUBROUTINE WriteDiagnostics_Fluid( myDGSEM )
    IMPLICIT NONE
    CLASS( Fluid ), INTENT(in) :: myDGSEM

    IF( myDGSEM % extComm % myRank == 0 )THEN
      WRITE(diagUnits(1),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % mass
      WRITE(diagUnits(2),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % KE
      WRITE(diagUnits(3),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % PE
      WRITE(diagUnits(4),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % heat
      WRITE(diagUnits(5),'(E15.5,2x,E15.5)') myDGSEM % simulationTime, myDGSEM % volume
    ENDIF

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
    INTEGER       :: thisRec, fUnit
    INTEGER       :: iEq, N
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

    PRINT*, 'S/R ReadPickup : Done.'

#ifdef HAVE_CUDA
    CALL myDGSEM % static % UpdateDevice( )
    CALL myDGSEM % state % UpdateDevice( )
    CALL myDGSEM % sourceTerms % UpdateDevice( )
    istat = cudaDeviceSynchronize( )

    tBlock = dim3(4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ), &
                  4*(ceiling( REAL(myDGSEM % params % polyDeg+1)/4 ) ) , &
                  myDGSEM % static % nEquations )
    grid = dim3(myDGSEM % static % nElements, 1, 1)  
    CALL CalculateFunctionsAtBoundaries_3D_CUDAKernel<<<grid, tBlock>>>( myDGSEM % static % solution_dev, &
                                                                         myDGSEM % static % boundarySolution_dev,  &
                                                                         myDGSEM % dgStorage % boundaryInterpolationMatrix_dev, &
                                                                         myDGSEM % params % polydeg_dev, myDGSEM % static % nEquations_dev,&
                                                                         myDGSEM % static % nElements_dev )
#else

    ! Interpolate the static state to the element boundaries
    CALL myDGSEM % static % Calculate_Solution_At_Boundaries( myDGSEM % dgStorage )

#endif

  END SUBROUTINE ReadPickup_Fluid
!
#ifdef HAVE_CUDA
! ============================================================================================================================ !
!------------------------------------------- CUDA Kernels Below -------------------------------------------------------------- !
! ============================================================================================================================ !
  ATTRIBUTES(Global) SUBROUTINE UpdateG3D_CUDAKernel( G3D, a, g, dt, solution, fluxDivergence, source, diffusiveFluxDivergence, Jac, N, nEq, nDiffEq, nElements )
  
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)       :: N, nEq, nDiffEq, nElements 
    REAL(prec), DEVICE, INTENT(inout) :: G3D(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)    :: a, g, dt
    REAL(prec), DEVICE, INTENT(inout) :: solution(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)    :: source(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)    :: fluxDivergence(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)    :: diffusiveFluxDivergence(0:N,0:N,0:N,1:nDiffEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)    :: Jac(0:N,0:N,0:N,1:nElements)
    ! Local
    INTEGER :: i, j, k, iEq, iEl
  
    iEl = blockIDx % x
    iEq = blockIDx % y
  
    i = threadIdx % x - 1
    j = threadIdx % y - 1
    k = threadIdx % z - 1

    IF( i <= N .AND. j <= N .AND. k <= N )THEN
  
      G3D(i,j,k,iEq,iEl)      = a*G3D(i,j,k,iEq,iEl) - fluxDivergence(i,j,k,iEq,iEl) + ( diffusiveFluxDivergence(i,j,k,iEq,iEl) )/Jac(i,j,k,iEl) + source(i,j,k,iEq,iEl)
  
      solution(i,j,k,iEq,iEl) = solution(i,j,k,iEq,iEl) + dt*g*G3D(i,j,k,iEq,iEl)

    ENDIF
  
  END SUBROUTINE UpdateG3D_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE CalculateSGSCoefficients_CUDAKernel( solution, static, smoothState, filterMat, sgsCoeffs, viscLengthScale, N, nEq, nElements  )
  
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)       :: N, nEq, nElements 
    REAL(prec), DEVICE, INTENT(in)    :: viscLengthScale
    REAL(prec), DEVICE, INTENT(in)    :: solution(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)    :: static(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(inout) :: smoothState(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)    :: filterMat(0:N,0:N)
    REAL(prec), DEVICE, INTENT(inout) :: sgsCoeffs(0:N,0:N,0:N,1:nEq-1,1:nElements)
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
    DO m = 1, nEq-1
      sgsCoeffs(i,j,k,m,iEl) = 0.09_prec*viscLengthScale*sqrt( sgsKE )
    ENDDO
  
  END SUBROUTINE CalculateSGSCoefficients_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE UpdateExternalSGSCoeffs_CUDAKernel( boundaryIDs, elementIDs, elementSides, procIDs, &
                                                                    externalsgsCoeffs, sgsCoeffsBsols, nHat, myRank, N, &
                                                                    nSGSEq, nBoundaryFaces, nFaces, nElements )
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: myRank, N, nSGSEq, nBoundaryFaces, nFaces, nElements
    INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nBoundaryFaces)
    INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: procIDs(1:nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(out) :: externalsgsCoeffs(0:N,0:N,1:nSGSEq,1:nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(in)  :: sgsCoeffsBsols(0:N,0:N,1:nSGSEq,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: nhat(1:3,0:N,0:N,1:6,1:nElements)
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
    
    IF( i <= N .AND. j <= N )THEN
    
      IF( p2 == myRank )THEN
        externalsgsCoeffs(i,j,iEq,iFace2) = sgsCoeffsBsols(i,j,iEq,s1,e1)
      ENDIF
    
    ENDIF
  
  END SUBROUTINE UpdateExternalSGSCoeffs_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE UpdateExternalState_CUDAKernel( boundaryIDs, elementIDs, elementSides, procIDs, &
                                                                externalState, stateBsols, prescribedState, nHat, myRank, cDrag, dragScale, N, &
                                                                nEq, nBoundaryFaces, nFaces, nElements)
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: myRank, N, nEq, nBoundaryFaces, nFaces, nElements
    INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nBoundaryFaces)
    INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: procIDs(1:nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(out) :: externalState(0:N,0:N,1:nEq,1:nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(in)  :: stateBsols(0:N,0:N,1:nEq,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: prescribedState(0:N,0:N,1:nEq,1:nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(in)  :: nhat(1:3,0:N,0:N,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: cDrag, dragScale
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
    
    IF( iFace <= nBoundaryFaces )THEN
    
      iFace2 = boundaryIDs( iFace ) ! Obtain the process-local face id for this boundary-face id
      e1     = elementIDs(1,iFace2)
      s1     = elementSides(1,iFace2)
      e2     = elementIDs(2,iFace2)
      p2     = procIDs( iFace )
    
      IF( i <= N .AND. j <= N .AND. p2 == myRank )THEN
    
        IF( e2 == PRESCRIBED )THEN
    
          DO iEq = 1, nEq
            externalState(i,j,iEq,iFace) = prescribedState(i,j,iEq,iFace)
          ENDDO
    
        ELSEIF( e2 == RADIATION )THEN
    
          DO iEq = 1, nEq
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
            (1.0_prec-Cdrag*dragScale*speed)
    
          ut = ( stateBsols(i,j,1,s1,e1)*tx + &
            stateBsols(i,j,2,s1,e1)*ty + &
            stateBsols(i,j,3,s1,e1)*tz )*&
            (1.0_prec-Cdrag*dragScale*speed)
    
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
 ATTRIBUTES(Global) SUBROUTINE CalculateFlux_CUDAKernel( solution, static, Ja, Jac, flux, N, nEq, nElements )
  
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)      :: N, nEq, nElements
    REAL(prec), DEVICE, INTENT(in)   :: solution(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)   :: static(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)   :: Ja(0:N,0:N,0:N,1:3,1:3,1:nElements)
    REAL(prec), DEVICE, INTENT(in)   :: Jac(0:N,0:N,0:N,1:nElements)
    REAL(prec), DEVICE, INTENT(out)  :: flux(1:3,0:N,0:N,0:N,1:nEq,1:nElements)
     ! Local
    INTEGER            :: i, j, k, row
    INTEGER            :: iEl, iEq
    REAL(prec)         :: tend, F
    
    iEq = blockIDx % x
    iEl = blockIDx % y
    
    i = threadIdx % x - 1
    j = threadIdx % y - 1
    k = threadIdx % z - 1

    flux(1,i,j,k,iEq,iEl) = 0.0_prec
    flux(2,i,j,k,iEq,iEl) = 0.0_prec
    flux(3,i,j,k,iEq,iEl) = 0.0_prec

    ! Here the flux tensor in physical space is calculated and rotated to give the
    ! contravariant flux tensor in the reference computational DOmain.

    !//////////////////////////////// Advection ///////////////////////////////////////!
    DO row = 1, 3

      F = solution(i,j,k,row,iEl)*( solution(i,j,k,iEq,iEl) + static(i,j,k,iEq,iEl) )/&
                                  ( solution(i,j,k,4,iEl) + static(i,j,k,4,iEl) )


      flux(1,i,j,k,iEq,iEl) = flux(1,i,j,k,iEq,iEl) + Ja(i,j,k,row,1,iEl)*F

      flux(2,i,j,k,iEq,iEl) = flux(2,i,j,k,iEq,iEl) + Ja(i,j,k,row,2,iEl)*F

      flux(3,i,j,k,iEq,iEl) = flux(3,i,j,k,iEq,iEl) + Ja(i,j,k,row,3,iEl)*F

    ENDDO

    ! //////////////////// Pressure (Momentum only) /////////////////////////// !

    IF( iEq <= 3 )THEN

      flux(1,i,j,k,iEq,iEl) = flux(1,i,j,k,iEq,iEl) + Ja(i,j,k,iEq,1,iEl)*solution(i,j,k,6,iEl)

      flux(2,i,j,k,iEq,iEl) = flux(2,i,j,k,iEq,iEl) + Ja(i,j,k,iEq,2,iEl)*solution(i,j,k,6,iEl)

      flux(3,i,j,k,iEq,iEl) = flux(3,i,j,k,iEq,iEl) + Ja(i,j,k,iEq,3,iEl)*solution(i,j,k,6,iEl)

    ENDIF
  
  END SUBROUTINE CalculateFlux_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE CalculateSourceTerms_CUDAKernel( solution, static, source, drag, fRotX, fRotY, fRotZ, g, N, nEq, nElements )
  
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: N, nEq, nElements
    REAL(prec), DEVICE, INTENT(in)  :: solution(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: static(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: drag(0:N,0:N,0:N,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: fRotX, fRotY, fRotZ, g
    REAL(prec), DEVICE, INTENT(out) :: source(0:N,0:N,0:N,1:nEq,1:nElements)
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

      source(i,j,k,1,iEl) = drag(i,j,k,iEl)*solution(i,j,k,1,iEl)*F-&
                            solution(i,j,k,3,iEl)*fRotY +&
                            solution(i,j,k,2,iEl)*fRotz
    
    ELSEIF( iEq == 2 )THEN

      source(i,j,k,2,iEl) = drag(i,j,k,iEl)*solution(i,j,k,2,iEl)*F -&
                            solution(i,j,k,1,iEl)*fRotZ +&
                            solution(i,j,k,3,iEl)*fRotX

    ELSEIF( iEq == 3 )THEN ! Add in the buoyancy acceleration

      source(i,j,k,3,iEl) = drag(i,j,k,iEl)*solution(i,j,k,3,iEl)*F -&
                            solution(i,j,k,2,iEl)*fRotX +&
                            solution(i,j,k,1,iEl)*fRotY -&
                            solution(i,j,k,4,iEl)*g

    ENDIF
  
  END SUBROUTINE CalculateSourceTerms_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE CalculateStressTensorFlux_CUDAKernel( solution, static, Ja, stressFlux, N, nEq, nDiffEq, nElements )
  
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: N, nEq, nDiffEq, nElements 
    REAL(prec), DEVICE, INTENT(in)  :: solution(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: static(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: Ja(0:N,0:N,0:N,1:3,1:3,1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: stressFlux(1:3,0:N,0:N,0:N,1:nDiffEq,1:nElements)
     ! Local
    INTEGER :: iEl, iEq, jEq, idir, i, j, k, m
    REAL(prec) :: strTens
    
    
    iEq = blockIDx % x
    iEl = blockIDx % y
    idir = blockIDx % z
    jEq = idir + (iEq-1)*3
    
    i = threadIDx % x-1
    j = threadIDx % y-1
    k = threadIDx % z-1
    
    IF( i <= N .AND. j <= N .AND. k <= N )THEN
      IF( iEq == 4 )THEN

        stressFlux(1,i,j,k,jEq,iEl) = Ja(i,j,k,idir,1,iEl)*solution(i,j,k,iEq,iEl)

        stressFlux(2,i,j,k,jEq,iEl) = Ja(i,j,k,idir,2,iEl)*solution(i,j,k,iEq,iEl)

        stressFlux(3,i,j,k,jEq,iEl) = Ja(i,j,k,idir,3,iEl)*solution(i,j,k,iEq,iEl)

      ELSE

        stressFlux(1,i,j,k,jEq,iEl) = Ja(i,j,k,idir,1,iEl)*solution(i,j,k,iEq,iEl)/&
          (solution(i,j,k,4,iEl)+static(i,j,k,4,iEl) )

        stressFlux(2,i,j,k,jEq,iEl) = Ja(i,j,k,idir,2,iEl)*solution(i,j,k,iEq,iEl)/&
          (solution(i,j,k,4,iEl)+static(i,j,k,4,iEl) )

        stressFlux(3,i,j,k,jEq,iEl) = Ja(i,j,k,idir,3,iEl)*solution(i,j,k,iEq,iEl)/&
          (solution(i,j,k,4,iEl)+static(i,j,k,4,iEl) )


      ENDIF

    ENDIF
  
  END SUBROUTINE CalculateStressTensorFlux_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE CalculateStressTensor_CUDAKernel( stressTensor, stressFluxDivergence, stressFlux, sgsCoeffs, state, static, Jac, Ja, N, nEq, nDiffEq, nSGSEq, nElements )
  
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: N, nEq, nDiffEq, nSGSEq, nElements 
    REAL(prec), DEVICE, INTENT(out) :: stressTensor(0:N,0:N,0:N,1:nDiffEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: stressFluxDivergence(0:N,0:N,0:N,1:nDiffEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: sgsCoeffs(0:N,0:N,0:N,1:nSGSEq,1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: stressFlux(1:3,0:N,0:N,0:N,1:nDiffEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: state(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: static(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: Jac(0:N,0:N,0:N,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: Ja(0:N,0:N,0:N,1:3,1:3,1:nElements)
    ! Local
    INTEGER :: iEl, iEq, jEq, idir, i, j, k
    REAL(prec) :: F

    iEq = blockIDx % x
    iEl = blockIDx % y

    i = threadIDx % x-1
    j = threadIDx % y-1
    k = threadIDx % z-1

    stressFlux(1:3,i,j,k,iEq,iEl) = 0.0_prec

    DO idir = 1, 3

      jEq = idir + (iEq-1)*3
      stressTensor(i,j,k,jEq,iEl) = stressFluxDivergence(i,j,k,jEq,iEl)/Jac(i,j,k,iEl)
  
      IF( iEq == 4 )THEN
  
        F = stressTensor(i,j,k,jEq,iEl)*&
            sgsCoeffs(i,j,k,iEq,iEl)
  
      ELSE
  
        F = stressTensor(i,j,k,jEq,iEl)*&
           (state(i,j,k,4,iEl)+&
            static(i,j,k,4,iEl))*&
           sgsCoeffs(i,j,k,iEq,iEl)
  
      ENDIF
  
      stressFlux(1,i,j,k,iEq,iEl) = stressFlux(1,i,j,k,iEq,iEl) + Ja(i,j,k,idir,1,iEl)*F
  
      stressFlux(2,i,j,k,iEq,iEl) = stressFlux(2,i,j,k,iEq,iEl) + Ja(i,j,k,idir,2,iEl)*F
  
      stressFlux(3,i,j,k,iEq,iEl) = stressFlux(3,i,j,k,iEq,iEl) + Ja(i,j,k,idir,3,iEl)*F

    ENDDO

  END SUBROUTINE CalculateStressTensor_CUDAKernel

  ATTRIBUTES(Global) SUBROUTINE UpdateExternalStress_CUDAKernel( boundaryIDs, elementIDs, elementSides, procIDs, &
                                                                 externalStress, stressBsols, prescribedStress, nHat, &
                                                                 myRank, N, nDiffEq, nBoundaryFaces, nFaces, nElements )
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: myRank, N, nDiffEq, nBoundaryFaces, nFaces, nElements
    INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nBoundaryFaces)
    INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: procIDs(1:nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(out) :: externalStress(0:N,0:N,1:nDiffEq,1:nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(in)  :: stressBsols(0:N,0:N,1:nDiffEq,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: prescribedStress(0:N,0:N,1:nDiffEq,1:nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(in)  :: nhat(1:3,0:N,0:N,1:6,1:nElements)
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
    
    IF( i <= N .AND. j <= N )THEN
    
      IF( p2 == myRank )THEN
        externalStress(i,j,iEq,iFace) = -stressBsols(i,j,iEq,s1,e1)
      ENDIF
    
    ENDIF
  
  END SUBROUTINE UpdateExternalStress_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE InternalStressFlux_CUDAKernel( elementIDs, elementSides, boundaryIDs, iMap, jMap, &
                                                               nHat, boundaryState, static, boundaryStress, sgsCoeffs, &
                                                               externalState, externalStress, boundaryFlux, viscLengthScale, &
                                                               N, nEq, nDiffEq, nSGSEq, nBoundaryFaces, nFaces, nElements )
  
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: N, nEq, nDiffEq, nSGSEq, nBoundaryFaces, nFaces, nElements 
    REAL(prec), DEVICE, INTENT(in)  :: viscLengthScale
    INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: iMap(0:N,0:N,1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: jMap(0:N,0:N,1:nFaces)
    REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:N,0:N,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryState(0:N,0:N,1:nEq,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryStress(0:N,0:N,1:nDiffEq,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: sgsCoeffs(0:N,0:N,1:nSGSEq,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: externalState(0:N,0:N,1:nEq,1:nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(in)  :: static(0:N,0:N,1:nEq,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: externalStress(0:N,0:N,1:nDiffEq,1:nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(out) :: boundaryFlux(0:N,0:N,1:nDiffEq,1:6,1:nElements)
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
        sgsCoeffs(i,j,iEq,s1,e1)*boundaryState(i,j,iEq,s1,e1))/viscLengthScale*norm
    
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
                                                               externalState, externalStress, boundaryFlux, viscLengthScale, &
                                                               N, nEq, nDiffEq, nSGSEq, nBoundaryFaces, nFaces, nElements  )
  
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: N, nEq, nDiffEq, nSGSEq, nBoundaryFaces, nFaces, nElements 
    REAL(prec), DEVICE, INTENT(in)  :: viscLengthScale
    INTEGER, DEVICE, INTENT(in)     :: elementIDs(1:2,1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: elementSides(1:2,1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: boundaryIDs(1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: iMap(0:N,0:N,1:nFaces)
    INTEGER, DEVICE, INTENT(in)     :: jMap(0:N,0:N,1:nFaces)
    REAL(prec), DEVICE, INTENT(in)  :: nHat(1:3,0:N,0:N,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryState(0:N,0:N,1:nEq,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryStress(0:N,0:N,1:nDiffEq,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: sgsCoeffs(0:N,0:N,1:nSGSEq,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: externalSGS(0:N,0:N,1:nSGSEq,nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(in)  :: externalState(0:N,0:N,1:nEq,1:nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(in)  :: static(0:N,0:N,1:nEq,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: externalStress(0:N,0:N,1:nDiffEq,1:nBoundaryFaces)
    REAL(prec), DEVICE, INTENT(out) :: boundaryFlux(0:N,0:N,1:nDiffEq,1:6,1:nElements)
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
          (externalState(ii,jj,iEq,bID)-boundaryState(i,j,iEq,s1,e1))/viscLengthScale*norm
    
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
          sgsCoeffs(i,j,iEq,s1,e1)*boundaryState(i,j,iEq,s1,e1))/viscLengthScale*norm
    
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
  ATTRIBUTES(Global) SUBROUTINE EquationOfState_CUDAKernel( solution, static, P0, R, hCapRatio, N, nEq, nElements )
    ! This routine calculates the anomalous pressure referenced to the static state.
    ! The pressure is calculated using the ideal gas law.
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)       :: N, nEq, nElements
    REAL(prec), DEVICE, INTENT(in)    :: P0, R, hCapRatio
    REAL(prec), DEVICE, INTENT(inout) :: solution(0:N,0:N,0:N,1:nEq,1:nElements)
    REAL(prec), DEVICE, INTENT(in)    :: static(0:N,0:N,0:N,1:nEq,1:nElements)
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
    solution(i,j,k,6,iEl) = P0*( rhoT*R/P0 )**hCapRatio - static(i,j,k,6,iEl)
  
  END SUBROUTINE EquationOfState_CUDAKernel
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
    
    
      IF( i <= polyDeg_dev .AND. j <= polyDeg_dev .AND. k <= polyDeg_dev )THEN
      
        fLocal(1:3,i,j,k) = f(1:3,i,j,k,iVar,iEl)
        
      ENDIF
      
      CALL syncthreads( )
      
      IF( i <= polyDeg_dev .AND. j <= polyDeg_dev .AND. k <= polyDeg_dev )THEN
      
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
                                      
      ENDIF

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
    
    
      IF( i <= polyDeg_dev .AND. j <= polyDeg_dev .AND. k <= polyDeg_dev )THEN
      
        fLocal(1:3,i,j,k) = f(1:3,i,j,k,iVar,iEl)
        
      ENDIF
      
      CALL syncthreads( )
      
      IF( i <= polyDeg_dev .AND. j <= polyDeg_dev .AND. k <= polyDeg_dev )THEN
      
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
                                      
      ENDIF

  END SUBROUTINE Stress_Mapped_DG_Divergence_3D_CUDAKernel
  ATTRIBUTES(Global) SUBROUTINE CalculateStateAtBoundaries_3D_CUDAKernel( f, fAtBoundaries, boundaryMatrix ) 
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
      
  END SUBROUTINE CalculateStateAtBoundaries_3D_CUDAKernel
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
  ATTRIBUTES(Global) SUBROUTINE CalculateStressAtBoundaries_3D_CUDAKernel( f, fAtBoundaries, boundaryMatrix ) 
    IMPLICIT NONE
    REAL(prec), DEVICE, INTENT(in)  :: f(0:polydeg_dev,0:polydeg_dev,0:polydeg_dev,1:nStress_dev,1:nEl_dev)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryMatrix(0:polydeg_dev,0:1)
    REAL(prec), DEVICE, INTENT(out) :: fAtBoundaries(0:polydeg_dev,0:polydeg_dev,1:nStress_dev,1:6,1:nEl_dev)
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
      
  END SUBROUTINE CalculateStressAtBoundaries_3D_CUDAKernel
#endif


END MODULE Fluid_Class


