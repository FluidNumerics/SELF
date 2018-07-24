! SELF_Fluids_Driver.f90
!
! Copyright 2018 Joseph Schoonover <joe@myFluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


PROGRAM SELF_Fluids_Driver

  USE ModelPrecision
  USE ModelParameters_Class
  USE Fluid_EquationParser_Class
  USE Fluid_Class

  IMPLICIT NONE

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !

  TYPE( Fluid )                :: myFluid
  TYPE( Fluid_EquationParser ) :: myFluidConditions
  LOGICAL                      :: setupSuccess
  LOGICAL                      :: initializeFromScratch
  LOGICAL                      :: pickupFileExists


! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !


  CALL Setup( )

  IF( setupSuccess )THEN

    CALL Initialize( )

    CALL MainLoop( )

    CALL Cleanup( )

  ENDIF

CONTAINS

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Setup( )
    IMPLICIT NONE

    CALL myFluid % Build( setupSuccess )

    IF( .NOT. setupSuccess )THEN
      RETURN
    ENDIF


    PRINT(MsgFMT), 'Setup Complete'

  END SUBROUTINE Setup

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !
  SUBROUTINE Initialize( )

    ! Attempt to read the fluid pickup file. If it doesn't exist, this routine
    ! returns FALSE.
    CALL myFluid % ReadPickup( pickupFileExists )

    ! If the pickup file doesn't exist, then the initial conditions are generated
    ! from the equation parser.
    IF( .NOT. pickupFileExists )THEN

      PRINT(MsgFMT), 'Pickup file not found.'
      PRINT(MsgFMT), 'Attempting initial condition generation from self.equations'

      CALL myFluidConditions % Build( 'self.equations' )
      CALL InitialConditions( )

      CALL myFluid % WritePickup( )
      CALL myFluid % WriteTecplot( )

    ENDIF

#ifdef HAVE_DIAGNOSTICS
    CALL myFluid % WriteDiagnostics( )
#endif

  END SUBROUTINE Initialize

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE InitialConditions( )
    ! Local
    INTEGER    :: i, j, k, iEl
    REAL(prec) :: x(1:3)
    REAL(prec) :: T, Tbar, u, v, w, rho, rhobar


    myFluid % state % solution = 0.0_prec
    !$OMP PARALLEL
    CALL myFluid % CalculateStaticState( ) !! CPU Kernel
    !$OMP END PARALLEL

    DO iEl = 1, myFluid % mesh % elements % nElements
      DO k = 0, myFluid % params % polyDeg
        DO j = 0, myFluid % params % polyDeg
          DO i = 0, myFluid % params % polyDeg

            x(1:3) = myFluid % mesh % elements % x(i,j,k,1:3,iEl)

            IF( myFluidConditions % calculate_density_from_T )THEN
               
              u = myFluidConditions % u % evaluate( x ) 
              v = myFluidConditions % v % evaluate( x ) 
              w = myFluidConditions % w % evaluate( x ) 
              T = myFluidConditions % t % evaluate( x ) ! Potential temperature anomaly

              Tbar = myFluid % static % solution(i,j,k,5,iEl)/myFluid % static % solution(i,j,k,4,iEl)

              myFluid % state % solution(i,j,k,4,iEl) = -myFluid % static % solution(i,j,k,4,iEl)*T/(Tbar + T)
              myFluid % state % solution(i,j,k,1,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*u
              myFluid % state % solution(i,j,k,2,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*v
              myFluid % state % solution(i,j,k,3,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*w

            ELSE

              u   = myFluidConditions % u % evaluate( x ) 
              v   = myFluidConditions % v % evaluate( x ) 
              w   = myFluidConditions % w % evaluate( x ) 
              rho = myFluidConditions % rho % evaluate( x ) 
              T   = myFluidConditions % t % evaluate( x ) ! Potential temperature anomaly


              myFluid % state % solution(i,j,k,4,iEl) = rho
              myFluid % state % solution(i,j,k,1,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*u
              myFluid % state % solution(i,j,k,2,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*v
              myFluid % state % solution(i,j,k,3,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*w
              myFluid % state % solution(i,j,k,5,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*T


              myFluid % sourceTerms % drag(i,j,k,iEl) = myFluidConditions % drag % evaluate( x )

            ENDIF

          ENDDO
        ENDDO
      ENDDO
    ENDDO

#ifdef HAVE_CUDA
    myFluid % state % solution_dev = myFluid % state % solution
    CALL myFluid % sourceTerms % Update_Device( )
#endif

    !$OMP PARALLEL
    CALL myFluid % EquationOfState( )
    !$OMP END PARALLEL

    CALL myFluid % UpdateExternalStaticState( )

  END SUBROUTINE InitialConditions

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Cleanup( )
    IMPLICIT NONE

    CALL myFluid % Trash( )

  END SUBROUTINE Cleanup

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE MainLoop( )
    IMPLICIT NONE
    INTEGER    :: iT
! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !
    DO iT = 1, myFluid % params % nDumps ! Loop over time-steps

      CALL myFluid % ForwardStepRK3( myFluid % params % nStepsPerDump ) ! Forward Step

#ifdef HAVE_CUDA
      myFluid % state % solution = myFluid % state % solution_dev ! Update the host from the GPU
#endif

      CALL myFluid % WritePickup( )
      CALL myFluid % WriteTecplot( )

#ifdef HAVE_DIAGNOSTICS
      CALL myFluid % Diagnostics( )
      CALL myFluid % WriteDiagnostics( )
#endif

    ENDDO


  END SUBROUTINE MainLoop


END PROGRAM SELF_Fluids_Driver

