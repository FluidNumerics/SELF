! Fluid_Driver.f90
!
! Copyright 2018 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


PROGRAM Fluid_Driver

  USE ModelPrecision
  USE ModelParameters_Class
  USE Fluid_Class

  IMPLICIT NONE

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !

  TYPE( Fluid )       :: myeu
  LOGICAL             :: setupSuccess


! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !


  CALL Setup( )

  IF( setupSuccess )THEN

    CALL MainLoop( )

    CALL Cleanup( )

  ENDIF

CONTAINS

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Setup( )
    IMPLICIT NONE

    CALL myeu % Build( setupSuccess )

    IF( .NOT. setupSuccess )THEN
      RETURN
    ENDIF

    CALL myeu % WriteDiagnostics( )

    PRINT(MsgFMT), 'Setup Complete'

  END SUBROUTINE Setup

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Cleanup( )
    IMPLICIT NONE

    CALL myeu % Trash( )

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
    !$OMP PARALLEL
    DO iT = 1, myeu % params % nDumps ! Loop over time-steps

      CALL myeu % ForwardStepRK3( myeu % params % nStepsPerDump ) ! Forward Step

#ifdef HAVE_CUDA
      myeu % state % solution = myeu % state % solution_dev ! Update the host from the GPU
#endif

      !$OMP MASTER
      CALL myeu % WritePickup( )
      CALL myeu % WriteTecplot( )
      !$OMP END MASTER

      CALL myeu % Diagnostics( )
      CALL myeu % WriteDiagnostics( )

    ENDDO
    !$OMP END PARALLEL


  END SUBROUTINE MainLoop


END PROGRAM Fluid_Driver

