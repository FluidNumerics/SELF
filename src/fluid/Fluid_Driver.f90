! Fluid_Driver.f90
! 
! Copyright 2017 Joseph Schoonover <schoonover.numerics@gmail.com>, SimonsTechnical
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

 
PROGRAM Fluid_Driver

! src/common/
USE ModelPrecision
#ifdef TIMING
USE Timing
#endif
! src/highend/euler/mpi-cuda/
USE FluidParams_Class
USE Fluid_Class

 IMPLICIT NONE

#ifdef INSITU_VIZ
INCLUDE "visitfortransimV2interface.inc"
#endif


 TYPE( Fluid ) :: myeu
#ifdef TIMING
 TYPE( MultiTimers ) :: timers
#endif
 INTEGER :: mpiErr, myRank, nProcs
 INTEGER :: iter0, nT, dFreq, iT, nDumps
 REAL(prec) :: tn, deltaT
 
      CALL Setup( )

#ifdef INSITU_VIZ
      CALL SetupLibSim( )
#endif INSITU_VIZ

#ifdef TIMING
      !$OMP MASTER
      IF( myRank == 0 )THEN
         CALL timers % Build( )
         CALL timers % AddTimer( 'MainLoop', 1 )
         CALL timers % StartTimer( 1 )
      ENDIF
      !$OMP END MASTER
#endif
      CALL MainLoop( )
#ifdef TIMING
      !$OMP MASTER
      IF( myRank == 0 )THEN
         CALL timers % StopTimer( 1 )
         CALL timers % Write_MultiTimers( )
         CALL timers % Trash( )
      ENDIF
      !$OMP END MASTER
#endif

      CALL Cleanup( )


CONTAINS

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Setup( )
     IMPLICIT NONE

#ifdef HAVE_MPI
      ! MPI Initialization
      CALL MPI_INIT( mpiErr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myRank, mpiErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nProcs, mpiErr )
      ! Sanity check
!      PRINT*, 'Fluid_Driver : Greetings from Process ', myRank, ' of ',nProcs
#else
      myRank = 0
      nProcs = 1
#endif
      CALL myeu % Build( myRank, nProcs )

      iter0  = myeu % params % iterInit
      nT     = myeu % params % nTimeSteps
      dFreq  = myeu % params % dumpFreq
      deltaT = myeu % params % dt
      nDumps = (nT)/dFreq

#ifdef HAVE_MPI
      CALL MPI_BARRIER( MPI_COMM_WORLD, mpiErr )
#endif
  END SUBROUTINE Setup

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !
#ifdef INSITU_VIZ
  SUBROUTINE SetupLibSim( )
     IMPLICIT NONE
     INTEGER :: vizErr

       vizErr = visitSetupEnv( ) 
       vizErr = visitInitializeSim( "SELF-Fluids",11, &
                                    "A simulation of compressible fluids",35, &
                                    "/home/joe/SELF-Fluids/bin/", 26, &
                                    VISIT_F77NULLSTRING, VISIT_F77NULLSTRINGLEN, &
                                    VISIT_F77NULLSTRING, VISIT_F77NULLSTRINGLEN, &
                                    VISIT_F77NULLSTRING, VISIT_F77NULLSTRINGLEN )

  END SUBROUTINE SetupLibSim
#endif
! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Cleanup( )
     IMPLICIT NONE

      CALL myeu % Trash( )

#ifdef HAVE_MPI      
      CALL MPI_FINALIZE( mpiErr )
#endif
  END SUBROUTINE Cleanup

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE MainLoop( )
     IMPLICIT NONE
#ifdef INSITU_VIZ
     INTEGER    :: visitState, visitResult, runFlag, blocking, iterate     
     REAL(prec) :: x(0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     1:myeu % mesh % nElems )
     REAL(prec) :: y(0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     1:myeu % mesh % nElems )
     REAL(prec) :: z(0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     1:myeu % mesh % nElems )
     REAL(prec) :: u(0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     1:myeu % mesh % nElems )
     REAL(prec) :: v(0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     1:myeu % mesh % nElems )
     REAL(prec) :: w(0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     1:myeu % mesh % nElems )
     REAL(prec) :: rho(0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     1:myeu % mesh % nElems )
     REAL(prec) :: potentialTemp(0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     1:myeu % mesh % nElems )
     REAL(prec) :: pressure(0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     0:myeu % params % nPlot, &
                     1:myeu % mesh % nElems )


      CALL myeu % ObtainPlottingMesh( x, y, z )

      runFlag = 1
      iterate = iter0
      !$OMP PARALLEL
      DO ! Loop indefinitely; really breaks when visitState < 0

         IF( runFlag == 1 )THEN
            blocking = 0
         ELSE
            blocking = 1
         ENDIF
 
         visitState = visitDetectInput( blocking, -1 )

         IF( visitState == 0 )THEN

            tn = real( iT, prec )*deltaT ! time at the beginning of this timestep          
            CALL myeu % ForwardStepRK3( tn, 1, myRank ) ! Forward Step
            iterate = iterate + 1
#ifdef HAVE_CUDA
            myeu % state % solution = myeu % state % solution_dev ! Update the host from the GPU
#endif
            CALL myeu % FluidStateAtPlottingPoints( u, v, w, potentialTemp, pressure )

          ELSEIF( visitState == 1 )THEN
 
            runFlag = 0
            visitResult = visitAttemptConnection( )
            IF( visitResult == 1 )THEN
               PRINT*, 'VisIt connected!'
            ELSE
               PRINT*, 'VisIt did not connect!'
            ENDIF

          ELSEIF( visitState == 2 )THEN
 
            runFlag = 0
            IF( visitProcessEngineCommand() == 0 )THEN
               visitResult = visitDisconnect( )
               runFlag = 1
            ENDIF

          ELSEIF( visitState < 0 )THEN

             EXIT

          ENDIF

         IF( MOD( iterate, dumpFreq ) == 0 )THEN
           !$OMP MASTER
           CALL myeu % WritePickup( iterate, myRank )
           !$OMP END MASTER
         ENDIF

      ENDDO
      !$OMP END PARALLEL
     

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !
#else
! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !
      !$OMP PARALLEL
      DO iT = iter0, iter0+nT-1, dFreq ! Loop over time-steps

         tn = real( iT, prec )*deltaT ! time at the beginning of this timestep          
         CALL myeu % ForwardStepRK3( tn, dFreq, myRank ) ! Forward Step
#ifdef HAVE_CUDA
         myeu % state % solution = myeu % state % solution_dev ! Update the host from the GPU
#endif


          ! In here, we will do the LIBSIM calls
! Turn off file I/O if timing is turned on
#ifndef TIMING
         !$OMP MASTER
         CALL myeu % WritePickup( iT+dFreq, myRank )
         CALL myeu % WriteTecplot( iT+dFreq, myeu % params % nPlot, myRank )
         !$OMP END MASTER
#endif

      ENDDO
      !$OMP END PARALLEL

#endif
     
  END SUBROUTINE MainLoop

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !
#ifdef INSITU_VIZ
  INTEGER FUNCTION VisitGetMetaData( handle )
     IMPLICIT NONE
     INTEGER :: handle
     
        VisitGetMetaData = VISIT_INVALID_HANDLE

  END FUNCTION VisitGetMetaData
#endif
END PROGRAM Fluid_Driver
