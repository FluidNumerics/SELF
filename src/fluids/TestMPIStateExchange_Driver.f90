! TestMPIStateExchange_Driver.f90
!
! Copyright 2017 Joseph Schoonover <schoonover.numerics@gmail.com>, SimonsTechnical
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


PROGRAM TestMPIStateExchange_Driver

! src/COMMON/
  USE ModelPrecision
  USE Timing
! src/highend/euler/mpi-cuda/
  USE CompressibleNavierStokesParams_CLASS
  USE CompressibleNavierStokes_CLASS

  IMPLICIT NONE


  TYPE( CompressibleNavierStokes ) :: myeu
  TYPE( MultiTimers ) :: timers
  INTEGER :: mpiErr, myRank, nProcs
  INTEGER :: iter0, nT, dFreq
  INTEGER :: iT, nDumps
  REAL(prec) :: tn, deltaT

#ifdef HAVE_MPI
  ! MPI Initialization
  CALL MPI_INIT( mpiErr )
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, myRank, mpiErr )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nProcs, mpiErr )
  ! Sanity check
  PRINT*, 'CompressibleNavierStokes_Driver : Greetings from Process ', myRank, ' of ',nProcs
#else
  myRank = 0
  nProcs = 1
#endif
  CALL myeu % Build( myRank, nProcs )

  ! ////////////////////////////////// Time Stepping /////////////////////////////////////////// !
  iter0  = myeu % params % iterInit
  nT     = myeu % params % nTimeSteps
  dFreq  = myeu % params % dumpFreq
  deltaT = myeu % params % dt
  nDumps = (nT)/dFreq

#ifdef HAVE_MPI
  CALL MPI_BARRIER( MPI_COMM_WORLD, mpiErr )
#endif

  CALL myeu % CalculateBoundarySolution( )      !
  CALL myeu % UpdateExternalState( 0.0_prec, myRank )

  IF( myRank == 0 )THEN
    CALL timers % Build( )
    CALL timers % AddTimer( 'MPI_StateExchange', 1 )
    CALL timers % StartTimer( 1 )
  ENDIF

  DO iT = 1, nT
    CALL myeu % MPI_StateExchange( myRank ) ! Forward Step
  ENDDO

  IF( myRank == 0 )THEN
    CALL timers % STOPTimer( 1 )
    CALL timers % Write_MultiTimers( )
    CALL timers % Trash( )

    PRINT*, myeu % externalState
  ENDIF



  CALL myeu % Trash( )
#ifdef HAVE_MPI
  CALL MPI_FINALIZE( mpiErr )
#endif

END PROGRAM TestMPIStateExchange_Driver
