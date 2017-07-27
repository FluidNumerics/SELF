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


 TYPE( Fluid ) :: myeu
#ifdef TIMING
 TYPE( MultiTimers ) :: timers
#endif
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
      PRINT*, 'Fluid_Driver : Greetings from Process ', myRank, ' of ',nProcs
#else
      myRank = 0
      nProcs = 1
#endif
      CALL myeu % Build( myRank, nProcs )
      
      !$OMP PARALLEL
      CALL myeu % GlobalTimeDerivative( 0.0_prec, myRank )
      !$OMP END PARALLEL
      
#ifdef HAVE_CUDA
      myeu % stressTensor % solution = myeu % stressTensor % solution_dev
#endif
     ! CALL myeu % WriteStressTensorTecplot( 0, myeu % params % nPlot, myRank )
    ! ////////////////////////////////// Time Stepping /////////////////////////////////////////// !
      iter0  = myeu % params % iterInit
      nT     = myeu % params % nTimeSteps
      dFreq  = myeu % params % dumpFreq
      deltaT = myeu % params % dt
      nDumps = (nT)/dFreq
      
#ifdef HAVE_MPI
      CALL MPI_BARRIER( MPI_COMM_WORLD, mpiErr )
#endif


! If loop level parallelism is turned on, the parallel regions are placed within each routine called underneath
! GlobalTimeDerivative rather than outside. This feature is for demonstration purposes only.
#ifndef LOOP_LEVEL
      !$OMP PARALLEL
#endif

#ifdef TIMING
      !$OMP MASTER
      IF( myRank == 0 )THEN
         CALL timers % Build( )
         CALL timers % AddTimer( 'Main', 1 )
         CALL timers % StartTimer( 1 )
      ENDIF
      !$OMP END MASTER
#endif
      DO iT = iter0, iter0+nT-1, dFreq ! Loop over time-steps

         tn = real( iT, prec )*deltaT ! time at the beginning of this timestep          
         CALL myeu % ForwardStepRK3( tn, dFreq, myRank ) ! Forward Step
#ifdef HAVE_CUDA
         myeu % state % solution = myeu % state % solution_dev ! Update the host from the GPU
#endif

! Turn off file I/O if timing is turned on
#ifndef TIMING
         !$OMP MASTER
         CALL myeu % WritePickup( iT+dFreq, myRank )
         CALL myeu % WriteTecplot( iT+dFreq, myeu % params % nPlot, myRank )
         !$OMP END MASTER
#endif

      ENDDO
#ifdef TIMING
      !$OMP MASTER
      IF( myRank == 0 )THEN
         CALL timers % StopTimer( 1 )
         CALL timers % Write_MultiTimers( )
         CALL timers % Trash( )
      ENDIF
      !$OMP END MASTER
#endif
#ifndef LOOP_LEVEL
      !$OMP END PARALLEL
#endif


      CALL myeu % Trash( )
#ifdef HAVE_MPI      
      CALL MPI_FINALIZE( mpiErr )
#endif

END PROGRAM Fluid_Driver
