! Fluid_PickupToTecplot.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

 
PROGRAM Fluid_PickupToTecplot

! src/common/
USE ModelPrecision
! src/highend/euler/
USE FluidParams_Class
USE Fluid_Class

 IMPLICIT NONE

 TYPE( Fluid ) :: myeu
 INTEGER :: mpiErr, myRank, nProcs
 INTEGER :: iter0, nT, dFreq
 INTEGER :: iT

#ifdef HAVE_MPI
      ! MPI Initialization
      CALL MPI_INIT( mpiErr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myRank, mpiErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nProcs, mpiErr )
      ! Sanity check
      PRINT*, 'Fluid_PickupToTecplot : Greetings from Process ', myRank, ' of ',nProcs
#else
      myRank = 0
      nProcs = 1
#endif
      CALL myeu % Build( myRank, nProcs )
 
    ! ////////////////////////////////// Time Stepping /////////////////////////////////////////// !
      iter0  = myeu % params % iterInit
      nT     = myeu % params % nTimeSteps
      dFreq  = myeu % params % dumpFreq
      
      DO iT = iter0, iter0+nT-1, dFreq ! Loop over time-steps

         CALL myeu % ReadPickup( iT+dFreq, myRank )
         CALL myeu % WriteTecplot( iT+dFreq, myeu % params % nPlot, myRank )

      ENDDO
       
      CALL myeu % Trash( )
#ifdef HAVE_MPI
      CALL MPI_FINALIZE( mpiErr )
#endif

END PROGRAM Fluid_PickupToTecplot
