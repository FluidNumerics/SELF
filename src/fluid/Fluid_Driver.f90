! Fluid_Driver.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
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

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> ! 

 TYPE( Fluid )       :: myeu
 INTEGER             :: mpiErr, myRank, nProcs
 LOGICAL             :: setupSuccess

#ifdef TIMING
 TYPE( MultiTimers ) :: timers
#endif

#ifdef DIAGNOSTICS
 INTEGER             :: diagUnits(1:nDiagnostics)
#endif

#ifdef INSITU_VIZ
 INTEGER             :: simulationCycle, simulationTime, runFlag
 COMMON/SIMSTATE/ runFlag, simulationCycle, simulationTime
#endif

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> ! 


      CALL Setup( )

      IF( setupSuccess )THEN
#ifdef INSITU_VIZ
         CALL SetupLibSim( )
#endif

         CALL MainLoop( )

         CALL Cleanup( )

      ENDIF

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
      CALL myeu % Build( myRank, nProcs, setupSuccess )
      IF( .NOT. setupSuccess )THEN
#ifdef HAVE_MPI
         CALL MPI_FINALIZE( mpiErr )
#endif
         RETURN
      ENDIF

#ifdef DIAGNOSTICS
      CALL myeu % OpenDiagnosticsFiles( diagUnits )
      CALL myeu % Diagnostics( )
      CALL myeu % WriteDiagnostics( diagUnits )
#endif


#ifdef TIMING
      IF( myRank == 0 )THEN
         CALL timers % Build( )
         CALL timers % AddTimer( 'ForwardStepRK3', 1 )
      ENDIF
#endif



  END SUBROUTINE Setup

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !
#ifdef INSITU_VIZ
  SUBROUTINE SetupLibSim( )
     IMPLICIT NONE
     INCLUDE "visitfortransimV2interface.inc"
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

#ifdef DIAGNOSTICS
      CALL myeu % CloseDiagnosticsFiles( diagUnits )
#endif

#ifdef TIMING
       IF( myRank == 0 )THEN
            CALL timers % Trash( )
       ENDIF
#endif

  END SUBROUTINE Cleanup

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE MainLoop( )
     IMPLICIT NONE
     INTEGER    :: iT
#ifdef INSITU_VIZ
     INCLUDE "visitfortransimV2interface.inc"
     INTEGER    :: visitState, visitResult, blocking, iterate     

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
     REAL(prec) :: density(0:myeu % params % nPlot, &
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
      simulationCycle = 0
      simulationTime  = myeu % simulationTime
      !$OMP PARALLEL
      DO ! Loop indefinitely; really breaks when visitState < 0

         IF( runFlag == 1 )THEN
            blocking = 0
         ELSE
            blocking = 1
         ENDIF
 
         visitState = visitDetectInput( blocking, -1 )

         IF( visitState == 0 )THEN

            CALL myeu % ForwardStepRK3( 1, myRank ) ! Forward Step
            simulationCycle = simulationCycle + 1
            simulationTime = myeu % simulationTime
#ifdef HAVE_CUDA
            myeu % state % solution = myeu % state % solution_dev ! Update the host from the GPU
#endif
            CALL myeu % FluidStateAtPlottingPoints( u, v, w, density, potentialTemp, pressure )

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

         IF( MOD( simulationCycle, dFreq ) == 0 )THEN
           !$OMP MASTER
           CALL myeu % WritePickup( myRank )
#ifdef DIAGNOSTICS
           CALL myeu % Diagnostics( ) 
           CALL myeu % WriteDiagnosticsFiles( diagUnits )
#endif
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
      DO iT = 1, myeu % params % nDumps ! Loop over time-steps


#ifdef TIMING
         !$OMP MASTER
         IF( myRank == 0 )THEN
            CALL timers % StartTimer( 1 )
         ENDIF
         !$OMP END MASTER
#endif


         CALL myeu % ForwardStepRK3( myeu % params % nStepsPerDump, myRank ) ! Forward Step

#ifdef TIMING
         !$OMP MASTER
         IF( myRank == 0 )THEN
            CALL timers % StopTimer( 1 )
         ENDIF
         !$OMP END MASTER
#endif


#ifdef HAVE_CUDA
         myeu % state % solution = myeu % state % solution_dev ! Update the host from the GPU
#endif


#ifdef DIAGNOSTICS
         CALL myeu % Diagnostics( ) 
         CALL myeu % WriteDiagnostics( diagUnits )
#endif
         !$OMP MASTER
         CALL myeu % WritePickup( myRank )
         CALL myeu % WriteTecplot( myRank )
         !$OMP END MASTER
#endif

      ENDDO
      !$OMP END PARALLEL

     
  END SUBROUTINE MainLoop


END PROGRAM Fluid_Driver

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !
#ifdef INSITU_VIZ
  SUBROUTINE VisitCommandCallBack (cmd, lcmd, intdata,   &
                                   floatdata, stringdata,&
                                   lstringdata )
     IMPLICIT NONE
     INCLUDE "visitfortransimV2interface.inc"
     CHARACTER(8) :: cmd, stringdata
     INTEGER      :: lcmd, lstringdata, intdata
     REAL         :: floatdata
     INTEGER      :: runflag, simcycle
     REAL         :: simtime

      IF( visitstrcmp(cmd, lcmd, "halt", 4) == 0)THEN
          runflag = 0
     ! ELSEIF(visitstrcmp(cmd, lcmd, "step", 4).eq.0) then
     !     CALL SimulateOneTimeStep( )
      ELSEIF(visitstrcmp(cmd, lcmd, "run", 3).eq.0) then
          runflag = 1
      ELSEIF(visitstrcmp(cmd, lcmd, "testcommand", 11).eq.0) then
          write (6,*) 'Received testcommand'
      endif

  END SUBROUTINE VisitCommandCallBack
!
  INTEGER FUNCTION VisitBroadCastIntFunction(value, sender)
     IMPLICIT NONE
     INTEGER :: value, sender
        
        VisitBroadCastIntFunction = 0

  END FUNCTION VisitBroadCastIntFunction
!
  INTEGER FUNCTION VisitBroadCastStringFunction(str, lstr, sender)
     IMPLICIT NONE
     CHARACTER(8) :: str
     INTEGER      :: lstr, sender

        VisitBroadCastStringFunction = 0

  END FUNCTION VisitBroadCastStringFunction
!
  SUBROUTINE VisitSlaveProcessCallBack( )
     IMPLICIT NONE
  END SUBROUTINE VisitSlaveProcessCallBack
!
  INTEGER FUNCTION VisitGetMetaData( )
     IMPLICIT NONE
     !INTEGER :: handle
     INCLUDE "visitfortransimV2interface.inc"
     INTEGER :: runFlag, simulationCycle, simulationTime
     COMMON /SIMSTATE/ runFlag, simulationCycle, simulationTime 
     INTEGER :: err, md, m1, cmd, vmd


        IF( VisitMDSimAlloc( md ) == VISIT_OKAY )THEN
           err = VisitMDSimSetCycleTime( md, simulationCycle, simulationTime )
           IF( runFlag == 1 )THEN
              err = VisitMDSimSetMode( md, VISIT_SIMMODE_RUNNING )
           ELSE
              err = VisitMDSimSetMode( md, VISIT_SIMMODE_STOPPED )
           ENDIF
        ENDIF

        VisitGetMetaData = md

        ! Set the Mesh Meta-Data
        IF( VisitMDMeshAlloc(m1) == VISIT_OKAY )THEN
           err = VisitMDMeshSetName( m1, "Mesh", 4 )
           err = VisitMDMeshSetMeshType( m1, VISIT_MESHTYPE_UNSTRUCTURED )
           err = VisitMDMeshSetTopologicalDim( m1, 3 )
           err = VisitMDMeshSetSpatialDim( m1, 3 )
           err = VisitMDMeshSetXUnits( m1, "m", 1 )
           err = VisitMDMeshSetYUnits( m1, "m", 1 )
           err = VisitMDMeshSetZUnits( m1, "m", 1 )
           err = VisitMDMeshSetXLabel( m1, "x", 1 )
           err = VisitMDMeshSetYLabel( m1, "y", 1 )
           err = VisitMDMeshSetZLabel( m1, "z", 1 )

           err = VisitMDSimAddMesh( md, m1 )
        ENDIF

        ! Add the density variable
        IF( VisitMDvarAlloc( vmd ) == VISIT_OKAY )THEN
           err = VisitMDvarSetName( vmd, "density", 7 )
           err = VisitMDvarSetMeshName( vmd, "Mesh", 4 )
           err = VisitMDvarSetType( vmd, VISIT_VARTYPE_SCALAR )
           err = VisitMDvarSetCentering( vmd, VISIT_VARCENTERING_NODE)
           err = VisitMDSimAddVariable( md, vmd )
        ENDIF

        ! Add meta-data for simulation control commands
        IF( VisitMDcmdAlloc( cmd ) == VISIT_OKAY )THEN
           err = VisitMDcmdSetName( cmd, "halt", 4 )
           err = VisitMDSimAddGenericCommand( md, cmd )
        ENDIF
        IF( VisitMDcmdAlloc( cmd ) == VISIT_OKAY )THEN
           err = VisitMDcmdSetName( cmd, "run", 3 )
           err = VisitMDSimAddGenericCommand( md, cmd )
        ENDIF
      !  IF( VisitMDcmdAlloc( cmd ) == VISIT_OKAY )THEN
      !     err = VisitMDcmdSetName( cmd, "step", 4 )
      !     err = VisitMDSimAddGenericCommand( md, cmd )
      !  ENDIF


       


  END FUNCTION VisitGetMetaData
!
  INTEGER FUNCTION VisitGetMesh(handle, domain, name, lname)
     IMPLICIT NONE
     CHARACTER(8) :: name
     INTEGER      :: handle, domain, lname
     INCLUDE "visitfortransimV2interface.inc"

       VisitGetMesh = VISIT_ERROR

  END FUNCTION VisitGetMesh
!
  INTEGER FUNCTION VisitGetMaterial(handle, domain, name, lname)
     IMPLICIT NONE
     CHARACTER(8) :: name
     INTEGER      :: handle, domain, lname
     INCLUDE "visitfortransimV2interface.inc"

       VisitGetMaterial = VISIT_ERROR

  END FUNCTION VisitGetMaterial
!
  INTEGER FUNCTION VisitGetDomainNesting(handle, domain, name, lname)
     IMPLICIT NONE
     CHARACTER(8) :: name
     INTEGER      :: handle, domain, lname
     INCLUDE "visitfortransimV2interface.inc"

       VisitGetDomainNesting = VISIT_ERROR

  END FUNCTION VisitGetDomainNesting
!
  INTEGER FUNCTION VisitGetDomainBounds(handle, domain, name, lname)
     IMPLICIT NONE
     CHARACTER(8) :: name
     INTEGER      :: handle, domain, lname
     INCLUDE "visitfortransimV2interface.inc"

       VisitGetDomainBounds = VISIT_ERROR

  END FUNCTION VisitGetDomainBounds
!
  INTEGER FUNCTION VisitGetDomainList(handle, domain, name, lname)
     IMPLICIT NONE
     CHARACTER(8) :: name
     INTEGER      :: handle, domain, lname
     INCLUDE "visitfortransimV2interface.inc"

       VisitGetDomainList = VISIT_ERROR

  END FUNCTION VisitGetDomainList
!
  INTEGER FUNCTION VisitGetCurve(handle, domain, name, lname)
     IMPLICIT NONE
     CHARACTER(8) :: name
     INTEGER      :: handle, domain, lname
     INCLUDE "visitfortransimV2interface.inc"

       VisitGetCurve = VISIT_ERROR

  END FUNCTION VisitGetCurve
!
  INTEGER FUNCTION VisitGetMixedVariable(handle, domain, name, lname)
     IMPLICIT NONE
     CHARACTER(8) :: name
     INTEGER      :: handle, domain, lname
     INCLUDE "visitfortransimV2interface.inc"

       VisitGetMixedVariable = VISIT_ERROR

  END FUNCTION VisitGetMixedVariable
!
  INTEGER FUNCTION VisitGetVariable(handle, domain, name, lname)
     IMPLICIT NONE
     CHARACTER(8) :: name
     INTEGER      :: handle, domain, lname
     INCLUDE "visitfortransimV2interface.inc"

       VisitGetVariable = VISIT_ERROR

  END FUNCTION VisitGetVariable
!
  INTEGER FUNCTION VisitActivateTimeStep(handle, domain, name, lname)
     IMPLICIT NONE
     CHARACTER(8) :: name
     INTEGER      ::  handle, domain, lname
     INCLUDE "visitfortransimV2interface.inc"

       VisitActivateTimeStep = VISIT_ERROR

  END FUNCTION VisitActivateTimeStep
!

#endif

!END PROGRAM Fluid_Driver
