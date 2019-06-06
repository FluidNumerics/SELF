! SELF_Fluids_Driver.f90
!
! Copyright 2018 Joseph Schoonover <joe@myFluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


PROGRAM SELF_Fluids_Driver

  USE ModelPrecision
  USE ModelParameters_Class
  USE HexMesh_Class
  USE Fluid_Class

  IMPLICIT NONE

#include "self_macros.h"

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
#undef __FUNC__
#define __FUNC__ "Main"
  TYPE( Fluid )                :: myFluid
  LOGICAL                      :: setupSuccess
  LOGICAL                      :: meshgenSuccess
  LOGICAL                      :: initializeFromScratch
  LOGICAL                      :: pickupFileExists
  LOGICAL                      :: run_MeshGenerator
  LOGICAL                      :: run_Initializer
  LOGICAL                      :: run_Integrator
  CHARACTER(500)               :: equationFile, paramFile


! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !


  INFO('Start')

  CALL Initialize_MPILayer()

  CALL GetCLIConf( )

  IF( setupSuccess )THEN

    IF( run_MeshGenerator )THEN
      CALL MeshGen( )
    ENDIF

    IF( run_Initializer )THEN
      CALL Initialize( )
    ENDIF
  
    IF( run_Integrator )THEN
      CALL Integrate( )
    ENDIF


  ENDIF

  CALL Finalize_MPILayer( )

  INFO('End')

CONTAINS

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE GetCLIConf( )
    
    ! Local
    INTEGER :: nArg, argID
    CHARACTER(500) :: argName
    LOGICAL :: helpNeeded, equationFileProvided, paramFileProvided

    helpNeeded           = .FALSE.
    run_MeshGenerator    = .FALSE.
    run_Initializer      = .FALSE.
    run_Integrator       = .FALSE.
    equationFileProvided = .FALSE.
    paramFileProvided    = .FALSE.

    paramFile = './runtime.params'
    equationFile = './self.equations'

    nArg = command_argument_count( )

    setupSuccess = .TRUE.
    DO argID = 1, nArg

      CALL get_command_argument( argID, argName )

      SELECT CASE( TRIM( argName ) )

        CASE( "meshgen" )

          run_MeshGenerator = .TRUE.
          setupSuccess      = .TRUE.

        CASE( "initialize" )

          run_Initializer = .TRUE.
          setupSuccess    = .TRUE.

       CASE( "integrate" )
          run_Integrator    = .TRUE.
          setupSuccess      = .TRUE.

        CASE( "help" )
          helpNeeded   = .TRUE.
          setupSuccess = .FALSE.

        CASE( "--param-file" )
          paramFileProvided = .TRUE.

        CASE( "--equation-file" )
          equationFileProvided = .TRUE.

        CASE DEFAULT

          run_MeshGenerator = .TRUE.
          run_Initializer   = .TRUE.
          run_Integrator    = .TRUE.
          setupSuccess      = .TRUE.

          IF( paramFileProvided )THEN

            paramFile = TRIM( argName )
            paramFileProvided = .FALSE.

          ENDIF

          IF( equationFileProvided )THEN

            equationFile = TRIM( argName )
            equationFileProvided = .FALSE.

          ENDIF

      END SELECT

    ENDDO

    IF( helpNeeded ) THEN

      PRINT*, 'SELF-Fluids (sfluid) Command Line Interface'      
      PRINT*, ' '
      PRINT*, ' A program for solving Compressible Navier-Stokes using the'
      PRINT*, ' Nodal Discontinuous Galerkin Spectral Element Method.'
      PRINT*, ' '
      PRINT*, '  sfluid [tool] [options]'      
      PRINT*, ' '
      PRINT*, ' [tool] can be :'
      PRINT*, ' '
      PRINT*, '   help'
      PRINT*, '     Display this help message'
      PRINT*, ' '
      PRINT*, '   meshgen'
      PRINT*, '     Run only the mesh generator to generate a structured mesh.'
      PRINT*, '     The structured mesh is built using nXElem, nYElem, and nZElem'
      PRINT*, '     specified in runtime.params. Further, domain decomposition '
      PRINT*, '     for the structured mesh is done by setting the number of'
      PRINT*, '     of processes in each direction (nProcX, nProcY, nProcZ)'
      PRINT*, '     Topography can be set using an equation like'
      PRINT*, ' '
      PRINT*, '           h = exp( -(x-500.0)^2/200.0 )'
      PRINT*, ' '
      PRINT*, '     in the self.equations file. This will result in a terrain-following'
      PRINT*, '     structured mesh.'
      PRINT*, ' '
      PRINT*, '     Future releases of SELF-Fluids will offer more complete support'
      PRINT*, '     for working with typical unstructured mesh formats. '
      PRINT*, ' '
      PRINT*, '   initialize'
      PRINT*, '     Generate model initial conditions but do not forward step the'
      PRINT*, '     model. The initial conditions are read in from the self.equations file. '
      PRINT*, ' '
      PRINT*, '   integrate'
      PRINT*, '     Forward step the model using provided mesh file and initial conditions.'
      PRINT*, '     If no mesh file or initial conditions are provided, then they are generated. '
      PRINT*, ' '
      PRINT*, ' '
      PRINT*, '  [options] can be :'
      PRINT*, ' '
      PRINT*, '    --param-file /path/to/param/file'
      PRINT*, '       Specifies the full path to a file with namelist settings for'
      PRINT*, '       the sfluid application. If not provided, runtime.params in  '
      PRINT*, '       your current directory is assumed.                          '
      PRINT*, ' '
      PRINT*, '   --equation-file /path/to/equation/file'
      PRINT*, '       Specifies the full path to an equation file for setting the '
      PRINT*, '       initial conditions, topography shape (for structured mesh), '
      PRINT*, '       and the drag field. If not provided, self.equations in your '
      PRINT*, '       current directory is assumed.                               '
      PRINT*, ' '
      PRINT*, ' '

      setupSuccess = .FALSE.
      RETURN
    ENDIF



  END SUBROUTINE GetCLIConf

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE MeshGen( )

      IF( myRank == 0 )THEN
        CALL Generate_SELFMesh(paramFile, equationFile, nProc)
      ENDIF

  END SUBROUTINE MeshGen

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Initialize( )
#undef __FUNC__
#define __FUNC__ "Initialize"

    INFO('Start')
    CALL myFluid % Build( equationFile, paramFile, setupSuccess )
    CALL myFluid % SetInitialConditions( )
    CALL myFluid % IO( )
    CALL myFluid % Trash( )
    INFO('End')

  END SUBROUTINE Initialize

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Integrate( )
#undef __FUNC__
#define __FUNC__ "Integrate"
    INTEGER    :: iT
    CHARACTER(50) :: msg

    INFO('Start')

    CALL myFluid % Build( equationFile, paramFile, setupSuccess )
    CALL myFluid % Read_from_HDF5( pickupFileExists ) 

    IF( pickupFileExists )THEN
      DO iT = 1, myFluid % params % nDumps ! Loop over time-steps

        WRITE(msg,'(E15.5)')myFluid % simulationTime
        msg = 'Starting time loop at t='//msg//'s'
        INFO(TRIM(msg))

        CALL myFluid % ForwardStepRK3( myFluid % params % nStepsPerDump ) ! Forward Step
        CALL myFluid % IO( )

      ENDDO
    ELSE
      INFO('Pickup file not found')
      INFO('Integration not executed')
    ENDIF

    CALL myFluid % Trash( )

    INFO('End')


  END SUBROUTINE Integrate


END PROGRAM SELF_Fluids_Driver

