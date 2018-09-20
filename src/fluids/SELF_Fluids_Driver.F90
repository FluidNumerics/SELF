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

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !

  TYPE( Fluid )                :: myFluid
  LOGICAL                      :: setupSuccess
  LOGICAL                      :: meshgenSuccess
  LOGICAL                      :: initializeFromScratch
  LOGICAL                      :: pickupFileExists
  LOGICAL                      :: run_MeshGenOnly, run_UpToInitOnly
  CHARACTER(500)               :: equationFile, paramFile


! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !


  CALL Setup( )

  IF( setupSuccess )THEN

    IF( run_MeshGenOnly )THEN

      CALL MeshGen( )
    
    ELSE

      CALL Initialize( )

      IF( .NOT. run_UpToInitOnly )THEN
        CALL MainLoop( )
      ENDIF

      CALL Cleanup( )

    ENDIF


  ENDIF

CONTAINS

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Setup( )
    IMPLICIT NONE
    ! Local
    INTEGER :: nArg, argID
    CHARACTER(500) :: argName
    LOGICAL :: helpNeeded, equationFileProvided, paramFileProvided

    helpNeeded           = .FALSE.
    run_MeshGenOnly      = .FALSE.
    run_UpToInitOnly     = .FALSE.
    equationFileProvided = .FALSE.
    paramFileProvided    = .FALSE.

    paramFile = './runtime.params'
    equationFile = './self.equations'

    nArg = command_argument_count( )

    DO argID = 1, nArg

      CALL get_command_argument( argID, argName )

      SELECT CASE( TRIM( argName ) )

        CASE( "meshgen" )

          run_MeshGenOnly  = .TRUE.
          run_UpToInitOnly = .FALSE.
          setupSuccess     = .TRUE.

        CASE( "init" )

          run_MeshGenOnly  = .FALSE.
          run_UpToInitOnly = .TRUE.
          setupSuccess     = .TRUE.

        CASE( "help" )
          helpNeeded   = .TRUE.
          setupSuccess = .FALSE.

        CASE( "--param-file" )
          paramFileProvided = .TRUE.

        CASE( "--equation-file" )
          equationFileProvided = .TRUE.

        CASE DEFAULT

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

    PRINT*, TRIM(equationFile)

    IF( helpNeeded ) THEN

      PRINT*, 'SELF-Fluids (sfluid) Command Line Tool'      
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
      PRINT*, '   init'
      PRINT*, '     Run up to the initial condition generation and do not forward'
      PRINT*, '     step the model. The initial conditions are read in from the '
      PRINT*, '     self.equations file. '
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

    IF( .NOT. run_MeshGenOnly )THEN
      CALL myFluid % Build( equationFile, paramFile, setupSuccess )
    ENDIF


  END SUBROUTINE Setup

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE MeshGen( )


      CALL myFluid % ExtComm % SetRanks( )

      IF( myFluid % ExtComm % myRank == 0 )THEN
        PRINT*, '  Generating structured mesh...'
        CALL StructuredMeshGenerator_3D( paramFile, meshGenSuccess )
        PRINT*, '  Done'
      ENDIF

      CALL myFluid % ExtComm % Finalize( )


  END SUBROUTINE MeshGen

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Initialize( )

    ! Attempt to read the fluid pickup file. If it doesn't exist, this routine
    ! returns FALSE.
    CALL myFluid % Read_from_HDF5( pickupFileExists ) 

    ! If the pickup file doesn't exist, then the initial conditions are generated
    ! from the equation parser.
    IF( .NOT. pickupFileExists )THEN

      PRINT(MsgFMT), 'Pickup file not found.'

    ENDIF

    IF( .NOT. pickupFileExists .OR. run_UpToInitOnly )THEN

      PRINT(MsgFMT), 'Attempting initial condition generation from self.equations'
      CALL myFluid % SetInitialConditions( )
      CALL myFluid % IO( )

    ENDIF

  END SUBROUTINE Initialize

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

      CALL myFluid % IO( )

    ENDDO


  END SUBROUTINE MainLoop


END PROGRAM SELF_Fluids_Driver

