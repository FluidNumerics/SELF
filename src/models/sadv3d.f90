PROGRAM sadv3d
        
USE SELF_Advection3D
USE SELF_Mesh
USE FEQParse
USE SELF_Constants
USE SELF_CLI

  TYPE( Advection3D ) :: model
  TYPE( CLI ) :: cliConf
  LOGICAL :: convergenceCheck
  REAL(prec) :: maxError

  CALL InitCLI( )

  ! Now that we have CLI variables,
  ! we need to determine what to do
  !
  ! When convergence check is requested, we want to advance the model multiple
  ! times, gradually increasing N. Additionally, we want to calculate the
  ! difference between the exact solution and the numerical solution to obtain
  ! the max(abs(error)).

  CALL cliConf % cliObj % get(val=convergenceCheck,switch='--convergence-check')

  IF( convergenceCheck ) THEN

    ! Do the initial run with N = --control-degree
    CALL model % InitWithCLI(cliConf)
    CALL model % ModelExecute( io=.FALSE. ) 
    CALL model % MaxSolutionError(maxError) 
    CALL model % Free()


  ELSE

    CALL model % InitWithCLI(cliConf)
    CALL model % ModelExecute() 
    CALL model % Free()

  ENDIF

CONTAINS

  SUBROUTINE InitCLI( ) 
    ! Local
    LOGICAL :: fileExists
    CHARACTER(LEN=self_FileNameLength) :: selfInstallDir

      CALL GET_ENVIRONMENT_VARIABLE('SELF_INSTALL_DIR', selfInstallDir)
      IF( TRIM(selfInstallDir) == "" )THEN
        PRINT*, "SELF_INSTALL_DIR environment variable unset."
        PRINT*, "Trying /opt/view/"
        selfInstallDir = "/opt/view"
      ENDIF

      INQUIRE(FILE=TRIM(selfInstallDir)//'/etc/sadv3d.json', EXIST=fileExists)
      IF( .NOT. fileExists )THEN
        PRINT*, TRIM(selfInstallDir)//'/etc/sadv3d.json'//' configuration file not found!'
        STOP "ERROR"
      ENDIF

      CALL cliConf % Init(TRIM(selfInstallDir)//'/etc/sadv3d.json')
      CALL cliConf % LoadFromCLI() 

  END SUBROUTINE InitCLI


END PROGRAM sadv3d
