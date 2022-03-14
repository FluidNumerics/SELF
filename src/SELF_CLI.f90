MODULE SELF_CLI
!
! Copyright 2021-2022 Fluid Numerics LLC
!
! ============================================================
!
!! Loads options for a self model from a json file
!! and creates a FLAP command line interface object.
!!
!! The provided json file must meet the following schema :
!!
!!   {
!!     "self_model": {
!!       "name": "",
!!       "version": "",
!!       "description":"",
!!       "license": "",
!!       "authors": "",
!!       "options": [
!!          {
!!           "type": "real | integer | string | logical",
!!           "value": "",
!!           "cli_long": "--this",
!!           "cli_short": "-t",
!!           "description": "",
!!   	    "action": "",
!!   	    "required": false | true,
!!           "choices": ""
!!          }
!!       ]
!!     }
!!   }
!!
!!

  USE SELF_Constants
  ! External Modules
  USE json_module
  USE FLAP
  USE ISO_FORTRAN_ENV

  TYPE,PUBLIC :: CLI
    TYPE(JSON_FILE) :: json
    TYPE(COMMAND_LINE_INTERFACE) :: cliObj
    CHARACTER(SELF_FILE_DEFAULT_LENGTH) :: config

    CONTAINS

      PROCEDURE, PUBLIC :: Init => Init_CLI
      PROCEDURE, PUBLIC :: Free => Free_CLI
      PROCEDURE, PRIVATE :: SetOptions_CLI
      PROCEDURE, PUBLIC :: LoadFromCLI

      GENERIC, PUBLIC :: Get_CLI => Get_CLI_int32, &
                                    Get_CLI_int64, &
                                    Get_CLI_real32, &
                                    Get_CLI_real64, &
                                    Get_CLI_logical

      PROCEDURE, PRIVATE :: Get_CLI_int32
      PROCEDURE, PRIVATE :: Get_CLI_int64
      PROCEDURE, PRIVATE :: Get_CLI_real32
      PROCEDURE, PRIVATE :: Get_CLI_real64
      PROCEDURE, PRIVATE :: Get_CLI_logical

  END TYPE CLI 

  INTEGER, PARAMETER :: SELF_JSON_DEFAULT_KEY_LENGTH = 200
  INTEGER, PARAMETER :: SELF_JSON_DEFAULT_VALUE_LENGTH = 200

  CONTAINS

  SUBROUTINE Init_CLI( this, configFile )
    !! Initialize an instance of the CLI class using a provided configuration file
    !! An example configuration file can be found in `self/etc/config_example.json`
    !! The initialization procedure initilizes a json_file class (from json-fortran)
    !! and a COMMAND_LINE_INTERFACE class (from FLAP). Further, the CLI options are
    !! set according to the self_model.options field in the  configuration file.
    IMPLICIT NONE
    CLASS(CLI), INTENT(out) :: this
    CHARACTER(*), INTENT(in) :: configFile
    !
    CHARACTER(LEN=:),ALLOCATABLE :: progname
    CHARACTER(LEN=:),ALLOCATABLE :: version
    CHARACTER(LEN=:),ALLOCATABLE :: description
    CHARACTER(LEN=:),ALLOCATABLE :: license
    CHARACTER(LEN=:),ALLOCATABLE :: authors
    LOGICAL :: found

      this % config = configFile

      CALL this % json % initialize(stop_on_error = .true., &
                                    comment_char = '#')

      CALL this % json % load_file(filename = TRIM(this % config))

!      CALL this % json % print_file()

      CALL this % json % get('self_model.name',progname,found)
      CALL this % json % get('self_model.version',version,found)
      CALL this % json % get('self_model.description',description,found)
      CALL this % json % get('self_model.license',license,found)
      CALL this % json % get('self_model.authors',authors,found)

      CALL this % cliObj % init(progname=TRIM(progname),&
                             version=TRIM(version),&
                             description=TRIM(description),&
                             license=TRIM(license),&
                             authors=TRIM(authors))

      CALL this % SetOptions_CLI()

      IF( ALLOCATED(progname) ) DEALLOCATE(progname)
      IF( ALLOCATED(version) ) DEALLOCATE(version)
      IF( ALLOCATED(description) ) DEALLOCATE(description)
      IF( ALLOCATED(license) ) DEALLOCATE(license)
      IF( ALLOCATED(authors) ) DEALLOCATE(authors)
       
  END SUBROUTINE Init_CLI

  SUBROUTINE Free_CLI( this )
    !! Frees the attributes of the CLI class and reset the config attribute
    !! to an empty string
    IMPLICIT NONE
    CLASS(CLI), INTENT(inout) :: this

      CALL this % json % destroy()
      this % config = ""

  END SUBROUTINE Free_CLI

  SUBROUTINE SetOptions_CLI( this )
    IMPLICIT NONE
    CLASS(CLI), INTENT(inout) :: this
    ! Local
    INTEGER :: nopts
    INTEGER :: i
    CHARACTER(4) :: arrayCount
    CHARACTER(SELF_JSON_DEFAULT_KEY_LENGTH) :: jsonKey
    CHARACTER(LEN=:),ALLOCATABLE :: cliObjLong
    CHARACTER(LEN=:),ALLOCATABLE :: cliObjShort
    CHARACTER(LEN=:),ALLOCATABLE :: cliObjDescription
    CHARACTER(LEN=:),ALLOCATABLE :: cliObjDefault
    CHARACTER(LEN=:),ALLOCATABLE :: cliObjType
    CHARACTER(LEN=:),ALLOCATABLE :: cliObjAction
    CHARACTER(LEN=:),ALLOCATABLE :: cliObjChoices
    LOGICAL :: cliObjRequired
    LOGICAL :: found

      CALL this % json % info('self_model.options',n_children=nopts)

      DO i = 1, nopts
        WRITE(arrayCount,"(I4)") i

        jsonKey = "self_model.options["//&
                  TRIM(arrayCount)//&
                  "].cli_long" 
        CALL this % json % get( TRIM(jsonKey), cliObjLong, found )

        jsonKey = "self_model.options["//&
                  TRIM(arrayCount)//&
                  "].cli_short" 
        CALL this % json % get( TRIM(jsonKey), cliObjShort, found )

        jsonKey = "self_model.options["//&
                  TRIM(arrayCount)//&
                  "].description" 
        CALL this % json % get( TRIM(jsonKey), cliObjDescription, found )

        jsonKey = "self_model.options["//&
                  TRIM(arrayCount)//&
                  "].value" 
        CALL this % json % get( TRIM(jsonKey), cliObjDefault, found )

        jsonKey = "self_model.options["//&
                  TRIM(arrayCount)//&
                  "].action" 
        CALL this % json % get( TRIM(jsonKey), cliObjAction, found )

        jsonKey = "self_model.options["//&
                  TRIM(arrayCount)//&
                  "].choices" 
        CALL this % json % get( TRIM(jsonKey), cliObjChoices, found )

        jsonKey = "self_model.options["//&
                  TRIM(arrayCount)//&
                  "].required" 
        CALL this % json % get( TRIM(jsonKey), cliObjRequired, found )

        IF( TRIM(cliObjAction) == "" )THEN

          IF( TRIM(cliObjChoices) == "" )THEN

            CALL this % cliObj % add( switch = TRIM(cliObjLong), &
                                   switch_ab = TRIM(cliObjShort), &
                                   help = TRIM(cliObjDescription)//NEW_LINE("A"), &
                                   def = TRIM(cliObjDefault), &
                                   required = cliObjRequired )

          ELSE

            CALL this % cliObj % add( switch = TRIM(cliObjLong), &
                                   switch_ab = TRIM(cliObjShort), &
                                   help = TRIM(cliObjDescription)//NEW_LINE("A"), &
                                   def = TRIM(cliObjDefault), &
                                   choices = TRIM(cliObjChoices), &
                                   required = cliObjRequired )

          ENDIF

        ELSE

          CALL this % cliObj % add( switch = TRIM(cliObjLong), &
                                 help = TRIM(cliObjDescription)//NEW_LINE("A"), &
                                 act = TRIM(cliObjAction), &
                                 def = TRIM(cliObjDefault), &
                                 required = cliObjRequired )

        ENDIF

        IF( ALLOCATED(cliObjLong) ) DEALLOCATE(cliObjLong)
        IF( ALLOCATED(cliObjShort) ) DEALLOCATE(cliObjShort)
        IF( ALLOCATED(cliObjDescription) ) DEALLOCATE(cliObjDescription)
        IF( ALLOCATED(cliObjDefault) ) DEALLOCATE(cliObjDefault)
        IF( ALLOCATED(cliObjType) ) DEALLOCATE(cliObjType)
        IF( ALLOCATED(cliObjAction) ) DEALLOCATE(cliObjAction)
        IF( ALLOCATED(cliObjChoices) ) DEALLOCATE(cliObjChoices)

      ENDDO

  END SUBROUTINE SetOptions_CLI

  SUBROUTINE LoadFromCLI( this )
    IMPLICIT NONE
    CLASS(CLI), INTENT(inout) :: this
    ! Local
    INTEGER :: nopts
    INTEGER :: i
    CHARACTER(4) :: arrayCount
    CHARACTER(SELF_JSON_DEFAULT_KEY_LENGTH) :: jsonKey
    CHARACTER(SELF_JSON_DEFAULT_VALUE_LENGTH) :: tmpVal
    CHARACTER(LEN=:),ALLOCATABLE :: cliObjLong
    LOGICAL :: found

      CALL this % json % info('self_model.options',n_children=nopts)

      DO i = 1, nopts
      
        ! Get the cli_long option
        WRITE(arrayCount,"(I4)") i
        jsonKey = "self_model.options["//&
                  TRIM(arrayCount)//&
                  "].cli_long" 
        CALL this % json % get( TRIM(jsonKey), cliObjLong, found )

        CALL this % cliObj % get(val=tmpVal,switch=TRIM(cliObjLong))

        ! Set the value for the option in the json object
        jsonKey = "self_model.options["//&
                  TRIM(arrayCount)//&
                  "].value" 

        CALL this % json % update( TRIM(jsonKey), TRIM(tmpVal), found )

        IF( ALLOCATED(cliObjLong) ) DEALLOCATE(cliObjLong)

      ENDDO

  END SUBROUTINE LoadFromCLI

  SUBROUTINE Get_CLI_int32( this, option, res )
    IMPLICIT NONE
    CLASS(CLI), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: option
    INTEGER(int32), INTENT(out) :: res 
    ! Local
    INTEGER :: error

      CALL this % cliObj % get(val=res,switch=TRIM(option),error=error)
      IF( error /= 0 )THEN
        PRINT*, "Configuration key not found : "//TRIM(option)
        STOP 1
      ENDIF

  END SUBROUTINE Get_CLI_int32

  SUBROUTINE Get_CLI_int64( this, option, res )
    IMPLICIT NONE
    CLASS(CLI), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: option
    INTEGER(int64), INTENT(out) :: res 
    ! Local
    INTEGER :: error

      CALL this % cliObj % get(val=res,switch=TRIM(option),error=error)
      IF( error /= 0 )THEN
        PRINT*, "Configuration key not found : "//TRIM(option)
        STOP 1
      ENDIF

  END SUBROUTINE Get_CLI_int64

  SUBROUTINE Get_CLI_real32( this, option, res )
    IMPLICIT NONE
    CLASS(CLI), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: option
    REAL(real32), INTENT(out) :: res 
    ! Local
    INTEGER :: error

      CALL this % cliObj % get(val=res,switch=TRIM(option),error=error)
      IF( error /= 0 )THEN
        PRINT*, "Configuration key not found : "//TRIM(option)
        STOP 1
      ENDIF

  END SUBROUTINE Get_CLI_real32

  SUBROUTINE Get_CLI_real64( this, option, res )
    IMPLICIT NONE
    CLASS(CLI), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: option
    REAL(real64), INTENT(out) :: res 
    ! Local
    INTEGER :: error

      CALL this % cliObj % get(val=res,switch=TRIM(option),error=error)
      IF( error /= 0 )THEN
        PRINT*, "Configuration key not found : "//TRIM(option)
        STOP 1
      ENDIF

  END SUBROUTINE Get_CLI_real64

  SUBROUTINE Get_CLI_logical( this, option, res )
    IMPLICIT NONE
    CLASS(CLI), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: option
    LOGICAL, INTENT(out) :: res 
    ! Local
    INTEGER :: error

      CALL this % cliObj % get(val=res,switch=TRIM(option),error=error)
      IF( error /= 0 )THEN
        PRINT*, "Configuration key not found : "//TRIM(option)
        STOP 1
      ENDIF

  END SUBROUTINE Get_CLI_logical


END MODULE SELF_CLI
