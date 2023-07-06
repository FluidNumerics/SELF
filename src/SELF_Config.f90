MODULE SELF_Config
   !
   ! Copyright 2021-2022 Fluid Numerics LLC
   !
   ! ============================================================
   !
   !! A class to load a schema for model configuration, and can be used
   !! to load concretized configurations. The base schema for all models
   !! is defined in ${SELF_PREFIX}/etc/schema/self.json
   !!
   !!
   !!

   !USE SELF_Constants
   USE SELF_CLI
   ! External Modules
   USE json_module
   USE ISO_FORTRAN_ENV

   IMPLICIT NONE

#include "SELF_Macros.h"

   INTEGER, PARAMETER :: SELF_FILE_DEFAULT_LENGTH=500

   TYPE,PUBLIC :: SELFConfig
      !TYPE(JSON_FILE) :: schema
      TYPE(JSON_FILE) :: concretization
      !CHARACTER(SELF_FILE_DEFAULT_LENGTH) :: schemaFile

   CONTAINS

      GENERIC, PUBLIC :: Init => Init_SELFConfig_FromFile, Init_SELFConfig_FromCLI
      PROCEDURE, PRIVATE :: Init_SELFConfig_FromFile
      PROCEDURE, PRIVATE :: Init_SELFConfig_FromCLI

      !GENERIC, PUBLIC :: LoadSchema => LoadSchema_SELFConfig_FromFile
      !PROCEDURE, PRIVATE :: LoadSchema_SELFConfig_FromFile

      GENERIC, PUBLIC :: LoadConcretization => LoadConcretization_SELFConfig_FromFile
      PROCEDURE, PRIVATE :: LoadConcretization_SELFConfig_FromFile

      PROCEDURE, PUBLIC :: Free => Free_SELFConfig

      GENERIC, PUBLIC :: Get => Get_SELFConfig_int32, &
         Get_SELFConfig_int64, &
         Get_SELFConfig_real32, &
         Get_SELFConfig_real64, &
         Get_SELFConfig_logical, &
         Get_SELFConfig_char

      PROCEDURE, PRIVATE :: Get_SELFConfig_int32
      PROCEDURE, PRIVATE :: Get_SELFConfig_int64
      PROCEDURE, PRIVATE :: Get_SELFConfig_real32
      PROCEDURE, PRIVATE :: Get_SELFConfig_real64
      PROCEDURE, PRIVATE :: Get_SELFConfig_logical
      PROCEDURE, PRIVATE :: Get_SELFConfig_char
      

   END TYPE SELFConfig

   INTEGER, PARAMETER :: SELF_JSON_DEFAULT_KEY_LENGTH = 200
   INTEGER, PARAMETER :: SELF_JSON_DEFAULT_VALUE_LENGTH = 200
   

CONTAINS

   SUBROUTINE Init_SELFConfig_FromFile( this, concretizationFile )
      IMPLICIT NONE
      CLASS(SELFConfig), INTENT(out) :: this
      CHARACTER(*), INTENT(in) :: concretizationFile

      !CALL this % LoadSchema( schemaFile )
      CALL this % LoadConcretization( concretizationFile )

   END SUBROUTINE Init_SELFConfig_FromFile

   SUBROUTINE Init_SELFConfig_FromCLI( this )
#undef __FUNC__
#define __FUNC__ "Init"
      IMPLICIT NONE
      CLASS(SELFConfig), INTENT(out) :: this
      ! Local
      CHARACTER(LEN=SELF_FILE_DEFAULT_LENGTH) :: concretizationFile 
      CHARACTER(LEN=200) :: SELF_PREFIX
      LOGICAL :: fileExists
      
      ! Set Default configuration file to
      ! ${SELF_PREFIX}/etc/schema/defaults/self.json
      CALL get_environment_variable("SELF_PREFIX", SELF_PREFIX)
      concretizationFile = TRIM(SELF_PREFIX)//"/etc/schema/defaults/self.json"

      IF ( CommandLineArgumentIsPresent(argument = "-i") )     THEN
         concretizationFile = StringValueForArgument(argument = "-i")
      END IF
      INFO("Using configuration file : "//TRIM(concretizationFile))
      INQUIRE(FILE=TRIM(concretizationFile), EXIST=fileExists )
      IF( fileExists )THEN
        CALL this % LoadConcretization( TRIM(concretizationFile) )
      ELSE
         ERROR("Configuration file does not exist : "//TRIM(concretizationFile))
         STOP 1
      ENDIF

   END SUBROUTINE Init_SELFConfig_FromCLI

   ! SUBROUTINE LoadSchema_SELFConfig_FromFile( this, schemaFile )
   !    !! Loads schema from file and stores in schema attribute
   !    IMPLICIT NONE
   !    CLASS(SELFConfig), INTENT(out) :: this
   !    CHARACTER(*), INTENT(in) :: schemaFile

   !    this % schemaFile = schemaFile
   !    CALL this % schema % initialize(stop_on_error = .true., &
   !       comment_char = '#')

   !    CALL this % schema % load_file(filename = TRIM(schemaFile))

   !    CALL this % schema % print_file()

   ! END SUBROUTINE LoadSchema_SELFConfig_FromFile

   SUBROUTINE LoadConcretization_SELFConfig_FromFile( this, concretizationFile )
      !! Loads a concretization and stores in concretization attributes
      IMPLICIT NONE
      CLASS(SELFConfig), INTENT(out) :: this
      CHARACTER(*), INTENT(in) :: concretizationFile

      CALL this % concretization % initialize(stop_on_error = .true., &
      comment_char = '#')

      CALL this % concretization % load_file(filename = TRIM(concretizationFile))

      CALL this % concretization % print_file()

   END SUBROUTINE LoadConcretization_SELFConfig_FromFile

   SUBROUTINE Free_SELFConfig( this )
      !! Frees the attributes of the SELFConfig class and reset the config attribute
      !! to an empty string
      IMPLICIT NONE
      CLASS(SELFConfig), INTENT(inout) :: this

      !CALL this % schema % destroy()
      CALL this % concretization % destroy()
     ! this % schemaFile = ""

   END SUBROUTINE Free_SELFConfig

   SUBROUTINE Get_SELFConfig_int32( this, jsonKey, res )
      IMPLICIT NONE
      CLASS(SELFConfig), INTENT(inout) :: this
      CHARACTER(*), INTENT(in) :: jsonKey
      INTEGER(int32), INTENT(out) :: res
      ! Local
      LOGICAL :: found

      CALL this % concretization % get( TRIM(jsonKey), res, found )
      IF( .NOT. found )THEN
         PRINT*, "JSON key not found : "//TRIM(jsonKey)
         STOP 1
      ENDIF

   END SUBROUTINE Get_SELFConfig_int32

   SUBROUTINE Get_SELFConfig_int64( this, jsonKey, res )
      IMPLICIT NONE
      CLASS(SELFConfig), INTENT(inout) :: this
      CHARACTER(*), INTENT(in) :: jsonKey
      INTEGER(int64), INTENT(out) :: res
      ! Local
      LOGICAL :: found
      INTEGER(int32) :: res32

      CALL this % concretization % get( TRIM(jsonKey), res32, found )
      IF( .NOT. found )THEN
         PRINT*, "JSON key not found : "//TRIM(jsonKey)
         STOP 1
      ENDIF
      res = INT( res32, kind=int64 )

   END SUBROUTINE Get_SELFConfig_int64

   SUBROUTINE Get_SELFConfig_real32( this, jsonKey, res )
      IMPLICIT NONE
      CLASS(SELFConfig), INTENT(inout) :: this
      CHARACTER(*), INTENT(in) :: jsonKey
      REAL(real32), INTENT(out) :: res
      ! Local
      LOGICAL :: found

      CALL this % concretization % get( TRIM(jsonKey), res, found )
      IF( .NOT. found )THEN
         PRINT*, "JSON key not found : "//TRIM(jsonKey)
         STOP 1
      ENDIF

   END SUBROUTINE Get_SELFConfig_real32

   SUBROUTINE Get_SELFConfig_real64( this, jsonKey, res )
      IMPLICIT NONE
      CLASS(SELFConfig), INTENT(inout) :: this
      CHARACTER(*), INTENT(in) :: jsonKey
      REAL(real64), INTENT(out) :: res
      ! Local
      LOGICAL :: found
      
      CALL this % concretization % get( TRIM(jsonKey), res, found )
      IF( .NOT. found )THEN
         PRINT*, "JSON key not found : "//TRIM(jsonKey)
         STOP 1
      ENDIF

   END SUBROUTINE Get_SELFConfig_real64

   SUBROUTINE Get_SELFConfig_logical( this, jsonKey, res )
      IMPLICIT NONE
      CLASS(SELFConfig), INTENT(inout) :: this
      CHARACTER(*), INTENT(in) :: jsonKey
      LOGICAL, INTENT(out) :: res
      ! Local
      LOGICAL :: found

      CALL this % concretization % get( TRIM(jsonKey), res, found )
      IF( .NOT. found )THEN
         PRINT*, "JSON key not found : "//TRIM(jsonKey)
         STOP 1
      ENDIF

   END SUBROUTINE Get_SELFConfig_logical

   SUBROUTINE Get_SELFConfig_char( this, jsonKey, res )
      IMPLICIT NONE
      CLASS(SELFConfig), INTENT(inout) :: this
      CHARACTER(*), INTENT(in) :: jsonKey
      CHARACTER(*), INTENT(out) :: res
      ! Local
      LOGICAL :: found
      CHARACTER(LEN=:),ALLOCATABLE :: resLoc

      CALL this % concretization % get( TRIM(jsonKey), resLoc, found )
      IF( .NOT. found )THEN
         PRINT*, "JSON key not found : "//TRIM(jsonKey)
         STOP 1
      ENDIF
      res = TRIM(resLoc)

   END SUBROUTINE Get_SELFConfig_char


END MODULE SELF_Config
