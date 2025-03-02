! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2024 Fluid Numerics LLC
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
!    the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUsLESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARIsLG IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_JSON_Config

  !USE SELF_Constants
  !USE SELF_CLI
  ! External Modules
  use json_module
  use iso_fortran_env

  implicit none

  integer,parameter :: SELF_FILE_DEFAULT_LENGTH = 500

  type,public :: SELFConfig
    !TYPE(JSON_FILE) :: schema
    type(JSON_FILE) :: concretization
    !CHARACTER(SELF_FILE_DEFAULT_LENGTH) :: schemaFile

  contains

    generic,public :: Init => Init_SELFConfig_FromFile !, Init_SELFConfig_FromCLI
    procedure,private :: Init_SELFConfig_FromFile
    !PROCEDURE, PRIVATE :: Init_SELFConfig_FromCLI

    !GENERIC, PUBLIC :: LoadSchema => LoadSchema_SELFConfig_FromFile
    !PROCEDURE, PRIVATE :: LoadSchema_SELFConfig_FromFile

    generic,public :: LoadConcretization => LoadConcretization_SELFConfig_FromFile
    procedure,private :: LoadConcretization_SELFConfig_FromFile

    procedure,public :: Free => Free_SELFConfig

    generic,public :: Get => Get_SELFConfig_int32, &
      Get_SELFConfig_int64, &
      Get_SELFConfig_real32, &
      Get_SELFConfig_real64, &
      Get_SELFConfig_logical, &
      Get_SELFConfig_char

    procedure,private :: Get_SELFConfig_int32
    procedure,private :: Get_SELFConfig_int64
    procedure,private :: Get_SELFConfig_real32
    procedure,private :: Get_SELFConfig_real64
    procedure,private :: Get_SELFConfig_logical
    procedure,private :: Get_SELFConfig_char

  endtype SELFConfig

  integer,parameter :: SELF_JSON_DEFAULT_KEY_LENGTH = 200
  integer,parameter :: SELF_JSON_DEFAULT_VALUE_LENGTH = 200

contains

  subroutine Init_SELFConfig_FromFile(this,concretizationFile)
    implicit none
    class(SELFConfig),intent(out) :: this
    character(*),intent(in) :: concretizationFile

    !CALL this % LoadSchema( schemaFile )
    call this%LoadConcretization(concretizationFile)

  endsubroutine Init_SELFConfig_FromFile

!     SUBROUTINE Init_SELFConfig_FromCLI( this )
!  #undef __FUNC__
!  #define __FUNC__ "Init"
!        IMPLICIT NONE
!        CLASS(SELFConfig), INTENT(out) :: this
!        ! Local
!        CHARACTER(LEN=SELF_FILE_DEFAULT_LENGTH) :: concretizationFile
!        CHARACTER(LEN=200) :: SELF_PREFIX
!        LOGICAL :: fileExists

!        ! Set Default configuration file to
!        ! ${SELF_PREFIX}/etc/schema/defaults/self.json
!        CALL get_environment_variable("SELF_PREFIX", SELF_PREFIX)
!        concretizationFile = TRIM(SELF_PREFIX)//"/etc/schema/defaults/self.json"

!        IF ( CommandLineArgumentIsPresent(argument = "-i") )     THEN
!           concretizationFile = StringValueForArgument(argument = "-i")
!        END IF
!        INFO("Using configuration file : "//TRIM(concretizationFile))
!        INQUIRE(FILE=TRIM(concretizationFile), EXIST=fileExists )
!        IF( fileExists )THEN
!          CALL this % LoadConcretization( TRIM(concretizationFile) )
!        ELSE
!           ERROR("Configuration file does not exist : "//TRIM(concretizationFile))
!           STOP 1
!        ENDIF

!     END SUBROUTINE Init_SELFConfig_FromCLI

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

  subroutine LoadConcretization_SELFConfig_FromFile(this,concretizationFile)
       !! Loads a concretization and stores in concretization attributes
    implicit none
    class(SELFConfig),intent(out) :: this
    character(*),intent(in) :: concretizationFile

    call this%concretization%initialize(stop_on_error=.true., &
                                        comment_char='#')

    call this%concretization%load_file(filename=trim(concretizationFile))

    !CALL get_environment_variable("SELF_DEBUG", SELF_DEBUG)
    !IF (SELF_DEBUG == 1)THEN
    !  CALL this % concretization % print_file()
    !ENDIF
  endsubroutine LoadConcretization_SELFConfig_FromFile

  subroutine Free_SELFConfig(this)
       !! Frees the attributes of the SELFConfig class and reset the config attribute
       !! to an empty string
    implicit none
    class(SELFConfig),intent(inout) :: this

    !CALL this % schema % destroy()
    call this%concretization%destroy()
    ! this % schemaFile = ""

  endsubroutine Free_SELFConfig

  subroutine Get_SELFConfig_int32(this,jsonKey,res)
    implicit none
    class(SELFConfig),intent(inout) :: this
    character(*),intent(in) :: jsonKey
    integer(int32),intent(out) :: res
    ! Local
    logical :: found

    call this%concretization%get(trim(jsonKey),res,found)
    if(.not. found) then
      print*,"JSON key not found : "//trim(jsonKey)
      stop 1
    endif

  endsubroutine Get_SELFConfig_int32

  subroutine Get_SELFConfig_int64(this,jsonKey,res)
    implicit none
    class(SELFConfig),intent(inout) :: this
    character(*),intent(in) :: jsonKey
    integer(int64),intent(out) :: res
    ! Local
    logical :: found
    integer(int32) :: res32

    call this%concretization%get(trim(jsonKey),res32,found)
    if(.not. found) then
      print*,"JSON key not found : "//trim(jsonKey)
      stop 1
    endif
    res = int(res32,kind=int64)

  endsubroutine Get_SELFConfig_int64

  subroutine Get_SELFConfig_real32(this,jsonKey,res)
    implicit none
    class(SELFConfig),intent(inout) :: this
    character(*),intent(in) :: jsonKey
    real(real32),intent(out) :: res
    ! Local
    logical :: found

    call this%concretization%get(trim(jsonKey),res,found)
    if(.not. found) then
      print*,"JSON key not found : "//trim(jsonKey)
      stop 1
    endif

  endsubroutine Get_SELFConfig_real32

  subroutine Get_SELFConfig_real64(this,jsonKey,res)
    implicit none
    class(SELFConfig),intent(inout) :: this
    character(*),intent(in) :: jsonKey
    real(real64),intent(out) :: res
    ! Local
    logical :: found

    call this%concretization%get(trim(jsonKey),res,found)
    if(.not. found) then
      print*,"JSON key not found : "//trim(jsonKey)
      stop 1
    endif

  endsubroutine Get_SELFConfig_real64

  subroutine Get_SELFConfig_logical(this,jsonKey,res)
    implicit none
    class(SELFConfig),intent(inout) :: this
    character(*),intent(in) :: jsonKey
    logical,intent(out) :: res
    ! Local
    logical :: found

    call this%concretization%get(trim(jsonKey),res,found)
    if(.not. found) then
      print*,"JSON key not found : "//trim(jsonKey)
      stop 1
    endif

  endsubroutine Get_SELFConfig_logical

  subroutine Get_SELFConfig_char(this,jsonKey,res)
    implicit none
    class(SELFConfig),intent(inout) :: this
    character(*),intent(in) :: jsonKey
    character(*),intent(out) :: res
    ! Local
    logical :: found
    character(LEN=:),allocatable :: resLoc

    call this%concretization%get(trim(jsonKey),resLoc,found)
    if(.not. found) then
      print*,"JSON key not found : "//trim(jsonKey)
      stop 1
    endif
    res = trim(resLoc)

  endsubroutine Get_SELFConfig_char

endmodule SELF_JSON_Config
