! CompareModelDataInstances_Driver.f90
! 
!  Copyright (2017), Joseph Schoonover, Cooperative Institute for Research in
!  Environmental Sciences, NOAA, (joseph.schoonover@noaa.gov)
!
!  Author : Joseph Schoonover
!  Creation Date : July 25, 2017
!
! ------------------------------------------------------------------------------------------------ !
!
! This program reads in model data instances from 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

PROGRAM CompareModelDataInstances

 USE ModelPrecision
 USE ModelDataInstances_Class

 IMPLICIT NONE

   TYPE(ModelDataInstances) :: mdi1, mdi2
   INTEGER                  :: nTotalObs, i
   CHARACTER(50)            :: dir1, dir2, baseName
   CHARACTER(8)             :: date
   CHARACTER(10)            :: time
   CHARACTER(5)             :: zone
   LOGICAL                  :: setupFail, fileExists
   


     ! Default settings
     dir1      = './'
     dir2      = './'
     baseName  = 'foo'
     nTotalObs = 0
     setupFail = .FALSE.

     ! Set up program from command line arguments
     CALL InitializeFromCommandLine( )

     CALL DATE_AND_TIME( DATE=date, TIME=time, ZONE=zone )
     PRINT '(3x,A,2x,A,2x,A)' , date, time, zone

     PRINT*, 'dir1 = '//TRIM(dir1)
     PRINT*, 'dir2 = '//TRIM(dir2)
     PRINT '(A,1x,I5)', 'nTotalObs = ', nTotalObs

     IF( setupFail )THEN
        PRINT*, '  cmdi setup failed!'
        STOP
     ENDIF

     

     ! Initialize the ModelDataInstances linked list, by pointing all of the
     ! pointers to Null
     CALL mdi1 % Build( ) ! mdi1 should refer to "new" version of data
     CALL mdi2 % Build( ) ! mdi2 should refer to "old" version of data

     DO i = 1, nTotalObs

        CALL mdi1 % Read_ModelDataInstances( TRIM(dir1)//'/'//TRIM(baseName), &
                                             i, fileExists ) 
        CALL mdi2 % Read_ModelDataInstances( TRIM(dir2)//'/'//TRIM(baseName), &
                                             i, fileExists ) 

        CALL mdi1 % CompareWith( mdi2 )

     ENDDO

     CALL mdi1 % Trash( )
     CALL mdi2 % Trash( )

CONTAINS

 SUBROUTINE InitializeFromCommandLine( )
   IMPLICIT NONE

   INTEGER       :: nArg, argID
   CHARACTER(50) :: argname
   LOGICAL       :: dir1Given, dir2Given, baseNameGiven, nGiven


     dir1Given     = .FALSE.
     dir2Given     = .FALSE.
     baseNameGiven = .FALSE.
     nGiven        = .FALSE.

     nArg = command_argument_count( )

     IF( nArg > 0 )THEN
 
        DO argID = 1, nArg
  
          CALL get_command_argument( argID, argName )

          SELECT CASE( TRIM(argName) )

             CASE("--dir1")
                dir1Given = .TRUE.
             CASE("--dir2")
                dir2Given = .TRUE.
             CASE("--basename")
                baseNameGiven = .TRUE.
             CASE("--n")
                nGiven = .TRUE.

             CASE DEFAULT

               IF( dir1Given )THEN

                  dir1      = TRIM(argName) ! Capture the directory name
                  dir1Given = .FALSE.

               ELSEIF( dir2Given )THEN

                  dir2      = TRIM(argName)
                  dir2Given = .FALSE.

               ELSEIF( baseNameGiven )THEN

                  baseName      = TRIM(argName)
                  baseNameGiven = .FALSE.

               ELSEIF( nGiven )THEN
               
                  READ( argName, '(I50)') nTotalObs
                  nGiven = .FALSE.
 
               ENDIF

          END SELECT 
        ENDDO

        
        INQUIRE( FILE = TRIM(dir1)//'/'//TRIM(baseName)//'.mdi.hdr', &
                 EXIST = fileExists )
        IF( .NOT.(fileExists) )THEN
          PRINT*, ' MDI Header file not found :'//TRIM(dir1)//TRIM(baseName)//'.mdi.hdr'
          setupFail = .TRUE.
        ENDIF

        INQUIRE( FILE = TRIM(dir2)//'/'//TRIM(baseName)//'.mdi.hdr', &
                 EXIST = fileExists )
        IF( .NOT.(fileExists) )THEN
          PRINT*, ' MDI Header file not found :'//TRIM(dir2)//TRIM(baseName)//'.mdi.hdr'
          setupFail = .TRUE.
        ENDIF

        IF( nTotalObs == 0 )THEN
          PRINT*, ' Number of observations need to be specified.'
          setupFail = .TRUE.
        ENDIF
        
     ELSE

        ! List possible options
        PRINT*, '  cmdi : Compare Model Data Instances'
        PRINT*, '    A tool for aiding in integration testing.'
        PRINT*, '--------------------------------------------------------------'
        PRINT*, '  Usage : cmdi --dir1 /path/to/mdifiles '
        PRINT*, '               --dir2 /path/to/other/mdifiles'
        PRINT*, '               --basename < .mdi and .mdi.hdr prefix >'
        PRINT*, '               --n < number of *.mdi files >'
        PRINT*, '--------------------------------------------------------------'
        PRINT*, ' This executable compares two sets of mdi header and binary '
        PRINT*, ' files and reports the differences of recorded model instances.'
        PRINT*, '--------------------------------------------------------------'
        STOP
        
     ENDIF   


 END SUBROUTINE InitializeFromCommandLine

END PROGRAM CompareModelDataInstances
