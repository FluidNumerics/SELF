! Copyright 2020 Fluid Numerics LLC
! All Rights Reserved.
!
! A module for obtain system configuration information

MODULE SysConf


CONTAINS


FUNCTION GetCPUModel_Linux() RESULT(cpuModel)
! Reads /proc/cpuinfo and parses for line with model_name to gather cpuModel
  IMPLICIT NONE
  CHARACTER(50) :: cpuModel
  ! Local
  INTEGER :: ios, separatorLoc
  CHARACTER(100) :: line
  CHARACTER(30)  :: key
  CHARACTER(50)  :: val

  cpuModel = 'unknown'

  OPEN(UNIT=100, FILE='/proc/cpuinfo', STATUS='OLD', ACTION='READ', IOSTAT=ios)

  DO WHILE (ios == 0)
    READ(100,'(A100)', IOSTAT=ios) line 
    
    IF( ios == 0 )THEN
      separatorLoc = INDEX( line, ":" )

      IF( separatorLoc /= 0 )THEN
        key = TRIM( line(1:separatorLoc-1) )
        val = TRIM( line(separatorLoc+1:) )
        IF( INDEX(TRIM(key),'model name') == 1 )THEN
          cpuModel = TRIM(val)
          EXIT
        ENDIF
      ENDIF

    ENDIF

  ENDDO
  CLOSE(UNIT=100)

END FUNCTION GetCPUModel_Linux


END MODULE SysConf
