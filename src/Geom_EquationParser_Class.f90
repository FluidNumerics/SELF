MODULE Geom_EquationParser_Class

USE ModelPrecision
USE CommonRoutines
USE EquationParser_Class

IMPLICIT NONE

  TYPE Geom_EquationParser
    TYPE( EquationParser ) :: topography

    CONTAINS

      PROCEDURE :: Build => Build_Geom_EquationParser

  END TYPE Geom_EquationParser

CONTAINS

  SUBROUTINE Build_Geom_EquationParser( geomEqs, icFile )
   CLASS( Geom_EquationParser ), INTENT(out) :: geomEqs
   CHARACTER(*), INTENT(IN)                   :: icFile
   ! Local
   INTEGER :: fUnit, ioErr
   CHARACTER(200) :: functionLine

     OPEN( UNIT=NewUnit(fUnit), &
           FILE=TRIM(icFile), &
           FORM='FORMATTED', &
           ACCESS='SEQUENTIAL', &
           ACTION='READ' )

     ioErr = 0
     DO WHILE( ioErr == 0 )

       READ( fUnit, '(A200)', ioStat=ioErr ) functionLine

       IF( functionLine(1:1) == 'h' )THEN

         geomEqs % topography = EquationParser( functionLine )

       ENDIF      

     ENDDO

     CLOSE( fUnit )

  END SUBROUTINE Build_Geom_EquationParser

END MODULE Geom_EquationParser_Class
