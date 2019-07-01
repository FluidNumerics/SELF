MODULE Geom_EquationParser_Class

USE ModelPrecision
USE CommonRoutines
USE EquationParser_Class

IMPLICIT NONE

  TYPE Geom_EquationParser
    TYPE( EquationParser ) :: topography
    INTEGER                :: boundaryConditionFlags(1:6)
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
   LOGICAL :: boundaryRead, bracketOpen

     geomEqs % boundaryConditionFlags(1:6) = NO_NORMAL_FLOW

     OPEN( UNIT=NewUnit(fUnit), &
           FILE=TRIM(icFile), &
           FORM='FORMATTED', &
           ACCESS='SEQUENTIAL', &
           ACTION='READ' )

     ioErr = 0
     boundaryRead = .FALSE.
     bracketOpen  = .FALSE.

     DO WHILE( ioErr == 0 )

       READ( fUnit, '(A200)', ioStat=ioErr ) functionLine

       IF( functionLine(1:8) == 'boundary' )THEN

         boundaryRead = .TRUE.
         CYCLE

       ENDIF

       IF( functionLine(1:8) == '{' )THEN

         bracketOpen = .TRUE.
         CYCLE

       ELSEIF( functionLine(1:1) == '}' )THEN
         
         bracketOpen = .FALSE.
         IF( boundaryRead )THEN 
           boundaryRead=.FALSE.
         ENDIF
         CYCLE

       ENDIF
       
       IF( boundaryRead .AND. bracketOpen )THEN
         IF( functionLine(1:1) == 'h' )THEN
  
           geomEqs % topography = EquationParser( functionLine )
  
         ELSEIF( functionLine(1:5) == 'south' )THEN
  
           geomEqs % boundaryConditionFlags(1) = Parse_BCLine( functionLine )
  
         ELSEIF( functionLine(1:4) == 'east' )THEN
  
           geomEqs % boundaryConditionFlags(2) = Parse_BCLine( functionLine )
  
         ELSEIF( functionLine(1:5) == 'north' )THEN
  
           geomEqs % boundaryConditionFlags(3) = Parse_BCLine( functionLine )
  
         ELSEIF( functionLine(1:4) == 'west' )THEN
  
           geomEqs % boundaryConditionFlags(4) = Parse_BCLine( functionLine )
  
         ELSEIF( functionLine(1:6) == 'bottom' )THEN
  
           geomEqs % boundaryConditionFlags(5) = Parse_BCLine( functionLine )
  
         ELSEIF( functionLine(1:3) == 'top' )THEN
  
           geomEqs % boundaryConditionFlags(6) = Parse_BCLine( functionLine )

         ELSEIF( functionLine(1:5) == 'sides' )THEN
 
           geomEqs % boundaryConditionFlags(1) = Parse_BCLine( functionLine )
           geomEqs % boundaryConditionFlags(2) = Parse_BCLine( functionLine )
           geomEqs % boundaryConditionFlags(3) = Parse_BCLine( functionLine )
           geomEqs % boundaryConditionFlags(4) = Parse_BCLine( functionLine )
  
         ENDIF      
       ENDIF

     ENDDO

     CLOSE( fUnit )

  END SUBROUTINE Build_Geom_EquationParser

  FUNCTION Parse_BCLine( bcLine ) RESULT( bcFlag )
    CHARACTER(*) :: bcLine
    INTEGER      :: bcflag
    ! Local
    INTEGER       :: nChar, rhsLoc, i
    CHARACTER(20) :: bcString

        nChar = LEN_TRIM( bcLine )
        rhsLoc = INDEX( bcLine, "=" ) + 1

        DO i = 1, 3
          IF( bcLine(rhsLoc:rhsLoc) == ' ' )THEN
            rhsLoc = rhsLoc + 1
          ELSE
            EXIT
          ENDIF
        ENDDO

        bcString = bcLine(rhsLoc:nChar) 

        IF( TRIM( UpperCase(bcString) ) == 'NO_NORMAL_FLOW' )THEN

          bcFlag = NO_NORMAL_FLOW

        ELSEIF( TRIM( UpperCase(bcString) ) == 'PRESCRIBED' )THEN

          bcFlag = PRESCRIBED

        ELSEIF( TRIM( UpperCase(bcString) ) == 'RADIATION' )THEN

          bcFlag = RADIATION

        ELSEIF( TRIM( UpperCase(bcString) ) == 'PERIODIC' )THEN

          bcFlag = PERIODIC

        ENDIF

  END FUNCTION Parse_BCLine

END MODULE Geom_EquationParser_Class
