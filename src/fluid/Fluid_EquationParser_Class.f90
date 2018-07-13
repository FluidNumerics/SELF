MODULE Fluid_EquationParser_Class

USE ModelPrecision
USE CommonRoutines
USE FortranParser, only : EquationParser

IMPLICIT NONE

  TYPE Fluid_EquationParser
    TYPE( EquationParser ) :: u
    TYPE( EquationParser ) :: v
    TYPE( EquationParser ) :: w
    TYPE( EquationParser ) :: t
    TYPE( EquationParser ) :: rho
    TYPE( EquationParser ) :: topography
    LOGICAL                :: calculate_density_from_T
    LOGICAL                :: topography_equation_provided

    CONTAINS

      PROCEDURE :: Build => Build_Fluid_EquationParser

  END TYPE Fluid_EquationParser

CONTAINS

  SUBROUTINE Build_Fluid_EquationParser( fluidEqs, icFile )
   CLASS( Fluid_EquationParser ), INTENT(out) :: fluidEqs
   CHARACTER(*), INTENT(IN)                   :: icFile
   ! Local
   INTEGER :: fUnit, ioErr
   CHARACTER(200) :: functionLine, actualFunction
   CHARACTER(10)  :: variables(3)
   CHARACTER(10)  :: variables_2d(2)

     variables    = ['x', 'y', 'z']
     variables_2d = ['x', 'y']

     fluidEqs % calculate_density_from_T = .TRUE.    
     fluidEqs % topography_equation_provided = .FALSE.    

     OPEN( UNIT=NewUnit(fUnit), &
           FILE=TRIM(icFile), &
           FORM='FORMATTED', &
           ACCESS='SEQUENTIAL', &
           ACTION='READ' )

     ioErr = 0
     DO WHILE( ioErr == 0 )

       READ( fUnit, '(A200)', ioStat=ioErr ) functionLine

       IF( functionLine(1:1) == 'u' )THEN

         WRITE( actualFunction, '("u = ",A200)' ) functionLine
         fluidEqs % u = EquationParser( actualFunction, variables )

       ELSEIF( functionLine(1:1) == 'b' )THEN

         WRITE( actualFunction, '("w = ",A200)' ) functionLine
         fluidEqs % w = EquationParser( actualFunction, variables )

       ELSEIF( functionLine(1:1) == 'w' )THEN

         WRITE( actualFunction, '("w = ",A200)' ) functionLine
         fluidEqs % w = EquationParser( actualFunction, variables )

       ELSEIF( functionLine(1:1) == 'r' )THEN

         WRITE( actualFunction, '("rho = ",A200)' ) functionLine
         fluidEqs % w = EquationParser( actualFunction, variables )
         fluidEqs % calculate_density_from_t = .FALSE.

       ELSEIF( functionLine(1:1) == 't' )THEN

         WRITE( actualFunction, '("t = ",A200)' ) functionLine
         fluidEqs % t = EquationParser( actualFunction, variables )

       ELSEIF( functionLine(1:1) == 'h' )THEN

         WRITE( actualFunction, '("h = ",A200)' ) functionLine
         fluidEqs % topography = EquationParser( actualFunction, variables_2d )
         fluidEqs % topography_equation_provided = .TRUE.
       ENDIF      

     ENDDO

  END SUBROUTINE Build_Fluid_EquationParser

END MODULE Fluid_EquationParser_Class
