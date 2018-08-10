MODULE Fluid_EquationParser_Class

USE ModelPrecision
USE CommonRoutines
USE EquationParser_Class

IMPLICIT NONE

  TYPE Fluid_EquationParser
    TYPE( EquationParser ) :: u
    TYPE( EquationParser ) :: v
    TYPE( EquationParser ) :: w
    TYPE( EquationParser ) :: t
    TYPE( EquationParser ) :: rho
    TYPE( EquationParser ) :: topography
    TYPE( EquationParser ) :: drag
    TYPE( EquationParser ) :: tracer
    TYPE( EquationParser ) :: staticTracer
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
   CHARACTER(200) :: functionLine


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

         fluidEqs % u = EquationParser( functionLine )

       ELSEIF( functionLine(1:1) == 'b' )THEN

         fluidEqs % w = EquationParser( functionLine )

       ELSEIF( functionLine(1:1) == 'w' )THEN

         fluidEqs % w = EquationParser( functionLine )

       ELSEIF( functionLine(1:1) == 'r' )THEN

         fluidEqs % w = EquationParser( functionLine )
         fluidEqs % calculate_density_from_t = .FALSE.

       ELSEIF( functionLine(1:1) == 't' )THEN

         fluidEqs % t = EquationParser( functionLine )

       ELSEIF( functionLine(1:1) == 'h' )THEN

         fluidEqs % topography = EquationParser( functionLine )
         fluidEqs % topography_equation_provided = .TRUE.

       ELSEIF( functionLine(1:1) == 'c' )THEN

         fluidEqs % drag = EquationParser( functionLine )

       ELSEIF( functionLine(1:1) == 's' )THEN

         fluidEqs % tracer = EquationParser( functionLine )

       ELSEIF( functionLine(1:1) == 'S' )THEN

         fluidEqs % staticTracer = EquationParser( functionLine )

       ENDIF      

     ENDDO

     CLOSE( fUnit )

  END SUBROUTINE Build_Fluid_EquationParser

END MODULE Fluid_EquationParser_Class
