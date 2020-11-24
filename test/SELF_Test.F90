MODULE SELF_Test

USE SELF_Constants
USE SELF_Memory
USE SELF_SupportRoutines
!USE FEQParse
!USE json_fortran

IMPLICIT NONE

  TYPE SELFTest
   INTEGER :: Nc
   INTEGER :: Nt
   INTEGER :: qType
   INTEGER :: Ne
   INTEGER :: Nf
   INTEGER :: nEq, nDim
   TYPE(EquationParser) :: f(:)
   TYPE(EquationParser) :: df(:)

   CONTAINS

     PROCEDURE, PUBLIC :: ParseCLI => ParseCLI_SELFTest
     PROCEDURE, PUBLIC :: LoadJSON => LoadJSON_SELFTest

     PROCEDURE, PUBLIC :: L2Error => L2Error_Scalar1D, &
                                     L2Error_Scalar2D, &
                                     L2Error_Scalar3D, &
                                     L2Error_Vector2D, &
                                     L2Error_Vector3D, &
                                     L2Error_Tensor2D, &
                                     L2Error_Tensor3D

     PROCEDURE, PUBLIC :: L2DError => L2DError_Scalar1D, &
                                      L2DError_Scalar2D, &
                                      L2DError_Scalar3D, &
                                      L2DError_Vector2D, &
                                      L2DError_Vector3D, &
                                      L2DError_Tensor2D, &
                                      L2DError_Tensor3D

  END TYPE SELFTest

CONTAINS


 SUBROUTINE ParseCLI_SELFTest(sTest)
   IMPLICIT NONE
   CLASS(SELFTest), INTENT(inout) :: sTest
   ! Local
   INTEGER :: nArg, argId
   CHARACTER(6) :: arg
   LOGICAL :: c, t, q, e, f

   nArg = command_argument_count()

   IF( nArg > 0 )THEN
     DO argId = 1, nArg
       CALL get_command_argument(argId, arg)

       SELECT CASE( TRIM(arg) )

         CASE("-c") ! Number of Control Points
           c = .TRUE.
         CASE("-t") ! Number of Target Points
           t = .TRUE.
         CASE("-q") ! Quadrature type
           q = .TRUE.
         CASE("-e") ! Number of Elements
           e = .TRUE.
         CASE("-f") ! Number of functions/nVar
           f = .TRUE.

         CASE DEFAULT

           IF(c)THEN
             READ(arg,*) sTest % Nc
             c = .FALSE. 
           ENDIF

           IF(t)THEN
             READ(arg,*) sTest % Nt
             t = .FALSE. 
           ENDIF

           IF(q)THEN
             READ(arg,*) sTest % qType
             q = .FALSE. 
           ENDIF

           IF(e)THEN
             READ(arg,*) sTest % Ne
             e = .FALSE. 
           ENDIF

           IF(f)THEN
             READ(arg,*) sTest % Nf
             f = .FALSE. 
           ENDIF

       END SELECT
     ENDDO
   ENDIF

 END SUBROUTINE ParseCLI_SELFTest

 SUBROUTINE LoadJSON_SELFTest( sTest, testType )
   IMPLICIT NONE
   CLASS(SELFTest), INTENT(inout) :: sTest
   CHARACTER(*), INTENT(in) :: testType
   ! Local
   TYPE(JSON_FILE) :: json
   TYPE(JSON_VALUE), POINTER :: objPointer, testPointer
   TYPE(JSON_CORE) :: jCore
   INTEGER(JSON_IK) :: iTest
   LOGICAL :: found
  
    CALL json % Initialize()
    CALL json % Load(filename = './tests.json')
    CALL json % info('scalar_1d', n_children=sTest % nEq)
    CALL json % get('scalar_1d', objPointer, found)
    CALL json % get_core(jCore)
  
    ALLOCATE( sTest % f(1:sTest % nEq),&
              sTest % dfdx(1:sTest % nEq) )
    

    IF( TRIM(testType) == 'Scalar1D' )THEN  
      DO iTest = 1, sTest % nEq
        ! Point to the i-th scalar_1d function for testing
        CALL jCore % get_child(objPointer, iTest, testPointer, found)
        IF( found )THEN
          ! Pull the function and derivative strings from the JSON
          CALL Get_Char_Obj(jCore, testPointer, 'name', funcName)
          CALL Get_Char_Obj(jCore, testPointer, 'function', funcValue)
          CALL Get_Char_Obj(jCore, testPointer, 'dfdx', funcDerivative)
          ! Create the exact function equation parsers
          f(iTest) = EquationParser(funcValue,(/'x'/))
          dfdx(iTest) = EquationParser(funcDerivative,(/'x'/))
        ENDIF
      ENDDO

   ! TO DO, Determine number of dimensions

   !
 END SUBROUTINE LoadJSON_SELFTest


END MODULE SELF_Test
