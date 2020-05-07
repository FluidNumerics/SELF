PROGRAM Lagrange_Test

USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
USE Lagrange_Class
USE EquationParser_Class

USE json_module

IMPLICIT NONE

  TYPE(JSON_FILE) :: json
  

    CALL Load_JSON(json)

    CALL Scalar_1D_Tests(json)

! Read interpolators.test.json
! For each scalar_1d
!   > Interpolate_1d + Validate
!   > ApplyInterpolationMatrix_1D + Validate
!   > CalculateDerivative_1D + Validate
!
! For each scalar_2d
!   > Interpolate_2D + Validate
!   > ApplyInterpolationMatrix_2D + Validate
!   > CalculateGradient_2D + Validate
!
! For each vector_2d
!   > CalculateDivergence_2D + Validate
!   > CalculateCurl_2D + Validate
!
! For each scalar_3d
!   > Interpolate_3D + Validate
!   > ApplyInterpolationMatrix_3D + Validate
!   > CalculateGradient_3D + Validate
!
! For each vector_3d
!   > CalculateDivergence_3D + Validate
!   > CalculateCurl_3D + Validate

    CALL json % Destroy()
    IF(json % Failed()) STOP 1

CONTAINS

  SUBROUTINE Load_JSON(json)
    IMPLICIT NONE
    TYPE(JSON_FILE), INTENT(out) :: json

      CALL json % Initialize()
      CALL json % Load(filename = './lagrange.test.json')
      !CALL json % Print()

  END SUBROUTINE Load_JSON
 
  SUBROUTINE Scalar_1D_Tests(json)
    IMPLICIT NONE
    TYPE(JSON_FILE), INTENT(inout) :: json
    ! Local
    TYPE(EquationParser) :: f, dfdx
    TYPE(Lagrange) :: interp
    TYPE(JSON_VALUE), POINTER :: objPointer, testPointer
    TYPE(JSON_CORE) :: jCore
    INTEGER(JSON_IK) :: nTests, iTest
    LOGICAL :: found
    CHARACTER(50) :: funcName, funcValue, funcDerivative 
    INTEGER, ALLOCATABLE :: polyDeg(:) 
    INTEGER :: nPlotPoints, i, pdeg
    REAL(prec), ALLOCATABLE :: interpNodes(:), qWeights(:), targetNodes(:)
    REAL(prec), ALLOCATABLE :: fNodal(:,:,:), fInterp(:,:,:), dfdxInterp(:,:,:),  fActual(:), dfdxActual(:)  
  
      CALL json % info('scalar_1d', n_children=nTests)
      CALL json % get('scalar_1d', objPointer, found)
      CALL json % get('polynomial_range', polyDeg, found)
      CALL json % get('n_plot_points', nPlotPoints, found)
      CALL json % get_core(jCore)

      ALLOCATE(targetNodes(0:nPlotPoints), fActual(0:nPlotPoints), fInterp(0:nPlotPoints,1,1))
      targetNodes = UniformPoints(-1.0_prec,1.0_prec,nPlotPoints)

      DO iTest = 1, nTests
 
        ! Point to the i-th scalar_1d function for testing
        CALL jCore % get_child(objPointer, iTest, testPointer, found)
        IF( found )THEN

          ! Pull the function and derivative strings from the JSON
          CALL Get_Char_Obj(jCore, testPointer, 'name', funcName)
          CALL Get_Char_Obj(jCore, testPointer, 'function', funcValue)
          CALL Get_Char_Obj(jCore, testPointer, 'derivative', funcDerivative)

          ! Create the exact function equation parsers
          f = EquationParser(funcValue)
          dfdx = EquationParser(funcDerivative)

          ! Create the exact function equation parsers
          DO i = 0, nPlotPoints
            fActual(i) = f % Evaluate( (// targetNodes(i), 0.0_prec, 0.0_prec //) )          
          ENDDO


          DO pdeg = polyDeg(1), polyDeg(2)

             ALLOCATE(interpNodes(0:pdeg), qWeights(0:pdeg), fNodal(0:pdeg,1,1), dfdxInterp(0:pdeg,1,1), dfdxActual(0:pdeg))

             CALL LegendreQuadrature( pdeg, interpNodes, qWeights, GAUSS )
             CALL interp % Build_Lagrange(pdeg, nTargetNodes, interpNodes, targetNodes)

             ! Evaluate the nodal values at quadrature points and the exact derivative
             DO i = 0, pdeg
               fNodal(i,1,1) = f % Evaluate( (// interpNodes(i), 0.0_prec, 0.0_prec //) )          
               dfdxActual(i) = dfdx % Evaluate( (// interpNodes(i), 0.0_prec, 0.0_prec //) )          
             ENDDO

             ! Interpolate the function to the targetNodes
             CALL interp % ApplyInterpolationMatrix_1D(fNodal, fInterp, 1, 1)  
             ! Estimate the derivative by applying the derivative matrix
             CALL interp % CalculateDerivative_1D(fNodal, dfdxInterp, 1, 1)

             ! Calculate Errors and report to screen

             DEALLOCATE(interpNodes, qWeights, fNodal, dfdxInterp, dfdxActual)

          ENDDO

        ELSE
          PRINT*, 'FAIL!'
          STOP 1
        ENDIF
        
        CALL jCore % destroy(testPointer)

      ENDDO
    
      DEALLOCATE(targetNodes, fActual, fInterp)
      CALL jCore % destroy(objPointer)
      CALL jCore % destroy()
      CALL json % destroy()

  END SUBROUTINE Scalar_1D_Tests

  SUBROUTINE Get_Char_Obj(jCore, testPointer, key, charValue)
    IMPLICIT NONE
    TYPE(JSON_CORE), INTENT(inout) :: jCore
    TYPE(JSON_VALUE), POINTER, INTENT(inout) :: testPointer
    CHARACTER(*), INTENT(in) :: key
    CHARACTER(*), INTENT(out) :: charValue
    ! Local
    TYPE(JSON_VALUE), POINTER :: p
    CHARACTER(KIND=JSON_CK, LEN=:), ALLOCATABLE :: var
    LOGICAL :: found

      CALL jCore % get_child(testPointer, TRIM(key), p, Found)
      IF( found )THEN
        CALL jCore % get(p, var)
        charValue = TRIM(var)
      ENDIF
      CALL jCore % destroy(p)

  END SUBROUTINE Get_Char_Obj

! SUBROUTINE Scalar_2D_Tests()
! END SUBROUTINE Scalar_2D_Tests
!
! SUBROUTINE Scalar_3D_Tests()
! END SUBROUTINE Scalar_3D_Tests
!
! SUBROUTINE Vector_2D_Tests()
! END SUBROUTINE Vector_2D_Tests
!
! SUBROUTINE Vector_3D_Tests()
! END SUBROUTINE Vector_3D_Tests

END PROGRAM Lagrange_Test
