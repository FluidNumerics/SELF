! Scalar1D.F90
! 
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self-fluids@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

SUBROUTINE Scalar1D_Tests( )
  IMPLICIT NONE
  ! Local
  TYPE(EquationParser), ALLOCATABLE :: f(:), dfdx(:)
  TYPE(JSON_FILE) :: json
  TYPE(JSON_VALUE), POINTER :: objPointer, testPointer
  TYPE(JSON_CORE) :: jCore
  INTEGER(JSON_IK) :: nFunctions, iTest
  LOGICAL :: found
  CHARACTER(50) :: funcName, funcValue, funcDerivative 
  INTEGER, ALLOCATABLE :: polyDeg(:) 
  INTEGER :: nPlotPoints, i, iEl, pdeg, nelements, nRepeats
  
    CALL json % Initialize()
    CALL json % Load(filename = './tests.json')
    CALL json % info('scalar_1d', n_children=nFunctions)
    CALL json % get('scalar_1d', objPointer, found)
    CALL json % get('polynomial_range', polyDeg, found)
    CALL json % get('n_plot_points', nPlotPoints, found)
    CALL json % get('element_dimensions', nElements, found)
    CALL json % get('n_timing_repeats', nRepeats, found)
    CALL json % get_core(jCore)
  
    ALLOCATE(f(1:nFunctions), dfdx(1:nFunctions))
  
    DO iTest = 1, nFunctions
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

    DO pdeg = polyDeg(1), polyDeg(2)
      CALL ScalarDerivative_1D_Test( f, dfdx, pdeg, nPlotPoints, nFunctions, nElements, nRepeats )
      CALL ScalarGridInterp_1D_Test( f, pdeg, nPlotPoints, nFunctions, nElements, nRepeats )
    ENDDO

    DEALLOCATE( f, dfdx, polyDeg )

END SUBROUTINE Scalar1D_Tests

SUBROUTINE FillScalarValues1D( scalar, f )
  IMPLICIT NONE
  TYPE(Scalar1D), INTENT(inout) :: scalar
  TYPE(EquationParser), INTENT(in) :: f(1:scalar % nVar)
  ! Local
  INTEGER :: iEl, iVar, i
  REAL(prec) :: xL, xR, x, dx

     dx = 2.0_prec/REAL( scalar % nElem, prec )
     DO iEl = 1, scalar % nElem
       DO iVar = 1, scalar % nVar
         DO i = 0, scalar % N
           xL = -1.0_prec + REAL((iEl-1),prec)*dx
           xR = xL + dx
           x = 0.5_prec*( xR*(scalar % interp % controlPoints % hostData(i)+1.0_prec) - &
                          xL*(scalar % interp % controlPoints % hostData(i)-1.0_prec) )
           scalar % interior % hostData(i,iVar,iEl) = f(iVar) % Evaluate( (/x/) )          
         ENDDO
       ENDDO
     ENDDO

END SUBROUTINE FillScalarValues1D

SUBROUTINE ScalarDerivative_1D_Test( f, dfdx, N, M, nVar, nElem, nRepeats )
  IMPLICIT NONE
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem
  INTEGER, INTENT(in) :: nRepeats
  TYPE(EquationParser), INTENT(in)  :: f(1:nVar)
  TYPE(EquationParser), INTENT(in)  :: dfdx(1:nVar)
  ! Local
  CHARACTER(LEN=50), DIMENSION(1:nVar) :: functionNames
  CHARACTER(LEN=50) :: cpuModel
  CHARACTER(LEN=50) :: gpuModel
  TYPE(Scalar1D) :: fScalar, dfdxInterp, dfdxActual, dfdxError
  REAL(prec) :: errors(1:nVar)
  REAL(prec) :: t1, t2, runtime
  INTEGER :: i, iVar

     cpuModel = GetCPUModel_Linux()
     gpuModel = 'None'

     DO iVar = 1, nVar
       functionNames(iVar) = f(iVar) % equation
     ENDDO

     CALL fScalar % Build( N, GAUSS, M, GAUSS, nVar, nElem ) 
     CALL dfdxInterp % Build( N, GAUSS, M, GAUSS, nVar, nElem ) 
     CALL dfdxActual % Build( N, GAUSS, M, GAUSS, nVar, nElem ) 
     CALL dfdxError % Build( N, GAUSS, M, GAUSS, nVar, nElem ) 

     CALL FillScalarValues1D( fScalar, f )
#ifdef GPU
     CALL fScalar % UpdateDevice( )
#endif
     CALL FillScalarValues1D( dfdxActual, dfdx )

     CALL fScalar % Derivative(dfdxInterp, .FALSE.)
     dfdxInterp % interior % hostData = dfdxInterp % interior % hostData*REAL( fScalar % nElem, prec ) 

     dfdxError = dfdxActual - dfdxInterp
     errors = dfdxError % AbsMaxInterior( )

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL fScalar % Derivative(dfdxInterp, .FALSE.)
     ENDDO
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec

     CALL Report_JSON(dfdxInterp % interp, &
                      'Derivative_1D', &
                      cpuModel, &
                      gpuModel, &
                      functionNames, &
                      runtime, &
                      errors, &
                      nVar, &
                      nElem)


#ifdef GPU
     CALL fScalar % Derivative(dfdxInterp, .TRUE.)
     CALL dfdxInterp % UpdateHost( )
     dfdxInterp % interior % hostData = dfdxInterp % interior % hostData*REAL( fScalar % nElem, prec ) 

     dfdxError = dfdxActual - dfdxInterp
     errors = dfdxError % AbsMaxInterior( )

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL fScalar % Derivative(dfdxInterp, .TRUE.)
     ENDDO
     CALL hipCheck(hipDeviceSynchronize()) 
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec


     CALL Report_JSON(dfdxInterp % interp, &
                      'Derivative_1D+gpuAccel', &
                      cpuModel, &
                      gpuModel, &
                      functionNames, &
                      runtime, &
                      errors, &
                      nVar, &
                      nElem)
     
#endif

     CALL fScalar % Trash( )
     CALL dfdxInterp % Trash( )
     CALL dfdxActual % Trash( )
     CALL dfdxError % Trash( )

END SUBROUTINE ScalarDerivative_1D_Test

SUBROUTINE ScalarGridInterp_1D_Test( f, N, M, nVar, nElem, nRepeats )
  IMPLICIT NONE
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem
  INTEGER, INTENT(in) :: nRepeats
  TYPE(EquationParser), INTENT(in)  :: f(1:nVar)
  ! Local
  CHARACTER(LEN=50), DIMENSION(1:nVar) :: functionNames
  CHARACTER(LEN=50) :: cpuModel
  CHARACTER(LEN=50) :: gpuModel
  TYPE(Scalar1D) :: fScalar, fInterp, fActual , fError
  REAL(prec) :: errors(1:nVar)
  REAL(prec) :: t1, t2, runtime
  INTEGER :: i, iVar

     cpuModel = GetCPUModel_Linux()
     gpuModel = 'None'

     DO iVar = 1, nVar
       functionNames(iVar) = f(iVar) % equation
     ENDDO

     CALL fScalar % Build( N, GAUSS, M, GAUSS, nVar, nElem ) 
     CALL fInterp % Build( M, GAUSS, M, GAUSS, nVar, nElem ) 
     CALL fActual % Build( M, GAUSS, M, GAUSS, nVar, nElem ) 
     CALL fError % Build( M, GAUSS, M, GAUSS, nVar, nElem ) 

     CALL FillScalarValues1D( fScalar, f )
#ifdef GPU
     CALL fScalar % UpdateDevice( )
#endif
     CALL FillScalarValues1D( fActual, f )

     CALL fScalar % GridInterp(fInterp, .FALSE.)

     fError = fActual - fInterp
     errors = fError % AbsMaxInterior( )

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL fScalar % GridInterp(fInterp, .FALSE.)
     ENDDO
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec

     CALL Report_JSON(fScalar % interp, &
                      'GridInterp_Scalar1D', &
                      cpuModel, &
                      gpuModel, &
                      functionNames, &
                      runtime, &
                      errors, &
                      nVar, &
                      nElem)


#ifdef GPU
     CALL fScalar % GridInterp(fInterp, .TRUE.)
     CALL fInterp % UpdateHost( )

     fError = fActual - fInterp
     errors = fError % AbsMaxInterior( )

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL fScalar % GridInterp(fInterp, .TRUE.)
     ENDDO
     CALL hipCheck(hipDeviceSynchronize()) 
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec


     CALL Report_JSON(fScalar % interp, &
                      'GridInterp_Scalar1D+gpuAccel', &
                      cpuModel, &
                      gpuModel, &
                      functionNames, &
                      runtime, &
                      errors, &
                      nVar, &
                      nElem)
     
#endif

     CALL fScalar % Trash( )
     CALL fInterp % Trash( )
     CALL fActual % Trash( )
     CALL fError % Trash( )

END SUBROUTINE ScalarGridInterp_1D_Test
