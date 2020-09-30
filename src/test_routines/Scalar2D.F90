! Scalar2D.F90
! 
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self-fluids@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

SUBROUTINE Scalar2D_Tests( )
  IMPLICIT NONE
  ! Local
  TYPE(EquationParser), ALLOCATABLE :: f(:), gradF(:,:)
  TYPE(JSON_FILE) :: json
  TYPE(JSON_VALUE), POINTER :: objPointer, testPointer
  TYPE(JSON_CORE) :: jCore
  INTEGER(JSON_IK) :: nFunctions, iTest
  LOGICAL :: found
  CHARACTER(50) :: funcName, funcValue, funcDx, funcDy
  INTEGER, ALLOCATABLE :: polyDeg(:) 
  INTEGER :: nPlotPoints, i, iEl, pdeg, nelements, nRepeats
  
    CALL json % Initialize()
    CALL json % Load(filename = './tests.json')
    CALL json % info('scalar_2d', n_children=nFunctions)
    CALL json % get('scalar_2d', objPointer, found)
    CALL json % get('polynomial_range', polyDeg, found)
    CALL json % get('n_plot_points', nPlotPoints, found)
    CALL json % get('element_dimensions', nElements, found)
    CALL json % get('n_timing_repeats', nRepeats, found)
    CALL json % get_core(jCore)
  
    ALLOCATE(f(1:nFunctions), gradF(1:2,1:nFunctions))
  
    DO iTest = 1, nFunctions
      ! Point to the i-th scalar_1d function for testing
      CALL jCore % get_child(objPointer, iTest, testPointer, found)
      IF( found )THEN
        ! Pull the function and derivative strings from the JSON
        CALL Get_Char_Obj(jCore, testPointer, 'name', funcName)
        CALL Get_Char_Obj(jCore, testPointer, 'function', funcValue)
        CALL Get_Char_Obj(jCore, testPointer, 'dfdx', funcDx)
        CALL Get_Char_Obj(jCore, testPointer, 'dfdy', funcDy)
        ! Create the exact function equation parsers
        f(iTest) = EquationParser(funcValue,(/'x','y'/))
        gradF(1,iTest) = EquationParser(funcDx,(/'x','y'/))
        gradF(2,iTest) = EquationParser(funcDy,(/'x','y'/))
      ENDIF
    ENDDO

    DO pdeg = polyDeg(1), polyDeg(2)
!      CALL ScalarDerivative_1D_Test( f, dfdx, dfdy, pdeg, nPlotPoints, nFunctions, nElements, nRepeats )
      CALL ScalarGridInterp_2D_Test( f, pdeg, nPlotPoints, nFunctions, nElements, nRepeats )
    ENDDO

    DEALLOCATE( f, gradF, polyDeg )

END SUBROUTINE Scalar2D_Tests

SUBROUTINE FillScalarValues2D( scalar, f, nEl )
  IMPLICIT NONE
  TYPE(Scalar2D), INTENT(inout) :: scalar
  TYPE(EquationParser), INTENT(in) :: f(1:scalar % nVar)
  INTEGER, INTENT(in) :: nEl ! number of elements in each direction
  ! Local
  INTEGER :: iEl, jEl, iVar, i, j, eid
  REAL(prec) :: xL, yL, xR, yR, x, y, ds

     ds = 2.0_prec/REAL( nEl, prec )
     DO jEl = 1, nEl
       yL = -1.0_prec + REAL((jEl-1),prec)*ds
       yR = yL + ds
       DO iEl = 1, nEl
         xL = -1.0_prec + REAL((iEl-1),prec)*ds
         xR = xL + ds

         eid = iEl + (jEl-1)*nEl

         DO iVar = 1, scalar % nVar

           DO j = 0, scalar % N
             y = 0.5_prec*( yR*(scalar % interp % controlPoints % hostData(j)+1.0_prec) - &
                            yL*(scalar % interp % controlPoints % hostData(j)-1.0_prec) )
             DO i = 0, scalar % N
               x = 0.5_prec*( xR*(scalar % interp % controlPoints % hostData(i)+1.0_prec) - &
                              xL*(scalar % interp % controlPoints % hostData(i)-1.0_prec) )

               scalar % interior % hostData(i,j,iVar,eid) = f(iVar) % Evaluate( (/x, y/) )          
             ENDDO
           ENDDO
         ENDDO
       ENDDO
     ENDDO

END SUBROUTINE FillScalarValues2D

SUBROUTINE ScalarGridInterp_2D_Test( f, N, M, nVar, nElem, nRepeats )
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
  TYPE(Scalar2D) :: fScalar, fInterp, fActual , fError
  REAL(prec) :: errors(1:nVar)
  REAL(prec) :: t1, t2, runtime
  INTEGER :: i, iVar

     cpuModel = GetCPUModel_Linux()
     gpuModel = 'None'

     DO iVar = 1, nVar
       functionNames(iVar) = f(iVar) % equation
     ENDDO

     CALL fScalar % Build( N, GAUSS, M, GAUSS, nVar, nElem*nElem ) 
     CALL fInterp % Build( M, GAUSS, M, GAUSS, nVar, nElem*nElem ) 
     CALL fActual % Build( M, GAUSS, M, GAUSS, nVar, nElem*nElem ) 
     CALL fError % Build( M, GAUSS, M, GAUSS, nVar, nElem*nElem ) 

     CALL FillScalarValues2D( fScalar, f, nElem )
#ifdef GPU
     CALL fScalar % UpdateDevice( )
#endif
     CALL FillScalarValues2D( fActual, f, nElem )

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
                      'GridInterp_Scalar2D', &
                      cpuModel, &
                      gpuModel, &
                      functionNames, &
                      runtime, &
                      errors, &
                      nVar, &
                      nElem*nElem)

#ifdef GPU
     CALL fScalar % GridInterp(fInterp, .TRUE.)
     CALL fInterp % UpdateHost( )
     CALL hipCheck(hipDeviceSynchronize()) 

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
                      'GridInterp_Scalar2D+gpuAccel', &
                      cpuModel, &
                      gpuModel, &
                      functionNames, &
                      runtime, &
                      errors, &
                      nVar, &
                      nElem*nElem)
     
#endif

     CALL fScalar % Trash( )
     CALL fInterp % Trash( )
     CALL fActual % Trash( )
     CALL fError % Trash( )

END SUBROUTINE ScalarGridInterp_2D_Test
