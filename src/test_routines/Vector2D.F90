! Vector2D.F90
! 
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self-fluids@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

SUBROUTINE Vector2D_Tests( )
  IMPLICIT NONE
  ! Local
  TYPE(EquationParser), ALLOCATABLE :: f(:,:), divF(:)
  TYPE(JSON_FILE) :: json
  TYPE(JSON_VALUE), POINTER :: objPointer, testPointer
  TYPE(JSON_CORE) :: jCore
  INTEGER(JSON_IK) :: nFunctions, iTest
  LOGICAL :: found
  CHARACTER(50) :: funcName, fx, fy, dF
  INTEGER, ALLOCATABLE :: polyDeg(:) 
  INTEGER :: nPlotPoints, i, iEl, pdeg, nelements, nRepeats
  
    CALL json % Initialize()
    CALL json % Load(filename = './tests.json')
    CALL json % info('vector_2d', n_children=nFunctions)
    CALL json % get('vector_2d', objPointer, found)
    CALL json % get('polynomial_range', polyDeg, found)
    CALL json % get('n_plot_points', nPlotPoints, found)
    CALL json % get('element_dimensions', nElements, found)
    CALL json % get('n_timing_repeats', nRepeats, found)
    CALL json % get_core(jCore)
  
    ALLOCATE(f(1:2,1:nFunctions), divF(1:nFunctions))
  
    DO iTest = 1, nFunctions
      ! Point to the i-th scalar_1d function for testing
      CALL jCore % get_child(objPointer, iTest, testPointer, found)
      IF( found )THEN
        ! Pull the function and derivative strings from the JSON
        CALL Get_Char_Obj(jCore, testPointer, 'name', funcName)
        CALL Get_Char_Obj(jCore, testPointer, 'fx', fx)
        CALL Get_Char_Obj(jCore, testPointer, 'fy', fy)
        CALL Get_Char_Obj(jCore, testPointer, 'divf', dF)
        ! Create the exact function equation parsers
        f(1,iTest) = EquationParser(fx,(/'x','y'/))
        f(2,iTest) = EquationParser(fy,(/'x','y'/))
        divF(iTest) = EquationParser(dF,(/'x','y'/))
      ENDIF
    ENDDO

    DO pdeg = polyDeg(1), polyDeg(2)
      CALL VectorGridInterp_2D_Test( f, pdeg, nPlotPoints, nFunctions, nElements, nRepeats )
      CALL VectorDivergence_2D_Test( f, divF, pdeg, nFunctions, nElements, nRepeats )
    ENDDO

    DEALLOCATE( f, divF, polyDeg )

END SUBROUTINE Vector2D_Tests

SUBROUTINE FillVectorValues2D( vector, f, nEl )
  IMPLICIT NONE
  TYPE(Vector2D), INTENT(inout) :: vector
  TYPE(EquationParser), INTENT(in) :: f(1:2,1:vector % nVar)
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

         DO iVar = 1, vector % nVar

           DO j = 0, vector % N
             y = 0.5_prec*( yR*(vector % interp % controlPoints % hostData(j)+1.0_prec) - &
                            yL*(vector % interp % controlPoints % hostData(j)-1.0_prec) )
             DO i = 0, vector % N
               x = 0.5_prec*( xR*(vector % interp % controlPoints % hostData(i)+1.0_prec) - &
                              xL*(vector % interp % controlPoints % hostData(i)-1.0_prec) )

               vector % interior % hostData(1,i,j,iVar,eid) = f(1,iVar) % Evaluate( (/x, y/) )
               vector % interior % hostData(2,i,j,iVar,eid) = f(2,iVar) % Evaluate( (/x, y/) )
             ENDDO
           ENDDO
         ENDDO
       ENDDO
     ENDDO

END SUBROUTINE FillVectorValues2D

SUBROUTINE FillDivergenceValues2D( divergence, divF, nEl )
  IMPLICIT NONE
  TYPE(Scalar2D), INTENT(inout) :: divergence
  TYPE(EquationParser), INTENT(in) :: divF(1:divergence % nVar)
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

         DO iVar = 1, divergence % nVar

           DO j = 0, divergence % N
             y = 0.5_prec*( yR*(divergence % interp % controlPoints % hostData(j)+1.0_prec) - &
                            yL*(divergence % interp % controlPoints % hostData(j)-1.0_prec) )
             DO i = 0, divergence % N
               x = 0.5_prec*( xR*(divergence % interp % controlPoints % hostData(i)+1.0_prec) - &
                              xL*(divergence % interp % controlPoints % hostData(i)-1.0_prec) )

               divergence % interior % hostData(i,j,iVar,eid) = divF(iVar) % Evaluate( (/x, y/) )
             ENDDO
           ENDDO
         ENDDO
       ENDDO
     ENDDO

END SUBROUTINE FillDivergenceValues2D

SUBROUTINE VectorGridInterp_2D_Test( f, N, M, nVar, nElem, nRepeats )
  IMPLICIT NONE
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem
  INTEGER, INTENT(in) :: nRepeats
  TYPE(EquationParser), INTENT(in)  :: f(1:2,1:nVar)
  ! Local
  CHARACTER(LEN=50), DIMENSION(1:nVar) :: functionNames
  CHARACTER(LEN=50) :: cpuModel
  CHARACTER(LEN=50) :: gpuModel
  TYPE(Vector2D) :: fVector, fInterp, fActual , fError
  REAL(prec) :: errors(1:nVar)
  REAL(prec) :: t1, t2, runtime
  INTEGER :: i, iVar

     cpuModel = GetCPUModel_Linux()
     gpuModel = 'None'

     DO iVar = 1, nVar
       functionNames(iVar) = TRIM(f(1,iVar) % equation)//', '//TRIM(f(2,iVar) % equation)
     ENDDO

     CALL fVector % Build( N, GAUSS, M, GAUSS, nVar, nElem*nElem ) 
     CALL fInterp % Build( M, GAUSS, M, GAUSS, nVar, nElem*nElem ) 
     CALL fActual % Build( M, GAUSS, M, GAUSS, nVar, nElem*nElem ) 
     CALL fError % Build( M, GAUSS, M, GAUSS, nVar, nElem*nElem ) 

     CALL FillVectorValues2D( fVector, f, nElem )
     CALL FillVectorValues2D( fActual, f, nElem )

     CALL fVector % GridInterp(fInterp, .FALSE.)

     fError = fActual - fInterp
     errors = fError % AbsMaxInterior( )

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL fVector % GridInterp(fInterp, .FALSE.)
     ENDDO
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec

     CALL Report_JSON(fVector % interp, &
                      'GridInterp_Vector2D', &
                      cpuModel, &
                      gpuModel, &
                      functionNames, &
                      runtime, &
                      errors, &
                      nVar, &
                      nElem*nElem)

#ifdef GPU
     CALL fVector % UpdateDevice( )
     CALL fVector % GridInterp(fInterp, .TRUE.)
     CALL fInterp % UpdateHost( )
     CALL hipCheck(hipDeviceSynchronize()) 

     fError = fActual - fInterp
     errors = fError % AbsMaxInterior( )

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL fVector % GridInterp(fInterp, .TRUE.)
     ENDDO
     CALL hipCheck(hipDeviceSynchronize()) 
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec


     CALL Report_JSON(fVector % interp, &
                      'GridInterp_Vector2D+gpuAccel', &
                      cpuModel, &
                      gpuModel, &
                      functionNames, &
                      runtime, &
                      errors, &
                      nVar, &
                      nElem*nElem)

#endif

     CALL fVector % Trash( )
     CALL fInterp % Trash( )
     CALL fActual % Trash( )
     CALL fError % Trash( )

END SUBROUTINE VectorGridInterp_2D_Test

SUBROUTINE VectorDivergence_2D_Test( f, divF, N, nVar, nElem, nRepeats )
  IMPLICIT NONE
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem
  INTEGER, INTENT(in) :: nRepeats
  TYPE(EquationParser), INTENT(in)  :: f(1:2,1:nVar)
  TYPE(EquationParser), INTENT(in)  :: divF(1:nVar)
  ! Local
  CHARACTER(LEN=50), DIMENSION(1:nVar) :: functionNames
  CHARACTER(LEN=50) :: cpuModel
  CHARACTER(LEN=50) :: gpuModel
  TYPE(Vector2D) :: fVector
  TYPE(Scalar2D) :: divFInterp, divFActual , divFError
  REAL(prec) :: errors(1:nVar)
  REAL(prec) :: t1, t2, runtime
  INTEGER :: i, iVar

     cpuModel = GetCPUModel_Linux()
     gpuModel = 'None'

     DO iVar = 1, nVar
       functionNames(iVar) = TRIM(f(1,iVar) % equation)//', '//TRIM(f(2,iVar) % equation)
     ENDDO

     CALL fVector % Build( N, GAUSS, N, GAUSS, nVar, nElem*nElem ) 
     CALL divFInterp % Build( N, GAUSS, N, GAUSS, nVar, nElem*nElem ) 
     CALL divFActual % Build( N, GAUSS, N, GAUSS, nVar, nElem*nElem ) 
     CALL divFError % Build( N, GAUSS, N, GAUSS, nVar, nElem*nElem ) 

     CALL FillVectorValues2D( fVector, f, nElem )
#ifdef GPU
     CALL fVector % UpdateDevice( )
#endif
     CALL FillDivergenceValues2D( divFActual, divF, nElem )

     CALL fVector % Divergence(divFInterp, .FALSE.)
     divFInterp % interior % hostData = divFInterp % interior % hostData*REAL( nElem, prec ) 

     divFError = divFActual - divFInterp
     errors = divFError % AbsMaxInterior( )

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL fVector % Divergence(divFInterp, .FALSE.)
     ENDDO
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec

     CALL Report_JSON(fVector % interp, &
                      'Divergence_Vector2D', &
                      cpuModel, &
                      gpuModel, &
                      functionNames, &
                      runtime, &
                      errors, &
                      nVar, &
                      nElem*nElem)

#ifdef GPU
     CALL fVector % Divergence(divFInterp, .TRUE.)
     CALL divFInterp % UpdateHost( )
     CALL hipCheck(hipDeviceSynchronize())
     divFInterp % interior % hostData = divFInterp % interior % hostData*REAL( nElem, prec ) 

     divFError = divFActual - divFInterp
     errors = divFError % AbsMaxInterior( )

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL fVector % Divergence(divFInterp, .TRUE.)
     ENDDO
     CALL hipCheck(hipDeviceSynchronize()) 
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec


     CALL Report_JSON(fVector % interp, &
                      'Divergence_Vector2D+gpuAccel', &
                      cpuModel, &
                      gpuModel, &
                      functionNames, &
                      runtime, &
                      errors, &
                      nVar, &
                      nElem*nElem)

#endif

     CALL fVector % Trash( )
     CALL divFInterp % Trash( )
     CALL divFActual % Trash( )
     CALL divFError % Trash( )

END SUBROUTINE VectorDivergence_2D_Test
