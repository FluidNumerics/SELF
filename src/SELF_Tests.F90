! SELF_Tests.F90
! 
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self-fluids@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE SELF_Tests

USE SysConf
USE FEQParse
USE json_module

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Data
USE SELF_Mesh
USE SELF_MappedData


USE ISO_C_BINDING
USE hipfort

IMPLICIT NONE


     ! GENERIC, PUBLIC :: Derivative_1D => Derivative_1D_cpu, Derivative_1D_gpu
     PUBLIC :: ScalarDerivative_1D_Test
     ! GENERIC, PUBLIC :: ScalarGridInterp_1D => ScalarGridInterp_1D_cpu, ScalarGridInterp_1D_gpu
!     PUBLIC :: ScalarGridInterp_1D_Test
!    ! GENERIC, PUBLIC :: ScalarBoundaryInterp_1D => ScalarBoundaryInterp_1D_cpu, ScalarBoundaryInterp_1D_gpu
!     PUBLIC :: Scalar_BoundaryInterp_1D_Test
!
!
!     ! GENERIC, PUBLIC :: ScalarGradient_2D => ScalarGradient_2D_cpu, ScalarGradient_2D_gpu
!     PUBLIC :: ScalarGradient_2D_Test
!     ! GENERIC, PUBLIC :: VectorGradient_2D => VectorGradient_2D_cpu, VectorGradient_2D_gpu
!     PUBLIC :: Vector_Gradient_2D_Test
!     ! GENERIC, PUBLIC :: VectorDivergence_2D => VectorDivergence_2D_cpu, VectorDivergence_2D_gpu
!     PUBLIC :: Vector_Divergence_2D_Test
!     ! GENERIC, PUBLIC :: VectorCurl_2D => VectorCurl_2D_cpu, VectorCurl_2D_gpu
!     PUBLIC :: Vector_Curl_2D_Test
!     ! GENERIC, PUBLIC :: ScalarGridInterp_2D => ScalarGridInterp_2D_cpu, ScalarGridInterp_2D_gpu
!     PUBLIC :: ScalarGridInterp_2D_Test
!     ! GENERIC, PUBLIC :: ScalarBoundaryInterp_2D => ScalarBoundaryInterp_2D_cpu, ScalarBoundaryInterp_2D_gpu
!     PUBLIC :: Scalar_BoundaryInterp_2D_Test
!     ! GENERIC, PUBLIC :: VectorGridInterp_2D => VectorGridInterp_2D_cpu, VectorGridInterp_2D_gpu
!     PUBLIC :: Vector_GridInterp_2D_Test
!     ! GENERIC, PUBLIC :: VectorBoundaryInterp_2D => VectorBoundaryInterp_2D_cpu, VectorBoundaryInterp_2D_gpu
!     PUBLIC :: Vector_BoundaryInterp_2D_Test
!     ! GENERIC, PUBLIC :: TensorGridInterp_2D => TensorGridInterp_2D_cpu, TensorGridInterp_2D_gpu
!     PUBLIC :: Tensor_GridInterp_2D_Test
!     ! GENERIC, PUBLIC :: TensorBoundaryInterp_2D => TensorBoundaryInterp_2D_cpu, TensorBoundaryInterp_2D_gpu
!     PUBLIC :: Tensor_BoundaryInterp_2D_Test
!
!     ! GENERIC, PUBLIC :: ScalarGradient_3D => ScalarGradient_3D_cpu, ScalarGradient_3D_gpu
!     PUBLIC :: Scalar_Gradient_3D_Test
!     ! GENERIC, PUBLIC :: VectorGradient_3D => VectorGradient_3D_cpu, VectorGradient_3D_gpu
!     PUBLIC :: Vector_Gradient_3D_Test
!     ! GENERIC, PUBLIC :: VectorDivergence_3D => VectorDivergence_3D_cpu, VectorDivergence_3D_gpu
!     PUBLIC :: Vector_Divergence_3D_Test
!     ! GENERIC, PUBLIC :: VectorCurl_3D => VectorCurl_3D_cpu, VectorCurl_3D_gpu
!     PUBLIC :: Vector_Curl_3D_Test
!     ! GENERIC, PUBLIC :: ScalarGridInterp_3D => ScalarGridInterp_3D_cpu, ScalarGridInterp_3D_gpu
!     PUBLIC :: Scalar_GridInterp_3D_Test
!     ! GENERIC, PUBLIC :: ScalarBoundaryInterp_3D => ScalarBoundaryInterp_3D_cpu, ScalarBoundaryInterp_3D_gpu
!     PUBLIC :: Scalar_BoundaryInterp_3D_Test
!     ! GENERIC, PUBLIC :: VectorGridInterp_3D => VectorGridInterp_3D_cpu, VectorGridInterp_3D_gpu
!     PUBLIC :: Vector_GridInterp_3D_Test
!     ! GENERIC, PUBLIC :: VectorBoundaryInterp_3D => VectorBoundaryInterp_3D_cpu, VectorBoundaryInterp_3D_gpu
!     PUBLIC :: Vector_BoundaryInterp_3D_Test
!     ! GENERIC, PUBLIC :: TensorGridInterp_3D => TensorGridInterp_3D_cpu, TensorGridInterp_3D_gpu
!     PUBLIC :: Tensor_GridInterp_3D_Test
!     ! GENERIC, PUBLIC :: TensorBoundaryInterp_3D => TensorBoundaryInterp_3D_cpu, TensorBoundaryInterp_3D_gpu
!     PUBLIC :: Tensor_BoundaryInterp_3D_Test

     PRIVATE :: Report_JSON

CONTAINS


 
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
    nullify(p)

END SUBROUTINE Get_Char_Obj

SUBROUTINE Report_JSON(myPoly, routineName, cpu, gpu, functions, runtime, errors, nFunctions, nElements)
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE
  TYPE(Lagrange), INTENT(in) :: myPoly 
  INTEGER, INTENT(in)              :: nFunctions
  INTEGER, INTENT(in)              :: nElements
  CHARACTER(*), INTENT(in)         :: cpu
  CHARACTER(*), INTENT(in)         :: gpu
  CHARACTER(*), INTENT(in)         :: routineName
  CHARACTER(LEN=50),DIMENSION(1:nFunctions), INTENT(in) :: functions
  REAL(prec), INTENT(in)           :: runtime
  REAL(prec), INTENT(in)           :: errors(1:nFunctions)
  ! Local
  TYPE(JSON_VALUE), POINTER :: p, res, conf, fArray
  TYPE(JSON_CORE) :: json
  CHARACTER(10) :: date
  CHARACTER(8) :: time
  CHARACTER(5) :: zone
  INTEGER :: i

    CALL DATE_AND_TIME(date, time, zone)
    DO i = 1, nFunctions
      CALL json % initialize()
      CALL json % create_object(p,'')

      CALL json % add(p, 'fortran_compiler', COMPILER_VERSION())
      CALL json % add(p, 'fortran_compiler_options', COMPILER_OPTIONS())
      CALL json % add(p, 'c_compiler', 'Unknown')
      CALL json % add(p, 'c_compiler_options', 'Unknown')
      CALL json % add(p, 'precision', prec)
      CALL json % add(p, 'cpu', TRIM(cpu))
      CALL json % add(p, 'gpu', TRIM(gpu))

      CALL json % add(p, 'date_time', TRIM(date(1:4))//'-'//TRIM(date(5:6))//'-'//TRIM(date(7:8))//'T'//&
                                        TRIM(time(1:2))//':'//TRIM(time(3:4))//':'//TRIM(time(5:6))//TRIM(zone) )
      CALL json % add(p, 'routine_name', TRIM(routineName))
      CALL json % add(p, 'polynomial_degree', myPoly % N)
      CALL json % add(p, 'n_function_batch', nFunctions)
      CALL json % add(p, 'n_elements', nElements)
      CALL json % add(p, 'functions', TRIM(functions(i)))
      CALL json % add(p, 'max_errors', errors(i))
      CALL json % add(p, 'runtime_ms', runtime)

      nullify(res)
      nullify(conf)

      CALL json % print(p)

      CALL json % destroy(res)
      IF( json % failed() ) STOP 1
    ENDDO

END SUBROUTINE Report_JSON

SUBROUTINE Scalar1DTests( )
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
    ENDDO

    DEALLOCATE( f, dfdx, polyDeg )

END SUBROUTINE Scalar1DTests

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
!
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

     CALL fScalar % Build( N, GAUSS, M, UNIFORM, nVar, nElem ) 
     CALL dfdxInterp % Build( N, GAUSS, M, UNIFORM, nVar, nElem ) 
     CALL dfdxActual % Build( N, GAUSS, M, UNIFORM, nVar, nElem ) 
     CALL dfdxError % Build( N, GAUSS, M, UNIFORM, nVar, nElem ) 

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

END MODULE SELF_Tests
