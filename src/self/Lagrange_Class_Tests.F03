MODULE Lagrange_Class_Tests

USE ModelPrecision
USE EquationParser_Class
USE Lagrange_Class
USE json_module

USE ISO_C_BINDING
USE hip_fortran

IMPLICIT NONE

   TYPE, EXTENDS(Lagrange) :: Lagrange_Tests

     CONTAINS

     ! GENERIC, PUBLIC :: Derivative_1D => Derivative_1D_cpu, Derivative_1D_gpu
     PROCEDURE :: Scalar_Derivative_1D_Test
     ! GENERIC, PUBLIC :: ScalarGridInterp_1D => ScalarGridInterp_1D_cpu, ScalarGridInterp_1D_gpu
!     PROCEDURE :: Scalar_GridInterp_1D_Test
!     ! GENERIC, PUBLIC :: ScalarBoundaryInterp_1D => ScalarBoundaryInterp_1D_cpu, ScalarBoundaryInterp_1D_gpu
!     PROCEDURE :: Scalar_BoundaryInterp_1D_Test
!
!
!     ! GENERIC, PUBLIC :: ScalarGradient_2D => ScalarGradient_2D_cpu, ScalarGradient_2D_gpu
!     PROCEDURE :: Scalar_Gradient_2D_Test
!     ! GENERIC, PUBLIC :: VectorGradient_2D => VectorGradient_2D_cpu, VectorGradient_2D_gpu
!     PROCEDURE :: Vector_Gradient_2D_Test
!     ! GENERIC, PUBLIC :: VectorDivergence_2D => VectorDivergence_2D_cpu, VectorDivergence_2D_gpu
!     PROCEDURE :: Vector_Divergence_2D_Test
!     ! GENERIC, PUBLIC :: VectorCurl_2D => VectorCurl_2D_cpu, VectorCurl_2D_gpu
!     PROCEDURE :: Vector_Curl_2D_Test
!     ! GENERIC, PUBLIC :: ScalarGridInterp_2D => ScalarGridInterp_2D_cpu, ScalarGridInterp_2D_gpu
!     PROCEDURE :: Scalar_GridInterp_2D_Test
!     ! GENERIC, PUBLIC :: ScalarBoundaryInterp_2D => ScalarBoundaryInterp_2D_cpu, ScalarBoundaryInterp_2D_gpu
!     PROCEDURE :: Scalar_BoundaryInterp_2D_Test
!     ! GENERIC, PUBLIC :: VectorGridInterp_2D => VectorGridInterp_2D_cpu, VectorGridInterp_2D_gpu
!     PROCEDURE :: Vector_GridInterp_2D_Test
!     ! GENERIC, PUBLIC :: VectorBoundaryInterp_2D => VectorBoundaryInterp_2D_cpu, VectorBoundaryInterp_2D_gpu
!     PROCEDURE :: Vector_BoundaryInterp_2D_Test
!     ! GENERIC, PUBLIC :: TensorGridInterp_2D => TensorGridInterp_2D_cpu, TensorGridInterp_2D_gpu
!     PROCEDURE :: Tensor_GridInterp_2D_Test
!     ! GENERIC, PUBLIC :: TensorBoundaryInterp_2D => TensorBoundaryInterp_2D_cpu, TensorBoundaryInterp_2D_gpu
!     PROCEDURE :: Tensor_BoundaryInterp_2D_Test
!
!     ! GENERIC, PUBLIC :: ScalarGradient_3D => ScalarGradient_3D_cpu, ScalarGradient_3D_gpu
!     PROCEDURE :: Scalar_Gradient_3D_Test
!     ! GENERIC, PUBLIC :: VectorGradient_3D => VectorGradient_3D_cpu, VectorGradient_3D_gpu
!     PROCEDURE :: Vector_Gradient_3D_Test
!     ! GENERIC, PUBLIC :: VectorDivergence_3D => VectorDivergence_3D_cpu, VectorDivergence_3D_gpu
!     PROCEDURE :: Vector_Divergence_3D_Test
!     ! GENERIC, PUBLIC :: VectorCurl_3D => VectorCurl_3D_cpu, VectorCurl_3D_gpu
!     PROCEDURE :: Vector_Curl_3D_Test
!     ! GENERIC, PUBLIC :: ScalarGridInterp_3D => ScalarGridInterp_3D_cpu, ScalarGridInterp_3D_gpu
!     PROCEDURE :: Scalar_GridInterp_3D_Test
!     ! GENERIC, PUBLIC :: ScalarBoundaryInterp_3D => ScalarBoundaryInterp_3D_cpu, ScalarBoundaryInterp_3D_gpu
!     PROCEDURE :: Scalar_BoundaryInterp_3D_Test
!     ! GENERIC, PUBLIC :: VectorGridInterp_3D => VectorGridInterp_3D_cpu, VectorGridInterp_3D_gpu
!     PROCEDURE :: Vector_GridInterp_3D_Test
!     ! GENERIC, PUBLIC :: VectorBoundaryInterp_3D => VectorBoundaryInterp_3D_cpu, VectorBoundaryInterp_3D_gpu
!     PROCEDURE :: Vector_BoundaryInterp_3D_Test
!     ! GENERIC, PUBLIC :: TensorGridInterp_3D => TensorGridInterp_3D_cpu, TensorGridInterp_3D_gpu
!     PROCEDURE :: Tensor_GridInterp_3D_Test
!     ! GENERIC, PUBLIC :: TensorBoundaryInterp_3D => TensorBoundaryInterp_3D_cpu, TensorBoundaryInterp_3D_gpu
!     PROCEDURE :: Tensor_BoundaryInterp_3D_Test

     PROCEDURE :: Report_JSON

   END TYPE Lagrange_Tests
  
   TYPE Phrase
     CHARACTER(512) :: val
   END TYPE Phrase

CONTAINS

SUBROUTINE Report_JSON(myPoly, routineName, functions, runtime, errors, nFunctions, nElements)
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE
  CLASS(Lagrange_Tests), INTENT(in) :: myPoly 
  INTEGER, INTENT(in)              :: nFunctions
  INTEGER, INTENT(in)              :: nElements
  CHARACTER(*), INTENT(in)         :: routineName
  CHARACTER(LEN=15),DIMENSION(1:nFunctions), INTENT(in) :: functions
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
    CALL json % initialize()
    CALL json % create_object(p,'')

    CALL json % create_object(conf,'build_conf')
    CALL json % add(p, conf)
    CALL json % add(conf, 'compiler', COMPILER_VERSION())
    CALL json % add(conf, 'compiler_options', COMPILER_OPTIONS())
    CALL json % add(conf, 'precision', prec)

    CALL json % create_object(res,'test_output')
    CALL json % add(p, res)

    CALL json % add(res, 'date_time', TRIM(date)//' '//TRIM(time)//' '//TRIM(zone) )
    CALL json % add(res, 'routine_name', TRIM(routineName))
    CALL json % add(res, 'polynomial_degree', myPoly % N)
    CALL json % add(res, 'n_variables', nFunctions)
    CALL json % add(res, 'n_elements', nElements)
    CALL json % add(res, 'functions', functions)
    CALL json % add(res, 'max_errors', errors)
    CALL json % add(res, 'runtime_ms', runtime)

    nullify(res)
    nullify(conf)

    CALL json % print(p)

    CALL json % destroy(res)
    IF( json % failed() ) STOP 1

END SUBROUTINE Report_JSON

SUBROUTINE Scalar_Derivative_1D_Test(interp, f, dfdx, nFunctions, nElements, nRepeats)
  IMPLICIT NONE
  CLASS(Lagrange_Tests), INTENT(in) :: interp
  TYPE(EquationParser), INTENT(in)  :: f(1:nFunctions)
  TYPE(EquationParser), INTENT(in)  :: dfdx(1:nFunctions)
  INTEGER, INTENT(in) :: nFunctions
  INTEGER, INTENT(in) :: nElements
  INTEGER, INTENT(in) :: nRepeats
  ! Local
  CHARACTER(LEN=15), DIMENSION(1:nFunctions) :: functionNames
  INTEGER :: iEl, iVar, i
  REAL(prec) :: xL, xR, x, dx, t1, t2, runtime
  REAL(prec), POINTER :: fNodal(:,:,:)
  REAL(prec), POINTER :: dfdxInterp(:,:,:)
  REAL(prec), POINTER :: dfdxActual(:,:,:)
  REAL(prec) :: errors(1:nFunctions)
#ifdef GPU
  TYPE(c_ptr) :: fNodal_dev, fInterp_dev, dfdxInterp_dev
#endif


     ALLOCATE(fNodal(0:interp % N, 1:nFunctions, 1:nElements),&
              dfdxInterp(0:interp % N, 1:nFunctions, 1:nElements),&
              dfdxActual(0:interp % N, 1:nFunctions, 1:nElements))
     DO iVar = 1, nFunctions
       functionNames(iVar) = f(iVar) % equation
     ENDDO

     dx = 2.0_prec/REAL(nElements,prec)
     ! Evaluate the exact functions and derivatives at the control points
     DO iEl = 1, nElements
       DO iVar = 1, nFunctions
         DO i = 0, interp % N
           xL = -1.0_prec + REAL((iEl-1),prec)*dx
           xR = xL + dx
           x = 0.5_prec*( xR*(interp % controlPoints(i)+1.0_prec) - xL*(interp % controlPoints(i)-1.0_prec) )
           fNodal(i,iVar,iEl) = f(iVar) % Evaluate( (/ x, 0.0_prec, 0.0_prec /) )          
           dfdxActual(i,iVar,iEl) = dfdx(iVar) % Evaluate( (/ x, 0.0_prec, 0.0_prec /) )
         ENDDO
       ENDDO
     ENDDO

     ! Derivative_1D
     CALL interp % Derivative_1D(fNodal, dfdxInterp, nFunctions, nElements)

     errors = 0.0_prec
     DO iEl = 1, nElements
       DO iVar = 1, nFunctions
         DO i = 0, interp % N
           errors(iVar) = MAX(ABS(dfdxInterp(i,iVar,iEl)*2.0_prec/dx - dfdxActual(i,iVar,iEl)),errors(iVar))
         ENDDO
       ENDDO
     ENDDO

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL interp % Derivative_1D(fNodal, dfdxInterp, nFunctions, nElements)
     ENDDO
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec
     CALL interp % Report_JSON('Derivative_1D', &
                               functionNames, &
                               runtime, &
                               errors, &
                               nFunctions, &
                               nElements)
#ifdef GPU
     CALL hipFortran(hipMalloc(fNodal_dev, SIZEOF(fNodal)))
     CALL hipFortran(hipMalloc(dfdxInterp_dev, SIZEOF(dfdxInterp)))
     CALL hipFortran(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))
     CALL interp % Derivative_1D(fNodal_dev, dfdxInterp_dev, nFunctions, nElements)
     CALL hipFortran(hipMemcpy(c_loc(dfdxInterp), dfdxInterp_dev, SIZEOF(dfdxInterp), hipMemcpyDeviceToHost))

     errors = 0.0_prec
     DO iEl = 1, nElements
       DO iVar = 1, nFunctions
         DO i = 0, interp % N
           errors(iVar) = MAX(ABS(dfdxInterp(i,iVar,iEl)*2.0_prec/dx - dfdxActual(i,iVar,iEl)),errors(iVar))
         ENDDO
       ENDDO
     ENDDO

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL interp % Derivative_1D(fNodal_dev, dfdxInterp_dev, nFunctions, nElements)
     ENDDO
     CALL hipFortran(hipDeviceSynchronize()) 
     CALL CPU_TIME(t2)

     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec
     CALL interp % Report_JSON('Derivative_1D - GPU', &
                               functionNames, &
                               runtime, &
                               errors, &
                               nFunctions, &
                               nElements)

     CALL hipFortran(hipFree(fNodal_dev))
     CALL hipFortran(hipFree(dfdxInterp_dev))
#endif
     DEALLOCATE(fNodal,dfdxInterp, dfdxActual)

END SUBROUTINE Scalar_Derivative_1D_Test
!
!SUBROUTINE Scalar_GridInterp_1D_Test(myPoly, f, nFunctions, nElements )
!  IMPLICIT NONE
!  CLASS(Lagrange_Tests), INTENT(in) :: myPoly
!  TYPE(EquationParser), INTENT(in)  :: f(1:nFunctions)
!  TYPE(EquationParser), INTENT(in)  :: dfdx(1:nFunctions)
!  INTEGER, INTENT(in) :: nFunctions
!  INTEGER, INTENT(in) :: nElements
!  ! Local
!  TYPE(Phrase) :: functionNames(1:nFunctions)
!
!
!     DO iVar = 1, nFunctions
!       ! Pull out the function Names
!     ENDDO
!
!     ! Evaluate the exact functions at the target nodes
!     DO iEl = 1, nElements
!       DO iVar = 1, nFunctions
!         DO i = 0, nPlotPoints
!           xL = -1.0_prec + REAL((iEl-1),prec)*dx
!           xR = xL + dx
!           x = 0.5_prec*( xR*(myPoly % targetNodes(i)+1.0_prec) - xL*(myPoly % targetNodes(i)-1.0_prec) )
!           fActual(i,iVar,iEl) = f(iVar) % Evaluate( (/ x, 0.0_prec, 0.0_prec /) )          
!         ENDDO
!       ENDDO
!     ENDDO
!
!     ! Evaluate the exact functions and derivatives at the control points
!     DO iEl = 1, nElements
!       DO iVar = 1, nFunctions
!         DO i = 0, pdeg
!           xL = -1.0_prec + REAL((iEl-1),prec)*dx
!           xR = xL + dx
!           x = 0.5_prec*( xR*(interp % controlPoints(i)+1.0_prec) - xL*(interp % controlPoints(i)-1.0_prec) )
!           fNodal(i,iVar,iEl) = f(iVar) % Evaluate( (/ x, 0.0_prec, 0.0_prec /) )          
!         ENDDO
!       ENDDO
!     ENDDO
!
!     CALL interp % ScalarGridInterp_1D(fNodal, fInterp, nFunctions, nElements)  
!
!     errors = 0.0_prec
!     DO iEl = 1, nElements
!       DO iVar = 1, nFunctions
!         DO i = 0, nPlotPoints
!           errors(iVar) = MAX(ABS(fInterp(i,iVar,iEl) -fActual(i,iVar,iEl)),fErr)
!         ENDDO
!       ENDDO
!     ENDDO
!
!     ! Estimate routine runtime
!     CALL CPU_TIME(t1)
!     DO i = 1, nRepeats
!       CALL interp % ScalarGridInterp_1D(fNodal, fInterp, nFunctions, nElements)
!     ENDDO
!     CALL CPU_TIME(t2)
!
!     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec
! 
!     CALL myPoly % Report_JSON('Scalar_GridInterp_1D', &
!                               functionNames, &
!                               runtime, &
!                               errors, &
!                               nFunctions, &
!                               nElements)
!
!#ifdef GPU
!     CALL hfMalloc(fNodal_dev, SIZEOF(fNodal))
!     CALL hipFortran(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))
!
!     CALL interp % ScalarGridInterp_1D(fNodal_dev, fInterp_dev, nFunctions, nElements)  
!
!     CALL hipFortran(hipMemcpy(c_loc(fInterp), fInterp_dev, SIZEOF(fInterp), hipMemcpyDeviceToHost))
!
!     fErr = 0.0_prec
!     DO iEl = 1, nElements
!       DO iVar = 1, nFunctions
!         DO i = 0, nPlotPoints
!           fErr(iVar) = MAX(ABS(fInterp(i,iVar,iEl) -fActual(i,iVar,iEl)),fErr)
!         ENDDO
!       ENDDO
!     ENDDO
!
!     ! Estimate routine runtime
!     CALL CPU_TIME(t1)
!     DO i = 1, nRepeats
!       CALL interp % ScalarGridInterp_1D(fNodal_dev, fInterp_dev, nFunctions, nElements)
!     ENDDO
!     CALL hipFortran(hipDeviceSynchronize()) 
!     CALL CPU_TIME(t2)
!
!     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec
! 
!     CALL myPoly % Report_JSON('Scalar_GridInterp_1D - GPU', &
!                               functionNames, &
!                               runtime, &
!                               errors, &
!                               nFunctions, &
!                               nElements)
!#endif
!
!END SUBROUTINE Scalar_GridInterp_1D_Test

END MODULE Lagrange_Class_Tests
