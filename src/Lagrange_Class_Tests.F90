! Lagrange_Class_Tests.F90
! 
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self-fluids@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE Lagrange_Class_Tests

USE ModelPrecision
USE SysConf
USE Lagrange_Class
USE FEQParse
USE json_module

USE ISO_C_BINDING
USE hip

IMPLICIT NONE

   TYPE, EXTENDS(Lagrange) :: Lagrange_Tests

     CONTAINS

     ! GENERIC, PUBLIC :: Derivative_1D => Derivative_1D_cpu, Derivative_1D_gpu
     PROCEDURE :: ScalarDerivative_1D_Test
     ! GENERIC, PUBLIC :: ScalarGridInterp_1D => ScalarGridInterp_1D_cpu, ScalarGridInterp_1D_gpu
     PROCEDURE :: ScalarGridInterp_1D_Test
!    ! GENERIC, PUBLIC :: ScalarBoundaryInterp_1D => ScalarBoundaryInterp_1D_cpu, ScalarBoundaryInterp_1D_gpu
!     PROCEDURE :: Scalar_BoundaryInterp_1D_Test
!
!
!     ! GENERIC, PUBLIC :: ScalarGradient_2D => ScalarGradient_2D_cpu, ScalarGradient_2D_gpu
     PROCEDURE :: ScalarGradient_2D_Test
!     ! GENERIC, PUBLIC :: VectorGradient_2D => VectorGradient_2D_cpu, VectorGradient_2D_gpu
!     PROCEDURE :: Vector_Gradient_2D_Test
!     ! GENERIC, PUBLIC :: VectorDivergence_2D => VectorDivergence_2D_cpu, VectorDivergence_2D_gpu
!     PROCEDURE :: Vector_Divergence_2D_Test
!     ! GENERIC, PUBLIC :: VectorCurl_2D => VectorCurl_2D_cpu, VectorCurl_2D_gpu
!     PROCEDURE :: Vector_Curl_2D_Test
!     ! GENERIC, PUBLIC :: ScalarGridInterp_2D => ScalarGridInterp_2D_cpu, ScalarGridInterp_2D_gpu
     PROCEDURE :: ScalarGridInterp_2D_Test
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

SUBROUTINE Report_JSON(myPoly, routineName, cpu, gpu, functions, runtime, errors, nFunctions, nElements)
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE
  CLASS(Lagrange_Tests), INTENT(in) :: myPoly 
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

FUNCTION FillScalarValues1D(f, cpoints, N, nFunctions, nElem) RESULT(fNodal)
  IMPLICIT NONE
  TYPE(EquationParser) :: f(1:nFunctions)
  INTEGER :: N, nFunctions, nElem
  REAL(prec) :: fNodal(0:N,1:nFunctions,1:nElem)
  REAL(prec) :: cpoints(0:N) 
  ! Local
  INTEGER :: iEl, iVar, i
  REAL(prec) :: xL, xR, x, dx

     dx = 2.0_prec/REAL(nElem,prec)
     DO iEl = 1, nElem
       DO iVar = 1, nFunctions
         DO i = 0, N
           xL = -1.0_prec + REAL((iEl-1),prec)*dx
           xR = xL + dx
           x = 0.5_prec*( xR*(cpoints(i)+1.0_prec) - xL*(cpoints(i)-1.0_prec) )
           fNodal(i,iVar,iEl) = f(iVar) % Evaluate( (/x/) )          
         ENDDO
       ENDDO
     ENDDO

END FUNCTION FillScalarValues1D

FUNCTION GetMaxErrorsScalar1D(fEst, fAct, N, nFunctions, nElem) RESULT(errors)
  IMPLICIT NONE
  INTEGER :: N, nFunctions, nElem
  REAL(prec) :: fEst(0:N,1:nFunctions,1:nElem)
  REAL(prec) :: fAct(0:N,1:nFunctions,1:nElem)
  REAL(prec) :: errors(1:nFunctions)
  ! Local
  INTEGER :: iEl, iVar, i

     errors = 0.0_prec
     DO iEl = 1, nElem
       DO iVar = 1, nFunctions
         DO i = 0, N
           errors(iVar) = MAX(ABS(fEst(i,iVar,iEl) - fAct(i,iVar,iEl)),errors(iVar))
         ENDDO
       ENDDO
     ENDDO

END FUNCTION GetMaxErrorsScalar1D

SUBROUTINE ScalarDerivative_1D_Test(interp, f, dfdx, nFunctions, nElements, nRepeats)
  IMPLICIT NONE
  CLASS(Lagrange_Tests), INTENT(in) :: interp
  TYPE(EquationParser), INTENT(in)  :: f(1:nFunctions)
  TYPE(EquationParser), INTENT(in)  :: dfdx(1:nFunctions)
  INTEGER, INTENT(in) :: nFunctions
  INTEGER, INTENT(in) :: nElements
  INTEGER, INTENT(in) :: nRepeats
  ! Local
  CHARACTER(LEN=50), DIMENSION(1:nFunctions) :: functionNames
  CHARACTER(LEN=50) :: cpuModel
  CHARACTER(LEN=50) :: gpuModel
  INTEGER :: iEl, iVar, i
  REAL(prec) :: t1, t2, runtime
  REAL(prec), POINTER :: fNodal(:,:,:)
  REAL(prec), POINTER :: dfdxInterp(:,:,:)
  REAL(prec), POINTER :: dfdxActual(:,:,:)
  REAL(prec) :: errors(1:nFunctions)
#ifdef GPU
  TYPE(c_ptr) :: fNodal_dev, fInterp_dev, dfdxInterp_dev
#endif

     cpuModel = GetCPUModel_Linux()
     gpuModel = 'None'

     ALLOCATE(fNodal(0:interp % N, 1:nFunctions, 1:nElements),&
              dfdxInterp(0:interp % N, 1:nFunctions, 1:nElements),&
              dfdxActual(0:interp % N, 1:nFunctions, 1:nElements))
     DO iVar = 1, nFunctions
       functionNames(iVar) = f(iVar) % equation
     ENDDO

     fNodal = FillScalarValues1D(f, interp % controlPoints, interp % N, nFunctions, nElements)
     dfdxActual = FillScalarValues1D(dfdx, interp % controlPoints, interp % N, nFunctions, nElements)
     
     ! Derivative_1D
     CALL interp % Derivative_1D(fNodal, dfdxInterp, nFunctions, nElements)
     ! Transform from computational coordinates to physical
     dfdxInterp = dfdxInterp*REAL(nElements,prec)

     errors = GetMaxErrorsScalar1D(dfdxInterp, dfdxActual, interp % N, nFunctions, nElements)

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL interp % Derivative_1D(fNodal, dfdxInterp, nFunctions, nElements)
     ENDDO
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec

     CALL interp % Report_JSON('Derivative_1D', &
                               cpuModel, &
                               gpuModel, &
                               functionNames, &
                               runtime, &
                               errors, &
                               nFunctions, &
                               nElements)
#ifdef GPU
     CALL hipCheck(hipMalloc(fNodal_dev, SIZEOF(fNodal)))
     CALL hipCheck(hipMalloc(dfdxInterp_dev, SIZEOF(dfdxInterp)))
     CALL hipCheck(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))
     CALL interp % Derivative_1D(fNodal_dev, dfdxInterp_dev, nFunctions, nElements)
     CALL hipCheck(hipMemcpy(c_loc(dfdxInterp), dfdxInterp_dev, SIZEOF(dfdxInterp), hipMemcpyDeviceToHost))
     ! Transform from computational coordinates to physical
     dfdxInterp = dfdxInterp*REAL(nElements,prec)

     errors = GetMaxErrorsScalar1D(dfdxInterp, dfdxActual, interp % N, nFunctions, nElements)

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL interp % Derivative_1D(fNodal_dev, dfdxInterp_dev, nFunctions, nElements)
     ENDDO
     CALL hipCheck(hipDeviceSynchronize()) 
     CALL CPU_TIME(t2)

     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec
     CALL interp % Report_JSON('Derivative_1D', &
                               cpuModel, &
                               gpuModel, &
                               functionNames, &
                               runtime, &
                               errors, &
                               nFunctions, &
                               nElements)

     CALL hipCheck(hipFree(fNodal_dev))
     CALL hipCheck(hipFree(dfdxInterp_dev))
#endif
     DEALLOCATE(fNodal,dfdxInterp, dfdxActual)

END SUBROUTINE ScalarDerivative_1D_Test
!
SUBROUTINE ScalarGridInterp_1D_Test(interp, f, nFunctions, nElements, nRepeats)
  IMPLICIT NONE
  CLASS(Lagrange_Tests), INTENT(in) :: interp
  TYPE(EquationParser), INTENT(in)  :: f(1:nFunctions)
  INTEGER, INTENT(in) :: nFunctions
  INTEGER, INTENT(in) :: nElements
  INTEGER, INTENT(in) :: nRepeats
  ! Local
  CHARACTER(LEN=50), DIMENSION(1:nFunctions) :: functionNames
  CHARACTER(LEN=50) :: cpuModel
  CHARACTER(LEN=50) :: gpuModel
  INTEGER :: iEl, iVar, i
  REAL(prec) :: xL, xR, x, dx, t1, t2, runtime
  REAL(prec), POINTER :: fNodal(:,:,:)
  REAL(prec), POINTER :: fInterp(:,:,:)
  REAL(prec), POINTER :: fActual(:,:,:)
  REAL(prec) :: errors(1:nFunctions)
#ifdef GPU
  TYPE(c_ptr) :: fNodal_dev, fInterp_dev
#endif

     cpuModel = GetCPUModel_Linux()
     gpuModel = 'None'

     ALLOCATE(fNodal(0:interp % N, 1:nFunctions, 1:nElements),&
              fInterp(0:interp % M, 1:nFunctions, 1:nElements),&
              fActual(0:interp % M, 1:nFunctions, 1:nElements))
     DO iVar = 1, nFunctions
       functionNames(iVar) = f(iVar) % equation
     ENDDO

     dx = 2.0_prec/REAL(nElements,prec)
     ! Evaluate the exact functions and derivatives at the control points
     DO iEl = 1, nElements
       xL = -1.0_prec + REAL((iEl-1),prec)*dx
       xR = xL + dx
       DO iVar = 1, nFunctions
         DO i = 0, interp % N
           x = 0.5_prec*( xR*(interp % controlPoints(i)+1.0_prec) - xL*(interp % controlPoints(i)-1.0_prec) )
           fNodal(i,iVar,iEl) = f(iVar) % Evaluate( (/ x, 0.0_prec, 0.0_prec /) )          
         ENDDO
         DO i = 0, interp % M
           x = 0.5_prec*( xR*(interp % targetPoints(i)+1.0_prec) - xL*(interp % targetPoints(i)-1.0_prec) )
           fActual(i,iVar,iEl) = f(iVar) % Evaluate( (/ x, 0.0_prec, 0.0_prec /) )          
         ENDDO
       ENDDO
     ENDDO

     CALL interp % ScalarGridInterp_1D(fNodal, fInterp, nFunctions, nElements)

     errors = 0.0_prec
     DO iEl = 1, nElements
       DO iVar = 1, nFunctions
         DO i = 0, interp % M
           errors(iVar) = MAX(ABS(fInterp(i,iVar,iEl) - fActual(i,iVar,iEl)),errors(iVar))
         ENDDO
       ENDDO
     ENDDO

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL interp % ScalarGridInterp_1D(fNodal, fInterp, nFunctions, nElements)
     ENDDO
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec
     CALL interp % Report_JSON('ScalarGridInterp_1D', &
                               cpuModel, &
                               gpuModel, &
                               functionNames, &
                               runtime, &
                               errors, &
                               nFunctions, &
                               nElements)
#ifdef GPU
     CALL hipCheck(hipMalloc(fNodal_dev, SIZEOF(fNodal)))
     CALL hipCheck(hipMalloc(fInterp_dev, SIZEOF(fInterp)))
     CALL hipCheck(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))
     CALL interp % ScalarGridInterp_1D(fNodal_dev, fInterp_dev, nFunctions, nElements)
     CALL hipCheck(hipMemcpy(c_loc(fInterp), fInterp_dev, SIZEOF(fInterp), hipMemcpyDeviceToHost))

     errors = 0.0_prec
     DO iEl = 1, nElements
       DO iVar = 1, nFunctions
         DO i = 0, interp % M
           errors(iVar) = MAX(ABS(fInterp(i,iVar,iEl) - fActual(i,iVar,iEl)),errors(iVar))
         ENDDO
       ENDDO
     ENDDO

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL interp % ScalarGridInterp_1D(fNodal_dev, fInterp_dev, nFunctions, nElements)
     ENDDO
     CALL hipCheck(hipDeviceSynchronize()) 
     CALL CPU_TIME(t2)

     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec
     CALL interp % Report_JSON('ScalarGridInterp_1D', &
                               cpuModel, &
                               gpuModel, &
                               functionNames, &
                               runtime, &
                               errors, &
                               nFunctions, &
                               nElements)

     CALL hipCheck(hipFree(fNodal_dev))
     CALL hipCheck(hipFree(fInterp_dev))
#endif
     DEALLOCATE(fNodal,fInterp, fActual)

END SUBROUTINE ScalarGridInterp_1D_Test

SUBROUTINE ScalarGradient_2D_Test(interp, f, dfdx, dfdy, nFunctions, nElements, nRepeats)
  IMPLICIT NONE
  CLASS(Lagrange_Tests), INTENT(in) :: interp
  TYPE(EquationParser), INTENT(in)  :: f(1:nFunctions)
  TYPE(EquationParser), INTENT(in)  :: dfdy(1:nFunctions)
  TYPE(EquationParser), INTENT(in)  :: dfdx(1:nFunctions)
  INTEGER, INTENT(in) :: nFunctions
  INTEGER, INTENT(in) :: nElements
  INTEGER, INTENT(in) :: nRepeats
  ! Local
  CHARACTER(LEN=50), DIMENSION(1:nFunctions) :: functionNames
  CHARACTER(LEN=50) :: cpuModel
  CHARACTER(LEN=50) :: gpuModel
  INTEGER :: iX, iY, iEl, iVar, i, j
  INTEGER :: nElem
  REAL(prec) :: xL, xR, yS, yN, x, y, dx, t1, t2, runtime
  REAL(prec), POINTER :: fNodal(:,:,:,:)
  REAL(prec), POINTER :: gradFActual(:,:,:,:,:)
  REAL(prec), POINTER :: gradF(:,:,:,:,:)
  REAL(prec) :: errors(1:nFunctions)
#ifdef GPU
  TYPE(c_ptr) :: fNodal_dev, gradF_dev
#endif

     cpuModel = GetCPUModel_Linux()
     gpuModel = 'None'

     nElem = nElements*nElements
     ALLOCATE(fNodal(0:interp % N, 0:interp % N, 1:nFunctions, 1:nElem),&
              gradFActual(1:2,0:interp % N, 0:interp % N, 1:nFunctions, 1:nElem),&
              gradF(1:2,0:interp % N, 0:interp % N, 1:nFunctions, 1:nElem))
     DO iVar = 1, nFunctions
       functionNames(iVar) = f(iVar) % equation
     ENDDO

     ! Evaluate the exact functions and derivatives at the control points
     dx = 2.0_prec/REAL(nElements,prec)
     DO iY = 1, nElements
       yS = -1.0_prec + REAL((iY-1),prec)*dx
       yN = yS + dx
       DO iX = 1, nElements
         xL = -1.0_prec + REAL((iX-1),prec)*dx
         xR = xL + dx
         iEl = iX + (iY-1)*nElements
         DO iVar = 1, nFunctions

           DO j = 0, interp % N
             y = 0.5_prec*( yN*(interp % controlPoints(j)+1.0_prec) - yS*(interp % controlPoints(j)-1.0_prec) )
             DO i = 0, interp % N
               x = 0.5_prec*( xR*(interp % controlPoints(i)+1.0_prec) - xL*(interp % controlPoints(i)-1.0_prec) )
               fNodal(i,j,iVar,iEl) = f(iVar) % Evaluate( (/ x, y, 0.0_prec /) )          
               gradFActual(1,i,j,iVar,iEl) = dfdx(iVar) % Evaluate( (/ x, y, 0.0_prec /) )          
               gradFActual(2,i,j,iVar,iEl) = dfdy(iVar) % Evaluate( (/ x, y, 0.0_prec /) )          
             ENDDO
           ENDDO

         ENDDO
       ENDDO
     ENDDO

     CALL interp % ScalarGradient_2D(fNodal, gradF, nFunctions, nElem)

     errors = 0.0_prec
     DO iEl = 1, nElem
       DO iVar = 1, nFunctions
         DO j = 0, interp % N
           DO i = 0, interp % N
             errors(iVar) = MAX(ABS(gradF(1,i,j,iVar,iEl)*2.0_prec/dx - gradFActual(1,i,j,iVar,iEl)),errors(iVar))
             errors(iVar) = MAX(ABS(gradF(2,i,j,iVar,iEl)*2.0_prec/dx - gradFActual(2,i,j,iVar,iEl)),errors(iVar))
           ENDDO
         ENDDO
       ENDDO
     ENDDO

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL interp % ScalarGradient_2D(fNodal, gradF, nFunctions, nElem)
     ENDDO
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec
     CALL interp % Report_JSON('ScalarGradient_2D', &
                               cpuModel, &
                               gpuModel, &
                               functionNames, &
                               runtime, &
                               errors, &
                               nFunctions, &
                               nElem)
#ifdef GPU
     CALL hipCheck(hipMalloc(fNodal_dev, SIZEOF(fNodal)))
     CALL hipCheck(hipMalloc(gradF_dev, SIZEOF(gradF)))
     CALL hipCheck(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))
     CALL interp % ScalarGradient_2D(fNodal_dev, gradF_dev, nFunctions, nElem)
     CALL hipCheck(hipMemcpy(c_loc(gradF), gradF_dev, SIZEOF(gradF), hipMemcpyDeviceToHost))

     errors = 0.0_prec
     DO iEl = 1, nElem
       DO iVar = 1, nFunctions
         DO j = 0, interp % N
           DO i = 0, interp % N
             errors(iVar) = MAX(ABS(gradF(1,i,j,iVar,iEl)*2.0_prec/dx - gradFActual(1,i,j,iVar,iEl)),errors(iVar))
             errors(iVar) = MAX(ABS(gradF(2,i,j,iVar,iEl)*2.0_prec/dx - gradFActual(2,i,j,iVar,iEl)),errors(iVar))
           ENDDO
         ENDDO
       ENDDO
     ENDDO

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL interp % ScalarGradient_2D(fNodal_dev, gradF_dev, nFunctions, nElem)
     ENDDO
     CALL hipCheck(hipDeviceSynchronize()) 
     CALL CPU_TIME(t2)

     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec
     CALL interp % Report_JSON('ScalarGradient_2D', &
                               cpuModel, &
                               gpuModel, &
                               functionNames, &
                               runtime, &
                               errors, &
                               nFunctions, &
                               nElem)

     CALL hipCheck(hipFree(fNodal_dev))
     CALL hipCheck(hipFree(gradF_dev))
#endif
     DEALLOCATE(fNodal,gradF,gradFActual)

END SUBROUTINE ScalarGradient_2D_Test

SUBROUTINE ScalarGridInterp_2D_Test(interp, f, nFunctions, nElements, nRepeats)
  IMPLICIT NONE
  CLASS(Lagrange_Tests), INTENT(in) :: interp
  TYPE(EquationParser), INTENT(in)  :: f(1:nFunctions)
  INTEGER, INTENT(in) :: nFunctions
  INTEGER, INTENT(in) :: nElements
  INTEGER, INTENT(in) :: nRepeats
  ! Local
  CHARACTER(LEN=50), DIMENSION(1:nFunctions) :: functionNames
  CHARACTER(LEN=50) :: cpuModel
  CHARACTER(LEN=50) :: gpuModel
  INTEGER :: iX, iY, iEl, iVar, i, j
  INTEGER :: nElem
  REAL(prec) :: xL, xR, yS, yN, x, y, dx, t1, t2, runtime
  REAL(prec), POINTER :: fNodal(:,:,:,:)
  REAL(prec), POINTER :: fInterp(:,:,:,:)
  REAL(prec), POINTER :: fActual(:,:,:,:)
  REAL(prec) :: errors(1:nFunctions)
#ifdef GPU
  TYPE(c_ptr) :: fNodal_dev, fInterp_dev
#endif

     cpuModel = GetCPUModel_Linux()
     gpuModel = 'None'
     nElem = nElements*nElements
     ALLOCATE(fNodal(0:interp % N, 0:interp % N, 1:nFunctions, 1:nElem),&
              fInterp(0:interp % M, 0:interp % M, 1:nFunctions, 1:nElem),&
              fActual(0:interp % M, 0:interp % M, 1:nFunctions, 1:nElem))
     DO iVar = 1, nFunctions
       functionNames(iVar) = f(iVar) % equation
     ENDDO

     ! Evaluate the exact functions at the control points and target points
     dx = 2.0_prec/REAL(nElements,prec)
     DO iY = 1, nElements
       yS = -1.0_prec + REAL((iY-1),prec)*dx
       yN = yS + dx
       DO iX = 1, nElements
         xL = -1.0_prec + REAL((iX-1),prec)*dx
         xR = xL + dx
         iEl = iX + (iY-1)*nElements
         DO iVar = 1, nFunctions

           ! Control point evaluation
           DO j = 0, interp % N
             y = 0.5_prec*( yN*(interp % controlPoints(j)+1.0_prec) - yS*(interp % controlPoints(j)-1.0_prec) )
             DO i = 0, interp % N
               x = 0.5_prec*( xR*(interp % controlPoints(i)+1.0_prec) - xL*(interp % controlPoints(i)-1.0_prec) )
               fNodal(i,j,iVar,iEl) = f(iVar) % Evaluate( (/ x, y, 0.0_prec /) )          
             ENDDO
           ENDDO

           ! Target point evaluation
           DO j = 0, interp % M
             y = 0.5_prec*( yN*(interp % targetPoints(j)+1.0_prec) - yS*(interp % targetPoints(j)-1.0_prec) )
             DO i = 0, interp % M
               x = 0.5_prec*( xR*(interp % targetPoints(i)+1.0_prec) - xL*(interp % targetPoints(i)-1.0_prec) )
               fActual(i,j,iVar,iEl) = f(iVar) % Evaluate( (/ x, y, 0.0_prec /) )          
             ENDDO
           ENDDO

         ENDDO
       ENDDO
     ENDDO

     CALL interp % ScalarGridInterp_2D(fNodal, fInterp, nFunctions, nElem)

     errors = 0.0_prec
     DO iEl = 1, nElem
       DO iVar = 1, nFunctions
         DO j = 0, interp % M
           DO i = 0, interp % M
             errors(iVar) = MAX(ABS(fInterp(i,j,iVar,iEl) - fActual(i,j,iVar,iEl)),errors(iVar))
           ENDDO
         ENDDO
       ENDDO
     ENDDO

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL interp % ScalarGridInterp_2D(fNodal, fInterp, nFunctions, nElem)
     ENDDO
     CALL CPU_TIME(t2)
     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec
     CALL interp % Report_JSON('ScalarGridInterp_2D', &
                               cpuModel, &
                               gpuModel, &
                               functionNames, &
                               runtime, &
                               errors, &
                               nFunctions, &
                               nElem)
#ifdef GPU
     CALL hipCheck(hipMalloc(fNodal_dev, SIZEOF(fNodal)))
     CALL hipCheck(hipMalloc(fInterp_dev, SIZEOF(fInterp)))
     CALL hipCheck(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))
     CALL interp % ScalarGridInterp_2D(fNodal_dev, fInterp_dev, nFunctions, nElem)
     CALL hipCheck(hipMemcpy(c_loc(fInterp), fInterp_dev, SIZEOF(fInterp), hipMemcpyDeviceToHost))

     errors = 0.0_prec
     DO iEl = 1, nElem
       DO iVar = 1, nFunctions
         DO j = 0, interp % M
           DO i = 0, interp % M
             errors(iVar) = MAX(ABS(fInterp(i,j,iVar,iEl) - fActual(i,j,iVar,iEl)),errors(iVar))
           ENDDO
         ENDDO
       ENDDO
     ENDDO

     ! Estimate routine runtime
     CALL CPU_TIME(t1)
     DO i = 1, nRepeats
       CALL interp % ScalarGridInterp_2D(fNodal_dev, fInterp_dev, nFunctions, nElem)
     ENDDO
     CALL hipCheck(hipDeviceSynchronize()) 
     CALL CPU_TIME(t2)

     runtime = (t2-t1)*REAL(nRepeats,prec)/1000.0_prec
     CALL interp % Report_JSON('ScalarGridInterp_2D', &
                               cpuModel, &
                               gpuModel, &
                               functionNames, &
                               runtime, &
                               errors, &
                               nFunctions, &
                               nElem)

     CALL hipCheck(hipFree(fNodal_dev))
     CALL hipCheck(hipFree(fInterp_dev))
#endif
     DEALLOCATE(fNodal,fInterp,fActual)

END SUBROUTINE ScalarGridInterp_2D_Test

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
!     CALL hipCheck(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))
!
!     CALL interp % ScalarGridInterp_1D(fNodal_dev, fInterp_dev, nFunctions, nElements)  
!
!     CALL hipCheck(hipMemcpy(c_loc(fInterp), fInterp_dev, SIZEOF(fInterp), hipMemcpyDeviceToHost))
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
!     CALL hipCheck(hipDeviceSynchronize()) 
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
