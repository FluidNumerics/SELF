! Lagrange_Test.F90
! 
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self-fluids@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

PROGRAM Lagrange_Test

USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
USE Lagrange_Class
USE Lagrange_Class_Tests
USE FEQParse
USE Quadrature

USE json_module

IMPLICIT NONE

    CALL Scalar_1D_Tests()

!    CALL Scalar_2D_Tests()
!
!    CALL Scalar_3D_Tests(json)

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

  SUBROUTINE Scalar_1D_Tests()
    USE ISO_FORTRAN_ENV
    USE ISO_C_BINDING
    IMPLICIT NONE
    ! Local
    TYPE(EquationParser), ALLOCATABLE :: f(:), dfdx(:)
    TYPE(Lagrange_Tests) :: interp
    TYPE(JSON_FILE) :: json
    TYPE(JSON_VALUE), POINTER :: objPointer, testPointer
    TYPE(JSON_CORE) :: jCore
    INTEGER(JSON_IK) :: nFunctions, iTest
    LOGICAL :: found
    CHARACTER(50) :: funcName, funcValue, funcDerivative 
    !INTEGER :: polyDeg(1:2) 
    INTEGER, ALLOCATABLE :: polyDeg(:) 
    INTEGER :: nPlotPoints, i, iEl, pdeg, nelements, nRepeats
    REAL(prec), ALLOCATABLE :: targetNodes(:)
  
      CALL json % Initialize()
      CALL json % Load(filename = './lagrange.test.json')
      CALL json % info('scalar_1d', n_children=nFunctions)
      CALL json % get('scalar_1d', objPointer, found)
      CALL json % get('polynomial_range', polyDeg, found)
      CALL json % get('n_plot_points', nPlotPoints, found)
      CALL json % get('element_dimensions', nElements, found)
      CALL json % get('n_timing_repeats', nRepeats, found)
      CALL json % get_core(jCore)

      ALLOCATE(targetNodes(0:nPlotPoints), f(1:nFunctions), dfdx(1:nFunctions))
      targetNodes = UniformPoints(-1.0_prec,1.0_prec,nPlotPoints)

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
        CALL interp % Build(pdeg, GAUSS, nPlotPoints, UNIFORM)
        !CALL interp % ScalarGridInterp_1D_Test(f, nFunctions, nElements, nRepeats)
        CALL interp % ScalarDerivative_1D_Test(f, dfdx, nFunctions, nElements, nRepeats)
        CALL interp % Trash()
      ENDDO

      DEALLOCATE(targetNodes, f, dfdx)

  END SUBROUTINE Scalar_1D_Tests

  SUBROUTINE Scalar_2D_Tests()
    USE ISO_FORTRAN_ENV
    USE ISO_C_BINDING
    IMPLICIT NONE
    ! Local
    TYPE(EquationParser), ALLOCATABLE :: f(:), dfdx(:), dfdy(:)
    TYPE(Lagrange_Tests) :: interp
    TYPE(JSON_FILE) :: json
    TYPE(JSON_VALUE), POINTER :: objPointer, testPointer
    TYPE(JSON_CORE) :: jCore
    INTEGER(JSON_IK) :: nFunctions, iTest
    LOGICAL :: found
    CHARACTER(50) :: funcName, funcValue, dfdxChar, dfdyChar 
    !INTEGER :: polyDeg(1:2) 
    INTEGER, ALLOCATABLE :: polyDeg(:) 
    INTEGER :: nPlotPoints, i, iEl, pdeg, nelements, nRepeats
    REAL(prec), ALLOCATABLE :: targetNodes(:)
  
      CALL json % Initialize()
      CALL json % Load(filename = './lagrange.test.json')
      CALL json % info('scalar_2d', n_children=nFunctions)
      CALL json % get('scalar_2d', objPointer, found)
      CALL json % get('polynomial_range', polyDeg, found)
      CALL json % get('n_plot_points', nPlotPoints, found)
      CALL json % get('element_dimensions', nElements, found)
      CALL json % get('n_timing_repeats', nRepeats, found)
      CALL json % get_core(jCore)

      ALLOCATE(targetNodes(0:nPlotPoints), f(1:nFunctions), dfdx(1:nFunctions), dfdy(1:nFunctions))
      targetNodes = UniformPoints(-1.0_prec,1.0_prec,nPlotPoints)

      DO iTest = 1, nFunctions
        ! Point to the i-th scalar_1d function for testing
        CALL jCore % get_child(objPointer, iTest, testPointer, found)
        IF( found )THEN
          ! Pull the function and derivative strings from the JSON
          CALL Get_Char_Obj(jCore, testPointer, 'name', funcName)
          CALL Get_Char_Obj(jCore, testPointer, 'function', funcValue)
          CALL Get_Char_Obj(jCore, testPointer, 'dfdx', dfdxChar)
          CALL Get_Char_Obj(jCore, testPointer, 'dfdy', dfdyChar)
          ! Create the exact function equation parsers
          f(iTest) = EquationParser(funcValue,(/'x','y'/))
          dfdx(iTest) = EquationParser(dfdxChar,(/'x','y'/))
          dfdy(iTest) = EquationParser(dfdyChar,(/'x','y'/))
        ENDIF
      ENDDO

      DO pdeg = polyDeg(1), polyDeg(2)
        CALL interp % Build(pdeg, GAUSS, nPlotPoints, UNIFORM)
        CALL interp % ScalarGridInterp_2D_Test(f, nFunctions, nElements, nRepeats)
        CALL interp % ScalarGradient_2D_Test(f, dfdx, dfdy, nFunctions, nElements, nRepeats)
        CALL interp % Trash()
      ENDDO

      DEALLOCATE(targetNodes, f, dfdx, dfdy)

  END SUBROUTINE Scalar_2D_Tests

!
!  SUBROUTINE Scalar_2D_Tests(json)
!    USE ISO_FORTRAN_ENV
!    USE ISO_C_BINDING
!    IMPLICIT NONE
!    TYPE(JSON_FILE), INTENT(inout) :: json
!    ! Local
!    TYPE(EquationParser) :: f, dfdx, dfdy
!    TYPE(Lagrange) :: interp
!    TYPE(JSON_VALUE), POINTER :: objPointer, testPointer
!    TYPE(JSON_CORE) :: jCore
!    INTEGER(JSON_IK) :: nTests, iTest
!    LOGICAL :: found
!    CHARACTER(50) :: funcName, funcValue, dfdx_str, dfdy_str 
!    INTEGER, ALLOCATABLE :: polyDeg(:) 
!    INTEGER :: nPlotPoints, i, j, iX, iY, iEl, pdeg, nelements, nRepeats
!    REAL(prec), ALLOCATABLE :: targetNodes(:)
!    REAL(prec), ALLOCATABLE, TARGET :: fNodal(:,:,:,:), fInterp(:,:,:,:), gradFInterp(:,:,:,:,:),  fActual(:,:,:,:), gradFActual(:,:,:,:,:)  
!#ifdef GPU
!    TYPE(c_ptr) :: fNodal_dev, fInterp_dev, gradFInterp_dev
!#endif
!    REAL(prec) :: fErr, t1, t2, dx, x, y, xL, xR, yL, yR
!    TYPE(JSON_CORE) :: json_output
!    TYPE(JSON_VALUE), POINTER :: p, resArray, res, buildconf
!  
!
!      ! Initialize output json
!      CALL json_output % Initialize()
!      CALL json_output % create_object(p,'')
!      ! Record system information
!      CALL json_output % create_object(buildconf,'build')
!      CALL json_output % add(p, buildconf)
!      CALL json_output % add(buildconf, 'compiler', COMPILER_VERSION())
!      CALL json_output % add(buildconf, 'compiler_options', COMPILER_OPTIONS())
!      CALL json_output % add(buildconf, 'precision', prec)
!
!      CALL json_output % create_array(resArray,'results')
!
!
!      CALL json % info('scalar_2d', n_children=nTests)
!      CALL json % get('scalar_2d', objPointer, found)
!      CALL json % get('polynomial_range', polyDeg, found)
!      CALL json % get('n_plot_points', nPlotPoints, found)
!      CALL json % get('element_dimensions', nElements, found)
!      CALL json % get('n_timing_repeats', nRepeats, found)
!      CALL json % get_core(jCore)
!
!      ALLOCATE(targetNodes(0:nPlotPoints), &
!               fActual(0:nPlotPoints,0:nPlotPoints,1,1:nElements*nElements), &
!               fInterp(0:nPlotPoints,0:nPlotPoints,1,1:nElements*nElements))
!#ifdef GPU      
!      CALL hfMalloc(fInterp_dev, SIZEOF(fInterp))
!#endif
!
!      targetNodes = UniformPoints(-1.0_prec,1.0_prec,nPlotPoints)
!
!      dx = 2.0_prec/REAL(nElements,prec)
!
!      DO iTest = 1, nTests
! 
!        ! Point to the i-th scalar_1d function for testing
!        CALL jCore % get_child(objPointer, iTest, testPointer, found)
!        IF( found )THEN
!
!         ! Pull the function and derivative strings from the JSON
!         CALL Get_Char_Obj(jCore, testPointer, 'name', funcName)
!         CALL Get_Char_Obj(jCore, testPointer, 'function', funcValue)
!         CALL Get_Char_Obj(jCore, testPointer, 'dfdx', dfdx_str)
!         CALL Get_Char_Obj(jCore, testPointer, 'dfdy', dfdy_str)
!
!          ! Create the exact function equation parsers
!          f = EquationParser(funcValue)
!          dfdx = EquationParser(dfdx_str)
!          dfdy = EquationParser(dfdy_str)
!
!          ! Create the exact function equation parsers
!          DO iY = 1, nElements
!            yL = -1.0_prec + REAL((iY-1),prec)*dx
!            yR = yL + dx
!            DO iX = 1, nElements
!              xL = -1.0_prec + REAL((iX-1),prec)*dx
!              xR = xL + dx
!              DO j = 0, nPlotPoints
!                y = 0.5_prec*( yR*(targetNodes(j)+1.0_prec) - yL*(targetNodes(j)-1.0_prec) )
!                DO i = 0, nPlotPoints
!                  x = 0.5_prec*( xR*(targetNodes(i)+1.0_prec) - xL*(targetNodes(i)-1.0_prec) )
!                  iEl = iX + nElements*(iY-1)
!                  fActual(i,j,1,iEl) = f % Evaluate( (/ x, y, 0.0_prec /) )          
!                ENDDO
!              ENDDO
!            ENDDO
!          ENDDO
!
!
!          DO pdeg = polyDeg(1), polyDeg(2)
!
!            ALLOCATE(fNodal(0:pdeg,0:pdeg,1,1:nElements*nElements), &
!                     gradFInterp(1:2,0:pdeg,0:pdeg,1,1:nElements*nElements), &
!                     gradFActual(1:2,0:pdeg,0:pdeg,1,1:nElements*nElements))
!
!            CALL interp % Build(pdeg, GAUSS, nPlotPoints, UNIFORM)
!
!            ! Evaluate the nodal values at quadrature points and the exact derivative
!            DO iY = 1, nElements
!              yL = -1.0_prec + REAL((iY-1),prec)*dx
!              yR = yL + dx
!              DO iX = 1, nElements
!                xL = -1.0_prec + REAL((iX-1),prec)*dx
!                xR = xL + dx
!                DO j = 0, pdeg
!                  y = 0.5_prec*( yR*(interp % controlPoints(j)+1.0_prec) - yL*(interp % controlPoints(j)-1.0_prec) )
!                  DO i = 0, pdeg
!                    x = 0.5_prec*( xR*(interp % controlPoints(i)+1.0_prec) - xL*(interp % controlPoints(i)-1.0_prec) )
!                    iEl = iX + nElements*(iY-1)
!                    fNodal(i,j,1,iEl) = f % Evaluate( (/ x, y, 0.0_prec /) )          
!                    gradFActual(1,i,j,1,iEl) = dfdx % Evaluate( (/ x, y, 0.0_prec /) )
!                    gradFActual(2,i,j,1,iEl) = dfdy % Evaluate( (/ x, y, 0.0_prec /) )
!                  ENDDO
!                ENDDO
!              ENDDO
!            ENDDO
!
!            ! ***** ScalarGridInterp_2D Testing ****** !
!            ! Interpolate the function to the targetNodes
!            CALL interp % ScalarGridInterp_2D(fNodal, fInterp, 1, nElements*nElements)  
!            fErr = 0.0_prec
!            DO iEl = 1, nElements*nElements
!              DO j = 0, nPlotPoints
!                DO i = 0, nPlotPoints
!                  fErr = MAX(ABS(fInterp(i,j,1,iEl) -fActual(i,j,1,iEl)),fErr)
!                ENDDO
!              ENDDO
!            ENDDO
!
!            ! Estimate routine runtime
!            CALL CPU_TIME(t1)
!            DO i = 1, nRepeats
!              CALL interp % ScalarGridInterp_2D(fNodal, fInterp, 1, nElements*nElements)
!            ENDDO
!            CALL CPU_TIME(t2)
!
!            CALL json_output % create_object(res,'')
!            CALL json_output % add(res, 'routine_name', 'ScalarGridInterp_2D')
!            CALL json_output % add(res, 'polynomial_degree', pdeg)
!            CALL json_output % add(res, 'n_variables', 1)
!            CALL json_output % add(res, 'n_elements', nElements*nElements)
!            CALL json_output % add(res, 'function', TRIM(funcValue))
!            CALL json_output % add(res, 'max_error', fErr)
!            CALL json_output % add(res, 'runtime_ms', (t2-t1)*REAL(nRepeats,prec)/1000.0_prec)
!            CALL json_output % add(resArray, res)
!            nullify(res)
!            ! ***** !
!#ifdef GPU
!            CALL hfMalloc(fNodal_dev, SIZEOF(fNodal))
!
!            ! ***** ScalarGridInterp_2D Testing ****** !
!            ! Interpolate the function to the targetNodes
!            CALL hipFortran(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))
!
!            CALL interp % ScalarGridInterp_2D(fNodal_dev, fInterp_dev, 1, nElements*nElements)  
!
!            CALL hipFortran(hipMemcpy(c_loc(fInterp), fInterp_dev, SIZEOF(fInterp), hipMemcpyDeviceToHost))
!
!            fErr = 0.0_prec
!            DO iEl = 1, nElements*nElements
!              DO j = 0, nPlotPoints
!                DO i = 0, nPlotPoints
!                  fErr = MAX(ABS(fInterp(i,j,1,iEl) -fActual(i,j,1,iEl)),fErr)
!                ENDDO
!              ENDDO
!            ENDDO
!
!            ! Estimate routine runtime
!            CALL CPU_TIME(t1)
!            DO i = 1, nRepeats
!              CALL interp % ScalarGridInterp_2D(fNodal_dev, fInterp_dev, 1, nElements*nElements)
!            ENDDO
!
!            CALL hipFortran(hipDeviceSynchronize()) 
!            CALL CPU_TIME(t2)
!
!            CALL json_output % create_object(res,'')
!            CALL json_output % add(res, 'routine_name', 'ScalarGridInterp_2D (GPU)')
!            CALL json_output % add(res, 'polynomial_degree', pdeg)
!            CALL json_output % add(res, 'n_variables', 1)
!            CALL json_output % add(res, 'n_elements', nElements*nElements)
!            CALL json_output % add(res, 'function', TRIM(funcValue))
!            CALL json_output % add(res, 'max_error', fErr)
!            CALL json_output % add(res, 'runtime_ms', (t2-t1)*REAL(nRepeats,prec)/1000.0_prec)
!            CALL json_output % add(resArray, res)
!            nullify(res)
!            ! ***** !
!#endif
!
!
!            ! ***** Derivative_2D Testing ****** !
!            ! Estimate the derivative by applying the derivative matrix
!            CALL interp % ScalarGradient_2D(fNodal, gradFInterp, 1, nElements*nElements)
!            fErr = 0.0_prec
!            DO iEl = 1, nElements*nElements
!              DO j = 0, pdeg
!                DO i = 0, pdeg
!                  fErr = MAX(ABS(gradFInterp(1,i,j,1,iEl)*2.0_prec/dx - gradFActual(1,i,j,1,iEl)),fErr)
!                  fErr = MAX(ABS(gradFInterp(2,i,j,1,iEl)*2.0_prec/dx - gradFActual(2,i,j,1,iEl)),fErr)
!                ENDDO
!              ENDDO
!            ENDDO
!
!            ! Estimate routine runtime
!            CALL CPU_TIME(t1)
!            DO i = 1, nRepeats
!              CALL interp % ScalarGradient_2D(fNodal, gradFInterp, 1, nElements*nElements)
!            ENDDO
!            CALL CPU_TIME(t2)
!
!            CALL json_output % create_object(res,'')
!            CALL json_output % add(res, 'routine_name', 'ScalarGradient_2D')
!            CALL json_output % add(res, 'polynomial_degree', pdeg)
!            CALL json_output % add(res, 'n_variables', 1)
!            CALL json_output % add(res, 'n_elements', nElements*nElements)
!            CALL json_output % add(res, 'function', TRIM(funcValue))
!            CALL json_output % add(res, 'max_error', fErr)
!            CALL json_output % add(res, 'runtime_ms', (t2-t1)*REAL(nRepeats,prec)/1000.0_prec)
!            CALL json_output % add(resArray, res)
!            nullify(res)
!            ! ***** !
!
!#ifdef GPU
!            CALL hfMalloc(gradFInterp_dev, SIZEOF(gradFInterp))
!
!            ! ***** Derivative_2D Testing ****** !
!            ! Estimate the derivative by applying the derivative matrix
!            CALL hipFortran(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))
!            CALL interp % ScalarGradient_2D(fNodal_dev, gradFInterp_dev, 1, nElements*nElements)
!            CALL hipFortran(hipMemcpy(c_loc(gradFInterp), gradFInterp_dev, SIZEOF(gradFInterp), hipMemcpyDeviceToHost))
!
!            fErr = 0.0_prec
!            DO iEl = 1, nElements*nElements
!              DO j = 0, pdeg
!                DO i = 0, pdeg
!                  fErr = MAX(ABS(gradFInterp(1,i,j,1,iEl)*2.0_prec/dx - gradFActual(1,i,j,1,iEl)),fErr)
!                  fErr = MAX(ABS(gradFInterp(2,i,j,1,iEl)*2.0_prec/dx - gradFActual(2,i,j,1,iEl)),fErr)
!                ENDDO
!              ENDDO
!            ENDDO
!
!            ! Estimate routine runtime
!            CALL CPU_TIME(t1)
!            DO i = 1, nRepeats
!              CALL interp % ScalarGradient_2D(fNodal_dev, gradFInterp_dev, 1, nElements*nElements)
!            ENDDO
!            CALL hipFortran(hipDeviceSynchronize()) 
!            CALL CPU_TIME(t2)
!
!            CALL json_output % create_object(res,'')
!            CALL json_output % add(res, 'routine_name', 'ScalarGradient_2D (GPU)')
!            CALL json_output % add(res, 'polynomial_degree', pdeg)
!            CALL json_output % add(res, 'n_variables', 1)
!            CALL json_output % add(res, 'n_elements', nElements*nElements)
!            CALL json_output % add(res, 'function', TRIM(funcValue))
!            CALL json_output % add(res, 'max_error', fErr)
!            CALL json_output % add(res, 'runtime_ms', (t2-t1)*REAL(nRepeats,prec)/1000.0_prec)
!            CALL json_output % add(resArray, res)
!            nullify(res)
!            ! ***** !
!#endif
!
!            DEALLOCATE(fNodal, gradFInterp, gradFActual)
!#ifdef GPU      
!            CALL hfFree(fNodal_dev)
!            CALL hfFree(gradFInterp_dev)
!#endif
!
!          ENDDO
!
!        ELSE
!          PRINT*, 'FAIL!'
!          STOP 1
!        ENDIF
!             
!        nullify(testPointer)
!        
!      ENDDO
!
!      CALL json_output % add(p, resArray) ! Add the results array to the payload
!      CALL json_output % print(p,'lagrange_scalar2d.results.json')
!
!      nullify(objPointer)
!      nullify(p)
!      CALL json_output % Destroy()
!      DEALLOCATE(targetNodes, fActual, fInterp)
!#ifdef GPU      
!      CALL hfFree(fInterp_dev)
!#endif
!
!  END SUBROUTINE Scalar_2D_Tests
!
!  SUBROUTINE Scalar_3D_Tests(json)
!    USE ISO_FORTRAN_ENV
!    USE ISO_C_BINDING
!    IMPLICIT NONE
!    TYPE(JSON_FILE), INTENT(inout) :: json
!    ! Local
!    TYPE(EquationParser) :: f, dfdx
!    TYPE(Lagrange) :: interp
!    TYPE(JSON_VALUE), POINTER :: objPointer, testPointer
!    TYPE(JSON_CORE) :: jCore
!    INTEGER(JSON_IK) :: nTests, iTest
!    LOGICAL :: found
!    CHARACTER(50) :: funcName, funcValue, funcDerivative 
!    INTEGER, ALLOCATABLE :: polyDeg(:) 
!    INTEGER :: nPlotPoints, i, j, k, iX, iY, iZ, iEl, pdeg, nelements, nRepeats
!    REAL(prec), ALLOCATABLE :: targetNodes(:)
!    REAL(prec), ALLOCATABLE, TARGET :: fNodal(:,:,:,:,:), fInterp(:,:,:,:,:), dfdxInterp(:,:,:,:,:),  fActual(:,:,:,:,:), dfdxActual(:,:,:,:,:)  
!#ifdef GPU
!    TYPE(c_ptr) :: fNodal_dev, fInterp_dev, dfdxInterp_dev
!#endif
!    REAL(prec) :: fErr, t1, t2, dx, x, y, z, xL, xR, yL, yR, zL, zR
!    TYPE(JSON_CORE) :: json_output
!    TYPE(JSON_VALUE), POINTER :: p, resArray, res, buildconf
!  
!
!      ! Initialize output json
!      CALL json_output % Initialize()
!      CALL json_output % create_object(p,'')
!      ! Record system information
!      CALL json_output % create_object(buildconf,'build')
!      CALL json_output % add(p, buildconf)
!      CALL json_output % add(buildconf, 'compiler', COMPILER_VERSION())
!      CALL json_output % add(buildconf, 'compiler_options', COMPILER_OPTIONS())
!      CALL json_output % add(buildconf, 'precision', prec)
!
!      CALL json_output % create_array(resArray,'results')
!
!
!      CALL json % info('scalar_3d', n_children=nTests)
!      CALL json % get('scalar_3d', objPointer, found)
!      CALL json % get('polynomial_range', polyDeg, found)
!      CALL json % get('n_plot_points', nPlotPoints, found)
!      CALL json % get('element_dimensions', nElements, found)
!      CALL json % get('n_timing_repeats', nRepeats, found)
!      CALL json % get_core(jCore)
!
!      ALLOCATE(targetNodes(0:nPlotPoints), &
!               fActual(0:nPlotPoints,0:nPlotPoints,0:nPlotPoints,1,1:nElements*nElements*nElements), &
!               fInterp(0:nPlotPoints,0:nPlotPoints,0:nPlotPoints,1,1:nElements*nElements*nElements))
!#ifdef GPU      
!      CALL hfMalloc(fInterp_dev, SIZEOF(fInterp))
!#endif
!
!      targetNodes = UniformPoints(-1.0_prec,1.0_prec,nPlotPoints)
!
!      dx = 2.0_prec/REAL(nElements,prec)
!
!      DO iTest = 1, nTests
! 
!        ! Point to the i-th scalar_1d function for testing
!        CALL jCore % get_child(objPointer, iTest, testPointer, found)
!        IF( found )THEN
!
!         ! Pull the function and derivative strings from the JSON
!         CALL Get_Char_Obj(jCore, testPointer, 'name', funcName)
!         CALL Get_Char_Obj(jCore, testPointer, 'function', funcValue)
!         CALL Get_Char_Obj(jCore, testPointer, 'derivative', funcDerivative)
!
!          ! Create the exact function equation parsers
!          f = EquationParser(funcValue)
!!          dfdx = EquationParser(funcDerivative)
!
!          ! Create the exact function equation parsers
!          DO iz = 1, nElements
!            zL = -1.0_prec + REAL((iz-1),prec)*dx
!            zR = zL + dx
!            DO iY = 1, nElements
!              yL = -1.0_prec + REAL((iY-1),prec)*dx
!              yR = yL + dx
!              DO iX = 1, nElements
!                xL = -1.0_prec + REAL((iX-1),prec)*dx
!                xR = xL + dx
!                DO k = 0, nPlotPoints
!                  z = 0.5_prec*( zR*(targetNodes(k)+1.0_prec) - zL*(targetNodes(k)-1.0_prec) )
!                  DO j = 0, nPlotPoints
!                    y = 0.5_prec*( yR*(targetNodes(j)+1.0_prec) - yL*(targetNodes(j)-1.0_prec) )
!                    DO i = 0, nPlotPoints
!                      x = 0.5_prec*( xR*(targetNodes(i)+1.0_prec) - xL*(targetNodes(i)-1.0_prec) )
!                      iEl = iX + nElements*(iY-1) + nElements*nElements*(iZ-1)
!                      fActual(i,j,k,1,iEl) = f % Evaluate( (/ x, y, z /) )          
!                    ENDDO
!                  ENDDO
!                ENDDO
!              ENDDO
!            ENDDO
!          ENDDO
!
!
!          DO pdeg = polyDeg(1), polyDeg(2)
!
!            ALLOCATE(fNodal(0:pdeg,0:pdeg,0:pdeg,1,1:nElements*nElements*nElements), &
!                     dfdxInterp(0:pdeg,0:pdeg,0:pdeg,1,1:nElements*nElements*nElements),&
!                     dfdxActual(0:pdeg,0:pdeg,0:pdeg,1,1:nElements*nElements*nElements))
!
!            CALL interp % Build(pdeg, GAUSS, nPlotPoints, UNIFORM)
!
!            ! Evaluate the nodal values at quadrature points and the exact derivative
!            DO iz = 1, nElements
!              zL = -1.0_prec + REAL((iz-1),prec)*dx
!              zR = zL + dx
!              DO iY = 1, nElements
!                yL = -1.0_prec + REAL((iY-1),prec)*dx
!                yR = yL + dx
!                DO iX = 1, nElements
!                  xL = -1.0_prec + REAL((iX-1),prec)*dx
!                  xR = xL + dx
!                  DO k = 0, pdeg
!                    z = 0.5_prec*( zR*(interp % controlPoints(k)+1.0_prec) - zL*(interp % controlPoints(k)-1.0_prec) )
!                    DO j = 0, pdeg
!                      y = 0.5_prec*( yR*(interp % controlPoints(j)+1.0_prec) - yL*(interp % controlPoints(j)-1.0_prec) )
!                      DO i = 0, pdeg
!                        x = 0.5_prec*( xR*(interp % controlPoints(i)+1.0_prec) - xL*(interp % controlPoints(i)-1.0_prec) )
!                        iEl = iX + nElements*(iY-1) + nElements*nElements*(iZ-1)
!                        fNodal(i,j,k,1,iEl) = f % Evaluate( (/ x, y, z /) )          
!                      ENDDO
!                    ENDDO
!                  ENDDO
!                ENDDO
!              ENDDO
!            ENDDO
!
!            ! ***** ScalarGridInterp_3D Testing ****** !
!            ! Interpolate the function to the targetNodes
!            CALL interp % ScalarGridInterp_3D(fNodal, fInterp, 1, nElements*nElements*nElements)  
!            fErr = 0.0_prec
!            DO iEl = 1, nElements*nElements*nElements
!              DO k = 0, nPlotPoints
!                DO j = 0, nPlotPoints
!                  DO i = 0, nPlotPoints
!                    fErr = MAX(ABS(fInterp(i,j,k,1,iEl) -fActual(i,j,k,1,iEl)),fErr)
!                  ENDDO
!                ENDDO
!              ENDDO
!            ENDDO
!
!            ! Estimate routine runtime
!            CALL CPU_TIME(t1)
!            DO i = 1, nRepeats
!              CALL interp % ScalarGridInterp_3D(fNodal, fInterp, 1, nElements*nElements*nElements)
!            ENDDO
!            CALL CPU_TIME(t2)
!
!            CALL json_output % create_object(res,'')
!            CALL json_output % add(res, 'routine_name', 'ScalarGridInterp_3D')
!            CALL json_output % add(res, 'polynomial_degree', pdeg)
!            CALL json_output % add(res, 'n_variables', 1)
!            CALL json_output % add(res, 'n_elements', nElements*nElements*nElements)
!            CALL json_output % add(res, 'function', TRIM(funcValue))
!            CALL json_output % add(res, 'max_error', fErr)
!            CALL json_output % add(res, 'runtime_ms', (t2-t1)*REAL(nRepeats,prec)/1000.0_prec)
!            CALL json_output % add(resArray, res)
!            nullify(res)
!            ! ***** !
!#ifdef GPU
!            CALL hfMalloc(fNodal_dev, SIZEOF(fNodal))
!
!            ! ***** ScalarGridInterp_3D Testing ****** !
!            ! Interpolate the function to the targetNodes
!            CALL hipFortran(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))
!
!            CALL interp % ScalarGridInterp_3D(fNodal_dev, fInterp_dev, 1, nElements*nElements*nElements)  
!
!            CALL hipFortran(hipMemcpy(c_loc(fInterp), fInterp_dev, SIZEOF(fInterp), hipMemcpyDeviceToHost))
!
!            fErr = 0.0_prec
!            DO iEl = 1, nElements*nElements
!              DO k = 0, nPlotPoints
!                DO j = 0, nPlotPoints
!                  DO i = 0, nPlotPoints
!                    fErr = MAX(ABS(fInterp(i,j,k,1,iEl) -fActual(i,j,k,1,iEl)),fErr)
!                  ENDDO
!                ENDDO
!              ENDDO
!            ENDDO
!
!            ! Estimate routine runtime
!            CALL CPU_TIME(t1)
!            DO i = 1, nRepeats
!              CALL interp % ScalarGridInterp_3D(fNodal_dev, fInterp_dev, 1, nElements*nElements*nElements)
!            ENDDO
!
!            CALL hipFortran(hipDeviceSynchronize()) 
!            CALL CPU_TIME(t2)
!
!            CALL json_output % create_object(res,'')
!            CALL json_output % add(res, 'routine_name', 'ScalarGridInterp_3D (GPU)')
!            CALL json_output % add(res, 'polynomial_degree', pdeg)
!            CALL json_output % add(res, 'n_variables', 1)
!            CALL json_output % add(res, 'n_elements', nElements*nElements*nElements)
!            CALL json_output % add(res, 'function', TRIM(funcValue))
!            CALL json_output % add(res, 'max_error', fErr)
!            CALL json_output % add(res, 'runtime_ms', (t2-t1)*REAL(nRepeats,prec)/1000.0_prec)
!            CALL json_output % add(resArray, res)
!            nullify(res)
!            ! ***** !
!#endif
!
!!
!!            ! ***** Derivative_3D Testing ****** !
!!            ! Estimate the derivative by applying the derivative matrix
!!            CALL interp % Derivative_3D(fNodal, dfdxInterp, 1, nElements)
!!            fErr = 0.0_prec
!!            DO iEl = 1, nElements
!!              DO i = 0, pdeg
!!                fErr = MAX(ABS(dfdxInterp(i,1,iEl)*2.0_prec/dx - dfdxActual(i,1,iEl)),fErr)
!!              ENDDO
!!            ENDDO
!!
!!            ! Estimate routine runtime
!!            CALL CPU_TIME(t1)
!!            DO i = 1, nRepeats
!!              CALL interp % Derivative_3D(fNodal, dfdxInterp, 1, nElements)
!!            ENDDO
!!            CALL CPU_TIME(t2)
!!
!!            CALL json_output % create_object(res,'')
!!            CALL json_output % add(res, 'routine_name', 'Derivative_3D')
!!            CALL json_output % add(res, 'polynomial_degree', pdeg)
!!            CALL json_output % add(res, 'n_variables', 1)
!!            CALL json_output % add(res, 'n_elements', nElements)
!!            CALL json_output % add(res, 'function', TRIM(funcValue))
!!            CALL json_output % add(res, 'max_error', fErr)
!!            CALL json_output % add(res, 'runtime_ms', (t2-t1)*REAL(nRepeats,prec)/1000.0_prec)
!!            CALL json_output % add(resArray, res)
!!            nullify(res)
!!            ! ***** !
!!
!!#ifdef GPU
!!            CALL hfMalloc(dfdxInterp_dev, SIZEOF(dfdxInterp))
!!
!!            ! ***** Derivative_3D Testing ****** !
!!            ! Estimate the derivative by applying the derivative matrix
!!            CALL hipFortran(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))
!!            CALL interp % Derivative_3D(fNodal_dev, dfdxInterp_dev, 1, nElements)
!!            CALL hipFortran(hipMemcpy(c_loc(dfdxInterp), dfdxInterp_dev, SIZEOF(dfdxInterp), hipMemcpyDeviceToHost))
!!
!!            fErr = 0.0_prec
!!            DO iEl = 1, nElements
!!              DO i = 0, pdeg
!!                fErr = MAX(ABS(dfdxInterp(i,1,iEl)*2.0_prec/dx - dfdxActual(i,1,iEl)),fErr)
!!              ENDDO
!!            ENDDO
!!
!!            ! Estimate routine runtime
!!            CALL CPU_TIME(t1)
!!            DO i = 1, nRepeats
!!              CALL interp % Derivative_3D(fNodal_dev, dfdxInterp_dev, 1, nElements)
!!            ENDDO
!!            CALL hipFortran(hipDeviceSynchronize()) 
!!            CALL CPU_TIME(t2)
!!
!!            CALL json_output % create_object(res,'')
!!            CALL json_output % add(res, 'routine_name', 'Derivative_3D (GPU)')
!!            CALL json_output % add(res, 'polynomial_degree', pdeg)
!!            CALL json_output % add(res, 'n_variables', 1)
!!            CALL json_output % add(res, 'n_elements', nElements)
!!            CALL json_output % add(res, 'function', TRIM(funcValue))
!!            CALL json_output % add(res, 'max_error', fErr)
!!            CALL json_output % add(res, 'runtime_ms', (t2-t1)*REAL(nRepeats,prec)/1000.0_prec)
!!            CALL json_output % add(resArray, res)
!!            nullify(res)
!!            ! ***** !
!!#endif
!
!            DEALLOCATE(fNodal, dfdxInterp, dfdxActual)
!#ifdef GPU      
!            CALL hfFree(fNodal_dev)
!!            CALL hfFree(dfdxInterp_dev)
!#endif
!
!          ENDDO
!
!        ELSE
!          PRINT*, 'FAIL!'
!          STOP 1
!        ENDIF
!             
!        nullify(testPointer)
!        
!      ENDDO
!
!      CALL json_output % add(p, resArray) ! Add the results array to the payload
!      CALL json_output % print(p,'lagrange_scalar3d.results.json')
!
!      nullify(objPointer)
!      nullify(p)
!      CALL json_output % Destroy()
!      DEALLOCATE(targetNodes, fActual, fInterp)
!#ifdef GPU      
!      CALL hfFree(fInterp_dev)
!#endif
!
!  END SUBROUTINE Scalar_3D_Tests
!
!
!
!! SUBROUTINE Scalar_2D_Tests()
!! END SUBROUTINE Scalar_2D_Tests
!!
!! SUBROUTINE Scalar_3D_Tests()
!! END SUBROUTINE Scalar_3D_Tests
!!
!! SUBROUTINE Vector_2D_Tests()
!! END SUBROUTINE Vector_2D_Tests
!!
!! SUBROUTINE Vector_3D_Tests()
!! END SUBROUTINE Vector_3D_Tests

END PROGRAM Lagrange_Test
