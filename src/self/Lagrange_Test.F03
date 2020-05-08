PROGRAM Lagrange_Test

USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
USE Lagrange_Class
USE EquationParser_Class
USE Quadrature

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

  END SUBROUTINE Load_JSON
 
  SUBROUTINE Scalar_1D_Tests(json)
    USE ISO_FORTRAN_ENV
    USE ISO_C_BINDING
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
    INTEGER :: nPlotPoints, i, iEl, pdeg, nelements, nRepeats
    REAL(prec), ALLOCATABLE :: interpNodes(:), qWeights(:), targetNodes(:)
    REAL(prec), ALLOCATABLE, TARGET :: fNodal(:,:,:), fInterp(:,:,:), dfdxInterp(:,:,:),  fActual(:,:,:), dfdxActual(:,:,:)  
#ifdef GPU
    TYPE(c_ptr) :: fNodal_dev, fInterp_dev, dfdxInterp_dev
#endif
    REAL(prec) :: fErr, t1, t2, dx, x, xL, xR
    TYPE(JSON_CORE) :: json_output
    TYPE(JSON_VALUE), POINTER :: p, resArray, res, buildconf
  

      ! Initialize output json
      CALL json_output % Initialize()
      CALL json_output % create_object(p,'')
      ! Record system information
      CALL json_output % create_object(buildconf,'build')
      CALL json_output % add(p, buildconf)
      CALL json_output % add(buildconf, 'compiler', COMPILER_VERSION())
      CALL json_output % add(buildconf, 'compiler_options', COMPILER_OPTIONS())
      CALL json_output % add(buildconf, 'precision', prec)

      CALL json_output % create_array(resArray,'results')


      CALL json % info('scalar_1d', n_children=nTests)
      CALL json % get('scalar_1d', objPointer, found)
      CALL json % get('polynomial_range', polyDeg, found)
      CALL json % get('n_plot_points', nPlotPoints, found)
      CALL json % get('element_dimensions', nElements, found)
      CALL json % get('n_timing_repeats', nRepeats, found)
      CALL json % get_core(jCore)

      ALLOCATE(targetNodes(0:nPlotPoints), fActual(0:nPlotPoints,1,1:nElements), fInterp(0:nPlotPoints,1,1:nElements))
#ifdef GPU      
      CALL hfMalloc(fInterp_dev, SIZEOF(fInterp))
#endif

      targetNodes = UniformPoints(-1.0_prec,1.0_prec,nPlotPoints)

      dx = 2.0_prec/REAL(nElements,prec)

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
          DO iEl = 1, nElements
            DO i = 0, nPlotPoints
              xL = -1.0_prec + REAL((iEl-1),prec)*dx
              xR = xL + dx
              x = 0.5_prec*( xR*(targetNodes(i)+1.0_prec) - xL*(targetNodes(i)-1.0_prec) )
              fActual(i,1,iEl) = f % Evaluate( (/ x, 0.0_prec, 0.0_prec /) )          
            ENDDO
          ENDDO


          DO pdeg = polyDeg(1), polyDeg(2)

            ALLOCATE(interpNodes(0:pdeg), qWeights(0:pdeg), fNodal(0:pdeg,1,1:nElements), dfdxInterp(0:pdeg,1,1:nElements), dfdxActual(0:pdeg,1,1:nElements))

            CALL LegendreQuadrature( pdeg, interpNodes, qWeights, GAUSS )
            CALL interp % Build(pdeg, nPlotPoints, interpNodes, targetNodes)

             ! Evaluate the nodal values at quadrature points and the exact derivative
            DO iEl = 1, nElements
              DO i = 0, pdeg
                xL = -1.0_prec + REAL((iEl-1),prec)*dx
                xR = xL + dx
                x = 0.5_prec*( xR*(interpNodes(i)+1.0_prec) - xL*(interpNodes(i)-1.0_prec) )
                fNodal(i,1,iEl) = f % Evaluate( (/ x, 0.0_prec, 0.0_prec /) )          
                dfdxActual(i,1,iEl) = dfdx % Evaluate( (/ x, 0.0_prec, 0.0_prec /) )
              ENDDO
            ENDDO
            ! ***** GridInterpolate_1D Testing ****** !
            ! Interpolate the function to the targetNodes
            CALL interp % GridInterpolate_1D(fNodal, fInterp, 1, nElements)  
            fErr = 0.0_prec
            DO iEl = 1, nElements
              DO i = 0, nPlotPoints
                fErr = MAX(ABS(fInterp(i,1,iEl) -fActual(i,1,iEl)),fErr)
              ENDDO
            ENDDO

            ! Estimate routine runtime
            CALL CPU_TIME(t1)
            DO i = 1, nRepeats
              CALL interp % GridInterpolate_1D(fNodal, fInterp, 1, nElements)
            ENDDO
            CALL CPU_TIME(t2)

            CALL json_output % create_object(res,'')
            CALL json_output % add(res, 'routine_name', 'GridInterpolate_1D')
            CALL json_output % add(res, 'polynomial_degree', pdeg)
            CALL json_output % add(res, 'n_variables', 1)
            CALL json_output % add(res, 'n_elements', nElements)
            CALL json_output % add(res, 'function', TRIM(funcValue))
            CALL json_output % add(res, 'max_error', fErr)
            CALL json_output % add(res, 'runtime_ms', (t2-t1)*REAL(nRepeats,prec)/1000.0_prec)
            CALL json_output % add(resArray, res)
            nullify(res)
            ! ***** !
#ifdef GPU
            CALL hfMalloc(fNodal_dev, SIZEOF(fNodal))

            ! ***** GridInterpolate_1D Testing ****** !
            ! Interpolate the function to the targetNodes
            CALL hipFortran(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))

            CALL interp % GridInterpolate_1D(fNodal_dev, fInterp_dev, 1, nElements)  

            CALL hipFortran(hipMemcpy(c_loc(fInterp), fInterp_dev, SIZEOF(fInterp), hipMemcpyDeviceToHost))

            fErr = 0.0_prec
            DO iEl = 1, nElements
              DO i = 0, nPlotPoints
                fErr = MAX(ABS(fInterp(i,1,iEl) -fActual(i,1,iEl)),fErr)
              ENDDO
            ENDDO

            ! Estimate routine runtime
            CALL CPU_TIME(t1)
            DO i = 1, nRepeats
              CALL interp % GridInterpolate_1D(fNodal_dev, fInterp_dev, 1, nElements)
            ENDDO
            CALL hipFortran(hipDeviceSynchronize()) 
            CALL CPU_TIME(t2)

            CALL json_output % create_object(res,'')
            CALL json_output % add(res, 'routine_name', 'GridInterpolate_1D (GPU)')
            CALL json_output % add(res, 'polynomial_degree', pdeg)
            CALL json_output % add(res, 'n_variables', 1)
            CALL json_output % add(res, 'n_elements', nElements)
            CALL json_output % add(res, 'function', TRIM(funcValue))
            CALL json_output % add(res, 'max_error', fErr)
            CALL json_output % add(res, 'runtime_ms', (t2-t1)*REAL(nRepeats,prec)/1000.0_prec)
            CALL json_output % add(resArray, res)
            nullify(res)
            ! ***** !
#endif


            ! ***** GridDerivative_1D Testing ****** !
            ! Estimate the derivative by applying the derivative matrix
            !CALL interp % GridDerivative_1D(fNodal, dfdxInterp, 1, nElements)
            fErr = 0.0_prec
           ! DO iEl = 1, nElements
           !   DO i = 0, pdeg
           !     fErr = MAX(ABS(dfdxInterp(i,1,iEl)*2.0_prec/dx - dfdxActual(i,1,iEl)),fErr)
           !   ENDDO
           ! ENDDO

            ! Estimate routine runtime
            CALL CPU_TIME(t1)
           ! DO i = 1, nRepeats
           !   CALL interp % GridDerivative_1D(fNodal, dfdxInterp, 1, nElements)
           ! ENDDO
            CALL CPU_TIME(t2)

            CALL json_output % create_object(res,'')
            CALL json_output % add(res, 'routine_name', 'GridDerivative_1D')
            CALL json_output % add(res, 'polynomial_degree', pdeg)
            CALL json_output % add(res, 'n_variables', 1)
            CALL json_output % add(res, 'n_elements', nElements)
            CALL json_output % add(res, 'function', TRIM(funcValue))
            CALL json_output % add(res, 'max_error', fErr)
            CALL json_output % add(res, 'runtime_ms', (t2-t1)*REAL(nRepeats,prec)/1000.0_prec)
            CALL json_output % add(resArray, res)
            nullify(res)
            ! ***** !

#ifdef GPU
            CALL hfMalloc(dfdxInterp_dev, SIZEOF(dfdxInterp))

            ! ***** GridDerivative_1D Testing ****** !
            ! Estimate the derivative by applying the derivative matrix
            CALL hipFortran(hipMemcpy(fNodal_dev, c_loc(fNodal), SIZEOF(fNodal), hipMemcpyHostToDevice))
           ! CALL interp % GridDerivative_1D(fNodal_dev, dfdxInterp_dev, 1, nElements)
            CALL hipFortran(hipMemcpy(c_loc(dfdxInterp), dfdxInterp_dev, SIZEOF(dfdxInterp), hipMemcpyDeviceToHost))

            fErr = 0.0_prec
           ! DO iEl = 1, nElements
           !   DO i = 0, pdeg
           !     fErr = MAX(ABS(dfdxInterp(i,1,iEl)*2.0_prec/dx - dfdxActual(i,1,iEl)),fErr)
           !   ENDDO
           ! ENDDO

            ! Estimate routine runtime
            CALL CPU_TIME(t1)
           ! DO i = 1, nRepeats
           !   CALL interp % GridDerivative_1D(fNodal_dev, dfdxInterp_dev, 1, nElements)
           ! ENDDO
           ! CALL hipFortran(hipDeviceSynchronize()) 
            CALL CPU_TIME(t2)

            CALL json_output % create_object(res,'')
            CALL json_output % add(res, 'routine_name', 'GridDerivative_1D (GPU)')
            CALL json_output % add(res, 'polynomial_degree', pdeg)
            CALL json_output % add(res, 'n_variables', 1)
            CALL json_output % add(res, 'n_elements', nElements)
            CALL json_output % add(res, 'function', TRIM(funcValue))
            CALL json_output % add(res, 'max_error', fErr)
            CALL json_output % add(res, 'runtime_ms', (t2-t1)*REAL(nRepeats,prec)/1000.0_prec)
            CALL json_output % add(resArray, res)
            nullify(res)
            ! ***** !
#endif

            DEALLOCATE(interpNodes, qWeights, fNodal, dfdxInterp, dfdxActual)
#ifdef GPU      
            CALL hfFree(fNodal_dev)
            CALL hfFree(dfdxInterp_dev)
#endif

          ENDDO

        ELSE
          PRINT*, 'FAIL!'
          STOP 1
        ENDIF
             
        nullify(testPointer)
        
      ENDDO

      CALL json_output % add(p, resArray) ! Add the results array to the payload
      CALL json_output % print(p,'lagrange.results.json')

      nullify(objPointer)
      nullify(p)
      CALL json_output % Destroy()
      DEALLOCATE(targetNodes, fActual, fInterp)
#ifdef GPU      
      CALL hfFree(fInterp_dev)
#endif

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
      nullify(p)

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
