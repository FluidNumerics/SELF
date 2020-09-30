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
  PRIVATE

  PUBLIC :: Scalar1D_Tests
  PUBLIC :: Scalar2D_Tests
  PUBLIC :: Scalar3D_Tests


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

#include "test_routines/Scalar1D.F90"
#include "test_routines/Scalar2D.F90"
#include "test_routines/Scalar3D.F90"

END MODULE SELF_Tests
