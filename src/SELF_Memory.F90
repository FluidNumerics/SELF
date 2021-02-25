MODULE SELF_Memory

  USE SELF_Constants

#ifdef GPU
  USE hipfort
  USE hipfort_check
#endif

  USE ISO_FORTRAN_ENV
  USE ISO_C_BINDING

  IMPLICIT NONE

  TYPE hfReal_r1
  !! Data type for storing one-dimensional real arrays on the host and the device
    REAL(prec),POINTER :: hostData(:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfReal_r1
    PROCEDURE,PUBLIC :: Free => Free_hfReal_r1

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfReal_r1
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfReal_r1
#endif

  END TYPE hfReal_r1

  TYPE hfReal_r2
  !! Data type for storing two-dimensional real arrays on the host and the device
    REAL(prec),POINTER :: hostData(:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfReal_r2
    PROCEDURE,PUBLIC :: Free => Free_hfReal_r2

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfReal_r2
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfReal_r2
#endif

  END TYPE hfReal_r2

  TYPE hfReal_r3
  !! Data type for storing three-dimensional real arrays on the host and the device
    REAL(prec),POINTER :: hostData(:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfReal_r3
    PROCEDURE,PUBLIC :: Free => Free_hfReal_r3

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfReal_r3
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfReal_r3
#endif

  END TYPE hfReal_r3

  TYPE hfReal_r4
  !! Data type for storing four-dimensional real arrays on the host and the device
    REAL(prec),POINTER :: hostData(:,:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfReal_r4
    PROCEDURE,PUBLIC :: Free => Free_hfReal_r4

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfReal_r4
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfReal_r4
#endif

  END TYPE hfReal_r4

  TYPE hfReal_r5
  !! Data type for storing five-dimensional real arrays on the host and the device
    REAL(prec),POINTER :: hostData(:,:,:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfReal_r5
    PROCEDURE,PUBLIC :: Free => Free_hfReal_r5

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfReal_r5
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfReal_r5
#endif

  END TYPE hfReal_r5

  TYPE hfReal_r6
  !! Data type for storing one-dimensional real arrays on the host and the device
    REAL(prec),POINTER :: hostData(:,:,:,:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfReal_r6
    PROCEDURE,PUBLIC :: Free => Free_hfReal_r6

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfReal_r6
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfReal_r6
#endif

  END TYPE hfReal_r6

  TYPE hfReal_r7
  !! Data type for storing seven-dimensional real arrays on the host and the device
    REAL(prec),POINTER :: hostData(:,:,:,:,:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfReal_r7
    PROCEDURE,PUBLIC :: Free => Free_hfReal_r7

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfReal_r7
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfReal_r7
#endif

  END TYPE hfReal_r7

  TYPE hfInt32_r1
  !! Data type for storing one-dimensional int32 arrays on the host and the device
    INTEGER(int32),POINTER :: hostData(:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt32_r1
    PROCEDURE,PUBLIC :: Free => Free_hfInt32_r1

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt32_r1
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt32_r1
#endif

  END TYPE hfInt32_r1

  TYPE hfInt32_r2
  !! Data type for storing two-dimensional int32 arrays on the host and the device
    INTEGER(int32),POINTER :: hostData(:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt32_r2
    PROCEDURE,PUBLIC :: Free => Free_hfInt32_r2

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt32_r2
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt32_r2
#endif

  END TYPE hfInt32_r2

  TYPE hfInt32_r3
  !! Data type for storing three-dimensional int32 arrays on the host and the device
    INTEGER(int32),POINTER :: hostData(:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt32_r3
    PROCEDURE,PUBLIC :: Free => Free_hfInt32_r3

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt32_r3
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt32_r3
#endif

  END TYPE hfInt32_r3

  TYPE hfInt32_r4
  !! Data type for storing four-dimensional int32 arrays on the host and the device
    INTEGER(int32),POINTER :: hostData(:,:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt32_r4
    PROCEDURE,PUBLIC :: Free => Free_hfInt32_r4

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt32_r4
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt32_r4
#endif

  END TYPE hfInt32_r4

  TYPE hfInt32_r5
  !! Data type for storing five-dimensional int32 arrays on the host and the device
    INTEGER(int32),POINTER :: hostData(:,:,:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt32_r5
    PROCEDURE,PUBLIC :: Free => Free_hfInt32_r5

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt32_r5
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt32_r5
#endif

  END TYPE hfInt32_r5

  TYPE hfInt32_r6
  !! Data type for storing six-dimensional int32 arrays on the host and the device
    INTEGER(int32),POINTER :: hostData(:,:,:,:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt32_r6
    PROCEDURE,PUBLIC :: Free => Free_hfInt32_r6

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt32_r6
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt32_r6
#endif

  END TYPE hfInt32_r6

  TYPE hfInt32_r7
  !! Data type for storing seven-dimensional int32 arrays on the host and the device
    INTEGER(int32),POINTER :: hostData(:,:,:,:,:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt32_r7
    PROCEDURE,PUBLIC :: Free => Free_hfInt32_r7

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt32_r7
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt32_r7
#endif

  END TYPE hfInt32_r7

  TYPE hfInt64_r1
  !! Data type for storing one-dimensional int64 arrays on the host and the device
    INTEGER(int64),POINTER :: hostData(:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt64_r1
    PROCEDURE,PUBLIC :: Free => Free_hfInt64_r1

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt64_r1
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt64_r1
#endif

  END TYPE hfInt64_r1

  TYPE hfInt64_r2
  !! Data type for storing two-dimensional int64 arrays on the host and the device
    INTEGER(int64),POINTER :: hostData(:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt64_r2
    PROCEDURE,PUBLIC :: Free => Free_hfInt64_r2

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt64_r2
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt64_r2
#endif

  END TYPE hfInt64_r2

  TYPE hfInt64_r3
  !! Data type for storing three-dimensional int64 arrays on the host and the device
    INTEGER(int64),POINTER :: hostData(:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt64_r3
    PROCEDURE,PUBLIC :: Free => Free_hfInt64_r3

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt64_r3
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt64_r3
#endif

  END TYPE hfInt64_r3

  TYPE hfInt64_r4
  !! Data type for storing four-dimensional int64 arrays on the host and the device
    INTEGER(int64),POINTER :: hostData(:,:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt64_r4
    PROCEDURE,PUBLIC :: Free => Free_hfInt64_r4

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt64_r4
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt64_r4
#endif

  END TYPE hfInt64_r4

  TYPE hfInt64_r5
  !! Data type for storing five-dimensional int64 arrays on the host and the device
    INTEGER(int64),POINTER :: hostData(:,:,:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt64_r5
    PROCEDURE,PUBLIC :: Free => Free_hfInt64_r5

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt64_r5
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt64_r5
#endif

  END TYPE hfInt64_r5

  TYPE hfInt64_r6
  !! Data type for storing six-dimensional int64 arrays on the host and the device
    INTEGER(int64),POINTER :: hostData(:,:,:,:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt64_r6
    PROCEDURE,PUBLIC :: Free => Free_hfInt64_r6

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt64_r6
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt64_r6
#endif

  END TYPE hfInt64_r6

  TYPE hfInt64_r7
  !! Data type for storing seven-dimensional int64 arrays on the host and the device
    INTEGER(int64),POINTER :: hostData(:,:,:,:,:,:,:)
    TYPE(c_ptr) :: deviceData

  CONTAINS

    PROCEDURE,PUBLIC :: Alloc => Alloc_hfInt64_r7
    PROCEDURE,PUBLIC :: Free => Free_hfInt64_r7

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_hfInt64_r7
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_hfInt64_r7
#endif

  END TYPE hfInt64_r7

CONTAINS

  SUBROUTINE Alloc_hfReal_r1(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfReal_r1),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound
    INTEGER,INTENT(in) :: upBound

    ALLOCATE (this % hostData(loBound:upBound))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfReal_r1

  SUBROUTINE Alloc_hfReal_r2(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfReal_r2),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:2)
    INTEGER,INTENT(in) :: upBound(1:2)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfReal_r2

  SUBROUTINE Alloc_hfReal_r3(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfReal_r3),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:3)
    INTEGER,INTENT(in) :: upBound(1:3)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfReal_r3

  SUBROUTINE Alloc_hfReal_r4(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfReal_r4),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:4)
    INTEGER,INTENT(in) :: upBound(1:4)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3), &
                              loBound(4):upBound(4)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfReal_r4

  SUBROUTINE Alloc_hfReal_r5(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfReal_r5),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:5)
    INTEGER,INTENT(in) :: upBound(1:5)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3), &
                              loBound(4):upBound(4), &
                              loBound(5):upBound(5)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfReal_r5

  SUBROUTINE Alloc_hfReal_r6(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfReal_r6),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:6)
    INTEGER,INTENT(in) :: upBound(1:6)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3), &
                              loBound(4):upBound(4), &
                              loBound(5):upBound(5), &
                              loBound(6):upBound(6)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfReal_r6

  SUBROUTINE Alloc_hfReal_r7(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfReal_r7),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:7)
    INTEGER,INTENT(in) :: upBound(1:7)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3), &
                              loBound(4):upBound(4), &
                              loBound(5):upBound(5), &
                              loBound(6):upBound(6), &
                              loBound(7):upBound(7)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfReal_r7

  SUBROUTINE Alloc_hfInt32_r1(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt32_r1),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound
    INTEGER,INTENT(in) :: upBound

    ALLOCATE (this % hostData(loBound:upBound))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt32_r1

  SUBROUTINE Alloc_hfInt32_r2(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt32_r2),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:2)
    INTEGER,INTENT(in) :: upBound(1:2)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt32_r2

  SUBROUTINE Alloc_hfInt32_r3(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt32_r3),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:3)
    INTEGER,INTENT(in) :: upBound(1:3)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt32_r3

  SUBROUTINE Alloc_hfInt32_r4(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt32_r4),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:4)
    INTEGER,INTENT(in) :: upBound(1:4)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3), &
                              loBound(4):upBound(4)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt32_r4

  SUBROUTINE Alloc_hfInt32_r5(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt32_r5),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:5)
    INTEGER,INTENT(in) :: upBound(1:5)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3), &
                              loBound(4):upBound(4), &
                              loBound(5):upBound(5)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt32_r5

  SUBROUTINE Alloc_hfInt32_r6(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt32_r6),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:6)
    INTEGER,INTENT(in) :: upBound(1:6)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3), &
                              loBound(4):upBound(4), &
                              loBound(5):upBound(5), &
                              loBound(6):upBound(6)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt32_r6

  SUBROUTINE Alloc_hfInt32_r7(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt32_r7),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:7)
    INTEGER,INTENT(in) :: upBound(1:7)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3), &
                              loBound(4):upBound(4), &
                              loBound(5):upBound(5), &
                              loBound(6):upBound(6), &
                              loBound(7):upBound(7)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt32_r7

  SUBROUTINE Alloc_hfInt64_r1(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt64_r1),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound
    INTEGER,INTENT(in) :: upBound

    ALLOCATE (this % hostData(loBound:upBound))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt64_r1

  SUBROUTINE Alloc_hfInt64_r2(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt64_r2),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:2)
    INTEGER,INTENT(in) :: upBound(1:2)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt64_r2

  SUBROUTINE Alloc_hfInt64_r3(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt64_r3),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:3)
    INTEGER,INTENT(in) :: upBound(1:3)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt64_r3

  SUBROUTINE Alloc_hfInt64_r4(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt64_r4),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:4)
    INTEGER,INTENT(in) :: upBound(1:4)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3), &
                              loBound(4):upBound(4)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt64_r4

  SUBROUTINE Alloc_hfInt64_r5(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt64_r5),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:5)
    INTEGER,INTENT(in) :: upBound(1:5)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3), &
                              loBound(4):upBound(4), &
                              loBound(5):upBound(5)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt64_r5

  SUBROUTINE Alloc_hfInt64_r6(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt64_r6),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:6)
    INTEGER,INTENT(in) :: upBound(1:6)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3), &
                              loBound(4):upBound(4), &
                              loBound(5):upBound(5), &
                              loBound(6):upBound(6)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt64_r6

  SUBROUTINE Alloc_hfInt64_r7(this,loBound,upBound)
    IMPLICIT NONE
    CLASS(hfInt64_r7),INTENT(out) :: this
    INTEGER,INTENT(in) :: loBound(1:7)
    INTEGER,INTENT(in) :: upBound(1:7)

    ALLOCATE (this % hostData(loBound(1):upBound(1), &
                              loBound(2):upBound(2), &
                              loBound(3):upBound(3), &
                              loBound(4):upBound(4), &
                              loBound(5):upBound(5), &
                              loBound(6):upBound(6), &
                              loBound(7):upBound(7)))

#ifdef GPU
    CALL hipCheck(hipMalloc(this % deviceData,SIZEOF(this % hostData)))
#endif

  END SUBROUTINE Alloc_hfInt64_r7

  SUBROUTINE Free_hfReal_r1(this)
    IMPLICIT NONE
    CLASS(hfReal_r1),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfReal_r1

  SUBROUTINE Free_hfReal_r2(this)
    IMPLICIT NONE
    CLASS(hfReal_r2),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfReal_r2

  SUBROUTINE Free_hfReal_r3(this)
    IMPLICIT NONE
    CLASS(hfReal_r3),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfReal_r3

  SUBROUTINE Free_hfReal_r4(this)
    IMPLICIT NONE
    CLASS(hfReal_r4),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfReal_r4

  SUBROUTINE Free_hfReal_r5(this)
    IMPLICIT NONE
    CLASS(hfReal_r5),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfReal_r5

  SUBROUTINE Free_hfReal_r6(this)
    IMPLICIT NONE
    CLASS(hfReal_r6),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfReal_r6

  SUBROUTINE Free_hfReal_r7(this)
    IMPLICIT NONE
    CLASS(hfReal_r7),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfReal_r7

  SUBROUTINE Free_hfInt32_r1(this)
    IMPLICIT NONE
    CLASS(hfInt32_r1),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt32_r1

  SUBROUTINE Free_hfInt32_r2(this)
    IMPLICIT NONE
    CLASS(hfInt32_r2),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt32_r2

  SUBROUTINE Free_hfInt32_r3(this)
    IMPLICIT NONE
    CLASS(hfInt32_r3),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt32_r3

  SUBROUTINE Free_hfInt32_r4(this)
    IMPLICIT NONE
    CLASS(hfInt32_r4),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt32_r4

  SUBROUTINE Free_hfInt32_r5(this)
    IMPLICIT NONE
    CLASS(hfInt32_r5),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt32_r5

  SUBROUTINE Free_hfInt32_r6(this)
    IMPLICIT NONE
    CLASS(hfInt32_r6),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt32_r6

  SUBROUTINE Free_hfInt32_r7(this)
    IMPLICIT NONE
    CLASS(hfInt32_r7),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt32_r7

  SUBROUTINE Free_hfInt64_r1(this)
    IMPLICIT NONE
    CLASS(hfInt64_r1),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt64_r1

  SUBROUTINE Free_hfInt64_r2(this)
    IMPLICIT NONE
    CLASS(hfInt64_r2),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt64_r2

  SUBROUTINE Free_hfInt64_r3(this)
    IMPLICIT NONE
    CLASS(hfInt64_r3),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt64_r3

  SUBROUTINE Free_hfInt64_r4(this)
    IMPLICIT NONE
    CLASS(hfInt64_r4),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt64_r4

  SUBROUTINE Free_hfInt64_r5(this)
    IMPLICIT NONE
    CLASS(hfInt64_r5),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt64_r5

  SUBROUTINE Free_hfInt64_r6(this)
    IMPLICIT NONE
    CLASS(hfInt64_r6),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt64_r6

  SUBROUTINE Free_hfInt64_r7(this)
    IMPLICIT NONE
    CLASS(hfInt64_r7),INTENT(inout) :: this

    DEALLOCATE (this % hostData)
#ifdef GPU
    CALL hipCheck(hipFree(this % deviceData))
#endif

  END SUBROUTINE Free_hfInt64_r7

#ifdef GPU
  SUBROUTINE UpdateHost_hfReal_r1(this)
    IMPLICIT NONE
    CLASS(hfReal_r1),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfReal_r1

  SUBROUTINE UpdateHost_hfReal_r2(this)
    IMPLICIT NONE
    CLASS(hfReal_r2),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfReal_r2

  SUBROUTINE UpdateHost_hfReal_r3(this)
    IMPLICIT NONE
    CLASS(hfReal_r3),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfReal_r3

  SUBROUTINE UpdateHost_hfReal_r4(this)
    IMPLICIT NONE
    CLASS(hfReal_r4),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfReal_r4

  SUBROUTINE UpdateHost_hfReal_r5(this)
    IMPLICIT NONE
    CLASS(hfReal_r5),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfReal_r5

  SUBROUTINE UpdateHost_hfReal_r6(this)
    IMPLICIT NONE
    CLASS(hfReal_r6),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfReal_r6

  SUBROUTINE UpdateHost_hfReal_r7(this)
    IMPLICIT NONE
    CLASS(hfReal_r7),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfReal_r7

  SUBROUTINE UpdateHost_hfInt32_r1(this)
    IMPLICIT NONE
    CLASS(hfInt32_r1),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt32_r1

  SUBROUTINE UpdateHost_hfInt32_r2(this)
    IMPLICIT NONE
    CLASS(hfInt32_r2),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt32_r2

  SUBROUTINE UpdateHost_hfInt32_r3(this)
    IMPLICIT NONE
    CLASS(hfInt32_r3),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt32_r3

  SUBROUTINE UpdateHost_hfInt32_r4(this)
    IMPLICIT NONE
    CLASS(hfInt32_r4),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt32_r4

  SUBROUTINE UpdateHost_hfInt32_r5(this)
    IMPLICIT NONE
    CLASS(hfInt32_r5),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt32_r5

  SUBROUTINE UpdateHost_hfInt32_r6(this)
    IMPLICIT NONE
    CLASS(hfInt32_r6),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt32_r6

  SUBROUTINE UpdateHost_hfInt32_r7(this)
    IMPLICIT NONE
    CLASS(hfInt32_r7),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt32_r7

  SUBROUTINE UpdateHost_hfInt64_r1(this)
    IMPLICIT NONE
    CLASS(hfInt64_r1),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt64_r1

  SUBROUTINE UpdateHost_hfInt64_r2(this)
    IMPLICIT NONE
    CLASS(hfInt64_r2),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt64_r2

  SUBROUTINE UpdateHost_hfInt64_r3(this)
    IMPLICIT NONE
    CLASS(hfInt64_r3),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt64_r3

  SUBROUTINE UpdateHost_hfInt64_r4(this)
    IMPLICIT NONE
    CLASS(hfInt64_r4),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt64_r4

  SUBROUTINE UpdateHost_hfInt64_r5(this)
    IMPLICIT NONE
    CLASS(hfInt64_r5),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt64_r5

  SUBROUTINE UpdateHost_hfInt64_r6(this)
    IMPLICIT NONE
    CLASS(hfInt64_r6),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt64_r6

  SUBROUTINE UpdateHost_hfInt64_r7(this)
    IMPLICIT NONE
    CLASS(hfInt64_r7),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(c_loc(this % hostData), &
                            this % deviceData, &
                            SIZEOF(this % hostData), &
                            hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_hfInt64_r7

  SUBROUTINE UpdateDevice_hfReal_r1(this)
    IMPLICIT NONE
    CLASS(hfReal_r1),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfReal_r1

  SUBROUTINE UpdateDevice_hfReal_r2(this)
    IMPLICIT NONE
    CLASS(hfReal_r2),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfReal_r2

  SUBROUTINE UpdateDevice_hfReal_r3(this)
    IMPLICIT NONE
    CLASS(hfReal_r3),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfReal_r3

  SUBROUTINE UpdateDevice_hfReal_r4(this)
    IMPLICIT NONE
    CLASS(hfReal_r4),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfReal_r4

  SUBROUTINE UpdateDevice_hfReal_r5(this)
    IMPLICIT NONE
    CLASS(hfReal_r5),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfReal_r5

  SUBROUTINE UpdateDevice_hfReal_r6(this)
    IMPLICIT NONE
    CLASS(hfReal_r6),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfReal_r6

  SUBROUTINE UpdateDevice_hfReal_r7(this)
    IMPLICIT NONE
    CLASS(hfReal_r7),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfReal_r7

  SUBROUTINE UpdateDevice_hfInt32_r1(this)
    IMPLICIT NONE
    CLASS(hfInt32_r1),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt32_r1

  SUBROUTINE UpdateDevice_hfInt32_r2(this)
    IMPLICIT NONE
    CLASS(hfInt32_r2),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt32_r2

  SUBROUTINE UpdateDevice_hfInt32_r3(this)
    IMPLICIT NONE
    CLASS(hfInt32_r3),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt32_r3

  SUBROUTINE UpdateDevice_hfInt32_r4(this)
    IMPLICIT NONE
    CLASS(hfInt32_r4),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt32_r4

  SUBROUTINE UpdateDevice_hfInt32_r5(this)
    IMPLICIT NONE
    CLASS(hfInt32_r5),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt32_r5

  SUBROUTINE UpdateDevice_hfInt32_r6(this)
    IMPLICIT NONE
    CLASS(hfInt32_r6),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt32_r6

  SUBROUTINE UpdateDevice_hfInt32_r7(this)
    IMPLICIT NONE
    CLASS(hfInt32_r7),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt32_r7

  SUBROUTINE UpdateDevice_hfInt64_r1(this)
    IMPLICIT NONE
    CLASS(hfInt64_r1),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt64_r1

  SUBROUTINE UpdateDevice_hfInt64_r2(this)
    IMPLICIT NONE
    CLASS(hfInt64_r2),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt64_r2

  SUBROUTINE UpdateDevice_hfInt64_r3(this)
    IMPLICIT NONE
    CLASS(hfInt64_r3),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt64_r3

  SUBROUTINE UpdateDevice_hfInt64_r4(this)
    IMPLICIT NONE
    CLASS(hfInt64_r4),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt64_r4

  SUBROUTINE UpdateDevice_hfInt64_r5(this)
    IMPLICIT NONE
    CLASS(hfInt64_r5),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt64_r5

  SUBROUTINE UpdateDevice_hfInt64_r6(this)
    IMPLICIT NONE
    CLASS(hfInt64_r6),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt64_r6

  SUBROUTINE UpdateDevice_hfInt64_r7(this)
    IMPLICIT NONE
    CLASS(hfInt64_r7),INTENT(inout) :: this

    CALL hipCheck(hipMemcpy(this % deviceData, &
                            c_loc(this % hostData), &
                            SIZEOF(this % hostData), &
                            hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_hfInt64_r7
#endif

END MODULE SELF_Memory
