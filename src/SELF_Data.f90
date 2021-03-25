! SELF_Data.F90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self-fluids@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE SELF_Data

  USE SELF_Constants
  USE SELF_Lagrange

  USE ISO_C_BINDING

  IMPLICIT NONE

! ---------------------- Scalars ---------------------- !
  TYPE,PUBLIC :: Scalar1D

    INTEGER :: N
    INTEGER :: M
    INTEGER :: nVar
    INTEGER :: nElem
    INTEGER :: controlType
    INTEGER :: targetType
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r3) :: interior
    TYPE(hfReal_r3) :: boundary

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Scalar1D
    PROCEDURE,PUBLIC :: Free => Free_Scalar1D

    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Scalar1D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Scalar1D

    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Scalar1D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Scalar1D

    GENERIC,PUBLIC :: Derivative => Derivative_Scalar1D
    PROCEDURE,PRIVATE :: Derivative_Scalar1D


    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Scalar1D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Scalar1D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Scalar1D
    GENERIC,PUBLIC :: OPERATOR(+) => Add_Scalar1D
    GENERIC,PUBLIC :: OPERATOR(-) => Subtract_Scalar1D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Scalar1D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Add_Scalar1D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Subtract_Scalar1D

  END TYPE Scalar1D

  TYPE,PUBLIC :: Scalar2D

    INTEGER :: N
    INTEGER :: M
    INTEGER :: nVar
    INTEGER :: nElem
    INTEGER :: controlType
    INTEGER :: targetType
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r4) :: interior
    TYPE(hfReal_r4) :: boundary

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Scalar2D
    PROCEDURE,PUBLIC :: Free => Free_Scalar2D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Scalar2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Scalar2D
    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Scalar2D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Scalar2D

    GENERIC,PUBLIC :: Gradient => Gradient_Scalar2D
    PROCEDURE,PRIVATE :: Gradient_Scalar2D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Scalar2D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Scalar2D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Scalar2D
    GENERIC,PUBLIC :: OPERATOR(+) => Add_Scalar2D
    GENERIC,PUBLIC :: OPERATOR(-) => Subtract_Scalar2D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Scalar2D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Add_Scalar2D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Subtract_Scalar2D

  END TYPE Scalar2D

  TYPE,PUBLIC :: Scalar3D

    INTEGER :: N
    INTEGER :: M
    INTEGER :: nVar
    INTEGER :: nElem
    INTEGER :: controlType
    INTEGER :: targetType
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r5) :: interior
    TYPE(hfReal_r5) :: boundary

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Scalar3D
    PROCEDURE,PUBLIC :: Free => Free_Scalar3D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Scalar3D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Scalar3D
    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Scalar3D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Scalar3D

    GENERIC,PUBLIC :: Gradient => Gradient_Scalar3D
    PROCEDURE,PRIVATE :: Gradient_Scalar3D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Scalar3D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Scalar3D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Scalar3D
    GENERIC,PUBLIC :: OPERATOR(+) => Add_Scalar3D
    GENERIC,PUBLIC :: OPERATOR(-) => Subtract_Scalar3D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Scalar3D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Add_Scalar3D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Subtract_Scalar3D

  END TYPE Scalar3D

! ---------------------- Vectors ---------------------- !

  TYPE,PUBLIC :: Vector2D

    INTEGER :: N
    INTEGER :: M
    INTEGER :: nVar
    INTEGER :: nElem
    INTEGER :: controlType
    INTEGER :: targetType
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r5) :: interior
    TYPE(hfReal_r5) :: boundary

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Vector2D
    PROCEDURE,PUBLIC :: Free => Free_Vector2D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Vector2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Vector2D
    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Vector2D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Vector2D

    GENERIC,PUBLIC :: Gradient => Gradient_Vector2D
    PROCEDURE,PRIVATE :: Gradient_Vector2D

    GENERIC,PUBLIC :: Divergence => Divergence_Vector2D
    PROCEDURE,PRIVATE :: Divergence_Vector2D

    GENERIC,PUBLIC :: Curl => Curl_Vector2D
    PROCEDURE,PRIVATE :: Curl_Vector2D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Vector2D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Vector2D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Vector2D
    GENERIC,PUBLIC :: OPERATOR(+) => Add_Vector2D
    GENERIC,PUBLIC :: OPERATOR(-) => Subtract_Vector2D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Vector2D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Add_Vector2D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Subtract_Vector2D

  END TYPE Vector2D

  TYPE,PUBLIC :: Vector3D

    INTEGER :: N
    INTEGER :: M
    INTEGER :: nVar
    INTEGER :: nElem
    INTEGER :: controlType
    INTEGER :: targetType
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r6) :: interior
    TYPE(hfReal_r6) :: boundary

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Vector3D
    PROCEDURE,PUBLIC :: Free => Free_Vector3D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Vector3D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Vector3D
    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Vector3D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Vector3D

    GENERIC,PUBLIC :: Gradient => Gradient_Vector3D
    PROCEDURE,PRIVATE :: Gradient_Vector3D

    GENERIC,PUBLIC :: Divergence => Divergence_Vector3D
    PROCEDURE,PRIVATE :: Divergence_Vector3D

    GENERIC,PUBLIC :: Curl => Curl_Vector3D
    PROCEDURE,PRIVATE :: Curl_Vector3D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Vector3D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Vector3D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Vector3D
    GENERIC,PUBLIC :: OPERATOR(+) => Add_Vector3D
    GENERIC,PUBLIC :: OPERATOR(-) => Subtract_Vector3D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Vector3D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Add_Vector3D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Subtract_Vector3D

  END TYPE Vector3D

! ---------------------- Tensors ---------------------- !

  TYPE,PUBLIC :: Tensor2D

    INTEGER :: N
    INTEGER :: M
    INTEGER :: nVar
    INTEGER :: nElem
    INTEGER :: controlType
    INTEGER :: targetType
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r6) :: interior
    TYPE(hfReal_r6) :: boundary

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Tensor2D
    PROCEDURE,PUBLIC :: Free => Free_Tensor2D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Tensor2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Tensor2D
    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Tensor2D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Tensor2D

    PROCEDURE,PUBLIC :: Determinant => Determinant_Tensor2D

    GENERIC,PUBLIC :: Divergence => Divergence_Tensor2D
    PROCEDURE,PRIVATE :: Divergence_Tensor2D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Tensor2D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Tensor2D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Tensor2D
    GENERIC,PUBLIC :: OPERATOR(+) => Add_Tensor2D
    GENERIC,PUBLIC :: OPERATOR(-) => Subtract_Tensor2D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Tensor2D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Add_Tensor2D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Subtract_Tensor2D

  END TYPE Tensor2D

  TYPE,PUBLIC :: Tensor3D

    INTEGER :: N
    INTEGER :: M
    INTEGER :: nVar
    INTEGER :: nElem
    INTEGER :: controlType
    INTEGER :: targetType
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r7) :: interior
    TYPE(hfReal_r7) :: boundary

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Tensor3D
    PROCEDURE,PUBLIC :: Free => Free_Tensor3D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Tensor3D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Tensor3D
    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Tensor3D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Tensor3D

    PROCEDURE,PUBLIC :: Determinant => Determinant_Tensor3D

    GENERIC,PUBLIC :: Divergence => Divergence_Tensor3D
    PROCEDURE,PRIVATE :: Divergence_Tensor3D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Tensor3D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Tensor3D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Tensor3D
    GENERIC,PUBLIC :: OPERATOR(+) => Add_Tensor3D
    GENERIC,PUBLIC :: OPERATOR(-) => Subtract_Tensor3D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Tensor3D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Add_Tensor3D
    PROCEDURE,PRIVATE,PASS(SELFa) :: Subtract_Tensor3D

  END TYPE Tensor3D

  INTEGER, PARAMETER :: selfStrongForm = 0
  INTEGER, PARAMETER :: selfWeakDGForm = 1
  INTEGER, PARAMETER :: selfWeakCGForm = 2

#ifdef GPU
  INTERFACE
    SUBROUTINE Determinant_Tensor2D_gpu_wrapper(tensor_dev,detTensor_dev,N,nVar,nEl) &
      bind(c,name="Determinant_Tensor2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: tensor_dev,detTensor_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE Determinant_Tensor2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE Determinant_Tensor3D_gpu_wrapper(tensor_dev,detTensor_dev,N,nVar,nEl) &
      bind(c,name="Determinant_Tensor3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: tensor_dev,detTensor_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE Determinant_Tensor3D_gpu_wrapper
  END INTERFACE
#endif

CONTAINS

! -- Scalar1D -- !

  SUBROUTINE Init_Scalar1D(SELFStorage,N,quadratureType,M,targetNodeType,nVar,nElem)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(out) :: SELFStorage
    INTEGER,INTENT(in) :: N
    INTEGER,INTENT(in) :: quadratureType
    INTEGER,INTENT(in) :: M
    INTEGER,INTENT(in) :: targetNodeType
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % M = M
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    SELFStorage % controlType = quadratureType
    SELFStorage % targetType = targetNodeType

    CALL SELFStorage % interp % Init(N=N, &
                                     controlNodeType=quadratureType, &
                                     M=M, &
                                     targetNodeType=targetNodeType)

    CALL SELFStorage % interior % Alloc(loBound=(/0,1,1/), &
                                        upBound=(/N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/1,1,1/), &
                                        upBound=(/nVar,2,nElem/))

  END SUBROUTINE Init_Scalar1D

  SUBROUTINE Free_Scalar1D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Free()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

  END SUBROUTINE Free_Scalar1D

  SUBROUTINE UpdateHost_Scalar1D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateHost()
    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
#endif

  END SUBROUTINE UpdateHost_Scalar1D

  SUBROUTINE UpdateDevice_Scalar1D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateDevice()
    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
#endif

  END SUBROUTINE UpdateDevice_Scalar1D

  SUBROUTINE BoundaryInterp_Scalar1D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarBoundaryInterp_1D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarBoundaryInterp_1D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Scalar1D

  SUBROUTINE GridInterp_Scalar1D(SELFStorage,SELFout,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(in) :: SELFStorage
    TYPE(Scalar1D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarGridInterp_1D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarGridInterp_1D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Scalar1D

  SUBROUTINE Derivative_Scalar1D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(in) :: SELFStorage
    TYPE(Scalar1D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % Derivative_1D(SELFStorage % interior % deviceData, &
                                                SELFout % interior % deviceData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % Derivative_1D(SELFStorage % interior % hostData, &
                                                SELFout % interior % hostData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem)
    ENDIF

  END SUBROUTINE Derivative_Scalar1D

  FUNCTION AbsMaxInterior_Scalar1D(scalar) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Scalar1D) :: scalar
    REAL(prec) :: absMax(1:scalar % nVar)
    ! Local
    INTEGER :: iEl,iVar,i

    absMax = 0.0_prec
    DO iEl = 1,scalar % nElem
      DO iVar = 1,scalar % nVar
        DO i = 0,scalar % N
          absMax(iVar) = MAX(ABS(scalar % interior % hostData(i,iVar,iEl)),absMax(iVar))
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Scalar1D

  FUNCTION AbsMaxBoundary_Scalar1D(scalar) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Scalar1D) :: scalar
    REAL(prec) :: absMax(1:scalar % nVar,1:2)
    ! Local
    INTEGER :: iEl,iVar,iSide

    absMax = 0.0_prec
    DO iEl = 1,scalar % nElem
      DO iSide = 1,2
        DO iVar = 1,scalar % nVar
          absMax(iVar,iSide) = MAX(ABS(scalar % boundary % hostData(iVar,iSide,iEl)),absMax(iVar,iSide))
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Scalar1D

  SUBROUTINE Equals_Scalar1D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(inout) :: SELFOut
    TYPE(Scalar1D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Scalar1D

  FUNCTION Add_Scalar1D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(in) :: SELFa,SELFb
    TYPE(Scalar1D):: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData + &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData + &
                                    SELFb % boundary % hostData

  END FUNCTION Add_Scalar1D

  FUNCTION Subtract_Scalar1D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(in) :: SELFa,SELFb
    TYPE(Scalar1D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData - &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData - &
                                    SELFb % boundary % hostData

  END FUNCTION Subtract_Scalar1D

! -- Scalar2D -- !

  SUBROUTINE Init_Scalar2D(SELFStorage,N,quadratureType,M,targetNodeType,nVar,nElem)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(out) :: SELFStorage
    INTEGER,INTENT(in) :: N
    INTEGER,INTENT(in) :: quadratureType
    INTEGER,INTENT(in) :: M
    INTEGER,INTENT(in) :: targetNodeType
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % M = M
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    SELFStorage % controlType = quadratureType
    SELFStorage % targetType = targetNodeType

    CALL SELFStorage % interp % Init(N=N, &
                                     controlNodeType=quadratureType, &
                                     M=M, &
                                     targetNodeType=targetNodeType)

    CALL SELFStorage % interior % Alloc(loBound=(/0,0,1,1/), &
                                        upBound=(/N,N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/0,1,1,1/), &
                                        upBound=(/N,nVar,4,nElem/))

  END SUBROUTINE Init_Scalar2D

  SUBROUTINE Free_Scalar2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Free()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

  END SUBROUTINE Free_Scalar2D

  SUBROUTINE UpdateHost_Scalar2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateHost()
    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
#endif

  END SUBROUTINE UpdateHost_Scalar2D

  SUBROUTINE UpdateDevice_Scalar2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateDevice()
    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
#endif

  END SUBROUTINE UpdateDevice_Scalar2D

  SUBROUTINE BoundaryInterp_Scalar2D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarBoundaryInterp_2D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarBoundaryInterp_2D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Scalar2D

  SUBROUTINE GridInterp_Scalar2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(in) :: SELFStorage
    TYPE(Scalar2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarGridInterp_2D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarGridInterp_2D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Scalar2D

  SUBROUTINE Gradient_Scalar2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(in) :: SELFStorage
    TYPE(Vector2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarGradient_2D(SELFStorage % interior % deviceData, &
                                                    SELFout % interior % deviceData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarGradient_2D(SELFStorage % interior % hostData, &
                                                    SELFout % interior % hostData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    END IF

  END SUBROUTINE Gradient_Scalar2D

  FUNCTION AbsMaxInterior_Scalar2D(scalar) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Scalar2D) :: scalar
    REAL(prec) :: absMax(1:scalar % nVar)
    ! Local
    INTEGER :: iEl,iVar,i,j

    absMax = 0.0_prec
    DO iEl = 1,scalar % nElem
      DO iVar = 1,scalar % nVar
        DO j = 0,scalar % N
          DO i = 0,scalar % N
            absMax(iVar) = MAX(ABS(scalar % interior % hostData(i,j,iVar,iEl)),absMax(iVar))
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Scalar2D

  FUNCTION AbsMaxBoundary_Scalar2D(scalar) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Scalar2D) :: scalar
    REAL(prec) :: absMax(1:scalar % nVar,1:4)
    ! Local
    INTEGER :: iEl,iVar,i,iSide

    absMax = 0.0_prec
    DO iEl = 1,scalar % nElem
      DO iSide = 1,4
        DO iVar = 1,scalar % nVar
          DO i = 0,scalar % N
            absMax(iVar,iSide) = MAX(ABS(scalar % boundary % hostData(i,iVar,iSide,iEl)),absMax(iVar,iSide))
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Scalar2D

  SUBROUTINE Equals_Scalar2D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(inout) :: SELFOut
    TYPE(Scalar2D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Scalar2D

  FUNCTION Add_Scalar2D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(in) :: SELFa,SELFb
    TYPE(Scalar2D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData + &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData + &
                                    SELFb % boundary % hostData

  END FUNCTION Add_Scalar2D

  FUNCTION Subtract_Scalar2D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(in) :: SELFa,SELFb
    TYPE(Scalar2D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData - &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData - &
                                    SELFb % boundary % hostData

  END FUNCTION Subtract_Scalar2D

! -- Scalar3D -- !

  SUBROUTINE Init_Scalar3D(SELFStorage,N,quadratureType,M,targetNodeType,nVar,nElem)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(out) :: SELFStorage
    INTEGER,INTENT(in) :: N
    INTEGER,INTENT(in) :: quadratureType
    INTEGER,INTENT(in) :: M
    INTEGER,INTENT(in) :: targetNodeType
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % M = M
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    SELFStorage % controlType = quadratureType
    SELFStorage % targetType = targetNodeType

    CALL SELFStorage % interp % Init(N=N, &
                                     controlNodeType=quadratureType, &
                                     M=M, &
                                     targetNodeType=targetNodeType)

    CALL SELFStorage % interior % Alloc(loBound=(/0,0,0,1,1/), &
                                        upBound=(/N,N,N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/0,0,1,1,1/), &
                                        upBound=(/N,N,nVar,6,nElem/))

  END SUBROUTINE Init_Scalar3D

  SUBROUTINE Free_Scalar3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Free()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

  END SUBROUTINE Free_Scalar3D

  SUBROUTINE UpdateHost_Scalar3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateHost()
    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
#endif

  END SUBROUTINE UpdateHost_Scalar3D

  SUBROUTINE UpdateDevice_Scalar3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateDevice()
    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
#endif

  END SUBROUTINE UpdateDevice_Scalar3D

  SUBROUTINE BoundaryInterp_Scalar3D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarBoundaryInterp_3D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarBoundaryInterp_3D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Scalar3D

  SUBROUTINE GridInterp_Scalar3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(in) :: SELFStorage
    TYPE(Scalar3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarGridInterp_3D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarGridInterp_3D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Scalar3D

  SUBROUTINE Gradient_Scalar3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(in) :: SELFStorage
    TYPE(Vector3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarGradient_3D(SELFStorage % interior % deviceData, &
                                                    SELFout % interior % deviceData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarGradient_3D(SELFStorage % interior % hostData, &
                                                    SELFout % interior % hostData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    END IF

  END SUBROUTINE Gradient_Scalar3D

  SUBROUTINE Equals_Scalar3D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(inout) :: SELFOut
    TYPE(Scalar3D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Scalar3D

  FUNCTION AbsMaxInterior_Scalar3D(scalar) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Scalar3D) :: scalar
    REAL(prec) :: absMax(1:scalar % nVar)
    ! Local
    INTEGER :: iEl,iVar,i,j,k

    absMax = 0.0_prec
    DO iEl = 1,scalar % nElem
      DO iVar = 1,scalar % nVar
        DO k = 0,scalar % N
          DO j = 0,scalar % N
            DO i = 0,scalar % N
              absMax(iVar) = MAX(ABS(scalar % interior % hostData(i,j,k,iVar,iEl)),absMax(iVar))
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Scalar3D

  FUNCTION AbsMaxBoundary_Scalar3D(scalar) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Scalar3D) :: scalar
    REAL(prec) :: absMax(1:scalar % nVar,1:6)
    ! Local
    INTEGER :: iEl,iVar,i,j,iSide

    absMax = 0.0_prec
    DO iEl = 1,scalar % nElem
      DO iSide = 1,6
        DO iVar = 1,scalar % nVar
          DO j = 0,scalar % N
            DO i = 0,scalar % N
              absMax(iVar,iSide) = MAX(ABS(scalar % boundary % hostData(i,j,iVar,iSide,iEl)),absMax(iVar,iSide))
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Scalar3D

  FUNCTION Add_Scalar3D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(in) :: SELFa,SELFb
    TYPE(Scalar3D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData + &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData + &
                                    SELFb % boundary % hostData

  END FUNCTION Add_Scalar3D

  FUNCTION Subtract_Scalar3D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(in) :: SELFa,SELFb
    TYPE(Scalar3D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData - &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData - &
                                    SELFb % boundary % hostData

  END FUNCTION Subtract_Scalar3D

! -- Vector2D -- !

  SUBROUTINE Init_Vector2D(SELFStorage,N,quadratureType,M,targetNodeType,nVar,nElem)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(out) :: SELFStorage
    INTEGER,INTENT(in) :: N
    INTEGER,INTENT(in) :: quadratureType
    INTEGER,INTENT(in) :: M
    INTEGER,INTENT(in) :: targetNodeType
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % M = M
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    SELFStorage % controlType = quadratureType
    SELFStorage % targetType = targetNodeType

    CALL SELFStorage % interp % Init(N=N, &
                                     controlNodeType=quadratureType, &
                                     M=M, &
                                     targetNodeType=targetNodeType)

    CALL SELFStorage % interior % Alloc(loBound=(/1,0,0,1,1/), &
                                        upBound=(/2,N,N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/1,0,1,1,1/), &
                                        upBound=(/2,N,nVar,4,nElem/))

  END SUBROUTINE Init_Vector2D

  SUBROUTINE Free_Vector2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Free()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

  END SUBROUTINE Free_Vector2D

  SUBROUTINE UpdateHost_Vector2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateHost()
    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
#endif

  END SUBROUTINE UpdateHost_Vector2D

  SUBROUTINE UpdateDevice_Vector2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateDevice()
    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
#endif

  END SUBROUTINE UpdateDevice_Vector2D

  SUBROUTINE BoundaryInterp_Vector2D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorBoundaryInterp_2D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorBoundaryInterp_2D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Vector2D

  SUBROUTINE GridInterp_Vector2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(in) :: SELFStorage
    TYPE(Vector2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorGridInterp_2D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorGridInterp_2D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Vector2D

  SUBROUTINE Gradient_Vector2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(in) :: SELFStorage
    TYPE(Tensor2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorGradient_2D(SELFStorage % interior % deviceData, &
                                                    SELFout % interior % deviceData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorGradient_2D(SELFStorage % interior % hostData, &
                                                    SELFout % interior % hostData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    END IF

  END SUBROUTINE Gradient_Vector2D

  SUBROUTINE Divergence_Vector2D(SELFStorage,SELFOut,dForm,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(in) :: SELFStorage
    TYPE(Scalar2D),INTENT(inout) :: SELFOut
    INTEGER,INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

    IF (dForm == selfWeakDGForm) THEN

      IF (gpuAccel) THEN
        CALL SELFStorage % interp % VectorDGDivergence_2D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFout % interior % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
      ELSE
        CALL SELFStorage % interp % VectorDGDivergence_2D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFout % interior % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
      END IF

    ELSE IF (dForm == selfStrongForm) THEN

      IF (gpuAccel) THEN
        CALL SELFStorage % interp % VectorDivergence_2D(SELFStorage % interior % deviceData, &
                                                        SELFout % interior % deviceData, &
                                                        SELFStorage % nVar, &
                                                        SELFStorage % nElem)
      ELSE
        CALL SELFStorage % interp % VectorDivergence_2D(SELFStorage % interior % hostData, &
                                                        SELFout % interior % hostData, &
                                                        SELFStorage % nVar, &
                                                        SELFStorage % nElem)
      END IF

    END IF

  END SUBROUTINE Divergence_Vector2D

  SUBROUTINE Curl_Vector2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(in) :: SELFStorage
    TYPE(Scalar2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorCurl_2D(SELFStorage % interior % deviceData, &
                                                SELFout % interior % deviceData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorCurl_2D(SELFStorage % interior % hostData, &
                                                SELFout % interior % hostData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem)
    END IF

  END SUBROUTINE Curl_Vector2D

  SUBROUTINE Equals_Vector2D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(inout) :: SELFOut
    TYPE(Vector2D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Vector2D

  FUNCTION AbsMaxInterior_Vector2D(vector) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Vector2D) :: vector
    REAL(prec) :: absMax(1:vector % nVar)
    ! Local
    INTEGER :: iEl,iVar,i,j,iDir

    absMax = 0.0_prec
    DO iEl = 1,vector % nElem
      DO iVar = 1,vector % nVar
        DO j = 0,vector % N
          DO i = 0,vector % N
            DO iDir = 1,2
              absMax(iVar) = MAX(ABS(vector % interior % hostData(iDir,i,j,iVar,iEl)),absMax(iVar))
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Vector2D

  FUNCTION AbsMaxBoundary_Vector2D(vector) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Vector2D) :: vector
    REAL(prec) :: absMax(1:vector % nVar,1:4)
    ! Local
    INTEGER :: iEl,iVar,i,iDir,iSide

    absMax = 0.0_prec
    DO iEl = 1,vector % nElem
      DO iSide = 1,4
        DO iVar = 1,vector % nVar
          DO i = 0,vector % N
            DO iDir = 1,2
              absMax(iVar,iSide) = MAX(ABS(vector % boundary % hostData(iDir,i,iVar,iSide,iEl)),absMax(iVar,iSide))
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Vector2D

  FUNCTION Add_Vector2D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(in) :: SELFa,SELFb
    TYPE(Vector2D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData + &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData + &
                                    SELFb % boundary % hostData

  END FUNCTION Add_Vector2D

  FUNCTION Subtract_Vector2D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(in) :: SELFa,SELFb
    TYPE(Vector2D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData - &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData - &
                                    SELFb % boundary % hostData

  END FUNCTION Subtract_Vector2D

! -- Vector3D -- !

  SUBROUTINE Init_Vector3D(SELFStorage,N,quadratureType,M,targetNodeType,nVar,nElem)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(out) :: SELFStorage
    INTEGER,INTENT(in) :: N
    INTEGER,INTENT(in) :: quadratureType
    INTEGER,INTENT(in) :: M
    INTEGER,INTENT(in) :: targetNodeType
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % M = M
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    SELFStorage % controlType = quadratureType
    SELFStorage % targetType = targetNodeType

    CALL SELFStorage % interp % Init(N=N, &
                                     controlNodeType=quadratureType, &
                                     M=M, &
                                     targetNodeType=targetNodeType)

    CALL SELFStorage % interior % Alloc(loBound=(/1,0,0,0,1,1/), &
                                        upBound=(/3,N,N,N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/1,0,0,1,1,1/), &
                                        upBound=(/3,N,N,nVar,6,nElem/))

  END SUBROUTINE Init_Vector3D

  SUBROUTINE Free_Vector3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Free()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

  END SUBROUTINE Free_Vector3D

  SUBROUTINE UpdateHost_Vector3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateHost()
    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
#endif

  END SUBROUTINE UpdateHost_Vector3D

  SUBROUTINE UpdateDevice_Vector3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateDevice()
    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
#endif

  END SUBROUTINE UpdateDevice_Vector3D

  SUBROUTINE BoundaryInterp_Vector3D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorBoundaryInterp_3D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorBoundaryInterp_3D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Vector3D

  SUBROUTINE GridInterp_Vector3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(in) :: SELFStorage
    TYPE(Vector3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorGridInterp_3D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorGridInterp_3D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Vector3D

  SUBROUTINE Gradient_Vector3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(in) :: SELFStorage
    TYPE(Tensor3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorGradient_3D(SELFStorage % interior % deviceData, &
                                                    SELFout % interior % deviceData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorGradient_3D(SELFStorage % interior % hostData, &
                                                    SELFout % interior % hostData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    END IF

  END SUBROUTINE Gradient_Vector3D

  SUBROUTINE Divergence_Vector3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(in) :: SELFStorage
    TYPE(Scalar3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorDivergence_3D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorDivergence_3D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE Divergence_Vector3D

  SUBROUTINE Curl_Vector3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(in) :: SELFStorage
    TYPE(Vector3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorCurl_3D(SELFStorage % interior % deviceData, &
                                                SELFout % interior % deviceData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorCurl_3D(SELFStorage % interior % hostData, &
                                                SELFout % interior % hostData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem)
    END IF

  END SUBROUTINE Curl_Vector3D

  FUNCTION AbsMaxInterior_Vector3D(vector) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Vector3D) :: vector
    REAL(prec) :: absMax(1:vector % nVar)
    ! Local
    INTEGER :: iEl,iVar,i,j,k,iDir

    absMax = 0.0_prec
    DO iEl = 1,vector % nElem
      DO iVar = 1,vector % nVar
        DO k = 0,vector % N
          DO j = 0,vector % N
            DO i = 0,vector % N
              DO iDir = 1,3
                absMax(iVar) = MAX(ABS(vector % interior % hostData(iDir,i,j,k,iVar,iEl)),absMax(iVar))
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Vector3D

  FUNCTION AbsMaxBoundary_Vector3D(vector) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Vector3D) :: vector
    REAL(prec) :: absMax(1:vector % nVar,1:6)
    ! Local
    INTEGER :: iEl,iVar,i,j,iSide,iDir

    absMax = 0.0_prec
    DO iEl = 1,vector % nElem
      DO iSide = 1,6
        DO iVar = 1,vector % nVar
          DO j = 0,vector % N
            DO i = 0,vector % N
              DO iDir = 1,3
                absMax(iVar,iSide) = MAX(ABS(vector % boundary % hostData(iDir,i,j,iVar,iSide,iEl)),absMax(iVar,iSide))
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Vector3D

  SUBROUTINE Equals_Vector3D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(inout) :: SELFOut
    TYPE(Vector3D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Vector3D

  FUNCTION Add_Vector3D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(in) :: SELFa,SELFb
    TYPE(Vector3D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData + &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData + &
                                    SELFb % boundary % hostData

  END FUNCTION Add_Vector3D

  FUNCTION Subtract_Vector3D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(in) :: SELFa,SELFb
    TYPE(Vector3D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData - &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData - &
                                    SELFb % boundary % hostData

  END FUNCTION Subtract_Vector3D

! -- Tensor2D -- !

  SUBROUTINE Init_Tensor2D(SELFStorage,N,quadratureType,M,targetNodeType,nVar,nElem)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(out) :: SELFStorage
    INTEGER,INTENT(in) :: N
    INTEGER,INTENT(in) :: quadratureType
    INTEGER,INTENT(in) :: M
    INTEGER,INTENT(in) :: targetNodeType
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % M = M
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    SELFStorage % controlType = quadratureType
    SELFStorage % targetType = targetNodeType

    CALL SELFStorage % interp % Init(N=N, &
                                     controlNodeType=quadratureType, &
                                     M=M, &
                                     targetNodeType=targetNodeType)

    CALL SELFStorage % interior % Alloc(loBound=(/1,1,0,0,1,1/), &
                                        upBound=(/2,2,N,N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/1,1,0,1,1,1/), &
                                        upBound=(/2,2,N,nVar,4,nElem/))

  END SUBROUTINE Init_Tensor2D

  SUBROUTINE Free_Tensor2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Free()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

  END SUBROUTINE Free_Tensor2D

  SUBROUTINE UpdateHost_Tensor2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateHost()
    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
#endif

  END SUBROUTINE UpdateHost_Tensor2D

  SUBROUTINE UpdateDevice_Tensor2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateDevice()
    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
#endif

  END SUBROUTINE UpdateDevice_Tensor2D

  SUBROUTINE BoundaryInterp_Tensor2D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % TensorBoundaryInterp_2D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % TensorBoundaryInterp_2D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Tensor2D

  SUBROUTINE GridInterp_Tensor2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(in) :: SELFStorage
    TYPE(Tensor2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % TensorGridInterp_2D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % TensorGridInterp_2D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Tensor2D

  SUBROUTINE Divergence_Tensor2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(in) :: SELFStorage
    TYPE(Vector2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % TensorDivergence_2D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % TensorDivergence_2D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE Divergence_Tensor2D

  SUBROUTINE Determinant_Tensor2D(SELFStorage,SELFout,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(in) :: SELFStorage
    TYPE(Scalar2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iEl,iVar,i,j

    IF (gpuAccel) THEN

      CALL Determinant_Tensor2D_gpu_wrapper(SELFStorage % interior % deviceData, &
                                            SELFOut % interior % deviceData, &
                                            SELFStorage % N, &
                                            SELFStorage % nVar, &
                                            SELFStorage % nElem)

    ELSE

      DO iEl = 1,SELFStorage % nElem
        DO iVar = 1,SELFStorage % nVar
          DO j = 0,SELFStorage % N
            DO i = 0,SELFStorage % N
  
              SELFOut % interior % hostData(i,j,iVar,iEl) = SELFStorage % interior % hostData(1,1,i,j,iVar,iEl)* &
                                                            SELFStorage % interior % hostData(2,2,i,j,iVar,iEl) - &
                                                            SELFStorage % interior % hostData(1,2,i,j,iVar,iEl)* &
                                                            SELFStorage % interior % hostData(2,1,i,j,iVar,iEl)
  
            END DO
          END DO
        END DO
      END DO

    ENDIF

  END SUBROUTINE Determinant_Tensor2D

  FUNCTION AbsMaxInterior_Tensor2D(tensor) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Tensor2D) :: tensor
    REAL(prec) :: absMax(1:tensor % nVar)
    ! Local
    INTEGER :: iEl,iVar,i,j,row,col

    absMax = 0.0_prec
    DO iEl = 1,tensor % nElem
      DO iVar = 1,tensor % nVar
        DO j = 0,tensor % N
          DO i = 0,tensor % N
            DO col = 1,2
              DO row = 1,2
                absMax(iVar) = MAX(ABS(tensor % interior % hostData(row,col,i,j,iVar,iEl)),absMax(iVar))
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Tensor2D

  FUNCTION AbsMaxBoundary_Tensor2D(tensor) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Tensor2D) :: tensor
    REAL(prec) :: absMax(1:tensor % nVar,1:4)
    ! Local
    INTEGER :: iEl,iVar,i,iSide,row,col

    absMax = 0.0_prec
    DO iEl = 1,tensor % nElem
      DO iSide = 1,4
        DO iVar = 1,tensor % nVar
          DO i = 0,tensor % N
            DO col = 1,2
              DO row = 1,2
                absMax(iVar,iSide) = MAX(ABS(tensor % boundary % hostData(row,col,i,iVar,iSide,iEl)),absMax(iVar,iSide))
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Tensor2D

  SUBROUTINE Equals_Tensor2D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(inout) :: SELFOut
    TYPE(Tensor2D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Tensor2D

  FUNCTION Add_Tensor2D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(in) :: SELFa,SELFb
    TYPE(Tensor2D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData + &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData + &
                                    SELFb % boundary % hostData

  END FUNCTION Add_Tensor2D

  FUNCTION Subtract_Tensor2D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(in) :: SELFa,SELFb
    TYPE(Tensor2D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData - &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData - &
                                    SELFb % boundary % hostData

  END FUNCTION Subtract_Tensor2D

! -- Tensor3D -- !

  SUBROUTINE Init_Tensor3D(SELFStorage,N,quadratureType,M,targetNodeType,nVar,nElem)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(out) :: SELFStorage
    INTEGER,INTENT(in) :: N
    INTEGER,INTENT(in) :: quadratureType
    INTEGER,INTENT(in) :: M
    INTEGER,INTENT(in) :: targetNodeType
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % M = M
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    SELFStorage % controlType = quadratureType
    SELFStorage % targetType = targetNodeType

    CALL SELFStorage % interp % Init(N=N, &
                                     controlNodeType=quadratureType, &
                                     M=M, &
                                     targetNodeType=targetNodeType)

    CALL SELFStorage % interior % Alloc(loBound=(/1,1,0,0,0,1,1/), &
                                        upBound=(/3,3,N,N,N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/1,1,0,0,1,1,1/), &
                                        upBound=(/3,3,N,N,nVar,6,nElem/))

  END SUBROUTINE Init_Tensor3D

  SUBROUTINE Free_Tensor3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Free()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

  END SUBROUTINE Free_Tensor3D

  SUBROUTINE UpdateHost_Tensor3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateHost()
    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
#endif

  END SUBROUTINE UpdateHost_Tensor3D

  SUBROUTINE UpdateDevice_Tensor3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(inout) :: SELFStorage

#ifdef GPU
    CALL SELFStorage % interp % UpdateDevice()
    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
#endif

  END SUBROUTINE UpdateDevice_Tensor3D

  SUBROUTINE BoundaryInterp_Tensor3D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % TensorBoundaryInterp_3D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % TensorBoundaryInterp_3D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Tensor3D

  SUBROUTINE GridInterp_Tensor3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(in) :: SELFStorage
    TYPE(Tensor3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % TensorGridInterp_3D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % TensorGridInterp_3D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Tensor3D

  SUBROUTINE Divergence_Tensor3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(in) :: SELFStorage
    TYPE(Vector3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % TensorDivergence_3D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % TensorDivergence_3D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE Divergence_Tensor3D

  SUBROUTINE Determinant_Tensor3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(in) :: SELFStorage
    TYPE(Scalar3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iEl,iVar,i,j,k

    IF (gpuAccel) THEN

      CALL Determinant_Tensor3D_gpu_wrapper(SELFStorage % interior % deviceData, &
                                            SELFOut % interior % deviceData, &
                                            SELFStorage % N, &
                                            SELFStorage % nVar, &
                                            SELFStorage % nElem)

    ELSE

      DO iEl = 1,SELFStorage % nElem
        DO iVar = 1,SELFStorage % nVar
          DO k = 0,SELFStorage % N
            DO j = 0,SELFStorage % N
              DO i = 0,SELFStorage % N

                SELFOut % interior % hostData(i,j,k,iVar,iEl) = &
                  SELFStorage % interior % hostData(1,1,i,j,k,iVar,iEl)* &
                  (SELFStorage % interior % hostData(2,2,i,j,k,iVar,iEl)* &
                   SELFStorage % interior % hostData(3,3,i,j,k,iVar,iEl) - &
                   SELFStorage % interior % hostData(2,3,i,j,k,iVar,iEl)* &
                   SELFStorage % interior % hostData(3,2,i,j,k,iVar,iEl)) - &
                  SELFStorage % interior % hostData(2,1,i,j,k,iVar,iEl)* &
                  (SELFStorage % interior % hostData(1,2,i,j,k,iVar,iEl)* &
                   SELFStorage % interior % hostData(3,3,i,j,k,iVar,iEl) - &
                   SELFStorage % interior % hostData(1,3,i,j,k,iVar,iEl)* &
                   SELFStorage % interior % hostData(3,2,i,j,k,iVar,iEl)) + &
                  SELFStorage % interior % hostData(3,1,i,j,k,iVar,iEl)* &
                  (SELFStorage % interior % hostData(1,2,i,j,k,iVar,iEl)* &
                   SELFStorage % interior % hostData(2,3,i,j,k,iVar,iEl) - &
                   SELFStorage % interior % hostData(1,3,i,j,k,iVar,iEl)* &
                   SELFStorage % interior % hostData(2,2,i,j,k,iVar,iEl))

              END DO
            END DO
          END DO
        END DO
      END DO

    ENDIF

  END SUBROUTINE Determinant_Tensor3D

  FUNCTION AbsMaxInterior_Tensor3D(tensor) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Tensor3D) :: tensor
    REAL(prec) :: absMax(1:tensor % nVar)
    ! Local
    INTEGER :: iEl,iVar,i,j,k,row,col

    absMax = 0.0_prec
    DO iEl = 1,tensor % nElem
      DO iVar = 1,tensor % nVar
        DO k = 0,tensor % N
          DO j = 0,tensor % N
            DO i = 0,tensor % N
              DO col = 1,3
                DO row = 1,3
                  absMax(iVar) = MAX(ABS(tensor % interior % hostData(row,col,i,j,k,iVar,iEl)),absMax(iVar))
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Tensor3D

  FUNCTION AbsMaxBoundary_Tensor3D(tensor) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Tensor3D) :: tensor
    REAL(prec) :: absMax(1:tensor % nVar,1:6)
    ! Local
    INTEGER :: iEl,iVar,i,j,iSide,row,col

    absMax = 0.0_prec
    DO iEl = 1,tensor % nElem
      DO iSide = 1,6
        DO iVar = 1,tensor % nVar
          DO j = 0,tensor % N
            DO i = 0,tensor % N
              DO col = 1,3
                DO row = 1,3
                  absMax(iVar,iSide) = MAX(ABS(tensor % boundary % hostData(row,col,i,j,iVar,iSide,iEl)),absMax(iVar,iSide))
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Tensor3D

  SUBROUTINE Equals_Tensor3D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(inout) :: SELFOut
    TYPE(Tensor3D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Tensor3D

  FUNCTION Add_Tensor3D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(in) :: SELFa,SELFb
    TYPE(Tensor3D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData + &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData + &
                                    SELFb % boundary % hostData

  END FUNCTION Add_Tensor3D

  FUNCTION Subtract_Tensor3D(SELFa,SELFb) RESULT(SELFOut)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(in) :: SELFa,SELFb
    TYPE(Tensor3D) :: SELFOut

    CALL SELFOut % Init(SELFa % N, &
                        SELFa % controlType, &
                        SELFa % M, &
                        SELFa % targetType, &
                        SELFa % nVar, &
                        SELFa % nElem)

    SELFOut % interior % hostData = SELFa % interior % hostData - &
                                    SELFb % interior % hostData

    SELFOut % boundary % hostData = SELFa % boundary % hostData - &
                                    SELFb % boundary % hostData

  END FUNCTION Subtract_Tensor3D

END MODULE SELF_Data
