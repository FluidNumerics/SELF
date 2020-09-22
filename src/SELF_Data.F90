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

USE hipfort
USE ISO_C_BINDING


IMPLICIT NONE

! ---------------------- Scalars ---------------------- !
  TYPE, PUBLIC :: Scalar1D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r3) :: interior
    TYPE(hfReal_r3) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_Scalar1D
      PROCEDURE, PUBLIC :: Trash => Trash_Scalar1D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_Scalar1D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_Scalar1D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_Scalar1D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_Scalar1D
      PROCEDURE, PUBLIC :: Derivative => Derivative_Scalar1D

  END TYPE Scalar1D

  TYPE, PUBLIC :: Scalar2D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r4) :: interior
    TYPE(hfReal_r4) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_Scalar2D
      PROCEDURE, PUBLIC :: Trash => Trash_Scalar2D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_Scalar2D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_Scalar2D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_Scalar2D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_Scalar2D

      GENERIC, PUBLIC :: Gradient => Gradient_Scalar2D
      PROCEDURE, PRIVATE :: Gradient_Scalar2D

  END TYPE Scalar2D

  TYPE, PUBLIC :: Scalar3D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r5) :: interior
    TYPE(hfReal_r5) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_Scalar3D
      PROCEDURE, PUBLIC :: Trash => Trash_Scalar3D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_Scalar3D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_Scalar3D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_Scalar3D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_Scalar3D

      GENERIC, PUBLIC :: Gradient => Gradient_Scalar3D
      PROCEDURE, PRIVATE :: Gradient_Scalar3D

  END TYPE Scalar3D

! ---------------------- Vectors ---------------------- !

  TYPE, PUBLIC :: Vector2D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r5) :: interior
    TYPE(hfReal_r5) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_Vector2D
      PROCEDURE, PUBLIC :: Trash => Trash_Vector2D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_Vector2D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_Vector2D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_Vector2D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_Vector2D

      GENERIC, PUBLIC :: Gradient => Gradient_Vector2D
      PROCEDURE, PRIVATE :: Gradient_Vector2D

      GENERIC, PUBLIC :: Divergence => Divergence_Vector2D
      PROCEDURE, PRIVATE :: Divergence_Vector2D

      GENERIC, PUBLIC :: Curl => Curl_Vector2D
      PROCEDURE, PRIVATE :: Curl_Vector2D


  END TYPE Vector2D

  TYPE, PUBLIC :: Vector3D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r6) :: interior
    TYPE(hfReal_r6) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_Vector3D
      PROCEDURE, PUBLIC :: Trash => Trash_Vector3D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_Vector3D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_Vector3D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_Vector3D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_Vector3D

      GENERIC, PUBLIC :: Gradient => Gradient_Vector3D
      PROCEDURE, PRIVATE :: Gradient_Vector3D

      GENERIC, PUBLIC :: Divergence => Divergence_Vector3D
      PROCEDURE, PRIVATE :: Divergence_Vector3D

      GENERIC, PUBLIC :: Curl => Curl_Vector3D
      PROCEDURE, PRIVATE :: Curl_Vector3D


  END TYPE Vector3D

! ---------------------- Tensors ---------------------- !

  TYPE, PUBLIC :: Tensor2D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r6) :: interior
    TYPE(hfReal_r6) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_Tensor2D
      PROCEDURE, PUBLIC :: Trash => Trash_Tensor2D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_Tensor2D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_Tensor2D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_Tensor2D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_Tensor2D 

      PROCEDURE, PUBLIC :: Determinant => Determinant_Tensor2D

      GENERIC, PUBLIC :: Divergence => Divergence_Tensor2D
      PROCEDURE, PRIVATE :: Divergence_Tensor2D

  END TYPE Tensor2D

  TYPE, PUBLIC :: Tensor3D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r7) :: interior
    TYPE(hfReal_r7) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_Tensor3D
      PROCEDURE, PUBLIC :: Trash => Trash_Tensor3D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_Tensor3D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_Tensor3D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_Tensor3D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_Tensor3D 

      PROCEDURE, PUBLIC :: Determinant => Determinant_Tensor3D

      GENERIC, PUBLIC :: Divergence => Divergence_Tensor3D
      PROCEDURE, PRIVATE :: Divergence_Tensor3D

  END TYPE Tensor3D

  INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE Equals_Scalar1D
    MODULE PROCEDURE Equals_Scalar2D
    MODULE PROCEDURE Equals_Scalar3D
    MODULE PROCEDURE Equals_Vector2D
    MODULE PROCEDURE Equals_Vector3D
    MODULE PROCEDURE Equals_Tensor2D
    MODULE PROCEDURE Equals_Tensor3D
  END INTERFACE

!  INTERFACE OPERATOR (+)
!    MODULE PROCEDURE Add_Scalar1D
!    MODULE PROCEDURE Add_Scalar2D
!    MODULE PROCEDURE Add_Scalar3D
!    MODULE PROCEDURE Add_Vector2D
!    MODULE PROCEDURE Add_Vector3D
!    MODULE PROCEDURE Add_Tensor2D
!  END INTERFACE
!
!  INTERFACE OPERATOR (-)
!    MODULE PROCEDURE Subtract_Scalar1D
!    MODULE PROCEDURE Subtract_Scalar2D
!    MODULE PROCEDURE Subtract_Scalar3D
!    MODULE PROCEDURE Subtract_Vector2D
!    MODULE PROCEDURE Subtract_Vector3D
!    MODULE PROCEDURE Subtract_Tensor2D
!  END INTERFACE

CONTAINS

! -- Scalar1D -- !

SUBROUTINE Build_Scalar1D( SELFStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(Scalar1D), INTENT(out) :: SELFStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem

    CALL SELFStorage % interp % Build(N = N, &
                                     controlNodeType = quadratureType, &
                                     M = M, &
                                     targetNodeType = targetNodeType)

    CALL SELFStorage % interior % Alloc(loBound = (/0, 1, 1/),&
                                       upBound = (/N, nVar, nElem/))

    CALL SELFStorage % boundary % Alloc(loBound = (/1, 1, 1/),&
                                       upBound = (/nVar, 2, nElem/))


END SUBROUTINE Build_Scalar1D

SUBROUTINE Trash_Scalar1D( SELFStorage ) 
  IMPLICIT NONE
  CLASS(Scalar1D), INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Trash()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

END SUBROUTINE Trash_Scalar1D

#ifdef GPU
SUBROUTINE UpdateHost_Scalar1D( SELFStorage )
  IMPLICIT NONE
  CLASS(Scalar1D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateHost()
     CALL SELFStorage % interior % UpdateHost()
     CALL SELFStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_Scalar1D

SUBROUTINE UpdateDevice_Scalar1D( SELFStorage )
  IMPLICIT NONE
  CLASS(Scalar1D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateDevice()
     CALL SELFStorage % interior % UpdateDevice()
     CALL SELFStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_Scalar1D
#endif

FUNCTION BoundaryInterp_Scalar1D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Scalar1D) :: SELFStorage
  TYPE(Scalar1D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % ScalarBoundaryInterp_1D( SELFStorage % interior % deviceData, &
                                                          SELFout % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % ScalarBoundaryInterp_1D( SELFStorage % interior % hostData, &
                                                          SELFout % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_Scalar1D

FUNCTION GridInterp_Scalar1D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Scalar1D) :: SELFStorage
  TYPE(Scalar1D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % ScalarGridInterp_1D( SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % ScalarGridInterp_1D( SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_Scalar1D

FUNCTION Derivative_Scalar1D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Scalar1D) :: SELFStorage
  TYPE(Scalar1D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % Derivative_1D( SELFStorage % interior % deviceData, &
                                                SELFout % interior % deviceData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % Derivative_1D( SELFStorage % interior % hostData, &
                                                SELFout % interior % hostData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem )  
    ENDIF

END FUNCTION Derivative_Scalar1D

SUBROUTINE Equals_Scalar1D( SELFOut, SELFStorage )
  IMPLICIT NONE
  TYPE(Scalar1D), INTENT(out) :: SELFOut
  TYPE(Scalar1D), INTENT(in) :: SELFStorage

    SELFOut % interior % hostData = SELFStorage % interior % hostData
    SELFOut % boundary % hostData = SELFStorage % boundary % hostData
  
END SUBROUTINE Equals_Scalar1D

! -- Scalar2D -- !

SUBROUTINE Build_Scalar2D( SELFStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(Scalar2D), INTENT(out) :: SELFStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem

    CALL SELFStorage % interp % Build( N = N, &
                                      controlNodeType = quadratureType, &
                                      M = M, &
                                      targetNodeType = targetNodeType )

    CALL SELFStorage % interior % Alloc(loBound = (/0, 0, 1, 1/),&
                                       upBound = (/N, N, nVar, nElem/))

    CALL SELFStorage % boundary % Alloc(loBound = (/0, 1, 1, 1/),&
                                       upBound = (/N, nVar, 4, nElem/))


END SUBROUTINE Build_Scalar2D

SUBROUTINE Trash_Scalar2D( SELFStorage ) 
  IMPLICIT NONE
  CLASS(Scalar2D), INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Trash()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

END SUBROUTINE Trash_Scalar2D

#ifdef GPU
SUBROUTINE UpdateHost_Scalar2D( SELFStorage )
  IMPLICIT NONE
  CLASS(Scalar2D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateHost()
     CALL SELFStorage % interior % UpdateHost()
     CALL SELFStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_Scalar2D

SUBROUTINE UpdateDevice_Scalar2D( SELFStorage )
  IMPLICIT NONE
  CLASS(Scalar2D), INTENT(inout) :: SELFStorage
          
     CALL SELFStorage % interp % UpdateDevice()
     CALL SELFStorage % interior % UpdateDevice()
     CALL SELFStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_Scalar2D
#endif

FUNCTION BoundaryInterp_Scalar2D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Scalar2D) :: SELFStorage
  TYPE(Scalar2D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % ScalarBoundaryInterp_2D( SELFStorage % interior % deviceData, &
                                                          SELFout % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % ScalarBoundaryInterp_2D( SELFStorage % interior % hostData, &
                                                          SELFout % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_Scalar2D

FUNCTION GridInterp_Scalar2D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Scalar2D) :: SELFStorage
  TYPE(Scalar2D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % ScalarGridInterp_2D( SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % ScalarGridInterp_2D( SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_Scalar2D

FUNCTION Gradient_Scalar2D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Scalar2D) :: SELFStorage
  TYPE(Vector2D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % ScalarGradient_2D( SELFStorage % interior % deviceData, &
                                                    SELFout % interior % deviceData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % ScalarGradient_2D( SELFStorage % interior % hostData, &
                                                    SELFout % interior % hostData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem )  
    ENDIF

END FUNCTION Gradient_Scalar2D

SUBROUTINE Equals_Scalar2D( SELFOut, SELFStorage )
  IMPLICIT NONE
  TYPE(Scalar2D), INTENT(out) :: SELFOut
  TYPE(Scalar2D), INTENT(in) :: SELFStorage

    SELFOut % interior % hostData = SELFStorage % interior % hostData
    SELFOut % boundary % hostData = SELFStorage % boundary % hostData
  
END SUBROUTINE Equals_Scalar2D

! -- Scalar3D -- !

SUBROUTINE Build_Scalar3D( SELFStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(Scalar3D), INTENT(out) :: SELFStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem

    CALL SELFStorage % interp % Build( N = N, &
                                      controlNodeType = quadratureType, &
                                      M = M, &
                                      targetNodeType = targetNodeType )

    CALL SELFStorage % interior % Alloc(loBound = (/0, 0, 0, 1, 1/),&
                                       upBound = (/N, N, N, nVar, nElem/))

    CALL SELFStorage % boundary % Alloc(loBound = (/0, 0, 1, 1, 1/),&
                                       upBound = (/N, N, nVar, 6, nElem/))

END SUBROUTINE Build_Scalar3D

SUBROUTINE Trash_Scalar3D( SELFStorage ) 
  IMPLICIT NONE
  CLASS(Scalar3D), INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Trash()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

END SUBROUTINE Trash_Scalar3D

#ifdef GPU
SUBROUTINE UpdateHost_Scalar3D( SELFStorage )
  IMPLICIT NONE
  CLASS(Scalar3D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateHost()
     CALL SELFStorage % interior % UpdateHost()
     CALL SELFStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_Scalar3D

SUBROUTINE UpdateDevice_Scalar3D( SELFStorage )
  IMPLICIT NONE
  CLASS(Scalar3D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateDevice()
     CALL SELFStorage % interior % UpdateDevice()
     CALL SELFStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_Scalar3D
#endif

FUNCTION BoundaryInterp_Scalar3D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Scalar3D) :: SELFStorage
  TYPE(Scalar3D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % ScalarBoundaryInterp_3D( SELFStorage % interior % deviceData, &
                                                          SELFout % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % ScalarBoundaryInterp_3D( SELFStorage % interior % hostData, &
                                                          SELFout % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_Scalar3D

FUNCTION GridInterp_Scalar3D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Scalar3D) :: SELFStorage
  TYPE(Scalar3D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % ScalarGridInterp_3D( SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % ScalarGridInterp_3D( SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_Scalar3D

FUNCTION Gradient_Scalar3D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Scalar3D) :: SELFStorage
  TYPE(Vector3D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % ScalarGradient_3D( SELFStorage % interior % deviceData, &
                                                    SELFout % interior % deviceData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % ScalarGradient_3D( SELFStorage % interior % hostData, &
                                                    SELFout % interior % hostData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem )  
    ENDIF

END FUNCTION Gradient_Scalar3D

SUBROUTINE Equals_Scalar3D( SELFOut, SELFStorage )
  IMPLICIT NONE
  TYPE(Scalar3D), INTENT(out) :: SELFOut
  TYPE(Scalar3D), INTENT(in) :: SELFStorage

    SELFOut % interior % hostData = SELFStorage % interior % hostData
    SELFOut % boundary % hostData = SELFStorage % boundary % hostData
  
END SUBROUTINE Equals_Scalar3D

! -- Vector2D -- !

SUBROUTINE Build_Vector2D( SELFStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(Vector2D), INTENT(out) :: SELFStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem

    CALL SELFStorage % interp % Build( N = N, &
                                      controlNodeType = quadratureType, &
                                      M = M, &
                                      targetNodeType = targetNodeType )

    CALL SELFStorage % interior % Alloc(loBound = (/1, 0, 0, 1, 1/),&
                                       upBound = (/2, N, N, nVar, nElem/))

    CALL SELFStorage % boundary % Alloc(loBound = (/1, 0, 1, 1, 1/),&
                                       upBound = (/2, N, nVar, 4, nElem/))

END SUBROUTINE Build_Vector2D

SUBROUTINE Trash_Vector2D( SELFStorage ) 
  IMPLICIT NONE
  CLASS(Vector2D), INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Trash()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

END SUBROUTINE Trash_Vector2D

#ifdef GPU
SUBROUTINE UpdateHost_Vector2D( SELFStorage )
  IMPLICIT NONE
  CLASS(Vector2D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateHost()
     CALL SELFStorage % interior % UpdateHost()
     CALL SELFStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_Vector2D

SUBROUTINE UpdateDevice_Vector2D( SELFStorage )
  IMPLICIT NONE
  CLASS(Vector2D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateDevice()
     CALL SELFStorage % interior % UpdateDevice()
     CALL SELFStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_Vector2D
#endif

FUNCTION BoundaryInterp_Vector2D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Vector2D) :: SELFStorage
  TYPE(Vector2D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % VectorBoundaryInterp_2D( SELFStorage % interior % deviceData, &
                                                          SELFout % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % VectorBoundaryInterp_2D( SELFStorage % interior % hostData, &
                                                          SELFout % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_Vector2D

FUNCTION GridInterp_Vector2D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Vector2D) :: SELFStorage
  TYPE(Vector2D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % VectorGridInterp_2D( SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % VectorGridInterp_2D( SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_Vector2D

FUNCTION Gradient_Vector2D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Vector2D) :: SELFStorage
  TYPE(Tensor2D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % VectorGradient_2D( SELFStorage % interior % deviceData, &
                                                    SELFout % interior % deviceData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % VectorGradient_2D( SELFStorage % interior % hostData, &
                                                    SELFout % interior % hostData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem )  
    ENDIF

END FUNCTION Gradient_Vector2D

FUNCTION Divergence_Vector2D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Vector2D) :: SELFStorage
  TYPE(Scalar2D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % VectorDivergence_2D( SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % VectorDivergence_2D( SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ENDIF

END FUNCTION Divergence_Vector2D

FUNCTION Curl_Vector2D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Vector2D) :: SELFStorage
  TYPE(Scalar2D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % VectorCurl_2D( SELFStorage % interior % deviceData, &
                                                SELFout % interior % deviceData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % VectorCurl_2D( SELFStorage % interior % hostData, &
                                                SELFout % interior % hostData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem )  
    ENDIF

END FUNCTION Curl_Vector2D

SUBROUTINE Equals_Vector2D( SELFOut, SELFStorage )
  IMPLICIT NONE
  TYPE(Vector2D), INTENT(out) :: SELFOut
  TYPE(Vector2D), INTENT(in) :: SELFStorage

    SELFOut % interior % hostData = SELFStorage % interior % hostData
    SELFOut % boundary % hostData = SELFStorage % boundary % hostData
  
END SUBROUTINE Equals_Vector2D

! -- Vector3D -- !

SUBROUTINE Build_Vector3D( SELFStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(Vector3D), INTENT(out) :: SELFStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem

    CALL SELFStorage % interp % Build( N = N, &
                                      controlNodeType = quadratureType, &
                                      M = M, &
                                      targetNodeType = targetNodeType )

    CALL SELFStorage % interior % Alloc(loBound = (/1, 0, 0, 0, 1, 1/),&
                                       upBound = (/3, N, N, N, nVar, nElem/))

    CALL SELFStorage % boundary % Alloc(loBound = (/1, 0, 0, 1, 1, 1/),&
                                       upBound = (/3, N, N, nVar, 6, nElem/))

END SUBROUTINE Build_Vector3D

SUBROUTINE Trash_Vector3D( SELFStorage ) 
  IMPLICIT NONE
  CLASS(Vector3D), INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Trash()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

END SUBROUTINE Trash_Vector3D

#ifdef GPU
SUBROUTINE UpdateHost_Vector3D( SELFStorage )
  IMPLICIT NONE
  CLASS(Vector3D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateHost()
     CALL SELFStorage % interior % UpdateHost()
     CALL SELFStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_Vector3D

SUBROUTINE UpdateDevice_Vector3D( SELFStorage )
  IMPLICIT NONE
  CLASS(Vector3D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateDevice()
     CALL SELFStorage % interior % UpdateDevice()
     CALL SELFStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_Vector3D
#endif

FUNCTION BoundaryInterp_Vector3D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Vector3D) :: SELFStorage
  TYPE(Vector3D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % VectorBoundaryInterp_3D( SELFStorage % interior % deviceData, &
                                                          SELFout % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % VectorBoundaryInterp_3D( SELFStorage % interior % hostData, &
                                                          SELFout % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_Vector3D

FUNCTION GridInterp_Vector3D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Vector3D) :: SELFStorage
  TYPE(Vector3D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % VectorGridInterp_3D( SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % VectorGridInterp_3D( SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_Vector3D

FUNCTION Gradient_Vector3D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Vector3D) :: SELFStorage
  TYPE(Tensor3D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % VectorGradient_3D( SELFStorage % interior % deviceData, &
                                                    SELFout % interior % deviceData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % VectorGradient_3D( SELFStorage % interior % hostData, &
                                                    SELFout % interior % hostData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem )  
    ENDIF

END FUNCTION Gradient_Vector3D

FUNCTION Divergence_Vector3D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Vector3D) :: SELFStorage
  TYPE(Scalar3D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % VectorDivergence_3D( SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % VectorDivergence_3D( SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ENDIF

END FUNCTION Divergence_Vector3D

FUNCTION Curl_Vector3D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Vector3D) :: SELFStorage
  TYPE(Vector3D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % VectorCurl_3D( SELFStorage % interior % deviceData, &
                                                SELFout % interior % deviceData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % VectorCurl_3D( SELFStorage % interior % hostData, &
                                                SELFout % interior % hostData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem )  
    ENDIF

END FUNCTION Curl_Vector3D

SUBROUTINE Equals_Vector3D( SELFOut, SELFStorage )
  IMPLICIT NONE
  TYPE(Vector3D), INTENT(out) :: SELFOut
  TYPE(Vector3D), INTENT(in) :: SELFStorage

    SELFOut % interior % hostData = SELFStorage % interior % hostData
    SELFOut % boundary % hostData = SELFStorage % boundary % hostData
  
END SUBROUTINE Equals_Vector3D

! -- Tensor2D -- !

SUBROUTINE Build_Tensor2D( SELFStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(Tensor2D), INTENT(out) :: SELFStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem

    CALL SELFStorage % interp % Build( N = N, &
                                      controlNodeType = quadratureType, &
                                      M = M, &
                                      targetNodeType = targetNodeType )

    CALL SELFStorage % interior % Alloc(loBound = (/1, 1, 0, 0, 1, 1/),&
                                       upBound = (/2, 2, N, N, nVar, nElem/))

    CALL SELFStorage % boundary % Alloc(loBound = (/1, 1, 0, 1, 1, 1/),&
                                       upBound = (/2, 2, N, nVar, 4, nElem/))

END SUBROUTINE Build_Tensor2D

SUBROUTINE Trash_Tensor2D( SELFStorage ) 
  IMPLICIT NONE
  CLASS(Tensor2D), INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Trash()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

END SUBROUTINE Trash_Tensor2D

#ifdef GPU
SUBROUTINE UpdateHost_Tensor2D( SELFStorage )
  IMPLICIT NONE
  CLASS(Tensor2D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateHost()
     CALL SELFStorage % interior % UpdateHost()
     CALL SELFStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_Tensor2D

SUBROUTINE UpdateDevice_Tensor2D( SELFStorage )
  IMPLICIT NONE
  CLASS(Tensor2D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateDevice()
     CALL SELFStorage % interior % UpdateDevice()
     CALL SELFStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_Tensor2D
#endif

FUNCTION BoundaryInterp_Tensor2D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Tensor2D) :: SELFStorage
  TYPE(Tensor2D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % TensorBoundaryInterp_2D( SELFStorage % interior % deviceData, &
                                                          SELFout % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % TensorBoundaryInterp_2D( SELFStorage % interior % hostData, &
                                                          SELFout % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_Tensor2D

FUNCTION GridInterp_Tensor2D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Tensor2D) :: SELFStorage
  TYPE(Tensor2D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % TensorGridInterp_2D( SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % TensorGridInterp_2D( SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_Tensor2D

FUNCTION Divergence_Tensor2D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Tensor2D) :: SELFStorage
  TYPE(Vector2D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % TensorDivergence_2D( SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % TensorDivergence_2D( SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ENDIF

END FUNCTION Divergence_Tensor2D

FUNCTION Determinant_Tensor2D( SELFStorage ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Tensor2D) :: SELFStorage
  TYPE(Scalar2D) :: SELFOut
  ! Local
  INTEGER :: iEl, iVar, i, j

    DO iEl = 1, SELFStorage % nElem
      DO iVar = 1, SELFStorage % nVar
          DO j = 0, SELFStorage % N
            DO i = 0, SELFStorage % N

              SELFOut % interior % hostData(i,j,iVar,iEl) = SELFStorage % interior % hostData(1,1,i,j,iVar,iEl)*&
                                                           SELFStorage % interior % hostData(2,2,i,j,iVar,iEl)-&
                                                           SELFStorage % interior % hostData(1,2,i,j,iVar,iEl)*&
                                                           SELFStorage % interior % hostData(2,1,i,j,iVar,iEl)
                                                     

            ENDDO
          ENDDO
        ENDDO
    ENDDO
    
END FUNCTION Determinant_Tensor2D

SUBROUTINE Equals_Tensor2D( SELFOut, SELFStorage )
  IMPLICIT NONE
  TYPE(Tensor2D), INTENT(out) :: SELFOut
  TYPE(Tensor2D), INTENT(in) :: SELFStorage

    SELFOut % interior % hostData = SELFStorage % interior % hostData
    SELFOut % boundary % hostData = SELFStorage % boundary % hostData
  
END SUBROUTINE Equals_Tensor2D

! -- Tensor3D -- !

SUBROUTINE Build_Tensor3D( SELFStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(Tensor3D), INTENT(out) :: SELFStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SELFStorage % N = N
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem

    CALL SELFStorage % interp % Build( N = N, &
                                      controlNodeType = quadratureType, &
                                      M = M, &
                                      targetNodeType = targetNodeType )

    CALL SELFStorage % interior % Alloc(loBound = (/1, 1, 0, 0, 0, 1, 1/),&
                                       upBound = (/3, 3, N, N, N, nVar, nElem/))

    CALL SELFStorage % boundary % Alloc(loBound = (/1, 1, 0, 0, 1, 1, 1/),&
                                       upBound = (/3, 3, N, N, nVar, 6, nElem/))

END SUBROUTINE Build_Tensor3D

SUBROUTINE Trash_Tensor3D( SELFStorage ) 
  IMPLICIT NONE
  CLASS(Tensor3D), INTENT(inout) :: SELFStorage

    SELFStorage % N = 0
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interp % Trash()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()

END SUBROUTINE Trash_Tensor3D

#ifdef GPU
SUBROUTINE UpdateHost_Tensor3D( SELFStorage )
  IMPLICIT NONE
  CLASS(Tensor3D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateHost()
     CALL SELFStorage % interior % UpdateHost()
     CALL SELFStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_Tensor3D

SUBROUTINE UpdateDevice_Tensor3D( SELFStorage )
  IMPLICIT NONE
  CLASS(Tensor3D), INTENT(inout) :: SELFStorage

     CALL SELFStorage % interp % UpdateDevice()
     CALL SELFStorage % interior % UpdateDevice()
     CALL SELFStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_Tensor3D
#endif

FUNCTION BoundaryInterp_Tensor3D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Tensor3D) :: SELFStorage
  TYPE(Tensor3D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % TensorBoundaryInterp_3D( SELFStorage % interior % deviceData, &
                                                          SELFout % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % TensorBoundaryInterp_3D( SELFStorage % interior % hostData, &
                                                          SELFout % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_Tensor3D

FUNCTION GridInterp_Tensor3D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Tensor3D) :: SELFStorage
  TYPE(Tensor3D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % TensorGridInterp_3D( SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % TensorGridInterp_3D( SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_Tensor3D

FUNCTION Divergence_Tensor3D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Tensor3D) :: SELFStorage
  TYPE(Vector3D) :: SELFOut
  LOGICAL :: gpuAccel

    IF( gpuAccel )THEN
      CALL SELFStorage % interp % TensorDivergence_3D( SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ELSE
      CALL SELFStorage % interp % TensorDivergence_3D( SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem )  
    ENDIF

END FUNCTION Divergence_Tensor3D

FUNCTION Determinant_Tensor3D( SELFStorage ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(Tensor3D) :: SELFStorage
  TYPE(Scalar3D) :: SELFOut
  ! Local
  INTEGER :: iEl, iVar, i, j, k

    DO iEl = 1, SELFStorage % nElem
      DO iVar = 1, SELFStorage % nVar
        DO k = 0, SELFStorage % N
          DO j = 0, SELFStorage % N
            DO i = 0, SELFStorage % N

              SELFOut % interior % hostData(i,j,k,iVar,iEl) = SELFStorage % interior % hostData(1,1,i,j,k,iVar,iEl)*&
                                                  ( SELFStorage % interior % hostData(2,2,i,j,k,iVar,iEl)*&
                                                    SELFStorage % interior % hostData(3,3,i,j,k,iVar,iEl)-&
                                                    SELFStorage % interior % hostData(2,3,i,j,k,iVar,iEl)*&
                                                    SELFStorage % interior % hostData(3,2,i,j,k,iVar,iEl) )-&
                                                 SELFStorage % interior % hostData(2,1,i,j,k,iVar,iEl)*&
                                                  ( SELFStorage % interior % hostData(1,2,i,j,k,iVar,iEl)*&
                                                    SELFStorage % interior % hostData(3,3,i,j,k,iVar,iEl)-&
                                                    SELFStorage % interior % hostData(1,3,i,j,k,iVar,iEl)*&
                                                    SELFStorage % interior % hostData(3,2,i,j,k,iVar,iEl) )+&
                                                 SELFStorage % interior % hostData(3,1,i,j,k,iVar,iEl)*&
                                                  ( SELFStorage % interior % hostData(1,2,i,j,k,iVar,iEl)*&
                                                    SELFStorage % interior % hostData(2,3,i,j,k,iVar,iEl)-&
                                                    SELFStorage % interior % hostData(1,3,i,j,k,iVar,iEl)*&
                                                    SELFStorage % interior % hostData(2,2,i,j,k,iVar,iEl) )
                                                     

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    
END FUNCTION Determinant_Tensor3D

SUBROUTINE Equals_Tensor3D( SELFOut, SELFStorage )
  IMPLICIT NONE
  TYPE(Tensor3D), INTENT(out) :: SELFOut
  TYPE(Tensor3D), INTENT(in) :: SELFStorage

    SELFOut % interior % hostData = SELFStorage % interior % hostData
    SELFOut % boundary % hostData = SELFStorage % boundary % hostData
  
END SUBROUTINE Equals_Tensor3D

END MODULE SELF_Data
