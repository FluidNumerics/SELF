! NodalSEMData.F90
! 
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self-fluids@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE NodalSEMData

USE SELFConstants
USE Lagrange_Class

USE hipfort
USE ISO_C_BINDING


IMPLICIT NONE

! ---------------------- Scalars ---------------------- !
  TYPE, PUBLIC :: SEMScalar1D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r3) :: interior
    TYPE(hfReal_r3) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_SEMScalar1D
      PROCEDURE, PUBLIC :: Trash => Trash_SEMScalar1D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMScalar1D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMScalar1D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_SEMScalar1D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMScalar1D
      PROCEDURE, PUBLIC :: Derivative => Derivative_SEMScalar1D

  END TYPE SEMScalar1D

  TYPE, PUBLIC :: SEMScalar2D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r4) :: interior
    TYPE(hfReal_r4) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_SEMScalar2D
      PROCEDURE, PUBLIC :: Trash => Trash_SEMScalar2D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMScalar2D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMScalar2D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_SEMScalar2D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMScalar2D
      PROCEDURE, PUBLIC :: Gradient => Gradient_SEMScalar2D

  END TYPE SEMScalar2D

  TYPE, PUBLIC :: SEMScalar3D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r5) :: interior
    TYPE(hfReal_r5) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_SEMScalar3D
      PROCEDURE, PUBLIC :: Trash => Trash_SEMScalar3D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMScalar3D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMScalar3D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_SEMScalar3D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMScalar3D
      PROCEDURE, PUBLIC :: Gradient => Gradient_SEMScalar3D

  END TYPE SEMScalar3D

! ---------------------- Vectors ---------------------- !

  TYPE, PUBLIC :: SEMVector2D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r5) :: interior
    TYPE(hfReal_r5) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_SEMVector2D
      PROCEDURE, PUBLIC :: Trash => Trash_SEMVector2D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMVector2D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMVector2D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_SEMVector2D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMVector2D
      PROCEDURE, PUBLIC :: Gradient => Gradient_SEMVector2D
      PROCEDURE, PUBLIC :: Divergence => Divergence_SEMVector2D
      PROCEDURE, PUBLIC :: Curl => Curl_SEMVector2D


  END TYPE SEMVector2D

  TYPE, PUBLIC :: SEMVector3D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r6) :: interior
    TYPE(hfReal_r6) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_SEMVector3D
      PROCEDURE, PUBLIC :: Trash => Trash_SEMVector3D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMVector3D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMVector3D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_SEMVector3D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMVector3D
      PROCEDURE, PUBLIC :: Gradient => Gradient_SEMVector3D
      PROCEDURE, PUBLIC :: Divergence => Divergence_SEMVector3D
      PROCEDURE, PUBLIC :: Curl => Curl_SEMVector3D


  END TYPE SEMVector3D
!
!! ---------------------- Tensors ---------------------- !
!
  TYPE, PUBLIC :: SEMTensor2D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r6) :: interior
    TYPE(hfReal_r6) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_SEMTensor2D
      PROCEDURE, PUBLIC :: Trash => Trash_SEMTensor2D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMTensor2D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMTensor2D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_SEMTensor2D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMTensor2D 

      PROCEDURE, PUBLIC :: Determinant => Determinant_SEMTensor2D

  END TYPE SEMTensor2D

  TYPE, PUBLIC :: SEMTensor3D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange) :: interp
    TYPE(hfReal_r7) :: interior
    TYPE(hfReal_r7) :: boundary

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_SEMTensor3D
      PROCEDURE, PUBLIC :: Trash => Trash_SEMTensor3D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMTensor3D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMTensor3D
#endif
      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_SEMTensor3D
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMTensor3D 

      PROCEDURE, PUBLIC :: Determinant => Determinant_SEMTensor3D

  END TYPE SEMTensor3D

  INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE Equals_SEMScalar1D
    MODULE PROCEDURE Equals_SEMScalar2D
    MODULE PROCEDURE Equals_SEMScalar3D
    MODULE PROCEDURE Equals_SEMVector2D
    MODULE PROCEDURE Equals_SEMVector3D
    MODULE PROCEDURE Equals_SEMTensor2D
    MODULE PROCEDURE Equals_SEMTensor3D
  END INTERFACE

!  INTERFACE OPERATOR (+)
!    MODULE PROCEDURE Add_SEMScalar1D
!    MODULE PROCEDURE Add_SEMScalar2D
!    MODULE PROCEDURE Add_SEMScalar3D
!    MODULE PROCEDURE Add_SEMVector2D
!    MODULE PROCEDURE Add_SEMVector3D
!    MODULE PROCEDURE Add_SEMTensor2D
!  END INTERFACE
!
!  INTERFACE OPERATOR (-)
!    MODULE PROCEDURE Subtract_SEMScalar1D
!    MODULE PROCEDURE Subtract_SEMScalar2D
!    MODULE PROCEDURE Subtract_SEMScalar3D
!    MODULE PROCEDURE Subtract_SEMVector2D
!    MODULE PROCEDURE Subtract_SEMVector3D
!    MODULE PROCEDURE Subtract_SEMTensor2D
!  END INTERFACE

CONTAINS

! -- SEMScalar1D -- !

SUBROUTINE Build_SEMScalar1D( SEMStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(SEMScalar1D), INTENT(out) :: SEMStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SEMStorage % N = N
    SEMStorage % nVar = nVar
    SEMStorage % nElem = nElem

    CALL SEMStorage % interp % Build(N = N, &
                                     controlNodeType = quadratureType, &
                                     M = M, &
                                     targetNodeType = targetNodeType)

    CALL SEMStorage % interior % Alloc(loBound = (/0, 1, 1/),&
                                       upBound = (/N, nVar, nElem/))

    CALL SEMStorage % boundary % Alloc(loBound = (/1, 1, 1/),&
                                       upBound = (/nVar, 2, nElem/))


END SUBROUTINE Build_SEMScalar1D

SUBROUTINE Trash_SEMScalar1D( SEMStorage ) 
  IMPLICIT NONE
  CLASS(SEMScalar1D), INTENT(inout) :: SEMStorage

    SEMStorage % N = 0
    SEMStorage % nVar = 0
    SEMStorage % nElem = 0
    CALL SEMStorage % interp % Trash()
    CALL SEMStorage % interior % Free()
    CALL SEMStorage % boundary % Free()

END SUBROUTINE Trash_SEMScalar1D

#ifdef GPU
SUBROUTINE UpdateHost_SEMScalar1D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMScalar1D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateHost()
     CALL SEMStorage % interior % UpdateHost()
     CALL SEMStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_SEMScalar1D

SUBROUTINE UpdateDevice_SEMScalar1D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMScalar1D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateDevice()
     CALL SEMStorage % interior % UpdateDevice()
     CALL SEMStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_SEMScalar1D
#endif

FUNCTION BoundaryInterp_SEMScalar1D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar1D) :: SEMStorage
  TYPE(SEMScalar1D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % ScalarBoundaryInterp_1D( SEMStorage % interior % deviceData, &
                                                          SEMout % boundary % deviceData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % ScalarBoundaryInterp_1D( SEMStorage % interior % hostData, &
                                                          SEMout % boundary % hostData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_SEMScalar1D

FUNCTION GridInterp_SEMScalar1D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar1D) :: SEMStorage
  TYPE(SEMScalar1D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % ScalarGridInterp_1D( SEMStorage % interior % deviceData, &
                                                      SEMout % interior % deviceData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % ScalarGridInterp_1D( SEMStorage % interior % hostData, &
                                                      SEMout % interior % hostData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_SEMScalar1D

FUNCTION Derivative_SEMScalar1D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar1D) :: SEMStorage
  TYPE(SEMScalar1D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % Derivative_1D( SEMStorage % interior % deviceData, &
                                                SEMout % interior % deviceData, &
                                                SEMStorage % nVar, &
                                                SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % Derivative_1D( SEMStorage % interior % hostData, &
                                                SEMout % interior % hostData, &
                                                SEMStorage % nVar, &
                                                SEMStorage % nElem )  
    ENDIF

END FUNCTION Derivative_SEMScalar1D

SUBROUTINE Equals_SEMScalar1D( SEMOut, SEMStorage )
  IMPLICIT NONE
  TYPE(SEMScalar1D), INTENT(out) :: SEMOut
  TYPE(SEMScalar1D), INTENT(in) :: SEMStorage

    SEMOut % interior % hostData = SEMStorage % interior % hostData
    SEMOut % boundary % hostData = SEMStorage % boundary % hostData
  
END SUBROUTINE Equals_SEMScalar1D

! -- SEMScalar2D -- !

SUBROUTINE Build_SEMScalar2D( SEMStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(SEMScalar2D), INTENT(out) :: SEMStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SEMStorage % N = N
    SEMStorage % nVar = nVar
    SEMStorage % nElem = nElem

    CALL SEMStorage % interp % Build( N = N, &
                                      controlNodeType = quadratureType, &
                                      M = M, &
                                      targetNodeType = targetNodeType )

    CALL SEMStorage % interior % Alloc(loBound = (/0, 0, 1, 1/),&
                                       upBound = (/N, N, nVar, nElem/))

    CALL SEMStorage % boundary % Alloc(loBound = (/0, 1, 1, 1/),&
                                       upBound = (/N, nVar, 4, nElem/))


END SUBROUTINE Build_SEMScalar2D

SUBROUTINE Trash_SEMScalar2D( SEMStorage ) 
  IMPLICIT NONE
  CLASS(SEMScalar2D), INTENT(inout) :: SEMStorage

    SEMStorage % N = 0
    SEMStorage % nVar = 0
    SEMStorage % nElem = 0
    CALL SEMStorage % interp % Trash()
    CALL SEMStorage % interior % Free()
    CALL SEMStorage % boundary % Free()

END SUBROUTINE Trash_SEMScalar2D

#ifdef GPU
SUBROUTINE UpdateHost_SEMScalar2D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMScalar2D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateHost()
     CALL SEMStorage % interior % UpdateHost()
     CALL SEMStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_SEMScalar2D

SUBROUTINE UpdateDevice_SEMScalar2D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMScalar2D), INTENT(inout) :: SEMStorage
          
     CALL SEMStorage % interp % UpdateDevice()
     CALL SEMStorage % interior % UpdateDevice()
     CALL SEMStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_SEMScalar2D
#endif

FUNCTION BoundaryInterp_SEMScalar2D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar2D) :: SEMStorage
  TYPE(SEMScalar2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % ScalarBoundaryInterp_2D( SEMStorage % interior % deviceData, &
                                                          SEMout % boundary % deviceData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % ScalarBoundaryInterp_2D( SEMStorage % interior % hostData, &
                                                          SEMout % boundary % hostData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_SEMScalar2D

FUNCTION GridInterp_SEMScalar2D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar2D) :: SEMStorage
  TYPE(SEMScalar2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % ScalarGridInterp_2D( SEMStorage % interior % deviceData, &
                                                      SEMout % interior % deviceData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % ScalarGridInterp_2D( SEMStorage % interior % hostData, &
                                                      SEMout % interior % hostData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_SEMScalar2D

FUNCTION Gradient_SEMScalar2D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar2D) :: SEMStorage
  TYPE(SEMVector2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % ScalarGradient_2D( SEMStorage % interior % deviceData, &
                                                    SEMout % interior % deviceData, &
                                                    SEMStorage % nVar, &
                                                    SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % ScalarGradient_2D( SEMStorage % interior % hostData, &
                                                    SEMout % interior % hostData, &
                                                    SEMStorage % nVar, &
                                                    SEMStorage % nElem )  
    ENDIF

END FUNCTION Gradient_SEMScalar2D

SUBROUTINE Equals_SEMScalar2D( SEMOut, SEMStorage )
  IMPLICIT NONE
  TYPE(SEMScalar2D), INTENT(out) :: SEMOut
  TYPE(SEMScalar2D), INTENT(in) :: SEMStorage

    SEMOut % interior % hostData = SEMStorage % interior % hostData
    SEMOut % boundary % hostData = SEMStorage % boundary % hostData
  
END SUBROUTINE Equals_SEMScalar2D

! -- SEMScalar3D -- !

SUBROUTINE Build_SEMScalar3D( SEMStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(SEMScalar3D), INTENT(out) :: SEMStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SEMStorage % N = N
    SEMStorage % nVar = nVar
    SEMStorage % nElem = nElem

    CALL SEMStorage % interp % Build( N = N, &
                                      controlNodeType = quadratureType, &
                                      M = M, &
                                      targetNodeType = targetNodeType )

    CALL SEMStorage % interior % Alloc(loBound = (/0, 0, 0, 1, 1/),&
                                       upBound = (/N, N, N, nVar, nElem/))

    CALL SEMStorage % boundary % Alloc(loBound = (/0, 0, 1, 1, 1/),&
                                       upBound = (/N, N, nVar, 6, nElem/))

END SUBROUTINE Build_SEMScalar3D

SUBROUTINE Trash_SEMScalar3D( SEMStorage ) 
  IMPLICIT NONE
  CLASS(SEMScalar3D), INTENT(inout) :: SEMStorage

    SEMStorage % N = 0
    SEMStorage % nVar = 0
    SEMStorage % nElem = 0
    CALL SEMStorage % interp % Trash()
    CALL SEMStorage % interior % Free()
    CALL SEMStorage % boundary % Free()

END SUBROUTINE Trash_SEMScalar3D

#ifdef GPU
SUBROUTINE UpdateHost_SEMScalar3D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMScalar3D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateHost()
     CALL SEMStorage % interior % UpdateHost()
     CALL SEMStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_SEMScalar3D

SUBROUTINE UpdateDevice_SEMScalar3D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMScalar3D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateDevice()
     CALL SEMStorage % interior % UpdateDevice()
     CALL SEMStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_SEMScalar3D
#endif

FUNCTION BoundaryInterp_SEMScalar3D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar3D) :: SEMStorage
  TYPE(SEMScalar3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % ScalarBoundaryInterp_3D( SEMStorage % interior % deviceData, &
                                                          SEMout % boundary % deviceData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % ScalarBoundaryInterp_3D( SEMStorage % interior % hostData, &
                                                          SEMout % boundary % hostData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_SEMScalar3D

FUNCTION GridInterp_SEMScalar3D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar3D) :: SEMStorage
  TYPE(SEMScalar3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % ScalarGridInterp_3D( SEMStorage % interior % deviceData, &
                                                      SEMout % interior % deviceData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % ScalarGridInterp_3D( SEMStorage % interior % hostData, &
                                                      SEMout % interior % hostData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_SEMScalar3D

FUNCTION Gradient_SEMScalar3D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar3D) :: SEMStorage
  TYPE(SEMVector3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % ScalarGradient_3D( SEMStorage % interior % deviceData, &
                                                    SEMout % interior % deviceData, &
                                                    SEMStorage % nVar, &
                                                    SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % ScalarGradient_3D( SEMStorage % interior % hostData, &
                                                    SEMout % interior % hostData, &
                                                    SEMStorage % nVar, &
                                                    SEMStorage % nElem )  
    ENDIF

END FUNCTION Gradient_SEMScalar3D

SUBROUTINE Equals_SEMScalar3D( SEMOut, SEMStorage )
  IMPLICIT NONE
  TYPE(SEMScalar3D), INTENT(out) :: SEMOut
  TYPE(SEMScalar3D), INTENT(in) :: SEMStorage

    SEMOut % interior % hostData = SEMStorage % interior % hostData
    SEMOut % boundary % hostData = SEMStorage % boundary % hostData
  
END SUBROUTINE Equals_SEMScalar3D

! -- SEMVector2D -- !

SUBROUTINE Build_SEMVector2D( SEMStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(SEMVector2D), INTENT(out) :: SEMStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SEMStorage % N = N
    SEMStorage % nVar = nVar
    SEMStorage % nElem = nElem

    CALL SEMStorage % interp % Build( N = N, &
                                      controlNodeType = quadratureType, &
                                      M = M, &
                                      targetNodeType = targetNodeType )

    CALL SEMStorage % interior % Alloc(loBound = (/1, 0, 0, 1, 1/),&
                                       upBound = (/2, N, N, nVar, nElem/))

    CALL SEMStorage % boundary % Alloc(loBound = (/1, 0, 1, 1, 1/),&
                                       upBound = (/2, N, nVar, 4, nElem/))

END SUBROUTINE Build_SEMVector2D

SUBROUTINE Trash_SEMVector2D( SEMStorage ) 
  IMPLICIT NONE
  CLASS(SEMVector2D), INTENT(inout) :: SEMStorage

    SEMStorage % N = 0
    SEMStorage % nVar = 0
    SEMStorage % nElem = 0
    CALL SEMStorage % interp % Trash()
    CALL SEMStorage % interior % Free()
    CALL SEMStorage % boundary % Free()

END SUBROUTINE Trash_SEMVector2D

#ifdef GPU
SUBROUTINE UpdateHost_SEMVector2D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMVector2D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateHost()
     CALL SEMStorage % interior % UpdateHost()
     CALL SEMStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_SEMVector2D

SUBROUTINE UpdateDevice_SEMVector2D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMVector2D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateDevice()
     CALL SEMStorage % interior % UpdateDevice()
     CALL SEMStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_SEMVector2D
#endif

FUNCTION BoundaryInterp_SEMVector2D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector2D) :: SEMStorage
  TYPE(SEMVector2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % VectorBoundaryInterp_2D( SEMStorage % interior % deviceData, &
                                                          SEMout % boundary % deviceData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % VectorBoundaryInterp_2D( SEMStorage % interior % hostData, &
                                                          SEMout % boundary % hostData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_SEMVector2D

FUNCTION GridInterp_SEMVector2D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector2D) :: SEMStorage
  TYPE(SEMVector2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % VectorGridInterp_2D( SEMStorage % interior % deviceData, &
                                                      SEMout % interior % deviceData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % VectorGridInterp_2D( SEMStorage % interior % hostData, &
                                                      SEMout % interior % hostData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_SEMVector2D

FUNCTION Gradient_SEMVector2D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector2D) :: SEMStorage
  TYPE(SEMTensor2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % VectorGradient_2D( SEMStorage % interior % deviceData, &
                                                    SEMout % interior % deviceData, &
                                                    SEMStorage % nVar, &
                                                    SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % VectorGradient_2D( SEMStorage % interior % hostData, &
                                                    SEMout % interior % hostData, &
                                                    SEMStorage % nVar, &
                                                    SEMStorage % nElem )  
    ENDIF

END FUNCTION Gradient_SEMVector2D

FUNCTION Divergence_SEMVector2D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector2D) :: SEMStorage
  TYPE(SEMScalar2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % VectorDivergence_2D( SEMStorage % interior % deviceData, &
                                                      SEMout % interior % deviceData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % VectorDivergence_2D( SEMStorage % interior % hostData, &
                                                      SEMout % interior % hostData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ENDIF

END FUNCTION Divergence_SEMVector2D

FUNCTION Curl_SEMVector2D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector2D) :: SEMStorage
  TYPE(SEMScalar2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % VectorCurl_2D( SEMStorage % interior % deviceData, &
                                                SEMout % interior % deviceData, &
                                                SEMStorage % nVar, &
                                                SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % VectorCurl_2D( SEMStorage % interior % hostData, &
                                                SEMout % interior % hostData, &
                                                SEMStorage % nVar, &
                                                SEMStorage % nElem )  
    ENDIF

END FUNCTION Curl_SEMVector2D

SUBROUTINE Equals_SEMVector2D( SEMOut, SEMStorage )
  IMPLICIT NONE
  TYPE(SEMVector2D), INTENT(out) :: SEMOut
  TYPE(SEMVector2D), INTENT(in) :: SEMStorage

    SEMOut % interior % hostData = SEMStorage % interior % hostData
    SEMOut % boundary % hostData = SEMStorage % boundary % hostData
  
END SUBROUTINE Equals_SEMVector2D

! -- SEMVector3D -- !

SUBROUTINE Build_SEMVector3D( SEMStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(SEMVector3D), INTENT(out) :: SEMStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SEMStorage % N = N
    SEMStorage % nVar = nVar
    SEMStorage % nElem = nElem

    CALL SEMStorage % interp % Build( N = N, &
                                      controlNodeType = quadratureType, &
                                      M = M, &
                                      targetNodeType = targetNodeType )

    CALL SEMStorage % interior % Alloc(loBound = (/1, 0, 0, 0, 1, 1/),&
                                       upBound = (/3, N, N, N, nVar, nElem/))

    CALL SEMStorage % boundary % Alloc(loBound = (/1, 0, 0, 1, 1, 1/),&
                                       upBound = (/3, N, N, nVar, 6, nElem/))

END SUBROUTINE Build_SEMVector3D

SUBROUTINE Trash_SEMVector3D( SEMStorage ) 
  IMPLICIT NONE
  CLASS(SEMVector3D), INTENT(inout) :: SEMStorage

    SEMStorage % N = 0
    SEMStorage % nVar = 0
    SEMStorage % nElem = 0
    CALL SEMStorage % interp % Trash()
    CALL SEMStorage % interior % Free()
    CALL SEMStorage % boundary % Free()

END SUBROUTINE Trash_SEMVector3D

#ifdef GPU
SUBROUTINE UpdateHost_SEMVector3D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMVector3D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateHost()
     CALL SEMStorage % interior % UpdateHost()
     CALL SEMStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_SEMVector3D

SUBROUTINE UpdateDevice_SEMVector3D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMVector3D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateDevice()
     CALL SEMStorage % interior % UpdateDevice()
     CALL SEMStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_SEMVector3D
#endif

FUNCTION BoundaryInterp_SEMVector3D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector3D) :: SEMStorage
  TYPE(SEMVector3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % VectorBoundaryInterp_3D( SEMStorage % interior % deviceData, &
                                                          SEMout % boundary % deviceData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % VectorBoundaryInterp_3D( SEMStorage % interior % hostData, &
                                                          SEMout % boundary % hostData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_SEMVector3D

FUNCTION GridInterp_SEMVector3D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector3D) :: SEMStorage
  TYPE(SEMVector3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % VectorGridInterp_3D( SEMStorage % interior % deviceData, &
                                                      SEMout % interior % deviceData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % VectorGridInterp_3D( SEMStorage % interior % hostData, &
                                                      SEMout % interior % hostData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_SEMVector3D

FUNCTION Gradient_SEMVector3D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector3D) :: SEMStorage
  TYPE(SEMTensor3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % VectorGradient_3D( SEMStorage % interior % deviceData, &
                                                    SEMout % interior % deviceData, &
                                                    SEMStorage % nVar, &
                                                    SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % VectorGradient_3D( SEMStorage % interior % hostData, &
                                                    SEMout % interior % hostData, &
                                                    SEMStorage % nVar, &
                                                    SEMStorage % nElem )  
    ENDIF

END FUNCTION Gradient_SEMVector3D

FUNCTION Divergence_SEMVector3D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector3D) :: SEMStorage
  TYPE(SEMScalar3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % VectorDivergence_3D( SEMStorage % interior % deviceData, &
                                                      SEMout % interior % deviceData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % VectorDivergence_3D( SEMStorage % interior % hostData, &
                                                      SEMout % interior % hostData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ENDIF

END FUNCTION Divergence_SEMVector3D

FUNCTION Curl_SEMVector3D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector3D) :: SEMStorage
  TYPE(SEMVector3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % VectorCurl_3D( SEMStorage % interior % deviceData, &
                                                SEMout % interior % deviceData, &
                                                SEMStorage % nVar, &
                                                SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % VectorCurl_3D( SEMStorage % interior % hostData, &
                                                SEMout % interior % hostData, &
                                                SEMStorage % nVar, &
                                                SEMStorage % nElem )  
    ENDIF

END FUNCTION Curl_SEMVector3D

SUBROUTINE Equals_SEMVector3D( SEMOut, SEMStorage )
  IMPLICIT NONE
  TYPE(SEMVector3D), INTENT(out) :: SEMOut
  TYPE(SEMVector3D), INTENT(in) :: SEMStorage

    SEMOut % interior % hostData = SEMStorage % interior % hostData
    SEMOut % boundary % hostData = SEMStorage % boundary % hostData
  
END SUBROUTINE Equals_SEMVector3D

! -- SEMTensor2D -- !

SUBROUTINE Build_SEMTensor2D( SEMStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(SEMTensor2D), INTENT(out) :: SEMStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SEMStorage % N = N
    SEMStorage % nVar = nVar
    SEMStorage % nElem = nElem

    CALL SEMStorage % interp % Build( N = N, &
                                      controlNodeType = quadratureType, &
                                      M = M, &
                                      targetNodeType = targetNodeType )

    CALL SEMStorage % interior % Alloc(loBound = (/1, 1, 0, 0, 1, 1/),&
                                       upBound = (/2, 2, N, N, nVar, nElem/))

    CALL SEMStorage % boundary % Alloc(loBound = (/1, 1, 0, 1, 1, 1/),&
                                       upBound = (/2, 2, N, nVar, 4, nElem/))

END SUBROUTINE Build_SEMTensor2D

SUBROUTINE Trash_SEMTensor2D( SEMStorage ) 
  IMPLICIT NONE
  CLASS(SEMTensor2D), INTENT(inout) :: SEMStorage

    SEMStorage % N = 0
    SEMStorage % nVar = 0
    SEMStorage % nElem = 0
    CALL SEMStorage % interp % Trash()
    CALL SEMStorage % interior % Free()
    CALL SEMStorage % boundary % Free()

END SUBROUTINE Trash_SEMTensor2D

#ifdef GPU
SUBROUTINE UpdateHost_SEMTensor2D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMTensor2D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateHost()
     CALL SEMStorage % interior % UpdateHost()
     CALL SEMStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_SEMTensor2D

SUBROUTINE UpdateDevice_SEMTensor2D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMTensor2D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateDevice()
     CALL SEMStorage % interior % UpdateDevice()
     CALL SEMStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_SEMTensor2D
#endif

FUNCTION BoundaryInterp_SEMTensor2D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMTensor2D) :: SEMStorage
  TYPE(SEMTensor2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % TensorBoundaryInterp_2D( SEMStorage % interior % deviceData, &
                                                          SEMout % boundary % deviceData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % TensorBoundaryInterp_2D( SEMStorage % interior % hostData, &
                                                          SEMout % boundary % hostData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_SEMTensor2D

FUNCTION GridInterp_SEMTensor2D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMTensor2D) :: SEMStorage
  TYPE(SEMTensor2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % TensorGridInterp_2D( SEMStorage % interior % deviceData, &
                                                      SEMout % interior % deviceData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % TensorGridInterp_2D( SEMStorage % interior % hostData, &
                                                      SEMout % interior % hostData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_SEMTensor2D

FUNCTION Determinant_SEMTensor2D( SEMStorage ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMTensor2D) :: SEMStorage
  TYPE(SEMScalar2D) :: SEMOut
  ! Local
  INTEGER :: iEl, iVar, i, j

    DO iEl = 1, SEMStorage % nElem
      DO iVar = 1, SEMStorage % nVar
          DO j = 0, SEMStorage % N
            DO i = 0, SEMStorage % N

              SEMOut % interior % hostData(i,j,iVar,iEl) = SEMStorage % interior % hostData(1,1,i,j,iVar,iEl)*&
                                                           SEMStorage % interior % hostData(2,2,i,j,iVar,iEl)-&
                                                           SEMStorage % interior % hostData(1,2,i,j,iVar,iEl)*&
                                                           SEMStorage % interior % hostData(2,1,i,j,iVar,iEl)
                                                     

            ENDDO
          ENDDO
        ENDDO
    ENDDO
    
END FUNCTION Determinant_SEMTensor2D

SUBROUTINE Equals_SEMTensor2D( SEMOut, SEMStorage )
  IMPLICIT NONE
  TYPE(SEMTensor2D), INTENT(out) :: SEMOut
  TYPE(SEMTensor2D), INTENT(in) :: SEMStorage

    SEMOut % interior % hostData = SEMStorage % interior % hostData
    SEMOut % boundary % hostData = SEMStorage % boundary % hostData
  
END SUBROUTINE Equals_SEMTensor2D

! -- SEMTensor3D -- !

SUBROUTINE Build_SEMTensor3D( SEMStorage, N, quadratureType, M, targetNodeType, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(SEMTensor3D), INTENT(out) :: SEMStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: quadratureType
  INTEGER, INTENT(in) :: M
  INTEGER, INTENT(in) :: targetNodeType
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SEMStorage % N = N
    SEMStorage % nVar = nVar
    SEMStorage % nElem = nElem

    CALL SEMStorage % interp % Build( N = N, &
                                      controlNodeType = quadratureType, &
                                      M = M, &
                                      targetNodeType = targetNodeType )

    CALL SEMStorage % interior % Alloc(loBound = (/1, 1, 0, 0, 0, 1, 1/),&
                                       upBound = (/3, 3, N, N, N, nVar, nElem/))

    CALL SEMStorage % boundary % Alloc(loBound = (/1, 1, 0, 0, 1, 1, 1/),&
                                       upBound = (/3, 3, N, N, nVar, 6, nElem/))

END SUBROUTINE Build_SEMTensor3D

SUBROUTINE Trash_SEMTensor3D( SEMStorage ) 
  IMPLICIT NONE
  CLASS(SEMTensor3D), INTENT(inout) :: SEMStorage

    SEMStorage % N = 0
    SEMStorage % nVar = 0
    SEMStorage % nElem = 0
    CALL SEMStorage % interp % Trash()
    CALL SEMStorage % interior % Free()
    CALL SEMStorage % boundary % Free()

END SUBROUTINE Trash_SEMTensor3D

#ifdef GPU
SUBROUTINE UpdateHost_SEMTensor3D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMTensor3D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateHost()
     CALL SEMStorage % interior % UpdateHost()
     CALL SEMStorage % boundary % UpdateHost()

END SUBROUTINE UpdateHost_SEMTensor3D

SUBROUTINE UpdateDevice_SEMTensor3D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMTensor3D), INTENT(inout) :: SEMStorage

     CALL SEMStorage % interp % UpdateDevice()
     CALL SEMStorage % interior % UpdateDevice()
     CALL SEMStorage % boundary % UpdateDevice()

END SUBROUTINE UpdateDevice_SEMTensor3D
#endif

FUNCTION BoundaryInterp_SEMTensor3D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMTensor3D) :: SEMStorage
  TYPE(SEMTensor3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % TensorBoundaryInterp_3D( SEMStorage % interior % deviceData, &
                                                          SEMout % boundary % deviceData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % TensorBoundaryInterp_3D( SEMStorage % interior % hostData, &
                                                          SEMout % boundary % hostData, &
                                                          SEMStorage % nVar, &
                                                          SEMStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_SEMTensor3D

FUNCTION GridInterp_SEMTensor3D( SEMStorage, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMTensor3D) :: SEMStorage
  TYPE(SEMTensor3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL SEMStorage % interp % TensorGridInterp_3D( SEMStorage % interior % deviceData, &
                                                      SEMout % interior % deviceData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ELSE
      CALL SEMStorage % interp % TensorGridInterp_3D( SEMStorage % interior % hostData, &
                                                      SEMout % interior % hostData, &
                                                      SEMStorage % nVar, &
                                                      SEMStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_SEMTensor3D

FUNCTION Determinant_SEMTensor3D( SEMStorage ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMTensor3D) :: SEMStorage
  TYPE(SEMScalar3D) :: SEMOut
  ! Local
  INTEGER :: iEl, iVar, i, j, k

    DO iEl = 1, SEMStorage % nElem
      DO iVar = 1, SEMStorage % nVar
        DO k = 0, SEMStorage % N
          DO j = 0, SEMStorage % N
            DO i = 0, SEMStorage % N

              SEMOut % interior % hostData(i,j,k,iVar,iEl) = SEMStorage % interior % hostData(1,1,i,j,k,iVar,iEl)*&
                                                  ( SEMStorage % interior % hostData(2,2,i,j,k,iVar,iEl)*&
                                                    SEMStorage % interior % hostData(3,3,i,j,k,iVar,iEl)-&
                                                    SEMStorage % interior % hostData(2,3,i,j,k,iVar,iEl)*&
                                                    SEMStorage % interior % hostData(3,2,i,j,k,iVar,iEl) )-&
                                                 SEMStorage % interior % hostData(2,1,i,j,k,iVar,iEl)*&
                                                  ( SEMStorage % interior % hostData(1,2,i,j,k,iVar,iEl)*&
                                                    SEMStorage % interior % hostData(3,3,i,j,k,iVar,iEl)-&
                                                    SEMStorage % interior % hostData(1,3,i,j,k,iVar,iEl)*&
                                                    SEMStorage % interior % hostData(3,2,i,j,k,iVar,iEl) )+&
                                                 SEMStorage % interior % hostData(3,1,i,j,k,iVar,iEl)*&
                                                  ( SEMStorage % interior % hostData(1,2,i,j,k,iVar,iEl)*&
                                                    SEMStorage % interior % hostData(2,3,i,j,k,iVar,iEl)-&
                                                    SEMStorage % interior % hostData(1,3,i,j,k,iVar,iEl)*&
                                                    SEMStorage % interior % hostData(2,2,i,j,k,iVar,iEl) )
                                                     

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    
END FUNCTION Determinant_SEMTensor3D

SUBROUTINE Equals_SEMTensor3D( SEMOut, SEMStorage )
  IMPLICIT NONE
  TYPE(SEMTensor3D), INTENT(out) :: SEMOut
  TYPE(SEMTensor3D), INTENT(in) :: SEMStorage

    SEMOut % interior % hostData = SEMStorage % interior % hostData
    SEMOut % boundary % hostData = SEMStorage % boundary % hostData
  
END SUBROUTINE Equals_SEMTensor3D

END MODULE NodalSEMData
