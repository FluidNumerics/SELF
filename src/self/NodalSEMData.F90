! NodalSEMData.F90
! 
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self-fluids@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE NodalSEMData

USE ModelPrecision
USE ConstantsDictionary
USE Lagrange_Class

USE hip
USE ISO_C_BINDING


IMPLICIT NONE

! ---------------------- Scalars ---------------------- !
  TYPE, PUBLIC :: SEMScalar1D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    REAL(prec), POINTER :: fInterior(:,:,:)
    REAL(prec), POINTER :: fBoundary(:,:,:)
    TYPE(c_ptr) :: fInterior_dev
    TYPE(c_ptr) :: fBoundary_dev

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
    REAL(prec), POINTER :: fInterior(:,:,:,:)
    REAL(prec), POINTER :: fBoundary(:,:,:,:)
    TYPE(c_ptr) :: fInterior_dev
    TYPE(c_ptr) :: fBoundary_dev

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
    REAL(prec), POINTER :: fInterior(:,:,:,:,:)
    REAL(prec), POINTER :: fBoundary(:,:,:,:,:)
    TYPE(c_ptr) :: fInterior_dev
    TYPE(c_ptr) :: fBoundary_dev

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
    REAL(prec), POINTER :: fInterior(:,:,:,:,:)
    REAL(prec), POINTER :: fBoundary(:,:,:,:,:)
    TYPE(c_ptr) :: fInterior_dev
    TYPE(c_ptr) :: fBoundary_dev

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_SEMVector2D
      PROCEDURE, PUBLIC :: Trash => Trash_SEMVector2D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMVector2D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMVector2D
#endif
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMVector2D
      PROCEDURE, PUBLIC :: Gradient => Gradient_SEMVector2D
      PROCEDURE, PUBLIC :: Divergence => Divergence_SEMVector2D
      PROCEDURE, PUBLIC :: Curl => Curl_SEMVector2D

  END TYPE SEMVector2D

  TYPE, PUBLIC :: SEMVector3D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    REAL(prec), POINTER :: fInterior(:,:,:,:,:,:)
    REAL(prec), POINTER :: fBoundary(:,:,:,:,:,:)
    TYPE(c_ptr) :: fInterior_dev
    TYPE(c_ptr) :: fBoundary_dev

    CONTAINS

!      PROCEDURE, PUBLIC :: Build => Build_SEMVector3D
!      PROCEDURE, PUBLIC :: Trash => Trash_SEMVector3D
!#ifdef GPU
!      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMVector3D
!      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMVector3D
!#endif
!      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMVector3D
!      PROCEDURE, PUBLIC :: Gradient => Gradient_SEMVector3D
!      PROCEDURE, PUBLIC :: Divergence => Divergence_SEMVector3D
!      PROCEDURE, PUBLIC :: Curl => Curl_SEMVector3D

  END TYPE SEMVector3D
!
!! ---------------------- Tensors ---------------------- !
!
  TYPE, PUBLIC :: SEMTensor2D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    REAL(prec), POINTER :: fInterior(:,:,:,:,:,:)
    REAL(prec), POINTER :: fBoundary(:,:,:,:,:,:)
    TYPE(c_ptr) :: fInterior_dev
    TYPE(c_ptr) :: fBoundary_dev

!    CONTAINS
!
!      PROCEDURE, PUBLIC :: Build
!      PROCEDURE, PUBLIC :: Trash
!
!#ifdef GPU
!      PROCEDURE, PUBLIC :: UpdateHost
!      PROCEDURE, PUBLIC :: UpdateDevice
!#endif
!
!      PROCEDURE, PUBLIC :: Set
!      PROCEDURE, PUBLIC :: SetInterior
!      PROCEDURE, PUBLIC :: SetBoundary
!      PROCEDURE, PUBLIC :: GetInterior
!      PROCEDURE, PUBLIC :: GetBoundary
!
!      PROCEDURE, PUBLIC :: GridInterp 
!      PROCEDURE, PUBLIC :: Divergence
!      PROCEDURE, PUBLIC :: Rotate 
!
  END TYPE SEMTensor2D

  TYPE, PUBLIC :: SEMTensor3D

    INTEGER :: N
    INTEGER :: nVar
    INTEGER :: nElem
    REAL(prec), POINTER :: fInterior(:,:,:,:,:,:,:)
    REAL(prec), POINTER :: fBoundary(:,:,:,:,:,:,:)
    TYPE(c_ptr) :: fInterior_dev
    TYPE(c_ptr) :: fBoundary_dev

!    CONTAINS
!
!      PROCEDURE, PUBLIC :: Build
!      PROCEDURE, PUBLIC :: Trash
!
!#ifdef GPU
!      PROCEDURE, PUBLIC :: UpdateHost
!      PROCEDURE, PUBLIC :: UpdateDevice
!#endif
!
!      PROCEDURE, PUBLIC :: Set
!      PROCEDURE, PUBLIC :: SetInterior
!      PROCEDURE, PUBLIC :: SetBoundary
!      PROCEDURE, PUBLIC :: GetInterior
!      PROCEDURE, PUBLIC :: GetBoundary
!
!      PROCEDURE, PUBLIC :: GridInterp 
!      PROCEDURE, PUBLIC :: Divergence
!      PROCEDURE, PUBLIC :: Rotate 
!
  END TYPE SEMTensor3D


CONTAINS

! -- SEMScalar1D -- !

SUBROUTINE Build_SEMScalar1D( SEMStorage, N, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(SEMScalar1D), INTENT(out) :: SEMStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SEMStorage % N = N
    SEMStorage % nVar = nVar
    SEMStorage % nElem = nElem

    ALLOCATE( SEMStorage % fInterior(0:N,1:nVar,1:nElem), &
              SEMStorage % fBoundary(1:nVar,1:2,1:nElem) )

#ifdef GPU
    CALL hipCheck(hipMalloc(SEMStorage % fInterior_dev, SIZEOF(SEMStorage % fInterior)))
    CALL hipCheck(hipMalloc(SEMStorage % fBoundary_dev, SIZEOF(SEMStorage % fBoundary)))
#endif


END SUBROUTINE Build_SEMScalar1D

SUBROUTINE Trash_SEMScalar1D( SEMStorage ) 
  IMPLICIT NONE
  CLASS(SEMScalar1D), INTENT(inout) :: SEMStorage

    DEALLOCATE( SEMStorage % fInterior, &
                SEMStorage % fBoundary )
#ifdef GPU
    CALL hipCheck(hipFree(SEMStorage % fInterior_dev))
    CALL hipCheck(hipFree(SEMStorage % fBoundary_dev))
#endif

END SUBROUTINE Trash_SEMScalar1D

#ifdef GPU
SUBROUTINE UpdateHost_SEMScalar1D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMScalar1D), INTENT(inout) :: SEMStorage

     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % fInterior), &
                             SEMStorage % fInterior_dev, &
                             SIZEOF(SEMStorage % fInterior), &
                             hipMemcpyDeviceToHost))

     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % fBoundary), &
                             SEMStorage % fBoundary_dev, &
                             SIZEOF(SEMStorage % fBoundary), &
                             hipMemcpyDeviceToHost))

END SUBROUTINE UpdateHost_SEMScalar1D

SUBROUTINE UpdateDevice_SEMScalar1D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMScalar1D), INTENT(inout) :: SEMStorage

     CALL hipCheck(hipMemcpy(SEMStorage % fInterior_dev, &
                             c_loc(SEMStorage % fInterior), &
                             SIZEOF(SEMStorage % fInterior), &
                             hipMemcpyHostToDevice))

     CALL hipCheck(hipMemcpy(SEMStorage % fBoundary_dev, &
                             c_loc(SEMStorage % fBoundary), &
                             SIZEOF(SEMStorage % fBoundary), &
                             hipMemcpyHostToDevice))

END SUBROUTINE UpdateDevice_SEMScalar1D
#endif

FUNCTION BoundaryInterp_SEMScalar1D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar1D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMScalar1D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % ScalarBoundaryInterp_1D( SEMStorage % fInterior_dev, &
                                         SEMout % fInterior_dev, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ELSE
      CALL interp % ScalarBoundaryInterp_1D( SEMStorage % fInterior, &
                                         SEMout % fInterior, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_SEMScalar1D

FUNCTION GridInterp_SEMScalar1D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar1D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMScalar1D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % ScalarGridInterp_1D( SEMStorage % fInterior_dev, &
                                         SEMout % fInterior_dev, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ELSE
      CALL interp % ScalarGridInterp_1D( SEMStorage % fInterior, &
                                         SEMout % fInterior, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_SEMScalar1D

FUNCTION Derivative_SEMScalar1D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar1D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMScalar1D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % Derivative_1D( SEMStorage % fInterior_dev, &
                                   SEMout % fInterior_dev, &
                                   SEMStorage % nVar, &
                                   SEMStorage % nElem )  
    ELSE
      CALL interp % Derivative_1D( SEMStorage % fInterior, &
                                   SEMout % fInterior, &
                                   SEMStorage % nVar, &
                                   SEMStorage % nElem )  
    ENDIF

END FUNCTION Derivative_SEMScalar1D

! -- SEMScalar2D -- !

SUBROUTINE Build_SEMScalar2D( SEMStorage, N, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(SEMScalar2D), INTENT(out) :: SEMStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SEMStorage % N = N
    SEMStorage % nVar = nVar
    SEMStorage % nElem = nElem

    ALLOCATE( SEMStorage % fInterior(0:N,0:N,1:nVar,1:nElem), &
              SEMStorage % fBoundary(0:N,1:nVar,1:4,1:nElem) )

#ifdef GPU
    CALL hipCheck(hipMalloc(SEMStorage % fInterior_dev, SIZEOF(SEMStorage % fInterior)))
    CALL hipCheck(hipMalloc(SEMStorage % fBoundary_dev, SIZEOF(SEMStorage % fBoundary)))
#endif


END SUBROUTINE Build_SEMScalar2D

SUBROUTINE Trash_SEMScalar2D( SEMStorage ) 
  IMPLICIT NONE
  CLASS(SEMScalar2D), INTENT(inout) :: SEMStorage

    DEALLOCATE( SEMStorage % fInterior, &
                SEMStorage % fBoundary )
#ifdef GPU
    CALL hipCheck(hipFree(SEMStorage % fInterior_dev))
    CALL hipCheck(hipFree(SEMStorage % fBoundary_dev))
#endif

END SUBROUTINE Trash_SEMScalar2D

#ifdef GPU
SUBROUTINE UpdateHost_SEMScalar2D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMScalar2D), INTENT(inout) :: SEMStorage

     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % fInterior), &
                             SEMStorage % fInterior_dev, &
                             SIZEOF(SEMStorage % fInterior), &
                             hipMemcpyDeviceToHost))

     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % fBoundary), &
                             SEMStorage % fBoundary_dev, &
                             SIZEOF(SEMStorage % fBoundary), &
                             hipMemcpyDeviceToHost))

END SUBROUTINE UpdateHost_SEMScalar2D

SUBROUTINE UpdateDevice_SEMScalar2D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMScalar2D), INTENT(inout) :: SEMStorage

     CALL hipCheck(hipMemcpy(SEMStorage % fInterior_dev, &
                             c_loc(SEMStorage % fInterior), &
                             SIZEOF(SEMStorage % fInterior), &
                             hipMemcpyHostToDevice))

     CALL hipCheck(hipMemcpy(SEMStorage % fBoundary_dev, &
                             c_loc(SEMStorage % fBoundary), &
                             SIZEOF(SEMStorage % fBoundary), &
                             hipMemcpyHostToDevice))

END SUBROUTINE UpdateDevice_SEMScalar2D
#endif

FUNCTION BoundaryInterp_SEMScalar2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar2D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMScalar2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % ScalarBoundaryInterp_2D( SEMStorage % fInterior_dev, &
                                         SEMout % fInterior_dev, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ELSE
      CALL interp % ScalarBoundaryInterp_2D( SEMStorage % fInterior, &
                                         SEMout % fInterior, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_SEMScalar2D

FUNCTION GridInterp_SEMScalar2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar2D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMScalar2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % ScalarGridInterp_2D( SEMStorage % fInterior_dev, &
                                         SEMout % fInterior_dev, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ELSE
      CALL interp % ScalarGridInterp_2D( SEMStorage % fInterior, &
                                         SEMout % fInterior, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_SEMScalar2D

FUNCTION Gradient_SEMScalar2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar2D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMVector2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % ScalarGradient_2D( SEMStorage % fInterior_dev, &
                                       SEMout % fInterior_dev, &
                                       SEMStorage % nVar, &
                                       SEMStorage % nElem )  
    ELSE
      CALL interp % ScalarGradient_2D( SEMStorage % fInterior, &
                                       SEMout % fInterior, &
                                       SEMStorage % nVar, &
                                       SEMStorage % nElem )  
    ENDIF

END FUNCTION Gradient_SEMScalar2D

! -- SEMScalar3D -- !

SUBROUTINE Build_SEMScalar3D( SEMStorage, N, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(SEMScalar3D), INTENT(out) :: SEMStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SEMStorage % N = N
    SEMStorage % nVar = nVar
    SEMStorage % nElem = nElem

    ALLOCATE( SEMStorage % fInterior(0:N,0:N,0:N,1:nVar,1:nElem), &
              SEMStorage % fBoundary(0:N,0:N,1:nVar,1:6,1:nElem) )

#ifdef GPU
    CALL hipCheck(hipMalloc(SEMStorage % fInterior_dev, SIZEOF(SEMStorage % fInterior)))
    CALL hipCheck(hipMalloc(SEMStorage % fBoundary_dev, SIZEOF(SEMStorage % fBoundary)))
#endif

END SUBROUTINE Build_SEMScalar3D

SUBROUTINE Trash_SEMScalar3D( SEMStorage ) 
  IMPLICIT NONE
  CLASS(SEMScalar3D), INTENT(inout) :: SEMStorage

    DEALLOCATE( SEMStorage % fInterior, &
                SEMStorage % fBoundary )
#ifdef GPU
    CALL hipCheck(hipFree(SEMStorage % fInterior_dev))
    CALL hipCheck(hipFree(SEMStorage % fBoundary_dev))
#endif

END SUBROUTINE Trash_SEMScalar3D

#ifdef GPU
SUBROUTINE UpdateHost_SEMScalar3D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMScalar3D), INTENT(inout) :: SEMStorage

     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % fInterior), &
                             SEMStorage % fInterior_dev, &
                             SIZEOF(SEMStorage % fInterior), &
                             hipMemcpyDeviceToHost))

     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % fBoundary), &
                             SEMStorage % fBoundary_dev, &
                             SIZEOF(SEMStorage % fBoundary), &
                             hipMemcpyDeviceToHost))

END SUBROUTINE UpdateHost_SEMScalar3D

SUBROUTINE UpdateDevice_SEMScalar3D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMScalar3D), INTENT(inout) :: SEMStorage

     CALL hipCheck(hipMemcpy(SEMStorage % fInterior_dev, &
                             c_loc(SEMStorage % fInterior), &
                             SIZEOF(SEMStorage % fInterior), &
                             hipMemcpyHostToDevice))

     CALL hipCheck(hipMemcpy(SEMStorage % fBoundary_dev, &
                             c_loc(SEMStorage % fBoundary), &
                             SIZEOF(SEMStorage % fBoundary), &
                             hipMemcpyHostToDevice))

END SUBROUTINE UpdateDevice_SEMScalar3D
#endif

FUNCTION BoundaryInterp_SEMScalar3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar3D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMScalar3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % ScalarBoundaryInterp_3D( SEMStorage % fInterior_dev, &
                                         SEMout % fInterior_dev, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ELSE
      CALL interp % ScalarBoundaryInterp_3D( SEMStorage % fInterior, &
                                         SEMout % fInterior, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_SEMScalar3D

FUNCTION GridInterp_SEMScalar3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar3D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMScalar3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % ScalarGridInterp_3D( SEMStorage % fInterior_dev, &
                                         SEMout % fInterior_dev, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ELSE
      CALL interp % ScalarGridInterp_3D( SEMStorage % fInterior, &
                                         SEMout % fInterior, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_SEMScalar3D

FUNCTION Gradient_SEMScalar3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMScalar3D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMVector3D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % ScalarGradient_3D( SEMStorage % fInterior_dev, &
                                       SEMout % fInterior_dev, &
                                       SEMStorage % nVar, &
                                       SEMStorage % nElem )  
    ELSE
      CALL interp % ScalarGradient_3D( SEMStorage % fInterior, &
                                       SEMout % fInterior, &
                                       SEMStorage % nVar, &
                                       SEMStorage % nElem )  
    ENDIF

END FUNCTION Gradient_SEMScalar3D

! -- SEMVector2D -- !

SUBROUTINE Build_SEMVector2D( SEMStorage, N, nVar, nElem ) 
  IMPLICIT NONE
  CLASS(SEMVector2D), INTENT(out) :: SEMStorage
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: nVar
  INTEGER, INTENT(in) :: nElem

    SEMStorage % N = N
    SEMStorage % nVar = nVar
    SEMStorage % nElem = nElem

    ALLOCATE( SEMStorage % fInterior(1:2,0:N,0:N,1:nVar,1:nElem), &
              SEMStorage % fBoundary(1:2,0:N,1:nVar,1:4,1:nElem) )

#ifdef GPU
    CALL hipCheck(hipMalloc(SEMStorage % fInterior_dev, SIZEOF(SEMStorage % fInterior)))
    CALL hipCheck(hipMalloc(SEMStorage % fBoundary_dev, SIZEOF(SEMStorage % fBoundary)))
#endif


END SUBROUTINE Build_SEMVector2D

SUBROUTINE Trash_SEMVector2D( SEMStorage ) 
  IMPLICIT NONE
  CLASS(SEMVector2D), INTENT(inout) :: SEMStorage

    DEALLOCATE( SEMStorage % fInterior, &
                SEMStorage % fBoundary )
#ifdef GPU
    CALL hipCheck(hipFree(SEMStorage % fInterior_dev))
    CALL hipCheck(hipFree(SEMStorage % fBoundary_dev))
#endif

END SUBROUTINE Trash_SEMVector2D

#ifdef GPU
SUBROUTINE UpdateHost_SEMVector2D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMVector2D), INTENT(inout) :: SEMStorage

     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % fInterior), &
                             SEMStorage % fInterior_dev, &
                             SIZEOF(SEMStorage % fInterior), &
                             hipMemcpyDeviceToHost))

     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % fBoundary), &
                             SEMStorage % fBoundary_dev, &
                             SIZEOF(SEMStorage % fBoundary), &
                             hipMemcpyDeviceToHost))

END SUBROUTINE UpdateHost_SEMVector2D

SUBROUTINE UpdateDevice_SEMVector2D( SEMStorage )
  IMPLICIT NONE
  CLASS(SEMVector2D), INTENT(inout) :: SEMStorage

     CALL hipCheck(hipMemcpy(SEMStorage % fInterior_dev, &
                             c_loc(SEMStorage % fInterior), &
                             SIZEOF(SEMStorage % fInterior), &
                             hipMemcpyHostToDevice))

     CALL hipCheck(hipMemcpy(SEMStorage % fBoundary_dev, &
                             c_loc(SEMStorage % fBoundary), &
                             SIZEOF(SEMStorage % fBoundary), &
                             hipMemcpyHostToDevice))

END SUBROUTINE UpdateDevice_SEMVector2D
#endif

FUNCTION BoundaryInterp_SEMVector2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector2D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMVector2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % VectorBoundaryInterp_2D( SEMStorage % fInterior_dev, &
                                         SEMout % fInterior_dev, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ELSE
      CALL interp % VectorBoundaryInterp_2D( SEMStorage % fInterior, &
                                         SEMout % fInterior, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ENDIF

END FUNCTION BoundaryInterp_SEMVector2D

FUNCTION GridInterp_SEMVector2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector2D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMVector2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % VectorGridInterp_2D( SEMStorage % fInterior_dev, &
                                         SEMout % fInterior_dev, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ELSE
      CALL interp % VectorGridInterp_2D( SEMStorage % fInterior, &
                                         SEMout % fInterior, &
                                         SEMStorage % nVar, &
                                         SEMStorage % nElem )  
    ENDIF

END FUNCTION GridInterp_SEMVector2D

FUNCTION Gradient_SEMVector2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector2D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMTensor2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % VectorGradient_2D( SEMStorage % fInterior_dev, &
                                       SEMout % fInterior_dev, &
                                       SEMStorage % nVar, &
                                       SEMStorage % nElem )  
    ELSE
      CALL interp % VectorGradient_2D( SEMStorage % fInterior, &
                                       SEMout % fInterior, &
                                       SEMStorage % nVar, &
                                       SEMStorage % nElem )  
    ENDIF

END FUNCTION Gradient_SEMVector2D

FUNCTION Divergence_SEMVector2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector2D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMScalar2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % VectorDivergence_2D( SEMStorage % fInterior_dev, &
                                       SEMout % fInterior_dev, &
                                       SEMStorage % nVar, &
                                       SEMStorage % nElem )  
    ELSE
      CALL interp % VectorDivergence_2D( SEMStorage % fInterior, &
                                       SEMout % fInterior, &
                                       SEMStorage % nVar, &
                                       SEMStorage % nElem )  
    ENDIF

END FUNCTION Divergence_SEMVector2D

FUNCTION Curl_SEMVector2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
  IMPLICIT NONE
  CLASS(SEMVector2D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMScalar2D) :: SEMOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
      CALL interp % VectorCurl_2D( SEMStorage % fInterior_dev, &
                                       SEMout % fInterior_dev, &
                                       SEMStorage % nVar, &
                                       SEMStorage % nElem )  
    ELSE
      CALL interp % VectorCurl_2D( SEMStorage % fInterior, &
                                       SEMout % fInterior, &
                                       SEMStorage % nVar, &
                                       SEMStorage % nElem )  
    ENDIF

END FUNCTION Curl_SEMVector2D

END MODULE NodalSEMData
