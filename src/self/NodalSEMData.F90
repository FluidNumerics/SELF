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
!
!      PROCEDURE, PUBLIC :: Set
!      PROCEDURE, PUBLIC :: SetInterior
!      PROCEDURE, PUBLIC :: SetBoundary
!      PROCEDURE, PUBLIC :: GetInterior
!      PROCEDURE, PUBLIC :: GetBoundary
!
      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMScalar1D
!      PROCEDURE, PUBLIC :: Derivative 

  END TYPE SEMScalar1D

!  TYPE, PUBLIC :: SEMScalar2D
!
!    INTEGER :: N
!    INTEGER :: nVar
!    INTEGER :: nElem
!    REAL(prec), POINTER :: fInterior(:,:,:,:)
!    REAL(prec), POINTER :: fBoundary(:,:,:,:)
!    TYPE(c_ptr) :: fInterior_dev
!    TYPE(c_ptr) :: fBoundary_dev
!
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
!      PROCEDURE, PUBLIC :: Gradient 
!
!  END TYPE SEMScalar2D
!
!  TYPE, PUBLIC :: SEMScalar3D
!
!    INTEGER :: N
!    INTEGER :: nVar
!    INTEGER :: nElem
!    REAL(prec), POINTER :: fInterior(:,:,:,:,:)
!    REAL(prec), POINTER :: fBoundary(:,:,:,:,:)
!    TYPE(c_ptr) :: fInterior_dev
!    TYPE(c_ptr) :: fBoundary_dev
!
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
!      PROCEDURE, PUBLIC :: Gradient 
!
!  END TYPE SEMScalar3D
!
!! ---------------------- Vectors ---------------------- !
!
!  TYPE, PUBLIC :: SEMVector2D
!
!    INTEGER :: N
!    INTEGER :: nVar
!    INTEGER :: nElem
!    REAL(prec), POINTER :: fInterior(:,:,:,:,:)
!    REAL(prec), POINTER :: fBoundary(:,:,:,:,:)
!    TYPE(c_ptr) :: fInterior_dev
!    TYPE(c_ptr) :: fBoundary_dev
!
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
!      PROCEDURE, PUBLIC :: Gradient 
!      PROCEDURE, PUBLIC :: Divergence
!      PROCEDURE, PUBLIC :: Curl
!
!  END TYPE SEMVector2D
!
!  TYPE, PUBLIC :: SEMVector3D
!
!    INTEGER :: N
!    INTEGER :: nVar
!    INTEGER :: nElem
!    REAL(prec), POINTER :: fInterior(:,:,:,:,:,:)
!    REAL(prec), POINTER :: fBoundary(:,:,:,:,:,:)
!    TYPE(c_ptr) :: fInterior_dev
!    TYPE(c_ptr) :: fBoundary_dev
!
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
!      PROCEDURE, PUBLIC :: Gradient 
!      PROCEDURE, PUBLIC :: Divergence
!      PROCEDURE, PUBLIC :: Curl
!
!  END TYPE SEMVector3D
!
!! ---------------------- Tensors ---------------------- !
!
!  TYPE, PUBLIC :: SEMTensor2D
!
!    INTEGER :: N
!    INTEGER :: nVar
!    INTEGER :: nElem
!    REAL(prec), POINTER :: fInterior(:,:,:,:,:,:)
!    REAL(prec), POINTER :: fBoundary(:,:,:,:,:,:)
!    TYPE(c_ptr) :: fInterior_dev
!    TYPE(c_ptr) :: fBoundary_dev
!
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
!  END TYPE SEMTensor2D
!
!  TYPE, PUBLIC :: SEMTensor3D
!
!    INTEGER :: N
!    INTEGER :: nVar
!    INTEGER :: nElem
!    REAL(prec), POINTER :: fInterior(:,:,:,:,:,:,:)
!    REAL(prec), POINTER :: fBoundary(:,:,:,:,:,:,:)
!    TYPE(c_ptr) :: fInterior_dev
!    TYPE(c_ptr) :: fBoundary_dev
!
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
!  END TYPE SEMTensor3D


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
              SEMStorage % fBoundary(1:2,1:nVar,1:nElem) )

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

END MODULE NodalSEMData
