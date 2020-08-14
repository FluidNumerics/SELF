! MappedNodalSEMData.F90
! 
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self-fluids@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE MappedNodalSEMData

USE ModelPrecision
USE ConstantsDictionary
USE Lagrange_Class

USE NodalSEMData
USE SEMMesh
!USE MappedSEM

USE hip
USE ISO_C_BINDING


IMPLICIT NONE

!! ---------------------- Scalars ---------------------- !
!  TYPE, EXTENDS(SEMScalar1D), PUBLIC :: MappedSEMScalar1D
!
!    CONTAINS
!
!      PROCEDURE, PUBLIC :: MappedDerivative => MappedDerivative_SEMScalar1D
!
!  END TYPE SEMScalar1D
!
  TYPE, EXTENDS(SEMScalar2D), PUBLIC :: MappedSEMScalar2D

    CONTAINS

      PROCEDURE, PUBLIC :: Gradient => Gradient_MappedSEMScalar2D

  END TYPE MappedSEMScalar2D
!
!  TYPE, PUBLIC :: SEMScalar3D
!
!    INTEGER :: N
!    INTEGER :: nVar
!    INTEGER :: nElem
!    REAL(prec), POINTER :: interior(:,:,:,:,:)
!    REAL(prec), POINTER :: boundary(:,:,:,:,:)
!    TYPE(c_ptr) :: interior_dev
!    TYPE(c_ptr) :: boundary_dev
!
!    CONTAINS
!
!      PROCEDURE, PUBLIC :: Build => Build_SEMScalar3D
!      PROCEDURE, PUBLIC :: Trash => Trash_SEMScalar3D
!#ifdef GPU
!      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMScalar3D
!      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMScalar3D
!#endif
!      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_SEMScalar3D
!      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMScalar3D
!      PROCEDURE, PUBLIC :: Gradient => Gradient_SEMScalar3D
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
!    REAL(prec), POINTER :: interior(:,:,:,:,:)
!    REAL(prec), POINTER :: boundary(:,:,:,:,:)
!    TYPE(c_ptr) :: interior_dev
!    TYPE(c_ptr) :: boundary_dev
!
!    CONTAINS
!
!      PROCEDURE, PUBLIC :: Build => Build_SEMVector2D
!      PROCEDURE, PUBLIC :: Trash => Trash_SEMVector2D
!#ifdef GPU
!      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMVector2D
!      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMVector2D
!#endif
!      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_SEMVector2D
!      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMVector2D
!      PROCEDURE, PUBLIC :: Gradient => Gradient_SEMVector2D
!      PROCEDURE, PUBLIC :: Divergence => Divergence_SEMVector2D
!      PROCEDURE, PUBLIC :: Curl => Curl_SEMVector2D
!
!
!  END TYPE SEMVector2D
!
!  TYPE, PUBLIC :: SEMVector3D
!
!    INTEGER :: N
!    INTEGER :: nVar
!    INTEGER :: nElem
!    REAL(prec), POINTER :: interior(:,:,:,:,:,:)
!    REAL(prec), POINTER :: boundary(:,:,:,:,:,:)
!    TYPE(c_ptr) :: interior_dev
!    TYPE(c_ptr) :: boundary_dev
!
!    CONTAINS
!
!      PROCEDURE, PUBLIC :: Build => Build_SEMVector3D
!      PROCEDURE, PUBLIC :: Trash => Trash_SEMVector3D
!#ifdef GPU
!      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMVector3D
!      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMVector3D
!#endif
!      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_SEMVector3D
!      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMVector3D
!      PROCEDURE, PUBLIC :: Gradient => Gradient_SEMVector3D
!      PROCEDURE, PUBLIC :: Divergence => Divergence_SEMVector3D
!      PROCEDURE, PUBLIC :: Curl => Curl_SEMVector3D
!
!
!  END TYPE SEMVector3D
!!
!!! ---------------------- Tensors ---------------------- !
!!
!  TYPE, PUBLIC :: SEMTensor2D
!
!    INTEGER :: N
!    INTEGER :: nVar
!    INTEGER :: nElem
!    REAL(prec), POINTER :: interior(:,:,:,:,:,:)
!    REAL(prec), POINTER :: boundary(:,:,:,:,:,:)
!    TYPE(c_ptr) :: interior_dev
!    TYPE(c_ptr) :: boundary_dev
!
!    CONTAINS
!
!      PROCEDURE, PUBLIC :: Build => Build_SEMTensor2D
!      PROCEDURE, PUBLIC :: Trash => Trash_SEMTensor2D
!#ifdef GPU
!      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMTensor2D
!      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMTensor2D
!#endif
!      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_SEMTensor2D
!      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMTensor2D 
!
!      PROCEDURE, PUBLIC :: Determinant => Determinant_SEMTensor2D
!
!  END TYPE SEMTensor2D
!
!  TYPE, PUBLIC :: SEMTensor3D
!
!    INTEGER :: N
!    INTEGER :: nVar
!    INTEGER :: nElem
!    REAL(prec), POINTER :: interior(:,:,:,:,:,:,:)
!    REAL(prec), POINTER :: boundary(:,:,:,:,:,:,:)
!    TYPE(c_ptr) :: interior_dev
!    TYPE(c_ptr) :: boundary_dev
!
!    CONTAINS
!
!      PROCEDURE, PUBLIC :: Build => Build_SEMTensor3D
!      PROCEDURE, PUBLIC :: Trash => Trash_SEMTensor3D
!#ifdef GPU
!      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMTensor3D
!      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMTensor3D
!#endif
!      PROCEDURE, PUBLIC :: BoundaryInterp => BoundaryInterp_SEMTensor3D
!      PROCEDURE, PUBLIC :: GridInterp => GridInterp_SEMTensor3D 
!
!      PROCEDURE, PUBLIC :: Determinant => Determinant_SEMTensor3D
!
!  END TYPE SEMTensor3D
!
!  INTERFACE ASSIGNMENT (=)
!    MODULE PROCEDURE Equals_SEMScalar1D
!    MODULE PROCEDURE Equals_SEMScalar2D
!    MODULE PROCEDURE Equals_SEMScalar3D
!    MODULE PROCEDURE Equals_SEMVector2D
!    MODULE PROCEDURE Equals_SEMVector3D
!    MODULE PROCEDURE Equals_SEMTensor2D
!    MODULE PROCEDURE Equals_SEMTensor3D
!  END INTERFACE
!
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
!
CONTAINS
!
!! -- SEMScalar1D -- !
!
FUNCTION MappedGradient_SEMScalar2D( SEMStorage, mesh ) RESULT( SEMOut )
  IMPLICIT NONE
  CLASS(SEMScalar2D) :: SEMStorage
  TYPE(Lagrange) :: interp
  TYPE(SEMVector2D) :: SEMOut
  !LOGICAL, OPTIONAL :: gpuAccel

  
 
  
END FUNCTION MappedGradient_SEMScalar2D
!SUBROUTINE Build_SEMScalar1D( SEMStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(SEMScalar1D), INTENT(out) :: SEMStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SEMStorage % N = N
!    SEMStorage % nVar = nVar
!    SEMStorage % nElem = nElem
!
!    ALLOCATE( SEMStorage % interior(0:N,1:nVar,1:nElem), &
!              SEMStorage % boundary(1:nVar,1:2,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SEMStorage % interior_dev, SIZEOF(SEMStorage % interior)))
!    CALL hipCheck(hipMalloc(SEMStorage % boundary_dev, SIZEOF(SEMStorage % boundary)))
!#endif
!
!
!END SUBROUTINE Build_SEMScalar1D
!
!SUBROUTINE Trash_SEMScalar1D( SEMStorage ) 
!  IMPLICIT NONE
!  CLASS(SEMScalar1D), INTENT(inout) :: SEMStorage
!
!    DEALLOCATE( SEMStorage % interior, &
!                SEMStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SEMStorage % interior_dev))
!    CALL hipCheck(hipFree(SEMStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_SEMScalar1D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_SEMScalar1D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMScalar1D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % interior), &
!                             SEMStorage % interior_dev, &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % boundary), &
!                             SEMStorage % boundary_dev, &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_SEMScalar1D
!
!SUBROUTINE UpdateDevice_SEMScalar1D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMScalar1D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(SEMStorage % interior_dev, &
!                             c_loc(SEMStorage % interior), &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SEMStorage % boundary_dev, &
!                             c_loc(SEMStorage % boundary), &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_SEMScalar1D
!#endif
!
!FUNCTION BoundaryInterp_SEMScalar1D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMScalar1D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMScalar1D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarBoundaryInterp_1D( SEMStorage % interior_dev, &
!                                         SEMout % boundary_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % ScalarBoundaryInterp_1D( SEMStorage % interior, &
!                                         SEMout % boundary, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_SEMScalar1D
!
!FUNCTION GridInterp_SEMScalar1D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMScalar1D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMScalar1D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarGridInterp_1D( SEMStorage % interior_dev, &
!                                         SEMout % interior_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % ScalarGridInterp_1D( SEMStorage % interior, &
!                                         SEMout % interior, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_SEMScalar1D
!
!FUNCTION Derivative_SEMScalar1D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMScalar1D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMScalar1D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % Derivative_1D( SEMStorage % interior_dev, &
!                                   SEMout % interior_dev, &
!                                   SEMStorage % nVar, &
!                                   SEMStorage % nElem )  
!    ELSE
!      CALL interp % Derivative_1D( SEMStorage % interior, &
!                                   SEMout % interior, &
!                                   SEMStorage % nVar, &
!                                   SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION Derivative_SEMScalar1D
!
!SUBROUTINE Equals_SEMScalar1D( SEMOut, SEMStorage )
!  IMPLICIT NONE
!  TYPE(SEMScalar1D), INTENT(out) :: SEMOut
!  TYPE(SEMScalar1D), INTENT(in) :: SEMStorage
!
!    SEMOut % interior = SEMStorage % interior
!    SEMOut % boundary = SEMStorage % boundary
!  
!END SUBROUTINE Equals_SEMScalar1D
!
!! -- SEMScalar2D -- !
!
!SUBROUTINE Build_SEMScalar2D( SEMStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(SEMScalar2D), INTENT(out) :: SEMStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SEMStorage % N = N
!    SEMStorage % nVar = nVar
!    SEMStorage % nElem = nElem
!
!    ALLOCATE( SEMStorage % interior(0:N,0:N,1:nVar,1:nElem), &
!              SEMStorage % boundary(0:N,1:nVar,1:4,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SEMStorage % interior_dev, SIZEOF(SEMStorage % interior)))
!    CALL hipCheck(hipMalloc(SEMStorage % boundary_dev, SIZEOF(SEMStorage % boundary)))
!#endif
!
!
!END SUBROUTINE Build_SEMScalar2D
!
!SUBROUTINE Trash_SEMScalar2D( SEMStorage ) 
!  IMPLICIT NONE
!  CLASS(SEMScalar2D), INTENT(inout) :: SEMStorage
!
!    DEALLOCATE( SEMStorage % interior, &
!                SEMStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SEMStorage % interior_dev))
!    CALL hipCheck(hipFree(SEMStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_SEMScalar2D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_SEMScalar2D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMScalar2D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % interior), &
!                             SEMStorage % interior_dev, &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % boundary), &
!                             SEMStorage % boundary_dev, &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_SEMScalar2D
!
!SUBROUTINE UpdateDevice_SEMScalar2D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMScalar2D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(SEMStorage % interior_dev, &
!                             c_loc(SEMStorage % interior), &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SEMStorage % boundary_dev, &
!                             c_loc(SEMStorage % boundary), &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_SEMScalar2D
!#endif
!
!FUNCTION BoundaryInterp_SEMScalar2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMScalar2D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMScalar2D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarBoundaryInterp_2D( SEMStorage % interior_dev, &
!                                         SEMout % boundary_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % ScalarBoundaryInterp_2D( SEMStorage % interior, &
!                                         SEMout % boundary, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_SEMScalar2D
!
!FUNCTION GridInterp_SEMScalar2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMScalar2D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMScalar2D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarGridInterp_2D( SEMStorage % interior_dev, &
!                                         SEMout % interior_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % ScalarGridInterp_2D( SEMStorage % interior, &
!                                         SEMout % interior, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_SEMScalar2D
!
!FUNCTION Gradient_SEMScalar2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMScalar2D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMVector2D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarGradient_2D( SEMStorage % interior_dev, &
!                                       SEMout % interior_dev, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ELSE
!      CALL interp % ScalarGradient_2D( SEMStorage % interior, &
!                                       SEMout % interior, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION Gradient_SEMScalar2D
!
!SUBROUTINE Equals_SEMScalar2D( SEMOut, SEMStorage )
!  IMPLICIT NONE
!  TYPE(SEMScalar2D), INTENT(out) :: SEMOut
!  TYPE(SEMScalar2D), INTENT(in) :: SEMStorage
!
!    SEMOut % interior = SEMStorage % interior
!    SEMOut % boundary = SEMStorage % boundary
!  
!END SUBROUTINE Equals_SEMScalar2D
!
!! -- SEMScalar3D -- !
!
!SUBROUTINE Build_SEMScalar3D( SEMStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(SEMScalar3D), INTENT(out) :: SEMStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SEMStorage % N = N
!    SEMStorage % nVar = nVar
!    SEMStorage % nElem = nElem
!
!    ALLOCATE( SEMStorage % interior(0:N,0:N,0:N,1:nVar,1:nElem), &
!              SEMStorage % boundary(0:N,0:N,1:nVar,1:6,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SEMStorage % interior_dev, SIZEOF(SEMStorage % interior)))
!    CALL hipCheck(hipMalloc(SEMStorage % boundary_dev, SIZEOF(SEMStorage % boundary)))
!#endif
!
!END SUBROUTINE Build_SEMScalar3D
!
!SUBROUTINE Trash_SEMScalar3D( SEMStorage ) 
!  IMPLICIT NONE
!  CLASS(SEMScalar3D), INTENT(inout) :: SEMStorage
!
!    DEALLOCATE( SEMStorage % interior, &
!                SEMStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SEMStorage % interior_dev))
!    CALL hipCheck(hipFree(SEMStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_SEMScalar3D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_SEMScalar3D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMScalar3D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % interior), &
!                             SEMStorage % interior_dev, &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % boundary), &
!                             SEMStorage % boundary_dev, &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_SEMScalar3D
!
!SUBROUTINE UpdateDevice_SEMScalar3D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMScalar3D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(SEMStorage % interior_dev, &
!                             c_loc(SEMStorage % interior), &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SEMStorage % boundary_dev, &
!                             c_loc(SEMStorage % boundary), &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_SEMScalar3D
!#endif
!
!FUNCTION BoundaryInterp_SEMScalar3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMScalar3D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMScalar3D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarBoundaryInterp_3D( SEMStorage % interior_dev, &
!                                         SEMout % boundary_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % ScalarBoundaryInterp_3D( SEMStorage % interior, &
!                                         SEMout % boundary, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_SEMScalar3D
!
!FUNCTION GridInterp_SEMScalar3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMScalar3D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMScalar3D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarGridInterp_3D( SEMStorage % interior_dev, &
!                                         SEMout % interior_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % ScalarGridInterp_3D( SEMStorage % interior, &
!                                         SEMout % interior, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_SEMScalar3D
!
!FUNCTION Gradient_SEMScalar3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMScalar3D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMVector3D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarGradient_3D( SEMStorage % interior_dev, &
!                                       SEMout % interior_dev, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ELSE
!      CALL interp % ScalarGradient_3D( SEMStorage % interior, &
!                                       SEMout % interior, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION Gradient_SEMScalar3D
!
!SUBROUTINE Equals_SEMScalar3D( SEMOut, SEMStorage )
!  IMPLICIT NONE
!  TYPE(SEMScalar3D), INTENT(out) :: SEMOut
!  TYPE(SEMScalar3D), INTENT(in) :: SEMStorage
!
!    SEMOut % interior = SEMStorage % interior
!    SEMOut % boundary = SEMStorage % boundary
!  
!END SUBROUTINE Equals_SEMScalar3D
!
!! -- SEMVector2D -- !
!
!SUBROUTINE Build_SEMVector2D( SEMStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(SEMVector2D), INTENT(out) :: SEMStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SEMStorage % N = N
!    SEMStorage % nVar = nVar
!    SEMStorage % nElem = nElem
!
!    ALLOCATE( SEMStorage % interior(1:2,0:N,0:N,1:nVar,1:nElem), &
!              SEMStorage % boundary(1:2,0:N,1:nVar,1:4,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SEMStorage % interior_dev, SIZEOF(SEMStorage % interior)))
!    CALL hipCheck(hipMalloc(SEMStorage % boundary_dev, SIZEOF(SEMStorage % boundary)))
!#endif
!
!
!END SUBROUTINE Build_SEMVector2D
!
!SUBROUTINE Trash_SEMVector2D( SEMStorage ) 
!  IMPLICIT NONE
!  CLASS(SEMVector2D), INTENT(inout) :: SEMStorage
!
!    DEALLOCATE( SEMStorage % interior, &
!                SEMStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SEMStorage % interior_dev))
!    CALL hipCheck(hipFree(SEMStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_SEMVector2D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_SEMVector2D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMVector2D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % interior), &
!                             SEMStorage % interior_dev, &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % boundary), &
!                             SEMStorage % boundary_dev, &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_SEMVector2D
!
!SUBROUTINE UpdateDevice_SEMVector2D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMVector2D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(SEMStorage % interior_dev, &
!                             c_loc(SEMStorage % interior), &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SEMStorage % boundary_dev, &
!                             c_loc(SEMStorage % boundary), &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_SEMVector2D
!#endif
!
!FUNCTION BoundaryInterp_SEMVector2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMVector2D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMVector2D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorBoundaryInterp_2D( SEMStorage % interior_dev, &
!                                         SEMout % boundary_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % VectorBoundaryInterp_2D( SEMStorage % interior, &
!                                         SEMout % boundary, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_SEMVector2D
!
!FUNCTION GridInterp_SEMVector2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMVector2D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMVector2D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorGridInterp_2D( SEMStorage % interior_dev, &
!                                         SEMout % interior_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % VectorGridInterp_2D( SEMStorage % interior, &
!                                         SEMout % interior, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_SEMVector2D
!
!FUNCTION Gradient_SEMVector2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMVector2D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMTensor2D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorGradient_2D( SEMStorage % interior_dev, &
!                                       SEMout % interior_dev, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ELSE
!      CALL interp % VectorGradient_2D( SEMStorage % interior, &
!                                       SEMout % interior, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION Gradient_SEMVector2D
!
!FUNCTION Divergence_SEMVector2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMVector2D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMScalar2D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorDivergence_2D( SEMStorage % interior_dev, &
!                                       SEMout % interior_dev, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ELSE
!      CALL interp % VectorDivergence_2D( SEMStorage % interior, &
!                                       SEMout % interior, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION Divergence_SEMVector2D
!
!FUNCTION Curl_SEMVector2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMVector2D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMScalar2D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorCurl_2D( SEMStorage % interior_dev, &
!                                       SEMout % interior_dev, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ELSE
!      CALL interp % VectorCurl_2D( SEMStorage % interior, &
!                                       SEMout % interior, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION Curl_SEMVector2D
!
!SUBROUTINE Equals_SEMVector2D( SEMOut, SEMStorage )
!  IMPLICIT NONE
!  TYPE(SEMVector2D), INTENT(out) :: SEMOut
!  TYPE(SEMVector2D), INTENT(in) :: SEMStorage
!
!    SEMOut % interior = SEMStorage % interior
!    SEMOut % boundary = SEMStorage % boundary
!  
!END SUBROUTINE Equals_SEMVector2D
!
!! -- SEMVector3D -- !
!
!SUBROUTINE Build_SEMVector3D( SEMStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(SEMVector3D), INTENT(out) :: SEMStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SEMStorage % N = N
!    SEMStorage % nVar = nVar
!    SEMStorage % nElem = nElem
!
!    ALLOCATE( SEMStorage % interior(1:3,0:N,0:N,0:N,1:nVar,1:nElem), &
!              SEMStorage % boundary(1:3,0:N,0:N,1:nVar,1:6,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SEMStorage % interior_dev, SIZEOF(SEMStorage % interior)))
!    CALL hipCheck(hipMalloc(SEMStorage % boundary_dev, SIZEOF(SEMStorage % boundary)))
!#endif
!
!
!END SUBROUTINE Build_SEMVector3D
!
!SUBROUTINE Trash_SEMVector3D( SEMStorage ) 
!  IMPLICIT NONE
!  CLASS(SEMVector3D), INTENT(inout) :: SEMStorage
!
!    DEALLOCATE( SEMStorage % interior, &
!                SEMStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SEMStorage % interior_dev))
!    CALL hipCheck(hipFree(SEMStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_SEMVector3D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_SEMVector3D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMVector3D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % interior), &
!                             SEMStorage % interior_dev, &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % boundary), &
!                             SEMStorage % boundary_dev, &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_SEMVector3D
!
!SUBROUTINE UpdateDevice_SEMVector3D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMVector3D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(SEMStorage % interior_dev, &
!                             c_loc(SEMStorage % interior), &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SEMStorage % boundary_dev, &
!                             c_loc(SEMStorage % boundary), &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_SEMVector3D
!#endif
!
!FUNCTION BoundaryInterp_SEMVector3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMVector3D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMVector3D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorBoundaryInterp_3D( SEMStorage % interior_dev, &
!                                         SEMout % boundary_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % VectorBoundaryInterp_3D( SEMStorage % interior, &
!                                         SEMout % boundary, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_SEMVector3D
!
!FUNCTION GridInterp_SEMVector3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMVector3D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMVector3D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorGridInterp_3D( SEMStorage % interior_dev, &
!                                         SEMout % interior_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % VectorGridInterp_3D( SEMStorage % interior, &
!                                         SEMout % interior, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_SEMVector3D
!
!FUNCTION Gradient_SEMVector3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMVector3D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMTensor3D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorGradient_3D( SEMStorage % interior_dev, &
!                                       SEMout % interior_dev, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ELSE
!      CALL interp % VectorGradient_3D( SEMStorage % interior, &
!                                       SEMout % interior, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION Gradient_SEMVector3D
!
!FUNCTION Divergence_SEMVector3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMVector3D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMScalar3D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorDivergence_3D( SEMStorage % interior_dev, &
!                                       SEMout % interior_dev, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ELSE
!      CALL interp % VectorDivergence_3D( SEMStorage % interior, &
!                                       SEMout % interior, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION Divergence_SEMVector3D
!
!FUNCTION Curl_SEMVector3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMVector3D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMVector3D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorCurl_3D( SEMStorage % interior_dev, &
!                                       SEMout % interior_dev, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ELSE
!      CALL interp % VectorCurl_3D( SEMStorage % interior, &
!                                       SEMout % interior, &
!                                       SEMStorage % nVar, &
!                                       SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION Curl_SEMVector3D
!
!SUBROUTINE Equals_SEMVector3D( SEMOut, SEMStorage )
!  IMPLICIT NONE
!  TYPE(SEMVector3D), INTENT(out) :: SEMOut
!  TYPE(SEMVector3D), INTENT(in) :: SEMStorage
!
!    SEMOut % interior = SEMStorage % interior
!    SEMOut % boundary = SEMStorage % boundary
!  
!END SUBROUTINE Equals_SEMVector3D
!
!! -- SEMTensor2D -- !
!
!SUBROUTINE Build_SEMTensor2D( SEMStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(SEMTensor2D), INTENT(out) :: SEMStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SEMStorage % N = N
!    SEMStorage % nVar = nVar
!    SEMStorage % nElem = nElem
!
!    ALLOCATE( SEMStorage % interior(1:2,1:2,0:N,0:N,1:nVar,1:nElem), &
!              SEMStorage % boundary(1:2,1:2,0:N,1:nVar,1:4,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SEMStorage % interior_dev, SIZEOF(SEMStorage % interior)))
!    CALL hipCheck(hipMalloc(SEMStorage % boundary_dev, SIZEOF(SEMStorage % boundary)))
!#endif
!
!
!END SUBROUTINE Build_SEMTensor2D
!
!SUBROUTINE Trash_SEMTensor2D( SEMStorage ) 
!  IMPLICIT NONE
!  CLASS(SEMTensor2D), INTENT(inout) :: SEMStorage
!
!    DEALLOCATE( SEMStorage % interior, &
!                SEMStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SEMStorage % interior_dev))
!    CALL hipCheck(hipFree(SEMStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_SEMTensor2D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_SEMTensor2D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMTensor2D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % interior), &
!                             SEMStorage % interior_dev, &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % boundary), &
!                             SEMStorage % boundary_dev, &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_SEMTensor2D
!
!SUBROUTINE UpdateDevice_SEMTensor2D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMTensor2D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(SEMStorage % interior_dev, &
!                             c_loc(SEMStorage % interior), &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SEMStorage % boundary_dev, &
!                             c_loc(SEMStorage % boundary), &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_SEMTensor2D
!#endif
!
!FUNCTION BoundaryInterp_SEMTensor2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMTensor2D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMTensor2D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % TensorBoundaryInterp_2D( SEMStorage % interior_dev, &
!                                         SEMout % boundary_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % TensorBoundaryInterp_2D( SEMStorage % interior, &
!                                         SEMout % boundary, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_SEMTensor2D
!
!FUNCTION GridInterp_SEMTensor2D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMTensor2D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMTensor2D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % TensorGridInterp_2D( SEMStorage % interior_dev, &
!                                         SEMout % interior_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % TensorGridInterp_2D( SEMStorage % interior, &
!                                         SEMout % interior, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_SEMTensor2D
!
!FUNCTION Determinant_SEMTensor2D( SEMStorage ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMTensor2D) :: SEMStorage
!  TYPE(SEMScalar2D) :: SEMOut
!  ! Local
!  INTEGER :: iEl, iVar, i, j
!
!    DO iEl = 1, SEMStorage % nElem
!      DO iVar = 1, SEMStorage % nVar
!          DO j = 0, SEMStorage % N
!            DO i = 0, SEMStorage % N
!
!              SEMOut % interior(i,j,iVar,iEl) = SEMStorage % interior(1,1,i,j,iVar,iEl)*&
!                                                SEMStorage % interior(2,2,i,j,iVar,iEl)-&
!                                                SEMStorage % interior(1,2,i,j,iVar,iEl)*&
!                                                SEMStorage % interior(2,1,i,j,iVar,iEl)
!                                                     
!
!            ENDDO
!          ENDDO
!        ENDDO
!    ENDDO
!    
!END FUNCTION Determinant_SEMTensor2D
!
!SUBROUTINE Equals_SEMTensor2D( SEMOut, SEMStorage )
!  IMPLICIT NONE
!  TYPE(SEMTensor2D), INTENT(out) :: SEMOut
!  TYPE(SEMTensor2D), INTENT(in) :: SEMStorage
!
!    SEMOut % interior = SEMStorage % interior
!    SEMOut % boundary = SEMStorage % boundary
!  
!END SUBROUTINE Equals_SEMTensor2D
!
!! -- SEMTensor2D -- !
!
!SUBROUTINE Build_SEMTensor3D( SEMStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(SEMTensor3D), INTENT(out) :: SEMStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SEMStorage % N = N
!    SEMStorage % nVar = nVar
!    SEMStorage % nElem = nElem
!
!    ALLOCATE( SEMStorage % interior(1:3,1:3,0:N,0:N,0:N,1:nVar,1:nElem), &
!              SEMStorage % boundary(1:3,1:3,0:N,0:N,1:nVar,1:6,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SEMStorage % interior_dev, SIZEOF(SEMStorage % interior)))
!    CALL hipCheck(hipMalloc(SEMStorage % boundary_dev, SIZEOF(SEMStorage % boundary)))
!#endif
!
!
!END SUBROUTINE Build_SEMTensor3D
!
!SUBROUTINE Trash_SEMTensor3D( SEMStorage ) 
!  IMPLICIT NONE
!  CLASS(SEMTensor3D), INTENT(inout) :: SEMStorage
!
!    DEALLOCATE( SEMStorage % interior, &
!                SEMStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SEMStorage % interior_dev))
!    CALL hipCheck(hipFree(SEMStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_SEMTensor3D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_SEMTensor3D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMTensor3D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % interior), &
!                             SEMStorage % interior_dev, &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SEMStorage % boundary), &
!                             SEMStorage % boundary_dev, &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_SEMTensor3D
!
!SUBROUTINE UpdateDevice_SEMTensor3D( SEMStorage )
!  IMPLICIT NONE
!  CLASS(SEMTensor3D), INTENT(inout) :: SEMStorage
!
!     CALL hipCheck(hipMemcpy(SEMStorage % interior_dev, &
!                             c_loc(SEMStorage % interior), &
!                             SIZEOF(SEMStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SEMStorage % boundary_dev, &
!                             c_loc(SEMStorage % boundary), &
!                             SIZEOF(SEMStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_SEMTensor3D
!#endif
!
!FUNCTION BoundaryInterp_SEMTensor3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMTensor3D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMTensor3D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % TensorBoundaryInterp_3D( SEMStorage % interior_dev, &
!                                         SEMout % boundary_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % TensorBoundaryInterp_3D( SEMStorage % interior, &
!                                         SEMout % boundary, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_SEMTensor3D
!
!FUNCTION GridInterp_SEMTensor3D( SEMStorage, interp, gpuAccel ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMTensor3D) :: SEMStorage
!  TYPE(Lagrange) :: interp
!  TYPE(SEMTensor3D) :: SEMOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % TensorGridInterp_3D( SEMStorage % interior_dev, &
!                                         SEMout % interior_dev, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ELSE
!      CALL interp % TensorGridInterp_3D( SEMStorage % interior, &
!                                         SEMout % interior, &
!                                         SEMStorage % nVar, &
!                                         SEMStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_SEMTensor3D
!
!FUNCTION Determinant_SEMTensor3D( SEMStorage ) RESULT( SEMout )
!  IMPLICIT NONE
!  CLASS(SEMTensor3D) :: SEMStorage
!  TYPE(SEMScalar3D) :: SEMOut
!  ! Local
!  INTEGER :: iEl, iVar, i, j, k
!
!    DO iEl = 1, SEMStorage % nElem
!      DO iVar = 1, SEMStorage % nVar
!        DO k = 0, SEMStorage % N
!          DO j = 0, SEMStorage % N
!            DO i = 0, SEMStorage % N
!
!              SEMOut % interior(i,j,k,iVar,iEl) = SEMStorage % interior(1,1,i,j,k,iVar,iEl)*&
!                                                  ( SEMStorage % interior(2,2,i,j,k,iVar,iEl)*&
!                                                    SEMStorage % interior(3,3,i,j,k,iVar,iEl)-&
!                                                    SEMStorage % interior(2,3,i,j,k,iVar,iEl)*&
!                                                    SEMStorage % interior(3,2,i,j,k,iVar,iEl) )-&
!                                                 SEMStorage % interior(2,1,i,j,k,iVar,iEl)*&
!                                                  ( SEMStorage % interior(1,2,i,j,k,iVar,iEl)*&
!                                                    SEMStorage % interior(3,3,i,j,k,iVar,iEl)-&
!                                                    SEMStorage % interior(1,3,i,j,k,iVar,iEl)*&
!                                                    SEMStorage % interior(3,2,i,j,k,iVar,iEl) )+&
!                                                 SEMStorage % interior(3,1,i,j,k,iVar,iEl)*&
!                                                  ( SEMStorage % interior(1,2,i,j,k,iVar,iEl)*&
!                                                    SEMStorage % interior(2,3,i,j,k,iVar,iEl)-&
!                                                    SEMStorage % interior(1,3,i,j,k,iVar,iEl)*&
!                                                    SEMStorage % interior(2,2,i,j,k,iVar,iEl) )
!                                                     
!
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO
!    ENDDO
!    
!END FUNCTION Determinant_SEMTensor3D
!
!SUBROUTINE Equals_SEMTensor3D( SEMOut, SEMStorage )
!  IMPLICIT NONE
!  TYPE(SEMTensor3D), INTENT(out) :: SEMOut
!  TYPE(SEMTensor3D), INTENT(in) :: SEMStorage
!
!    SEMOut % interior = SEMStorage % interior
!    SEMOut % boundary = SEMStorage % boundary
!  
!END SUBROUTINE Equals_SEMTensor3D
!
!END MODULE NodalSEMData
