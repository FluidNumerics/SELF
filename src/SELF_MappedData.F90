! SELF_MappedData.F90
! 
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self-fluids@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE SELF_MappedData

USE SELF_Constants
USE SELF_Lagrange

USE SELF_Data
USE SELF_Mesh

USE hipfort
USE ISO_C_BINDING


IMPLICIT NONE

!! ---------------------- Scalars ---------------------- !
!  TYPE, EXTENDS(Scalar1D), PUBLIC :: MappedScalar1D
!
!    CONTAINS
!
!      PROCEDURE, PUBLIC :: MappedDerivative => MappedDerivative_Scalar1D
!
!  END TYPE Scalar1D
!
  TYPE, EXTENDS(Scalar2D), PUBLIC :: MappedScalar2D

    CONTAINS

      PROCEDURE, PUBLIC :: MappedGradient => MappedGradient_MappedScalar2D

  END TYPE MappedScalar2D

!  TYPE, EXTENDS(Scalar3D), PUBLIC :: MappedScalar3D
!
!    CONTAINS
!
!      PROCEDURE, PUBLIC :: Gradient => MappedGradient_MappedScalar3D
!
!  END TYPE MappedScalar3D
!
!  TYPE, EXTENDS(Vector2D), PUBLIC :: MappedVector2D
!
!    CONTAINS
!
!      PROCEDURE, PUBLIC :: ContravariantProjection => ContravariantProjection_MappedVector2D
!      PROCEDURE, PUBLIC :: CovariantProjection => CovariantProjection_MappedVector2D
!      PROCEDURE, PUBLIC :: Gradient => Gradient_MappedVector2D
!      PROCEDURE, PUBLIC :: Divergence => Divergence_MappedVector2D
!      PROCEDURE, PUBLIC :: Curl => Curl_MappedVector2D
!
!  END TYPE MappedVector2D
!
!  TYPE, EXTENDS(Vector3D), PUBLIC :: MappedVector3D
!
!    CONTAINS
!
!      PROCEDURE, PUBLIC :: ContravariantProjection => ContravariantProjection_MappedVector3D
!      PROCEDURE, PUBLIC :: CovariantProjection => CovariantProjection_MappedVector3D
!      PROCEDURE, PUBLIC :: Gradient => Gradient_MappedVector3D
!      PROCEDURE, PUBLIC :: Divergence => Divergence_MappedVector3D
!      PROCEDURE, PUBLIC :: Curl => Curl_MappedVector3D
!
!  END TYPE MappedVector3D
!
!  TYPE, EXTENDS(Tensor2D), PUBLIC :: MappedTensor2D
!
!    CONTAINS
!
!  END TYPE MappedTensor2D
!
!  TYPE, EXTENDS(Tensor3D), PUBLIC :: MappedTensor3D
!
!    CONTAINS
!
!  END TYPE MappedTensor3D

CONTAINS

FUNCTION MappedGradient_Scalar2D( SELFStorage, mesh ) RESULT( SELFOut )
  IMPLICIT NONE
  CLASS(MappedScalar2D) :: SELFStorage
  TYPE(Mesh2D) :: mesh
  TYPE(MappedVector2D) :: SELFOut
  
    SELFOut = SELFStorage % Gradient( 
  
END FUNCTION MappedGradient_Scalar2D
!SUBROUTINE Build_Scalar1D( SELFStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(Scalar1D), INTENT(out) :: SELFStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SELFStorage % N = N
!    SELFStorage % nVar = nVar
!    SELFStorage % nElem = nElem
!
!    ALLOCATE( SELFStorage % interior(0:N,1:nVar,1:nElem), &
!              SELFStorage % boundary(1:nVar,1:2,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SELFStorage % interior_dev, SIZEOF(SELFStorage % interior)))
!    CALL hipCheck(hipMalloc(SELFStorage % boundary_dev, SIZEOF(SELFStorage % boundary)))
!#endif
!
!
!END SUBROUTINE Build_Scalar1D
!
!SUBROUTINE Trash_Scalar1D( SELFStorage ) 
!  IMPLICIT NONE
!  CLASS(Scalar1D), INTENT(inout) :: SELFStorage
!
!    DEALLOCATE( SELFStorage % interior, &
!                SELFStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SELFStorage % interior_dev))
!    CALL hipCheck(hipFree(SELFStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_Scalar1D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_Scalar1D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Scalar1D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % interior), &
!                             SELFStorage % interior_dev, &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % boundary), &
!                             SELFStorage % boundary_dev, &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_Scalar1D
!
!SUBROUTINE UpdateDevice_Scalar1D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Scalar1D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(SELFStorage % interior_dev, &
!                             c_loc(SELFStorage % interior), &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SELFStorage % boundary_dev, &
!                             c_loc(SELFStorage % boundary), &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_Scalar1D
!#endif
!
!FUNCTION BoundaryInterp_Scalar1D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Scalar1D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Scalar1D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarBoundaryInterp_1D( SELFStorage % interior_dev, &
!                                         SELFout % boundary_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % ScalarBoundaryInterp_1D( SELFStorage % interior, &
!                                         SELFout % boundary, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_Scalar1D
!
!FUNCTION GridInterp_Scalar1D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Scalar1D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Scalar1D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarGridInterp_1D( SELFStorage % interior_dev, &
!                                         SELFout % interior_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % ScalarGridInterp_1D( SELFStorage % interior, &
!                                         SELFout % interior, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_Scalar1D
!
!FUNCTION Derivative_Scalar1D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Scalar1D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Scalar1D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % Derivative_1D( SELFStorage % interior_dev, &
!                                   SELFout % interior_dev, &
!                                   SELFStorage % nVar, &
!                                   SELFStorage % nElem )  
!    ELSE
!      CALL interp % Derivative_1D( SELFStorage % interior, &
!                                   SELFout % interior, &
!                                   SELFStorage % nVar, &
!                                   SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION Derivative_Scalar1D
!
!SUBROUTINE Equals_Scalar1D( SELFOut, SELFStorage )
!  IMPLICIT NONE
!  TYPE(Scalar1D), INTENT(out) :: SELFOut
!  TYPE(Scalar1D), INTENT(in) :: SELFStorage
!
!    SELFOut % interior = SELFStorage % interior
!    SELFOut % boundary = SELFStorage % boundary
!  
!END SUBROUTINE Equals_Scalar1D
!
!! -- Scalar2D -- !
!
!SUBROUTINE Build_Scalar2D( SELFStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(Scalar2D), INTENT(out) :: SELFStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SELFStorage % N = N
!    SELFStorage % nVar = nVar
!    SELFStorage % nElem = nElem
!
!    ALLOCATE( SELFStorage % interior(0:N,0:N,1:nVar,1:nElem), &
!              SELFStorage % boundary(0:N,1:nVar,1:4,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SELFStorage % interior_dev, SIZEOF(SELFStorage % interior)))
!    CALL hipCheck(hipMalloc(SELFStorage % boundary_dev, SIZEOF(SELFStorage % boundary)))
!#endif
!
!
!END SUBROUTINE Build_Scalar2D
!
!SUBROUTINE Trash_Scalar2D( SELFStorage ) 
!  IMPLICIT NONE
!  CLASS(Scalar2D), INTENT(inout) :: SELFStorage
!
!    DEALLOCATE( SELFStorage % interior, &
!                SELFStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SELFStorage % interior_dev))
!    CALL hipCheck(hipFree(SELFStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_Scalar2D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_Scalar2D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Scalar2D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % interior), &
!                             SELFStorage % interior_dev, &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % boundary), &
!                             SELFStorage % boundary_dev, &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_Scalar2D
!
!SUBROUTINE UpdateDevice_Scalar2D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Scalar2D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(SELFStorage % interior_dev, &
!                             c_loc(SELFStorage % interior), &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SELFStorage % boundary_dev, &
!                             c_loc(SELFStorage % boundary), &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_Scalar2D
!#endif
!
!FUNCTION BoundaryInterp_Scalar2D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Scalar2D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Scalar2D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarBoundaryInterp_2D( SELFStorage % interior_dev, &
!                                         SELFout % boundary_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % ScalarBoundaryInterp_2D( SELFStorage % interior, &
!                                         SELFout % boundary, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_Scalar2D
!
!FUNCTION GridInterp_Scalar2D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Scalar2D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Scalar2D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarGridInterp_2D( SELFStorage % interior_dev, &
!                                         SELFout % interior_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % ScalarGridInterp_2D( SELFStorage % interior, &
!                                         SELFout % interior, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_Scalar2D
!
!FUNCTION Gradient_Scalar2D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Scalar2D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Vector2D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarGradient_2D( SELFStorage % interior_dev, &
!                                       SELFout % interior_dev, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ELSE
!      CALL interp % ScalarGradient_2D( SELFStorage % interior, &
!                                       SELFout % interior, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION Gradient_Scalar2D
!
!SUBROUTINE Equals_Scalar2D( SELFOut, SELFStorage )
!  IMPLICIT NONE
!  TYPE(Scalar2D), INTENT(out) :: SELFOut
!  TYPE(Scalar2D), INTENT(in) :: SELFStorage
!
!    SELFOut % interior = SELFStorage % interior
!    SELFOut % boundary = SELFStorage % boundary
!  
!END SUBROUTINE Equals_Scalar2D
!
!! -- Scalar3D -- !
!
!SUBROUTINE Build_Scalar3D( SELFStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(Scalar3D), INTENT(out) :: SELFStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SELFStorage % N = N
!    SELFStorage % nVar = nVar
!    SELFStorage % nElem = nElem
!
!    ALLOCATE( SELFStorage % interior(0:N,0:N,0:N,1:nVar,1:nElem), &
!              SELFStorage % boundary(0:N,0:N,1:nVar,1:6,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SELFStorage % interior_dev, SIZEOF(SELFStorage % interior)))
!    CALL hipCheck(hipMalloc(SELFStorage % boundary_dev, SIZEOF(SELFStorage % boundary)))
!#endif
!
!END SUBROUTINE Build_Scalar3D
!
!SUBROUTINE Trash_Scalar3D( SELFStorage ) 
!  IMPLICIT NONE
!  CLASS(Scalar3D), INTENT(inout) :: SELFStorage
!
!    DEALLOCATE( SELFStorage % interior, &
!                SELFStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SELFStorage % interior_dev))
!    CALL hipCheck(hipFree(SELFStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_Scalar3D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_Scalar3D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Scalar3D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % interior), &
!                             SELFStorage % interior_dev, &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % boundary), &
!                             SELFStorage % boundary_dev, &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_Scalar3D
!
!SUBROUTINE UpdateDevice_Scalar3D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Scalar3D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(SELFStorage % interior_dev, &
!                             c_loc(SELFStorage % interior), &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SELFStorage % boundary_dev, &
!                             c_loc(SELFStorage % boundary), &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_Scalar3D
!#endif
!
!FUNCTION BoundaryInterp_Scalar3D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Scalar3D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Scalar3D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarBoundaryInterp_3D( SELFStorage % interior_dev, &
!                                         SELFout % boundary_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % ScalarBoundaryInterp_3D( SELFStorage % interior, &
!                                         SELFout % boundary, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_Scalar3D
!
!FUNCTION GridInterp_Scalar3D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Scalar3D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Scalar3D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarGridInterp_3D( SELFStorage % interior_dev, &
!                                         SELFout % interior_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % ScalarGridInterp_3D( SELFStorage % interior, &
!                                         SELFout % interior, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_Scalar3D
!
!FUNCTION Gradient_Scalar3D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Scalar3D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Vector3D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % ScalarGradient_3D( SELFStorage % interior_dev, &
!                                       SELFout % interior_dev, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ELSE
!      CALL interp % ScalarGradient_3D( SELFStorage % interior, &
!                                       SELFout % interior, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION Gradient_Scalar3D
!
!SUBROUTINE Equals_Scalar3D( SELFOut, SELFStorage )
!  IMPLICIT NONE
!  TYPE(Scalar3D), INTENT(out) :: SELFOut
!  TYPE(Scalar3D), INTENT(in) :: SELFStorage
!
!    SELFOut % interior = SELFStorage % interior
!    SELFOut % boundary = SELFStorage % boundary
!  
!END SUBROUTINE Equals_Scalar3D
!
!! -- Vector2D -- !
!
!SUBROUTINE Build_Vector2D( SELFStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(Vector2D), INTENT(out) :: SELFStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SELFStorage % N = N
!    SELFStorage % nVar = nVar
!    SELFStorage % nElem = nElem
!
!    ALLOCATE( SELFStorage % interior(1:2,0:N,0:N,1:nVar,1:nElem), &
!              SELFStorage % boundary(1:2,0:N,1:nVar,1:4,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SELFStorage % interior_dev, SIZEOF(SELFStorage % interior)))
!    CALL hipCheck(hipMalloc(SELFStorage % boundary_dev, SIZEOF(SELFStorage % boundary)))
!#endif
!
!
!END SUBROUTINE Build_Vector2D
!
!SUBROUTINE Trash_Vector2D( SELFStorage ) 
!  IMPLICIT NONE
!  CLASS(Vector2D), INTENT(inout) :: SELFStorage
!
!    DEALLOCATE( SELFStorage % interior, &
!                SELFStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SELFStorage % interior_dev))
!    CALL hipCheck(hipFree(SELFStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_Vector2D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_Vector2D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Vector2D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % interior), &
!                             SELFStorage % interior_dev, &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % boundary), &
!                             SELFStorage % boundary_dev, &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_Vector2D
!
!SUBROUTINE UpdateDevice_Vector2D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Vector2D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(SELFStorage % interior_dev, &
!                             c_loc(SELFStorage % interior), &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SELFStorage % boundary_dev, &
!                             c_loc(SELFStorage % boundary), &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_Vector2D
!#endif
!
!FUNCTION BoundaryInterp_Vector2D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Vector2D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Vector2D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorBoundaryInterp_2D( SELFStorage % interior_dev, &
!                                         SELFout % boundary_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % VectorBoundaryInterp_2D( SELFStorage % interior, &
!                                         SELFout % boundary, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_Vector2D
!
!FUNCTION GridInterp_Vector2D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Vector2D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Vector2D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorGridInterp_2D( SELFStorage % interior_dev, &
!                                         SELFout % interior_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % VectorGridInterp_2D( SELFStorage % interior, &
!                                         SELFout % interior, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_Vector2D
!
!FUNCTION Gradient_Vector2D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Vector2D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Tensor2D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorGradient_2D( SELFStorage % interior_dev, &
!                                       SELFout % interior_dev, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ELSE
!      CALL interp % VectorGradient_2D( SELFStorage % interior, &
!                                       SELFout % interior, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION Gradient_Vector2D
!
!FUNCTION Divergence_Vector2D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Vector2D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Scalar2D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorDivergence_2D( SELFStorage % interior_dev, &
!                                       SELFout % interior_dev, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ELSE
!      CALL interp % VectorDivergence_2D( SELFStorage % interior, &
!                                       SELFout % interior, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION Divergence_Vector2D
!
!FUNCTION Curl_Vector2D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Vector2D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Scalar2D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorCurl_2D( SELFStorage % interior_dev, &
!                                       SELFout % interior_dev, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ELSE
!      CALL interp % VectorCurl_2D( SELFStorage % interior, &
!                                       SELFout % interior, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION Curl_Vector2D
!
!SUBROUTINE Equals_Vector2D( SELFOut, SELFStorage )
!  IMPLICIT NONE
!  TYPE(Vector2D), INTENT(out) :: SELFOut
!  TYPE(Vector2D), INTENT(in) :: SELFStorage
!
!    SELFOut % interior = SELFStorage % interior
!    SELFOut % boundary = SELFStorage % boundary
!  
!END SUBROUTINE Equals_Vector2D
!
!! -- Vector3D -- !
!
!SUBROUTINE Build_Vector3D( SELFStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(Vector3D), INTENT(out) :: SELFStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SELFStorage % N = N
!    SELFStorage % nVar = nVar
!    SELFStorage % nElem = nElem
!
!    ALLOCATE( SELFStorage % interior(1:3,0:N,0:N,0:N,1:nVar,1:nElem), &
!              SELFStorage % boundary(1:3,0:N,0:N,1:nVar,1:6,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SELFStorage % interior_dev, SIZEOF(SELFStorage % interior)))
!    CALL hipCheck(hipMalloc(SELFStorage % boundary_dev, SIZEOF(SELFStorage % boundary)))
!#endif
!
!
!END SUBROUTINE Build_Vector3D
!
!SUBROUTINE Trash_Vector3D( SELFStorage ) 
!  IMPLICIT NONE
!  CLASS(Vector3D), INTENT(inout) :: SELFStorage
!
!    DEALLOCATE( SELFStorage % interior, &
!                SELFStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SELFStorage % interior_dev))
!    CALL hipCheck(hipFree(SELFStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_Vector3D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_Vector3D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Vector3D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % interior), &
!                             SELFStorage % interior_dev, &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % boundary), &
!                             SELFStorage % boundary_dev, &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_Vector3D
!
!SUBROUTINE UpdateDevice_Vector3D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Vector3D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(SELFStorage % interior_dev, &
!                             c_loc(SELFStorage % interior), &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SELFStorage % boundary_dev, &
!                             c_loc(SELFStorage % boundary), &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_Vector3D
!#endif
!
!FUNCTION BoundaryInterp_Vector3D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Vector3D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Vector3D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorBoundaryInterp_3D( SELFStorage % interior_dev, &
!                                         SELFout % boundary_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % VectorBoundaryInterp_3D( SELFStorage % interior, &
!                                         SELFout % boundary, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_Vector3D
!
!FUNCTION GridInterp_Vector3D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Vector3D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Vector3D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorGridInterp_3D( SELFStorage % interior_dev, &
!                                         SELFout % interior_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % VectorGridInterp_3D( SELFStorage % interior, &
!                                         SELFout % interior, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_Vector3D
!
!FUNCTION Gradient_Vector3D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Vector3D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Tensor3D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorGradient_3D( SELFStorage % interior_dev, &
!                                       SELFout % interior_dev, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ELSE
!      CALL interp % VectorGradient_3D( SELFStorage % interior, &
!                                       SELFout % interior, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION Gradient_Vector3D
!
!FUNCTION Divergence_Vector3D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Vector3D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Scalar3D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorDivergence_3D( SELFStorage % interior_dev, &
!                                       SELFout % interior_dev, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ELSE
!      CALL interp % VectorDivergence_3D( SELFStorage % interior, &
!                                       SELFout % interior, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION Divergence_Vector3D
!
!FUNCTION Curl_Vector3D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Vector3D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Vector3D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % VectorCurl_3D( SELFStorage % interior_dev, &
!                                       SELFout % interior_dev, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ELSE
!      CALL interp % VectorCurl_3D( SELFStorage % interior, &
!                                       SELFout % interior, &
!                                       SELFStorage % nVar, &
!                                       SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION Curl_Vector3D
!
!SUBROUTINE Equals_Vector3D( SELFOut, SELFStorage )
!  IMPLICIT NONE
!  TYPE(Vector3D), INTENT(out) :: SELFOut
!  TYPE(Vector3D), INTENT(in) :: SELFStorage
!
!    SELFOut % interior = SELFStorage % interior
!    SELFOut % boundary = SELFStorage % boundary
!  
!END SUBROUTINE Equals_Vector3D
!
!! -- Tensor2D -- !
!
!SUBROUTINE Build_Tensor2D( SELFStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(Tensor2D), INTENT(out) :: SELFStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SELFStorage % N = N
!    SELFStorage % nVar = nVar
!    SELFStorage % nElem = nElem
!
!    ALLOCATE( SELFStorage % interior(1:2,1:2,0:N,0:N,1:nVar,1:nElem), &
!              SELFStorage % boundary(1:2,1:2,0:N,1:nVar,1:4,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SELFStorage % interior_dev, SIZEOF(SELFStorage % interior)))
!    CALL hipCheck(hipMalloc(SELFStorage % boundary_dev, SIZEOF(SELFStorage % boundary)))
!#endif
!
!
!END SUBROUTINE Build_Tensor2D
!
!SUBROUTINE Trash_Tensor2D( SELFStorage ) 
!  IMPLICIT NONE
!  CLASS(Tensor2D), INTENT(inout) :: SELFStorage
!
!    DEALLOCATE( SELFStorage % interior, &
!                SELFStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SELFStorage % interior_dev))
!    CALL hipCheck(hipFree(SELFStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_Tensor2D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_Tensor2D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Tensor2D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % interior), &
!                             SELFStorage % interior_dev, &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % boundary), &
!                             SELFStorage % boundary_dev, &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_Tensor2D
!
!SUBROUTINE UpdateDevice_Tensor2D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Tensor2D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(SELFStorage % interior_dev, &
!                             c_loc(SELFStorage % interior), &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SELFStorage % boundary_dev, &
!                             c_loc(SELFStorage % boundary), &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_Tensor2D
!#endif
!
!FUNCTION BoundaryInterp_Tensor2D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Tensor2D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Tensor2D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % TensorBoundaryInterp_2D( SELFStorage % interior_dev, &
!                                         SELFout % boundary_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % TensorBoundaryInterp_2D( SELFStorage % interior, &
!                                         SELFout % boundary, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_Tensor2D
!
!FUNCTION GridInterp_Tensor2D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Tensor2D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Tensor2D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % TensorGridInterp_2D( SELFStorage % interior_dev, &
!                                         SELFout % interior_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % TensorGridInterp_2D( SELFStorage % interior, &
!                                         SELFout % interior, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_Tensor2D
!
!FUNCTION Determinant_Tensor2D( SELFStorage ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Tensor2D) :: SELFStorage
!  TYPE(Scalar2D) :: SELFOut
!  ! Local
!  INTEGER :: iEl, iVar, i, j
!
!    DO iEl = 1, SELFStorage % nElem
!      DO iVar = 1, SELFStorage % nVar
!          DO j = 0, SELFStorage % N
!            DO i = 0, SELFStorage % N
!
!              SELFOut % interior(i,j,iVar,iEl) = SELFStorage % interior(1,1,i,j,iVar,iEl)*&
!                                                SELFStorage % interior(2,2,i,j,iVar,iEl)-&
!                                                SELFStorage % interior(1,2,i,j,iVar,iEl)*&
!                                                SELFStorage % interior(2,1,i,j,iVar,iEl)
!                                                     
!
!            ENDDO
!          ENDDO
!        ENDDO
!    ENDDO
!    
!END FUNCTION Determinant_Tensor2D
!
!SUBROUTINE Equals_Tensor2D( SELFOut, SELFStorage )
!  IMPLICIT NONE
!  TYPE(Tensor2D), INTENT(out) :: SELFOut
!  TYPE(Tensor2D), INTENT(in) :: SELFStorage
!
!    SELFOut % interior = SELFStorage % interior
!    SELFOut % boundary = SELFStorage % boundary
!  
!END SUBROUTINE Equals_Tensor2D
!
!! -- Tensor2D -- !
!
!SUBROUTINE Build_Tensor3D( SELFStorage, N, nVar, nElem ) 
!  IMPLICIT NONE
!  CLASS(Tensor3D), INTENT(out) :: SELFStorage
!  INTEGER, INTENT(in) :: N
!  INTEGER, INTENT(in) :: nVar
!  INTEGER, INTENT(in) :: nElem
!
!    SELFStorage % N = N
!    SELFStorage % nVar = nVar
!    SELFStorage % nElem = nElem
!
!    ALLOCATE( SELFStorage % interior(1:3,1:3,0:N,0:N,0:N,1:nVar,1:nElem), &
!              SELFStorage % boundary(1:3,1:3,0:N,0:N,1:nVar,1:6,1:nElem) )
!
!#ifdef GPU
!    CALL hipCheck(hipMalloc(SELFStorage % interior_dev, SIZEOF(SELFStorage % interior)))
!    CALL hipCheck(hipMalloc(SELFStorage % boundary_dev, SIZEOF(SELFStorage % boundary)))
!#endif
!
!
!END SUBROUTINE Build_Tensor3D
!
!SUBROUTINE Trash_Tensor3D( SELFStorage ) 
!  IMPLICIT NONE
!  CLASS(Tensor3D), INTENT(inout) :: SELFStorage
!
!    DEALLOCATE( SELFStorage % interior, &
!                SELFStorage % boundary )
!#ifdef GPU
!    CALL hipCheck(hipFree(SELFStorage % interior_dev))
!    CALL hipCheck(hipFree(SELFStorage % boundary_dev))
!#endif
!
!END SUBROUTINE Trash_Tensor3D
!
!#ifdef GPU
!SUBROUTINE UpdateHost_Tensor3D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Tensor3D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % interior), &
!                             SELFStorage % interior_dev, &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyDeviceToHost))
!
!     CALL hipCheck(hipMemcpy(c_loc(SELFStorage % boundary), &
!                             SELFStorage % boundary_dev, &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyDeviceToHost))
!
!END SUBROUTINE UpdateHost_Tensor3D
!
!SUBROUTINE UpdateDevice_Tensor3D( SELFStorage )
!  IMPLICIT NONE
!  CLASS(Tensor3D), INTENT(inout) :: SELFStorage
!
!     CALL hipCheck(hipMemcpy(SELFStorage % interior_dev, &
!                             c_loc(SELFStorage % interior), &
!                             SIZEOF(SELFStorage % interior), &
!                             hipMemcpyHostToDevice))
!
!     CALL hipCheck(hipMemcpy(SELFStorage % boundary_dev, &
!                             c_loc(SELFStorage % boundary), &
!                             SIZEOF(SELFStorage % boundary), &
!                             hipMemcpyHostToDevice))
!
!END SUBROUTINE UpdateDevice_Tensor3D
!#endif
!
!FUNCTION BoundaryInterp_Tensor3D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Tensor3D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Tensor3D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % TensorBoundaryInterp_3D( SELFStorage % interior_dev, &
!                                         SELFout % boundary_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % TensorBoundaryInterp_3D( SELFStorage % interior, &
!                                         SELFout % boundary, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION BoundaryInterp_Tensor3D
!
!FUNCTION GridInterp_Tensor3D( SELFStorage, interp, gpuAccel ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Tensor3D) :: SELFStorage
!  TYPE(Lagrange) :: interp
!  TYPE(Tensor3D) :: SELFOut
!  LOGICAL, OPTIONAL :: gpuAccel
!
!    IF( PRESENT(gpuAccel) )THEN
!      CALL interp % TensorGridInterp_3D( SELFStorage % interior_dev, &
!                                         SELFout % interior_dev, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ELSE
!      CALL interp % TensorGridInterp_3D( SELFStorage % interior, &
!                                         SELFout % interior, &
!                                         SELFStorage % nVar, &
!                                         SELFStorage % nElem )  
!    ENDIF
!
!END FUNCTION GridInterp_Tensor3D
!
!FUNCTION Determinant_Tensor3D( SELFStorage ) RESULT( SELFout )
!  IMPLICIT NONE
!  CLASS(Tensor3D) :: SELFStorage
!  TYPE(Scalar3D) :: SELFOut
!  ! Local
!  INTEGER :: iEl, iVar, i, j, k
!
!    DO iEl = 1, SELFStorage % nElem
!      DO iVar = 1, SELFStorage % nVar
!        DO k = 0, SELFStorage % N
!          DO j = 0, SELFStorage % N
!            DO i = 0, SELFStorage % N
!
!              SELFOut % interior(i,j,k,iVar,iEl) = SELFStorage % interior(1,1,i,j,k,iVar,iEl)*&
!                                                  ( SELFStorage % interior(2,2,i,j,k,iVar,iEl)*&
!                                                    SELFStorage % interior(3,3,i,j,k,iVar,iEl)-&
!                                                    SELFStorage % interior(2,3,i,j,k,iVar,iEl)*&
!                                                    SELFStorage % interior(3,2,i,j,k,iVar,iEl) )-&
!                                                 SELFStorage % interior(2,1,i,j,k,iVar,iEl)*&
!                                                  ( SELFStorage % interior(1,2,i,j,k,iVar,iEl)*&
!                                                    SELFStorage % interior(3,3,i,j,k,iVar,iEl)-&
!                                                    SELFStorage % interior(1,3,i,j,k,iVar,iEl)*&
!                                                    SELFStorage % interior(3,2,i,j,k,iVar,iEl) )+&
!                                                 SELFStorage % interior(3,1,i,j,k,iVar,iEl)*&
!                                                  ( SELFStorage % interior(1,2,i,j,k,iVar,iEl)*&
!                                                    SELFStorage % interior(2,3,i,j,k,iVar,iEl)-&
!                                                    SELFStorage % interior(1,3,i,j,k,iVar,iEl)*&
!                                                    SELFStorage % interior(2,2,i,j,k,iVar,iEl) )
!                                                     
!
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO
!    ENDDO
!    
!END FUNCTION Determinant_Tensor3D
!
!SUBROUTINE Equals_Tensor3D( SELFOut, SELFStorage )
!  IMPLICIT NONE
!  TYPE(Tensor3D), INTENT(out) :: SELFOut
!  TYPE(Tensor3D), INTENT(in) :: SELFStorage
!
!    SELFOut % interior = SELFStorage % interior
!    SELFOut % boundary = SELFStorage % boundary
!  
!END SUBROUTINE Equals_Tensor3D
!
!END MODULE NodalSELFData
