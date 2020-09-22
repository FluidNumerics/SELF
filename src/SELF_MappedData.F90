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

  TYPE, EXTENDS(Scalar2D), PUBLIC :: MappedScalar2D

    CONTAINS

      GENERIC, PUBLIC :: Gradient => Gradient_MappedScalar2D
      PROCEDURE, PRIVATE :: Gradient_MappedScalar2D
      PROCEDURE, PRIVATE :: ContravariantWeight => ContravariantWeight_MappedScalar2D_cpu

  END TYPE MappedScalar2D

  TYPE, EXTENDS(Scalar3D), PUBLIC :: MappedScalar3D

    CONTAINS

      GENERIC, PUBLIC :: Gradient => Gradient_MappedScalar3D
      PROCEDURE, PRIVATE :: Gradient_MappedScalar3D
      PROCEDURE, PRIVATE :: ContravariantWeight => ContravariantWeight_MappedScalar3D_cpu

  END TYPE MappedScalar3D

  TYPE, EXTENDS(Vector2D), PUBLIC :: MappedVector2D

    CONTAINS

      PROCEDURE, PUBLIC :: ContravariantProjection => ContravariantProjection_MappedVector2D
!      PROCEDURE, PUBLIC :: Divergence => Divergence_MappedVector2D
!      PROCEDURE, PUBLIC :: Curl => Curl_MappedVector2D

  END TYPE MappedVector2D

  TYPE, EXTENDS(Vector3D), PUBLIC :: MappedVector3D

    CONTAINS

      PROCEDURE, PUBLIC :: ContravariantProjection => ContravariantProjection_MappedVector3D
!      PROCEDURE, PUBLIC :: Divergence => Divergence_MappedVector3D
!      PROCEDURE, PUBLIC :: Curl => Curl_MappedVector3D

  END TYPE MappedVector3D

  TYPE, EXTENDS(Tensor2D), PUBLIC :: MappedTensor2D

    CONTAINS

      PROCEDURE, PUBLIC :: CompDivergence => CompDivergence_MappedTensor2D

  END TYPE MappedTensor2D

  TYPE, EXTENDS(Tensor3D), PUBLIC :: MappedTensor3D

    CONTAINS

      PROCEDURE, PUBLIC :: CompDivergence => CompDivergence_MappedTensor3D

  END TYPE MappedTensor3D

CONTAINS


! ---------------------- Scalars ---------------------- !

FUNCTION Gradient_MappedScalar2D( SELFStorage, workTensor, mesh ) RESULT( SELFOut )
  ! Strong Form Operator
  !
  ! Calculates the gradient of a scalar 2D function using the conservative form of the
  ! mapped gradient operator 
  !
  ! \grad_{phys}( f ) =  (1 / J)*(\partial / \partial \xi_i ( J\vec{a}_i f )
  ! 
  ! where the sum over i is implied.
  IMPLICIT NONE
  CLASS(MappedScalar2D) :: SELFStorage
  TYPE(MappedTensor2D) :: workTensor
  TYPE(Mesh2D) :: mesh
  TYPE(MappedVector2D) :: SELFOut

    CALL SELFStorage % ContravariantWeight( workTensor, mesh ) 
    SELFOut = workTensor % CompDivergence( ) 

END FUNCTION Gradient_MappedScalar2D

SUBROUTINE ContravariantWeight_MappedScalar2D_cpu( SELFStorage, workTensor, mesh)
  IMPLICIT NONE
  CLASS(MappedScalar2D), INTENT(in) :: SELFStorage
  TYPE(MappedTensor2D), INTENT(inout) :: workTensor
  TYPE(Mesh2D), INTENT(in) :: mesh
  ! Local
  INTEGER :: i, j, iVar, iEl

    DO iEl = 1, SELFStorage % nElem
      DO iVar = 1, SELFStorage % nVar
        DO j = 0, SELFStorage % N
          DO i = 0, SELFStorage % N

             workTensor % interior % hostData(1,1,i,j,iVar,iEl) = mesh % geometry % dsdx % interior % &
                                                                      hostData(1,1,i,j,iVar,iEl)*SELFStorage % interior % hostData(i,j,iVar,iEl) 

             workTensor % interior % hostData(2,1,i,j,iVar,iEl) = mesh % geometry % dsdx % interior % & 
                                                                      hostData(1,2,i,j,iVar,iEl)*SELFStorage % interior % hostData(i,j,iVar,iEl) 

             workTensor % interior % hostData(1,2,i,j,iVar,iEl) = mesh % geometry % dsdx % interior % &
                                                                      hostData(2,1,i,j,iVar,iEl)*SELFStorage % interior % hostData(i,j,iVar,iEl) 

             workTensor % interior % hostData(2,2,i,j,iVar,iEl) = mesh % geometry % dsdx % interior % &
                                                                      hostData(2,2,i,j,iVar,iEl)*SELFStorage % interior % hostData(i,j,iVar,iEl) 

          ENDDO
        ENDDO
      ENDDO
    ENDDO

END SUBROUTINE ContravariantWeight_MappedScalar2D_cpu

FUNCTION Gradient_MappedScalar3D( SELFStorage, workTensor, mesh ) RESULT( SELFOut )
  ! Strong Form Operator
  !
  ! Calculates the gradient of a scalar 3D function using the conservative form of the
  ! mapped gradient operator 
  !
  ! \grad_{phys}( f ) =  (1 / J)*(\partial / \partial \xi_i ( J\vec{a}_i f )
  ! 
  ! where the sum over i is implied.
  IMPLICIT NONE
  CLASS(MappedScalar3D) :: SELFStorage
  TYPE(MappedTensor3D) :: workTensor
  TYPE(Mesh3D) :: mesh
  TYPE(MappedVector3D) :: SELFOut

    CALL SELFStorage % ContravariantWeight( workTensor, mesh ) 
    SELFOut = workTensor % CompDivergence( ) 

END FUNCTION Gradient_MappedScalar3D

SUBROUTINE ContravariantWeight_MappedScalar3D_cpu( SELFStorage, workTensor, mesh)
  IMPLICIT NONE
  CLASS(MappedScalar3D), INTENT(in) :: SELFStorage
  TYPE(MappedTensor3D), INTENT(inout) :: workTensor
  TYPE(Mesh3D), INTENT(in) :: mesh
  ! Local
  INTEGER :: i, j, k, iVar, iEl

    DO iEl = 1, SELFStorage % nElem
      DO iVar = 1, SELFStorage % nVar
        DO k = 0, SELFStorage % N
          DO j = 0, SELFStorage % N
            DO i = 0, SELFStorage % N

               ! Get the x-component of the Jacobian weighted contravariant basis vectors multipled by the scalar
               workTensor % interior % hostData(1,1,i,j,k,iVar,iEl) = mesh % geometry % dsdx % interior % &
                                                                        hostData(1,1,i,j,k,iVar,iEl)*SELFStorage % interior % hostData(i,j,k,iVar,iEl) 

               workTensor % interior % hostData(2,1,i,j,k,iVar,iEl) = mesh % geometry % dsdx % interior % & 
                                                                        hostData(1,2,i,j,k,iVar,iEl)*SELFStorage % interior % hostData(i,j,k,iVar,iEl) 

               workTensor % interior % hostData(3,1,i,j,k,iVar,iEl) = mesh % geometry % dsdx % interior % & 
                                                                        hostData(1,3,i,j,k,iVar,iEl)*SELFStorage % interior % hostData(i,j,k,iVar,iEl) 

               ! Get the y-component of the Jacobian weighted contravariant basis vectors multipled by the scalar
               workTensor % interior % hostData(1,2,i,j,k,iVar,iEl) = mesh % geometry % dsdx % interior % &
                                                                        hostData(2,1,i,j,k,iVar,iEl)*SELFStorage % interior % hostData(i,j,k,iVar,iEl) 

               workTensor % interior % hostData(2,2,i,j,k,iVar,iEl) = mesh % geometry % dsdx % interior % &
                                                                        hostData(2,2,i,j,k,iVar,iEl)*SELFStorage % interior % hostData(i,j,k,iVar,iEl) 

               workTensor % interior % hostData(3,2,i,j,k,iVar,iEl) = mesh % geometry % dsdx % interior % &
                                                                        hostData(2,3,i,j,k,iVar,iEl)*SELFStorage % interior % hostData(i,j,k,iVar,iEl) 

               ! Get the z-component of the Jacobian weighted contravariant basis vectors multipled by the scalar
               workTensor % interior % hostData(1,3,i,j,k,iVar,iEl) = mesh % geometry % dsdx % interior % &
                                                                        hostData(3,1,i,j,k,iVar,iEl)*SELFStorage % interior % hostData(i,j,k,iVar,iEl) 

               workTensor % interior % hostData(2,3,i,j,k,iVar,iEl) = mesh % geometry % dsdx % interior % &
                                                                        hostData(3,2,i,j,k,iVar,iEl)*SELFStorage % interior % hostData(i,j,k,iVar,iEl) 

               workTensor % interior % hostData(3,3,i,j,k,iVar,iEl) = mesh % geometry % dsdx % interior % &
                                                                        hostData(3,3,i,j,k,iVar,iEl)*SELFStorage % interior % hostData(i,j,k,iVar,iEl) 

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

END SUBROUTINE ContravariantWeight_MappedScalar3D_cpu

! ---------------------- Vectors ---------------------- !

SUBROUTINE ContravariantProjection_MappedVector2D( physVector, compVector, mesh )
  ! Takes a vector that has physical space coordinate directions (x,y,z) and projects the vector
  ! into the the contravariant basis vector directions. Keep in mind that the contravariant basis
  ! vectors are really the Jacobian weighted contravariant basis vectors
  IMPLICIT NONE
  CLASS(MappedVector2D), INTENT(in) :: physVector
  TYPE(Vector2D), INTENT(inout) :: compVector
  TYPE(Mesh2D), INTENT(in) :: mesh
  ! Local
  INTEGER :: i, j, iVar, iEl

  ! Assume that tensor(j,i) is vector i, component j => dot product is done along first dimension to project onto computational
    DO iEl = 1, physVector % nElem
      DO iVar = 1, physVector % nVar
        DO j = 0, physVector % N
          DO i = 0, physVector % N
            
             compVector % interior % hostData(1,i,j,iVar,iEl) = mesh % geometry % dsdx % interior % hostData(1,1,i,j,1,iEl)*&
                                                                physVector % interior % hostData(1,i,j,iVar,iEl)+&
                                                                mesh % geometry % dsdx % interior % hostData(2,1,i,j,1,iEl)*&
                                                                physVector % interior % hostData(2,i,j,iVar,iEl)

             compVector % interior % hostData(2,i,j,iVar,iEl) = mesh % geometry % dsdx % interior % hostData(1,2,i,j,1,iEl)*&
                                                                physVector % interior % hostData(1,i,j,iVar,iEl)+&
                                                                mesh % geometry % dsdx % interior % hostData(2,2,i,j,1,iEl)*&
                                                                physVector % interior % hostData(2,i,j,iVar,iEl)

          ENDDO
        ENDDO
      ENDDO
    ENDDO

END SUBROUTINE ContravariantProjection_MappedVector2D

SUBROUTINE ContravariantProjection_MappedVector3D( physVector, compVector, mesh )
  ! Takes a vector that has physical space coordinate directions (x,y,z) and projects the vector
  ! into the the contravariant basis vector directions. Keep in mind that the contravariant basis
  ! vectors are really the Jacobian weighted contravariant basis vectors
  IMPLICIT NONE
  CLASS(MappedVector3D), INTENT(in) :: physVector
  TYPE(MappedVector3D), INTENT(inout) :: compVector
  TYPE(Mesh3D), INTENT(in) :: mesh
  ! Local
  INTEGER :: i, j, k, iVar, iEl, iDir

  ! Assume that tensor(j,i) is vector i, component j => dot product is done along first dimension to project onto computational
  ! space
    DO iEl = 1, physVector % nElem
      DO iVar = 1, physVector % nVar
        DO k = 0, physVector % N
          DO j = 0, physVector % N
            DO i = 0, physVector % N

               compVector % interior % hostData(1,i,j,k,iVar,iEl) = mesh % geometry % dsdx % interior % hostData(1,1,i,j,k,1,iEl)*&
                                                                    physVector % interior % hostData(1,i,j,k,iVar,iEl)+&
                                                                    mesh % geometry % dsdx % interior % hostData(2,1,i,j,k,1,iEl)*&
                                                                    physVector % interior % hostData(2,i,j,k,iVar,iEl)+&
                                                                    mesh % geometry % dsdx % interior % hostData(3,1,i,j,k,1,iEl)*&
                                                                    physVector % interior % hostData(3,i,j,k,iVar,iEl)

               compVector % interior % hostData(2,i,j,k,iVar,iEl) = mesh % geometry % dsdx % interior % hostData(1,2,i,j,k,1,iEl)*&
                                                                    physVector % interior % hostData(1,i,j,k,iVar,iEl)+&
                                                                    mesh % geometry % dsdx % interior % hostData(2,2,i,j,k,1,iEl)*&
                                                                    physVector % interior % hostData(2,i,j,k,iVar,iEl)+&
                                                                    mesh % geometry % dsdx % interior % hostData(3,2,i,j,k,1,iEl)*&
                                                                    physVector % interior % hostData(3,i,j,k,iVar,iEl)

               compVector % interior % hostData(3,i,j,k,iVar,iEl) = mesh % geometry % dsdx % interior % hostData(1,3,i,j,k,1,iEl)*&
                                                                    physVector % interior % hostData(1,i,j,k,iVar,iEl)+&
                                                                    mesh % geometry % dsdx % interior % hostData(2,3,i,j,k,1,iEl)*&
                                                                    physVector % interior % hostData(2,i,j,k,iVar,iEl)+&
                                                                    mesh % geometry % dsdx % interior % hostData(3,3,i,j,k,1,iEl)*&
                                                                    physVector % interior % hostData(3,i,j,k,iVar,iEl)

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

END SUBROUTINE ContravariantProjection_MappedVector3D

! ---------------------- Tensors ---------------------- !

FUNCTION CompDivergence_MappedTensor2D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(MappedTensor2D) :: SELFStorage
  TYPE(MappedVector2D) :: SELFOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
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

END FUNCTION CompDivergence_MappedTensor2D

FUNCTION CompDivergence_MappedTensor3D( SELFStorage, gpuAccel ) RESULT( SELFout )
  IMPLICIT NONE
  CLASS(MappedTensor3D) :: SELFStorage
  TYPE(MappedVector3D) :: SELFOut
  LOGICAL, OPTIONAL :: gpuAccel

    IF( PRESENT(gpuAccel) )THEN
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

END FUNCTION CompDivergence_MappedTensor3D

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
END MODULE SELF_MappedData
