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

!USE ISO_C_BINDING


IMPLICIT NONE

  TYPE, EXTENDS(Scalar2D), PUBLIC :: MappedScalar2D

    CONTAINS

      GENERIC, PUBLIC :: Gradient => Gradient_MappedScalar2D
      PROCEDURE, PRIVATE :: Gradient_MappedScalar2D
      PROCEDURE, PRIVATE :: ContravariantWeight => ContravariantWeight_MappedScalar2D

  END TYPE MappedScalar2D

  TYPE, EXTENDS(Scalar3D), PUBLIC :: MappedScalar3D

    CONTAINS

      GENERIC, PUBLIC :: Gradient => Gradient_MappedScalar3D
      PROCEDURE, PRIVATE :: Gradient_MappedScalar3D
      PROCEDURE, PRIVATE :: ContravariantWeight => ContravariantWeight_MappedScalar3D

  END TYPE MappedScalar3D

  TYPE, EXTENDS(Vector2D), PUBLIC :: MappedVector2D

    CONTAINS

      GENERIC, PUBLIC :: ContravariantProjection => ContravariantProjection_MappedVector2D
      PROCEDURE, PRIVATE :: ContravariantProjection_MappedVector2D
!      PROCEDURE, PUBLIC :: Divergence => Divergence_MappedVector2D
!      PROCEDURE, PUBLIC :: Curl => Curl_MappedVector2D

  END TYPE MappedVector2D

  TYPE, EXTENDS(Vector3D), PUBLIC :: MappedVector3D

    CONTAINS

      GENERIC, PUBLIC :: ContravariantProjection => ContravariantProjection_MappedVector3D
      PROCEDURE, PRIVATE :: ContravariantProjection_MappedVector3D
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

  INTERFACE
    SUBROUTINE ContravariantWeight_MappedScalar2D_gpu_wrapper(scalar, workTensor, dsdx, N, nVar, nEl) bind(c,name="ContravariantWeight_MappedScalar2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar, workTensor, dsdx
      INTEGER, VALUE :: N, nVar, nEl
    END SUBROUTINE ContravariantWeight_MappedScalar2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantWeight_MappedScalar3D_gpu_wrapper(scalar, workTensor, dsdx, N, nVar, nEl) bind(c,name="ContravariantWeight_MappedScalar3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar, workTensor, dsdx
      INTEGER, VALUE :: N, nVar, nEl
    END SUBROUTINE ContravariantWeight_MappedScalar3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantProjection_MappedVector2D_gpu_wrapper(physVector, compVector, dsdx, N, nVar, nEl) bind(c,name="ContravariantProjection_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: physVector, compVector, dsdx
      INTEGER, VALUE :: N, nVar, nEl
    END SUBROUTINE ContravariantProjection_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantProjection_MappedVector3D_gpu_wrapper(physVector, compVector, dsdx, N, nVar, nEl) bind(c,name="ContravariantProjection_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: physVector, compVector, dsdx
      INTEGER, VALUE :: N, nVar, nEl
    END SUBROUTINE ContravariantProjection_MappedVector3D_gpu_wrapper
  END INTERFACE

CONTAINS


! ---------------------- Scalars ---------------------- !

  FUNCTION Gradient_MappedScalar2D( SELFStorage, workTensor, mesh, gpuAccel ) RESULT( SELFOut )
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
    LOGICAL :: gpuAccel
  
      CALL SELFStorage % ContravariantWeight( workTensor, mesh, gpuAccel ) 
      SELFOut = workTensor % CompDivergence( gpuAccel ) 
  
  END FUNCTION Gradient_MappedScalar2D
  
  SUBROUTINE ContravariantWeight_MappedScalar2D(SELFStorage, workTensor, mesh, gpuAccel)
    IMPLICIT NONE
    CLASS(MappedScalar2D), INTENT(in) :: SELFStorage
    TYPE(MappedTensor2D), INTENT(inout) :: workTensor
    TYPE(Mesh2D), INTENT(in) :: mesh
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: i, j, iVar, iEl
  
      IF( gpuAccel )THEN

        CALL ContravariantWeight_MappedScalar2D_gpu_wrapper( SELFStorage % interior % deviceData, &
                                                             workTensor % interior % deviceData, &
                                                             mesh % geometry % dsdx % interior % deviceData, &
                                                             SELFStorage % N, &
                                                             SELFStorage % nVar, &
                                                             SELFStorage % nElem )

      ELSE

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

      ENDIF
  
  END SUBROUTINE ContravariantWeight_MappedScalar2D
  
  FUNCTION Gradient_MappedScalar3D( SELFStorage, workTensor, mesh, gpuAccel ) RESULT( SELFOut )
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
    LOGICAL :: gpuAccel
  
      CALL SELFStorage % ContravariantWeight( workTensor, mesh, gpuAccel ) 
      SELFOut = workTensor % CompDivergence( gpuAccel ) 
  
  END FUNCTION Gradient_MappedScalar3D
  
  SUBROUTINE ContravariantWeight_MappedScalar3D( SELFStorage, workTensor, mesh, gpuAccel )
    IMPLICIT NONE
    CLASS(MappedScalar3D), INTENT(in) :: SELFStorage
    TYPE(MappedTensor3D), INTENT(inout) :: workTensor
    TYPE(Mesh3D), INTENT(in) :: mesh
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: i, j, k, iVar, iEl
  
      IF( gpuAccel )THEN

        CALL ContravariantWeight_MappedScalar3D_gpu_wrapper( SELFStorage % interior % deviceData, &
                                                             workTensor % interior % deviceData, &
                                                             mesh % geometry % dsdx % interior % deviceData, &
                                                             SELFStorage % N, &
                                                             SELFStorage % nVar, &
                                                             SELFStorage % nElem )
      ELSE

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

      ENDIF
  
  END SUBROUTINE ContravariantWeight_MappedScalar3D
  
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
  
  END FUNCTION CompDivergence_MappedTensor2D
  
  FUNCTION CompDivergence_MappedTensor3D( SELFStorage, gpuAccel ) RESULT( SELFout )
    IMPLICIT NONE
    CLASS(MappedTensor3D) :: SELFStorage
    TYPE(MappedVector3D) :: SELFOut
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
  
  END FUNCTION CompDivergence_MappedTensor3D
  
END MODULE SELF_MappedData
