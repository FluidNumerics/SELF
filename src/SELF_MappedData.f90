! SELF_MappedData.F90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE SELF_MappedData

  USE SELF_Constants
  USE SELF_Lagrange

  USE SELF_Data
  USE SELF_Mesh
  USE SELF_Geometry

  USE ISO_C_BINDING

  IMPLICIT NONE

#include "SELF_Macros.h"

  TYPE,EXTENDS(Scalar1D),PUBLIC :: MappedScalar1D

  CONTAINS
    GENERIC,PUBLIC :: Derivative => Derivative_MappedScalar1D
    PROCEDURE,PRIVATE :: Derivative_MappedScalar1D
    PROCEDURE,PRIVATE :: JacobianWeight => JacobianWeight_MappedScalar1D

  END TYPE MappedScalar1D

  TYPE,EXTENDS(Scalar2D),PUBLIC :: MappedScalar2D

  CONTAINS

    PROCEDURE,PUBLIC :: SideExchange => SideExchange_MappedScalar2D 
    PROCEDURE,PUBLIC :: BassiRebaySides => BassiRebaySides_MappedScalar2D 

    GENERIC,PUBLIC :: Gradient => Gradient_MappedScalar2D
    PROCEDURE,PRIVATE :: Gradient_MappedScalar2D
    PROCEDURE,PRIVATE :: ContravariantWeight => ContravariantWeight_MappedScalar2D
    PROCEDURE,PRIVATE :: JacobianWeight => JacobianWeight_MappedScalar2D

  END TYPE MappedScalar2D

  TYPE,EXTENDS(Scalar3D),PUBLIC :: MappedScalar3D

  CONTAINS

    PROCEDURE,PUBLIC :: SideExchange => SideExchange_MappedScalar3D
    PROCEDURE,PUBLIC :: BassiRebaySides => BassiRebaySides_MappedScalar3D
    GENERIC,PUBLIC :: Gradient => Gradient_MappedScalar3D
    PROCEDURE,PRIVATE :: Gradient_MappedScalar3D
    PROCEDURE,PRIVATE :: ContravariantWeight => ContravariantWeight_MappedScalar3D
    PROCEDURE,PRIVATE :: JacobianWeight => JacobianWeight_MappedScalar3D

  END TYPE MappedScalar3D

  TYPE,EXTENDS(Vector2D),PUBLIC :: MappedVector2D

  CONTAINS

    PROCEDURE,PUBLIC :: SideExchange => SideExchange_MappedVector2D
    PROCEDURE,PUBLIC :: BassiRebaySides => BassiRebaySides_MappedVector2D

    GENERIC,PUBLIC :: Divergence => Divergence_MappedVector2D
    GENERIC,PUBLIC :: Gradient => Gradient_MappedVector2D
!    GENERIC,PUBLIC :: Curl => Curl_MappedVector2D

    PROCEDURE,PRIVATE :: Divergence_MappedVector2D
    PROCEDURE,PRIVATE :: Gradient_MappedVector2D
!    PROCEDURE,PRIVATE :: Curl_MappedVector2D
    PROCEDURE,PRIVATE :: ContravariantProjection => ContravariantProjection_MappedVector2D
    PROCEDURE,PRIVATE :: JacobianWeight => JacobianWeight_MappedVector2D
    PROCEDURE,PRIVATE :: MapToScalar => MapToScalar_MappedVector2D
    PROCEDURE,PRIVATE :: MapToTensor => MapToTensor_MappedVector2D

  END TYPE MappedVector2D

  TYPE,EXTENDS(Vector3D),PUBLIC :: MappedVector3D

  CONTAINS

    PROCEDURE,PUBLIC :: SideExchange => SideExchange_MappedVector3D
    PROCEDURE,PUBLIC :: BassiRebaySides => BassiRebaySides_MappedVector3D 

    GENERIC,PUBLIC :: Divergence => Divergence_MappedVector3D
!    GENERIC,PUBLIC :: Curl => Curl_MappedVector3D
    GENERIC,PUBLIC :: Gradient => Gradient_MappedVector3D
    PROCEDURE,PRIVATE :: Divergence_MappedVector3D
!    PROCEDURE,PRIVATE :: Curl_MappedVector3D
    PROCEDURE,PRIVATE :: Gradient_MappedVector3D
    PROCEDURE,PRIVATE :: ContravariantProjection => ContravariantProjection_MappedVector3D
    PROCEDURE,PRIVATE :: JacobianWeight => JacobianWeight_MappedVector3D
    PROCEDURE,PRIVATE :: MapToScalar => MapToScalar_MappedVector3D
    PROCEDURE,PRIVATE :: MapToTensor => MapToTensor_MappedVector3D

  END TYPE MappedVector3D

  TYPE,EXTENDS(Tensor2D),PUBLIC :: MappedTensor2D

  CONTAINS

    PROCEDURE,PUBLIC :: SideExchange => SideExchange_MappedTensor2D
    PROCEDURE,PUBLIC :: BassiRebaySides => BassiRebaySides_MappedTensor2D

!    GENERIC,PUBLIC :: Divergence => Divergence_MappedTensor2D
!    PROCEDURE,PRIVATE :: Divergence_MappedTensor2D
    PROCEDURE,PRIVATE :: JacobianWeight => JacobianWeight_MappedTensor2D

  END TYPE MappedTensor2D

  TYPE,EXTENDS(Tensor3D),PUBLIC :: MappedTensor3D

  CONTAINS

    PROCEDURE,PUBLIC :: SideExchange => SideExchange_MappedTensor3D 
    PROCEDURE,PUBLIC :: BassiRebaySides => BassiRebaySides_MappedTensor3D


!    GENERIC,PUBLIC :: Divergence => Divergence_MappedTensor3D
!    PROCEDURE,PRIVATE :: Divergence_MappedTensor3D
    PROCEDURE,PRIVATE :: JacobianWeight => JacobianWeight_MappedTensor3D

  END TYPE MappedTensor3D

#ifdef GPU
  INTERFACE
    SUBROUTINE JacobianWeight_MappedScalar1D_gpu_wrapper(scalar,dxds,N,nVar,nEl) &
      bind(c,name="JacobianWeight_MappedScalar1D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,dxds
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE JacobianWeight_MappedScalar1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE JacobianWeight_MappedScalar2D_gpu_wrapper(scalar,jacobian,N,nVar,nEl) &
      bind(c,name="JacobianWeight_MappedScalar2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,jacobian
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE JacobianWeight_MappedScalar2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE JacobianWeight_MappedScalar3D_gpu_wrapper(scalar,jacobian,N,nVar,nEl) &
      bind(c,name="JacobianWeight_MappedScalar3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,jacobian
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE JacobianWeight_MappedScalar3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantWeight_MappedScalar2D_gpu_wrapper(scalar,workTensor,dsdx,N,nVar,nEl) &
      bind(c,name="ContravariantWeight_MappedScalar2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,workTensor,dsdx
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ContravariantWeight_MappedScalar2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantWeightBoundary_MappedScalar2D_gpu_wrapper(scalar,workTensor,dsdx,N,nVar,nEl) &
      bind(c,name="ContravariantWeightBoundary_MappedScalar2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,workTensor,dsdx
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ContravariantWeightBoundary_MappedScalar2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantWeight_MappedScalar3D_gpu_wrapper(scalar,workTensor,dsdx,N,nVar,nEl) &
      bind(c,name="ContravariantWeight_MappedScalar3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,workTensor,dsdx
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ContravariantWeight_MappedScalar3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantWeightBoundary_MappedScalar3D_gpu_wrapper(scalar,workTensor,dsdx,N,nVar,nEl) &
      bind(c,name="ContravariantWeightBoundary_MappedScalar3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,workTensor,dsdx
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ContravariantWeightBoundary_MappedScalar3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantProjection_MappedVector2D_gpu_wrapper(physVector,compVector,dsdx,N,nVar,nEl) &
      bind(c,name="ContravariantProjection_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: physVector,compVector,dsdx
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ContravariantProjection_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantProjectionBoundary_MappedVector2D_gpu_wrapper(physVector,compVector,dsdx,N,nVar,nEl) &
      bind(c,name="ContravariantProjectionBoundary_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: physVector,compVector,dsdx
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ContravariantProjectionBoundary_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE JacobianWeight_MappedVector2D_gpu_wrapper(vector,jacobian,N,nVar,nEl) &
      bind(c,name="JacobianWeight_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: vector,jacobian
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE JacobianWeight_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantProjection_MappedVector3D_gpu_wrapper(physVector,compVector,dsdx,N,nVar,nEl) &
      bind(c,name="ContravariantProjection_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: physVector,compVector,dsdx
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ContravariantProjection_MappedVector3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantProjectionBoundary_MappedVector3D_gpu_wrapper(physVector,compVector,dsdx,N,nVar,nEl) &
      bind(c,name="ContravariantProjectionBoundary_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: physVector,compVector,dsdx
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ContravariantProjectionBoundary_MappedVector3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE JacobianWeight_MappedVector3D_gpu_wrapper(vector,jacobian,N,nVar,nEl) &
      bind(c,name="JacobianWeight_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: vector,jacobian
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE JacobianWeight_MappedVector3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE JacobianWeight_MappedTensor2D_gpu_wrapper(tensor,jacobian,N,nVar,nEl) &
      bind(c,name="JacobianWeight_MappedTensor2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: tensor,jacobian
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE JacobianWeight_MappedTensor2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE JacobianWeight_MappedTensor3D_gpu_wrapper(tensor,jacobian,N,nVar,nEl) &
      bind(c,name="JacobianWeight_MappedTensor3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: tensor,jacobian
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE JacobianWeight_MappedTensor3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE MapToScalar_MappedVector2D_gpu_wrapper(scalar,vector,N,nVar,nEl) &
      bind(c,name="MapToScalar_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,vector
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE MapToScalar_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE MapToScalarBoundary_MappedVector2D_gpu_wrapper(scalar,vector,N,nVar,nEl) &
      bind(c,name="MapToScalarBoundary_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,vector
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE MapToScalarBoundary_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE MapToTensor_MappedVector2D_gpu_wrapper(tensor,vector,N,nVar,nEl) &
      bind(c,name="MapToTensor_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: tensor,vector
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE MapToTensor_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE MapToTensorBoundary_MappedVector2D_gpu_wrapper(tensor,vector,N,nVar,nEl) &
      bind(c,name="MapToTensorBoundary_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: tensor,vector
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE MapToTensorBoundary_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE MapToScalar_MappedVector3D_gpu_wrapper(scalar,vector,N,nVar,nEl) &
      bind(c,name="MapToScalar_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,vector
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE MapToScalar_MappedVector3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE MapToScalarBoundary_MappedVector3D_gpu_wrapper(scalar,vector,N,nVar,nEl) &
      bind(c,name="MapToScalarBoundary_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,vector
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE MapToScalarBoundary_MappedVector3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE MapToTensor_MappedVector3D_gpu_wrapper(tensor,vector,N,nVar,nEl) &
      bind(c,name="MapToTensor_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: tensor,vector
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE MapToTensor_MappedVector3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE MapToTensorBoundary_MappedVector3D_gpu_wrapper(tensor,vector,N,nVar,nEl) &
      bind(c,name="MapToTensorBoundary_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: tensor,vector
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE MapToTensorBoundary_MappedVector3D_gpu_wrapper
  END INTERFACE
#endif

CONTAINS

! ---------------------- Scalars ---------------------- !

  SUBROUTINE Derivative_MappedScalar1D(scalar,geometry,dF,dForm,gpuAccel)
    ! Strong Form Operator
    !
    ! Calculates the gradient of a scalar 1D function using the conservative form of the
    ! mapped gradient operator
    !
    ! df/dx =  d\xi/dx( df/d\xi )
    !
    ! where the sum over i is implied.
    IMPLICIT NONE
    CLASS(MappedScalar1D),INTENT(in) :: scalar
    TYPE(Geometry1D),INTENT(in) :: geometry
    TYPE(MappedScalar1D),INTENT(inout) :: dF
    INTEGER,INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

    IF (dForm == selfWeakDGForm) THEN

      IF (gpuAccel) THEN
        CALL scalar % interp % DGDerivative_1D(scalar % interior % deviceData, &
                                               scalar % boundary % deviceData, &
                                               df % interior % deviceData, &
                                               scalar % nVar, &
                                               scalar % nElem)
      ELSE
        CALL scalar % interp % DGDerivative_1D(scalar % interior % hostData, &
                                               scalar % boundary % hostData, &
                                               df % interior % hostData, &
                                               scalar % nVar, &
                                               scalar % nElem)
      END IF

    ELSEIF (dForm == selfStrongForm) THEN

      IF (gpuAccel) THEN
        CALL scalar % interp % Derivative_1D(scalar % interior % deviceData, &
                                             df % interior % deviceData, &
                                             scalar % nVar, &
                                             scalar % nElem)
      ELSE
        CALL scalar % interp % Derivative_1D(scalar % interior % hostData, &
                                             df % interior % hostData, &
                                             scalar % nVar, &
                                             scalar % nElem)
      END IF

    ENDIF

    CALL df % JacobianWeight(geometry,gpuAccel)

  END SUBROUTINE Derivative_MappedScalar1D

  SUBROUTINE JacobianWeight_MappedScalar1D(scalar,geometry,gpuAccel)
#undef __FUNC__
#define __FUNC__ "JacobianWeight_MappedScalar1D"
  ! Applies the inverse jacobian
    IMPLICIT NONE
    CLASS(MappedScalar1D),INTENT(inout) :: scalar
    TYPE(Geometry1D),INTENT(in) :: geometry
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iEl, iVar, i
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN
#ifdef GPU
      CALL JacobianWeight_MappedScalar1D_gpu_wrapper(scalar % interior % deviceData, &
                                                     geometry % dxds % interior % deviceData, &
                                                     scalar % N, &
                                                     scalar % nVar, &
                                                     scalar % nElem)
#else
      msg = "GPU Acceleration is not currently enabled in SELF."
      WARNING(msg)
#endif
    ELSE

      DO iEl = 1,scalar % nElem
        DO iVar = 1,scalar % nVar
          DO i = 0,scalar % N
            scalar % interior % hostData(i,iVar,iEl) = scalar % interior % hostData(i,iVar,iEl)/&
                                                       geometry % dxds % interior % hostData(i,1,iEl)
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE JacobianWeight_MappedScalar1D

  SUBROUTINE SideExchange_MappedScalar2D(scalar,mesh,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedScalar2D),INTENT(inout) :: scalar
    TYPE(Mesh2D),INTENT(in) :: mesh
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1, e2, s1, s2, sid 
    INTEGER :: flip, bcid, globalSideId
    INTEGER :: i1, i2, ivar

    IF(gpuAccel)THEN

      ! TO DO ! 
      PRINT*, 'Woopsie! GPU Acceleration not implemented yet for SideExchange'

    ELSE

      DO e1 = 1, mesh % nElem
        s1 = 0
        DO sid = mesh % elemInfo % hostData(3,e1)+1, mesh % elemInfo % hostData(4,e1) ! Loop over local sides 
          s1 = s1 + 1 ! Increment local side ID 
          globalSideId = mesh % sideInfo % hostData(2,sid)
          e2 = mesh % sideInfo % hostData(3,sid)
          s2 = mesh % sideInfo % hostData(4,sid)/10
          flip = mesh % sideInfo % hostData(4,sid)-s2*10
          bcid = mesh % sideInfo % hostData(5,sid)

          IF(bcid /= 0)THEN   

            IF(flip == 1)THEN 
          
              DO ivar = 1, scalar % nvar
                DO i1 = 0, scalar % N
                  scalar % extBoundary % hostData(i1,ivar,s1,e1) = &
                      scalar % boundary % hostData(i1,ivar,s2,e2)
                ENDDO
              ENDDO

            ELSEIF(flip == 2)THEN

              DO ivar = 1, scalar % nvar
                DO i1 = 0, scalar % N
                  i2 = scalar % N - i1
                  scalar % extBoundary % hostData(i1,ivar,s1,e1) = &
                     scalar % boundary % hostData(i2,ivar,s2,e2)
                ENDDO
              ENDDO

            ENDIF

          ENDIF

        ENDDO
      ENDDO

    END IF
    
  END SUBROUTINE SideExchange_MappedScalar2D

  SUBROUTINE BassiRebaySides_MappedScalar2D(scalar,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedScalar2D),INTENT(inout) :: scalar
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iel
    INTEGER :: iside
    INTEGER :: ivar
    INTEGER :: i

    IF(gpuAccel)THEN

      ! TO DO ! 
      PRINT*, 'Woopsie! GPU Acceleration not implemented yet for BassiRebay'

    ELSE

      DO iel = 1, scalar % nElem
        DO iside = 1, 4
          DO ivar = 1, scalar % nVar
            DO i = 0, scalar % N
              scalar % boundary % hostData(i,ivar,iside,iel) = 0.5_prec*(&
                scalar % boundary % hostData(i,ivar,iside,iel) +&
                scalar % extBoundary % hostData(i,ivar,iside,iel))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE BassiRebaySides_MappedScalar2D

  SUBROUTINE Gradient_MappedScalar2D(scalar,workTensor,geometry,gradF,dForm,gpuAccel)
    ! Strong Form Operator - (Conservative Form)
    !
    ! Calculates the gradient of a scalar 2D function using the conservative form of the
    ! mapped gradient operator
    !
    ! \grad_{phys}( f ) =  (1 / J)*(\partial / \partial \xi_i ( J\vec{a}_i f )
    !
    ! where the sum over i is implied.
    IMPLICIT NONE
    CLASS(MappedScalar2D),INTENT(in) :: scalar
    TYPE(MappedTensor2D),INTENT(inout) :: workTensor
    TYPE(SEMQuad),INTENT(in) :: geometry
    TYPE(MappedVector2D),INTENT(inout) :: gradF
    INTEGER,INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

    CALL scalar % ContravariantWeight(geometry,workTensor,gpuAccel)

    IF (dForm == selfWeakDGForm) THEN

      IF (gpuAccel) THEN
        CALL workTensor % interp % TensorDGDivergence_2D(workTensor % interior % deviceData, &
                                                         workTensor % boundary % deviceData, &
                                                         gradF % interior % deviceData, &
                                                         workTensor % nVar, &
                                                         workTensor % nElem)
      ELSE
        CALL workTensor % interp % TensorDGDivergence_2D(workTensor % interior % hostData, &
                                                         workTensor % boundary % hostData, &
                                                         gradF % interior % hostData, &
                                                         workTensor % nVar, &
                                                         workTensor % nElem)
      END IF

    ELSEIF (dForm == selfStrongForm) THEN

      IF (gpuAccel) THEN
        CALL workTensor % interp % TensorDivergence_2D(workTensor % interior % deviceData, &
                                                       gradF % interior % deviceData, &
                                                       workTensor % nVar, &
                                                       workTensor % nElem)
      ELSE
        CALL workTensor % interp % TensorDivergence_2D(workTensor % interior % hostData, &
                                                       gradF % interior % hostData, &
                                                       workTensor % nVar, &
                                                       workTensor % nElem)
      END IF

    ENDIF

    CALL gradF % JacobianWeight(geometry,gpuAccel)

  END SUBROUTINE Gradient_MappedScalar2D

  SUBROUTINE ContravariantWeight_MappedScalar2D(scalar,geometry,workTensor,gpuAccel)
#undef __FUNC__
#define __FUNC__ "ContravariantWeight_MappedScalar2D"
    IMPLICIT NONE
    CLASS(MappedScalar2D),INTENT(in) :: scalar
    TYPE(SEMQuad),INTENT(in) :: geometry
    LOGICAL,INTENT(in):: gpuAccel
    TYPE(MappedTensor2D),INTENT(inout) :: workTensor
    ! Local
    INTEGER :: i,j,iVar,iEl,iside
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN
#ifdef GPU
      CALL ContravariantWeight_MappedScalar2D_gpu_wrapper(scalar % interior % deviceData, &
                                                          workTensor % interior % deviceData, &
                                                          geometry % dsdx % interior % deviceData, &
                                                          scalar % N, &
                                                          scalar % nVar, &
                                                          scalar % nElem)

      CALL ContravariantWeightBoundary_MappedScalar2D_gpu_wrapper(scalar % boundary % deviceData, &
                                                          workTensor % boundary % deviceData, &
                                                          geometry % dsdx % boundary % deviceData, &
                                                          scalar % N, &
                                                          scalar % nVar, &
                                                          scalar % nElem)
#else
      msg = "GPU Acceleration is not currently enabled in SELF."
      WARNING(msg)
#endif
    ELSE

      DO iEl = 1,scalar % nElem
        DO iVar = 1,scalar % nVar
          DO j = 0,scalar % N
            DO i = 0,scalar % N

              workTensor % interior % hostData(1,1,i,j,iVar,iEl) = geometry % dsdx % interior % &
                                                                   hostData(1,1,i,j,1,iEl)* &
                                                                   scalar % interior % hostData(i,j,iVar,iEl)

              workTensor % interior % hostData(2,1,i,j,iVar,iEl) = geometry % dsdx % interior % &
                                                                   hostData(1,2,i,j,1,iEl)* &
                                                                   scalar % interior % hostData(i,j,iVar,iEl)

              workTensor % interior % hostData(1,2,i,j,iVar,iEl) = geometry % dsdx % interior % &
                                                                   hostData(2,1,i,j,1,iEl)* &
                                                                   scalar % interior % hostData(i,j,iVar,iEl)

              workTensor % interior % hostData(2,2,i,j,iVar,iEl) = geometry % dsdx % interior % &
                                                                   hostData(2,2,i,j,1,iEl)* &
                                                                   scalar % interior % hostData(i,j,iVar,iEl)

            END DO
          END DO
        END DO
      END DO

      ! Boundary Terms
      DO iEl = 1,scalar % nElem
        DO iside = 1,4
          DO iVar = 1,scalar % nVar
            DO j = 0,scalar % N
              workTensor % boundary % hostData(1,1,j,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                   hostData(1,1,j,1,iside,iEl)* &
                                                                   scalar % boundary % hostData(j,iVar,iside,iEl)

              workTensor % boundary % hostData(2,1,j,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                   hostData(1,2,j,1,iside,iEl)* &
                                                                   scalar % boundary % hostData(j,iVar,iside,iEl)

              workTensor % boundary % hostData(1,2,j,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                   hostData(2,1,j,1,iside,iEl)* &
                                                                   scalar % boundary % hostData(j,iVar,iside,iEl)

              workTensor % boundary % hostData(2,2,j,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                   hostData(2,2,j,1,iside,iEl)* &
                                                                   scalar % boundary % hostData(j,iVar,iside,iEl)

            ENDDO
          ENDDO
        ENDDO
      ENDDO

    END IF

  END SUBROUTINE ContravariantWeight_MappedScalar2D

  SUBROUTINE JacobianWeight_MappedScalar2D(scalar,geometry,gpuAccel)
#undef __FUNC__
#define __FUNC__ "JacobianWeight_MappedScalar2D"
  ! Applies the inverse jacobian
    IMPLICIT NONE
    CLASS(MappedScalar2D),INTENT(inout) :: scalar
    TYPE(SEMQuad),INTENT(in) :: geometry
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iEl,iVar,i,j
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN
#ifdef GPU
      CALL JacobianWeight_MappedScalar2D_gpu_wrapper(scalar % interior % deviceData, &
                                                     geometry % J % interior % deviceData, &
                                                     scalar % N, &
                                                     scalar % nVar, &
                                                     scalar % nElem)
#else
      msg = "GPU Acceleration is not currently enabled in SELF."
      WARNING(msg)
#endif
    ELSE

      DO iEl = 1,scalar % nElem
        DO iVar = 1,scalar % nVar
          DO j = 0,scalar % N
            DO i = 0,scalar % N
              scalar % interior % hostData(i,j,iVar,iEl) = scalar % interior % hostData(i,j,iVar,iEl)/&
                                                           geometry % J % interior % hostData(i,j,1,iEl)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE JacobianWeight_MappedScalar2D

  ! SideExchange_MappedScalar3D is used to populate scalar % extBoundary
  ! by finding neighboring elements that share a side and copying the neighboring
  ! elements solution % boundary data.

  SUBROUTINE SideExchange_MappedScalar3D(scalar,mesh,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedScalar3D),INTENT(inout) :: scalar
    TYPE(Mesh3D),INTENT(in) :: mesh
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1, e2, s1, s2, sid 
    INTEGER :: flip, bcid, globalSideId
    INTEGER :: i1, i2, j1, j2, ivar

    IF(gpuAccel)THEN

      ! TO DO ! 
      PRINT*, 'Woopsie! GPU Acceleration not implemented yet for SideExchange'

    ELSE

      DO e1 = 1, mesh % nElem
        s1 = 0
        DO sid = mesh % elemInfo % hostData(3,e1)+1, mesh % elemInfo % hostData(4,e1) ! Loop over local sides 
          s1 = s1 + 1 ! Increment local side ID 
          globalSideId = mesh % sideInfo % hostData(2,sid)
          e2 = mesh % sideInfo % hostData(3,sid)
          s2 = mesh % sideInfo % hostData(4,sid)/10
          flip = mesh % sideInfo % hostData(4,sid)-s2*10
          bcid = mesh % sideInfo % hostData(5,sid)

          IF(bcid /= 0)THEN   

            IF(flip == 1)THEN 
          
              DO ivar = 1, scalar % nvar
                DO j1 = 0, scalar % N
                  DO i1 = 0, scalar % N
                    scalar % extBoundary % hostData(i1,j1,ivar,s1,e1) = &
                      scalar % boundary % hostData(i1,j1,ivar,s2,e2)
                  ENDDO
                ENDDO
              ENDDO

            ELSEIF(flip == 2)THEN

              DO ivar = 1, scalar % nvar
                DO j1 = 0, scalar % N
                  DO i1 = 0, scalar % N
                    i2 = scalar % N - j1
                    j2 = i1
                    scalar % extBoundary % hostData(i1,j1,ivar,s1,e1) = &
                      scalar % boundary % hostData(i2,j2,ivar,s2,e2)
                  ENDDO
                ENDDO
              ENDDO

            ELSEIF(flip == 3)THEN
                    
              DO ivar = 1, scalar % nvar
                DO j1 = 0, scalar % N
                  DO i1 = 0, scalar % N
                    i2 = scalar % N - i1
                    j2 = scalar % N - j1 
                    scalar % extBoundary % hostData(i1,j1,ivar,s1,e1) = &
                      scalar % boundary % hostData(i2,j2,ivar,s2,e2)
                  ENDDO
                ENDDO
              ENDDO
          
            ELSEIF(flip == 4)THEN
                    
              DO ivar = 1, scalar % nvar
                DO j1 = 0, scalar % N
                  DO i1 = 0, scalar % N
                    i2 = j1
                    j2 = scalar % N - i1
                    scalar % extBoundary % hostData(i1,j1,ivar,s1,e1) = &
                      scalar % boundary % hostData(i2,j2,ivar,s2,e2)
                  ENDDO
                ENDDO
              ENDDO
          
            ENDIF

          ENDIF

        ENDDO
      ENDDO

    END IF
    
  END SUBROUTINE SideExchange_MappedScalar3D

  SUBROUTINE BassiRebaySides_MappedScalar3D(scalar,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedScalar3D),INTENT(inout) :: scalar
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iel
    INTEGER :: iside
    INTEGER :: ivar
    INTEGER :: i, j

    IF(gpuAccel)THEN

      ! TO DO ! 
      PRINT*, 'Woopsie! GPU Acceleration not implemented yet for BassiRebay'

    ELSE

      DO iel = 1, scalar % nElem
        DO iside = 1, 6
          DO ivar = 1, scalar % nVar
            DO j = 0, scalar % N
              DO i = 0, scalar % N
                scalar % boundary % hostData(i,j,ivar,iside,iel) = 0.5_prec*(&
                  scalar % boundary % hostData(i,j,ivar,iside,iel) +&
                  scalar % extBoundary % hostData(i,j,ivar,iside,iel))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE BassiRebaySides_MappedScalar3D

  SUBROUTINE Gradient_MappedScalar3D(scalar,workTensor,geometry,gradF,dForm,gpuAccel)
    ! Strong Form Operator
    !
    ! Calculates the gradient of a scalar 3D function using the conservative form of the
    ! mapped gradient operator
    !
    ! \grad_{phys}( f ) =  (1 / J)*(\partial / \partial \xi_i ( J\vec{a}_i f )
    !
    ! where the sum over i is implied.
    IMPLICIT NONE
    CLASS(MappedScalar3D),INTENT(in) :: scalar
    TYPE(MappedTensor3D),INTENT(inout) :: workTensor
    TYPE(SEMHex),INTENT(in) :: geometry
    TYPE(MappedVector3D),INTENT(inout) :: gradF
    INTEGER,INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

    CALL scalar % ContravariantWeight(geometry,workTensor,gpuAccel)
    IF (dForm == selfWeakDGForm) THEN

      IF (gpuAccel) THEN
        CALL workTensor % interp % TensorDGDivergence_3D(workTensor % interior % deviceData, &
                                                         workTensor % boundary % deviceData, &
                                                         gradF % interior % deviceData, &
                                                         workTensor % nVar, &
                                                         workTensor % nElem)
      ELSE
        CALL workTensor % interp % TensorDGDivergence_3D(workTensor % interior % hostData, &
                                                         workTensor % boundary % hostData, &
                                                         gradF % interior % hostData, &
                                                         workTensor % nVar, &
                                                         workTensor % nElem)
      END IF

    ELSE IF (dForm == selfStrongForm) THEN

      IF (gpuAccel) THEN
        CALL workTensor % interp % TensorDivergence_3D(workTensor % interior % deviceData, &
                                                       gradF % interior % deviceData, &
                                                       workTensor % nVar, &
                                                       workTensor % nElem)
      ELSE
        CALL workTensor % interp % TensorDivergence_3D(workTensor % interior % hostData, &
                                                       gradF % interior % hostData, &
                                                       workTensor % nVar, &
                                                       workTensor % nElem)
      END IF

    END IF
    CALL gradF % JacobianWeight(geometry,gpuAccel)

  END SUBROUTINE Gradient_MappedScalar3D

  SUBROUTINE ContravariantWeight_MappedScalar3D(scalar,geometry,workTensor,gpuAccel)
#undef __FUNC__
#define __FUNC__ "ContravariantWeight_MappedScalar3D"
    IMPLICIT NONE
    CLASS(MappedScalar3D),INTENT(in) :: scalar
    TYPE(SEMHex),INTENT(in) :: geometry
    LOGICAL,INTENT(in) :: gpuAccel
    TYPE(MappedTensor3D),INTENT(inout) :: workTensor
    ! Local
    INTEGER :: i,j,k,iVar,iEl,iside
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN
#ifdef GPU
      CALL ContravariantWeight_MappedScalar3D_gpu_wrapper(scalar % interior % deviceData, &
                                                          workTensor % interior % deviceData, &
                                                          geometry % dsdx % interior % deviceData, &
                                                          scalar % N, &
                                                          scalar % nVar, &
                                                          scalar % nElem)

      CALL ContravariantWeightBoundary_MappedScalar3D_gpu_wrapper(scalar % boundary % deviceData, &
                                                          workTensor % boundary % deviceData, &
                                                          geometry % dsdx % boundary % deviceData, &
                                                          scalar % N, &
                                                          scalar % nVar, &
                                                          scalar % nElem)
#else
      msg = "GPU Acceleration is not currently enabled in SELF."
      WARNING(msg)
#endif
    ELSE

      DO iEl = 1,scalar % nElem
        DO iVar = 1,scalar % nVar
          DO k = 0,scalar % N
            DO j = 0,scalar % N
              DO i = 0,scalar % N

                ! Get the x-component of the Jacobian weighted
                ! contravariant basis vectors multipled by the scalar
                workTensor % interior % hostData(1,1,i,j,k,iVar,iEl) = geometry % dsdx % interior % &
                                                                       hostData(1,1,i,j,k,1,iEl)* &
                                                                       scalar % interior % hostData(i,j,k,iVar,iEl)

                workTensor % interior % hostData(2,1,i,j,k,iVar,iEl) = geometry % dsdx % interior % &
                                                                       hostData(1,2,i,j,k,1,iEl)* &
                                                                       scalar % interior % hostData(i,j,k,iVar,iEl)

                workTensor % interior % hostData(3,1,i,j,k,iVar,iEl) = geometry % dsdx % interior % &
                                                                       hostData(1,3,i,j,k,1,iEl)* &
                                                                       scalar % interior % hostData(i,j,k,iVar,iEl)

                ! Get the y-component of the Jacobian weighted
                ! contravariant basis vectors multipled by the scalar
                workTensor % interior % hostData(1,2,i,j,k,iVar,iEl) = geometry % dsdx % interior % &
                                                                       hostData(2,1,i,j,k,1,iEl)* &
                                                                       scalar % interior % hostData(i,j,k,iVar,iEl)

                workTensor % interior % hostData(2,2,i,j,k,iVar,iEl) = geometry % dsdx % interior % &
                                                                       hostData(2,2,i,j,k,1,iEl)* &
                                                                       scalar % interior % hostData(i,j,k,iVar,iEl)

                workTensor % interior % hostData(3,2,i,j,k,iVar,iEl) = geometry % dsdx % interior % &
                                                                       hostData(2,3,i,j,k,1,iEl)* &
                                                                       scalar % interior % hostData(i,j,k,iVar,iEl)

                ! Get the z-component of the Jacobian weighted
                ! contravariant basis vectors multipled by the scalar
                workTensor % interior % hostData(1,3,i,j,k,iVar,iEl) = geometry % dsdx % interior % &
                                                                       hostData(3,1,i,j,k,1,iEl)* &
                                                                       scalar % interior % hostData(i,j,k,iVar,iEl)

                workTensor % interior % hostData(2,3,i,j,k,iVar,iEl) = geometry % dsdx % interior % &
                                                                       hostData(3,2,i,j,k,1,iEl)* &
                                                                       scalar % interior % hostData(i,j,k,iVar,iEl)

                workTensor % interior % hostData(3,3,i,j,k,iVar,iEl) = geometry % dsdx % interior % &
                                                                       hostData(3,3,i,j,k,1,iEl)* &
                                                                       scalar % interior % hostData(i,j,k,iVar,iEl)

              END DO
            END DO
          END DO
        END DO
      END DO

      ! Boundary Term
      DO iEl = 1,scalar % nElem
        DO iside = 1,6
          DO iVar = 1,scalar % nVar
            DO k = 0,scalar % N
              DO j = 0,scalar % N
                ! Get the x-component of the Jacobian weighted
                ! contravariant basis vectors multipled by the scalar
                workTensor % boundary % hostData(1,1,j,k,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                       hostData(1,1,j,k,1,iside,iEl)* &
                                                                       scalar % boundary % hostData(j,k,iVar,iside,iEl)

                workTensor % boundary % hostData(2,1,j,k,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                       hostData(1,2,j,k,1,iside,iEl)* &
                                                                       scalar % boundary % hostData(j,k,iVar,iside,iEl)

                workTensor % boundary % hostData(3,1,j,k,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                       hostData(1,3,j,k,1,iside,iEl)* &
                                                                       scalar % boundary % hostData(j,k,iVar,iside,iEl)

                ! Get the y-component of the Jacobian weighted
                ! contravariant basis vectors multipled by the scalar
                workTensor % boundary % hostData(1,2,j,k,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                       hostData(2,1,j,k,1,iside,iEl)* &
                                                                       scalar % boundary % hostData(j,k,iVar,iside,iEl)

                workTensor % boundary % hostData(2,2,j,k,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                       hostData(2,2,j,k,1,iside,iEl)* &
                                                                       scalar % boundary % hostData(j,k,iVar,iside,iEl)

                workTensor % boundary % hostData(3,2,j,k,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                       hostData(2,3,j,k,1,iside,iEl)* &
                                                                       scalar % boundary % hostData(j,k,iVar,iside,iEl)

                ! Get the z-component of the Jacobian weighted
                ! contravariant basis vectors multipled by the scalar
                workTensor % boundary % hostData(1,3,j,k,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                       hostData(3,1,j,k,1,iside,iEl)* &
                                                                       scalar % boundary % hostData(j,k,iVar,iside,iEl)

                workTensor % boundary % hostData(2,3,j,k,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                       hostData(3,2,j,k,1,iside,iEl)* &
                                                                       scalar % boundary % hostData(j,k,iVar,iside,iEl)

                workTensor % boundary % hostData(3,3,j,k,iVar,iside,iEl) = geometry % dsdx % boundary % &
                                                                       hostData(3,3,j,k,1,iside,iEl)* &
                                                                       scalar % boundary % hostData(j,k,iVar,iside,iEl)


              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    END IF

  END SUBROUTINE ContravariantWeight_MappedScalar3D

  SUBROUTINE JacobianWeight_MappedScalar3D(scalar,geometry,gpuAccel)
#undef __FUNC__
#define __FUNC__ "JacobianWeight_MappedScalar3D"
  ! Applies the inverse jacobian
    IMPLICIT NONE
    CLASS(MappedScalar3D),INTENT(inout) :: scalar
    TYPE(SEMHex),INTENT(in) :: geometry
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iEl,iVar,i,j,k
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN
#ifdef GPU
      CALL JacobianWeight_MappedScalar3D_gpu_wrapper(scalar % interior % deviceData, &
                                                     geometry % J % interior % deviceData, &
                                                     scalar % N, &
                                                     scalar % nVar, &
                                                     scalar % nElem)
#else
      msg = "GPU Acceleration is not currently enabled in SELF."
      WARNING(msg)
#endif

    ELSE

      DO iEl = 1,scalar % nElem
        DO iVar = 1,scalar % nVar
          DO k = 0,scalar % N
            DO j = 0,scalar % N
              DO i = 0,scalar % N
                scalar % interior % hostData(i,j,k,iVar,iEl) = scalar % interior % hostData(i,j,k,iVar,iEl)/&
                                                               geometry % J % interior % hostData(i,j,k,1,iEl)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE JacobianWeight_MappedScalar3D

  ! ---------------------- Vectors ---------------------- !
  ! SideExchange_MappedVectorvector2D is used to populate vector % extBoundary
  ! by finding neighboring elements that share a side and copying the neighboring
  ! elements solution % boundary data.

  SUBROUTINE SideExchange_MappedVector2D(vector,mesh,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(inout) :: vector
    TYPE(Mesh2D),INTENT(in) :: mesh
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1, e2, s1, s2, sid 
    INTEGER :: flip, bcid, globalSideId
    INTEGER :: i1, i2, ivar

    IF(gpuAccel)THEN

      ! TO DO ! 
      PRINT*, 'Woopsie! GPU Acceleration not implemented yet for SideExchange'

    ELSE

      DO e1 = 1, mesh % nElem
        s1 = 0
        DO sid = mesh % elemInfo % hostData(3,e1)+1, mesh % elemInfo % hostData(4,e1) ! Loop over local sides 
          s1 = s1 + 1 ! Increment local side ID 
          globalSideId = mesh % sideInfo % hostData(2,sid)
          e2 = mesh % sideInfo % hostData(3,sid)
          s2 = mesh % sideInfo % hostData(4,sid)/10
          flip = mesh % sideInfo % hostData(4,sid)-s2*10
          bcid = mesh % sideInfo % hostData(5,sid)

          IF(bcid /= 0)THEN   

            IF(flip == 1)THEN 
          
              DO ivar = 1, vector % nvar
                DO i1 = 0, vector % N
                  vector % extBoundary % hostData(1:2,i1,ivar,s1,e1) = &
                      vector % boundary % hostData(1:2,i1,ivar,s2,e2)
                ENDDO
              ENDDO

            ELSEIF(flip == 2)THEN

              DO ivar = 1, vector % nvar
                DO i1 = 0, vector % N
                  i2 = vector % N - i1
                  vector % extBoundary % hostData(1:2,i1,ivar,s1,e1) = &
                      vector % boundary % hostData(1:2,i2,ivar,s2,e2)
                ENDDO
              ENDDO

            ENDIF

          ENDIF

        ENDDO
      ENDDO

    END IF
    
  END SUBROUTINE SideExchange_MappedVector2D

  SUBROUTINE BassiRebaySides_MappedVector2D(vector,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(inout) :: vector
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iel
    INTEGER :: iside
    INTEGER :: ivar
    INTEGER :: i

    IF(gpuAccel)THEN

      ! TO DO ! 
      PRINT*, 'Woopsie! GPU Acceleration not implemented yet for BassiRebay'

    ELSE

      DO iel = 1, vector % nElem
        DO iside = 1, 4
          DO ivar = 1, vector % nVar
            DO i = 0, vector % N
              vector % boundary % hostData(1:2,i,ivar,iside,iel) = 0.5_prec*(&
                vector % boundary % hostData(1:2,i,ivar,iside,iel) +&
                vector % extBoundary % hostData(1:2,i,ivar,iside,iel))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF
           
  END SUBROUTINE BassiRebaySides_MappedVector2D

  SUBROUTINE Divergence_MappedVector2D(physVector,compVector,geometry,divVector,dForm,gpuAccel)
    ! Strong Form Operator
    !
    ! DG Weak Form Operator
    !
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(in) :: physVector
    TYPE(MappedVector2D),INTENT(inout) :: compVector
    TYPE(SEMQuad),INTENT(in) :: geometry
    TYPE(MappedScalar2D),INTENT(inout) :: divVector
    INTEGER,INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

    CALL physVector % ContravariantProjection(geometry,compVector,gpuAccel)
    IF (dForm == selfWeakDGForm) THEN

      IF (gpuAccel) THEN
        CALL compVector % interp % VectorDGDivergence_2D(compVector % interior % deviceData, &
                                                         compVector % boundary % deviceData, &
                                                         divVector % interior % deviceData, &
                                                         compVector % nvar, &
                                                         compVector % nelem)
      ELSE
        CALL compVector % interp % VectorDGDivergence_2D(compVector % interior % hostData, &
                                                         compVector % boundary % hostData, &
                                                         divVector % interior % hostData, &
                                                         compVector % nvar, &
                                                         compVector % nelem)
      ENDIF

    ELSE IF (dForm == selfStrongForm) THEN

      IF (gpuAccel) THEN
        CALL compVector % interp % VectorDivergence_2D(compVector % interior % deviceData, &
                                                       divVector % interior % deviceData, &
                                                       compVector % nvar, &
                                                       compVector % nelem)
      ELSE
        CALL compVector % interp % VectorDivergence_2D(compVector % interior % hostData, &
                                                       divVector % interior % hostData, &
                                                       compVector % nvar, &
                                                       compVector % nelem)
      ENDIF

    END IF

    CALL divVector % JacobianWeight(geometry,gpuAccel)

  END SUBROUTINE Divergence_MappedVector2D

  SUBROUTINE Gradient_MappedVector2D(vector,workScalar,workVector,workTensor,geometry,gradF,dForm,gpuAccel)
    ! Strong Form Operator - (Conservative Form)
    !
    ! Calculates the gradient of a scalar 2D function using the conservative form of the
    ! mapped gradient operator
    !
    ! \grad_{phys}( f ) =  (1 / J)*(\partial / \partial \xi_i ( J\vec{a}_i f )
    !
    ! where the sum over i is implied.
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(in) :: vector
    TYPE(MappedScalar2D),INTENT(inout) :: workScalar ! (scalar) nvar = 2*nvar
    TYPE(MappedVector2D),INTENT(inout) :: workVector ! (scalar) nvar = 2*nvar
    TYPE(MappedTensor2D),INTENT(inout) :: workTensor ! (tensor) nvar = 2*nvar
    TYPE(SEMQuad),INTENT(in) :: geometry
    TYPE(MappedTensor2D),INTENT(inout) :: gradF
    INTEGER,INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

    CALL vector % MapToScalar(workScalar,gpuAccel)

    CALL workScalar % ContravariantWeight(geometry,workTensor,gpuAccel)

    IF (dForm == selfWeakDGForm) THEN

      IF (gpuAccel) THEN
        CALL workTensor % interp % TensorDGDivergence_2D(workTensor % interior % deviceData, &
                                                         workTensor % boundary % deviceData, &
                                                         workVector % interior % deviceData, &
                                                         workTensor % nVar, &
                                                         workTensor % nElem)
      ELSE
        CALL workTensor % interp % TensorDGDivergence_2D(workTensor % interior % hostData, &
                                                         workTensor % boundary % hostData, &
                                                         workVector % interior % hostData, &
                                                         workTensor % nVar, &
                                                         workTensor % nElem)
      END IF

    ELSE IF (dForm == selfStrongForm) THEN

      IF (gpuAccel) THEN
        CALL workTensor % interp % TensorDivergence_2D(workTensor % interior % deviceData, &
                                                       workVector % interior % deviceData, &
                                                       workTensor % nVar, &
                                                       workTensor % nElem)
      ELSE
        CALL workTensor % interp % TensorDivergence_2D(workTensor % interior % hostData, &
                                                       workVector % interior % hostData, &
                                                       workTensor % nVar, &
                                                       workTensor % nElem)
      END IF

    END IF

    CALL workVector % MapToTensor(gradF,gpuAccel)

    CALL gradF % JacobianWeight(geometry,gpuAccel)

  END SUBROUTINE Gradient_MappedVector2D

  SUBROUTINE MapToScalar_MappedVector2D(vector,scalar,gpuAccel)
#undef __FUNC__
#define __FUNC__ "MapToScalar_MappedVector2D"
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(in) :: vector
    TYPE(MappedScalar2D),INTENT(inout) :: scalar
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: row,i,j,ivar,jvar,iel,iside
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN
#ifdef GPU
      CALL MapToScalar_MappedVector2D_gpu_wrapper( scalar % interior % deviceData,&
                                                   vector % interior % deviceData,&
                                                   vector % N, &
                                                   vector % nVar, &
                                                   vector % nelem )

      CALL MapToScalarBoundary_MappedVector2D_gpu_wrapper( scalar % boundary % deviceData,&
                                                   vector % boundary % deviceData,&
                                                   vector % N, &
                                                   vector % nVar, &
                                                   vector % nelem )
#else
      msg = "GPU Acceleration is not currently enabled in SELF."
      WARNING(msg)
#endif
    ELSE
      DO iel = 1,vector % nelem
        DO ivar = 1,vector % nvar
          DO j = 0,vector % N
            DO i = 0,vector % N
              DO row = 1,2
                jvar = row + 2*(ivar-1)
                scalar % interior % hostData(i,j,jvar,iel) = vector % interior % hostData(row,i,j,ivar,iel)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      ! Boundary Terms
      DO iel = 1,vector % nelem
        DO iside = 1,4
          DO ivar = 1,vector % nvar
            DO j = 0,vector % N
              DO row = 1,2
                jvar = row + 2*(ivar-1)
                scalar % boundary % hostData(j,jvar,iside,iel) = vector % boundary % hostData(row,j,ivar,iside,iel)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE MapToScalar_MappedVector2D

  SUBROUTINE MapToTensor_MappedVector2D(vector,tensor,gpuAccel)
#undef __FUNC__
#define __FUNC__ "MapToTensor_MappedVector2D"
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(in) :: vector
    TYPE(MappedTensor2D),INTENT(inout) :: tensor
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: row,col,i,j,ivar,jvar,iel,iside
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN
#ifdef GPU
      CALL MapToTensor_MappedVector2D_gpu_wrapper( tensor % interior % deviceData,&
                                                   vector % interior % deviceData,&
                                                   tensor % N, &
                                                   tensor % nVar, &
                                                   tensor % nelem )

      CALL MapToTensorBoundary_MappedVector2D_gpu_wrapper( tensor % boundary % deviceData,&
                                                   vector % boundary % deviceData,&
                                                   tensor % N, &
                                                   tensor % nVar, &
                                                   tensor % nelem )
#else
      msg = "GPU Acceleration is not currently enabled in SELF."
      WARNING(msg)
#endif
    ELSE
      DO iel = 1,tensor % nelem
        DO ivar = 1,tensor % nvar
          DO j = 0,tensor % N
            DO i = 0,tensor % N
              DO col = 1,2
                DO row = 1,2
                  jvar = row + 2*(ivar-1)
                  tensor % interior % hostData(row,col,i,j,ivar,iel) = vector % interior % hostData(col,i,j,jvar,iel)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      ! Boundary Terms
      DO iel = 1,tensor % nelem
        DO iside = 1,4
          DO ivar = 1,tensor % nvar
            DO j = 0,tensor % N
              DO col = 1,2
                DO row = 1,2
                  jvar = row + 2*(ivar-1)
                  tensor % boundary % hostData(row,col,j,ivar,iside,iel) = vector % boundary % hostData(col,j,jvar,iside,iel)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF
  END SUBROUTINE MapToTensor_MappedVector2D

  SUBROUTINE ContravariantProjection_MappedVector2D(physVector,geometry,compVector,gpuAccel)
#undef __FUNC__
#define __FUNC__ "ContravariantProjection_MappedVector2D"
    ! Takes a vector that has physical space coordinate directions (x,y,z) and projects the vector
    ! into the the contravariant basis vector directions. Keep in mind that the contravariant basis
    ! vectors are really the Jacobian weighted contravariant basis vectors
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(in) :: physVector
    TYPE(SEMQuad),INTENT(in) :: geometry
    LOGICAL,INTENT(in) :: gpuAccel
    TYPE(MappedVector2D),INTENT(inout) :: compVector
    ! Local
    INTEGER :: i,j,ivar,iel,iside
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN
#ifdef GPU
      CALL ContravariantProjection_MappedVector2D_gpu_wrapper(physVector % interior % deviceData, &
                                                              compVector % interior % deviceData, &
                                                              geometry % dsdx % interior % deviceData, &
                                                              physVector % N, &
                                                              physVector % nVar, &
                                                              physVector % nElem)

      CALL ContravariantProjectionBoundary_MappedVector2D_gpu_wrapper(physVector % boundary % deviceData, &
                                                              compVector % boundary % deviceData, &
                                                              geometry % dsdx % boundary % deviceData, &
                                                              physVector % N, &
                                                              physVector % nVar, &
                                                              physVector % nElem)
#else
      msg = "GPU Acceleration is not currently enabled in SELF."
      WARNING(msg)
#endif

    ELSE
      ! Assume that tensor(j,i) is vector i, component j
      ! => dot product is done along first dimension
      ! to project onto computational space
      DO iel = 1,physVector % nElem
        DO ivar = 1,physVector % nVar
          DO j = 0,physVector % N
            DO i = 0,physVector % N

              compVector % interior % hostData(1,i,j,ivar,iel) = &
                geometry % dsdx % interior % hostData(1,1,i,j,1,iel)* &
                physVector % interior % hostData(1,i,j,ivar,iel) + &
                geometry % dsdx % interior % hostData(2,1,i,j,1,iel)* &
                physVector % interior % hostData(2,i,j,ivar,iel)

              compVector % interior % hostData(2,i,j,ivar,iel) = &
                geometry % dsdx % interior % hostData(1,2,i,j,1,iel)* &
                physVector % interior % hostData(1,i,j,ivar,iel) + &
                geometry % dsdx % interior % hostData(2,2,i,j,1,iel)* &
                physVector % interior % hostData(2,i,j,ivar,iel)

            END DO
          END DO
        END DO
      END DO

      ! Boundary Terms
      DO iel = 1,physVector % nElem
        DO iside = 1,4
          DO ivar = 1,physVector % nVar
            DO j = 0,physVector % N
              compVector % boundary % hostData(1,j,ivar,iside,iel) = &
                geometry % dsdx % boundary % hostData(1,1,j,1,iside,iel)* &
                physVector % boundary % hostData(1,j,ivar,iside,iel) + &
                geometry % dsdx % boundary % hostData(2,1,j,1,iside,iel)* &
                physVector % boundary % hostData(2,j,ivar,iside,iel)

              compVector % boundary % hostData(2,j,ivar,iside,iel) = &
                geometry % dsdx % boundary % hostData(1,2,j,1,iside,iel)* &
                physVector % boundary % hostData(1,j,ivar,iside,iel) + &
                geometry % dsdx % boundary % hostData(2,2,j,1,iside,iel)* &
                physVector % boundary % hostData(2,j,ivar,iside,iel)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    END IF

  END SUBROUTINE ContravariantProjection_MappedVector2D

  SUBROUTINE JacobianWeight_MappedVector2D(vector,geometry,gpuAccel)
#undef __FUNC__
#define __FUNC__ "JacobianWeight_MappedVector2D"
  ! Applies the inverse jacobian
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(inout) :: vector
    TYPE(SEMQuad),INTENT(in) :: geometry
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iEl,iVar,i,j
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN

#ifdef GPU
      CALL JacobianWeight_MappedVector2D_gpu_wrapper(vector % interior % deviceData, &
                                                     geometry % J % interior % deviceData, &
                                                     vector % N, &
                                                     vector % nVar, &
                                                     vector % nElem)
#else
      msg = "GPU Acceleration is not currently enabled in SELF."
      WARNING(msg)
#endif

    ELSE

      DO iEl = 1,vector % nElem
        DO iVar = 1,vector % nVar
          DO j = 0,vector % N
            DO i = 0,vector % N
              vector % interior % hostData(1,i,j,iVar,iEl) = vector % interior % hostData(1,i,j,iVar,iEl)/&
                                                             geometry % J % interior % hostData(i,j,1,iEl)
              vector % interior % hostData(2,i,j,iVar,iEl) = vector % interior % hostData(2,i,j,iVar,iEl)/&
                                                             geometry % J % interior % hostData(i,j,1,iEl)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE JacobianWeight_MappedVector2D

! SideExchange_MappedVector3D is used to populate vector % extBoundary
  ! by finding neighboring elements that share a side and copying the neighboring
  ! elements solution % boundary data.

  SUBROUTINE SideExchange_MappedVector3D(vector,mesh,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(inout) :: vector
    TYPE(Mesh3D),INTENT(in) :: mesh
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1, e2, s1, s2, sid 
    INTEGER :: flip, bcid, globalSideId
    INTEGER :: i1, i2, j1, j2, ivar

    IF(gpuAccel)THEN

      ! TO DO ! 
      PRINT*, 'Woopsie! GPU Acceleration not implemented yet for SideExchange'

    ELSE

      DO e1 = 1, mesh % nElem
        s1 = 0
        DO sid = mesh % elemInfo % hostData(3,e1)+1, mesh % elemInfo % hostData(4,e1) ! Loop over local sides 
          s1 = s1 + 1 ! Increment local side ID 
          globalSideId = mesh % sideInfo % hostData(2,sid)
          e2 = mesh % sideInfo % hostData(3,sid)
          s2 = mesh % sideInfo % hostData(4,sid)/10
          flip = mesh % sideInfo % hostData(4,sid)-s2*10
          bcid = mesh % sideInfo % hostData(5,sid)

          IF(bcid /= 0)THEN   

            IF(flip == 1)THEN 
          
              DO ivar = 1, vector % nvar
                DO j1 = 0, vector % N
                  DO i1 = 0, vector % N
                    vector % extBoundary % hostData(1:3,i1,j1,ivar,s1,e1) = &
                      vector % boundary % hostData(1:3,i1,j1,ivar,s2,e2)
                  ENDDO
                ENDDO
              ENDDO

            ELSEIF(flip == 2)THEN

              DO ivar = 1, vector % nvar
                DO j1 = 0, vector % N
                  DO i1 = 0, vector % N
                    i2 = vector % N - j1
                    j2 = i1
                    vector % extBoundary % hostData(1:3,i1,j1,ivar,s1,e1) = &
                      vector % boundary % hostData(1:3,i2,j2,ivar,s2,e2)
                  ENDDO
                ENDDO
              ENDDO

            ELSEIF(flip == 3)THEN
                    
              DO ivar = 1, vector % nvar
                DO j1 = 0, vector % N
                  DO i1 = 0, vector % N
                    i2 = vector % N - i1
                    j2 = vector % N - j1 
                    vector % extBoundary % hostData(1:3,i1,j1,ivar,s1,e1) = &
                      vector % boundary % hostData(1:3,i2,j2,ivar,s2,e2)
                  ENDDO
                ENDDO
              ENDDO
          
            ELSEIF(flip == 4)THEN
                    
              DO ivar = 1, vector % nvar
                DO j1 = 0, vector % N
                  DO i1 = 0, vector % N
                    i2 = j1
                    j2 = vector % N - i1
                    vector % extBoundary % hostData(1:3,i1,j1,ivar,s1,e1) = &
                      vector % boundary % hostData(1:3,i2,j2,ivar,s2,e2)
                  ENDDO
                ENDDO
              ENDDO
          
            ENDIF

          ENDIF

        ENDDO
      ENDDO

    END IF
    
  END SUBROUTINE SideExchange_MappedVector3D

  SUBROUTINE BassiRebaySides_MappedVector3D(vector,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(inout) :: vector
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iel
    INTEGER :: iside
    INTEGER :: ivar
    INTEGER :: i, j

    IF(gpuAccel)THEN

      ! TO DO ! 
      PRINT*, 'Woopsie! GPU Acceleration not implemented yet for BassiRebay'

    ELSE

      DO iel = 1, vector % nElem
        DO iside = 1, 6
          DO ivar = 1, vector % nVar
            DO j = 0, vector % N
              DO i = 0, vector % N
                vector % boundary % hostData(1:3,i,j,ivar,iside,iel) = 0.5_prec*(&
                  vector % boundary % hostData(1:3,i,j,ivar,iside,iel) +&
                  vector % extBoundary % hostData(1:3,i,j,ivar,iside,iel))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF
           
  END SUBROUTINE BassiRebaySides_MappedVector3D

  SUBROUTINE Divergence_MappedVector3D(physVector,compVector,geometry,divVector,dForm,gpuAccel)
    !
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(in) :: physVector
    TYPE(MappedVector3D),INTENT(inout) :: compVector
    TYPE(SEMHex),INTENT(in) :: geometry
    TYPE(MappedScalar3D),INTENT(inout) :: divVector
    INTEGER,INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

    CALL physVector % ContravariantProjection(geometry,compVector,gpuAccel)

    IF (dForm == selfWeakDGForm) THEN

      IF (gpuAccel) THEN
        CALL compVector % interp % VectorDGDivergence_3D(compVector % interior % deviceData, &
                                                         compVector % boundary % deviceData, &
                                                         divVector % interior % deviceData, &
                                                         compVector % nvar, &
                                                         compVector % nelem)
      ELSE
        CALL compVector % interp % VectorDGDivergence_3D(compVector % interior % hostData, &
                                                         compVector % boundary % hostData, &
                                                         divVector % interior % hostData, &
                                                         compVector % nvar, &
                                                         compVector % nelem)
      ENDIF

    ELSE IF (dForm == selfStrongForm) THEN

      IF (gpuAccel) THEN
        CALL compVector % interp % VectorDivergence_3D(compVector % interior % deviceData, &
                                                       divVector % interior % deviceData, &
                                                       compVector % nvar, &
                                                       compVector % nelem)
      ELSE
        CALL compVector % interp % VectorDivergence_3D(compVector % interior % hostData, &
                                                       divVector % interior % hostData, &
                                                       compVector % nvar, &
                                                       compVector % nelem)
      ENDIF

    END IF

    CALL divVector % JacobianWeight(geometry,gpuAccel)

  END SUBROUTINE Divergence_MappedVector3D

  SUBROUTINE Gradient_MappedVector3D(vector,workScalar,workVector,workTensor,geometry,gradF,dForm,gpuAccel)
    ! Strong Form Operator - (Conservative Form)
    !
    ! Calculates the gradient of a scalar 3D function using the conservative form of the
    ! mapped gradient operator
    !
    ! \grad_{phys}( f ) =  (1 / J)*(\partial / \partial \xi_i ( J\vec{a}_i f )
    !
    ! where the sum over i is implied.
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(in) :: vector
    TYPE(MappedScalar3D),INTENT(inout) :: workScalar ! (scalar) nvar = 3*nvar
    TYPE(MappedVector3D),INTENT(inout) :: workVector ! (vector) nvar = 3*nvar
    TYPE(MappedTensor3D),INTENT(inout) :: workTensor ! (tensor) nvar = 3*nvar
    TYPE(SEMHex),INTENT(in) :: geometry
    TYPE(MappedTensor3D),INTENT(inout) :: gradF
    INTEGER,INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

    CALL vector % MapToScalar(workScalar,gpuAccel)

    CALL workScalar % ContravariantWeight(geometry,workTensor,gpuAccel)

    IF (dForm == selfWeakDGForm) THEN

      IF (gpuAccel) THEN
        CALL workTensor % interp % TensorDGDivergence_3D(workTensor % interior % deviceData, &
                                                         workTensor % boundary % deviceData, &
                                                         workVector % interior % deviceData, &
                                                         workTensor % nVar, &
                                                         workTensor % nElem)
      ELSE
        CALL workTensor % interp % TensorDGDivergence_3D(workTensor % interior % hostData, &
                                                         workTensor % boundary % hostData, &
                                                         workVector % interior % hostData, &
                                                         workTensor % nVar, &
                                                         workTensor % nElem)
      END IF

    ELSE IF (dForm == selfStrongForm) THEN

      IF (gpuAccel) THEN
        CALL workTensor % interp % TensorDivergence_3D(workTensor % interior % deviceData, &
                                                       workVector % interior % deviceData, &
                                                       workTensor % nVar, &
                                                       workTensor % nElem)
      ELSE
        CALL workTensor % interp % TensorDivergence_3D(workTensor % interior % hostData, &
                                                       workVector % interior % hostData, &
                                                       workTensor % nVar, &
                                                       workTensor % nElem)
      END IF

    END IF

    CALL workVector % MapToTensor(gradF,gpuAccel)

    CALL gradF % JacobianWeight(geometry,gpuAccel)

  END SUBROUTINE Gradient_MappedVector3D

  SUBROUTINE MapToScalar_MappedVector3D(vector,scalar,gpuAccel)
#undef __FUNC__
#define __FUNC__ "MapToScalar_MappedVector3D"
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(in) :: vector
    TYPE(MappedScalar3D),INTENT(inout) :: scalar
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: row,i,j,k,ivar,jvar,iel,iside
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN
#ifdef GPU
      CALL MapToScalar_MappedVector3D_gpu_wrapper( scalar % interior % deviceData,&
                                                   vector % interior % deviceData,&
                                                   vector % N, &
                                                   vector % nVar, &
                                                   vector % nelem )

      CALL MapToScalarBoundary_MappedVector3D_gpu_wrapper( scalar % boundary % deviceData,&
                                                   vector % boundary % deviceData,&
                                                   vector % N, &
                                                   vector % nVar, &
                                                   vector % nelem )
#else
      msg = "GPU Acceleration is not currently enabled in SELF."
      WARNING(msg)
#endif
    ELSE
      DO iel = 1,vector % nelem
        DO ivar = 1,vector % nvar
          DO k = 0,vector % N
            DO j = 0,vector % N
              DO i = 0,vector % N
                DO row = 1,3
                  jvar = row + 3*(ivar-1)
                  scalar % interior % hostData(i,j,k,jvar,iel) = vector % interior % hostData(row,i,j,k,ivar,iel)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      ! Boundary Terms
      DO iel = 1,vector % nelem
        DO iside = 1,6
          DO ivar = 1,vector % nvar
            DO k = 0,vector % N
              DO j = 0,vector % N
                DO row = 1,3
                  jvar = row + 3*(ivar-1)
                  scalar % boundary % hostData(j,k,jvar,iside,iel) = vector % boundary % hostData(row,j,k,ivar,iside,iel)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF
  END SUBROUTINE MapToScalar_MappedVector3D

  SUBROUTINE MapToTensor_MappedVector3D(vector,tensor,gpuAccel)
#undef __FUNC__
#define __FUNC__ "MapToTensor_MappedVector3D"
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(in) :: vector
    TYPE(MappedTensor3D),INTENT(inout) :: tensor
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: row,col,i,j,k,ivar,jvar,iel,iside
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN
#ifdef GPU
      CALL MapToTensor_MappedVector3D_gpu_wrapper( tensor % interior % deviceData,&
                                                   vector % interior % deviceData,&
                                                   tensor % N, &
                                                   tensor % nVar, &
                                                   tensor % nelem )

      CALL MapToTensorBoundary_MappedVector3D_gpu_wrapper( tensor % boundary % deviceData,&
                                                   vector % boundary % deviceData,&
                                                   tensor % N, &
                                                   tensor % nVar, &
                                                   tensor % nelem )
#else
      msg = "GPU Acceleration is not currently enabled in SELF."
      WARNING(msg)
#endif
    ELSE
      DO iel = 1,tensor % nelem
        DO ivar = 1,tensor % nvar
          DO k = 0,tensor % N
            DO j = 0,tensor % N
              DO i = 0,tensor % N
                DO col = 1,3
                  DO row = 1,3
                    jvar = row + 3*(ivar-1)
                    tensor % interior % hostData(row,col,i,j,k,ivar,iel) = vector % interior % hostData(col,i,j,k,jvar,iel)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      ! Boundary Terms
      DO iel = 1,tensor % nelem
        DO iside = 1,6
          DO ivar = 1,tensor % nvar
            DO k = 0,tensor % N
              DO j = 0,tensor % N
                DO col = 1,3
                  DO row = 1,3
                    jvar = row + 3*(ivar-1)
                    tensor % boundary % hostData(row,col,j,k,ivar,iside,iel) = vector % boundary % hostData(col,j,k,jvar,iside,iel)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF
  END SUBROUTINE MapToTensor_MappedVector3D

  SUBROUTINE ContravariantProjection_MappedVector3D(physVector,geometry,compVector,gpuAccel)
#undef __FUNC__
#define __FUNC__ "ContravariantProjection_MappedVector3D"
    ! Takes a vector that has physical space coordinate directions (x,y,z) and projects the vector
    ! into the the contravariant basis vector directions. Keep in mind that the contravariant basis
    ! vectors are really the Jacobian weighted contravariant basis vectors
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(in) :: physVector
    TYPE(MappedVector3D),INTENT(inout) :: compVector
    TYPE(SEMHex),INTENT(in) :: geometry
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: i,j,k,iVar,iEl,iDir,iside
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN
#ifdef GPU
      CALL ContravariantProjection_MappedVector3D_gpu_wrapper(physVector % interior % deviceData, &
                                                              compVector % interior % deviceData, &
                                                              geometry % dsdx % interior % deviceData, &
                                                              physVector % N, &
                                                              physVector % nVar, &
                                                              physVector % nElem)

      CALL ContravariantProjectionBoundary_MappedVector3D_gpu_wrapper(physVector % boundary % deviceData, &
                                                              compVector % boundary % deviceData, &
                                                              geometry % dsdx % boundary % deviceData, &
                                                              physVector % N, &
                                                              physVector % nVar, &
                                                              physVector % nElem)
#else
     msg = "GPU Acceleration currently not enabled in SELF"
     WARNING(msg)
#endif

    ELSE
      ! Assume that tensor(j,i) is vector i, component j
      ! => dot product is done along first dimension to
      ! project onto computational space
      DO iEl = 1,physVector % nElem
        DO iVar = 1,physVector % nVar
          DO k = 0,physVector % N
            DO j = 0,physVector % N
              DO i = 0,physVector % N

                compVector % interior % hostData(1,i,j,k,iVar,iEl) = &
                        geometry % dsdx % interior % hostData(1,1,i,j,k,1,iEl)* &
                        physVector % interior % hostData(1,i,j,k,iVar,iEl) + &
                        geometry % dsdx % interior % hostData(2,1,i,j,k,1,iEl)* &
                        physVector % interior % hostData(2,i,j,k,iVar,iEl) + &
                        geometry % dsdx % interior % hostData(3,1,i,j,k,1,iEl)* &
                        physVector % interior % hostData(3,i,j,k,iVar,iEl)

                compVector % interior % hostData(2,i,j,k,iVar,iEl) = &
                        geometry % dsdx % interior % hostData(1,2,i,j,k,1,iEl)* &
                        physVector % interior % hostData(1,i,j,k,iVar,iEl) + &
                        geometry % dsdx % interior % hostData(2,2,i,j,k,1,iEl)* &
                        physVector % interior % hostData(2,i,j,k,iVar,iEl) + &
                        geometry % dsdx % interior % hostData(3,2,i,j,k,1,iEl)* &
                        physVector % interior % hostData(3,i,j,k,iVar,iEl)

                compVector % interior % hostData(3,i,j,k,iVar,iEl) = &
                        geometry % dsdx % interior % hostData(1,3,i,j,k,1,iEl)* &
                        physVector % interior % hostData(1,i,j,k,iVar,iEl) + &
                        geometry % dsdx % interior % hostData(2,3,i,j,k,1,iEl)* &
                        physVector % interior % hostData(2,i,j,k,iVar,iEl) + &
                        geometry % dsdx % interior % hostData(3,3,i,j,k,1,iEl)* &
                        physVector % interior % hostData(3,i,j,k,iVar,iEl)

              END DO
            END DO
          END DO
        END DO
      END DO

      ! Boundary Terms
      DO iEl = 1,physVector % nElem
        DO iside = 1,6
          DO iVar = 1,physVector % nVar
            DO k = 0,physVector % N
              DO j = 0,physVector % N
                compVector % boundary % hostData(1,j,k,iVar,iside,iEl) = &
                        geometry % dsdx % boundary % hostData(1,1,j,k,1,iside,iEl)* &
                        physVector % boundary % hostData(1,j,k,iVar,iside,iEl) + &
                        geometry % dsdx % boundary % hostData(2,1,j,k,1,iside,iEl)* &
                        physVector % boundary % hostData(2,j,k,iVar,iside,iEl) + &
                        geometry % dsdx % boundary % hostData(3,1,j,k,1,iside,iEl)* &
                        physVector % boundary % hostData(3,j,k,iVar,iside,iEl)

                compVector % boundary % hostData(2,j,k,iVar,iside,iEl) = &
                        geometry % dsdx % boundary % hostData(1,2,j,k,1,iside,iEl)* &
                        physVector % boundary % hostData(1,j,k,iVar,iside,iEl) + &
                        geometry % dsdx % boundary % hostData(2,2,j,k,1,iside,iEl)* &
                        physVector % boundary % hostData(2,j,k,iVar,iside,iEl) + &
                        geometry % dsdx % boundary % hostData(3,2,j,k,1,iside,iEl)* &
                        physVector % boundary % hostData(3,j,k,iVar,iside,iEl)

                compVector % boundary % hostData(3,j,k,iVar,iside,iEl) = &
                        geometry % dsdx % boundary % hostData(1,3,j,k,1,iside,iEl)* &
                        physVector % boundary % hostData(1,j,k,iVar,iside,iEl) + &
                        geometry % dsdx % boundary % hostData(2,3,j,k,1,iside,iEl)* &
                        physVector % boundary % hostData(2,j,k,iVar,iside,iEl) + &
                        geometry % dsdx % boundary % hostData(3,3,j,k,1,iside,iEl)* &
                        physVector % boundary % hostData(3,j,k,iVar,iside,iEl)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      
    END IF

  END SUBROUTINE ContravariantProjection_MappedVector3D

  SUBROUTINE JacobianWeight_MappedVector3D(vector,geometry,gpuAccel)
#undef __FUNC__
#define __FUNC__ "JacobianWeight_MappedVector3D"
  ! Applies the inverse jacobian
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(inout) :: vector
    TYPE(SEMHex),INTENT(in) :: geometry
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iEl,iVar,i,j,k
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN

#ifdef GPU
      CALL JacobianWeight_MappedVector3D_gpu_wrapper(vector % interior % deviceData, &
                                                     geometry % J % interior % deviceData, &
                                                     vector % N, &
                                                     vector % nVar, &
                                                     vector % nElem)
#else
      msg = "GPU Acceleration currently not enabled in SELF."
      WARNING(msg)
#endif

    ELSE

      DO iEl = 1,vector % nElem
        DO iVar = 1,vector % nVar
          DO k = 0,vector % N
            DO j = 0,vector % N
              DO i = 0,vector % N
                vector % interior % hostData(1,i,j,k,iVar,iEl) = vector % interior % hostData(1,i,j,k,iVar,iEl)/&
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
                vector % interior % hostData(2,i,j,k,iVar,iEl) = vector % interior % hostData(2,i,j,k,iVar,iEl)/&
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
                vector % interior % hostData(3,i,j,k,iVar,iEl) = vector % interior % hostData(3,i,j,k,iVar,iEl)/&
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE JacobianWeight_MappedVector3D

  ! ---------------------- Tensors ---------------------- !
  ! SideExchange_MappedTensor2D is used to populate tensor % extBoundary
  ! by finding neighboring elements that share a side and copying the neighboring
  ! elements solution % boundary data.

  SUBROUTINE SideExchange_MappedTensor2D(tensor,mesh,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedTensor2D),INTENT(inout) :: tensor
    TYPE(Mesh2D),INTENT(in) :: mesh
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1, e2, s1, s2, sid 
    INTEGER :: flip, bcid, globalSideId
    INTEGER :: i1, i2, ivar

    IF(gpuAccel)THEN

      ! TO DO ! 
      PRINT*, 'Woopsie! GPU Acceleration not implemented yet for SideExchange'

    ELSE

      DO e1 = 1, mesh % nElem
        s1 = 0
        DO sid = mesh % elemInfo % hostData(3,e1)+1, mesh % elemInfo % hostData(4,e1) ! Loop over local sides 
          s1 = s1 + 1 ! Increment local side ID 
          globalSideId = mesh % sideInfo % hostData(2,sid)
          e2 = mesh % sideInfo % hostData(3,sid)
          s2 = mesh % sideInfo % hostData(4,sid)/10
          flip = mesh % sideInfo % hostData(4,sid)-s2*10
          bcid = mesh % sideInfo % hostData(5,sid)

          IF(bcid /= 0)THEN   

            IF(flip == 1)THEN 
          
              DO ivar = 1, tensor % nvar
                DO i1 = 0, tensor % N
                  tensor % extBoundary % hostData(1:2,1:2,i1,ivar,s1,e1) = &
                      tensor % boundary % hostData(1:2,1:2,i1,ivar,s2,e2)
                ENDDO
              ENDDO

            ELSEIF(flip == 2)THEN

              DO ivar = 1, tensor % nvar
                DO i1 = 0, tensor % N
                  i2 = tensor % N - i1
                  tensor % extBoundary % hostData(1:2,1:2,i1,ivar,s1,e1) = &
                      tensor % boundary % hostData(1:2,1:2,i2,ivar,s2,e2)
                ENDDO
              ENDDO

            ENDIF

          ENDIF

        ENDDO
      ENDDO

    END IF
    
  END SUBROUTINE SideExchange_MappedTensor2D

  SUBROUTINE BassiRebaySides_MappedTensor2D(tensor,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedTensor2D),INTENT(inout) :: tensor
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iel
    INTEGER :: iside
    INTEGER :: ivar
    INTEGER :: i

    IF(gpuAccel)THEN

      ! TO DO ! 
      PRINT*, 'Woopsie! GPU Acceleration not implemented yet for BassiRebay'

    ELSE

      DO iel = 1, tensor % nElem
        DO iside = 1, 4
          DO ivar = 1, tensor % nVar
            DO i = 0, tensor % N
              tensor % boundary % hostData(1:2,1:2,i,ivar,iside,iel) = 0.5_prec*(&
                tensor % boundary % hostData(1:2,1:2,i,ivar,iside,iel) +&
                tensor % extBoundary % hostData(1:2,1:2,i,ivar,iside,iel))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF
           
  END SUBROUTINE BassiRebaySides_MappedTensor2D

  SUBROUTINE JacobianWeight_MappedTensor2D(tensor,geometry,gpuAccel)
#undef __FUNC__
#define __FUNC__ "JacobianWeight_MappedTensor2D"
  ! Applies the inverse jacobian
    IMPLICIT NONE
    CLASS(MappedTensor2D),INTENT(inout) :: tensor
    TYPE(SEMQuad),INTENT(in) :: geometry
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iEl,iVar,i,j
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN

#ifdef GPU
      CALL JacobianWeight_MappedTensor2D_gpu_wrapper(tensor % interior % deviceData, &
                                                     geometry % J % interior % deviceData, &
                                                     tensor % N, &
                                                     tensor % nVar, &
                                                     tensor % nElem)
#else
      msg = "GPU Acceleration currently not enabled in SELF."
      WARNING(msg)
#endif

    ELSE

      DO iEl = 1,tensor % nElem
        DO iVar = 1,tensor % nVar
          DO j = 0,tensor % N
            DO i = 0,tensor % N
              tensor % interior % hostData(1,1,i,j,iVar,iEl) = tensor % interior % hostData(1,1,i,j,iVar,iEl)/&
                                                             geometry % J % interior % hostData(i,j,1,iEl)
              tensor % interior % hostData(2,1,i,j,iVar,iEl) = tensor % interior % hostData(2,1,i,j,iVar,iEl)/&
                                                             geometry % J % interior % hostData(i,j,1,iEl)
              tensor % interior % hostData(1,2,i,j,iVar,iEl) = tensor % interior % hostData(1,2,i,j,iVar,iEl)/&
                                                             geometry % J % interior % hostData(i,j,1,iEl)
              tensor % interior % hostData(2,2,i,j,iVar,iEl) = tensor % interior % hostData(2,2,i,j,iVar,iEl)/&
                                                             geometry % J % interior % hostData(i,j,1,iEl)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE JacobianWeight_MappedTensor2D

  ! SideExchange_MappedVector3D is used to populate vector % extBoundary
  ! by finding neighboring elements that share a side and copying the neighboring
  ! elements solution % boundary data.

  SUBROUTINE SideExchange_MappedTensor3D(tensor,mesh,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedTensor3D),INTENT(inout) :: tensor
    TYPE(Mesh3D),INTENT(in) :: mesh
    LOGICAL, INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1, e2, s1, s2, sid 
    INTEGER :: flip, bcid, globalSideId
    INTEGER :: i1, i2, j1, j2, ivar

    IF(gpuAccel)THEN

      ! TO DO ! 
      PRINT*, 'Woopsie! GPU Acceleration not implemented yet for SideExchange'

    ELSE

      DO e1 = 1, mesh % nElem
        s1 = 0
        DO sid = mesh % elemInfo % hostData(3,e1)+1, mesh % elemInfo % hostData(4,e1) ! Loop over local sides 
          s1 = s1 + 1 ! Increment local side ID 
          globalSideId = mesh % sideInfo % hostData(2,sid)
          e2 = mesh % sideInfo % hostData(3,sid)
          s2 = mesh % sideInfo % hostData(4,sid)/10
          flip = mesh % sideInfo % hostData(4,sid)-s2*10
          bcid = mesh % sideInfo % hostData(5,sid)

          IF(bcid /= 0)THEN   

            IF(flip == 1)THEN 
          
              DO ivar = 1, tensor % nvar
                DO j1 = 0, tensor % N
                  DO i1 = 0, tensor % N
                    tensor % extBoundary % hostData(1:3,1:3,i1,j1,ivar,s1,e1) = &
                      tensor % boundary % hostData(1:3,1:3,i1,j1,ivar,s2,e2)
                  ENDDO
                ENDDO
              ENDDO

            ELSEIF(flip == 2)THEN

              DO ivar = 1, tensor % nvar
                DO j1 = 0, tensor % N
                  DO i1 = 0, tensor % N
                    i2 = tensor % N - j1
                    j2 = i1
                    tensor % extBoundary % hostData(1:3,1:3,i1,j1,ivar,s1,e1) = &
                      tensor % boundary % hostData(1:3,1:3,i2,j2,ivar,s2,e2)
                  ENDDO
                ENDDO
              ENDDO

            ELSEIF(flip == 3)THEN
                    
              DO ivar = 1, tensor % nvar
                DO j1 = 0, tensor % N
                  DO i1 = 0, tensor % N
                    i2 = tensor % N - i1
                    j2 = tensor % N - j1 
                    tensor % extBoundary % hostData(1:3,1:3,i1,j1,ivar,s1,e1) = &
                      tensor % boundary % hostData(1:3,1:3,i2,j2,ivar,s2,e2)
                  ENDDO
                ENDDO
              ENDDO
          
            ELSEIF(flip == 4)THEN
                    
              DO ivar = 1, tensor % nvar
                DO j1 = 0, tensor % N
                  DO i1 = 0, tensor % N
                    i2 = j1
                    j2 = tensor % N - i1
                    tensor % extBoundary % hostData(1:3,1:3,i1,j1,ivar,s1,e1) = &
                      tensor % boundary % hostData(1:3,1:3,i2,j2,ivar,s2,e2)
                  ENDDO
                ENDDO
              ENDDO
          
            ENDIF

          ENDIF

        ENDDO
      ENDDO

    END IF
    
  END SUBROUTINE SideExchange_MappedTensor3D

  SUBROUTINE BassiRebaySides_MappedTensor3D(tensor,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedTensor3D),INTENT(inout) :: tensor
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iel
    INTEGER :: iside
    INTEGER :: ivar
    INTEGER :: i, j

    IF(gpuAccel)THEN

      ! TO DO ! 
      PRINT*, 'Woopsie! GPU Acceleration not implemented yet for BassiRebay'

    ELSE

      DO iel = 1, tensor % nElem
        DO iside = 1, 6
          DO ivar = 1, tensor % nVar
            DO j = 0, tensor % N
              DO i = 0, tensor % N
                tensor % boundary % hostData(1:3,1:3,i,j,ivar,iside,iel) = 0.5_prec*(&
                  tensor % boundary % hostData(1:3,1:3,i,j,ivar,iside,iel) +&
                  tensor % extBoundary % hostData(1:3,1:3,i,j,ivar,iside,iel))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF
           
  END SUBROUTINE BassiRebaySides_MappedTensor3D

  SUBROUTINE JacobianWeight_MappedTensor3D(tensor,geometry,gpuAccel)
#undef __FUNC__
#define __FUNC__ "JacobianWeight_MappedTensor3D"
  ! Applies the inverse jacobian
    IMPLICIT NONE
    CLASS(MappedTensor3D),INTENT(inout) :: tensor
    TYPE(SEMHex),INTENT(in) :: geometry
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iEl,iVar,i,j,k
    CHARACTER(100) :: msg

    IF (gpuAccel) THEN

#ifdef GPU
      CALL JacobianWeight_MappedTensor3D_gpu_wrapper(tensor % interior % deviceData, &
                                                     geometry % J % interior % deviceData, &
                                                     tensor % N, &
                                                     tensor % nVar, &
                                                     tensor % nElem)
#else
      msg = "GPU Acceleration currently not enabled in SELF."
      WARNING(msg)
#endif

    ELSE

      DO iEl = 1,tensor % nElem
        DO iVar = 1,tensor % nVar
          DO k = 0,tensor % N
            DO j = 0,tensor % N
              DO i = 0,tensor % N
                tensor % interior % hostData(1,1,i,j,k,iVar,iEl) = tensor % interior % hostData(1,1,i,j,k,iVar,iEl)/&
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
                tensor % interior % hostData(2,1,i,j,k,iVar,iEl) = tensor % interior % hostData(2,1,i,j,k,iVar,iEl)/&
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
                tensor % interior % hostData(3,1,i,j,k,iVar,iEl) = tensor % interior % hostData(3,1,i,j,k,iVar,iEl)/&
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
                tensor % interior % hostData(1,2,i,j,k,iVar,iEl) = tensor % interior % hostData(1,2,i,j,k,iVar,iEl)/&
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
                tensor % interior % hostData(2,2,i,j,k,iVar,iEl) = tensor % interior % hostData(2,2,i,j,k,iVar,iEl)/&
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
                tensor % interior % hostData(3,2,i,j,k,iVar,iEl) = tensor % interior % hostData(3,2,i,j,k,iVar,iEl)/&
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
                tensor % interior % hostData(1,3,i,j,k,iVar,iEl) = tensor % interior % hostData(1,3,i,j,k,iVar,iEl)/&
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
                tensor % interior % hostData(2,3,i,j,k,iVar,iEl) = tensor % interior % hostData(2,3,i,j,k,iVar,iEl)/&
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
                tensor % interior % hostData(3,3,i,j,k,iVar,iEl) = tensor % interior % hostData(3,3,i,j,k,iVar,iEl)/&
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE JacobianWeight_MappedTensor3D

END MODULE SELF_MappedData
