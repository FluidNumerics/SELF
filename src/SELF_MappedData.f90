! SELF_MappedData.F90
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE SELF_MappedData

  USE SELF_Constants
  USE SELF_Lagrange
  USE SELF_Data
  USE SELF_Mesh
  USE SELF_Geometry
  USE SELF_HDF5
  USE HDF5

  USE FEQParse

  USE ISO_C_BINDING

  IMPLICIT NONE

#include "SELF_Macros.h"

  TYPE,EXTENDS(Scalar1D),PUBLIC :: MappedScalar1D

  CONTAINS
    PROCEDURE,PUBLIC :: SideExchange => SideExchange_MappedScalar1D
    PROCEDURE,PUBLIC :: BassiRebaySides => BassiRebaySides_MappedScalar1D
    GENERIC,PUBLIC :: Derivative => Derivative_MappedScalar1D
    PROCEDURE,PRIVATE :: Derivative_MappedScalar1D
    PROCEDURE,PUBLIC :: JacobianWeight => JacobianWeight_MappedScalar1D

    PROCEDURE,PUBLIC :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar1D

  END TYPE MappedScalar1D

  TYPE,EXTENDS(Scalar2D),PUBLIC :: MappedScalar2D

  CONTAINS

    PROCEDURE,PUBLIC :: SideExchange => SideExchange_MappedScalar2D
    PROCEDURE,PUBLIC :: BassiRebaySides => BassiRebaySides_MappedScalar2D

    GENERIC,PUBLIC :: Gradient => Gradient_MappedScalar2D
    PROCEDURE,PRIVATE :: Gradient_MappedScalar2D
    PROCEDURE,PRIVATE :: GradientSF_MappedScalar2D ! Strong-Form Gradient
    PROCEDURE,PRIVATE :: GradientBR_MappedScalar2D ! Bassi-Rebay Gradient

    PROCEDURE,PUBLIC :: JacobianWeight => JacobianWeight_MappedScalar2D

    PROCEDURE,PRIVATE :: MPIExchangeAsync => MPIExchangeAsync_MappedScalar2D
    PROCEDURE,PRIVATE :: ApplyFlip => ApplyFlip_MappedScalar2D

    PROCEDURE,PUBLIC :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar2D

    PROCEDURE,PUBLIC :: Integral => Integral_MappedScalar2D

  END TYPE MappedScalar2D

  TYPE,EXTENDS(Scalar3D),PUBLIC :: MappedScalar3D

  CONTAINS

    PROCEDURE,PUBLIC :: SideExchange => SideExchange_MappedScalar3D
    PROCEDURE,PUBLIC :: BassiRebaySides => BassiRebaySides_MappedScalar3D
    
    GENERIC,PUBLIC :: Gradient => Gradient_MappedScalar3D
    PROCEDURE,PRIVATE :: Gradient_MappedScalar3D
    PROCEDURE,PRIVATE :: GradientSF_MappedScalar3D ! Strong-Form Gradient
    PROCEDURE,PRIVATE :: GradientBR_MappedScalar3D ! Bassi-Rebay Gradient

    PROCEDURE,PUBLIC :: JacobianWeight => JacobianWeight_MappedScalar3D

    PROCEDURE,PRIVATE :: MPIExchangeAsync => MPIExchangeAsync_MappedScalar3D
    PROCEDURE,PRIVATE :: ApplyFlip => ApplyFlip_MappedScalar3D

    PROCEDURE,PUBLIC :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar3D

  END TYPE MappedScalar3D

  TYPE,EXTENDS(Vector2D),PUBLIC :: MappedVector2D

  CONTAINS

    PROCEDURE,PUBLIC :: SideExchange => SideExchange_MappedVector2D
    PROCEDURE,PUBLIC :: BassiRebaySides => BassiRebaySides_MappedVector2D

    GENERIC,PUBLIC :: Divergence => Divergence_MappedVector2D

    PROCEDURE,PRIVATE :: Divergence_MappedVector2D
    PROCEDURE,PUBLIC :: ContravariantProjection => ContravariantProjection_MappedVector2D
    PROCEDURE,PUBLIC :: JacobianWeight => JacobianWeight_MappedVector2D

    PROCEDURE,PRIVATE :: MPIExchangeAsync => MPIExchangeAsync_MappedVector2D
    PROCEDURE,PRIVATE :: ApplyFlip => ApplyFlip_MappedVector2D

    PROCEDURE,PUBLIC :: SetInteriorFromEquation => SetInteriorFromEquation_MappedVector2D


  END TYPE MappedVector2D

  TYPE,EXTENDS(Vector3D),PUBLIC :: MappedVector3D

  CONTAINS

    PROCEDURE,PUBLIC :: SideExchange => SideExchange_MappedVector3D
    PROCEDURE,PUBLIC :: BassiRebaySides => BassiRebaySides_MappedVector3D

    GENERIC,PUBLIC :: Divergence => Divergence_MappedVector3D
    PROCEDURE,PRIVATE :: Divergence_MappedVector3D

    PROCEDURE,PUBLIC :: ContravariantProjection => ContravariantProjection_MappedVector3D
    PROCEDURE,PUBLIC :: JacobianWeight => JacobianWeight_MappedVector3D

    PROCEDURE,PRIVATE :: MPIExchangeAsync => MPIExchangeAsync_MappedVector3D
    PROCEDURE,PRIVATE :: ApplyFlip => ApplyFlip_MappedVector3D

    PROCEDURE,PUBLIC :: SetInteriorFromEquation => SetInteriorFromEquation_MappedVector3D

  END TYPE MappedVector3D

  ! TYPE,EXTENDS(P2Vector2D),PUBLIC :: MappedP2Vector2D

  ! CONTAINS

  !   PROCEDURE,PUBLIC :: SideExchange => SideExchange_MappedP2Vector2D

  !   GENERIC,PUBLIC :: Divergence => Divergence_MappedP2Vector2D

  !   PROCEDURE,PRIVATE :: Divergence_MappedP2Vector2D
  !   PROCEDURE,PUBLIC :: ContravariantProjection => ContravariantProjection_MappedP2Vector2D
  !   !PROCEDURE,PUBLIC :: JacobianWeight => JacobianWeight_MappedP2Vector2D

  !   PROCEDURE,PRIVATE :: MPIExchangeAsync => MPIExchangeAsync_MappedP2Vector2D
  !   PROCEDURE,PRIVATE :: ApplyFlip => ApplyFlip_MappedP2Vector2D

  ! END TYPE MappedP2Vector2D

  INTERFACE
    SUBROUTINE GradientBR_MappedScalar2D_gpu_wrapper(scalar,avgBoundary,dsdx,jacobian,nHat,nScale,&
                    gradF,dgMatrix,bMatrix,qWeights,N,nVar,nEl) &
      bind(c,name="GradientBR_MappedScalar2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,avgBoundary,dsdx,jacobian,nHat,nScale,gradF,dgMatrix,bMatrix,qWeights
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE GradientBR_MappedScalar2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE GradientSF_MappedScalar2D_gpu_wrapper(scalar,dsdx,jacobian,gradF,dMatrix,N,nVar,nEl) &
      bind(c,name="GradientSF_MappedScalar2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,dsdx,jacobian,gradF,dMatrix
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE GradientSF_MappedScalar2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE GradientBR_MappedScalar3D_gpu_wrapper(scalar,avgBoundary,dsdx,jacobian,nHat,nScale,&
                    gradF,dgMatrix,bMatrix,qWeights,N,nVar,nEl) &
      bind(c,name="GradientBR_MappedScalar3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,avgBoundary,dsdx,jacobian,nHat,nScale,gradF,dgMatrix,bMatrix,qWeights
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE GradientBR_MappedScalar3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE GradientSF_MappedScalar3D_gpu_wrapper(scalar,dsdx,jacobian,gradF,dMatrix,N,nVar,nEl) &
      bind(c,name="GradientSF_MappedScalar3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,dsdx,jacobian,gradF,dMatrix
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE GradientSF_MappedScalar3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE JacobianWeight_MappedScalar1D_gpu_wrapper(scalar,dxds,N,nVar,nEl) &
      bind(c,name="JacobianWeight_MappedScalar1D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,dxds
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE JacobianWeight_MappedScalar1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE JacobianWeight_MappedScalar2D_gpu_wrapper(scalar,jacobian,N,nVar,nEl) &
      bind(c,name="JacobianWeight_MappedScalar2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,jacobian
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE JacobianWeight_MappedScalar2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE JacobianWeight_MappedScalar3D_gpu_wrapper(scalar,jacobian,N,nVar,nEl) &
      bind(c,name="JacobianWeight_MappedScalar3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: scalar,jacobian
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE JacobianWeight_MappedScalar3D_gpu_wrapper
  END INTERFACE

!  INTERFACE
!    SUBROUTINE ContravariantWeight_MappedScalar2D_gpu_wrapper(scalar,workTensor,dsdx,N,nVar,nEl) &
!      bind(c,name="ContravariantWeight_MappedScalar2D_gpu_wrapper")
!      USE ISO_C_BINDING
!      IMPLICIT NONE
!      TYPE(c_ptr) :: scalar,workTensor,dsdx
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!    END SUBROUTINE ContravariantWeight_MappedScalar2D_gpu_wrapper
!  END INTERFACE

!  INTERFACE
!    SUBROUTINE ContravariantWeightBoundary_MappedScalar2D_gpu_wrapper(scalar,workTensor,dsdx,N,nVar,nEl) &
!      bind(c,name="ContravariantWeightBoundary_MappedScalar2D_gpu_wrapper")
!      USE ISO_C_BINDING
!      IMPLICIT NONE
!      TYPE(c_ptr) :: scalar,workTensor,dsdx
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!    END SUBROUTINE ContravariantWeightBoundary_MappedScalar2D_gpu_wrapper
!  END INTERFACE
!
!  INTERFACE
!    SUBROUTINE ContravariantWeight_MappedScalar3D_gpu_wrapper(scalar,workTensor,dsdx,N,nVar,nEl) &
!      bind(c,name="ContravariantWeight_MappedScalar3D_gpu_wrapper")
!      USE ISO_C_BINDING
!      IMPLICIT NONE
!      TYPE(c_ptr) :: scalar,workTensor,dsdx
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!    END SUBROUTINE ContravariantWeight_MappedScalar3D_gpu_wrapper
!  END INTERFACE
!
!  INTERFACE
!    SUBROUTINE ContravariantWeightBoundary_MappedScalar3D_gpu_wrapper(scalar,workTensor,dsdx,N,nVar,nEl) &
!      bind(c,name="ContravariantWeightBoundary_MappedScalar3D_gpu_wrapper")
!      USE ISO_C_BINDING
!      IMPLICIT NONE
!      TYPE(c_ptr) :: scalar,workTensor,dsdx
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!    END SUBROUTINE ContravariantWeightBoundary_MappedScalar3D_gpu_wrapper
!  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantProjection_MappedVector2D_gpu_wrapper(vector,dsdx,N,nVar,nEl) &
      bind(c,name="ContravariantProjection_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: vector,dsdx
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ContravariantProjection_MappedVector2D_gpu_wrapper
  END INTERFACE

  ! INTERFACE
  !   SUBROUTINE ContravariantProjection_MappedP2Vector2D_gpu_wrapper(vector,physical,dsdx,N,nVar,nEl) &
  !     bind(c,name="ContravariantProjection_MappedP2Vector2D_gpu_wrapper")
  !     USE ISO_C_BINDING
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: vector,physical,dsdx
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE ContravariantProjection_MappedP2Vector2D_gpu_wrapper
  ! END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantProjectionBoundary_MappedVector2D_gpu_wrapper(physVector,compVector,dsdx,N,nVar,nEl) &
      bind(c,name="ContravariantProjectionBoundary_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: physVector,compVector,dsdx
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ContravariantProjectionBoundary_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE JacobianWeight_MappedVector2D_gpu_wrapper(vector,jacobian,N,nVar,nEl) &
      bind(c,name="JacobianWeight_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: vector,jacobian
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE JacobianWeight_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantProjection_MappedVector3D_gpu_wrapper(vector,dsdx,N,nVar,nEl) &
      bind(c,name="ContravariantProjection_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: vector,dsdx
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ContravariantProjection_MappedVector3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ContravariantProjectionBoundary_MappedVector3D_gpu_wrapper(physVector,compVector,dsdx,N,nVar,nEl) &
      bind(c,name="ContravariantProjectionBoundary_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: physVector,compVector,dsdx
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ContravariantProjectionBoundary_MappedVector3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE JacobianWeight_MappedVector3D_gpu_wrapper(vector,jacobian,N,nVar,nEl) &
      bind(c,name="JacobianWeight_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: vector,jacobian
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE JacobianWeight_MappedVector3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE SideExchange_MappedScalar2D_gpu_wrapper(extBoundary,boundary, &
                                                       sideInfo,elemToRank,rankId,offset,N,nVar,nEl) &
      bind(c,name="SideExchange_MappedScalar2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: extBoundary,boundary,sideInfo,elemToRank
      INTEGER(C_INT),VALUE :: rankId,offset,N,nVar,nEl
    END SUBROUTINE SideExchange_MappedScalar2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE SideExchange_MappedVector2D_gpu_wrapper(extBoundary,boundary, &
                                                       sideInfo,elemToRank,rankId,offset,N,nVar,nEl) &
      bind(c,name="SideExchange_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: extBoundary,boundary,sideInfo,elemToRank
      INTEGER(C_INT),VALUE :: rankId,offset,N,nVar,nEl
    END SUBROUTINE SideExchange_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE SideExchange_MappedScalar3D_gpu_wrapper(extBoundary,boundary, &
                                                       sideInfo,elemToRank,rankId,offset,N,nVar,nEl) &
      bind(c,name="SideExchange_MappedScalar3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: extBoundary,boundary,sideInfo,elemToRank
      INTEGER(C_INT),VALUE :: rankId,offset,N,nVar,nEl
    END SUBROUTINE SideExchange_MappedScalar3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE SideExchange_MappedVector3D_gpu_wrapper(extBoundary,boundary, &
                                                       sideInfo,elemToRank,rankId,offset,N,nVar,nEl) &
      bind(c,name="SideExchange_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: extBoundary,boundary,sideInfo,elemToRank
      INTEGER(C_INT),VALUE :: rankId,offset,N,nVar,nEl
    END SUBROUTINE SideExchange_MappedVector3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE BassiRebaySides_MappedScalar2D_gpu_wrapper(avgBoundary,boundary,extBoundary,N,nVar,nEl) &
      bind(c,name="BassiRebaySides_MappedScalar2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: extBoundary,boundary,avgBoundary
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE BassiRebaySides_MappedScalar2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE BassiRebaySides_MappedVector2D_gpu_wrapper(extBoundary,boundary,N,nVar,nEl) &
      bind(c,name="BassiRebaySides_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: extBoundary,boundary
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE BassiRebaySides_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE BassiRebaySides_MappedScalar3D_gpu_wrapper(avgBoundary,boundary,extBoundary,N,nVar,nEl) &
      bind(c,name="BassiRebaySides_MappedScalar3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: avgBoundary,extBoundary,boundary
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE BassiRebaySides_MappedScalar3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE BassiRebaySides_MappedVector3D_gpu_wrapper(extBoundary,boundary,N,nVar,nEl) &
      bind(c,name="BassiRebaySides_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: extBoundary,boundary
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE BassiRebaySides_MappedVector3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ApplyFlip_MappedScalar2D_gpu_wrapper(extBoundary,selfSideInfo,elemToRank,rankId,N,nVar,nEl) &
      bind(c,name="ApplyFlip_MappedScalar2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: selfSideInfo,elemToRank,extBoundary
      INTEGER(C_INT),VALUE :: rankId,N,nVar,nEl
    END SUBROUTINE ApplyFlip_MappedScalar2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ApplyFlip_MappedVector2D_gpu_wrapper(extBoundary,selfSideInfo,elemToRank,rankId,N,nVar,nEl) &
      bind(c,name="ApplyFlip_MappedVector2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: selfSideInfo,elemToRank,extBoundary
      INTEGER(C_INT),VALUE :: rankId,N,nVar,nEl
    END SUBROUTINE ApplyFlip_MappedVector2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ApplyFlip_MappedScalar3D_gpu_wrapper(extBoundary,selfSideInfo,elemToRank,rankId,N,nVar,nEl) &
      bind(c,name="ApplyFlip_MappedScalar3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: selfSideInfo,elemToRank,extBoundary
      INTEGER(C_INT),VALUE :: rankId,N,nVar,nEl
    END SUBROUTINE ApplyFlip_MappedScalar3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ApplyFlip_MappedVector3D_gpu_wrapper(extBoundary,selfSideInfo,elemToRank,rankId,N,nVar,nEl) &
      bind(c,name="ApplyFlip_MappedVector3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: selfSideInfo,elemToRank,extBoundary
      INTEGER(C_INT),VALUE :: rankId,N,nVar,nEl
    END SUBROUTINE ApplyFlip_MappedVector3D_gpu_wrapper
  END INTERFACE

CONTAINS

! ---------------------- Scalars ---------------------- !

  SUBROUTINE SetInteriorFromEquation_MappedScalar1D( scalar, geometry, time )
    !!  Sets the scalar % interior attribute using the eqn attribute,
    !!  geometry (for physical positions), and provided simulation time. 
    IMPLICIT NONE
    CLASS(MappedScalar1D), INTENT(inout) :: scalar
    TYPE(Geometry1D), INTENT(in) :: geometry
    REAL(prec), INTENT(in) :: time
    ! Local
    INTEGER :: i, iEl, iVar
    REAL(prec) :: x

    ! TO DO : Check if scalar % eqn is set before proceeding

    DO iEl = 1,scalar % nElem
      DO iVar = 1, scalar % nVar
        DO i = 0, scalar % interp % N

          ! Get the mesh positions
          x = geometry % x % interior % hostData(i,1,iEl)

          scalar % interior % hostData(i,iVar,iEl) = &
            scalar % eqn(iVar) % Evaluate((/x, 0.0_prec, 0.0_prec, time/))

        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE SetInteriorFromEquation_MappedScalar1D

  SUBROUTINE SideExchange_MappedScalar1D(scalar,mesh,decomp,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedScalar1D),INTENT(inout) :: scalar
    TYPE(Mesh1D),INTENT(in) :: mesh
    TYPE(MPILayer),INTENT(inout) :: decomp
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1,e2,s1,s2,e2Global
    INTEGER :: flip,bcid
    INTEGER :: i1,i2,ivar
    INTEGER :: neighborRank
    INTEGER :: rankId, offset

      rankId = decomp % rankId
      offset = decomp % offsetElem % hostData(rankId)

      if( gpuAccel )then
        call scalar % boundary % updatehost()
      endif

      DO e1 = 1,mesh % nElem
        
        IF( e1 == 1 )THEN

          s1 = 2
          e2 = e1 + 1
          s2 = 1
          !neighborRank = decomp % elemToRank % hostData(e2Global)
          DO ivar = 1,scalar % nvar
            scalar % extBoundary % hostData(ivar,s1,e1) = scalar % boundary % hostData(ivar,s2,e2)
          ENDDO

        ELSEIF( e1 == mesh % nElem )THEN

          s1 = 1
          e2 = e1 - 1
          s2 = 2
          !neighborRank = decomp % elemToRank % hostData(e2Global)
          DO ivar = 1,scalar % nvar
            scalar % extBoundary % hostData(ivar,s1,e1) = scalar % boundary % hostData(ivar,s2,e2)
          ENDDO

        ELSE

          s1 = 1
          e2 = e1 - 1
          s2 = 2
          !neighborRank = decomp % elemToRank % hostData(e2Global)
          DO ivar = 1,scalar % nvar
            scalar % extBoundary % hostData(ivar,s1,e1) = scalar % boundary % hostData(ivar,s2,e2)
          ENDDO

          s1 = 2
          e2 = e1 + 1
          s2 = 1
          !neighborRank = decomp % elemToRank % hostData(e2Global)
          DO ivar = 1,scalar % nvar
            scalar % extBoundary % hostData(ivar,s1,e1) = scalar % boundary % hostData(ivar,s2,e2)
          ENDDO

        ENDIF

      ENDDO

      if( gpuAccel )then
        call scalar % extBoundary % updatedevice()
      endif


  END SUBROUTINE SideExchange_MappedScalar1D

  SUBROUTINE BassiRebaySides_MappedScalar1D(scalar,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedScalar1D),INTENT(inout) :: scalar
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iel
    INTEGER :: iside
    INTEGER :: ivar
    INTEGER :: i

      if( gpuAccel )then
        call scalar % boundary % updatehost()
        call scalar % extBoundary % updatehost()
      endif

      DO iel = 1,scalar % nElem
          DO ivar = 1,scalar % nVar

              ! Left side - we account for the -\hat{x} normal
              scalar % avgBoundary % hostData(ivar,1,iel) = -0.5_prec*( &
                                                               scalar % boundary % hostData(ivar,1,iel) + &
                                                               scalar % extBoundary % hostData(ivar,1,iel))

              ! Right side - we account for the +\hat{x} normal
              scalar % avgBoundary % hostData(ivar,2,iel) = 0.5_prec*( &
                                                               scalar % boundary % hostData(ivar,2,iel) + &
                                                               scalar % extBoundary % hostData(ivar,2,iel))
          END DO
      END DO

      if( gpuAccel )then
        call scalar % avgBoundary % updateDevice()
      endif

  END SUBROUTINE BassiRebaySides_MappedScalar1D

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

    ELSEIF (dForm == selfWeakBRForm) THEN

        IF (gpuAccel) THEN
          CALL scalar % interp % DGDerivative_1D(scalar % interior % deviceData, &
                                                 scalar % avgBoundary % deviceData, &
                                                 df % interior % deviceData, &
                                                 scalar % nVar, &
                                                 scalar % nElem)
        ELSE
          CALL scalar % interp % DGDerivative_1D(scalar % interior % hostData, &
                                                 scalar % avgBoundary % hostData, &
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

    END IF

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
    INTEGER :: iEl,iVar,i

    IF (gpuAccel) THEN

      CALL JacobianWeight_MappedScalar1D_gpu_wrapper(scalar % interior % deviceData, &
                                                     geometry % dxds % interior % deviceData, &
                                                     scalar % interp % N, &
                                                     scalar % nVar, &
                                                     scalar % nElem)
    ELSE

      DO iEl = 1,scalar % nElem
        DO iVar = 1,scalar % nVar
          DO i = 0,scalar % interp % N
            scalar % interior % hostData(i,iVar,iEl) = scalar % interior % hostData(i,iVar,iEl)/ &
                                                       geometry % dxds % interior % hostData(i,1,iEl)
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE JacobianWeight_MappedScalar1D

  SUBROUTINE SetInteriorFromEquation_MappedScalar2D( scalar, geometry, time )
  !!  Sets the scalar % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time. 
    IMPLICIT NONE
    CLASS(MappedScalar2D), INTENT(inout) :: scalar
    TYPE(SEMQuad), INTENT(in) :: geometry
    REAL(prec), INTENT(in) :: time
    ! Local
    INTEGER :: i, j, iEl, iVar
    REAL(prec) :: x
    REAL(prec) :: y


    DO iEl = 1,scalar % nElem
      DO iVar = 1, scalar % nVar
        DO j = 0, scalar % interp % N
          DO i = 0, scalar % interp % N

            ! Get the mesh positions
            x = geometry % x % interior % hostData(1,i,j,1,iEl)
            y = geometry % x % interior % hostData(2,i,j,1,iEl)

            scalar % interior % hostData(i,j,iVar,iEl) = &
              scalar % eqn(iVar) % Evaluate((/x, y, 0.0_prec, time/))

          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE SetInteriorFromEquation_MappedScalar2D

  SUBROUTINE SideExchange_MappedScalar2D(scalar,mesh,decomp,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedScalar2D),INTENT(inout) :: scalar
    TYPE(Mesh2D),INTENT(in) :: mesh
    TYPE(MPILayer),INTENT(inout) :: decomp
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1,e2,s1,s2,e2Global
    INTEGER :: flip,bcid
    INTEGER :: i1,i2,ivar
    INTEGER :: neighborRank
    INTEGER :: rankId, offset

    rankId = decomp % rankId
    offset = decomp % offsetElem % hostData(rankId)

    IF (gpuAccel) THEN

      CALL scalar % boundary % UpdateHost()
      CALL scalar % MPIExchangeAsync(decomp,mesh,resetCount=.TRUE.)
      CALL decomp % FinalizeMPIExchangeAsync()
      CALL scalar % extBoundary % UpdateDevice()

      CALL SideExchange_MappedScalar2D_gpu_wrapper(scalar % extBoundary % deviceData, &
                                                   scalar % boundary % deviceData, &
                                                   mesh % sideInfo % deviceData, &
                                                   decomp % elemToRank % deviceData, &
                                                   decomp % rankId, &
                                                   offset, &
                                                   scalar % interp % N, &
                                                   scalar % nvar, &
                                                   scalar % nElem)

    ELSE

      CALL scalar % MPIExchangeAsync(decomp,mesh,resetCount=.TRUE.)
      DO e1 = 1,mesh % nElem
        DO s1 = 1,4
          e2Global = mesh % sideInfo % hostData(3,s1,e1)
          e2 = e2Global - offset
          s2 = mesh % sideInfo % hostData(4,s1,e1)/10
          flip = mesh % sideInfo % hostData(4,s1,e1) - s2*10
          bcid = mesh % sideInfo % hostData(5,s1,e1)

          IF (bcid == 0) THEN

            neighborRank = decomp % elemToRank % hostData(e2Global)

            IF (neighborRank == decomp % rankId) THEN

              IF (flip == 0) THEN

                DO ivar = 1,scalar % nvar
                  DO i1 = 0,scalar % interp % N
                    scalar % extBoundary % hostData(i1,ivar,s1,e1) = &
                      scalar % boundary % hostData(i1,ivar,s2,e2)
                  END DO
                END DO

              ELSEIF (flip == 1) THEN

                DO ivar = 1,scalar % nvar
                  DO i1 = 0,scalar % interp % N
                    i2 = scalar % interp % N - i1
                    scalar % extBoundary % hostData(i1,ivar,s1,e1) = &
                      scalar % boundary % hostData(i2,ivar,s2,e2)
                  END DO
                END DO

              END IF

            END IF

          END IF

        END DO
      END DO
      CALL decomp % FinalizeMPIExchangeAsync()

    END IF

    CALL scalar % ApplyFlip(decomp,mesh,gpuAccel)

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

    IF (gpuAccel) THEN

      CALL BassiRebaySides_MappedScalar2D_gpu_wrapper(scalar % avgBoundary % deviceData, &
                                                      scalar % boundary % deviceData, &
                                                      scalar % extBoundary % deviceData, &
                                                      scalar % interp % N, &
                                                      scalar % nvar, &
                                                      scalar % nElem)

    ELSE

      DO iel = 1,scalar % nElem
        DO iside = 1,4
          DO ivar = 1,scalar % nVar
            DO i = 0,scalar % interp % N
              scalar % avgBoundary % hostData(i,ivar,iside,iel) = 0.5_prec*( &
                                                               scalar % boundary % hostData(i,ivar,iside,iel) + &
                                                               scalar % extBoundary % hostData(i,ivar,iside,iel))
            END DO
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE BassiRebaySides_MappedScalar2D

  SUBROUTINE Gradient_MappedScalar2D(scalar,geometry,gradF,dForm,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedScalar2D),INTENT(inout) :: scalar
    TYPE(SEMQuad),INTENT(in) :: geometry
    TYPE(MappedVector2D),INTENT(inout) :: gradF
    INTEGER, INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

      if( dForm == selfStrongForm )then
        call scalar % GradientSF_MappedScalar2D( geometry, gradF, gpuAccel ) 
      elseif( dForm == selfWeakBRForm )then
        call scalar % GradientBR_MappedScalar2D( geometry, gradF, gpuAccel ) 
      endif

  END SUBROUTINE Gradient_MappedScalar2D

  SUBROUTINE GradientBR_MappedScalar2D(scalar,geometry,gradF,gpuAccel)
    !! Calculates the gradient of a scalar 2D function using a bassi-rebay method
    !!
    !! This method will call the BassiRebaySides method, which assumes the SideExchange
    !! has already been completed, to update the avgBoundary attribute.
    !!
    IMPLICIT NONE
    CLASS(MappedScalar2D),INTENT(inout) :: scalar
    TYPE(SEMQuad),INTENT(in) :: geometry
    TYPE(MappedVector2D),INTENT(inout) :: gradF
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl
    REAL(prec) :: gFx, gFy
    REAL(prec) :: f1, f2
     
    CALL scalar % BassiRebaySides( gpuAccel )

    IF( gpuAccel )THEN

      CALL GradientBR_MappedScalar2D_gpu_wrapper(scalar % interior % deviceData, &
                                                 scalar % avgBoundary % deviceData, &
                                                 geometry % dsdx % interior % deviceData, &
                                                 geometry % J % interior % deviceData, &
                                                 geometry % nHat % boundary % deviceData, &
                                                 geometry % nScale % boundary % deviceData, &
                                                 gradF % interior % deviceData, &
                                                 scalar % interp % dgMatrix % deviceData, &
                                                 scalar % interp % bMatrix % deviceData, &
                                                 scalar % interp % qWeights % deviceData, &
                                                 scalar % interp % N, &
                                                 scalar % nVar, &
                                                 scalar % nElem)
    ELSE

     DO iEl = 1, scalar % nElem
       DO iVar = 1, scalar % nVar
         DO j = 0, scalar % interp % N
           DO i = 0, scalar % interp % N

             gFx = 0.0_prec
             gFy = 0.0_prec
             DO ii = 0, scalar % interp % N

               f1 = scalar % interior % hostData(ii,j,iVar,iEl)*&
                      geometry % dsdx % interior % hostData(1,1,ii,j,1,iEl)

               f2 = scalar % interior % hostData(i,ii,iVar,iEl)*&
                      geometry % dsdx % interior % hostData(1,2,i,ii,1,iEl)

               gFx = gFx + scalar % interp % dgMatrix % hostData(ii,i)*f1 +&
                               scalar % interp % dgMatrix % hostData(ii,j)*f2

               f1 = scalar % interior % hostData(ii,j,iVar,iEl)*&
                      geometry % dsdx % interior % hostData(2,1,ii,j,1,iEl)

               f2 = scalar % interior % hostData(i,ii,iVar,iEl)*&
                      geometry % dsdx % interior % hostData(2,2,i,ii,1,iEl)

               gFy = gFy + scalar % interp % dgMatrix % hostData(ii,i)*f1 +&
                               scalar % interp % dgMatrix % hostData(ii,j)*f2

             END DO

             ! Boundary Contribution
             f1 = scalar % avgBoundary % hostData(j,iVar,2,iEl)*&
                     geometry % nHat % boundary % hostData(1,j,1,2,iEl)*&
                     geometry % nScale % boundary % hostData(j,1,2,iEl) ! East

             f2 = scalar % avgBoundary % hostData(j,iVar,4,iEl)*&
                     geometry % nHat % boundary % hostData(1,j,1,4,iEl)*&
                     geometry % nScale % boundary % hostData(j,1,4,iEl) ! West

             gFx = gFx + (f1*scalar % interp % bMatrix % hostData(i,1) + &
                          f2*scalar % interp % bMatrix % hostData(i,0))/ &
                     scalar % interp % qWeights % hostData(i)

             f1 = scalar % avgBoundary % hostData(i,iVar,3,iEl)*&
                     geometry % nHat % boundary % hostData(1,i,1,3,iEl)*&
                     geometry % nScale % boundary % hostData(i,1,3,iEl) ! North

             f2 = scalar % avgBoundary % hostData(i,iVar,1,iEl)*&
                     geometry % nHat % boundary % hostData(1,i,1,1,iEl)*&
                     geometry % nScale % boundary % hostData(i,1,1,iEl) ! South

             gFx = gFx + (f1*scalar % interp % bMatrix % hostData(j,1) + &
                          f2*scalar % interp % bMatrix % hostData(j,0))/ &
                     scalar % interp % qWeights % hostData(j)

             f1 = scalar % avgBoundary % hostData(j,iVar,2,iEl)*&
                     geometry % nHat % boundary % hostData(2,j,1,2,iEl)*&
                     geometry % nScale % boundary % hostData(j,1,2,iEl) ! East

             f2 = scalar % avgBoundary % hostData(j,iVar,4,iEl)*&
                     geometry % nHat % boundary % hostData(2,j,1,4,iEl)*&
                     geometry % nScale % boundary % hostData(j,1,4,iEl) ! West

             gFy = gFy + (f1*scalar % interp % bMatrix % hostData(i,1) + &
                          f2*scalar % interp % bMatrix % hostData(i,0))/ &
                     scalar % interp % qWeights % hostData(i)

             f1 = scalar % avgBoundary % hostData(i,iVar,3,iEl)*&
                     geometry % nHat % boundary % hostData(2,i,1,3,iEl)*&
                     geometry % nScale % boundary % hostData(i,1,3,iEl) ! North

             f2 = scalar % avgBoundary % hostData(i,iVar,1,iEl)*&
                     geometry % nHat % boundary % hostData(2,i,1,1,iEl)*&
                     geometry % nScale % boundary % hostData(i,1,1,iEl) ! South

             gFy = gFy + (f1*scalar % interp % bMatrix % hostData(j,1) + &
                          f2*scalar % interp % bMatrix % hostData(j,0))/ &
                     scalar % interp % qWeights % hostData(j)

             gradF % interior % hostData(1,i,j,iVar,iEl) = gFx/geometry % J % interior % hostData(i,j,1,iEl)
             gradF % interior % hostData(2,i,j,iVar,iEl) = gFy/geometry % J % interior % hostData(i,j,1,iEl)

           END DO
         END DO
       END DO
     END DO

    ENDIF

  END SUBROUTINE GradientBR_MappedScalar2D

  SUBROUTINE GradientSF_MappedScalar2D(scalar,geometry,gradF,gpuAccel)
    !! Calculates the gradient of a scalar 2D function using the conservative form of the
    !! mapped gradient operator
    !!
    !! \grad_{phys}( f ) =  (1 / J)*(\partial / \partial \xi_i ( J\vec{a}_i f )
    !!
    !! where the sum over i is implied.
    IMPLICIT NONE
    CLASS(MappedScalar2D),INTENT(in) :: scalar
    TYPE(SEMQuad),INTENT(in) :: geometry
    TYPE(MappedVector2D),INTENT(inout) :: gradF
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl
    REAL(prec) :: gFx, gFy
    REAL(prec) :: f1, f2

    IF (gpuAccel) THEN
      CALL GradientSF_MappedScalar2D_gpu_wrapper(scalar % interior % deviceData, &
                                               geometry % dsdx % interior % deviceData, &
                                               geometry % J % interior % deviceData, &
                                               gradF % interior % deviceData, &
                                               scalar % interp % dMatrix % deviceData, &
                                               scalar % interp % N, &
                                               scalar % nVar, &
                                               scalar % nElem)

    ELSE


      DO iEl = 1, scalar % nElem
        DO iVar = 1, scalar % nVar
          DO j = 0, scalar % interp % N
            DO i = 0, scalar % interp % N

              gFx = 0.0_prec
              gFy = 0.0_prec
              DO ii = 0, scalar % interp % N

                f1 = scalar % interior % hostData(ii,j,iVar,iEl)*&
                       geometry % dsdx % interior % hostData(1,1,ii,j,1,iEl)

                f2 = scalar % interior % hostData(i,ii,iVar,iEl)*&
                       geometry % dsdx % interior % hostData(1,2,i,ii,1,iEl)

                gFx = gFx + scalar % interp % dMatrix % hostData(ii,i)*f1 +&
                                scalar % interp % dMatrix % hostData(ii,j)*f2

                f1 = scalar % interior % hostData(ii,j,iVar,iEl)*&
                       geometry % dsdx % interior % hostData(2,1,ii,j,1,iEl)

                f2 = scalar % interior % hostData(i,ii,iVar,iEl)*&
                       geometry % dsdx % interior % hostData(2,2,i,ii,1,iEl)

                gFy = gFy + scalar % interp % dMatrix % hostData(ii,i)*f1 +&
                                scalar % interp % dMatrix % hostData(ii,j)*f2

              END DO

              gradF % interior % hostData(1,i,j,iVar,iEl) = gFx/geometry % J % interior % hostData(i,j,1,iEl)
              gradF % interior % hostData(2,i,j,iVar,iEl) = gFy/geometry % J % interior % hostData(i,j,1,iEl)

            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE GradientSF_MappedScalar2D

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

    IF (gpuAccel) THEN

      CALL JacobianWeight_MappedScalar2D_gpu_wrapper(scalar % interior % deviceData, &
                                                     geometry % J % interior % deviceData, &
                                                     scalar % interp % N, &
                                                     scalar % nVar, &
                                                     scalar % nElem)
    ELSE

      DO iEl = 1,scalar % nElem
        DO iVar = 1,scalar % nVar
          DO j = 0,scalar % interp % N
            DO i = 0,scalar % interp % N
              scalar % interior % hostData(i,j,iVar,iEl) = scalar % interior % hostData(i,j,iVar,iEl)/ &
                                                           geometry % J % interior % hostData(i,j,1,iEl)
            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE JacobianWeight_MappedScalar2D

  FUNCTION Integral_MappedScalar2D(this, geometry, decomp, gpuAccel) RESULT( fRes )
    !! Calculates the area integral the scalar over all of the geometry.
    !! Global reduction is done across all MPI ranks when the domain
    !! decomposition indicates MPI is enabled. 
    IMPLICIT NONE
    CLASS(MappedScalar2D) :: this
    TYPE(SEMQuad) :: geometry
    TYPE(MPILayer) :: decomp
    LOGICAL :: gpuAccel
    REAL(prec) :: fRes
    ! Local
    INTEGER :: i, j, iEl
    REAL(prec) :: wi, wj, fint, Jacobian, f

      IF( gpuAccel ) THEN
        CALL this % interior % UpdateHost()
      ENDIF

      fint = 0.0_prec

      DO iEl = 1, geometry % x % nElem
        DO j = 0, geometry % x % interp % N
          DO i = 0, geometry % x % interp % N

            ! Coordinate mapping Jacobian
            Jacobian = geometry % J % interior % hostData(i,j,1,iEl)

            ! Quadrature weights
            wi = geometry % x % interp % qWeights % hostData(i) 
            wj = geometry % x % interp % qWeights % hostData(j)
            
            f = this % interior % hostData(i,j,4,iEl)
            
            fint = fint + f*wi*wj*Jacobian
          
          ENDDO
        ENDDO
      ENDDO

      CALL decomp % GlobalReduce( fint, fRes )

  END FUNCTION Integral_MappedScalar2D

  ! SideExchange_MappedScalar3D is used to populate scalar % extBoundary
  ! by finding neighboring elements that share a side and copying the neighboring
  ! elements solution % boundary data.

  SUBROUTINE SetInteriorFromEquation_MappedScalar3D( scalar, geometry, time )
  !!  Sets the scalar % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time. 
    IMPLICIT NONE
    CLASS(MappedScalar3D), INTENT(inout) :: scalar
    TYPE(SEMHex), INTENT(in) :: geometry
    REAL(prec), INTENT(in) :: time
    ! Local
    INTEGER :: i, j, k, iEl, iVar
    REAL(prec) :: x
    REAL(prec) :: y
    REAL(prec) :: z


    DO iEl = 1,scalar % nElem
      DO iVar = 1, scalar % nVar
        DO k = 0, scalar % interp % N
          DO j = 0, scalar % interp % N
            DO i = 0, scalar % interp % N

              ! Get the mesh positions
              x = geometry % x % interior % hostData(1,i,j,k,1,iEl)
              y = geometry % x % interior % hostData(2,i,j,k,1,iEl)
              z = geometry % x % interior % hostData(3,i,j,k,1,iEl)

              scalar % interior % hostData(i,j,k,iVar,iEl) = &
                scalar % eqn(iVar) % Evaluate((/x, y, z, time/))

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE SetInteriorFromEquation_MappedScalar3D

  SUBROUTINE SideExchange_MappedScalar3D(scalar,mesh,decomp,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedScalar3D),INTENT(inout) :: scalar
    TYPE(Mesh3D),INTENT(in) :: mesh
    TYPE(MPILayer),INTENT(inout) :: decomp
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1,e2,s1,s2,e2Global
    INTEGER :: flip,bcid
    INTEGER :: neighborRank
    INTEGER :: i1,i2,j1,j2,ivar
    INTEGER :: rankId, offset

      rankId = decomp % rankId
      offset = decomp % offsetElem % hostData(rankId)

    IF (gpuAccel) THEN

      CALL scalar % boundary % UpdateHost()
      CALL scalar % MPIExchangeAsync(decomp,mesh,resetCount=.TRUE.)
      CALL decomp % FinalizeMPIExchangeAsync()
      CALL scalar % extBoundary % UpdateDevice()

      CALL SideExchange_MappedScalar3D_gpu_wrapper(scalar % extBoundary % deviceData, &
                                                   scalar % boundary % deviceData, &
                                                   mesh % sideInfo % deviceData, &
                                                   decomp % elemToRank % deviceData, &
                                                   decomp % rankId, &
                                                   offset, &
                                                   scalar % interp % N, &
                                                   scalar % nvar, &
                                                   scalar % nElem)

    ELSE

      CALL scalar % MPIExchangeAsync(decomp,mesh,resetCount=.TRUE.)

      DO e1 = 1,mesh % nElem
        DO s1 = 1,6
          e2Global = mesh % sideInfo % hostData(3,s1,e1)
          e2 = e2Global - offset
          s2 = mesh % sideInfo % hostData(4,s1,e1)/10
          flip = mesh % sideInfo % hostData(4,s1,e1) - s2*10
          bcid = mesh % sideInfo % hostData(5,s1,e1)
          

          IF (bcid == 0) THEN

            neighborRank = decomp % elemToRank % hostData(e2Global)

            IF (neighborRank == decomp % rankId) THEN

              IF (flip == 0) THEN ! Orientation matches on both sides of the face

                DO ivar = 1,scalar % nvar
                  DO j1 = 0,scalar % interp % N
                    DO i1 = 0,scalar % interp % N
                      scalar % extBoundary % hostData(i1,j1,ivar,s1,e1) = &
                        scalar % boundary % hostData(i1,j1,ivar,s2,e2)
                    END DO
                  END DO
                END DO

              ELSEIF (flip == 1) THEN

                DO ivar = 1,scalar % nvar
                  DO j1 = 0,scalar % interp % N
                    DO i1 = 0,scalar % interp % N

                      i2 = j1
                      j2 = scalar % interp % N - i1
                      scalar % extBoundary % hostData(i1,j1,ivar,s1,e1) = &
                        scalar % boundary % hostData(i2,j2,ivar,s2,e2)

                    END DO
                  END DO
                END DO

              ELSEIF (flip == 2) THEN

                DO ivar = 1,scalar % nvar
                  DO j1 = 0,scalar % interp % N
                    DO i1 = 0,scalar % interp % N
                      i2 = scalar % interp % N - i1
                      j2 = scalar % interp % N - j1
                      scalar % extBoundary % hostData(i1,j1,ivar,s1,e1) = &
                        scalar % boundary % hostData(i2,j2,ivar,s2,e2)
                    END DO
                  END DO
                END DO

              ELSEIF (flip == 3) THEN

                DO ivar = 1,scalar % nvar
                  DO j1 = 0,scalar % interp % N
                    DO i1 = 0,scalar % interp % N
                      i2 = scalar % interp % N - j1
                      j2 = i1
                      scalar % extBoundary % hostData(i1,j1,ivar,s1,e1) = &
                        scalar % boundary % hostData(i2,j2,ivar,s2,e2)
                    END DO
                  END DO
                END DO

              ELSEIF (flip == 4) THEN

                DO ivar = 1,scalar % nvar
                  DO j1 = 0,scalar % interp % N
                    DO i1 = 0,scalar % interp % N
                      i2 = j1
                      j2 = i1
                      scalar % extBoundary % hostData(i1,j1,ivar,s1,e1) = &
                        scalar % boundary % hostData(i2,j2,ivar,s2,e2)
                    END DO
                  END DO
                END DO

              END IF

            END IF

          END IF

        END DO
      END DO

      CALL decomp % FinalizeMPIExchangeAsync()

    END IF

    CALL scalar % ApplyFlip(decomp,mesh,gpuAccel)

  END SUBROUTINE SideExchange_MappedScalar3D

  SUBROUTINE BassiRebaySides_MappedScalar3D(scalar,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedScalar3D),INTENT(inout) :: scalar
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iel
    INTEGER :: iside
    INTEGER :: ivar
    INTEGER :: i,j

    IF (gpuAccel) THEN

      CALL BassiRebaySides_MappedScalar3D_gpu_wrapper(scalar % avgBoundary % deviceData, &
                                                      scalar % boundary % deviceData, &
                                                      scalar % extBoundary % deviceData, &
                                                      scalar % interp % N, &
                                                      scalar % nvar, &
                                                      scalar % nElem)

    ELSE

      DO iel = 1,scalar % nElem
        DO iside = 1,6
          DO ivar = 1,scalar % nVar
            DO j = 0,scalar % interp % N
              DO i = 0,scalar % interp % N
                scalar % avgBoundary % hostData(i,j,ivar,iside,iel) = 0.5_prec*( &
                                                                   scalar % boundary % hostData(i,j,ivar,iside,iel) + &
                                                                   scalar % extBoundary % hostData(i,j,ivar,iside,iel))
              END DO
            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE BassiRebaySides_MappedScalar3D

  SUBROUTINE Gradient_MappedScalar3D(scalar,geometry,gradF,dForm,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedScalar3D),INTENT(inout) :: scalar
    TYPE(SEMHex),INTENT(in) :: geometry
    TYPE(MappedVector3D),INTENT(inout) :: gradF
    INTEGER, INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

      if( dForm == selfStrongForm )then
        call scalar % GradientSF_MappedScalar3D( geometry, gradF, gpuAccel ) 
      elseif( dForm == selfWeakBRForm )then
        call scalar % GradientBR_MappedScalar3D( geometry, gradF, gpuAccel ) 
      endif

  END SUBROUTINE Gradient_MappedScalar3D

  SUBROUTINE GradientBR_MappedScalar3D(scalar,geometry,gradF,gpuAccel)
    !! Calculates the gradient of a scalar 3D function using a bassi-rebay method
    !!
    !! This method will call the BassiRebaySides method, which assumes the SideExchange
    !! has already been completed, to update the avgBoundary attribute.
    !!
    IMPLICIT NONE
    CLASS(MappedScalar3D),INTENT(inout) :: scalar
    TYPE(SEMHex),INTENT(in) :: geometry
    TYPE(MappedVector3D),INTENT(inout) :: gradF
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER    :: i,j,k,ii,iVar,iEl
    REAL(prec) :: gFx, gFy, gFz
    REAL(prec) :: f1, f2, f3
     
    CALL scalar % BassiRebaySides( gpuAccel )

    IF( gpuAccel )THEN

      CALL GradientBR_MappedScalar3D_gpu_wrapper(scalar % interior % deviceData, &
                                                 scalar % avgBoundary % deviceData, &
                                                 geometry % dsdx % interior % deviceData, &
                                                 geometry % J % interior % deviceData, &
                                                 geometry % nHat % boundary % deviceData, &
                                                 geometry % nScale % boundary % deviceData, &
                                                 gradF % interior % deviceData, &
                                                 scalar % interp % dgMatrix % deviceData, &
                                                 scalar % interp % bMatrix % deviceData, &
                                                 scalar % interp % qWeights % deviceData, &
                                                 scalar % interp % N, &
                                                 scalar % nVar, &
                                                 scalar % nElem)
    ELSE

     DO iEl = 1, scalar % nElem
       DO iVar = 1, scalar % nVar
         DO k = 0, scalar % interp % N
           DO j = 0, scalar % interp % N
             DO i = 0, scalar % interp % N

               gFx = 0.0_prec
               gFy = 0.0_prec
               gFz = 0.0_prec
               DO ii = 0, scalar % interp % N

                 f1 = scalar % interior % hostData(ii,j,k,iVar,iEl)*&
                        geometry % dsdx % interior % hostData(1,1,ii,j,k,1,iEl)

                 f2 = scalar % interior % hostData(i,ii,k,iVar,iEl)*&
                        geometry % dsdx % interior % hostData(1,2,i,ii,k,1,iEl)

                 f3 = scalar % interior % hostData(i,j,ii,iVar,iEl)*&
                        geometry % dsdx % interior % hostData(1,3,i,j,ii,1,iEl)

                 gFx = gFx + scalar % interp % dgMatrix % hostData(ii,i)*f1 +&
                             scalar % interp % dgMatrix % hostData(ii,j)*f2 +&
                             scalar % interp % dgMatrix % hostData(ii,k)*f3

                 f1 = scalar % interior % hostData(ii,j,k,iVar,iEl)*&
                        geometry % dsdx % interior % hostData(2,1,ii,j,k,1,iEl)

                 f2 = scalar % interior % hostData(i,ii,k,iVar,iEl)*&
                        geometry % dsdx % interior % hostData(2,2,i,ii,k,1,iEl)

                 f3 = scalar % interior % hostData(i,j,ii,iVar,iEl)*&
                        geometry % dsdx % interior % hostData(2,3,i,j,ii,1,iEl)

                 gFy = gFy + scalar % interp % dgMatrix % hostData(ii,i)*f1 +&
                             scalar % interp % dgMatrix % hostData(ii,j)*f2 +&
                             scalar % interp % dgMatrix % hostData(ii,k)*f3

                 f1 = scalar % interior % hostData(ii,j,k,iVar,iEl)*&
                        geometry % dsdx % interior % hostData(3,1,ii,j,k,1,iEl)

                 f2 = scalar % interior % hostData(i,ii,k,iVar,iEl)*&
                        geometry % dsdx % interior % hostData(3,2,i,ii,k,1,iEl)

                 f3 = scalar % interior % hostData(i,j,ii,iVar,iEl)*&
                        geometry % dsdx % interior % hostData(3,3,i,j,ii,1,iEl)

                 gFz = gFz + scalar % interp % dgMatrix % hostData(ii,i)*f1 +&
                             scalar % interp % dgMatrix % hostData(ii,j)*f2 +&
                             scalar % interp % dgMatrix % hostData(ii,k)*f3

               END DO

               ! Boundary Contribution
               f1 = scalar % avgBoundary % hostData(j,k,iVar,2,iEl)*&
                       geometry % nHat % boundary % hostData(1,j,k,1,2,iEl)*&
                       geometry % nScale % boundary % hostData(j,k,1,2,iEl) ! East

               f2 = scalar % avgBoundary % hostData(j,k,iVar,4,iEl)*&
                       geometry % nHat % boundary % hostData(1,j,k,1,4,iEl)*&
                       geometry % nScale % boundary % hostData(j,k,1,4,iEl) ! West

               gFx = gFx + (f1*scalar % interp % bMatrix % hostData(i,1) + &
                            f2*scalar % interp % bMatrix % hostData(i,0))/ &
                       scalar % interp % qWeights % hostData(i)

               f1 = scalar % avgBoundary % hostData(i,k,iVar,3,iEl)*&
                       geometry % nHat % boundary % hostData(1,i,k,1,3,iEl)*&
                       geometry % nScale % boundary % hostData(i,k,1,3,iEl) ! North

               f2 = scalar % avgBoundary % hostData(i,k,iVar,1,iEl)*&
                       geometry % nHat % boundary % hostData(1,i,k,1,1,iEl)*&
                       geometry % nScale % boundary % hostData(i,k,1,1,iEl) ! South

               gFx = gFx + (f1*scalar % interp % bMatrix % hostData(j,1) + &
                            f2*scalar % interp % bMatrix % hostData(j,0))/ &
                       scalar % interp % qWeights % hostData(j)

               f1 = scalar % avgBoundary % hostData(i,j,iVar,6,iEl)*&
                       geometry % nHat % boundary % hostData(1,i,j,1,6,iEl)*&
                       geometry % nScale % boundary % hostData(i,j,1,6,iEl) ! Top

               f2 = scalar % avgBoundary % hostData(i,j,iVar,5,iEl)*&
                       geometry % nHat % boundary % hostData(1,i,j,1,5,iEl)*&
                       geometry % nScale % boundary % hostData(i,j,1,5,iEl) ! Bottom

               gFx = gFx + (f1*scalar % interp % bMatrix % hostData(k,1) + &
                            f2*scalar % interp % bMatrix % hostData(k,0))/ &
                       scalar % interp % qWeights % hostData(k)

               f1 = scalar % avgBoundary % hostData(j,k,iVar,2,iEl)*&
                       geometry % nHat % boundary % hostData(2,j,k,1,2,iEl)*&
                       geometry % nScale % boundary % hostData(j,k,1,2,iEl) ! East

               f2 = scalar % avgBoundary % hostData(j,k,iVar,4,iEl)*&
                       geometry % nHat % boundary % hostData(2,j,k,1,4,iEl)*&
                       geometry % nScale % boundary % hostData(j,k,1,4,iEl) ! West

               gFy = gFy + (f1*scalar % interp % bMatrix % hostData(i,1) + &
                            f2*scalar % interp % bMatrix % hostData(i,0))/ &
                       scalar % interp % qWeights % hostData(i)

               f1 = scalar % avgBoundary % hostData(i,k,iVar,3,iEl)*&
                       geometry % nHat % boundary % hostData(2,i,k,1,3,iEl)*&
                       geometry % nScale % boundary % hostData(i,k,1,3,iEl) ! North

               f2 = scalar % avgBoundary % hostData(i,k,iVar,1,iEl)*&
                       geometry % nHat % boundary % hostData(2,i,k,1,1,iEl)*&
                       geometry % nScale % boundary % hostData(i,k,1,1,iEl) ! South

               gFy = gFy + (f1*scalar % interp % bMatrix % hostData(j,1) + &
                            f2*scalar % interp % bMatrix % hostData(j,0))/ &
                       scalar % interp % qWeights % hostData(j)

               f1 = scalar % avgBoundary % hostData(i,j,iVar,6,iEl)*&
                       geometry % nHat % boundary % hostData(2,i,j,1,6,iEl)*&
                       geometry % nScale % boundary % hostData(i,j,1,6,iEl) ! Top

               f2 = scalar % avgBoundary % hostData(i,j,iVar,5,iEl)*&
                       geometry % nHat % boundary % hostData(2,i,j,1,5,iEl)*&
                       geometry % nScale % boundary % hostData(i,j,1,5,iEl) ! Bottom

               gFy = gFy + (f1*scalar % interp % bMatrix % hostData(k,1) + &
                            f2*scalar % interp % bMatrix % hostData(k,0))/ &
                       scalar % interp % qWeights % hostData(k)

               f1 = scalar % avgBoundary % hostData(j,k,iVar,2,iEl)*&
                       geometry % nHat % boundary % hostData(3,j,k,1,2,iEl)*&
                       geometry % nScale % boundary % hostData(j,k,1,2,iEl) ! East

               f2 = scalar % avgBoundary % hostData(j,k,iVar,4,iEl)*&
                       geometry % nHat % boundary % hostData(3,j,k,1,4,iEl)*&
                       geometry % nScale % boundary % hostData(j,k,1,4,iEl) ! West

               gFz = gFz + (f1*scalar % interp % bMatrix % hostData(i,1) + &
                            f2*scalar % interp % bMatrix % hostData(i,0))/ &
                       scalar % interp % qWeights % hostData(i)

               f1 = scalar % avgBoundary % hostData(i,k,iVar,3,iEl)*&
                       geometry % nHat % boundary % hostData(3,i,k,1,3,iEl)*&
                       geometry % nScale % boundary % hostData(i,k,1,3,iEl) ! North

               f2 = scalar % avgBoundary % hostData(i,k,iVar,1,iEl)*&
                       geometry % nHat % boundary % hostData(3,i,k,1,1,iEl)*&
                       geometry % nScale % boundary % hostData(i,k,1,1,iEl) ! South

               gFz = gFz + (f1*scalar % interp % bMatrix % hostData(j,1) + &
                            f2*scalar % interp % bMatrix % hostData(j,0))/ &
                       scalar % interp % qWeights % hostData(j)

               f1 = scalar % avgBoundary % hostData(i,j,iVar,6,iEl)*&
                       geometry % nHat % boundary % hostData(3,i,j,1,6,iEl)*&
                       geometry % nScale % boundary % hostData(i,j,1,6,iEl) ! Top

               f2 = scalar % avgBoundary % hostData(i,j,iVar,5,iEl)*&
                       geometry % nHat % boundary % hostData(3,i,j,1,5,iEl)*&
                       geometry % nScale % boundary % hostData(i,j,1,5,iEl) ! Bottom

               gFz = gFz + (f1*scalar % interp % bMatrix % hostData(k,1) + &
                            f2*scalar % interp % bMatrix % hostData(k,0))/ &
                       scalar % interp % qWeights % hostData(k)

               gradF % interior % hostData(1,i,j,k,iVar,iEl) = gFx/geometry % J % interior % hostData(i,j,k,1,iEl)
               gradF % interior % hostData(2,i,j,k,iVar,iEl) = gFy/geometry % J % interior % hostData(i,j,k,1,iEl)
               gradF % interior % hostData(3,i,j,k,iVar,iEl) = gFz/geometry % J % interior % hostData(i,j,k,1,iEl)

             END DO
           END DO
         END DO
       END DO
     END DO

    ENDIF

  END SUBROUTINE GradientBR_MappedScalar3D

  SUBROUTINE GradientSF_MappedScalar3D(scalar,geometry,gradF,gpuAccel)
    !! Calculates the gradient of a scalar 3D function using the conservative form of the
    !! mapped gradient operator
    !!
    !! \grad_{phys}( f ) =  (1 / J)*(\partial / \partial \xi_i ( J\vec{a}_i f )
    !!
    !! where the sum over i is implied.
    IMPLICIT NONE
    CLASS(MappedScalar3D),INTENT(in) :: scalar
    TYPE(SEMHex),INTENT(in) :: geometry
    TYPE(MappedVector3D),INTENT(inout) :: gradF
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER    :: i,j,k,ii,iVar,iEl
    REAL(prec) :: gFx, gFy, gFz
    REAL(prec) :: f1, f2, f3

    IF (gpuAccel) THEN
      CALL GradientSF_MappedScalar3D_gpu_wrapper(scalar % interior % deviceData, &
                                               geometry % dsdx % interior % deviceData, &
                                               geometry % J % interior % deviceData, &
                                               gradF % interior % deviceData, &
                                               scalar % interp % dMatrix % deviceData, &
                                               scalar % interp % N, &
                                               scalar % nVar, &
                                               scalar % nElem)

    ELSE


      DO iEl = 1, scalar % nElem
        DO iVar = 1, scalar % nVar
          DO k = 0, scalar % interp % N
            DO j = 0, scalar % interp % N
              DO i = 0, scalar % interp % N

                gFx = 0.0_prec
                gFy = 0.0_prec
                gFz = 0.0_prec
                DO ii = 0, scalar % interp % N

                  f1 = scalar % interior % hostData(ii,j,k,iVar,iEl)*&
                         geometry % dsdx % interior % hostData(1,1,ii,j,k,1,iEl)

                  f2 = scalar % interior % hostData(i,ii,k,iVar,iEl)*&
                         geometry % dsdx % interior % hostData(1,2,i,ii,k,1,iEl)

                  f3 = scalar % interior % hostData(i,j,ii,iVar,iEl)*&
                         geometry % dsdx % interior % hostData(1,3,i,j,ii,1,iEl)

                  gFx = gFx + scalar % interp % dMatrix % hostData(ii,i)*f1 +&
                              scalar % interp % dMatrix % hostData(ii,j)*f2 +&
                              scalar % interp % dMatrix % hostData(ii,k)*f3

                  f1 = scalar % interior % hostData(ii,j,k,iVar,iEl)*&
                         geometry % dsdx % interior % hostData(2,1,ii,j,k,1,iEl)

                  f2 = scalar % interior % hostData(i,ii,k,iVar,iEl)*&
                         geometry % dsdx % interior % hostData(2,2,i,ii,k,1,iEl)

                  f3 = scalar % interior % hostData(i,j,ii,iVar,iEl)*&
                         geometry % dsdx % interior % hostData(2,2,i,j,ii,1,iEl)

                  gFy = gFy + scalar % interp % dMatrix % hostData(ii,i)*f1 +&
                              scalar % interp % dMatrix % hostData(ii,j)*f2 +&
                              scalar % interp % dMatrix % hostData(ii,k)*f3

                  f1 = scalar % interior % hostData(ii,j,k,iVar,iEl)*&
                         geometry % dsdx % interior % hostData(3,1,ii,j,k,1,iEl)

                  f2 = scalar % interior % hostData(i,ii,k,iVar,iEl)*&
                         geometry % dsdx % interior % hostData(3,2,i,ii,k,1,iEl)

                  f3 = scalar % interior % hostData(i,j,ii,iVar,iEl)*&
                         geometry % dsdx % interior % hostData(3,2,i,j,ii,1,iEl)

                  gFz = gFz + scalar % interp % dMatrix % hostData(ii,i)*f1 +&
                              scalar % interp % dMatrix % hostData(ii,j)*f2 +&
                              scalar % interp % dMatrix % hostData(ii,k)*f3

                END DO

                gradF % interior % hostData(1,i,j,k,iVar,iEl) = gFx/geometry % J % interior % hostData(i,j,k,1,iEl)
                gradF % interior % hostData(2,i,j,k,iVar,iEl) = gFy/geometry % J % interior % hostData(i,j,k,1,iEl)
                gradF % interior % hostData(3,i,j,k,iVar,iEl) = gFz/geometry % J % interior % hostData(i,j,k,1,iEl)

              END DO
            END DO
          ENDDO
        END DO
      END DO

    END IF

  END SUBROUTINE GradientSF_MappedScalar3D

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

    IF (gpuAccel) THEN

      CALL JacobianWeight_MappedScalar3D_gpu_wrapper(scalar % interior % deviceData, &
                                                     geometry % J % interior % deviceData, &
                                                     scalar % interp % N, &
                                                     scalar % nVar, &
                                                     scalar % nElem)

    ELSE

      DO iEl = 1,scalar % nElem
        DO iVar = 1,scalar % nVar
          DO k = 0,scalar % interp % N
            DO j = 0,scalar % interp % N
              DO i = 0,scalar % interp % N
                scalar % interior % hostData(i,j,k,iVar,iEl) = scalar % interior % hostData(i,j,k,iVar,iEl)/ &
                                                               geometry % J % interior % hostData(i,j,k,1,iEl)
              END DO
            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE JacobianWeight_MappedScalar3D

  ! ---------------------- Vectors ---------------------- !

  SUBROUTINE SetInteriorFromEquation_MappedVector2D( vector, geometry, time )
  !!  Sets the scalar % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time. 
    IMPLICIT NONE
    CLASS(MappedVector2D), INTENT(inout) :: vector
    TYPE(SEMQuad), INTENT(in) :: geometry
    REAL(prec), INTENT(in) :: time
    ! Local
    INTEGER :: i, j, iEl, iVar
    REAL(prec) :: x
    REAL(prec) :: y


    DO iEl = 1,vector % nElem
      DO iVar = 1, vector % nVar
        DO j = 0, vector % interp % N
          DO i = 0, vector % interp % N

            ! Get the mesh positions
            x = geometry % x % interior % hostData(1,i,j,1,iEl)
            y = geometry % x % interior % hostData(2,i,j,1,iEl)

            vector % interior % hostData(1,i,j,iVar,iEl) = &
              vector % eqn(1+2*(iVar-1)) % Evaluate((/x, y, 0.0_prec, time/))

            vector % interior % hostData(2,i,j,iVar,iEl) = &
              vector % eqn(2+2*(iVar-1)) % Evaluate((/x, y, 0.0_prec, time/))

          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE SetInteriorFromEquation_MappedVector2D

  SUBROUTINE SideExchange_MappedVector2D(vector,mesh,decomp,gpuAccel)
  !! SideExchange_MappedVectorvector2D is used to populate vector % extBoundary
  !! by finding neighboring elements that share a side and copying the neighboring
  !! elements solution % boundary data.
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(inout) :: vector
    TYPE(Mesh2D),INTENT(in) :: mesh
    TYPE(MPILayer),INTENT(inout) :: decomp
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1,e2,s1,s2,e2Global
    INTEGER :: flip,bcid
    INTEGER :: neighborRank
    INTEGER :: i1,i2,ivar
    INTEGER :: rankId, offset

      rankId = decomp % rankId
      offset = decomp % offsetElem % hostData(rankId)

    IF (gpuAccel) THEN

      CALL vector % boundary % UpdateHost()
      CALL vector % MPIExchangeAsync(decomp,mesh,resetCount=.TRUE.)
      CALL decomp % FinalizeMPIExchangeAsync()
      CALL vector % extBoundary % UpdateDevice()

      CALL SideExchange_MappedVector2D_gpu_wrapper(vector % extBoundary % deviceData, &
                                                   vector % boundary % deviceData, &
                                                   mesh % sideInfo % deviceData, &
                                                   decomp % elemToRank % deviceData, &
                                                   decomp % rankId, &
                                                   offset, &
                                                   vector % interp % N, &
                                                   vector % nvar, &
                                                   vector % nElem)

    ELSE

      CALL vector % MPIExchangeAsync(decomp,mesh,resetCount=.TRUE.)

      DO e1 = 1,mesh % nElem
        DO s1 = 1,4
          e2Global = mesh % sideInfo % hostData(3,s1,e1)
          e2 = e2Global - offset
          s2 = mesh % sideInfo % hostData(4,s1,e1)/10
          flip = mesh % sideInfo % hostData(4,s1,e1) - s2*10
          bcid = mesh % sideInfo % hostData(5,s1,e1)

          IF (bcid == 0) THEN

            neighborRank = decomp % elemToRank % hostData(e2Global)

            IF (neighborRank == decomp % rankId) THEN

              IF (flip == 0) THEN

                DO ivar = 1,vector % nvar
                  DO i1 = 0,vector % interp % N
                    vector % extBoundary % hostData(1:2,i1,ivar,s1,e1) = &
                      vector % boundary % hostData(1:2,i1,ivar,s2,e2)
                  END DO
                END DO

              ELSEIF (flip == 1) THEN

                DO ivar = 1,vector % nvar
                  DO i1 = 0,vector % interp % N
                    i2 = vector % interp % N - i1
                    vector % extBoundary % hostData(1:2,i1,ivar,s1,e1) = &
                      vector % boundary % hostData(1:2,i2,ivar,s2,e2)
                  END DO
                END DO

              END IF

            END IF

          END IF

        END DO
      END DO

      CALL decomp % FinalizeMPIExchangeAsync()

    END IF

    CALL vector % ApplyFlip(decomp,mesh,gpuAccel)

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

    IF (gpuAccel) THEN

      CALL BassiRebaySides_MappedVector2D_gpu_wrapper(vector % extBoundary % deviceData, &
                                                      vector % boundary % deviceData, &
                                                      vector % interp % N, &
                                                      vector % nvar, &
                                                      vector % nElem)

    ELSE

      DO iel = 1,vector % nElem
        DO iside = 1,4
          DO ivar = 1,vector % nVar
            DO i = 0,vector % interp % N
              vector % boundary % hostData(1:2,i,ivar,iside,iel) = 0.5_prec*( &
                                                                  vector % boundary % hostData(1:2,i,ivar,iside,iel) + &
                                                                  vector % extBoundary % hostData(1:2,i,ivar,iside,iel))
            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE BassiRebaySides_MappedVector2D

  SUBROUTINE Divergence_MappedVector2D(compVector,geometry,divVector,dForm,gpuAccel)
    ! Strong Form Operator
    !
    ! DG Weak Form Operator
    !
    ! Assumes vector has been projected to computational coordinates
    !
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(in) :: compVector
    TYPE(SEMQuad),INTENT(in) :: geometry
    TYPE(MappedScalar2D),INTENT(inout) :: divVector
    INTEGER,INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

    IF (dForm == selfWeakDGForm) THEN

      IF (gpuAccel) THEN
        CALL compVector % interp % VectorDGDivergence_2D(compVector % interior % deviceData, &
                                                         compVector % boundaryNormal % deviceData, &
                                                         divVector % interior % deviceData, &
                                                         compVector % nvar, &
                                                         compVector % nelem)
      ELSE
        CALL compVector % interp % VectorDGDivergence_2D(compVector % interior % hostData, &
                                                         compVector % boundaryNormal % hostData, &
                                                         divVector % interior % hostData, &
                                                         compVector % nvar, &
                                                         compVector % nelem)
      END IF

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
      END IF

    END IF

    CALL divVector % JacobianWeight(geometry,gpuAccel)

  END SUBROUTINE Divergence_MappedVector2D

  SUBROUTINE ContravariantProjection_MappedVector2D(vector,geometry,gpuAccel)
#undef __FUNC__
#define __FUNC__ "ContravariantProjection_MappedVector2D"
    ! Takes a vector that has physical space coordinate directions (x,y,z) and projects the vector
    ! into the the contravariant basis vector directions. Keep in mind that the contravariant basis
    ! vectors are really the Jacobian weighted contravariant basis vectors
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(inout) :: vector
    TYPE(SEMQuad),INTENT(in) :: geometry
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: i,j,ivar,iel
    REAL(prec) :: Fx, Fy

    IF (gpuAccel) THEN

      CALL ContravariantProjection_MappedVector2D_gpu_wrapper(vector % interior % deviceData, &
                                                              geometry % dsdx % interior % deviceData, &
                                                              vector % interp % N, &
                                                              vector % nVar, &
                                                              vector % nElem)

    ELSE
      ! Assume that tensor(j,i) is vector i, component j
      ! => dot product is done along first dimension
      ! to project onto computational space
      DO iel = 1,vector % nElem
        DO ivar = 1,vector % nVar
          DO j = 0,vector % interp % N
            DO i = 0,vector % interp % N

              Fx = vector % interior % hostData(1,i,j,ivar,iel)
              Fy = vector % interior % hostData(2,i,j,ivar,iel)

              vector % interior % hostData(1,i,j,ivar,iel) = &
                geometry % dsdx % interior % hostData(1,1,i,j,1,iel)*Fx + &
                geometry % dsdx % interior % hostData(2,1,i,j,1,iel)*Fy

              vector % interior % hostData(2,i,j,ivar,iel) = &
                geometry % dsdx % interior % hostData(1,2,i,j,1,iel)*Fx + &
                geometry % dsdx % interior % hostData(2,2,i,j,1,iel)*Fy

            END DO
          END DO
        END DO
      END DO

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

    IF (gpuAccel) THEN

      CALL JacobianWeight_MappedVector2D_gpu_wrapper(vector % interior % deviceData, &
                                                     geometry % J % interior % deviceData, &
                                                     vector % interp % N, &
                                                     vector % nVar, &
                                                     vector % nElem)
    ELSE

      DO iEl = 1,vector % nElem
        DO iVar = 1,vector % nVar
          DO j = 0,vector % interp % N
            DO i = 0,vector % interp % N
              vector % interior % hostData(1,i,j,iVar,iEl) = vector % interior % hostData(1,i,j,iVar,iEl)/ &
                                                             geometry % J % interior % hostData(i,j,1,iEl)
              vector % interior % hostData(2,i,j,iVar,iEl) = vector % interior % hostData(2,i,j,iVar,iEl)/ &
                                                             geometry % J % interior % hostData(i,j,1,iEl)
            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE JacobianWeight_MappedVector2D

  SUBROUTINE SetInteriorFromEquation_MappedVector3D( vector, geometry, time )
  !!  Sets the scalar % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time. 
    IMPLICIT NONE
    CLASS(MappedVector3D), INTENT(inout) :: vector
    TYPE(SEMHex), INTENT(in) :: geometry
    REAL(prec), INTENT(in) :: time
    ! Local
    INTEGER :: i, j, k, iEl, iVar
    REAL(prec) :: x
    REAL(prec) :: y
    REAL(prec) :: z


    DO iEl = 1,vector % nElem
      DO iVar = 1, vector % nVar
        DO k = 0, vector % interp % N
          DO j = 0, vector % interp % N
            DO i = 0, vector % interp % N

              ! Get the mesh positions
              x = geometry % x % interior % hostData(1,i,j,k,1,iEl)
              y = geometry % x % interior % hostData(2,i,j,k,1,iEl)
              z = geometry % x % interior % hostData(3,i,j,k,1,iEl)

              vector % interior % hostData(1,i,j,k,iVar,iEl) = &
                vector % eqn(1+3*(iVar-1)) % Evaluate((/x, y, z, time/))

              vector % interior % hostData(2,i,j,k,iVar,iEl) = &
                vector % eqn(2+3*(iVar-1)) % Evaluate((/x, y, z, time/))

              vector % interior % hostData(3,i,j,k,iVar,iEl) = &
                vector % eqn(3+3*(iVar-1)) % Evaluate((/x, y, z, time/))

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE SetInteriorFromEquation_MappedVector3D

  ! SideExchange_MappedVector3D is used to populate vector % extBoundary
  ! by finding neighboring elements that share a side and copying the neighboring
  ! elements solution % boundary data.

  SUBROUTINE SideExchange_MappedVector3D(vector,mesh,decomp,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(inout) :: vector
    TYPE(Mesh3D),INTENT(in) :: mesh
    TYPE(MPILayer),INTENT(inout) :: decomp
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1,e2,s1,s2,e2Global
    INTEGER :: flip,bcid
    INTEGER :: neighborRank
    INTEGER :: i1,i2,j1,j2,ivar
    INTEGER :: rankId, offset

      rankId = decomp % rankId
      offset = decomp % offsetElem % hostData(rankId)

    IF (gpuAccel) THEN

      CALL vector % boundary % UpdateHost()
      CALL vector % MPIExchangeAsync(decomp,mesh,resetCount=.TRUE.)
      CALL decomp % FinalizeMPIExchangeAsync()
      CALL vector % extBoundary % UpdateDevice()

      CALL SideExchange_MappedVector3D_gpu_wrapper(vector % extBoundary % deviceData, &
                                                   vector % boundary % deviceData, &
                                                   mesh % sideInfo % deviceData, &
                                                   decomp % elemToRank % deviceData, &
                                                   decomp % rankId, &
                                                   offset, &
                                                   vector % interp % N, &
                                                   vector % nvar, &
                                                   vector % nElem)

    ELSE

      CALL vector % MPIExchangeAsync(decomp,mesh,resetCount=.TRUE.)

      DO e1 = 1,mesh % nElem
        DO s1 = 1,6
          e2Global = mesh % sideInfo % hostData(3,s1,e1)
          e2 = e2Global - offset
          s2 = mesh % sideInfo % hostData(4,s1,e1)/10
          flip = mesh % sideInfo % hostData(4,s1,e1) - s2*10
          bcid = mesh % sideInfo % hostData(5,s1,e1)

          IF (bcid == 0) THEN ! Interior

            neighborRank = decomp % elemToRank % hostData(e2Global)

            IF (neighborRank == decomp % rankId) THEN

              IF (flip == 0) THEN

                DO ivar = 1,vector % nvar
                  DO j1 = 0,vector % interp % N
                    DO i1 = 0,vector % interp % N
                      vector % extBoundary % hostData(1:3,i1,j1,ivar,s1,e1) = &
                        vector % boundary % hostData(1:3,i1,j1,ivar,s2,e2)
                    END DO
                  END DO
                END DO

              ELSEIF (flip == 1) THEN

                DO ivar = 1,vector % nvar
                  DO j1 = 0,vector % interp % N
                    DO i1 = 0,vector % interp % N
                      i2 = j1
                      j2 = vector % interp % N - i1
                      vector % extBoundary % hostData(1:3,i1,j1,ivar,s1,e1) = &
                        vector % boundary % hostData(1:3,i2,j2,ivar,s2,e2)
                    END DO
                  END DO
                END DO

              ELSEIF (flip == 2) THEN

                DO ivar = 1,vector % nvar
                  DO j1 = 0,vector % interp % N
                    DO i1 = 0,vector % interp % N
                      i2 = vector % interp % N - i1
                      j2 = vector % interp % N - j1
                      vector % extBoundary % hostData(1:3,i1,j1,ivar,s1,e1) = &
                        vector % boundary % hostData(1:3,i2,j2,ivar,s2,e2)
                    END DO
                  END DO
                END DO

              ELSEIF (flip == 3) THEN

                DO ivar = 1,vector % nvar
                  DO j1 = 0,vector % interp % N
                    DO i1 = 0,vector % interp % N
                      i2 = vector % interp % N - j1
                      j2 = i1
                      vector % extBoundary % hostData(1:3,i1,j1,ivar,s1,e1) = &
                        vector % boundary % hostData(1:3,i2,j2,ivar,s2,e2)
                    END DO
                  END DO
                END DO

              ELSEIF (flip == 4) THEN

                DO ivar = 1,vector % nvar
                  DO j1 = 0,vector % interp % N
                    DO i1 = 0,vector % interp % N
                      i2 = j1
                      j2 = i1
                      vector % extBoundary % hostData(1:3,i1,j1,ivar,s1,e1) = &
                        vector % boundary % hostData(1:3,i2,j2,ivar,s2,e2)
                    END DO
                  END DO
                END DO

              END IF

            END IF

          END IF

        END DO
      END DO

      CALL decomp % FinalizeMPIExchangeAsync()

    END IF

    CALL vector % ApplyFlip(decomp,mesh,gpuAccel)

  END SUBROUTINE SideExchange_MappedVector3D

  SUBROUTINE BassiRebaySides_MappedVector3D(vector,gpuAccel)
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(inout) :: vector
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iel
    INTEGER :: iside
    INTEGER :: ivar
    INTEGER :: i,j

    IF (gpuAccel) THEN

      CALL BassiRebaySides_MappedVector3D_gpu_wrapper(vector % extBoundary % deviceData, &
                                                      vector % boundary % deviceData, &
                                                      vector % interp % N, &
                                                      vector % nvar, &
                                                      vector % nElem)
    ELSE

      DO iel = 1,vector % nElem
        DO iside = 1,6
          DO ivar = 1,vector % nVar
            DO j = 0,vector % interp % N
              DO i = 0,vector % interp % N
                vector % boundary % hostData(1:3,i,j,ivar,iside,iel) = 0.5_prec*( &
                                                                vector % boundary % hostData(1:3,i,j,ivar,iside,iel) + &
                                                                vector % extBoundary % hostData(1:3,i,j,ivar,iside,iel))
              END DO
            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE BassiRebaySides_MappedVector3D

  SUBROUTINE Divergence_MappedVector3D(compVector,geometry,divVector,dForm,gpuAccel)
    !
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(in) :: compVector
    TYPE(SEMHex),INTENT(in) :: geometry
    TYPE(MappedScalar3D),INTENT(inout) :: divVector
    INTEGER,INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

    IF (dForm == selfWeakDGForm) THEN

      IF (gpuAccel) THEN
        CALL compVector % interp % VectorDGDivergence_3D(compVector % interior % deviceData, &
                                                         compVector % boundaryNormal % deviceData, &
                                                         divVector % interior % deviceData, &
                                                         compVector % nvar, &
                                                         compVector % nelem)
      ELSE
        CALL compVector % interp % VectorDGDivergence_3D(compVector % interior % hostData, &
                                                         compVector % boundaryNormal % hostData, &
                                                         divVector % interior % hostData, &
                                                         compVector % nvar, &
                                                         compVector % nelem)
      END IF

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
      END IF

    END IF

    CALL divVector % JacobianWeight(geometry,gpuAccel)

  END SUBROUTINE Divergence_MappedVector3D

  SUBROUTINE ContravariantProjection_MappedVector3D(vector,geometry,gpuAccel)
#undef __FUNC__
#define __FUNC__ "ContravariantProjection_MappedVector3D"
    ! Takes a vector that has physical space coordinate directions (x,y,z) and projects the vector
    ! into the the contravariant basis vector directions. Keep in mind that the contravariant basis
    ! vectors are really the Jacobian weighted contravariant basis vectors
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(inout) :: vector
    TYPE(SEMHex),INTENT(in) :: geometry
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: i,j,k,iVar,iEl
    REAL(prec) :: Fx, Fy, Fz

    IF (gpuAccel) THEN

      CALL ContravariantProjection_MappedVector3D_gpu_wrapper(vector % interior % deviceData, &
                                                              geometry % dsdx % interior % deviceData, &
                                                              vector % interp % N, &
                                                              vector % nVar, &
                                                              vector % nElem)

    ELSE
      ! Assume that tensor(j,i) is vector i, component j
      ! => dot product is done along first dimension to
      ! project onto computational space
      DO iEl = 1,vector % nElem
        DO iVar = 1,vector % nVar
          DO k = 0,vector % interp % N
            DO j = 0,vector % interp % N
              DO i = 0,vector % interp % N

                Fx = vector % interior % hostData(1,i,j,k,iVar,iEl)
                Fy = vector % interior % hostData(2,i,j,k,iVar,iEl)
                Fz = vector % interior % hostData(3,i,j,k,iVar,iEl)

                vector % interior % hostData(1,i,j,k,iVar,iEl) = &
                  geometry % dsdx % interior % hostData(1,1,i,j,k,1,iEl)*Fx + &
                  geometry % dsdx % interior % hostData(2,1,i,j,k,1,iEl)*Fy + &
                  geometry % dsdx % interior % hostData(3,1,i,j,k,1,iEl)*Fz

                vector % interior % hostData(2,i,j,k,iVar,iEl) = &
                  geometry % dsdx % interior % hostData(1,2,i,j,k,1,iEl)*Fx + &
                  geometry % dsdx % interior % hostData(2,2,i,j,k,1,iEl)*Fy + &
                  geometry % dsdx % interior % hostData(3,2,i,j,k,1,iEl)*Fz

                vector % interior % hostData(3,i,j,k,iVar,iEl) = &
                  geometry % dsdx % interior % hostData(1,3,i,j,k,1,iEl)*Fx + &
                  geometry % dsdx % interior % hostData(2,3,i,j,k,1,iEl)*Fy + &
                  geometry % dsdx % interior % hostData(3,3,i,j,k,1,iEl)*Fz

              END DO
            END DO
          END DO
        END DO
      END DO

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

    IF (gpuAccel) THEN

      CALL JacobianWeight_MappedVector3D_gpu_wrapper(vector % interior % deviceData, &
                                                     geometry % J % interior % deviceData, &
                                                     vector % interp % N, &
                                                     vector % nVar, &
                                                     vector % nElem)
    ELSE

      DO iEl = 1,vector % nElem
        DO iVar = 1,vector % nVar
          DO k = 0,vector % interp % N
            DO j = 0,vector % interp % N
              DO i = 0,vector % interp % N
                vector % interior % hostData(1,i,j,k,iVar,iEl) = vector % interior % hostData(1,i,j,k,iVar,iEl)/ &
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
                vector % interior % hostData(2,i,j,k,iVar,iEl) = vector % interior % hostData(2,i,j,k,iVar,iEl)/ &
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
                vector % interior % hostData(3,i,j,k,iVar,iEl) = vector % interior % hostData(3,i,j,k,iVar,iEl)/ &
                                                                 geometry % J % interior % hostData(i,j,k,1,iEl)
              END DO
            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE JacobianWeight_MappedVector3D

  ! --- MPI Routines --- !

  SUBROUTINE MPIExchangeAsync_MappedScalar2D(scalar,mpiHandler,mesh,resetCount)
    IMPLICIT NONE
    CLASS(MappedScalar2D),INTENT(inout) :: scalar
    TYPE(MPILayer),INTENT(inout) :: mpiHandler
    TYPE(Mesh2D),INTENT(in) :: mesh
    LOGICAL,INTENT(in) :: resetCount
    ! Local
    INTEGER :: e1,s1,e2,s2
    INTEGER :: globalSideId,r2
    INTEGER :: iError
    INTEGER :: msgCount

    IF (mpiHandler % mpiEnabled) THEN
      IF (resetCount) THEN
        msgCount = 0
      ELSE
        msgCount = mpiHandler % msgCount
      END IF

      DO e1 = 1,scalar % nElem
        DO s1 = 1,4

          e2 = mesh % sideInfo % hostData(3,s1,e1) ! Neighbor Element
          IF( e2 > 0 )THEN
            r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

            IF (r2 /= mpiHandler % rankId) THEN

              s2 = mesh % sideInfo % hostData(4,s1,e1)/10
              globalSideId = ABS(mesh % sideInfo % hostdata(2,s1,e1))

              msgCount = msgCount + 1
              CALL MPI_IRECV(scalar % extBoundary % hostData(:,:,s1,e1), &
                             (scalar % interp % N + 1)*scalar % nVar, &
                             mpiHandler % mpiPrec, &
                             r2,globalSideId, &
                             mpiHandler % mpiComm, &
                             mpiHandler % requests(msgCount),iError)

              msgCount = msgCount + 1
              CALL MPI_ISEND(scalar % boundary % hostData(:,:,s1,e1), &
                             (scalar % interp % N + 1)*scalar % nVar, &
                             mpiHandler % mpiPrec, &
                             r2,globalSideId, &
                             mpiHandler % mpiComm, &
                             mpiHandler % requests(msgCount),iError)
            END IF
          END IF

        END DO
      END DO

      mpiHandler % msgCount = msgCount
    END IF

  END SUBROUTINE MPIExchangeAsync_MappedScalar2D
!
  SUBROUTINE ApplyFlip_MappedScalar2D(scalar,mpiHandler,mesh,gpuAccel)
    ! Apply side flips to sides where MPI exchanges took place.
    IMPLICIT NONE
    CLASS(MappedScalar2D),INTENT(inout) :: scalar
    TYPE(MPILayer),INTENT(inout) :: mpiHandler
    TYPE(Mesh2D),INTENT(in) :: mesh
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1,s1,e2,s2
    INTEGER :: i,i2
    INTEGER :: r2,flip,ivar
    INTEGER :: globalSideId
    INTEGER :: bcid
    REAL(prec) :: extBuff(0:scalar % interp % N)

    IF (mpiHandler % mpiEnabled) THEN
      IF (gpuAccel) THEN

        CALL ApplyFlip_MappedScalar2D_gpu_wrapper(scalar % extBoundary % deviceData, &
                                                  mesh % sideInfo % deviceData, &
                                                  mpiHandler % elemToRank % deviceData, &
                                                  mpiHandler % rankId, &
                                                  scalar % interp % N, &
                                                  scalar % nVar, &
                                                  scalar % nElem)
      ELSE
        DO e1 = 1,scalar % nElem
          DO s1 = 1,4

            e2 = mesh % sideInfo % hostData(3,s1,e1) ! Neighbor Element
            bcid = mesh % sideInfo % hostData(5,s1,e1)
            IF (bcid == 0) THEN ! Interior Element
              r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

              IF (r2 /= mpiHandler % rankId) THEN

                s2 = mesh % sideInfo % hostData(4,s1,e1)/10
                flip = mesh % sideInfo % hostData(4,s1,e1) - s2*10
                globalSideId = mesh % sideInfo % hostdata(2,s1,e1)

                ! Need to update extBoundary with flip applied
                IF (flip == 1) THEN

                  DO ivar = 1,scalar % nvar
                    DO i = 0,scalar % interp % N
                      i2 = scalar % interp % N - i
                      extBuff(i) = scalar % extBoundary % hostData(i2,ivar,s1,e1)
                    END DO
                    DO i = 0,scalar % interp % N
                      scalar % extBoundary % hostData(i,ivar,s1,e1) = extBuff(i)
                    END DO
                  END DO

                END IF
              END IF

            ENDIF

          END DO
        END DO
      END IF
    END IF

  END SUBROUTINE ApplyFlip_MappedScalar2D

  SUBROUTINE MPIExchangeAsync_MappedVector2D(vector,mpiHandler,mesh,resetCount)
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(inout) :: vector
    TYPE(MPILayer),INTENT(inout) :: mpiHandler
    TYPE(Mesh2D),INTENT(in) :: mesh
    LOGICAL,INTENT(in) :: resetCount
    ! Local
    INTEGER :: e1,s1,e2,s2
    INTEGER :: globalSideId,r2
    INTEGER :: iError
    INTEGER :: msgCount

    IF (mpiHandler % mpiEnabled) THEN
      IF (resetCount) THEN
        msgCount = 0
      ELSE
        msgCount = mpiHandler % msgCount
      END IF

      DO e1 = 1,vector % nElem
        DO s1 = 1,4

          e2 = mesh % sideInfo % hostData(3,s1,e1) ! Neighbor Element
          IF( e2 > 0 )THEN
            r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

            IF (r2 /= mpiHandler % rankId) THEN

              s2 = mesh % sideInfo % hostData(4,s1,e1)/10
              globalSideId = ABS(mesh % sideInfo % hostdata(2,s1,e1))

              msgCount = msgCount + 1
              CALL MPI_IRECV(vector % extBoundary % hostData(:,:,:,s1,e1), &
                             2*(vector % interp % N + 1)*vector % nVar, &
                             mpiHandler % mpiPrec, &
                             r2,globalSideId, &
                             mpiHandler % mpiComm, &
                             mpiHandler % requests(msgCount),iError)

              msgCount = msgCount + 1
              CALL MPI_ISEND(vector % boundary % hostData(:,:,:,s1,e1), &
                             2*(vector % interp % N + 1)*vector % nVar, &
                             mpiHandler % mpiPrec, &
                             r2,globalSideId, &
                             mpiHandler % mpiComm, &
                             mpiHandler % requests(msgCount),iError)

            END IF
          ENDIF

        END DO
      END DO

      mpiHandler % msgCount = msgCount
    END IF

  END SUBROUTINE MPIExchangeAsync_MappedVector2D

  SUBROUTINE ApplyFlip_MappedVector2D(vector,mpiHandler,mesh,gpuAccel)
    ! Apply side flips to sides where MPI exchanges took place.
    IMPLICIT NONE
    CLASS(MappedVector2D),INTENT(inout) :: vector
    TYPE(MPILayer),INTENT(inout) :: mpiHandler
    TYPE(Mesh2D),INTENT(in) :: mesh
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1,s1,e2,s2
    INTEGER :: i,i2
    INTEGER :: r2,flip,ivar
    INTEGER :: globalSideId
    INTEGER :: bcid
    REAL(prec) :: extBuff(1:2,0:vector % interp % N)

    IF (mpiHandler % mpiEnabled) THEN
      IF (gpuAccel) THEN

        CALL ApplyFlip_MappedVector2D_gpu_wrapper(vector % extBoundary % deviceData, &
                                                  mesh % sideInfo % deviceData, &
                                                  mpiHandler % elemToRank % deviceData, &
                                                  mpiHandler % rankId, &
                                                  vector % interp % N, &
                                                  vector % nVar, &
                                                  vector % nElem)
      ELSE
        DO e1 = 1,vector % nElem
          DO s1 = 1,4

            e2 = mesh % sideInfo % hostData(3,s1,e1) ! Neighbor Element
            bcid = mesh % sideInfo % hostData(5,s1,e1)
            IF (bcid == 0) THEN ! Interior Element
              r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

              IF (r2 /= mpiHandler % rankId) THEN

                s2 = mesh % sideInfo % hostData(4,s1,e1)/10
                flip = mesh % sideInfo % hostData(4,s1,e1) - s2*10
                globalSideId = mesh % sideInfo % hostdata(2,s1,e1)

                ! Need to update extBoundary with flip applied
                IF (flip == 1) THEN

                  DO ivar = 1,vector % nvar
                    DO i = 0,vector % interp % N
                      i2 = vector % interp % N - i
                      extBuff(1:2,i) = vector % extBoundary % hostData(1:2,i2,ivar,s1,e1)
                    END DO
                    DO i = 0,vector % interp % N
                      vector % extBoundary % hostData(1:2,i,ivar,s1,e1) = extBuff(1:2,i)
                    END DO
                  END DO

                END IF
              END IF
            ENDIF

          END DO
        END DO
      END IF
    END IF

  END SUBROUTINE ApplyFlip_MappedVector2D

  SUBROUTINE MPIExchangeAsync_MappedScalar3D(scalar,mpiHandler,mesh,resetCount)
    IMPLICIT NONE
    CLASS(MappedScalar3D),INTENT(inout) :: scalar
    TYPE(MPILayer),INTENT(inout) :: mpiHandler
    TYPE(Mesh3D),INTENT(in) :: mesh
    LOGICAL,INTENT(in) :: resetCount
    ! Local
    INTEGER :: e1,s1,e2,s2
    INTEGER :: globalSideId,r2
    INTEGER :: iError
    INTEGER :: msgCount

    IF (mpiHandler % mpiEnabled) THEN
      IF (resetCount) THEN
        msgCount = 0
      ELSE
        msgCount = mpiHandler % msgCount
      END IF

      DO e1 = 1,scalar % nElem
        DO s1 = 1,6

          e2 = mesh % sideInfo % hostData(3,s1,e1) ! Neighbor Element
          IF( e2 > 0 )THEN
            r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

            IF (r2 /= mpiHandler % rankId) THEN

              s2 = mesh % sideInfo % hostData(4,s1,e1)/10
              globalSideId = ABS(mesh % sideInfo % hostdata(2,s1,e1))

              msgCount = msgCount + 1
              CALL MPI_IRECV(scalar % extBoundary % hostData(:,:,:,s1,e1), &
                             (scalar % interp % N + 1)*(scalar % interp % N + 1)*scalar % nVar, &
                             mpiHandler % mpiPrec, &
                             r2,globalSideId, &
                             mpiHandler % mpiComm, &
                             mpiHandler % requests(msgCount),iError)

              msgCount = msgCount + 1
              CALL MPI_ISEND(scalar % boundary % hostData(:,:,:,s1,e1), &
                             (scalar % interp % N + 1)*(scalar % interp % N + 1)*scalar % nVar, &
                             mpiHandler % mpiPrec, &
                             r2,globalSideId, &
                             mpiHandler % mpiComm, &
                             mpiHandler % requests(msgCount),iError)

            END IF

          ENDIF

        END DO
      END DO

      mpiHandler % msgCount = msgCount
    END IF

  END SUBROUTINE MPIExchangeAsync_MappedScalar3D
!
  SUBROUTINE ApplyFlip_MappedScalar3D(scalar,mpiHandler,mesh,gpuAccel)
    ! Apply side flips to sides where MPI exchanges took place.
    IMPLICIT NONE
    CLASS(MappedScalar3D),INTENT(inout) :: scalar
    TYPE(MPILayer),INTENT(inout) :: mpiHandler
    TYPE(Mesh3D),INTENT(in) :: mesh
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1,s1,e2,s2
    INTEGER :: i,i2,j,j2
    INTEGER :: r2,flip,ivar
    INTEGER :: globalSideId
    INTEGER :: bcid
    REAL(prec) :: extBuff(0:scalar % interp % N,0:scalar % interp % N)

    IF (mpiHandler % mpiEnabled) THEN
      IF (gpuAccel) THEN

        CALL ApplyFlip_MappedScalar3D_gpu_wrapper(scalar % extBoundary % deviceData, &
                                                  mesh % sideInfo % deviceData, &
                                                  mpiHandler % elemToRank % deviceData, &
                                                  mpiHandler % rankId, &
                                                  scalar % interp % N, &
                                                  scalar % nVar, &
                                                  scalar % nElem)
      ELSE
        DO e1 = 1,scalar % nElem
          DO s1 = 1,6

            e2 = mesh % sideInfo % hostData(3,s1,e1) ! Neighbor Element
            bcid = mesh % sideInfo % hostData(5,s1,e1)
            IF( bcid == 0 )THEN ! Interior Element
              r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

              IF (r2 /= mpiHandler % rankId) THEN

                s2 = mesh % sideInfo % hostData(4,s1,e1)/10
                flip = mesh % sideInfo % hostData(4,s1,e1) - s2*10
                globalSideId = mesh % sideInfo % hostdata(2,s1,e1)

                ! Need to update extBoundary with flip applied
                IF (flip == 1) THEN

                  DO ivar = 1,scalar % nvar
                    DO j = 0,scalar % interp % N
                      DO i = 0,scalar % interp % N
                        i2 = j
                        j2 = scalar % interp % N-i
                        extBuff(i,j) = scalar % extBoundary % hostData(i2,j2,ivar,s1,e1)
                      END DO
                    END DO
                    DO j = 0,scalar % interp % N
                      DO i = 0,scalar % interp % N
                        scalar % extBoundary % hostData(i,j,ivar,s1,e1) = extBuff(i,j)
                      END DO
                    END DO
                  END DO

                ELSEIF (flip == 2) THEN

                  DO ivar = 1,scalar % nvar
                    DO j = 0,scalar % interp % N
                      DO i = 0,scalar % interp % N
                        i2 = scalar % interp % N - i
                        j2 = scalar % interp % N - j
                        extBuff(i,j) = scalar % extBoundary % hostData(i2,j2,ivar,s1,e1)
                      END DO
                    END DO
                    DO j = 0,scalar % interp % N
                      DO i = 0,scalar % interp % N
                        scalar % extBoundary % hostData(i,j,ivar,s1,e1) = extBuff(i,j)
                      END DO
                    END DO
                  END DO

                ELSEIF (flip == 3) THEN

                  DO ivar = 1,scalar % nvar
                    DO j = 0,scalar % interp % N
                      DO i = 0,scalar % interp % N
                        i2 = scalar % interp % N-j
                        j2 = i
                        extBuff(i,j) = scalar % extBoundary % hostData(i2,j2,ivar,s1,e1)
                      END DO
                    END DO
                    DO j = 0,scalar % interp % N
                      DO i = 0,scalar % interp % N
                        scalar % extBoundary % hostData(i,j,ivar,s1,e1) = extBuff(i,j)
                      END DO
                    END DO
                  END DO

                ELSEIF (flip == 4) THEN

                  DO ivar = 1,scalar % nvar
                    DO j = 0,scalar % interp % N
                      DO i = 0,scalar % interp % N
                        i2 = j
                        j2 = i
                        extBuff(i,j) = scalar % extBoundary % hostData(i2,j2,ivar,s1,e1)
                      END DO
                    END DO
                    DO j = 0,scalar % interp % N
                      DO i = 0,scalar % interp % N
                        scalar % extBoundary % hostData(i,j,ivar,s1,e1) = extBuff(i,j)
                      END DO
                    END DO
                  END DO

                END IF
              END IF

            ENDIF

          END DO
        END DO
      END IF
    END IF

  END SUBROUTINE ApplyFlip_MappedScalar3D

  SUBROUTINE MPIExchangeAsync_MappedVector3D(vector,mpiHandler,mesh,resetCount)
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(inout) :: vector
    TYPE(MPILayer),INTENT(inout) :: mpiHandler
    TYPE(Mesh3D),INTENT(in) :: mesh
    LOGICAL,INTENT(in) :: resetCount
    ! Local
    INTEGER :: e1,s1,e2,s2
    INTEGER :: globalSideId,r2
    INTEGER :: iError
    INTEGER :: msgCount

    IF (mpiHandler % mpiEnabled) THEN
      IF (resetCount) THEN
        msgCount = 0
      ELSE
        msgCount = mpiHandler % msgCount
      END IF

      DO e1 = 1,vector % nElem
        DO s1 = 1,6

          e2 = mesh % sideInfo % hostData(3,s1,e1) ! Neighbor Element
          r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

          IF (r2 /= mpiHandler % rankId) THEN

            s2 = mesh % sideInfo % hostData(4,s1,e1)/10
            globalSideId = ABS(mesh % sideInfo % hostdata(2,s1,e1))

            msgCount = msgCount + 1
            CALL MPI_IRECV(vector % extBoundary % hostData(:,:,:,:,s1,e1), &
                           3*(vector % interp % N + 1)*(vector % interp % N + 1)*vector % nVar, &
                           mpiHandler % mpiPrec, &
                           r2,globalSideId, &
                           mpiHandler % mpiComm, &
                           mpiHandler % requests(msgCount),iError)

            msgCount = msgCount + 1
            CALL MPI_ISEND(vector % boundary % hostData(:,:,:,:,s1,e1), &
                           3*(vector % interp % N + 1)*(vector % interp % N + 1)*vector % nVar, &
                           mpiHandler % mpiPrec, &
                           r2,globalSideId, &
                           mpiHandler % mpiComm, &
                           mpiHandler % requests(msgCount),iError)
          END IF

        END DO
      END DO

      mpiHandler % msgCount = msgCount
    END IF

  END SUBROUTINE MPIExchangeAsync_MappedVector3D

  SUBROUTINE ApplyFlip_MappedVector3D(vector,mpiHandler,mesh,gpuAccel)
    ! Apply side flips to sides where MPI exchanges took place.
    IMPLICIT NONE
    CLASS(MappedVector3D),INTENT(inout) :: vector
    TYPE(MPILayer),INTENT(inout) :: mpiHandler
    TYPE(Mesh3D),INTENT(in) :: mesh
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: e1,s1,e2,s2
    INTEGER :: i,i2,j,j2
    INTEGER :: r2,flip,ivar
    INTEGER :: globalSideId
    INTEGER :: bcid     
    REAL(prec) :: extBuff(1:3,0:vector % interp % N,0:vector % interp % N)

    IF (mpiHandler % mpiEnabled) THEN
      IF (gpuAccel) THEN

        CALL ApplyFlip_MappedVector3D_gpu_wrapper(vector % extBoundary % deviceData, &
                                                  mesh % sideInfo % deviceData, &
                                                  mpiHandler % elemToRank % deviceData, &
                                                  mpiHandler % rankId, &
                                                  vector % interp % N, &
                                                  vector % nVar, &
                                                  vector % nElem)
      ELSE
        DO e1 = 1,vector % nElem
          DO s1 = 1,6

            e2 = mesh % sideInfo % hostData(3,s1,e1) ! Neighbor Element
            bcid = mesh % sideInfo % hostData(5,s1,e1)
            IF (bcid == 0) THEN ! Interior Element
              r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

              IF (r2 /= mpiHandler % rankId) THEN

                s2 = mesh % sideInfo % hostData(4,s1,e1)/10
                flip = mesh % sideInfo % hostData(4,s1,e1) - s2*10
                globalSideId = mesh % sideInfo % hostdata(2,s1,e1)

                IF (flip == 1) THEN

                  DO ivar = 1,vector % nvar
                    DO j = 0,vector % interp % N
                      DO i = 0,vector % interp % N
                        i2 = j
                        j2 = vector % interp % N - i
                        extBuff(1:3,i,j) = vector % extBoundary % hostData(1:3,i2,j2,ivar,s1,e1)
                      END DO
                    END DO
                    DO j = 0,vector % interp % N
                      DO i = 0,vector % interp % N
                        vector % extBoundary % hostData(1:3,i,j,ivar,s1,e1) = extBuff(1:3,i,j)
                      END DO
                    END DO
                  END DO

                ELSEIF (flip == 2) THEN

                  DO ivar = 1,vector % nvar
                    DO j = 0,vector % interp % N
                      DO i = 0,vector % interp % N
                        i2 = vector % interp % N - i
                        j2 = vector % interp % N - j
                        extBuff(1:3,i,j) = vector % extBoundary % hostData(1:3,i2,j2,ivar,s1,e1)
                      END DO
                    END DO
                    DO j = 0,vector % interp % N
                      DO i = 0,vector % interp % N
                        vector % extBoundary % hostData(1:3,i,j,ivar,s1,e1) = extBuff(1:3,i,j)
                      END DO
                    END DO
                  END DO

                ELSEIF (flip == 3) THEN

                  DO ivar = 1,vector % nvar
                    DO j = 0,vector % interp % N
                      DO i = 0,vector % interp % N
                        i2 = vector % interp % N - j
                        j2 = i
                        extBuff(1:3,i,j) = vector % extBoundary % hostData(1:3,i2,j2,ivar,s1,e1)
                      END DO
                    END DO
                    DO j = 0,vector % interp % N
                      DO i = 0,vector % interp % N
                        vector % extBoundary % hostData(1:3,i,j,ivar,s1,e1) = extBuff(1:3,i,j)
                      END DO
                    END DO
                  END DO

                ELSEIF (flip == 4) THEN

                  DO ivar = 1,vector % nvar
                    DO j = 0,vector % interp % N
                      DO i = 0,vector % interp % N
                        i2 = j
                        j2 = i
                        extBuff(1:3,i,j) = vector % extBoundary % hostData(1:3,i2,j2,ivar,s1,e1)
                      END DO
                    END DO
                    DO j = 0,vector % interp % N
                      DO i = 0,vector % interp % N
                        vector % extBoundary % hostData(1:3,i,j,ivar,s1,e1) = extBuff(1:3,i,j)
                      END DO
                    END DO
                  END DO

                END IF
              END IF
            ENDIF
          END DO
        END DO
      END IF
    END IF

  END SUBROUTINE ApplyFlip_MappedVector3D

!   ! ---------------------- Two point Vectors ---------------------- !

!   SUBROUTINE SideExchange_MappedP2Vector2D(vector,mesh,decomp,gpuAccel)
!   !! SideExchange_MappedP2Vectorvector2D is used to populate vector % extBoundary
!   !! by finding neighboring elements that share a side and copying the neighboring
!   !! elements solution % boundary data.
!     IMPLICIT NONE
!     CLASS(MappedP2Vector2D),INTENT(inout) :: vector
!     TYPE(Mesh2D),INTENT(in) :: mesh
!     TYPE(MPILayer),INTENT(inout) :: decomp
!     LOGICAL,INTENT(in) :: gpuAccel
!     ! Local
!     INTEGER :: e1,e2,s1,s2,e2Global
!     INTEGER :: flip,bcid
!     INTEGER :: neighborRank
!     INTEGER :: i1,i2,ivar
!     INTEGER :: rankId, offset

!       rankId = decomp % rankId
!       offset = decomp % offsetElem % hostData(rankId)

!     IF (gpuAccel) THEN

!       CALL vector % boundary % UpdateHost()
!       CALL vector % MPIExchangeAsync(decomp,mesh,resetCount=.TRUE.)
!       CALL decomp % FinalizeMPIExchangeAsync()
!       CALL vector % extBoundary % UpdateDevice()

!       CALL SideExchange_MappedVector2D_gpu_wrapper(vector % extBoundary % deviceData, &
!                                                    vector % boundary % deviceData, &
!                                                    mesh % sideInfo % deviceData, &
!                                                    decomp % elemToRank % deviceData, &
!                                                    decomp % rankId, &
!                                                    offset, &
!                                                    vector % interp % N, &
!                                                    vector % nvar, &
!                                                    vector % nElem)

!     ELSE

!       CALL vector % MPIExchangeAsync(decomp,mesh,resetCount=.TRUE.)

!       DO e1 = 1,mesh % nElem
!         DO s1 = 1,4
!           e2Global = mesh % sideInfo % hostData(3,s1,e1)
!           e2 = e2Global - offset
!           s2 = mesh % sideInfo % hostData(4,s1,e1)/10
!           flip = mesh % sideInfo % hostData(4,s1,e1) - s2*10
!           bcid = mesh % sideInfo % hostData(5,s1,e1)

!           IF (bcid == 0) THEN

!             neighborRank = decomp % elemToRank % hostData(e2Global)

!             IF (neighborRank == decomp % rankId) THEN

!               IF (flip == 0) THEN

!                 DO ivar = 1,vector % nvar
!                   DO i1 = 0,vector % interp % N
!                     vector % extBoundary % hostData(1:2,i1,ivar,s1,e1) = &
!                       vector % boundary % hostData(1:2,i1,ivar,s2,e2)
!                   END DO
!                 END DO

!               ELSEIF (flip == 1) THEN

!                 DO ivar = 1,vector % nvar
!                   DO i1 = 0,vector % interp % N
!                     i2 = vector % interp % N - i1
!                     vector % extBoundary % hostData(1:2,i1,ivar,s1,e1) = &
!                       vector % boundary % hostData(1:2,i2,ivar,s2,e2)
!                   END DO
!                 END DO

!               END IF

!             END IF

!           END IF

!         END DO
!       END DO

!       CALL decomp % FinalizeMPIExchangeAsync()

!     END IF

!     CALL vector % ApplyFlip(decomp,mesh,gpuAccel)

!   END SUBROUTINE SideExchange_MappedP2Vector2D

!   SUBROUTINE Divergence_MappedP2Vector2D(compVector,geometry,divVector,dForm,gpuAccel)
!     ! Strong Form Operator
!     !
!     ! DG Weak Form Operator
!     !
!     ! Assumes vector has been projected to computational coordinates
!     !
!     IMPLICIT NONE
!     CLASS(MappedP2Vector2D),INTENT(in) :: compVector
!     TYPE(SEMQuad),INTENT(in) :: geometry
!     TYPE(MappedScalar2D),INTENT(inout) :: divVector
!     INTEGER,INTENT(in) :: dForm
!     LOGICAL,INTENT(in) :: gpuAccel

!     IF (dForm == selfWeakDGForm) THEN

!       IF (gpuAccel) THEN
!         CALL compVector % interp % P2VectorDGDivergence_2D(compVector % interior % deviceData, &
!                                                          compVector % boundaryNormal % deviceData, &
!                                                          divVector % interior % deviceData, &
!                                                          compVector % nvar, &
!                                                          compVector % nelem)
!       ELSE
!         CALL compVector % interp % P2VectorDGDivergence_2D(compVector % interior % hostData, &
!                                                          compVector % boundaryNormal % hostData, &
!                                                          divVector % interior % hostData, &
!                                                          compVector % nvar, &
!                                                          compVector % nelem)
!       END IF

!     ELSE IF (dForm == selfStrongForm) THEN

!       IF (gpuAccel) THEN
!         CALL compVector % interp % P2VectorDivergence_2D(compVector % interior % deviceData, &
!                                                        divVector % interior % deviceData, &
!                                                        compVector % nvar, &
!                                                        compVector % nelem)
!       ELSE
!         CALL compVector % interp % P2VectorDivergence_2D(compVector % interior % hostData, &
!                                                        divVector % interior % hostData, &
!                                                        compVector % nvar, &
!                                                        compVector % nelem)
!       END IF

!     END IF

!     CALL divVector % JacobianWeight(geometry,gpuAccel)

!   END SUBROUTINE Divergence_MappedP2Vector2D

!   SUBROUTINE ContravariantProjection_MappedP2Vector2D(vector,geometry,gpuAccel)
! #undef __FUNC__
! #define __FUNC__ "ContravariantProjection_MappedP2Vector2D"
!     ! Takes a vector that has physical space coordinate directions (x,y,z) and projects the vector
!     ! into the the contravariant basis vector directions. Keep in mind that the contravariant basis
!     ! vectors are really the Jacobian weighted contravariant basis vectors
!     IMPLICIT NONE
!     CLASS(MappedP2Vector2D),INTENT(inout) :: vector
!     TYPE(SEMQuad),INTENT(in) :: geometry
!     LOGICAL,INTENT(in) :: gpuAccel
!     ! Local
!     INTEGER :: i,j,n,ivar,iel
!     REAL(prec) :: Fx, Fy

!     IF (gpuAccel) THEN

!       CALL ContravariantProjection_MappedP2Vector2D_gpu_wrapper(vector % interior % deviceData, &
!                                                               vector % physical % deviceData, &
!                                                               geometry % dsdx % interior % deviceData, &
!                                                               vector % interp % N, &
!                                                               vector % nVar, &
!                                                               vector % nElem)

!     ELSE
!       ! Assume that tensor(j,i) is vector i, component j
!       ! => dot product is done along first dimension
!       ! to project onto computational space
!       DO iel = 1,vector % nElem
!         DO ivar = 1,vector % nVar
!           DO j = 0,vector % interp % N
!             DO i = 0,vector % interp % N

!               ! From Winters et al. 2020, Kopriva and Gassner 2014, and Kopriva et al. 2019,  we use two point averaging of the
!               ! metric terms for dealiasing
!               ! > See pages 60-62 of "Construction of Modern Robust Nodal Discontinuous Galerkin Spectral Element Methods for 
!               !   the Compressible Navier-Stokes Equations", Winters et al. 2020
!               DO n = 0, vector % interp % N

!                 ! I think we need another attribute here, where
!                 ! two point values are stored for each Fx, Fy
!                 ! for each computational dimension
!                 ! Fx_{(i,n),j}, Fx_{i,(j,n)}
!                 ! Fy_{(i,n),j}, Fy_{i,(j,n)}

!                 ! Fx_{(i,n),j}
!                 Fx = vector % physical % hostData(1,1,n,i,j,ivar,iel)
!                 ! Fy_{(i,n),j}
!                 Fy = vector % physical % hostData(2,1,n,i,j,ivar,iel)

!                 vector % interior % hostData(1,n,i,j,ivar,iel) = &
!                   0.5_prec*( geometry % dsdx % interior % hostData(1,1,i,j,1,iel) + &
!                              geometry % dsdx % interior % hostData(1,1,n,j,1,iel) )*Fx + &
!                   0.5_prec*( geometry % dsdx % interior % hostData(2,1,i,j,1,iel) + &
!                              geometry % dsdx % interior % hostData(2,1,n,j,1,iel) )*Fy

!                 ! Fx_{i,(j,n)}
!                 Fx = vector % physical % hostData(1,2,n,i,j,ivar,iel)
!                 ! Fy_{i,(j,n)}
!                 Fy = vector % physical % hostData(2,2,n,i,j,ivar,iel)
!                 vector % interior % hostData(2,n,i,j,ivar,iel) = &
!                   0.5_prec*( geometry % dsdx % interior % hostData(1,2,i,j,1,iel) + &
!                              geometry % dsdx % interior % hostData(1,2,i,n,1,iel) )*Fx + &
!                   0.5_prec*( geometry % dsdx % interior % hostData(2,2,i,j,1,iel) + &
!                              geometry % dsdx % interior % hostData(2,2,i,n,1,iel) )*Fy

!               ENDDO

!             END DO
!           END DO
!         END DO
!       END DO

!     END IF

!   END SUBROUTINE ContravariantProjection_MappedP2Vector2D

!   SUBROUTINE MPIExchangeAsync_MappedP2Vector2D(vector,mpiHandler,mesh,resetCount)
!     IMPLICIT NONE
!     CLASS(MappedP2Vector2D),INTENT(inout) :: vector
!     TYPE(MPILayer),INTENT(inout) :: mpiHandler
!     TYPE(Mesh2D),INTENT(in) :: mesh
!     LOGICAL,INTENT(in) :: resetCount
!     ! Local
!     INTEGER :: e1,s1,e2,s2
!     INTEGER :: globalSideId,r2
!     INTEGER :: iError
!     INTEGER :: msgCount

!     IF (mpiHandler % mpiEnabled) THEN
!       IF (resetCount) THEN
!         msgCount = 0
!       ELSE
!         msgCount = mpiHandler % msgCount
!       END IF

!       DO e1 = 1,vector % nElem
!         DO s1 = 1,4

!           e2 = mesh % sideInfo % hostData(3,s1,e1) ! Neighbor Element
!           IF( e2 > 0 )THEN
!             r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

!             IF (r2 /= mpiHandler % rankId) THEN

!               s2 = mesh % sideInfo % hostData(4,s1,e1)/10
!               globalSideId = ABS(mesh % sideInfo % hostdata(2,s1,e1))

!               msgCount = msgCount + 1
!               CALL MPI_IRECV(vector % extBoundary % hostData(:,:,:,s1,e1), &
!                              2*(vector % interp % N + 1)*vector % nVar, &
!                              mpiHandler % mpiPrec, &
!                              r2,globalSideId, &
!                              mpiHandler % mpiComm, &
!                              mpiHandler % requests(msgCount),iError)

!               msgCount = msgCount + 1
!               CALL MPI_ISEND(vector % boundary % hostData(:,:,:,s1,e1), &
!                              2*(vector % interp % N + 1)*vector % nVar, &
!                              mpiHandler % mpiPrec, &
!                              r2,globalSideId, &
!                              mpiHandler % mpiComm, &
!                              mpiHandler % requests(msgCount),iError)

!             END IF
!           ENDIF

!         END DO
!       END DO

!       mpiHandler % msgCount = msgCount
!     END IF

!   END SUBROUTINE MPIExchangeAsync_MappedP2Vector2D

!   SUBROUTINE ApplyFlip_MappedP2Vector2D(vector,mpiHandler,mesh,gpuAccel)
!     ! Apply side flips to sides where MPI exchanges took place.
!     IMPLICIT NONE
!     CLASS(MappedP2Vector2D),INTENT(inout) :: vector
!     TYPE(MPILayer),INTENT(inout) :: mpiHandler
!     TYPE(Mesh2D),INTENT(in) :: mesh
!     LOGICAL,INTENT(in) :: gpuAccel
!     ! Local
!     INTEGER :: e1,s1,e2,s2
!     INTEGER :: i,i2
!     INTEGER :: r2,flip,ivar
!     INTEGER :: globalSideId
!     INTEGER :: bcid
!     REAL(prec) :: extBuff(1:2,0:vector % interp % N)

!     IF (mpiHandler % mpiEnabled) THEN
!       IF (gpuAccel) THEN

!         ! Since the boundary data for a p2 vector and a vector are identical,
!         ! we can reuse the applyFlip method for MappedVector here
!         CALL ApplyFlip_MappedVector2D_gpu_wrapper(vector % extBoundary % deviceData, &
!                                                   mesh % sideInfo % deviceData, &
!                                                   mpiHandler % elemToRank % deviceData, &
!                                                   mpiHandler % rankId, &
!                                                   vector % interp % N, &
!                                                   vector % nVar, &
!                                                   vector % nElem)
!       ELSE
!         DO e1 = 1,vector % nElem
!           DO s1 = 1,4

!             e2 = mesh % sideInfo % hostData(3,s1,e1) ! Neighbor Element
!             bcid = mesh % sideInfo % hostData(5,s1,e1)
!             IF (bcid == 0) THEN ! Interior Element
!               r2 = mpiHandler % elemToRank % hostData(e2) ! Neighbor Rank

!               IF (r2 /= mpiHandler % rankId) THEN

!                 s2 = mesh % sideInfo % hostData(4,s1,e1)/10
!                 flip = mesh % sideInfo % hostData(4,s1,e1) - s2*10
!                 globalSideId = mesh % sideInfo % hostdata(2,s1,e1)

!                 ! Need to update extBoundary with flip applied
!                 IF (flip == 1) THEN

!                   DO ivar = 1,vector % nvar
!                     DO i = 0,vector % interp % N
!                       i2 = vector % interp % N - i
!                       extBuff(1:2,i) = vector % extBoundary % hostData(1:2,i2,ivar,s1,e1)
!                     END DO
!                     DO i = 0,vector % interp % N
!                       vector % extBoundary % hostData(1:2,i,ivar,s1,e1) = extBuff(1:2,i)
!                     END DO
!                   END DO

!                 END IF
!               END IF
!             ENDIF

!           END DO
!         END DO
!       END IF
!     END IF

!   END SUBROUTINE ApplyFlip_MappedP2Vector2D

END MODULE SELF_MappedData
