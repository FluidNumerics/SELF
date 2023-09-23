! SELF_Data.F90
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE SELF_Data

  USE SELF_Constants
  USE SELF_Lagrange
  USE SELF_Metadata
  USE FEQParse

  USE ISO_C_BINDING

  IMPLICIT NONE

#include "SELF_Macros.h"

  TYPE,PUBLIC :: SELF_DataObj
  !! The SELF_DataObj class is a base class for all data objects in SELF.
  !! A data object in SELF is a multidimensional array of data, represented
  !! on both host and device, that is associated with an interpolant, metadata,
  !! and (optionally) an equation string.
  !! Type extensions of the SELF_DataObj include scalars, vectors, and tensors
  !! in 1-D, 2-D, and 3-D using the storage patterns that are expected for
  !! derivative and interpolation operations defined in SELF_Lagrange.f90
  !! Additionally, each extended type has the necessary attributes to store
  !! information on element interiors and element boundaries, both of which
  !! are commonly used for spectral element solvers.

    INTEGER :: nVar
    INTEGER :: nElem
    TYPE(Lagrange),POINTER :: interp
    TYPE(Metadata),ALLOCATABLE :: meta(:)
    TYPE(EquationParser),ALLOCATABLE :: eqn(:)

    CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_DataObj
    PROCEDURE,PUBLIC :: Free => Free_DataObj

    ! Procedures for setting metadata for 
    PROCEDURE,PUBLIC :: SetName => SetName_DataObj
    PROCEDURE,PUBLIC :: SetDescription => SetDescription_DataObj
    PROCEDURE,PUBLIC :: SetUnits => SetUnits_DataObj
    GENERIC,PUBLIC :: SetEquation => SetEquation_DataObj
    PROCEDURE,PRIVATE :: SetEquation_DataObj

  END TYPE SELF_DataObj

! ---------------------- Scalars ---------------------- !
  TYPE,EXTENDS(SELF_DataObj),PUBLIC :: Scalar1D

    TYPE(hfReal_r3) :: interior
    TYPE(hfReal_r3) :: boundary
    TYPE(hfReal_r3) :: extBoundary
    TYPE(hfReal_r3) :: avgBoundary
    TYPE(hfReal_r3) :: jumpBoundary

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Scalar1D
    PROCEDURE,PUBLIC :: Free => Free_Scalar1D

    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Scalar1D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Scalar1D

    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Scalar1D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Scalar1D

    GENERIC,PUBLIC :: Derivative => Derivative_Scalar1D
    PROCEDURE,PRIVATE :: Derivative_Scalar1D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Scalar1D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Scalar1D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Scalar1D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Scalar1D

    GENERIC,PUBLIC :: WriteHDF5 => WriteHDF5_MPI_Scalar1D, WriteHDF5_Scalar1D
    PROCEDURE, PRIVATE :: WriteHDF5_MPI_Scalar1D
    PROCEDURE, PRIVATE :: WriteHDF5_Scalar1D


  END TYPE Scalar1D

  TYPE,EXTENDS(SELF_DataObj),PUBLIC :: Scalar2D

    TYPE(hfReal_r4) :: interior
    TYPE(hfReal_r4) :: boundary
    TYPE(hfReal_r4) :: extBoundary
    TYPE(hfReal_r4) :: avgBoundary
    TYPE(hfReal_r4) :: jumpBoundary

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Scalar2D
    PROCEDURE,PUBLIC :: Free => Free_Scalar2D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Scalar2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Scalar2D
    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Scalar2D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Scalar2D

    GENERIC,PUBLIC :: Gradient => Gradient_Scalar2D
    PROCEDURE,PRIVATE :: Gradient_Scalar2D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Scalar2D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Scalar2D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Scalar2D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Scalar2D

    GENERIC,PUBLIC :: WriteHDF5 => WriteHDF5_MPI_Scalar2D, WriteHDF5_Scalar2D
    PROCEDURE, PRIVATE :: WriteHDF5_MPI_Scalar2D
    PROCEDURE, PRIVATE :: WriteHDF5_Scalar2D

  END TYPE Scalar2D

  TYPE,EXTENDS(SELF_DataObj),PUBLIC :: Scalar3D

    TYPE(hfReal_r5) :: interior
    TYPE(hfReal_r5) :: boundary
    TYPE(hfReal_r5) :: extBoundary
    TYPE(hfReal_r5) :: avgBoundary
    TYPE(hfReal_r5) :: jumpBoundary

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Scalar3D
    PROCEDURE,PUBLIC :: Free => Free_Scalar3D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Scalar3D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Scalar3D
    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Scalar3D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Scalar3D

    GENERIC,PUBLIC :: Gradient => Gradient_Scalar3D
    PROCEDURE,PRIVATE :: Gradient_Scalar3D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Scalar3D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Scalar3D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Scalar3D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Scalar3D

  END TYPE Scalar3D

! ---------------------- Vectors ---------------------- !

  TYPE,EXTENDS(SELF_DataObj),PUBLIC :: Vector2D

    TYPE(hfReal_r5) :: interior
    TYPE(hfReal_r5) :: boundary
    TYPE(hfReal_r5) :: extBoundary
    TYPE(hfReal_r4) :: boundaryNormal

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Vector2D
    PROCEDURE,PUBLIC :: Free => Free_Vector2D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Vector2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Vector2D
    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Vector2D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Vector2D

    GENERIC,PUBLIC :: Gradient => Gradient_Vector2D
    PROCEDURE,PRIVATE :: Gradient_Vector2D

    GENERIC,PUBLIC :: Divergence => Divergence_Vector2D
    PROCEDURE,PRIVATE :: Divergence_Vector2D

    GENERIC,PUBLIC :: Curl => Curl_Vector2D
    PROCEDURE,PRIVATE :: Curl_Vector2D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Vector2D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Vector2D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Vector2D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Vector2D

    GENERIC,PUBLIC :: SetEquation => SetEquation_Vector2D
    PROCEDURE,PRIVATE :: SetEquation_Vector2D

    GENERIC,PUBLIC :: WriteHDF5 => WriteHDF5_MPI_Vector2D, WriteHDF5_Vector2D
    PROCEDURE, PRIVATE :: WriteHDF5_MPI_Vector2D
    PROCEDURE, PRIVATE :: WriteHDF5_Vector2D

  END TYPE Vector2D

  TYPE,EXTENDS(SELF_DataObj),PUBLIC :: Vector3D

    TYPE(hfReal_r6) :: interior
    TYPE(hfReal_r6) :: boundary
    TYPE(hfReal_r6) :: extBoundary
    TYPE(hfReal_r5) :: boundaryNormal

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Vector3D
    PROCEDURE,PUBLIC :: Free => Free_Vector3D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Vector3D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Vector3D
    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Vector3D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Vector3D

    GENERIC,PUBLIC :: Gradient => Gradient_Vector3D
    PROCEDURE,PRIVATE :: Gradient_Vector3D

    GENERIC,PUBLIC :: Divergence => Divergence_Vector3D
    PROCEDURE,PRIVATE :: Divergence_Vector3D

    GENERIC,PUBLIC :: Curl => Curl_Vector3D
    PROCEDURE,PRIVATE :: Curl_Vector3D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Vector3D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Vector3D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Vector3D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Vector3D

    GENERIC,PUBLIC :: SetEquation => SetEquation_Vector3D
    PROCEDURE,PRIVATE :: SetEquation_Vector3D

  END TYPE Vector3D

! ----------------- Two-point Vectors ----------------- !

  TYPE,EXTENDS(SELF_DataObj),PUBLIC :: P2Vector2D

    TYPE(hfReal_r6) :: interior
    TYPE(hfReal_r7) :: physical
    TYPE(hfReal_r5) :: boundary
    TYPE(hfReal_r5) :: extBoundary
    TYPE(hfReal_r4) :: boundaryNormal

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_P2Vector2D
    PROCEDURE,PUBLIC :: Free => Free_P2Vector2D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_P2Vector2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_P2Vector2D

    GENERIC,PUBLIC :: Divergence => Divergence_P2Vector2D
    PROCEDURE,PRIVATE :: Divergence_P2Vector2D

  END TYPE P2Vector2D
! ---------------------- Tensors ---------------------- !

  TYPE,EXTENDS(SELF_DataObj),PUBLIC :: Tensor2D

    TYPE(hfReal_r6) :: interior
    TYPE(hfReal_r6) :: boundary
    TYPE(hfReal_r6) :: extBoundary

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Tensor2D
    PROCEDURE,PUBLIC :: Free => Free_Tensor2D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Tensor2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Tensor2D
    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Tensor2D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Tensor2D

    PROCEDURE,PUBLIC :: Determinant => Determinant_Tensor2D

    GENERIC,PUBLIC :: Divergence => Divergence_Tensor2D
    PROCEDURE,PRIVATE :: Divergence_Tensor2D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Tensor2D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Tensor2D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Tensor2D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Tensor2D

    GENERIC,PUBLIC :: SetEquation => SetEquation_Tensor2D
    PROCEDURE,PRIVATE :: SetEquation_Tensor2D

  END TYPE Tensor2D

  TYPE,EXTENDS(SELF_DataObj),PUBLIC :: Tensor3D

    TYPE(hfReal_r7) :: interior
    TYPE(hfReal_r7) :: boundary
    TYPE(hfReal_r7) :: extBoundary

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Tensor3D
    PROCEDURE,PUBLIC :: Free => Free_Tensor3D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Tensor3D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Tensor3D
    PROCEDURE,PUBLIC :: BoundaryInterp => BoundaryInterp_Tensor3D
    PROCEDURE,PUBLIC :: GridInterp => GridInterp_Tensor3D

    PROCEDURE,PUBLIC :: Determinant => Determinant_Tensor3D

    GENERIC,PUBLIC :: Divergence => Divergence_Tensor3D
    PROCEDURE,PRIVATE :: Divergence_Tensor3D

    PROCEDURE,PUBLIC :: AbsMaxInterior => AbsMaxInterior_Tensor3D
    PROCEDURE,PUBLIC :: AbsMaxBoundary => AbsMaxBoundary_Tensor3D

    GENERIC,PUBLIC :: ASSIGNMENT(=) => Equals_Tensor3D
    PROCEDURE,PRIVATE,PASS(SELFOut) :: Equals_Tensor3D

    GENERIC,PUBLIC :: SetEquation => SetEquation_Tensor3D
    PROCEDURE,PRIVATE :: SetEquation_Tensor3D

  END TYPE Tensor3D

  INTEGER,PARAMETER :: selfStrongForm = 0
  INTEGER,PARAMETER :: selfWeakDGForm = 1
  INTEGER,PARAMETER :: selfWeakCGForm = 2
  INTEGER,PARAMETER :: selfWeakBRForm = 3

  INTERFACE
    SUBROUTINE Determinant_Tensor2D_gpu_wrapper(tensor_dev,detTensor_dev,N,nVar,nEl) &
      bind(c,name="Determinant_Tensor2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: tensor_dev,detTensor_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE Determinant_Tensor2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE Determinant_Tensor3D_gpu_wrapper(tensor_dev,detTensor_dev,N,nVar,nEl) &
      bind(c,name="Determinant_Tensor3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: tensor_dev,detTensor_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE Determinant_Tensor3D_gpu_wrapper
  END INTERFACE

CONTAINS

! -- DataObj -- !
  SUBROUTINE Init_DataObj(SELFStorage,interp,nVar,nElem)
    IMPLICIT NONE
    CLASS(SELF_DataObj),INTENT(out)  :: SELFStorage
    TYPE(Lagrange),INTENT(in),TARGET :: interp
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem

    SELFStorage % nElem = nElem
    SELFStorage % nVar = nVar
    SELFStorage % interp => interp
    ALLOCATE( SELFStorage % meta(1:nVar) )
    ALLOCATE( SELFStorage % eqn(1:nVar) )

  END SUBROUTINE Init_DataObj

  SUBROUTINE Free_DataObj(SELFStorage)
    IMPLICIT NONE
    CLASS(SELF_DataObj),INTENT(inout) :: SELFStorage

    SELFStorage % interp => NULL()
    DEALLOCATE( SELFStorage % meta )
    DEALLOCATE( SELFStorage % eqn )

  END SUBROUTINE Free_DataObj

  SUBROUTINE SetName_DataObj(SELFStorage,ivar,name)
    !! Set the name of the `ivar-th` variable
    IMPLICIT NONE
    CLASS(SELF_DataObj),INTENT(inout) :: SELFStorage
    INTEGER,INTENT(in) :: ivar
    CHARACTER(*),INTENT(in) :: name

    CALL SELFStorage % meta(ivar) % SetName(name) 

  END SUBROUTINE SetName_DataObj

  SUBROUTINE SetDescription_DataObj(SELFStorage,ivar,description)
    !! Set the description of the `ivar-th` variable
    IMPLICIT NONE
    CLASS(SELF_DataObj),INTENT(inout) :: SELFStorage
    INTEGER,INTENT(in) :: ivar
    CHARACTER(*),INTENT(in) :: description

    CALL SELFStorage % meta(ivar) % SetDescription(description) 

  END SUBROUTINE SetDescription_DataObj

  SUBROUTINE SetUnits_DataObj(SELFStorage,ivar,units)
    !! Set the units of the `ivar-th` variable
    IMPLICIT NONE
    CLASS(SELF_DataObj),INTENT(inout) :: SELFStorage
    INTEGER,INTENT(in) :: ivar
    CHARACTER(*),INTENT(in) :: units

    CALL SELFStorage % meta(ivar) % SetUnits(units) 

  END SUBROUTINE SetUnits_DataObj

  SUBROUTINE SetEquation_DataObj(SELFStorage,ivar,eqnChar)
    !! Sets the equation parser for the `ivar-th` variable
    IMPLICIT NONE
    CLASS(SELF_DataObj),INTENT(inout) :: SELFStorage
    INTEGER,INTENT(in) :: ivar
    CHARACTER(*),INTENT(in) :: eqnChar

    SELFStorage % eqn(ivar) = EquationParser( TRIM(eqnChar), &
                                              (/'x','y','z','t'/) )

  END SUBROUTINE SetEquation_DataObj

! -- Scalar1D -- !

  SUBROUTINE Init_Scalar1D(SELFStorage,interp,nVar,nElem)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(out) :: SELFStorage
    TYPE(Lagrange),INTENT(in),TARGET :: interp
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem

    SELFStorage % interp => interp
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem

    CALL SELFStorage % interior % Alloc(loBound=(/0,1,1/), &
                                        upBound=(/interp % N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/1,1,1/), &
                                        upBound=(/nVar,2,nElem/))

    CALL SELFStorage % extBoundary % Alloc(loBound=(/1,1,1/), &
                                           upBound=(/nVar,2,nElem/))

    CALL SELFStorage % avgBoundary % Alloc(loBound=(/1,1,1/), &
                                           upBound=(/nVar,2,nElem/))

    CALL SELFStorage % jumpBoundary % Alloc(loBound=(/1,1,1/), &
                                           upBound=(/nVar,2,nElem/))

    ALLOCATE( SELFStorage % meta(1:nVar) )
    ALLOCATE( SELFStorage % eqn(1:nVar) )

  END SUBROUTINE Init_Scalar1D

  SUBROUTINE Free_Scalar1D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(inout) :: SELFStorage

    SELFStorage % interp => NULL()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()
    CALL SELFStorage % extBoundary % Free()
    CALL SELFStorage % avgBoundary % Free()
    CALL SELFStorage % jumpBoundary % Free()
    DEALLOCATE( SELFStorage % meta )
    DEALLOCATE( SELFStorage % eqn )

  END SUBROUTINE Free_Scalar1D

  SUBROUTINE UpdateHost_Scalar1D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
    CALL SELFStorage % extBoundary % UpdateHost()
    CALL SELFStorage % avgBoundary % UpdateHost()
    CALL SELFStorage % jumpBoundary % UpdateHost()

  END SUBROUTINE UpdateHost_Scalar1D

  SUBROUTINE UpdateDevice_Scalar1D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
    CALL SELFStorage % extBoundary % UpdateDevice()
    CALL SELFStorage % avgBoundary % UpdateDevice()
    CALL SELFStorage % jumpBoundary % UpdateDevice()

  END SUBROUTINE UpdateDevice_Scalar1D

  SUBROUTINE BoundaryInterp_Scalar1D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarBoundaryInterp_1D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarBoundaryInterp_1D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Scalar1D

  SUBROUTINE GridInterp_Scalar1D(SELFStorage,SELFout,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(in) :: SELFStorage
    TYPE(Scalar1D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarGridInterp_1D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarGridInterp_1D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Scalar1D

  SUBROUTINE Derivative_Scalar1D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(in) :: SELFStorage
    TYPE(Scalar1D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % Derivative_1D(SELFStorage % interior % deviceData, &
                                                SELFout % interior % deviceData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % Derivative_1D(SELFStorage % interior % hostData, &
                                                SELFout % interior % hostData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem)
    END IF

  END SUBROUTINE Derivative_Scalar1D

  FUNCTION AbsMaxInterior_Scalar1D(scalar) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Scalar1D) :: scalar
    REAL(prec) :: absMax(1:scalar % nVar)
    ! Local
    INTEGER :: iEl,iVar,i

    absMax = 0.0_prec
    DO iEl = 1,scalar % nElem
      DO iVar = 1,scalar % nVar
        DO i = 0,scalar % interp % N
          absMax(iVar) = MAX(ABS(scalar % interior % hostData(i,iVar,iEl)),absMax(iVar))
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Scalar1D

  FUNCTION AbsMaxBoundary_Scalar1D(scalar) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Scalar1D) :: scalar
    REAL(prec) :: absMax(1:scalar % nVar,1:2)
    ! Local
    INTEGER :: iEl,iVar,iSide

    absMax = 0.0_prec
    DO iEl = 1,scalar % nElem
      DO iSide = 1,2
        DO iVar = 1,scalar % nVar
          absMax(iVar,iSide) = MAX(ABS(scalar % boundary % hostData(iVar,iSide,iEl)),absMax(iVar,iSide))
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Scalar1D

  SUBROUTINE Equals_Scalar1D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Scalar1D),INTENT(inout) :: SELFOut
    TYPE(Scalar1D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Scalar1D

  SUBROUTINE WriteHDF5_MPI_Scalar1D(this,fileId,group,elemoffset,nglobalelem)
    IMPLICIT NONE
    CLASS(Scalar1D), INTENT(in) :: this
    CHARACTER(*), INTENT(in) :: group
    INTEGER(HID_T), INTENT(in) :: fileId
    INTEGER, INTENT(in) :: elemoffset
    INTEGER, INTENT(in) :: nglobalelem
    ! Local
    INTEGER(HID_T) :: offset(1:3)
    INTEGER(HID_T) :: bOffset(1:3)
    INTEGER(HID_T) :: globalDims(1:3)
    INTEGER(HID_T) :: bGlobalDims(1:3)
    INTEGER :: ivar

      offset(1:3) = (/0,0,elemoffset/)
      globalDims(1:3) = (/this % interp % N + 1, &
                          this % nVar, &
                          nGlobalElem/)

      ! Offsets and dimensions for element boundary data
      bOffset(1:3) = (/0,0,elemoffset/)
      bGlobalDims(1:3) = (/this % nVar, &
                           2, &
                           nGlobalElem/)
  
      CALL CreateGroup_HDF5(fileId,TRIM(group))

      DO ivar = 1, this % nVar
        CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
      ENDDO 

      CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
                           this % interior,offset,globalDims)

      CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
                           this % boundary,bOffset,bGlobalDims)


  END SUBROUTINE WriteHDF5_MPI_Scalar1D

  SUBROUTINE WriteHDF5_Scalar1D(this,fileId,group)
    IMPLICIT NONE
    CLASS(Scalar1D), INTENT(in) :: this
    INTEGER(HID_T), INTENT(in) :: fileId
    CHARACTER(*), INTENT(in) :: group
    ! Local
    INTEGER :: ivar

      CALL CreateGroup_HDF5(fileId,TRIM(group))

      DO ivar = 1, this % nVar
        CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
      ENDDO 

      CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
                           this % interior)

      CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
                           this % boundary)

  END SUBROUTINE WriteHDF5_Scalar1D

! -- Scalar2D -- !

  SUBROUTINE Init_Scalar2D(SELFStorage,interp,nVar,nElem)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(out) :: SELFStorage
    TYPE(Lagrange),INTENT(in),TARGET :: interp
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem

    SELFStorage % interp => interp
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem

    CALL SELFStorage % interior % Alloc(loBound=(/0,0,1,1/), &
                                        upBound=(/interp % N,interp % N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/0,1,1,1/), &
                                        upBound=(/interp % N,nVar,4,nElem/))

    CALL SELFStorage % extBoundary % Alloc(loBound=(/0,1,1,1/), &
                                           upBound=(/interp % N,nVar,4,nElem/))

    CALL SELFStorage % avgBoundary % Alloc(loBound=(/0,1,1,1/), &
                                           upBound=(/interp % N,nVar,4,nElem/))

    CALL SELFStorage % jumpBoundary % Alloc(loBound=(/0,1,1,1/), &
                                           upBound=(/interp % N,nVar,4,nElem/))

    ALLOCATE( SELFStorage % meta(1:nVar) )
    ALLOCATE( SELFStorage % eqn(1:nVar) )
 

  END SUBROUTINE Init_Scalar2D

  SUBROUTINE Free_Scalar2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(inout) :: SELFStorage

    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    SELFStorage % interp => NULL()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()
    CALL SELFStorage % extBoundary % Free()
    CALL SELFStorage % avgBoundary % Free()
    CALL SELFStorage % jumpBoundary % Free()

    DEALLOCATE( SELFStorage % meta )
    DEALLOCATE( SELFStorage % eqn )

  END SUBROUTINE Free_Scalar2D

  SUBROUTINE UpdateHost_Scalar2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
    CALL SELFStorage % extBoundary % UpdateHost()
    CALL SELFStorage % avgBoundary % UpdateHost()
    CALL SELFStorage % jumpBoundary % UpdateHost()

  END SUBROUTINE UpdateHost_Scalar2D

  SUBROUTINE UpdateDevice_Scalar2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
    CALL SELFStorage % extBoundary % UpdateDevice()
    CALL SELFStorage % avgBoundary % UpdateDevice()
    CALL SELFStorage % jumpBoundary % UpdateDevice()

  END SUBROUTINE UpdateDevice_Scalar2D

  SUBROUTINE BoundaryInterp_Scalar2D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarBoundaryInterp_2D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarBoundaryInterp_2D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Scalar2D

  SUBROUTINE GridInterp_Scalar2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(in) :: SELFStorage
    TYPE(Scalar2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarGridInterp_2D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarGridInterp_2D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Scalar2D

  SUBROUTINE Gradient_Scalar2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(in) :: SELFStorage
    TYPE(Vector2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarGradient_2D(SELFStorage % interior % deviceData, &
                                                    SELFout % interior % deviceData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarGradient_2D(SELFStorage % interior % hostData, &
                                                    SELFout % interior % hostData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    END IF

  END SUBROUTINE Gradient_Scalar2D

  FUNCTION AbsMaxInterior_Scalar2D(scalar) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Scalar2D) :: scalar
    REAL(prec) :: absMax(1:scalar % nVar)
    ! Local
    INTEGER :: iEl,iVar,i,j

    absMax = 0.0_prec
    DO iEl = 1,scalar % nElem
      DO iVar = 1,scalar % nVar
        DO j = 0,scalar % interp % N
          DO i = 0,scalar % interp % N
            absMax(iVar) = MAX(ABS(scalar % interior % hostData(i,j,iVar,iEl)),absMax(iVar))
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Scalar2D

  FUNCTION AbsMaxBoundary_Scalar2D(scalar) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Scalar2D) :: scalar
    REAL(prec) :: absMax(1:scalar % nVar,1:4)
    ! Local
    INTEGER :: iEl,iVar,i,iSide

    absMax = 0.0_prec
    DO iEl = 1,scalar % nElem
      DO iSide = 1,4
        DO iVar = 1,scalar % nVar
          DO i = 0,scalar % interp % N
            absMax(iVar,iSide) = MAX(ABS(scalar % boundary % hostData(i,iVar,iSide,iEl)),absMax(iVar,iSide))
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Scalar2D

  SUBROUTINE Equals_Scalar2D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Scalar2D),INTENT(inout) :: SELFOut
    TYPE(Scalar2D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Scalar2D

  SUBROUTINE WriteHDF5_MPI_Scalar2D(this,fileId,group,elemoffset,nglobalelem)
    IMPLICIT NONE
    CLASS(Scalar2D), INTENT(in) :: this
    CHARACTER(*), INTENT(in) :: group
    INTEGER(HID_T), INTENT(in) :: fileId
    INTEGER, INTENT(in) :: elemoffset
    INTEGER, INTENT(in) :: nglobalelem
    ! Local
    INTEGER(HID_T) :: offset(1:4)
    INTEGER(HID_T) :: bOffset(1:4)
    INTEGER(HID_T) :: globalDims(1:4)
    INTEGER(HID_T) :: bGlobalDims(1:4)
    INTEGER :: ivar

      offset(1:4) = (/0,0,0,elemoffset/)
      globalDims(1:4) = (/this % interp % N + 1, &
                          this % interp % N + 1, &
                          this % nVar, &
                          nglobalelem/)

      ! Offsets and dimensions for element boundary data
      bOffset(1:4) = (/0,0,0,elemoffset/)
      bGlobalDims(1:4) = (/this % interp % N + 1, &
                           this % nVar, &
                           4, &
                           nglobalelem/)
  
      CALL CreateGroup_HDF5(fileId,TRIM(group))

      DO ivar = 1, this % nVar
        CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
      ENDDO 

      CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
                           this % interior,offset,globalDims)

      CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
                           this % boundary,bOffset,bGlobalDims)


  END SUBROUTINE WriteHDF5_MPI_Scalar2D

  SUBROUTINE WriteHDF5_Scalar2D(this,fileId,group)
    IMPLICIT NONE
    CLASS(Scalar2D), INTENT(in) :: this
    INTEGER(HID_T), INTENT(in) :: fileId
    CHARACTER(*), INTENT(in) :: group   
    ! Local
    INTEGER :: ivar

      CALL CreateGroup_HDF5(fileId,TRIM(group))

      DO ivar = 1, this % nVar
        CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
      ENDDO 

      CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
                           this % interior)

      CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
                           this % boundary)

  END SUBROUTINE WriteHDF5_Scalar2D

! -- Scalar3D -- !

  SUBROUTINE Init_Scalar3D(SELFStorage,interp,nVar,nElem)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(out) :: SELFStorage
    TYPE(Lagrange),TARGET,INTENT(in) :: interp
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem
    ! Local
    INTEGER :: N

    SELFStorage % interp => interp
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    N = interp % N

    CALL SELFStorage % interior % Alloc(loBound=(/0,0,0,1,1/), &
                                        upBound=(/N,N,N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/0,0,1,1,1/), &
                                        upBound=(/N,N,nVar,6,nElem/))

    CALL SELFStorage % extBoundary % Alloc(loBound=(/0,0,1,1,1/), &
                                           upBound=(/N,N,nVar,6,nElem/))

    CALL SELFStorage % avgBoundary % Alloc(loBound=(/0,0,1,1,1/), &
                                           upBound=(/N,N,nVar,6,nElem/))

    CALL SELFStorage % jumpBoundary % Alloc(loBound=(/0,0,1,1,1/), &
                                           upBound=(/N,N,nVar,6,nElem/))

    ALLOCATE( SELFStorage % meta(1:nVar) )
    ALLOCATE( SELFStorage % eqn(1:nVar) )

  END SUBROUTINE Init_Scalar3D

  SUBROUTINE Free_Scalar3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(inout) :: SELFStorage

    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    SELFStorage % interp => NULL()
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()
    CALL SELFStorage % extBoundary % Free()
    CALL SELFStorage % avgBoundary % Free()
    CALL SELFStorage % jumpBoundary % Free()

    DEALLOCATE( SELFStorage % meta )
    DEALLOCATE( SELFStorage % eqn )

  END SUBROUTINE Free_Scalar3D

  SUBROUTINE UpdateHost_Scalar3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
    CALL SELFStorage % extBoundary % UpdateHost()
    CALL SELFStorage % avgBoundary % UpdateHost()
    CALL SELFStorage % jumpBoundary % UpdateHost()

  END SUBROUTINE UpdateHost_Scalar3D

  SUBROUTINE UpdateDevice_Scalar3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
    CALL SELFStorage % extBoundary % UpdateDevice()
    CALL SELFStorage % avgBoundary % UpdateDevice()
    CALL SELFStorage % jumpBoundary % UpdateDevice()

  END SUBROUTINE UpdateDevice_Scalar3D

  SUBROUTINE BoundaryInterp_Scalar3D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarBoundaryInterp_3D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarBoundaryInterp_3D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Scalar3D

  SUBROUTINE GridInterp_Scalar3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(in) :: SELFStorage
    TYPE(Scalar3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarGridInterp_3D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarGridInterp_3D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Scalar3D

  SUBROUTINE Gradient_Scalar3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(in) :: SELFStorage
    TYPE(Vector3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % ScalarGradient_3D(SELFStorage % interior % deviceData, &
                                                    SELFout % interior % deviceData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % ScalarGradient_3D(SELFStorage % interior % hostData, &
                                                    SELFout % interior % hostData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    END IF

  END SUBROUTINE Gradient_Scalar3D

  SUBROUTINE Equals_Scalar3D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Scalar3D),INTENT(inout) :: SELFOut
    TYPE(Scalar3D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Scalar3D

  FUNCTION AbsMaxInterior_Scalar3D(scalar) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Scalar3D) :: scalar
    REAL(prec) :: absMax(1:scalar % nVar)
    ! Local
    INTEGER :: iEl,iVar,i,j,k

    absMax = 0.0_prec
    DO iEl = 1,scalar % nElem
      DO iVar = 1,scalar % nVar
        DO k = 0,scalar % interp % N
          DO j = 0,scalar % interp % N
            DO i = 0,scalar % interp % N
              absMax(iVar) = MAX(ABS(scalar % interior % hostData(i,j,k,iVar,iEl)),absMax(iVar))
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Scalar3D

  FUNCTION AbsMaxBoundary_Scalar3D(scalar) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Scalar3D) :: scalar
    REAL(prec) :: absMax(1:scalar % nVar,1:6)
    ! Local
    INTEGER :: iEl,iVar,i,j,iSide

    absMax = 0.0_prec
    DO iEl = 1,scalar % nElem
      DO iSide = 1,6
        DO iVar = 1,scalar % nVar
          DO j = 0,scalar % interp % N
            DO i = 0,scalar % interp % N
              absMax(iVar,iSide) = MAX(ABS(scalar % boundary % hostData(i,j,iVar,iSide,iEl)),absMax(iVar,iSide))
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Scalar3D

! -- Vector2D -- !

  SUBROUTINE Init_Vector2D(SELFStorage,interp,nVar,nElem)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(out) :: SELFStorage
    TYPE(Lagrange),TARGET,INTENT(in) :: interp
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem
    ! Local
    INTEGER :: N

    SELFStorage % interp => interp
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    N = interp % N


    CALL SELFStorage % interior % Alloc(loBound=(/1,0,0,1,1/), &
                                        upBound=(/2,N,N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/1,0,1,1,1/), &
                                        upBound=(/2,N,nVar,4,nElem/))

    CALL SELFStorage % boundaryNormal % Alloc(loBound=(/0,1,1,1/), &
                                        upBound=(/N,nVar,4,nElem/))

    CALL SELFStorage % extBoundary % Alloc(loBound=(/1,0,1,1,1/), &
                                           upBound=(/2,N,nVar,4,nElem/))

    ALLOCATE( SELFStorage % meta(1:nVar) )
    ALLOCATE( SELFStorage % eqn(1:2*nVar) )

  END SUBROUTINE Init_Vector2D

  SUBROUTINE Free_Vector2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(inout) :: SELFStorage

    SELFStorage % interp => NULL()
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()
    CALL SELFStorage % boundaryNormal % Free()
    CALL SELFStorage % extBoundary % Free()

    DEALLOCATE( SELFStorage % meta )
    DEALLOCATE( SELFStorage % eqn )

  END SUBROUTINE Free_Vector2D

  SUBROUTINE SetEquation_Vector2D(SELFStorage,idir,ivar,eqnChar)
    !! Sets the equation parser for the `idir` direction and `ivar-th` variable
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(inout) :: SELFStorage
    INTEGER,INTENT(in) :: idir,ivar
    CHARACTER(*),INTENT(in) :: eqnChar

    SELFStorage % eqn(idir+2*(ivar-1)) = EquationParser( TRIM(eqnChar), &
                                              (/'x','y','z','t'/) )

  END SUBROUTINE SetEquation_Vector2D

  SUBROUTINE UpdateHost_Vector2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
    CALL SELFStorage % boundaryNormal % UpdateHost()
    CALL SELFStorage % extBoundary % UpdateHost()

  END SUBROUTINE UpdateHost_Vector2D

  SUBROUTINE UpdateDevice_Vector2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
    CALL SELFStorage % boundaryNormal % UpdateDevice()
    CALL SELFStorage % extBoundary % UpdateDevice()

  END SUBROUTINE UpdateDevice_Vector2D

  SUBROUTINE BoundaryInterp_Vector2D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorBoundaryInterp_2D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorBoundaryInterp_2D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Vector2D

  SUBROUTINE GridInterp_Vector2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(in) :: SELFStorage
    TYPE(Vector2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorGridInterp_2D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorGridInterp_2D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Vector2D

  SUBROUTINE Gradient_Vector2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(in) :: SELFStorage
    TYPE(Tensor2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorGradient_2D(SELFStorage % interior % deviceData, &
                                                    SELFout % interior % deviceData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorGradient_2D(SELFStorage % interior % hostData, &
                                                    SELFout % interior % hostData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    END IF

  END SUBROUTINE Gradient_Vector2D

  SUBROUTINE Divergence_Vector2D(SELFStorage,SELFOut,dForm,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(in) :: SELFStorage
    TYPE(Scalar2D),INTENT(inout) :: SELFOut
    INTEGER,INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

    IF (dForm == selfWeakDGForm) THEN

      IF (gpuAccel) THEN
        CALL SELFStorage % interp % VectorDGDivergence_2D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundaryNormal % deviceData, &
                                                          SELFout % interior % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
      ELSE
        CALL SELFStorage % interp % VectorDGDivergence_2D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundaryNormal % hostData, &
                                                          SELFout % interior % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
      END IF

    ELSE IF (dForm == selfStrongForm) THEN

      IF (gpuAccel) THEN
        CALL SELFStorage % interp % VectorDivergence_2D(SELFStorage % interior % deviceData, &
                                                        SELFout % interior % deviceData, &
                                                        SELFStorage % nVar, &
                                                        SELFStorage % nElem)
      ELSE
        CALL SELFStorage % interp % VectorDivergence_2D(SELFStorage % interior % hostData, &
                                                        SELFout % interior % hostData, &
                                                        SELFStorage % nVar, &
                                                        SELFStorage % nElem)
      END IF

    END IF

  END SUBROUTINE Divergence_Vector2D

  SUBROUTINE Curl_Vector2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(in) :: SELFStorage
    TYPE(Scalar2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorCurl_2D(SELFStorage % interior % deviceData, &
                                                SELFout % interior % deviceData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorCurl_2D(SELFStorage % interior % hostData, &
                                                SELFout % interior % hostData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem)
    END IF

  END SUBROUTINE Curl_Vector2D

  SUBROUTINE Equals_Vector2D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Vector2D),INTENT(inout) :: SELFOut
    TYPE(Vector2D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Vector2D

  FUNCTION AbsMaxInterior_Vector2D(vector) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Vector2D) :: vector
    REAL(prec) :: absMax(1:vector % nVar)
    ! Local
    INTEGER :: iEl,iVar,i,j,iDir

    absMax = 0.0_prec
    DO iEl = 1,vector % nElem
      DO iVar = 1,vector % nVar
        DO j = 0,vector % interp % N
          DO i = 0,vector % interp % N
            DO iDir = 1,2
              absMax(iVar) = MAX(ABS(vector % interior % hostData(iDir,i,j,iVar,iEl)),absMax(iVar))
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Vector2D

  FUNCTION AbsMaxBoundary_Vector2D(vector) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Vector2D) :: vector
    REAL(prec) :: absMax(1:vector % nVar,1:4)
    ! Local
    INTEGER :: iEl,iVar,i,iDir,iSide

    absMax = 0.0_prec
    DO iEl = 1,vector % nElem
      DO iSide = 1,4
        DO iVar = 1,vector % nVar
          DO i = 0,vector % interp % N
            DO iDir = 1,2
              absMax(iVar,iSide) = MAX(ABS(vector % boundary % hostData(iDir,i,iVar,iSide,iEl)),absMax(iVar,iSide))
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Vector2D

  SUBROUTINE WriteHDF5_MPI_Vector2D(this,fileId,group,elemoffset,nglobalelem)
    IMPLICIT NONE
    CLASS(Vector2D), INTENT(in) :: this
    CHARACTER(*), INTENT(in) :: group
    INTEGER(HID_T), INTENT(in) :: fileId
    INTEGER, INTENT(in) :: elemoffset
    INTEGER, INTENT(in) :: nglobalelem
    ! Local
    INTEGER(HID_T) :: offset(1:5)
    INTEGER(HID_T) :: bOffset(1:5)
    INTEGER(HID_T) :: globalDims(1:5)
    INTEGER(HID_T) :: bGlobalDims(1:5)
    INTEGER :: ivar

      offset(1:5) = (/0,0,0,0,elemoffset/)
      globalDims(1:5) = (/2, &
                          this % interp % N + 1, &
                          this % interp % N + 1, &
                          this % nVar, &
                          nglobalelem/)

      ! Offsets and dimensions for element boundary data
      bOffset(1:5) = (/0,0,0,0,elemoffset/)
      bGlobalDims(1:5) = (/2, &
                           this % interp % N + 1, &
                           this % nVar, &
                           4, &
                           nglobalelem/)
  
      CALL CreateGroup_HDF5(fileId,TRIM(group))

      DO ivar = 1, this % nVar
        CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
      ENDDO 

      CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
                           this % interior,offset,globalDims)

      CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
                           this % boundary,bOffset,bGlobalDims)


  END SUBROUTINE WriteHDF5_MPI_Vector2D

  SUBROUTINE WriteHDF5_Vector2D(this,fileId,group)
    IMPLICIT NONE
    CLASS(Vector2D), INTENT(in) :: this
    INTEGER(HID_T), INTENT(in) :: fileId
    CHARACTER(*), INTENT(in) :: group
    ! Local
    INTEGER :: ivar

      CALL CreateGroup_HDF5(fileId,TRIM(group))

      DO ivar = 1, this % nVar
        CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
      ENDDO 

      CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
                           this % interior)

      CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
                           this % boundary)

  END SUBROUTINE WriteHDF5_Vector2D

! -- Vector3D -- !

  SUBROUTINE Init_Vector3D(SELFStorage,interp,nVar,nElem)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(out) :: SELFStorage
    TYPE(Lagrange),TARGET,INTENT(in) :: interp
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem
    ! Local
    INTEGER :: N

    SELFStorage % interp => interp
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    N = interp % N

    CALL SELFStorage % interior % Alloc(loBound=(/1,0,0,0,1,1/), &
                                        upBound=(/3,N,N,N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/1,0,0,1,1,1/), &
                                        upBound=(/3,N,N,nVar,6,nElem/))

    CALL SELFStorage % boundaryNormal % Alloc(loBound=(/0,0,1,1,1/), &
                                        upBound=(/N,N,nVar,6,nElem/))

    CALL SELFStorage % extBoundary % Alloc(loBound=(/1,0,0,1,1,1/), &
                                           upBound=(/3,N,N,nVar,6,nElem/))

    ALLOCATE( SELFStorage % meta(1:nVar) )
    ALLOCATE( SELFStorage % eqn(1:3*nVar) )

  END SUBROUTINE Init_Vector3D

  SUBROUTINE Free_Vector3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(inout) :: SELFStorage

    SELFStorage % interp => NULL()
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()
    CALL SELFStorage % boundaryNormal % Free()
    CALL SELFStorage % extBoundary % Free()

    DEALLOCATE( SELFStorage % meta )
    DEALLOCATE( SELFStorage % eqn )

  END SUBROUTINE Free_Vector3D

  SUBROUTINE SetEquation_Vector3D(SELFStorage,idir,ivar,eqnChar)
    !! Sets the equation parser for the `idir` direction and `ivar-th` variable
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(inout) :: SELFStorage
    INTEGER,INTENT(in) :: idir,ivar
    CHARACTER(*),INTENT(in) :: eqnChar

    SELFStorage % eqn(idir+3*(ivar-1)) = EquationParser( TRIM(eqnChar), &
                                              (/'x','y','z','t'/) )

  END SUBROUTINE SetEquation_Vector3D

  SUBROUTINE UpdateHost_Vector3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
    CALL SELFStorage % boundaryNormal % UpdateHost()
    CALL SELFStorage % extBoundary % UpdateHost()

  END SUBROUTINE UpdateHost_Vector3D

  SUBROUTINE UpdateDevice_Vector3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
    CALL SELFStorage % boundaryNormal % UpdateDevice()
    CALL SELFStorage % extBoundary % UpdateDevice()

  END SUBROUTINE UpdateDevice_Vector3D

  SUBROUTINE BoundaryInterp_Vector3D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorBoundaryInterp_3D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorBoundaryInterp_3D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Vector3D

  SUBROUTINE GridInterp_Vector3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(in) :: SELFStorage
    TYPE(Vector3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorGridInterp_3D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorGridInterp_3D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Vector3D

  SUBROUTINE Gradient_Vector3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(in) :: SELFStorage
    TYPE(Tensor3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorGradient_3D(SELFStorage % interior % deviceData, &
                                                    SELFout % interior % deviceData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorGradient_3D(SELFStorage % interior % hostData, &
                                                    SELFout % interior % hostData, &
                                                    SELFStorage % nVar, &
                                                    SELFStorage % nElem)
    END IF

  END SUBROUTINE Gradient_Vector3D

  SUBROUTINE Divergence_Vector3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(in) :: SELFStorage
    TYPE(Scalar3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorDivergence_3D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorDivergence_3D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE Divergence_Vector3D

  SUBROUTINE Curl_Vector3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(in) :: SELFStorage
    TYPE(Vector3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % VectorCurl_3D(SELFStorage % interior % deviceData, &
                                                SELFout % interior % deviceData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % VectorCurl_3D(SELFStorage % interior % hostData, &
                                                SELFout % interior % hostData, &
                                                SELFStorage % nVar, &
                                                SELFStorage % nElem)
    END IF

  END SUBROUTINE Curl_Vector3D

  FUNCTION AbsMaxInterior_Vector3D(vector) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Vector3D) :: vector
    REAL(prec) :: absMax(1:vector % nVar)
    ! Local
    INTEGER :: iEl,iVar,i,j,k,iDir

    absMax = 0.0_prec
    DO iEl = 1,vector % nElem
      DO iVar = 1,vector % nVar
        DO k = 0,vector % interp % N
          DO j = 0,vector % interp % N
            DO i = 0,vector % interp % N
              DO iDir = 1,3
                absMax(iVar) = MAX(ABS(vector % interior % hostData(iDir,i,j,k,iVar,iEl)),absMax(iVar))
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Vector3D

  FUNCTION AbsMaxBoundary_Vector3D(vector) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Vector3D) :: vector
    REAL(prec) :: absMax(1:vector % nVar,1:6)
    ! Local
    INTEGER :: iEl,iVar,i,j,iSide,iDir

    absMax = 0.0_prec
    DO iEl = 1,vector % nElem
      DO iSide = 1,6
        DO iVar = 1,vector % nVar
          DO j = 0,vector % interp % N
            DO i = 0,vector % interp % N
              DO iDir = 1,3
                absMax(iVar,iSide) = MAX(ABS(vector % boundary % hostData(iDir,i,j,iVar,iSide,iEl)),absMax(iVar,iSide))
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Vector3D

  SUBROUTINE Equals_Vector3D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Vector3D),INTENT(inout) :: SELFOut
    TYPE(Vector3D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Vector3D

! -- P2Vector2D -- !

  SUBROUTINE Init_P2Vector2D(SELFStorage,interp,nVar,nElem)
    IMPLICIT NONE
    CLASS(P2Vector2D),INTENT(out) :: SELFStorage
    TYPE(Lagrange),TARGET,INTENT(in) :: interp
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem
    ! Local
    INTEGER :: N

    SELFStorage % interp => interp
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    N = interp % N


    CALL SELFStorage % interior % Alloc(loBound=(/1,0,0,0,1,1/), &
                                        upBound=(/2,N,N,N,nVar,nElem/))

    CALL SELFStorage % physical % Alloc(loBound=(/1,1,0,0,0,1,1/), &
                                        upBound=(/2,2,N,N,N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/1,0,1,1,1/), &
                                        upBound=(/2,N,nVar,4,nElem/))

    CALL SELFStorage % boundaryNormal % Alloc(loBound=(/0,1,1,1/), &
                                        upBound=(/N,nVar,4,nElem/))

    CALL SELFStorage % extBoundary % Alloc(loBound=(/1,0,1,1,1/), &
                                           upBound=(/2,N,nVar,4,nElem/))

    ALLOCATE( SELFStorage % meta(1:nVar) )
    ALLOCATE( SELFStorage % eqn(1:2*nVar) )

  END SUBROUTINE Init_P2Vector2D

  SUBROUTINE Free_P2Vector2D(SELFStorage)
    IMPLICIT NONE
    CLASS(P2Vector2D),INTENT(inout) :: SELFStorage

    SELFStorage % interp => NULL()
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % physical % Free()
    CALL SELFStorage % boundary % Free()
    CALL SELFStorage % boundaryNormal % Free()
    CALL SELFStorage % extBoundary % Free()

    DEALLOCATE( SELFStorage % meta )
    DEALLOCATE( SELFStorage % eqn )

  END SUBROUTINE Free_P2Vector2D

  SUBROUTINE UpdateHost_P2Vector2D(SELFStorage)
    IMPLICIT NONE
    CLASS(P2Vector2D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % physical % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
    CALL SELFStorage % boundaryNormal % UpdateHost()
    CALL SELFStorage % extBoundary % UpdateHost()

  END SUBROUTINE UpdateHost_P2Vector2D

  SUBROUTINE UpdateDevice_P2Vector2D(SELFStorage)
    IMPLICIT NONE
    CLASS(P2Vector2D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % physical % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
    CALL SELFStorage % boundaryNormal % UpdateDevice()
    CALL SELFStorage % extBoundary % UpdateDevice()

  END SUBROUTINE UpdateDevice_P2Vector2D

  SUBROUTINE Divergence_P2Vector2D(SELFStorage,SELFOut,dForm,gpuAccel)
    IMPLICIT NONE
    CLASS(P2Vector2D),INTENT(in) :: SELFStorage
    TYPE(Scalar2D),INTENT(inout) :: SELFOut
    INTEGER,INTENT(in) :: dForm
    LOGICAL,INTENT(in) :: gpuAccel

    IF (dForm == selfWeakDGForm) THEN

      IF (gpuAccel) THEN
        CALL SELFStorage % interp % P2VectorDGDivergence_2D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundaryNormal % deviceData, &
                                                          SELFout % interior % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
      ELSE
        CALL SELFStorage % interp % P2VectorDGDivergence_2D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundaryNormal % hostData, &
                                                          SELFout % interior % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
      END IF

    ELSE IF (dForm == selfStrongForm) THEN

      IF (gpuAccel) THEN
        CALL SELFStorage % interp % P2VectorDivergence_2D(SELFStorage % interior % deviceData, &
                                                        SELFout % interior % deviceData, &
                                                        SELFStorage % nVar, &
                                                        SELFStorage % nElem)
      ELSE
        CALL SELFStorage % interp % P2VectorDivergence_2D(SELFStorage % interior % hostData, &
                                                        SELFout % interior % hostData, &
                                                        SELFStorage % nVar, &
                                                        SELFStorage % nElem)
      END IF

    END IF

  END SUBROUTINE Divergence_P2Vector2D

! -- Tensor2D -- !

  SUBROUTINE Init_Tensor2D(SELFStorage,interp,nVar,nElem)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(out) :: SELFStorage
    TYPE(Lagrange),TARGET,INTENT(in) :: interp
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem
    ! Local
    INTEGER :: N

    SELFStorage % interp => interp
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    N = interp % N

    CALL SELFStorage % interior % Alloc(loBound=(/1,1,0,0,1,1/), &
                                        upBound=(/2,2,N,N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/1,1,0,1,1,1/), &
                                        upBound=(/2,2,N,nVar,4,nElem/))

    CALL SELFStorage % extBoundary % Alloc(loBound=(/1,1,0,1,1,1/), &
                                           upBound=(/2,2,N,nVar,4,nElem/))

    ALLOCATE( SELFStorage % meta(1:nVar) )
    ALLOCATE( SELFStorage % eqn(1:4*nVar) )

  END SUBROUTINE Init_Tensor2D

  SUBROUTINE Free_Tensor2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(inout) :: SELFStorage

    SELFStorage % interp => NULL()
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()
    CALL SELFStorage % extBoundary % Free()

    DEALLOCATE( SELFStorage % meta )
    DEALLOCATE( SELFStorage % eqn )

  END SUBROUTINE Free_Tensor2D

  SUBROUTINE SetEquation_Tensor2D(SELFStorage,row,col,ivar,eqnChar)
    !! Sets the equation parser for row, col  of the ivar-th tensor
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(inout) :: SELFStorage
    INTEGER,INTENT(in) :: row,col,ivar
    CHARACTER(*),INTENT(in) :: eqnChar
    ! Local
    INTEGER :: ind

    ind = row+2*(col-1+2*(ivar-1))
    SELFStorage % eqn(ind) = EquationParser( TRIM(eqnChar), &
                                              (/'x','y','z','t'/) )

  END SUBROUTINE SetEquation_Tensor2D

  SUBROUTINE UpdateHost_Tensor2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
    CALL SELFStorage % extBoundary % UpdateHost()

  END SUBROUTINE UpdateHost_Tensor2D

  SUBROUTINE UpdateDevice_Tensor2D(SELFStorage)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
    CALL SELFStorage % extBoundary % UpdateDevice()

  END SUBROUTINE UpdateDevice_Tensor2D

  SUBROUTINE BoundaryInterp_Tensor2D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % TensorBoundaryInterp_2D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % TensorBoundaryInterp_2D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Tensor2D

  SUBROUTINE GridInterp_Tensor2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(in) :: SELFStorage
    TYPE(Tensor2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % TensorGridInterp_2D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % TensorGridInterp_2D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Tensor2D

  SUBROUTINE Divergence_Tensor2D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(in) :: SELFStorage
    TYPE(Vector2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % TensorDivergence_2D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % TensorDivergence_2D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE Divergence_Tensor2D

  SUBROUTINE Determinant_Tensor2D(SELFStorage,SELFout,gpuAccel)
#undef __FUNC__
#define __FUNC__ "Determinant_Tensor2D"
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(in) :: SELFStorage
    TYPE(Scalar2D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iEl,iVar,i,j

    IF (gpuAccel) THEN

      CALL Determinant_Tensor2D_gpu_wrapper(SELFStorage % interior % deviceData, &
                                            SELFOut % interior % deviceData, &
                                            SELFStorage % interp % N, &
                                            SELFStorage % nVar, &
                                            SELFStorage % nElem)

    ELSE

      DO iEl = 1,SELFStorage % nElem
        DO iVar = 1,SELFStorage % nVar
          DO j = 0,SELFStorage % interp % N
            DO i = 0,SELFStorage % interp % N

              SELFOut % interior % hostData(i,j,iVar,iEl) = SELFStorage % interior % hostData(1,1,i,j,iVar,iEl)* &
                                                            SELFStorage % interior % hostData(2,2,i,j,iVar,iEl) - &
                                                            SELFStorage % interior % hostData(1,2,i,j,iVar,iEl)* &
                                                            SELFStorage % interior % hostData(2,1,i,j,iVar,iEl)

            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE Determinant_Tensor2D

  FUNCTION AbsMaxInterior_Tensor2D(tensor) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Tensor2D) :: tensor
    REAL(prec) :: absMax(1:tensor % nVar)
    ! Local
    INTEGER :: iEl,iVar,i,j,row,col

    absMax = 0.0_prec
    DO iEl = 1,tensor % nElem
      DO iVar = 1,tensor % nVar
        DO j = 0,tensor % interp % N
          DO i = 0,tensor % interp % N
            DO col = 1,2
              DO row = 1,2
                absMax(iVar) = MAX(ABS(tensor % interior % hostData(row,col,i,j,iVar,iEl)),absMax(iVar))
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Tensor2D

  FUNCTION AbsMaxBoundary_Tensor2D(tensor) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Tensor2D) :: tensor
    REAL(prec) :: absMax(1:tensor % nVar,1:4)
    ! Local
    INTEGER :: iEl,iVar,i,iSide,row,col

    absMax = 0.0_prec
    DO iEl = 1,tensor % nElem
      DO iSide = 1,4
        DO iVar = 1,tensor % nVar
          DO i = 0,tensor % interp % N
            DO col = 1,2
              DO row = 1,2
                absMax(iVar,iSide) = MAX(ABS(tensor % boundary % hostData(row,col,i,iVar,iSide,iEl)),absMax(iVar,iSide))
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Tensor2D

  SUBROUTINE Equals_Tensor2D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Tensor2D),INTENT(inout) :: SELFOut
    TYPE(Tensor2D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Tensor2D

! -- Tensor3D -- !

  SUBROUTINE Init_Tensor3D(SELFStorage,interp,nVar,nElem)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(out) :: SELFStorage
    TYPE(Lagrange),TARGET,INTENT(in) :: interp
    INTEGER,INTENT(in) :: nVar
    INTEGER,INTENT(in) :: nElem
    ! Local
    INTEGER :: N

    SELFStorage % interp => interp
    SELFStorage % nVar = nVar
    SELFStorage % nElem = nElem
    N = interp % N

    CALL SELFStorage % interior % Alloc(loBound=(/1,1,0,0,0,1,1/), &
                                        upBound=(/3,3,N,N,N,nVar,nElem/))

    CALL SELFStorage % boundary % Alloc(loBound=(/1,1,0,0,1,1,1/), &
                                        upBound=(/3,3,N,N,nVar,6,nElem/))

    CALL SELFStorage % extBoundary % Alloc(loBound=(/1,1,0,0,1,1,1/), &
                                           upBound=(/3,3,N,N,nVar,6,nElem/))

    ALLOCATE( SELFStorage % meta(1:nVar) )
    ALLOCATE( SELFStorage % eqn(1:9*nVar) )

  END SUBROUTINE Init_Tensor3D

  SUBROUTINE Free_Tensor3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(inout) :: SELFStorage

    SELFStorage % interp => NULL()
    SELFStorage % nVar = 0
    SELFStorage % nElem = 0
    CALL SELFStorage % interior % Free()
    CALL SELFStorage % boundary % Free()
    CALL SELFStorage % extBoundary % Free()
    
    DEALLOCATE( SELFStorage % meta )
    DEALLOCATE( SELFStorage % eqn )

  END SUBROUTINE Free_Tensor3D

  SUBROUTINE SetEquation_Tensor3D(SELFStorage,row,col,ivar,eqnChar)
    !! Sets the equation parser for row, col  of the ivar-th tensor
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(inout) :: SELFStorage
    INTEGER,INTENT(in) :: row,col,ivar
    CHARACTER(*),INTENT(in) :: eqnChar
    ! Local
    INTEGER :: ind

    ind = row+3*(col-1+3*(ivar-1))
    SELFStorage % eqn(ind) = EquationParser( TRIM(eqnChar), &
                                              (/'x','y','z','t'/) )

  END SUBROUTINE SetEquation_Tensor3D

  SUBROUTINE UpdateHost_Tensor3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateHost()
    CALL SELFStorage % boundary % UpdateHost()
    CALL SELFStorage % extBoundary % UpdateHost()

  END SUBROUTINE UpdateHost_Tensor3D

  SUBROUTINE UpdateDevice_Tensor3D(SELFStorage)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(inout) :: SELFStorage

    CALL SELFStorage % interior % UpdateDevice()
    CALL SELFStorage % boundary % UpdateDevice()
    CALL SELFStorage % extBoundary % UpdateHost()

  END SUBROUTINE UpdateDevice_Tensor3D

  SUBROUTINE BoundaryInterp_Tensor3D(SELFStorage,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(inout) :: SELFStorage
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % TensorBoundaryInterp_3D(SELFStorage % interior % deviceData, &
                                                          SELFStorage % boundary % deviceData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % TensorBoundaryInterp_3D(SELFStorage % interior % hostData, &
                                                          SELFStorage % boundary % hostData, &
                                                          SELFStorage % nVar, &
                                                          SELFStorage % nElem)
    END IF

  END SUBROUTINE BoundaryInterp_Tensor3D

  SUBROUTINE GridInterp_Tensor3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(in) :: SELFStorage
    TYPE(Tensor3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % TensorGridInterp_3D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % TensorGridInterp_3D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE GridInterp_Tensor3D

  SUBROUTINE Divergence_Tensor3D(SELFStorage,SELFOut,gpuAccel)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(in) :: SELFStorage
    TYPE(Vector3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel

    IF (gpuAccel) THEN
      CALL SELFStorage % interp % TensorDivergence_3D(SELFStorage % interior % deviceData, &
                                                      SELFout % interior % deviceData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    ELSE
      CALL SELFStorage % interp % TensorDivergence_3D(SELFStorage % interior % hostData, &
                                                      SELFout % interior % hostData, &
                                                      SELFStorage % nVar, &
                                                      SELFStorage % nElem)
    END IF

  END SUBROUTINE Divergence_Tensor3D

  SUBROUTINE Determinant_Tensor3D(SELFStorage,SELFOut,gpuAccel)
#undef __FUNC__
#define __FUNC__ "Determinant_Tensor3D"
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(in) :: SELFStorage
    TYPE(Scalar3D),INTENT(inout) :: SELFOut
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: iEl,iVar,i,j,k

    IF (gpuAccel) THEN

      CALL Determinant_Tensor3D_gpu_wrapper(SELFStorage % interior % deviceData, &
                                            SELFOut % interior % deviceData, &
                                            SELFStorage % interp % N, &
                                            SELFStorage % nVar, &
                                            SELFStorage % nElem)

    ELSE

      DO iEl = 1,SELFStorage % nElem
        DO iVar = 1,SELFStorage % nVar
          DO k = 0,SELFStorage % interp % N
            DO j = 0,SELFStorage % interp % N
              DO i = 0,SELFStorage % interp % N

                SELFOut % interior % hostData(i,j,k,iVar,iEl) = &
                  SELFStorage % interior % hostData(1,1,i,j,k,iVar,iEl)* &
                  (SELFStorage % interior % hostData(2,2,i,j,k,iVar,iEl)* &
                   SELFStorage % interior % hostData(3,3,i,j,k,iVar,iEl) - &
                   SELFStorage % interior % hostData(2,3,i,j,k,iVar,iEl)* &
                   SELFStorage % interior % hostData(3,2,i,j,k,iVar,iEl)) - &
                  SELFStorage % interior % hostData(2,1,i,j,k,iVar,iEl)* &
                  (SELFStorage % interior % hostData(1,2,i,j,k,iVar,iEl)* &
                   SELFStorage % interior % hostData(3,3,i,j,k,iVar,iEl) - &
                   SELFStorage % interior % hostData(1,3,i,j,k,iVar,iEl)* &
                   SELFStorage % interior % hostData(3,2,i,j,k,iVar,iEl)) + &
                  SELFStorage % interior % hostData(3,1,i,j,k,iVar,iEl)* &
                  (SELFStorage % interior % hostData(1,2,i,j,k,iVar,iEl)* &
                   SELFStorage % interior % hostData(2,3,i,j,k,iVar,iEl) - &
                   SELFStorage % interior % hostData(1,3,i,j,k,iVar,iEl)* &
                   SELFStorage % interior % hostData(2,2,i,j,k,iVar,iEl))

              END DO
            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE Determinant_Tensor3D

  FUNCTION AbsMaxInterior_Tensor3D(tensor) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Tensor3D) :: tensor
    REAL(prec) :: absMax(1:tensor % nVar)
    ! Local
    INTEGER :: iEl,iVar,i,j,k,row,col

    absMax = 0.0_prec
    DO iEl = 1,tensor % nElem
      DO iVar = 1,tensor % nVar
        DO k = 0,tensor % interp % N
          DO j = 0,tensor % interp % N
            DO i = 0,tensor % interp % N
              DO col = 1,3
                DO row = 1,3
                  absMax(iVar) = MAX(ABS(tensor % interior % hostData(row,col,i,j,k,iVar,iEl)),absMax(iVar))
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxInterior_Tensor3D

  FUNCTION AbsMaxBoundary_Tensor3D(tensor) RESULT(absMax)
    IMPLICIT NONE
    CLASS(Tensor3D) :: tensor
    REAL(prec) :: absMax(1:tensor % nVar,1:6)
    ! Local
    INTEGER :: iEl,iVar,i,j,iSide,row,col

    absMax = 0.0_prec
    DO iEl = 1,tensor % nElem
      DO iSide = 1,6
        DO iVar = 1,tensor % nVar
          DO j = 0,tensor % interp % N
            DO i = 0,tensor % interp % N
              DO col = 1,3
                DO row = 1,3
              absMax(iVar,iSide) = MAX(ABS(tensor % boundary % hostData(row,col,i,j,iVar,iSide,iEl)),absMax(iVar,iSide))
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END FUNCTION AbsMaxBoundary_Tensor3D

  SUBROUTINE Equals_Tensor3D(SELFOut,SELFin)
    IMPLICIT NONE
    CLASS(Tensor3D),INTENT(inout) :: SELFOut
    TYPE(Tensor3D),INTENT(in) :: SELFin

    SELFOut % interior % hostData = SELFin % interior % hostData
    SELFOut % boundary % hostData = SELFin % boundary % hostData

  END SUBROUTINE Equals_Tensor3D

END MODULE SELF_Data
