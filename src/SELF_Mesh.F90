MODULE SELF_Mesh

  USE SELF_Constants
  USE SELF_Lagrange
  USE SELF_Data
  USE SELF_SupportRoutines

  USE hipfort
  USE ISO_C_BINDING

  IMPLICIT NONE

#include "SELF_Macros.h"
! ========================================================================= !
! Node, Edge, Face, Element and Connectivity Standard
! ========================================================================= !
!
! To define the element corner nodes, the side order and side connectivity,
! we follow the standard from CGNS SIDS (CFD General Notation System,
! Standard Interface Data Structures, http: //cgns.sourceforge.net/ ).
!
! Computational coordinate directions are defined as follows
!
! xi1 direction points from "West" (xi1=-1) to "East" (xi1=1)
! xi2 direction points from "South" (xi2=-1) to "North" (xi2=1)
! xi3 direction points from "Bottom" (xi3=-1) to "Top" (xi3=1)
!
! 2-D Hexahedreal Element sides are defined as
!
! Side 1 = South  (xi2 = -1) = [CN1, CN2]
! Side 2 = East   (xi1 = 1) = [CN2, CN3]
! Side 3 = North  (xi2 = 1) = [CN4, CN3]
! Side 4 = West   (xi1 = -1) = [CN1, CN4]
!
! 3-D Hexahedreal Element sides are defined as
!
! Side 1 = Bottom (xi3 = -1) = [CN1, CN4, CN3, CN2]
! Side 2 = South  (xi2 = -1) = [CN1, CN2, CN6, CN5]
! Side 3 = East   (xi1 = 1) = [CN2, CN3, CN7, CN6]
! Side 4 = North  (xi2 = 1) = [CN3, CN4, CN8, CN7]
! Side 5 = West   (xi1 = -1) = [CN1, CN5, CN8, CN4]
! Side 6 = Top    (xi3 = 1) = [CN5, CN6, CN7, CN8]
!
! In 2-D, corner nodes are order counter-clockwise (looking in the -xi3 direction).
!
! CornerNode 1 = South-West = (-1,-1)
! CornerNode 2 = South-East = ( 1,-1)
! CornerNode 3 = North-East = ( 1, 1)
! CornerNode 4 = North-West = (-1, 1)
!
! In 3-D, corner nodes are order counter-clockwise (looking in the -xi3 direction) from
! bottom to top.
!
! CornerNode 1 = Bottom-South-West = (-1,-1,-1)
! CornerNode 2 = Bottom-South-East = ( 1,-1,-1)
! CornerNode 3 = Bottom-North-East = ( 1, 1,-1)
! CornerNode 4 = Bottom-North-West = (-1, 1,-1)
! CornerNode 5 = Top-South-West = (-1,-1, 1)
! CornerNode 6 = Top-South-East = ( 1,-1, 1)
! CornerNode 7 = Top-North-East = ( 1, 1, 1)
! CornerNode 8 = Top-North-West = (-1, 1, 1)
!
! ========================================================================= !

  INTEGER,PARAMETER :: selfMinNodalValence2D = 4
  INTEGER,PARAMETER :: selfMinNodalValence3D = 8
  INTEGER,PARAMETER :: selfMaxNodalValence2D = 6
  INTEGER,PARAMETER :: selfMaxNodalValence3D = 10
  ! Side Ordering
  INTEGER,PARAMETER :: selfSide2D_South = 1
  INTEGER,PARAMETER :: selfSide2D_East = 2
  INTEGER,PARAMETER :: selfSide2D_North = 3
  INTEGER,PARAMETER :: selfSide2D_West = 4
  INTEGER,PARAMETER :: selfSide3D_Bottom = 1
  INTEGER,PARAMETER :: selfSide3D_South = 2
  INTEGER,PARAMETER :: selfSide3D_East = 3
  INTEGER,PARAMETER :: selfSide3D_North = 4
  INTEGER,PARAMETER :: selfSide3D_West = 5
  INTEGER,PARAMETER :: selfSide3D_Top = 6

  TYPE,PUBLIC :: Geometry1D
    INTEGER :: cqType ! Control Quadrature Type
    INTEGER :: tqType ! Target Quadrature Type
    INTEGER :: nElem
    TYPE(Scalar1D) :: x ! Physical Positions
    TYPE(Scalar1D) :: dxds ! Conversion from computational to physical space

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Geometry1D
    PROCEDURE,PUBLIC :: Free => Free_Geometry1D
#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Geometry1D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Geometry1D
#endif
    PROCEDURE,PUBLIC :: CalculateMetricTerms => CalculateMetricTerms_Geometry1D

  END TYPE Geometry1D

  TYPE,PUBLIC :: Geometry2D
    INTEGER :: cqType ! Control Quadrature Type
    INTEGER :: tqType ! Target Quadrature Type
    INTEGER :: nElem
    TYPE(Vector2D) :: x ! Physical positions
    TYPE(Tensor2D) :: dxds ! Covariant basis vectors
    TYPE(Tensor2D) :: dsdx ! Contavariant basis vectors
    TYPE(Scalar2D) :: J ! Jacobian of the transformation

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Geometry2D
    PROCEDURE,PUBLIC :: Free => Free_Geometry2D
#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Geometry2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Geometry2D
#endif
    PROCEDURE,PUBLIC :: CalculateMetricTerms => CalculateMetricTerms_Geometry2D
    PROCEDURE,PRIVATE :: CalculateContravariantBasis => CalculateContravariantBasis_Geometry2D

  END TYPE Geometry2D

  TYPE,PUBLIC :: Geometry3D
    INTEGER :: cqType ! Control Quadrature Type
    INTEGER :: tqType ! Target Quadrature Type
    INTEGER :: nElem
    TYPE(Vector3D) :: x ! Physical positions
    TYPE(Tensor3D) :: dxds ! Covariant basis vectors
    TYPE(Tensor3D) :: dsdx ! Contavariant basis vectors
    TYPE(Scalar3D) :: J ! Jacobian of the transformation

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Geometry3D
    PROCEDURE,PUBLIC :: Free => Free_Geometry3D
#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Geometry3D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Geometry3D
#endif
    PROCEDURE,PUBLIC :: CalculateMetricTerms => CalculateMetricTerms_Geometry3D
    PROCEDURE,PRIVATE :: CalculateContravariantBasis => CalculateContravariantBasis_Geometry3D

  END TYPE Geometry3D

  TYPE,PUBLIC :: Mesh1D
    INTEGER :: nGeo
    INTEGER :: nElem
    INTEGER :: nNodes
    INTEGER :: nUniqueNodes
    INTEGER :: nBCs
    TYPE(Geometry1D) :: geometry
    TYPE(hfInt32_r2) :: elemInfo
    TYPE(hfReal_r1)  :: nodeCoords
    TYPE(hfInt32_r2) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS
    PROCEDURE,PUBLIC :: Init => Init_Mesh1D
    PROCEDURE,PUBLIC :: Free => Free_Mesh1D
#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Mesh1D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh1D
#endif
    PROCEDURE,PUBLIC :: UniformBlockMesh => UniformBlockMesh_Mesh1D

  END TYPE Mesh1D

  ! Mesh format is set up as the HOPr format
  ! See https://hopr-project.org/externals/MeshFormat.pdf

  TYPE,PUBLIC :: Mesh2D
    INTEGER :: nGeo
    INTEGER :: nElem
    INTEGER :: nSides
    INTEGER :: nNodes
    INTEGER :: nUniqueSides
    INTEGER :: nUniqueNodes
    INTEGER :: nBCs
    TYPE(Geometry2D) :: geometry
    TYPE(hfInt32_r2) :: elemInfo
    TYPE(hfInt32_r2) :: sideInfo
    TYPE(hfReal_r2)  :: nodeCoords
    TYPE(hfInt32_r1) :: globalNodeIDs
    TYPE(hfInt32_r2) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS
    PROCEDURE,PUBLIC :: Init => Init_Mesh2D
    PROCEDURE,PUBLIC :: Free => Free_Mesh2D
#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Mesh2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh2D
#endif

  END TYPE Mesh2D

  TYPE,PUBLIC :: Mesh3D
    INTEGER :: nGeo
    INTEGER :: nElem
    INTEGER :: nSides
    INTEGER :: nNodes
    INTEGER :: nUniqueSides
    INTEGER :: nUniqueNodes
    INTEGER :: nBCs
    TYPE(Geometry3D) :: geometry
    TYPE(hfInt32_r2) :: elemInfo
    TYPE(hfInt32_r2) :: sideInfo
    TYPE(hfReal_r2)  :: nodeCoords
    TYPE(hfInt32_r1) :: globalNodeIDs
    TYPE(hfInt32_r2) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Mesh3D
    PROCEDURE,PUBLIC :: Free => Free_Mesh3D
#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Mesh3D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh3D
#endif

!      PROCEDURE, PUBLIC :: LoadHOPRMesh => LoadHOPRMesh_3D
!
!      PROCEDURE, PUBLIC :: Read_CGNSMesh
!      PROCEDURE, PUBLIC :: Read_UCDMesh
!      PROCEDURE, PUBLIC :: Read_TrellisUCDMesh
!
!      PROCEDURE, PUBLIC :: GenerateConnectivity => GenerateConnectivity_Mesh3D

  END TYPE Mesh3D

CONTAINS

  SUBROUTINE Init_Mesh1D(myMesh,cqType,tqType,nControlPoints,nTargetPoints, &
                         nElem,nNodes,nUniqueNodes,nBCs)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(out) :: myMesh
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: nControlPoints
    INTEGER,INTENT(in) :: nTargetPoints
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nNodes
    INTEGER,INTENT(in) :: nUniqueNodes
    INTEGER,INTENT(in) :: nBCs

    myMesh % nElem = nElem
    myMesh % nNodes = nNodes
    myMesh % nUniqueNodes = nUniqueNodes
    myMesh % nBCs = nBCs

    CALL myMesh % geometry % Init(cqType,tqType,nControlPoints,nTargetPoints,nElem)

    CALL myMesh % elemInfo % Alloc(loBound=(/1,1/), &
                                   upBound=(/4,nElem/))

    CALL myMesh % nodeCoords % Alloc(loBound=1, &
                                     upBound=nNodes)

    CALL myMesh % BCType % Alloc(loBound=(/1,1/), &
                                 upBound=(/4,nBCs/))

    ALLOCATE (myMesh % BCNames(1:nBCs))

  END SUBROUTINE Init_Mesh1D

  SUBROUTINE Free_Mesh1D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(inout) :: myMesh

    myMesh % nElem = 0
    myMesh % nNodes = 0
    myMesh % nUniqueNodes = 0
    myMesh % nBCs = 0

    CALL myMesh % geometry % Free()
    CALL myMesh % elemInfo % Free()
    CALL myMesh % nodeCoords % Free()
    CALL myMesh % BCType % Free()

    DEALLOCATE (myMesh % BCNames)

  END SUBROUTINE Free_Mesh1D
#ifdef GPU
  SUBROUTINE UpdateHost_Mesh1D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(inout) :: myMesh

    CALL myMesh % geometry % UpdateHost()
    CALL myMesh % elemInfo % UpdateHost()
    CALL myMesh % nodeCoords % UpdateHost()
    CALL myMesh % BCType % UpdateHost()

  END SUBROUTINE UpdateHost_Mesh1D

  SUBROUTINE UpdateDevice_Mesh1D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(inout) :: myMesh

    CALL myMesh % geometry % UpdateDevice()
    CALL myMesh % elemInfo % UpdateDevice()
    CALL myMesh % nodeCoords % UpdateDevice()
    CALL myMesh % BCType % UpdateDevice()

  END SUBROUTINE UpdateDevice_Mesh1D
#endif

  SUBROUTINE UniformBlockMesh_Mesh1D(myMesh,cqType,tqType,nControlPoints,nTargetPoints,nElem,x0,x1)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(out) :: myMesh
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: nControlPoints
    INTEGER,INTENT(in) :: nTargetPoints
    INTEGER,INTENT(in) :: nElem
    REAL(prec),INTENT(in) :: x0
    REAL(prec),INTENT(in) :: x1
    ! Local
    INTEGER :: iel,nl,nr,i
    REAL(prec) :: xl,xr,s

    CALL myMesh % Init(cqType,tqType,nControlPoints, &
                       nTargetPoints,nElem,nElem + 1,nElem + 1,2)

    ! Set the element boundaries
    myMesh % nodeCoords % hostData = UniformPoints(x0,x1,nElem)

    ! Set the element information
    DO iel = 1,nElem
      myMesh % eleminfo % hostData(1,iel) = 1 ! Element Type
      myMesh % eleminfo % hostData(2,iel) = 1 ! Element Zone
      myMesh % eleminfo % hostData(3,iel) = iel ! Node Index Start
      myMesh % eleminfo % hostData(4,iel) = iel + 1 ! Node Index End
    END DO

    ! Set the element internal mesh locations
    DO iel = 1,nElem
      nl = myMesh % eleminfo % hostData(3,iel)
      nr = myMesh % eleminfo % hostData(4,iel)
      xl = myMesh % nodeCoords % hostData(nl)
      xr = myMesh % nodeCoords % hostData(nr)
      DO i = 0,nControlPoints
        s = myMesh % geometry % x % interp % controlPoints % hostData(i)
        ! Transfinite interpolation (Linear)
        myMesh % geometry % x % interior % hostData(i,1,iel) = 0.5_prec* &
                                                               ((s + 1.0_prec)*xr &
                                                                - (s - 1.0_prec)*xl)
      END DO
    END DO

    CALL myMesh % geometry % x % BoundaryInterp(gpuAccel=.FALSE.)
    CALL myMesh % geometry % CalculateMetricTerms()

#ifdef GPU
    CALL myMesh % UpdateDevice()
#endif

  END SUBROUTINE UniformBlockMesh_Mesh1D

  SUBROUTINE Init_Geometry1D(myGeom,cqType,tqType,nControlPoints,nTargetPoints,nElem)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(out) :: myGeom
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: nControlPoints
    INTEGER,INTENT(in) :: nTargetPoints
    INTEGER,INTENT(in) :: nElem

    myGeom % cqType = cqType
    myGeom % tqType = tqType
    myGeom % nElem = nElem

    CALL myGeom % x % Init(N=nControlPoints, &
                           quadratureType=cqType, &
                           M=nTargetPoints, &
                           targetNodeType=tqType, &
                           nVar=1, &
                           nElem=nElem)

    CALL myGeom % dxds % Init(N=nControlPoints, &
                              quadratureType=cqType, &
                              M=nTargetPoints, &
                              targetNodeType=tqType, &
                              nVar=1, &
                              nElem=nElem)

  END SUBROUTINE Init_Geometry1D

  SUBROUTINE Free_Geometry1D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(inout) :: myGeom

    CALL myGeom % x % Free()
    CALL myGeom % dxds % Free()

  END SUBROUTINE Free_Geometry1D
#ifdef GPU
  SUBROUTINE UpdateHost_Geometry1D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(inout) :: myGeom

    CALL myGeom % x % UpdateHost()
    CALL myGeom % dxds % UpdateHost()

  END SUBROUTINE UpdateHost_Geometry1D

  SUBROUTINE UpdateDevice_Geometry1D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(inout) :: myGeom

    CALL myGeom % x % UpdateDevice()
    CALL myGeom % dxds % UpdateDevice()

  END SUBROUTINE UpdateDevice_Geometry1D
#endif

  SUBROUTINE CalculateMetricTerms_Geometry1D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(inout) :: myGeom

    CALL myGeom % x % Derivative(myGeom % dxds,gpuAccel=.FALSE.)
    CALL myGeom % dxds % BoundaryInterp(gpuAccel=.FALSE.)

#ifdef GPU
    CALL myGeom % UpdateDevice()
#endif

  END SUBROUTINE CalculateMetricTerms_Geometry1D

  SUBROUTINE Init_Mesh2D(myMesh,cqType,tqType,nControlPoints,nTargetPoints, &
                         nElem,nSides,nNodes,nUniqueSides,nUniqueNodes,nBCs)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(out) :: myMesh
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: nControlPoints
    INTEGER,INTENT(in) :: nTargetPoints
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nSides
    INTEGER,INTENT(in) :: nNodes
    INTEGER,INTENT(in) :: nUniqueSides
    INTEGER,INTENT(in) :: nUniqueNodes
    INTEGER,INTENT(in) :: nBCs

    myMesh % nElem = nElem
    myMesh % nSides = nSides
    myMesh % nNodes = nNodes
    myMesh % nUniqueSides = nUniqueSides
    myMesh % nUniqueNodes = nUniqueNodes
    myMesh % nBCs = nBCs

    CALL myMesh % geometry % Init(cqType,tqType,nControlPoints,nTargetPoints,nElem)

    CALL myMesh % elemInfo % Alloc(loBound=(/1,1/), &
                                   upBound=(/6,nElem/))

    CALL myMesh % sideInfo % Alloc(loBound=(/1,1/), &
                                   upBound=(/5,nSides/))

    CALL myMesh % nodeCoords % Alloc(loBound=(/1,1/), &
                                     upBound=(/3,nNodes/))

    CALL myMesh % globalNodeIDs % Alloc(loBound=1, &
                                        upBound=nNodes)

    CALL myMesh % BCType % Alloc(loBound=(/1,1/), &
                                 upBound=(/4,nBCs/))

    ALLOCATE (myMesh % BCNames(1:nBCs))

  END SUBROUTINE Init_Mesh2D

  SUBROUTINE Free_Mesh2D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh

    myMesh % nElem = 0
    myMesh % nSides = 0
    myMesh % nNodes = 0
    myMesh % nUniqueSides = 0
    myMesh % nUniqueNodes = 0
    myMesh % nBCs = 0

    CALL myMesh % elemInfo % Free()
    CALL myMesh % sideInfo % Free()
    CALL myMesh % nodeCoords % Free()
    CALL myMesh % globalNodeIDs % Free()
    CALL myMesh % BCType % Free()

    DEALLOCATE (myMesh % BCNames)

  END SUBROUTINE Free_Mesh2D
#ifdef GPU
  SUBROUTINE UpdateHost_Mesh2D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh

    CALL myMesh % geometry % UpdateHost()
    CALL myMesh % elemInfo % UpdateHost()
    CALL myMesh % sideInfo % UpdateHost()
    CALL myMesh % nodeCoords % UpdateHost()
    CALL myMesh % globalNodeIDs % UpdateHost()
    CALL myMesh % BCType % UpdateHost()

  END SUBROUTINE UpdateHost_Mesh2D

  SUBROUTINE UpdateDevice_Mesh2D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh

    CALL myMesh % geometry % UpdateDevice()
    CALL myMesh % elemInfo % UpdateDevice()
    CALL myMesh % sideInfo % UpdateDevice()
    CALL myMesh % nodeCoords % UpdateDevice()
    CALL myMesh % globalNodeIDs % UpdateDevice()
    CALL myMesh % BCType % UpdateDevice()

  END SUBROUTINE UpdateDevice_Mesh2D
#endif
  SUBROUTINE Init_Geometry2D(myGeom,cqType,tqType,nControlPoints,nTargetPoints,nElem)
    IMPLICIT NONE
    CLASS(Geometry2D),INTENT(out) :: myGeom
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: nControlPoints
    INTEGER,INTENT(in) :: nTargetPoints
    INTEGER,INTENT(in) :: nElem

    myGeom % cqType = cqType
    myGeom % tqType = tqType
    myGeom % nElem = nElem

    CALL myGeom % x % Init(N=nControlPoints, &
                           quadratureType=cqType, &
                           M=nTargetPoints, &
                           targetNodeType=tqType, &
                           nVar=1, &
                           nElem=nElem)

    CALL myGeom % dxds % Init(N=nControlPoints, &
                              quadratureType=cqType, &
                              M=nTargetPoints, &
                              targetNodeType=tqType, &
                              nVar=1, &
                              nElem=nElem)

    CALL myGeom % dsdx % Init(N=nControlPoints, &
                              quadratureType=cqType, &
                              M=nTargetPoints, &
                              targetNodeType=tqType, &
                              nVar=1, &
                              nElem=nElem)

    CALL myGeom % J % Init(N=nControlPoints, &
                           quadratureType=cqType, &
                           M=nTargetPoints, &
                           targetNodeType=tqType, &
                           nVar=1, &
                           nElem=nElem)

  END SUBROUTINE Init_Geometry2D

  SUBROUTINE Free_Geometry2D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry2D),INTENT(inout) :: myGeom

    CALL myGeom % x % Free()
    CALL myGeom % dxds % Free()
    CALL myGeom % dsdx % Free()
    CALL myGeom % J % Free()

  END SUBROUTINE Free_Geometry2D
#ifdef GPU
  SUBROUTINE UpdateHost_Geometry2D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry2D),INTENT(inout) :: myGeom

    CALL myGeom % x % UpdateHost()
    CALL myGeom % dxds % UpdateHost()
    CALL myGeom % dsdx % UpdateHost()
    CALL myGeom % J % UpdateHost()

  END SUBROUTINE UpdateHost_Geometry2D

  SUBROUTINE UpdateDevice_Geometry2D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry2D),INTENT(inout) :: myGeom

    CALL myGeom % x % UpdateDevice()
    CALL myGeom % dxds % UpdateDevice()
    CALL myGeom % dsdx % UpdateDevice()
    CALL myGeom % J % UpdateDevice()

  END SUBROUTINE UpdateDevice_Geometry2D
#endif
  SUBROUTINE CalculateContravariantBasis_Geometry2D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry2D),INTENT(inout) :: myGeom
    ! Local
    INTEGER :: iEl,i,j,k
    REAL(prec) :: fac

    ! Now calculate the contravariant basis vectors
    ! In this convention, dsdx(j,i) is contravariant vector i, component j
    ! To project onto contravariant vector i, dot vector along the first dimension
    DO iEl = 1,myGeom % nElem
      DO j = 0,myGeom % dxds % N
        DO i = 0,myGeom % dxds % N

          myGeom % dsdx % interior % hostData(1,1,i,j,1,iEl) = myGeom % dxds % interior % hostData(2,2,i,j,1,iEl)
          myGeom % dsdx % interior % hostData(1,2,i,j,1,iEl) = -myGeom % dxds % interior % hostData(1,2,i,j,1,iEl)
          myGeom % dsdx % interior % hostData(2,1,i,j,1,iEl) = -myGeom % dxds % interior % hostData(2,1,i,j,1,iEl)
          myGeom % dsdx % interior % hostData(2,2,i,j,1,iEl) = myGeom % dxds % interior % hostData(1,1,i,j,1,iEl)

        END DO
      END DO
    END DO

    ! Interpolate the contravariant tensor to the boundaries
    CALL myGeom % dsdx % BoundaryInterp(gpuAccel=.FALSE.)

    ! Now, modify the sign of dsdx so that
    ! myGeom % dsdx % boundary is equal to the outward pointing normal vector
    DO iEl = 1,myGeom % nElem
      DO k = 1,4
        DO i = 0,myGeom % J % N
          IF (k == selfSide2D_East .OR. k == selfSide2D_North) THEN
            fac = SIGN(1.0_prec,myGeom % J % boundary % hostData(i,1,k,iEl))
          ELSE
            fac = -SIGN(1.0_prec,myGeom % J % boundary % hostData(i,1,k,iEl))
          END IF

          myGeom % dsdx % boundary % hostData(1:2,1:2,i,1,k,iEl) = &
            fac*myGeom % dsdx % boundary % hostData(1:2,1:2,i,1,k,iEl)

        END DO
      END DO
    END DO

  END SUBROUTINE CalculateContravariantBasis_Geometry2D

  SUBROUTINE CalculateMetricTerms_Geometry2D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry2D),INTENT(inout) :: myGeom

    CALL myGeom % x % Gradient(myGeom % dxds,gpuAccel=.FALSE.)
    CALL myGeom % dxds % BoundaryInterp(gpuAccel=.FALSE.)

    ! Calculate the Jacobian = determinant of the covariant matrix at each point
    CALL myGeom % dxds % Determinant(myGeom % J)
    CALL myGeom % J % BoundaryInterp(gpuAccel=.FALSE.)

    CALL myGeom % CalculateContravariantBasis()

#ifdef GPU
    CALL myGeom % UpdateDevice()
#endif

  END SUBROUTINE CalculateMetricTerms_Geometry2D

  SUBROUTINE Init_Mesh3D(myMesh,cqType,tqType,nControlPoints,nTargetPoints, &
                         nElem,nSides,nNodes,nUniqueSides,nUniqueNodes,nBCs)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(out) :: myMesh
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: nControlPoints
    INTEGER,INTENT(in) :: nTargetPoints
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nSides
    INTEGER,INTENT(in) :: nNodes
    INTEGER,INTENT(in) :: nUniqueSides
    INTEGER,INTENT(in) :: nUniqueNodes
    INTEGER,INTENT(in) :: nBCs

    myMesh % nElem = nElem
    myMesh % nSides = nSides
    myMesh % nNodes = nNodes
    myMesh % nUniqueSides = nUniqueSides
    myMesh % nUniqueNodes = nUniqueNodes
    myMesh % nBCs = nBCs

    CALL myMesh % geometry % Init(cqType,tqType,nControlPoints,nTargetPoints,nElem)

    CALL myMesh % elemInfo % Alloc(loBound=(/1,1/), &
                                   upBound=(/6,nElem/))

    CALL myMesh % sideInfo % Alloc(loBound=(/1,1/), &
                                   upBound=(/5,nSides/))

    CALL myMesh % nodeCoords % Alloc(loBound=(/1,1/), &
                                     upBound=(/3,nNodes/))

    CALL myMesh % globalNodeIDs % Alloc(loBound=1, &
                                        upBound=nNodes)

    CALL myMesh % BCType % Alloc(loBound=(/1,1/), &
                                 upBound=(/4,nBCs/))

    ALLOCATE (myMesh % BCNames(1:nBCs))

  END SUBROUTINE Init_Mesh3D

  SUBROUTINE Free_Mesh3D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh

    myMesh % nElem = 0
    myMesh % nSides = 0
    myMesh % nNodes = 0
    myMesh % nUniqueSides = 0
    myMesh % nUniqueNodes = 0
    myMesh % nBCs = 0

    CALL myMesh % elemInfo % Free()
    CALL myMesh % sideInfo % Free()
    CALL myMesh % nodeCoords % Free()
    CALL myMesh % globalNodeIDs % Free()
    CALL myMesh % BCType % Free()

    DEALLOCATE (myMesh % BCNames)

  END SUBROUTINE Free_Mesh3D
#ifdef GPU
  SUBROUTINE UpdateHost_Mesh3D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh

    CALL myMesh % geometry % UpdateHost()
    CALL myMesh % elemInfo % UpdateHost()
    CALL myMesh % sideInfo % UpdateHost()
    CALL myMesh % nodeCoords % UpdateHost()
    CALL myMesh % globalNodeIDs % UpdateHost()
    CALL myMesh % BCType % UpdateHost()

  END SUBROUTINE UpdateHost_Mesh3D

  SUBROUTINE UpdateDevice_Mesh3D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh

    CALL myMesh % geometry % UpdateDevice()
    CALL myMesh % elemInfo % UpdateDevice()
    CALL myMesh % sideInfo % UpdateDevice()
    CALL myMesh % nodeCoords % UpdateDevice()
    CALL myMesh % globalNodeIDs % UpdateDevice()
    CALL myMesh % BCType % UpdateDevice()

  END SUBROUTINE UpdateDevice_Mesh3D
#endif
  SUBROUTINE Init_Geometry3D(myGeom,cqType,tqType,nControlPoints,nTargetPoints,nElem)
    IMPLICIT NONE
    CLASS(Geometry3D),INTENT(out) :: myGeom
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: nControlPoints
    INTEGER,INTENT(in) :: nTargetPoints
    INTEGER,INTENT(in) :: nElem

    myGeom % cqType = cqType
    myGeom % tqType = tqType
    myGeom % nElem = nElem

    CALL myGeom % x % Init(N=nControlPoints, &
                           quadratureType=cqType, &
                           M=nTargetPoints, &
                           targetNodeType=tqType, &
                           nVar=1, &
                           nElem=nElem)

    CALL myGeom % dxds % Init(N=nControlPoints, &
                              quadratureType=cqType, &
                              M=nTargetPoints, &
                              targetNodeType=tqType, &
                              nVar=1, &
                              nElem=nElem)

    CALL myGeom % dsdx % Init(N=nControlPoints, &
                              quadratureType=cqType, &
                              M=nTargetPoints, &
                              targetNodeType=tqType, &
                              nVar=1, &
                              nElem=nElem)

    CALL myGeom % J % Init(N=nControlPoints, &
                           quadratureType=cqType, &
                           M=nTargetPoints, &
                           targetNodeType=tqType, &
                           nVar=1, &
                           nElem=nElem)

  END SUBROUTINE Init_Geometry3D

  SUBROUTINE Free_Geometry3D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry3D),INTENT(inout) :: myGeom

    CALL myGeom % x % Free()
    CALL myGeom % dxds % Free()
    CALL myGeom % dsdx % Free()
    CALL myGeom % J % Free()

  END SUBROUTINE Free_Geometry3D
#ifdef GPU
  SUBROUTINE UpdateHost_Geometry3D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry3D),INTENT(inout) :: myGeom

    CALL myGeom % x % UpdateHost()
    CALL myGeom % dxds % UpdateHost()
    CALL myGeom % dsdx % UpdateHost()
    CALL myGeom % J % UpdateHost()

  END SUBROUTINE UpdateHost_Geometry3D

  SUBROUTINE UpdateDevice_Geometry3D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry3D),INTENT(inout) :: myGeom

    CALL myGeom % x % UpdateDevice()
    CALL myGeom % dxds % UpdateDevice()
    CALL myGeom % dsdx % UpdateDevice()
    CALL myGeom % J % UpdateDevice()

  END SUBROUTINE UpdateDevice_Geometry3D
#endif

  SUBROUTINE CalculateContravariantBasis_Geometry3D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry3D),INTENT(inout) :: myGeom
    ! Local
    INTEGER :: iEl,i,j,k
    REAL(prec) :: fac

    ! Now calculate the contravariant basis vectors
    ! In this convention, dsdx(j,i) is contravariant vector i, component j
    ! To project onto contravariant vector i, dot vector along the first dimension
    DO iEl = 1,myGeom % nElem
      DO k = 0,myGeom % dxds % N
        DO j = 0,myGeom % dxds % N
          DO i = 0,myGeom % dxds % N

            myGeom % dsdx % interior % hostData(1,1,i,j,k,1,iEl) = &
              myGeom % dxds % interior % hostData(2,2,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(2,3,i,j,k,1,iEl) - &
              myGeom % dxds % interior % hostData(3,2,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(2,3,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(1,2,i,j,k,1,iEl) = &
              myGeom % dxds % interior % hostData(3,2,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(1,3,i,j,k,1,iEl) - &
              myGeom % dxds % interior % hostData(1,2,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(3,3,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(1,3,i,j,k,1,iEl) = &
              myGeom % dxds % interior % hostData(1,2,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(2,3,i,j,k,1,iEl) - &
              myGeom % dxds % interior % hostData(2,2,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(1,3,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(2,1,i,j,k,1,iEl) = &
              myGeom % dxds % interior % hostData(2,3,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(3,1,i,j,k,1,iEl) - &
              myGeom % dxds % interior % hostData(3,3,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(2,1,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(2,2,i,j,k,1,iEl) = &
              myGeom % dxds % interior % hostData(3,3,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(1,1,i,j,k,1,iEl) - &
              myGeom % dxds % interior % hostData(1,3,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(3,1,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(2,3,i,j,k,1,iEl) = &
              myGeom % dxds % interior % hostData(1,3,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(2,1,i,j,k,1,iEl) - &
              myGeom % dxds % interior % hostData(2,3,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(1,1,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(3,1,i,j,k,1,iEl) = &
              myGeom % dxds % interior % hostData(2,1,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(3,2,i,j,k,1,iEl) - &
              myGeom % dxds % interior % hostData(3,1,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(2,2,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(3,2,i,j,k,1,iEl) = &
              myGeom % dxds % interior % hostData(3,1,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(1,2,i,j,k,1,iEl) - &
              myGeom % dxds % interior % hostData(1,1,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(3,2,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(3,3,i,j,k,1,iEl) = &
              myGeom % dxds % interior % hostData(1,1,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(2,2,i,j,k,1,iEl) - &
              myGeom % dxds % interior % hostData(2,1,i,j,k,1,iEl)* &
              myGeom % dxds % interior % hostData(1,2,i,j,k,1,iEl)

          END DO
        END DO
      END DO
    END DO

    ! Interpolate the contravariant tensor to the boundaries
    CALL myGeom % dsdx % BoundaryInterp(gpuAccel=.FALSE.)

    ! Now, modify the sign of dsdx so that
    ! myGeom % dsdx % boundary is equal to the outward pointing normal vector

    DO iEl = 1,myGeom % nElem
      DO k = 1,6
        DO j = 0,myGeom % J % N
          DO i = 0,myGeom % J % N
            IF (k == selfSide3D_Top .OR. k == selfSide3D_East .OR. k == selfSide3D_North) THEN
              fac = SIGN(1.0_prec,myGeom % J % boundary % hostData(i,j,1,k,iEl))
            ELSE
              fac = -SIGN(1.0_prec,myGeom % J % boundary % hostData(i,j,1,k,iEl))
            END IF

            myGeom % dsdx % boundary % hostData(1:3,1:3,i,j,1,k,iEl) = &
              fac*myGeom % dsdx % boundary % hostData(1:3,1:3,i,j,1,k,iEl)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE CalculateContravariantBasis_Geometry3D

  SUBROUTINE CalculateMetricTerms_Geometry3D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry3D),INTENT(inout) :: myGeom

    CALL myGeom % x % Gradient(myGeom % dxds,gpuAccel=.FALSE.)
    CALL myGeom % dxds % BoundaryInterp(gpuAccel=.FALSE.)

    ! Calculate the Jacobian = determinant of the covariant matrix at each point
    CALL myGeom % dxds % Determinant(myGeom % J)
    CALL myGeom % J % BoundaryInterp(gpuAccel=.FALSE.)

    CALL myGeom % CalculateContravariantBasis()
#ifdef GPU
    CALL myGeom % UpdateDevice()
#endif

  END SUBROUTINE CalculateMetricTerms_Geometry3D

END MODULE SELF_Mesh
