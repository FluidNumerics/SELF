!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_Mesh

  USE SELF_Constants
  USE SELF_HIP
  USE SELF_Lagrange
  USE SELF_Data
  USE SELF_SupportRoutines
  USE SELF_HDF5
  
  ! External Libs !
  USE HDF5

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
!
! Notes:
!  * cornerNode attributes have not been implemented yet
!
!  * For line segments, quads, and hexes, SELF uses Legendre-Gauss-Lobatto quadrature
!
!
! Connectivity information
!
!  sideInfo(1:5,iSide,iEl)
!
!    1 - Side Type
!    2 - Global Side ID
!    3 - Neighbor Element ID
!    4 - 10*( neighbor local side )  + flip
!    5 - Boundary Condition ID
!
!
! ========================================================================= !

  ! Element Types - From Table 4.1 of https://www.hopr-project.org/externals/Meshformat.pdf
  INTEGER,PARAMETER :: selfLineLinear = 1
  INTEGER,PARAMETER :: selfLineNonlinear = 2
  INTEGER,PARAMETER :: selfTriangleLinear = 3
  INTEGER,PARAMETER :: selfQuadLinear = 4
  INTEGER,PARAMETER :: selfQuadBilinear = 14
  INTEGER,PARAMETER :: selfTriangleNonlinear = 23
  INTEGER,PARAMETER :: selfQuadNonlinear = 24
  INTEGER,PARAMETER :: selfTetrahedronLinear = 104
  INTEGER,PARAMETER :: selfPyramidLinear = 105
  INTEGER,PARAMETER :: selfPrismLinear = 106
  INTEGER,PARAMETER :: selfHexahedronLinear = 108
  INTEGER,PARAMETER :: selfPyramidBilinear = 115
  INTEGER,PARAMETER :: selfPrismBilinear = 116
  INTEGER,PARAMETER :: selfHexahedronBilinear = 118
  INTEGER,PARAMETER :: selfTetrahedronNonlinear = 204
  INTEGER,PARAMETER :: selfPyramidNonlinear = 205
  INTEGER,PARAMETER :: selfPrismNonlinear = 206
  INTEGER,PARAMETER :: selfHexahedronNonlinear = 208

  !
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
  !
  INTEGER,PARAMETER :: self_BCDefault = 1
  INTEGER,PARAMETER :: self_nBCsDefault = 5

  !==============================================!
  ! --------------- File Types------------------ !
  !==============================================!
  INTEGER, PARAMETER :: SELF_MESH_ISM_V2_2D = 1
  INTEGER, PARAMETER :: SELF_MESH_ISM_V2_3D = 2
  INTEGER, PARAMETER :: SELF_MESH_HOPR_2D = 3
  INTEGER, PARAMETER :: SELF_MESH_HOPR_3D = 4

  TYPE MeshSpec
    CHARACTER(self_FileNameLength) :: filename
    INTEGER :: fileType

    LOGICAL :: blockMesh
    INTEGER :: blockMesh_nGeo
    INTEGER :: blockMesh_nElemX
    INTEGER :: blockMesh_nElemY
    INTEGER :: blockMesh_nElemZ
    REAL(prec) :: blockMesh_x0,blockMesh_x1
    REAL(prec) :: blockMesh_y0,blockMesh_y1
    REAL(prec) :: blockMesh_z0,blockMesh_z1

  END TYPE MeshSpec

  TYPE MPILayer
    LOGICAL :: mpiEnabled
    INTEGER :: mpiComm
    INTEGER :: mpiPrec
    INTEGER :: rankId
    INTEGER :: nRanks
    INTEGER :: nElem
    INTEGER :: maxMsg
    INTEGER :: msgCount
    TYPE(hfInt32_r1) :: elemToRank
    TYPE(hfInt32_r1) :: offSetElem
    INTEGER, ALLOCATABLE :: requests(:)
    INTEGER, ALLOCATABLE :: stats(:,:)

  CONTAINS

    PROCEDURE :: Init => Init_MPILayer
    PROCEDURE :: Free => Free_MPILayer
    PROCEDURE :: Finalize => Finalize_MPILayer

    PROCEDURE :: GenerateDecomposition => GenerateDecomposition_MPILayer
    PROCEDURE :: SetElemToRank
    PROCEDURE :: SetMaxMsg

    PROCEDURE,PUBLIC :: FinalizeMPIExchangeAsync

    GENERIC, PUBLIC :: GlobalReduce => GlobalReduce_RealScalar
    PROCEDURE, PRIVATE :: GlobalReduce_RealScalar

  END TYPE MPILayer

  TYPE :: SEMMesh
    INTEGER :: nGeo
    INTEGER :: nElem
    INTEGER :: nGlobalElem
    INTEGER :: nNodes
    INTEGER :: nSides
    INTEGER :: nCornerNodes
    INTEGER :: nUniqueNodes
    INTEGER :: nUniqueSides
    INTEGER :: nBCs
    INTEGER :: quadrature

  END TYPE SEMMesh

  TYPE,EXTENDS(SEMMesh) :: Mesh1D
    TYPE(hfInt32_r2) :: elemInfo
    TYPE(hfReal_r1) :: nodeCoords
    TYPE(hfInt32_r1) :: globalNodeIDs
    TYPE(hfInt32_r2) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS
    PROCEDURE,PUBLIC :: Init => Init_Mesh1D
    PROCEDURE,PUBLIC :: Free => Free_Mesh1D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Mesh1D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh1D
    PROCEDURE,PUBLIC :: UniformBlockMesh => UniformBlockMesh_Mesh1D

    PROCEDURE,PUBLIC  :: Write_Mesh => Write_Mesh1D

  END TYPE Mesh1D

  ! Mesh format is set up similar to the HOPr format
  ! See https://hopr-project.org/externals/MeshFormat.pdf

  TYPE,EXTENDS(SEMMesh) :: Mesh2D
    TYPE(hfInt32_r3) :: sideInfo
    TYPE(hfReal_r4) :: nodeCoords
    TYPE(hfInt32_r2) :: elemInfo
    TYPE(hfInt32_r3) :: globalNodeIDs
    TYPE(hfInt32_r2) :: CGNSCornerMap
    TYPE(hfInt32_r2) :: CGNSSideMap
    TYPE(hfInt32_r2) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS
    PROCEDURE,PUBLIC :: Init => Init_Mesh2D
    PROCEDURE,PUBLIC :: Free => Free_Mesh2D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Mesh2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh2D

    PROCEDURE,PUBLIC :: ResetBoundaryConditionType => ResetBoundaryConditionType_Mesh2D

    PROCEDURE,PUBLIC :: Read_HOPr => Read_HOPr_Mesh2D

    PROCEDURE,PUBLIC :: Write_Mesh => Write_Mesh2D

    PROCEDURE,PRIVATE :: RecalculateFlip => RecalculateFlip_Mesh2D

  END TYPE Mesh2D

  TYPE,EXTENDS(SEMMesh) :: Mesh3D
    TYPE(hfInt32_r3) :: sideInfo
    TYPE(hfReal_r5) :: nodeCoords
    TYPE(hfInt32_r2) :: elemInfo
    TYPE(hfInt32_r4) :: globalNodeIDs
    TYPE(hfInt32_r2) :: CGNSCornerMap
    TYPE(hfInt32_r2) :: CGNSSideMap
    TYPE(hfInt32_r2) :: sideMap
    TYPE(hfInt32_r2) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Mesh3D
    PROCEDURE,PUBLIC :: Free => Free_Mesh3D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Mesh3D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh3D

    PROCEDURE,PUBLIC :: Read_HOPr => Read_HOPr_Mesh3D

    PROCEDURE,PUBLIC :: ResetBoundaryConditionType => ResetBoundaryConditionType_Mesh3D

    PROCEDURE,PUBLIC :: Write_Mesh => Write_Mesh3D

    PROCEDURE,PRIVATE :: RecalculateFlip => RecalculateFlip_Mesh3D

  END TYPE Mesh3D

CONTAINS

  SUBROUTINE Init_Mesh1D(myMesh,nGeo,nElem,nNodes,nBCs)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(out) :: myMesh
    INTEGER,INTENT(in) :: nGeo
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nNodes
    INTEGER,INTENT(in) :: nBCs

    myMesh % nGeo = nGeo
    myMesh % nElem = nElem
    myMesh % nGlobalElem = nElem
    myMesh % nNodes = nNodes
    myMesh % nCornerNodes = nElem*2
    myMesh % nUniqueNodes = 0
    myMesh % nBCs = nBCs

    CALL myMesh % elemInfo % Alloc(loBound=(/1,1/), &
                                        upBound=(/4,nElem/))

    CALL myMesh % nodeCoords % Alloc(loBound=1, &
                                          upBound=nNodes)

    CALL myMesh % globalNodeIDs % Alloc(loBound=1, &
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
    myMesh % nCornerNodes = 0
    myMesh % nUniqueNodes = 0
    myMesh % nBCs = 0

    CALL myMesh % elemInfo % Free()
    CALL myMesh % nodeCoords % Free()
    CALL myMesh % globalNodeIDs % Free()
    CALL myMesh % BCType % Free()

    DEALLOCATE (myMesh % BCNames)

  END SUBROUTINE Free_Mesh1D

  SUBROUTINE UpdateHost_Mesh1D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(inout) :: myMesh

    CALL myMesh % elemInfo % UpdateHost()
    CALL myMesh % nodeCoords % UpdateHost()
    CALL myMesh % globalNodeIDs % UpdateHost()
    CALL myMesh % BCType % UpdateHost()

  END SUBROUTINE UpdateHost_Mesh1D

  SUBROUTINE UpdateDevice_Mesh1D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(inout) :: myMesh

    CALL myMesh % elemInfo % UpdateDevice()
    CALL myMesh % nodeCoords % UpdateDevice()
    CALL myMesh % globalNodeIDs % UpdateDevice()
    CALL myMesh % BCType % UpdateDevice()

  END SUBROUTINE UpdateDevice_Mesh1D

  SUBROUTINE UniformBlockMesh_Mesh1D(myMesh,nGeo,nElem,x)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(out) :: myMesh
    INTEGER,INTENT(in) :: nGeo
    INTEGER,INTENT(in) :: nElem
    REAL(prec),INTENT(in) :: x(1:2)
    ! Local
    INTEGER :: iel
    INTEGER :: nid,nNodes
    INTEGER :: i
    REAL(prec) :: xU(1:nElem + 1)
    TYPE(Lagrange), TARGET :: linearInterp
    TYPE(Lagrange), TARGET :: nGeoInterp
    TYPE(Scalar1D) :: xLinear
    TYPE(Scalar1D) :: xGeo

    nNodes = nElem*(nGeo + 1)
    CALL myMesh % Init(nGeo,nElem,nNodes,2)
    myMesh % quadrature = GAUSS_LOBATTO

    ! Set the hopr_nodeCoords
    xU = UniformPoints(x(1),x(2),1,nElem + 1)

    CALL linearInterp % Init(1,GAUSS_LOBATTO,&
            nGeo,GAUSS_LOBATTO)

    CALL nGeoInterp % Init(nGeo,GAUSS_LOBATTO,&
            nGeo,GAUSS_LOBATTO)

    ! Create a linear interpolant to interpolate to nGeo grid
    CALL xLinear % Init(linearInterp,1,nElem)
    CALL xGeo % Init(nGeoInterp,1,nElem)

    DO iel = 1,nElem
      xLinear % interior % hostData(0:1,1,iel) = xU(iel:iel + 1)
    END DO

    CALL xLinear % GridInterp(xGeo,.FALSE.)

    ! Set the element information
    nid = 1
    DO iel = 1,nElem
      myMesh % elemInfo % hostData(1,iel) = selfLineLinear ! Element Type
      myMesh % elemInfo % hostData(2,iel) = 1 ! Element Zone
      myMesh % elemInfo % hostData(3,iel) = nid ! Node Index Start
      DO i = 0,nGeo
        myMesh % nodeCoords % hostData(nid) = xGeo % interior % hostData(i,1,iel)
        nid = nid + 1
      END DO
      myMesh % elemInfo % hostData(4,iel) = nid - 1 ! Node Index End
    END DO


    CALL myMesh % UpdateDevice()

    CALL xLinear % Free()
    CALL xGeo % Free()
    CALL linearInterp % Free()
    CALL nGeoInterp % Free()

  END SUBROUTINE UniformBlockMesh_Mesh1D

  SUBROUTINE Write_Mesh1D(myMesh,meshFile)
    ! Writes mesh output in HOPR format (serial IO only)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(inout) :: myMesh
    CHARACTER(*),INTENT(in) :: meshFile
    ! Local
    INTEGER(HID_T) :: fileId

    CALL Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId)

    CALL WriteAttribute_HDF5(fileId,'nElems',myMesh % nElem)
    CALL WriteAttribute_HDF5(fileId,'Ngeo',myMesh % nGeo)
    CALL WriteAttribute_HDF5(fileId,'nBCs',myMesh % nBCs)

    CALL WriteArray_HDF5(fileId,'BCType',myMesh % bcType)

    ! Read local subarray of ElemInfo
    CALL WriteArray_HDF5(fileId,'ElemInfo',myMesh % elemInfo)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    CALL WriteArray_HDF5(fileId,'NodeCoords',myMesh % nodeCoords)
    CALL WriteArray_HDF5(fileId,'GlobalNodeIDs',myMesh % globalNodeIDs)

    CALL Close_HDF5(fileID)

  END SUBROUTINE Write_Mesh1D

  SUBROUTINE Init_Mesh2D(myMesh,nGeo,nElem,nSides,nNodes,nBCs)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(out) :: myMesh
    INTEGER,INTENT(in) :: nGeo
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nSides
    INTEGER,INTENT(in) :: nNodes
    INTEGER,INTENT(in) :: nBCs
    ! Local
    INTEGER :: i,j,l

    myMesh % nGeo = nGeo
    myMesh % nElem = nElem
    myMesh % nGlobalElem = nElem
    myMesh % nNodes = nNodes
    myMesh % nSides = nSides
    myMesh % nCornerNodes = 0
    myMesh % nUniqueNodes = 0
    myMesh % nUniqueSides = 0
    myMesh % nBCs = nBCs

    CALL myMesh % elemInfo % Alloc(loBound=(/1,1/), &
                                        upBound=(/6,nElem/))

    CALL myMesh % sideInfo % Alloc(loBound=(/1,1,1/), &
                                        upBound=(/5,4,nElem/))

    CALL myMesh % nodeCoords % Alloc(loBound=(/1,0,0,1/), &
                                          upBound=(/2,nGeo,nGeo,nElem/))

    CALL myMesh % globalNodeIDs % Alloc(loBound=(/0,0,1/), &
                                        upBound=(/nGeo,nGeo,nElem/))

    CALL myMesh % CGNSCornerMap % Alloc(loBound=(/1,1/), &
                                        upBound=(/2,4/))

    CALL myMesh % CGNSSideMap % Alloc(loBound=(/1,1/), &
                                      upBound=(/2,4/))

    CALL myMesh % BCType % Alloc(loBound=(/1,1/), &
                                 upBound=(/4,nBCs/))

    ALLOCATE (myMesh % BCNames(1:nBCs))

    ! Create lookup tables to assist with connectivity generation
    myMesh % CGNSCornerMap % hostData(1:2,1) = (/0, 0/)
    myMesh % CGNSCornerMap % hostData(1:2,2) = (/nGeo, 0/)
    myMesh % CGNSCornerMap % hostData(1:2,3) = (/nGeo, nGeo/)
    myMesh % CGNSCornerMap % hostData(1:2,4) = (/0, nGeo/)

    ! Maps from local corner node id to CGNS side
    myMesh % CGNSSideMap % hostData(1:2,1) = (/1,2/)
    myMesh % CGNSSideMap % hostData(1:2,2) = (/2,3/)
    myMesh % CGNSSideMap % hostData(1:2,3) = (/4,3/)
    myMesh % CGNSSideMap % hostData(1:2,4) = (/1,4/)

  END SUBROUTINE Init_Mesh2D

  SUBROUTINE Free_Mesh2D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh

    myMesh % nElem = 0
    myMesh % nNodes = 0
    myMesh % nSides = 0
    myMesh % nCornerNodes = 0
    myMesh % nUniqueSides = 0
    myMesh % nUniqueNodes = 0
    myMesh % nBCs = 0

    CALL myMesh % elemInfo % Free()
    CALL myMesh % sideInfo % Free()
    CALL myMesh % nodeCoords % Free()
    CALL myMesh % CGNSCornerMap % Free()
    CALL myMesh % globalNodeIDs % Free()
    CALL myMesh % BCType % Free()

    DEALLOCATE (myMesh % BCNames)

  END SUBROUTINE Free_Mesh2D

  SUBROUTINE UpdateHost_Mesh2D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh

    CALL myMesh % elemInfo % UpdateHost()
    CALL myMesh % sideInfo % UpdateHost()
    CALL myMesh % nodeCoords % UpdateHost()
    CALL myMesh % globalNodeIDs % UpdateHost()
    CALL myMesh % BCType % UpdateHost()

  END SUBROUTINE UpdateHost_Mesh2D

  SUBROUTINE UpdateDevice_Mesh2D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh

    CALL myMesh % elemInfo % UpdateDevice()
    CALL myMesh % sideInfo % UpdateDevice()
    CALL myMesh % nodeCoords % UpdateDevice()
    CALL myMesh % globalNodeIDs % UpdateDevice()
    CALL myMesh % BCType % UpdateDevice()

  END SUBROUTINE UpdateDevice_Mesh2D

  SUBROUTINE ResetBoundaryConditionType_Mesh2D(myMesh,bcid)
    !! This method can be used to reset all of the boundary elements
    !! boundary condition type to the desired value.
    !!
    !! Note that ALL physical boundaries will be set to have this boundary 
    !! condition
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh
    INTEGER, INTENT(in) :: bcid
    ! Local
    INTEGER :: iSide,iEl,e2      

    DO iEl = 1, myMesh % nElem
      DO iSide = 1, 4

        e2 = myMesh % sideInfo % hostData(3,iSide,iEl)

        IF( e2 == 0 )THEN
          myMesh % sideInfo % hostData(5,iSide,iEl) = bcid
        ENDIF

      ENDDO
    ENDDO

  END SUBROUTINE ResetBoundaryConditionType_Mesh2D 

  SUBROUTINE Read_HOPr_Mesh2D(myMesh,meshFile,decomp)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    ! Adapted for 2D Mesh : Note that HOPR does not have 2D mesh output.
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(out) :: myMesh
    CHARACTER(*),INTENT(in) :: meshFile
    TYPE(MPILayer),INTENT(inout) :: decomp
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: offset(1:2),gOffset(1)
    INTEGER :: nGlobalElem
    INTEGER :: firstElem
    INTEGER :: firstNode
    INTEGER :: firstSide
    INTEGER :: nLocalElems
    INTEGER :: nLocalNodes3D
    INTEGER :: nLocalSides3D
    INTEGER :: nUniqueSides3D
    INTEGER :: nLocalNodes2D
    INTEGER :: nLocalSides2D
    INTEGER :: nUniqueSides2D
    INTEGER :: nGeo,nBCs
    INTEGER :: eid, lsid, iSide
    INTEGER :: i, j, nid
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfInt32_r2) :: hopr_sideInfo
    TYPE(hfReal_r2) :: hopr_nodeCoords
    TYPE(hfInt32_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: bcType


    IF ( decomp % mpiEnabled )THEN
      CALL Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId,decomp % mpiComm)
    ELSE
      CALL Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId)
    ENDIF

    CALL ReadAttribute_HDF5(fileId,'nElems',nGlobalElem)
    CALL ReadAttribute_HDF5(fileId,'Ngeo',nGeo)
    CALL ReadAttribute_HDF5(fileId,'nBCs',nBCs)
    CALL ReadAttribute_HDF5(fileId,'nUniqueSides',nUniqueSides3D)

    ! Read BCType
    CALL bcType % Alloc(loBound=(/1,1/), &
                        upBound=(/4,nBCs/))
    IF ( decomp % mpiEnabled )THEN
      offset(:) = 0
      CALL ReadArray_HDF5(fileId,'BCType',bcType,offset)
    ELSE
      CALL ReadArray_HDF5(fileId,'BCType',bcType)
    ENDIF

    ! Read local subarray of ElemInfo
    CALL decomp % GenerateDecomposition(nGlobalElem,nUniqueSides3D)

    firstElem = decomp % offsetElem % hostData(decomp % rankId) + 1
    nLocalElems = decomp % offsetElem % hostData(decomp % rankId + 1) - &
                  decomp % offsetElem % hostData(decomp % rankId)

    ! Allocate Space for hopr_elemInfo!
    CALL hopr_elemInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/6,nLocalElems/))

    IF ( decomp % mpiEnabled )THEN
      offset = (/0,firstElem - 1/)
      CALL ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo,offset)
    ELSE
      CALL ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo)
    ENDIF

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode = hopr_elemInfo % hostData(5,1) + 1
    nLocalNodes3D = hopr_elemInfo % hostData(6,nLocalElems) - hopr_elemInfo % hostData(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !
    CALL hopr_nodeCoords % Alloc(loBound=(/1,1/), &
                                 upBound=(/3,nLocalNodes3D/))

    CALL hopr_globalNodeIDs % Alloc(loBound=1, &
                                    upBound=nLocalNodes3D)

    IF ( decomp % mpiEnabled )THEN
      offset = (/0,firstNode - 1/)
      CALL ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords,offset)
      gOffset = (/firstNode - 1/)
      CALL ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs,gOffset)
    ELSE
      CALL ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords)
      CALL ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs)
    ENDIF

    ! Read local subarray of SideInfo
    firstSide = hopr_elemInfo % hostData(3,1) + 1
    nLocalSides3D = hopr_elemInfo % hostData(4,nLocalElems) - hopr_elemInfo % hostData(3,1)

    ! Allocate space for hopr_sideInfo
    CALL hopr_sideInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/5,nLocalSides3D/))
    IF ( decomp % mpiEnabled )THEN
      offset = (/0,firstSide - 1/)
      CALL ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo,offset)
    ELSE
      CALL ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo)
    ENDIF

    CALL Close_HDF5(fileID)
    ! ---- Done reading 3-D Mesh information ---- !

    ! Now we need to convert from 3-D to 2-D !
    nLocalSides2D = nLocalSides3D - 2*nGlobalElem
    nUniqueSides2D = nUniqueSides3D - 2*nGlobalElem ! Remove the "top" and "bottom" faces
    nLocalNodes2D = nLocalNodes2D - nGlobalElem*nGeo*(nGeo+1)**2 ! Remove the third dimension

    CALL myMesh % Init(nGeo,nLocalElems,nLocalSides2D,nLocalNodes2D,nBCs)

    ! Copy data from local arrays into myMesh
    !  elemInfo(1:6,iEl)
    !    1 - Element Type
    !    2 - Zone
    !    3 - offset index for side array (not needed when all quads are assumed)
    !    4 - last index for side array (not needed when all quads are assumed)
    !    5 - offset index for node array (not needed when all quads are assumed)
    !    6 - last index for node array (not needed when all quads are assumed)
    myMesh % elemInfo % hostData = hopr_elemInfo % hostData
    myMesh % quadrature = UNIFORM  ! HOPr uses uniformly spaced points

    ! Grab the node coordinates (x and y only) from the "bottom" layer of the extruded mesh
    DO eid = 1, myMesh % nElem
      DO j = 0, nGeo
        DO i = 0, nGeo
          nid = i+1 + (nGeo+1)*(j + (nGeo+1)*((nGeo+1)*(eid-1)))
          myMesh % nodeCoords % hostData(1:2,i,j,eid) = hopr_nodeCoords % hostData(1:2,nid)
          myMesh % globalNodeIDs % hostData(i,j,eid) = hopr_globalNodeIDs % hostData(nid)
        ENDDO
      ENDDO
    ENDDO

    ! Grab the south, west, north, and south sides of the elements 
    !  sideInfo(1:5,iSide,iEl)
    !
    !    1 - Side Type (currently unused in SELF)
    !    2 - Global Side ID (Used for message passing. Don't need to change)
    !    3 - Neighbor Element ID (Can stay the same)
    !    4 - 10*( neighbor local side )  + flip (Need to recalculate flip)
    !    5 - Boundary Condition ID (Can stay the same)
    DO eid = 1,myMesh % nElem
      DO lsid = 1,4
        ! Calculate the 3-D side ID from the 2-D local side id and element ID
        iSide = lsid + 1 + 6*(eid-1)
        myMesh % sideInfo % hostData(1:5,lsid,eid) = hopr_sideInfo % hostData(1:5,iSide)
        ! Adjust the secondary side index for 2-D
        myMesh % sideInfo % hostData(4,lsid,eid) = myMesh % sideInfo % hostData(4,lsid,eid)-10
      ENDDO
    ENDDO

    CALL myMesh % RecalculateFlip()

    CALL myMesh % UpdateDevice()

    CALL hopr_elemInfo % Free()
    CALL hopr_nodeCoords % Free()
    CALL hopr_globalNodeIDs % Free()
    CALL hopr_sideInfo % Free()

  END SUBROUTINE Read_HOPr_Mesh2D

  SUBROUTINE RecalculateFlip_Mesh2D(myMesh,decomp)  
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh
    TYPE(MPILayer),INTENT(inout),OPTIONAL :: decomp
    ! Local
    INTEGER :: e1
    INTEGER :: s1
    INTEGER :: e2
    INTEGER :: e2Global
    INTEGER :: s2
    INTEGER :: flip
    INTEGER :: bcid
    INTEGER :: lnid1(1:2)
    INTEGER :: lnid2(1:2)
    INTEGER :: nid1(1:2,1:4,1:myMesh % nElem)
    INTEGER :: nid2(1:2,1:4,1:myMesh % nElem)
    INTEGER :: nloc1(1:2)
    INTEGER :: nloc2(1:2)
    INTEGER :: n1
    INTEGER :: n1Global
    INTEGER :: n2
    INTEGER :: n2Global
    INTEGER :: c1
    INTEGER :: c2
    INTEGER :: i, j
    INTEGER :: l
    INTEGER :: nShifts
    INTEGER :: neighborRank
    INTEGER :: rankId
    INTEGER :: offset
    INTEGER :: msgCount
    INTEGER :: globalSideId
    INTEGER, ALLOCATABLE :: requests(:)
    INTEGER, ALLOCATABLE :: stats(:,:)
    INTEGER :: iError
    LOGICAL :: theyMatch


    ALLOCATE(requests(1:myMesh % nSides*2))
    ALLOCATE(stats(MPI_STATUS_SIZE,1:myMesh % nSides*2))

    IF (PRESENT(decomp)) THEN
      rankId = decomp % rankId
      offset = decomp % offsetElem % hostData(rankId)
    ELSE
      rankId = 0
      offset = 0
    ENDIF

    msgCount = 0
    DO e1 = 1,myMesh % nElem
      DO s1 = 1,4

        e2Global = myMesh % sideInfo % hostData(3,s1,e1)
        e2 = e2Global - offset
        s2 = myMesh % sideInfo % hostData(4,s1,e1)/10
        flip = myMesh % sideInfo % hostData(4,s1,e1) - s2*10
        bcid = myMesh % sideInfo % hostData(5,s1,e1)

        IF (bcid == 0) THEN

          IF (PRESENT(decomp)) THEN
            neighborRank = decomp % elemToRank % hostData(e2Global)
          ELSE
            neighborRank = 0
          ENDIF

          IF (neighborRank == rankId) THEN

            lnid1 = myMesh % CGNSSideMap % hostData(1:2,s1) ! local CGNS corner node ids for element 1 side
            lnid2 = myMesh % CGNSSideMap % hostData(1:2,s2) ! local CGNS corner node ids for element 2 side

            DO l = 1, 2

              i = myMesh % CGNSCornerMap % hostData(1,lnid1(l))
              j = myMesh % CGNSCornerMap % hostData(2,lnid1(l))
              nid1(l,s1,e1) = myMesh % globalNodeIDs % hostData(i,j,e1)

              i = myMesh % CGNSCornerMap % hostData(1,lnid2(l))
              j = myMesh % CGNSCornerMap % hostData(2,lnid2(l))
              nid2(l,s1,e1) = myMesh % globalNodeIDs % hostData(i,j,e2)

            ENDDO

          ELSE ! In this case, we need to exchange

            globalSideId = ABS(myMesh % sideInfo % hostdata(2,s1,e1))

            lnid1 = myMesh % CGNSSideMap % hostData(1:2,s1) ! local CGNS corner node ids for element 1 side

            DO l = 1, 2

              i = myMesh % CGNSCornerMap % hostData(1,lnid1(l))
              j = myMesh % CGNSCornerMap % hostData(2,lnid1(l))
              nid1(l,s1,e1) = myMesh % globalNodeIDs % hostData(i,j,e1)

              msgCount = msgCount + 1
              CALL MPI_IRECV(nid2(l,s1,e1), &
                             1, &
                             MPI_INTEGER, &
                             neighborRank,globalSideId, &
                             decomp % mpiComm, &
                             requests(msgCount),iError)
  
              ! Send nid1(l) from this rank to nid2(l) on the other rank
              msgCount = msgCount + 1
              CALL MPI_ISEND(nid1(l,s1,e1), &
                             1, &
                             MPI_INTEGER, &
                             neighborRank,globalSideId, &
                             decomp % mpiComm, &
                             requests(msgCount),iError)
  

            ENDDO

          ENDIF ! MPI or not

        ENDIF ! If not physical boundary

      ENDDO
    ENDDO

    IF (PRESENT(decomp) .AND. msgCount > 0) THEN
      CALL MPI_WaitAll(msgCount, &
                       requests(1:msgCount), &
                       stats(1:MPI_STATUS_SIZE,1:msgCount), &
                       iError)
    ENDIF

    DO e1 = 1,myMesh % nElem
      DO s1 = 1,4

        s2 = myMesh % sideInfo % hostData(4,s1,e1)/10
        bcid = myMesh % sideInfo % hostData(5,s1,e1)
        nloc1(1:2) = nid1(1:2,s1,e1)
        nloc2(1:2) = nid2(1:2,s1,e1)

        IF (bcid == 0) THEN
          theyMatch = CompareArray( nloc1, nloc2, 2 )

          IF( theyMatch )THEN
            myMesh % sideInfo % hostData(4,s1,e1) = 10*s2
          ELSE
            myMesh % sideInfo % hostData(4,s1,e1) = 10*s2+1
          ENDIF


        ENDIF

      ENDDO
    ENDDO

    DEALLOCATE(requests)
    DEALLOCATE(stats)

  END SUBROUTINE RecalculateFlip_Mesh2D

  SUBROUTINE Write_Mesh2D(myMesh,meshFile)
    ! Writes mesh output in HOPR format (serial only)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh
    CHARACTER(*),INTENT(in) :: meshFile
    ! Local
    INTEGER(HID_T) :: fileId

    CALL Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId)
    CALL WriteAttribute_HDF5(fileId,'nElems',myMesh % nElem)
    CALL WriteAttribute_HDF5(fileId,'Ngeo',myMesh % nGeo)
    CALL WriteAttribute_HDF5(fileId,'nBCs',myMesh % nBCs)

    CALL WriteArray_HDF5(fileId,'BCType',myMesh % bcType)

    ! Write local subarray of ElemInfo
    CALL WriteArray_HDF5(fileId,'ElemInfo',myMesh % elemInfo)

    ! Write local subarray of NodeCoords and GlobalNodeIDs
    CALL WriteArray_HDF5(fileId,'NodeCoords',myMesh % nodeCoords)
    CALL WriteArray_HDF5(fileId,'GlobalNodeIDs',myMesh % globalNodeIDs)

    ! Write local subarray of SideInfo
    CALL WriteArray_HDF5(fileId,'SideInfo',myMesh % sideInfo)

    CALL Close_HDF5(fileID)

  END SUBROUTINE Write_Mesh2D

  SUBROUTINE Init_Mesh3D(myMesh,nGeo,nElem,nSides,nNodes,nBCs)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(out) :: myMesh
    INTEGER,INTENT(in) :: nGeo
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nSides
    INTEGER,INTENT(in) :: nNodes
    INTEGER,INTENT(in) :: nBCs
    ! Local
    INTEGER :: i,j,k,l

    myMesh % nElem = nElem
    myMesh % nGlobalElem = nElem
    myMesh % nGeo = nGeo
    myMesh % nSides = nSides
    myMesh % nNodes = nNodes
    myMesh % nCornerNodes = 0
    myMesh % nUniqueSides = 0
    myMesh % nUniqueNodes = 0
    myMesh % nBCs = nBCs

    CALL myMesh % elemInfo % Alloc(loBound=(/1,1/), &
                                   upBound=(/6,nElem/))

    CALL myMesh % sideInfo % Alloc(loBound=(/1,1,1/), &
                                   upBound=(/5,6,nElem/))

    CALL myMesh % nodeCoords % Alloc(loBound=(/1,0,0,0,1/), &
                                     upBound=(/3,nGeo,nGeo,nGeo,nElem/))

    CALL myMesh % globalNodeIDs % Alloc(loBound=(/0,0,0,1/), &
                                        upBound=(/nGeo,nGeo,nGeo,nElem/))

    CALL myMesh % CGNSCornerMap % Alloc(loBound=(/1,1/), &
                                        upBound=(/3,8/))

    CALL myMesh % sideMap % Alloc(loBound=(/1,1/), &
                                  upBound=(/4,6/))

    CALL myMesh % CGNSSideMap % Alloc(loBound=(/1,1/), &
                                      upBound=(/4,6/))

    CALL myMesh % BCType % Alloc(loBound=(/1,1/), &
                                 upBound=(/4,nBCs/))

    ALLOCATE (myMesh % BCNames(1:nBCs))

    ! Create lookup tables to assist with connectivity generation
    myMesh % CGNSCornerMap % hostData(1:3,1) = (/0,0,0/) ! Bottom-South-West
    myMesh % CGNSCornerMap % hostData(1:3,2) = (/nGeo,0,0/) ! Bottom-South-East
    myMesh % CGNSCornerMap % hostData(1:3,3) = (/nGeo,nGeo,0/)! Bottom-North-East
    myMesh % CGNSCornerMap % hostData(1:3,4) = (/0,nGeo,0/)! Bottom-North-West
    myMesh % CGNSCornerMap % hostData(1:3,5) = (/0,0,nGeo/) ! Top-South-West
    myMesh % CGNSCornerMap % hostData(1:3,6) = (/nGeo,0,nGeo/) ! Top-South-East
    myMesh % CGNSCornerMap % hostData(1:3,7) = (/nGeo,nGeo,nGeo/)! Top-North-East
    myMesh % CGNSCornerMap % hostData(1:3,8) = (/0,nGeo,nGeo/)! Top-North-West

    ! Maps from local corner node id to CGNS side
    myMesh % CGNSSideMap % hostData(1:4,1) = (/1,4,3,2/)
    myMesh % CGNSSideMap % hostData(1:4,2) = (/1,2,6,5/)
    myMesh % CGNSSideMap % hostData(1:4,3) = (/2,3,7,6/)
    myMesh % CGNSSideMap % hostData(1:4,4) = (/3,4,8,7/)
    myMesh % CGNSSideMap % hostData(1:4,5) = (/1,5,8,4/)
    myMesh % CGNSSideMap % hostData(1:4,6) = (/5,6,7,8/)

    myMesh % sideMap % hostData(1:4,1) = (/1,2,3,4/) ! Bottom
    myMesh % sideMap % hostData(1:4,2) = (/1,2,6,5/) ! South
    myMesh % sideMap % hostData(1:4,3) = (/2,3,7,6/) ! East
    myMesh % sideMap % hostData(1:4,4) = (/4,3,7,8/) ! North
    myMesh % sideMap % hostData(1:4,5) = (/1,4,8,5/) ! West
    myMesh % sideMap % hostData(1:4,6) = (/5,6,7,8/) ! Top

  END SUBROUTINE Init_Mesh3D

  SUBROUTINE Free_Mesh3D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh

    myMesh % nElem = 0
    myMesh % nSides = 0
    myMesh % nNodes = 0
    myMesh % nCornerNodes = 0
    myMesh % nUniqueSides = 0
    myMesh % nUniqueNodes = 0
    myMesh % nBCs = 0

    CALL myMesh % elemInfo % Free()
    CALL myMesh % sideInfo % Free()
    CALL myMesh % nodeCoords % Free()
    CALL myMesh % CGNSCornerMap % Free()
    CALL myMesh % globalNodeIDs % Free()
    CALL myMesh % sideMap % Free()
    CALL myMesh % BCType % Free()

    DEALLOCATE (myMesh % BCNames)

  END SUBROUTINE Free_Mesh3D
  SUBROUTINE UpdateHost_Mesh3D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh

    CALL myMesh % elemInfo % UpdateHost()
    CALL myMesh % sideInfo % UpdateHost()
    CALL myMesh % nodeCoords % UpdateHost()
    CALL myMesh % globalNodeIDs % UpdateHost()
    CALL myMesh % BCType % UpdateHost()

  END SUBROUTINE UpdateHost_Mesh3D

  SUBROUTINE UpdateDevice_Mesh3D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh

    CALL myMesh % elemInfo % UpdateDevice()
    CALL myMesh % sideInfo % UpdateDevice()
    CALL myMesh % nodeCoords % UpdateDevice()
    CALL myMesh % globalNodeIDs % UpdateDevice()
    CALL myMesh % BCType % UpdateDevice()

  END SUBROUTINE UpdateDevice_Mesh3D

  SUBROUTINE ResetBoundaryConditionType_Mesh3D(myMesh,bcid)
    !! This method can be used to reset all of the boundary elements
    !! boundary condition type to the desired value.
    !!
    !! Note that ALL physical boundaries will be set to have this boundary 
    !! condition
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh
    INTEGER, INTENT(in) :: bcid
    ! Local
    INTEGER :: iSide,iEl,e2      

    DO iEl = 1, myMesh % nElem
      DO iSide = 1, 6

        e2 = myMesh % sideInfo % hostData(3,iSide,iEl)

        IF( e2 == 0 )THEN
          myMesh % sideInfo % hostData(5,iSide,iEl) = bcid
        ENDIF

      ENDDO
    ENDDO

  END SUBROUTINE ResetBoundaryConditionType_Mesh3D 

  SUBROUTINE RecalculateFlip_Mesh3D(myMesh,decomp)  
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh
    TYPE(MPILayer),INTENT(inout),OPTIONAL :: decomp
    ! Local
    INTEGER :: e1
    INTEGER :: s1
    INTEGER :: e2
    INTEGER :: e2Global
    INTEGER :: s2
    INTEGER :: flip
    INTEGER :: bcid
    INTEGER :: lnid1(1:4)
    INTEGER :: lnid2(1:4)
    INTEGER :: nid1(1:4,1:6,1:myMesh % nElem)
    INTEGER :: nid2(1:4,1:6,1:myMesh % nElem)
    INTEGER :: nloc1(1:4)
    INTEGER :: nloc2(1:4)
    INTEGER :: n1
    INTEGER :: n1Global
    INTEGER :: n2
    INTEGER :: n2Global
    INTEGER :: c1
    INTEGER :: c2
    INTEGER :: i,j,k
    INTEGER :: l
    INTEGER :: nShifts
    INTEGER :: neighborRank
    INTEGER :: rankId
    INTEGER :: offset
    INTEGER :: msgCount
    INTEGER :: globalSideId
    INTEGER, ALLOCATABLE :: requests(:)
    INTEGER, ALLOCATABLE :: stats(:,:)
    INTEGER :: iError
    LOGICAL :: theyMatch


    ALLOCATE(requests(1:myMesh % nSides*2))
    ALLOCATE(stats(MPI_STATUS_SIZE,1:myMesh % nSides*2))

    IF (PRESENT(decomp)) THEN
      rankId = decomp % rankId
      offset = decomp % offsetElem % hostData(rankId)
    ELSE
      rankId = 0
      offset = 0
    ENDIF

    msgCount = 0
    DO e1 = 1,myMesh % nElem
      DO s1 = 1,6

        e2Global = myMesh % sideInfo % hostData(3,s1,e1)
        e2 = e2Global - offset
        s2 = myMesh % sideInfo % hostData(4,s1,e1)/10
        flip = myMesh % sideInfo % hostData(4,s1,e1) - s2*10
        bcid = myMesh % sideInfo % hostData(5,s1,e1)

        IF (bcid == 0) THEN

          IF (PRESENT(decomp)) THEN
            neighborRank = decomp % elemToRank % hostData(e2Global)
          ELSE
            neighborRank = 0
          ENDIF

          IF (neighborRank == rankId) THEN

            lnid1 = myMesh % sideMap % hostData(1:4,s1) ! local CGNS corner node ids for element 1 side
            lnid2 = myMesh % sideMap % hostData(1:4,s2) ! local CGNS corner node ids for element 2 side

            DO l = 1, 4
              
              i = myMesh % CGNSCornerMap % hostData(1,lnid1(l))
              j = myMesh % CGNSCornerMap % hostData(2,lnid1(l))
              k = myMesh % CGNSCornerMap % hostData(3,lnid1(l))
              nid1(l,s1,e1) = myMesh % globalNodeIDs % hostData(i,j,k,e1)

              i = myMesh % CGNSCornerMap % hostData(1,lnid2(l))
              j = myMesh % CGNSCornerMap % hostData(2,lnid2(l))
              k = myMesh % CGNSCornerMap % hostData(3,lnid2(l))
              nid2(l,s1,e1) = myMesh % globalNodeIDs % hostData(i,j,k,e1)

            ENDDO

          ELSE ! In this case, we need to exchange

            globalSideId = ABS(myMesh % sideInfo % hostdata(2,s1,e1))

            lnid1 = myMesh % sideMap % hostData(1:4,s1) ! local CGNS corner node ids for element 1 side

            DO l = 1, 4
              
              i = myMesh % CGNSCornerMap % hostData(1,lnid1(l))
              j = myMesh % CGNSCornerMap % hostData(2,lnid1(l))
              k = myMesh % CGNSCornerMap % hostData(3,lnid1(l))
              nid1(l,s1,e1) = myMesh % globalNodeIDs % hostData(i,j,k,e1)

              ! Receive nid2(l) on this rank from  nid1(l) on the other rank
              msgCount = msgCount + 1
              CALL MPI_IRECV(nid2(l,s1,e1), &
                             1, &
                             MPI_INTEGER, &
                             neighborRank,globalSideId, &
                             decomp % mpiComm, &
                             requests(msgCount),iError)
  
              ! Send nid1(l) from this rank to nid2(l) on the other rank
              msgCount = msgCount + 1
              CALL MPI_ISEND(nid1(l,s1,e1), &
                             1, &
                             MPI_INTEGER, &
                             neighborRank,globalSideId, &
                             decomp % mpiComm, &
                             requests(msgCount),iError)
  
            ENDDO

          ENDIF ! MPI or not

        ENDIF ! If not physical boundary

      ENDDO
    ENDDO

    IF (PRESENT(decomp) .AND. msgCount > 0) THEN
      CALL MPI_WaitAll(msgCount, &
                       requests(1:msgCount), &
                       stats(1:MPI_STATUS_SIZE,1:msgCount), &
                       iError)
    ENDIF

    DO e1 = 1,myMesh % nElem
      DO s1 = 1,6

        s2 = myMesh % sideInfo % hostData(4,s1,e1)/10
        bcid = myMesh % sideInfo % hostData(5,s1,e1)
        nloc1(1:4) = nid1(1:4,s1,e1)
        nloc2(1:4) = nid2(1:4,s1,e1)

        IF (bcid == 0) THEN
          nShifts = 0
          theyMatch = .FALSE.

          DO i = 1, 4

            theyMatch = CompareArray( nloc1, nloc2, 4 )

            IF( theyMatch )THEN
              EXIT
            ELSE
              nShifts = nShifts + 1
              CALL ForwardShift( nloc1, 4 )
            ENDIF

          ENDDO

          myMesh % sideInfo % hostData(4,s1,e1) = 10*s2+nShifts

        ENDIF

      ENDDO
    ENDDO

    DEALLOCATE(requests)
    DEALLOCATE(stats)

  END SUBROUTINE RecalculateFlip_Mesh3D

  SUBROUTINE Read_HOPr_Mesh3D(myMesh,meshFile,decomp)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(out) :: myMesh
    CHARACTER(*),INTENT(in) :: meshFile
    TYPE(MPILayer),INTENT(inout) :: decomp
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: offset(1:2),gOffset(1)
    INTEGER :: nGlobalElem
    INTEGER :: firstElem
    INTEGER :: firstNode
    INTEGER :: firstSide
    INTEGER :: nLocalElems
    INTEGER :: nLocalNodes
    INTEGER :: nLocalSides
    INTEGER :: nUniqueSides
    INTEGER :: nGeo,nBCs
    INTEGER :: eid, lsid, iSide
    INTEGER :: i, j, k, nid
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfInt32_r2) :: hopr_sideInfo
    TYPE(hfReal_r2) :: hopr_nodeCoords
    TYPE(hfInt32_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: bcType


    IF ( decomp % mpiEnabled )THEN
      CALL Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId,decomp % mpiComm)
    ELSE
      CALL Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId)
    ENDIF

    CALL ReadAttribute_HDF5(fileId,'nElems',nGlobalElem)
    CALL ReadAttribute_HDF5(fileId,'Ngeo',nGeo)
    CALL ReadAttribute_HDF5(fileId,'nBCs',nBCs)
    CALL ReadAttribute_HDF5(fileId,'nUniqueSides',nUniqueSides)

    ! Read BCType
    CALL bcType % Alloc(loBound=(/1,1/), &
                        upBound=(/4,nBCs/))
    IF ( decomp % mpiEnabled )THEN
      offset(:) = 0
      CALL ReadArray_HDF5(fileId,'BCType',bcType,offset)
    ELSE
      CALL ReadArray_HDF5(fileId,'BCType',bcType)
    ENDIF

    ! Read local subarray of ElemInfo
    CALL decomp % GenerateDecomposition(nGlobalElem,nUniqueSides)

    firstElem = decomp % offsetElem % hostData(decomp % rankId) + 1
    nLocalElems = decomp % offsetElem % hostData(decomp % rankId + 1) - &
                  decomp % offsetElem % hostData(decomp % rankId)

    ! Allocate Space for hopr_elemInfo!
    CALL hopr_elemInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/6,nLocalElems/))

    IF ( decomp % mpiEnabled )THEN
      offset = (/0,firstElem - 1/)
      CALL ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo,offset)
    ELSE
      CALL ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo)
    ENDIF

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode = hopr_elemInfo % hostData(5,1) + 1
    nLocalNodes = hopr_elemInfo % hostData(6,nLocalElems) - hopr_elemInfo % hostData(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !
    CALL hopr_nodeCoords % Alloc(loBound=(/1,1/), &
                                 upBound=(/3,nLocalNodes/))

    CALL hopr_globalNodeIDs % Alloc(loBound=1, &
                                    upBound=nLocalNodes)

    IF ( decomp % mpiEnabled )THEN
      offset = (/0,firstNode - 1/)
      CALL ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords,offset)
      gOffset = (/firstNode - 1/)
      CALL ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs,gOffset)
    ELSE
      CALL ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords)
      CALL ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs)
    ENDIF

    ! Read local subarray of SideInfo
    firstSide = hopr_elemInfo % hostData(3,1) + 1
    nLocalSides = hopr_elemInfo % hostData(4,nLocalElems) - hopr_elemInfo % hostData(3,1)

    ! Allocate space for hopr_sideInfo
    CALL hopr_sideInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/5,nLocalSides/))
    IF ( decomp % mpiEnabled )THEN
      offset = (/0,firstSide - 1/)
      CALL ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo,offset)
    ELSE
      CALL ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo)
    ENDIF

    CALL Close_HDF5(fileID)
    ! ---- Done reading 3-D Mesh information ---- !
    ! Load hopr data into mesh data structure

    CALL myMesh % Init(nGeo,nLocalElems,nLocalSides,nLocalNodes,nBCs)

    ! Copy data from local arrays into myMesh
    myMesh % elemInfo % hostData = hopr_elemInfo % hostData
    myMesh % nUniqueSides = nUniqueSides
    myMesh % quadrature = UNIFORM

    ! Grab the node coordinates
    DO eid = 1, myMesh % nElem
      DO k = 0, nGeo
        DO j = 0, nGeo
          DO i = 0, nGeo
            nid = i+1 + (nGeo+1)*(j + (nGeo+1)*(k + (nGeo+1)*(eid-1)))
            myMesh % nodeCoords % hostData(1:3,i,j,k,eid) = hopr_nodeCoords % hostData(1:3,nid)
            myMesh % globalNodeIDs % hostData(i,j,k,eid) = hopr_globalNodeIDs % hostData(nid)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    iSide = 0 
    DO eid = 1,myMesh % nElem
      DO lsid = 1,6
        iSide = iSide + 1
        myMesh % sideInfo % hostData(1:5,lsid,eid) = hopr_sideInfo % hostData(1:5,iSide)
      ENDDO
    ENDDO

    CALL myMesh % RecalculateFlip()

    CALL myMesh % UpdateDevice()

    CALL hopr_elemInfo % Free()
    CALL hopr_nodeCoords % Free()
    CALL hopr_globalNodeIDs % Free()
    CALL hopr_sideInfo % Free()

  END SUBROUTINE Read_HOPr_Mesh3D

  SUBROUTINE Write_Mesh3D(myMesh,meshFile)
    ! Writes mesh output in HOPR format (serial only)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh
    CHARACTER(*),INTENT(in) :: meshFile
    ! Local
    INTEGER(HID_T) :: fileId

    CALL Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId)

    CALL WriteAttribute_HDF5(fileId,'nElems',myMesh % nElem)
    CALL WriteAttribute_HDF5(fileId,'Ngeo',myMesh % nGeo)
    CALL WriteAttribute_HDF5(fileId,'nBCs',myMesh % nBCs)

    CALL WriteArray_HDF5(fileId,'BCType',myMesh % bcType)
    CALL WriteArray_HDF5(fileId,'ElemInfo',myMesh % elemInfo)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    CALL WriteArray_HDF5(fileId,'NodeCoords',myMesh % nodeCoords)
    CALL WriteArray_HDF5(fileId,'GlobalNodeIDs',myMesh % globalNodeIDs)

    ! Read local subarray of SideInfo
    CALL WriteArray_HDF5(fileId,'SideInfo',myMesh % sideInfo)

    CALL Close_HDF5(fileID)

  END SUBROUTINE Write_Mesh3D

!  FUNCTION TransfiniteInterpolation_2D( interp, bCurves, iEl, a, b ) RESULT( P )
!    ! TransfiniteInterpolation
!    !  Takes in the six surfaces (south, east, north, west, bottom, top) and evaluates the
!    !  bidirectional mapping at xi^1 = a, xi^2 = b, xi^3 = c. The boundary of the computational
!    !  coordinate system is assumed to be at +/- 1 in each direction.
!    !
!    ! =============================================================================================== !
!    ! DECLARATIONS
!    IMPLICIT NONE
!    REAL(prec) :: bCurves(1:2,1:4)
!    INTEGER :: i,j
!    REAL(prec) :: P(1:2)
!    ! LOCAL
!    REAL(prec)  :: P1(1:2), P2(1:2), P2(1:2)
!    REAL(prec)  :: sSurf(1:2), nSurf(1:2), eSurf(1:2), wSurf(1:2)
!    REAL(prec)  :: l1(1:2), l2(1:2), l2(1:2)
!    REAL(prec)  :: ref(1:2)
!    INTEGER     :: i, j, iSurf
!
!    ref = (/ -1.0_prec, 1.0_prec /)
!
!    ! Transfinite interpolation with linear blending USEs linear lagrange interpolating polynomials
!    ! to blend the bounding surfaces.
!    ! The linear blending weights in the first computational direction are calculated.
!
!    l1 = LinearBlend(a)
!    l2 = LinearBlend(b)
!    l3 = LinearBlend(c)
!
!    ! The bounding surfaces need to be evaluated at the provided computational coordinates
!
!    wSurf = bCurves(1:2,4)
!    eSurf = bCurves(1:2,2)
!    sSurf = bCurves(1:2,1)
!    nSurf = bCurves(1:2,3)
!
!    ! P1 CONTAINS the interpolation in the first computational coordinate
!    ! The first computational coordinate is assumed to vary between the (computational) east and
!    ! west boundaries.
!
!    P1 = l1(1)*wSurf + l1(2)*eSurf
!
!    ! P2 CONTAINS the interpolation in the second computational coordinate
!    ! The second computational coordinate is assumed to vary between the (computational) south and
!    ! north boundaries.
!
!    P2 = l2(1)*sSurf + l2(2)*nSurf
!
!    ! P3 CONTAINS the interpolation in the first computational coordinate
!    ! The first computational coordinate is assumed to vary between the (computational) bottom and
!    ! top boundaries.
!
!    P3 = l3(1)*bSurf + l3(2)*tSurf
!
!    DO i = 1, 2
!
!      ! Now we need to compute the tensor product of the first and second computational direction
!      ! interpolants and subtract from P1.
!
!      wSurf = boundingsurfaces % Evaluate( (/ref(i), c/), west + (iEl-1)*6 )
!      eSurf = boundingsurfaces % Evaluate( (/ref(i), c/), east + (iEl-1)*6 )
!      P1 = P1 - l2(i)*( wSurf*l1(1) + eSurf*l1(2) )
!
!      ! Now we need to compute the tensor product of the first and third computational direction
!      ! interpolants and subtract from P1.
!
!      wSurf = boundingsurfaces % Evaluate( (/b, ref(i)/), west + (iEl-1)*6 )
!      eSurf = boundingsurfaces % Evaluate( (/b, ref(i)/), east + (iEl-1)*6 )
!
!      P1 = P1 - l3(i)*( wSurf*l1(1) + eSurf*l1(2) )
!
!      ! Now we need to compute the tensor product of the second and third computational direction
!      ! interpolants and subtract from P2.
!
!      sSurf = boundingsurfaces % Evaluate( (/a, ref(i)/), south + (iEl-1)*6 )
!      nSurf = boundingsurfaces % Evaluate( (/a, ref(i)/), north + (iEl-1)*6 )
!
!      P2 = P2 - l3(i)*( sSurf*l2(1) + nSurf*l2(2) )
!
!    ENDDO
!
!    ! Next, the compounded tensor product is computed and added to P3.
!    DO j = 1,2
!      DO i = 1,2
!
!        wSurf = boundingsurfaces % Evaluate( (/ref(i), ref(j)/), west + (iEl-1)*6 )
!        eSurf = boundingsurfaces % Evaluate( (/ref(i), ref(j)/), east + (iEl-1)*6 )
!        P3 = P3 + l2(i)*l3(j)*( wSurf*l1(1) + eSurf*l1(2) )
!
!      ENDDO
!    ENDDO
!
!    !Finally, the sum the interpolants is computed to yield the computational coordinate
!    P = P1 + P2 + P3
!
!  END FUNCTION TransfiniteInterpolation_2D

!  FUNCTION TransfiniteInterpolation( boundingSurfaces, iEl, a, b, c ) RESULT( P )
!    ! TransfiniteInterpolation
!    !  Takes in the six surfaces (south, east, north, west, bottom, top) and evaluates the
!    !  bidirectional mapping at xi^1 = a, xi^2 = b, xi^3 = c. The boundary of the computational
!    !  coordinate system is assumed to be at +/- 1 in each direction.
!    !
!    ! =============================================================================================== !
!    ! DECLARATIONS
!    IMPLICIT NONE
!    TYPE( Surfaces )  :: boundingSurfaces
!    INTEGER           :: iEl
!    REAL(prec)       :: a, b, c
!    REAL(prec)       :: P(1:3)
!    ! LOCAL
!    REAL(prec)  :: P1(1:3), P2(1:3), P3(1:3)
!    REAL(prec)  :: sSurf(1:3), nSurf(1:3), eSurf(1:3), wSurf(1:3), bSurf(1:3), tSurf(1:3)
!    REAL(prec)  :: l1(1:2), l2(1:2), l3(1:2)
!    REAL(prec)  :: ref(1:2)
!    INTEGER     :: i, j, iSurf
!
!    ref = (/ -1.0_prec, 1.0_prec /)
!
!    ! Transfinite interpolation with linear blending USEs linear lagrange interpolating polynomials
!    ! to blend the bounding surfaces.
!    ! The linear blending weights in the first computational direction are calculated.
!
!    l1 = LinearBlend( a )
!    l2 = LinearBlend( b )
!    l3 = LinearBlend( c )
!
!    ! The bounding surfaces need to be evaluated at the provided computational coordinates
!
!    wSurf = boundingSurfaces % Evaluate( (/b, c/), west + (iEl-1)*6 )   ! west
!    eSurf = boundingSurfaces % Evaluate( (/b, c/), east + (iEl-1)*6 )   ! east
!    sSurf = boundingSurfaces % Evaluate( (/a, c/), south + (iEl-1)*6 )  ! south
!    nSurf = boundingSurfaces % Evaluate( (/a, c/), north + (iEl-1)*6 )  ! north
!    bSurf = boundingSurfaces % Evaluate( (/a, b/), bottom + (iEl-1)*6 ) ! bottom
!    tSurf = boundingSurfaces % Evaluate( (/a, b/), top + (iEl-1)*6 )    ! top
!
!    ! P1 CONTAINS the interpolation in the first computational coordinate
!    ! The first computational coordinate is assumed to vary between the (computational) east and
!    ! west boundaries.
!
!    P1 = l1(1)*wSurf + l1(2)*eSurf
!
!    ! P2 CONTAINS the interpolation in the second computational coordinate
!    ! The second computational coordinate is assumed to vary between the (computational) south and
!    ! north boundaries.
!
!    P2 = l2(1)*sSurf + l2(2)*nSurf
!
!    ! P3 CONTAINS the interpolation in the first computational coordinate
!    ! The first computational coordinate is assumed to vary between the (computational) bottom and
!    ! top boundaries.
!
!    P3 = l3(1)*bSurf + l3(2)*tSurf
!
!    DO i = 1, 2
!
!      ! Now we need to compute the tensor product of the first and second computational direction
!      ! interpolants and subtract from P1.
!
!      wSurf = boundingsurfaces % Evaluate( (/ref(i), c/), west + (iEl-1)*6 )
!      eSurf = boundingsurfaces % Evaluate( (/ref(i), c/), east + (iEl-1)*6 )
!      P1 = P1 - l2(i)*( wSurf*l1(1) + eSurf*l1(2) )
!
!      ! Now we need to compute the tensor product of the first and third computational direction
!      ! interpolants and subtract from P1.
!
!      wSurf = boundingsurfaces % Evaluate( (/b, ref(i)/), west + (iEl-1)*6 )
!      eSurf = boundingsurfaces % Evaluate( (/b, ref(i)/), east + (iEl-1)*6 )
!
!      P1 = P1 - l3(i)*( wSurf*l1(1) + eSurf*l1(2) )
!
!      ! Now we need to compute the tensor product of the second and third computational direction
!      ! interpolants and subtract from P2.
!
!      sSurf = boundingsurfaces % Evaluate( (/a, ref(i)/), south + (iEl-1)*6 )
!      nSurf = boundingsurfaces % Evaluate( (/a, ref(i)/), north + (iEl-1)*6 )
!
!      P2 = P2 - l3(i)*( sSurf*l2(1) + nSurf*l2(2) )
!
!    ENDDO
!
!    ! Next, the compounded tensor product is computed and added to P3.
!    DO j = 1,2
!      DO i = 1,2
!
!        wSurf = boundingsurfaces % Evaluate( (/ref(i), ref(j)/), west + (iEl-1)*6 )
!        eSurf = boundingsurfaces % Evaluate( (/ref(i), ref(j)/), east + (iEl-1)*6 )
!        P3 = P3 + l2(i)*l3(j)*( wSurf*l1(1) + eSurf*l1(2) )
!
!      ENDDO
!    ENDDO
!
!    !Finally, the sum the interpolants is computed to yield the computational coordinate
!    P = P1 + P2 + P3
!
!  END FUNCTION TransfiniteInterpolation
  FUNCTION Unidirectional(valLeft,valRight,a) RESULT(P)
    !
    ! =============================================================================================== !
    ! DECLARATIONS
    IMPLICIT NONE
    REAL(prec) :: valLeft(1:3),valRight(1:3)
    REAL(prec) :: a
    REAL(prec) :: P(1:3)

    P = 0.5_prec*((1.0_prec - a)*valLeft + (1.0_prec + a)*valRight)

  END FUNCTION Unidirectional
  FUNCTION LinearBlend(a) RESULT(weights)

    IMPLICIT NONE
    REAL(prec) :: a
    REAL(prec) :: weights(1:2)

    weights(1) = 0.5_prec*(1.0_prec - a)
    weights(2) = 0.5_prec*(1.0_prec + a)

  END FUNCTION LinearBlend

  SUBROUTINE Init_MPILayer(this,enableMPI)
#undef __FUNC__
#define __FUNC__ "Init_MPILayer"
    IMPLICIT NONE
    CLASS(MPILayer),INTENT(out) :: this
    LOGICAL,INTENT(in) :: enableMPI
    ! Local
    INTEGER       :: ierror
    CHARACTER(50) :: msg
    INTEGER :: nGPU, gpuID
    CHARACTER(2) :: msg2

    this % mpiComm = 0
    this % mpiPrec = prec
    this % rankId = 0
    this % nRanks = 1
    this % nElem = 0
    this % mpiEnabled = enableMPI

    IF (enableMPI) THEN
      this % mpiComm = MPI_COMM_WORLD
      CALL MPI_INIT(ierror)
      CALL MPI_COMM_RANK(this % mpiComm,this % rankId,ierror)
      CALL MPI_COMM_SIZE(this % mpiComm,this % nRanks,ierror)
    END IF

    IF (prec == real32) THEN
      this % mpiPrec = MPI_FLOAT
    ELSE
      this % mpiPrec = MPI_DOUBLE
    END IF

    CALL this % offSetElem % Alloc(0,this % nRanks)

    WRITE (msg,'(I5)') this % rankId
    msg = "Greetings from rank "//TRIM(msg)//"."
    INFO(TRIM(msg))

    IF( GPUAvailable() )THEN
      ! Get the number of GPUs per node
      CALL hipCheck(hipGetDeviceCount(nGPU))

      ! Assume that we have the 1 GPU per rank
      ! implying that nMPIRanksPerNode = nGPU
      ! Assume that all nodes have the same number of GPUs per node
      gpuID = MOD(this % rankId, nGPU)

      CALL hipCheck(hipSetDevice(gpuID))
      WRITE (msg,'(I5)') this % rankId
      WRITE (msg2,'(I2)') gpuID
      msg = "Rank "//TRIM(msg)//": Setting device to GPU"//TRIM(msg2)
      INFO(TRIM(msg))

    ENDIF

  END SUBROUTINE Init_MPILayer

  SUBROUTINE Free_MPILayer(this)
    IMPLICIT NONE
    CLASS(MPILayer),INTENT(inout) :: this

    IF (ASSOCIATED(this % offSetElem % hostData)) THEN
      CALL this % offSetElem % Free()
    ENDIF
    IF (ASSOCIATED(this % elemToRank % hostData)) THEN
      CALL this % elemToRank % Free()
    ENDIF

    DEALLOCATE( this % requests )
    DEALLOCATE( this % stats )

  END SUBROUTINE Free_MPILayer

  SUBROUTINE Finalize_MPILayer(this)
#undef __FUNC__
#define __FUNC__ "Finalize_MPILayer"
    IMPLICIT NONE
    CLASS(MPILayer),INTENT(inout) :: this
    ! Local
    INTEGER       :: ierror
    CHARACTER(30) :: msg

    IF (this % mpiEnabled) THEN
      WRITE (msg,'(I5)') this % rankId
      msg = "Goodbye from rank "//TRIM(msg)//"."
      INFO(TRIM(msg))
      CALL MPI_FINALIZE(ierror)
    ENDIF

  END SUBROUTINE Finalize_MPILayer

  SUBROUTINE GenerateDecomposition_MPILayer(this,nGlobalElem,maxMsg)
#undef __FUNC__
#define __FUNC__ "GenerateDecomposition_MPILayer"
    IMPLICIT NONE
    CLASS(MPILayer),INTENT(inout) :: this
    INTEGER,INTENT(in) :: nGlobalElem
    INTEGER,INTENT(in) :: maxMsg
    ! Local
    INTEGER :: maxMsgLoc
    CHARACTER(50) :: msg
    CHARACTER(5) :: msg2

    CALL this % setElemToRank(nGlobalElem)
    CALL this % SetMaxMsg(maxMsg)

    WRITE (msg,'(I5)') this % rankId
    WRITE (msg2,'(I5)') this % offSetElem % hostData(this % rankId + 1)-&
                        this % offSetElem % hostData(this % rankId)
    msg = "Rank "//TRIM(msg)//": nElem = "//TRIM(msg2)
    INFO(TRIM(msg))

  END SUBROUTINE GenerateDecomposition_MPILayer

  SUBROUTINE SetMaxMsg(this,maxMsg)
    IMPLICIT NONE
    CLASS(MPILayer),INTENT(inout) :: this
    INTEGER,INTENT(in) :: maxMsg

    IF (ALLOCATED(this % requests) ) DEALLOCATE( this % requests )
    IF (ALLOCATED(this % stats) ) DEALLOCATE( this % stats )

    ALLOCATE( this % requests(1:maxMsg) )
    ALLOCATE( this % stats(MPI_STATUS_SIZE,1:maxMsg) )
    this % maxMsg = maxMsg

  END SUBROUTINE SetMaxMsg

  SUBROUTINE SetElemToRank(this,nElem)
    IMPLICIT NONE
    CLASS(MPILayer),INTENT(inout) :: this
    INTEGER,INTENT(in) :: nElem
    ! Local
    INTEGER :: iel

    this % nElem = nElem

    CALL this % elemToRank % Alloc(1,nElem)

    CALL DomainDecomp(nElem, &
                      this % nRanks, &
                      this % offSetElem % hostData)

    DO iel = 1,nElem
      CALL ElemToRank(this % nRanks, &
                      this % offSetElem % hostData, &
                      iel, &
                      this % elemToRank % hostData(iel))
    END DO

    CALL this % offSetElem % UpdateDevice()
    CALL this % elemToRank % UpdateDevice()

  END SUBROUTINE SetElemToRank

  SUBROUTINE DomainDecomp(nElems,nDomains,offSetElem)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 4
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nElems
    INTEGER,INTENT(in) :: nDomains
    INTEGER,INTENT(out) :: offsetElem(0:nDomains)
    ! Local
    INTEGER :: nLocalElems
    INTEGER :: remainElems
    INTEGER :: iDom

    nLocalElems = nElems/nDomains
    remainElems = nElems - nLocalElems*nDomains
    DO iDom = 0,nDomains - 1
      offSetElem(iDom) = iDom*nLocalElems + MIN(iDom,remainElems)
    END DO
    offSetElem(nDomains) = nElems

  END SUBROUTINE DomainDecomp

  SUBROUTINE ElemToRank(nDomains,offsetElem,elemID,domain)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 7
    !   "Find domain containing element index"
    !
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nDomains
    INTEGER,INTENT(in) :: offsetElem(0:nDomains)
    INTEGER,INTENT(in) :: elemID
    INTEGER,INTENT(out) :: domain
    ! Local
    INTEGER :: maxSteps
    INTEGER :: low,up,mid
    INTEGER :: i

    domain = 0
    maxSteps = INT(LOG10(REAL(nDomains))/LOG10(2.0)) + 1
    low = 0
    up = nDomains - 1

    IF (offsetElem(low) < elemID .AND. elemID <= offsetElem(low + 1)) THEN
      domain = low
    ELSEIF (offsetElem(up) < elemID .AND. elemID <= offsetElem(up + 1)) THEN
      domain = up
    ELSE
      DO i = 1,maxSteps
        mid = (up - low)/2 + low
        IF (offsetElem(mid) < elemID .AND. elemID <= offsetElem(mid + 1)) THEN
          domain = mid
          RETURN
        ELSEIF (elemID > offsetElem(mid + 1)) THEN
          low = mid + 1
        ELSE
          up = mid
        END IF
      END DO
    END IF

  END SUBROUTINE ElemToRank

  SUBROUTINE FinalizeMPIExchangeAsync(mpiHandler)
    CLASS(MPILayer),INTENT(inout) :: mpiHandler
    ! Local
    INTEGER :: ierror

    IF( mpiHandler % mpiEnabled )THEN
    CALL MPI_WaitAll(mpiHandler % msgCount, &
                     mpiHandler % requests(1:mpiHandler % msgCount), &
                     mpiHandler % stats(1:MPI_STATUS_SIZE,1:mpiHandler % msgCount), &
                     iError)
    ENDIF

  END SUBROUTINE FinalizeMPIExchangeAsync

  SUBROUTINE GlobalReduce_RealScalar(mpiHandler, sendBuf, recvBuf)
    CLASS(MPILayer), INTENT(in) :: mpiHandler
    REAL(prec), INTENT(in) :: sendBuf
    REAL(prec), INTENT(out) :: recvBuf
    ! Local
    INTEGER :: iError

      IF (mpiHandler % mpiEnabled) THEN
        CALL MPI_ALLREDUCE(sendBuf, &
                           recvBuf, &
                           1, &
                           mpiHandler % mpiPrec, &
                           MPI_SUM, &
                           mpiHandler % mpiComm, &
                           iError)
      ELSE
        recvBuf = sendBuf
      ENDIF

  END SUBROUTINE GlobalReduce_RealScalar

END MODULE SELF_Mesh
