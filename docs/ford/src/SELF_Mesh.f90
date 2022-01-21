!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_Mesh

  USE SELF_Constants
  USE SELF_HashTable
  USE SELF_Lagrange
  USE SELF_MPI
  USE SELF_Data
  USE SELF_SupportRoutines
  USE SELF_HDF5
  ! External Libs !
  USE HDF5
  USE MPI

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

  TYPE,PUBLIC :: Mesh1D
    INTEGER :: nGeo
    INTEGER :: nGlobalElem
    INTEGER :: nElem
    INTEGER :: nNodes
    INTEGER :: nCornerNodes
    INTEGER :: nUniqueNodes
    INTEGER :: nBCs
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfReal_r1) :: hopr_nodeCoords
    TYPE(hfInt32_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS
    PROCEDURE,PUBLIC :: Init => Init_Mesh1D
    PROCEDURE,PUBLIC :: Free => Free_Mesh1D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Mesh1D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh1D
    PROCEDURE,PUBLIC :: UniformBlockMesh => UniformBlockMesh_Mesh1D

    GENERIC,PUBLIC  :: Read_HOPr => Read_HOPr_Mesh1D_serial,Read_HOPr_Mesh1D_parallel
    PROCEDURE,PRIVATE :: Read_HOPr_Mesh1D_serial,Read_HOPr_Mesh1D_parallel
    PROCEDURE,PUBLIC  :: Write_HOPr => Write_HOPr_Mesh1D

  END TYPE Mesh1D

  ! Mesh format is set up as the HOPr format
  ! See https://hopr-project.org/externals/MeshFormat.pdf

  TYPE,PUBLIC :: Mesh2D
    INTEGER :: nGeo
    INTEGER :: nGlobalElem
    INTEGER :: nElem
    INTEGER :: nNodes
    INTEGER :: nSides
    INTEGER :: nCornerNodes
    INTEGER :: nUniqueNodes
    INTEGER :: nUniqueSides
    INTEGER :: nBCs
    TYPE(hfInt32_r3) :: self_sideInfo
    TYPE(hfReal_r3) :: self_nodeCoords
    TYPE(hfReal_r2) :: hohq_cornerNodes
    TYPE(hfInt32_r2) :: hohq_elemInfo
    TYPE(hfInt32_r2) :: hohq_sideInfo
    TYPE(hfReal_r4) :: hohq_sideCurves
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfInt32_r2) :: hopr_sideInfo
    TYPE(hfReal_r2) :: hopr_nodeCoords
    TYPE(hfInt32_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r1) :: hopr_CGNSCornerMap
    TYPE(hfInt32_r2) :: hopr_CGNSSideMap
    TYPE(hfInt32_r2) :: hopr_curveNodeMap
    TYPE(hfInt32_r2) :: hopr_curveNodeMapInv
    TYPE(hfInt32_r2) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS
    PROCEDURE,PUBLIC :: Init => Init_Mesh2D
    PROCEDURE,PUBLIC :: Free => Free_Mesh2D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Mesh2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh2D
    PROCEDURE,PUBLIC :: UniformBlockMesh => UniformBlockMesh_Mesh2D

    GENERIC,PUBLIC :: Load => Load_Mesh2D_serial,Load_Mesh2D_parallel
    PROCEDURE,PRIVATE :: Load_Mesh2D_serial,Load_Mesh2D_parallel
    GENERIC,PUBLIC :: Read_HOPr => Read_HOPr_Mesh2D_serial,Read_HOPr_Mesh2D_parallel
    PROCEDURE,PRIVATE :: Read_HOPr_Mesh2D_serial,Read_HOPr_Mesh2D_parallel
    PROCEDURE,PUBLIC :: Write_HOPr => Write_HOPr_Mesh2D

    PROCEDURE,PUBLIC :: Read_ISMv2 => Read_ISMv2_Mesh2D

    PROCEDURE,PRIVATE :: GenerateConnectivity => GenerateConnectivity_Mesh2D

  END TYPE Mesh2D

  TYPE,PUBLIC :: Mesh3D
    INTEGER :: nGeo
    INTEGER :: nElem
    INTEGER :: nGlobalElem
    INTEGER :: nNodes
    INTEGER :: nSides
    INTEGER :: nCornerNodes
    INTEGER :: nUniqueNodes
    INTEGER :: nUniqueSides
    INTEGER :: nBCs
    TYPE(hfInt32_r3) :: self_sideInfo
    TYPE(hfReal_r3) :: self_nodeCoords
    TYPE(hfReal_r2) :: hohq_cornerNodes
    TYPE(hfInt32_r2) :: hohq_elemInfo
    TYPE(hfInt32_r2) :: hohq_sideInfo
    TYPE(hfReal_r5) :: hohq_sideSurfaces
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfInt32_r2) :: hopr_sideInfo
    TYPE(hfReal_r2) :: hopr_nodeCoords
    TYPE(hfInt32_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r1) :: hopr_CGNSCornerMap
    TYPE(hfInt32_r2) :: hopr_CGNSSideMap
    TYPE(hfInt32_r2) :: self_sideMap
    TYPE(hfInt32_r2) :: hopr_curveNodeMap
    TYPE(hfInt32_r3) :: hopr_curveNodeMapInv
    TYPE(hfInt32_r2) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Mesh3D
    PROCEDURE,PUBLIC :: Free => Free_Mesh3D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Mesh3D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh3D
    PROCEDURE,PUBLIC :: UniformBlockMesh => UniformBlockMesh_Mesh3D

    GENERIC,PUBLIC :: Load => Load_Mesh3D_serial,Load_Mesh3D_parallel
    PROCEDURE,PRIVATE :: Load_Mesh3D_serial,Load_Mesh3D_parallel
    GENERIC,PUBLIC :: Read_HOPr => Read_HOPr_Mesh3D_serial,Read_HOPr_Mesh3D_parallel
    PROCEDURE,PRIVATE :: Read_HOPr_Mesh3D_serial,Read_HOPr_Mesh3D_parallel
    PROCEDURE,PUBLIC :: Write_HOPr => Write_HOPr_Mesh3D

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

    CALL myMesh % hopr_elemInfo % Alloc(loBound=(/1,1/), &
                                        upBound=(/4,nElem/))

    CALL myMesh % hopr_nodeCoords % Alloc(loBound=1, &
                                          upBound=nNodes)

    CALL myMesh % hopr_globalNodeIDs % Alloc(loBound=1, &
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

    CALL myMesh % hopr_elemInfo % Free()
    CALL myMesh % hopr_nodeCoords % Free()
    CALL myMesh % hopr_globalNodeIDs % Free()
    CALL myMesh % BCType % Free()

    DEALLOCATE (myMesh % BCNames)

  END SUBROUTINE Free_Mesh1D

  SUBROUTINE UpdateHost_Mesh1D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(inout) :: myMesh

    CALL myMesh % hopr_elemInfo % UpdateHost()
    CALL myMesh % hopr_nodeCoords % UpdateHost()
    CALL myMesh % hopr_globalNodeIDs % UpdateHost()
    CALL myMesh % BCType % UpdateHost()

  END SUBROUTINE UpdateHost_Mesh1D

  SUBROUTINE UpdateDevice_Mesh1D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(inout) :: myMesh

    CALL myMesh % hopr_elemInfo % UpdateDevice()
    CALL myMesh % hopr_nodeCoords % UpdateDevice()
    CALL myMesh % hopr_globalNodeIDs % UpdateDevice()
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
    TYPE(Scalar1D) :: xLinear
    TYPE(Scalar1D) :: xGeo

    nNodes = nElem*(nGeo + 1)
    CALL myMesh % Init(nGeo,nElem,nNodes,2)

    ! Set the hopr_nodeCoords
    xU = UniformPoints(x(1),x(2),1,nElem + 1)

    ! Create a linear interpolant to interpolate to nGeo grid
    CALL xLinear % Init(1,GAUSS_LOBATTO, &
                        nGeo,GAUSS_LOBATTO, &
                        1,nElem)

    CALL xGeo % Init(nGeo,GAUSS_LOBATTO, &
                     nGeo,GAUSS_LOBATTO, &
                     1,nElem)
    DO iel = 1,nElem
      xLinear % interior % hostData(0:1,1,iel) = xU(iel:iel + 1)
    END DO

    CALL xLinear % GridInterp(xGeo,.FALSE.)

    ! Set the element information
    nid = 1
    DO iel = 1,nElem
      myMesh % hopr_elemInfo % hostData(1,iel) = selfLineLinear ! Element Type
      myMesh % hopr_elemInfo % hostData(2,iel) = 1 ! Element Zone
      myMesh % hopr_elemInfo % hostData(3,iel) = nid ! Node Index Start
      DO i = 0,nGeo
        myMesh % hopr_nodeCoords % hostData(nid) = xGeo % interior % hostData(i,1,iel)
        nid = nid + 1
      END DO
      myMesh % hopr_elemInfo % hostData(4,iel) = nid - 1 ! Node Index End
    END DO

    CALL myMesh % UpdateDevice()

    CALL xLinear % Free()
    CALL xGeo % Free()

  END SUBROUTINE UniformBlockMesh_Mesh1D

  SUBROUTINE Read_HOPr_Mesh1D_serial(myMesh,meshFile)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    ! Adapted for 1D Mesh : Note that HOPR does not have 1D mesh output.
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(out) :: myMesh
    CHARACTER(*),INTENT(in) :: meshFile
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER :: nGlobalElem
    INTEGER :: nLocalNodes
    INTEGER :: nGeo,nBCs
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfReal_r1) :: hopr_nodeCoords
    TYPE(hfInt32_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: bcType

    CALL Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId)

    CALL ReadAttribute_HDF5(fileId,'nElems',nGlobalElem)
    CALL ReadAttribute_HDF5(fileId,'Ngeo',nGeo)
    CALL ReadAttribute_HDF5(fileId,'nBCs',nBCs)

    ! Read BCType
    CALL bcType % Alloc(loBound=(/1,1/), &
                        upBound=(/4,nBCs/))
    CALL ReadArray_HDF5(fileId,'BCType',bcType)

    ! Read local subarray of ElemInfo
    ! Allocate Space for hopr_elemInfo!
    CALL hopr_elemInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/6,nGlobalElem/))

    CALL ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    nLocalNodes = hopr_elemInfo % hostData(6,nGlobalElem) - hopr_elemInfo % hostData(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !
    CALL hopr_nodeCoords % Alloc(loBound=1, &
                                 upBound=nLocalNodes)

    CALL hopr_globalNodeIDs % Alloc(loBound=1, &
                                    upBound=nLocalNodes)

    CALL ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords)
    CALL ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs)

    CALL Close_HDF5(fileID)

    CALL myMesh % Init(nGeo,nGlobalElem,nLocalNodes,nBCs)

    ! Copy data from local arrays into myMesh
    myMesh % hopr_elemInfo = hopr_elemInfo
    myMesh % hopr_nodeCoords = hopr_nodeCoords
    myMesh % hopr_globalNodeIDs = hopr_globalNodeIDs

    CALL myMesh % UpdateDevice()

    CALL hopr_elemInfo % Free()
    CALL hopr_nodeCoords % Free()
    CALL hopr_globalNodeIDs % Free()

  END SUBROUTINE Read_HOPr_Mesh1D_serial

  SUBROUTINE Read_HOPr_Mesh1D_parallel(myMesh,meshFile,decomp)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    ! Adapted for 1D Mesh : Note that HOPR does not have 1D mesh output.
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(out) :: myMesh
    CHARACTER(*),INTENT(in) :: meshFile
    TYPE(MPILayer),INTENT(inout) :: decomp
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: offset(1:2),gOffset(1)
    INTEGER :: nGlobalElem
    INTEGER :: firstElem,nLocalElems
    INTEGER :: firstNode,nLocalNodes
    INTEGER :: nGeo,nBCs
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfReal_r1) :: hopr_nodeCoords
    TYPE(hfInt32_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: bcType

    CALL Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId,decomp % mpiComm)

    CALL ReadAttribute_HDF5(fileId,'nElems',nGlobalElem)
    CALL ReadAttribute_HDF5(fileId,'Ngeo',nGeo)
    CALL ReadAttribute_HDF5(fileId,'nBCs',nBCs)

    CALL decomp % SetElemToRank(nGlobalElem)

    ! Read BCType
    CALL bcType % Alloc(loBound=(/1,1/), &
                        upBound=(/4,nBCs/))
    offset(:) = 0
    CALL ReadArray_HDF5(fileId,'BCType',bcType,offset)

    ! Read local subarray of ElemInfo
    firstElem = decomp % offsetElem % hostData(decomp % rankId) + 1
    nLocalElems = decomp % offsetElem % hostData(decomp % rankId + 1) - &
                  decomp % offsetElem % hostData(decomp % rankId)

    ! Allocate Space for hopr_elemInfo!
    CALL hopr_elemInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/6,nLocalElems/))

    offset = (/0,firstElem - 1/)
    CALL ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo,offset)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode = hopr_elemInfo % hostData(5,1) + 1
    nLocalNodes = hopr_elemInfo % hostData(6,nLocalElems) - hopr_elemInfo % hostData(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !
    CALL hopr_nodeCoords % Alloc(loBound=1, &
                                 upBound=nLocalNodes)

    CALL hopr_globalNodeIDs % Alloc(loBound=1, &
                                    upBound=nLocalNodes)

    gOffset = (/firstNode - 1/)
    CALL ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords,gOffset)
    CALL ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs,gOffset)

    CALL Close_HDF5(fileID)

    CALL myMesh % Init(nGeo,nLocalElems,nLocalNodes,nBCs)

    ! Copy data from local arrays into myMesh
    myMesh % hopr_elemInfo = hopr_elemInfo
    myMesh % hopr_nodeCoords = hopr_nodeCoords
    myMesh % hopr_globalNodeIDs = hopr_globalNodeIDs

    CALL myMesh % UpdateDevice()

    CALL hopr_elemInfo % Free()
    CALL hopr_nodeCoords % Free()
    CALL hopr_globalNodeIDs % Free()

  END SUBROUTINE Read_HOPr_Mesh1D_parallel

  SUBROUTINE Write_HOPr_Mesh1D(myMesh,meshFile)
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
    CALL WriteArray_HDF5(fileId,'ElemInfo',myMesh % hopr_elemInfo)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    CALL WriteArray_HDF5(fileId,'NodeCoords',myMesh % hopr_nodeCoords)
    CALL WriteArray_HDF5(fileId,'GlobalNodeIDs',myMesh % hopr_globalNodeIDs)

    CALL Close_HDF5(fileID)

  END SUBROUTINE Write_HOPr_Mesh1D

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

    CALL myMesh % hopr_elemInfo % Alloc(loBound=(/1,1/), &
                                        upBound=(/6,nElem/))

    CALL myMesh % self_sideInfo % Alloc(loBound=(/1,1,1/), &
                                        upBound=(/5,4,nElem/))

    CALL myMesh % hopr_sideInfo % Alloc(loBound=(/1,1/), &
                                        upBound=(/5,nSides/))

    CALL myMesh % hopr_nodeCoords % Alloc(loBound=(/1,1/), &
                                          upBound=(/2,nNodes/))

    CALL myMesh % self_nodeCoords % Alloc(loBound=(/1,1,1/), &
                                          upBound=(/2, (nGeo + 1)**2,nElem/))

    CALL myMesh % hopr_globalNodeIDs % Alloc(loBound=1, &
                                             upBound=nNodes)

    CALL myMesh % hopr_CGNSCornerMap % Alloc(loBound=1, &
                                             upBound=4)

    CALL myMesh % hopr_CGNSSideMap % Alloc(loBound=(/1,1/), &
                                           upBound=(/2,4/))

    CALL myMesh % hopr_curveNodeMap % Alloc(loBound=(/1,1/), &
                                            upBound=(/2, (nGeo + 1)**2/))

    CALL myMesh % hopr_curveNodeMapInv % Alloc(loBound=(/0,0/), &
                                               upBound=(/nGeo,nGeo/))

    CALL myMesh % BCType % Alloc(loBound=(/1,1/), &
                                 upBound=(/4,nBCs/))

    ALLOCATE (myMesh % BCNames(1:nBCs))

    ! Create lookup tables to assist with connectivity generation
    myMesh % hopr_CGNSCornerMap % hostData(1) = 1
    myMesh % hopr_CGNSCornerMap % hostData(2) = nGeo + 1
    myMesh % hopr_CGNSCornerMap % hostData(3) = (nGeo + 1)**2
    myMesh % hopr_CGNSCornerMap % hostData(4) = nGeo*(nGeo + 1) + 1

    l = 0
    DO j = 0,nGeo
      DO i = 0,nGeo
        l = l + 1
        myMesh % hopr_curveNodeMap % hostData(1:2,l) = (/i,j/)
        myMesh % hopr_curveNodeMapInv % hostData(i,j) = l
      END DO
    END DO

    ! Maps from local corner node id to CGNS side
    myMesh % hopr_CGNSSideMap % hostData(1:2,1) = (/1,2/)
    myMesh % hopr_CGNSSideMap % hostData(1:2,2) = (/2,3/)
    myMesh % hopr_CGNSSideMap % hostData(1:2,3) = (/4,3/)
    myMesh % hopr_CGNSSideMap % hostData(1:2,4) = (/1,4/)

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

    CALL myMesh % hopr_elemInfo % Free()
    CALL myMesh % hopr_sideInfo % Free()
    CALL myMesh % self_sideInfo % Free()
    CALL myMesh % hopr_nodeCoords % Free()
    CALL myMesh % self_nodeCoords % Free()
    CALL myMesh % hopr_CGNSCornerMap % Free()
    CALL myMesh % hopr_globalNodeIDs % Free()
    CALL myMesh % BCType % Free()

    DEALLOCATE (myMesh % BCNames)

  END SUBROUTINE Free_Mesh2D

  SUBROUTINE UpdateHost_Mesh2D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh

    CALL myMesh % hopr_elemInfo % UpdateHost()
    CALL myMesh % hopr_sideInfo % UpdateHost()
    CALL myMesh % self_sideInfo % UpdateHost()
    CALL myMesh % hopr_nodeCoords % UpdateHost()
    CALL myMesh % self_nodeCoords % UpdateHost()
    CALL myMesh % hopr_globalNodeIDs % UpdateHost()
    CALL myMesh % BCType % UpdateHost()

  END SUBROUTINE UpdateHost_Mesh2D

  SUBROUTINE UpdateDevice_Mesh2D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh

    CALL myMesh % hopr_elemInfo % UpdateDevice()
    CALL myMesh % hopr_sideInfo % UpdateDevice()
    CALL myMesh % self_sideInfo % UpdateDevice()
    CALL myMesh % hopr_nodeCoords % UpdateDevice()
    CALL myMesh % self_nodeCoords % UpdateDevice()
    CALL myMesh % hopr_globalNodeIDs % UpdateDevice()
    CALL myMesh % BCType % UpdateDevice()

  END SUBROUTINE UpdateDevice_Mesh2D

  SUBROUTINE UniformBlockMesh_Mesh2D(myMesh,nGeo,nElem,x)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(out) :: myMesh
    INTEGER,INTENT(in) :: nGeo
    INTEGER,INTENT(in) :: nElem(1:2)
    REAL(prec),INTENT(in) :: x(1:4)
    ! Local
    INTEGER :: iel,jel,nEl,elid
    INTEGER :: sid,nid,nNodes
    INTEGER :: nSides
    INTEGER :: i,j
    REAL(prec) :: xU(1:nElem(1) + 1)
    REAL(prec) :: yU(1:nElem(2) + 1)
    TYPE(Vector2D) :: xLinear
    TYPE(Vector2D) :: xGeo

    nEl = nElem(1)*nElem(2)
    nNodes = nEl*(nGeo + 1)*(nGeo + 1)
    nSides = nEl*4
    CALL myMesh % Init(nGeo,nEl,nSides,nNodes,1)
    myMesh % nUniqueSides = (nElem(1) + 1)*nElem(2) + (nElem(2) + 1)*nElem(1)

    ! Set the hopr_nodeCoords
    xU = UniformPoints(x(1),x(2),1,nElem(1) + 1)
    yU = UniformPoints(x(3),x(4),1,nElem(2) + 1)

    ! Create a linear interpolant to interpolate to nGeo grid
    CALL xLinear % Init(1,GAUSS_LOBATTO, &
                        nGeo,GAUSS_LOBATTO, &
                        1,nEl)

    CALL xGeo % Init(nGeo,GAUSS_LOBATTO, &
                     nGeo,GAUSS_LOBATTO, &
                     1,nEl)
    elid = 1
    DO jel = 1,nElem(2)
      DO iel = 1,nElem(1)
        ! x component
        xLinear % interior % hostData(1,0:1,0,1,elid) = xU(iel:iel + 1)
        xLinear % interior % hostData(1,0:1,1,1,elid) = xU(iel:iel + 1)
        ! y component
        xLinear % interior % hostData(2,0,0:1,1,elid) = yU(jel:jel + 1)
        xLinear % interior % hostData(2,1,0:1,1,elid) = yU(jel:jel + 1)
        ! Incremenent the element ID
        elid = elid + 1
      END DO
    END DO

    CALL xLinear % GridInterp(xGeo,.FALSE.)

    ! Set the element information
    nid = 1
    sid = 0
    elid = 1
    DO jel = 1,nElem(2)
      DO iel = 1,nElem(1)
        myMesh % hopr_elemInfo % hostData(1,elid) = selfQuadLinear ! Element Type
        myMesh % hopr_elemInfo % hostData(2,elid) = 1 ! Element Zone
        myMesh % hopr_elemInfo % hostData(3,elid) = sid ! Side Index Start
        sid = sid + 4
        myMesh % hopr_elemInfo % hostData(4,elid) = sid ! Side Index End
        myMesh % hopr_elemInfo % hostData(5,elid) = nid - 1 ! Node Index Start
        DO j = 0,nGeo
          DO i = 0,nGeo
            myMesh % hopr_nodeCoords % hostData(1:2,nid) = xGeo % interior % hostData(1:2,i,j,1,elid)
            nid = nid + 1
          END DO
        END DO
        myMesh % hopr_elemInfo % hostData(6,elid) = nid ! Node Index End
        elid = elid + 1
      END DO
    END DO

    CALL myMesh % GenerateConnectivity()

    CALL myMesh % UpdateDevice()

    CALL xLinear % Free()
    CALL xGeo % Free()

  END SUBROUTINE UniformBlockMesh_Mesh2D

  SUBROUTINE Load_Mesh2D_serial(myMesh,myMeshSpec)
#undef __FUNC__
#define __FUNC__ "Load_Mesh2D"
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(out) :: myMesh
    TYPE(MeshSpec),INTENT(in) :: myMeshSpec

    IF (myMeshSpec % blockMesh) THEN

      CALL myMesh % UniformBlockMesh(myMeshSpec % blockMesh_nGeo, &
                                     (/myMeshSpec % blockMesh_nElemX,myMeshSpec % blockMesh_nElemY/), &
                                     (/myMeshSpec % blockMesh_x0,myMeshSpec % blockMesh_x1, &
                                       myMeshSpec % blockMesh_y0,myMeshSpec % blockMesh_y1/))
    ELSE

      IF (myMeshSpec % fileType == SELF_MESH_ISM_V2_2D)THEN

        CALL myMesh % Read_ISMv2(myMeshSpec % filename)

      ELSEIF (myMeshSpec % fileType == SELF_MESH_HOPR_2D)THEN

        CALL myMesh % Read_HOPr(myMeshSpec % filename)

      ENDIF


    END IF

    CALL myMesh % UpdateDevice()

  END SUBROUTINE Load_Mesh2D_serial

  SUBROUTINE Load_Mesh2D_parallel(myMesh,myMeshSpec,decomp)
#undef __FUNC__
#define __FUNC__ "Load_Mesh2D"
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(out) :: myMesh
    TYPE(MeshSpec),INTENT(in) :: myMeshSpec
    TYPE(MPILayer),INTENT(inout) :: decomp

    IF (myMeshSpec % blockMesh) THEN

      IF (decomp % nRanks > 1) THEN ! Error out
        ERROR("Block Mesh only supported in serial")
        STOP ! TO DO : Safe exit
      ELSE
        CALL myMesh % UniformBlockMesh(myMeshSpec % blockMesh_nGeo, &
                                       (/myMeshSpec % blockMesh_nElemX,myMeshSpec % blockMesh_nElemY/), &
                                       (/myMeshSpec % blockMesh_x0,myMeshSpec % blockMesh_x1, &
                                         myMeshSpec % blockMesh_y0,myMeshSpec % blockMesh_y1/))
      END IF

    ELSE

      CALL myMesh % Read_HOPr(myMeshSpec % filename,decomp)

    END IF

    CALL myMesh % UpdateDevice()

  END SUBROUTINE Load_Mesh2D_parallel

  SUBROUTINE Read_ISMv2_Mesh2D(myMesh,meshFile)
#undef __FUNC__
#define __FUNC__ "Read_ISMv2_Mesh2D"
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(out) :: myMesh
    CHARACTER(*),INTENT(in) :: meshFile
    ! Local
    INTEGER :: nNodes,nSides,nElem,nGeo
    INTEGER :: lnid,nid,usid,lsid,sid,eid
    INTEGER :: fUnit
    INTEGER :: i,j
    CHARACTER(100) :: line
    CHARACTER(500) :: msg
    REAL(prec) :: x(1:3),x0(1:2),x1(1:2)
    REAL(prec) :: wSurf(1:2),eSurf(1:2),nSurf(1:2),sSurf(1:2)
    REAL(prec) :: P1(1:2),P2(1:2)
    REAL(prec) :: l1(1:2),l2(1:2)
    REAL(prec) :: se(1:2),sw(1:2),ne(1:2),nw(1:2)
    INTEGER :: bCurveFlag(1:4)
    TYPE(Lagrange) :: interp

    OPEN (UNIT=NEWUNIT(fUnit), &
          FILE=TRIM(meshFile), &
          FORM='FORMATTED', &
          STATUS='OLD', &
          ACCESS='SEQUENTIAL')

    READ (fUnit,*) line
    IF (TRIM(line) /= 'ISM-V2') THEN
      msg = 'Unrecognized file format : '//TRIM(line)
      ERROR(msg)
      STOP
    END IF

    ! Number of Nodes, Number of Edges (sides; unique!), number of elements, polynomial degree of elements
    READ (fUnit,*) nNodes,nSides,nElem,nGeo

    ! HOHQMesh reports interpolant data on Chebyshev Lobatto points
    ! We want data to be interpolated to Legendre Gauss Lobatto points
    CALL interp % Init(nGeo,CHEBYSHEV_GAUSS_LOBATTO,nGeo,GAUSS_LOBATTO)

    ! When we initialize the mesh, we set nNodes=nElem*4*(nGeo+1)**2 and
    ! nSides = nElem*4 since we still use `nNodes` and `nSides`
    ! in the input correspond to the HOPR definitions of these
    ! variables.
    ! `nSides` in HOHQMesh corresponds to nUniqueSides in HOPR and SELF
    ! `nNodes` in HOHQMesh corresponds to nCornerNodes (unique) in HOPR and SELF
    CALL myMesh % Init(nGeo,nElem,nElem*4,nElem*(nGeo + 1)**2,self_nBCsDefault)
    myMesh % nUniqueSides = nSides
    myMesh % nCornerNodes = nNodes

    CALL myMesh % hohq_cornerNodes % Alloc((/1,1/), &
                                           (/2,nNodes/))
    CALL myMesh % hohq_sideInfo % Alloc((/1,1/), &
                                        (/6,nSides/))
    CALL myMesh % hohq_elemInfo % Alloc((/1,1/), &
                                        (/4,nElem/))
    CALL myMesh % hohq_sideCurves % Alloc((/1,0,1,1/), &
                                          (/2,nGeo,4,nElem/))

    DO nid = 1,myMesh % nCornerNodes
      READ (fUnit,*) x
      myMesh % hohq_cornerNodes % hostData(1:2,nid) = x(1:2)
    END DO

    DO usid = 1,myMesh % nUniqueSides
      READ (fUnit,*) myMesh % hohq_sideInfo % hostData(1:6,usid)
    END DO

    DO eid = 1,myMesh % nElem
      READ (fUnit,*) myMesh % hohq_elemInfo % hostData(1:4,eid)
      READ (fUnit,*) bCurveFlag(1:4)
      DO lsid = 1,4
        IF (bCurveFlag(lsid) == 1) THEN
          DO i = 0,nGeo
            READ (fUnit,*) x
            myMesh % hohq_sideCurves % hostData(1:2,i,lsid,eid) = x(1:2)
          END DO
        ELSE

          ! For non-polynomial sides, create the side curve through interpolation between corner nodes
          lnid = myMesh % hopr_CGNSSideMap % hostData(1,lsid)
          nid = myMesh % hohq_elemInfo % hostData(lnid,eid)
          x0(1:2) = myMesh % hohq_cornerNodes % hostData(1:2,nid)

          lnid = myMesh % hopr_CGNSSideMap % hostData(2,lsid)
          nid = myMesh % hohq_elemInfo % hostData(lnid,eid)
          x1(1:2) = myMesh % hohq_cornerNodes % hostData(1:2,nid)

          DO i = 0,nGeo
            myMesh % hohq_sideCurves % hostData(1:2,i,lsid,eid) = 0.5_prec*( &
                                                           x0(1:2)*(1.0_prec - interp % controlPoints % hostData(i)) + &
                                                              x1(1:2)*(interp % controlPoints % hostData(i) + 1.0_prec))
          END DO

        END IF
      END DO
      READ (fUnit,*) line
      ! TO DO : Parse line for boundary conditions
    END DO

    CLOSE (fUnit)

    ! Fill in hopr_elemInfo
    ! To fill in the element info, we apply the assumption of a mesh with all quadrilateral elements
    ! with uniform polynomial order
    sid = 1
    nid = 1
    DO eid = 1,myMesh % nElem
      myMesh % hopr_elemInfo % hostData(1,eid) = selfQuadLinear ! Element Type
      myMesh % hopr_elemInfo % hostData(2,eid) = 1 ! Element Zone
      myMesh % hopr_elemInfo % hostData(3,eid) = sid ! Side Index Start
      sid = sid + 4
      myMesh % hopr_elemInfo % hostData(4,eid) = sid ! Side Index End
      myMesh % hopr_elemInfo % hostData(5,eid) = nid - 1 ! Node Index Start
      nid = nid + (nGeo + 1)**2
      myMesh % hopr_elemInfo % hostData(6,eid) = nid ! Node Index End
    END DO

    ! Generate the self_nodeCoords through transfinite interpolation with linear blending
    nid = 1
    DO eid = 1,myMesh % nElem
      lnid = 1
      ! Evaluate for corner points of mapping. This requires computational coordinates
      ! that include -1 and 1 in each computational direction (e.g. Gauss Lobatto)
      sw = myMesh % hohq_sideCurves % hostData(1:2,0,1,eid)
      se = myMesh % hohq_sideCurves % hostData(1:2,nGeo,1,eid)
      ne = myMesh % hohq_sideCurves % hostData(1:2,nGeo,3,eid)
      nw = myMesh % hohq_sideCurves % hostData(1:2,0,3,eid)
      DO j = 0,nGeo
        wSurf = myMesh % hohq_sideCurves % hostData(1:2,j,4,eid)
        eSurf = myMesh % hohq_sideCurves % hostData(1:2,j,2,eid)
        l2 = LinearBlend(interp % controlPoints % hostData(j))
        DO i = 0,nGeo
          sSurf = myMesh % hohq_sideCurves % hostData(1:2,i,1,eid)
          nSurf = myMesh % hohq_sideCurves % hostData(1:2,i,3,eid)
          l1 = LinearBlend(interp % controlPoints % hostData(i))

          P1 = l1(1)*wSurf + l1(2)*eSurf
          P2 = l2(1)*sSurf + l2(2)*nSurf

          ! Apply transfinite interpolation with linear blending
          myMesh % self_nodeCoords % hostData(1:2,lnid,eid) = P1 + P2 - &
                                                              l1(1)*(l2(1)*sw + l2(2)*nw) - l1(2)*(l2(1)*se + l2(2)*ne)

          myMesh % hopr_nodeCoords % hostData(1:2,nid) = myMesh % self_nodeCoords % hostData(1:2,lnid,eid)

          ! Increment node ids
          nid = nid + 1
          lnid = lnid + 1

        END DO
      END DO
    END DO

    CALL myMesh % GenerateConnectivity()

    CALL interp % Free()

  END SUBROUTINE Read_ISMv2_Mesh2D

  SUBROUTINE GenerateConnectivity_Mesh2D(myMesh)
#undef __FUNC__
#define __FUNC__ "GenerateConnectivity_Mesh2D"
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh
    ! Local
    INTEGER :: nid,unid,cnid,lnid,rnid
    INTEGER :: eid,sid,lsid,usid,gn1,gn2
    INTEGER :: key1,key2,e2,s2,e2gn1
    INTEGER :: nUniqueSides,flip
    INTEGER :: side(1:2,1:myMesh % nUniqueSides)
    TYPE(HashTable) :: sideTable

    ! Set the globalNodeIDs
    !  > Nodes in the nodeCoords list are possibly repeated since elements may share common sides. When the sides are shared
    !    the nodeCoords have the same x,y values (to machine precision)

    myMesh % nUniqueNodes = 1
    myMesh % hopr_globalNodeIDs % hostData(1) = 1
    DO nid = 2,myMesh % nNodes

      unid = myMesh % nUniqueNodes + 1
      DO rnid = 1,nid - 1
        IF (AlmostEqual(myMesh % hopr_nodeCoords % hostData(1,nid),myMesh % hopr_nodeCoords % hostData(1,rnid)) .AND. &
            AlmostEqual(myMesh % hopr_nodeCoords % hostData(2,nid),myMesh % hopr_nodeCoords % hostData(2,rnid))) THEN

          unid = myMesh % hopr_globalNodeIDs % hostData(rnid)
          EXIT

        END IF
      END DO

      myMesh % hopr_globalNodeIDs % hostData(nid) = unid
      IF (unid == myMesh % nUniqueNodes + 1) THEN
        myMesh % nUniqueNodes = myMesh % nUniqueNodes + 1
      END IF

    END DO

    CALL sideTable % Init(myMesh % nUniqueNodes)

    ! Set the sideInfo
    sid = 1
    nUniqueSides = 0
    DO eid = 1,myMesh % nElem
      DO lsid = 1,4

        cnid = myMesh % hopr_CGNSSideMap % hostData(1,lsid) ! Local Corner ID
        lnid = myMesh % hopr_CGNSCornerMap % hostData(cnid) ! Reference to Local Node ID
        nid = myMesh % hopr_elemInfo % hostData(5,eid) + lnid ! Add the offSetIndNODE to get the hopr node id
        gn1 = myMesh % hopr_globalNodeIDs % hostData(nid) ! Get the global Node ID for n1

        cnid = myMesh % hopr_CGNSSideMap % hostData(2,lsid) ! Local Corner ID
        lnid = myMesh % hopr_CGNSCornerMap % hostData(cnid) ! Reference to Local Node ID
        nid = myMesh % hopr_elemInfo % hostData(5,eid) + lnid ! Add the offSetIndNODE to get the hopr node id
        gn2 = myMesh % hopr_globalNodeIDs % hostData(nid) ! Get the global Node ID for n1

        ! Fill side info for eid
        myMesh % self_sideInfo % hostData(1,lsid,eid) = selfLineNonlinear ! SideType
        myMesh % self_sideInfo % hostData(2,lsid,eid) = 0 ! Global Side ID
        myMesh % self_sideInfo % hostData(3,lsid,eid) = 0 ! Neighbor Element ID
        myMesh % self_sideInfo % hostData(4,lsid,eid) = 0 ! Encoding for neighbor side ID and flip (10*s2+flip)
        myMesh % self_sideInfo % hostData(5,lsid,eid) = 1 ! Boundary Condition ID

        key1 = MIN(gn1,gn2)
        key2 = MAX(gn1,gn2)

        IF (sideTable % ContainsKeys(key1,key2)) THEN

          ! Get e2, s2, and flip
          CALL sideTable % GetDataForKeys(usid,key1,key2)

          e2 = side(1,usid)
          s2 = side(2,usid)

          ! Calculate flip
          ! > Get the starting global node ID for the other element
          cnid = myMesh % hopr_CGNSSideMap % hostData(1,s2) ! Local Corner ID
          lnid = myMesh % hopr_CGNSCornerMap % hostData(cnid) ! Reference to Local Node ID
          nid = myMesh % hopr_elemInfo % hostData(5,e2) + lnid ! Add the offSetIndNODE to get the hopr node id
          e2gn1 = myMesh % hopr_globalNodeIDs % hostData(nid) ! Get the global Node ID for n1

          flip = 0
          IF (e2gn1 /= gn1) THEN
            flip = 1
          END IF

          ! Populate information for this element
          myMesh % self_sideInfo % hostData(2,lsid,eid) = usid ! Global Side ID
          myMesh % self_sideInfo % hostData(3,lsid,eid) = e2 ! Neighbor Element ID
          myMesh % self_sideInfo % hostData(4,lsid,eid) = 10*s2 + flip ! Neighbor Element ID
          myMesh % self_sideInfo % hostData(5,lsid,eid) = 0 ! boundary condition id

          ! Population information for the other element
          myMesh % self_sideInfo % hostData(2,s2,e2) = usid ! Global Side ID
          myMesh % self_sideInfo % hostData(3,s2,e2) = eid ! Neighbor Element ID
          myMesh % self_sideInfo % hostData(4,s2,e2) = 10*lsid + flip ! Neighbor Element ID
          myMesh % self_sideInfo % hostData(5,s2,e2) = 0 ! boundary condition id

        ELSE

          nUniqueSides = nUniqueSides + 1
          side(1,nUniqueSides) = eid ! Store the element ID
          side(2,nUniqueSides) = lsid ! Store the local side ID
          ! Add the side to the hash table

          CALL sideTable % AddDataForKeys(nUniqueSides,key1,key2)

        END IF

        sid = sid + 1

      END DO
    END DO

    CALL sideTable % Free()

  END SUBROUTINE GenerateConnectivity_Mesh2D

  SUBROUTINE Read_HOPr_Mesh2D_serial(myMesh,meshFile)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    ! Adapted for 2D Mesh : Note that HOPR does not have 2D mesh output.
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(out) :: myMesh
    CHARACTER(*),INTENT(in) :: meshFile
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER :: nGlobalElem
    INTEGER :: nLocalNodes
    INTEGER :: nLocalSides
    INTEGER :: nUniqueSides
    INTEGER :: nGeo,nBCs
    INTEGER :: eid, lsid, iSide
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfInt32_r2) :: hopr_sideInfo
    TYPE(hfReal_r2) :: hopr_nodeCoords
    TYPE(hfInt32_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: bcType

    CALL Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId)

    CALL ReadAttribute_HDF5(fileId,'nElems',nGlobalElem)
    CALL ReadAttribute_HDF5(fileId,'Ngeo',nGeo)
    CALL ReadAttribute_HDF5(fileId,'nBCs',nBCs)
    CALL ReadAttribute_HDF5(fileId,'nUniqueSides',nUniqueSides)

    ! Read BCType
    CALL bcType % Alloc(loBound=(/1,1/), &
                        upBound=(/4,nBCs/))
    CALL ReadArray_HDF5(fileId,'BCType',bcType)

    ! Allocate Space for hopr_elemInfo!
    CALL hopr_elemInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/6,nGlobalElem/))

    CALL ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    nLocalNodes = hopr_elemInfo % hostData(6,nGlobalElem) - hopr_elemInfo % hostData(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !
    CALL hopr_nodeCoords % Alloc(loBound=(/1,1/), &
                                 upBound=(/2,nLocalNodes/))

    CALL hopr_globalNodeIDs % Alloc(loBound=1, &
                                    upBound=nLocalNodes)

    CALL ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords)
    CALL ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs)

    ! Read local subarray of SideInfo
    nLocalSides = hopr_elemInfo % hostData(4,nGlobalElem) - hopr_elemInfo % hostData(3,1)

    ! Allocate space for hopr_sideInfo
    CALL hopr_sideInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/5,nLocalSides/))

    CALL ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo)

    CALL Close_HDF5(fileID)

    CALL myMesh % Init(nGeo,nGlobalElem,nLocalSides,nLocalNodes,nBCs)

    ! Copy data from local arrays into myMesh
    myMesh % hopr_elemInfo = hopr_elemInfo
    myMesh % hopr_nodeCoords = hopr_nodeCoords
    myMesh % hopr_globalNodeIDs = hopr_globalNodeIDs
    myMesh % hopr_sideInfo = hopr_sideInfo
    myMesh % nUniqueSides = nUniqueSides

    iSide = 0 
    DO eid = 1,myMesh % nElem
      DO lsid = 1,4
        iSide = iSide + 1
        myMesh % self_sideInfo % hostData(1:5,lsid,eid) = myMesh % hopr_sideInfo % hostData(1:5,iSide)
      ENDDO
    ENDDO

    CALL myMesh % UpdateDevice()

    CALL hopr_elemInfo % Free()
    CALL hopr_nodeCoords % Free()
    CALL hopr_globalNodeIDs % Free()
    CALL hopr_sideInfo % Free()

  END SUBROUTINE Read_HOPr_Mesh2D_serial

  SUBROUTINE Read_HOPr_Mesh2D_parallel(myMesh,meshFile,decomp)
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
    INTEGER :: firstElem,nLocalElems
    INTEGER :: firstNode,nLocalNodes
    INTEGER :: firstSide,nLocalSides
    INTEGER :: nUniqueSides
    INTEGER :: nGeo,nBCs
    INTEGER :: eid, lsid, iSide
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfInt32_r2) :: hopr_sideInfo
    TYPE(hfReal_r2) :: hopr_nodeCoords
    TYPE(hfInt32_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: bcType

    CALL Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId,decomp % mpiComm)

    CALL ReadAttribute_HDF5(fileId,'nElems',nGlobalElem)
    CALL ReadAttribute_HDF5(fileId,'Ngeo',nGeo)
    CALL ReadAttribute_HDF5(fileId,'nBCs',nBCs)
    CALL ReadAttribute_HDF5(fileId,'nUniqueSides',nUniqueSides)

    ! Read BCType
    CALL bcType % Alloc(loBound=(/1,1/), &
                        upBound=(/4,nBCs/))
    offset(:) = 0
    CALL ReadArray_HDF5(fileId,'BCType',bcType,offset)

    ! Read local subarray of ElemInfo
    CALL decomp % SetElemToRank(nGlobalElem)
    firstElem = decomp % offsetElem % hostData(decomp % rankId) + 1
    nLocalElems = decomp % offsetElem % hostData(decomp % rankId + 1) - &
                  decomp % offsetElem % hostData(decomp % rankId)

    ! Allocate Space for hopr_elemInfo!
    CALL hopr_elemInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/6,nLocalElems/))

    offset = (/0,firstElem - 1/)
    CALL ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo,offset)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode = hopr_elemInfo % hostData(5,1) + 1
    nLocalNodes = hopr_elemInfo % hostData(6,nLocalElems) - hopr_elemInfo % hostData(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !
    CALL hopr_nodeCoords % Alloc(loBound=(/1,1/), &
                                 upBound=(/2,nLocalNodes/))

    CALL hopr_globalNodeIDs % Alloc(loBound=1, &
                                    upBound=nLocalNodes)

    offset = (/0,firstNode - 1/)
    CALL ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords,offset)
    gOffset = (/firstNode - 1/)
    CALL ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs,gOffset)

    ! Read local subarray of SideInfo
    firstSide = hopr_elemInfo % hostData(3,1) + 1
    nLocalSides = hopr_elemInfo % hostData(4,nLocalElems) - hopr_elemInfo % hostData(3,1)

    ! Allocate space for hopr_sideInfo
    CALL hopr_sideInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/5,nLocalSides/))
    offset = (/0,firstSide - 1/)
    CALL ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo,offset)

    CALL Close_HDF5(fileID)

    CALL myMesh % Init(nGeo,nLocalElems,nLocalSides,nLocalNodes,nBCs)

    ! Copy data from local arrays into myMesh
    myMesh % hopr_elemInfo = hopr_elemInfo
    myMesh % hopr_nodeCoords = hopr_nodeCoords
    myMesh % hopr_globalNodeIDs = hopr_globalNodeIDs
    myMesh % hopr_sideInfo = hopr_sideInfo
    myMesh % nUniqueSides = nUniqueSides
    myMesh % nGlobalElem = nGlobalElem

    iSide = 0 
    DO eid = 1,myMesh % nElem
      DO lsid = 1,4
        iSide = iSide + 1
        myMesh % self_sideInfo % hostData(1:5,lsid,eid) = myMesh % hopr_sideInfo % hostData(1:5,iSide)
      ENDDO
    ENDDO

    CALL myMesh % UpdateDevice()

    CALL hopr_elemInfo % Free()
    CALL hopr_nodeCoords % Free()
    CALL hopr_globalNodeIDs % Free()
    CALL hopr_sideInfo % Free()

  END SUBROUTINE Read_HOPr_Mesh2D_parallel

  SUBROUTINE Write_HOPr_Mesh2D(myMesh,meshFile)
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
    CALL WriteArray_HDF5(fileId,'ElemInfo',myMesh % hopr_elemInfo)

    ! Write local subarray of NodeCoords and GlobalNodeIDs
    CALL WriteArray_HDF5(fileId,'NodeCoords',myMesh % hopr_nodeCoords)
    CALL WriteArray_HDF5(fileId,'GlobalNodeIDs',myMesh % hopr_globalNodeIDs)

    ! Write local subarray of SideInfo
    CALL WriteArray_HDF5(fileId,'SideInfo',myMesh % hopr_sideInfo)

    CALL Close_HDF5(fileID)

  END SUBROUTINE Write_HOPr_Mesh2D

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

    CALL myMesh % hopr_elemInfo % Alloc(loBound=(/1,1/), &
                                        upBound=(/6,nElem/))

    CALL myMesh % hopr_sideInfo % Alloc(loBound=(/1,1/), &
                                        upBound=(/5,nSides/))

    CALL myMesh % self_sideInfo % Alloc(loBound=(/1,1,1/), &
                                        upBound=(/5,6,nElem/))

    CALL myMesh % hopr_nodeCoords % Alloc(loBound=(/1,1/), &
                                          upBound=(/3,nNodes/))

    CALL myMesh % self_nodeCoords % Alloc(loBound=(/1,1,1/), &
                                          upBound=(/3, (nGeo + 1)**3,nElem/))

    CALL myMesh % hopr_globalNodeIDs % Alloc(loBound=1, &
                                             upBound=nNodes)

    CALL myMesh % hopr_CGNSCornerMap % Alloc(loBound=1, &
                                             upBound=8)

    CALL myMesh % self_sideMap % Alloc(loBound=(/1,1/), &
                                       upBound=(/4,6/))

    CALL myMesh % hopr_CGNSSideMap % Alloc(loBound=(/1,1/), &
                                           upBound=(/4,6/))

    CALL myMesh % hopr_curveNodeMap % Alloc(loBound=(/1,1/), &
                                            upBound=(/3, (nGeo + 1)**3/))

    CALL myMesh % hopr_curveNodeMapInv % Alloc(loBound=(/0,0,0/), &
                                               upBound=(/nGeo,nGeo,nGeo/))

    CALL myMesh % BCType % Alloc(loBound=(/1,1/), &
                                 upBound=(/4,nBCs/))

    ALLOCATE (myMesh % BCNames(1:nBCs))

    ! Create lookup tables to assist with connectivity generation
    myMesh % hopr_CGNSCornerMap % hostData(1) = 1
    myMesh % hopr_CGNSCornerMap % hostData(2) = nGeo + 1
    myMesh % hopr_CGNSCornerMap % hostData(3) = (nGeo + 1)**2
    myMesh % hopr_CGNSCornerMap % hostData(4) = nGeo*(nGeo + 1) + 1
    myMesh % hopr_CGNSCornerMap % hostData(5) = nGeo*(nGeo + 1)**2 + 1
    myMesh % hopr_CGNSCornerMap % hostData(6) = nGeo*(nGeo + 1)**2 + (nGeo + 1)
    myMesh % hopr_CGNSCornerMap % hostData(7) = (nGeo + 1)**3
    myMesh % hopr_CGNSCornerMap % hostData(8) = nGeo*(nGeo + 1)*(nGeo + 2) + 1

    l = 0
    DO k = 0,nGeo
      DO j = 0,nGeo
        DO i = 0,nGeo
          l = l + 1
          myMesh % hopr_curveNodeMap % hostData(1:3,l) = (/i,j,k/)
          myMesh % hopr_curveNodeMapInv % hostData(i,j,k) = l
        END DO
      END DO
    END DO

    ! Maps from local corner node id to CGNS side
    myMesh % hopr_CGNSSideMap % hostData(1:4,1) = (/1,4,3,2/)
    myMesh % hopr_CGNSSideMap % hostData(1:4,2) = (/1,2,6,5/)
    myMesh % hopr_CGNSSideMap % hostData(1:4,3) = (/2,3,7,6/)
    myMesh % hopr_CGNSSideMap % hostData(1:4,4) = (/3,4,8,7/)
    myMesh % hopr_CGNSSideMap % hostData(1:4,5) = (/1,5,8,4/)
    myMesh % hopr_CGNSSideMap % hostData(1:4,6) = (/5,6,7,8/)

    myMesh % self_sideMap % hostData(1:4,1) = (/1,2,3,4/) ! Bottom
    myMesh % self_sideMap % hostData(1:4,2) = (/1,2,6,5/) ! South
    myMesh % self_sideMap % hostData(1:4,3) = (/2,3,7,6/) ! East
    myMesh % self_sideMap % hostData(1:4,4) = (/4,3,7,8/) ! North
    myMesh % self_sideMap % hostData(1:4,5) = (/1,4,8,5/) ! West
    myMesh % self_sideMap % hostData(1:4,6) = (/5,6,7,8/) ! Top

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

    CALL myMesh % hopr_elemInfo % Free()
    CALL myMesh % hopr_sideInfo % Free()
    CALL myMesh % self_sideInfo % Free()
    CALL myMesh % hopr_nodeCoords % Free()
    CALL myMesh % self_nodeCoords % Free()
    CALL myMesh % hopr_CGNSCornerMap % Free()
    CALL myMesh % hopr_globalNodeIDs % Free()
    CALL myMesh % self_sideMap % Free()
    CALL myMesh % BCType % Free()

    DEALLOCATE (myMesh % BCNames)

  END SUBROUTINE Free_Mesh3D
  SUBROUTINE UpdateHost_Mesh3D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh

    CALL myMesh % hopr_elemInfo % UpdateHost()
    CALL myMesh % hopr_sideInfo % UpdateHost()
    CALL myMesh % self_sideInfo % UpdateHost()
    CALL myMesh % hopr_nodeCoords % UpdateHost()
    CALL myMesh % self_nodeCoords % UpdateHost()
    CALL myMesh % hopr_globalNodeIDs % UpdateHost()
    CALL myMesh % BCType % UpdateHost()

  END SUBROUTINE UpdateHost_Mesh3D

  SUBROUTINE UpdateDevice_Mesh3D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh

    CALL myMesh % hopr_elemInfo % UpdateDevice()
    CALL myMesh % hopr_sideInfo % UpdateDevice()
    CALL myMesh % self_sideInfo % UpdateDevice()
    CALL myMesh % hopr_nodeCoords % UpdateDevice()
    CALL myMesh % self_nodeCoords % UpdateDevice()
    CALL myMesh % hopr_globalNodeIDs % UpdateDevice()
    CALL myMesh % BCType % UpdateDevice()

  END SUBROUTINE UpdateDevice_Mesh3D

  SUBROUTINE UniformBlockMesh_Mesh3D(myMesh,nGeo,nElem,x)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(out) :: myMesh
    INTEGER,INTENT(in) :: nGeo
    INTEGER,INTENT(in) :: nElem(1:3)
    REAL(prec),INTENT(in) :: x(1:6)
    ! Local
    INTEGER :: iel,jel,kel,nEl,elid
    INTEGER :: sid,nid,nNodes
    INTEGER :: nSides
    INTEGER :: bcid
    INTEGER :: i,j,k
    INTEGER :: nbeid, nbsid, usid, lsid
    REAL(prec) :: xU(1:nElem(1) + 1)
    REAL(prec) :: yU(1:nElem(2) + 1)
    REAL(prec) :: zU(1:nElem(3) + 1)
    TYPE(Vector3D) :: xLinear
    TYPE(Vector3D) :: xGeo

    nEl = nElem(1)*nElem(2)*nElem(3)
    nNodes = nEl*(nGeo + 1)*(nGeo + 1)*(nGeo + 1)
    nSides = nEl*6
    CALL myMesh % Init(nGeo,nEl,nSides,nNodes,1)
    myMesh % nUniqueSides = (nElem(1) + 1)*nElem(2)*nElem(3) + &
             (nElem(2) + 1)*nElem(1)*nElem(3) + &
             (nElem(3) + 1)*nElem(1)*nElem(2)

    ! Set the hopr_nodeCoords
    xU = UniformPoints(x(1),x(2),1,nElem(1) + 1)
    yU = UniformPoints(x(3),x(4),1,nElem(2) + 1)
    zU = UniformPoints(x(5),x(6),1,nElem(3) + 1)

    ! Create a linear interpolant to interpolate to nGeo grid
    CALL xLinear % Init(1,CHEBYSHEV_GAUSS_LOBATTO, &
                        nGeo,CHEBYSHEV_GAUSS_LOBATTO, &
                        1,nEl)

    CALL xGeo % Init(nGeo,CHEBYSHEV_GAUSS_LOBATTO, &
                     nGeo,CHEBYSHEV_GAUSS_LOBATTO, &
                     1,nEl)
    elid = 1
    DO kel = 1,nElem(3)
      DO jel = 1,nElem(2)
        DO iel = 1,nElem(1)
          ! x component
          xLinear % interior % hostData(1,0:1,0,0,1,elid) = xU(iel:iel + 1)
          xLinear % interior % hostData(1,0:1,1,0,1,elid) = xU(iel:iel + 1)
          xLinear % interior % hostData(1,0:1,0,1,1,elid) = xU(iel:iel + 1)
          xLinear % interior % hostData(1,0:1,1,1,1,elid) = xU(iel:iel + 1)
          ! y component
          xLinear % interior % hostData(2,0,0:1,0,1,elid) = yU(jel:jel + 1)
          xLinear % interior % hostData(2,1,0:1,0,1,elid) = yU(jel:jel + 1)
          xLinear % interior % hostData(2,0,0:1,1,1,elid) = yU(jel:jel + 1)
          xLinear % interior % hostData(2,1,0:1,1,1,elid) = yU(jel:jel + 1)
          ! z component
          xLinear % interior % hostData(3,0,0,0:1,1,elid) = zU(kel:kel + 1)
          xLinear % interior % hostData(3,1,0,0:1,1,elid) = zU(kel:kel + 1)
          xLinear % interior % hostData(3,0,1,0:1,1,elid) = zU(kel:kel + 1)
          xLinear % interior % hostData(3,1,1,0:1,1,elid) = zU(kel:kel + 1)
          ! Incremenent the element ID
          elid = elid + 1
        END DO
      END DO
    END DO

    CALL xLinear % GridInterp(xGeo,.FALSE.)

    ! Set the element information
    nid = 1
    sid = 0
    elid = 1
    DO kel = 1,nElem(3)
      DO jel = 1,nElem(2)
        DO iel = 1,nElem(1)
          myMesh % hopr_elemInfo % hostData(1,elid) = selfQuadLinear ! Element Type
          myMesh % hopr_elemInfo % hostData(2,elid) = 1 ! Element Zone
          myMesh % hopr_elemInfo % hostData(3,elid) = sid ! Side Index Start
          sid = sid + 6
          myMesh % hopr_elemInfo % hostData(4,elid) = sid ! Side Index End
          myMesh % hopr_elemInfo % hostData(5,elid) = nid - 1 ! Node Index Start
          DO k = 0,nGeo
            DO j = 0,nGeo
              DO i = 0,nGeo
                myMesh % hopr_nodeCoords % hostData(1:3,nid) = xGeo % interior % hostData(1:3,i,j,k,1,elid)
                nid = nid + 1
              END DO
            END DO
          END DO
          myMesh % hopr_elemInfo % hostData(6,elid) = nid ! Node Index End
          elid = elid + 1
        END DO
      END DO
    END DO

    ! Set up the side info
    elid = 1
    sid = 0
    usid = 0
    DO kel = 1,nElem(3)
      DO jel = 1,nElem(2)
        DO iel = 1,nElem(1)

          ! Bottom Face ! Local Side = 1
          IF (kel == 1) THEN

             sid = sid + 1
             usid = usid + 1 ! unique side id
             myMesh % hopr_sideInfo % hostData(1,sid) = selfQuadLinear ! Side type set to linear quad
             myMesh % hopr_sideInfo % hostData(2,sid) = usid ! Unique side id
             myMesh % hopr_sideInfo % hostData(3,sid) = 0 ! Neighbor Element ID (0=boundary)
             myMesh % hopr_sideInfo % hostData(4,sid) = 0 ! 10*nbLocalSide + flip
             myMesh % hopr_sideInfo % hostData(5,sid) = self_BCDefault ! Boundary condition ID

          ELSE

             sid = sid + 1
             nbeid = iel + nElem(1)*(jel-1 + nElem(2)*(kel-2)) ! Get the element id for the element below
             nbsid = myMesh % hopr_elemInfo % hostData(3,nbeid) + selfSide3D_Top-1 ! Get sid for the top face of the element below
             usid = myMesh % hopr_sideInfo % hostData(2,nbsid) ! Get the unique side address 
             myMesh % hopr_sideInfo % hostData(1,sid) = selfQuadLinear ! Side type set to linear quad
             myMesh % hopr_sideInfo % hostData(2,sid) = usid ! Unique side id
             myMesh % hopr_sideInfo % hostData(3,sid) = nbeid ! Neighbor Element ID (0=boundary)
             ! flip == 1 indicates internal geometry of each neighboring element
             ! has the same orientation
             myMesh % hopr_sideInfo % hostData(4,sid) = 10*selfSide3D_Top + 1 ! 10*nbLocalSide + flip
             myMesh % hopr_sideInfo % hostData(5,sid) = 0 ! Boundary condition ID

          ENDIF

          ! South Face ! Local Side = 2
          IF (jel == 1) THEN

             sid = sid + 1
             usid = usid + 1 ! unique side id
             myMesh % hopr_sideInfo % hostData(1,sid) = selfQuadLinear ! Side type set to linear quad
             myMesh % hopr_sideInfo % hostData(2,sid) = usid ! Unique side id
             myMesh % hopr_sideInfo % hostData(3,sid) = 0 ! Neighbor Element ID (0=boundary)
             myMesh % hopr_sideInfo % hostData(4,sid) = 0 ! 10*nbLocalSide + flip
             myMesh % hopr_sideInfo % hostData(5,sid) = self_BCDefault ! Boundary condition ID

          ELSE

             sid = sid + 1
             nbeid = iel + nElem(1)*(jel-2 + nElem(2)*(kel-1)) ! Get the element id for the element to the south
             nbsid = myMesh % hopr_elemInfo % hostData(3,nbeid) + selfSide3D_North-1 ! Get sid for the north face of the element to  the south
             usid = myMesh % hopr_sideInfo % hostData(2,nbsid) ! Get the unique side address 
             myMesh % hopr_sideInfo % hostData(1,sid) = selfQuadLinear ! Side type set to linear quad
             myMesh % hopr_sideInfo % hostData(2,sid) = usid ! Unique side id
             myMesh % hopr_sideInfo % hostData(3,sid) = nbeid ! Neighbor Element ID (0=boundary)
             myMesh % hopr_sideInfo % hostData(4,sid) = 10*selfSide3D_North ! 10*nbLocalSide + flip
             myMesh % hopr_sideInfo % hostData(5,sid) = 0 ! Boundary condition ID

          ENDIF

          ! East Face ! Local Side = 3
          sid = sid + 1
          IF (iel == nElem(1)) THEN
            nbeid = 0
            bcid = self_BCDefault
          ELSE   
            nbeid = iel + 1 + nElem(1)*(jel-1 + nElem(2)*(kel-1)) ! Get the element id for the element to the east
            bcid = 0
          ENDIF
          usid = usid + 1 ! unique side id
          myMesh % hopr_sideInfo % hostData(1,sid) = selfQuadLinear ! Side type set to linear quad
          myMesh % hopr_sideInfo % hostData(2,sid) = usid ! Unique side id
          myMesh % hopr_sideInfo % hostData(3,sid) = nbeid ! Neighbor Element ID (0=boundary)
          myMesh % hopr_sideInfo % hostData(4,sid) = 10*selfSide3D_West ! 10*nbLocalSide + flip
          myMesh % hopr_sideInfo % hostData(5,sid) = bcid ! Boundary condition ID


          ! North Face ! Local Side = 4
          sid = sid + 1
          IF (jel == nElem(2)) THEN
            nbeid = 0
            bcid = self_BCDefault
          ELSE   
            nbeid = iel + nElem(1)*(jel + nElem(2)*(kel-1)) ! Get the element id for the element to the north
            bcid = 0
          ENDIF
          usid = usid + 1 ! unique side id
          myMesh % hopr_sideInfo % hostData(1,sid) = selfQuadLinear ! Side type set to linear quad
          myMesh % hopr_sideInfo % hostData(2,sid) = usid ! Unique side id
          myMesh % hopr_sideInfo % hostData(3,sid) = nbeid ! Neighbor Element ID (0=boundary)
          myMesh % hopr_sideInfo % hostData(4,sid) = 10*selfSide3D_South ! 10*nbLocalSide + flip
          myMesh % hopr_sideInfo % hostData(5,sid) = bcid ! Boundary condition ID

          ! West Face ! Local Side = 5
          IF (iel == 1) THEN

             sid = sid + 1
             usid = usid + 1 ! unique side id
             myMesh % hopr_sideInfo % hostData(1,sid) = selfQuadLinear ! Side type set to linear quad
             myMesh % hopr_sideInfo % hostData(2,sid) = usid ! Unique side id
             myMesh % hopr_sideInfo % hostData(3,sid) = 0 ! Neighbor Element ID (0=boundary)
             myMesh % hopr_sideInfo % hostData(4,sid) = 0 ! 10*nbLocalSide + flip
             myMesh % hopr_sideInfo % hostData(5,sid) = self_BCDefault ! Boundary condition ID

          ELSE

             sid = sid + 1
             nbeid = iel - 1 + nElem(1)*(jel-1 + nElem(2)*(kel-1)) ! Get the element id for the element to the west
             nbsid = myMesh % hopr_elemInfo % hostData(3,nbeid) + selfSide3D_East-1 ! Get sid for the east face of the element to the west
             usid = myMesh % hopr_sideInfo % hostData(2,nbsid) ! Get the unique side address 
             myMesh % hopr_sideInfo % hostData(1,sid) = selfQuadLinear ! Side type set to linear quad
             myMesh % hopr_sideInfo % hostData(2,sid) = usid ! Unique side id
             myMesh % hopr_sideInfo % hostData(3,sid) = nbeid ! Neighbor Element ID (0=boundary)
             myMesh % hopr_sideInfo % hostData(4,sid) = 10*selfSide3D_East ! 10*nbLocalSide + flip
             myMesh % hopr_sideInfo % hostData(5,sid) = 0 ! Boundary condition ID

          ENDIF

          ! Top Face ! Local Side = 6
          sid = sid + 1
          IF (kel == nElem(3)) THEN
            nbeid = 0
            bcid = self_BCDefault
          ELSE   
            nbeid = iel + nElem(1)*(jel-1 + nElem(2)*(kel)) ! Get the element id for the element above
            bcid = 0
          ENDIF
          usid = usid + 1 ! unique side id
          myMesh % hopr_sideInfo % hostData(1,sid) = selfQuadLinear ! Side type set to linear quad
          myMesh % hopr_sideInfo % hostData(2,sid) = usid ! Unique side id
          myMesh % hopr_sideInfo % hostData(3,sid) = nbeid ! Neighbor Element ID (0=boundary)
          myMesh % hopr_sideInfo % hostData(4,sid) = 10*selfSide3D_Bottom ! 10*nbLocalSide + flip
          myMesh % hopr_sideInfo % hostData(5,sid) = bcid ! Boundary condition ID


        ENDDO
      ENDDO
    ENDDO

    elid = 0
    sid = 0
    DO kel = 1,nElem(3)
      DO jel = 1,nElem(2)
        DO iel = 1,nElem(1)
          elid = elid + 1
          DO lsid = 1, 6
            sid = sid + 1
            myMesh % self_sideInfo % hostData(1:5,lsid,elid) = myMesh % hopr_sideInfo % hostData(1:5,sid)
          ENDDO
        ENDDO
      ENDDO
    ENDDO


    CALL myMesh % UpdateDevice()

    CALL xLinear % Free()
    CALL xGeo % Free()

  END SUBROUTINE UniformBlockMesh_Mesh3D

  SUBROUTINE Load_Mesh3D_serial(myMesh,myMeshSpec)
#undef __FUNC__
#define __FUNC__ "Load_Mesh3D"
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(out) :: myMesh
    TYPE(MeshSpec),INTENT(in) :: myMeshSpec

    IF (myMeshSpec % blockMesh) THEN

      CALL myMesh % UniformBlockMesh(myMeshSpec % blockMesh_nGeo, &
                        (/myMeshSpec % blockMesh_nElemX,myMeshSpec % blockMesh_nElemY,myMeshSpec % blockMesh_nElemZ/), &
                                     (/myMeshSpec % blockMesh_x0,myMeshSpec % blockMesh_x1, &
                                       myMeshSpec % blockMesh_y0,myMeshSpec % blockMesh_y1, &
                                       myMeshSpec % blockMesh_z0,myMeshSpec % blockMesh_z1/))

    ELSE

      CALL myMesh % Read_HOPr(myMeshSpec % filename)

    END IF

  END SUBROUTINE Load_Mesh3D_serial

  SUBROUTINE Load_Mesh3D_parallel(myMesh,myMeshSpec,decomp)
#undef __FUNC__
#define __FUNC__ "Load_Mesh3D"
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(out) :: myMesh
    TYPE(MeshSpec),INTENT(in) :: myMeshSpec
    TYPE(MPILayer),INTENT(inout) :: decomp

    IF (myMeshSpec % blockMesh) THEN

      IF (decomp % nRanks > 1) THEN ! Error out
        ERROR("Block Mesh only supported in serial")
        STOP ! TO DO : Safe exit for serial and parallel
      ELSE
        CALL myMesh % UniformBlockMesh(myMeshSpec % blockMesh_nGeo, &
                        (/myMeshSpec % blockMesh_nElemX,myMeshSpec % blockMesh_nElemY,myMeshSpec % blockMesh_nElemZ/), &
                                       (/myMeshSpec % blockMesh_x0,myMeshSpec % blockMesh_x1, &
                                         myMeshSpec % blockMesh_y0,myMeshSpec % blockMesh_y1, &
                                         myMeshSpec % blockMesh_z0,myMeshSpec % blockMesh_z1/))
      END IF

    ELSE

      CALL myMesh % Read_HOPr(myMeshSpec % filename,decomp)

    END IF

  END SUBROUTINE Load_Mesh3D_parallel

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
    INTEGER :: i
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

        e2Global = myMesh % self_sideInfo % hostData(3,s1,e1)
        e2 = e2Global - offset
        s2 = myMesh % self_sideInfo % hostData(4,s1,e1)/10
        flip = myMesh % self_sideInfo % hostData(4,s1,e1) - s2*10
        bcid = myMesh % self_sideInfo % hostData(5,s1,e1)

        IF (bcid == 0) THEN

          IF (PRESENT(decomp)) THEN
            neighborRank = decomp % elemToRank % hostData(e2Global)
          ELSE
            neighborRank = 0
          ENDIF

          IF (neighborRank == rankId) THEN

            ! With 8 nodes per element, and the nodes provided in order, we can also shift the node indices
            n1Global = myMesh % hopr_elemInfo % hostData(5,e1) ! Starting node index for element 1
            n1 = n1Global - 8*offset

            n2Global = myMesh % hopr_elemInfo % hostData(5,e2) ! Starting node index for element 2
            n2 = n2Global - 8*offset

            lnid1 = myMesh % self_sideMap % hostData(1:4,s1) ! local CGNS corner node ids for element 1 side
            lnid2 = myMesh % self_sideMap % hostData(1:4,s2) ! local CGNS corner node ids for element 2 side

            DO l = 1, 4
              
              c1 = myMesh % hopr_CGNSCornerMap % hostData(lnid1(l)) ! Get the local HOPR node id for element 1
              c2 = myMesh % hopr_CGNSCornerMap % hostData(lnid2(l)) ! Get the local HOPR node id for element 2
              nid1(l,s1,e1) = myMesh % hopr_globalNodeIDs % hostData(n1+c1) ! Global node IDs for element 1 side
              nid2(l,s1,e1) = myMesh % hopr_globalNodeIDs % hostData(n2+c2) ! Global node IDs for element 2 side

            ENDDO

          ELSE ! In this case, we need to exchange

            globalSideId = ABS(myMesh % self_sideInfo % hostdata(2,s1,e1))

            n1Global = myMesh % hopr_elemInfo % hostData(5,e1) ! Starting node index for element 1
            n1 = n1Global - 8*offset

            lnid1 = myMesh % self_sideMap % hostData(1:4,s1) ! local CGNS corner node ids for element 1 side

            DO l = 1, 4
              
              c1 = myMesh % hopr_CGNSCornerMap % hostData(lnid1(l)) ! Get the local HOPR node id for element 1
              nid1(l,s1,e1) = myMesh % hopr_globalNodeIDs % hostData(n1+c1) ! Global node IDs for element 1 side

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

        s2 = myMesh % self_sideInfo % hostData(4,s1,e1)/10
        bcid = myMesh % self_sideInfo % hostData(5,s1,e1)
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

          myMesh % self_sideInfo % hostData(4,s1,e1) = 10*s2+nShifts

        ENDIF

      ENDDO
    ENDDO

    DEALLOCATE(requests)
    DEALLOCATE(stats)

  END SUBROUTINE RecalculateFlip_Mesh3D

  SUBROUTINE Read_HOPr_Mesh3D_serial(myMesh,meshFile)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(out) :: myMesh
    CHARACTER(*),INTENT(in) :: meshFile
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER :: nGlobalElem
    INTEGER :: nLocalNodes
    INTEGER :: nLocalSides
    INTEGER :: nUniqueSides
    INTEGER :: nGeo,nBCs
    INTEGER :: eid, lsid, iSide
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfInt32_r2) :: hopr_sideInfo
    TYPE(hfReal_r2) :: hopr_nodeCoords
    TYPE(hfInt32_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: bcType

    CALL Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId)

    CALL ReadAttribute_HDF5(fileId,'nElems',nGlobalElem)
    CALL ReadAttribute_HDF5(fileId,'Ngeo',nGeo)
    CALL ReadAttribute_HDF5(fileId,'nBCs',nBCs)
    CALL ReadAttribute_HDF5(fileId,'nUniqueSides',nUniqueSides)

    PRINT*, 'Reading mesh file '//TRIM(meshFile)
    PRINT*, 'Number of Elements : ', nGlobalElem
    PRINT*, 'Input polynomial degree : ', nGeo
    PRINT*, 'Number of BC types : ', nBCs

    ! Read BCType
    CALL bcType % Alloc(loBound=(/1,1/), &
                        upBound=(/4,nBCs/))
    CALL ReadArray_HDF5(fileId,'BCType',bcType)

    ! Allocate Space for hopr_elemInfo!
    CALL hopr_elemInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/6,nGlobalElem/))

    CALL ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    nLocalNodes = hopr_elemInfo % hostData(6,nGlobalElem) - hopr_elemInfo % hostData(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !
    CALL hopr_nodeCoords % Alloc(loBound=(/1,1/), &
                                 upBound=(/3,nLocalNodes/))

    CALL hopr_globalNodeIDs % Alloc(loBound=1, &
                                    upBound=nLocalNodes)

    CALL ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords)
    CALL ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs)
    ! Read local subarray of SideInfo
    nLocalSides = hopr_elemInfo % hostData(4,nGlobalElem) - hopr_elemInfo % hostData(3,1)
    PRINT*, 'Number of local sides : ', nLocalSides

    ! Allocate space for hopr_sideInfo
    CALL hopr_sideInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/5,nLocalSides/))

    CALL ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo)

    CALL Close_HDF5(fileID)

    CALL myMesh % Init(nGeo,nGlobalElem,nLocalSides,nLocalNodes,nBCs)

    ! Copy data from local arrays into myMesh
    myMesh % hopr_elemInfo % hostData = hopr_elemInfo % hostData
    myMesh % hopr_nodeCoords % hostData = hopr_nodeCoords % hostData
    myMesh % hopr_globalNodeIDs % hostData = hopr_globalNodeIDs % hostData
    myMesh % hopr_sideInfo % hostData = hopr_sideInfo % hostData
    myMesh % nUniqueSides = nUniqueSides

    iSide = 0 
    DO eid = 1,myMesh % nElem
      DO lsid = 1,6
        iSide = iSide + 1
        myMesh % self_sideInfo % hostData(1:5,lsid,eid) = myMesh % hopr_sideInfo % hostData(1:5,iSide)
      ENDDO
    ENDDO

    CALL myMesh % RecalculateFlip()

    CALL myMesh % UpdateDevice()

    CALL hopr_elemInfo % Free()
    CALL hopr_nodeCoords % Free()
    CALL hopr_globalNodeIDs % Free()
    CALL hopr_sideInfo % Free()

  END SUBROUTINE Read_HOPr_Mesh3D_serial

  SUBROUTINE Read_HOPr_Mesh3D_parallel(myMesh,meshFile,decomp)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(out) :: myMesh
    CHARACTER(*),INTENT(in) :: meshFile
    TYPE(MPILayer),INTENT(inout) :: decomp
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: offset(1:2),gOffset(1)
    INTEGER :: nGlobalElem
    INTEGER :: firstElem,nLocalElems
    INTEGER :: firstNode,nLocalNodes
    INTEGER :: firstSide,nLocalSides
    INTEGER :: nUniqueSides
    INTEGER :: nGeo,nBCs
    INTEGER :: eid, lsid, iSide
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfInt32_r2) :: hopr_sideInfo
    TYPE(hfReal_r2) :: hopr_nodeCoords
    TYPE(hfInt32_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: bcType

    PRINT*, 'Reading mesh file '//TRIM(meshFile)
    CALL Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId,decomp % mpiComm)

    CALL ReadAttribute_HDF5(fileId,'nElems',nGlobalElem)
    CALL ReadAttribute_HDF5(fileId,'Ngeo',nGeo)
    CALL ReadAttribute_HDF5(fileId,'nBCs',nBCs)
    CALL ReadAttribute_HDF5(fileId,'nUniqueSides',nUniqueSides)

    PRINT*, 'Number of Global Elements : ', nGlobalElem
    PRINT*, 'Input polynomial degree : ', nGeo
    PRINT*, 'Number of BC types : ', nBCs

    ! Read BCType
    CALL bcType % Alloc(loBound=(/1,1/), &
                        upBound=(/4,nBCs/))
    offset(:) = 0
    CALL ReadArray_HDF5(fileId,'BCType',bcType,offset)

    ! Read local subarray of ElemInfo
    CALL decomp % SetElemToRank(nGlobalElem)
    firstElem = decomp % offsetElem % hostData(decomp % rankId) + 1
    nLocalElems = decomp % offsetElem % hostData(decomp % rankId+1)-& 
                  decomp % offsetElem % hostData(decomp % rankId) 

    PRINT*, 'Number of Local Elements : ', nLocalElems

    ! Allocate Space for hopr_elemInfo!
    CALL hopr_elemInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/6,nLocalElems/))

    offset = (/0,firstElem - 1/)
    CALL ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo,offset)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode = hopr_elemInfo % hostData(5,1) + 1
    nLocalNodes = hopr_elemInfo % hostData(6,nLocalElems) - hopr_elemInfo % hostData(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !
    CALL hopr_nodeCoords % Alloc(loBound=(/1,1/), &
                                 upBound=(/3,nLocalNodes/))

    CALL hopr_globalNodeIDs % Alloc(loBound=1, &
                                    upBound=nLocalNodes)

    offset = (/0,firstNode - 1/)
    CALL ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords,offset)
    gOffset = (/firstNode - 1/)
    CALL ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs,gOffset)

    ! Read local subarray of SideInfo
    firstSide = hopr_elemInfo % hostData(3,1) + 1
    nLocalSides = hopr_elemInfo % hostData(4,nLocalElems) - hopr_elemInfo % hostData(3,1)
    PRINT*, 'Number of local sides : ', nLocalSides

    ! Allocate space for hopr_sideInfo
    CALL hopr_sideInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/5,nLocalSides/))

    offset = (/0,firstSide - 1/)
    CALL ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo,offset)

    CALL Close_HDF5(fileID)

    CALL myMesh % Init(nGeo,nLocalElems,nLocalSides,nLocalNodes,nBCs)

    ! Copy data from local arrays into myMesh
    myMesh % hopr_elemInfo % hostData = hopr_elemInfo % hostData
    myMesh % hopr_nodeCoords % hostData = hopr_nodeCoords % hostData
    myMesh % hopr_globalNodeIDs % hostData = hopr_globalNodeIDs % hostData
    myMesh % hopr_sideInfo % hostData = hopr_sideInfo % hostData
    myMesh % nUniqueSides = nUniqueSides
    myMesh % nGlobalElem = nGlobalElem

    iSide = 0 
    DO eid = 1,myMesh % nElem
      DO lsid = 1,6
        iSide = iSide + 1
        myMesh % self_sideInfo % hostData(1:5,lsid,eid) = myMesh % hopr_sideInfo % hostData(1:5,iSide)
      ENDDO
    ENDDO

    CALL myMesh % RecalculateFlip(decomp)

    CALL myMesh % UpdateDevice()

    CALL hopr_elemInfo % Free()
    CALL hopr_nodeCoords % Free()
    CALL hopr_globalNodeIDs % Free()
    CALL hopr_sideInfo % Free()

  END SUBROUTINE Read_HOPr_Mesh3D_parallel

  SUBROUTINE Write_HOPr_Mesh3D(myMesh,meshFile)
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
    CALL WriteArray_HDF5(fileId,'ElemInfo',myMesh % hopr_elemInfo)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    CALL WriteArray_HDF5(fileId,'NodeCoords',myMesh % hopr_nodeCoords)
    CALL WriteArray_HDF5(fileId,'GlobalNodeIDs',myMesh % hopr_globalNodeIDs)

    ! Read local subarray of SideInfo
    CALL WriteArray_HDF5(fileId,'SideInfo',myMesh % hopr_sideInfo)

    CALL Close_HDF5(fileID)

  END SUBROUTINE Write_HOPr_Mesh3D

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

END MODULE SELF_Mesh
