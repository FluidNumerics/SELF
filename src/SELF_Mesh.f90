!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_Mesh

  USE SELF_Constants
  USE SELF_Lagrange
  USE SELF_MPI
  USE SELF_Data
  USE SELF_SupportRoutines
  USE SELF_HDF5
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
!  * For line segments, quads, and hexes, Gauss-Lobatto quadrature is required
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
  INTEGER,PARAMETER :: self_nBCsDefault = 5

  TYPE MeshSpec
    CHARACTER(self_FileNameLength) :: hoprFile

    LOGICAL :: blockMesh
    INTEGER :: blockMesh_nGeo
    INTEGER :: blockMesh_nElemX
    INTEGER :: blockMesh_nElemY
    INTEGER :: blockMesh_nElemZ
    REAL(prec) :: blockMesh_x0, blockMesh_x1
    REAL(prec) :: blockMesh_y0, blockMesh_y1
    REAL(prec) :: blockMesh_z0, blockMesh_z1

  END TYPE MeshSpec

  TYPE,PUBLIC :: Mesh1D
    INTEGER :: nGeo
    INTEGER :: nElem
    INTEGER :: nNodes
    INTEGER :: nCornerNodes
    INTEGER :: nUniqueNodes
    INTEGER :: nBCs
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfReal_r1) :: hopr_nodeCoords
    TYPE(hfReal_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS
    PROCEDURE,PUBLIC :: Init => Init_Mesh1D
    PROCEDURE,PUBLIC :: Free => Free_Mesh1D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Mesh1D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh1D
    PROCEDURE,PUBLIC :: UniformBlockMesh => UniformBlockMesh_Mesh1D

    PROCEDURE,PUBLIC :: Read_HOPr => Read_HOPr_Mesh1D
    PROCEDURE,PUBLIC :: Write_HOPr => Write_HOPr_Mesh1D

  END TYPE Mesh1D

  ! Mesh format is set up as the HOPr format
  ! See https://hopr-project.org/externals/MeshFormat.pdf

  TYPE,PUBLIC :: Mesh2D
    INTEGER :: nGeo
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
    TYPE(hfReal_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r1) :: hopr_CGNSCornerMap
    TYPE(hfInt32_r2) :: hopr_CGNSSideMap
    TYPE(hfInt32_r2) :: hopr_curveNodeMap
    TYPE(hfInt32_r2) :: hopr_hopr_curveNodeMapInv
    TYPE(hfInt32_r2) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS
    PROCEDURE,PUBLIC :: Init => Init_Mesh2D
    PROCEDURE,PUBLIC :: Free => Free_Mesh2D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Mesh2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh2D
    PROCEDURE,PUBLIC :: UniformBlockMesh => UniformBlockMesh_Mesh2D

    PROCEDURE,PUBLIC :: Load => Load_Mesh2D
    PROCEDURE,PUBLIC :: Read_HOPr => Read_HOPr_Mesh2D
    PROCEDURE,PUBLIC :: Write_HOPr => Write_HOPr_Mesh2D

    PROCEDURE,PUBLIC :: Read_ISMv2 => Read_ISMv2_Mesh2D

  END TYPE Mesh2D

  TYPE,PUBLIC :: Mesh3D
    INTEGER :: nGeo
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
    TYPE(hfReal_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r1) :: hopr_CGNSCornerMap
    TYPE(hfInt32_r2) :: hopr_CGNSSideMap
    TYPE(hfInt32_r2) :: hopr_curveNodeMap
    TYPE(hfInt32_r3) :: hopr_hopr_curveNodeMapInv
    TYPE(hfInt32_r2) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Mesh3D
    PROCEDURE,PUBLIC :: Free => Free_Mesh3D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Mesh3D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh3D
    PROCEDURE,PUBLIC :: UniformBlockMesh => UniformBlockMesh_Mesh3D

    PROCEDURE,PUBLIC :: Load => Load_Mesh3D
    PROCEDURE,PUBLIC :: Read_HOPr => Read_HOPr_Mesh3D
    PROCEDURE,PUBLIC :: Write_HOPr => Write_HOPr_Mesh3D

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

#ifdef GPU
    CALL myMesh % hopr_elemInfo % UpdateHost()
    CALL myMesh % hopr_nodeCoords % UpdateHost()
    CALL myMesh % hopr_globalNodeIDs % UpdateHost()
    CALL myMesh % BCType % UpdateHost()
#endif

  END SUBROUTINE UpdateHost_Mesh1D

  SUBROUTINE UpdateDevice_Mesh1D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(inout) :: myMesh

#ifdef GPU
    CALL myMesh % hopr_elemInfo % UpdateDevice()
    CALL myMesh % hopr_nodeCoords % UpdateDevice()
    CALL myMesh % hopr_globalNodeIDs % UpdateDevice()
    CALL myMesh % BCType % UpdateDevice()
#endif

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

  SUBROUTINE Read_HOPr_Mesh1D( myMesh, meshFile, nRanks, myRank, mpiComm )
  ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
  ! Adapted for 1D Mesh : Note that HOPR does not have 1D mesh output.
    IMPLICIT NONE
    CLASS(Mesh1D), INTENT(out) :: myMesh
    CHARACTER(*), INTENT(in) :: meshFile
    INTEGER, INTENT(in) :: nRanks
    INTEGER, INTENT(in) :: myRank
    INTEGER, OPTIONAL, INTENT(in) :: mpiComm
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: offset(1:2), gOffset(1) 
    INTEGER :: nGlobalElem
    INTEGER :: offSetElem(0:nRanks)
    INTEGER :: firstElem, nLocalElems
    INTEGER :: firstNode, nLocalNodes
    INTEGER :: nGeo, nBCs
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfReal_r1) :: hopr_nodeCoords
    TYPE(hfReal_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: bcType

    IF( PRESENT(mpiComm) )THEN
      CALL Open_HDF5(meshFile, H5F_ACC_RDWR_F, fileId, mpiComm)
    ELSE
      CALL Open_HDF5(meshFile, H5F_ACC_RDWR_F, fileId)
    ENDIF

    CALL ReadAttribute_HDF5(fileId, 'nElems', nGlobalElem)
    CALL ReadAttribute_HDF5(fileId, 'Ngeo', nGeo)
    CALL ReadAttribute_HDF5(fileId, 'nBCs', nBCs)

    CALL DomainDecomp(nGlobalElem,nRanks,offSetElem)

    ! Read BCType
    CALL bcType % Alloc(loBound=(/1,1/), &
                        upBound=(/4,nBCs/))
    offset(:) = 0
    CALL ReadArray_HDF5(fileId, 'BCType', offset, bcType)

    ! Read local subarray of ElemInfo
    firstElem = offsetElem(myRank)+1

    nLocalElems = offsetElem(myRank+1)-offsetElem(myRank)

    ! Allocate Space for hopr_elemInfo!
    CALL hopr_elemInfo % Alloc(loBound=(/1,1/), &
                          upBound=(/6,nLocalElems/))

    offset = (/0, firstElem-1/)
    CALL ReadArray_HDF5(fileId, 'ElemInfo', offset, hopr_elemInfo  )

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode= hopr_elemInfo % hostData(5,1)+1

    nLocalNodes = hopr_elemInfo % hostData(6,nLocalElems)-hopr_elemInfo % hostData(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !

    CALL hopr_nodeCoords % Alloc(loBound=1, &
                            upBound=nLocalNodes)

    CALL hopr_globalNodeIDs % Alloc(loBound=1, &
                               upBound=nLocalNodes)

    gOffset = (/firstNode-1/)
    CALL ReadArray_HDF5(fileId, 'NodeCoords', gOffset, hopr_nodeCoords  )
    CALL ReadArray_HDF5(fileId, 'GlobalNodeIDs', gOffset, hopr_globalNodeIDs  )

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

  END SUBROUTINE Read_HOPr_Mesh1D

  SUBROUTINE Write_HOPr_Mesh1D( myMesh, meshFile )
  ! Writes mesh output in HOPR format
    IMPLICIT NONE
    CLASS(Mesh1D), INTENT(inout) :: myMesh
    CHARACTER(*), INTENT(in) :: meshFile
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: offset(1:2), gOffset(1) 
    INTEGER :: nGlobalElem
    INTEGER :: firstElem, nLocalElems
    INTEGER :: firstNode, nLocalNodes
    INTEGER :: firstSide, nLocalSides
    INTEGER :: nGeo, nBCs

 
    CALL Open_HDF5(meshFile, H5F_ACC_RDWR_F, fileId)
 
    CALL WriteAttribute_HDF5(fileId, 'nElems', myMesh % nElem)
    CALL WriteAttribute_HDF5(fileId, 'Ngeo', myMesh % nGeo )
    CALL WriteAttribute_HDF5(fileId, 'nBCs', myMesh % nBCs)

    offset(:) = 0
    CALL WriteArray_HDF5(fileId, 'BCType', offset, myMesh % bcType)

    ! Read local subarray of ElemInfo
    firstElem = 1

    nLocalElems = myMesh % nElem

    offset = (/0, firstElem-1/)
    CALL WriteArray_HDF5(fileId, 'ElemInfo', offset, myMesh % hopr_elemInfo  )

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode= myMesh % hopr_elemInfo % hostData(5,1)+1

    nLocalNodes = myMesh % hopr_elemInfo % hostData(6,nLocalElems)-myMesh % hopr_elemInfo % hostData(5,1)

    offset = (/0, firstNode-1/)
    CALL WriteArray_HDF5(fileId, 'NodeCoords', offset, myMesh % hopr_nodeCoords  )
    gOffset = (/firstNode-1/)
    CALL WriteArray_HDF5(fileId, 'GlobalNodeIDs', gOffset, myMesh % hopr_globalNodeIDs  )

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
                                     upBound=(/2,(nGeo+1)**2,nElem/))

    CALL myMesh % hopr_globalNodeIDs % Alloc(loBound=1, &
                                        upBound=nNodes)

    CALL myMesh % hopr_CGNSCornerMap % Alloc(loBound=1, &
                                        upBound=4)

    CALL myMesh % hopr_CGNSSideMap % Alloc(loBound=(/1,1/), &
                                      upBound=(/2,4/))

    CALL myMesh % hopr_curveNodeMap % Alloc(loBound=(/1,1/),&
                                       upBound=(/2,(nGeo+1)**2/))

    CALL myMesh % hopr_hopr_curveNodeMapInv % Alloc(loBound=(/0,0/), &
                                          upBound=(/nGeo,nGeo/))

    CALL myMesh % BCType % Alloc(loBound=(/1,1/), &
                                 upBound=(/4,nBCs/))

    ALLOCATE (myMesh % BCNames(1:nBCs))

    ! Create lookup tables to assist with connectivity generation
    myMesh % hopr_CGNSCornerMap % hostData(1) = 1
    myMesh % hopr_CGNSCornerMap % hostData(2) = nGeo+1
    myMesh % hopr_CGNSCornerMap % hostData(3) = (nGeo+1)**2
    myMesh % hopr_CGNSCornerMap % hostData(4) = nGeo*(nGeo+1)+1

    DO j = 0, nGeo
      DO i = 0, nGeo
        l = l+1
        myMesh % hopr_curveNodeMap % hostData(1:2,l) = (/i,j/)
        myMesh % hopr_hopr_curveNodeMapInv % hostData(i,j) = l
      ENDDO
    ENDDO

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

#ifdef GPU
    CALL myMesh % hopr_elemInfo % UpdateHost()
    CALL myMesh % hopr_sideInfo % UpdateHost()
    CALL myMesh % self_sideInfo % UpdateHost()
    CALL myMesh % hopr_nodeCoords % UpdateHost()
    CALL myMesh % self_nodeCoords % UpdateHost()
    CALL myMesh % hopr_globalNodeIDs % UpdateHost()
    CALL myMesh % BCType % UpdateHost()
#endif

  END SUBROUTINE UpdateHost_Mesh2D

  SUBROUTINE UpdateDevice_Mesh2D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: myMesh

#ifdef GPU
    CALL myMesh % hopr_elemInfo % UpdateDevice()
    CALL myMesh % hopr_sideInfo % UpdateDevice()
    CALL myMesh % self_sideInfo % UpdateDevice() 
    CALL myMesh % hopr_nodeCoords % UpdateDevice()
    CALL myMesh % self_nodeCoords % UpdateDevice()
    CALL myMesh % hopr_globalNodeIDs % UpdateDevice()
    CALL myMesh % BCType % UpdateDevice()
#endif

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
    nNodes = nel*(nGeo + 1)*(nGeo + 1)
    nSides = (nElem(1) + 1)*nElem(2) + (nElem(2) + 1)*nElem(1)
    CALL myMesh % Init(nGeo,nEl,nSides,nNodes,1)

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
        myMesh % hopr_elemInfo % hostData(1,iel) = selfQuadLinear ! Element Type
        myMesh % hopr_elemInfo % hostData(2,iel) = 1 ! Element Zone
        myMesh % hopr_elemInfo % hostData(3,iel) = sid ! Side Index Start
        sid = sid+4
        myMesh % hopr_elemInfo % hostData(4,iel) = sid ! Side Index End
        myMesh % hopr_elemInfo % hostData(5,iel) = nid-1 ! Node Index Start
        DO j = 0,nGeo
          DO i = 0,nGeo
            myMesh % hopr_nodeCoords % hostData(1:2,nid) = xGeo % interior % hostData(1:2,i,j,1,elid)
            nid = nid + 1
          END DO
        END DO
        myMesh % hopr_elemInfo % hostData(6,iel) = nid ! Node Index End
        elid = elid + 1
      END DO
    END DO

    ! Set the SideInfo


    CALL myMesh % UpdateDevice()

    CALL xLinear % Free()
    CALL xGeo % Free()

  END SUBROUTINE UniformBlockMesh_Mesh2D
  
  SUBROUTINE Load_Mesh2D(myMesh,myMeshSpec,decomp)
#undef __FUNC__
#define __FUNC__ "Load_Mesh2D"
    IMPLICIT NONE
    CLASS(Mesh2D), INTENT(out) :: myMesh
    TYPE(MeshSpec), INTENT(in) :: myMeshSpec
    TYPE(MPILayer), INTENT(inout) :: decomp

      IF(myMeshSpec % blockMesh)THEN

        IF(decomp % nRanks > 1)THEN ! Error out
          ERROR("Block Mesh only supported in serial")
          STOP ! TO DO : Safe exit for serial and parallel
        ELSE
          CALL myMesh % UniformBlockMesh(myMeshSpec % blockMesh_nGeo,&
                (/myMeshSpec % blockMesh_nElemX,myMeshSpec % blockMesh_nElemY/),&
                (/myMeshSpec % blockMesh_x0,myMeshSpec % blockMesh_x1,&
                  myMeshSpec % blockMesh_y0,myMeshSpec % blockMesh_y1/))
        ENDIF

      ELSE

        IF(decomp % mpiEnabled)THEN
          CALL myMesh % Read_HOPr(myMeshSpec % hoprFile,decomp)
        ELSE
          CALL myMesh % Read_HOPr(myMeshSpec % hoprFile)
        ENDIF

      ENDIF
    
  END SUBROUTINE Load_Mesh2D

  SUBROUTINE Read_ISMv2_Mesh2D( myMesh, meshFile, decomp )
#undef __FUNC__
#define __FUNC__ "Read_ISMv2_Mesh2D"
    IMPLICIT NONE
    CLASS(Mesh2D), INTENT(out) :: myMesh
    CHARACTER(*), INTENT(in) :: meshFile
    TYPE(MPILayer), OPTIONAL, INTENT(inout) :: decomp
    ! Local
    INTEGER :: nNodes, nSides, nElem, nGeo
    INTEGER :: nRanks, myRank
    INTEGER :: unid, lnid, nid, usid, lsid, sid, eid
    INTEGER :: fUnit
    INTEGER :: i, j
    CHARACTER(100) :: line
    CHARACTER(500) :: msg
    REAL(prec) :: x(1:3), x0(1:2), x1(1:2)
    INTEGER :: bCurveFlag(1:4)
    TYPE(Lagrange) :: interp

      OPEN(UNIT=NEWUNIT(fUnit),&
           FILE=TRIM(meshFile),&
           FORM='FORMATTED',&
           STATUS='OLD',&
           ACCESS='SEQUENTIAL')

      READ(fUnit,*) line
      IF( TRIM(line) /= ' ISM-v2' )THEN
        msg='Unrecognized file format : '//TRIM(line)
        ERROR(msg)
        STOP
      ENDIF

      READ(fUnit,*) nNodes, nSides, nElem, nGeo

      ! TO DO : Need to check if HOHQMesh places boundary curves on Gauss Lobatto points or uniform points
      CALL interp % Init(nGeo,GAUSS_LOBATTO,nGeo,GAUSS_LOBATTO)

      ! When we initialize the mesh, we set nNodes=nElem*4*(nGeo+1)**2 and
      ! nSides = nElem*4 since we still use `nNodes` and `nSides`
      ! in the input correspond to the HOPR definitions of these
      ! variables. 
      ! `nSides` in HOHQMesh corresponds to nUniqueSides in HOPR and SELF 
      ! `nNodes` in HOHQMesh corresponds to nCornerNodes (unique) in HOPR and SELF
      CALL myMesh % Init(nGeo,nElem,nElem*4,nElem*4*(nGeo+1)**2,self_nBCsDefault)
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

      DO nid = 1, myMesh % nCornerNodes
        READ(fUnit,*) x
        myMesh % hohq_cornerNodes % hostData(1:2,nid) = x(1:2)
      ENDDO

      DO usid = 1, myMesh % nUniqueSides
        READ(fUnit,*) myMesh % hohq_sideInfo % hostData(1:6,usid)
      ENDDO

      DO eid = 1, myMesh % nElem
        READ(fUnit,*) myMesh % hohq_elemInfo % hostData(1:4,eid)
        READ(fUnit,*) bCurveFlag(1:4)
        DO lsid = 1, 4
          IF(bCurveFlag(lsid) == 1)THEN
            DO i = 0, nGeo
              READ(fUnit,*) x
              myMesh % hohq_sideCurves % hostData(1:2,i,lsid,eid) = x(1:2)
            ENDDO
          ELSE

            ! For non-polynomial sides, create the side curve through interpolation between corner nodes
            lnid = myMesh % hopr_CGNSSideMap % hostData(1,lsid)
            nid = myMesh % hohq_elemInfo % hostData(lnid,eid)
            x0(1:2) = myMesh % hohq_cornerNodes % hostData(1:2,nid)

            lnid = myMesh % hopr_CGNSSideMap % hostData(2,lsid)
            nid = myMesh % hohq_elemInfo % hostData(lnid,eid)
            x1(1:2) = myMesh % hohq_cornerNodes % hostData(1:2,nid)

            DO i = 0, nGeo
              myMesh % hohq_sideCurves % hostData(1:2,i,lsid,eid) = 0.5_prec*(&
                       x0(1:2)*(1.0_prec-interp % controlPoints % hostData(i))+&
                       x1(1:2)*(interp % controlPoints % hostData(i)+1.0_prec))
            ENDDO

          ENDIF
        ENDDO
        READ(fUnit,*) line
        ! TO DO : Parse line for boundary conditions
      ENDDO

      CLOSE(fUnit)

      ! Generate the self_nodeCoords through transfinite interpolation with linear blending
      !DO eid = 1, myMesh % nElem


      !ENDDO

      CALL interp % Free()

      ! TO DO : Fill in self_* attributes from hohq_ attributes 
      !   hohq_sideInfo(1:6)
      !       1. start node ID (corner node ID)
      !       2. end node ID (corner node ID)
      !       3. primary element ID
      !       4. secondary element ID
      !       5. local side ID for primary element
      !       6. local side ID for secondary element
      !
      !   hohq_elemInfo(1:4)
      !       > Global Corner Node IDs, locally ordered by CGNS mapping
      !
      !   To create self_* attributes, we need to 
      !       1. Use hohq_sideCurves to generate self_nodeCoords (transfinite interpolation with linear blending)
      !       2. Copy self_nodeCoords to hopr_nodeCoords (array reshape)
      !       3. Find unique hopr_nodeCoords to set hopr_globalNodeIDs
      !       4. Use hopr_globalNodeIDs and hopr_CGNSSideMap to set hopr_sideInfo (uniqe side identifcation and connectivity
      !       mapping)
      !       
      
  END SUBROUTINE Read_ISMv2_Mesh2D

  SUBROUTINE Read_HOPr_Mesh2D( myMesh, meshFile, decomp )
  ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
  ! Adapted for 2D Mesh : Note that HOPR does not have 2D mesh output.
    IMPLICIT NONE
    CLASS(Mesh2D), INTENT(out) :: myMesh
    CHARACTER(*), INTENT(in) :: meshFile
    TYPE(MPILayer), OPTIONAL, INTENT(inout) :: decomp
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: offset(1:2), gOffset(1) 
    INTEGER :: nGlobalElem
    INTEGER :: firstElem, nLocalElems
    INTEGER :: firstNode, nLocalNodes
    INTEGER :: firstSide, nLocalSides
    INTEGER :: nGeo, nBCs
    INTEGER :: nRanks, myRank
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfInt32_r2) :: hopr_sideInfo
    TYPE(hfReal_r2) :: hopr_nodeCoords
    TYPE(hfReal_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: bcType

    IF( PRESENT(decomp) )THEN
      CALL Open_HDF5(meshFile, H5F_ACC_RDWR_F, fileId, decomp % mpiComm)
      nRanks = decomp % nRanks
      myRank = decomp % rankId
    ELSE
      CALL Open_HDF5(meshFile, H5F_ACC_RDWR_F, fileId)
      nRanks = 1
      myRank = 0
    ENDIF

    CALL ReadAttribute_HDF5(fileId, 'nElems', nGlobalElem)
    CALL ReadAttribute_HDF5(fileId, 'Ngeo', nGeo)
    CALL ReadAttribute_HDF5(fileId, 'nBCs', nBCs)

    ! Read BCType
    CALL bcType % Alloc(loBound=(/1,1/), &
                        upBound=(/4,nBCs/))
    offset(:) = 0
    CALL ReadArray_HDF5(fileId, 'BCType', offset, bcType)

    ! Read local subarray of ElemInfo
    IF( PRESENT(decomp) )THEN
      CALL decomp % SetElemToRank(nGlobalElem)
      firstElem = decomp % offsetElem % hostData(myRank)+1
      nLocalElems = decomp % offsetElem % hostData(myRank+1)-decomp % offsetElem % hostData(myRank)
    ELSE
      firstElem = 1
      nLocalElems = nGlobalElem
    ENDIF

    ! Allocate Space for hopr_elemInfo!
    CALL hopr_elemInfo % Alloc(loBound=(/1,1/), &
                          upBound=(/6,nLocalElems/))

    offset = (/0, firstElem-1/)
    CALL ReadArray_HDF5(fileId, 'ElemInfo', offset, hopr_elemInfo  )

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode= hopr_elemInfo % hostData(5,1)+1

    nLocalNodes = hopr_elemInfo % hostData(6,nLocalElems)-hopr_elemInfo % hostData(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !

    CALL hopr_nodeCoords % Alloc(loBound=(/1,1/), &
                            upBound=(/2,nLocalNodes/))

    CALL hopr_globalNodeIDs % Alloc(loBound=1, &
                               upBound=nLocalNodes)

    offset = (/0, firstNode-1/)
    CALL ReadArray_HDF5(fileId, 'NodeCoords', offset, hopr_nodeCoords  )
    gOffset = (/firstNode-1/)
    CALL ReadArray_HDF5(fileId, 'GlobalNodeIDs', gOffset, hopr_globalNodeIDs  )

    ! Read local subarray of SideInfo
    firstSide= hopr_elemInfo % hostData(3,1)+1

    nLocalSides = hopr_elemInfo % hostData(4,nLocalElems)-hopr_elemInfo % hostData(3,1)

    ! Allocate space for hopr_sideInfo
    CALL hopr_sideInfo % Alloc(loBound=(/1,1/), &
                               upBound=(/5,nLocalSides/))

    offset = (/0, firstSide-1/)
    CALL ReadArray_HDF5(fileId, 'SideInfo', offset, hopr_sideInfo  )

    CALL Close_HDF5(fileID)

    CALL myMesh % Init(nGeo,nLocalElems,nLocalSides,nLocalNodes,nBCs)

    ! Copy data from local arrays into myMesh
    myMesh % hopr_elemInfo = hopr_elemInfo
    myMesh % hopr_nodeCoords = hopr_nodeCoords
    myMesh % hopr_globalNodeIDs = hopr_globalNodeIDs
    myMesh % hopr_sideInfo = hopr_sideInfo

    ! TO DO : Set SELF Side Info from HOPR side info

    CALL myMesh % UpdateDevice()

    CALL hopr_elemInfo % Free()
    CALL hopr_nodeCoords % Free()
    CALL hopr_globalNodeIDs % Free()
    CALL hopr_sideInfo % Free()

  END SUBROUTINE Read_HOPr_Mesh2D

  SUBROUTINE Write_HOPr_Mesh2D( myMesh, meshFile)
  ! Writes mesh output in HOPR format
    IMPLICIT NONE
    CLASS(Mesh2D), INTENT(inout) :: myMesh
    CHARACTER(*), INTENT(in) :: meshFile
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: offset(1:2), gOffset(1) 
    INTEGER :: nGlobalElem
    INTEGER :: firstElem, nLocalElems
    INTEGER :: firstNode, nLocalNodes
    INTEGER :: firstSide, nLocalSides
    INTEGER :: nGeo, nBCs

    CALL Open_HDF5(meshFile, H5F_ACC_RDWR_F, fileId)
    CALL WriteAttribute_HDF5(fileId, 'nElems', myMesh % nElem)
    CALL WriteAttribute_HDF5(fileId, 'Ngeo', myMesh % nGeo)
    CALL WriteAttribute_HDF5(fileId, 'nBCs', myMesh % nBCs)

    offset(:) = 0
    CALL WriteArray_HDF5(fileId, 'BCType', offset, myMesh % bcType)

    ! Write local subarray of ElemInfo
    firstElem = 1
    nLocalElems = myMesh % nElem

    offset = (/0, firstElem-1/)
    CALL WriteArray_HDF5(fileId, 'ElemInfo', offset, myMesh % hopr_elemInfo  )

    ! Write local subarray of NodeCoords and GlobalNodeIDs
    firstNode = myMesh % hopr_elemInfo % hostData(5,1)+1

    nLocalNodes = myMesh % hopr_elemInfo % hostData(6,nLocalElems)-myMesh % hopr_elemInfo % hostData(5,1)

    offset = (/0, firstNode-1/)
    CALL WriteArray_HDF5(fileId, 'NodeCoords', offset, myMesh % hopr_nodeCoords  )
    gOffset = (/firstNode-1/)
    CALL WriteArray_HDF5(fileId, 'GlobalNodeIDs', gOffset, myMesh % hopr_globalNodeIDs  )

    ! Read local subarray of SideInfo
    firstSide= myMesh % hopr_elemInfo % hostData(3,1)+1

    nLocalSides = myMesh % hopr_elemInfo % hostData(4,nLocalElems)-myMesh % hopr_elemInfo % hostData(3,1)

    offset = (/0, firstSide-1/)
    CALL WriteArray_HDF5(fileId, 'SideInfo', offset, myMesh % hopr_sideInfo  )

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
                                     upBound=(/3,(nGeo+1)**3,nElem/))

    CALL myMesh % hopr_globalNodeIDs % Alloc(loBound=1, &
                                        upBound=nNodes)

    CALL myMesh % hopr_CGNSCornerMap % Alloc(loBound=1, &
                                        upBound=8)

    CALL myMesh % hopr_CGNSSideMap % Alloc(loBound=(/1,1/), &
                                      upBound=(/4,6/))

    CALL myMesh % hopr_curveNodeMap % Alloc(loBound=(/1,1/), &
                                       upBound=(/3,(nGeo+1)**3/))

    CALL myMesh % hopr_hopr_curveNodeMapInv % Alloc(loBound=(/0,0,0/), &
                                          upBound=(/nGeo,nGeo,nGeo/))

    CALL myMesh % BCType % Alloc(loBound=(/1,1/), &
                                 upBound=(/4,nBCs/))

    ALLOCATE (myMesh % BCNames(1:nBCs))

    ! Create lookup tables to assist with connectivity generation
    myMesh % hopr_CGNSCornerMap % hostData(1) = 1
    myMesh % hopr_CGNSCornerMap % hostData(2) = nGeo+1
    myMesh % hopr_CGNSCornerMap % hostData(3) = (nGeo+1)**2
    myMesh % hopr_CGNSCornerMap % hostData(4) = nGeo*(nGeo+1)+1
    myMesh % hopr_CGNSCornerMap % hostData(5) = nGeo*(nGeo+1)**2+1
    myMesh % hopr_CGNSCornerMap % hostData(6) = nGeo*(nGeo+1)**2+(nGeo+1)
    myMesh % hopr_CGNSCornerMap % hostData(7) = (nGeo+1)**3
    myMesh % hopr_CGNSCornerMap % hostData(8) = nGeo*(nGeo+1)*(nGeo+2)+1

    DO k = 0, nGeo
      DO j = 0, nGeo
        DO i = 0, nGeo
          l = l+1
          myMesh % hopr_curveNodeMap % hostData(1:3,l) = (/i,j,k/)
          myMesh % hopr_hopr_curveNodeMapInv % hostData(i,j,k) = l
        ENDDO
      ENDDO
    ENDDO

    ! Maps from local corner node id to CGNS side
    myMesh % hopr_CGNSSideMap % hostData(1:4,1) = (/1,4,3,2/)
    myMesh % hopr_CGNSSideMap % hostData(1:4,2) = (/1,2,6,5/)
    myMesh % hopr_CGNSSideMap % hostData(1:4,3) = (/2,3,7,6/)
    myMesh % hopr_CGNSSideMap % hostData(1:4,4) = (/3,4,8,7/)
    myMesh % hopr_CGNSSideMap % hostData(1:4,5) = (/1,5,8,4/)
    myMesh % hopr_CGNSSideMap % hostData(1:4,6) = (/5,6,7,8/)

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
    CALL myMesh % BCType % Free()

    DEALLOCATE (myMesh % BCNames)

  END SUBROUTINE Free_Mesh3D
  SUBROUTINE UpdateHost_Mesh3D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh

#ifdef GPU
    CALL myMesh % hopr_elemInfo % UpdateHost()
    CALL myMesh % hopr_sideInfo % UpdateHost()
    CALL myMesh % self_sideInfo % UpdateHost()
    CALL myMesh % hopr_nodeCoords % UpdateHost()
    CALL myMesh % self_nodeCoords % UpdateHost()
    CALL myMesh % hopr_globalNodeIDs % UpdateHost()
    CALL myMesh % BCType % UpdateHost()
#endif

  END SUBROUTINE UpdateHost_Mesh3D

  SUBROUTINE UpdateDevice_Mesh3D(myMesh)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: myMesh

#ifdef GPU
    CALL myMesh % hopr_elemInfo % UpdateDevice()
    CALL myMesh % hopr_sideInfo % UpdateDevice()
    CALL myMesh % self_sideInfo % UpdateDevice()
    CALL myMesh % hopr_nodeCoords % UpdateDevice()
    CALL myMesh % self_nodeCoords % UpdateDevice()
    CALL myMesh % hopr_globalNodeIDs % UpdateDevice()
    CALL myMesh % BCType % UpdateDevice()
#endif

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
    INTEGER :: i,j,k
    REAL(prec) :: xU(1:nElem(1) + 1)
    REAL(prec) :: yU(1:nElem(2) + 1)
    REAL(prec) :: zU(1:nElem(3) + 1)
    TYPE(Vector3D) :: xLinear
    TYPE(Vector3D) :: xGeo

    nEl = nElem(1)*nElem(2)*nElem(3)
    nNodes = nel*(nGeo + 1)*(nGeo + 1)*(nGeo + 1)
    nSides = (nElem(1) + 1)*nElem(2)*nElem(3) + &
             (nElem(2) + 1)*nElem(1)*nElem(3) + &
             (nElem(3) + 1)*nElem(1)*nElem(2)
    CALL myMesh % Init(nGeo,nEl,nSides,nNodes,1)

    ! Set the hopr_nodeCoords
    xU = UniformPoints(x(1),x(2),1,nElem(1) + 1)
    yU = UniformPoints(x(3),x(4),1,nElem(2) + 1)
    zU = UniformPoints(x(5),x(6),1,nElem(3) + 1)

    ! Create a linear interpolant to interpolate to nGeo grid
    CALL xLinear % Init(1,GAUSS_LOBATTO, &
                        nGeo,GAUSS_LOBATTO, &
                        1,nEl)

    CALL xGeo % Init(nGeo,GAUSS_LOBATTO, &
                     nGeo,GAUSS_LOBATTO, &
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
          myMesh % hopr_elemInfo % hostData(1,iel) = selfQuadLinear ! Element Type
          myMesh % hopr_elemInfo % hostData(2,iel) = 1 ! Element Zone
          myMesh % hopr_elemInfo % hostData(3,iel) = sid ! Side Index Start
          sid = sid + 6
          myMesh % hopr_elemInfo % hostData(4,iel) = sid ! Side Index End
          myMesh % hopr_elemInfo % hostData(5,iel) = nid-1 ! Node Index Start
          DO k = 0,nGeo
            DO j = 0,nGeo
              DO i = 0,nGeo
                myMesh % hopr_nodeCoords % hostData(1:3,nid) = xGeo % interior % hostData(1:3,i,j,k,1,elid)
                nid = nid + 1
              END DO
            END DO
          END DO
          myMesh % hopr_elemInfo % hostData(6,iel) = nid ! Node Index End
          elid = elid + 1
        END DO
      END DO
    END DO

    ! TO DO: Add Side information !

    CALL myMesh % UpdateDevice()

    CALL xLinear % Free()
    CALL xGeo % Free()

  END SUBROUTINE UniformBlockMesh_Mesh3D

  SUBROUTINE Load_Mesh3D(myMesh,myMeshSpec,decomp)
#undef __FUNC__
#define __FUNC__ "Load_Mesh3D"
    IMPLICIT NONE
    CLASS(Mesh3D), INTENT(out) :: myMesh
    TYPE(MeshSpec), INTENT(in) :: myMeshSpec
    TYPE(MPILayer), INTENT(inout) :: decomp

      IF(myMeshSpec % blockMesh)THEN

        IF(decomp % nRanks > 1)THEN ! Error out
          ERROR("Block Mesh only supported in serial")
          STOP ! TO DO : Safe exit for serial and parallel
        ELSE
          CALL myMesh % UniformBlockMesh(myMeshSpec % blockMesh_nGeo,&
                (/myMeshSpec % blockMesh_nElemX,myMeshSpec % blockMesh_nElemY,myMeshSpec % blockMesh_nElemZ/),&
                (/myMeshSpec % blockMesh_x0,myMeshSpec % blockMesh_x1,&
                  myMeshSpec % blockMesh_y0,myMeshSpec % blockMesh_y1,&
                  myMeshSpec % blockMesh_z0,myMeshSpec % blockMesh_z1/))
        ENDIF

      ELSE

        IF(decomp % mpiEnabled)THEN
          CALL myMesh % Read_HOPr(myMeshSpec % hoprFile,decomp)
        ELSE
          CALL myMesh % Read_HOPr(myMeshSpec % hoprFile)
        ENDIF

      ENDIF
    
  END SUBROUTINE Load_Mesh3D

  SUBROUTINE Read_HOPr_Mesh3D( myMesh, meshFile, decomp )
  ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    IMPLICIT NONE
    CLASS(Mesh3D), INTENT(out) :: myMesh
    CHARACTER(*), INTENT(in) :: meshFile
    TYPE(MPILayer), OPTIONAL, INTENT(inout) :: decomp
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: offset(1:2), gOffset(1) 
    INTEGER :: nGlobalElem
    INTEGER :: firstElem, nLocalElems
    INTEGER :: firstNode, nLocalNodes
    INTEGER :: firstSide, nLocalSides
    INTEGER :: nGeo, nBCs
    INTEGER :: nRanks, myRank
    TYPE(hfInt32_r2) :: hopr_elemInfo
    TYPE(hfInt32_r2) :: hopr_sideInfo
    TYPE(hfReal_r2) :: hopr_nodeCoords
    TYPE(hfReal_r1) :: hopr_globalNodeIDs
    TYPE(hfInt32_r2) :: bcType

    IF( PRESENT(decomp) )THEN
      CALL Open_HDF5(meshFile, H5F_ACC_RDWR_F, fileId, decomp % mpiComm)
      nRanks = decomp % nRanks
      myRank = decomp % rankId
    ELSE
      CALL Open_HDF5(meshFile, H5F_ACC_RDWR_F, fileId)
      nRanks = 1
      myRank = 0
    ENDIF

    CALL ReadAttribute_HDF5(fileId, 'nElems', nGlobalElem)
    CALL ReadAttribute_HDF5(fileId, 'Ngeo', nGeo)
    CALL ReadAttribute_HDF5(fileId, 'nBCs', nBCs)

    ! Read BCType
    CALL bcType % Alloc(loBound=(/1,1/), &
                        upBound=(/4,nBCs/))
    offset(:) = 0
    CALL ReadArray_HDF5(fileId, 'BCType', offset, bcType)

    ! Read local subarray of ElemInfo
    IF( PRESENT(decomp) )THEN
      CALL decomp % SetElemToRank(nGlobalElem)
      firstElem = decomp % offsetElem % hostData(myRank)+1
      nLocalElems = decomp % offsetElem % hostData(myRank+1)-decomp % offsetElem % hostData(myRank)
    ELSE
      firstElem = 1
      nLocalElems = nGlobalElem
    ENDIF

    ! Allocate Space for hopr_elemInfo!
    CALL hopr_elemInfo % Alloc(loBound=(/1,1/), &
                          upBound=(/6,nLocalElems/))

    offset = (/0, firstElem-1/)
    CALL ReadArray_HDF5(fileId, 'ElemInfo', offset, hopr_elemInfo  )

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode= hopr_elemInfo % hostData(5,1)+1

    nLocalNodes = hopr_elemInfo % hostData(6,nLocalElems)-hopr_elemInfo % hostData(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !

    CALL hopr_nodeCoords % Alloc(loBound=(/1,1/), &
                            upBound=(/3,nLocalNodes/))

    CALL hopr_globalNodeIDs % Alloc(loBound=1, &
                               upBound=nLocalNodes)

    offset = (/0, firstNode-1/)
    CALL ReadArray_HDF5(fileId, 'NodeCoords', offset, hopr_nodeCoords  )
    gOffset = (/firstNode-1/)
    CALL ReadArray_HDF5(fileId, 'GlobalNodeIDs', gOffset, hopr_globalNodeIDs  )

    ! Read local subarray of SideInfo
    firstSide= hopr_elemInfo % hostData(3,1)+1

    nLocalSides = hopr_elemInfo % hostData(4,nLocalElems)-hopr_elemInfo % hostData(3,1)

    ! Allocate space for hopr_sideInfo
    CALL hopr_sideInfo % Alloc(loBound=(/1,1/), &
                          upBound=(/5,nLocalSides/))

    offset = (/0, firstSide-1/)
    CALL ReadArray_HDF5(fileId, 'SideInfo', offset, hopr_sideInfo  )

    CALL Close_HDF5(fileID)

    CALL myMesh % Init(nGeo,nLocalElems,nLocalSides,nLocalNodes,nBCs)

    ! Copy data from local arrays into myMesh
    myMesh % hopr_elemInfo = hopr_elemInfo
    myMesh % hopr_nodeCoords = hopr_nodeCoords
    myMesh % hopr_globalNodeIDs = hopr_globalNodeIDs
    myMesh % hopr_sideInfo = hopr_sideInfo

    ! TO DO : Copy SELF side info from HOPR side info

    CALL myMesh % UpdateDevice()

    CALL hopr_elemInfo % Free()
    CALL hopr_nodeCoords % Free()
    CALL hopr_globalNodeIDs % Free()
    CALL hopr_sideInfo % Free()

  END SUBROUTINE Read_HOPr_Mesh3D

  SUBROUTINE Write_HOPr_Mesh3D( myMesh, meshFile )
  ! Writes mesh output in HOPR format
    IMPLICIT NONE
    CLASS(Mesh3D), INTENT(inout) :: myMesh
    CHARACTER(*), INTENT(in) :: meshFile
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: offset(1:2), gOffset(1) 
    INTEGER :: nGlobalElem
    INTEGER :: firstElem, nLocalElems
    INTEGER :: firstNode, nLocalNodes
    INTEGER :: firstSide, nLocalSides
    INTEGER :: nGeo, nBCs


    CALL Open_HDF5(meshFile, H5F_ACC_RDWR_F, fileId)

    CALL WriteAttribute_HDF5(fileId, 'nElems', myMesh % nElem)
    CALL WriteAttribute_HDF5(fileId, 'Ngeo', myMesh % nGeo)
    CALL WriteAttribute_HDF5(fileId, 'nBCs', myMesh % nBCs)

    offset(:) = 0
    CALL WriteArray_HDF5(fileId, 'BCType', offset, myMesh % bcType)

    nLocalElems = myMesh % nElem

    offset = (/0, 0/)
    CALL WriteArray_HDF5(fileId, 'ElemInfo', offset, myMesh % hopr_elemInfo  )

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode = myMesh % hopr_elemInfo % hostData(5,1)+1

    nLocalNodes = myMesh % hopr_elemInfo % hostData(6,nLocalElems)-myMesh % hopr_elemInfo % hostData(5,1)

    offset = (/0, firstNode-1/)
    CALL WriteArray_HDF5(fileId, 'NodeCoords', offset, myMesh % hopr_nodeCoords  )
    gOffset = (/firstNode-1/)
    CALL WriteArray_HDF5(fileId, 'GlobalNodeIDs', gOffset, myMesh % hopr_globalNodeIDs  )

    ! Read local subarray of SideInfo
    firstSide = myMesh % hopr_elemInfo % hostData(3,1)+1

    nLocalSides = myMesh % hopr_elemInfo % hostData(4,nLocalElems)-myMesh % hopr_elemInfo % hostData(3,1)

    offset = (/0, firstSide-1/)
    CALL WriteArray_HDF5(fileId, 'SideInfo', offset, myMesh % hopr_sideInfo  )

    CALL Close_HDF5(fileID)

  END SUBROUTINE Write_HOPr_Mesh3D


END MODULE SELF_Mesh
