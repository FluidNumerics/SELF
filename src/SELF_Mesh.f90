!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_Mesh

  USE SELF_Constants
  !USE hipfort
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
    integer, pointer, dimension(:) :: elemToRank
    integer, pointer, dimension(:) :: offSetElem
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
    integer, pointer, dimension(:,:) :: elemInfo
    real(prec), pointer, dimension(:) :: nodeCoords
    integer, pointer, dimension(:) :: globalNodeIDs
    integer, pointer, dimension(:,:) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS
    PROCEDURE,PUBLIC :: Init => Init_Mesh1D
    PROCEDURE,PUBLIC :: Free => Free_Mesh1D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh1D
    PROCEDURE,PUBLIC :: UniformBlockMesh => UniformBlockMesh_Mesh1D

    PROCEDURE,PUBLIC  :: Write_Mesh => Write_Mesh1D

  END TYPE Mesh1D

  ! Mesh format is set up similar to the HOPr format
  ! See https://hopr-project.org/externals/MeshFormat.pdf

  TYPE,EXTENDS(SEMMesh) :: Mesh2D
    integer, pointer, dimension(:,:,:) :: sideInfo
    real(prec), pointer, dimension(:,:,:,:) :: nodeCoords
    integer, pointer, dimension(:,:) :: elemInfo
    integer, pointer, dimension(:,:,:) :: globalNodeIDs
    integer, pointer, dimension(:,:) :: CGNSCornerMap
    integer, pointer, dimension(:,:) :: CGNSSideMap
    integer, pointer, dimension(:,:) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS
    PROCEDURE,PUBLIC :: Init => Init_Mesh2D
    PROCEDURE,PUBLIC :: Free => Free_Mesh2D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh2D

    PROCEDURE,PUBLIC :: ResetBoundaryConditionType => ResetBoundaryConditionType_Mesh2D

    PROCEDURE,PUBLIC :: Read_HOPr => Read_HOPr_Mesh2D

    PROCEDURE,PUBLIC :: Write_Mesh => Write_Mesh2D

    PROCEDURE,PRIVATE :: RecalculateFlip => RecalculateFlip_Mesh2D

  END TYPE Mesh2D

  TYPE,EXTENDS(SEMMesh) :: Mesh3D
    integer, pointer, dimension(:,:,:) :: sideInfo
    real(prec), pointer, dimension(:,:,:,:,:) :: nodeCoords
    integer, pointer, dimension(:,:) :: elemInfo
    integer, pointer, dimension(:,:,:,:) :: globalNodeIDs
    integer, pointer, dimension(:,:) :: CGNSCornerMap
    integer, pointer, dimension(:,:) :: sideMap
    integer, pointer, dimension(:,:) :: CGNSSideMap
    integer, pointer, dimension(:,:) :: BCType
    CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Mesh3D
    PROCEDURE,PUBLIC :: Free => Free_Mesh3D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Mesh3D

    PROCEDURE,PUBLIC :: Read_HOPr => Read_HOPr_Mesh3D

    PROCEDURE,PUBLIC :: ResetBoundaryConditionType => ResetBoundaryConditionType_Mesh3D

    PROCEDURE,PUBLIC :: Write_Mesh => Write_Mesh3D

    PROCEDURE,PRIVATE :: RecalculateFlip => RecalculateFlip_Mesh3D

  END TYPE Mesh3D

CONTAINS

  SUBROUTINE Init_Mesh1D(this,nGeo,nElem,nNodes,nBCs)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(out) :: this
    INTEGER,INTENT(in) :: nGeo
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nNodes
    INTEGER,INTENT(in) :: nBCs

    this % nGeo = nGeo
    this % nElem = nElem
    this % nGlobalElem = nElem
    this % nNodes = nNodes
    this % nCornerNodes = nElem*2
    this % nUniqueNodes = 0
    this % nBCs = nBCs

    call hipcheck(hipMallocManaged(this % elemInfo,4,nElem,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % nodeCoords,nNodes,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % globalNodeIDs,nNodes,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % BCType,4,nBCs,hipMemAttachGlobal))

    ALLOCATE (this % BCNames(1:nBCs))

  END SUBROUTINE Init_Mesh1D

  SUBROUTINE Free_Mesh1D(this)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(inout) :: this

    this % nElem = 0
    this % nNodes = 0
    this % nCornerNodes = 0
    this % nUniqueNodes = 0
    this % nBCs = 0
    call hipcheck(hipFree(this % elemInfo))
    call hipcheck(hipFree(this % nodeCoords))
    call hipcheck(hipFree(this % globalNodeIDs))
    call hipcheck(hipFree(this % BCType))
    DEALLOCATE (this % BCNames)

  END SUBROUTINE Free_Mesh1D

  SUBROUTINE UpdateDevice_Mesh1D(this)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(inout) :: this

    call hipcheck(hipMemPrefetchAsync(c_loc(this % elemInfo),sizeof(this % elemInfo),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % nodeCoords),sizeof(this % nodeCoords),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % globalNodeIDs),sizeof(this % globalNodeIDs),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % BCType),sizeof(this % BCType),0,c_null_ptr))

  END SUBROUTINE UpdateDevice_Mesh1D

  SUBROUTINE UniformBlockMesh_Mesh1D(this,nGeo,nElem,x)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(out) :: this
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
    CALL this % Init(nGeo,nElem,nNodes,2)
    this % quadrature = GAUSS_LOBATTO

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
      xLinear % interior(1:2,1,iel) = xU(iel:iel + 1)
    END DO

    CALL xLinear % GridInterp(xGeo,.FALSE.)

    ! Set the element information
    nid = 1
    DO iel = 1,nElem
      this % elemInfo(1,iel) = selfLineLinear ! Element Type
      this % elemInfo(2,iel) = 1 ! Element Zone
      this % elemInfo(3,iel) = nid ! Node Index Start
      DO i = 1,nGeo+1
        this % nodeCoords(nid) = xGeo % interior(i,1,iel)
        nid = nid + 1
      END DO
      this % elemInfo(4,iel) = nid - 1 ! Node Index End
    END DO


    CALL this % UpdateDevice()

    CALL xLinear % Free()
    CALL xGeo % Free()
    CALL linearInterp % Free()
    CALL nGeoInterp % Free()

  END SUBROUTINE UniformBlockMesh_Mesh1D

  SUBROUTINE Write_Mesh1D(this,meshFile)
    ! Writes mesh output in HOPR format (serial IO only)
    IMPLICIT NONE
    CLASS(Mesh1D),INTENT(inout) :: this
    CHARACTER(*),INTENT(in) :: meshFile
    ! Local
    INTEGER(HID_T) :: fileId

    CALL Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId)

    CALL WriteAttribute_HDF5(fileId,'nElems',this % nElem)
    CALL WriteAttribute_HDF5(fileId,'Ngeo',this % nGeo)
    CALL WriteAttribute_HDF5(fileId,'nBCs',this % nBCs)

    CALL WriteArray_HDF5(fileId,'BCType',this % bcType)

    ! Read local subarray of ElemInfo
    CALL WriteArray_HDF5(fileId,'ElemInfo',this % elemInfo)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    CALL WriteArray_HDF5(fileId,'NodeCoords',this % nodeCoords)
    CALL WriteArray_HDF5(fileId,'GlobalNodeIDs',this % globalNodeIDs)

    CALL Close_HDF5(fileID)

  END SUBROUTINE Write_Mesh1D

  SUBROUTINE Init_Mesh2D(this,nGeo,nElem,nSides,nNodes,nBCs)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(out) :: this
    INTEGER,INTENT(in) :: nGeo
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nSides
    INTEGER,INTENT(in) :: nNodes
    INTEGER,INTENT(in) :: nBCs
    ! Local
    INTEGER :: i,j,l

    this % nGeo = nGeo
    this % nElem = nElem
    this % nGlobalElem = nElem
    this % nNodes = nNodes
    this % nSides = nSides
    this % nCornerNodes = 0
    this % nUniqueNodes = 0
    this % nUniqueSides = 0
    this % nBCs = nBCs

    call hipcheck(hipMallocManaged(this % elemInfo,6,nElem,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % sideInfo,5,4,nElem,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % nodeCoords,2,nGeo+1,nGeo+1,nElem,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % globalNodeIDs,nGeo+1,nGeo+1,nElem,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % CGNSCornerMap,2,4,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % CGNSSideMap,2,4,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % BCType,4,nBCs,hipMemAttachGlobal))

    ALLOCATE (this % BCNames(1:nBCs))

    ! Create lookup tables to assist with connectivity generation
    this % CGNSCornerMap(1:2,1) = (/1, 1/)
    this % CGNSCornerMap(1:2,2) = (/nGeo+1, 1/)
    this % CGNSCornerMap(1:2,3) = (/nGeo+1, nGeo+1/)
    this % CGNSCornerMap(1:2,4) = (/1, nGeo+1/)

    ! Maps from local corner node id to CGNS side
    this % CGNSSideMap(1:2,1) = (/1,2/)
    this % CGNSSideMap(1:2,2) = (/2,3/)
    this % CGNSSideMap(1:2,3) = (/4,3/)
    this % CGNSSideMap(1:2,4) = (/1,4/)

  END SUBROUTINE Init_Mesh2D

  SUBROUTINE Free_Mesh2D(this)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: this

    this % nElem = 0
    this % nNodes = 0
    this % nSides = 0
    this % nCornerNodes = 0
    this % nUniqueSides = 0
    this % nUniqueNodes = 0
    this % nBCs = 0

    call hipcheck(hipFree(this % elemInfo))
    call hipcheck(hipFree(this % sideInfo))
    call hipcheck(hipFree(this % nodeCoords))
    call hipcheck(hipFree(this % globalNodeIDs))
    call hipcheck(hipFree(this % CGNSCornerMap))
    call hipcheck(hipFree(this % CGNSSideMap))
    call hipcheck(hipFree(this % BCType))
    DEALLOCATE (this % BCNames)

  END SUBROUTINE Free_Mesh2D

  SUBROUTINE UpdateDevice_Mesh2D(this)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: this

    call hipcheck(hipMemPrefetchAsync(c_loc(this % elemInfo),sizeof(this % elemInfo),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % sideInfo),sizeof(this % sideInfo),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % nodeCoords),sizeof(this % nodeCoords),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % globalNodeIDs),sizeof(this % globalNodeIDs),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % BCType),sizeof(this % BCType),0,c_null_ptr))

  END SUBROUTINE UpdateDevice_Mesh2D

  SUBROUTINE ResetBoundaryConditionType_Mesh2D(this,bcid)
    !! This method can be used to reset all of the boundary elements
    !! boundary condition type to the desired value.
    !!
    !! Note that ALL physical boundaries will be set to have this boundary 
    !! condition
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: this
    INTEGER, INTENT(in) :: bcid
    ! Local
    INTEGER :: iSide,iEl,e2      

    DO iEl = 1, this % nElem
      DO iSide = 1, 4

        e2 = this % sideInfo(3,iSide,iEl)

        IF( e2 == 0 )THEN
          this % sideInfo(5,iSide,iEl) = bcid
        ENDIF

      ENDDO
    ENDDO

  END SUBROUTINE ResetBoundaryConditionType_Mesh2D 

  SUBROUTINE Read_HOPr_Mesh2D(this,meshFile,decomp)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    ! Adapted for 2D Mesh : Note that HOPR does not have 2D mesh output.
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(out) :: this
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
    integer, dimension(:,:), allocatable :: hopr_elemInfo
    integer, dimension(:,:), allocatable :: hopr_sideInfo
    real(prec), dimension(:,:), allocatable :: hopr_nodeCoords
    integer, dimension(:), allocatable :: hopr_globalNodeIDs
    integer, dimension(:,:), allocatable :: bcType


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
    allocate(bcTypes(1:4,1:nBCS))

    IF ( decomp % mpiEnabled )THEN
      offset(:) = 0
      CALL ReadArray_HDF5(fileId,'BCType',bcType,offset)
    ELSE
      CALL ReadArray_HDF5(fileId,'BCType',bcType)
    ENDIF

    ! Read local subarray of ElemInfo
    CALL decomp % GenerateDecomposition(nGlobalElem,nUniqueSides3D)

    firstElem = decomp % offsetElem(decomp % rankId) + 1
    nLocalElems = decomp % offsetElem(decomp % rankId + 1) - &
                  decomp % offsetElem(decomp % rankId)

    ! Allocate Space for hopr_elemInfo!
    allocate(hopr_elemInfo(1:6,1:nLocalElems))

    IF ( decomp % mpiEnabled )THEN
      offset = (/0,firstElem - 1/)
      CALL ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo,offset)
    ELSE
      CALL ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo)
    ENDIF

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode = hopr_elemInfo(5,1) + 1
    nLocalNodes3D = hopr_elemInfo(6,nLocalElems) - hopr_elemInfo(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !
    allocate( hopr_nodeCoords(1:3,nLocalNodes3D), hopr_globalNodeIDs(1:nLocalNodes3D))

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
    firstSide = hopr_elemInfo(3,1) + 1
    nLocalSides3D = hopr_elemInfo(4,nLocalElems) - hopr_elemInfo(3,1)

    ! Allocate space for hopr_sideInfo
    allocate( hopr_sideInfo(1:5,1:nLocalSides3D) )
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

    CALL this % Init(nGeo,nLocalElems,nLocalSides2D,nLocalNodes2D,nBCs)

    ! Copy data from local arrays into this
    !  elemInfo(1:6,iEl)
    !    1 - Element Type
    !    2 - Zone
    !    3 - offset index for side array (not needed when all quads are assumed)
    !    4 - last index for side array (not needed when all quads are assumed)
    !    5 - offset index for node array (not needed when all quads are assumed)
    !    6 - last index for node array (not needed when all quads are assumed)
    this % elemInfo = hopr_elemInfo
    this % quadrature = UNIFORM  ! HOPr uses uniformly spaced points

    ! Grab the node coordinates (x and y only) from the "bottom" layer of the extruded mesh
    DO eid = 1, this % nElem
      DO j = 1,nGeo+1
        DO i = 1,nGeo+1
          nid = i+1 + (nGeo+1)*(j + (nGeo+1)*((nGeo+1)*(eid-1)))
          this % nodeCoords(1:2,i,j,eid) = hopr_nodeCoords(1:2,nid)
          this % globalNodeIDs(i,j,eid) = hopr_globalNodeIDs(nid)
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
    DO eid = 1,this % nElem
      DO lsid = 1,4
        ! Calculate the 3-D side ID from the 2-D local side id and element ID
        iSide = lsid + 1 + 6*(eid-1)
        this % sideInfo(1:5,lsid,eid) = hopr_sideInfo(1:5,iSide)
        ! Adjust the secondary side index for 2-D
        this % sideInfo(4,lsid,eid) = this % sideInfo(4,lsid,eid)-10
      ENDDO
    ENDDO

    CALL this % RecalculateFlip()

    CALL this % UpdateDevice()

    deallocate( hopr_elemInfo, hopr_nodeCoords, hopr_globalNodeIDs, hopr_sideInfo)

  END SUBROUTINE Read_HOPr_Mesh2D

  SUBROUTINE RecalculateFlip_Mesh2D(this,decomp)  
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: this
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
    INTEGER :: nid1(1:2,1:4,1:this % nElem)
    INTEGER :: nid2(1:2,1:4,1:this % nElem)
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

    ALLOCATE(requests(1:this % nSides*2))
    ALLOCATE(stats(MPI_STATUS_SIZE,1:this % nSides*2))

    IF (PRESENT(decomp)) THEN
      rankId = decomp % rankId
      offset = decomp % offsetElem(rankId)
    ELSE
      rankId = 0
      offset = 0
    ENDIF

    msgCount = 0
    DO e1 = 1,this % nElem
      DO s1 = 1,4

        e2Global = this % sideInfo(3,s1,e1)
        e2 = e2Global - offset
        s2 = this % sideInfo(4,s1,e1)/10
        flip = this % sideInfo(4,s1,e1) - s2*10
        bcid = this % sideInfo(5,s1,e1)

        IF (bcid == 0) THEN

          IF (PRESENT(decomp)) THEN
            neighborRank = decomp % elemToRank(e2Global)
          ELSE
            neighborRank = 0
          ENDIF

          IF (neighborRank == rankId) THEN

            lnid1 = this % CGNSSideMap(1:2,s1) ! local CGNS corner node ids for element 1 side
            lnid2 = this % CGNSSideMap(1:2,s2) ! local CGNS corner node ids for element 2 side

            DO l = 1, 2

              i = this % CGNSCornerMap(1,lnid1(l))
              j = this % CGNSCornerMap(2,lnid1(l))
              nid1(l,s1,e1) = this % globalNodeIDs(i,j,e1)

              i = this % CGNSCornerMap(1,lnid2(l))
              j = this % CGNSCornerMap(2,lnid2(l))
              nid2(l,s1,e1) = this % globalNodeIDs(i,j,e2)

            ENDDO

          ELSE ! In this case, we need to exchange

            globalSideId = ABS(this % sideInfo(2,s1,e1))

            lnid1 = this % CGNSSideMap(1:2,s1) ! local CGNS corner node ids for element 1 side

            DO l = 1, 2

              i = this % CGNSCornerMap(1,lnid1(l))
              j = this % CGNSCornerMap(2,lnid1(l))
              nid1(l,s1,e1) = this % globalNodeIDs(i,j,e1)

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

    DO e1 = 1,this % nElem
      DO s1 = 1,4

        s2 = this % sideInfo(4,s1,e1)/10
        bcid = this % sideInfo(5,s1,e1)
        nloc1(1:2) = nid1(1:2,s1,e1)
        nloc2(1:2) = nid2(1:2,s1,e1)

        IF (bcid == 0) THEN
          theyMatch = CompareArray( nloc1, nloc2, 2 )

          IF( theyMatch )THEN
            this % sideInfo(4,s1,e1) = 10*s2
          ELSE
            this % sideInfo(4,s1,e1) = 10*s2+1
          ENDIF


        ENDIF

      ENDDO
    ENDDO

    DEALLOCATE(requests)
    DEALLOCATE(stats)

  END SUBROUTINE RecalculateFlip_Mesh2D

  SUBROUTINE Write_Mesh2D(this,meshFile)
    ! Writes mesh output in HOPR format (serial only)
    IMPLICIT NONE
    CLASS(Mesh2D),INTENT(inout) :: this
    CHARACTER(*),INTENT(in) :: meshFile
    ! Local
    INTEGER(HID_T) :: fileId

    CALL Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId)
    CALL WriteAttribute_HDF5(fileId,'nElems',this % nElem)
    CALL WriteAttribute_HDF5(fileId,'Ngeo',this % nGeo)
    CALL WriteAttribute_HDF5(fileId,'nBCs',this % nBCs)

    CALL WriteArray_HDF5(fileId,'BCType',this % bcType)

    ! Write local subarray of ElemInfo
    CALL WriteArray_HDF5(fileId,'ElemInfo',this % elemInfo)

    ! Write local subarray of NodeCoords and GlobalNodeIDs
    CALL WriteArray_HDF5(fileId,'NodeCoords',this % nodeCoords)
    CALL WriteArray_HDF5(fileId,'GlobalNodeIDs',this % globalNodeIDs)

    ! Write local subarray of SideInfo
    CALL WriteArray_HDF5(fileId,'SideInfo',this % sideInfo)

    CALL Close_HDF5(fileID)

  END SUBROUTINE Write_Mesh2D

  SUBROUTINE Init_Mesh3D(this,nGeo,nElem,nSides,nNodes,nBCs)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(out) :: this
    INTEGER,INTENT(in) :: nGeo
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nSides
    INTEGER,INTENT(in) :: nNodes
    INTEGER,INTENT(in) :: nBCs
    ! Local
    INTEGER :: i,j,k,l

    this % nElem = nElem
    this % nGlobalElem = nElem
    this % nGeo = nGeo
    this % nSides = nSides
    this % nNodes = nNodes
    this % nCornerNodes = 0
    this % nUniqueSides = 0
    this % nUniqueNodes = 0
    this % nBCs = nBCs

    call hipcheck(hipMallocManaged(this % elemInfo,6,nElem,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % sideInfo,5,6,nElem,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % nodeCoords,3,nGeo+1,nGeo+1,nGeo+1,nElem,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % globalNodeIDs,nGeo+1,nGeo+1,nGeo+1,nElem,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % CGNSCornerMap,3,8,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % CGNSSideMap,4,6,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % sideMap,4,6,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % BCType,4,nBCs,hipMemAttachGlobal))

    ALLOCATE (this % BCNames(1:nBCs))

    ! Create lookup tables to assist with connectivity generation
    this % CGNSCornerMap(1:3,1) = (/1,1,1/) ! Bottom-South-West
    this % CGNSCornerMap(1:3,2) = (/nGeo+1,1,1/) ! Bottom-South-East
    this % CGNSCornerMap(1:3,3) = (/nGeo+1,nGeo+1,1/)! Bottom-North-East
    this % CGNSCornerMap(1:3,4) = (/1,nGeo+1,1/)! Bottom-North-West
    this % CGNSCornerMap(1:3,5) = (/1,1,nGeo+1/) ! Top-South-West
    this % CGNSCornerMap(1:3,6) = (/nGeo+1,1,nGeo+1/) ! Top-South-East
    this % CGNSCornerMap(1:3,7) = (/nGeo+1,nGeo+1,nGeo+1/)! Top-North-East
    this % CGNSCornerMap(1:3,8) = (/1,nGeo+1,nGeo+1/)! Top-North-West

    ! Maps from local corner node id to CGNS side
    this % CGNSSideMap(1:4,1) = (/1,4,3,2/)
    this % CGNSSideMap(1:4,2) = (/1,2,6,5/)
    this % CGNSSideMap(1:4,3) = (/2,3,7,6/)
    this % CGNSSideMap(1:4,4) = (/3,4,8,7/)
    this % CGNSSideMap(1:4,5) = (/1,5,8,4/)
    this % CGNSSideMap(1:4,6) = (/5,6,7,8/)

    this % sideMap(1:4,1) = (/1,2,3,4/) ! Bottom
    this % sideMap(1:4,2) = (/1,2,6,5/) ! South
    this % sideMap(1:4,3) = (/2,3,7,6/) ! East
    this % sideMap(1:4,4) = (/4,3,7,8/) ! North
    this % sideMap(1:4,5) = (/1,4,8,5/) ! West
    this % sideMap(1:4,6) = (/5,6,7,8/) ! Top

  END SUBROUTINE Init_Mesh3D

  SUBROUTINE Free_Mesh3D(this)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: this

    this % nElem = 0
    this % nSides = 0
    this % nNodes = 0
    this % nCornerNodes = 0
    this % nUniqueSides = 0
    this % nUniqueNodes = 0
    this % nBCs = 0

    call hipcheck(hipFree(this % elemInfo))
    call hipcheck(hipFree(this % sideInfo))
    call hipcheck(hipFree(this % nodeCoords))
    call hipcheck(hipFree(this % globalNodeIDs))
    call hipcheck(hipFree(this % CGNSCornerMap))
    call hipcheck(hipFree(this % sideMap))
    call hipcheck(hipFree(this % CGNSSideMap))
    call hipcheck(hipFree(this % BCType))

    DEALLOCATE (this % BCNames)

  END SUBROUTINE Free_Mesh3D

  SUBROUTINE UpdateDevice_Mesh3D(this)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: this

    call hipcheck(hipMemPrefetchAsync(c_loc(this % elemInfo),sizeof(this % elemInfo),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % sideInfo),sizeof(this % sideInfo),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % nodeCoords),sizeof(this % nodeCoords),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % globalNodeIDs),sizeof(this % globalNodeIDs),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % BCType),sizeof(this % BCType),0,c_null_ptr))

  END SUBROUTINE UpdateDevice_Mesh3D

  SUBROUTINE ResetBoundaryConditionType_Mesh3D(this,bcid)
    !! This method can be used to reset all of the boundary elements
    !! boundary condition type to the desired value.
    !!
    !! Note that ALL physical boundaries will be set to have this boundary 
    !! condition
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: this
    INTEGER, INTENT(in) :: bcid
    ! Local
    INTEGER :: iSide,iEl,e2      

    DO iEl = 1, this % nElem
      DO iSide = 1, 6

        e2 = this % sideInfo(3,iSide,iEl)

        IF( e2 == 0 )THEN
          this % sideInfo(5,iSide,iEl) = bcid
        ENDIF

      ENDDO
    ENDDO

  END SUBROUTINE ResetBoundaryConditionType_Mesh3D 

  SUBROUTINE RecalculateFlip_Mesh3D(this,decomp)  
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: this
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
    INTEGER :: nid1(1:4,1:6,1:this % nElem)
    INTEGER :: nid2(1:4,1:6,1:this % nElem)
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


    ALLOCATE(requests(1:this % nSides*2))
    ALLOCATE(stats(MPI_STATUS_SIZE,1:this % nSides*2))

    IF (PRESENT(decomp)) THEN
      rankId = decomp % rankId
      offset = decomp % offsetElem(rankId)
    ELSE
      rankId = 0
      offset = 0
    ENDIF

    msgCount = 0
    DO e1 = 1,this % nElem
      DO s1 = 1,6

        e2Global = this % sideInfo(3,s1,e1)
        e2 = e2Global - offset
        s2 = this % sideInfo(4,s1,e1)/10
        flip = this % sideInfo(4,s1,e1) - s2*10
        bcid = this % sideInfo(5,s1,e1)

        IF (bcid == 0) THEN

          IF (PRESENT(decomp)) THEN
            neighborRank = decomp % elemToRank(e2Global)
          ELSE
            neighborRank = 0
          ENDIF

          IF (neighborRank == rankId) THEN

            lnid1 = this % sideMap(1:4,s1) ! local CGNS corner node ids for element 1 side
            lnid2 = this % sideMap(1:4,s2) ! local CGNS corner node ids for element 2 side

            DO l = 1, 4
              
              i = this % CGNSCornerMap(1,lnid1(l))
              j = this % CGNSCornerMap(2,lnid1(l))
              k = this % CGNSCornerMap(3,lnid1(l))
              nid1(l,s1,e1) = this % globalNodeIDs(i,j,k,e1)

              i = this % CGNSCornerMap(1,lnid2(l))
              j = this % CGNSCornerMap(2,lnid2(l))
              k = this % CGNSCornerMap(3,lnid2(l))
              nid2(l,s1,e1) = this % globalNodeIDs(i,j,k,e1)

            ENDDO

          ELSE ! In this case, we need to exchange

            globalSideId = ABS(this % sideInfo(2,s1,e1))

            lnid1 = this % sideMap(1:4,s1) ! local CGNS corner node ids for element 1 side

            DO l = 1, 4
              
              i = this % CGNSCornerMap(1,lnid1(l))
              j = this % CGNSCornerMap(2,lnid1(l))
              k = this % CGNSCornerMap(3,lnid1(l))
              nid1(l,s1,e1) = this % globalNodeIDs(i,j,k,e1)

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

    DO e1 = 1,this % nElem
      DO s1 = 1,6

        s2 = this % sideInfo(4,s1,e1)/10
        bcid = this % sideInfo(5,s1,e1)
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

          this % sideInfo(4,s1,e1) = 10*s2+nShifts

        ENDIF

      ENDDO
    ENDDO

    DEALLOCATE(requests)
    DEALLOCATE(stats)

  END SUBROUTINE RecalculateFlip_Mesh3D

  SUBROUTINE Read_HOPr_Mesh3D(this,meshFile,decomp)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(out) :: this
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
    integer, dimension(:,:), allocatable :: hopr_elemInfo
    integer, dimension(:,:), allocatable :: hopr_sideInfo
    real(prec), dimension(:,:), allocatable :: hopr_nodeCoords
    integer, dimension(:), allocatable :: hopr_globalNodeIDs
    integer, dimension(:,:), allocatable :: bcType


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

    firstElem = decomp % offsetElem(decomp % rankId) + 1
    nLocalElems = decomp % offsetElem(decomp % rankId + 1) - &
                  decomp % offsetElem(decomp % rankId)

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
    firstNode = hopr_elemInfo(5,1) + 1
    nLocalNodes = hopr_elemInfo(6,nLocalElems) - hopr_elemInfo(5,1)

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
    firstSide = hopr_elemInfo(3,1) + 1
    nLocalSides = hopr_elemInfo(4,nLocalElems) - hopr_elemInfo(3,1)

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

    CALL this % Init(nGeo,nLocalElems,nLocalSides,nLocalNodes,nBCs)

    ! Copy data from local arrays into this
    this % elemInfo = hopr_elemInfo
    this % nUniqueSides = nUniqueSides
    this % quadrature = UNIFORM

    ! Grab the node coordinates
    DO eid = 1, this % nElem
      DO k = 1,nGeo+1
        DO j = 1,nGeo+1
          DO i = 1,nGeo+1
            nid = i+1 + (nGeo+1)*(j + (nGeo+1)*(k + (nGeo+1)*(eid-1)))
            this % nodeCoords(1:3,i,j,k,eid) = hopr_nodeCoords(1:3,nid)
            this % globalNodeIDs(i,j,k,eid) = hopr_globalNodeIDs(nid)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    iSide = 0 
    DO eid = 1,this % nElem
      DO lsid = 1,6
        iSide = iSide + 1
        this % sideInfo(1:5,lsid,eid) = hopr_sideInfo(1:5,iSide)
      ENDDO
    ENDDO

    CALL this % RecalculateFlip()

    CALL this % UpdateDevice()

    CALL hopr_elemInfo % Free()
    CALL hopr_nodeCoords % Free()
    CALL hopr_globalNodeIDs % Free()
    CALL hopr_sideInfo % Free()

  END SUBROUTINE Read_HOPr_Mesh3D

  SUBROUTINE Write_Mesh3D(this,meshFile)
    ! Writes mesh output in HOPR format (serial only)
    IMPLICIT NONE
    CLASS(Mesh3D),INTENT(inout) :: this
    CHARACTER(*),INTENT(in) :: meshFile
    ! Local
    INTEGER(HID_T) :: fileId

    CALL Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId)

    CALL WriteAttribute_HDF5(fileId,'nElems',this % nElem)
    CALL WriteAttribute_HDF5(fileId,'Ngeo',this % nGeo)
    CALL WriteAttribute_HDF5(fileId,'nBCs',this % nBCs)

    CALL WriteArray_HDF5(fileId,'BCType',this % bcType)
    CALL WriteArray_HDF5(fileId,'ElemInfo',this % elemInfo)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    CALL WriteArray_HDF5(fileId,'NodeCoords',this % nodeCoords)
    CALL WriteArray_HDF5(fileId,'GlobalNodeIDs',this % globalNodeIDs)

    ! Read local subarray of SideInfo
    CALL WriteArray_HDF5(fileId,'SideInfo',this % sideInfo)

    CALL Close_HDF5(fileID)

  END SUBROUTINE Write_Mesh3D

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

    IF (ASSOCIATED(this % offSetElem)) THEN
      CALL this % offSetElem % Free()
    ENDIF
    IF (ASSOCIATED(this % elemToRank)) THEN
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
    WRITE (msg2,'(I5)') this % offSetElem(this % rankId + 1)-&
                        this % offSetElem(this % rankId)
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
                      this % offSetElem)

    DO iel = 1,nElem
      CALL ElemToRank(this % nRanks, &
                      this % offSetElem, &
                      iel, &
                      this % elemToRank(iel))
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
