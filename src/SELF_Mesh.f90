!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
module SELF_Mesh

  use SELF_Constants
  !USE hipfort
  use SELF_Lagrange
  use SELF_Data
  use SELF_SupportRoutines
  use SELF_HDF5

  ! External Libs !
  use HDF5

  use iso_c_binding

  implicit none

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
  integer,parameter :: selfLineLinear = 1
  integer,parameter :: selfLineNonlinear = 2
  integer,parameter :: selfTriangleLinear = 3
  integer,parameter :: selfQuadLinear = 4
  integer,parameter :: selfQuadBilinear = 14
  integer,parameter :: selfTriangleNonlinear = 23
  integer,parameter :: selfQuadNonlinear = 24
  integer,parameter :: selfTetrahedronLinear = 104
  integer,parameter :: selfPyramidLinear = 105
  integer,parameter :: selfPrismLinear = 106
  integer,parameter :: selfHexahedronLinear = 108
  integer,parameter :: selfPyramidBilinear = 115
  integer,parameter :: selfPrismBilinear = 116
  integer,parameter :: selfHexahedronBilinear = 118
  integer,parameter :: selfTetrahedronNonlinear = 204
  integer,parameter :: selfPyramidNonlinear = 205
  integer,parameter :: selfPrismNonlinear = 206
  integer,parameter :: selfHexahedronNonlinear = 208

  !
  integer,parameter :: selfMinNodalValence2D = 4
  integer,parameter :: selfMinNodalValence3D = 8
  integer,parameter :: selfMaxNodalValence2D = 6
  integer,parameter :: selfMaxNodalValence3D = 10

  ! Side Ordering
  integer,parameter :: selfSide2D_South = 1
  integer,parameter :: selfSide2D_East = 2
  integer,parameter :: selfSide2D_North = 3
  integer,parameter :: selfSide2D_West = 4
  integer,parameter :: selfSide3D_Bottom = 1
  integer,parameter :: selfSide3D_South = 2
  integer,parameter :: selfSide3D_East = 3
  integer,parameter :: selfSide3D_North = 4
  integer,parameter :: selfSide3D_West = 5
  integer,parameter :: selfSide3D_Top = 6
  !
  integer,parameter :: self_BCDefault = 1
  integer,parameter :: self_nBCsDefault = 5

  !==============================================!
  ! --------------- File Types------------------ !
  !==============================================!
  integer,parameter :: SELF_MESH_ISM_V2_2D = 1
  integer,parameter :: SELF_MESH_ISM_V2_3D = 2
  integer,parameter :: SELF_MESH_HOPR_2D = 3
  integer,parameter :: SELF_MESH_HOPR_3D = 4

  type MeshSpec
    character(self_FileNameLength) :: filename
    integer :: fileType

    logical :: blockMesh
    integer :: blockMesh_nGeo
    integer :: blockMesh_nElemX
    integer :: blockMesh_nElemY
    integer :: blockMesh_nElemZ
    real(prec) :: blockMesh_x0,blockMesh_x1
    real(prec) :: blockMesh_y0,blockMesh_y1
    real(prec) :: blockMesh_z0,blockMesh_z1

  endtype MeshSpec

  type MPILayer
    logical :: mpiEnabled
    integer :: mpiComm
    integer :: mpiPrec
    integer :: rankId
    integer :: nRanks
    integer :: nElem
    integer :: maxMsg
    integer :: msgCount
    integer,pointer,dimension(:) :: elemToRank
    integer,pointer,dimension(:) :: offSetElem
    integer,allocatable :: requests(:)
    integer,allocatable :: stats(:,:)

  contains

    procedure :: Init => Init_MPILayer
    procedure :: Free => Free_MPILayer
    procedure :: Finalize => Finalize_MPILayer

    procedure :: GenerateDecomposition => GenerateDecomposition_MPILayer
    procedure :: SetElemToRank
    procedure :: SetMaxMsg

    procedure,public :: FinalizeMPIExchangeAsync

    generic,public :: GlobalReduce => GlobalReduce_RealScalar
    procedure,private :: GlobalReduce_RealScalar

  endtype MPILayer

  type :: SEMMesh
    integer :: nGeo
    integer :: nElem
    integer :: nGlobalElem
    integer :: nNodes
    integer :: nSides
    integer :: nCornerNodes
    integer :: nUniqueNodes
    integer :: nUniqueSides
    integer :: nBCs
    integer :: quadrature

  endtype SEMMesh

  type,extends(SEMMesh) :: Mesh1D
    integer,pointer,dimension(:,:) :: elemInfo
    real(prec),pointer,dimension(:) :: nodeCoords
    integer,pointer,dimension(:) :: globalNodeIDs
    integer,pointer,dimension(:,:) :: BCType
    character(LEN=255),allocatable :: BCNames(:)

  contains
    procedure,public :: Init => Init_Mesh1D
    procedure,public :: Free => Free_Mesh1D
    procedure,public :: UniformBlockMesh => UniformBlockMesh_Mesh1D

    procedure,public  :: Write_Mesh => Write_Mesh1D

  endtype Mesh1D

  ! Mesh format is set up similar to the HOPr format
  ! See https://hopr-project.org/externals/MeshFormat.pdf

  type,extends(SEMMesh) :: Mesh2D
    integer,pointer,dimension(:,:,:) :: sideInfo
    real(prec),pointer,dimension(:,:,:,:) :: nodeCoords
    integer,pointer,dimension(:,:) :: elemInfo
    integer,pointer,dimension(:,:,:) :: globalNodeIDs
    integer,pointer,dimension(:,:) :: CGNSCornerMap
    integer,pointer,dimension(:,:) :: CGNSSideMap
    integer,pointer,dimension(:,:) :: BCType
    character(LEN=255),allocatable :: BCNames(:)

  contains
    procedure,public :: Init => Init_Mesh2D
    procedure,public :: Free => Free_Mesh2D

    procedure,public :: ResetBoundaryConditionType => ResetBoundaryConditionType_Mesh2D

    procedure,public :: Read_HOPr => Read_HOPr_Mesh2D

    procedure,public :: Write_Mesh => Write_Mesh2D

    procedure,private :: RecalculateFlip => RecalculateFlip_Mesh2D

  endtype Mesh2D

  type,extends(SEMMesh) :: Mesh3D
    integer,pointer,dimension(:,:,:) :: sideInfo
    real(prec),pointer,dimension(:,:,:,:,:) :: nodeCoords
    integer,pointer,dimension(:,:) :: elemInfo
    integer,pointer,dimension(:,:,:,:) :: globalNodeIDs
    integer,pointer,dimension(:,:) :: CGNSCornerMap
    integer,pointer,dimension(:,:) :: sideMap
    integer,pointer,dimension(:,:) :: CGNSSideMap
    integer,pointer,dimension(:,:) :: BCType
    character(LEN=255),allocatable :: BCNames(:)

  contains

    procedure,public :: Init => Init_Mesh3D
    procedure,public :: Free => Free_Mesh3D

    procedure,public :: Read_HOPr => Read_HOPr_Mesh3D

    procedure,public :: ResetBoundaryConditionType => ResetBoundaryConditionType_Mesh3D

    procedure,public :: Write_Mesh => Write_Mesh3D

    procedure,private :: RecalculateFlip => RecalculateFlip_Mesh3D

  endtype Mesh3D

contains

  subroutine Init_Mesh1D(this,nGeo,nElem,nNodes,nBCs)
    implicit none
    class(Mesh1D),intent(out) :: this
    integer,intent(in) :: nGeo
    integer,intent(in) :: nElem
    integer,intent(in) :: nNodes
    integer,intent(in) :: nBCs

    this%nGeo = nGeo
    this%nElem = nElem
    this%nGlobalElem = nElem
    this%nNodes = nNodes
    this%nCornerNodes = nElem*2
    this%nUniqueNodes = 0
    this%nBCs = nBCs

    allocate(this%elemInfo(1:4,1:nElem))
    allocate(this%nodeCoords(1:nNodes))
    allocate(this%globalNodeIDs(1:nNodes))
    allocate(this%BCType(1:4,1:nBCs))

    !$omp target enter data map(alloc: this % elemInfo)
    !$omp target enter data map(alloc: this % nodeCoords)
    !$omp target enter data map(alloc: this % globalNodeIDs)
    !$omp target enter data map(alloc: this % BCType)

    allocate(this%BCNames(1:nBCs))

  endsubroutine Init_Mesh1D

  subroutine Free_Mesh1D(this)
    implicit none
    class(Mesh1D),intent(inout) :: this

    this%nElem = 0
    this%nNodes = 0
    this%nCornerNodes = 0
    this%nUniqueNodes = 0
    this%nBCs = 0
    deallocate(this%elemInfo)
    deallocate(this%nodeCoords)
    deallocate(this%globalNodeIDs)
    deallocate(this%BCType)
    !$omp target exit data map(delete: this % elemInfo)
    !$omp target exit data map(delete: this % nodeCoords)
    !$omp target exit data map(delete: this % globalNodeIDs)
    !$omp target exit data map(delete: this % BCType)
    deallocate(this%BCNames)

  endsubroutine Free_Mesh1D

  subroutine UniformBlockMesh_Mesh1D(this,nGeo,nElem,x)
    implicit none
    class(Mesh1D),intent(out) :: this
    integer,intent(in) :: nGeo
    integer,intent(in) :: nElem
    real(prec),intent(in) :: x(1:2)
    ! Local
    integer :: iel
    integer :: nid,nNodes
    integer :: i
    real(prec) :: xU(1:nElem+1)
    type(Lagrange),target :: linearInterp
    type(Lagrange),target :: nGeoInterp
    type(Scalar1D) :: xLinear
    type(Scalar1D) :: xGeo

    nNodes = nElem*(nGeo+1)
    call this%Init(nGeo,nElem,nNodes,2)
    this%quadrature = GAUSS_LOBATTO

    ! Set the hopr_nodeCoords
    xU = UniformPoints(x(1),x(2),1,nElem+1)

    call linearInterp%Init(1,GAUSS_LOBATTO, &
                           nGeo,GAUSS_LOBATTO)

    call nGeoInterp%Init(nGeo,GAUSS_LOBATTO, &
                         nGeo,GAUSS_LOBATTO)

    ! Create a linear interpolant to interpolate to nGeo grid
    call xLinear%Init(linearInterp,1,nElem)
    call xGeo%Init(nGeoInterp,1,nElem)

    do iel = 1,nElem
      xLinear%interior(1:2,iel,1) = xU(iel:iel+1)
    enddo

    call xLinear%GridInterp(xGeo)

    ! Set the element information
    nid = 1
    do iel = 1,nElem
      this%elemInfo(1,iel) = selfLineLinear ! Element Type
      this%elemInfo(2,iel) = 1 ! Element Zone
      this%elemInfo(3,iel) = nid ! Node Index Start
      do i = 1,nGeo+1
        this%nodeCoords(nid) = xGeo%interior(i,iel,1)
        nid = nid+1
      enddo
      this%elemInfo(4,iel) = nid-1 ! Node Index End
    enddo

    call xLinear%Free()
    call xGeo%Free()
    call linearInterp%Free()
    call nGeoInterp%Free()

  endsubroutine UniformBlockMesh_Mesh1D

  subroutine Write_Mesh1D(this,meshFile)
    ! Writes mesh output in HOPR format (serial IO only)
    implicit none
    class(Mesh1D),intent(inout) :: this
    character(*),intent(in) :: meshFile
    ! Local
    integer(HID_T) :: fileId

    call Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId)

    call WriteAttribute_HDF5(fileId,'nElems',this%nElem)
    call WriteAttribute_HDF5(fileId,'Ngeo',this%nGeo)
    call WriteAttribute_HDF5(fileId,'nBCs',this%nBCs)

    call WriteArray_HDF5(fileId,'BCType',this%bcType)

    ! Read local subarray of ElemInfo
    call WriteArray_HDF5(fileId,'ElemInfo',this%elemInfo)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    call WriteArray_HDF5(fileId,'NodeCoords',this%nodeCoords)
    call WriteArray_HDF5(fileId,'GlobalNodeIDs',this%globalNodeIDs)

    call Close_HDF5(fileID)

  endsubroutine Write_Mesh1D

  subroutine Init_Mesh2D(this,nGeo,nElem,nSides,nNodes,nBCs)
    implicit none
    class(Mesh2D),intent(out) :: this
    integer,intent(in) :: nGeo
    integer,intent(in) :: nElem
    integer,intent(in) :: nSides
    integer,intent(in) :: nNodes
    integer,intent(in) :: nBCs
    ! Local
    integer :: i,j,l

    this%nGeo = nGeo
    this%nElem = nElem
    this%nGlobalElem = nElem
    this%nNodes = nNodes
    this%nSides = nSides
    this%nCornerNodes = 0
    this%nUniqueNodes = 0
    this%nUniqueSides = 0
    this%nBCs = nBCs

    allocate(this%elemInfo(1:6,1:nElem))
    allocate(this%sideInfo(1:5,1:4,1:nElem))
    allocate(this%nodeCoords(1:2,1:nGeo+1,1:nGeo+1,1:nElem))
    allocate(this%globalNodeIDs(1:nGeo+1,1:nGeo+1,1:nElem))
    allocate(this%CGNSCornerMap(1:2,1:4))
    allocate(this%CGNSSideMap(1:2,1:4))
    allocate(this%BCType(1:4,1:nBCs))

    !$omp target enter data map(alloc: this % elemInfo)
    !$omp target enter data map(alloc: this % sideInfo)
    !$omp target enter data map(alloc: this % nodeCoords)
    !$omp target enter data map(alloc: this % globalNodeIDs)
    !$omp target enter data map(alloc: this % CGNSCornerMap)
    !$omp target enter data map(alloc: this % CGNSSideMap)
    !$omp target enter data map(alloc: this % BCType)

    allocate(this%BCNames(1:nBCs))

    ! Create lookup tables to assist with connectivity generation
    this%CGNSCornerMap(1:2,1) = (/1,1/)
    this%CGNSCornerMap(1:2,2) = (/nGeo+1,1/)
    this%CGNSCornerMap(1:2,3) = (/nGeo+1,nGeo+1/)
    this%CGNSCornerMap(1:2,4) = (/1,nGeo+1/)

    ! Maps from local corner node id to CGNS side
    this%CGNSSideMap(1:2,1) = (/1,2/)
    this%CGNSSideMap(1:2,2) = (/2,3/)
    this%CGNSSideMap(1:2,3) = (/4,3/)
    this%CGNSSideMap(1:2,4) = (/1,4/)

    !$omp target update to(this % CGNSCornerMap)
    !$omp target update to(this % CGNSSideMap)

  endsubroutine Init_Mesh2D

  subroutine Free_Mesh2D(this)
    implicit none
    class(Mesh2D),intent(inout) :: this

    this%nElem = 0
    this%nNodes = 0
    this%nSides = 0
    this%nCornerNodes = 0
    this%nUniqueSides = 0
    this%nUniqueNodes = 0
    this%nBCs = 0

    deallocate(this%elemInfo)
    deallocate(this%sideInfo)
    deallocate(this%nodeCoords)
    deallocate(this%globalNodeIDs)
    deallocate(this%CGNSCornerMap)
    deallocate(this%CGNSSideMap)
    deallocate(this%BCType)
    deallocate(this%BCNames)

    !$omp target exit data map(delete: this % elemInfo)
    !$omp target exit data map(delete: this % sideInfo)
    !$omp target exit data map(delete: this % nodeCoords)
    !$omp target exit data map(delete: this % globalNodeIDs)
    !$omp target exit data map(delete: this % CGNSCornerMap)
    !$omp target exit data map(delete: this % CGNSSideMap)
    !$omp target exit data map(delete: this % BCType)

  endsubroutine Free_Mesh2D

  subroutine ResetBoundaryConditionType_Mesh2D(this,bcid)
    !! This method can be used to reset all of the boundary elements
    !! boundary condition type to the desired value.
    !!
    !! Note that ALL physical boundaries will be set to have this boundary
    !! condition
    implicit none
    class(Mesh2D),intent(inout) :: this
    integer,intent(in) :: bcid
    ! Local
    integer :: iSide,iEl,e2

    do iEl = 1,this%nElem
      do iSide = 1,4

        e2 = this%sideInfo(3,iSide,iEl)

        if(e2 == 0) then
          this%sideInfo(5,iSide,iEl) = bcid
        endif

      enddo
    enddo

  endsubroutine ResetBoundaryConditionType_Mesh2D

  subroutine Read_HOPr_Mesh2D(this,meshFile,decomp)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    ! Adapted for 2D Mesh : Note that HOPR does not have 2D mesh output.
    implicit none
    class(Mesh2D),intent(out) :: this
    character(*),intent(in) :: meshFile
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer(HID_T) :: fileId
    integer(HID_T) :: offset(1:2),gOffset(1)
    integer :: nGlobalElem
    integer :: firstElem
    integer :: firstNode
    integer :: firstSide
    integer :: nLocalElems
    integer :: nLocalNodes3D
    integer :: nLocalSides3D
    integer :: nUniqueSides3D
    integer :: nLocalNodes2D
    integer :: nLocalSides2D
    integer :: nUniqueSides2D
    integer :: nGeo,nBCs
    integer :: eid,lsid,iSide
    integer :: i,j,nid
    integer,dimension(:,:),allocatable :: hopr_elemInfo
    integer,dimension(:,:),allocatable :: hopr_sideInfo
    real(prec),dimension(:,:),allocatable :: hopr_nodeCoords
    integer,dimension(:),allocatable :: hopr_globalNodeIDs
    integer,dimension(:,:),allocatable :: bcType

    if(decomp%mpiEnabled) then
      call Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId,decomp%mpiComm)
    else
      call Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId)
    endif

    call ReadAttribute_HDF5(fileId,'nElems',nGlobalElem)
    call ReadAttribute_HDF5(fileId,'Ngeo',nGeo)
    call ReadAttribute_HDF5(fileId,'nBCs',nBCs)
    call ReadAttribute_HDF5(fileId,'nUniqueSides',nUniqueSides3D)

    ! Read BCType
    allocate(bcType(1:4,1:nBCS))

    if(decomp%mpiEnabled) then
      offset(:) = 0
      call ReadArray_HDF5(fileId,'BCType',bcType,offset)
    else
      call ReadArray_HDF5(fileId,'BCType',bcType)
    endif

    ! Read local subarray of ElemInfo
    call decomp%GenerateDecomposition(nGlobalElem,nUniqueSides3D)

    firstElem = decomp%offsetElem(decomp%rankId+1)+1
    nLocalElems = decomp%offsetElem(decomp%rankId+2)- &
                  decomp%offsetElem(decomp%rankId+1)

    ! Allocate Space for hopr_elemInfo!
    allocate(hopr_elemInfo(1:6,1:nLocalElems))

    if(decomp%mpiEnabled) then
      offset = (/0,firstElem-1/)
      call ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo,offset)
    else
      call ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo)
    endif

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode = hopr_elemInfo(5,1)+1
    nLocalNodes3D = hopr_elemInfo(6,nLocalElems)-hopr_elemInfo(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !
    allocate(hopr_nodeCoords(1:3,nLocalNodes3D),hopr_globalNodeIDs(1:nLocalNodes3D))

    if(decomp%mpiEnabled) then
      offset = (/0,firstNode-1/)
      call ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords,offset)
      gOffset = (/firstNode-1/)
      call ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs,gOffset)
    else
      call ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords)
      call ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs)
    endif

    ! Read local subarray of SideInfo
    firstSide = hopr_elemInfo(3,1)+1
    nLocalSides3D = hopr_elemInfo(4,nLocalElems)-hopr_elemInfo(3,1)

    ! Allocate space for hopr_sideInfo
    allocate(hopr_sideInfo(1:5,1:nLocalSides3D))
    if(decomp%mpiEnabled) then
      offset = (/0,firstSide-1/)
      call ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo,offset)
    else
      call ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo)
    endif

    call Close_HDF5(fileID)
    ! ---- Done reading 3-D Mesh information ---- !

    ! Now we need to convert from 3-D to 2-D !
    nLocalSides2D = nLocalSides3D-2*nGlobalElem
    nUniqueSides2D = nUniqueSides3D-2*nGlobalElem ! Remove the "top" and "bottom" faces
    nLocalNodes2D = nLocalNodes2D-nGlobalElem*nGeo*(nGeo+1)**2 ! Remove the third dimension

    call this%Init(nGeo,nLocalElems,nLocalSides2D,nLocalNodes2D,nBCs)

    ! Copy data from local arrays into this
    !  elemInfo(1:6,iEl)
    !    1 - Element Type
    !    2 - Zone
    !    3 - offset index for side array (not needed when all quads are assumed)
    !    4 - last index for side array (not needed when all quads are assumed)
    !    5 - offset index for node array (not needed when all quads are assumed)
    !    6 - last index for node array (not needed when all quads are assumed)
    this%elemInfo = hopr_elemInfo
    this%quadrature = UNIFORM ! HOPr uses uniformly spaced points

    ! Grab the node coordinates (x and y only) from the "bottom" layer of the extruded mesh
    do eid = 1,this%nElem
      do j = 1,nGeo+1
        do i = 1,nGeo+1
          nid = i+(nGeo+1)*(j-1+(nGeo+1)*((nGeo+1)*(eid-1)))
          this%nodeCoords(1:2,i,j,eid) = hopr_nodeCoords(1:2,nid)
          this%globalNodeIDs(i,j,eid) = hopr_globalNodeIDs(nid)
        enddo
      enddo
    enddo

    ! Grab the south, west, north, and south sides of the elements
    !  sideInfo(1:5,iSide,iEl)
    !
    !    1 - Side Type (currently unused in SELF)
    !    2 - Global Side ID (Used for message passing. Don't need to change)
    !    3 - Neighbor Element ID (Can stay the same)
    !    4 - 10*( neighbor local side )  + flip (Need to recalculate flip)
    !    5 - Boundary Condition ID (Can stay the same)
    do eid = 1,this%nElem
      do lsid = 1,4
        ! Calculate the 3-D side ID from the 2-D local side id and element ID
        iSide = lsid+1+6*(eid-1)
        this%sideInfo(1:5,lsid,eid) = hopr_sideInfo(1:5,iSide)
        ! Adjust the secondary side index for 2-D
        this%sideInfo(4,lsid,eid) = this%sideInfo(4,lsid,eid)-10
      enddo
    enddo

    call this%RecalculateFlip()

    deallocate(hopr_elemInfo,hopr_nodeCoords,hopr_globalNodeIDs,hopr_sideInfo)

  endsubroutine Read_HOPr_Mesh2D

  subroutine RecalculateFlip_Mesh2D(this,decomp)
    implicit none
    class(Mesh2D),intent(inout) :: this
    type(MPILayer),intent(inout),optional :: decomp
    ! Local
    integer :: e1
    integer :: s1
    integer :: e2
    integer :: e2Global
    integer :: s2
    integer :: flip
    integer :: bcid
    integer :: lnid1(1:2)
    integer :: lnid2(1:2)
    integer :: nid1(1:2,1:4,1:this%nElem)
    integer :: nid2(1:2,1:4,1:this%nElem)
    integer :: nloc1(1:2)
    integer :: nloc2(1:2)
    integer :: n1
    integer :: n1Global
    integer :: n2
    integer :: n2Global
    integer :: c1
    integer :: c2
    integer :: i,j
    integer :: l
    integer :: nShifts
    integer :: neighborRank
    integer :: rankId
    integer :: offset
    integer :: msgCount
    integer :: globalSideId
    integer,allocatable :: requests(:)
    integer,allocatable :: stats(:,:)
    integer :: iError
    logical :: theyMatch

    allocate(requests(1:this%nSides*2))
    allocate(stats(MPI_STATUS_SIZE,1:this%nSides*2))

    if(present(decomp)) then
      rankId = decomp%rankId
      offset = decomp%offsetElem(rankId+1)
    else
      rankId = 0
      offset = 0
    endif

    msgCount = 0
    do e1 = 1,this%nElem
      do s1 = 1,4

        e2Global = this%sideInfo(3,s1,e1)
        e2 = e2Global-offset
        s2 = this%sideInfo(4,s1,e1)/10
        flip = this%sideInfo(4,s1,e1)-s2*10
        bcid = this%sideInfo(5,s1,e1)

        if(bcid == 0) then

          if(present(decomp)) then
            neighborRank = decomp%elemToRank(e2Global)
          else
            neighborRank = 0
          endif

          if(neighborRank == rankId) then

            lnid1 = this%CGNSSideMap(1:2,s1) ! local CGNS corner node ids for element 1 side
            lnid2 = this%CGNSSideMap(1:2,s2) ! local CGNS corner node ids for element 2 side

            do l = 1,2

              i = this%CGNSCornerMap(1,lnid1(l))
              j = this%CGNSCornerMap(2,lnid1(l))
              nid1(l,s1,e1) = this%globalNodeIDs(i,j,e1)

              i = this%CGNSCornerMap(1,lnid2(l))
              j = this%CGNSCornerMap(2,lnid2(l))
              nid2(l,s1,e1) = this%globalNodeIDs(i,j,e2)

            enddo

          else ! In this case, we need to exchange

            globalSideId = abs(this%sideInfo(2,s1,e1))

            lnid1 = this%CGNSSideMap(1:2,s1) ! local CGNS corner node ids for element 1 side

            do l = 1,2

              i = this%CGNSCornerMap(1,lnid1(l))
              j = this%CGNSCornerMap(2,lnid1(l))
              nid1(l,s1,e1) = this%globalNodeIDs(i,j,e1)

              msgCount = msgCount+1
              call MPI_IRECV(nid2(l,s1,e1), &
                             1, &
                             MPI_INTEGER, &
                             neighborRank,globalSideId, &
                             decomp%mpiComm, &
                             requests(msgCount),iError)

              ! Send nid1(l) from this rank to nid2(l) on the other rank
              msgCount = msgCount+1
              call MPI_ISEND(nid1(l,s1,e1), &
                             1, &
                             MPI_INTEGER, &
                             neighborRank,globalSideId, &
                             decomp%mpiComm, &
                             requests(msgCount),iError)

            enddo

          endif ! MPI or not

        endif ! If not physical boundary

      enddo
    enddo

    if(present(decomp) .and. msgCount > 0) then
      call MPI_WaitAll(msgCount, &
                       requests(1:msgCount), &
                       stats(1:MPI_STATUS_SIZE,1:msgCount), &
                       iError)
    endif

    do e1 = 1,this%nElem
      do s1 = 1,4

        s2 = this%sideInfo(4,s1,e1)/10
        bcid = this%sideInfo(5,s1,e1)
        nloc1(1:2) = nid1(1:2,s1,e1)
        nloc2(1:2) = nid2(1:2,s1,e1)

        if(bcid == 0) then
          theyMatch = CompareArray(nloc1,nloc2,2)

          if(theyMatch) then
            this%sideInfo(4,s1,e1) = 10*s2
          else
            this%sideInfo(4,s1,e1) = 10*s2+1
          endif

        endif

      enddo
    enddo

    deallocate(requests)
    deallocate(stats)

  endsubroutine RecalculateFlip_Mesh2D

  subroutine Write_Mesh2D(this,meshFile)
    ! Writes mesh output in HOPR format (serial only)
    implicit none
    class(Mesh2D),intent(inout) :: this
    character(*),intent(in) :: meshFile
    ! Local
    integer(HID_T) :: fileId

    call Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId)
    call WriteAttribute_HDF5(fileId,'nElems',this%nElem)
    call WriteAttribute_HDF5(fileId,'Ngeo',this%nGeo)
    call WriteAttribute_HDF5(fileId,'nBCs',this%nBCs)

    call WriteArray_HDF5(fileId,'BCType',this%bcType)

    ! Write local subarray of ElemInfo
    call WriteArray_HDF5(fileId,'ElemInfo',this%elemInfo)

    ! Write local subarray of NodeCoords and GlobalNodeIDs
    call WriteArray_HDF5(fileId,'NodeCoords',this%nodeCoords)
    call WriteArray_HDF5(fileId,'GlobalNodeIDs',this%globalNodeIDs)

    ! Write local subarray of SideInfo
    call WriteArray_HDF5(fileId,'SideInfo',this%sideInfo)

    call Close_HDF5(fileID)

  endsubroutine Write_Mesh2D

  subroutine Init_Mesh3D(this,nGeo,nElem,nSides,nNodes,nBCs)
    implicit none
    class(Mesh3D),intent(out) :: this
    integer,intent(in) :: nGeo
    integer,intent(in) :: nElem
    integer,intent(in) :: nSides
    integer,intent(in) :: nNodes
    integer,intent(in) :: nBCs
    ! Local
    integer :: i,j,k,l

    this%nElem = nElem
    this%nGlobalElem = nElem
    this%nGeo = nGeo
    this%nSides = nSides
    this%nNodes = nNodes
    this%nCornerNodes = 0
    this%nUniqueSides = 0
    this%nUniqueNodes = 0
    this%nBCs = nBCs

    allocate(this%elemInfo(1:6,1:nElem))
    allocate(this%sideInfo(1:5,1:6,1:nElem))
    allocate(this%nodeCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nElem))
    allocate(this%globalNodeIDs(1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nElem))
    allocate(this%CGNSCornerMap(1:3,1:8))
    allocate(this%CGNSSideMap(1:4,1:6))
    allocate(this%sideMap(1:4,1:6))
    allocate(this%BCType(1:4,1:nBCs))

    !$omp target enter data map(alloc: this % elemInfo)
    !$omp target enter data map(alloc: this % sideInfo)
    !$omp target enter data map(alloc: this % nodeCoords)
    !$omp target enter data map(alloc: this % globalNodeIDs)
    !$omp target enter data map(alloc: this % CGNSCornerMap)
    !$omp target enter data map(alloc: this % CGNSSideMap)
    !$omp target enter data map(alloc: this % sideMap)
    !$omp target enter data map(alloc: this % BCType)

    allocate(this%BCNames(1:nBCs))

    ! Create lookup tables to assist with connectivity generation
    this%CGNSCornerMap(1:3,1) = (/1,1,1/) ! Bottom-South-West
    this%CGNSCornerMap(1:3,2) = (/nGeo+1,1,1/) ! Bottom-South-East
    this%CGNSCornerMap(1:3,3) = (/nGeo+1,nGeo+1,1/) ! Bottom-North-East
    this%CGNSCornerMap(1:3,4) = (/1,nGeo+1,1/) ! Bottom-North-West
    this%CGNSCornerMap(1:3,5) = (/1,1,nGeo+1/) ! Top-South-West
    this%CGNSCornerMap(1:3,6) = (/nGeo+1,1,nGeo+1/) ! Top-South-East
    this%CGNSCornerMap(1:3,7) = (/nGeo+1,nGeo+1,nGeo+1/) ! Top-North-East
    this%CGNSCornerMap(1:3,8) = (/1,nGeo+1,nGeo+1/) ! Top-North-West

    ! Maps from local corner node id to CGNS side
    this%CGNSSideMap(1:4,1) = (/1,4,3,2/)
    this%CGNSSideMap(1:4,2) = (/1,2,6,5/)
    this%CGNSSideMap(1:4,3) = (/2,3,7,6/)
    this%CGNSSideMap(1:4,4) = (/3,4,8,7/)
    this%CGNSSideMap(1:4,5) = (/1,5,8,4/)
    this%CGNSSideMap(1:4,6) = (/5,6,7,8/)

    this%sideMap(1:4,1) = (/1,2,3,4/) ! Bottom
    this%sideMap(1:4,2) = (/1,2,6,5/) ! South
    this%sideMap(1:4,3) = (/2,3,7,6/) ! East
    this%sideMap(1:4,4) = (/4,3,7,8/) ! North
    this%sideMap(1:4,5) = (/1,4,8,5/) ! West
    this%sideMap(1:4,6) = (/5,6,7,8/) ! Top

    !$omp target update to(this % CGNSCornerMap)
    !$omp target update to(this % CGNSSideMap)
    !$omp target update to(this % sideMap)

  endsubroutine Init_Mesh3D

  subroutine Free_Mesh3D(this)
    implicit none
    class(Mesh3D),intent(inout) :: this

    this%nElem = 0
    this%nSides = 0
    this%nNodes = 0
    this%nCornerNodes = 0
    this%nUniqueSides = 0
    this%nUniqueNodes = 0
    this%nBCs = 0

    deallocate(this%elemInfo)
    deallocate(this%sideInfo)
    deallocate(this%nodeCoords)
    deallocate(this%globalNodeIDs)
    deallocate(this%CGNSCornerMap)
    deallocate(this%sideMap)
    deallocate(this%CGNSSideMap)
    deallocate(this%BCType)

    !$omp target exit data map(delete: this % elemInfo)
    !$omp target exit data map(delete: this % sideInfo)
    !$omp target exit data map(delete: this % nodeCoords)
    !$omp target exit data map(delete: this % globalNodeIDs)
    !$omp target exit data map(delete: this % CGNSCornerMap)
    !$omp target exit data map(delete: this % CGNSSideMap)
    !$omp target exit data map(delete: this % sideMap)
    !$omp target exit data map(delete: this % BCType)

    deallocate(this%BCNames)

  endsubroutine Free_Mesh3D

  subroutine ResetBoundaryConditionType_Mesh3D(this,bcid)
    !! This method can be used to reset all of the boundary elements
    !! boundary condition type to the desired value.
    !!
    !! Note that ALL physical boundaries will be set to have this boundary
    !! condition
    implicit none
    class(Mesh3D),intent(inout) :: this
    integer,intent(in) :: bcid
    ! Local
    integer :: iSide,iEl,e2

    do iEl = 1,this%nElem
      do iSide = 1,6

        e2 = this%sideInfo(3,iSide,iEl)

        if(e2 == 0) then
          this%sideInfo(5,iSide,iEl) = bcid
        endif

      enddo
    enddo

  endsubroutine ResetBoundaryConditionType_Mesh3D

  subroutine RecalculateFlip_Mesh3D(this,decomp)
    implicit none
    class(Mesh3D),intent(inout) :: this
    type(MPILayer),intent(inout),optional :: decomp
    ! Local
    integer :: e1
    integer :: s1
    integer :: e2
    integer :: e2Global
    integer :: s2
    integer :: flip
    integer :: bcid
    integer :: lnid1(1:4)
    integer :: lnid2(1:4)
    integer :: nid1(1:4,1:6,1:this%nElem)
    integer :: nid2(1:4,1:6,1:this%nElem)
    integer :: nloc1(1:4)
    integer :: nloc2(1:4)
    integer :: n1
    integer :: n1Global
    integer :: n2
    integer :: n2Global
    integer :: c1
    integer :: c2
    integer :: i,j,k
    integer :: l
    integer :: nShifts
    integer :: neighborRank
    integer :: rankId
    integer :: offset
    integer :: msgCount
    integer :: globalSideId
    integer,allocatable :: requests(:)
    integer,allocatable :: stats(:,:)
    integer :: iError
    logical :: theyMatch

    allocate(requests(1:this%nSides*2))
    allocate(stats(MPI_STATUS_SIZE,1:this%nSides*2))

    if(present(decomp)) then
      rankId = decomp%rankId
      offset = decomp%offsetElem(rankId+1)
    else
      rankId = 0
      offset = 0
    endif

    msgCount = 0
    do e1 = 1,this%nElem
      do s1 = 1,6

        e2Global = this%sideInfo(3,s1,e1)
        e2 = e2Global-offset
        s2 = this%sideInfo(4,s1,e1)/10
        flip = this%sideInfo(4,s1,e1)-s2*10
        bcid = this%sideInfo(5,s1,e1)

        if(bcid == 0) then

          if(present(decomp)) then
            neighborRank = decomp%elemToRank(e2Global)
          else
            neighborRank = 0
          endif

          if(neighborRank == rankId) then

            lnid1 = this%sideMap(1:4,s1) ! local CGNS corner node ids for element 1 side
            lnid2 = this%sideMap(1:4,s2) ! local CGNS corner node ids for element 2 side

            do l = 1,4

              i = this%CGNSCornerMap(1,lnid1(l))
              j = this%CGNSCornerMap(2,lnid1(l))
              k = this%CGNSCornerMap(3,lnid1(l))
              nid1(l,s1,e1) = this%globalNodeIDs(i,j,k,e1)

              i = this%CGNSCornerMap(1,lnid2(l))
              j = this%CGNSCornerMap(2,lnid2(l))
              k = this%CGNSCornerMap(3,lnid2(l))
              nid2(l,s1,e1) = this%globalNodeIDs(i,j,k,e1)

            enddo

          else ! In this case, we need to exchange

            globalSideId = abs(this%sideInfo(2,s1,e1))

            lnid1 = this%sideMap(1:4,s1) ! local CGNS corner node ids for element 1 side

            do l = 1,4

              i = this%CGNSCornerMap(1,lnid1(l))
              j = this%CGNSCornerMap(2,lnid1(l))
              k = this%CGNSCornerMap(3,lnid1(l))
              nid1(l,s1,e1) = this%globalNodeIDs(i,j,k,e1)

              ! Receive nid2(l) on this rank from  nid1(l) on the other rank
              msgCount = msgCount+1
              call MPI_IRECV(nid2(l,s1,e1), &
                             1, &
                             MPI_INTEGER, &
                             neighborRank,globalSideId, &
                             decomp%mpiComm, &
                             requests(msgCount),iError)

              ! Send nid1(l) from this rank to nid2(l) on the other rank
              msgCount = msgCount+1
              call MPI_ISEND(nid1(l,s1,e1), &
                             1, &
                             MPI_INTEGER, &
                             neighborRank,globalSideId, &
                             decomp%mpiComm, &
                             requests(msgCount),iError)

            enddo

          endif ! MPI or not

        endif ! If not physical boundary

      enddo
    enddo

    if(present(decomp) .and. msgCount > 0) then
      call MPI_WaitAll(msgCount, &
                       requests(1:msgCount), &
                       stats(1:MPI_STATUS_SIZE,1:msgCount), &
                       iError)
    endif

    do e1 = 1,this%nElem
      do s1 = 1,6

        s2 = this%sideInfo(4,s1,e1)/10
        bcid = this%sideInfo(5,s1,e1)
        nloc1(1:4) = nid1(1:4,s1,e1)
        nloc2(1:4) = nid2(1:4,s1,e1)

        if(bcid == 0) then
          nShifts = 0
          theyMatch = .false.

          do i = 1,4

            theyMatch = CompareArray(nloc1,nloc2,4)

            if(theyMatch) then
              exit
            else
              nShifts = nShifts+1
              call ForwardShift(nloc1,4)
            endif

          enddo

          this%sideInfo(4,s1,e1) = 10*s2+nShifts

        endif

      enddo
    enddo

    deallocate(requests)
    deallocate(stats)

  endsubroutine RecalculateFlip_Mesh3D

  subroutine Read_HOPr_Mesh3D(this,meshFile,decomp)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    implicit none
    class(Mesh3D),intent(out) :: this
    character(*),intent(in) :: meshFile
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer(HID_T) :: fileId
    integer(HID_T) :: offset(1:2),gOffset(1)
    integer :: nGlobalElem
    integer :: firstElem
    integer :: firstNode
    integer :: firstSide
    integer :: nLocalElems
    integer :: nLocalNodes
    integer :: nLocalSides
    integer :: nUniqueSides
    integer :: nGeo,nBCs
    integer :: eid,lsid,iSide
    integer :: i,j,k,nid
    integer,dimension(:,:),allocatable :: hopr_elemInfo
    integer,dimension(:,:),allocatable :: hopr_sideInfo
    real(prec),dimension(:,:),allocatable :: hopr_nodeCoords
    integer,dimension(:),allocatable :: hopr_globalNodeIDs
    integer,dimension(:,:),allocatable :: bcType

    if(decomp%mpiEnabled) then
      call Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId,decomp%mpiComm)
    else
      call Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId)
    endif

    call ReadAttribute_HDF5(fileId,'nElems',nGlobalElem)
    call ReadAttribute_HDF5(fileId,'Ngeo',nGeo)
    call ReadAttribute_HDF5(fileId,'nBCs',nBCs)
    call ReadAttribute_HDF5(fileId,'nUniqueSides',nUniqueSides)

    ! Read BCType
    allocate(bcType(1:4,1:nBCs))
    if(decomp%mpiEnabled) then
      offset(:) = 0
      call ReadArray_HDF5(fileId,'BCType',bcType,offset)
    else
      call ReadArray_HDF5(fileId,'BCType',bcType)
    endif

    ! Read local subarray of ElemInfo
    call decomp%GenerateDecomposition(nGlobalElem,nUniqueSides)

    firstElem = decomp%offsetElem(decomp%rankId+1)+1
    nLocalElems = decomp%offsetElem(decomp%rankId+2)- &
                  decomp%offsetElem(decomp%rankId+1)

    ! Allocate Space for hopr_elemInfo!
    allocate(hopr_elemInfo(1:6,1:nLocalElems))
    if(decomp%mpiEnabled) then
      offset = (/0,firstElem-1/)
      call ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo,offset)
    else
      call ReadArray_HDF5(fileId,'ElemInfo',hopr_elemInfo)
    endif

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    firstNode = hopr_elemInfo(5,1)+1
    nLocalNodes = hopr_elemInfo(6,nLocalElems)-hopr_elemInfo(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !
    allocate(hopr_nodeCoords(1:3,1:nLocalNodes),hopr_globalNodeIDs(1:nLocalNodes))

    if(decomp%mpiEnabled) then
      offset = (/0,firstNode-1/)
      call ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords,offset)
      gOffset = (/firstNode-1/)
      call ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs,gOffset)
    else
      call ReadArray_HDF5(fileId,'NodeCoords',hopr_nodeCoords)
      call ReadArray_HDF5(fileId,'GlobalNodeIDs',hopr_globalNodeIDs)
    endif

    ! Read local subarray of SideInfo
    firstSide = hopr_elemInfo(3,1)+1
    nLocalSides = hopr_elemInfo(4,nLocalElems)-hopr_elemInfo(3,1)

    ! Allocate space for hopr_sideInfo
    allocate(hopr_sideInfo(1:5,1:nLocalSides))

    if(decomp%mpiEnabled) then
      offset = (/0,firstSide-1/)
      call ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo,offset)
    else
      call ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo)
    endif

    call Close_HDF5(fileID)
    ! ---- Done reading 3-D Mesh information ---- !
    ! Load hopr data into mesh data structure

    call this%Init(nGeo,nLocalElems,nLocalSides,nLocalNodes,nBCs)

    ! Copy data from local arrays into this
    this%elemInfo = hopr_elemInfo
    this%nUniqueSides = nUniqueSides
    this%quadrature = UNIFORM

    ! Grab the node coordinates
    do eid = 1,this%nElem
      do k = 1,nGeo+1
        do j = 1,nGeo+1
          do i = 1,nGeo+1
            nid = i+(nGeo+1)*(j-1+(nGeo+1)*(k-1+(nGeo+1)*(eid-1)))
            this%nodeCoords(1:3,i,j,k,eid) = hopr_nodeCoords(1:3,nid)
            this%globalNodeIDs(i,j,k,eid) = hopr_globalNodeIDs(nid)
          enddo
        enddo
      enddo
    enddo

    iSide = 0
    do eid = 1,this%nElem
      do lsid = 1,6
        iSide = iSide+1
        this%sideInfo(1:5,lsid,eid) = hopr_sideInfo(1:5,iSide)
      enddo
    enddo

    call this%RecalculateFlip()

    deallocate(hopr_elemInfo,hopr_nodeCoords,hopr_globalNodeIDs,hopr_sideInfo)

  endsubroutine Read_HOPr_Mesh3D

  subroutine Write_Mesh3D(this,meshFile)
    ! Writes mesh output in HOPR format (serial only)
    implicit none
    class(Mesh3D),intent(inout) :: this
    character(*),intent(in) :: meshFile
    ! Local
    integer(HID_T) :: fileId

    call Open_HDF5(meshFile,H5F_ACC_RDWR_F,fileId)

    call WriteAttribute_HDF5(fileId,'nElems',this%nElem)
    call WriteAttribute_HDF5(fileId,'Ngeo',this%nGeo)
    call WriteAttribute_HDF5(fileId,'nBCs',this%nBCs)

    call WriteArray_HDF5(fileId,'BCType',this%bcType)
    call WriteArray_HDF5(fileId,'ElemInfo',this%elemInfo)

    ! Read local subarray of NodeCoords and GlobalNodeIDs
    call WriteArray_HDF5(fileId,'NodeCoords',this%nodeCoords)
    call WriteArray_HDF5(fileId,'GlobalNodeIDs',this%globalNodeIDs)

    ! Read local subarray of SideInfo
    call WriteArray_HDF5(fileId,'SideInfo',this%sideInfo)

    call Close_HDF5(fileID)

  endsubroutine Write_Mesh3D

  subroutine Init_MPILayer(this,enableMPI)
#undef __FUNC__
#define __FUNC__ "Init_MPILayer"
    implicit none
    class(MPILayer),intent(out) :: this
    logical,intent(in) :: enableMPI
    ! Local
    integer       :: ierror
    character(50) :: msg
    integer :: nGPU,gpuID
    character(2) :: msg2

    this%mpiComm = 0
    this%mpiPrec = prec
    this%rankId = 0
    this%nRanks = 1
    this%nElem = 0
    this%mpiEnabled = enableMPI

    if(enableMPI) then
      this%mpiComm = MPI_COMM_WORLD
      call MPI_INIT(ierror)
      call MPI_COMM_RANK(this%mpiComm,this%rankId,ierror)
      call MPI_COMM_SIZE(this%mpiComm,this%nRanks,ierror)
    endif

    if(prec == real32) then
      this%mpiPrec = MPI_FLOAT
    else
      this%mpiPrec = MPI_DOUBLE
    endif

    allocate(this%offsetElem(1:this%nRanks+1))
    !$omp target enter data map(alloc: this % offsetElem)

    write(msg,'(I5)') this%rankId
    msg = "Greetings from rank "//trim(msg)//"."
    INFO(trim(msg))

    ! ! Get the number of GPUs per node
    ! CALL hipCheck(hipGetDeviceCount(nGPU))

    ! ! Assume that we have the 1 GPU per rank
    ! ! implying that nMPIRanksPerNode = nGPU
    ! ! Assume that all nodes have the same number of GPUs per node
    ! gpuID = MOD(this % rankId, nGPU)

    ! CALL hipCheck(hipSetDevice(gpuID))
    ! WRITE (msg,'(I5)') this % rankId
    !WRITE (msg2,'(I2)') gpuID
    !msg = "Rank "//TRIM(msg)//": Setting device to GPU"//TRIM(msg2)
    !NFO(TRIM(msg))

  endsubroutine Init_MPILayer

  subroutine Free_MPILayer(this)
    implicit none
    class(MPILayer),intent(inout) :: this

    if(associated(this%offSetElem)) then
      deallocate(this%offSetElem)
      !$omp target exit data map(delete: this % offsetElem)
    endif
    if(associated(this%elemToRank)) then
      deallocate(this%elemToRank)
      !$omp target exit data map(delete: this % elemToRank)
    endif

    deallocate(this%requests)
    deallocate(this%stats)

  endsubroutine Free_MPILayer

  subroutine Finalize_MPILayer(this)
#undef __FUNC__
#define __FUNC__ "Finalize_MPILayer"
    implicit none
    class(MPILayer),intent(inout) :: this
    ! Local
    integer       :: ierror
    character(30) :: msg

    if(this%mpiEnabled) then
      write(msg,'(I5)') this%rankId
      msg = "Goodbye from rank "//trim(msg)//"."
      INFO(trim(msg))
      call MPI_FINALIZE(ierror)
    endif

  endsubroutine Finalize_MPILayer

  subroutine GenerateDecomposition_MPILayer(this,nGlobalElem,maxMsg)
#undef __FUNC__
#define __FUNC__ "GenerateDecomposition_MPILayer"
    implicit none
    class(MPILayer),intent(inout) :: this
    integer,intent(in) :: nGlobalElem
    integer,intent(in) :: maxMsg
    ! Local
    integer :: maxMsgLoc
    character(50) :: msg
    character(5) :: msg2

    call this%setElemToRank(nGlobalElem)
    call this%SetMaxMsg(maxMsg)

    write(msg,'(I5)') this%rankId
    write(msg2,'(I5)') this%offSetElem(this%rankId+2)- &
      this%offSetElem(this%rankId+1)
    msg = "Rank "//trim(msg)//": nElem = "//trim(msg2)
    INFO(trim(msg))

  endsubroutine GenerateDecomposition_MPILayer

  subroutine SetMaxMsg(this,maxMsg)
    implicit none
    class(MPILayer),intent(inout) :: this
    integer,intent(in) :: maxMsg

    if(allocated(this%requests)) deallocate(this%requests)
    if(allocated(this%stats)) deallocate(this%stats)

    allocate(this%requests(1:maxMsg))
    allocate(this%stats(MPI_STATUS_SIZE,1:maxMsg))
    this%maxMsg = maxMsg

  endsubroutine SetMaxMsg

  subroutine SetElemToRank(this,nElem)
    implicit none
    class(MPILayer),intent(inout) :: this
    integer,intent(in) :: nElem
    ! Local
    integer :: iel

    this%nElem = nElem

    allocate(this%elemToRank(1:nelem))
    !$omp target enter data map(alloc:this % elemToRank)

    call DomainDecomp(nElem, &
                      this%nRanks, &
                      this%offSetElem)

    do iel = 1,nElem
      call ElemToRank(this%nRanks, &
                      this%offSetElem, &
                      iel, &
                      this%elemToRank(iel))
    enddo

    !$omp target update to(this % offsetElem)
    !$omp target update to(this % elemToRank)

  endsubroutine SetElemToRank

  subroutine DomainDecomp(nElems,nDomains,offSetElem)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 4
    implicit none
    integer,intent(in) :: nElems
    integer,intent(in) :: nDomains
    integer,intent(out) :: offsetElem(0:nDomains)
    ! Local
    integer :: nLocalElems
    integer :: remainElems
    integer :: iDom

    nLocalElems = nElems/nDomains
    remainElems = nElems-nLocalElems*nDomains
    do iDom = 0,nDomains-1
      offSetElem(iDom) = iDom*nLocalElems+min(iDom,remainElems)
    enddo
    offSetElem(nDomains) = nElems

  endsubroutine DomainDecomp

  subroutine ElemToRank(nDomains,offsetElem,elemID,domain)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 7
    !   "Find domain containing element index"
    !
    implicit none
    integer,intent(in) :: nDomains
    integer,intent(in) :: offsetElem(0:nDomains)
    integer,intent(in) :: elemID
    integer,intent(out) :: domain
    ! Local
    integer :: maxSteps
    integer :: low,up,mid
    integer :: i

    domain = 0
    maxSteps = int(log10(real(nDomains))/log10(2.0))+1
    low = 0
    up = nDomains-1

    if(offsetElem(low) < elemID .and. elemID <= offsetElem(low+1)) then
      domain = low
    elseif(offsetElem(up) < elemID .and. elemID <= offsetElem(up+1)) then
      domain = up
    else
      do i = 1,maxSteps
        mid = (up-low)/2+low
        if(offsetElem(mid) < elemID .and. elemID <= offsetElem(mid+1)) then
          domain = mid
          return
        elseif(elemID > offsetElem(mid+1)) then
          low = mid+1
        else
          up = mid
        endif
      enddo
    endif

  endsubroutine ElemToRank

  subroutine FinalizeMPIExchangeAsync(mpiHandler)
    class(MPILayer),intent(inout) :: mpiHandler
    ! Local
    integer :: ierror

    if(mpiHandler%mpiEnabled) then
      call MPI_WaitAll(mpiHandler%msgCount, &
                       mpiHandler%requests(1:mpiHandler%msgCount), &
                       mpiHandler%stats(1:MPI_STATUS_SIZE,1:mpiHandler%msgCount), &
                       iError)
    endif

  endsubroutine FinalizeMPIExchangeAsync

  subroutine GlobalReduce_RealScalar(mpiHandler,sendBuf,recvBuf)
    class(MPILayer),intent(in) :: mpiHandler
    real(prec),intent(in) :: sendBuf
    real(prec),intent(out) :: recvBuf
    ! Local
    integer :: iError

    if(mpiHandler%mpiEnabled) then
      call MPI_ALLREDUCE(sendBuf, &
                         recvBuf, &
                         1, &
                         mpiHandler%mpiPrec, &
                         MPI_SUM, &
                         mpiHandler%mpiComm, &
                         iError)
    else
      recvBuf = sendBuf
    endif

  endsubroutine GlobalReduce_RealScalar

endmodule SELF_Mesh
