! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2024 Fluid Numerics LLC
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
!    the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Mesh_2D_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_SupportRoutines
  use SELF_HDF5
  use SELF_Mesh
  use SELF_DomainDecomposition

  ! External Libs !
  use HDF5

  use iso_c_binding

  implicit none

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
!
! 2-D Hexahedreal Element sides are defined as
!
! Side 1 = South  (xi2 = -1) = [CN1, CN2]
! Side 2 = East   (xi1 = 1) = [CN2, CN3]
! Side 3 = North  (xi2 = 1) = [CN4, CN3]
! Side 4 = West   (xi1 = -1) = [CN1, CN4]
!
! In 2-D, corner nodes are order counter-clockwise (looking in the -xi3 direction).
!
! CornerNode 1 = South-West = (-1,-1)
! CornerNode 2 = South-East = ( 1,-1)
! CornerNode 3 = North-East = ( 1, 1)
! CornerNode 4 = North-West = (-1, 1)
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

  ! Side Ordering
  integer,parameter :: selfSide2D_South = 1
  integer,parameter :: selfSide2D_East = 2
  integer,parameter :: selfSide2D_North = 3
  integer,parameter :: selfSide2D_West = 4

  ! Mesh format is set up similar to the HOPr format
  ! See https://hopr-project.org/externals/MeshFormat.pdf

  type,extends(SEMMesh) :: Mesh2D_t
    integer,pointer,dimension(:,:,:) :: sideInfo
    real(prec),pointer,dimension(:,:,:,:) :: nodeCoords
    integer,pointer,dimension(:,:) :: elemInfo
    integer,pointer,dimension(:,:,:) :: globalNodeIDs
    integer,pointer,dimension(:,:) :: CGNSCornerMap
    integer,pointer,dimension(:,:) :: CGNSSideMap
    integer,pointer,dimension(:,:) :: BCType
    character(LEN=255),allocatable :: BCNames(:)

  contains
    procedure,public :: Init => Init_Mesh2D_t
    procedure,public :: Free => Free_Mesh2D_t
    procedure,public :: UpdateDevice => UpdateDevice_Mesh2D_t

    generic,public :: StructuredMesh => UniformStructuredMesh_Mesh2D_t
    procedure,private :: UniformStructuredMesh_Mesh2D_t
    procedure,public :: ResetBoundaryConditionType => ResetBoundaryConditionType_Mesh2D_t

    procedure,public :: Read_HOPr => Read_HOPr_Mesh2D_t

    procedure,public :: Write_Mesh => Write_Mesh2D_t

    procedure,public :: RecalculateFlip => RecalculateFlip_Mesh2D_t

  endtype Mesh2D_t

contains

  subroutine Init_Mesh2D_t(this,nGeo,nElem,nSides,nNodes,nBCs)
    implicit none
    class(Mesh2D_t),intent(inout) :: this
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

  endsubroutine Init_Mesh2D_t

  subroutine Free_Mesh2D_t(this)
    implicit none
    class(Mesh2D_t),intent(inout) :: this

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
    call this%decomp%Free()

  endsubroutine Free_Mesh2D_t

  subroutine UpdateDevice_Mesh2D_t(this)
    implicit none
    class(Mesh2D_t),intent(inout) :: this

    return

  endsubroutine UpdateDevice_Mesh2D_t

  subroutine ResetBoundaryConditionType_Mesh2D_t(this,bcid)
    !! This method can be used to reset all of the boundary elements
    !! boundary condition type to the desired value.
    !!
    !! Note that ALL physical boundaries will be set to have this boundary
    !! condition
    implicit none
    class(Mesh2D_t),intent(inout) :: this
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

    call this%UpdateDevice()

  endsubroutine ResetBoundaryConditionType_Mesh2D_t

  subroutine UniformStructuredMesh_Mesh2D_t(this,nxPerTile,nyPerTile,nTileX,nTileY,dx,dy,bcids)
  !!
  !! Create a structured mesh and store it in SELF's unstructured mesh format.
  !! The mesh is created in tiles of size (tnx,tny). Tiling is used to determine
  !! the element ordering.
  !!
  !!
  !!  Input
  !!    - this : Fresh/empty mesh2d_t object
  !!    - nxPerTile : The number of elements in the x direction within a tile
  !!    - nyPerTile : The number of elements in the y direction within a tile
  !!    - nTileX : The number of tiles in the x direction
  !!    - nTileY : The number of tiles in the y direction
  !!    - dx : Element width in the x-direction
  !!    - dy : Element width in the y-direction
  !!    - bcids(1:4) : Boundary condition flags for the south, east, north, and west sides of the domain
  !!    - enableDomainDecomposition : Boolean to determine if domain decomposition is used.
  !!
  !!  Output
  !!    - this : mesh2d_t object with vertices, edges, and element information
  !!
  !! Total number of elements in the x-direction is nX = nxPerTile*nTileX
  !! Total number of elements in the y-direction is nY = nyPerTile*nTileY
  !!
  !! Length of the domain in the x-direction is Lx = dx*nX
  !! Length of the domain in the y-direction is Ly = dy*nY
  !!
    implicit none
    class(Mesh2D_t),intent(out) :: this
    integer,intent(in) :: nxPerTile
    integer,intent(in) :: nyPerTile
    integer,intent(in) :: nTileX
    integer,intent(in) :: nTileY
    real(prec),intent(in) :: dx
    real(prec),intent(in) :: dy
    integer,intent(in) :: bcids(1:4)
    ! Local
    integer :: nX,nY,nGeo,nBCs
    integer :: nGlobalElem
    integer :: nUniqueSides
    integer :: nUniqueNodes
    integer :: nLocalElems
    integer :: nLocalSides
    integer :: nLocalNodes
    real(prec),allocatable :: nodeCoords(:,:,:,:)
    integer,allocatable :: globalNodeIDs(:,:,:)
    integer,allocatable :: sideInfo(:,:,:)
    integer :: i,j,ti,tj
    integer :: ix,iy,iel
    integer :: ni,nj
    integer :: e1,e2
    integer :: nedges

    call this%decomp%init()

    nX = nTileX*nxPerTile
    nY = nTileY*nyPerTile
    nGeo = 1 ! Force the geometry to be linear
    nBCs = 4 ! Force the number of boundary conditions to 4

    nGlobalElem = nX*nY
    nUniqueSides = (nX+1)*nY+(nY+1)*nX
    nUniqueNodes = (nX+1)*(nY+1)

    allocate(nodeCoords(1:2,1:nGeo+1,1:nGeo+1,1:nGlobalElem))
    allocate(globalNodeIDs(1:nGeo+1,1:nGeo+1,1:nGlobalElem))
    allocate(sideInfo(1:5,1:4,1:nGlobalElem))

    do tj = 1,nTileY
      do ti = 1,nTileX
        do j = 1,nyPerTile
          iy = j+nyPerTile*(tj-1)
          do i = 1,nxPerTile
            iel = i+nxPerTile*(j-1+nyPerTile*(ti-1+nTilex*(tj-1)))
            ix = i+nxPerTile*(ti-1) ! nxpertile + nxpertile*(nTileX-1) = nxperTile*nTilex = 1
            do nj = 1,nGeo+1
              do ni = 1,nGeo+1
                nodeCoords(1,ni,nj,iel) = real(ni-1+ix-1,prec)*dx
                nodeCoords(2,ni,nj,iel) = real(nj-1+iy-1,prec)*dy
                globalNodeIDs(ni,nj,iel) = ni-1+i+(nxPerTile+1)*( &
                                           nj-1+j-1+(nyPerTile+1)*( &
                                           ti-1+nTileX*(tj-1)))
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    ! Fill in edge information
    !  sideInfo(1:5,iSide,iEl)
    !    1 - Side Type (currently unused in SELF)
    !    2 - Global Side ID (Used for message passing. Don't need to change)
    !    3 - Neighbor Element ID (Can stay the same)
    !    4 - 10*( neighbor local side )  + flip (Need to recalculate flip)
    !    5 - Boundary Condition ID (Can stay the same)
    nedges = 0
    do tj = 1,nTileY
      do ti = 1,nTileX
        do j = 1,nyPerTile
          do i = 1,nxPerTile
            iel = i+nxPerTile*(j-1+nyPerTile*(ti-1+nTilex*(tj-1)))

            ! south, iside=1
            ! Get the corner node ids for this edge
            ! sideInfo(2,1,iel) = (nc1+nc2)*(nc1+nc2+1)/2 + nc2
            if(j == 1) then ! southern most part of the tile
              if(tj == 1) then ! southern most tile
                nedges = nedges+1
                sideinfo(2,1,iel) = nedges
                sideinfo(3,1,iel) = 0 ! Neigbor element (null, boundary condition)
                sideinfo(4,1,iel) = 0 ! Neighbor side id (null, boundary condition)
                sideinfo(5,1,iel) = bcids(1) ! Boundary condition id; set from the user input
              else ! interior tile, but souther most edge of the tile
                e2 = i+nxPerTile*(nyPerTile-1+nyPerTile*(ti-1+nTilex*(tj-2))) ! Neigbor element, northernmost element, in tile to the south
                sideinfo(2,1,iel) = sideInfo(2,3,e2) ! Copy the edge id from neighbor's north edge
                sideinfo(3,1,iel) = e2
                sideinfo(4,1,iel) = 10*3 ! Neighbor side id - neighbor to the south, north side (3)
                sideinfo(5,1,iel) = 0 ! Boundary condition id; (null, interior edge)
              endif
            else ! interior to the tile
              e2 = i+nxPerTile*(j-2+nyPerTile*(ti-1+nTilex*(tj-1))) ! Neigbor element, inside same tile, to the south
              sideinfo(2,1,iel) = sideInfo(2,3,e2) ! Copy the edge id from neighbor's north edge
              sideinfo(3,1,iel) = e2
              sideinfo(4,1,iel) = 10*3 ! Neighbor side id - neighbor to the south, north side (3)
              sideinfo(5,1,iel) = 0 ! Boundary condition id; (null, interior edge)
            endif

            ! east, iside=2
            ! Get the corner node ids for this edge
            ! East edges are always new edges, due to the way we are traversing the grid
            nedges = nedges+1
            sideinfo(2,2,iel) = nedges
            if(i == nxPerTile) then ! eastern most part of the tile
              if(ti == nTileX) then ! eastern most tile
                sideinfo(3,2,iel) = 0 ! Neigbor element (null, boundary condition)
                sideinfo(4,2,iel) = 0 ! Neighbor side id (null, boundary condition)
                sideinfo(5,2,iel) = bcids(2) ! Boundary condition id; eastern boundary set from the user input
              else ! interior tile, but eastern most edge of the tile
                sideinfo(3,2,iel) = 1+nxPerTile*(j-1+nyPerTile*(ti+nTilex*(tj-1))) ! Neigbor element, westernnmost element, in tile to the east
                sideinfo(4,2,iel) = 10*4 ! Neighbor side id - neighbor to the east, west side (4)
                sideinfo(5,2,iel) = 0 ! Boundary condition id; (null, interior edge)
              endif
            else ! interior to the tile
              sideinfo(3,2,iel) = i+1+nxPerTile*(j-1+nyPerTile*(ti-1+nTilex*(tj-1))) ! Neigbor element, inside same tile, to the east
              sideinfo(4,2,iel) = 10*4 ! Neighbor side id - neighbor to the east, west side (4)
              sideinfo(5,2,iel) = 0 ! Boundary condition id; (null, interior edge)
            endif

            ! north, iside=3
            ! Get the corner node ids for this edge
            ! East edges are always new edges, due to the way we are traversing the grid
            nedges = nedges+1
            sideinfo(2,3,iel) = nedges
            if(j == nyPerTile) then ! northern most part of the tile
              if(tj == nTileY) then ! northern most tile
                sideinfo(3,3,iel) = 0 ! Neigbor element (null, boundary condition)
                sideinfo(4,3,iel) = 0 ! Neighbor side id (null, boundary condition)
                sideinfo(5,3,iel) = bcids(3) ! Boundary condition id; set from the user input
              else ! interior tile, but northern most edge of the tile
                sideinfo(3,3,iel) = i+nxPerTile*(nyPerTile*(ti-1+nTilex*(tj))) ! Neigbor element, southernmost element in tile to the north
                sideinfo(4,3,iel) = 10*1 ! Neighbor side id - neighbor to the north, south side (1)
                sideinfo(5,3,iel) = 0 ! Boundary condition id; (null, interior edge)
              endif
            else ! interior to the tile
              sideinfo(3,3,iel) = i+nxPerTile*(j+nyPerTile*(ti-1+nTilex*(tj-1))) ! Neigbor element, inside same tile, to the north
              sideinfo(4,3,iel) = 10*1 ! Neighbor side id - neighbor to the north, south side (1)
              sideinfo(5,3,iel) = 0 ! Boundary condition id; (null, interior edge)
            endif

            ! west, iside=4
            ! Get the corner node ids for this edge
            ! n1 = globalNodeIds(this%CGNSCornerMap(1,1),this%CGNSCornerMap(2,1),iel)
            ! n2 = globalNodeIds(this%CGNSCornerMap(1,4),this%CGNSCornerMap(2,4),iel)
            ! nc1 = min(n1,n2)
            ! nc2 = max(n1,n2)
            ! sideInfo(2,1,iel) = (nc1+nc2)*(nc1+nc2+1)/2 + nc2
            if(i == 1) then ! western most part of the tile
              if(ti == 1) then ! western most tile
                nedges = nedges+1
                sideinfo(2,4,iel) = nedges
                sideinfo(3,4,iel) = 0 ! Neigbor element (null, boundary condition)
                sideinfo(4,4,iel) = 0 ! Neighbor side id (null, boundary condition)
                sideinfo(5,4,iel) = bcids(4) ! Boundary condition id; eastern boundary set from the user input
              else ! interior tile, but western most edge of the tile
                e2 = nxPerTile+nxPerTile*(j-1+nyPerTile*(ti-2+nTilex*(tj-1))) ! Neigbor element, easternnmost element in tile to the west
                sideinfo(3,4,iel) = sideInfo(2,2,e2) ! Copy the edge id from neighbor's east edge
                sideinfo(3,4,iel) = e2
                sideinfo(4,4,iel) = 10*2 ! Neighbor side id - neighbor to the west, east side (2)
                sideinfo(5,4,iel) = 0 ! Boundary condition id; (null, interior edge)
              endif
            else ! interior to the tile
              e2 = i-1+nxPerTile*(j-1+nyPerTile*(ti-1+nTilex*(tj-1))) ! Neigbor element, inside same tile, to the west
              sideinfo(3,4,iel) = sideInfo(2,2,e2) ! Copy the edge id from neighbor's east edge
              sideinfo(3,4,iel) = e2
              sideinfo(4,4,iel) = 10*2 ! Neighbor side id - neighbor to the west, east side (2)
              sideinfo(5,4,iel) = 0 ! Boundary condition id; (null, interior edge)
            endif

          enddo
        enddo
      enddo
    enddo

    call this%decomp%GenerateDecomposition(nGlobalElem,nUniqueSides)

    e1 = this%decomp%offsetElem(this%decomp%rankId+1)+1
    e2 = this%decomp%offsetElem(this%decomp%rankId+2)
    nLocalElems = e2-e1+1

    nLocalSides = nLocalElems*4
    nLocalNodes = nLocalElems*4
    call this%Init(nGeo,nLocalElems,nLocalSides,nLocalNodes,nBCs)
    this%nUniqueSides = nUniqueSides
    this%quadrature = UNIFORM

    this%nodeCoords(1:2,1:nGeo+1,1:nGeo+1,1:nLocalElems) = nodeCoords(1:2,1:nGeo+1,1:nGeo+1,e1:e2)
    this%globalNodeIDs(1:nGeo+1,1:nGeo+1,1:nLocalElems) = globalNodeIDs(1:nGeo+1,1:nGeo+1,e1:e2)
    this%sideInfo(1:5,1:4,1:nLocalElems) = sideInfo(1:5,1:4,e1:e2)

    deallocate(nodeCoords)
    deallocate(globalNodeIDs)
    deallocate(sideInfo)

    call this%UpdateDevice()

  endsubroutine UniformStructuredMesh_Mesh2D_t

  subroutine Read_HOPr_Mesh2D_t(this,meshFile)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    ! Adapted for 2D Mesh : Note that HOPR does not have 2D mesh output.
    implicit none
    class(Mesh2D_t),intent(out) :: this
    character(*),intent(in) :: meshFile
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

    call this%decomp%init()

    print*,__FILE__//' : Reading HOPr mesh from'//trim(meshfile)
    if(this%decomp%mpiEnabled) then
      call Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId,this%decomp%mpiComm)
    else
      call Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId)
    endif

    print*,__FILE__//' : Loading mesh attributes'
    call ReadAttribute_HDF5(fileId,'nElems',nGlobalElem)
    call ReadAttribute_HDF5(fileId,'Ngeo',nGeo)
    call ReadAttribute_HDF5(fileId,'nBCs',nBCs)
    call ReadAttribute_HDF5(fileId,'nUniqueSides',nUniqueSides3D)
    print*,__FILE__//' : N Global Elements = ',nGlobalElem
    print*,__FILE__//' : Mesh geometry degree = ',nGeo
    print*,__FILE__//' : N Boundary conditions = ',nBCs
    print*,__FILE__//' : N Unique Sides (3D) = ',nUniqueSides3D

    ! Read BCType
    allocate(bcType(1:4,1:nBCS))

    if(this%decomp%mpiEnabled) then
      offset(:) = 0
      call ReadArray_HDF5(fileId,'BCType',bcType,offset)
    else
      call ReadArray_HDF5(fileId,'BCType',bcType)
    endif

    ! Read local subarray of ElemInfo
    print*,__FILE__//' : Generating Domain Decomposition'
    call this%decomp%GenerateDecomposition(nGlobalElem,nUniqueSides3D)

    firstElem = this%decomp%offsetElem(this%decomp%rankId+1)+1
    nLocalElems = this%decomp%offsetElem(this%decomp%rankId+2)- &
                  this%decomp%offsetElem(this%decomp%rankId+1)

    print*,__FILE__//' : Rank ',this%decomp%rankId+1,' : element offset = ',firstElem
    print*,__FILE__//' : Rank ',this%decomp%rankId+1,' : n_elements = ',nLocalElems

    ! Allocate Space for hopr_elemInfo!
    allocate(hopr_elemInfo(1:6,1:nLocalElems))

    if(this%decomp%mpiEnabled) then
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

    if(this%decomp%mpiEnabled) then
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
    if(this%decomp%mpiEnabled) then
      offset = (/0,firstSide-1/)
      print*,__FILE__//' : Rank ',this%decomp%rankId+1,' Reading side information'
      call ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo,offset)
    else
      call ReadArray_HDF5(fileId,'SideInfo',hopr_sideInfo)
    endif

    call Close_HDF5(fileID)
    ! ---- Done reading 3-D Mesh information ---- !

    ! Now we need to convert from 3-D to 2-D !
    nLocalSides2D = nLocalSides3D-2*nLocalElems
    nUniqueSides2D = nUniqueSides3D-2*nGlobalElem ! Remove the "top" and "bottom" faces
    nLocalNodes2D = nLocalNodes2D-nLocalElems*nGeo*(nGeo+1)**2 ! Remove the third dimension

    print*,__FILE__//' : Rank ',this%decomp%rankId+1,' Allocating memory for mesh'
    print*,__FILE__//' : Rank ',this%decomp%rankId+1,' n local sides  : ',nLocalSides2D
    call this%Init(nGeo,nLocalElems,nLocalSides2D,nLocalNodes2D,nBCs)
    this%nUniqueSides = nUniqueSides2D ! Store the number of sides in the global mesh

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

    call this%UpdateDevice()

  endsubroutine Read_HOPr_Mesh2D_t

  subroutine RecalculateFlip_Mesh2D_t(this)
    implicit none
    class(Mesh2D_t),intent(inout) :: this
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
    integer :: tag
    logical :: theyMatch

    allocate(requests(1:this%nSides*2))
    allocate(stats(MPI_STATUS_SIZE,1:this%nSides*2))

    if(this%decomp%mpiEnabled) then
      rankId = this%decomp%rankId
      offset = this%decomp%offsetElem(rankId+1)
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

        if(e2Global > 0) then

          if(this%decomp%mpiEnabled) then
            neighborRank = this%decomp%elemToRank(e2Global)
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

              tag = l+2*globalSideId
              msgCount = msgCount+1
              call MPI_IRECV(nid2(l,s1,e1), &
                             1, &
                             MPI_INTEGER, &
                             neighborRank,tag, &
                             this%decomp%mpiComm, &
                             requests(msgCount),iError)

              ! Send nid1(l) from this rank to nid2(l) on the other rank
              msgCount = msgCount+1
              call MPI_ISEND(nid1(l,s1,e1), &
                             1, &
                             MPI_INTEGER, &
                             neighborRank,tag, &
                             this%decomp%mpiComm, &
                             requests(msgCount),iError)

            enddo

          endif ! MPI or not

        endif ! If not physical boundary

      enddo
    enddo

    if(this%decomp%mpiEnabled .and. msgCount > 0) then
      call MPI_WaitAll(msgCount, &
                       requests(1:msgCount), &
                       stats(1:MPI_STATUS_SIZE,1:msgCount), &
                       iError)
    endif

    do e1 = 1,this%nElem
      do s1 = 1,4
        e2Global = this%sideInfo(3,s1,e1)
        s2 = this%sideInfo(4,s1,e1)/10
        nloc1(1:2) = nid1(1:2,s1,e1)
        nloc2(1:2) = nid2(1:2,s1,e1)

        if(e2Global > 0) then
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

  endsubroutine RecalculateFlip_Mesh2D_t

  subroutine Write_Mesh2D_t(this,meshFile)
    ! Writes mesh output in HOPR format (serial only)
    implicit none
    class(Mesh2D_t),intent(inout) :: this
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

  endsubroutine Write_Mesh2D_t

endmodule SELF_Mesh_2D_t
