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
  use SELF_Quadrature
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

  ! Length of a material-name string stored in the mesh
  integer,parameter :: SELF_MESH_MATNAME_LENGTH = 64

  type,extends(SEMMesh) :: Mesh2D_t
    integer,pointer,dimension(:,:,:) :: sideInfo
    real(prec),pointer,dimension(:,:,:,:) :: nodeCoords
    integer,pointer,dimension(:,:) :: elemInfo
    integer,pointer,dimension(:,:,:) :: globalNodeIDs
    integer,pointer,dimension(:,:) :: CGNSCornerMap
    integer,pointer,dimension(:,:) :: CGNSSideMap
    integer,pointer,dimension(:,:) :: BCType
    character(LEN=255),allocatable :: BCNames(:)
    ! 2:1 nonconforming (mortar) interface bookkeeping.
    !
    ! A mortar interface joins one "big" element edge to the two half-edges of its 2:1
    ! refined neighbors ("small" sides). Sides that participate in a mortar carry
    ! sideInfo(3) = 0 (no conforming neighbor) and sideInfo(5) = 0 (no boundary condition),
    ! so all conforming-side machinery (SideExchange, RecalculateFlip, boundary-condition
    ! mapping) skips them; sideInfo(1) is set to the mortar index for reference.
    !
    ! mortarInfo(1:8,1:nMortars) is replicated on all ranks (element ids are global):
    !   1 - big element id ; 2 - big local side id
    !   3 - small element id on the sub-edge covering big-edge coordinate [-1,0]
    !   4 - 10*(small local side) + flip for that sub-edge
    !   5 - small element id on the sub-edge covering big-edge coordinate [0,1]
    !   6 - 10*(small local side) + flip for that sub-edge
    !   7 - global side id of sub-edge 1 (used as an MPI message tag)
    !   8 - global side id of sub-edge 2 (used as an MPI message tag)
    ! flip = 0 when the small side's edge coordinate runs in the same direction as the big
    ! side's edge coordinate, flip = 1 when reversed (same convention as sideInfo(4)).
    integer :: nMortars = 0
    integer,pointer,dimension(:,:) :: mortarInfo => null()
    ! Material tracking: every element has an integer material id
    ! indexing into materialNames. Single-material readers (HOPr,
    ! structured, ISM, ISM-v2) leave nMaterials = 1 with the name
    ! "default". The ISM-MM reader populates the table with the
    ! material strings from the .mesh file.
    integer :: nMaterials = 0
    integer,allocatable :: elemMaterial(:)
    character(LEN=SELF_MESH_MATNAME_LENGTH),allocatable :: materialNames(:)

  contains
    procedure,public :: Init => Init_Mesh2D_t
    procedure,public :: Free => Free_Mesh2D_t
    procedure,public :: UpdateDevice => UpdateDevice_Mesh2D_t

    generic,public :: StructuredMesh => UniformStructuredMesh_Mesh2D_t
    procedure,private :: UniformStructuredMesh_Mesh2D_t
    procedure,public :: SimpleMortarMesh => SimpleMortarMesh_Mesh2D_t
    procedure,public :: DoubleMortarMesh => DoubleMortarMesh_Mesh2D_t
    procedure,public :: ResetBoundaryConditionType => ResetBoundaryConditionType_Mesh2D_t

    procedure,public :: Read_HOPr => Read_HOPr_Mesh2D_t
    procedure,public :: Read_HOHQMesh => Read_HOHQMesh_Mesh2D_t

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

    ! Default material table: a single "default" material covers all
    ! elements. Readers that carry material information overwrite
    ! these allocations with the per-file table.
    this%nMaterials = 1
    allocate(this%elemMaterial(1:nElem))
    allocate(this%materialNames(1:1))
    this%elemMaterial = 1
    this%materialNames(1) = "default"

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
    if(allocated(this%elemMaterial)) deallocate(this%elemMaterial)
    if(allocated(this%materialNames)) deallocate(this%materialNames)
    this%nMaterials = 0
    if(associated(this%mortarInfo)) deallocate(this%mortarInfo)
    this%mortarInfo => null()
    this%nMortars = 0
    call this%decomp%Free()

  endsubroutine Free_Mesh2D_t

  subroutine UpdateDevice_Mesh2D_t(this)
    implicit none
    class(Mesh2D_t),intent(inout) :: this
    if(.false.) this%nElem = this%nElem ! CPU stub; suppress unused-dummy-argument warning
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

  subroutine SimpleMortarMesh_Mesh2D_t(this,dx,bcids)
  !!
  !! Create the smallest 2:1 nonconforming (mortar) mesh: one "big" element of size
  !! 2*dx x 2*dx whose east edge is shared with the west edges of two "small" dx x dx
  !! elements. The mesh is conforming everywhere except at the single mortar interface.
  !!
  !!      y = 2*dx  _______________ ______
  !!               |               |  e3  |
  !!               |               |______|
  !!               |      e1       |  e2  |
  !!               |_______________|______|
  !!             x = 0           2*dx    3*dx
  !!
  !!  Input
  !!    - this : Fresh/empty mesh2d_t object
  !!    - dx : Edge length of the small elements; the big element has edge length 2*dx
  !!    - bcids(1:4) : Boundary condition flags for the south, east, north, and west
  !!                   sides of the domain
  !!
  !!  Output
  !!    - this : mesh2d_t object with three elements and one mortar interface
  !!
  !! This mesh is primarily intended for testing and demonstration of the mortar
  !! interface support. Element geometry is bilinear (nGeo=1). With domain
  !! decomposition on two or more ranks, the mortar interface straddles the rank
  !! boundary (elements 1,2 | 3 for two ranks).
  !!
    implicit none
    class(Mesh2D_t),intent(out) :: this
    real(prec),intent(in) :: dx
    integer,intent(in) :: bcids(1:4)
    ! Local
    integer,parameter :: nGlobalElem = 3
    integer,parameter :: nUniqueSides = 10
    real(prec) :: nodeCoords(1:2,1:2,1:2,1:nGlobalElem)
    integer :: globalNodeIDs(1:2,1:2,1:nGlobalElem)
    integer :: sideInfo(1:5,1:4,1:nGlobalElem)
    integer :: e1,e2
    integer :: nLocalElems
    integer :: nGeo,nBCs

    call this%decomp%init()

    nGeo = 1 ! Bilinear element geometry
    nBCs = 4

    ! Element 1 (big) : [0,2dx] x [0,2dx]
    nodeCoords(1:2,1,1,1) = [0.0_prec,0.0_prec]
    nodeCoords(1:2,2,1,1) = [2.0_prec*dx,0.0_prec]
    nodeCoords(1:2,1,2,1) = [0.0_prec,2.0_prec*dx]
    nodeCoords(1:2,2,2,1) = [2.0_prec*dx,2.0_prec*dx]
    globalNodeIDs(1,1,1) = 1
    globalNodeIDs(2,1,1) = 2
    globalNodeIDs(1,2,1) = 6
    globalNodeIDs(2,2,1) = 7

    ! Element 2 (small, lower) : [2dx,3dx] x [0,dx]
    nodeCoords(1:2,1,1,2) = [2.0_prec*dx,0.0_prec]
    nodeCoords(1:2,2,1,2) = [3.0_prec*dx,0.0_prec]
    nodeCoords(1:2,1,2,2) = [2.0_prec*dx,dx]
    nodeCoords(1:2,2,2,2) = [3.0_prec*dx,dx]
    globalNodeIDs(1,1,2) = 2
    globalNodeIDs(2,1,2) = 3
    globalNodeIDs(1,2,2) = 4
    globalNodeIDs(2,2,2) = 5

    ! Element 3 (small, upper) : [2dx,3dx] x [dx,2dx]
    nodeCoords(1:2,1,1,3) = [2.0_prec*dx,dx]
    nodeCoords(1:2,2,1,3) = [3.0_prec*dx,dx]
    nodeCoords(1:2,1,2,3) = [2.0_prec*dx,2.0_prec*dx]
    nodeCoords(1:2,2,2,3) = [3.0_prec*dx,2.0_prec*dx]
    globalNodeIDs(1,1,3) = 4
    globalNodeIDs(2,1,3) = 5
    globalNodeIDs(1,2,3) = 7
    globalNodeIDs(2,2,3) = 8

    ! Side connectivity. Global side ids:
    !  1: e1-S, 2: mortar sub-edge 1 (e1-E lower / e2-W), 3: mortar sub-edge 2
    !  (e1-E upper / e3-W), 4: e1-N, 5: e1-W, 6: e2-S, 7: e2-E, 8: e2-N/e3-S,
    !  9: e3-E, 10: e3-N
    sideInfo = 0

    ! Element 1 (big)
    sideInfo(2,1,1) = 1
    sideInfo(5,1,1) = bcids(1) ! south -> domain south
    sideInfo(1,2,1) = 1 ! east -> mortar 1, big side
    sideInfo(2,2,1) = 2
    sideInfo(2,3,1) = 4
    sideInfo(5,3,1) = bcids(3) ! north -> domain north
    sideInfo(2,4,1) = 5
    sideInfo(5,4,1) = bcids(4) ! west -> domain west

    ! Element 2 (small, lower)
    sideInfo(2,1,2) = 6
    sideInfo(5,1,2) = bcids(1) ! south -> domain south
    sideInfo(2,2,2) = 7
    sideInfo(5,2,2) = bcids(2) ! east -> domain east
    sideInfo(2,3,2) = 8 ! north -> conforming interior side shared with element 3
    sideInfo(3,3,2) = 3
    sideInfo(4,3,2) = 10*1 ! neighbor's south side, flip 0
    sideInfo(1,4,2) = 1 ! west -> mortar 1, small side on sub-edge 1
    sideInfo(2,4,2) = 2

    ! Element 3 (small, upper)
    sideInfo(2,1,3) = 8 ! south -> conforming interior side shared with element 2
    sideInfo(3,1,3) = 2
    sideInfo(4,1,3) = 10*3 ! neighbor's north side, flip 0
    sideInfo(2,2,3) = 9
    sideInfo(5,2,3) = bcids(2) ! east -> domain east
    sideInfo(2,3,3) = 10
    sideInfo(5,3,3) = bcids(3) ! north -> domain north
    sideInfo(1,4,3) = 1 ! west -> mortar 1, small side on sub-edge 2
    sideInfo(2,4,3) = 3

    ! Domain decomposition. The message count upper bound is oversized relative to
    ! nUniqueSides to accommodate the per-variable (and per-direction) mortar and
    ! conforming side messages on this small mesh.
    call this%decomp%GenerateDecomposition(nGlobalElem,64*nUniqueSides)

    e1 = this%decomp%offsetElem(this%decomp%rankId+1)+1
    e2 = this%decomp%offsetElem(this%decomp%rankId+2)
    nLocalElems = e2-e1+1

    call this%Init(nGeo,nLocalElems,nLocalElems*4,nLocalElems*4,nBCs)
    this%nUniqueSides = nUniqueSides
    this%quadrature = UNIFORM
    this%BCType = 0
    this%elemInfo = 0

    this%nodeCoords(1:2,1:2,1:2,1:nLocalElems) = nodeCoords(1:2,1:2,1:2,e1:e2)
    this%globalNodeIDs(1:2,1:2,1:nLocalElems) = globalNodeIDs(1:2,1:2,e1:e2)
    this%sideInfo(1:5,1:4,1:nLocalElems) = sideInfo(1:5,1:4,e1:e2)

    ! The mortar table is replicated on all ranks; element ids are global
    this%nMortars = 1
    allocate(this%mortarInfo(1:8,1:1))
    this%mortarInfo(1:8,1) = [1,2, & ! big element, big local side (east)
                              2,10*4, & ! sub-edge 1 : element 2, west side, flip 0
                              3,10*4, & ! sub-edge 2 : element 3, west side, flip 0
                              2,3] ! global side ids for the two sub-edges

    call this%UpdateDevice()

  endsubroutine SimpleMortarMesh_Mesh2D_t

  subroutine DoubleMortarMesh_Mesh2D_t(this,dx,bcids)
  !!
  !! Create a six-element mesh with two 2:1 mortar interfaces, used to validate the
  !! mortar machinery in configurations the SimpleMortarMesh cannot reach:
  !!
  !!      y = 4*dx  ______________ ______
  !!               |              | e6 R |   R : element 6 is rotated 180 degrees, so
  !!               |     e2       |______|       its side facing e2 has flip = 1, and
  !!               |              |  e5  |       its side facing e5 is a conforming
  !!      y = 2*dx |______________|______|       interior side with flip = 1
  !!               |              |  e4  |
  !!               |     e1       |______|
  !!               |              |  e3  |
  !!               |______________|______|
  !!             x = 0          2*dx   3*dx
  !!
  !!  - Two mortar interfaces (element ordering exercises the mortar-index strides
  !!    in the exchange buffers).
  !!  - Mortar 2's second sub-edge has flip = 1 (element 6's edge coordinate runs
  !!    opposite the big side's), exercising every trace-reorientation branch.
  !!  - With two ranks (elements 1-3 | 4-6), mortar 1 splits big/small across ranks
  !!    while mortar 2 places BOTH small elements remote from the big element's
  !!    rank, so the big-side trace is received independently for each sub-edge.
  !!
  !!  Input
  !!    - this : Fresh/empty mesh2d_t object
  !!    - dx : Edge length of the small elements; big elements have edge length 2*dx
  !!    - bcids(1:4) : Boundary condition flags for the south, east, north, and west
  !!                   sides of the domain
  !!
    implicit none
    class(Mesh2D_t),intent(out) :: this
    real(prec),intent(in) :: dx
    integer,intent(in) :: bcids(1:4)
    ! Local
    integer,parameter :: nGlobalElem = 6
    integer,parameter :: nUniqueSides = 18
    real(prec) :: nodeCoords(1:2,1:2,1:2,1:nGlobalElem)
    integer :: globalNodeIDs(1:2,1:2,1:nGlobalElem)
    integer :: sideInfo(1:5,1:4,1:nGlobalElem)
    integer :: e1,e2
    integer :: nLocalElems
    integer :: nGeo,nBCs

    call this%decomp%init()

    nGeo = 1 ! Bilinear element geometry
    nBCs = 4

    ! Element 1 (big) : [0,2dx] x [0,2dx]
    nodeCoords(1:2,1,1,1) = [0.0_prec,0.0_prec]
    nodeCoords(1:2,2,1,1) = [2.0_prec*dx,0.0_prec]
    nodeCoords(1:2,1,2,1) = [0.0_prec,2.0_prec*dx]
    nodeCoords(1:2,2,2,1) = [2.0_prec*dx,2.0_prec*dx]
    globalNodeIDs(1:2,1,1) = [1,2]
    globalNodeIDs(1:2,2,1) = [6,7]

    ! Element 2 (big) : [0,2dx] x [2dx,4dx]
    nodeCoords(1:2,1,1,2) = [0.0_prec,2.0_prec*dx]
    nodeCoords(1:2,2,1,2) = [2.0_prec*dx,2.0_prec*dx]
    nodeCoords(1:2,1,2,2) = [0.0_prec,4.0_prec*dx]
    nodeCoords(1:2,2,2,2) = [2.0_prec*dx,4.0_prec*dx]
    globalNodeIDs(1:2,1,2) = [6,7]
    globalNodeIDs(1:2,2,2) = [11,12]

    ! Element 3 (small) : [2dx,3dx] x [0,dx]
    nodeCoords(1:2,1,1,3) = [2.0_prec*dx,0.0_prec]
    nodeCoords(1:2,2,1,3) = [3.0_prec*dx,0.0_prec]
    nodeCoords(1:2,1,2,3) = [2.0_prec*dx,dx]
    nodeCoords(1:2,2,2,3) = [3.0_prec*dx,dx]
    globalNodeIDs(1:2,1,3) = [2,3]
    globalNodeIDs(1:2,2,3) = [4,5]

    ! Element 4 (small) : [2dx,3dx] x [dx,2dx]
    nodeCoords(1:2,1,1,4) = [2.0_prec*dx,dx]
    nodeCoords(1:2,2,1,4) = [3.0_prec*dx,dx]
    nodeCoords(1:2,1,2,4) = [2.0_prec*dx,2.0_prec*dx]
    nodeCoords(1:2,2,2,4) = [3.0_prec*dx,2.0_prec*dx]
    globalNodeIDs(1:2,1,4) = [4,5]
    globalNodeIDs(1:2,2,4) = [7,8]

    ! Element 5 (small) : [2dx,3dx] x [2dx,3dx]
    nodeCoords(1:2,1,1,5) = [2.0_prec*dx,2.0_prec*dx]
    nodeCoords(1:2,2,1,5) = [3.0_prec*dx,2.0_prec*dx]
    nodeCoords(1:2,1,2,5) = [2.0_prec*dx,3.0_prec*dx]
    nodeCoords(1:2,2,2,5) = [3.0_prec*dx,3.0_prec*dx]
    globalNodeIDs(1:2,1,5) = [7,8]
    globalNodeIDs(1:2,2,5) = [9,10]

    ! Element 6 (small, rotated 180 degrees) : [2dx,3dx] x [3dx,4dx]
    ! Corner 1 sits at the physical northeast corner, so the element's local
    ! coordinate directions run opposite the physical x and y directions. The
    ! corner ordering remains counterclockwise (positive Jacobian).
    nodeCoords(1:2,1,1,6) = [3.0_prec*dx,4.0_prec*dx]
    nodeCoords(1:2,2,1,6) = [2.0_prec*dx,4.0_prec*dx]
    nodeCoords(1:2,1,2,6) = [3.0_prec*dx,3.0_prec*dx]
    nodeCoords(1:2,2,2,6) = [2.0_prec*dx,3.0_prec*dx]
    globalNodeIDs(1:2,1,6) = [13,12]
    globalNodeIDs(1:2,2,6) = [10,9]

    ! Side connectivity. Global side ids:
    !  1: e1-S, 2: mortar 1 sub-edge 1 (e1-E lower / e3-W), 3: mortar 1 sub-edge 2
    !  (e1-E upper / e4-W), 4: e1-N/e2-S, 5: e1-W, 6: mortar 2 sub-edge 1
    !  (e2-E lower / e5-W), 7: mortar 2 sub-edge 2 (e2-E upper / e6 local side 2,
    !  flip 1), 8: e2-N, 9: e2-W, 10: e3-S, 11: e3-E, 12: e3-N/e4-S, 13: e4-E,
    !  14: e4-N/e5-S, 15: e5-E, 16: e5-N/e6 local side 3 (flip 1), 17: e6 local
    !  side 4 (domain east), 18: e6 local side 1 (domain north)
    sideInfo = 0

    ! Element 1 (big, lower)
    sideInfo(2,1,1) = 1
    sideInfo(5,1,1) = bcids(1) ! south -> domain south
    sideInfo(1,2,1) = 1 ! east -> mortar 1, big side
    sideInfo(2,2,1) = 2
    sideInfo(2,3,1) = 4 ! north -> conforming interior side shared with element 2
    sideInfo(3,3,1) = 2
    sideInfo(4,3,1) = 10*1 ! neighbor's south side, flip 0
    sideInfo(2,4,1) = 5
    sideInfo(5,4,1) = bcids(4) ! west -> domain west

    ! Element 2 (big, upper)
    sideInfo(2,1,2) = 4 ! south -> conforming interior side shared with element 1
    sideInfo(3,1,2) = 1
    sideInfo(4,1,2) = 10*3 ! neighbor's north side, flip 0
    sideInfo(1,2,2) = 2 ! east -> mortar 2, big side
    sideInfo(2,2,2) = 6
    sideInfo(2,3,2) = 8
    sideInfo(5,3,2) = bcids(3) ! north -> domain north
    sideInfo(2,4,2) = 9
    sideInfo(5,4,2) = bcids(4) ! west -> domain west

    ! Element 3 (small, lower-right)
    sideInfo(2,1,3) = 10
    sideInfo(5,1,3) = bcids(1) ! south -> domain south
    sideInfo(2,2,3) = 11
    sideInfo(5,2,3) = bcids(2) ! east -> domain east
    sideInfo(2,3,3) = 12 ! north -> conforming interior side shared with element 4
    sideInfo(3,3,3) = 4
    sideInfo(4,3,3) = 10*1 ! neighbor's south side, flip 0
    sideInfo(1,4,3) = 1 ! west -> mortar 1, small side on sub-edge 1
    sideInfo(2,4,3) = 2

    ! Element 4 (small)
    sideInfo(2,1,4) = 12 ! south -> conforming interior side shared with element 3
    sideInfo(3,1,4) = 3
    sideInfo(4,1,4) = 10*3 ! neighbor's north side, flip 0
    sideInfo(2,2,4) = 13
    sideInfo(5,2,4) = bcids(2) ! east -> domain east
    sideInfo(2,3,4) = 14 ! north -> conforming interior side shared with element 5
    sideInfo(3,3,4) = 5
    sideInfo(4,3,4) = 10*1 ! neighbor's south side, flip 0
    sideInfo(1,4,4) = 1 ! west -> mortar 1, small side on sub-edge 2
    sideInfo(2,4,4) = 3

    ! Element 5 (small)
    sideInfo(2,1,5) = 14 ! south -> conforming interior side shared with element 4
    sideInfo(3,1,5) = 4
    sideInfo(4,1,5) = 10*3 ! neighbor's north side, flip 0
    sideInfo(2,2,5) = 15
    sideInfo(5,2,5) = bcids(2) ! east -> domain east
    sideInfo(2,3,5) = 16 ! north -> conforming side shared with element 6; element
    sideInfo(3,3,5) = 6 ! 6's edge coordinate runs opposite this element's
    sideInfo(4,3,5) = 10*3+1 ! neighbor's local north side, flip 1
    sideInfo(1,4,5) = 2 ! west -> mortar 2, small side on sub-edge 1
    sideInfo(2,4,5) = 6

    ! Element 6 (small, rotated). Local sides map to physical directions as:
    ! side 1 -> domain north, side 2 -> mortar 2 (faces element 2), side 3 ->
    ! conforming side shared with element 5, side 4 -> domain east.
    sideInfo(2,1,6) = 18
    sideInfo(5,1,6) = bcids(3) ! local south -> domain north
    sideInfo(1,2,6) = 2 ! local east -> mortar 2, small side on sub-edge 2
    sideInfo(2,2,6) = 7
    sideInfo(2,3,6) = 16 ! local north -> conforming side shared with element 5
    sideInfo(3,3,6) = 5
    sideInfo(4,3,6) = 10*3+1 ! neighbor's north side, flip 1
    sideInfo(2,4,6) = 17
    sideInfo(5,4,6) = bcids(2) ! local west -> domain east

    ! Domain decomposition; oversized message bound as in SimpleMortarMesh
    call this%decomp%GenerateDecomposition(nGlobalElem,64*nUniqueSides)

    e1 = this%decomp%offsetElem(this%decomp%rankId+1)+1
    e2 = this%decomp%offsetElem(this%decomp%rankId+2)
    nLocalElems = e2-e1+1

    call this%Init(nGeo,nLocalElems,nLocalElems*4,nLocalElems*4,nBCs)
    this%nUniqueSides = nUniqueSides
    this%quadrature = UNIFORM
    this%BCType = 0
    this%elemInfo = 0

    this%nodeCoords(1:2,1:2,1:2,1:nLocalElems) = nodeCoords(1:2,1:2,1:2,e1:e2)
    this%globalNodeIDs(1:2,1:2,1:nLocalElems) = globalNodeIDs(1:2,1:2,e1:e2)
    this%sideInfo(1:5,1:4,1:nLocalElems) = sideInfo(1:5,1:4,e1:e2)

    ! The mortar table is replicated on all ranks; element ids are global
    this%nMortars = 2
    allocate(this%mortarInfo(1:8,1:2))
    this%mortarInfo(1:8,1) = [1,2, & ! big element 1, east side
                              3,10*4, & ! sub-edge 1 : element 3, west side, flip 0
                              4,10*4, & ! sub-edge 2 : element 4, west side, flip 0
                              2,3] ! global side ids for the two sub-edges
    this%mortarInfo(1:8,2) = [2,2, & ! big element 2, east side
                              5,10*4, & ! sub-edge 1 : element 5, west side, flip 0
                              6,10*2+1, & ! sub-edge 2 : element 6, local east side, flip 1
                              6,7] ! global side ids for the two sub-edges

    call this%UpdateDevice()

  endsubroutine DoubleMortarMesh_Mesh2D_t

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

  subroutine Read_HOHQMesh_Mesh2D_t(this,meshFile)
    !! Reader for HOHQMesh text mesh files in the ISM and ISM-MM
    !! formats. The format is auto-detected from the first line:
    !!   * Line equal to "ISM-MM" (trimmed) => ISM-MM with per-element
    !!     material name strings and a 4-int count line that includes
    !!     an unused nEdges field (the ISM-MM writer in HOHQMesh does
    !!     NOT emit an edge block).
    !!   * Anything else is treated as plain ISM: the first line is
    !!     itself the count line "nNodes nElems polyOrder" and there
    !!     are no per-element material names.
    !!
    !! See HOHQMesh `WriteISMMeshFile` in `Source/IO/MeshOutputMethods.f90`
    !! for the canonical writer. Curve sample points are written at
    !! Chebyshev-Gauss-Lobatto nodes, so the resulting mesh has
    !! `quadrature = CHEBYSHEV_GAUSS_LOBATTO`.
    !!
    !! Element-side connectivity (sideInfo(3:4,...)) is not present in
    !! either of these formats, so neighbors are reconstructed here by
    !! matching pairs of corner-node IDs across elements. Boundary
    !! names from the .mesh file are collected into this%BCNames, and
    !! sideInfo(5,...) carries the corresponding 1-based index (or 0
    !! for interior faces). Material names (ISM-MM) populate
    !! this%materialNames and this%elemMaterial.
    implicit none
    class(Mesh2D_t),intent(out) :: this
    character(*),intent(in) :: meshFile
    ! Local
    integer :: iUnit
    integer :: ios
    integer :: nNodesFile,nEdgesFile,nElemFile,polyOrder
    integer :: nGeo,nLocal
    integer :: i,j,k,e,iSide
    integer :: cornerIDs(1:4)
    integer :: bCurveFlag(1:4)
    integer :: cN1,cN2
    integer :: ePair,sPair
    integer :: bcIdx
    integer :: matIdx
    integer :: nBCsLocal
    integer :: nMatsLocal
    integer :: nCornerPairs
    integer :: hashKey,bucket,probe
    integer :: hashSize
    integer,allocatable :: hashHead(:)
    integer,allocatable :: hashNext(:)
    integer,allocatable :: pairN1(:),pairN2(:)
    integer,allocatable :: pairElem(:),pairSide(:)
    integer,allocatable :: ismCorners(:,:) ! 4 x nElem
    integer,allocatable :: ismBCid(:,:) ! 4 x nElem
    integer,allocatable :: ismFlag(:,:) ! 4 x nElem
    integer,allocatable :: ismMat(:) ! nElem
    real(prec),allocatable :: nodeXY(:,:) ! 2 x nNodesFile
    real(prec),allocatable :: edgeCurve(:,:,:,:) ! 2, nGeo+1, 4, nElem
    real(prec) :: xyz(1:3)
    character(LEN=255) :: header
    character(LEN=SELF_MESH_MATNAME_LENGTH) :: matName
    character(LEN=255) :: bdyNames(1:4)
    character(LEN=255),allocatable :: BCNamesLocal(:)
    character(LEN=SELF_MESH_MATNAME_LENGTH),allocatable :: matNamesLocal(:)
    logical :: isISM_MM

    call this%decomp%init()

    open(newunit=iUnit,file=trim(meshFile),status='old',action='read', &
         form='formatted',iostat=ios)
    if(ios /= 0) then
      print*,__FILE__//' : Failed to open '//trim(meshFile)
      stop 1
    endif

    print*,__FILE__//' : Reading HOHQMesh mesh from '//trim(meshFile)

    ! ---- 1. Detect format from the first line ----
    read(iUnit,'(A)') header
    if(trim(adjustl(header)) == "ISM-MM") then
      isISM_MM = .true.
      read(iUnit,*) nNodesFile,nEdgesFile,nElemFile,polyOrder
    else
      isISM_MM = .false.
      ! The header line we just read IS the count line for plain ISM.
      read(header,*) nNodesFile,nElemFile,polyOrder
      nEdgesFile = 0
    endif

    nGeo = polyOrder
    print*,__FILE__//' : Format = ',merge("ISM-MM","ISM   ",isISM_MM)
    print*,__FILE__//' : nNodes = ',nNodesFile,' nElem = ',nElemFile, &
      ' polyOrder = ',polyOrder

    ! ---- 2. Read all node coordinates (drop the z component for 2D) ----
    allocate(nodeXY(1:2,1:nNodesFile))
    do i = 1,nNodesFile
      read(iUnit,*) xyz(1:3)
      nodeXY(1:2,i) = xyz(1:2)
    enddo

    ! Neither ISM nor ISM-MM writes an edge block, even though ISM-MM's
    ! count line declares nEdges. The HOHQMesh writer only emits edges
    ! for the ISM-V2 variant.

    ! ---- 3. Per-element block ----
    allocate(ismCorners(1:4,1:nElemFile))
    allocate(ismFlag(1:4,1:nElemFile))
    allocate(ismBCid(1:4,1:nElemFile))
    allocate(ismMat(1:nElemFile))
    allocate(edgeCurve(1:2,1:nGeo+1,1:4,1:nElemFile))
    edgeCurve = 0.0_prec

    ! Boundary-name table built incrementally
    nBCsLocal = 0
    allocate(BCNamesLocal(1:16))
    BCNamesLocal = ""

    ! Material-name table built incrementally
    nMatsLocal = 0
    allocate(matNamesLocal(1:8))
    matNamesLocal = ""

    do e = 1,nElemFile

      if(isISM_MM) then
        read(iUnit,*) cornerIDs(1:4),matName
        ! lookup/insert material name
        matIdx = 0
        do k = 1,nMatsLocal
          if(trim(matNamesLocal(k)) == trim(matName)) then
            matIdx = k; exit
          endif
        enddo
        if(matIdx == 0) then
          nMatsLocal = nMatsLocal+1
          if(nMatsLocal > size(matNamesLocal)) call grow_string_table(matNamesLocal)
          matNamesLocal(nMatsLocal) = matName
          matIdx = nMatsLocal
        endif
        ismMat(e) = matIdx
      else
        read(iUnit,*) cornerIDs(1:4)
        ismMat(e) = 1
      endif
      ismCorners(1:4,e) = cornerIDs

      read(iUnit,*) bCurveFlag(1:4)
      ismFlag(1:4,e) = bCurveFlag

      do k = 1,4
        if(bCurveFlag(k) == 1) then
          do j = 1,nGeo+1
            read(iUnit,*) xyz(1:3)
            edgeCurve(1:2,j,k,e) = xyz(1:2)
          enddo
        endif
      enddo

      read(iUnit,*) bdyNames(1:4)
      do k = 1,4
        if(trim(adjustl(bdyNames(k))) == "---") then
          ismBCid(k,e) = 0
        else
          ! lookup/insert bdy name
          bcIdx = 0
          do i = 1,nBCsLocal
            if(trim(BCNamesLocal(i)) == trim(adjustl(bdyNames(k)))) then
              bcIdx = i; exit
            endif
          enddo
          if(bcIdx == 0) then
            nBCsLocal = nBCsLocal+1
            if(nBCsLocal > size(BCNamesLocal)) call grow_bc_table(BCNamesLocal)
            BCNamesLocal(nBCsLocal) = trim(adjustl(bdyNames(k)))
            bcIdx = nBCsLocal
          endif
          ismBCid(k,e) = bcIdx
        endif
      enddo
    enddo

    close(iUnit)

    ! At least one BC slot must exist so that Init allocates BCNames/BCType
    nBCsLocal = max(nBCsLocal,1)

    ! Set up the domain decomposition arrays (elemToRank, offsetElem,
    ! request/stat slots) using the file's element count. This is
    ! required for SideExchange even in the serial single-rank case.
    call this%decomp%GenerateDecomposition(nElemFile,4*nElemFile)

    ! ---- 4. Allocate SELF Mesh2D_t and populate ----
    call this%Init(nGeo,nElemFile,4*nElemFile,nElemFile*(nGeo+1)**2,nBCsLocal)
    this%nUniqueSides = 0 ! filled below after pair matching
    this%quadrature = CHEBYSHEV_GAUSS_LOBATTO ! HOHQMesh curve samples are at CGL points
    this%BCType = 0
    do i = 1,nBCsLocal
      if(BCNamesLocal(i) /= "") then
        this%BCNames(i) = BCNamesLocal(i)
      else
        this%BCNames(i) = "unused"
      endif
    enddo

    ! Replace the default single-material table with the parsed one
    deallocate(this%materialNames)
    this%nMaterials = max(nMatsLocal,1)
    allocate(this%materialNames(1:this%nMaterials))
    if(nMatsLocal == 0) then
      this%materialNames(1) = "default"
    else
      this%materialNames(1:nMatsLocal) = matNamesLocal(1:nMatsLocal)
    endif
    this%elemMaterial = ismMat

    ! Place corner nodes and run transfinite interpolation per element
    do e = 1,nElemFile
      call build_nodeCoords_for_element(this,e,nGeo,cornerIDs=ismCorners(:,e), &
                                        flag=ismFlag(:,e),edgeCurve=edgeCurve(:,:,:,e), &
                                        nodeXY=nodeXY)
      ! Synthesize globalNodeIDs: stamp the four corners with their file
      ! IDs and leave interior IDs as 0 (interior nodes are private to
      ! the element under our high-order tensor product layout).
      this%globalNodeIDs(:,:,e) = 0
      this%globalNodeIDs(1,1,e) = ismCorners(1,e)
      this%globalNodeIDs(nGeo+1,1,e) = ismCorners(2,e)
      this%globalNodeIDs(nGeo+1,nGeo+1,e) = ismCorners(3,e)
      this%globalNodeIDs(1,nGeo+1,e) = ismCorners(4,e)

      ! Pack elemInfo with simple placeholders; SELF's 2D path does not
      ! depend on the HOPR-style offset fields when the mesh comes from
      ! a non-HOPR reader.
      this%elemInfo(1,e) = 0
      this%elemInfo(2,e) = ismMat(e) ! Zone = material id
      this%elemInfo(3,e) = 4*(e-1)
      this%elemInfo(4,e) = 4*e
      this%elemInfo(5,e) = (nGeo+1)**2*(e-1)
      this%elemInfo(6,e) = (nGeo+1)**2*e
    enddo

    ! ---- 5. Build sideInfo via corner-pair matching ----
    nCornerPairs = 4*nElemFile
    allocate(pairN1(1:nCornerPairs),pairN2(1:nCornerPairs))
    allocate(pairElem(1:nCornerPairs),pairSide(1:nCornerPairs))
    do e = 1,nElemFile
      do iSide = 1,4
        call corner_pair_for_side(ismCorners(:,e),iSide,cN1,cN2)
        i = iSide+4*(e-1)
        pairN1(i) = min(cN1,cN2)
        pairN2(i) = max(cN1,cN2)
        pairElem(i) = e
        pairSide(i) = iSide
      enddo
    enddo

    ! Hash chain by first node id
    hashSize = max(nNodesFile,nCornerPairs)+1
    allocate(hashHead(0:hashSize-1),hashNext(1:nCornerPairs))
    hashHead = 0
    hashNext = 0
    do i = 1,nCornerPairs
      hashKey = pairN1(i)
      bucket = modulo(hashKey,hashSize)
      hashNext(i) = hashHead(bucket)
      hashHead(bucket) = i
    enddo

    this%sideInfo = 0
    this%nUniqueSides = 0
    do e = 1,nElemFile
      do iSide = 1,4
        i = iSide+4*(e-1)
        ePair = 0; sPair = 0
        bucket = modulo(pairN1(i),hashSize)
        probe = hashHead(bucket)
        do while(probe /= 0)
          if(probe /= i .and. pairN1(probe) == pairN1(i) .and. &
             pairN2(probe) == pairN2(i)) then
            ePair = pairElem(probe)
            sPair = pairSide(probe)
            exit
          endif
          probe = hashNext(probe)
        enddo
        this%sideInfo(3,iSide,e) = ePair
        this%sideInfo(4,iSide,e) = 10*sPair ! flip resolved by RecalculateFlip
        this%sideInfo(5,iSide,e) = ismBCid(iSide,e)
        ! Allocate a globalSideID (count each shared side once)
        if(ePair == 0 .or. e < ePair) then
          this%nUniqueSides = this%nUniqueSides+1
          this%sideInfo(2,iSide,e) = this%nUniqueSides
        else
          this%sideInfo(2,iSide,e) = this%sideInfo(2,sPair,ePair)
        endif
      enddo
    enddo

    deallocate(hashHead,hashNext,pairN1,pairN2,pairElem,pairSide)
    deallocate(nodeXY,ismCorners,ismFlag,ismBCid,ismMat,edgeCurve)
    deallocate(BCNamesLocal,matNamesLocal)

    call this%RecalculateFlip()
    call this%UpdateDevice()

  contains

    subroutine grow_bc_table(tbl)
      character(LEN=255),allocatable,intent(inout) :: tbl(:)
      character(LEN=255),allocatable :: tmp(:)
      integer :: oldSize
      oldSize = size(tbl)
      allocate(tmp(1:2*oldSize))
      tmp(1:oldSize) = tbl(1:oldSize)
      tmp(oldSize+1:) = ""
      call move_alloc(tmp,tbl)
    endsubroutine grow_bc_table

    subroutine grow_string_table(tbl)
      character(LEN=SELF_MESH_MATNAME_LENGTH),allocatable,intent(inout) :: tbl(:)
      character(LEN=SELF_MESH_MATNAME_LENGTH),allocatable :: tmp(:)
      integer :: oldSize
      oldSize = size(tbl)
      allocate(tmp(1:2*oldSize))
      tmp(1:oldSize) = tbl(1:oldSize)
      tmp(oldSize+1:) = ""
      call move_alloc(tmp,tbl)
    endsubroutine grow_string_table

  endsubroutine Read_HOHQMesh_Mesh2D_t

  subroutine corner_pair_for_side(corners,iSide,cN1,cN2)
    !! Returns the two corner-node IDs delimiting a given local side
    !! using SELF's 2D side convention:
    !!   Side 1 South = [CN1, CN2]
    !!   Side 2 East  = [CN2, CN3]
    !!   Side 3 North = [CN4, CN3]
    !!   Side 4 West  = [CN1, CN4]
    implicit none
    integer,intent(in) :: corners(1:4)
    integer,intent(in) :: iSide
    integer,intent(out) :: cN1,cN2
    select case(iSide)
    case(1); cN1 = corners(1); cN2 = corners(2)
    case(2); cN1 = corners(2); cN2 = corners(3)
    case(3); cN1 = corners(4); cN2 = corners(3)
    case(4); cN1 = corners(1); cN2 = corners(4)
    case default; cN1 = 0; cN2 = 0
    endselect
  endsubroutine corner_pair_for_side

  subroutine build_nodeCoords_for_element(mesh,e,nGeo,cornerIDs,flag,edgeCurve,nodeXY)
    !! Fill mesh%nodeCoords(:,:,:,e) for one element. Place the four
    !! corner nodes from the file's node table, then perform linear
    !! (nGeo=1) placement or transfinite interpolation (nGeo>1) using
    !! the (possibly curved) edge curves. For straight sides the edge
    !! curve is filled in here by linear interpolation between corners
    !! at Chebyshev-Gauss-Lobatto parametric coordinates.
    implicit none
    class(Mesh2D_t),intent(inout) :: mesh
    integer,intent(in) :: e
    integer,intent(in) :: nGeo
    integer,intent(in) :: cornerIDs(1:4)
    integer,intent(in) :: flag(1:4)
    real(prec),intent(in) :: edgeCurve(1:2,1:nGeo+1,1:4)
    real(prec),intent(in) :: nodeXY(:,:)
    ! Local
    real(prec) :: P(1:2,1:4)
    real(prec) :: gS(1:2,1:nGeo+1),gE(1:2,1:nGeo+1)
    real(prec) :: gN(1:2,1:nGeo+1),gW(1:2,1:nGeo+1)
    real(prec) :: xi(0:nGeo),wts(0:nGeo)
    real(prec) :: u(1:nGeo+1),t
    integer :: i,j

    P(1:2,1) = nodeXY(1:2,cornerIDs(1))
    P(1:2,2) = nodeXY(1:2,cornerIDs(2))
    P(1:2,3) = nodeXY(1:2,cornerIDs(3))
    P(1:2,4) = nodeXY(1:2,cornerIDs(4))

    if(nGeo == 1) then
      ! Pure bilinear element with the 4 corners
      mesh%nodeCoords(1:2,1,1,e) = P(1:2,1)
      mesh%nodeCoords(1:2,nGeo+1,1,e) = P(1:2,2)
      mesh%nodeCoords(1:2,nGeo+1,nGeo+1,e) = P(1:2,3)
      mesh%nodeCoords(1:2,1,nGeo+1,e) = P(1:2,4)
      return
    endif

    ! Chebyshev-Gauss-Lobatto parametric coordinates on [-1,1] and
    ! their [0,1] image used in the TFI blending weights. We use
    ! SELF's quadrature routine to keep node placement consistent
    ! with the rest of the code.
    call ChebyshevQuadrature(nGeo,xi,wts,CHEBYSHEV_GAUSS_LOBATTO)
    do i = 1,nGeo+1
      u(i) = 0.5_prec*(xi(i-1)+1.0_prec)
    enddo

    ! Edge curves. HOHQMesh's `edgeMap(2,4)` is `(1,2 ; 2,3 ; 4,3 ; 1,4)`,
    ! so its sides 1/2/3/4 already run in SELF's (+ξ, +η, +ξ, +η)
    ! directions. No reversal is needed.
    do i = 1,nGeo+1
      if(flag(1) == 1) then
        gS(1:2,i) = edgeCurve(1:2,i,1)
      else
        t = u(i)
        gS(1:2,i) = (1.0_prec-t)*P(1:2,1)+t*P(1:2,2)
      endif

      if(flag(2) == 1) then
        gE(1:2,i) = edgeCurve(1:2,i,2)
      else
        t = u(i)
        gE(1:2,i) = (1.0_prec-t)*P(1:2,2)+t*P(1:2,3)
      endif

      if(flag(3) == 1) then
        gN(1:2,i) = edgeCurve(1:2,i,3)
      else
        t = u(i)
        gN(1:2,i) = (1.0_prec-t)*P(1:2,4)+t*P(1:2,3)
      endif

      if(flag(4) == 1) then
        gW(1:2,i) = edgeCurve(1:2,i,4)
      else
        t = u(i)
        gW(1:2,i) = (1.0_prec-t)*P(1:2,1)+t*P(1:2,4)
      endif
    enddo

    ! Transfinite (Coons-patch) interpolation. u,v in [0,1].
    do j = 1,nGeo+1
      do i = 1,nGeo+1
        mesh%nodeCoords(1:2,i,j,e) = &
          (1.0_prec-u(j))*gS(1:2,i)+u(j)*gN(1:2,i)+ &
          (1.0_prec-u(i))*gW(1:2,j)+u(i)*gE(1:2,j)- &
          ((1.0_prec-u(i))*(1.0_prec-u(j))*P(1:2,1)+ &
           u(i)*(1.0_prec-u(j))*P(1:2,2)+ &
           u(i)*u(j)*P(1:2,3)+ &
           (1.0_prec-u(i))*u(j)*P(1:2,4))
      enddo
    enddo
  endsubroutine build_nodeCoords_for_element

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
    integer :: i,j
    integer :: l
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
