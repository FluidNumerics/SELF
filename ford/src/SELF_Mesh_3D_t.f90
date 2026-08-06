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

module SELF_Mesh_3D_t

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
! xi3 direction points from "Bottom" (xi3=-1) to "Top" (xi3=1)
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

  ! Side Ordering
  integer,parameter :: selfSide3D_Bottom = 1
  integer,parameter :: selfSide3D_South = 2
  integer,parameter :: selfSide3D_East = 3
  integer,parameter :: selfSide3D_North = 4
  integer,parameter :: selfSide3D_West = 5
  integer,parameter :: selfSide3D_Top = 6

  type,extends(SEMMesh) :: Mesh3D_t
    integer,pointer,dimension(:,:,:) :: sideInfo
    real(prec),pointer,dimension(:,:,:,:,:) :: nodeCoords
    integer,pointer,dimension(:,:) :: elemInfo
    integer,pointer,dimension(:,:,:,:) :: globalNodeIDs
    integer,pointer,dimension(:,:) :: CGNSCornerMap
    integer,pointer,dimension(:,:) :: sideMap
    integer,pointer,dimension(:,:) :: CGNSSideMap
    integer,pointer,dimension(:,:) :: BCType
    character(LEN=255),allocatable :: BCNames(:)
    ! Material tracking: every element has an integer material id
    ! indexing into materialNames. Single-material constructors and
    ! readers (HOPr, structured, periodic) leave nMaterials = 1 with
    ! the name "default". The HOHQMesh ISM-MM reader populates the
    ! table with the material strings from the .mesh file.
    integer :: nMaterials = 0
    integer,allocatable :: elemMaterial(:)
    character(LEN=SELF_MESH_MATNAME_LENGTH),allocatable :: materialNames(:)

  contains

    procedure,public :: Init => Init_Mesh3D_t
    procedure,public :: Free => Free_Mesh3D_t
    procedure,public :: UpdateDevice => UpdateDevice_Mesh3D_t

    generic,public :: StructuredMesh => UniformStructuredMesh_Mesh3D_t
    procedure,private :: UniformStructuredMesh_Mesh3D_t

    generic,public :: PeriodicStructuredMesh => UniformPeriodicMesh_Mesh3D_t
    procedure,private :: UniformPeriodicMesh_Mesh3D_t

    procedure,public :: Read_HOPr => Read_HOPr_Mesh3D_t
    procedure,public :: Read_HOHQMesh => Read_HOHQMesh_Mesh3D_t

    procedure,public :: ResetBoundaryConditionType => ResetBoundaryConditionType_Mesh3D_t

    procedure,public :: Write_Mesh => Write_Mesh3D_t

    procedure,public :: RecalculateFlip => RecalculateFlip_Mesh3D_t

  endtype Mesh3D_t

  integer,private :: CGNStoSELFflip(1:6,1:6,1:4)

  ! This table maps the primary side, secondary side, and CGNS flip values
  ! to indexing flips that are used in SELF.
  ! This table is used after reading in HOPr mesh information in "RecalculateFlip"
  ! SELF's flip indices correspond to the following scenarios
  !
  ! 0    i2 = i1     j2 = j1
  ! 1    i2 = N-i1   j2 = j1
  ! 2    i2 = N-i1   j2 = N-j1
  ! 3    i2 = i1     j2 = N-j1
  ! 4    i2 = j1     j2 = i1
  ! 5    i2 = N-j1   j2 = i1
  ! 6    i2 = N-j1   j2 = N-i1
  ! 7    i2 = j1     j2 = N-i1
  !
  data CGNStoSELFflip/ &
    4,0,0,1,4,0, &
    0,4,4,5,0,4, &
    0,4,4,5,0,4, &
    1,7,7,6,1,7, &
    4,0,0,1,4,0, &
    0,4,4,5,0,4, &
    3,5,5,4,3,5, &
    7,1,1,0,7,1, &
    7,1,1,0,7,1, &
    4,0,0,1,4,0, &
    3,5,5,4,3,5, &
    7,1,1,0,7,1, &
    6,2,2,3,6,2, &
    2,6,6,7,2,6, &
    2,6,6,7,2,6, &
    3,5,5,4,3,5, &
    6,2,2,3,6,2, &
    2,6,6,7,2,6, &
    1,7,7,6,1,7, &
    5,3,3,2,5,3, &
    5,3,3,2,5,3, &
    6,2,2,3,6,2, &
    1,7,7,6,1,7, &
    5,3,3,2,5,3/

contains

  subroutine Init_Mesh3D_t(this,nGeo,nElem,nSides,nNodes,nBCs)
    implicit none
    class(Mesh3D_t),intent(inout) :: this
    integer,intent(in) :: nGeo
    integer,intent(in) :: nElem
    integer,intent(in) :: nSides
    integer,intent(in) :: nNodes
    integer,intent(in) :: nBCs
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

    ! Sidemap traverses each face so that the normal
    ! formed by the right hand rule is the coordinate
    ! positive pointing normal. For east,north,and top
    ! this is an outward facing normal.
    ! For bottom, south, and west, the normal is inward
    ! facing.
    this%sideMap(1:4,1) = (/1,2,3,4/) ! Bottom
    this%sideMap(1:4,2) = (/1,2,6,5/) ! South
    this%sideMap(1:4,3) = (/2,3,7,6/) ! East
    this%sideMap(1:4,4) = (/4,3,7,8/) ! North
    this%sideMap(1:4,5) = (/1,4,8,5/) ! West
    this%sideMap(1:4,6) = (/5,6,7,8/) ! Top

  endsubroutine Init_Mesh3D_t

  subroutine Free_Mesh3D_t(this)
    implicit none
    class(Mesh3D_t),intent(inout) :: this

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

    deallocate(this%BCNames)
    if(allocated(this%elemMaterial)) deallocate(this%elemMaterial)
    if(allocated(this%materialNames)) deallocate(this%materialNames)
    this%nMaterials = 0
    call this%decomp%Free()

  endsubroutine Free_Mesh3D_t

  subroutine UpdateDevice_Mesh3D_t(this)
    implicit none
    class(Mesh3D_t),intent(inout) :: this
    if(.false.) this%nElem = this%nElem ! CPU stub; suppress unused-dummy-argument warning
  endsubroutine UpdateDevice_Mesh3D_t

  subroutine ResetBoundaryConditionType_Mesh3D_t(this,bcid)
    !! This method can be used to reset all of the boundary elements
    !! boundary condition type to the desired value.
    !!
    !! Note that ALL physical boundaries will be set to have this boundary
    !! condition
    implicit none
    class(Mesh3D_t),intent(inout) :: this
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

    call this%UpdateDevice()

  endsubroutine ResetBoundaryConditionType_Mesh3D_t

  subroutine RecalculateFlip_Mesh3D_t(this)
    implicit none
    class(Mesh3D_t),intent(inout) :: this
    ! Local
    integer :: e1
    integer :: s1
    integer :: e2
    integer :: s2
    integer :: cgnsFlip,selfFlip

    do e1 = 1,this%nElem
      do s1 = 1,6

        e2 = this%sideInfo(3,s1,e1)
        s2 = this%sideInfo(4,s1,e1)/10
        cgnsFlip = this%sideInfo(4,s1,e1)-s2*10

        if(e2 /= 0) then

          selfFlip = CGNStoSELFflip(s2,s1,cgnsFlip)
          this%sideInfo(4,s1,e1) = 10*s2+selfFlip

        endif

      enddo
    enddo

  endsubroutine RecalculateFlip_Mesh3D_t

  pure function elementid(i,j,k,ti,tj,tk,nxpertile,nypertile,nzpertile, &
                          ntilex,ntiley,ntilez) result(eid)
    integer,intent(in) :: i,j,k
    integer,intent(in) :: ti,tj,tk
    integer,intent(in) :: nxpertile,nypertile,nzpertile
    integer,intent(in) :: ntilex,ntiley,ntilez
    integer :: eid

    eid = i+nxpertile*(j-1+nypertile*(k-1+nzpertile*( &
                                      ti-1+ntilex*(tj-1+ntiley*(tk-1)))))
    if(.false.) eid = eid+ntilez ! suppress unused-dummy-argument warning

  endfunction elementid

  pure function gid2eid(gx,gy,gz,nxpertile,nypertile,nzpertile, &
                        ntilex,ntiley,ntilez) result(eid)
  !! Map a global element position (gx,gy,gz), with gx in [1,nX] etc., to the
  !! tile-ordered element id used throughout the structured mesh. This is the
  !! inverse of the (i,ti) decomposition: gx = i + nxpertile*(ti-1).
    integer,intent(in) :: gx,gy,gz
    integer,intent(in) :: nxpertile,nypertile,nzpertile
    integer,intent(in) :: ntilex,ntiley,ntilez
    integer :: eid
    integer :: i,j,k,ti,tj,tk

    ti = (gx-1)/nxpertile+1; i = gx-(ti-1)*nxpertile
    tj = (gy-1)/nypertile+1; j = gy-(tj-1)*nypertile
    tk = (gz-1)/nzpertile+1; k = gz-(tk-1)*nzpertile

    eid = elementid(i,j,k,ti,tj,tk,nxpertile,nypertile,nzpertile, &
                    ntilex,ntiley,ntilez)

  endfunction gid2eid

  subroutine UniformStructuredMesh_Mesh3D_t(this,nxPerTile,nyPerTile,nzPerTile, &
                                            nTileX,nTileY,nTileZ,dx,dy,dz,bcids,comm)
  !!
  !! Create a structured mesh and store it in SELF's unstructured mesh format.
  !! The mesh is created in tiles of size (tnx,tny,tnz). Tiling is used to determine
  !! the element ordering.
  !!
  !!
  !!  Input
  !!    - this : Fresh/empty mesh2d_t object
  !!    - nxPerTile : The number of elements in the x direction within a tile
  !!    - nyPerTile : The number of elements in the y direction within a tile
  !!    - nzPerTile : The number of elements in the z direction within a tile
  !!    - nTileX : The number of tiles in the x direction
  !!    - nTileY : The number of tiles in the y direction
  !!    - nTileZ : The number of tiles in the z direction
  !!    - dx : Element width in the x-direction
  !!    - dy : Element width in the y-direction
  !!    - dz : Element width in the z-direction
  !!    - bcids(1:6) : Boundary condition flags for the south, east, north, and west sides of the domain
  !!    - comm (optional) : Externally managed MPI communicator (Fortran integer handle). When
  !!      provided, the caller owns MPI initialization/finalization.
  !!
  !!  Output
  !!    - this : mesh2d_t object with vertices, faces, and element information
  !!
  !! Total number of elements in the x-direction is nX = nxPerTile*nTileX
  !! Total number of elements in the y-direction is nY = nyPerTile*nTileY
  !!
  !! Length of the domain in the x-direction is Lx = dx*nX
  !! Length of the domain in the y-direction is Ly = dy*nY
  !!
    implicit none
    class(Mesh3D_t),intent(out) :: this
    integer,intent(in) :: nxPerTile
    integer,intent(in) :: nyPerTile
    integer,intent(in) :: nzPerTile
    integer,intent(in) :: nTileX
    integer,intent(in) :: nTileY
    integer,intent(in) :: nTileZ
    real(prec),intent(in) :: dx
    real(prec),intent(in) :: dy
    real(prec),intent(in) :: dz
    integer,intent(in) :: bcids(1:6)
    integer,intent(in),optional :: comm
    ! Local
    integer :: nX,nY,nZ,nGeo,nBCs
    integer :: nGlobalElem
    integer :: nUniqueSides
    integer :: nUniqueNodes
    integer :: nLocalElems
    integer :: nLocalSides
    integer :: nLocalNodes
    real(prec),allocatable :: nodeCoords(:,:,:,:,:)
    integer,allocatable :: globalNodeIDs(:,:,:,:)
    integer,allocatable :: sideInfo(:,:,:)
    integer :: i,j,k,ti,tj,tk
    integer :: ix,iy,iz,iel
    integer :: ni,nj,nk
    integer :: e1,e2,s1,s2
    integer :: nfaces

    call this%decomp%init(comm)

    nX = nTileX*nxPerTile
    nY = nTileY*nyPerTile
    nZ = nTileZ*nzPerTile
    nGeo = 1 ! Force the geometry to be linear
    nBCs = 6 ! Force the number of boundary conditions to 4

    nGlobalElem = nX*nY*nZ
    nUniqueSides = (nX+1)*nY*nZ+(nY+1)*nX*nZ+(nZ+1)*nX*nY
    nUniqueNodes = (nX+1)*(nY+1)*(nZ+1)

    allocate(nodeCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nGlobalElem))
    allocate(globalNodeIDs(1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nGlobalElem))
    allocate(sideInfo(1:5,1:6,1:nGlobalElem))

    do tk = 1,nTileZ
      do tj = 1,nTileY
        do ti = 1,nTileX
          do k = 1,nzPerTile
            iz = k+nzPerTile*(tk-1)
            do j = 1,nyPerTile
              iy = j+nyPerTile*(tj-1)
              do i = 1,nxPerTile

                iel = elementid(i,j,k,ti,tj,tk, &
                                nxpertile,nypertile,nzpertile, &
                                ntilex,ntiley,ntilez)
                ix = i+nxPerTile*(ti-1)

                do nk = 1,nGeo+1
                  do nj = 1,nGeo+1
                    do ni = 1,nGeo+1
                      nodeCoords(1,ni,nj,nk,iel) = real(ni-1+ix-1,prec)*dx
                      nodeCoords(2,ni,nj,nk,iel) = real(nj-1+iy-1,prec)*dy
                      nodeCoords(3,ni,nj,nk,iel) = real(nk-1+iz-1,prec)*dz
                      globalNodeIDs(ni,nj,nk,iel) = ni-1+i+(nxPerTile+1)*( &
                                                    nj-1+j-1+(nyPerTile+1)*( &
                                                    nk-1+k-1+(nzPerTile+1)*( &
                                                    (ti-1+nTileX*( &
                                                     tj-1+nTileY*(tk-1))))))
                    enddo
                  enddo
                enddo

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    ! Fill in face information
    !  sideInfo(1:5,iSide,iEl)
    !    1 - Side Type (currently unused in SELF)
    !    2 - Global Side ID (Used for message passing)
    !    3 - Neighbor Element ID
    !    4 - 10*( neighbor local side )  + flip
    !    5 - Boundary Condition ID
    nfaces = 0
    do tk = 1,nTileZ
      do tj = 1,nTileY
        do ti = 1,nTileX
          do k = 1,nzPerTile
            do j = 1,nyPerTile
              do i = 1,nxPerTile

                iel = elementid(i,j,k,ti,tj,tk, &
                                nxpertile,nypertile,nzpertile, &
                                ntilex,ntiley,ntilez)
                ! bottom, iside=1
                s1 = 1
                s2 = 6
                if(k == 1) then ! bottom most part of the tile
                  if(tk == 1) then ! bottom most tile
                    nfaces = nfaces+1
                    sideinfo(2,s1,iel) = nfaces
                    sideinfo(3,s1,iel) = 0 ! Neigbor element (null, boundary condition)
                    sideinfo(4,s1,iel) = 0 ! Neighbor side id (null, boundary condition)
                    sideinfo(5,s1,iel) = bcids(s1) ! Boundary condition id; set from the user input
                  else ! interior tile
                    !neighbor element is the top most element in the tile beneath
                    e2 = elementid(i,j,nzpertile,ti,tj,tk-1, &
                                   nxpertile,nypertile,nzpertile, &
                                   ntilex,ntiley,ntilez)

                    sideinfo(2,s1,iel) = sideInfo(2,s2,e2) ! Copy the face id from neighbor
                    sideinfo(3,s1,iel) = e2
                    sideinfo(4,s1,iel) = 10*s2 ! Neighbor side id
                    sideinfo(5,s1,iel) = 0 ! Boundary condition id; (null, interior face)
                  endif
                else ! interior to the tile
                  !neighbor element is in the same tile, but beneath
                  e2 = elementid(i,j,k-1,ti,tj,tk, &
                                 nxpertile,nypertile,nzpertile, &
                                 ntilex,ntiley,ntilez)

                  sideinfo(2,s1,iel) = sideInfo(2,s2,e2) ! Copy the face id from neighbor
                  sideinfo(3,s1,iel) = e2
                  sideinfo(4,s1,iel) = 10*s2 ! Neighbor side id
                  sideinfo(5,s1,iel) = 0 ! Boundary condition id; (null, interior face)
                endif

                ! south, iside=2
                s1 = 2
                s2 = 4 ! Neighbor side is north (4)
                if(j == 1) then ! southern  most part of the tile
                  if(tj == 1) then ! southern most tile
                    nfaces = nfaces+1
                    sideinfo(2,s1,iel) = nfaces
                    sideinfo(3,s1,iel) = 0 ! Neigbor element (null, boundary condition)
                    sideinfo(4,s1,iel) = 0 ! Neighbor side id (null, boundary condition)
                    sideinfo(5,s1,iel) = bcids(s1) ! Boundary condition id; eastern boundary set from the user input
                  else ! interior tile
                    !neighbor element is northernmost element in the tile to the south
                    e2 = elementid(i,nypertile,k,ti,tj-1,tk, &
                                   nxpertile,nypertile,nzpertile, &
                                   ntilex,ntiley,ntilez)

                    sideinfo(2,s1,iel) = sideInfo(2,s2,e2) ! Copy the face id from neighbor
                    sideinfo(3,s1,iel) = e2 ! Neigbor element
                    sideinfo(4,s1,iel) = 10*s2 ! Neighbor side id
                    sideinfo(5,s1,iel) = 0 ! Boundary condition id; (null, interior face)
                  endif
                else ! interior to the tile
                  !neighbor element is in the same tile, to the south
                  e2 = elementid(i,j-1,k,ti,tj,tk, &
                                 nxpertile,nypertile,nzpertile, &
                                 ntilex,ntiley,ntilez)

                  sideinfo(2,s1,iel) = sideInfo(2,s2,e2) ! Copy the face id from neighbor
                  sideinfo(3,s1,iel) = e2 ! Neigbor element
                  sideinfo(4,s1,iel) = 10*s2 ! Neighbor side id
                  sideinfo(5,s1,iel) = 0 ! Boundary condition id; (null, interior face)
                endif

                ! east, iside=3
                s1 = 3
                s2 = 5 ! neighbor side id is west (5)
                ! East faces are always new faces, due to the way we are traversing the grid
                nfaces = nfaces+1
                sideinfo(2,s1,iel) = nfaces
                if(i == nxPerTile) then ! eastern most part of the tile
                  if(ti == nTileX) then ! eastern most tile
                    sideinfo(3,s1,iel) = 0 ! Neigbor element (null, boundary condition)
                    sideinfo(4,s1,iel) = 0 ! Neighbor side id (null, boundary condition)
                    sideinfo(5,s1,iel) = bcids(s1) ! Boundary condition id;
                  else ! interior tile
                    !neighbor element is westernmost element in tile to the east
                    e2 = elementid(1,j,k,ti+1,tj,tk, &
                                   nxpertile,nypertile,nzpertile, &
                                   ntilex,ntiley,ntilez)
                    sideinfo(3,s1,iel) = e2 ! Neigbor element
                    sideinfo(4,s1,iel) = 10*s2 ! Neighbor side id
                    sideinfo(5,s1,iel) = 0 ! Boundary condition id; (null, interior face)
                  endif
                else ! interior to the tile
                  !neighbor element is in the same tile, to the east
                  e2 = elementid(i+1,j,k,ti,tj,tk, &
                                 nxpertile,nypertile,nzpertile, &
                                 ntilex,ntiley,ntilez)
                  sideinfo(3,s1,iel) = e2 ! Neigbor element
                  sideinfo(4,s1,iel) = 10*s2 ! Neighbor side id
                  sideinfo(5,s1,iel) = 0 ! Boundary condition id; (null, interior face)
                endif

                ! north, iside=4
                s1 = 4
                s2 = 2 ! neighbor side is south (2)
                ! North faces are always new faces, due to the way we are traversing the grid
                nfaces = nfaces+1
                sideinfo(2,s1,iel) = nfaces
                if(j == nyPerTile) then ! northern most part of the tile
                  if(tj == nTileY) then ! northern most tile
                    sideinfo(3,s1,iel) = 0 ! Neigbor element (null, boundary condition)
                    sideinfo(4,s1,iel) = 0 ! Neighbor side id (null, boundary condition)
                    sideinfo(5,s1,iel) = bcids(s1) ! Boundary condition id; set from the user input
                  else ! interior tile, but northern most face of the tile
                    !neighbor element is the southernmost element in the tile to the north
                    e2 = elementid(i,1,k,ti,tj+1,tk, &
                                   nxpertile,nypertile,nzpertile, &
                                   ntilex,ntiley,ntilez)
                    sideinfo(3,s1,iel) = e2 ! Neigbor element
                    sideinfo(4,s1,iel) = 10*s2 ! Neighbor side id
                    sideinfo(5,s1,iel) = 0 ! Boundary condition id; (null, interior face)
                  endif
                else ! interior to the tile
                  !neighbor element is the tile to the north
                  e2 = elementid(i,j+1,k,ti,tj,tk, &
                                 nxpertile,nypertile,nzpertile, &
                                 ntilex,ntiley,ntilez)
                  sideinfo(3,s1,iel) = e2 ! Neigbor element
                  sideinfo(4,s1,iel) = 10*s2 ! Neighbor side id
                  sideinfo(5,s1,iel) = 0 ! Boundary condition id; (null, interior face)
                endif

                ! west, iside=5
                s1 = 5
                s2 = 3 ! neighbor side id is east (3)
                if(i == 1) then ! western most part of the tile
                  if(ti == 1) then ! western most tile
                    nfaces = nfaces+1
                    sideinfo(2,s1,iel) = nfaces
                    sideinfo(3,s1,iel) = 0 ! Neigbor element (null, boundary condition)
                    sideinfo(4,s1,iel) = 0 ! Neighbor side id (null, boundary condition)
                    sideinfo(5,s1,iel) = bcids(s1) ! Boundary condition id
                  else ! interior tile, but western most face of the tile
                    !neighbor element is the easternmost element in the tile to the west
                    e2 = elementid(nxperTile,j,k,ti-1,tj,tk, &
                                   nxpertile,nypertile,nzpertile, &
                                   ntilex,ntiley,ntilez)

                    sideinfo(2,s1,iel) = sideInfo(2,s2,e2) ! Copy the face id from neighbor's east face
                    sideinfo(3,s1,iel) = e2
                    sideinfo(4,s1,iel) = 10*s2 ! Neighbor side id - neighbor to the west, east side (2)
                    sideinfo(5,s1,iel) = 0 ! Boundary condition id; (null, interior face)
                  endif
                else ! interior to the tile
                  !neighbor element is the element to the west in the same tile
                  e2 = elementid(i-1,j,k,ti,tj,tk, &
                                 nxpertile,nypertile,nzpertile, &
                                 ntilex,ntiley,ntilez)

                  sideinfo(2,s1,iel) = sideInfo(2,s2,e2) ! Copy the face id from neighbor's east face
                  sideinfo(3,s1,iel) = e2
                  sideinfo(4,s1,iel) = 10*s2 ! Neighbor side id - neighbor to the west, east side (2)
                  sideinfo(5,s1,iel) = 0 ! Boundary condition id; (null, interior face)
                endif

                ! top, iside=6
                s1 = 6
                s2 = 1 ! neighbor side is bottom (1)
                ! Top faces are always new faces, due to the way we are traversing the grid
                nfaces = nfaces+1
                sideinfo(2,s1,iel) = nfaces
                if(k == nzPerTile) then ! top most part of the tile
                  if(tk == nTileZ) then ! top most tile
                    sideinfo(3,s1,iel) = 0 ! Neigbor element (null, boundary condition)
                    sideinfo(4,s1,iel) = 0 ! Neighbor side id (null, boundary condition)
                    sideinfo(5,s1,iel) = bcids(s1) ! Boundary condition id; set from the user input
                  else ! interior tile, but top most face of the tile
                    !neighbor element is the bottom-most element in the tile above
                    e2 = elementid(i,j,1,ti,tj,tk+1, &
                                   nxpertile,nypertile,nzpertile, &
                                   ntilex,ntiley,ntilez)
                    sideinfo(3,s1,iel) = e2 ! Neigbor element
                    sideinfo(4,s1,iel) = 10*s2 ! Neighbor side id
                    sideinfo(5,s1,iel) = 0 ! Boundary condition id; (null, interior face)
                  endif
                else ! interior to the tile
                  !neighbor element is the tile above
                  e2 = elementid(i,j,k+1,ti,tj,tk, &
                                 nxpertile,nypertile,nzpertile, &
                                 ntilex,ntiley,ntilez)
                  sideinfo(3,s1,iel) = e2 ! Neigbor element, inside same tile, to the north
                  sideinfo(4,s1,iel) = 10*s2 ! Neighbor side id - neighbor to the north, south side (1)
                  sideinfo(5,s1,iel) = 0 ! Boundary condition id; (null, interior face)
                endif

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    call this%decomp%GenerateDecomposition(nGlobalElem,nUniqueSides)

    e1 = this%decomp%offsetElem(this%decomp%rankId+1)+1
    e2 = this%decomp%offsetElem(this%decomp%rankId+2)
    nLocalElems = e2-e1+1

    nLocalSides = nLocalElems*6
    nLocalNodes = nLocalElems*8
    call this%Init(nGeo,nLocalElems,nLocalSides,nLocalNodes,nBCs)
    this%nUniqueSides = nUniqueSides
    this%quadrature = UNIFORM

    this%nodeCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nLocalElems) = nodeCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,e1:e2)
    this%globalNodeIDs(1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nLocalElems) = globalNodeIDs(1:nGeo+1,1:nGeo+1,1:nGeo+1,e1:e2)
    this%sideInfo(1:5,1:6,1:nLocalElems) = sideInfo(1:5,1:6,e1:e2)

    deallocate(nodeCoords)
    deallocate(globalNodeIDs)
    deallocate(sideInfo)

    call this%UpdateDevice()

  endsubroutine UniformStructuredMesh_Mesh3D_t

  subroutine UniformPeriodicMesh_Mesh3D_t(this,nxPerTile,nyPerTile,nzPerTile, &
                                          nTileX,nTileY,nTileZ,dx,dy,dz,comm)
  !!
  !! Create a fully triply-periodic structured hexahedral mesh and store it in
  !! SELF's unstructured mesh format. Element geometry and ordering are identical
  !! to UniformStructuredMesh; the only difference is connectivity. The faces on
  !! the six domain boundaries are wired as interior faces whose neighbor is the
  !! element on the opposite side of the domain. This realises the triply
  !! periodic box T^3 = [0,Lx] x [0,Ly] x [0,Lz] required by, e.g., the
  !! Arnold-Beltrami-Childress (ABC) flow benchmark.
  !!
  !!  Input
  !!    - this : Fresh/empty mesh3d_t object
  !!    - nxPerTile,nyPerTile,nzPerTile : Elements per tile in each direction
  !!    - nTileX,nTileY,nTileZ : Number of tiles in each direction
  !!    - dx,dy,dz : Element widths in each direction
  !!
  !! Total elements:  nX = nxPerTile*nTileX, etc.
  !! Domain lengths:  Lx = dx*nX, etc.
  !!
  !! Connectivity notes
  !! ------------------
  !! Local side numbering is 1=bottom, 2=south, 3=east, 4=north, 5=west, 6=top.
  !! Because the mesh is uniform and axis-aligned, each periodic face pair has
  !! the same orientation as the corresponding interior face pair generated by
  !! UniformStructuredMesh, so the side "flip" index is 0 throughout (the
  !! identity face mapping in SideExchange). This matches the convention used by
  !! UniformStructuredMesh, which assigns flip 0 to every interior face and never
  !! calls RecalculateFlip.
  !!
  !! Global side IDs are assigned analytically so that the two faces of every
  !! periodic pair share a single ID. With nXYZ = nX*nY*nZ:
  !!   x-faces (east/west)   : 1          .. nXYZ
  !!   y-faces (south/north) : nXYZ+1     .. 2*nXYZ
  !!   z-faces (bottom/top)  : 2*nXYZ+1   .. 3*nXYZ
  !! The east face of the element at global position (gx,gy,gz) and the west face
  !! of its +x neighbor both reference xface(gx,gy,gz); the periodic wrap
  !! identifies gx=nX (east of the last element) with the west face of the first.
  !! There are therefore exactly 3*nElem unique sides and no domain boundaries.
  !!
    implicit none
    class(Mesh3D_t),intent(out) :: this
    integer,intent(in) :: nxPerTile
    integer,intent(in) :: nyPerTile
    integer,intent(in) :: nzPerTile
    integer,intent(in) :: nTileX
    integer,intent(in) :: nTileY
    integer,intent(in) :: nTileZ
    real(prec),intent(in) :: dx
    real(prec),intent(in) :: dy
    real(prec),intent(in) :: dz
    integer,intent(in),optional :: comm
    ! Local
    integer :: nX,nY,nZ,nGeo,nBCs
    integer :: nGlobalElem
    integer :: nUniqueSides
    integer :: nUniqueNodes
    integer :: nLocalElems
    integer :: nLocalSides
    integer :: nLocalNodes
    integer :: nXYZ
    real(prec),allocatable :: nodeCoords(:,:,:,:,:)
    integer,allocatable :: globalNodeIDs(:,:,:,:)
    integer,allocatable :: sideInfo(:,:,:)
    integer :: i,j,k,ti,tj,tk
    integer :: ix,iy,iz,iel
    integer :: ni,nj,nk
    integer :: e1,e2
    integer :: gx,gy,gz,gxm,gxp,gym,gyp,gzm,gzp

    call this%decomp%init(comm)

    nX = nTileX*nxPerTile
    nY = nTileY*nyPerTile
    nZ = nTileZ*nzPerTile
    nGeo = 1 ! Force the geometry to be linear
    nBCs = 1 ! No domain boundaries; allocate a single (unused) slot
    nXYZ = nX*nY*nZ

    nGlobalElem = nX*nY*nZ
    nUniqueSides = 3*nGlobalElem ! periodic: nXYZ x-faces + nXYZ y-faces + nXYZ z-faces
    nUniqueNodes = (nX+1)*(nY+1)*(nZ+1)

    allocate(nodeCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nGlobalElem))
    allocate(globalNodeIDs(1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nGlobalElem))
    allocate(sideInfo(1:5,1:6,1:nGlobalElem))

    ! Node coordinates and global node IDs (identical to UniformStructuredMesh)
    do tk = 1,nTileZ
      do tj = 1,nTileY
        do ti = 1,nTileX
          do k = 1,nzPerTile
            iz = k+nzPerTile*(tk-1)
            do j = 1,nyPerTile
              iy = j+nyPerTile*(tj-1)
              do i = 1,nxPerTile

                iel = elementid(i,j,k,ti,tj,tk, &
                                nxpertile,nypertile,nzpertile, &
                                ntilex,ntiley,ntilez)
                ix = i+nxPerTile*(ti-1)

                do nk = 1,nGeo+1
                  do nj = 1,nGeo+1
                    do ni = 1,nGeo+1
                      nodeCoords(1,ni,nj,nk,iel) = real(ni-1+ix-1,prec)*dx
                      nodeCoords(2,ni,nj,nk,iel) = real(nj-1+iy-1,prec)*dy
                      nodeCoords(3,ni,nj,nk,iel) = real(nk-1+iz-1,prec)*dz
                      globalNodeIDs(ni,nj,nk,iel) = ni-1+i+(nxPerTile+1)*( &
                                                    nj-1+j-1+(nyPerTile+1)*( &
                                                    nk-1+k-1+(nzPerTile+1)*( &
                                                    (ti-1+nTileX*( &
                                                     tj-1+nTileY*(tk-1))))))
                    enddo
                  enddo
                enddo

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    ! Fill in face information with full triple periodicity.
    !  sideInfo(1:5,iSide,iEl)
    !    1 - Side Type (currently unused in SELF)
    !    2 - Global Side ID (shared by both faces of a periodic pair)
    !    3 - Neighbor Element ID (always > 0; periodic mesh has no boundaries)
    !    4 - 10*( neighbor local side ) + flip (flip is 0 for axis-aligned wrap)
    !    5 - Boundary Condition ID (always 0; interior face)
    do gz = 1,nZ
      do gy = 1,nY
        do gx = 1,nX

          iel = gid2eid(gx,gy,gz,nxPerTile,nyPerTile,nzPerTile, &
                        nTileX,nTileY,nTileZ)

          ! +/- neighbor positions with periodic wrap (1..n)
          gxp = modulo(gx,nX)+1; gxm = modulo(gx-2,nX)+1
          gyp = modulo(gy,nY)+1; gym = modulo(gy-2,nY)+1
          gzp = modulo(gz,nZ)+1; gzm = modulo(gz-2,nZ)+1

          ! bottom (s1=1) <-> top (s2=6) of the -z neighbor
          sideInfo(1,1,iel) = 0
          sideInfo(2,1,iel) = 2*nXYZ+gzm+nZ*((gx-1)+nX*(gy-1)) ! zface(gx,gy,gzm)
          sideInfo(3,1,iel) = gid2eid(gx,gy,gzm,nxPerTile,nyPerTile,nzPerTile, &
                                      nTileX,nTileY,nTileZ)
          sideInfo(4,1,iel) = 10*6
          sideInfo(5,1,iel) = 0

          ! south (s1=2) <-> north (s2=4) of the -y neighbor
          sideInfo(1,2,iel) = 0
          sideInfo(2,2,iel) = nXYZ+gym+nY*((gx-1)+nX*(gz-1)) ! yface(gx,gym,gz)
          sideInfo(3,2,iel) = gid2eid(gx,gym,gz,nxPerTile,nyPerTile,nzPerTile, &
                                      nTileX,nTileY,nTileZ)
          sideInfo(4,2,iel) = 10*4
          sideInfo(5,2,iel) = 0

          ! east (s1=3) <-> west (s2=5) of the +x neighbor
          sideInfo(1,3,iel) = 0
          sideInfo(2,3,iel) = gx+nX*((gy-1)+nY*(gz-1)) ! xface(gx,gy,gz)
          sideInfo(3,3,iel) = gid2eid(gxp,gy,gz,nxPerTile,nyPerTile,nzPerTile, &
                                      nTileX,nTileY,nTileZ)
          sideInfo(4,3,iel) = 10*5
          sideInfo(5,3,iel) = 0

          ! north (s1=4) <-> south (s2=2) of the +y neighbor
          sideInfo(1,4,iel) = 0
          sideInfo(2,4,iel) = nXYZ+gy+nY*((gx-1)+nX*(gz-1)) ! yface(gx,gy,gz)
          sideInfo(3,4,iel) = gid2eid(gx,gyp,gz,nxPerTile,nyPerTile,nzPerTile, &
                                      nTileX,nTileY,nTileZ)
          sideInfo(4,4,iel) = 10*2
          sideInfo(5,4,iel) = 0

          ! west (s1=5) <-> east (s2=3) of the -x neighbor
          sideInfo(1,5,iel) = 0
          sideInfo(2,5,iel) = gxm+nX*((gy-1)+nY*(gz-1)) ! xface(gxm,gy,gz)
          sideInfo(3,5,iel) = gid2eid(gxm,gy,gz,nxPerTile,nyPerTile,nzPerTile, &
                                      nTileX,nTileY,nTileZ)
          sideInfo(4,5,iel) = 10*3
          sideInfo(5,5,iel) = 0

          ! top (s1=6) <-> bottom (s2=1) of the +z neighbor
          sideInfo(1,6,iel) = 0
          sideInfo(2,6,iel) = 2*nXYZ+gz+nZ*((gx-1)+nX*(gy-1)) ! zface(gx,gy,gz)
          sideInfo(3,6,iel) = gid2eid(gx,gy,gzp,nxPerTile,nyPerTile,nzPerTile, &
                                      nTileX,nTileY,nTileZ)
          sideInfo(4,6,iel) = 10*1
          sideInfo(5,6,iel) = 0

        enddo
      enddo
    enddo

    call this%decomp%GenerateDecomposition(nGlobalElem,nUniqueSides)

    e1 = this%decomp%offsetElem(this%decomp%rankId+1)+1
    e2 = this%decomp%offsetElem(this%decomp%rankId+2)
    nLocalElems = e2-e1+1

    nLocalSides = nLocalElems*6
    nLocalNodes = nLocalElems*8
    call this%Init(nGeo,nLocalElems,nLocalSides,nLocalNodes,nBCs)
    this%nUniqueSides = nUniqueSides
    this%quadrature = UNIFORM

    this%nodeCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nLocalElems) = nodeCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,e1:e2)
    this%globalNodeIDs(1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nLocalElems) = globalNodeIDs(1:nGeo+1,1:nGeo+1,1:nGeo+1,e1:e2)
    this%sideInfo(1:5,1:6,1:nLocalElems) = sideInfo(1:5,1:6,e1:e2)

    deallocate(nodeCoords)
    deallocate(globalNodeIDs)
    deallocate(sideInfo)

    call this%UpdateDevice()

  endsubroutine UniformPeriodicMesh_Mesh3D_t

  subroutine Read_HOPr_Mesh3D_t(this,meshFile,comm)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    implicit none
    class(Mesh3D_t),intent(out) :: this
    character(*),intent(in) :: meshFile
    integer,intent(in),optional :: comm
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

    call this%decomp%init(comm)

    if(this%decomp%mpiEnabled) then
      call Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId,this%decomp%mpiComm)
    else
      call Open_HDF5(meshFile,H5F_ACC_RDONLY_F,fileId)
    endif

    call ReadAttribute_HDF5(fileId,'nElems',nGlobalElem)
    call ReadAttribute_HDF5(fileId,'Ngeo',nGeo)
    call ReadAttribute_HDF5(fileId,'nBCs',nBCs)
    call ReadAttribute_HDF5(fileId,'nUniqueSides',nUniqueSides)

    ! Read BCType
    allocate(bcType(1:4,1:nBCs))
    if(this%decomp%mpiEnabled) then
      offset(:) = 0
      call ReadArray_HDF5(fileId,'BCType',bcType,offset)
    else
      call ReadArray_HDF5(fileId,'BCType',bcType)
    endif

    ! Read local subarray of ElemInfo
    call this%decomp%GenerateDecomposition(nGlobalElem,nUniqueSides)

    firstElem = this%decomp%offsetElem(this%decomp%rankId+1)+1
    nLocalElems = this%decomp%offsetElem(this%decomp%rankId+2)- &
                  this%decomp%offsetElem(this%decomp%rankId+1)

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
    nLocalNodes = hopr_elemInfo(6,nLocalElems)-hopr_elemInfo(5,1)

    ! Allocate Space for hopr_nodeCoords and hopr_globalNodeIDs !
    allocate(hopr_nodeCoords(1:3,1:nLocalNodes),hopr_globalNodeIDs(1:nLocalNodes))

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
    nLocalSides = hopr_elemInfo(4,nLocalElems)-hopr_elemInfo(3,1)

    ! Allocate space for hopr_sideInfo
    allocate(hopr_sideInfo(1:5,1:nLocalSides))

    if(this%decomp%mpiEnabled) then
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

    call this%UpdateDevice()

  endsubroutine Read_HOPr_Mesh3D_t

  subroutine Read_HOHQMesh_Mesh3D_t(this,meshFile,comm)
    !! Reader for HOHQMesh 3-D (hexahedral) text mesh files in the ISM
    !! and ISM-MM formats, as written by HOHQMesh's `WriteISMHexMeshFile`
    !! (`Source/3DSource/Mesh3DOutputMethods.f90`). Unlike the 2-D
    !! writer, the 3-D writer emits NO format header line: the first
    !! line is always the count line "nNodes nElems polyOrder". The two
    !! variants differ only in the per-element corner-node line:
    !!   * ISM    : 8 corner-node ids
    !!   * ISM-MM : 8 corner-node ids followed by a material-name string
    !! The variant is auto-detected from the presence of the 9th token.
    !!
    !! Each element block contains, in order: the corner-node line, a
    !! line of 6 boundary-face flags, a (polyOrder+1)^2 block of face
    !! points (x,y,z per line, inner index first) for every flagged
    !! face, and a line of 6 boundary-condition names ("---" marks an
    !! interior face). Face points are sampled at Chebyshev-Gauss-
    !! Lobatto points, so the resulting mesh has
    !! `quadrature = CHEBYSHEV_GAUSS_LOBATTO`.
    !!
    !! HOHQMesh numbers hex corners with nodes 1-4 on the bottom face
    !! (counter-clockwise) and 5-8 above them, which matches SELF's
    !! CGNS corner convention exactly. HOHQMesh face numbering
    !! (1=south, 2=north, 3=bottom, 4=east, 5=top, 6=west; see
    !! `FaceFromVolume`) is remapped to SELF's side ordering
    !! (1=bottom, 2=south, 3=east, 4=north, 5=west, 6=top). Face-point
    !! grids are written with the two on-face volume axes in natural
    !! order, which coincides with SELF's boundary index convention,
    !! so no reorientation of the face data is required.
    !!
    !! Element interior nodes are reconstructed by transfinite (Coons)
    !! interpolation of the six face grids; unflagged faces are
    !! bilinear patches of their corner nodes. Element-face
    !! connectivity is not present in the format, so neighbors are
    !! reconstructed by matching the sorted corner-node ids of each
    !! face across elements, and the side "flip" is computed by
    !! matching the corner-node orderings of the paired faces.
    !! Boundary names populate this%BCNames and sideInfo(5,...) carries
    !! the 1-based index into that table (0 for interior faces).
    !! Material names (ISM-MM) populate this%materialNames and
    !! this%elemMaterial; plain ISM meshes keep the single "default"
    !! material.
    implicit none
    class(Mesh3D_t),intent(out) :: this
    character(*),intent(in) :: meshFile
    integer,intent(in),optional :: comm
    ! Local
    integer :: iUnit
    integer :: ios
    integer :: nNodesFile,nElemFile,polyOrder
    integer :: nGeo,ng1
    integer :: i,j,k,e,f,s,iSide,l,p
    integer :: cornerIDs(1:8)
    integer :: probe8(1:8)
    integer :: hohqFlag(1:6)
    integer :: quadA(1:4),quadB(1:4)
    integer :: sortedA(1:4)
    integer :: ePair,sPair,flip
    integer :: bcIdx
    integer :: matIdx
    integer :: nBCsLocal
    integer :: nMatsLocal
    integer :: nFaceRecords
    integer :: hashKey,bucket,probe
    integer :: hashSize
    integer,allocatable :: hashHead(:)
    integer,allocatable :: hashNext(:)
    integer,allocatable :: pairKey(:,:) ! sorted corner ids, 4 x 6*nElem
    integer,allocatable :: pairElem(:),pairSide(:)
    integer,allocatable :: ismCorners(:,:) ! 8 x nElem
    integer,allocatable :: ismBCid(:,:) ! 6 x nElem (SELF side ordering)
    integer,allocatable :: ismFlag(:,:) ! 6 x nElem (SELF side ordering)
    integer,allocatable :: ismMat(:) ! nElem
    real(prec),allocatable :: nodeXYZ(:,:) ! 3 x nNodesFile
    real(prec),allocatable :: faceCurve(:,:,:,:,:) ! 3, ng1, ng1, 6, nElem
    real(prec) :: xyz(1:3)
    character(LEN=512) :: lineBuf
    character(LEN=SELF_MESH_MATNAME_LENGTH) :: matName
    character(LEN=255) :: bdyNames(1:6)
    character(LEN=255),allocatable :: BCNamesLocal(:)
    character(LEN=SELF_MESH_MATNAME_LENGTH),allocatable :: matNamesLocal(:)
    logical :: isISM_MM
    ! Map from HOHQMesh hex face id (1=south, 2=north, 3=bottom,
    ! 4=east, 5=top, 6=west) to SELF local side id
    integer,parameter :: hohqToSelfSide(1:6) = [2,4,1,3,6,5]
    ! Corner permutations realised by SELF's eight face flips: for a
    ! side pair (s1,s2) with flip p, corner l of side s1 (in face-index
    ! order, see sideMap) coincides with corner flipPerm(l,p) of side
    ! s2. Flips 0-3 are rotations/reflections keeping the face axes,
    ! flips 4-7 transpose them (see ApplyFlip/SideExchange kernels).
    integer,parameter :: flipPerm(1:4,0:7) = reshape( &
                         [1,2,3,4, &
                          2,1,4,3, &
                          3,4,1,2, &
                          4,3,2,1, &
                          1,4,3,2, &
                          2,3,4,1, &
                          3,2,1,4, &
                          4,1,2,3],[4,8])

    call this%decomp%init(comm)

    open(newunit=iUnit,file=trim(meshFile),status='old',action='read', &
         form='formatted',iostat=ios)
    if(ios /= 0) then
      print*,__FILE__//' : Failed to open '//trim(meshFile)
      stop 1
    endif

    print*,__FILE__//' : Reading HOHQMesh mesh from '//trim(meshFile)

    ! ---- 1. Count line (3-D ISM files carry no format header) ----
    read(iUnit,*) nNodesFile,nElemFile,polyOrder

    nGeo = polyOrder
    ng1 = nGeo+1
    print*,__FILE__//' : nNodes = ',nNodesFile,' nElem = ',nElemFile, &
      ' polyOrder = ',polyOrder

    ! ---- 2. Read all node coordinates ----
    allocate(nodeXYZ(1:3,1:nNodesFile))
    do i = 1,nNodesFile
      read(iUnit,*) xyz(1:3)
      nodeXYZ(1:3,i) = xyz(1:3)
    enddo

    ! ---- 3. Per-element block ----
    allocate(ismCorners(1:8,1:nElemFile))
    allocate(ismFlag(1:6,1:nElemFile))
    allocate(ismBCid(1:6,1:nElemFile))
    allocate(ismMat(1:nElemFile))
    allocate(faceCurve(1:3,1:ng1,1:ng1,1:6,1:nElemFile))
    faceCurve = 0.0_prec

    ! Boundary-name table built incrementally
    nBCsLocal = 0
    allocate(BCNamesLocal(1:16))
    BCNamesLocal = ""

    ! Material-name table built incrementally
    nMatsLocal = 0
    allocate(matNamesLocal(1:8))
    matNamesLocal = ""

    isISM_MM = .false.
    nFaceRecords = 0

    do e = 1,nElemFile

      ! Corner-node line; the trailing material-name token (ISM-MM)
      ! is detected by probing for a 9th list item after the 8
      ! integer tokens.
      read(iUnit,'(A)') lineBuf
      read(lineBuf,*) cornerIDs(1:8)
      read(lineBuf,*,iostat=ios) probe8,matName
      if(ios /= 0) matName = ""
      if(e == 1) then
        isISM_MM = (matName /= "")
        print*,__FILE__//' : Format = ',merge("ISM-MM","ISM   ",isISM_MM)
      endif

      if(matName /= "") then
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
        ismMat(e) = 1
      endif
      ismCorners(1:8,e) = cornerIDs

      ! Face flags and face-point grids, in HOHQMesh face order
      read(iUnit,*) hohqFlag(1:6)
      do f = 1,6
        s = hohqToSelfSide(f)
        ismFlag(s,e) = hohqFlag(f)
        if(hohqFlag(f) == 1) then
          nFaceRecords = nFaceRecords+1
          do j = 1,ng1
            do i = 1,ng1
              read(iUnit,*) xyz(1:3)
              faceCurve(1:3,i,j,s,e) = xyz(1:3)
            enddo
          enddo
        endif
      enddo

      ! Boundary-condition names, in HOHQMesh face order
      read(iUnit,*) bdyNames(1:6)
      do f = 1,6
        s = hohqToSelfSide(f)
        if(trim(adjustl(bdyNames(f))) == "---") then
          ismBCid(s,e) = 0
        else
          ! lookup/insert bdy name
          bcIdx = 0
          do i = 1,nBCsLocal
            if(trim(BCNamesLocal(i)) == trim(adjustl(bdyNames(f)))) then
              bcIdx = i; exit
            endif
          enddo
          if(bcIdx == 0) then
            nBCsLocal = nBCsLocal+1
            if(nBCsLocal > size(BCNamesLocal)) call grow_bc_table(BCNamesLocal)
            BCNamesLocal(nBCsLocal) = trim(adjustl(bdyNames(f)))
            bcIdx = nBCsLocal
          endif
          ismBCid(s,e) = bcIdx
        endif
      enddo
    enddo

    close(iUnit)

    print*,__FILE__//' : n curved faces = ',nFaceRecords, &
      ' n materials = ',nMatsLocal,' n boundary names = ',nBCsLocal

    ! At least one BC slot must exist so that Init allocates BCNames/BCType
    nBCsLocal = max(nBCsLocal,1)

    ! Set up the domain decomposition arrays (elemToRank, offsetElem,
    ! request/stat slots) using the file's element count. This is
    ! required for SideExchange even in the serial single-rank case.
    call this%decomp%GenerateDecomposition(nElemFile,6*nElemFile)

    ! ---- 4. Allocate SELF Mesh3D_t and populate ----
    call this%Init(nGeo,nElemFile,6*nElemFile,nElemFile*ng1**3,nBCsLocal)
    this%nUniqueSides = 0 ! filled below after face matching
    this%quadrature = CHEBYSHEV_GAUSS_LOBATTO ! HOHQMesh face samples are at CGL points
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
      call build_nodeCoords_for_hex(this,e,nGeo,ismCorners,ismFlag,faceCurve,nodeXYZ)
      ! Synthesize globalNodeIDs: stamp the eight corners with their
      ! file IDs and leave interior IDs as 0 (interior nodes are
      ! private to the element under our tensor product layout).
      this%globalNodeIDs(:,:,:,e) = 0
      do l = 1,8
        i = this%CGNSCornerMap(1,l)
        j = this%CGNSCornerMap(2,l)
        k = this%CGNSCornerMap(3,l)
        this%globalNodeIDs(i,j,k,e) = ismCorners(l,e)
      enddo

      ! Pack elemInfo with simple placeholders; SELF's 3D path does not
      ! depend on the HOPR-style offset fields when the mesh comes from
      ! a non-HOPR reader.
      this%elemInfo(1,e) = 0
      this%elemInfo(2,e) = ismMat(e) ! Zone = material id
      this%elemInfo(3,e) = 6*(e-1)
      this%elemInfo(4,e) = 6*e
      this%elemInfo(5,e) = ng1**3*(e-1)
      this%elemInfo(6,e) = ng1**3*e
    enddo

    ! ---- 5. Build sideInfo via face corner matching ----
    ! Each face is keyed by its four sorted corner-node ids; two local
    ! faces with the same key are the two sides of an interior face.
    allocate(pairKey(1:4,1:6*nElemFile))
    allocate(pairElem(1:6*nElemFile),pairSide(1:6*nElemFile))
    do e = 1,nElemFile
      do iSide = 1,6
        do l = 1,4
          quadA(l) = ismCorners(this%sideMap(l,iSide),e)
        enddo
        call sort4(quadA,sortedA)
        i = iSide+6*(e-1)
        pairKey(1:4,i) = sortedA
        pairElem(i) = e
        pairSide(i) = iSide
      enddo
    enddo

    ! Hash chain by smallest corner node id
    hashSize = max(nNodesFile,6*nElemFile)+1
    allocate(hashHead(0:hashSize-1),hashNext(1:6*nElemFile))
    hashHead = 0
    hashNext = 0
    do i = 1,6*nElemFile
      hashKey = pairKey(1,i)
      bucket = modulo(hashKey,hashSize)
      hashNext(i) = hashHead(bucket)
      hashHead(bucket) = i
    enddo

    this%sideInfo = 0
    this%nUniqueSides = 0
    do e = 1,nElemFile
      do iSide = 1,6
        i = iSide+6*(e-1)
        ePair = 0; sPair = 0
        bucket = modulo(pairKey(1,i),hashSize)
        probe = hashHead(bucket)
        do while(probe /= 0)
          if(probe /= i) then
            if(all(pairKey(1:4,probe) == pairKey(1:4,i))) then
              ePair = pairElem(probe)
              sPair = pairSide(probe)
              exit
            endif
          endif
          probe = hashNext(probe)
        enddo

        flip = 0
        if(ePair /= 0) then
          ! Determine the flip by matching the corner-node orderings of
          ! the paired faces (face-index corner order per sideMap).
          do l = 1,4
            quadA(l) = ismCorners(this%sideMap(l,iSide),e)
            quadB(l) = ismCorners(this%sideMap(l,sPair),ePair)
          enddo
          flip = -1
          do p = 0,7
            if(quadB(flipPerm(1,p)) == quadA(1) .and. &
               quadB(flipPerm(2,p)) == quadA(2) .and. &
               quadB(flipPerm(3,p)) == quadA(3) .and. &
               quadB(flipPerm(4,p)) == quadA(4)) then
              flip = p
              exit
            endif
          enddo
          if(flip < 0) then
            print*,__FILE__//' : Inconsistent face corner ordering between elements ', &
              e,' and ',ePair
            stop 1
          endif
        endif

        this%sideInfo(3,iSide,e) = ePair
        this%sideInfo(4,iSide,e) = 10*sPair+flip
        this%sideInfo(5,iSide,e) = ismBCid(iSide,e)
        ! Allocate a globalSideID (count each shared face once)
        if(ePair == 0 .or. e < ePair) then
          this%nUniqueSides = this%nUniqueSides+1
          this%sideInfo(2,iSide,e) = this%nUniqueSides
        else
          this%sideInfo(2,iSide,e) = this%sideInfo(2,sPair,ePair)
        endif
      enddo
    enddo

    deallocate(hashHead,hashNext,pairKey,pairElem,pairSide)
    deallocate(nodeXYZ,ismCorners,ismFlag,ismBCid,ismMat,faceCurve)
    deallocate(BCNamesLocal,matNamesLocal)

    call this%UpdateDevice()

  contains

    subroutine sort4(a,b)
      !! Ascending insertion sort of four integers
      integer,intent(in) :: a(1:4)
      integer,intent(out) :: b(1:4)
      integer :: m,n,tmp
      b = a
      do m = 2,4
        tmp = b(m)
        n = m-1
        do while(n >= 1)
          if(b(n) <= tmp) exit
          b(n+1) = b(n)
          n = n-1
        enddo
        b(n+1) = tmp
      enddo
    endsubroutine sort4

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

  endsubroutine Read_HOHQMesh_Mesh3D_t

  subroutine build_nodeCoords_for_hex(mesh,e,nGeo,allCorners,allFlags,allFaces,nodeXYZ)
    !! Fill mesh%nodeCoords(:,:,:,:,e) for one hexahedral element by
    !! transfinite (Coons) interpolation of its six face grids. Faces
    !! flagged in the mesh file use the file's face-point grids; the
    !! remaining faces are bilinear patches of their four corner nodes
    !! evaluated at Chebyshev-Gauss-Lobatto parametric coordinates.
    !! Edge curves are extracted from the face grids (preferring a
    !! flagged face when an edge borders one flagged and one bilinear
    !! face) and the standard Boolean-sum formula
    !!   x = Pxi + Peta + Pzeta - Pxi*Peta - Pxi*Pzeta - Peta*Pzeta
    !!       + Pxi*Peta*Pzeta
    !! combines face, edge, and corner contributions. For an element
    !! with all-straight faces this reduces to trilinear interpolation
    !! of the eight corners.
    implicit none
    class(Mesh3D_t),intent(inout) :: mesh
    integer,intent(in) :: e
    integer,intent(in) :: nGeo
    integer,intent(in) :: allCorners(:,:) ! 8 x nElem corner-node ids
    integer,intent(in) :: allFlags(:,:) ! 6 x nElem face flags (SELF side order)
    real(prec),intent(in) :: allFaces(:,:,:,:,:) ! 3 x nGeo+1 x nGeo+1 x 6 x nElem
    real(prec),intent(in) :: nodeXYZ(:,:)
    ! Local
    integer :: ng1
    integer :: i,j,k,s,l
    integer :: cornerIDs(1:8)
    integer :: flag(1:6)
    real(prec) :: P(1:3,1:8)
    real(prec) :: face(1:3,1:nGeo+1,1:nGeo+1,1:6)
    real(prec) :: exEdge(1:3,1:nGeo+1,1:4) ! xi-directed edges: bs, bn, ts, tn
    real(prec) :: eyEdge(1:3,1:nGeo+1,1:4) ! eta-directed edges: bw, be, tw, te
    real(prec) :: ezEdge(1:3,1:nGeo+1,1:4) ! zeta-directed edges: sw, se, nw, ne
    real(prec) :: xi(0:nGeo),wts(0:nGeo)
    real(prec) :: u(1:nGeo+1)
    real(prec) :: a,b,c
    real(prec) :: xfp(1:3),xep(1:3),xcp(1:3)

    ng1 = nGeo+1
    cornerIDs = allCorners(1:8,e)
    flag = allFlags(1:6,e)

    do l = 1,8
      P(1:3,l) = nodeXYZ(1:3,cornerIDs(l))
    enddo

    if(nGeo == 1) then
      ! Pure trilinear element with the 8 corners
      do l = 1,8
        i = mesh%CGNSCornerMap(1,l)
        j = mesh%CGNSCornerMap(2,l)
        k = mesh%CGNSCornerMap(3,l)
        mesh%nodeCoords(1:3,i,j,k,e) = P(1:3,l)
      enddo
      return
    endif

    ! Chebyshev-Gauss-Lobatto parametric coordinates on [-1,1] and
    ! their [0,1] image used in the blending weights. We use SELF's
    ! quadrature routine to keep node placement consistent with the
    ! rest of the code.
    call ChebyshevQuadrature(nGeo,xi,wts,CHEBYSHEV_GAUSS_LOBATTO)
    do i = 1,ng1
      u(i) = 0.5_prec*(xi(i-1)+1.0_prec)
    enddo

    ! Face grids: curved faces from the file, straight faces as
    ! bilinear patches of their four corners (face-index corner order
    ! per sideMap; all faces use the natural on-face volume axes).
    do s = 1,6
      if(flag(s) == 1) then
        face(1:3,1:ng1,1:ng1,s) = allFaces(1:3,1:ng1,1:ng1,s,e)
      else
        do j = 1,ng1
          do i = 1,ng1
            a = u(i)
            b = u(j)
            face(1:3,i,j,s) = (1.0_prec-a)*(1.0_prec-b)*P(1:3,mesh%sideMap(1,s))+ &
                              a*(1.0_prec-b)*P(1:3,mesh%sideMap(2,s))+ &
                              a*b*P(1:3,mesh%sideMap(3,s))+ &
                              (1.0_prec-a)*b*P(1:3,mesh%sideMap(4,s))
          enddo
        enddo
      endif
    enddo

    ! Edge curves extracted from the face grids. Each edge borders two
    ! faces; a flagged (curved) face is preferred so that straight-face
    ! bilinear patches never override curved edge data. When both
    ! bordering faces are straight, the two restrictions coincide (the
    ! straight line between the edge's corner nodes). SELF side ids:
    ! 1=bottom, 2=south, 3=east, 4=north, 5=west, 6=top.
    do i = 1,ng1
      ! xi-directed edges
      if(flag(1) == 1 .or. flag(2) /= 1) then ! bottom-south: bottom or south
        exEdge(1:3,i,1) = face(1:3,i,1,1)
      else
        exEdge(1:3,i,1) = face(1:3,i,1,2)
      endif
      if(flag(1) == 1 .or. flag(4) /= 1) then ! bottom-north: bottom or north
        exEdge(1:3,i,2) = face(1:3,i,ng1,1)
      else
        exEdge(1:3,i,2) = face(1:3,i,1,4)
      endif
      if(flag(6) == 1 .or. flag(2) /= 1) then ! top-south: top or south
        exEdge(1:3,i,3) = face(1:3,i,1,6)
      else
        exEdge(1:3,i,3) = face(1:3,i,ng1,2)
      endif
      if(flag(6) == 1 .or. flag(4) /= 1) then ! top-north: top or north
        exEdge(1:3,i,4) = face(1:3,i,ng1,6)
      else
        exEdge(1:3,i,4) = face(1:3,i,ng1,4)
      endif
      ! eta-directed edges
      if(flag(1) == 1 .or. flag(5) /= 1) then ! bottom-west: bottom or west
        eyEdge(1:3,i,1) = face(1:3,1,i,1)
      else
        eyEdge(1:3,i,1) = face(1:3,i,1,5)
      endif
      if(flag(1) == 1 .or. flag(3) /= 1) then ! bottom-east: bottom or east
        eyEdge(1:3,i,2) = face(1:3,ng1,i,1)
      else
        eyEdge(1:3,i,2) = face(1:3,i,1,3)
      endif
      if(flag(6) == 1 .or. flag(5) /= 1) then ! top-west: top or west
        eyEdge(1:3,i,3) = face(1:3,1,i,6)
      else
        eyEdge(1:3,i,3) = face(1:3,i,ng1,5)
      endif
      if(flag(6) == 1 .or. flag(3) /= 1) then ! top-east: top or east
        eyEdge(1:3,i,4) = face(1:3,ng1,i,6)
      else
        eyEdge(1:3,i,4) = face(1:3,i,ng1,3)
      endif
      ! zeta-directed edges
      if(flag(2) == 1 .or. flag(5) /= 1) then ! south-west: south or west
        ezEdge(1:3,i,1) = face(1:3,1,i,2)
      else
        ezEdge(1:3,i,1) = face(1:3,1,i,5)
      endif
      if(flag(2) == 1 .or. flag(3) /= 1) then ! south-east: south or east
        ezEdge(1:3,i,2) = face(1:3,ng1,i,2)
      else
        ezEdge(1:3,i,2) = face(1:3,1,i,3)
      endif
      if(flag(4) == 1 .or. flag(5) /= 1) then ! north-west: north or west
        ezEdge(1:3,i,3) = face(1:3,1,i,4)
      else
        ezEdge(1:3,i,3) = face(1:3,ng1,i,5)
      endif
      if(flag(4) == 1 .or. flag(3) /= 1) then ! north-east: north or east
        ezEdge(1:3,i,4) = face(1:3,ng1,i,4)
      else
        ezEdge(1:3,i,4) = face(1:3,ng1,i,3)
      endif
    enddo

    ! Boolean-sum transfinite interpolation. u,v,w in [0,1] are the
    ! blending parameters along xi, eta, zeta.
    do k = 1,ng1
      c = u(k)
      do j = 1,ng1
        b = u(j)
        do i = 1,ng1
          a = u(i)

          ! Face contributions (Pxi + Peta + Pzeta)
          xfp = (1.0_prec-a)*face(1:3,j,k,5)+a*face(1:3,j,k,3)+ &
                (1.0_prec-b)*face(1:3,i,k,2)+b*face(1:3,i,k,4)+ &
                (1.0_prec-c)*face(1:3,i,j,1)+c*face(1:3,i,j,6)

          ! Edge corrections (Pxi*Peta + Pxi*Pzeta + Peta*Pzeta)
          xep = (1.0_prec-a)*(1.0_prec-b)*ezEdge(1:3,k,1)+ &
                a*(1.0_prec-b)*ezEdge(1:3,k,2)+ &
                a*b*ezEdge(1:3,k,4)+ &
                (1.0_prec-a)*b*ezEdge(1:3,k,3)+ &
                (1.0_prec-a)*(1.0_prec-c)*eyEdge(1:3,j,1)+ &
                a*(1.0_prec-c)*eyEdge(1:3,j,2)+ &
                a*c*eyEdge(1:3,j,4)+ &
                (1.0_prec-a)*c*eyEdge(1:3,j,3)+ &
                (1.0_prec-b)*(1.0_prec-c)*exEdge(1:3,i,1)+ &
                b*(1.0_prec-c)*exEdge(1:3,i,2)+ &
                b*c*exEdge(1:3,i,4)+ &
                (1.0_prec-b)*c*exEdge(1:3,i,3)

          ! Corner contributions (Pxi*Peta*Pzeta)
          xcp = (1.0_prec-a)*(1.0_prec-b)*(1.0_prec-c)*P(1:3,1)+ &
                a*(1.0_prec-b)*(1.0_prec-c)*P(1:3,2)+ &
                a*b*(1.0_prec-c)*P(1:3,3)+ &
                (1.0_prec-a)*b*(1.0_prec-c)*P(1:3,4)+ &
                (1.0_prec-a)*(1.0_prec-b)*c*P(1:3,5)+ &
                a*(1.0_prec-b)*c*P(1:3,6)+ &
                a*b*c*P(1:3,7)+ &
                (1.0_prec-a)*b*c*P(1:3,8)

          mesh%nodeCoords(1:3,i,j,k,e) = xfp-xep+xcp

        enddo
      enddo
    enddo

  endsubroutine build_nodeCoords_for_hex

  subroutine Write_Mesh3D_t(this,meshFile)
    ! Writes mesh output in HOPR format (serial only)
    implicit none
    class(Mesh3D_t),intent(inout) :: this
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

  endsubroutine Write_Mesh3D_t

endmodule SELF_Mesh_3D_t
