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
  use SELF_SupportRoutines
  use SELF_HDF5
  use SELF_Mesh
  use SELF_DomainDecomposition

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

  contains

    procedure,public :: Init => Init_Mesh3D_t
    procedure,public :: Free => Free_Mesh3D_t
    procedure,public :: UpdateDevice => UpdateDevice_Mesh3D_t

    generic,public :: StructuredMesh => UniformStructuredMesh_Mesh3D_t
    procedure,private :: UniformStructuredMesh_Mesh3D_t

    procedure,public :: Read_HOPr => Read_HOPr_Mesh3D_t

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
    call this%decomp%Free()

  endsubroutine Free_Mesh3D_t

  subroutine UpdateDevice_Mesh3D_t(this)
    implicit none
    class(Mesh3D_t),intent(inout) :: this

    return

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

  endfunction elementid

  subroutine UniformStructuredMesh_Mesh3D_t(this,nxPerTile,nyPerTile,nzPerTile, &
                                            nTileX,nTileY,nTileZ,dx,dy,dz,bcids)
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
  !!    - enableDomainDecomposition : Boolean to determine if domain decomposition is used.
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

    call this%decomp%init()

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

  subroutine Read_HOPr_Mesh3D_t(this,meshFile)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    implicit none
    class(Mesh3D_t),intent(out) :: this
    character(*),intent(in) :: meshFile
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

    call this%decomp%init()

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
