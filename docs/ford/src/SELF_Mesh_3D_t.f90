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

  subroutine Read_HOPr_Mesh3D_t(this,meshFile,enableDomainDecomposition)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    implicit none
    class(Mesh3D_t),intent(out) :: this
    character(*),intent(in) :: meshFile
    logical,intent(in),optional :: enableDomainDecomposition
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

    if(present(enableDomainDecomposition)) then
      call this%decomp%init(enableDomainDecomposition)
    else
      call this%decomp%init(.false.)
    endif

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
