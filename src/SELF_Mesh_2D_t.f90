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
  use SELF_MPI

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

    procedure,public :: ResetBoundaryConditionType => ResetBoundaryConditionType_Mesh2D_t

    procedure,public :: Read_HOPr => Read_HOPr_Mesh2D_t

    procedure,public :: Write_Mesh => Write_Mesh2D_t

    procedure,public :: RecalculateFlip => RecalculateFlip_Mesh2D_t

  endtype Mesh2D_t

contains

  subroutine Init_Mesh2D_t(this,nGeo,nElem,nSides,nNodes,nBCs)
    implicit none
    class(Mesh2D_t),intent(out) :: this
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

  endsubroutine Free_Mesh2D_t

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

  endsubroutine ResetBoundaryConditionType_Mesh2D_t

  subroutine Read_HOPr_Mesh2D_t(this,meshFile,decomp)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    ! Adapted for 2D Mesh : Note that HOPR does not have 2D mesh output.
    implicit none
    class(Mesh2D_t),intent(out) :: this
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

  endsubroutine Read_HOPr_Mesh2D_t

  subroutine RecalculateFlip_Mesh2D_t(this,decomp)
    implicit none
    class(Mesh2D_t),intent(inout) :: this
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
