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

module SELF_Mesh_2D

  use SELF_Mesh_2D_t
  use SELF_GPU
  use iso_c_binding

  implicit none

  type,extends(Mesh2D_t) :: Mesh2D
    type(c_ptr) :: sideInfo_gpu

  contains
    procedure,public :: Init => Init_Mesh2D
    procedure,public :: Free => Free_Mesh2D
    procedure,public :: Read_HOPr => Read_HOPr_Mesh2D

  endtype Mesh2D

contains

  subroutine Init_Mesh2D(this,nGeo,nElem,nSides,nNodes,nBCs)
    implicit none
    class(Mesh2D),intent(inout) :: this
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

    call gpuCheck(hipMalloc(this%sideInfo_gpu,sizeof(this%sideInfo)))

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
    call this%decomp%Free()

    call gpuCheck(hipFree(this%sideInfo_gpu))

  endsubroutine Free_Mesh2D

  subroutine Read_HOPr_Mesh2D(this,meshFile,enableDomainDecomposition)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 6
    ! Adapted for 2D Mesh : Note that HOPR does not have 2D mesh output.
    implicit none
    class(Mesh2D),intent(out) :: this
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

    if(present(enableDomainDecomposition)) then
      call this%decomp%init(enableDomainDecomposition)
    else
      call this%decomp%init(.false.)
    endif

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
    call gpuCheck(hipMemcpy(this%sideInfo_gpu,c_loc(this%sideInfo),sizeof(this%sideInfo),hipMemcpyHostToDevice))

  endsubroutine Read_HOPr_Mesh2D

endmodule SELF_Mesh_2D
