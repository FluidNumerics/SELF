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

module SELF_Mesh_3D

  use SELF_Mesh_3D_t
  use SELF_GPU
  use iso_c_binding

  implicit none

  type,extends(Mesh3D_t) :: Mesh3D
    type(c_ptr) :: sideInfo_gpu

  contains
    procedure,public :: Init => Init_Mesh3D
    procedure,public :: Free => Free_Mesh3D
    procedure,public :: UpdateDevice => UpdateDevice_Mesh3D

  endtype Mesh3D

contains

  subroutine Init_Mesh3D(this,nGeo,nElem,nSides,nNodes,nBCs)
    implicit none
    class(Mesh3D),intent(inout) :: this
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

    call gpuCheck(hipMalloc(this%sideInfo_gpu,sizeof(this%sideInfo)))

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
    deallocate(this%BCNames)
    call this%decomp%Free()

    call gpuCheck(hipFree(this%sideInfo_gpu))

  endsubroutine Free_Mesh3D

  subroutine UpdateDevice_Mesh3D(this)
    implicit none
    class(Mesh3D),intent(inout) :: this

    call gpuCheck(hipMemcpy(this%sideInfo_gpu,c_loc(this%sideInfo),sizeof(this%sideInfo),hipMemcpyHostToDevice))

  endsubroutine UpdateDevice_Mesh3D

endmodule SELF_Mesh_3D
