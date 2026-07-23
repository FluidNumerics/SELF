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

module SELF_Mesh_1D

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Data
  use SELF_Scalar_1D
  use SELF_SupportRoutines
  use SELF_HDF5
  use SELF_Mesh

  ! External Libs !
  use HDF5

  use iso_c_binding

  implicit none

  type,extends(SEMMesh) :: Mesh1D
    integer,pointer,dimension(:,:) :: elemInfo
    real(prec),pointer,dimension(:) :: nodeCoords
    integer,pointer,dimension(:) :: globalNodeIDs
    integer,pointer,dimension(:,:) :: BCType
    character(LEN=255),allocatable :: BCNames(:)
    integer,dimension(2) :: bcid = 0 ! Boundary conditions for the left and right endpoints

  contains
    procedure,public :: Init => Init_Mesh1D
    procedure,public :: Free => Free_Mesh1D
    generic,public :: StructuredMesh => UniformBlockMesh_Mesh1D
    procedure,private ::  UniformBlockMesh_Mesh1D
    procedure,public :: ResetBoundaryConditionType => ResetBoundaryConditionType_Mesh1D

    procedure,public  :: Write_Mesh => Write_Mesh1D

  endtype Mesh1D

contains

  subroutine Init_Mesh1D(this,nElem,nNodes,nBCs)
    implicit none
    class(Mesh1D),intent(out) :: this
    integer,intent(in) :: nElem
    integer,intent(in) :: nNodes
    integer,intent(in) :: nBCs

    this%nGeo = 1
    this%nElem = nElem
    this%nGlobalElem = nElem
    this%nNodes = nNodes
    this%nCornerNodes = nElem*2
    this%nUniqueNodes = 0
    this%nBCs = nBCs
    this%bcid = 0

    allocate(this%elemInfo(1:4,1:nElem))
    allocate(this%nodeCoords(1:nNodes))
    allocate(this%globalNodeIDs(1:nNodes))
    allocate(this%BCType(1:4,1:nBCs))

    allocate(this%BCNames(1:nBCs))
    call this%decomp%Init()

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
    deallocate(this%BCNames)
    call this%decomp%Free()

  endsubroutine Free_Mesh1D

  subroutine UniformBlockMesh_Mesh1D(this,nElem,x)
    implicit none
    class(Mesh1D),intent(out) :: this
    integer,intent(in) :: nElem
    real(prec),intent(in) :: x(1:2)
    ! Local
    integer :: iel,ngeo
    integer :: nid,nNodes
    integer :: i
    real(prec) :: xU(1:nElem+1)
    type(Lagrange),target :: linearInterp
    type(Lagrange),target :: nGeoInterp
    type(Scalar1D) :: xLinear
    type(Scalar1D) :: xGeo

    ngeo = 1

    nNodes = nElem*(nGeo+1)
    call this%Init(nElem,nNodes,2)
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

    call xLinear%GridInterp(xGeo%interior)

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

  subroutine ResetBoundaryConditionType_Mesh1D(this,leftbc,rightbc)
    !! This method can be used to reset all of the boundary elements
    !! boundary condition type to the desired value.
    !!
    !! Note that ALL physical boundaries will be set to have this boundary
    !! condition
    implicit none
    class(Mesh1D),intent(inout) :: this
    integer,intent(in) ::leftbc,rightbc

    this%bcid(1) = leftbc
    this%bcid(2) = rightbc

  endsubroutine ResetBoundaryConditionType_Mesh1D

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

endmodule SELF_Mesh_1D
