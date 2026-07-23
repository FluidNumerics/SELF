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

module SELF_DomainDecomposition

  use SELF_DomainDecomposition_t
  use mpi
  use iso_c_binding
  use iso_fortran_env,only:int64

  implicit none

  type,extends(DomainDecomposition_t) :: DomainDecomposition
    type(c_ptr) :: elemToRank_gpu

    ! Aggregated MPI halo exchange tables, shared by every field on the mesh.
    ! All locally-owned sides with a neighbor element on another rank are
    ! enumerated once, grouped by neighbor rank and sorted by global side id
    ! within each group, so that the two ranks sharing an interface enumerate
    ! its sides in the same order and packed message buffers line up without
    ! any per-message metadata. Field classes allocate their own packed
    ! send/receive buffers (sized by their variable count) and share these
    ! tables for the pack/unpack kernels and message posting.
    logical :: halo_built = .false.
    integer :: halo_nnbr = 0 ! Number of neighboring ranks
    integer :: halo_nsides = 0 ! Total number of sides exchanged with other ranks
    integer,allocatable :: halo_rank(:) ! (1:halo_nnbr) neighbor rank ids
    integer,allocatable :: halo_offset(:) ! (1:halo_nnbr+1) side-list offsets per neighbor
    type(c_ptr) :: halo_sides_gpu = c_null_ptr ! (element,side,flip) triplets, device copy

  contains

    procedure :: Init => Init_DomainDecomposition
    procedure :: Free => Free_DomainDecomposition

    procedure :: SetElemToRank => SetElemToRank_DomainDecomposition

    procedure :: BuildHaloExchange => BuildHaloExchange_DomainDecomposition

  endtype DomainDecomposition

contains

  subroutine Init_DomainDecomposition(this)
    implicit none
    class(DomainDecomposition),intent(inout) :: this
    ! Local
    integer       :: ierror
    integer(c_int) :: num_devices,hip_err,device_id

    this%mpiComm = 0
    this%mpiPrec = prec
    this%rankId = 0
    this%nRanks = 1
    this%nElem = 0
    this%mpiEnabled = .false.

    this%mpiComm = MPI_COMM_WORLD
    print*,__FILE__," : Initializing MPI"
    call mpi_init(ierror)
    call mpi_comm_rank(this%mpiComm,this%rankId,ierror)
    call mpi_comm_size(this%mpiComm,this%nRanks,ierror)
    print*,__FILE__," : Rank ",this%rankId+1,"/",this%nRanks," checking in."

    if(this%nRanks > 1) then
      this%mpiEnabled = .true.
    else
      print*,__FILE__," : No domain decomposition used."
    endif

    if(prec == real32) then
      this%mpiPrec = MPI_FLOAT
    else
      this%mpiPrec = MPI_DOUBLE
    endif

    allocate(this%offsetElem(1:this%nRanks+1))

    hip_err = hipGetDeviceCount(num_devices)
    if(hip_err /= 0) then
      print*,'Failed to get device count on rank',this%rankId
      call MPI_Abort(MPI_COMM_WORLD,hip_err,ierror)
    endif

    ! Assign GPU device ID based on MPI rank
    device_id = modulo(this%rankId,num_devices) ! Assumes that mpi ranks are packed sequentially on a node until the node is filled up.
    hip_err = hipSetDevice(device_id)
    print*,__FILE__," : Rank ",this%rankId+1," assigned to device ",device_id
    if(hip_err /= 0) then
      print*,'Failed to set device for rank',this%rankId,'to device',device_id
      call MPI_Abort(MPI_COMM_WORLD,hip_err,ierror)
    endif

    this%initialized = .true.

  endsubroutine Init_DomainDecomposition
  subroutine Free_DomainDecomposition(this)
    implicit none
    class(DomainDecomposition),intent(inout) :: this
    ! Local
    integer :: ierror

    if(associated(this%offSetElem)) then
      deallocate(this%offSetElem)
    endif
    if(associated(this%elemToRank)) then
      deallocate(this%elemToRank)
      call gpuCheck(hipFree(this%elemToRank_gpu))
    endif

    if(allocated(this%requests)) deallocate(this%requests)
    if(allocated(this%stats)) deallocate(this%stats)

    if(allocated(this%halo_rank)) deallocate(this%halo_rank)
    if(allocated(this%halo_offset)) deallocate(this%halo_offset)
    if(c_associated(this%halo_sides_gpu)) then
      call gpuCheck(hipFree(this%halo_sides_gpu))
      this%halo_sides_gpu = c_null_ptr
    endif
    this%halo_built = .false.

    print*,__FILE__," : Rank ",this%rankId+1,"/",this%nRanks," checking out."
    call MPI_FINALIZE(ierror)

  endsubroutine Free_DomainDecomposition

  subroutine SetElemToRank_DomainDecomposition(this,nElem)
    implicit none
    class(DomainDecomposition),intent(inout) :: this
    integer,intent(in) :: nElem
    ! Local
    integer :: iel

    this%nElem = nElem

    allocate(this%elemToRank(1:nelem))
    call gpuCheck(hipMalloc(this%elemToRank_gpu,sizeof(this%elemToRank)))

    call DomainDecomp(nElem, &
                      this%nRanks, &
                      this%offSetElem)

    do iel = 1,nElem
      call ElemToRank(this%nRanks, &
                      this%offSetElem, &
                      iel, &
                      this%elemToRank(iel))
    enddo
    call gpuCheck(hipMemcpy(this%elemToRank_gpu,c_loc(this%elemToRank),sizeof(this%elemToRank),hipMemcpyHostToDevice))

  endsubroutine SetElemToRank_DomainDecomposition

  subroutine BuildHaloExchange_DomainDecomposition(this,sideInfo,nElem,nSidesPerElem)
  !! Build the aggregated halo exchange tables for the mesh this domain
  !! decomposition belongs to.
  !!
  !! Enumerates every locally-owned side whose neighbor element resides on
  !! another rank, groups the sides by neighbor rank, and sorts each group by
  !! global side id. Both ranks sharing an interface enumerate the interface
  !! sides in the same order (the global side id is shared), so a packed send
  !! buffer on one rank and the packed receive buffer on the other line up
  !! without any per-message metadata. The (element,side,flip) triplets are
  !! copied to the device for the field pack/unpack kernels.
  !!
  !!  Input
  !!    - sideInfo : the mesh sideInfo array, (1:5, 1:nSidesPerElem, 1:nElem)
  !!    - nElem : number of local elements
  !!    - nSidesPerElem : 4 for 2-D (quadrilateral), 6 for 3-D (hexahedral)
    implicit none
    class(DomainDecomposition),intent(inout) :: this
    integer,intent(in) :: nElem
    integer,intent(in) :: nSidesPerElem
    integer,intent(in) :: sideInfo(1:5,1:nSidesPerElem,1:nElem)
    ! Local
    integer :: e1,s1,e2,r2,n,nhalo
    integer,allocatable :: elems(:),sides(:),flips(:),ranks(:)
    integer,allocatable,target :: halosides(:)
    integer(int64),allocatable :: keys(:)

    nhalo = 0
    do e1 = 1,nElem
      do s1 = 1,nSidesPerElem
        e2 = sideInfo(3,s1,e1)
        if(e2 > 0) then
          if(this%elemToRank(e2) /= this%rankId) then
            nhalo = nhalo+1
          endif
        endif
      enddo
    enddo

    this%halo_nsides = nhalo
    this%halo_nnbr = 0
    if(nhalo == 0) then
      this%halo_built = .true.
      return
    endif

    allocate(elems(1:nhalo),sides(1:nhalo),flips(1:nhalo), &
             ranks(1:nhalo),keys(1:nhalo))

    n = 0
    do e1 = 1,nElem
      do s1 = 1,nSidesPerElem
        e2 = sideInfo(3,s1,e1)
        if(e2 > 0) then
          r2 = this%elemToRank(e2)
          if(r2 /= this%rankId) then
            n = n+1
            elems(n) = e1
            sides(n) = s1
            flips(n) = sideInfo(4,s1,e1)-10*(sideInfo(4,s1,e1)/10)
            ranks(n) = r2
            ! Composite sort key: neighbor rank (major), global side id (minor).
            keys(n) = int(r2,int64)*2147483648_int64+ &
                      int(abs(sideInfo(2,s1,e1)),int64)
          endif
        endif
      enddo
    enddo

    call HaloKeySort(keys,elems,sides,flips,ranks,nhalo)

    ! Count the neighboring ranks and record the segment offsets
    this%halo_nnbr = 1
    do n = 2,nhalo
      if(ranks(n) /= ranks(n-1)) this%halo_nnbr = this%halo_nnbr+1
    enddo

    allocate(this%halo_rank(1:this%halo_nnbr))
    allocate(this%halo_offset(1:this%halo_nnbr+1))

    this%halo_nnbr = 1
    this%halo_rank(1) = ranks(1)
    this%halo_offset(1) = 0
    do n = 2,nhalo
      if(ranks(n) /= ranks(n-1)) then
        this%halo_nnbr = this%halo_nnbr+1
        this%halo_rank(this%halo_nnbr) = ranks(n)
        this%halo_offset(this%halo_nnbr) = n-1
      endif
    enddo
    this%halo_offset(this%halo_nnbr+1) = nhalo

    ! Flatten the (element,side,flip) triplets with 0-based element and side
    ! indices for the device pack/unpack kernels
    allocate(halosides(1:3*nhalo))
    do n = 1,nhalo
      halosides(3*n-2) = elems(n)-1
      halosides(3*n-1) = sides(n)-1
      halosides(3*n) = flips(n)
    enddo

    call gpuCheck(hipMalloc(this%halo_sides_gpu,sizeof(halosides)))
    call gpuCheck(hipMemcpy(this%halo_sides_gpu,c_loc(halosides), &
                            sizeof(halosides),hipMemcpyHostToDevice))

    deallocate(elems,sides,flips,ranks,keys,halosides)

    this%halo_built = .true.

  endsubroutine BuildHaloExchange_DomainDecomposition

  subroutine HaloKeySort(keys,elems,sides,flips,ranks,n)
  !! Heap sort of the halo side list by ascending 64-bit key, permuting the
  !! companion arrays alongside the keys. O(n log n), in place; the relative
  !! order of equal keys is irrelevant because global side ids are unique.
    implicit none
    integer,intent(in) :: n
    integer(int64),intent(inout) :: keys(1:n)
    integer,intent(inout) :: elems(1:n),sides(1:n),flips(1:n),ranks(1:n)
    ! Local
    integer :: i,last

    do i = n/2,1,-1
      call HaloSiftDown(keys,elems,sides,flips,ranks,i,n)
    enddo
    do last = n,2,-1
      call HaloSwap(keys,elems,sides,flips,ranks,1,last)
      call HaloSiftDown(keys,elems,sides,flips,ranks,1,last-1)
    enddo

  endsubroutine HaloKeySort

  subroutine HaloSiftDown(keys,elems,sides,flips,ranks,start,last)
    implicit none
    integer(int64),intent(inout) :: keys(:)
    integer,intent(inout) :: elems(:),sides(:),flips(:),ranks(:)
    integer,intent(in) :: start,last
    ! Local
    integer :: root,child

    root = start
    do while(2*root <= last)
      child = 2*root
      if(child < last) then
        if(keys(child) < keys(child+1)) child = child+1
      endif
      if(keys(root) < keys(child)) then
        call HaloSwap(keys,elems,sides,flips,ranks,root,child)
        root = child
      else
        return
      endif
    enddo

  endsubroutine HaloSiftDown

  subroutine HaloSwap(keys,elems,sides,flips,ranks,i,j)
    implicit none
    integer(int64),intent(inout) :: keys(:)
    integer,intent(inout) :: elems(:),sides(:),flips(:),ranks(:)
    integer,intent(in) :: i,j
    ! Local
    integer(int64) :: k
    integer :: t

    k = keys(i); keys(i) = keys(j); keys(j) = k
    t = elems(i); elems(i) = elems(j); elems(j) = t
    t = sides(i); sides(i) = sides(j); sides(j) = t
    t = flips(i); flips(i) = flips(j); flips(j) = t
    t = ranks(i); ranks(i) = ranks(j); ranks(j) = t

  endsubroutine HaloSwap

endmodule SELF_DomainDecomposition
