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

module SELF_DomainDecomposition_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_SupportRoutines
  use mpi
  use iso_c_binding

  implicit none

  type DomainDecomposition_t
    logical :: mpiEnabled = .false.
    logical :: initialized = .false.
    logical :: ownsMpi = .false. ! true when SELF called mpi_init and is responsible for MPI_Finalize
    integer :: mpiComm = MPI_COMM_NULL
    integer :: mpiPrec
    integer :: rankId
    integer :: nRanks
    integer :: nElem
    integer :: maxMsg
    integer :: msgCount
    integer,pointer,dimension(:) :: elemToRank
    integer,pointer,dimension(:) :: offSetElem
    integer,allocatable :: requests(:)
    integer,allocatable :: stats(:,:)

  contains

    procedure :: Init => Init_DomainDecomposition_t
    procedure :: Free => Free_DomainDecomposition_t

    procedure :: GenerateDecomposition => GenerateDecomposition_DomainDecomposition_t
    procedure :: SetElemToRank => SetElemToRank_DomainDecomposition_t
    procedure :: ReserveMessages => ReserveMessages_DomainDecomposition_t
    procedure :: CountRemoteSides => CountRemoteSides_DomainDecomposition_t
    procedure :: CountRemoteMortarSides => CountRemoteMortarSides_DomainDecomposition_t

    procedure,public :: FinalizeMPIExchangeAsync

  endtype DomainDecomposition_t

  ! Process-wide bookkeeping for the MPI lifecycle.
  !
  ! Several decompositions can be live at once -- a process may hold more than
  ! one mesh, which is the normal case for a long-lived Python driver. SELF may
  ! only finalize MPI that it initialized itself, and only once the last
  ! decomposition has been freed; otherwise the first mesh to be torn down
  ! pulls MPI out from under the meshes that are still in use.
  integer,private :: nLiveDecomps = 0
  logical,private :: selfInitializedMpi = .false.

contains

  subroutine AcquireMPI(mpiComm,ownsMpi,comm)
  !! Resolve the communicator a new decomposition will run on and register it
  !! with the process-wide lifecycle bookkeeping.
  !!
  !! With `comm` present the caller (e.g. mpi4py) already owns MPI and SELF
  !! calls neither MPI_Init nor MPI_Finalize. Without it SELF falls back to
  !! MPI_COMM_WORLD, initializing MPI only if nobody else has yet.
  !!
  !! Shared by the CPU and GPU DomainDecomposition Init implementations so the
  !! two cannot drift apart.
    implicit none
    integer,intent(out) :: mpiComm
    logical,intent(out) :: ownsMpi
    integer,intent(in),optional :: comm
    ! Local
    integer :: ierror
    logical :: mpiIsInitialized

    call MPI_Initialized(mpiIsInitialized,ierror)

    if(present(comm)) then
      ! The caller (e.g. mpi4py) owns the MPI lifecycle and provides the communicator.
      if(.not. mpiIsInitialized) then
        error stop __FILE__//" : A communicator was provided but MPI is not initialized."// &
          " Initialize MPI (e.g. via mpi4py or MPI_Init) before creating a mesh with an external communicator."
      endif
      mpiComm = comm
    else
      mpiComm = MPI_COMM_WORLD
      if(.not. mpiIsInitialized) then
        print*,__FILE__," : Initializing MPI"
        call mpi_init(ierror)
        selfInitializedMpi = .true.
      endif
    endif

    nLiveDecomps = nLiveDecomps+1
    ownsMpi = selfInitializedMpi

  endsubroutine AcquireMPI

  subroutine ReleaseMPI(ownsMpi,rankId,nRanks)
  !! Retire one decomposition from the process-wide lifecycle bookkeeping and,
  !! if it was the last one and SELF initialized MPI, finalize.
  !!
  !! Shared by the CPU and GPU DomainDecomposition Free implementations.
    implicit none
    logical,intent(inout) :: ownsMpi
    integer,intent(in) :: rankId
    integer,intent(in) :: nRanks
    ! Local
    integer :: ierror
    logical :: mpiIsFinalized

    ownsMpi = .false.
    nLiveDecomps = max(nLiveDecomps-1,0)

    ! Other meshes are still live; they need MPI to stay up.
    if(nLiveDecomps > 0) return

    call MPI_Finalized(mpiIsFinalized,ierror)
    if(selfInitializedMpi .and. .not. mpiIsFinalized) then
      print*,__FILE__," : Rank ",rankId+1,"/",nRanks," checking out."
      call MPI_FINALIZE(ierror)
    endif
    selfInitializedMpi = .false.

  endsubroutine ReleaseMPI

  subroutine Init_DomainDecomposition_t(this,comm)
    implicit none
    class(DomainDecomposition_t),intent(inout) :: this
    integer,intent(in),optional :: comm
    ! Local
    integer       :: ierror

    this%mpiComm = MPI_COMM_NULL
    this%mpiPrec = prec
    this%rankId = 0
    this%nRanks = 1
    this%nElem = 0
    this%mpiEnabled = .false.
    this%ownsMpi = .false.

    call AcquireMPI(this%mpiComm,this%ownsMpi,comm)

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

    this%initialized = .true.

  endsubroutine Init_DomainDecomposition_t

  subroutine Free_DomainDecomposition_t(this)
    implicit none
    class(DomainDecomposition_t),intent(inout) :: this

    if(associated(this%offSetElem)) then
      deallocate(this%offSetElem)
    endif
    if(associated(this%elemToRank)) then
      deallocate(this%elemToRank)
    endif

    if(allocated(this%requests)) deallocate(this%requests)
    if(allocated(this%stats)) deallocate(this%stats)

    ! Guard against a second Free, and against freeing a decomposition that was
    ! never initialized, both of which would corrupt the live-decomposition count.
    if(this%initialized) then
      call ReleaseMPI(this%ownsMpi,this%rankId,this%nRanks)
      this%initialized = .false.
    endif

  endsubroutine Free_DomainDecomposition_t

  subroutine GenerateDecomposition_DomainDecomposition_t(this,nGlobalElem,maxMsg)
    implicit none
    class(DomainDecomposition_t),intent(inout) :: this
    integer,intent(in) :: nGlobalElem
    integer,intent(in) :: maxMsg

    call this%setElemToRank(nGlobalElem)
    if(allocated(this%requests)) deallocate(this%requests)
    if(allocated(this%stats)) deallocate(this%stats)

    allocate(this%requests(1:maxMsg))
    allocate(this%stats(MPI_STATUS_SIZE,1:maxMsg))
    this%maxMsg = maxMsg

    print*,__FILE__//" : Rank ",this%rankId+1," : n_elements = ", &
      this%offSetElem(this%rankId+2)-this%offSetElem(this%rankId+1)

  endsubroutine GenerateDecomposition_DomainDecomposition_t

  subroutine ReserveMessages_DomainDecomposition_t(this,nMsg)
  !! Ensures the non-blocking request and status scratch can hold `nMsg` messages
  !! in flight, growing it when it cannot.
  !!
  !! GenerateDecomposition sizes the scratch from the mesh alone (one entry per
  !! unique side), but an exchange posts a send and a receive per rank-remote side
  !! *per variable* -- and, for vectors, per direction. The message count is
  !! therefore a multiple of the side count that the mesh cannot know, and a
  !! solution with enough variables overruns the scratch. Each exchange counts the
  !! sides it will message (CountRemoteSides / CountRemoteMortarSides) and reserves
  !! from that count before posting; the call is a comparison once the scratch is
  !! large enough.
  !!
  !! Reserving from the rank-remote side count rather than from the local side
  !! count matters at scale: the remote sides are the partition boundary, while the
  !! local sides are the whole subdomain, and each entry carries an
  !! MPI_STATUS_SIZE row of `stats` with it.
  !!
  !! No message may be in flight when this is called: the scratch is reallocated,
  !! not grown in place. Every exchange posts and then waits (through
  !! FinalizeMPIExchangeAsync) before another can start, so a reservation at the
  !! head of an exchange is safe.
    implicit none
    class(DomainDecomposition_t),intent(inout) :: this
    integer,intent(in) :: nMsg

    if(nMsg <= this%maxMsg) return

    if(allocated(this%requests)) deallocate(this%requests)
    if(allocated(this%stats)) deallocate(this%stats)

    allocate(this%requests(1:nMsg))
    allocate(this%stats(MPI_STATUS_SIZE,1:nMsg))
    this%maxMsg = nMsg

  endsubroutine ReserveMessages_DomainDecomposition_t

  integer function CountRemoteSides_DomainDecomposition_t(this,sideInfo,nElem,nSidesPerElem) result(nRemote)
  !! Counts the conforming sides of the local elements whose neighbor element is
  !! owned by another rank -- the sides an exchange sends and receives over MPI.
  !!
  !! The scan is over the side table alone (a few integers per side), which is
  !! negligible next to packing and sending the traces those sides carry, and it
  !! stays correct when the mesh is re-partitioned under adaptive refinement.
  !!
  !!  Input
  !!    - sideInfo : the mesh sideInfo array, (1:5, 1:nSidesPerElem, 1:nElem)
  !!    - nElem : number of local elements
  !!    - nSidesPerElem : 4 for 2-D (quadrilateral), 6 for 3-D (hexahedral)
    implicit none
    class(DomainDecomposition_t),intent(in) :: this
    integer,intent(in) :: nElem
    integer,intent(in) :: nSidesPerElem
    integer,intent(in) :: sideInfo(1:5,1:nSidesPerElem,1:nElem)
    ! Local
    integer :: e1,s1,e2

    nRemote = 0
    do e1 = 1,nElem
      do s1 = 1,nSidesPerElem
        e2 = sideInfo(3,s1,e1) ! Neighbor element (global id); 0 on a boundary or mortar side
        if(e2 > 0) then
          if(this%elemToRank(e2) /= this%rankId) nRemote = nRemote+1
        endif
      enddo
    enddo

  endfunction CountRemoteSides_DomainDecomposition_t

  integer function CountRemoteMortarSides_DomainDecomposition_t(this,mortarInfo,nMortars,nSub) result(nRemote)
  !! Counts the mortar sub-sides this rank exchanges over MPI: those where exactly
  !! one of the big side and the sub-side is owned by this rank. A mortar entirely
  !! on one rank is served from memory, and one entirely on other ranks is not this
  !! rank's to carry.
  !!
  !!  Input
  !!    - mortarInfo : the mesh mortarInfo array, null on a mesh with no mortars.
  !!      Sub-side k has its small element id at row 2*k+1, in both the 2-D
  !!      (1:8,:) and 3-D (1:14,:) layouts.
  !!    - nMortars : number of mortar interfaces in the (replicated) mortar table
  !!    - nSub : sub-sides per mortar, 2 in 2-D (sub-edges) and 4 in 3-D (sub-faces)
    implicit none
    class(DomainDecomposition_t),intent(in) :: this
    integer,pointer,intent(in) :: mortarInfo(:,:)
    integer,intent(in) :: nMortars
    integer,intent(in) :: nSub
    ! Local
    integer :: m,k
    logical :: bigIsLocal,smallIsLocal

    nRemote = 0
    if(.not. associated(mortarInfo)) return

    do m = 1,nMortars
      bigIsLocal = this%elemToRank(mortarInfo(1,m)) == this%rankId
      do k = 1,nSub
        smallIsLocal = this%elemToRank(mortarInfo(2*k+1,m)) == this%rankId
        if(bigIsLocal .neqv. smallIsLocal) nRemote = nRemote+1
      enddo
    enddo

  endfunction CountRemoteMortarSides_DomainDecomposition_t

  subroutine SetElemToRank_DomainDecomposition_t(this,nElem)
    implicit none
    class(DomainDecomposition_t),intent(inout) :: this
    integer,intent(in) :: nElem
    ! Local
    integer :: iel

    this%nElem = nElem

    allocate(this%elemToRank(1:nelem))

    call DomainDecomp(nElem, &
                      this%nRanks, &
                      this%offSetElem)

    do iel = 1,nElem
      call ElemToRank(this%nRanks, &
                      this%offSetElem, &
                      iel, &
                      this%elemToRank(iel))
    enddo

  endsubroutine SetElemToRank_DomainDecomposition_t

  subroutine DomainDecomp(nElems,nDomains,offSetElem)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 4
    implicit none
    integer,intent(in) :: nElems
    integer,intent(in) :: nDomains
    integer,intent(out) :: offsetElem(0:nDomains)
    ! Local
    integer :: nLocalElems
    integer :: remainElems
    integer :: iDom

    nLocalElems = nElems/nDomains
    remainElems = nElems-nLocalElems*nDomains
    do iDom = 0,nDomains-1
      offSetElem(iDom) = iDom*nLocalElems+min(iDom,remainElems)
    enddo
    offSetElem(nDomains) = nElems

  endsubroutine DomainDecomp

  subroutine ElemToRank(nDomains,offsetElem,elemID,domain)
    ! From https://www.hopr-project.org/externals/Meshformat.pdf, Algorithm 7
    !   "Find domain containing element index"
    !
    implicit none
    integer,intent(in) :: nDomains
    integer,intent(in) :: offsetElem(0:nDomains)
    integer,intent(in) :: elemID
    integer,intent(out) :: domain
    ! Local
    integer :: maxSteps
    integer :: low,up,mid
    integer :: i

    domain = 0
    maxSteps = int(log10(real(nDomains))/log10(2.0))+1
    low = 0
    up = nDomains-1

    if(offsetElem(low) < elemID .and. elemID <= offsetElem(low+1)) then
      domain = low
    elseif(offsetElem(up) < elemID .and. elemID <= offsetElem(up+1)) then
      domain = up
    else
      do i = 1,maxSteps
        mid = (up-low)/2+low
        if(offsetElem(mid) < elemID .and. elemID <= offsetElem(mid+1)) then
          domain = mid
          return
        elseif(elemID > offsetElem(mid+1)) then
          low = mid+1
        else
          up = mid
        endif
      enddo
    endif

  endsubroutine ElemToRank

  subroutine FinalizeMPIExchangeAsync(mpiHandler)
    implicit none
    class(DomainDecomposition_t),intent(inout) :: mpiHandler
    ! Local
    integer :: ierror
    integer :: msgCount

    if(mpiHandler%mpiEnabled) then
      msgCount = mpiHandler%msgCount
      call MPI_WaitAll(msgCount, &
                       mpiHandler%requests(1:msgCount), &
                       mpiHandler%stats(1:MPI_STATUS_SIZE,1:msgCount), &
                       iError)
    endif

  endsubroutine FinalizeMPIExchangeAsync

endmodule SELF_DomainDecomposition_t
