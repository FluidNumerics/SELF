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
    integer :: mpiComm
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

    procedure,public :: FinalizeMPIExchangeAsync

  endtype DomainDecomposition_t

contains

  subroutine Init_DomainDecomposition_t(this)
    implicit none
    class(DomainDecomposition_t),intent(inout) :: this
    ! Local
    integer       :: ierror

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

    this%initialized = .true.

    this%initialized = .true.

  endsubroutine Init_DomainDecomposition_t

  subroutine Free_DomainDecomposition_t(this)
    implicit none
    class(DomainDecomposition_t),intent(inout) :: this
    ! Local
    integer :: ierror

    if(associated(this%offSetElem)) then
      deallocate(this%offSetElem)
    endif
    if(associated(this%elemToRank)) then
      deallocate(this%elemToRank)
    endif

    if(allocated(this%requests)) deallocate(this%requests)
    if(allocated(this%stats)) deallocate(this%stats)

    !if(this%mpiEnabled) then
    print*,__FILE__," : Rank ",this%rankId+1,"/",this%nRanks," checking out."
    call MPI_FINALIZE(ierror)
    !endif

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
