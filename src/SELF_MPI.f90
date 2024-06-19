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

module SELF_MPI

  use SELF_Constants
  use SELF_Lagrange
  use SELF_SupportRoutines

  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type MPILayer
    logical :: mpiEnabled
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

    procedure :: Init => Init_MPILayer
    procedure :: Free => Free_MPILayer
    procedure :: Finalize => Finalize_MPILayer

    procedure :: GenerateDecomposition => GenerateDecomposition_MPILayer
    procedure :: SetElemToRank
    procedure :: SetMaxMsg

    procedure,public :: FinalizeMPIExchangeAsync

    generic,public :: GlobalReduce => GlobalReduce_RealScalar
    procedure,private :: GlobalReduce_RealScalar

  endtype MPILayer

contains

  subroutine Init_MPILayer(this,enableMPI)
#undef __FUNC__
#define __FUNC__ "Init_MPILayer"
    implicit none
    class(MPILayer),intent(out) :: this
    logical,intent(in) :: enableMPI
    ! Local
    integer       :: ierror
    character(50) :: msg
    character(2) :: msg2

    this%mpiComm = 0
    this%mpiPrec = prec
    this%rankId = 0
    this%nRanks = 1
    this%nElem = 0
    this%mpiEnabled = enableMPI

    if(enableMPI) then
      this%mpiComm = MPI_COMM_WORLD
      call MPI_INIT(ierror)
      call MPI_COMM_RANK(this%mpiComm,this%rankId,ierror)
      call MPI_COMM_SIZE(this%mpiComm,this%nRanks,ierror)
    endif

    if(prec == real32) then
      this%mpiPrec = MPI_FLOAT
    else
      this%mpiPrec = MPI_DOUBLE
    endif

    allocate(this%offsetElem(1:this%nRanks+1))
    !$omp target enter data map(alloc: this % offsetElem)

    write(msg,'(I5)') this%rankId
    msg = "Greetings from rank "//trim(msg)//"."
    INFO(trim(msg))

  endsubroutine Init_MPILayer

  subroutine Free_MPILayer(this)
    implicit none
    class(MPILayer),intent(inout) :: this

    if(associated(this%offSetElem)) then
      deallocate(this%offSetElem)
      !$omp target exit data map(delete: this % offsetElem)
    endif
    if(associated(this%elemToRank)) then
      deallocate(this%elemToRank)
      !$omp target exit data map(delete: this % elemToRank)
    endif

    deallocate(this%requests)
    deallocate(this%stats)

  endsubroutine Free_MPILayer

  subroutine Finalize_MPILayer(this)
#undef __FUNC__
#define __FUNC__ "Finalize_MPILayer"
    implicit none
    class(MPILayer),intent(inout) :: this
    ! Local
    integer       :: ierror
    character(30) :: msg

    if(this%mpiEnabled) then
      write(msg,'(I5)') this%rankId
      msg = "Goodbye from rank "//trim(msg)//"."
      INFO(trim(msg))
      call MPI_FINALIZE(ierror)
    endif

  endsubroutine Finalize_MPILayer

  subroutine GenerateDecomposition_MPILayer(this,nGlobalElem,maxMsg)
#undef __FUNC__
#define __FUNC__ "GenerateDecomposition_MPILayer"
    implicit none
    class(MPILayer),intent(inout) :: this
    integer,intent(in) :: nGlobalElem
    integer,intent(in) :: maxMsg
    ! Local
    character(50) :: msg
    character(5) :: msg2

    call this%setElemToRank(nGlobalElem)
    call this%SetMaxMsg(maxMsg)

    write(msg,'(I5)') this%rankId
    write(msg2,'(I5)') this%offSetElem(this%rankId+2)- &
      this%offSetElem(this%rankId+1)
    msg = "Rank "//trim(msg)//": nElem = "//trim(msg2)
    INFO(trim(msg))

  endsubroutine GenerateDecomposition_MPILayer

  subroutine SetMaxMsg(this,maxMsg)
    implicit none
    class(MPILayer),intent(inout) :: this
    integer,intent(in) :: maxMsg

    if(allocated(this%requests)) deallocate(this%requests)
    if(allocated(this%stats)) deallocate(this%stats)

    allocate(this%requests(1:maxMsg))
    allocate(this%stats(MPI_STATUS_SIZE,1:maxMsg))
    this%maxMsg = maxMsg

  endsubroutine SetMaxMsg

  subroutine SetElemToRank(this,nElem)
    implicit none
    class(MPILayer),intent(inout) :: this
    integer,intent(in) :: nElem
    ! Local
    integer :: iel

    this%nElem = nElem

    allocate(this%elemToRank(1:nelem))
    !$omp target enter data map(alloc:this % elemToRank)

    call DomainDecomp(nElem, &
                      this%nRanks, &
                      this%offSetElem)

    do iel = 1,nElem
      call ElemToRank(this%nRanks, &
                      this%offSetElem, &
                      iel, &
                      this%elemToRank(iel))
    enddo

    !$omp target update to(this % offsetElem)
    !$omp target update to(this % elemToRank)

  endsubroutine SetElemToRank

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
    class(MPILayer),intent(inout) :: mpiHandler
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

  subroutine GlobalReduce_RealScalar(mpiHandler,sendBuf,recvBuf)
    class(MPILayer),intent(in) :: mpiHandler
    real(prec),intent(in) :: sendBuf
    real(prec),intent(out) :: recvBuf
    ! Local
    integer :: iError

    if(mpiHandler%mpiEnabled) then
      call MPI_ALLREDUCE(sendBuf, &
                         recvBuf, &
                         1, &
                         mpiHandler%mpiPrec, &
                         MPI_SUM, &
                         mpiHandler%mpiComm, &
                         iError)
    else
      recvBuf = sendBuf
    endif

  endsubroutine GlobalReduce_RealScalar

endmodule SELF_MPI
