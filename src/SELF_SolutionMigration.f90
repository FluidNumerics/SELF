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

module SELF_SolutionMigration
!! Point-to-point migration of the pre-regrid solution into a rank's old-element WINDOW: the
!! memory-space-agnostic core of AMR Stage-5 v2, shared by every backend and both dimensions.
!!
!! An adapting epoch repartitions the mesh, so the rank that will own a new element and the rank
!! that owned its old source need not be the same. Both partitions are contiguous ranges of the
!! same space-filling-curve leaf order and the transfer plan is rank-replicated, so each rank can
!! compute the contiguous window of old elements its own new range references (PlanWindows, which
!! needs the plan and therefore stays with the AMR controllers) and then receive exactly that
!! window - rather than allgathering the whole old field, which is what v1 did.
!!
!! Two properties make the schedule cheap, and both are worth understanding before editing this
!! module:
!!
!!   1. Routing needs NO communication. Send and receive runs are the same intersection
!!      (OwnedRun) evaluated from the same replicated offsetElem and window tables on both ends
!!      of every pair, so the two schedules match by construction: no handshake, no count
!!      exchange, no collective. Reaching for MPI_Alltoallv here would pay for a negotiation that
!!      replication has already answered.
!!
!!   2. NOTHING needs packing. A run of elements at a fixed variable is contiguous in both the
!!      source field and the window, so each message reads and writes the real storage directly.
!!      That costs one message per (peer,variable) instead of one per peer, and it removes a pack
!!      buffer, an unpack loop and two whole-field copies.
!!
!! Everything here works on FLAT buffers, indexed as (perElem, nElem, nvar) in Fortran order with
!! perElem = (N+1)**2 in 2-D and (N+1)**3 in 3-D. That is what lets one copy of the schedule
!! serve both dimensions and both backends: the GPU backend hands these routines Fortran
!! descriptors over DEVICE allocations and MPI receives straight into device memory (the idiom of
!! SideExchangeStart in SELF_MappedScalar_3D.f90; GPU-aware MPI is already an unconditional
!! requirement of any multi-rank GPU run in SELF), while the portable backend hands them host
!! allocatables.
!!
!! The buffer arguments are assumed-size, so a caller with a rank-4 or rank-5 window passes the
!! whole array and relies on sequence association. Every such actual argument in SELF is
!! contiguous (explicit-shape dummies, allocatables, or `pointer,contiguous`), so no copy-in
!! temporary can appear - which matters because MPI is handed the addresses.
!!
!! Runs once per adapting epoch, between time steps - never inside the time-stepping loop.
  use iso_fortran_env,only:int64
  use SELF_Constants
  use SELF_DomainDecomposition
  use mpi

  implicit none

contains

  subroutine OwnedRun(offsetElem,r,first,last,a,b)
    !! The run of elements that rank r-1 owns and that also lies in [first,last]: a..b, empty when
    !! b < a. offsetElem is a decomposition's contiguous ownership table, so rank r-1 owns
    !! offsetElem(r)+1 .. offsetElem(r+1).
    !!
    !! Both halves of the migration schedule are this one intersection, evaluated from replicated
    !! tables: my receive run from peer r-1 is (my window) n (r-1's old range), and my send run to
    !! peer r-1 is (r-1's window) n (my old range). Factored out so that the sender and the
    !! receiver cannot drift apart, and so a serial test can exercise the arithmetic the exchange
    !! actually uses rather than a copy of it.
    implicit none
    integer,intent(in) :: offsetElem(:)
    integer,intent(in) :: r !! 1-based table index; the rank is r-1
    integer,intent(in) :: first
    integer,intent(in) :: last
    integer,intent(out) :: a
    integer,intent(out) :: b

    a = max(first,offsetElem(r)+1)
    b = min(last,offsetElem(r+1))

  endsubroutine OwnedRun

  subroutine PostOldWindowExchange(decomp,perElem,nvar,nLocalOld,uLocal,winFirst,winLast, &
                                   wFirst,wLast,uWin,requests,msgCount, &
                                   nBytesRecv,nBytesSent,nElemRemote)
    !! Post every message of this rank's migration schedule and return without waiting: the
    !! receives that fill the remote runs of its window, then the sends that serve the runs of its
    !! old range that other ranks' windows need. Receives are posted before sends, as everywhere
    !! else in SELF.
    !!
    !! Neither buffer is touched other than by MPI, so this is also the device-resident form: a
    !! GPU backend passes descriptors over device memory and the received runs land there
    !! directly. The caller supplies the request array (sized for the worst case of every rank
    !! being a peer) and pairs this with FinishOldWindowExchange.
    !!
    !! The whole window is written every epoch, not just the elements the plan reads: the local
    !! part plus the per-peer runs tile [wFirst,wLast] exactly once, because the old partition
    !! tiles the old element list. So a reused window buffer carries nothing stale from a previous
    !! epoch, and the verification path may compare all of it.
    implicit none
    type(DomainDecomposition),intent(in) :: decomp
    integer,intent(in) :: perElem !! reals per element per variable: (N+1)**2 in 2-D, (N+1)**3 in 3-D
    integer,intent(in) :: nvar
    integer,intent(in) :: nLocalOld !! elements this rank owned before the epoch (uLocal's stride)
    real(prec),intent(in) :: uLocal(*) !! (perElem,nLocalOld,nvar), first element = this rank's first
    integer,intent(in) :: winFirst(1:decomp%nRanks)
    integer,intent(in) :: winLast(1:decomp%nRanks)
    integer,intent(in) :: wFirst !! this rank's window, normalized (wFirst > wLast if empty)
    integer,intent(in) :: wLast
    real(prec),intent(inout) :: uWin(*) !! (perElem,nWinElem,nvar), first element = global old wFirst
    !! intent(inout), not out: the receives write it through MPI rather than through the dummy,
    !! and an intent(out) dummy would license a compiler to treat it as undefined on entry.
    integer,intent(inout) :: requests(:) !! at least 2*nvar*decomp%nRanks entries
    integer,intent(inout) :: msgCount !! 0 on entry; the number of posted messages on exit
    integer(int64),intent(inout) :: nBytesRecv
    integer(int64),intent(inout) :: nBytesSent
    integer(int64),intent(inout) :: nElemRemote
    ! Local
    integer :: r,iv,a,b,cnt,ierror
    integer :: myFirst,nWinElem,nbyte
    ! int64: off spans the whole buffer in reals (perElem*nWinElem*nvar), which a
    ! high-order 3-D window can push past the 32-bit range even though each message's cnt
    ! stays small. Fortran's own subscript arithmetic used to carry this at pointer width;
    ! computing it explicitly is what makes the kind matter.
    integer(int64) :: off

    myFirst = decomp%offsetElem(decomp%rankId+1)+1
    nWinElem = max(wLast-wFirst+1,0)
    ! storage_size, not the kind value prec: the two coincide for the kinds SELF uses, but only
    ! the intrinsic actually reports a width, and these counters are quoted as bytes.
    nbyte = storage_size(1.0_prec)/8

    ! Receives first: the runs of my window that other ranks own.
    do r = 1,decomp%nRanks
      if(r-1 == decomp%rankId) cycle
      call OwnedRun(decomp%offsetElem,r,wFirst,wLast,a,b)
      if(b < a) cycle
      cnt = perElem*(b-a+1)
      nBytesRecv = nBytesRecv+int(cnt,int64)*nvar*nbyte
      nElemRemote = nElemRemote+int(b-a+1,int64)
      do iv = 1,nvar
        off = int(perElem,int64)*(int(a-wFirst,int64)+int(nWinElem,int64)*(iv-1))
        msgCount = msgCount+1
        call MPI_IRECV(uWin(off+1),cnt,decomp%mpiPrec, &
                       r-1,iv,decomp%mpiComm,requests(msgCount),ierror)
      enddo
    enddo

    ! Sends: the runs of my old range that other ranks' windows need. A rank that owns no new
    ! elements still reaches this loop, because its peers may need old elements it owns.
    do r = 1,decomp%nRanks
      if(r-1 == decomp%rankId) cycle
      call OwnedRun(decomp%offsetElem,decomp%rankId+1,winFirst(r),winLast(r),a,b)
      if(b < a) cycle
      cnt = perElem*(b-a+1)
      nBytesSent = nBytesSent+int(cnt,int64)*nvar*nbyte
      do iv = 1,nvar
        off = int(perElem,int64)*(int(a-myFirst,int64)+int(nLocalOld,int64)*(iv-1))
        msgCount = msgCount+1
        call MPI_ISEND(uLocal(off+1),cnt,decomp%mpiPrec, &
                       r-1,iv,decomp%mpiComm,requests(msgCount),ierror)
      enddo
    enddo

  endsubroutine PostOldWindowExchange

  subroutine FinishOldWindowExchange(requests,msgCount)
    !! Wait on the messages PostOldWindowExchange posted. A no-op when nothing was posted, which
    !! is the common case for a symmetric refinement: under a contiguous space-filling-curve
    !! repartition each rank's new range is sourced from its own old range, so the window is
    !! entirely local and not a single byte moves.
    implicit none
    integer,intent(inout) :: requests(:)
    integer,intent(in) :: msgCount
    ! Local
    integer :: ierror
    integer,allocatable :: stats(:,:)

    if(msgCount > 0) then
      allocate(stats(MPI_STATUS_SIZE,1:msgCount))
      call MPI_WAITALL(msgCount,requests(1:msgCount), &
                       stats(1:MPI_STATUS_SIZE,1:msgCount),ierror)
      deallocate(stats)
    endif

  endsubroutine FinishOldWindowExchange

  subroutine ExchangeOldWindowFlat(decomp,perElem,nvar,nLocalOld,uLocal,winFirst,winLast, &
                                   wFirst,wLast,uWin,nBytesRecv,nBytesSent,nElemRemote)
    !! The complete host-memory migration: post the messages, copy the part of the window this
    !! rank already owns while they are in flight, then wait. This is the portable backend's whole
    !! implementation; the GPU backend calls PostOldWindowExchange / FinishOldWindowExchange
    !! itself so that it can serve the local part with a device-to-device copy instead.
    implicit none
    type(DomainDecomposition),intent(in) :: decomp
    integer,intent(in) :: perElem
    integer,intent(in) :: nvar
    integer,intent(in) :: nLocalOld
    real(prec),intent(in) :: uLocal(*)
    integer,intent(in) :: winFirst(1:decomp%nRanks)
    integer,intent(in) :: winLast(1:decomp%nRanks)
    integer,intent(in) :: wFirst
    integer,intent(in) :: wLast
    real(prec),intent(inout) :: uWin(*)
    integer(int64),intent(inout) :: nBytesRecv
    integer(int64),intent(inout) :: nBytesSent
    integer(int64),intent(inout) :: nElemRemote
    ! Local
    integer :: a,b,e,iv,p,msgCount
    integer :: myFirst,nWinElem
    integer(int64) :: dst,src !! int64 for the reason given in PostOldWindowExchange
    integer,allocatable :: requests(:)

    myFirst = decomp%offsetElem(decomp%rankId+1)+1
    nWinElem = max(wLast-wFirst+1,0)

    allocate(requests(1:2*nvar*decomp%nRanks))
    msgCount = 0
    call PostOldWindowExchange(decomp,perElem,nvar,nLocalOld,uLocal,winFirst,winLast, &
                               wFirst,wLast,uWin,requests,msgCount, &
                               nBytesRecv,nBytesSent,nElemRemote)

    ! The part of my window I already own, copied while the messages are in flight. Element order
    ! within a variable is the storage order of both buffers, so this is a pure byte move and the
    ! migrated window is bit-identical to the v1 allgathered field.
    call OwnedRun(decomp%offsetElem,decomp%rankId+1,wFirst,wLast,a,b)
    do iv = 1,nvar
      do e = a,b
        dst = int(perElem,int64)*(int(e-wFirst,int64)+int(nWinElem,int64)*(iv-1))
        src = int(perElem,int64)*(int(e-myFirst,int64)+int(nLocalOld,int64)*(iv-1))
        do p = 1,perElem
          uWin(dst+p) = uLocal(src+p)
        enddo
      enddo
    enddo

    call FinishOldWindowExchange(requests,msgCount)
    deallocate(requests)

  endsubroutine ExchangeOldWindowFlat

endmodule SELF_SolutionMigration
