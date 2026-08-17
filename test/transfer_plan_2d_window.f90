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

program test

  implicit none
  integer :: exit_code

  exit_code = transfer_plan_2d_window()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function transfer_plan_2d_window() result(r)
    !! Unit test for the routing arithmetic behind the point-to-point AMR solution migration
    !! (Stage-5 v2): PlanWindows, the peer-overlap decomposition ExchangeOldWindow performs, and
    !! the windowed apply ApplyTransferPlanWindow. Serial by construction - the whole point of
    !! the design is that routing is a LOCAL computation on the rank-replicated plan, so every
    !! rank's decisions can be reproduced and checked in one process, for partitions that a
    !! 2-rank MPI test could never produce.
    !!
    !! One adaptation epoch of a 2 x 2 base forest produces a plan containing all three
    !! transfer kinds (COPY, PROLONG at depths 1 and 2, RESTRICT of a complete four-child
    !! family), so the checks below cover the family case, whose four children may straddle a
    !! partition boundary.
    !!
    !! For each of a table of (old partition, new partition) pairs, and for every rank in it:
    !!
    !!  (a) Coverage: every old element the rank's new range references lies inside the window
    !!      PlanWindows returned. This is what correctness requires of the window.
    !!  (b) Tightness: both window ends are themselves referenced, so the window is the exact
    !!      hull and the migration moves no more than it must. This is what the performance
    !!      claim requires.
    !!  (c) Empty ranges: a rank owning no new elements gets the empty sentinel
    !!      (winFirst > winLast).
    !!  (d) Overlap decomposition: the per-peer contiguous runs ExchangeOldWindow derives from
    !!      the old offsetElem table tile the window exactly once - no gap, no double cover, no
    !!      off-by-one - and every referenced element is claimed by exactly one owner.
    !!  (e) Windowed apply: applying the plan to the rank's new range from a window slice
    !!      reproduces, BIT FOR BIT, what applying it from the whole global old field produces.
    !!      Exact equality is the bar, not a tolerance: it is the same operands fed to the same
    !!      operators in the same order, which is the guarantee the migration rests on.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_QuadTreeMesh_2D
    use SELF_AdaptiveMesh_2D
    use SELF_TransferPlan_2D
    use SELF_AMRController_2D,only:PlanWindows,OwnedRun

    implicit none

    integer,parameter :: controlDegree = 3
    integer,parameter :: nvar = 2
    integer,parameter :: nCases = 5
    integer,parameter :: maxRanks = 5
    type(Lagrange),target :: interp
    type(Mesh2D),target :: baseMesh
    type(QuadTreeMesh2D) :: forest
    type(TransferPlan2D) :: plan
    integer :: bcids(1:4)
    integer :: i,j,e,c,iv,rank,icase,nOld,nNew,Np
    integer :: nRanks,eFirst,eLast,wFirst,wLast,a,b,src
    integer :: nPeer,nCovered,maxPeers,nAllRemote,nEmptyWindow
    integer,allocatable :: flag(:),oldLeaf(:)
    integer,allocatable :: winFirst(:),winLast(:)
    integer,allocatable :: claimed(:)
    integer :: oldOff(0:maxRanks),newOff(0:maxRanks)
    real(prec),allocatable :: uOld(:,:,:,:),uRef(:,:,:,:),uWin(:,:,:,:),uLoc(:,:,:,:)
    real(prec),allocatable :: uAsm(:,:,:,:)

    r = 0
    maxPeers = 0
    nAllRemote = 0
    nEmptyWindow = 0
    Np = controlDegree+1
    bcids(1:4) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]

    call baseMesh%StructuredMesh(2,2,1,1,0.5_prec,0.5_prec,bcids)
    call interp%Init(N=controlDegree,controlNodeType=GAUSS,M=controlDegree, &
                     targetNodeType=UNIFORM)

    ! ---- Pre-epoch state: refine root 1, so the epoch has a family available to coarsen ----
    call forest%Init(baseMesh)
    allocate(flag(1:forest%nLeaves))
    flag = QUADTREE_KEEP
    flag(1) = QUADTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)

    nOld = forest%nLeaves
    allocate(oldLeaf(1:nOld))
    oldLeaf(1:nOld) = forest%leaf(1:nOld)

    ! ---- One adaptation epoch: coarsen root 1's family, refine root 2, re-refine root 1, ----
    ! ---- refine one grandchild, then 2:1-balance. Yields COPY, PROLONG d1/d2, RESTRICT.  ----
    allocate(flag(1:nOld))
    flag = QUADTREE_KEEP
    flag(1:4) = QUADTREE_COARSEN
    flag(5) = QUADTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)
    call forest%RefineNode(1)
    call forest%RefineNode(forest%child(3,2))
    call forest%RebuildLeaves()
    call forest%Balance2to1()

    call BuildTransferPlan(forest,nOld,oldLeaf,plan)
    nNew = plan%nNew
    print*,"nOld, nNew :",nOld,nNew

    ! ---- Synthetic old field: whole numbers, so every transfer result is reproducible ----
    ! exactly and a bitwise comparison is meaningful.
    allocate(uOld(1:Np,1:Np,1:nOld,1:nvar))
    do iv = 1,nvar
      do e = 1,nOld
        do j = 1,Np
          do i = 1,Np
            uOld(i,j,e,iv) = real(100000*iv+1000*e+10*i+j,prec)
          enddo
        enddo
      enddo
    enddo

    allocate(winFirst(1:maxRanks),winLast(1:maxRanks))
    allocate(claimed(1:nOld))

    do icase = 1,nCases

      ! Partition tables, in the offsetElem convention: rank q (0-based) owns elements
      ! off(q)+1 .. off(q+1). Deliberately including cases DomainDecomp would never generate,
      ! to exercise the routing rather than the balance policy.
      select case(icase)
      case(1) ! one rank: the whole field, the degenerate case
        nRanks = 1
        oldOff(0:1) = [0,nOld]
        newOff(0:1) = [0,nNew]
      case(2) ! four near-equal ranks, as DomainDecomp would split them
        nRanks = 4
        oldOff(0:4) = [0,nOld/4,nOld/2,(3*nOld)/4,nOld]
        newOff(0:4) = [0,nNew/4,nNew/2,(3*nNew)/4,nNew]
      case(3) ! heavy skew: the new partition bears no resemblance to the old one, so rank 0's
        ! window spans every old owner - the multi-peer routing a 2-rank run cannot produce
        nRanks = 5
        oldOff(0:5) = [0,nOld/5,(2*nOld)/5,(3*nOld)/5,(4*nOld)/5,nOld]
        newOff(0:5) = [0,nNew-3,nNew-2,nNew-1,nNew,nNew]
      case(4) ! empty new ranges (ranks 0, 2 and 4 own nothing) and an empty old range
        nRanks = 5
        oldOff(0:5) = [0,0,nOld/2,nOld/2,nOld,nOld]
        newOff(0:5) = [0,0,nNew/2,nNew/2,nNew,nNew]
      case(5) ! reversed loading: rank 0 takes almost everything new while owning almost
        ! nothing old, so its window is nearly all remote
        nRanks = 2
        oldOff(0:2) = [0,1,nOld]
        newOff(0:2) = [0,nNew-1,nNew]
      endselect

      call PlanWindows(plan,nRanks,newOff(0:nRanks),winFirst(1:nRanks),winLast(1:nRanks))

      do rank = 0,nRanks-1
        eFirst = newOff(rank)+1
        eLast = newOff(rank+1)
        wFirst = winFirst(rank+1)
        wLast = winLast(rank+1)

        ! ---- (c) empty new range -> empty window sentinel ----
        if(eLast < eFirst) then
          if(wFirst <= wLast) then
            print*,"FAIL: case",icase," rank",rank, &
              " owns no new elements but got a non-empty window",wFirst,wLast
            r = 1
          endif
          ! The exchange normalizes an empty window to (1,0). Every peer overlap must then come
          ! out empty, so such a rank posts no receives at all - it only serves its peers.
          nPeer = 0
          do i = 1,nRanks
            call OwnedRun(oldOff(0:nRanks),i,1,0,a,b)
            if(b >= a) nPeer = nPeer+1
          enddo
          if(nPeer /= 0) then
            print*,"FAIL: case",icase," rank",rank, &
              " would post",nPeer," receives for an empty window"
            r = 1
          endif
          nEmptyWindow = nEmptyWindow+1
          cycle
        endif

        if(wFirst < 1 .or. wLast > nOld .or. wLast < wFirst) then
          print*,"FAIL: case",icase," rank",rank," window out of range",wFirst,wLast
          r = 1
          cycle
        endif

        ! ---- (a) coverage and (b) tightness ----
        claimed(1:nOld) = 0
        do e = eFirst,eLast
          if(plan%sourceKind(e) == SELF_TRANSFER_RESTRICT) then
            do c = 1,4
              src = plan%family(c,e)
              if(src < wFirst .or. src > wLast) then
                print*,"FAIL: case",icase," rank",rank," family source",src, &
                  " outside window",wFirst,wLast
                r = 1
              else
                claimed(src) = 1
              endif
            enddo
          else
            src = plan%sourceElem(e)
            if(src < wFirst .or. src > wLast) then
              print*,"FAIL: case",icase," rank",rank," source",src, &
                " outside window",wFirst,wLast
              r = 1
            else
              claimed(src) = 1
            endif
          endif
        enddo
        if(claimed(wFirst) == 0 .or. claimed(wLast) == 0) then
          print*,"FAIL: case",icase," rank",rank," window",wFirst,wLast, &
            " is not tight (an end is never referenced)"
          r = 1
        endif

        ! ---- (d) the peer-overlap decomposition tiles the window exactly once ----
        ! Driving OwnedRun, the routine ExchangeOldWindow itself uses for both its receive runs
        ! and its send runs, so this checks that code rather than a copy of its formula.
        nPeer = 0
        nCovered = 0
        do i = 1,nRanks
          call OwnedRun(oldOff(0:nRanks),i,wFirst,wLast,a,b)
          if(b < a) cycle
          nCovered = nCovered+(b-a+1)
          if(i-1 /= rank) nPeer = nPeer+1
        enddo
        if(nCovered /= wLast-wFirst+1) then
          print*,"FAIL: case",icase," rank",rank," peer runs cover",nCovered, &
            " of",wLast-wFirst+1," window elements"
          r = 1
        endif
        maxPeers = max(maxPeers,nPeer)
        call OwnedRun(oldOff(0:nRanks),rank+1,wFirst,wLast,a,b)
        if(b < a) nAllRemote = nAllRemote+1

        ! ---- (e) windowed apply is bit-identical to the whole-field apply ----
        allocate(uRef(1:Np,1:Np,1:eLast-eFirst+1,1:nvar))
        allocate(uWin(1:Np,1:Np,1:eLast-eFirst+1,1:nvar))
        call ApplyTransferPlanRange(plan,interp,nvar,uOld,eFirst,eLast,uRef)
        call ApplyTransferPlanWindow(plan,interp,nvar,uOld(:,:,wFirst:wLast,:), &
                                     wFirst,wLast,eFirst,eLast,uWin)
        if(any(uWin /= uRef)) then
          print*,"FAIL: case",icase," rank",rank, &
            " windowed apply differs from the whole-field apply; max|d| =", &
            maxval(abs(uWin-uRef))
          r = 1
        endif
        deallocate(uRef,uWin)

      enddo
    enddo

    ! ---- The table above must actually have reached the interesting configurations ----
    if(maxPeers < 3) then
      print*,"FAIL: no case produced a window spanning 3 or more peers (maxPeers =",maxPeers,")"
      r = 1
    endif
    if(nAllRemote == 0) then
      print*,"FAIL: no case produced a rank whose window is entirely remote"
      r = 1
    endif
    if(nEmptyWindow == 0) then
      print*,"FAIL: no case produced a rank owning no new elements"
      r = 1
    endif
    print*,"max peers, all-remote windows, empty windows :",maxPeers,nAllRemote,nEmptyWindow

    ! ---- A window assembled the way ExchangeOldWindow assembles it transfers identically ----
    ! Same routing as case 3 rank 0, but the window is filled run by run from per-rank local
    ! slices rather than sliced out of a global array, which is what the exchange does.
    nRanks = 5
    oldOff(0:5) = [0,nOld/5,(2*nOld)/5,(3*nOld)/5,(4*nOld)/5,nOld]
    newOff(0:5) = [0,nNew-3,nNew-2,nNew-1,nNew,nNew]
    call PlanWindows(plan,nRanks,newOff(0:nRanks),winFirst(1:nRanks),winLast(1:nRanks))
    wFirst = winFirst(1)
    wLast = winLast(1)
    eFirst = 1
    eLast = newOff(1)
    allocate(uWin(1:Np,1:Np,wFirst:wLast,1:nvar))
    uWin = 0.0_prec
    do i = 1,nRanks
      a = max(wFirst,oldOff(i-1)+1)
      b = min(wLast,oldOff(i))
      if(b < a) cycle
      ! rank i-1's own storage: its old elements, indexed from 1
      allocate(uLoc(1:Np,1:Np,1:oldOff(i)-oldOff(i-1),1:nvar))
      uLoc(:,:,:,:) = uOld(:,:,oldOff(i-1)+1:oldOff(i),:)
      uWin(:,:,a:b,:) = uLoc(:,:,a-oldOff(i-1):b-oldOff(i-1),:)
      deallocate(uLoc)
    enddo
    allocate(uRef(1:Np,1:Np,1:eLast-eFirst+1,1:nvar))
    allocate(uAsm(1:Np,1:Np,1:eLast-eFirst+1,1:nvar))
    call ApplyTransferPlanRange(plan,interp,nvar,uOld,eFirst,eLast,uRef)
    call ApplyTransferPlanWindow(plan,interp,nvar,uWin,wFirst,wLast,eFirst,eLast,uAsm)
    if(any(uAsm /= uRef)) then
      print*,"FAIL: run-by-run assembled window differs from the whole-field apply"
      r = 1
    endif
    deallocate(uRef,uAsm,uWin)

    if(r == 0) print*,"TRANSFER PLAN 2D WINDOW CHECKS PASSED"

    call plan%Free()
    call forest%Free()
    deallocate(uOld,winFirst,winLast,claimed,oldLeaf)
    call interp%Free()
    call baseMesh%Free()

  endfunction transfer_plan_2d_window

endprogram test
