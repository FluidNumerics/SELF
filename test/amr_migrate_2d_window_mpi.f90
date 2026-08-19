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

  exit_code = amr_migrate_2d_window_mpi()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function amr_migrate_2d_window_mpi() result(r)
    !! MPI regression for the point-to-point AMR solution migration (Stage-5 v2): it exercises
    !! ExchangeOldWindow itself, with old elements that genuinely change rank.
    !!
    !! Why this test exists and the AMR soundwave MPI runs are not enough. Under a symmetric
    !! refinement the contiguous space-filling-curve repartition keeps each rank's new range
    !! sourced from its own old range, so the window is entirely local, no message is posted at
    !! all, and the local copy alone produces the right answer. Such a run would still pass
    !! SELF_AMR_MIGRATE_VERIFY, so it proves nothing about the send/receive schedule. Here the
    !! adaptation is deliberately ASYMMETRIC - only the first few leaves of the forest refine -
    !! so the new partition boundaries land inside what used to be another rank's old range and
    !! elements must actually move. The test asserts that they did.
    !!
    !! The old field is analytic and encodes each element's GLOBAL old index in its values, so a
    !! misrouted element is not merely "different", it names the element that arrived instead.
    !!
    !! Assertions:
    !!  (a) at least one rank received at least one old element from a peer (otherwise the test
    !!      has stopped testing what it is for, and the migration counters are checked as well);
    !!  (b) every element of the received window carries exactly the analytic values of that
    !!      global old element - bit for bit, on every rank;
    !!  (c) the window covers everything the rank's new range references (this is what
    !!      ApplyTransferPlanWindow's own guard would catch, asserted here directly);
    !!  (d) applying the plan from the migrated window reproduces, bit for bit, the result of
    !!      applying it from the whole global old field.
    !!
    !! Runs on any rank count >= 1; with one rank it degenerates to the local-copy path, which is
    !! still a valid (if weaker) check, and (a) is then skipped.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_QuadTreeMesh_2D
    use SELF_AdaptiveMesh_2D
    use SELF_TransferPlan_2D
    use SELF_AMRController_2D,only:PlanWindows,ExchangeOldWindow, &
                                    InitForestFromDecomposedMesh
    use iso_fortran_env,only:int64
    use mpi

    implicit none

    integer,parameter :: controlDegree = 3
    integer,parameter :: nvar = 2
    integer,parameter :: nRefine = 5 !! leaves 1..nRefine refine: deliberately asymmetric
    type(Lagrange),target :: interp
    type(Mesh2D),target :: baseMesh
    type(Mesh2D),pointer :: newMesh
    type(QuadTreeMesh2D) :: forest
    type(TransferPlan2D) :: plan
    integer :: bcids(1:4)
    integer :: i,j,e,c,iv,nOld,nNew,Np,ierror
    integer :: nR,rankId,myFirst,myLast,nLocalOld
    integer :: eFirst,eLast,wFirst,wLast,src
    integer(int64) :: nBytesRecv,nBytesSent,nElemRemote,nRemoteMax
    integer,allocatable :: flag(:),oldLeaf(:),winFirst(:),winLast(:)
    real(prec),allocatable :: uLocal(:,:,:,:),uWin(:,:,:,:)
    real(prec),allocatable :: uAll(:,:,:,:),uRef(:,:,:,:),uNew(:,:,:,:)

    r = 0
    Np = controlDegree+1
    bcids(1:4) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]

    ! 8 x 8 base mesh: 64 elements, enough that a 2- or 4-rank split leaves several elements per
    ! rank and the refined band sits inside the first rank's range only.
    call baseMesh%StructuredMesh(8,8,1,1,0.125_prec,0.125_prec,bcids)
    call interp%Init(N=controlDegree,controlNodeType=GAUSS,M=controlDegree, &
                     targetNodeType=UNIFORM)

    nR = baseMesh%decomp%nRanks
    rankId = baseMesh%decomp%rankId
    myFirst = baseMesh%decomp%offsetElem(rankId+1)+1
    myLast = baseMesh%decomp%offsetElem(rankId+2)
    nLocalOld = myLast-myFirst+1

    ! ---- Rank-replicated forest over the decomposed base mesh ----
    if(nR > 1) then
      call InitForestFromDecomposedMesh(forest,baseMesh)
    else
      call forest%Init(baseMesh)
    endif

    nOld = forest%nLeaves
    if(nOld /= 64) then
      print*,"FAIL: unexpected old leaf count",nOld
      r = 1
    endif
    allocate(oldLeaf(1:nOld))
    oldLeaf(1:nOld) = forest%leaf(1:nOld)

    ! ---- One deliberately asymmetric epoch: refine only the first nRefine leaves ----
    allocate(flag(1:nOld))
    flag = QUADTREE_KEEP
    flag(1:nRefine) = QUADTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)
    call forest%Balance2to1()
    call BuildTransferPlan(forest,nOld,oldLeaf,plan)
    nNew = plan%nNew
    if(rankId == 0) print*,"nOld, nNew :",nOld,nNew

    ! ---- The emitted mesh supplies the new partition, exactly as the controller uses it ----
    allocate(newMesh)
    call EmitMesh(forest,baseMesh,newMesh)
    eFirst = newMesh%decomp%offsetElem(newMesh%decomp%rankId+1)+1
    eLast = newMesh%decomp%offsetElem(newMesh%decomp%rankId+2)

    allocate(winFirst(1:nR),winLast(1:nR))
    call PlanWindows(plan,nR,newMesh%decomp%offsetElem,winFirst,winLast)
    wFirst = winFirst(rankId+1)
    wLast = winLast(rankId+1)
    if(eLast < eFirst) then
      wFirst = 1
      wLast = 0
    endif

    ! ---- Analytic old field: the value names the global old element it belongs to ----
    allocate(uLocal(1:Np,1:Np,1:nLocalOld,1:nvar))
    do iv = 1,nvar
      do e = 1,nLocalOld
        do j = 1,Np
          do i = 1,Np
            uLocal(i,j,e,iv) = OldValue(e+myFirst-1,i,j,iv)
          enddo
        enddo
      enddo
    enddo

    ! ---- Migrate ----
    nBytesRecv = 0
    nBytesSent = 0
    nElemRemote = 0
    allocate(uWin(1:Np,1:Np,wFirst:wLast,1:nvar))
    call ExchangeOldWindow(baseMesh%decomp,Np,nvar,nLocalOld,uLocal,winFirst,winLast, &
                           wFirst,wLast,uWin,nBytesRecv,nBytesSent,nElemRemote)

    ! ---- (a) elements really did cross ranks ----
    if(nR > 1) then
      call mpi_allreduce(nElemRemote,nRemoteMax,1,MPI_INTEGER8,MPI_MAX, &
                         baseMesh%decomp%mpiComm,ierror)
      if(nRemoteMax <= 0) then
        print*,"FAIL: no rank received any old element from a peer, so the point-to-point"// &
          " schedule was never exercised"
        r = 1
      endif
      if(nElemRemote > 0 .and. nBytesRecv <= 0) then
        print*,"FAIL: rank",rankId," received elements but counted no bytes"
        r = 1
      endif
      print*,"rank",rankId," window",wFirst,wLast," remote elements",nElemRemote, &
        " bytes recv/sent",nBytesRecv,nBytesSent
    endif

    ! ---- (b) every window element carries its own analytic values ----
    do iv = 1,nvar
      do e = wFirst,wLast
        do j = 1,Np
          do i = 1,Np
            if(uWin(i,j,e,iv) /= OldValue(e,i,j,iv)) then
              print*,"FAIL: rank",rankId," window element",e," variable",iv, &
                " node",i,j," holds",uWin(i,j,e,iv)," expected",OldValue(e,i,j,iv)
              r = 1
            endif
          enddo
        enddo
      enddo
    enddo

    ! ---- (c) the window covers every old element the rank's new range references ----
    do e = eFirst,eLast
      if(plan%sourceKind(e) == SELF_TRANSFER_RESTRICT) then
        do c = 1,4
          src = plan%family(c,e)
          if(src < wFirst .or. src > wLast) then
            print*,"FAIL: rank",rankId," window",wFirst,wLast," misses family source",src
            r = 1
          endif
        enddo
      else
        src = plan%sourceElem(e)
        if(src < wFirst .or. src > wLast) then
          print*,"FAIL: rank",rankId," window",wFirst,wLast," misses source",src
          r = 1
        endif
      endif
    enddo

    ! ---- (d) transferring from the window matches transferring from the global field ----
    if(eLast >= eFirst) then
      allocate(uAll(1:Np,1:Np,1:nOld,1:nvar))
      do iv = 1,nvar
        do e = 1,nOld
          do j = 1,Np
            do i = 1,Np
              uAll(i,j,e,iv) = OldValue(e,i,j,iv)
            enddo
          enddo
        enddo
      enddo
      allocate(uRef(1:Np,1:Np,1:eLast-eFirst+1,1:nvar))
      allocate(uNew(1:Np,1:Np,1:eLast-eFirst+1,1:nvar))
      call ApplyTransferPlanRange(plan,interp,nvar,uAll,eFirst,eLast,uRef)
      call ApplyTransferPlanWindow(plan,interp,nvar,uWin,wFirst,wLast,eFirst,eLast,uNew)
      if(any(uNew /= uRef)) then
        print*,"FAIL: rank",rankId," transfer from the migrated window differs from the"// &
          " whole-field transfer; max|d| =",maxval(abs(uNew-uRef))
        r = 1
      endif
      deallocate(uAll,uRef,uNew)
    endif

    if(r == 0 .and. rankId == 0) print*,"AMR MIGRATE 2D WINDOW CHECKS PASSED"

    deallocate(uLocal,uWin,winFirst,winLast,oldLeaf)
    call plan%Free()
    call newMesh%Free()
    deallocate(newMesh)
    call forest%Free()
    call interp%Free()
    call baseMesh%Free()

  endfunction amr_migrate_2d_window_mpi

  pure real(prec) function OldValue(eGlobal,i,j,iv)
    !! Analytic old-field value: distinct for every (global element, node, variable) and exactly
    !! representable, so a comparison can be exact and a mismatch names the element that arrived.
    use SELF_Constants
    implicit none
    integer,intent(in) :: eGlobal,i,j,iv

    OldValue = real(1000000*iv+1000*eGlobal+10*i+j,prec)

  endfunction OldValue

endprogram test
