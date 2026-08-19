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

  use mpi
  implicit none
  integer :: exit_code,ierror

  call mpi_init(ierror)
  exit_code = amr_migrate_2d_device_apply_mpi()
  call mpi_finalize(ierror)
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function amr_migrate_2d_device_apply_mpi() result(r)
    !! Drives the MULTI-RANK 2-D solution migration and windowed apply through the model, which is
    !! the only configuration in which the device path is reachable, and compares the result
    !! against an independent whole-field reference.
    !!
    !! What this covers that amr_migrate_2d_window_mpi does not. That test calls ExchangeOldWindow
    !! directly, so it always migrates in HOST memory whatever the backend. On a GPU build the
    !! model's MigrateOldWindow instead assembles the window in DEVICE memory: the rank's own run
    !! by a device-to-device copy, the peers' runs received by MPI_Irecv straight into device
    !! memory, and the plan applied to it by TransferSolution_2D_gpu with the oldFirst0 offset.
    !! Nothing else exercises that. In particular the AMR soundwave MPI regressions do not: under a
    !! symmetric refinement the contiguous space-filling-curve repartition leaves every rank's new
    !! range sourced from its own old range, so not a single byte moves and only the local copy
    !! runs. Hence the deliberately ASYMMETRIC refinement below, inherited from
    !! amr_migrate_2d_window_mpi.
    !!
    !! Assertions, and the bar each is held to:
    !!
    !!   (a) some rank received old elements from a peer - otherwise the run proves nothing about
    !!       the messages, whatever else passes.
    !!   (b) some rank's window base is NONTRIVIAL: wFirst > 1 and wFirst /= its own first old
    !!       element. Without this the offset is untested even at four ranks, because a window that
    !!       begins exactly where the rank's own old range begins is indistinguishable from
    !!       indexing the local buffer from 1 - which is what a dropped oldFirst0 does.
    !!   (c) the migrated window itself, BITWISE against the analytic field. Migration is pure data
    !!       movement - a device-to-device copy, MPI byte transfers, and the download used to read
    !!       it here - so no arithmetic touches these values and exactness is available. A mismatch
    !!       names the element that arrived wrong.
    !!   (d) the transferred solution, against the portable whole-field apply. Tolerance on a GPU
    !!       build, because the kernel's multiply-accumulates are contracted into FMAs and its
    !!       descent prolongs only onto the child on the recorded path; BITWISE otherwise, where
    !!       both sides are the same code and the stronger bar is free.
    !!   (e) the byte counters agree with what moved: a rank that received elements must have
    !!       counted bytes. The related case of a rank owning NO new elements that still has to
    !!       serve peers is REPORTED, not required - whether the partition produces one depends on
    !!       the mesh and rank count, and asserting it here would make the test's premise depend on
    !!       a coincidence. The arithmetic for it is covered without MPI by
    !!       transfer_plan_{2,3}d_window's case table, which constructs empty new ranges directly.
    !!
    !! Runs at any rank count >= 1; it is registered at 2 and at 4, because two ranks cannot
    !! produce a window that straddles more than one peer.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_QuadTreeMesh_2D
    use SELF_AdaptiveMesh_2D
    use SELF_TransferPlan_2D
    use SELF_LinearEuler2D
    use SELF_AMRController_2D,only:PlanWindows,InitForestFromDecomposedMesh
    use iso_fortran_env,only:int64
    use mpi

    implicit none

    integer,parameter :: controlDegree = 3
    integer,parameter :: nRefine = 5 !! leaves 1..nRefine refine: deliberately asymmetric
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 1.0e-10_prec
#else
    real(prec),parameter :: tolerance = 1.0e-3_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh2D),target :: baseMesh
    type(Mesh2D),pointer :: newMesh
    type(SEMQuad),target :: oldGeom
    type(SEMQuad),pointer :: newGeom
    type(QuadTreeMesh2D) :: forest
    type(TransferPlan2D) :: plan
    type(LinearEuler2D) :: modelobj
    integer :: bcids(1:4)
    integer :: i,j,e,iv,nOld,nNew,Np,nvar,ierror
    integer :: nR,rankId,myFirst,myLast,nLocalOld,nLocal
    integer :: eFirst,eLast,wFirst,wLast,nWinElem
    integer :: offsetFlag,offsetFlagMax
    integer(int64) :: nBytesRecv,nBytesSent,nElemRemote,nRemoteMax
    integer,allocatable :: flag(:),oldLeaf(:),winFirst(:),winLast(:)
    real(prec),allocatable :: uWin(:,:,:,:),uAll(:,:,:,:),uRef(:,:,:,:)
    real(prec) :: err,fieldScale

    r = 0
    Np = controlDegree+1
    bcids(1:4) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]

    ! 8 x 8 base mesh: 64 elements, enough that a 2- or 4-rank split leaves several elements per
    ! rank and the refined band sits inside the first rank's range only - which is what makes the
    ! repartition asymmetric and forces real traffic.
    call baseMesh%StructuredMesh(8,8,1,1,0.125_prec,0.125_prec,bcids)
    call interp%Init(N=controlDegree,controlNodeType=GAUSS,M=controlDegree, &
                     targetNodeType=UNIFORM)

    nR = baseMesh%decomp%nRanks
    rankId = baseMesh%decomp%rankId
    myFirst = baseMesh%decomp%offsetElem(rankId+1)+1
    myLast = baseMesh%decomp%offsetElem(rankId+2)
    nLocalOld = myLast-myFirst+1

    if(nLocalOld /= baseMesh%nElem) then
      print*,"FAIL: rank",rankId," local old count",nLocalOld, &
        " disagrees with the mesh's element count",baseMesh%nElem
      r = 1
    endif

    call oldGeom%Init(interp,baseMesh%nElem)
    call oldGeom%GenerateFromMesh(baseMesh)

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

    ! ---- The live model on the old mesh, carrying the analytic field ----
    call modelobj%Init(baseMesh,oldGeom)
    modelobj%prescribed_bcs_enabled = .false.
    modelobj%tecplot_enabled = .false.
    nvar = modelobj%nvar

    do iv = 1,nvar
      do e = 1,nLocalOld
        do j = 1,Np
          do i = 1,Np
            modelobj%solution%interior(i,j,e,iv) = OldValue(e+myFirst-1,i,j,iv)
          enddo
        enddo
      enddo
    enddo
    call modelobj%solution%UpdateDevice()

    ! ---- One deliberately asymmetric epoch ----
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
    nWinElem = max(wLast-wFirst+1,0)

    print*,"rank",rankId," old [",myFirst,",",myLast,"] new [",eFirst,",",eLast, &
      "] window [",wFirst,",",wLast,"]"

    ! ---- Migrate, through the model: device-resident on a GPU build ----
    nBytesRecv = 0
    nBytesSent = 0
    nElemRemote = 0
    call modelobj%MigrateOldWindow(winFirst,winLast,wFirst,wLast, &
                                   nBytesRecv,nBytesSent,nElemRemote)

    ! ---- (a) elements really did cross ranks ----
    if(nR > 1) then
      call mpi_allreduce(nElemRemote,nRemoteMax,1,MPI_INTEGER8,MPI_MAX, &
                         baseMesh%decomp%mpiComm,ierror)
      if(nRemoteMax <= 0) then
        print*,"FAIL: no rank received any old element from a peer, so the device receive"// &
          " path was never exercised"
        r = 1
      endif
    endif

    ! ---- (b) the window base is nontrivial somewhere ----
    ! A window starting at 1, or starting exactly at the rank's own first old element, cannot
    ! distinguish a correct oldFirst0 from a dropped one. Require that some rank has neither.
    offsetFlag = 0
    if(nWinElem > 0 .and. wFirst > 1 .and. wFirst /= myFirst .and. nWinElem < nOld) then
      offsetFlag = 1
    endif
    if(nR > 1) then
      call mpi_allreduce(offsetFlag,offsetFlagMax,1,MPI_INTEGER,MPI_MAX, &
                         baseMesh%decomp%mpiComm,ierror)
    else
      offsetFlagMax = offsetFlag
    endif
    if(nR > 1 .and. offsetFlagMax == 0) then
      print*,"FAIL: no rank has a window base past 1 and away from its own first old element,"// &
        " so a dropped oldFirst0 offset would be invisible to this run"
      r = 1
    endif

    ! ---- (c) the migrated window carries its own analytic values, BITWISE ----
    if(nWinElem > 0) then
      allocate(uWin(1:Np,1:Np,1:nWinElem,1:nvar))
      call modelobj%DownloadOldWindow(wFirst,wLast,uWin)
      do iv = 1,nvar
        do e = 1,nWinElem
          do j = 1,Np
            do i = 1,Np
              if(uWin(i,j,e,iv) /= OldValue(e+wFirst-1,i,j,iv)) then
                print*,"FAIL: rank",rankId," window element",e+wFirst-1," variable",iv, &
                  " node",i,j," holds",uWin(i,j,e,iv)," expected",OldValue(e+wFirst-1,i,j,iv)
                r = 1
              endif
            enddo
          enddo
        enddo
      enddo
      deallocate(uWin)
    endif

    ! ---- (e) counters are consistent with what moved ----
    if(nElemRemote > 0 .and. nBytesRecv <= 0) then
      print*,"FAIL: rank",rankId," received elements but counted no bytes"
      r = 1
    endif
    if(nWinElem == 0 .and. nR > 1) then
      ! A rank owning no new elements still has to serve peers whose windows reach into its old
      ! range. Report rather than require: whether it happens depends on the partition.
      print*,"rank",rankId," owns no new elements; bytes sent =",nBytesSent
    endif
    print*,"rank",rankId," remote elements",nElemRemote," bytes recv/sent",nBytesRecv,nBytesSent

    ! ---- Regrid and apply ----
    allocate(newGeom)
    call newGeom%Init(interp,newMesh%nElem)
    call newGeom%GenerateFromMesh(newMesh)
    call modelobj%Regrid(newMesh,newGeom)

    if(eLast >= eFirst) then
      call modelobj%ApplyTransferPlan(plan,interp,eFirst,eLast)
      ! On a GPU build the result is left on the device by contract.
      call modelobj%solution%UpdateHost()

      ! ---- (d) against the portable whole-field apply, from an independently built old field ----
      ! uAll is the analytic field over ALL old elements, so the reference never passes through the
      ! migration and cannot inherit a routing error from it.
      nLocal = eLast-eFirst+1
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
      allocate(uRef(1:Np,1:Np,1:nLocal,1:nvar))
      call ApplyTransferPlanRange(plan,interp,nvar,uAll,eFirst,eLast,uRef)

      err = 0.0_prec
      fieldScale = 0.0_prec
      do iv = 1,nvar
        do e = 1,nLocal
          do j = 1,Np
            do i = 1,Np
              err = max(err,abs(modelobj%solution%interior(i,j,e,iv)-uRef(i,j,e,iv)))
              fieldScale = max(fieldScale,abs(uRef(i,j,e,iv)))
            enddo
          enddo
        enddo
      enddo
      print*,"rank",rankId," max |windowed apply - whole-field apply| =",err, &
        ", field scale =",fieldScale
#ifdef ENABLE_GPU
      if(err > tolerance*max(fieldScale,1.0_prec)) then
        print*,"FAIL: rank",rankId," the device windowed apply differs from the portable"// &
          " whole-field apply beyond round-off"
        r = 1
      endif
#else
      if(err > 0.0_prec) then
        print*,"FAIL: rank",rankId," the windowed apply differs from the whole-field apply,"// &
          " and on a CPU build these are the same code so they must agree exactly"
        r = 1
      endif
#endif
      deallocate(uAll,uRef)
    endif

    if(r == 0 .and. rankId == 0) print*,"AMR MIGRATE 2D DEVICE APPLY CHECKS PASSED"

    deallocate(winFirst,winLast,oldLeaf)
    call plan%Free()
    call modelobj%Free()
    call newGeom%Free()
    deallocate(newGeom)
    call oldGeom%Free()
    call newMesh%Free()
    deallocate(newMesh)
    call forest%Free()
    call interp%Free()
    call baseMesh%Free()

  endfunction amr_migrate_2d_device_apply_mpi

  pure real(prec) function OldValue(eGlobal,i,j,iv)
    !! Analytic old-field value: distinct for every (global element, node, variable) and exactly
    !! representable, so the window comparison can be exact and a mismatch names the element that
    !! arrived. Scaled down from the sibling test's magnitudes so that the transfer operators,
    !! which form weighted sums of these, stay comfortably inside the mantissa.
    use SELF_Constants
    implicit none
    integer,intent(in) :: eGlobal,i,j,iv

    OldValue = real(1000*iv+100*eGlobal+10*i+j,prec)/1024.0_prec

  endfunction OldValue

endprogram test
