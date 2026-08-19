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

  exit_code = solution_transfer_2d_device_offset()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function solution_transfer_2d_device_offset() result(r)
    !! Calls TransferSolution_2D_gpu DIRECTLY with a window of the old field and a nonzero
    !! oldFirst0, and compares against ApplyTransferPlanWindow over the same window.
    !!
    !! Why go under the model to do this. The kernel's window offset is only reachable through the
    !! model on more than one rank: a rank's window base is nonzero because some OTHER rank owned
    !! the old elements it reads, and at one rank the hull of all new elements is always [1,nOld].
    !! The launcher, though, takes oldFirst0 as an ordinary argument, so the arithmetic can be
    !! exercised serially by handing it a window directly. That isolates an indexing error in the
    !! kernel - a dropped or mis-signed rebase, or nOld read as the global count instead of the
    !! window stride - from an error in the migration that fills the window, which
    !! amr_migrate_2d_device_apply_mpi covers end to end but cannot separate.
    !!
    !! GPU-gated, unavoidably: on a CPU build there is no launcher to call and no second
    !! implementation to compare against, so the test reports that it was skipped rather than
    !! pretending to check something. Every other new test in this change is deliberately not
    !! gated; this one is the exception because the kernel IS the subject.
    use SELF_Constants
    implicit none

#ifdef ENABLE_GPU
    r = RunOffsetCase()
#else
    r = 0
    print*,"SKIP: solution_transfer_2d_device_offset needs a GPU build (the kernel is the subject)"
#endif

  endfunction solution_transfer_2d_device_offset

#ifdef ENABLE_GPU
  integer function RunOffsetCase() result(r)
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_QuadTreeMesh_2D
    use SELF_AdaptiveMesh_2D
    use SELF_TransferPlan_2D
    use SELF_GPU
    use SELF_GPU_enums
    use SELF_GPUInterfaces
    use SELF_AMRController_2D,only:PlanWindows
    use iso_c_binding

    implicit none

    integer,parameter :: controlDegree = 3
    integer,parameter :: nvar = 3
    real(prec),parameter :: tolerance = 100.0_prec*sqrt(epsilon(1.0_prec))
    type(Lagrange),target :: interp
    type(Mesh2D),target :: baseMesh
    type(QuadTreeMesh2D) :: forest
    type(TransferPlan2D),target :: plan
    integer :: bcids(1:4)
    integer :: i,j,e,iv,Np,nOld,nNew,nLocal,nWinElem,pathStride
    integer :: eFirst,eLast,wFirst,wLast
    integer :: newOffset(1:2),winFirst(1:1),winLast(1:1)
    integer,allocatable :: flag(:),oldLeaf(:)
    real(prec),allocatable,target :: uWin(:,:,:,:),uRef(:,:,:,:),uDev(:,:,:,:)
    real(prec) :: err,fieldScale
    integer(c_size_t) :: nb
    type(c_ptr) :: uWin_gpu,uNew_gpu
    type(c_ptr) :: kind_gpu,elem_gpu,family_gpu,depth_gpu,path_gpu

    r = 0
    Np = controlDegree+1
    bcids(1:4) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]

    call baseMesh%StructuredMesh(2,2,1,1,0.5_prec,0.5_prec,bcids)
    call interp%Init(N=controlDegree,controlNodeType=GAUSS,M=controlDegree, &
                     targetNodeType=UNIFORM)
    ! No explicit upload: on a GPU build Init_Lagrange allocates and fills the device copies of
    ! the mortar matrices itself, which is what mortarR_gpu/mortarP_gpu below rely on.

    if(baseMesh%decomp%nRanks > 1) then
      print*,"FAIL: solution_transfer_2d_device_offset is a single-rank test"
      r = 1
      return
    endif

    ! ---- The transfer_plan_2d epoch: one plan carrying COPY, PROLONG at depth >= 2, RESTRICT ----
    call forest%Init(baseMesh)
    allocate(flag(1:forest%nLeaves))
    flag = QUADTREE_KEEP
    flag(1) = QUADTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)
    nOld = forest%nLeaves
    allocate(oldLeaf(1:nOld))
    oldLeaf(1:nOld) = forest%leaf(1:nOld)

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

    ! ---- A sub-range of new elements, and the hull of just its sources ----
    ! Nothing here has to satisfy the model's whole-range contract, because the launcher is called
    ! directly: it is told the new field's element stride separately from the range it fills.
    eFirst = nNew/2+1
    eLast = nNew
    newOffset(1) = eFirst-1
    newOffset(2) = eLast
    call PlanWindows(plan,1,newOffset,winFirst,winLast)
    wFirst = winFirst(1)
    wLast = winLast(1)
    nWinElem = wLast-wFirst+1
    nLocal = eLast-eFirst+1

    print*,"new [",eFirst,",",eLast,"] of",nNew,"; window [",wFirst,",",wLast,"] of",nOld

    ! Guard the premise: without all three, the test cannot tell a correct rebase from a dropped
    ! one or from nOld being read as the global old count.
    if(wFirst <= 1) then
      print*,"FAIL: the window starts at 1, so oldFirst0 is zero and a dropped rebase would be ", &
        "invisible; this test cannot detect it"
      r = 1
    endif
    if(nWinElem >= nOld) then
      print*,"FAIL: the window spans the whole old field, so a wrong nOld stride would be ", &
        "invisible; this test cannot detect it"
      r = 1
    endif
    if(nWinElem < 2) then
      print*,"FAIL: the window is a single element, so the stride is untested"
      r = 1
    endif
    if(r /= 0) return

    ! ---- A window of old-field values ----
    ! One decimal place per index, so the map from (element,node,variable) to value is a
    ! BIJECTION: a wrong source element cannot coincide with a legitimate value. Divided by a
    ! power of two so every value is exactly representable and the difference below reflects
    ! only the operators, not the encoding.
    allocate(uWin(1:Np,1:Np,1:nWinElem,1:nvar))
    do iv = 1,nvar
      do e = 1,nWinElem
        do j = 1,Np
          do i = 1,Np
            uWin(i,j,e,iv) = real(1000000*iv+1000*(e+wFirst-1)+10*i+j,prec)/1024.0_prec
          enddo
        enddo
      enddo
    enddo

    ! ---- The portable reference over the same window ----
    allocate(uRef(1:Np,1:Np,1:nLocal,1:nvar))
    call ApplyTransferPlanWindow(plan,interp,nvar,uWin,wFirst,wLast,eFirst,eLast,uRef)

    ! ---- The kernel, given the window and its base ----
    pathStride = size(plan%path,1)
    nb = int(plan%nNew,c_size_t)*c_int
    call gpuCheck(hipMalloc(kind_gpu,nb))
    call gpuCheck(hipMalloc(elem_gpu,nb))
    call gpuCheck(hipMalloc(depth_gpu,nb))
    call gpuCheck(hipMalloc(family_gpu,4_c_size_t*nb))
    call gpuCheck(hipMalloc(path_gpu,int(pathStride,c_size_t)*nb))
    call gpuCheck(hipMemcpy(kind_gpu,c_loc(plan%sourceKind),nb,hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(elem_gpu,c_loc(plan%sourceElem),nb,hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(depth_gpu,c_loc(plan%depth),nb,hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(family_gpu,c_loc(plan%family),4_c_size_t*nb,hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(path_gpu,c_loc(plan%path), &
                            int(pathStride,c_size_t)*nb,hipMemcpyHostToDevice))

    nb = int(Np*Np,c_size_t)*int(nWinElem,c_size_t)*int(nvar,c_size_t)*prec
    call gpuCheck(hipMalloc(uWin_gpu,nb))
    call gpuCheck(hipMemcpy(uWin_gpu,c_loc(uWin),nb,hipMemcpyHostToDevice))

    allocate(uDev(1:Np,1:Np,1:nLocal,1:nvar))
    uDev = 0.0_prec
    nb = int(Np*Np,c_size_t)*int(nLocal,c_size_t)*int(nvar,c_size_t)*prec
    call gpuCheck(hipMalloc(uNew_gpu,nb))
    call gpuCheck(hipMemcpy(uNew_gpu,c_loc(uDev),nb,hipMemcpyHostToDevice))

    ! nOld is the WINDOW's element stride and oldFirst0 its 0-based global base; nNew here is the
    ! new field's own stride, which is the sub-range length because that is how uDev is shaped.
    call TransferSolution_2D_gpu(uWin_gpu,uNew_gpu,kind_gpu,elem_gpu,family_gpu, &
                                 depth_gpu,path_gpu,interp%mortarR_gpu,interp%mortarP_gpu, &
                                 pathStride,eFirst-1,wFirst-1,interp%N,nvar, &
                                 nWinElem,nLocal,nLocal)
    call gpuCheck(hipMemcpy(c_loc(uDev),uNew_gpu,nb,hipMemcpyDeviceToHost))

    err = 0.0_prec
    fieldScale = 0.0_prec
    do iv = 1,nvar
      do e = 1,nLocal
        do j = 1,Np
          do i = 1,Np
            err = max(err,abs(uDev(i,j,e,iv)-uRef(i,j,e,iv)))
            fieldScale = max(fieldScale,abs(uRef(i,j,e,iv)))
          enddo
        enddo
      enddo
    enddo
    print*,"max |kernel - windowed reference| =",err,", field scale =",fieldScale
    ! Tolerance, not bitwise: the device compiler contracts the kernel's multiply-accumulates into
    ! FMAs. A dropped or mis-signed oldFirst0 would read a different element entirely, and the
    ! analytic values differ between neighbouring elements by ~0.1 - six orders of magnitude above
    ! this bar - so the tolerance cannot hide an indexing error.
    if(err > tolerance*max(fieldScale,1.0_prec)) then
      print*,"FAIL: the kernel's windowed output does not match the portable windowed reference"
      r = 1
    endif

    call gpuCheck(hipFree(uWin_gpu))
    call gpuCheck(hipFree(uNew_gpu))
    call gpuCheck(hipFree(kind_gpu))
    call gpuCheck(hipFree(elem_gpu))
    call gpuCheck(hipFree(depth_gpu))
    call gpuCheck(hipFree(family_gpu))
    call gpuCheck(hipFree(path_gpu))
    deallocate(uWin,uRef,uDev,oldLeaf)
    call plan%Free()
    call forest%Free()
    call interp%Free()
    call baseMesh%Free()

    if(r == 0) print*,"PASS: solution_transfer_2d_device_offset"

  endfunction RunOffsetCase
#endif

endprogram test
