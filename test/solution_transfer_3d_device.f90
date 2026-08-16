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

  exit_code = solution_transfer_3d_device()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function solution_transfer_3d_device() result(r)
    !! Pins the model-level 3-D AMR solution transfer - StageSolutionForTransfer -> Regrid ->
    !! ApplyTransferPlan - against the portable reference ApplyTransferPlanRange, value by value.
    !!
    !! Why this exists. On a GPU build the model overrides both steps: the solution is staged
    !! device-to-device and the plan is applied in TransferSolution_3D_gpu, so the transferred
    !! field never crosses the host link. Nothing else in the suite compares that kernel's output
    !! against the host implementation; the AMR regressions only assert invariants of the result
    !! (conservation of the Jacobian-weighted integral, entropy non-growth), which a kernel with
    !! a transposed direction or a mis-indexed child could still satisfy. This test asserts the
    !! values themselves.
    !!
    !! It is not GPU-gated. On a CPU build both sides run the same reference code and the
    !! comparison is exact, which keeps the harness itself honest; on a GPU build the left-hand
    !! side is the kernel. Agreement is asserted to the DOUBLE_PRECISION tolerance rather than
    !! bitwise, because the device compiler contracts the kernel's multiply-accumulates into FMAs.
    !!
    !! The epoch is the one transfer_plan_3d builds, chosen because a single plan then carries all
    !! three transfer kinds and prolongation depths 0..2:
    !!
    !!   COPY      - roots left untouched
    !!   RESTRICT  - a coarsened family, immediately re-refined, so depth > 0 as well
    !!   PROLONG   - a refined leaf whose child is refined again (depth 2), plus a
    !!               2:1-balancing ripple into an untouched root (depth 1)
    !!
    !! The driving field is deliberately non-polynomial and DIFFERENT IN EVERY VARIABLE. A
    !! polynomial of degree <= N is reproduced exactly by both prolongation and restriction, so it
    !! would hide an error in the operators themselves; and LinearEuler3D carries nvar = 6, of
    !! which the last two (sound speed and background density) are static background fields, so a
    !! kernel that transferred only the prognostic variables would corrupt the medium invisibly on
    !! a uniform-medium test.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_OctreeMesh_3D
    use SELF_AdaptiveMesh_3D
    use SELF_TransferPlan_3D
    use SELF_LinearEuler3D

    implicit none

    integer,parameter :: controlDegree = 4
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 1.0e-10_prec
#else
    real(prec),parameter :: tolerance = 1.0e-3_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh3D),target :: baseMesh,oldMesh,newMesh
    type(SEMHex),target :: oldGeom,newGeom
    type(OctreeMesh3D) :: forest
    type(TransferPlan3D) :: plan
    type(LinearEuler3D) :: modelobj
    integer :: bcids(1:6)
    integer :: i,j,k,e,iv,Np,nOld,nNew,nvar
    integer :: nCopy,nProlong,nRestrict,nDeep
    integer,allocatable :: flag(:),oldLeaf(:)
    real(prec),allocatable :: uOld(:,:,:,:,:),uRef(:,:,:,:,:)
    real(prec) :: x,y,z,err,scale

    r = 0
    bcids(1:6) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED, &
                  SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]

    Np = controlDegree+1

    ! 2 x 2 x 2 element structured base mesh; element ordering is x-fastest, then y, then z.
    call baseMesh%StructuredMesh(2,2,2,1,1,1,0.5_prec,0.5_prec,0.5_prec,bcids)
    call interp%Init(N=controlDegree,controlNodeType=GAUSS,M=controlDegree, &
                     targetNodeType=UNIFORM)

    ! ---- Pre-epoch state: refine root 1 so the epoch has a family available to coarsen ----
    call forest%Init(baseMesh)
    allocate(flag(1:forest%nLeaves))
    flag = OCTREE_KEEP
    flag(1) = OCTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)

    call EmitMesh(forest,baseMesh,oldMesh)
    call oldGeom%Init(interp,oldMesh%nElem)
    call oldGeom%GenerateFromMesh(oldMesh)
    nOld = forest%nLeaves
    allocate(oldLeaf(1:nOld))
    oldLeaf(1:nOld) = forest%leaf(1:nOld)

    ! Single-rank test: the plan's global new-leaf range is the rank-local range, so eFirst = 1.
    ! (The multi-rank transfer takes the uGlobal path, which is the host reference on every
    ! backend and is covered by lineareuler3d_amr_soundwave_mpi.)
    if(oldMesh%decomp%nRanks > 1) then
      print*,"FAIL: solution_transfer_3d_device is a single-rank test"
      r = 1
      return
    endif

    ! ---- The live model on the old mesh ----
    call modelobj%Init(oldMesh,oldGeom)
    modelobj%prescribed_bcs_enabled = .false.
    modelobj%tecplot_enabled = .false.
    nvar = modelobj%nvar

    ! Smooth but non-polynomial, and distinct per variable so a variable-indexing error in the
    ! kernel cannot cancel. The wavenumbers are irrational multiples of the element size, so no
    ! variable is representable exactly in the degree-N nodal basis.
    do iv = 1,nvar
      do e = 1,oldMesh%nElem
        do k = 1,Np
          do j = 1,Np
            do i = 1,Np
              x = oldGeom%x%interior(i,j,k,e,1,1)
              y = oldGeom%x%interior(i,j,k,e,1,2)
              z = oldGeom%x%interior(i,j,k,e,1,3)
              modelobj%solution%interior(i,j,k,e,iv) = &
                real(iv,prec)+sin(7.3_prec*x+real(iv,prec))* &
                cos(5.1_prec*y-real(iv,prec))*sin(3.7_prec*z+0.5_prec*real(iv,prec))
            enddo
          enddo
        enddo
      enddo
    enddo
    call modelobj%solution%UpdateDevice()

    ! The reference input: the old field exactly as staged. On a GPU build the device staging
    ! copies interior_gpu, which UpdateDevice above has just made equal to this array.
    allocate(uOld(1:Np,1:Np,1:Np,1:nOld,1:nvar))
    uOld(1:Np,1:Np,1:Np,1:nOld,1:nvar) = &
      modelobj%solution%interior(1:Np,1:Np,1:Np,1:nOld,1:nvar)

    ! ---- One adaptation epoch, built to contain all three transfer kinds ----
    allocate(flag(1:nOld))
    flag = OCTREE_KEEP
    flag(1:8) = OCTREE_COARSEN
    flag(9) = OCTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)
    call forest%RefineNode(1)
    call forest%RefineNode(forest%child(3,2))
    call forest%RebuildLeaves()
    call forest%Balance2to1()

    call BuildTransferPlan(forest,nOld,oldLeaf,plan)
    nNew = plan%nNew

    ! Guard the premise: if the plan ever stopped covering all three kinds and a multi-level
    ! descent, this test would still pass while checking much less than it claims to.
    nCopy = count(plan%sourceKind(1:nNew) == SELF_TRANSFER_COPY)
    nProlong = count(plan%sourceKind(1:nNew) == SELF_TRANSFER_PROLONG)
    nRestrict = count(plan%sourceKind(1:nNew) == SELF_TRANSFER_RESTRICT)
    nDeep = count(plan%depth(1:nNew) >= 2)
    print*,"plan: copy",nCopy,", prolong",nProlong,", restrict",nRestrict,", depth>=2",nDeep
    if(nCopy == 0 .or. nProlong == 0 .or. nRestrict == 0 .or. nDeep == 0) then
      print*,"FAIL: the epoch no longer exercises every transfer kind and a depth >= 2 descent"
      r = 1
    endif

    call EmitMesh(forest,baseMesh,newMesh)
    call newGeom%Init(interp,newMesh%nElem)
    call newGeom%GenerateFromMesh(newMesh)

    ! ---- The transfer under test ----
    call modelobj%StageSolutionForTransfer()
    call modelobj%Regrid(newMesh,newGeom)
    call modelobj%ApplyTransferPlan(plan,interp,1,nNew)
    ! On a GPU build the result is left on the device (the host mirror is stale by contract), so
    ! it must be pulled back before it can be read.
    call modelobj%solution%UpdateHost()

    ! ---- The portable reference, from the same staged input ----
    allocate(uRef(1:Np,1:Np,1:Np,1:nNew,1:nvar))
    call ApplyTransferPlanRange(plan,interp,nvar,uOld,1,nNew,uRef)

    err = 0.0_prec
    scale = 0.0_prec
    do iv = 1,nvar
      do e = 1,nNew
        do k = 1,Np
          do j = 1,Np
            do i = 1,Np
              err = max(err,abs(modelobj%solution%interior(i,j,k,e,iv)-uRef(i,j,k,e,iv)))
              scale = max(scale,abs(uRef(i,j,k,e,iv)))
            enddo
          enddo
        enddo
      enddo
    enddo
    print*,"max |device - host| =",err,", field scale =",scale
    if(err > tolerance*max(scale,1.0_prec)) then
      print*,"FAIL: the transferred solution does not match the portable reference"
      r = 1
    endif

    ! ---- The staging buffer must be reusable across epochs ----
    ! A second stage/apply on the settled mesh exercises the buffer-reuse path (no growth, so no
    ! reallocation) and must reproduce the field it was staged from, element for element: with
    ! no forest mutation between the two calls, every plan entry is a COPY of itself.
    call BuildTransferPlan(forest,nNew,forest%leaf(1:nNew),plan)
    if(any(plan%sourceKind(1:plan%nNew) /= SELF_TRANSFER_COPY)) then
      print*,"FAIL: an unmutated forest did not produce an all-COPY plan"
      r = 1
    endif
    uRef(1:Np,1:Np,1:Np,1:nNew,1:nvar) = &
      modelobj%solution%interior(1:Np,1:Np,1:Np,1:nNew,1:nvar)
    call modelobj%StageSolutionForTransfer()
    call modelobj%Regrid(newMesh,newGeom)
    call modelobj%ApplyTransferPlan(plan,interp,1,plan%nNew)
    call modelobj%solution%UpdateHost()

    err = 0.0_prec
    do iv = 1,nvar
      do e = 1,nNew
        do k = 1,Np
          do j = 1,Np
            do i = 1,Np
              err = max(err,abs(modelobj%solution%interior(i,j,k,e,iv)-uRef(i,j,k,e,iv)))
            enddo
          enddo
        enddo
      enddo
    enddo
    print*,"max |round trip - original| on the re-staged epoch =",err
    if(err > 0.0_prec) then
      print*,"FAIL: an all-COPY transfer did not reproduce the field exactly"
      r = 1
    endif

    call plan%Free()
    call modelobj%Free()
    call newGeom%Free()
    call oldGeom%Free()
    call newMesh%Free()
    call oldMesh%Free()
    call baseMesh%Free()
    call forest%Free()
    call interp%Free()
    deallocate(uOld,uRef,oldLeaf)

    if(r == 0) print*,"PASS: solution_transfer_3d_device"

  endfunction solution_transfer_3d_device

endprogram test
