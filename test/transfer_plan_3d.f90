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

  exit_code = transfer_plan_3d()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function transfer_plan_3d() result(r)
    !! Validates the AMR adaptation-epoch transfer plan (old-leaf -> new-leaf solution mapping)
    !! end to end against emitted meshes and real geometry:
    !!
    !!   1. Classification: one epoch that simultaneously coarsens a family (which is then
    !!      re-refined, exercising RESTRICT with depth > 0), refines a leaf whose child is
    !!      refined again (PROLONG with depth 2, the Balance2to1-ripple shape), lets 2:1
    !!      balancing ripple into an untouched root (PROLONG depth 1), and keeps five roots
    !!      untouched (COPY). All three transfer kinds and depths 0..2 appear, with the
    !!      expected multiplicities.
    !!   2. Exactness: for a trilinear field (degree <= N in each direction) every transfer path
    !!      - copy, multi-level prolongation, and restriction of children that carry the parent
    !!      polynomial exactly - reproduces the analytic field at the new mesh's nodes to
    !!      roundoff.
    !!   3. Conservation: for a smooth non-polynomial field, the Jacobian-weighted global
    !!      integral int u dV is identical on the old and new emitted meshes to roundoff
    !!      (prolongation is exact interpolation; restriction is the conservative L2
    !!      projection; the base mesh is affine so child Jacobians are exact eighths).
    !!   4. Reversibility: refine-everything then coarsen-everything across two epochs returns
    !!      arbitrary nodal data to its original values to roundoff (sum_k P_k R_k = I).
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_OctreeMesh_3D
    use SELF_AdaptiveMesh_3D
    use SELF_TransferPlan_3D

    implicit none

    integer,parameter :: controlDegree = 4
    integer,parameter :: nvar = 2
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 1.0e-10_prec
#else
    real(prec),parameter :: tolerance = 1.0e-3_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh3D),target :: baseMesh,oldMesh,newMesh
    type(SEMHex),target :: oldGeom,newGeom
    type(OctreeMesh3D) :: forest
    type(TransferPlan3D) :: plan,planA,planB
    integer :: bcids(1:6)
    integer :: i,j,k,e,c,nOld,nNew
    integer :: nCopy,nProlong1,nProlong2,nRestrict
    integer,allocatable :: flag(:),oldLeaf(:)
    real(prec),allocatable :: uOld(:,:,:,:,:),uNew(:,:,:,:,:)
    real(prec),allocatable :: w0(:,:,:,:,:),wf(:,:,:,:,:),wb(:,:,:,:,:)
    real(prec) :: x,y,z,expected,err,intOld,intNew

    r = 0
    bcids(1:6) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED, &
                  SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]

    ! 2 x 2 x 2 element structured base mesh; element ordering is x-fastest, then y, then z:
    ! elements 1,2 are the south row of the bottom layer, 3,4 its north row, 5-8 the top layer.
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

    ! Old mesh/geometry and the snapshot that opens the epoch. Old element index = position in
    ! the leaf list: 1-8 are root 1's children, 9-15 are roots 2,3,...,8.
    call EmitMesh(forest,baseMesh,oldMesh)
    call oldGeom%Init(interp,oldMesh%nElem)
    call oldGeom%GenerateFromMesh(oldMesh)
    nOld = forest%nLeaves
    allocate(oldLeaf(1:nOld))
    oldLeaf(1:nOld) = forest%leaf(1:nOld)

    ! ---- One adaptation epoch ----
    ! Coarsen root 1's family and refine root 2 ...
    allocate(flag(1:nOld))
    flag = OCTREE_KEEP
    flag(1:8) = OCTREE_COARSEN
    flag(9) = OCTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)
    ! ... re-refine the just-coarsened root 1 (RESTRICT + prolong, the coarsen-then-balance
    ! shape) and refine root 2's child 3 (its upper-x/upper-y/lower-z octant) again (PROLONG
    ! depth 2) ...
    call forest%RefineNode(1)
    call forest%RefineNode(forest%child(3,2))
    call forest%RebuildLeaves()
    ! ... and 2:1-balance, which must ripple into root 4 (north of root 2's level-2 leaves;
    ! the other faces of root 2's refined child are domain boundaries or sibling faces).
    call forest%Balance2to1()
    if(forest%MaxLevelJump() > 1) then
      print*,"FAIL: epoch forest is not balanced"
      r = 1
    endif

    call BuildTransferPlan(forest,nOld,oldLeaf,plan)

    ! ---- 1. Classification ----
    ! Expected: roots 3,5,6,7,8 survive (5 COPY); root 1's eight new children restrict its
    ! coarsened family (8 RESTRICT, depth 1, family = old elements 1-8); root 2's seven
    ! remaining children and root 4's eight balance-ripple children prolong one level (15
    ! PROLONG depth 1); root 2's child-3 grandchildren prolong two levels (8 PROLONG depth 2).
    ! 36 new leaves in total.
    nNew = plan%nNew
    if(nNew /= forest%nLeaves .or. nNew /= 36) then
      print*,"FAIL: unexpected new leaf count",nNew
      r = 1
    endif
    nCopy = 0
    nProlong1 = 0
    nProlong2 = 0
    nRestrict = 0
    do e = 1,nNew
      select case(plan%sourceKind(e))
      case(SELF_TRANSFER_COPY)
        nCopy = nCopy+1
      case(SELF_TRANSFER_PROLONG)
        if(plan%depth(e) == 1) nProlong1 = nProlong1+1
        if(plan%depth(e) == 2) nProlong2 = nProlong2+1
      case(SELF_TRANSFER_RESTRICT)
        nRestrict = nRestrict+1
        if(plan%depth(e) /= 1) then
          print*,"FAIL: restricted leaf has depth",plan%depth(e)
          r = 1
        endif
        do c = 1,8
          if(plan%family(c,e) /= c) then
            print*,"FAIL: restrict family mismatch",plan%family(1:8,e)
            r = 1
          endif
        enddo
      endselect
    enddo
    print*,"copy, prolong(d1), prolong(d2), restrict :",nCopy,nProlong1,nProlong2,nRestrict
    if(nCopy /= 5 .or. nProlong1 /= 15 .or. nProlong2 /= 8 .or. nRestrict /= 8) then
      print*,"FAIL: transfer plan classification"
      r = 1
    endif
    if(plan%maxDepth /= 2) then
      print*,"FAIL: plan maxDepth",plan%maxDepth
      r = 1
    endif

    ! ---- New mesh/geometry for value and conservation checks ----
    call EmitMesh(forest,baseMesh,newMesh)
    call newGeom%Init(interp,newMesh%nElem)
    call newGeom%GenerateFromMesh(newMesh)

    ! ---- 2. Exactness on a trilinear field ----
    allocate(uOld(1:controlDegree+1,1:controlDegree+1,1:controlDegree+1,1:nOld,1:nvar))
    allocate(uNew(1:controlDegree+1,1:controlDegree+1,1:controlDegree+1,1:nNew,1:nvar))
    do e = 1,nOld
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            x = oldGeom%x%interior(i,j,k,e,1,1)
            y = oldGeom%x%interior(i,j,k,e,1,2)
            z = oldGeom%x%interior(i,j,k,e,1,3)
            uOld(i,j,k,e,1) = 2.0_prec+3.0_prec*x-5.0_prec*y+4.0_prec*z+0.25_prec*x*y*z
            uOld(i,j,k,e,2) = -1.0_prec+0.5_prec*x+y-2.0_prec*z+x*y-y*z
          enddo
        enddo
      enddo
    enddo
    call ApplyTransferPlan(plan,interp,nvar,uOld,uNew)
    err = 0.0_prec
    do e = 1,nNew
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            x = newGeom%x%interior(i,j,k,e,1,1)
            y = newGeom%x%interior(i,j,k,e,1,2)
            z = newGeom%x%interior(i,j,k,e,1,3)
            expected = 2.0_prec+3.0_prec*x-5.0_prec*y+4.0_prec*z+0.25_prec*x*y*z
            err = max(err,abs(uNew(i,j,k,e,1)-expected))
            expected = -1.0_prec+0.5_prec*x+y-2.0_prec*z+x*y-y*z
            err = max(err,abs(uNew(i,j,k,e,2)-expected))
          enddo
        enddo
      enddo
    enddo
    print*,"trilinear-field transfer max error :",err
    if(err > tolerance) then
      print*,"FAIL: transfer is not exact for a trilinear field"
      r = 1
    endif

    ! ---- 3. Conservation of int u dV for a non-polynomial field ----
    do e = 1,nOld
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            x = oldGeom%x%interior(i,j,k,e,1,1)
            y = oldGeom%x%interior(i,j,k,e,1,2)
            z = oldGeom%x%interior(i,j,k,e,1,3)
            uOld(i,j,k,e,1) = sin(1.7_prec*x)*cos(1.3_prec*y)*cos(1.1_prec*z)+0.5_prec*x*y*z
            uOld(i,j,k,e,2) = cos(2.1_prec*x)*sin(0.9_prec*y)*sin(1.2_prec*z)-x+y-z
          enddo
        enddo
      enddo
    enddo
    call ApplyTransferPlan(plan,interp,nvar,uOld,uNew)
    intOld = 0.0_prec
    do e = 1,nOld
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            intOld = intOld+uOld(i,j,k,e,1)*oldGeom%J%interior(i,j,k,e,1)* &
                     interp%qWeights(i)*interp%qWeights(j)*interp%qWeights(k)
          enddo
        enddo
      enddo
    enddo
    intNew = 0.0_prec
    do e = 1,nNew
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            intNew = intNew+uNew(i,j,k,e,1)*newGeom%J%interior(i,j,k,e,1)* &
                     interp%qWeights(i)*interp%qWeights(j)*interp%qWeights(k)
          enddo
        enddo
      enddo
    enddo
    print*,"integral old, new :",intOld,intNew
    if(abs(intOld-intNew) > tolerance*max(abs(intOld),1.0_prec)) then
      print*,"FAIL: transfer is not conservative"
      r = 1
    endif
    deallocate(uOld,uNew)
    call plan%Free()
    call forest%Free()

    ! ---- 4. Reversibility across two epochs (refine everything, coarsen everything) ----
    call forest%Init(baseMesh)
    nOld = forest%nLeaves
    deallocate(oldLeaf)
    allocate(oldLeaf(1:nOld))
    oldLeaf(1:nOld) = forest%leaf(1:nOld)
    allocate(flag(1:nOld))
    flag = OCTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)
    call BuildTransferPlan(forest,nOld,oldLeaf,planA)

    allocate(w0(1:controlDegree+1,1:controlDegree+1,1:controlDegree+1,1:nOld,1:nvar))
    allocate(wf(1:controlDegree+1,1:controlDegree+1,1:controlDegree+1,1:planA%nNew,1:nvar))
    allocate(wb(1:controlDegree+1,1:controlDegree+1,1:controlDegree+1,1:nOld,1:nvar))
    do e = 1,nOld ! arbitrary (non-polynomial) nodal data
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            w0(i,j,k,e,1) = sin(real(17*i+5*j+3*k,prec))+real(e,prec)
            w0(i,j,k,e,2) = cos(real(3*i-7*j+11*k,prec))-2.0_prec*real(e,prec)
          enddo
        enddo
      enddo
    enddo
    call ApplyTransferPlan(planA,interp,nvar,w0,wf)

    deallocate(oldLeaf)
    allocate(oldLeaf(1:forest%nLeaves))
    oldLeaf(1:forest%nLeaves) = forest%leaf(1:forest%nLeaves)
    allocate(flag(1:forest%nLeaves))
    flag = OCTREE_COARSEN
    call forest%AdaptFromFlags(flag)
    call BuildTransferPlan(forest,size(flag),oldLeaf,planB)
    deallocate(flag)
    call ApplyTransferPlan(planB,interp,nvar,wf,wb)

    err = 0.0_prec
    do e = 1,nOld
      err = max(err,maxval(abs(wb(:,:,:,e,:)-w0(:,:,:,e,:))))
    enddo
    print*,"refine-coarsen round-trip max error :",err
    if(err > tolerance) then
      print*,"FAIL: refine-then-coarsen is not the identity"
      r = 1
    endif

    if(r == 0) print*,"TRANSFER PLAN CHECKS PASSED"

    deallocate(w0,wf,wb,oldLeaf)
    call planA%Free()
    call planB%Free()
    call forest%Free()
    call oldGeom%Free()
    call newGeom%Free()
    call interp%Free()
    call oldMesh%Free()
    call newMesh%Free()
    call baseMesh%Free()

  endfunction transfer_plan_3d

endprogram test
