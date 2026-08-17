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

module SELF_TransferPlan_3D
!! Old-leaf -> new-leaf solution transfer plan for one adaptation epoch of the 3-D oct-forest:
!! the driver layer of AMR Stage 3 that connects the element-local transfer operators
!! (SELF_SolutionTransfer_3D) to the forest mutation (SELF_OctreeMesh_3D).
!!
!! One adaptation epoch mutates the forest between two leaf configurations:
!!
!!     nOld = forest%nLeaves ; oldLeaf = forest%leaf(1:nOld)   ! snapshot BEFORE mutating
!!     call forest%AdaptFromFlags(flag)                        ! at most one call, then
!!     call forest%Balance2to1()                               ! any number of further
!!                                                             ! refinements (RefineNode too)
!!     call BuildTransferPlan(forest,nOld,oldLeaf,plan)        ! AFTER the last mutation
!!
!! Solution element index = leaf-list position (the element ordering EmitMesh produces), so the
!! plan records, for every new leaf, where its data comes from in the old element ordering:
!!
!!   SELF_TRANSFER_COPY     : the leaf survived unchanged; copy old element sourceElem.
!!   SELF_TRANSFER_PROLONG  : the leaf descends from old leaf sourceElem; interpolate the old
!!                            element's degree-N polynomial down the octree path (one step per
!!                            level, so depth > 1 handles a fresh child that 2:1 balancing
!!                            refined again in the same epoch). Exact, no loss.
!!   SELF_TRANSFER_RESTRICT : the leaf is (an ancestor of) a coarsened family; L2-project the
!!                            eight old children family(1:8) onto their parent, then prolong
!!                            down depth >= 0 further steps (depth > 0 occurs when a
!!                            just-coarsened parent is immediately re-refined by 2:1
!!                            balancing). Conservative.
!!
!! The reconstruction is possible after the fact because forest node ids are stable: refinement
!! only appends nodes and coarsening only detaches children, whose level/parent/octant entries
!! persist. Each new leaf therefore ascends its parent chain until it meets either an old leaf
!! or a node holding a complete eight-child old-leaf family; anything else means the snapshot
!! does not describe the epoch that produced the forest, and the builder stops with an error.
!!
!! ApplyTransferPlan executes the plan on nodal data u(1:N+1,1:N+1,1:N+1,1:nElem,1:nVar) (the
!! layout of MappedScalar3D %interior; units are those of the transferred fields) and inherits
!! the Stage-3 operator identities: prolongation is exact polynomial interpolation and
!! restriction is the conservative L2 projection, so the Jacobian-weighted cell integral of
!! every variable is conserved and a refine-then-coarsen round trip is the identity to
!! roundoff. Transfer runs once per adaptation epoch, between time steps; it is not a per-step
!! hot path.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_OctreeMesh_3D
  use SELF_SolutionTransfer_3D

  implicit none

  integer,parameter :: SELF_TRANSFER_COPY = 0
  integer,parameter :: SELF_TRANSFER_PROLONG = 1
  integer,parameter :: SELF_TRANSFER_RESTRICT = 2

  type :: TransferPlan3D
    integer :: nOld = 0 !! old (pre-epoch) element count
    integer :: nNew = 0 !! new element count (= forest%nLeaves at build time)
    integer :: maxDepth = 0 !! largest prolongation depth in the plan
    integer,allocatable :: sourceKind(:) !! (nNew) SELF_TRANSFER_COPY / PROLONG / RESTRICT
    integer,allocatable :: sourceElem(:) !! (nNew) old element index (COPY/PROLONG; 0 otherwise)
    integer,allocatable :: family(:,:) !! (8,nNew) old element index per child octant (RESTRICT)
    integer,allocatable :: depth(:) !! (nNew) number of prolongation steps below the source
    integer,allocatable :: path(:,:) !! (>=maxDepth,nNew) octant taken at each step, top-down

  contains
    procedure,public :: Free => Free_TransferPlan3D

  endtype TransferPlan3D

contains

  subroutine Free_TransferPlan3D(this)
    implicit none
    class(TransferPlan3D),intent(inout) :: this

    if(allocated(this%sourceKind)) deallocate(this%sourceKind)
    if(allocated(this%sourceElem)) deallocate(this%sourceElem)
    if(allocated(this%family)) deallocate(this%family)
    if(allocated(this%depth)) deallocate(this%depth)
    if(allocated(this%path)) deallocate(this%path)
    this%nOld = 0
    this%nNew = 0
    this%maxDepth = 0

  endsubroutine Free_TransferPlan3D

  subroutine BuildTransferPlan(forest,nOld,oldLeaf,plan)
    !! Build the old->new transfer plan for the adaptation epoch that took the forest from the
    !! leaf configuration (nOld, oldLeaf) - a snapshot of (forest%nLeaves, forest%leaf) taken
    !! before mutating - to its current leaf configuration. See the module documentation for the
    !! allowed mutations within one epoch. The plan is valid until the forest is mutated again.
    implicit none
    type(OctreeMesh3D),intent(in) :: forest
    integer,intent(in) :: nOld
    integer,intent(in) :: oldLeaf(1:nOld)
    type(TransferPlan3D),intent(out) :: plan
    ! Local
    integer :: i,li,walk,d,mx,steps
    logical :: found
    integer,allocatable :: oldElemOfNode(:)
    integer,allocatable :: famElem(:,:)
    integer,allocatable :: rev(:)

    do i = 1,nOld
      if(oldLeaf(i) < 1 .or. oldLeaf(i) > forest%nNodes) then
        print*,__FILE__,':',__LINE__, &
          ' : Error : oldLeaf snapshot contains a node id outside the forest node pool.'
        stop 1
      endif
    enddo

    ! Invert the snapshot: old element index by node id, and each old family's element indices
    ! by (child octant, parent node id). Coarsening zeroes the parent's child pointers but
    ! never the children's parent/octant entries, so famElem recovers detached families.
    allocate(oldElemOfNode(1:forest%nNodes))
    oldElemOfNode = 0
    allocate(famElem(1:8,1:forest%nNodes))
    famElem = 0
    do i = 1,nOld
      oldElemOfNode(oldLeaf(i)) = i
      if(forest%parent(oldLeaf(i)) > 0) then
        famElem(forest%octant(oldLeaf(i)),forest%parent(oldLeaf(i))) = i
      endif
    enddo

    plan%nOld = nOld
    plan%nNew = forest%nLeaves
    mx = max(forest%MaxLevel(),1) ! prolongation depth is bounded by the deepest leaf level
    allocate(plan%sourceKind(1:plan%nNew))
    allocate(plan%sourceElem(1:plan%nNew))
    allocate(plan%family(1:8,1:plan%nNew))
    allocate(plan%depth(1:plan%nNew))
    allocate(plan%path(1:mx,1:plan%nNew))
    plan%sourceKind = SELF_TRANSFER_COPY
    plan%sourceElem = 0
    plan%family = 0
    plan%depth = 0
    plan%path = 0
    plan%maxDepth = 0

    allocate(rev(1:mx))

    do li = 1,plan%nNew

      ! Ascend from the new leaf until reaching its data source; collect the octant taken at
      ! each level (bottom-up in rev, reversed into plan%path top-down). The chain has at most
      ! level+1 nodes, ending at the root.
      walk = forest%leaf(li)
      d = 0
      found = .false.
      do steps = 0,forest%level(forest%leaf(li))
        if(oldElemOfNode(walk) > 0) then
          if(d == 0) then
            plan%sourceKind(li) = SELF_TRANSFER_COPY
          else
            plan%sourceKind(li) = SELF_TRANSFER_PROLONG
          endif
          plan%sourceElem(li) = oldElemOfNode(walk)
          found = .true.
          exit
        elseif(famElem(1,walk) > 0 .and. famElem(2,walk) > 0 .and. &
               famElem(3,walk) > 0 .and. famElem(4,walk) > 0 .and. &
               famElem(5,walk) > 0 .and. famElem(6,walk) > 0 .and. &
               famElem(7,walk) > 0 .and. famElem(8,walk) > 0) then
          plan%sourceKind(li) = SELF_TRANSFER_RESTRICT
          plan%family(1:8,li) = famElem(1:8,walk)
          found = .true.
          exit
        endif
        if(forest%parent(walk) == 0) exit
        d = d+1
        rev(d) = forest%octant(walk)
        walk = forest%parent(walk)
      enddo

      if(.not. found) then
        print*,__FILE__,':',__LINE__, &
          ' : Error : new leaf has no old-leaf ancestor or coarsened old family; the snapshot'// &
          ' does not describe one adaptation epoch of this forest.'
        stop 1
      endif

      plan%depth(li) = d
      do i = 1,d
        plan%path(i,li) = rev(d-i+1)
      enddo
      plan%maxDepth = max(plan%maxDepth,d)

    enddo

    deallocate(oldElemOfNode,famElem,rev)

  endsubroutine BuildTransferPlan

  subroutine ApplyTransferPlan(plan,interp,nVar,uOld,uNew)
    !! Execute a transfer plan on nodal element data: uNew(:,:,:,li,:) receives old element
    !! data copied, prolonged (exact interpolation), or restricted (conservative L2 projection)
    !! according to plan entry li. interp must be the solution interpolant the data lives on
    !! (its mortar operators drive the transfer). Runs once per adaptation epoch - not a
    !! per-time-step hot path - so clarity is preferred over fused loops here.
    implicit none
    type(TransferPlan3D),intent(in) :: plan
    type(Lagrange),intent(in) :: interp
    integer,intent(in) :: nVar
    real(prec),intent(in) :: uOld(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:plan%nOld,1:nVar)
    real(prec),intent(out) :: uNew(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:plan%nNew,1:nVar)

    call ApplyTransferPlanRange(plan,interp,nVar,uOld,1,plan%nNew,uNew)

  endsubroutine ApplyTransferPlan

  subroutine ApplyTransferPlanRange(plan,interp,nVar,uOld,eFirst,eLast,uNew)
    !! Execute the contiguous sub-range eFirst..eLast of a transfer plan: uNew(:,:,:,k,:)
    !! receives the data of new element eFirst+k-1. uOld is the full (global) old field; the
    !! output is only the requested slice. This is the gather-then-slice (AMR Stage-5 v1)
    !! decomposed-mesh entry point: each rank passes its own contiguous range of the new
    !! element ordering and a gathered global old solution, and fills exactly its rank-local
    !! storage. ApplyTransferPlanWindow is the point-to-point (v2) entry point, which needs
    !! only the sub-range of the old field the requested new range actually references; this
    !! routine is that one over the whole old field, so the two are bit-identical by
    !! construction rather than by inspection.
    implicit none
    type(TransferPlan3D),intent(in) :: plan
    type(Lagrange),intent(in) :: interp
    integer,intent(in) :: nVar
    real(prec),intent(in) :: uOld(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:plan%nOld,1:nVar)
    integer,intent(in) :: eFirst
    integer,intent(in) :: eLast
    real(prec),intent(out) :: uNew(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:eLast-eFirst+1,1:nVar)

    call ApplyTransferPlanWindow(plan,interp,nVar,uOld,1,plan%nOld,eFirst,eLast,uNew)

  endsubroutine ApplyTransferPlanRange

  subroutine ApplyTransferPlanWindow(plan,interp,nVar,uOld,oldFirst,oldLast,eFirst,eLast,uNew)
    !! Execute the contiguous sub-range eFirst..eLast of a transfer plan against a WINDOW of
    !! the old field: uOld holds old elements oldFirst..oldLast only, indexed in the global old
    !! element numbering that plan%sourceElem and plan%family use, rather than the whole global
    !! field. uNew(:,:,:,k,:) receives the data of new element eFirst+k-1.
    !!
    !! This is the point-to-point migration (AMR Stage-5 v2) entry point. Because both the old
    !! and the new partition are contiguous ranges of the same space-filling-curve leaf order,
    !! the old elements a rank's new range references form a contiguous window that a rank
    !! computes locally from the (rank-replicated) plan - see PlanWindows in the AMR
    !! controllers - and receives point-to-point instead of allgathering the global field.
    !!
    !! Units and layout are those of MappedScalar3D %interior. The transfer operators are
    !! unchanged, so the operator identities hold as documented in the module header: exact
    !! prolongation, conservative L2 restriction.
    implicit none
    type(TransferPlan3D),intent(in) :: plan
    type(Lagrange),intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: oldFirst
    integer,intent(in) :: oldLast
    real(prec),intent(in) :: uOld(1:interp%N+1,1:interp%N+1,1:interp%N+1,oldFirst:oldLast,1:nVar)
    integer,intent(in) :: eFirst
    integer,intent(in) :: eLast
    real(prec),intent(out) :: uNew(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:eLast-eFirst+1,1:nVar)
    ! Local
    integer :: li,lo,c,step,Np
    real(prec),allocatable :: buf(:,:,:,:),kids(:,:,:,:,:),fam(:,:,:,:,:)

    if(.not. allocated(plan%sourceKind)) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : ApplyTransferPlan called with an unbuilt TransferPlan3D.'
      stop 1
    endif
    if(eFirst < 1 .or. eLast > plan%nNew .or. eLast < eFirst) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : ApplyTransferPlanRange called with a range outside 1..nNew.'
      stop 1
    endif
    if(oldFirst < 1 .or. oldLast > plan%nOld .or. oldLast < oldFirst) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : ApplyTransferPlanWindow called with an old-element window outside'// &
        ' 1..nOld.'
      stop 1
    endif

    Np = interp%N+1
    allocate(buf(1:Np,1:Np,1:Np,1:nVar))
    allocate(kids(1:Np,1:Np,1:Np,1:nVar,1:8))
    allocate(fam(1:Np,1:Np,1:Np,1:nVar,1:8))

    do li = eFirst,eLast
      lo = li-eFirst+1

      ! Every old element this entry reads must lie inside the supplied window. A window that
      ! does not cover the entry is a routing error (the window was computed from a different
      ! plan or a different partition), so report it rather than read outside the array. On the
      ! whole-field path the condition can never fire; the cost is a handful of integer
      ! comparisons per element, once per adaptation epoch.
      if(plan%sourceKind(li) == SELF_TRANSFER_RESTRICT) then
        do c = 1,8
          if(plan%family(c,li) < oldFirst .or. plan%family(c,li) > oldLast) then
            print*,__FILE__,':',__LINE__, &
              ' : Error : new element',li,' restricts from old element',plan%family(c,li), &
              ' outside the supplied old-element window',oldFirst,oldLast
            stop 1
          endif
        enddo
      elseif(plan%sourceElem(li) < oldFirst .or. plan%sourceElem(li) > oldLast) then
        print*,__FILE__,':',__LINE__, &
          ' : Error : new element',li,' sources old element',plan%sourceElem(li), &
          ' outside the supplied old-element window',oldFirst,oldLast
        stop 1
      endif

      if(plan%sourceKind(li) == SELF_TRANSFER_COPY) then
        uNew(1:Np,1:Np,1:Np,lo,1:nVar) = uOld(1:Np,1:Np,1:Np,plan%sourceElem(li),1:nVar)
        cycle
      endif

      if(plan%sourceKind(li) == SELF_TRANSFER_RESTRICT) then
        do c = 1,8
          fam(1:Np,1:Np,1:Np,1:nVar,c) = uOld(1:Np,1:Np,1:Np,plan%family(c,li),1:nVar)
        enddo
        call RestrictFromChildren(interp,nVar,fam,buf)
      else ! SELF_TRANSFER_PROLONG
        buf(1:Np,1:Np,1:Np,1:nVar) = uOld(1:Np,1:Np,1:Np,plan%sourceElem(li),1:nVar)
      endif

      do step = 1,plan%depth(li)
        call ProlongToChildren(interp,nVar,buf,kids)
        buf(1:Np,1:Np,1:Np,1:nVar) = kids(1:Np,1:Np,1:Np,1:nVar,plan%path(step,li))
      enddo

      uNew(1:Np,1:Np,1:Np,lo,1:nVar) = buf(1:Np,1:Np,1:Np,1:nVar)

    enddo

    deallocate(buf,kids,fam)

  endsubroutine ApplyTransferPlanWindow

endmodule SELF_TransferPlan_3D
