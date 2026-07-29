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

module SELF_TransferPlan_2D
!! Old-leaf -> new-leaf solution transfer plan for one adaptation epoch of the 2-D quad-forest:
!! the driver layer of AMR Stage 3 that connects the element-local transfer operators
!! (SELF_SolutionTransfer_2D) to the forest mutation (SELF_QuadTreeMesh_2D).
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
!!                            element's degree-N polynomial down the quadtree path (one step per
!!                            level, so depth > 1 handles a fresh child that 2:1 balancing
!!                            refined again in the same epoch). Exact, no loss.
!!   SELF_TRANSFER_RESTRICT : the leaf is (an ancestor of) a coarsened family; L2-project the
!!                            four old children family(1:4) onto their parent, then prolong down
!!                            depth >= 0 further steps (depth > 0 occurs when a just-coarsened
!!                            parent is immediately re-refined by 2:1 balancing). Conservative.
!!
!! The reconstruction is possible after the fact because forest node ids are stable: refinement
!! only appends nodes and coarsening only detaches children, whose level/parent/quadrant entries
!! persist. Each new leaf therefore ascends its parent chain until it meets either an old leaf
!! or a node holding a complete four-child old-leaf family; anything else means the snapshot
!! does not describe the epoch that produced the forest, and the builder stops with an error.
!!
!! ApplyTransferPlan executes the plan on nodal data u(1:N+1,1:N+1,1:nElem,1:nVar) (the layout
!! of MappedScalar2D %interior; units are those of the transferred fields) and inherits the
!! Stage-3 operator identities: prolongation is exact polynomial interpolation and restriction
!! is the conservative L2 projection, so the Jacobian-weighted cell integral of every variable
!! is conserved and a refine-then-coarsen round trip is the identity to roundoff. Transfer runs
!! once per adaptation epoch, between time steps; it is not a per-step hot path.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_QuadTreeMesh_2D
  use SELF_SolutionTransfer_2D

  implicit none

  integer,parameter :: SELF_TRANSFER_COPY = 0
  integer,parameter :: SELF_TRANSFER_PROLONG = 1
  integer,parameter :: SELF_TRANSFER_RESTRICT = 2

  type :: TransferPlan2D
    integer :: nOld = 0 !! old (pre-epoch) element count
    integer :: nNew = 0 !! new element count (= forest%nLeaves at build time)
    integer :: maxDepth = 0 !! largest prolongation depth in the plan
    integer,allocatable :: sourceKind(:) !! (nNew) SELF_TRANSFER_COPY / PROLONG / RESTRICT
    integer,allocatable :: sourceElem(:) !! (nNew) old element index (COPY/PROLONG; 0 otherwise)
    integer,allocatable :: family(:,:) !! (4,nNew) old element index per child quadrant (RESTRICT)
    integer,allocatable :: depth(:) !! (nNew) number of prolongation steps below the source
    integer,allocatable :: path(:,:) !! (>=maxDepth,nNew) quadrant taken at each step, top-down

  contains
    procedure,public :: Free => Free_TransferPlan2D

  endtype TransferPlan2D

contains

  subroutine Free_TransferPlan2D(this)
    implicit none
    class(TransferPlan2D),intent(inout) :: this

    if(allocated(this%sourceKind)) deallocate(this%sourceKind)
    if(allocated(this%sourceElem)) deallocate(this%sourceElem)
    if(allocated(this%family)) deallocate(this%family)
    if(allocated(this%depth)) deallocate(this%depth)
    if(allocated(this%path)) deallocate(this%path)
    this%nOld = 0
    this%nNew = 0
    this%maxDepth = 0

  endsubroutine Free_TransferPlan2D

  subroutine BuildTransferPlan(forest,nOld,oldLeaf,plan)
    !! Build the old->new transfer plan for the adaptation epoch that took the forest from the
    !! leaf configuration (nOld, oldLeaf) - a snapshot of (forest%nLeaves, forest%leaf) taken
    !! before mutating - to its current leaf configuration. See the module documentation for the
    !! allowed mutations within one epoch. The plan is valid until the forest is mutated again.
    implicit none
    type(QuadTreeMesh2D),intent(in) :: forest
    integer,intent(in) :: nOld
    integer,intent(in) :: oldLeaf(1:nOld)
    type(TransferPlan2D),intent(out) :: plan
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
    ! by (child quadrant, parent node id). Coarsening zeroes the parent's child pointers but
    ! never the children's parent/quadrant entries, so famElem recovers detached families.
    allocate(oldElemOfNode(1:forest%nNodes))
    oldElemOfNode = 0
    allocate(famElem(1:4,1:forest%nNodes))
    famElem = 0
    do i = 1,nOld
      oldElemOfNode(oldLeaf(i)) = i
      if(forest%parent(oldLeaf(i)) > 0) then
        famElem(forest%quadrant(oldLeaf(i)),forest%parent(oldLeaf(i))) = i
      endif
    enddo

    plan%nOld = nOld
    plan%nNew = forest%nLeaves
    mx = max(forest%MaxLevel(),1) ! prolongation depth is bounded by the deepest leaf level
    allocate(plan%sourceKind(1:plan%nNew))
    allocate(plan%sourceElem(1:plan%nNew))
    allocate(plan%family(1:4,1:plan%nNew))
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

      ! Ascend from the new leaf until reaching its data source; collect the quadrant taken at
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
               famElem(3,walk) > 0 .and. famElem(4,walk) > 0) then
          plan%sourceKind(li) = SELF_TRANSFER_RESTRICT
          plan%family(1:4,li) = famElem(1:4,walk)
          found = .true.
          exit
        endif
        if(forest%parent(walk) == 0) exit
        d = d+1
        rev(d) = forest%quadrant(walk)
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
    !! Execute a transfer plan on nodal element data: uNew(:,:,li,:) receives old element data
    !! copied, prolonged (exact interpolation), or restricted (conservative L2 projection)
    !! according to plan entry li. interp must be the solution interpolant the data lives on
    !! (its mortar operators drive the transfer). Runs once per adaptation epoch - not a
    !! per-time-step hot path - so clarity is preferred over fused loops here.
    implicit none
    type(TransferPlan2D),intent(in) :: plan
    type(Lagrange),intent(in) :: interp
    integer,intent(in) :: nVar
    real(prec),intent(in) :: uOld(1:interp%N+1,1:interp%N+1,1:plan%nOld,1:nVar)
    real(prec),intent(out) :: uNew(1:interp%N+1,1:interp%N+1,1:plan%nNew,1:nVar)
    ! Local
    integer :: li,c,step,Np
    real(prec),allocatable :: buf(:,:,:),kids(:,:,:,:),fam(:,:,:,:)

    if(.not. allocated(plan%sourceKind)) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : ApplyTransferPlan called with an unbuilt TransferPlan2D.'
      stop 1
    endif

    Np = interp%N+1
    allocate(buf(1:Np,1:Np,1:nVar))
    allocate(kids(1:Np,1:Np,1:nVar,1:4))
    allocate(fam(1:Np,1:Np,1:nVar,1:4))

    do li = 1,plan%nNew

      if(plan%sourceKind(li) == SELF_TRANSFER_COPY) then
        uNew(1:Np,1:Np,li,1:nVar) = uOld(1:Np,1:Np,plan%sourceElem(li),1:nVar)
        cycle
      endif

      if(plan%sourceKind(li) == SELF_TRANSFER_RESTRICT) then
        do c = 1,4
          fam(1:Np,1:Np,1:nVar,c) = uOld(1:Np,1:Np,plan%family(c,li),1:nVar)
        enddo
        call RestrictFromChildren(interp,nVar,fam,buf)
      else ! SELF_TRANSFER_PROLONG
        buf(1:Np,1:Np,1:nVar) = uOld(1:Np,1:Np,plan%sourceElem(li),1:nVar)
      endif

      do step = 1,plan%depth(li)
        call ProlongToChildren(interp,nVar,buf,kids)
        buf(1:Np,1:Np,1:nVar) = kids(1:Np,1:Np,1:nVar,plan%path(step,li))
      enddo

      uNew(1:Np,1:Np,li,1:nVar) = buf(1:Np,1:Np,1:nVar)

    enddo

    deallocate(buf,kids,fam)

  endsubroutine ApplyTransferPlan

endmodule SELF_TransferPlan_2D
