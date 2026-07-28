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

module SELF_QuadTreeMesh_2D
!! Forest-of-quadtrees data structure for adaptive h-refinement of 2-D quadrilateral meshes
!! (AMR Stage 2b). Each base-mesh element is the root of a quadtree; refinement replaces a leaf
!! with four children (one per reference sub-quadrant, SELF child ordering 1=SW,2=SE,3=NE,4=NW),
!! and coarsening merges a family of four leaf siblings back into their parent. The forest tracks
!! refinement levels and the active leaf set; the physical geometry of any leaf is produced by
!! repeated exact isoparametric subdivision of its root (SELF_RefinementPrimitives_2D).
!!
!! This module owns the adaptive *mesh mutation* only. Turning an adaptively refined (and hence
!! generally nonconforming) forest into a solver-ready Mesh2D_t - face-neighbour queries, 2:1
!! balancing, hanging-node/mortar generation - is AMR Stage 4, and dynamic MPI re-partitioning is
!! Stage 5. The intended driver loop each adaptation step is:
!!
!!   indicator%Estimate(solution,ivar)      ! AMR Stage 1 : per-leaf refine/keep/coarsen flags
!!   forest%AdaptFromFlags(indicator%flag)  ! AMR Stage 2b : mutate the forest (this module)
!!   ... Stage 4 : balance + emit Mesh2D_t + mortars ; Stage 3 : transfer the solution ...
!!
!! Coarsening does not currently reclaim the storage of removed child nodes (the node arrays grow
!! monotonically); the active leaf set is always recovered by traversal from the roots, so
!! orphaned nodes are simply never revisited. Storage compaction is a future optimisation.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_RefinementPrimitives_2D

  implicit none

  ! Per-leaf adaptation flags (chosen to match SELF_RefinementIndicator_2D's SELF_AMR_* codes).
  integer,parameter :: QUADTREE_COARSEN = -1
  integer,parameter :: QUADTREE_KEEP = 0
  integer,parameter :: QUADTREE_REFINE = 1

  type :: QuadTreeMesh2D
    ! ---- Base-mesh (root) reference data ----
    integer :: nGeo = 0
    integer :: quadrature = 0
    integer :: nRoots = 0
    real(prec),allocatable :: rootCoords(:,:,:,:) ! (2,nGeo+1,nGeo+1,nRoots)
    ! Base-mesh face connectivity at the roots (conforming base assumed). rootNbr(s,r) is the
    ! base neighbour element of root r across local side s (0 = physical boundary), which is also
    ! that neighbour's root node id; rootNbrSide / rootFlip decode base sideInfo(4,s,r).
    integer,allocatable :: rootNbr(:,:) ! (4,nRoots)
    integer,allocatable :: rootNbrSide(:,:) ! (4,nRoots)
    integer,allocatable :: rootFlip(:,:) ! (4,nRoots)

    ! ---- Forest node storage (roots occupy node ids 1:nRoots) ----
    integer :: nNodes = 0
    integer :: capacity = 0
    integer,allocatable :: level(:) ! refinement level (0 at roots)
    integer,allocatable :: parent(:) ! parent node id (0 at roots)
    integer,allocatable :: quadrant(:) ! child index within parent (1..4; 0 at roots)
    integer,allocatable :: rootElem(:) ! base element this node descends from
    integer,allocatable :: child(:,:) ! (4,capacity) child node ids (0 => leaf)

    ! ---- Active leaf set (rebuilt after every mutation) ----
    integer :: nLeaves = 0
    integer,allocatable :: leaf(:) ! (nLeaves) leaf node ids, root-major DFS order

  contains
    procedure,public :: Init => Init_QuadTreeMesh2D
    procedure,public :: Free => Free_QuadTreeMesh2D
    procedure,public :: RefineNode => RefineNode_QuadTreeMesh2D
    procedure,public :: AdaptFromFlags => AdaptFromFlags_QuadTreeMesh2D
    procedure,public :: RebuildLeaves => RebuildLeaves_QuadTreeMesh2D
    procedure,public :: LeafCoords => LeafCoords_QuadTreeMesh2D
    procedure,public :: MaxLevel => MaxLevel_QuadTreeMesh2D
    procedure,public :: FaceNeighbor => FaceNeighbor_QuadTreeMesh2D
    procedure,public :: Balance2to1 => Balance2to1_QuadTreeMesh2D
    procedure,public :: MaxLevelJump => MaxLevelJump_QuadTreeMesh2D
    procedure,private :: EnsureCapacity => EnsureCapacity_QuadTreeMesh2D
  endtype QuadTreeMesh2D

contains

  subroutine Init_QuadTreeMesh2D(this,mesh)
    !! Initialize the forest with one root per base-mesh element (all leaves at level 0). The
    !! base geometry (node coordinates) is copied so leaf geometry can be regenerated after any
    !! amount of refinement without holding a reference to the mesh.
    implicit none
    class(QuadTreeMesh2D),intent(out) :: this
    type(Mesh2D),intent(in) :: mesh
    ! Local
    integer :: r,s

    this%nGeo = mesh%nGeo
    this%quadrature = mesh%quadrature
    this%nRoots = mesh%nElem

    allocate(this%rootCoords(1:2,1:mesh%nGeo+1,1:mesh%nGeo+1,1:mesh%nElem))
    this%rootCoords(1:2,1:mesh%nGeo+1,1:mesh%nGeo+1,1:mesh%nElem) = &
      mesh%nodeCoords(1:2,1:mesh%nGeo+1,1:mesh%nGeo+1,1:mesh%nElem)

    ! Root face connectivity from the (conforming) base mesh sideInfo.
    allocate(this%rootNbr(1:4,1:mesh%nElem))
    allocate(this%rootNbrSide(1:4,1:mesh%nElem))
    allocate(this%rootFlip(1:4,1:mesh%nElem))
    do r = 1,mesh%nElem
      do s = 1,4
        this%rootNbr(s,r) = mesh%sideInfo(3,s,r)
        this%rootNbrSide(s,r) = mesh%sideInfo(4,s,r)/10
        this%rootFlip(s,r) = mod(mesh%sideInfo(4,s,r),10)
      enddo
    enddo

    ! Roots are the first nRoots nodes.
    this%capacity = max(4*this%nRoots,16)
    this%nNodes = this%nRoots
    allocate(this%level(1:this%capacity))
    allocate(this%parent(1:this%capacity))
    allocate(this%quadrant(1:this%capacity))
    allocate(this%rootElem(1:this%capacity))
    allocate(this%child(1:4,1:this%capacity))
    this%level = 0
    this%parent = 0
    this%quadrant = 0
    this%rootElem = 0
    this%child = 0
    do r = 1,this%nRoots
      this%rootElem(r) = r
    enddo

    call this%RebuildLeaves()

  endsubroutine Init_QuadTreeMesh2D

  subroutine Free_QuadTreeMesh2D(this)
    implicit none
    class(QuadTreeMesh2D),intent(inout) :: this

    if(allocated(this%rootCoords)) deallocate(this%rootCoords)
    if(allocated(this%rootNbr)) deallocate(this%rootNbr)
    if(allocated(this%rootNbrSide)) deallocate(this%rootNbrSide)
    if(allocated(this%rootFlip)) deallocate(this%rootFlip)
    if(allocated(this%level)) deallocate(this%level)
    if(allocated(this%parent)) deallocate(this%parent)
    if(allocated(this%quadrant)) deallocate(this%quadrant)
    if(allocated(this%rootElem)) deallocate(this%rootElem)
    if(allocated(this%child)) deallocate(this%child)
    if(allocated(this%leaf)) deallocate(this%leaf)
    this%nGeo = 0
    this%nRoots = 0
    this%nNodes = 0
    this%capacity = 0
    this%nLeaves = 0

  endsubroutine Free_QuadTreeMesh2D

  subroutine EnsureCapacity_QuadTreeMesh2D(this,need)
    !! Grow the node arrays (amortized doubling) so at least `need` nodes fit.
    implicit none
    class(QuadTreeMesh2D),intent(inout) :: this
    integer,intent(in) :: need
    ! Local
    integer :: newCap
    integer,allocatable :: itmp(:),ctmp(:,:)

    if(need <= this%capacity) return

    newCap = this%capacity
    do while(newCap < need)
      newCap = 2*newCap
    enddo

    allocate(itmp(1:newCap))
    itmp = 0; itmp(1:this%nNodes) = this%level(1:this%nNodes); call move_alloc(itmp,this%level)
    allocate(itmp(1:newCap))
    itmp = 0; itmp(1:this%nNodes) = this%parent(1:this%nNodes); call move_alloc(itmp,this%parent)
    allocate(itmp(1:newCap))
    itmp = 0; itmp(1:this%nNodes) = this%quadrant(1:this%nNodes); call move_alloc(itmp,this%quadrant)
    allocate(itmp(1:newCap))
    itmp = 0; itmp(1:this%nNodes) = this%rootElem(1:this%nNodes); call move_alloc(itmp,this%rootElem)
    allocate(ctmp(1:4,1:newCap))
    ctmp = 0; ctmp(1:4,1:this%nNodes) = this%child(1:4,1:this%nNodes); call move_alloc(ctmp,this%child)

    this%capacity = newCap

  endsubroutine EnsureCapacity_QuadTreeMesh2D

  subroutine RefineNode_QuadTreeMesh2D(this,node)
    !! Subdivide one leaf node into four children. Does nothing (with a warning) if the node is
    !! already refined. The caller is responsible for rebuilding the leaf list afterwards (or
    !! calling AdaptFromFlags, which does so).
    implicit none
    class(QuadTreeMesh2D),intent(inout) :: this
    integer,intent(in) :: node
    ! Local
    integer :: c,newid

    if(this%child(1,node) /= 0) then
      print*,__FILE__,':',__LINE__,' : Warning : RefineNode called on an already-refined node.'
      return
    endif

    call this%EnsureCapacity(this%nNodes+4)
    do c = 1,4
      newid = this%nNodes+c
      this%level(newid) = this%level(node)+1
      this%parent(newid) = node
      this%quadrant(newid) = c
      this%rootElem(newid) = this%rootElem(node)
      this%child(1:4,newid) = 0
    enddo
    this%child(1,node) = this%nNodes+1
    this%child(2,node) = this%nNodes+2
    this%child(3,node) = this%nNodes+3
    this%child(4,node) = this%nNodes+4
    this%nNodes = this%nNodes+4

  endsubroutine RefineNode_QuadTreeMesh2D

  subroutine AdaptFromFlags_QuadTreeMesh2D(this,flag)
    !! Mutate the forest from a per-leaf flag array (indexed over the current leaves in this%leaf
    !! order, e.g. the flag array produced by SELF_RefinementIndicator_2D):
    !!
    !!   flag(i) = QUADTREE_REFINE  (+1) -> subdivide leaf i
    !!   flag(i) = QUADTREE_COARSEN (-1) -> merge leaf i's family if ALL four siblings are leaves
    !!                                      and ALL are flagged COARSEN (standard coarsening rule)
    !!   flag(i) = QUADTREE_KEEP     (0) -> unchanged
    !!
    !! Refinement and coarsening are both resolved against the pre-adaptation leaf snapshot, so
    !! the two operations never interfere. The leaf list is rebuilt on return.
    implicit none
    class(QuadTreeMesh2D),intent(inout) :: this
    integer,intent(in) :: flag(1:this%nLeaves)
    ! Local
    integer :: i,c,p,node,nSnap
    integer,allocatable :: flagOfNode(:)
    integer,allocatable :: refineList(:)
    logical :: family

    nSnap = this%nLeaves

    ! Map the flag onto node ids so coarsening can test whole families by node.
    allocate(flagOfNode(1:this%nNodes))
    flagOfNode = QUADTREE_KEEP
    allocate(refineList(1:nSnap))
    refineList = 0
    do i = 1,nSnap
      flagOfNode(this%leaf(i)) = flag(i)
      if(flag(i) == QUADTREE_REFINE) refineList(i) = this%leaf(i)
    enddo

    ! ---- Coarsening : merge families all of whose four leaf children are flagged COARSEN ----
    do i = 1,nSnap
      if(flag(i) /= QUADTREE_COARSEN) cycle
      node = this%leaf(i)
      p = this%parent(node)
      if(p == 0) cycle ! a root cannot be coarsened
      ! Only act once per family (when processing its first child), and only if every child is a
      ! leaf flagged COARSEN.
      if(this%child(1,p) /= node) cycle
      family = .true.
      do c = 1,4
        if(this%child(1,this%child(c,p)) /= 0) family = .false. ! child not a leaf
        if(flagOfNode(this%child(c,p)) /= QUADTREE_COARSEN) family = .false.
      enddo
      if(family) this%child(1:4,p) = 0 ! detach children -> p becomes a leaf again
    enddo

    ! ---- Refinement : subdivide flagged leaves that are still leaves ----
    do i = 1,nSnap
      if(refineList(i) == 0) cycle
      node = refineList(i)
      if(this%child(1,node) == 0) call this%RefineNode(node)
    enddo

    deallocate(flagOfNode,refineList)

    call this%RebuildLeaves()

  endsubroutine AdaptFromFlags_QuadTreeMesh2D

  subroutine RebuildLeaves_QuadTreeMesh2D(this)
    !! Recompute the active leaf set by depth-first traversal from the roots. Traversal (rather
    !! than a scan of all nodes) is what makes orphaned nodes left behind by coarsening invisible.
    implicit none
    class(QuadTreeMesh2D),intent(inout) :: this
    ! Local
    integer :: r,n

    ! First pass: count leaves.
    n = 0
    do r = 1,this%nRoots
      call count_leaves(r,n)
    enddo
    this%nLeaves = n
    if(allocated(this%leaf)) deallocate(this%leaf)
    allocate(this%leaf(1:max(n,1)))
    this%leaf = 0

    ! Second pass: collect leaves in root-major DFS order.
    n = 0
    do r = 1,this%nRoots
      call collect_leaves(r,n)
    enddo

  contains

    recursive subroutine count_leaves(node,cnt)
      integer,intent(in) :: node
      integer,intent(inout) :: cnt
      integer :: k
      if(this%child(1,node) == 0) then
        cnt = cnt+1
      else
        do k = 1,4
          call count_leaves(this%child(k,node),cnt)
        enddo
      endif
    endsubroutine count_leaves

    recursive subroutine collect_leaves(node,idx)
      integer,intent(in) :: node
      integer,intent(inout) :: idx
      integer :: k
      if(this%child(1,node) == 0) then
        idx = idx+1
        this%leaf(idx) = node
      else
        do k = 1,4
          call collect_leaves(this%child(k,node),idx)
        enddo
      endif
    endsubroutine collect_leaves

  endsubroutine RebuildLeaves_QuadTreeMesh2D

  function MaxLevel_QuadTreeMesh2D(this) result(mx)
    !! Highest refinement level among the active leaves.
    implicit none
    class(QuadTreeMesh2D),intent(in) :: this
    integer :: mx
    integer :: i

    mx = 0
    do i = 1,this%nLeaves
      mx = max(mx,this%level(this%leaf(i)))
    enddo

  endfunction MaxLevel_QuadTreeMesh2D

  subroutine LeafCoords_QuadTreeMesh2D(this,leafIndex,geomInterp,coords)
    !! Physical geometry-node coordinates of leaf `leafIndex`, produced by repeated exact
    !! isoparametric subdivision of its root element along the quadtree path. geomInterp must be a
    !! degree-nGeo Lagrange interpolant on the mesh's geometry (quadrature) nodes - the same
    !! interpolant SELF_MeshRefinement_2D builds. Level 0 leaves return the root geometry directly.
    implicit none
    class(QuadTreeMesh2D),intent(in) :: this
    integer,intent(in) :: leafIndex
    type(Lagrange),intent(in) :: geomInterp
    real(prec),intent(out) :: coords(1:2,1:this%nGeo+1,1:this%nGeo+1)
    ! Local
    integer :: node,lvl,step
    integer,allocatable :: path(:)
    real(prec),allocatable :: cur(:,:,:),kids(:,:,:,:)

    node = this%leaf(leafIndex)
    lvl = this%level(node)

    allocate(cur(1:2,1:this%nGeo+1,1:this%nGeo+1))

    if(lvl == 0) then
      cur(1:2,:,:) = this%rootCoords(1:2,:,:,this%rootElem(node))
      coords = cur
      deallocate(cur)
      return
    endif

    ! Path of quadrant indices from root (step 1) down to the leaf (step lvl).
    allocate(path(1:lvl))
    do step = lvl,1,-1
      path(step) = this%quadrant(node)
      node = this%parent(node)
    enddo
    ! `node` is now the root.
    cur(1:2,:,:) = this%rootCoords(1:2,:,:,this%rootElem(this%leaf(leafIndex)))

    allocate(kids(1:2,1:this%nGeo+1,1:this%nGeo+1,1:4))
    do step = 1,lvl
      call SubdivideNodeCoords(geomInterp,this%nGeo,cur,kids)
      cur(1:2,:,:) = kids(1:2,:,:,path(step))
    enddo
    coords = cur

    deallocate(cur,kids,path)

  endsubroutine LeafCoords_QuadTreeMesh2D

  recursive subroutine FaceNeighbor_QuadTreeMesh2D(this,node,s,nbr,ns,nf)
    !! Find the equal-or-larger face neighbour of `node` across its local side s using the
    !! classic quadtree ascend/descend search. Returns:
    !!   nbr - neighbour node id (0 = physical domain boundary). It is either a LEAF at any level
    !!         <= level(node), or an INTERNAL node at exactly level(node) (meaning the shared face
    !!         is subdivided on the neighbour side, i.e. finer neighbours exist).
    !!   ns  - the neighbour's local side that faces `node`.
    !!   nf  - the flip between the two shared edges (0 same direction, 1 reversed), inherited
    !!         from the base-mesh face where the search crosses a root boundary.
    !! With this, a 2:1 hanging face is exactly "nbr is a leaf with level(nbr) = level(node)-1",
    !! and finer neighbours are exactly "nbr is internal".
    implicit none
    class(QuadTreeMesh2D),intent(in) :: this
    integer,intent(in) :: node,s
    integer,intent(out) :: nbr,ns,nf
    ! Local
    integer :: p,c,pnbr,ps,pf,t,tq

    if(this%level(node) == 0) then
      ! Cross a base-mesh face: neighbour root, its side, and the base flip.
      nbr = this%rootNbr(s,this%rootElem(node))
      ns = this%rootNbrSide(s,this%rootElem(node))
      nf = this%rootFlip(s,this%rootElem(node))
      return
    endif

    p = this%parent(node)
    c = this%quadrant(node)

    if(qt_internal(s,c)) then
      ! Neighbour is the sibling on the other side of an interior face of the parent.
      nbr = this%child(qt_reflect(s,c),p)
      ns = qt_opposite(s)
      nf = 0
      return
    endif

    ! Otherwise ascend: find the parent's neighbour across the same side, then descend.
    call FaceNeighbor_QuadTreeMesh2D(this,p,s,pnbr,ps,pf)
    if(pnbr == 0) then
      nbr = 0; ns = 0; nf = 0
      return
    endif

    if(this%child(1,pnbr) == 0) then
      ! Parent's neighbour is a leaf (equal or larger than the parent) -> our larger neighbour.
      nbr = pnbr; ns = ps; nf = pf
      return
    endif

    ! Parent's neighbour is internal (at level(node)-1): descend one level to the child that
    ! borders `node`, matching sub-positions across the face through the flip.
    t = qt_subpos(c,s)
    if(pf == 0) then
      tq = t
    else
      tq = 3-t
    endif
    nbr = this%child(childOfSide(tq,ps),pnbr)
    ns = ps
    nf = pf

  endsubroutine FaceNeighbor_QuadTreeMesh2D

  subroutine Balance2to1_QuadTreeMesh2D(this)
    !! Enforce the 2:1 balance condition: no leaf face may separate elements differing by more
    !! than one refinement level. Iterates to a fixed point - in each sweep, any leaf whose
    !! equal-or-larger neighbour is a leaf two or more levels coarser triggers refinement of that
    !! coarser neighbour; refinement can ripple, so sweeps repeat until nothing changes. The leaf
    !! set is rebuilt on return.
    implicit none
    class(QuadTreeMesh2D),intent(inout) :: this
    ! Local
    integer :: li,s,node,nbr,ns,nf,nSnap
    integer,allocatable :: snap(:)
    logical :: changed

    do
      changed = .false.
      nSnap = this%nLeaves
      allocate(snap(1:nSnap))
      snap(1:nSnap) = this%leaf(1:nSnap)

      do li = 1,nSnap
        node = snap(li)
        if(this%child(1,node) /= 0) cycle ! refined earlier in this sweep
        do s = 1,4
          call this%FaceNeighbor(node,s,nbr,ns,nf)
          if(nbr /= 0) then
            if(this%child(1,nbr) == 0 .and. this%level(nbr) <= this%level(node)-2) then
              call this%RefineNode(nbr)
              changed = .true.
            endif
          endif
        enddo
      enddo

      deallocate(snap)
      call this%RebuildLeaves()
      if(.not. changed) exit
    enddo

  endsubroutine Balance2to1_QuadTreeMesh2D

  function MaxLevelJump_QuadTreeMesh2D(this) result(mx)
    !! Largest refinement-level difference across any leaf face (0 on a conforming or uniformly
    !! refined forest, 1 on a 2:1-balanced adaptive forest). Because FaceNeighbor returns the
    !! equal-or-larger neighbour, every level difference is observed from the finer leaf as a
    !! coarser leaf neighbour; internal (finer) neighbours contribute nothing from this side.
    implicit none
    class(QuadTreeMesh2D),intent(in) :: this
    integer :: mx
    ! Local
    integer :: li,s,node,nbr,ns,nf

    mx = 0
    do li = 1,this%nLeaves
      node = this%leaf(li)
      do s = 1,4
        call this%FaceNeighbor(node,s,nbr,ns,nf)
        if(nbr /= 0) then
          if(this%child(1,nbr) == 0) mx = max(mx,this%level(node)-this%level(nbr))
        endif
      enddo
    enddo

  endfunction MaxLevelJump_QuadTreeMesh2D

  ! -------------------------------------------------------------------------------------------- !
  ! Quadtree face-adjacency helpers (SELF child ordering 1=SW,2=SE,3=NE,4=NW; sides 1=S,2=E,3=N,
  ! 4=W). See SELF_RefinementPrimitives_2D for the quadrant/side conventions and childOfSide.
  ! -------------------------------------------------------------------------------------------- !

  pure function qt_opposite(s) result(o)
    !! The local side directly across an element from side s.
    implicit none
    integer,intent(in) :: s
    integer :: o
    integer,parameter :: opp(1:4) = [3,4,1,2]
    o = opp(s)
  endfunction qt_opposite

  pure function qt_internal(s,c) result(isInternal)
    !! .true. if child c's side s is interior to its parent (its neighbour across s is a sibling).
    implicit none
    integer,intent(in) :: s,c
    logical :: isInternal
    select case(s)
    case(1); isInternal = (c == 3 .or. c == 4) ! South interior for the top children
    case(2); isInternal = (c == 1 .or. c == 4) ! East interior for the left children
    case(3); isInternal = (c == 1 .or. c == 2) ! North interior for the bottom children
    case(4); isInternal = (c == 2 .or. c == 3) ! West interior for the right children
    case default; isInternal = .false.
    endselect
  endfunction qt_internal

  pure function qt_reflect(s,c) result(rc)
    !! Sibling child index obtained by reflecting c across side s (swap the y-half for the
    !! horizontal faces S/N, the x-half for the vertical faces E/W).
    implicit none
    integer,intent(in) :: s,c
    integer :: rc
    integer,parameter :: vref(1:4) = [4,3,2,1] ! swap ay (s = 1,3)
    integer,parameter :: href(1:4) = [2,1,4,3] ! swap ax (s = 2,4)
    if(s == 1 .or. s == 3) then
      rc = vref(c)
    else
      rc = href(c)
    endif
  endfunction qt_reflect

  pure function qt_subpos(c,s) result(t)
    !! Sub-position (1 or 2, in the positive direction of side s) of child c along side s;
    !! 0 if c does not touch side s.
    implicit none
    integer,intent(in) :: c,s
    integer :: t
    if(childOfSide(1,s) == c) then
      t = 1
    elseif(childOfSide(2,s) == c) then
      t = 2
    else
      t = 0
    endif
  endfunction qt_subpos

endmodule SELF_QuadTreeMesh_2D
