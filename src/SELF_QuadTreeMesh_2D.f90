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
    integer :: r

    this%nGeo = mesh%nGeo
    this%quadrature = mesh%quadrature
    this%nRoots = mesh%nElem

    allocate(this%rootCoords(1:2,1:mesh%nGeo+1,1:mesh%nGeo+1,1:mesh%nElem))
    this%rootCoords(1:2,1:mesh%nGeo+1,1:mesh%nGeo+1,1:mesh%nElem) = &
      mesh%nodeCoords(1:2,1:mesh%nGeo+1,1:mesh%nGeo+1,1:mesh%nElem)

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

endmodule SELF_QuadTreeMesh_2D
