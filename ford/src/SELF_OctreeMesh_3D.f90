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

module SELF_OctreeMesh_3D
!! Forest-of-octrees data structure for adaptive h-refinement of 3-D hexahedral meshes,
!! the direct analogue of SELF_QuadTreeMesh_2D. Each base-mesh element is the root of an
!! octree; refinement replaces a leaf with eight children (one per reference octant, SELF
!! child ordering = CGNS corner order), and coarsening merges a family of eight leaf
!! siblings back into their parent. The forest tracks refinement levels and the active
!! leaf set; the physical geometry of any leaf is produced by repeated exact
!! isoparametric subdivision of its root (SELF_RefinementPrimitives_3D).
!!
!! This module owns the adaptive *mesh mutation* only. Turning an adaptively refined (and
!! hence generally nonconforming) forest into a solver-ready Mesh3D_t - face-neighbour
!! queries, 2:1 balancing, hanging-face/mortar generation - is SELF_AdaptiveMesh_3D.
!! The intended driver loop each adaptation step matches the 2-D one:
!!
!!   indicator%Estimate(solution,ivar)      ! per-leaf refine/keep/coarsen flags
!!   forest%AdaptFromFlags(indicator%flag)  ! mutate the forest (this module)
!!   forest%Balance2to1() ; EmitMesh(...)   ! balance + emit Mesh3D_t + mortars
!!
!! 2:1 balance is enforced across faces only : hanging edges and corners carry no
!! interface data in the DG discretization (face mortars carry all of it), exactly as the
!! 2-D quadtree tolerates corner-level jumps.
!!
!! Coarsening does not currently reclaim the storage of removed child nodes (the node
!! arrays grow monotonically); the active leaf set is always recovered by traversal from
!! the roots, so orphaned nodes are simply never revisited. Storage compaction is a
!! future optimisation.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_3D
  use SELF_RefinementPrimitives_3D

  implicit none

  ! Per-leaf adaptation flags (chosen to match the refinement indicator's SELF_AMR_* codes).
  integer,parameter :: OCTREE_COARSEN = -1
  integer,parameter :: OCTREE_KEEP = 0
  integer,parameter :: OCTREE_REFINE = 1

  type :: OctreeMesh3D
    ! ---- Base-mesh (root) reference data ----
    integer :: nGeo = 0
    integer :: quadrature = 0
    integer :: nRoots = 0
    real(prec),allocatable :: rootCoords(:,:,:,:,:) ! (3,nGeo+1,nGeo+1,nGeo+1,nRoots)
    ! Base-mesh face connectivity at the roots (conforming base assumed). rootNbr(s,r) is
    ! the base neighbour element of root r across local face s (0 = physical boundary),
    ! which is also that neighbour's root node id; rootNbrSide / rootFlip decode base
    ! sideInfo(4,s,r).
    integer,allocatable :: rootNbr(:,:) ! (6,nRoots)
    integer,allocatable :: rootNbrSide(:,:) ! (6,nRoots)
    integer,allocatable :: rootFlip(:,:) ! (6,nRoots)
    integer,allocatable :: rootBC(:,:) ! (6,nRoots) base boundary-condition id per face
    integer,allocatable :: rootMaterial(:) ! (nRoots) base material id per root element

    ! ---- Forest node storage (roots occupy node ids 1:nRoots) ----
    integer :: nNodes = 0
    integer :: capacity = 0
    integer,allocatable :: level(:) ! refinement level (0 at roots)
    integer,allocatable :: parent(:) ! parent node id (0 at roots)
    integer,allocatable :: octant(:) ! child index within parent (1..8; 0 at roots)
    integer,allocatable :: rootElem(:) ! base element this node descends from
    integer,allocatable :: child(:,:) ! (8,capacity) child node ids (0 => leaf)

    ! ---- Active leaf set (rebuilt after every mutation) ----
    integer :: nLeaves = 0
    integer,allocatable :: leaf(:) ! (nLeaves) leaf node ids, root-major DFS order

  contains
    procedure,public :: Init => Init_OctreeMesh3D
    procedure,public :: InitGlobal => InitGlobal_OctreeMesh3D
    procedure,public :: Free => Free_OctreeMesh3D
    procedure,public :: RefineNode => RefineNode_OctreeMesh3D
    procedure,public :: AdaptFromFlags => AdaptFromFlags_OctreeMesh3D
    procedure,public :: RebuildLeaves => RebuildLeaves_OctreeMesh3D
    procedure,public :: LeafCoords => LeafCoords_OctreeMesh3D
    procedure,public :: MaxLevel => MaxLevel_OctreeMesh3D
    procedure,public :: FaceNeighbor => FaceNeighbor_OctreeMesh3D
    procedure,public :: Balance2to1 => Balance2to1_OctreeMesh3D
    procedure,public :: MaxLevelJump => MaxLevelJump_OctreeMesh3D
    procedure,private :: EnsureCapacity => EnsureCapacity_OctreeMesh3D
  endtype OctreeMesh3D

contains

  subroutine Init_OctreeMesh3D(this,mesh)
    !! Initialize the forest with one root per base-mesh element (all leaves at level 0).
    !! The base geometry (node coordinates) is copied so leaf geometry can be regenerated
    !! after any amount of refinement without holding a reference to the mesh. Requires a
    !! single-rank mesh (a decomposed mesh only stores its local elements); a
    !! rank-replicated forest over a decomposed base is built by gathering the global
    !! tables and calling InitGlobal.
    implicit none
    class(OctreeMesh3D),intent(out) :: this
    type(Mesh3D),intent(in) :: mesh
    ! Local
    integer :: r,s
    integer,allocatable :: nbr(:,:),nbrSide(:,:),flip(:,:),bc(:,:),mat(:)

    if(mesh%decomp%nRanks > 1) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : OctreeMesh3D%Init requires a single-rank mesh; gather the global base'// &
        ' tables and call InitGlobal for a decomposed base mesh.'
      stop 1
    endif

    allocate(nbr(1:6,1:mesh%nElem),nbrSide(1:6,1:mesh%nElem))
    allocate(flip(1:6,1:mesh%nElem),bc(1:6,1:mesh%nElem))
    allocate(mat(1:mesh%nElem))
    do r = 1,mesh%nElem
      do s = 1,6
        nbr(s,r) = mesh%sideInfo(3,s,r)
        nbrSide(s,r) = mesh%sideInfo(4,s,r)/10
        flip(s,r) = mod(mesh%sideInfo(4,s,r),10)
        bc(s,r) = mesh%sideInfo(5,s,r)
      enddo
      mat(r) = mesh%elemMaterial(r)
    enddo

    call this%InitGlobal(mesh%nElem,mesh%nGeo,mesh%quadrature, &
                         mesh%nodeCoords,nbr,nbrSide,flip,bc,mat)

    deallocate(nbr,nbrSide,flip,bc,mat)

  endsubroutine Init_OctreeMesh3D

  subroutine InitGlobal_OctreeMesh3D(this,nRoots,nGeo,quadrature,rootCoords, &
                                     rootNbr,rootNbrSide,rootFlip,rootBC,rootMaterial)
    !! Initialize the forest directly from GLOBAL base-mesh tables (one root per global
    !! base element, all leaves at level 0). This is the initialization path for a
    !! rank-replicated forest over a decomposed base mesh: every rank passes the same
    !! gathered tables and holds an identical forest. rootNbr carries global element ids
    !! (0 = physical boundary), rootNbrSide/rootFlip decode the base sideInfo(4) pairing,
    !! rootBC the base boundary-condition id per face, and rootMaterial the base material
    !! id per element.
    implicit none
    class(OctreeMesh3D),intent(out) :: this
    integer,intent(in) :: nRoots
    integer,intent(in) :: nGeo
    integer,intent(in) :: quadrature
    real(prec),intent(in) :: rootCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nRoots)
    integer,intent(in) :: rootNbr(1:6,1:nRoots)
    integer,intent(in) :: rootNbrSide(1:6,1:nRoots)
    integer,intent(in) :: rootFlip(1:6,1:nRoots)
    integer,intent(in) :: rootBC(1:6,1:nRoots)
    integer,intent(in) :: rootMaterial(1:nRoots)
    ! Local
    integer :: r

    this%nGeo = nGeo
    this%quadrature = quadrature
    this%nRoots = nRoots

    allocate(this%rootCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nRoots))
    this%rootCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nRoots) = &
      rootCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nRoots)

    allocate(this%rootNbr(1:6,1:nRoots))
    allocate(this%rootNbrSide(1:6,1:nRoots))
    allocate(this%rootFlip(1:6,1:nRoots))
    allocate(this%rootBC(1:6,1:nRoots))
    allocate(this%rootMaterial(1:nRoots))
    this%rootNbr(1:6,1:nRoots) = rootNbr(1:6,1:nRoots)
    this%rootNbrSide(1:6,1:nRoots) = rootNbrSide(1:6,1:nRoots)
    this%rootFlip(1:6,1:nRoots) = rootFlip(1:6,1:nRoots)
    this%rootBC(1:6,1:nRoots) = rootBC(1:6,1:nRoots)
    this%rootMaterial(1:nRoots) = rootMaterial(1:nRoots)

    ! Roots are the first nRoots nodes.
    this%capacity = max(8*this%nRoots,16)
    this%nNodes = this%nRoots
    allocate(this%level(1:this%capacity))
    allocate(this%parent(1:this%capacity))
    allocate(this%octant(1:this%capacity))
    allocate(this%rootElem(1:this%capacity))
    allocate(this%child(1:8,1:this%capacity))
    this%level = 0
    this%parent = 0
    this%octant = 0
    this%rootElem = 0
    this%child = 0
    do r = 1,this%nRoots
      this%rootElem(r) = r
    enddo

    call this%RebuildLeaves()

  endsubroutine InitGlobal_OctreeMesh3D

  subroutine Free_OctreeMesh3D(this)
    implicit none
    class(OctreeMesh3D),intent(inout) :: this

    if(allocated(this%rootCoords)) deallocate(this%rootCoords)
    if(allocated(this%rootNbr)) deallocate(this%rootNbr)
    if(allocated(this%rootNbrSide)) deallocate(this%rootNbrSide)
    if(allocated(this%rootFlip)) deallocate(this%rootFlip)
    if(allocated(this%rootBC)) deallocate(this%rootBC)
    if(allocated(this%rootMaterial)) deallocate(this%rootMaterial)
    if(allocated(this%level)) deallocate(this%level)
    if(allocated(this%parent)) deallocate(this%parent)
    if(allocated(this%octant)) deallocate(this%octant)
    if(allocated(this%rootElem)) deallocate(this%rootElem)
    if(allocated(this%child)) deallocate(this%child)
    if(allocated(this%leaf)) deallocate(this%leaf)
    this%nGeo = 0
    this%nRoots = 0
    this%nNodes = 0
    this%capacity = 0
    this%nLeaves = 0

  endsubroutine Free_OctreeMesh3D

  subroutine EnsureCapacity_OctreeMesh3D(this,need)
    !! Grow the node arrays (amortized doubling) so at least `need` nodes fit.
    implicit none
    class(OctreeMesh3D),intent(inout) :: this
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
    itmp = 0
    itmp(1:this%nNodes) = this%level(1:this%nNodes)
    call move_alloc(itmp,this%level)
    allocate(itmp(1:newCap))
    itmp = 0
    itmp(1:this%nNodes) = this%parent(1:this%nNodes)
    call move_alloc(itmp,this%parent)
    allocate(itmp(1:newCap))
    itmp = 0
    itmp(1:this%nNodes) = this%octant(1:this%nNodes)
    call move_alloc(itmp,this%octant)
    allocate(itmp(1:newCap))
    itmp = 0
    itmp(1:this%nNodes) = this%rootElem(1:this%nNodes)
    call move_alloc(itmp,this%rootElem)
    allocate(ctmp(1:8,1:newCap))
    ctmp = 0
    ctmp(1:8,1:this%nNodes) = this%child(1:8,1:this%nNodes)
    call move_alloc(ctmp,this%child)

    this%capacity = newCap

  endsubroutine EnsureCapacity_OctreeMesh3D

  subroutine RefineNode_OctreeMesh3D(this,node)
    !! Subdivide one leaf node into eight children. Does nothing (with a warning) if the
    !! node is already refined. The caller is responsible for rebuilding the leaf list
    !! afterwards (or calling AdaptFromFlags, which does so).
    implicit none
    class(OctreeMesh3D),intent(inout) :: this
    integer,intent(in) :: node
    ! Local
    integer :: c,newid

    if(this%child(1,node) /= 0) then
      print*,__FILE__,':',__LINE__,' : Warning : RefineNode called on an already-refined node.'
      return
    endif

    call this%EnsureCapacity(this%nNodes+8)
    do c = 1,8
      newid = this%nNodes+c
      this%level(newid) = this%level(node)+1
      this%parent(newid) = node
      this%octant(newid) = c
      this%rootElem(newid) = this%rootElem(node)
      this%child(1:8,newid) = 0
      this%child(c,node) = newid
    enddo
    this%nNodes = this%nNodes+8

  endsubroutine RefineNode_OctreeMesh3D

  subroutine AdaptFromFlags_OctreeMesh3D(this,flag)
    !! Mutate the forest from a per-leaf flag array (indexed over the current leaves in
    !! this%leaf order, e.g. the flag array produced by the refinement indicator):
    !!
    !!   flag(i) = OCTREE_REFINE  (+1) -> subdivide leaf i
    !!   flag(i) = OCTREE_COARSEN (-1) -> merge leaf i's family if ALL eight siblings are
    !!                                    leaves and ALL are flagged COARSEN
    !!   flag(i) = OCTREE_KEEP     (0) -> unchanged
    !!
    !! Refinement and coarsening are both resolved against the pre-adaptation leaf
    !! snapshot, so the two operations never interfere. The leaf list is rebuilt on
    !! return.
    implicit none
    class(OctreeMesh3D),intent(inout) :: this
    integer,intent(in) :: flag(1:this%nLeaves)
    ! Local
    integer :: i,c,p,node,nSnap
    integer,allocatable :: flagOfNode(:)
    integer,allocatable :: refineList(:)
    logical :: family

    nSnap = this%nLeaves

    ! Map the flag onto node ids so coarsening can test whole families by node.
    allocate(flagOfNode(1:this%nNodes))
    flagOfNode = OCTREE_KEEP
    allocate(refineList(1:nSnap))
    refineList = 0
    do i = 1,nSnap
      flagOfNode(this%leaf(i)) = flag(i)
      if(flag(i) == OCTREE_REFINE) refineList(i) = this%leaf(i)
    enddo

    ! ---- Coarsening : merge families all of whose eight leaf children are flagged COARSEN ----
    do i = 1,nSnap
      if(flag(i) /= OCTREE_COARSEN) cycle
      node = this%leaf(i)
      p = this%parent(node)
      if(p == 0) cycle ! a root cannot be coarsened
      ! Only act once per family (when processing its first child), and only if every
      ! child is a leaf flagged COARSEN.
      if(this%child(1,p) /= node) cycle
      family = .true.
      do c = 1,8
        if(this%child(1,this%child(c,p)) /= 0) family = .false. ! child not a leaf
        if(flagOfNode(this%child(c,p)) /= OCTREE_COARSEN) family = .false.
      enddo
      if(family) this%child(1:8,p) = 0 ! detach children -> p becomes a leaf again
    enddo

    ! ---- Refinement : subdivide flagged leaves that are still leaves ----
    do i = 1,nSnap
      if(refineList(i) == 0) cycle
      node = refineList(i)
      if(this%child(1,node) == 0) call this%RefineNode(node)
    enddo

    deallocate(flagOfNode,refineList)

    call this%RebuildLeaves()

  endsubroutine AdaptFromFlags_OctreeMesh3D

  subroutine RebuildLeaves_OctreeMesh3D(this)
    !! Recompute the active leaf set by depth-first traversal from the roots. Traversal
    !! (rather than a scan of all nodes) is what makes orphaned nodes left behind by
    !! coarsening invisible. The root-major DFS order with children visited 1..8 is the
    !! Morton (Z-order) curve within each root, which is what makes the contiguous block
    !! decomposition of the emitted element list a locality-preserving SFC partition.
    implicit none
    class(OctreeMesh3D),intent(inout) :: this
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
        do k = 1,8
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
        do k = 1,8
          call collect_leaves(this%child(k,node),idx)
        enddo
      endif
    endsubroutine collect_leaves

  endsubroutine RebuildLeaves_OctreeMesh3D

  function MaxLevel_OctreeMesh3D(this) result(mx)
    !! Highest refinement level among the active leaves.
    implicit none
    class(OctreeMesh3D),intent(in) :: this
    integer :: mx
    integer :: i

    mx = 0
    do i = 1,this%nLeaves
      mx = max(mx,this%level(this%leaf(i)))
    enddo

  endfunction MaxLevel_OctreeMesh3D

  subroutine LeafCoords_OctreeMesh3D(this,leafIndex,geomInterp,coords)
    !! Physical geometry-node coordinates of leaf `leafIndex`, produced by repeated exact
    !! isoparametric subdivision of its root element along the octree path. geomInterp
    !! must be a degree-nGeo Lagrange interpolant on the mesh's geometry (quadrature)
    !! nodes - the same interpolant SELF_MeshRefinement_3D builds. Level 0 leaves return
    !! the root geometry directly. Pure function of (root coords, level, octant path), so
    !! regenerated leaf geometry is bit-identical across epochs.
    implicit none
    class(OctreeMesh3D),intent(in) :: this
    integer,intent(in) :: leafIndex
    type(Lagrange),intent(in) :: geomInterp
    real(prec),intent(out) :: coords(1:3,1:this%nGeo+1,1:this%nGeo+1,1:this%nGeo+1)
    ! Local
    integer :: node,lvl,step
    integer,allocatable :: path(:)
    real(prec),allocatable :: cur(:,:,:,:),kids(:,:,:,:,:)

    node = this%leaf(leafIndex)
    lvl = this%level(node)

    allocate(cur(1:3,1:this%nGeo+1,1:this%nGeo+1,1:this%nGeo+1))

    if(lvl == 0) then
      cur(1:3,:,:,:) = this%rootCoords(1:3,:,:,:,this%rootElem(node))
      coords = cur
      deallocate(cur)
      return
    endif

    ! Path of octant indices from root (step 1) down to the leaf (step lvl).
    allocate(path(1:lvl))
    do step = lvl,1,-1
      path(step) = this%octant(node)
      node = this%parent(node)
    enddo
    ! `node` is now the root.
    cur(1:3,:,:,:) = this%rootCoords(1:3,:,:,:,this%rootElem(this%leaf(leafIndex)))

    allocate(kids(1:3,1:this%nGeo+1,1:this%nGeo+1,1:this%nGeo+1,1:8))
    do step = 1,lvl
      call SubdivideNodeCoords(geomInterp,this%nGeo,cur,kids)
      cur(1:3,:,:,:) = kids(1:3,:,:,:,path(step))
    enddo
    coords = cur

    deallocate(cur,kids,path)

  endsubroutine LeafCoords_OctreeMesh3D

  recursive subroutine FaceNeighbor_OctreeMesh3D(this,node,s,nbr,ns,nf)
    !! Find the equal-or-larger face neighbour of `node` across its local face s using the
    !! classic octree ascend/descend search. Returns:
    !!   nbr - neighbour node id (0 = physical domain boundary). It is either a LEAF at
    !!         any level <= level(node), or an INTERNAL node at exactly level(node)
    !!         (meaning the shared face is subdivided on the neighbour side, i.e. finer
    !!         neighbours exist).
    !!   ns  - the neighbour's local face that faces `node`.
    !!   nf  - the flip between the two shared faces (0..7, sideInfo(4) convention),
    !!         inherited from the base-mesh face where the search crosses a root boundary.
    !! With this, a 2:1 hanging face is exactly "nbr is a leaf with level(nbr) =
    !! level(node)-1", and finer neighbours are exactly "nbr is internal".
    implicit none
    class(OctreeMesh3D),intent(in) :: this
    integer,intent(in) :: node,s
    integer,intent(out) :: nbr,ns,nf
    ! Local
    integer :: p,c,pnbr,ps,pf,t,tq

    if(this%level(node) == 0) then
      ! Cross a base-mesh face: neighbour root, its face, and the base flip.
      nbr = this%rootNbr(s,this%rootElem(node))
      ns = this%rootNbrSide(s,this%rootElem(node))
      nf = this%rootFlip(s,this%rootElem(node))
      return
    endif

    p = this%parent(node)
    c = this%octant(node)

    if(oc_internal(s,c)) then
      ! Neighbour is the sibling on the other side of an interior face of the parent.
      nbr = this%child(oc_reflect(s,c),p)
      ns = oc_opposite(s)
      nf = 0
      return
    endif

    ! Otherwise ascend: find the parent's neighbour across the same face, then descend.
    call FaceNeighbor_OctreeMesh3D(this,p,s,pnbr,ps,pf)
    if(pnbr == 0) then
      nbr = 0
      ns = 0
      nf = 0
      return
    endif

    if(this%child(1,pnbr) == 0) then
      ! Parent's neighbour is a leaf (equal or larger than the parent) -> our larger
      ! neighbour.
      nbr = pnbr
      ns = ps
      nf = pf
      return
    endif

    ! Parent's neighbour is internal (at level(node)-1): descend one level to the child
    ! that borders `node`, matching face quadrants across the face through the flip.
    t = oc_subpos(c,s)
    tq = faceQuadPerm(t,pf)
    nbr = this%child(childOfFace(tq,ps),pnbr)
    ns = ps
    nf = pf

  endsubroutine FaceNeighbor_OctreeMesh3D

  subroutine Balance2to1_OctreeMesh3D(this)
    !! Enforce the 2:1 face balance condition: no leaf face may separate elements
    !! differing by more than one refinement level. Iterates to a fixed point - in each
    !! sweep, any leaf whose equal-or-larger neighbour is a leaf two or more levels
    !! coarser triggers refinement of that coarser neighbour; refinement can ripple, so
    !! sweeps repeat until nothing changes. The leaf set is rebuilt on return.
    implicit none
    class(OctreeMesh3D),intent(inout) :: this
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
        do s = 1,6
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

  endsubroutine Balance2to1_OctreeMesh3D

  function MaxLevelJump_OctreeMesh3D(this) result(mx)
    !! Largest refinement-level difference across any leaf face (0 on a conforming or
    !! uniformly refined forest, 1 on a 2:1-balanced adaptive forest). Because
    !! FaceNeighbor returns the equal-or-larger neighbour, every level difference is
    !! observed from the finer leaf as a coarser leaf neighbour; internal (finer)
    !! neighbours contribute nothing from this side.
    implicit none
    class(OctreeMesh3D),intent(in) :: this
    integer :: mx
    ! Local
    integer :: li,s,node,nbr,ns,nf

    mx = 0
    do li = 1,this%nLeaves
      node = this%leaf(li)
      do s = 1,6
        call this%FaceNeighbor(node,s,nbr,ns,nf)
        if(nbr /= 0) then
          if(this%child(1,nbr) == 0) mx = max(mx,this%level(node)-this%level(nbr))
        endif
      enddo
    enddo

  endfunction MaxLevelJump_OctreeMesh3D

  ! -------------------------------------------------------------------------------------------- !
  ! Octree face-adjacency helpers (SELF child ordering = CGNS corner order; faces 1=Bottom,
  ! 2=South, 3=East, 4=North, 5=West, 6=Top). See SELF_RefinementPrimitives_3D for the
  ! octant/face conventions, childOfFace, and faceQuadPerm.
  ! -------------------------------------------------------------------------------------------- !

  pure function oc_opposite(s) result(o)
    !! The local face directly across an element from face s.
    implicit none
    integer,intent(in) :: s
    integer :: o
    integer,parameter :: opp(1:6) = [6,4,5,2,3,1]
    o = opp(s)
  endfunction oc_opposite

  pure function oc_internal(s,c) result(isInternal)
    !! .true. if child c's face s is interior to its parent (its neighbour across s is a
    !! sibling): the child sits in the half-cube away from parent face s.
    implicit none
    integer,intent(in) :: s,c
    logical :: isInternal
    select case(s)
    case(1) ! Bottom (xi3 = -1) interior for the top-half children
      isInternal = (octAzc(c) == 1)
    case(2) ! South (xi2 = -1) interior for the north-half children
      isInternal = (octAyc(c) == 1)
    case(3) ! East (xi1 = +1) interior for the west-half children
      isInternal = (octAxc(c) == 0)
    case(4) ! North (xi2 = +1) interior for the south-half children
      isInternal = (octAyc(c) == 0)
    case(5) ! West (xi1 = -1) interior for the east-half children
      isInternal = (octAxc(c) == 1)
    case default ! Top (xi3 = +1) interior for the bottom-half children
      isInternal = (octAzc(c) == 0)
    endselect
  endfunction oc_internal

  pure function oc_reflect(s,c) result(rc)
    !! Sibling child index obtained by reflecting c across face s (swap the half-index of
    !! the direction normal to s).
    implicit none
    integer,intent(in) :: s,c
    integer :: rc
    integer,parameter :: xref(1:8) = [2,1,4,3,6,5,8,7] ! swap ax (s = 3,5)
    integer,parameter :: yref(1:8) = [4,3,2,1,8,7,6,5] ! swap ay (s = 2,4)
    integer,parameter :: zref(1:8) = [5,6,7,8,1,2,3,4] ! swap az (s = 1,6)
    if(s == 1 .or. s == 6) then
      rc = zref(c)
    elseif(s == 2 .or. s == 4) then
      rc = yref(c)
    else
      rc = xref(c)
    endif
  endfunction oc_reflect

  pure function oc_subpos(c,s) result(t)
    !! Face quadrant (1..4, in face s's trace coordinates) of child c on face s;
    !! 0 if c does not touch face s.
    implicit none
    integer,intent(in) :: c,s
    integer :: t
    integer :: q
    t = 0
    do q = 1,4
      if(childOfFace(q,s) == c) then
        t = q
        return
      endif
    enddo
  endfunction oc_subpos

endmodule SELF_OctreeMesh_3D
