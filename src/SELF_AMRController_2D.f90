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

module SELF_AMRController_2D
!! Adaptive-mesh-refinement controller for 2-D DG models: the driver that closes the serial
!! AMR loop around a live, time-stepping model. One call to Adapt performs one adaptation
!! epoch between time steps:
!!
!!   1. estimate  - the Legendre modal-decay indicator flags each element refine/keep/coarsen
!!                  from the current solution (driving variable ivar, e.g. pressure);
!!   2. cap       - refine flags on leaves already at maxLevel are demoted to keep, bounding
!!                  the finest resolution (and the stable time step) a priori;
!!   3. halo      - refine flags spread to face neighbours for nHalo passes, so a feature
!!                  moving at speed c stays inside the refined band provided the adaptation
!!                  cadence satisfies k*dt*c <= nHalo * h_fine;
!!   4. mutate    - AdaptFromFlags + Balance2to1 update the quad-forest; if the leaf set is
!!                  unchanged the epoch is a no-op and the model is untouched;
!!   5. transfer  - BuildTransferPlan maps old leaves to new (copy / exact prolongation /
!!                  conservative restriction), EmitMesh produces the solver-ready
!!                  nonconforming mesh, and a new SEMQuad geometry is generated;
!!   6. regrid    - model%Regrid rebinds the model's storage and boundary conditions to the
!!                  new mesh, preserving its time state and parameters; the transferred
!!                  solution is applied and uploaded to the device.
!!
!! The controller owns the forest, the indicator, and the meshes/geometries it emits
!! (double-buffered: the previous pair is freed after the model is rebound to the new pair).
!! The base mesh and geometry the model was initialized with belong to the caller and are
!! never freed here, but the base mesh must outlive the controller: it supplies the
!! boundary-condition metadata and the communicator for every emitted mesh. After
!! controller%Free the model's mesh/geometry pointers are dangling, so free or stop using the
!! model first.
!!
!! MPI (AMR Stage 5): the forest is rank-replicated. At Init the global base-mesh tables are
!! allgathered so every rank builds an identical forest; each epoch the rank-local indicator
!! flags are allgathered (one small collective) so every rank applies identical mutations and
!! computes identical transfer plans and global connectivity. EmitMesh re-decomposes the new
!! leaf list into contiguous (space-filling-curve) ranges, so repartitioning and load balance
!! are implicit in every epoch. Solution migration gathers the old rank-local solutions into a
!! global field and applies the plan to the new rank-local range only - simple and correct at
!! single-node scale; a point-to-point exchange is a drop-in replacement behind the same
!! interface. All collectives run between time steps at the adaptation cadence; nothing is
!! added to the time-stepping loop.
!!
!! Because refinement halves the element scale per level, an explicit-stability time step
!! chosen for the base mesh must shrink with the finest active level;
!! RecommendedTimeStep(dtBase) = dtBase / 2**MaxLevel gives the level-based bound to pass to
!! ForwardStep after each epoch.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use SELF_DGModel2D_t
  use SELF_QuadTreeMesh_2D
  use SELF_AdaptiveMesh_2D
  use SELF_RefinementIndicator_2D
  use SELF_TransferPlan_2D
  use SELF_DomainDecomposition
  use mpi

  implicit none

  type :: AMRController2D
    type(QuadTreeMesh2D) :: forest
    type(RefinementIndicator2D) :: indicator
    type(Mesh2D),pointer :: baseMesh => null() !! caller-owned; metadata source for EmitMesh
    type(Mesh2D),pointer :: activeMesh => null() !! mesh the model currently runs on
    type(SEMQuad),pointer :: activeGeom => null() !! geometry the model currently runs on
    type(Lagrange),pointer :: interp => null() !! the model's solution interpolant
    logical :: ownsActive = .false. !! whether activeMesh/activeGeom were emitted (vs caller's)
    integer :: ivar = SELF_AMR_ALLVARS !! driving variable for the indicator
    integer :: maxLevel = 1 !! refinement-level cap
    integer :: nHalo = 1 !! refine-flag halo-expansion passes
    real(prec) :: refineThreshold = 0.0_prec
    real(prec) :: coarsenThreshold = 0.0_prec

  contains
    procedure,public :: Init => Init_AMRController2D
    procedure,public :: Free => Free_AMRController2D
    procedure,public :: Adapt => Adapt_AMRController2D
    procedure,public :: RecommendedTimeStep => RecommendedTimeStep_AMRController2D

  endtype AMRController2D

contains

  subroutine Init_AMRController2D(this,model,refineThreshold,coarsenThreshold,ivar, &
                                  maxLevel,nHalo)
    !! Attach the controller to an initialized model. The model's current mesh becomes the
    !! forest's base mesh (level 0); its geometry interpolant drives the indicator and the
    !! solution transfer. Thresholds are the sigma = log10 modal-energy-ratio cut-offs of the
    !! refinement indicator (refineThreshold > coarsenThreshold; see
    !! SELF_RefinementIndicator_2D). ivar is the driving solution variable (or
    !! SELF_AMR_ALLVARS). maxLevel >= 0 caps the refinement depth; nHalo >= 0 sets the
    !! refine-flag halo width in elements.
    implicit none
    class(AMRController2D),intent(out) :: this
    class(DGModel2D_t),intent(in) :: model
    real(prec),intent(in) :: refineThreshold
    real(prec),intent(in) :: coarsenThreshold
    integer,intent(in) :: ivar
    integer,intent(in) :: maxLevel
    integer,intent(in) :: nHalo

    if(.not. associated(model%mesh)) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : AMRController2D%Init requires an initialized model.'
      stop 1
    endif
    if(maxLevel < 0 .or. nHalo < 0) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : AMRController2D%Init requires maxLevel >= 0 and nHalo >= 0.'
      stop 1
    endif

    this%baseMesh => model%mesh
    this%activeMesh => model%mesh
    this%activeGeom => model%geometry
    this%interp => model%geometry%x%interp
    this%ownsActive = .false.
    this%ivar = ivar
    this%maxLevel = maxLevel
    this%nHalo = nHalo
    this%refineThreshold = refineThreshold
    this%coarsenThreshold = coarsenThreshold

    ! Rank-replicated forest: on one rank, straight from the mesh; on several, from the
    ! allgathered global base tables (every rank builds the identical forest).
    if(model%mesh%decomp%nRanks > 1) then
      call InitForestFromDecomposedMesh(this%forest,model%mesh)
    else
      call this%forest%Init(model%mesh)
    endif
    call this%indicator%Init(this%interp,model%mesh%nElem,refineThreshold,coarsenThreshold)

  endsubroutine Init_AMRController2D

  subroutine InitForestFromDecomposedMesh(forest,mesh)
    !! Build a rank-replicated forest over a decomposed base mesh: allgather the global
    !! node coordinates, side table, and material ids (by the decomposition's contiguous
    !! element ownership) and initialize the forest from the global tables. Collective over
    !! the mesh communicator; runs once, at controller initialization.
    implicit none
    type(QuadTreeMesh2D),intent(out) :: forest
    type(Mesh2D),intent(in) :: mesh
    ! Local
    integer :: nG,nGeo,r,s
    real(prec),allocatable :: coordsG(:,:,:,:)
    integer,allocatable :: siG(:,:,:),matG(:)
    integer,allocatable :: nbr(:,:),nbrSide(:,:),flip(:,:),bc(:,:)

    nG = mesh%decomp%nElem ! global element count
    nGeo = mesh%nGeo

    allocate(coordsG(1:2,1:nGeo+1,1:nGeo+1,1:nG))
    allocate(siG(1:5,1:4,1:nG))
    allocate(matG(1:nG))
    call AllgatherPerElemReals(mesh%decomp,2*(nGeo+1)*(nGeo+1), &
                               mesh%nodeCoords,coordsG)
    call AllgatherPerElemInts(mesh%decomp,20,mesh%sideInfo,siG)
    call AllgatherPerElemInts(mesh%decomp,1,mesh%elemMaterial,matG)

    ! Decode the global side table (sideInfo(3) already carries global element ids).
    allocate(nbr(1:4,1:nG),nbrSide(1:4,1:nG),flip(1:4,1:nG),bc(1:4,1:nG))
    do r = 1,nG
      do s = 1,4
        nbr(s,r) = siG(3,s,r)
        nbrSide(s,r) = siG(4,s,r)/10
        flip(s,r) = mod(siG(4,s,r),10)
        bc(s,r) = siG(5,s,r)
      enddo
    enddo

    call forest%InitGlobal(nG,nGeo,mesh%quadrature,coordsG,nbr,nbrSide,flip,bc,matG)

    deallocate(coordsG,siG,matG,nbr,nbrSide,flip,bc)

  endsubroutine InitForestFromDecomposedMesh

  subroutine AllgatherPerElemInts(decomp,perElem,localArr,globalArr)
    !! Allgather an integer array with perElem entries per element from the decomposition's
    !! contiguous rank-local element ranges into the global element ordering.
    implicit none
    type(DomainDecomposition),intent(in) :: decomp
    integer,intent(in) :: perElem
    integer,intent(in) :: localArr(*)
    integer,intent(out) :: globalArr(*)
    ! Local
    integer :: r,ierror
    integer,allocatable :: counts(:),displs(:)

    allocate(counts(1:decomp%nRanks),displs(1:decomp%nRanks))
    do r = 1,decomp%nRanks
      counts(r) = perElem*(decomp%offsetElem(r+1)-decomp%offsetElem(r))
      displs(r) = perElem*decomp%offsetElem(r)
    enddo
    call mpi_allgatherv(localArr,counts(decomp%rankId+1),MPI_INTEGER, &
                        globalArr,counts,displs,MPI_INTEGER, &
                        decomp%mpiComm,ierror)
    deallocate(counts,displs)

  endsubroutine AllgatherPerElemInts

  subroutine AllgatherPerElemReals(decomp,perElem,localArr,globalArr)
    !! Allgather a real(prec) array with perElem entries per element from the decomposition's
    !! contiguous rank-local element ranges into the global element ordering.
    implicit none
    type(DomainDecomposition),intent(in) :: decomp
    integer,intent(in) :: perElem
    real(prec),intent(in) :: localArr(*)
    real(prec),intent(out) :: globalArr(*)
    ! Local
    integer :: r,ierror
    integer,allocatable :: counts(:),displs(:)

    allocate(counts(1:decomp%nRanks),displs(1:decomp%nRanks))
    do r = 1,decomp%nRanks
      counts(r) = perElem*(decomp%offsetElem(r+1)-decomp%offsetElem(r))
      displs(r) = perElem*decomp%offsetElem(r)
    enddo
    call mpi_allgatherv(localArr,counts(decomp%rankId+1),decomp%mpiPrec, &
                        globalArr,counts,displs,decomp%mpiPrec, &
                        decomp%mpiComm,ierror)
    deallocate(counts,displs)

  endsubroutine AllgatherPerElemReals

  subroutine Free_AMRController2D(this)
    !! Release the forest, the indicator, and any controller-emitted mesh/geometry. The
    !! caller-owned base mesh/geometry are untouched. A model still pointing at a
    !! controller-emitted mesh must not be used after this call.
    implicit none
    class(AMRController2D),intent(inout) :: this

    call this%forest%Free()
    call this%indicator%Free()
    if(this%ownsActive) then
      call this%activeGeom%Free()
      call this%activeMesh%Free()
      deallocate(this%activeGeom)
      deallocate(this%activeMesh)
    endif
    this%baseMesh => null()
    this%activeMesh => null()
    this%activeGeom => null()
    this%interp => null()
    this%ownsActive = .false.

  endsubroutine Free_AMRController2D

  subroutine Adapt_AMRController2D(this,model,adapted)
    !! Perform one adaptation epoch on the model (see the module documentation). On return,
    !! adapted reports whether the mesh changed; when it did, the model is already rebound to
    !! the new mesh with the solution transferred (conservatively), and the caller should
    !! re-evaluate its time step (RecommendedTimeStep) before the next ForwardStep. When the
    !! leaf set is unchanged the model is untouched.
    !!
    !! Where the transferred solution lives on return: on a GPU build the transfer is performed
    !! on the device (Stage 6a) and the result is left there, so solution%interior (the host
    !! mirror) is STALE afterwards. This matches the rest of the time loop, where the device is
    !! authoritative and a caller that wants host data calls solution%UpdateHost() first - as
    !! Write_DGModel2D_t does before writing a snapshot. Before the device transfer existed the
    !! mirror happened to be fresh here because the transfer ran on the host; do not rely on
    !! that. On CPU builds host and device are the same storage and the question does not arise.
    implicit none
    class(AMRController2D),intent(inout) :: this
    class(DGModel2D_t),intent(inout) :: model
    logical,intent(out) :: adapted
    ! Local
    integer :: li,s,pass,node,nbr,ns,nf,nOld,Np,changed,eFirst,eLast,iv
    integer,allocatable :: flag(:),spread(:)
    integer,allocatable :: oldLeaf(:)
    integer,allocatable :: leafIdx(:)
    type(TransferPlan2D) :: plan
    type(Mesh2D),pointer :: newMesh
    type(SEMQuad),pointer :: newGeom
    real(prec),allocatable :: uOld(:,:,:,:)

    adapted = .false.

    if(.not. associated(model%mesh,this%activeMesh)) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : the model is not running on this controller''s active mesh.'
      stop 1
    endif

    ! ---- 1. Indicator flags from the current solution ----
    ! The indicator is rank-local; the (replicated) forest needs the global per-leaf flags, so
    ! on nRanks > 1 they are allgathered by the active decomposition's element ranges. From
    ! here on every rank applies identical mutations to its identical forest copy.
    call this%indicator%Estimate(model%solution,this%ivar)
    nOld = this%forest%nLeaves
    allocate(flag(1:nOld))
    if(model%mesh%decomp%nRanks > 1) then
      if(model%mesh%decomp%nElem /= nOld) then
        print*,__FILE__,':',__LINE__, &
          ' : Error : the active decomposition does not span the forest leaf list.'
        stop 1
      endif
      call AllgatherPerElemInts(model%mesh%decomp,1,this%indicator%flag,flag)
    else
      flag(1:nOld) = this%indicator%flag(1:nOld)
    endif

    ! ---- 2. Cap refinement at maxLevel ----
    do li = 1,nOld
      if(flag(li) == SELF_AMR_REFINE .and. &
         this%forest%level(this%forest%leaf(li)) >= this%maxLevel) then
        flag(li) = SELF_AMR_KEEP
      endif
    enddo

    ! ---- 3. Halo expansion: spread refine flags to face neighbours ----
    ! A leaf neighbour of a refine-flagged leaf is also flagged (up to the level cap) so the
    ! refined band extends nHalo elements beyond where the indicator fires; internal (finer)
    ! neighbours are already refined and need nothing.
    allocate(leafIdx(1:this%forest%nNodes))
    leafIdx = 0
    do li = 1,nOld
      leafIdx(this%forest%leaf(li)) = li
    enddo
    allocate(spread(1:nOld))
    do pass = 1,this%nHalo
      spread(1:nOld) = flag(1:nOld)
      do li = 1,nOld
        if(flag(li) /= SELF_AMR_REFINE) cycle
        node = this%forest%leaf(li)
        do s = 1,4
          call this%forest%FaceNeighbor(node,s,nbr,ns,nf)
          if(nbr == 0) cycle ! physical boundary
          if(this%forest%child(1,nbr) /= 0) cycle ! finer neighbour, already refined
          if(this%forest%level(nbr) >= this%maxLevel) cycle ! at the cap
          spread(leafIdx(nbr)) = SELF_AMR_REFINE
        enddo
      enddo
      flag(1:nOld) = spread(1:nOld)
    enddo
    deallocate(spread,leafIdx)

    ! ---- 4. Mutate the forest; detect a no-op epoch ----
    allocate(oldLeaf(1:nOld))
    oldLeaf(1:nOld) = this%forest%leaf(1:nOld)

    call this%forest%AdaptFromFlags(flag)
    deallocate(flag)
    call this%forest%Balance2to1()

    changed = 1
    if(this%forest%nLeaves == nOld) then
      changed = 0
      do li = 1,nOld
        if(this%forest%leaf(li) /= oldLeaf(li)) changed = 1
      enddo
    endif
    if(changed == 0) then
      deallocate(oldLeaf)
      return
    endif

    ! ---- 5. Transfer plan, emitted mesh, and geometry ----
    call BuildTransferPlan(this%forest,nOld,oldLeaf,plan)
    deallocate(oldLeaf)

    allocate(newMesh)
    call EmitMesh(this%forest,this%baseMesh,newMesh)
    allocate(newGeom)
    call newGeom%Init(this%interp,newMesh%nElem)
    call newGeom%GenerateFromMesh(newMesh)

    ! ---- 6. Regrid the model and transfer (migrate) the solution ----
    ! The solution is staged before Regrid (which releases the storage it lives in) and
    ! transferred onto the new mesh afterwards. Both steps are type-bound and backend-specific:
    ! the portable implementation stages on the host and runs ApplyTransferPlanRange, while the
    ! GPU backend stages device-to-device and applies the plan in a kernel, so an adapting run
    ! on one GPU moves no solution data across the host link at all (Stage 6a).
    !
    ! On several ranks the old field must first be assembled globally, because each rank then
    ! fills exactly its new contiguous element range and elements that changed ranks are
    ! migrated by construction (Stage-5 v1 migration). That allgather is a host operation, so
    ! the multi-rank path stays on the portable host transfer: a device transfer only pays off
    ! there once migration is point-to-point (Stage-5 v2).
    eFirst = -1 ! set below, once newMesh's decomposition is known
    if(model%mesh%decomp%nRanks > 1) then
      Np = this%interp%N+1
      allocate(uOld(1:Np,1:Np,1:nOld,1:model%nvar))
      call model%solution%UpdateHost()
      do iv = 1,model%nvar
        call AllgatherPerElemReals(model%mesh%decomp,Np*Np, &
                                   model%solution%interior(:,:,:,iv),uOld(:,:,:,iv))
      enddo

      call model%Regrid(newMesh,newGeom)

      eFirst = newMesh%decomp%offsetElem(newMesh%decomp%rankId+1)+1
      eLast = newMesh%decomp%offsetElem(newMesh%decomp%rankId+2)
      call model%ApplyTransferPlan(plan,this%interp,eFirst,eLast,uOld)
      deallocate(uOld)
    else
      call model%StageSolutionForTransfer()

      call model%Regrid(newMesh,newGeom)

      eFirst = newMesh%decomp%offsetElem(newMesh%decomp%rankId+1)+1
      eLast = newMesh%decomp%offsetElem(newMesh%decomp%rankId+2)
      call model%ApplyTransferPlan(plan,this%interp,eFirst,eLast)
    endif
    call plan%Free()

    ! ---- 7. Retire the previous mesh/geometry and re-size the indicator ----
    if(this%ownsActive) then
      call this%activeGeom%Free()
      call this%activeMesh%Free()
      deallocate(this%activeGeom)
      deallocate(this%activeMesh)
    endif
    this%activeMesh => newMesh
    this%activeGeom => newGeom
    this%ownsActive = .true.

    call this%indicator%Free()
    call this%indicator%Init(this%interp,newMesh%nElem, &
                             this%refineThreshold,this%coarsenThreshold)

    adapted = .true.

  endsubroutine Adapt_AMRController2D

  function RecommendedTimeStep_AMRController2D(this,dtBase) result(dt)
    !! Level-based explicit-stability time step: refinement halves the element scale per
    !! level, so a time step dtBase that is stable on the base (level-0) mesh scales to
    !! dtBase / 2**MaxLevel on the current forest. Deterministic and exact for the quadtree
    !! (child elements are exact half-scale subdivisions); no geometry reduction is needed.
    implicit none
    class(AMRController2D),intent(in) :: this
    real(prec),intent(in) :: dtBase
    real(prec) :: dt

    dt = dtBase/real(2**this%forest%MaxLevel(),prec)

  endfunction RecommendedTimeStep_AMRController2D

endmodule SELF_AMRController_2D
