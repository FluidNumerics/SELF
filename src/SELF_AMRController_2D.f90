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
!! boundary-condition and material metadata for every emitted mesh. After controller%Free the
!! model's mesh/geometry pointers are dangling, so free or stop using the model first.
!!
!! Serial only, matching the rest of the AMR stack (MPI repartitioning is AMR Stage 5).
!! Adapt runs between time steps at the caller's cadence; it is not a per-step hot path.
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
    if(model%mesh%decomp%nRanks > 1) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : AMRController2D is serial-only pending AMR Stage 5 (MPI repartitioning).'
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

    call this%forest%Init(model%mesh)
    call this%indicator%Init(this%interp,model%mesh%nElem,refineThreshold,coarsenThreshold)

  endsubroutine Init_AMRController2D

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
    !! the new mesh with the solution transferred (conservatively) and uploaded to the device,
    !! and the caller should re-evaluate its time step (RecommendedTimeStep) before the next
    !! ForwardStep. When the leaf set is unchanged the model is untouched.
    implicit none
    class(AMRController2D),intent(inout) :: this
    class(DGModel2D_t),intent(inout) :: model
    logical,intent(out) :: adapted
    ! Local
    integer :: li,s,pass,node,nbr,ns,nf,nOld,Np,changed
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
    call this%indicator%Estimate(model%solution,this%ivar)
    nOld = this%forest%nLeaves
    allocate(flag(1:nOld))
    flag(1:nOld) = this%indicator%flag(1:nOld)

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

    ! ---- 6. Regrid the model and transfer the solution ----
    Np = this%interp%N+1
    allocate(uOld(1:Np,1:Np,1:nOld,1:model%nvar))
    call model%solution%UpdateHost()
    uOld(1:Np,1:Np,1:nOld,1:model%nvar) = &
      model%solution%interior(1:Np,1:Np,1:nOld,1:model%nvar)

    call model%Regrid(newMesh,newGeom)

    call ApplyTransferPlan(plan,this%interp,model%nvar, &
                           uOld,model%solution%interior)
    call model%solution%UpdateDevice()
    deallocate(uOld)
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
