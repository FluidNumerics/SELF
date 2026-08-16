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

module SELF_AMRController_3D
!! Adaptive-mesh-refinement controller for 3-D DG models, the direct analogue of
!! SELF_AMRController_2D: the driver that closes the AMR loop around a live,
!! time-stepping model. One call to Adapt performs one adaptation epoch between time
!! steps:
!!
!!   1. estimate  - the Legendre modal-decay indicator flags each element refine/keep/coarsen
!!                  from the current solution (driving variable ivar, e.g. pressure);
!!   2. cap       - refine flags on leaves already at maxLevel are demoted to keep, bounding
!!                  the finest resolution (and the stable time step) a priori;
!!   3. halo      - refine flags spread to face neighbours for nHalo passes, so a feature
!!                  moving at speed c stays inside the refined band provided the adaptation
!!                  cadence satisfies k*dt*c <= nHalo * h_fine;
!!   4. mutate    - AdaptFromFlags + Balance2to1 update the octree forest; if the leaf set is
!!                  unchanged the epoch is a no-op and the model is untouched;
!!   5. transfer  - BuildTransferPlan maps old leaves to new (copy / exact prolongation /
!!                  conservative restriction), EmitMesh produces the solver-ready
!!                  nonconforming mesh, and a new SEMHex geometry is generated;
!!   6. regrid    - model%Regrid rebinds the model's storage and boundary conditions to the
!!                  new mesh, preserving its time state and parameters; the transferred
!!                  solution is applied and uploaded to the device.
!!
!! The controller owns the forest, the indicator, and the meshes/geometries it emits. The
!! mesh is rebuilt each epoch and the previous one freed after the model is rebound.
!! Geometry instead lives in two long-lived buffers that alternate, so it is resized rather
!! than reallocated and the previous epoch's geometry remains available while the new one is
!! built - geometry reuse matters even more in 3-D than 2-D because the curl-invariant
!! metric-term computation is substantially more expensive. The base mesh and geometry the
!! model was initialized with belong to the caller and are never freed here, but the base
!! mesh must outlive the controller: it supplies the boundary-condition metadata and the
!! communicator for every emitted mesh. After controller%Free the model's mesh/geometry
!! pointers are dangling, so free or stop using the model first.
!!
!! MPI: the forest is rank-replicated. At Init the global base-mesh tables are allgathered
!! so every rank builds an identical forest; each epoch the rank-local indicator flags are
!! allgathered (one small collective) so every rank applies identical mutations and computes
!! identical transfer plans and global connectivity. EmitMesh re-decomposes the new leaf
!! list into contiguous (space-filling-curve) ranges, so repartitioning and load balance are
!! implicit in every epoch. Solution migration gathers the old rank-local solutions into a
!! global field and applies the plan to the new rank-local range only - simple and correct
!! at single-node scale; a point-to-point exchange is a drop-in replacement behind the same
!! interface. All collectives run between time steps at the adaptation cadence; nothing is
!! added to the time-stepping loop.
!!
!! Because refinement halves the element scale per level, an explicit-stability time step
!! chosen for the base mesh must shrink with the finest active level;
!! RecommendedTimeStep(dtBase) = dtBase / 2**MaxLevel gives the level-based bound to pass to
!! ForwardStep after each epoch.

  use iso_fortran_env,only:int64
  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_DGModel3D_t
  use SELF_OctreeMesh_3D
  use SELF_AdaptiveMesh_3D
  use SELF_RefinementIndicator_3D
  use SELF_TransferPlan_3D
  use SELF_DomainDecomposition
  use mpi

  implicit none

  !! Debug switches for the incremental geometry path, resolved once from the environment on
  !! first use and shared with the 2-D controller's conventions:
  !!
  !!   SELF_AMR_GEOM_NO_REUSE=1  regenerate every element (the reuse predicate never fires)
  !!   SELF_AMR_GEOM_VERIFY=1    additionally generate the full geometry the old way and compare
  !!                             element by element, reporting the first quantity that differs
  logical,save :: geomDebugResolved = .false.
  logical,save :: geomNoReuse = .false.
  logical,save :: geomVerify = .false.
  logical,save :: geomFull = .false. !! SELF_AMR_GEOM_FULL=1: bypass the incremental path

  type :: AMRController3D
    type(OctreeMesh3D) :: forest
    type(RefinementIndicator3D) :: indicator
    type(Mesh3D),pointer :: baseMesh => null() !! caller-owned; metadata source for EmitMesh
    type(Mesh3D),pointer :: activeMesh => null() !! mesh the model currently runs on
    type(SEMHex),pointer :: activeGeom => null() !! geometry the model currently runs on
    type(Lagrange),pointer :: interp => null() !! the model's solution interpolant
    logical :: ownsActive = .false. !! whether activeMesh was emitted by us (vs the caller's)
    !! Two long-lived geometry buffers, alternated each epoch. Geometry is resized in place
    !! rather than allocated and freed per epoch; alternating means the PREVIOUS epoch's
    !! geometry is still intact while the new one is filled, which the incremental reuse path
    !! needs. geomSlot records which buffer activeGeom currently is, or 0 while it is still the
    !! caller's geometry.
    type(SEMHex),pointer :: geomA => null()
    type(SEMHex),pointer :: geomB => null()
    integer :: geomSlot = 0
    !! Scratch geometry holding just the elements an epoch actually changed, generated compacted
    !! and then scattered into place. Persistent and resized, so it allocates only when an epoch
    !! changes more elements than any epoch before it.
    type(SEMHex),pointer :: genGeom => null()
    !! Cumulative geometry-reuse accounting, so a run can report what fraction of elements
    !! avoided regeneration rather than leaving the payoff to be estimated.
    integer(int64) :: nGeomReused = 0
    integer(int64) :: nGeomGenerated = 0
    integer :: ivar = SELF_AMR_ALLVARS !! driving variable for the indicator
    integer :: maxLevel = 1 !! refinement-level cap
    integer :: nHalo = 1 !! refine-flag halo-expansion passes
    real(prec) :: refineThreshold = 0.0_prec
    real(prec) :: coarsenThreshold = 0.0_prec
    ! Amplitude-gate settings forwarded to the indicator. They are held here because the
    ! indicator is Freed and re-Init-ed (intent(out)) at the end of every adapting epoch, which
    ! resets its own copies; the controller re-applies these afterwards.
    real(prec) :: relativeEnergyFloor = SELF_AMR_DEFAULT_RELFLOOR
    real(prec) :: significantEnergyFloor = SELF_AMR_DEFAULT_RELFLOOR
    real(prec),allocatable :: energyWeights(:)

  contains
    procedure,public :: Init => Init_AMRController3D
    procedure,public :: Free => Free_AMRController3D
    procedure,public :: Adapt => Adapt_AMRController3D
    procedure,public :: RecommendedTimeStep => RecommendedTimeStep_AMRController3D
    procedure,private :: ApplyIndicatorSettings => ApplyIndicatorSettings_AMRController3D
    procedure,private :: NextGeomBuffer => NextGeomBuffer_AMRController3D
    procedure,private :: BuildGeometry => BuildGeometry_AMRController3D

  endtype AMRController3D

contains

  subroutine Init_AMRController3D(this,model,refineThreshold,coarsenThreshold,ivar, &
                                  maxLevel,nHalo,relativeEnergyFloor,energyWeights, &
                                  significantEnergyFloor)
    !! Attach the controller to an initialized model. The model's current mesh becomes the
    !! forest's base mesh (level 0); its geometry interpolant drives the indicator and the
    !! solution transfer. Thresholds are the sigma = log10 modal-energy-ratio cut-offs of the
    !! refinement indicator (refineThreshold > coarsenThreshold; see
    !! SELF_RefinementIndicator_3D). ivar is the driving solution variable (or
    !! SELF_AMR_ALLVARS). maxLevel >= 0 caps the refinement depth; nHalo >= 0 sets the
    !! refine-flag halo width in elements.
    !!
    !! relativeEnergyFloor tunes the indicator's amplitude gate: an element whose energy is
    !! below this fraction of the field's peak element energy is treated as quiescent (hence
    !! resolved) regardless of its modal shape, which is what allows the mesh to be released
    !! behind a passing front. It is an energy fraction, so the square of the corresponding
    !! amplitude fraction; 0 disables the gate. energyWeights optionally weights the
    !! per-variable energies that form the gate, e.g. from the model's entropy function - see
    !! RefinementIndicator3D%SetEnergyWeights. significantEnergyFloor opens a hysteresis band
    !! on the energy axis: elements between the two floors hold their mesh (SELF_AMR_KEEP)
    !! instead of flipping across a single hard cut. All three default to the indicator's own
    !! defaults.
    implicit none
    class(AMRController3D),intent(out) :: this
    class(DGModel3D_t),intent(in) :: model
    real(prec),intent(in) :: refineThreshold
    real(prec),intent(in) :: coarsenThreshold
    integer,intent(in) :: ivar
    integer,intent(in) :: maxLevel
    integer,intent(in) :: nHalo
    real(prec),intent(in),optional :: relativeEnergyFloor
    real(prec),intent(in),optional :: energyWeights(:)
    real(prec),intent(in),optional :: significantEnergyFloor

    if(.not. associated(model%mesh)) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : AMRController3D%Init requires an initialized model.'
      stop 1
    endif
    if(maxLevel < 0 .or. nHalo < 0) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : AMRController3D%Init requires maxLevel >= 0 and nHalo >= 0.'
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
    if(present(relativeEnergyFloor)) then
      this%relativeEnergyFloor = relativeEnergyFloor
      this%significantEnergyFloor = relativeEnergyFloor
    endif
    if(present(significantEnergyFloor)) this%significantEnergyFloor = significantEnergyFloor
    if(present(energyWeights)) then
      allocate(this%energyWeights(1:size(energyWeights)))
      this%energyWeights = energyWeights
    endif

    ! Rank-replicated forest: on one rank, straight from the mesh; on several, from the
    ! allgathered global base tables (every rank builds the identical forest).
    if(model%mesh%decomp%nRanks > 1) then
      call InitForestFromDecomposedMesh(this%forest,model%mesh)
    else
      call this%forest%Init(model%mesh)
    endif
    call this%indicator%Init(this%interp,model%mesh%nElem,refineThreshold,coarsenThreshold)
    call this%ApplyIndicatorSettings()

  endsubroutine Init_AMRController3D

  subroutine ApplyIndicatorSettings_AMRController3D(this)
    !! Push the controller-owned amplitude-gate settings onto the indicator. Called after every
    !! indicator Init, because Init is intent(out) and therefore resets them to the indicator's
    !! defaults. The setters carry the validation.
    implicit none
    class(AMRController3D),intent(inout) :: this

    call this%indicator%SetRelativeEnergyFloor(this%relativeEnergyFloor, &
                                               significantEnergyFloor=this%significantEnergyFloor)
    if(allocated(this%energyWeights)) then
      call this%indicator%SetEnergyWeights(this%energyWeights)
    endif

  endsubroutine ApplyIndicatorSettings_AMRController3D

  subroutine InitForestFromDecomposedMesh(forest,mesh)
    !! Build a rank-replicated forest over a decomposed base mesh: allgather the global
    !! node coordinates, side table, and material ids (by the decomposition's contiguous
    !! element ownership) and initialize the forest from the global tables. Collective over
    !! the mesh communicator; runs once, at controller initialization.
    implicit none
    type(OctreeMesh3D),intent(out) :: forest
    type(Mesh3D),intent(in) :: mesh
    ! Local
    integer :: nG,nGeo,r,s
    real(prec),allocatable :: coordsG(:,:,:,:,:)
    integer,allocatable :: siG(:,:,:),matG(:)
    integer,allocatable :: nbr(:,:),nbrSide(:,:),flip(:,:),bc(:,:)

    nG = mesh%decomp%nElem ! global element count
    nGeo = mesh%nGeo

    allocate(coordsG(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nG))
    allocate(siG(1:5,1:6,1:nG))
    allocate(matG(1:nG))
    call AllgatherPerElemReals(mesh%decomp,3*(nGeo+1)*(nGeo+1)*(nGeo+1), &
                               mesh%nodeCoords,coordsG)
    call AllgatherPerElemInts(mesh%decomp,30,mesh%sideInfo,siG)
    call AllgatherPerElemInts(mesh%decomp,1,mesh%elemMaterial,matG)

    ! Decode the global side table (sideInfo(3) already carries global element ids).
    allocate(nbr(1:6,1:nG),nbrSide(1:6,1:nG),flip(1:6,1:nG),bc(1:6,1:nG))
    do r = 1,nG
      do s = 1,6
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

  subroutine BuildGeometry_AMRController3D(this,newMesh,plan,newGeom,nReused)
    !! Fill newGeom for the emitted mesh, reusing the previous epoch's geometry for every
    !! element that did not change and generating only the rest. nReused reports how many
    !! elements were copied rather than computed.
    !!
    !! Why the reuse is exact. The transfer plan marks a new leaf SELF_TRANSFER_COPY only when
    !! the walk depth is zero, i.e. when the new element and old element sourceElem(li) are the
    !! SAME forest node (see BuildTransferPlan). A leaf's mesh node coordinates come from
    !! LeafCoords, which is a pure function of the root element's coordinates, the leaf's level
    !! and its octant path, evaluated in a fixed order; root coordinates are never mutated and
    !! node ids, levels and octants are stable across forest mutations. So an unchanged leaf's
    !! coordinates are bit-identical between epochs, and because per-element geometry
    !! generation touches only that element's own coordinates, its whole geometry block is too.
    !!
    !! Multi-rank: sourceElem indexes the GLOBAL old element list while each rank holds only
    !! its own slice, so reuse additionally requires the source to be locally owned. That is a
    !! range test against the OLD decomposition, which is still valid here because Regrid has
    !! not run yet. On one rank every COPY element qualifies and the test is always true; on
    !! several ranks whatever migrated is simply regenerated. No communication is added either
    !! way.
    implicit none
    class(AMRController3D),intent(inout) :: this
    type(Mesh3D),intent(in) :: newMesh
    type(TransferPlan3D),intent(in) :: plan
    type(SEMHex),intent(inout) :: newGeom
    integer,intent(out) :: nReused
    ! Local
    integer :: li,gi,src,nLocal,nGen,oldFirst,oldLast,eFirst,nGeo,k
    integer,allocatable :: srcIdx(:),dstIdx(:),genIdx(:)
    real(prec),allocatable :: genCoords(:,:,:,:,:)

    nLocal = newMesh%nElem
    nGeo = newMesh%nGeo
    eFirst = newMesh%decomp%offsetElem(newMesh%decomp%rankId+1)+1

    ! Rank-local range of the OLD element list, i.e. what activeGeom actually holds.
    oldFirst = this%activeMesh%decomp%offsetElem(this%activeMesh%decomp%rankId+1)+1
    oldLast = this%activeMesh%decomp%offsetElem(this%activeMesh%decomp%rankId+2)

    allocate(srcIdx(1:nLocal),dstIdx(1:nLocal),genIdx(1:nLocal))

    call ResolveGeomDebug()

    ! Diagnostic bypass: full regeneration on the target buffer while keeping the persistent
    ! alternating buffers, to tell a defect in the incremental assembly apart from one in the
    ! buffer reuse itself.
    if(geomFull) then
      call newGeom%GenerateFromMesh(newMesh)
      nReused = 0
      deallocate(srcIdx,dstIdx,genIdx)
      return
    endif

    nReused = 0
    nGen = 0
    do li = 1,nLocal
      gi = eFirst+li-1 ! this element's index in the plan's global new-leaf arrays
      src = plan%sourceElem(gi)
      if(.not. geomNoReuse .and. plan%sourceKind(gi) == SELF_TRANSFER_COPY .and. &
         src >= oldFirst .and. src <= oldLast) then
        nReused = nReused+1
        srcIdx(nReused) = src-oldFirst+1 ! rank-local index into activeGeom
        dstIdx(nReused) = li
      else
        nGen = nGen+1
        genIdx(nGen) = li
      endif
    enddo

    ! Generate the changed elements, compacted, so the generation loops run over nGen elements
    ! instead of all of them. Their geometry is then scattered into place.
    if(nGen > 0) then
      allocate(genCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nGen))
      do k = 1,nGen
        genCoords(1:3,:,:,:,k) = newMesh%nodeCoords(1:3,:,:,:,genIdx(k))
      enddo

      if(.not. associated(this%genGeom)) allocate(this%genGeom)
      if(this%genGeom%nElem == 0) then
        call this%genGeom%Init(this%interp,nGen)
      else
        call this%genGeom%Resize(this%interp,nGen)
      endif
      call this%genGeom%GenerateFromNodeCoords(genCoords,nGeo,newMesh%quadrature,nGen)
      call this%genGeom%x%UpdateDevice()
      call this%genGeom%x%BoundaryInterp()
      call this%genGeom%x%UpdateHost()
      call this%genGeom%CalculateMetricTerms()

      do k = 1,nGen
        srcIdx(nReused+k) = k
        dstIdx(nReused+k) = genIdx(k)
      enddo
      call newGeom%CopyElements(this%genGeom,srcIdx(nReused+1:),dstIdx(nReused+1:),nGen)
      deallocate(genCoords)
    endif

    ! Carry the unchanged elements across from the previous epoch's geometry.
    if(nReused > 0) then
      call newGeom%CopyElements(this%activeGeom,srcIdx,dstIdx,nReused)
    endif

    call newGeom%UploadGeometry()

    if(geomVerify) call VerifyGeometry(this,newMesh,newGeom)

    deallocate(srcIdx,dstIdx,genIdx)

  endsubroutine BuildGeometry_AMRController3D

  subroutine ResolveGeomDebug()
    !! Read the incremental-geometry debug switches once.
    implicit none
    character(8) :: envstr
    integer :: envstat

    if(geomDebugResolved) return
    geomDebugResolved = .true.

    call get_environment_variable("SELF_AMR_GEOM_NO_REUSE",envstr,status=envstat)
    if(envstat == 0) geomNoReuse = (trim(envstr) == "1")
    call get_environment_variable("SELF_AMR_GEOM_VERIFY",envstr,status=envstat)
    if(envstat == 0) geomVerify = (trim(envstr) == "1")
    call get_environment_variable("SELF_AMR_GEOM_FULL",envstr,status=envstat)
    if(envstat == 0) geomFull = (trim(envstr) == "1")

    if(geomNoReuse) print*,"SELF_AMR_GEOM_NO_REUSE: geometry reuse disabled"
    if(geomVerify) print*,"SELF_AMR_GEOM_VERIFY: cross-checking incremental geometry"
    if(geomFull) print*,"SELF_AMR_GEOM_FULL: incremental geometry bypassed"

  endsubroutine ResolveGeomDebug

  subroutine VerifyGeometry(this,newMesh,newGeom)
    !! Generate the geometry the original way and report, per quantity, the largest
    !! discrepancy against the incrementally built one. Diagnostic only; gated by
    !! SELF_AMR_GEOM_VERIFY.
    !!
    !! Whole-array maxima rather than per-element reporting: naming the quantity that diverges
    !! is what localizes the defect, and it keeps this to F2008 array intrinsics.
    implicit none
    class(AMRController3D),intent(inout) :: this
    type(Mesh3D),intent(in) :: newMesh
    type(SEMHex),intent(in) :: newGeom
    ! Local
    type(SEMHex) :: ref

    call ref%Init(this%interp,newMesh%nElem)
    call ref%GenerateFromMesh(newMesh)

    print*,"GEOM_VERIFY nElem            =",newMesh%nElem
    print*,"GEOM_VERIFY max|d x%interior|     =", &
      maxval(abs(newGeom%x%interior-ref%x%interior))
    print*,"GEOM_VERIFY max|d x%boundary|     =", &
      maxval(abs(newGeom%x%boundary-ref%x%boundary))
    print*,"GEOM_VERIFY max|d J%interior|     =", &
      maxval(abs(newGeom%J%interior-ref%J%interior))
    print*,"GEOM_VERIFY max|d J%boundary|     =", &
      maxval(abs(newGeom%J%boundary-ref%J%boundary))
    print*,"GEOM_VERIFY max|d dsdx%interior|  =", &
      maxval(abs(newGeom%dsdx%interior-ref%dsdx%interior))
    print*,"GEOM_VERIFY max|d dsdx%boundary|  =", &
      maxval(abs(newGeom%dsdx%boundary-ref%dsdx%boundary))
    print*,"GEOM_VERIFY max|d nHat%boundary|  =", &
      maxval(abs(newGeom%nHat%boundary-ref%nHat%boundary))
    print*,"GEOM_VERIFY max|d nScale%boundary|=", &
      maxval(abs(newGeom%nScale%boundary-ref%nScale%boundary))

    call ref%Free()

  endsubroutine VerifyGeometry

  subroutine NextGeomBuffer_AMRController3D(this,nElem,geom,slot)
    !! Select the geometry buffer to fill this epoch: whichever of the two is not currently
    !! active, so the previous epoch's geometry stays intact and readable.
    !!
    !! Each buffer is Init-ed the first time it is used and Resized on every subsequent epoch,
    !! so after the first two adaptations geometry performs no allocation at all - which is the
    !! point of the change. On the very first adaptation activeGeom is still the caller's
    !! geometry (geomSlot == 0); that object belongs to the caller and is never freed or reused
    !! here.
    implicit none
    class(AMRController3D),intent(inout) :: this
    integer,intent(in) :: nElem
    type(SEMHex),pointer,intent(out) :: geom
    integer,intent(out) :: slot

    if(this%geomSlot == 1) then
      slot = 2
    else
      slot = 1
    endif

    if(slot == 1) then
      if(.not. associated(this%geomA)) allocate(this%geomA)
      geom => this%geomA
    else
      if(.not. associated(this%geomB)) allocate(this%geomB)
      geom => this%geomB
    endif

    if(geom%nElem == 0) then
      call geom%Init(this%interp,nElem)
    else
      call geom%Resize(this%interp,nElem)
    endif

  endsubroutine NextGeomBuffer_AMRController3D

  subroutine Free_AMRController3D(this)
    !! Release the forest, the indicator, and any controller-emitted mesh/geometry. The
    !! caller-owned base mesh/geometry are untouched. A model still pointing at a
    !! controller-emitted mesh must not be used after this call.
    implicit none
    class(AMRController3D),intent(inout) :: this

    call this%forest%Free()
    call this%indicator%Free()
    if(this%ownsActive) then
      call this%activeMesh%Free()
      deallocate(this%activeMesh)
    endif
    if(associated(this%geomA)) then
      call this%geomA%Free()
      deallocate(this%geomA)
    endif
    if(associated(this%geomB)) then
      call this%geomB%Free()
      deallocate(this%geomB)
    endif
    if(associated(this%genGeom)) then
      call this%genGeom%Free()
      deallocate(this%genGeom)
    endif
    this%baseMesh => null()
    this%activeMesh => null()
    this%activeGeom => null()
    this%geomA => null()
    this%geomB => null()
    this%genGeom => null()
    this%geomSlot = 0
    this%nGeomReused = 0
    this%nGeomGenerated = 0
    this%interp => null()
    this%ownsActive = .false.
    if(allocated(this%energyWeights)) deallocate(this%energyWeights)

  endsubroutine Free_AMRController3D

  subroutine Adapt_AMRController3D(this,model,adapted)
    !! Perform one adaptation epoch on the model (see the module documentation). On return,
    !! adapted reports whether the mesh changed; when it did, the model is already rebound to
    !! the new mesh with the solution transferred (conservatively), and the caller should
    !! re-evaluate its time step (RecommendedTimeStep) before the next ForwardStep. When the
    !! leaf set is unchanged the model is untouched.
    !!
    !! Where the transferred solution lives on return: on a GPU build the transfer is performed
    !! on the device (the 3-D analogue of the 2-D Stage 6a) and the result is left there, so
    !! solution%interior (the host mirror) is STALE afterwards. This matches the rest of the
    !! time loop, where the device is authoritative and a caller that wants host data calls
    !! solution%UpdateHost() first - as Write_DGModel3D_t does before writing a snapshot. Before
    !! the device transfer existed the mirror happened to be fresh here because the transfer ran
    !! on the host; do not rely on that. On CPU builds host and device are the same storage and
    !! the question does not arise. The multi-rank path still transfers on the host (see step 6),
    !! so both mirrors are fresh there.
    implicit none
    class(AMRController3D),intent(inout) :: this
    class(DGModel3D_t),intent(inout) :: model
    logical,intent(out) :: adapted
    ! Local
    integer :: li,s,pass,node,nbr,ns,nf,nOld,Np,changed,eFirst,eLast,iv
    integer,allocatable :: flag(:),spread(:)
    integer,allocatable :: oldLeaf(:)
    integer,allocatable :: leafIdx(:)
    type(TransferPlan3D) :: plan
    type(Mesh3D),pointer :: newMesh
    type(SEMHex),pointer :: newGeom
    integer :: newSlot
    integer :: nReused
    real(prec),allocatable :: uOld(:,:,:,:,:)

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
    ! The indicator's amplitude gate normalizes each element's energy by the peak element
    ! energy, so on several ranks that peak must be global: otherwise the gate - and hence the
    ! flags and the adapted mesh - would depend on the decomposition. That is one extra small
    ! collective per epoch, alongside the flag allgather below, and none inside the
    ! time-stepping loop.
    if(model%mesh%decomp%mpiEnabled) then
      call this%indicator%Estimate(model%solution,this%ivar, &
                                   comm=model%mesh%decomp%mpiComm)
    else
      call this%indicator%Estimate(model%solution,this%ivar)
    endif
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
        do s = 1,6
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

    ! Take the geometry buffer that is NOT currently active, so the previous epoch's geometry
    ! stays readable. Each buffer is Init-ed once and resized thereafter.
    call this%NextGeomBuffer(newMesh%nElem,newGeom,newSlot)
    call this%BuildGeometry(newMesh,plan,newGeom,nReused)
    this%nGeomReused = this%nGeomReused+int(nReused,int64)
    this%nGeomGenerated = this%nGeomGenerated+int(newMesh%nElem-nReused,int64)

    ! ---- 6. Regrid the model and transfer (migrate) the solution ----
    ! The solution is staged before Regrid (which releases the storage it lives in) and
    ! transferred onto the new mesh afterwards. Both steps are type-bound and backend-specific:
    ! the portable implementation stages on the host and runs ApplyTransferPlanRange, while the
    ! GPU backend stages device-to-device and applies the plan in a kernel. That moves the
    ! per-element tensor-product interpolation onto the device - which is where the measured
    ! saving comes from - and incidentally leaves no solution data crossing the host link.
    !
    ! On several ranks the old field must first be assembled globally, because each rank then
    ! fills exactly its new contiguous element range and elements that changed ranks are
    ! migrated by construction (the 2-D Stage-5 v1 migration). That allgather is a host
    ! operation, so the multi-rank path stays on the portable host transfer: a device transfer
    ! only pays off there once migration is point-to-point (Stage-5 v2).
    eFirst = -1 ! set below, once newMesh's decomposition is known
    if(model%mesh%decomp%nRanks > 1) then
      Np = this%interp%N+1
      allocate(uOld(1:Np,1:Np,1:Np,1:nOld,1:model%nvar))
      call model%solution%UpdateHost()
      do iv = 1,model%nvar
        call AllgatherPerElemReals(model%mesh%decomp,Np*Np*Np, &
                                   model%solution%interior(:,:,:,:,iv),uOld(:,:,:,:,iv))
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
    ! The geometry buffers are retained and reused; only the mesh is still rebuilt per epoch.
    if(this%ownsActive) then
      call this%activeMesh%Free()
      deallocate(this%activeMesh)
    endif
    this%activeMesh => newMesh
    this%activeGeom => newGeom
    this%geomSlot = newSlot
    this%ownsActive = .true.

    call this%indicator%Free()
    call this%indicator%Init(this%interp,newMesh%nElem, &
                             this%refineThreshold,this%coarsenThreshold)
    ! Only controller-owned indicator settings survive an epoch. A driver that called
    ! indicator%SetEnergyScale (or SetEnergyWeights) directly must re-apply it after any epoch
    ! that reported adapted = .true.
    call this%ApplyIndicatorSettings()

    adapted = .true.

  endsubroutine Adapt_AMRController3D

  function RecommendedTimeStep_AMRController3D(this,dtBase) result(dt)
    !! Level-based explicit-stability time step: refinement halves the element scale per
    !! level, so a time step dtBase that is stable on the base (level-0) mesh scales to
    !! dtBase / 2**MaxLevel on the current forest. Deterministic and exact for the octree
    !! (child elements are exact half-scale subdivisions); no geometry reduction is needed.
    implicit none
    class(AMRController3D),intent(in) :: this
    real(prec),intent(in) :: dtBase
    real(prec) :: dt

    dt = dtBase/real(2**this%forest%MaxLevel(),prec)

  endfunction RecommendedTimeStep_AMRController3D

endmodule SELF_AMRController_3D
