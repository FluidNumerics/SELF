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

module SELF_LocalTimeStepping_2D
!! Local time stepping (LTS) for the 2-D adaptive solver - AMR Stage 7.
!!
!! Under the level-based global time step, every element is advanced at the rate the FINEST
!! level requires: `dt = dtBase / 2**maxLevel` (AMRController2D%RecommendedTimeStep). The
!! coarse bulk, which is nearly all of the mesh in a typical wavefront-tracking run, is then
!! over-resolved in time by a factor of `2**maxLevel`. LTS removes that: elements at
!! refinement level ell are advanced at `dtBase / 2**ell`, so each element is stepped at the
!! rate its own size demands.
!!
!! ## Why the level interfaces are the whole problem
!!
!! EmitMesh only ever emits a conforming side between two leaves of the SAME level; every
!! level jump is a 2:1 mortar (SELF_AdaptiveMesh_2D). Therefore:
!!
!!    all coupling between different time-step sizes flows through mesh%mortarInfo.
!!
!! Within a level, SideExchange behaves exactly as it always has - the neighbours of a
!! level-ell element across conforming sides are themselves at level ell, and are being
!! advanced by the same substep. Only the mortars need new machinery.
!!
!! ## The scheme : Berger-Oliger recursion with conservative refluxing
!!
!! One macro step of size dt0 = dtBase is
!!
!!   StepLevel(0, dt0):  advance level 0 by dt0
!!                       StepLevel(1, dt0/2) ; StepLevel(1, dt0/2)
!!                       reflux the level-0/1 interfaces
!!
!! and StepLevel recurses to maxLevel. At the end of a macro step every level has reached
!! t0 + dt0, so entropy diagnostics, file IO and AMR adaptation see a globally synchronized
!! state and need no changes.
!!
!! Coupling in each direction:
!!
!!   coarse -> fine.  The coarse element has already finished its step when the fine
!!     substeps run, so its stored edge trace is at the wrong time. Before each fine RK
!!     stage the big side's trace is overwritten with a time interpolation of a short
!!     history (`bigTrace`), and the ordinary MortarExchange then restricts it onto the
!!     small sides exactly as before - including over MPI, since MortarExchange already
!!     ships the big trace to the small sides' ranks. Quadratic interpolation through
!!     t0-dtC, t0, t0+dtC is used once the history is deep enough (third order, matching
!!     RK3); linear at the start of a run and after each adaptation epoch.
!!
!!   fine -> coarse.  During the coarse stages the fine traces are still at t0 - the
!!     standard Berger-Oliger lag - so the coarse element temporarily sees a lagged
!!     interface flux. Both what it used and what the fine level actually used are
!!     accumulated into flux registers, and Reflux then replaces the former with the
!!     latter. Refluxing is what makes the interface conservative; it also restores the
!!     accuracy of the interface flux, so the lag does not degrade the coarse solution.
!!
!! ## Why the flux registers are exact
!!
!! For a Williamson 2N low-storage scheme (`S_m = a_m S_{m-1} + R_m`, `U_m = U_{m-1} +
!! g_m dt S_m`) a full step telescopes into a plain weighted sum of the stage tendencies,
!!
!!     U^{n+1} = U^n + dt ( c_1 R_1 + c_2 R_2 + c_3 R_3 )
!!     c_3 = g_3,  c_2 = g_2 + g_3 a_3,  c_1 = g_1 + g_2 a_2 + g_3 a_3 a_2
!!
!! (`rk3StageWeight` below; the weights sum to one, which is the consistency condition).
!! The DG surface term enters R linearly through flux%boundaryNormal, so accumulating
!! `c_m * dt * boundaryNormal` over the stages gives the exact time-integrated interface
!! flux of the step, and adding the surface term of a register DIFFERENCE afterwards
!! reproduces exactly the update that a different interface flux would have produced.
!!
!! One honest caveat: the coarse element's INTERMEDIATE stage values still saw the lagged
!! flux. For a linear system - LinearEuler2D and LinearShallowWater2D, the models AMR
!! targets - the correction is therefore exact. For a nonlinear system it is the standard
!! AMR refluxing approximation: conservation is exact either way, the stage states are not.
!!
!! ## Scope and guards (this stage)
!!
!!   - RK3 only. The recursion and the registers are written around the 2N weights above.
!!   - CPU only. The GPU operators flatten the element dimension into one degree-of-freedom
!!     count and have no element-subset kernel; passing a subset to them is a hard error
!!     (RejectSubset in SELF_Data), which is what stops a GPU build silently single-rating.
!!   - Hyperbolic only (gradient_enabled = .false.). Parabolic terms across a level
!!     interface need the gradient traces interpolated in time as well.
!!
!! See docs/Learning/AdaptiveMeshRefinement.md section 5.4 for how this fits the AMR staging.

  use SELF_Constants
  use SELF_Mesh_2D
  use SELF_Model
  use SELF_DGModel2D

  implicit none

  !! Register selectors for AccumulateRegister.
  integer,parameter :: LTS_REG_COARSE = 1 !! what the interface's own (big) level used
  integer,parameter :: LTS_REG_FINE = 2 !! what the level one finer used, projected

  !! Per-level index lists. All three are GLOBAL mortar indices or RANK-LOCAL element
  !! indices, matching what the subset-aware operators expect.
  type :: LTSLevel2D
    !! Rank-local elements at this level.
    integer,pointer,contiguous :: elems(:) => null()
    !! Mortars this level touches at all, as the big side or as a small side. This is the
    !! subset handed to MortarExchange while the level is being advanced.
    integer,pointer,contiguous :: mortars(:) => null()
    !! Mortars where this level is the BIG side (mortarLevel == this level). These are the
    !! interfaces this level must reflux after its children have caught up.
    integer,pointer,contiguous :: mortarsBig(:) => null()
    !! Mortars where this level is a SMALL side (mortarLevel == this level - 1). These are
    !! the interfaces whose big-side trace must be interpolated in time for this level, and
    !! whose flux this level contributes to the parent's register.
    integer,pointer,contiguous :: mortarsSmall(:) => null()
    !! Number of valid slots in the big-side trace history of this level's mortarsBig:
    !! 0 before the first step, 1 after seeding, 2 (linear), 3 (quadratic).
    integer :: historyDepth = 0
    !! Start time of the step whose end value sits in history slot 3, i.e. slot 2 was
    !! captured at windowStart and slot 3 at windowStart + this level's dt. A child level
    !! interpolates its parent's trace against this window.
    real(prec) :: windowStart = 0.0_prec
  endtype LTSLevel2D

  type :: LTSSchedule2D
    logical :: built = .false.
    integer :: maxLevel = 0
    integer :: nMortars = 0
    integer :: N = 0
    integer :: nVar = 0
    type(LTSLevel2D),allocatable :: level(:) ! 0:maxLevel

    !! Flux registers, in the BIG side's trace space, indexed by global mortar id.
    !! regCoarse is the time-integrated interface flux the big side actually used over its
    !! own step; regFine is the time-integrated flux the two fine substeps used, already
    !! projected onto the big side by MortarFluxCollect. Reflux applies the difference.
    real(prec),allocatable :: regCoarse(:,:,:) ! (N+1, nVar, nMortars)
    real(prec),allocatable :: regFine(:,:,:)
    !! Big-side edge trace history for time interpolation, slots (t-dt, t, t+dt) in
    !! 1,2,3 with 3 the most recent. Only the rank owning the big element fills its entry.
    real(prec),allocatable :: bigTrace(:,:,:,:) ! (N+1, nVar, 3, nMortars)

  contains
    procedure,public :: Build => Build_LTSSchedule2D
    procedure,public :: Free => Free_LTSSchedule2D
    procedure,public :: Matches => Matches_LTSSchedule2D
  endtype LTSSchedule2D

contains

  pure function rk3StageWeight(m) result(c)
    !! Effective quadrature weight of stage m in a full Williamson 2N RK3 step:
    !! U^{n+1} = U^n + dt * sum_m c_m R_m. Derived once from rk3_a and rk3_g rather than
    !! hard-coded, so the registers cannot drift away from the integrator's coefficients.
    implicit none
    integer,intent(in) :: m
    real(prec) :: c

    select case(m)
    case(1)
      c = rk3_g(1)+rk3_g(2)*rk3_a(2)+rk3_g(3)*rk3_a(3)*rk3_a(2)
    case(2)
      c = rk3_g(2)+rk3_g(3)*rk3_a(3)
    case default
      c = rk3_g(3)
    endselect

  endfunction rk3StageWeight

! ---------------------------------------------------------------------------------------- !
!  Schedule construction
! ---------------------------------------------------------------------------------------- !

  subroutine Build_LTSSchedule2D(this,mesh,N,nVar)
    !! Derive the per-level element and mortar lists from mesh%elemLevel and
    !! mesh%mortarLevel, and size the registers and trace history.
    !!
    !! The element lists are rank-local; the mortar lists are global and identical on every
    !! rank, because mortarInfo and mortarLevel are replicated. That is what lets every rank
    !! walk the same recursion and post matching messages even for levels on which it owns
    !! no elements.
    implicit none
    class(LTSSchedule2D),intent(inout) :: this
    type(Mesh2D),intent(in) :: mesh
    integer,intent(in) :: N
    integer,intent(in) :: nVar
    ! Local
    integer :: L,e,m,n1,n2,n3,n4,lvl

    call this%Free()

    this%maxLevel = mesh%maxElemLevel
    this%nMortars = mesh%nMortars
    this%N = N
    this%nVar = nVar

    allocate(this%level(0:this%maxLevel))

    do L = 0,this%maxLevel

      ! Pass 1 : count. Pass 2 : fill. Element ids ascend, so each list is sorted, which
      ! keeps the restricted loops in the same relative order as the unrestricted ones.
      n1 = 0
      do e = 1,mesh%nElem
        if(mesh%elemLevel(e) == L) n1 = n1+1
      enddo
      allocate(this%level(L)%elems(1:n1))
      n1 = 0
      do e = 1,mesh%nElem
        if(mesh%elemLevel(e) == L) then
          n1 = n1+1
          this%level(L)%elems(n1) = e
        endif
      enddo

      n2 = 0
      n3 = 0
      n4 = 0
      do m = 1,mesh%nMortars
        lvl = mesh%mortarLevel(m)
        if(lvl == L) then
          n2 = n2+1
          n3 = n3+1
        elseif(lvl == L-1) then
          n2 = n2+1
          n4 = n4+1
        endif
      enddo
      allocate(this%level(L)%mortars(1:n2))
      allocate(this%level(L)%mortarsBig(1:n3))
      allocate(this%level(L)%mortarsSmall(1:n4))
      n2 = 0
      n3 = 0
      n4 = 0
      do m = 1,mesh%nMortars
        lvl = mesh%mortarLevel(m)
        if(lvl == L) then
          n2 = n2+1
          this%level(L)%mortars(n2) = m
          n3 = n3+1
          this%level(L)%mortarsBig(n3) = m
        elseif(lvl == L-1) then
          n2 = n2+1
          this%level(L)%mortars(n2) = m
          n4 = n4+1
          this%level(L)%mortarsSmall(n4) = m
        endif
      enddo

      this%level(L)%historyDepth = 0

    enddo

    if(this%nMortars > 0) then
      allocate(this%regCoarse(1:N+1,1:nVar,1:this%nMortars))
      allocate(this%regFine(1:N+1,1:nVar,1:this%nMortars))
      allocate(this%bigTrace(1:N+1,1:nVar,1:3,1:this%nMortars))
      this%regCoarse = 0.0_prec
      this%regFine = 0.0_prec
      this%bigTrace = 0.0_prec
    endif

    this%built = .true.

  endsubroutine Build_LTSSchedule2D

  logical function Matches_LTSSchedule2D(this,mesh,N,nVar) result(ok)
    !! Is this schedule still the one for the given mesh? Used to rebuild lazily after an
    !! adaptation epoch has replaced the mesh under the model.
    implicit none
    class(LTSSchedule2D),intent(in) :: this
    type(Mesh2D),intent(in) :: mesh
    integer,intent(in) :: N
    integer,intent(in) :: nVar

    ok = this%built .and. &
         this%maxLevel == mesh%maxElemLevel .and. &
         this%nMortars == mesh%nMortars .and. &
         this%N == N .and. &
         this%nVar == nVar

    if(ok) then
      ! Element counts are the cheap discriminator between two meshes that happen to share
      ! a level cap and mortar count.
      ok = (size(this%level(0)%elems)+CountAbove(this,0)) == mesh%nElem
    endif

  endfunction Matches_LTSSchedule2D

  integer function CountAbove(sched,L) result(n)
    !! Total number of rank-local elements at levels above L.
    implicit none
    class(LTSSchedule2D),intent(in) :: sched
    integer,intent(in) :: L
    ! Local
    integer :: k

    n = 0
    do k = L+1,sched%maxLevel
      n = n+size(sched%level(k)%elems)
    enddo

  endfunction CountAbove

  subroutine Free_LTSSchedule2D(this)
    implicit none
    class(LTSSchedule2D),intent(inout) :: this
    ! Local
    integer :: L

    if(allocated(this%level)) then
      do L = lbound(this%level,1),ubound(this%level,1)
        if(associated(this%level(L)%elems)) deallocate(this%level(L)%elems)
        if(associated(this%level(L)%mortars)) deallocate(this%level(L)%mortars)
        if(associated(this%level(L)%mortarsBig)) deallocate(this%level(L)%mortarsBig)
        if(associated(this%level(L)%mortarsSmall)) deallocate(this%level(L)%mortarsSmall)
        this%level(L)%elems => null()
        this%level(L)%mortars => null()
        this%level(L)%mortarsBig => null()
        this%level(L)%mortarsSmall => null()
      enddo
      deallocate(this%level)
    endif

    if(allocated(this%regCoarse)) deallocate(this%regCoarse)
    if(allocated(this%regFine)) deallocate(this%regFine)
    if(allocated(this%bigTrace)) deallocate(this%bigTrace)

    this%built = .false.
    this%maxLevel = 0
    this%nMortars = 0
    this%N = 0
    this%nVar = 0

  endsubroutine Free_LTSSchedule2D

! ---------------------------------------------------------------------------------------- !
!  The recursion
! ---------------------------------------------------------------------------------------- !

  subroutine ForwardStepLTS(model,sched,tn,dtBase,ioInterval)
    !! Advance the model to tn with local time stepping, writing output every ioInterval.
    !!
    !! dtBase is the time step for the BASE (level-0) mesh - not the level-derived global dt
    !! that AMRController2D%RecommendedTimeStep returns. Level ell is advanced at
    !! dtBase / 2**ell, so the stability constraint each level sees is its own.
    !!
    !! This deliberately mirrors ForwardStep_Model (entropy, metrics, model and tecplot
    !! output on the same cadence) rather than hooking into SetTimeIntegrator, so that local
    !! time stepping is opt-in and the single-rate path is untouched.
    implicit none
    class(DGModel2D_t),intent(inout) :: model
    type(LTSSchedule2D),intent(inout) :: sched
    real(prec),intent(in) :: tn
    real(prec),intent(in) :: dtBase
    real(prec),intent(in) :: ioInterval
    ! Local
    integer :: i,nIO
    real(prec) :: targetTime,tNext,tStep

    call GuardLTS(model)

    if(.not. sched%Matches(model%mesh,model%solution%interp%N,model%solution%nVar)) then
      call sched%Build(model%mesh,model%solution%interp%N,model%solution%nVar)
    endif

    model%dt = dtBase
    targetTime = tn
    nIO = int((targetTime-model%t)/ioInterval)

    do i = 1,nIO

      tNext = model%t+ioInterval

      do while(model%t < tNext)
        call model%PreStepHook()
        ! tStep is a COPY of model%t on purpose. StepLevel writes model%t as it walks the
        ! RK stages, so passing model%t directly would alias its own t0 argument and the
        ! stage times would shift under it.
        tStep = model%t
        call StepLevel(model,sched,0,min(dtBase,tNext-model%t),tStep)
        call model%PostStepHook()
      enddo

      model%t = tNext

      call model%CalculateEntropy()
      call model%ReportEntropy()
      call model%ReportMetrics()

      call model%WriteModel()
      if(model%tecplot_enabled) then
        call model%WriteTecplot()
      endif
      call model%IncrementIOCounter()

    enddo

  endsubroutine ForwardStepLTS

  subroutine GuardLTS(model)
    !! Reject the configurations this stage does not implement, loudly and up front rather
    !! than as a wrong answer. The CPU-only restriction is additionally enforced deeper down
    !! by RejectSubset, which fires as soon as a GPU operator is handed a subset.
    implicit none
    class(DGModel2D_t),intent(in) :: model
    ! Local
    integer :: m

    if(model%gradient_enabled) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : local time stepping does not support gradient-enabled (parabolic) '// &
        'models; the gradient traces across a level interface are not time-interpolated.'
      stop 1
    endif

    if(model%mesh%nMortars > 0) then
      if(.not. allocated(model%mesh%mortarLevel)) then
        print*,__FILE__,':',__LINE__, &
          ' : Error : local time stepping needs mesh%mortarLevel, which this mesh does '// &
          'not carry. Emit the mesh from a quadtree forest (EmitMesh) or use one of the '// &
          'built-in mortar meshes.'
        stop 1
      endif
      do m = 1,model%mesh%nMortars
        if(model%mesh%mortarLevel(m) < 0 .or. &
           model%mesh%mortarLevel(m) > model%mesh%maxElemLevel-1) then
          print*,__FILE__,':',__LINE__, &
            ' : Error : mortar ',m,' has big-side level ',model%mesh%mortarLevel(m), &
            ' which is not one level below the mesh level cap ',model%mesh%maxElemLevel, &
            '. The forest must be 2:1 balanced.'
          stop 1
        endif
      enddo
    endif

  endsubroutine GuardLTS

  recursive subroutine StepLevel(model,sched,L,dt,t0)
    !! Advance every element at refinement level L by dt, then recurse twice on level L+1 at
    !! dt/2 and reflux the level-L/L+1 interfaces.
    implicit none
    class(DGModel2D_t),intent(inout) :: model
    type(LTSSchedule2D),intent(inout) :: sched
    integer,intent(in) :: L
    real(prec),intent(in) :: dt
    real(prec),intent(in) :: t0
    ! Local
    integer :: m
    real(prec) :: dtSave,tBase

    ! Copy t0 before touching model%t. A caller may legitimately hand us its own t0 dummy
    ! (the recursive calls below do), which is ultimately anchored to model%t; reading t0
    ! after model%t has been advanced would then read the advanced value.
    tBase = t0
    dtSave = model%dt

    ! Registers for the interfaces this level owns are re-opened for this step: regCoarse
    ! fills during the stages below, regFine during the level-(L+1) substeps that follow.
    call ZeroRegisters(sched,L)

    ! Seed the big-side trace history the first time this level is stepped, so that the
    ! children have a t0 value to interpolate from.
    if(sched%level(L)%historyDepth == 0 .and. size(sched%level(L)%mortarsBig) > 0) then
      call model%solution%BoundaryInterp(sched%level(L)%elems)
      call CaptureBigTrace(model,sched,L)
      sched%level(L)%historyDepth = 1
    endif
    call ShiftBigTrace(sched,L)

    ! ---- Advance level L with the low-storage RK3, stage times exactly as in
    !      LowStorageRK3_timeIntegrator ----
    model%dt = dt
    model%t = tBase
    model%activeElem => sched%level(L)%elems
    model%activeMortar => sched%level(L)%mortars

    do m = 1,3
      ! The big sides of this level's mortarsSmall belong to level L-1, which has already
      ! finished its step; give them their trace at the current stage time.
      call WriteInterpolatedBigTrace(model,sched,L,model%t,dt)

      call model%CalculateTendency()

      ! flux%boundaryNormal now holds this stage's mortar surface integrand, already
      ! projected from the small sides by MortarFluxCollect. Bank it in both directions:
      ! as "what the big side used" for the interfaces this level owns, and as "what the
      ! fine level used" for the interfaces its parent owns.
      call AccumulateRegister(model,sched,L,LTS_REG_COARSE,rk3StageWeight(m)*dt)
      call AccumulateRegister(model,sched,L,LTS_REG_FINE,rk3StageWeight(m)*dt)

      call model%UpdateGRK3(m)
      model%t = tBase+rk3_b(m)*dt
    enddo

    model%t = tBase+dt
    model%activeElem => null()
    model%activeMortar => null()
    sched%level(L)%windowStart = tBase

    ! Record the level's new big-side trace at t0+dt, completing the interpolation window
    ! [t0, t0+dt] the children are about to step across.
    call model%solution%BoundaryInterp(sched%level(L)%elems)
    call CaptureBigTrace(model,sched,L)
    sched%level(L)%historyDepth = min(sched%level(L)%historyDepth+1,3)

    ! ---- Children, then reflux ----
    if(L < sched%maxLevel) then
      call StepLevel(model,sched,L+1,0.5_prec*dt,tBase)
      call StepLevel(model,sched,L+1,0.5_prec*dt,tBase+0.5_prec*dt)
      call Reflux(model,sched,L)
      ! The reflux changed the big elements' interiors, so the trace captured above is
      ! stale; refresh it before the children read it again.
      call model%solution%BoundaryInterp(sched%level(L)%elems)
      call CaptureBigTrace(model,sched,L)
      model%t = tBase+dt
    endif

    model%dt = dtSave

  endsubroutine StepLevel

! ---------------------------------------------------------------------------------------- !
!  Big-side trace history and time interpolation
! ---------------------------------------------------------------------------------------- !

  subroutine ShiftBigTrace(sched,L)
    !! Age the history of this level's big-side traces by one step: (t-dt, t, t+dt) becomes
    !! the window ending at the value about to be overwritten.
    implicit none
    type(LTSSchedule2D),intent(inout) :: sched
    integer,intent(in) :: L
    ! Local
    integer :: k,m

    do k = 1,size(sched%level(L)%mortarsBig)
      m = sched%level(L)%mortarsBig(k)
      sched%bigTrace(:,:,1,m) = sched%bigTrace(:,:,2,m)
      sched%bigTrace(:,:,2,m) = sched%bigTrace(:,:,3,m)
    enddo

  endsubroutine ShiftBigTrace

  subroutine CaptureBigTrace(model,sched,L)
    !! Copy the current big-side edge trace of every mortar this level owns into the newest
    !! history slot. Only the rank that owns the big element has the trace, and only that
    !! rank needs it: MortarExchange is what ships it to the small sides' ranks.
    implicit none
    class(DGModel2D_t),intent(in) :: model
    type(LTSSchedule2D),intent(inout) :: sched
    integer,intent(in) :: L
    ! Local
    integer :: k,m,eB,sB,offset,rankId

    rankId = model%mesh%decomp%rankId
    offset = model%mesh%decomp%offsetElem(rankId+1)

    do k = 1,size(sched%level(L)%mortarsBig)
      m = sched%level(L)%mortarsBig(k)
      eB = model%mesh%mortarInfo(1,m)
      if(model%mesh%decomp%elemToRank(eB) /= rankId) cycle
      sB = model%mesh%mortarInfo(2,m)
      sched%bigTrace(1:sched%N+1,1:sched%nVar,3,m) = &
        model%solution%boundary(1:sched%N+1,sB,eB-offset,1:sched%nVar)
    enddo

  endsubroutine CaptureBigTrace

  subroutine WriteInterpolatedBigTrace(model,sched,L,t,dtChild)
    !! For every mortar in which level L is a small side, overwrite the big element's stored
    !! edge trace with the parent's trace interpolated to time t, so that the MortarExchange
    !! inside CalculateTendency restricts a time-accurate value onto the small sides.
    !!
    !! This is safe to do to a live solution object: the big element belongs to level L-1,
    !! which is not being advanced now, and its trace is recomputed by BoundaryInterp the
    !! next time its own level steps.
    !!
    !! The interpolation is quadratic through the history slots (t0-dtP, t0, t0+dtP) once
    !! three are valid - third order, matching RK3 - and linear on (t0, t0+dtP) before then,
    !! which is the case on the first parent step of a run and just after an adaptation
    !! epoch has reset the history.
    implicit none
    class(DGModel2D_t),intent(inout) :: model
    type(LTSSchedule2D),intent(inout) :: sched
    integer,intent(in) :: L
    real(prec),intent(in) :: t
    real(prec),intent(in) :: dtChild
    ! Local
    integer :: k,m,eB,sB,offset,rankId,i,ivar,depth
    real(prec) :: dtParent,tParent0,s,w1,w2,w3

    if(L == 0) return
    if(size(sched%level(L)%mortarsSmall) == 0) return

    depth = sched%level(L-1)%historyDepth
    if(depth < 2) return ! nothing to interpolate between yet; leave the stored trace alone

    rankId = model%mesh%decomp%rankId
    offset = model%mesh%decomp%offsetElem(rankId+1)

    ! The parent's step is twice this level's, and its history window ends at the parent's
    ! current time; slot 2 sits one parent step back from slot 3.
    dtParent = 2.0_prec*dtChild
    tParent0 = sched%level(L-1)%windowStart
    s = (t-tParent0)/dtParent

    if(depth >= 3) then
      ! Lagrange through s = -1, 0, 1 (history slots 1, 2, 3)
      w1 = 0.5_prec*s*(s-1.0_prec)
      w2 = 1.0_prec-s*s
      w3 = 0.5_prec*s*(s+1.0_prec)
    else
      w1 = 0.0_prec
      w2 = 1.0_prec-s
      w3 = s
    endif

    do k = 1,size(sched%level(L)%mortarsSmall)
      m = sched%level(L)%mortarsSmall(k)
      eB = model%mesh%mortarInfo(1,m)
      if(model%mesh%decomp%elemToRank(eB) /= rankId) cycle
      sB = model%mesh%mortarInfo(2,m)
      do ivar = 1,sched%nVar
        do i = 1,sched%N+1
          model%solution%boundary(i,sB,eB-offset,ivar) = &
            w1*sched%bigTrace(i,ivar,1,m)+ &
            w2*sched%bigTrace(i,ivar,2,m)+ &
            w3*sched%bigTrace(i,ivar,3,m)
        enddo
      enddo
    enddo

  endsubroutine WriteInterpolatedBigTrace

! ---------------------------------------------------------------------------------------- !
!  Flux registers and refluxing
! ---------------------------------------------------------------------------------------- !

  subroutine ZeroRegisters(sched,L)
    !! Re-open both registers for the interfaces level L owns.
    implicit none
    type(LTSSchedule2D),intent(inout) :: sched
    integer,intent(in) :: L
    ! Local
    integer :: k,m

    if(sched%nMortars == 0) return

    do k = 1,size(sched%level(L)%mortarsBig)
      m = sched%level(L)%mortarsBig(k)
      sched%regCoarse(:,:,m) = 0.0_prec
      sched%regFine(:,:,m) = 0.0_prec
    enddo

  endsubroutine ZeroRegisters

  subroutine AccumulateRegister(model,sched,L,which,weight)
    !! Add weight * flux%boundaryNormal(big side) to one of the registers, over the mortars
    !! that level L owns (which = LTS_REG_COARSE) or in which it is a small side
    !! (which = LTS_REG_FINE).
    !!
    !! By the time this is called, MortarFluxCollect has already replaced the big side's
    !! integrand with the L2 projection of the two small sides' integrands, so a single read
    !! serves both roles: for the level that owns the interface it is the flux that level
    !! used; for the level below it, it is the exact projection of what the fine sides used.
    implicit none
    class(DGModel2D_t),intent(in) :: model
    type(LTSSchedule2D),intent(inout) :: sched
    integer,intent(in) :: L
    integer,intent(in) :: which
    real(prec),intent(in) :: weight
    ! Local
    integer :: k,m,eB,sB,offset,rankId,i,ivar,nm

    if(sched%nMortars == 0) return

    if(which == LTS_REG_COARSE) then
      nm = size(sched%level(L)%mortarsBig)
    else
      nm = size(sched%level(L)%mortarsSmall)
    endif
    if(nm == 0) return

    rankId = model%mesh%decomp%rankId
    offset = model%mesh%decomp%offsetElem(rankId+1)

    do k = 1,nm
      if(which == LTS_REG_COARSE) then
        m = sched%level(L)%mortarsBig(k)
      else
        m = sched%level(L)%mortarsSmall(k)
      endif
      eB = model%mesh%mortarInfo(1,m)
      if(model%mesh%decomp%elemToRank(eB) /= rankId) cycle
      sB = model%mesh%mortarInfo(2,m)
      do ivar = 1,sched%nVar
        do i = 1,sched%N+1
          if(which == LTS_REG_COARSE) then
            sched%regCoarse(i,ivar,m) = sched%regCoarse(i,ivar,m)+ &
                                        weight*model%flux%boundaryNormal(i,sB,eB-offset,ivar)
          else
            sched%regFine(i,ivar,m) = sched%regFine(i,ivar,m)+ &
                                      weight*model%flux%boundaryNormal(i,sB,eB-offset,ivar)
          endif
        enddo
      enddo
    enddo

  endsubroutine AccumulateRegister

  subroutine Reflux(model,sched,L)
    !! Replace, on every interface level L owns, the time-integrated flux the big side used
    !! with the one the fine level actually used.
    !!
    !! The DG surface term of MappedDGDivergence contributes, for a change Delta in the
    !! boundaryNormal of side sB,
    !!
    !!    dF(i,j) += B_sB(i,j;Delta) / J(i,j)
    !!
    !! with B given per side by the bMatrix/qWeights combinations below, and dSdt = source -
    !! dF. Over a full RK3 step the solution therefore picked up -B(regCoarse)/J from this
    !! interface where it should have picked up -B(regFine)/J, so the correction is
    !!
    !!    U(i,j) += B(regCoarse - regFine)(i,j) / J(i,j)
    !!
    !! which is exactly conservative because regFine is, by construction, minus the sum of
    !! the small sides' own time-integrated surface integrals (see MortarFluxCollect).
    implicit none
    class(DGModel2D_t),intent(inout) :: model
    type(LTSSchedule2D),intent(in) :: sched
    integer,intent(in) :: L
    ! Local
    integer :: k,m,eB,sB,offset,rankId,i,j,ivar,Np
    real(prec) :: d,surf

    rankId = model%mesh%decomp%rankId
    offset = model%mesh%decomp%offsetElem(rankId+1)
    Np = sched%N+1

    do k = 1,size(sched%level(L)%mortarsBig)
      m = sched%level(L)%mortarsBig(k)
      eB = model%mesh%mortarInfo(1,m)
      if(model%mesh%decomp%elemToRank(eB) /= rankId) cycle
      sB = model%mesh%mortarInfo(2,m)

      do ivar = 1,model%nstepped
        do j = 1,Np
          do i = 1,Np

            select case(sB)
            case(1) ! South : the trace runs along i
              d = sched%regCoarse(i,ivar,m)-sched%regFine(i,ivar,m)
              surf = model%solution%interp%bMatrix(j,1)*d/ &
                     model%solution%interp%qWeights(j)
            case(2) ! East : the trace runs along j
              d = sched%regCoarse(j,ivar,m)-sched%regFine(j,ivar,m)
              surf = model%solution%interp%bMatrix(i,2)*d/ &
                     model%solution%interp%qWeights(i)
            case(3) ! North : the trace runs along i
              d = sched%regCoarse(i,ivar,m)-sched%regFine(i,ivar,m)
              surf = model%solution%interp%bMatrix(j,2)*d/ &
                     model%solution%interp%qWeights(j)
            case default ! West : the trace runs along j
              d = sched%regCoarse(j,ivar,m)-sched%regFine(j,ivar,m)
              surf = model%solution%interp%bMatrix(i,1)*d/ &
                     model%solution%interp%qWeights(i)
            endselect

            model%solution%interior(i,j,eB-offset,ivar) = &
              model%solution%interior(i,j,eB-offset,ivar)+ &
              surf/model%geometry%J%interior(i,j,eB-offset,1)

          enddo
        enddo
      enddo

    enddo

  endsubroutine Reflux

endmodule SELF_LocalTimeStepping_2D
