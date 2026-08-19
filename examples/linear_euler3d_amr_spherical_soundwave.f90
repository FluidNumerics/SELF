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

program LinearEuler3D_AMR_SphericalSoundWave
!! A spherical acoustic pulse in water propagating through a 1 m^3 domain on a dynamically
!! adapting (octree / face-mortar) mesh - the demonstration problem for 3-D AMR with the linear
!! Euler model, and the benchmark vehicle for the 3-D adaptation cost.
!!
!! It is the 3-D counterpart of linear_euler2d_amr_ultrasound_pointsource, and deliberately
!! keeps that example's structure, environment overrides and BENCH_ reporting so the two are
!! directly comparable. The physics differ only in dimension.
!!
!! Physical configuration
!! ----------------------
!! Water at rest: sound speed c0 = 1500 m/s, background density rho0 = 1000 kg/m^3. A Gaussian
!! pressure pulse of half-width Lr is released at the domain centre and radiates as a spherical
!! shell whose width is ~4*Lr, i.e. a dominant wavelength lambda = 4*Lr and dominant frequency
!! f0 = c0/lambda. All six domain boundaries carry impedance-matched radiation (non-reflecting
!! outflow) conditions.
!!
!! Discretization and adaptivity
!! -----------------------------
!! The base mesh is 4 x 4 x 4 elements (h0 = 250 mm) with degree N = 4. At the default depth the
!! base mesh is marginal rather than hopeless for the pulse (10 points per wavelength once
!! refined, 5 at level 0); raising maxLevel halves Lr with it, so the case becomes progressively
!! more under-resolved at level 0 and the refined band correspondingly more localized. The AMR
!! controller is driven by the Legendre modal-decay indicator on the pressure field (variable 4
!! in 3-D), with a one-element refine-flag halo. The time step follows the level-based bound
!! dtBase / 2**MaxLevel per epoch.
!!
!! The base mesh is smaller and the degree lower than in the 2-D example on purpose: a level of
!! refinement multiplies the element count by 8 in 3-D rather than 4, and every element carries
!! (N+1)^3 nodes rather than (N+1)^2, so an equivalently-sized 3-D case is roughly two orders of
!! magnitude more work. The defaults are sized to run in a few seconds at -O3 so the example
!! stays inside the CI budget (examples run in the 90-minute coverage jobs); refinement depth is
!! the knob for benchmarking, via SELF_AMR_LE3D_MAXLEVEL.
!!
!! Environment overrides (all optional; the defaults are the CI-sized configuration):
!!
!!   SELF_AMR_LE3D_EPOCHS    number of adaptation epochs (default 2).
!!   SELF_AMR_LE3D_MAXLEVEL  refinement-level cap (default 1). Lr and epochLength track it
!!                           automatically, so this alone selects the frequency band, and it is
!!                           the knob for a depth sweep.
!!   SELF_AMR_LE3D_LR        source half-width in m, overriding the maxLevel-derived value;
!!                           f0 = c0/(4*Lr).
!!   SELF_AMR_LE3D_EPOCHLEN  adaptation cadence in s, overriding stepsPerEpoch * dt.
!!   SELF_AMR_LE3D_RELFLOOR  relative energy floor of the indicator's amplitude gate
!!                           (default SELF_AMR_DEFAULT_RELFLOOR); 0 disables the gate.
!!   SELF_AMR_LE3D_SIGNIFFLOOR  upper edge of the gate's hysteresis band (default: the floor
!!                           itself, i.e. a single hard cut).
!!
!! The run prints a CONFIG_* block (resolved resolution, f0, points per wavelength) and a BENCH_*
!! block splitting epoch wall time between time integration and adaptation, so a sweep over
!! maxLevel is interpretable and machine-parseable. Epochs in which the leaf set does not change
!! cost only an indicator pass, so the adaptation bucket is reported both in total and per
!! ADAPTING epoch (BENCH_nAdaptEpochs, BENCH_tAdaptPerAdaptation_s): the per-adaptation figure is
!! the one to compare between builds, because the number of adapting epochs is not guaranteed
!! identical across backends. BENCH_elemTransferred is its denominator - the elements the
!! transfer actually produced, summed over adapting epochs.
!!
!! Snapshots are written to the working directory as self-describing HDF5 (solution + geometry
!! per file), one per epoch plus the initial condition. Because every snapshot carries its own
!! geometry, the changing mesh needs no special handling downstream.
!!
!! The run doubles as an integration check: it verifies that refinement occurred, that the
!! acoustic energy stays finite and non-increasing (upwind flux + radiation boundaries are
!! dissipative; the solution transfer is conservative), and that the solution is NaN-free.

  use iso_fortran_env,only:int64,real64,output_unit
  use self_data
  use self_lineareuler3d
  use self_mesh_3d
  use SELF_Geometry_3D
  use SELF_AMRController_3D

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 4
  integer,parameter :: targetDegree = 4
  integer,parameter :: defaultMaxLevel = 1
  real(prec),parameter :: Lx = 1.0_prec ! domain extent (m)
  integer,parameter :: nElemX = 4 ! base-mesh elements per direction
  real(prec),parameter :: c0 = 1500.0_prec ! sound speed in water (m/s)
  real(prec),parameter :: rho0 = 1000.0_prec ! background density (kg/m^3)
  real(prec),parameter :: defaultLr = 6.25e-2_prec ! pulse half-width (m) at defaultMaxLevel
  real(prec),parameter :: rhoprime = 1.0e-4_prec ! density-anomaly amplitude (kg/m^3)
  real(prec),parameter :: dtBase = 4.0e-6_prec ! stable on the level-0 mesh (s)
  ! Timesteps per adaptation epoch. Because dt = dtBase/2**maxLevel and hFinest = h0/2**maxLevel
  ! scale together, the wavefront crosses one finest-level element in hFinest/(c0*dt) = 42 steps
  ! at any depth. The cadence is set to 0.6 of that, so nHalo = 1 still holds the front inside
  ! the refined band across an epoch that happens not to adapt - the 2-D example buys the same
  ! margin with nHalo = 2 instead, which is far more expensive in 3-D (26 face/edge/corner
  ! neighbours per element rather than 8). epochLength is derived from this rather than fixed, so
  ! it tracks maxLevel automatically.
  integer,parameter :: stepsPerEpoch = 25
  ! CI-sized default. Two epochs still exercise the full loop: initial static refinement to the
  ! level cap, two dynamic adaptations, and the entropy/NaN checks.
  integer,parameter :: defaultEpochs = 2
  ! Slack on the entropy non-growth check. The dissipation argument is exact in exact arithmetic,
  ! so the tolerance is pure round-off accumulation over the run; single precision needs far more
  ! of it. Same idiom (and same values) as test/lineareuler3d_amr_soundwave.f90.
#ifdef DOUBLE_PRECISION
  real(prec),parameter :: entropyTol = 1.0e-10_prec
#else
  real(prec),parameter :: entropyTol = 1.0e-3_prec
#endif

  type(LinearEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  type(AMRController3D) :: controller
  integer :: bcids(1:6)
  integer :: nEpochs,epoch,i,envstat
  integer :: nAdaptEpochs
  integer(int64) :: nElemTransferred
  integer :: maxLevel
  logical :: adapted
  real(prec) :: e0,ef,dt
  real(prec) :: Lr,LrRamp,epochLength,hFinest,lambda,f0,ppw
  real(prec) :: relativeEnergyFloor,significantEnergyFloor
  character(32) :: envstr
  ! Wall-clock instrumentation of the epoch loop. system_clock with an integer(int64) count is a
  ! monotonic wall clock; ForwardStep's own TIMER macro (src/SELF_Macros.h) expands to cpu_time
  ! in non-multithreaded builds and so cannot be used for GPU timing.
  integer(int64) :: c0clk,c1clk,c2clk,crate
  integer(int64) :: nStepsTotal
  real(real64) :: tFwd,tAdapt,tSim

  ! Number of adaptation epochs: CI-sized by default, override for a benchmark sweep.
  nEpochs = defaultEpochs
  call get_environment_variable("SELF_AMR_LE3D_EPOCHS",envstr,status=envstat)
  if(envstat == 0) then
    read(envstr,*) nEpochs
  endif

  ! ---- Refinement depth and the quantities that must track it ----
  ! The base mesh spacing is Lx/nElemX = 250 mm and each level halves it, so the finest spacing
  ! is 250 mm / 2**maxLevel.
  !
  ! Two other quantities must move with maxLevel, or extra depth buys nothing:
  !
  !   Lr  - the source half-width sets the dominant frequency, f0 = c0/(4 Lr). Holding Lr fixed
  !         while raising maxLevel just refines the same pulse until it is resolved and then
  !         stops, so Lr is halved per level. This holds points-per-wavelength constant and
  !         raises f0 geometrically with depth.
  !
  !   epochLength - dt = dtBase/2**maxLevel, so a fixed epochLength would put 2**k times more
  !         timesteps in an epoch and the wavefront would leave the refined patch mid-epoch. It
  !         is therefore derived as stepsPerEpoch * dt.
  !
  ! At maxLevel = defaultMaxLevel both reduce to the documented default constants.
  maxLevel = defaultMaxLevel
  call get_environment_variable("SELF_AMR_LE3D_MAXLEVEL",envstr,status=envstat)
  if(envstat == 0) then
    read(envstr,*) maxLevel
  endif

  Lr = defaultLr*2.0_prec**(defaultMaxLevel-maxLevel)
  call get_environment_variable("SELF_AMR_LE3D_LR",envstr,status=envstat)
  if(envstat == 0) then
    read(envstr,*) Lr
  endif

  epochLength = real(stepsPerEpoch,prec)*dtBase/2.0_prec**maxLevel
  call get_environment_variable("SELF_AMR_LE3D_EPOCHLEN",envstr,status=envstat)
  if(envstat == 0) then
    read(envstr,*) epochLength
  endif

  ! Amplitude gate of the refinement indicator (see SELF_RefinementIndicator_3D): elements below
  ! this fraction of the field's peak element ENERGY are treated as quiescent and released, which
  ! is what lets the mesh coarsen behind the passing shell. 0 restores the pre-gate behaviour, in
  ! which the refined region grows to the whole swept volume.
  relativeEnergyFloor = SELF_AMR_DEFAULT_RELFLOOR
  call get_environment_variable("SELF_AMR_LE3D_RELFLOOR",envstr,status=envstat)
  if(envstat == 0) then
    read(envstr,*) relativeEnergyFloor
  endif

  ! Upper edge of the amplitude gate's hysteresis band; defaults to the floor itself, i.e. a
  ! single hard cut.
  significantEnergyFloor = relativeEnergyFloor
  call get_environment_variable("SELF_AMR_LE3D_SIGNIFFLOOR",envstr,status=envstat)
  if(envstat == 0) then
    read(envstr,*) significantEnergyFloor
  endif

  ! Resolution actually being requested, reported so a sweep over maxLevel is interpretable.
  hFinest = Lx/real(nElemX,prec)/2.0_prec**maxLevel
  lambda = 4.0_prec*Lr ! dominant wavelength of the pulse
  f0 = c0/lambda
  ppw = real(controlDegree+1,prec)*lambda/hFinest
  write(output_unit,'(A,I0)') "CONFIG_maxLevel = ",maxLevel
  write(output_unit,'(A,ES16.7E3)') "CONFIG_hFinest_m = ",hFinest
  write(output_unit,'(A,ES16.7E3)') "CONFIG_Lr_m = ",Lr
  write(output_unit,'(A,ES16.7E3)') "CONFIG_lambda_m = ",lambda
  write(output_unit,'(A,ES16.7E3)') "CONFIG_f0_Hz = ",f0
  write(output_unit,'(A,ES16.7E3)') "CONFIG_pointsPerWavelength = ",ppw
  write(output_unit,'(A,ES16.7E3)') "CONFIG_epochLength_s = ",epochLength
  write(output_unit,'(A,ES16.7E3)') "CONFIG_relativeEnergyFloor = ",relativeEnergyFloor
  write(output_unit,'(A,ES16.7E3)') "CONFIG_significantEnergyFloor = ",significantEnergyFloor

  bcids(1:6) = [SELF_BC_RADIATION, & ! bottom
                SELF_BC_RADIATION, & ! south
                SELF_BC_RADIATION, & ! east
                SELF_BC_RADIATION, & ! north
                SELF_BC_RADIATION, & ! west
                SELF_BC_RADIATION] ! top

  call mesh%StructuredMesh(nElemX,nElemX,nElemX,1,1,1, &
                           Lx/real(nElemX,prec),Lx/real(nElemX,prec),Lx/real(nElemX,prec),bcids)
  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%prescribed_bcs_enabled = .false.
  modelobj%tecplot_enabled = .false.
  modelobj%rho0 = rho0
  modelobj%c = c0
  call modelobj%SetTimeIntegrator(integrator)

  ! Pressure (variable 4 in 3-D) drives the indicator; thresholds are the recommended starting
  ! values from the AMR documentation.
  call controller%Init(modelobj,refineThreshold=-3.0_prec,coarsenThreshold=-8.0_prec, &
                       ivar=4,maxLevel=maxLevel,nHalo=1, &
                       relativeEnergyFloor=relativeEnergyFloor, &
                       significantEnergyFloor=significantEnergyFloor)

  ! ---- Initial condition + static refinement to the pulse ----
  ! After each mesh change the pulse is re-evaluated analytically on the new mesh so the
  ! indicator sees the true field, not its coarse-mesh interpolant; the loop converges once the
  ! pulse is resolved (or the level cap is reached).
  !
  ! Cold start: indicator-driven refinement can only chase structure the current mesh can
  ! actually see. The base-mesh nodal spacing is Lx/(nElemX*controlDegree) = 62.5 mm - the average
  ! gap between the N+1 Gauss nodes, which span N intervals, so it is the N convention rather than
  ! the N+1 one CONFIG_pointsPerWavelength counts with - so a pulse
  ! narrower than that samples as ~0 on the level-0 mesh, the indicator finds nothing, and
  ! refinement never starts. The fix is a width ramp - begin with a pulse the base mesh resolves
  ! and halve its width as each level appears, so there is always resolvable structure to chase,
  ! until the target Lr is reached. Only engaged for deeper-than-default targets, so default
  ! behaviour is unchanged.
  if(maxLevel > defaultMaxLevel) then
    LrRamp = defaultLr*2.0_prec**defaultMaxLevel ! widest: resolvable on the level-0 mesh
    call modelobj%SphericalSoundWave(rhoprime,LrRamp,0.5_prec*Lx,0.5_prec*Lx,0.5_prec*Lx)
    do i = 1,2*maxLevel+4
      call controller%Adapt(modelobj,adapted)
      LrRamp = max(Lr,defaultLr*2.0_prec**(defaultMaxLevel-controller%forest%MaxLevel()))
      call modelobj%SphericalSoundWave(rhoprime,LrRamp,0.5_prec*Lx,0.5_prec*Lx,0.5_prec*Lx)
      if(.not. adapted .and. LrRamp <= Lr) exit
    enddo
    print*,"cold-start ramp finished: Lr =",LrRamp,"(target",Lr,")"
  else
    call modelobj%SphericalSoundWave(rhoprime,Lr,0.5_prec*Lx,0.5_prec*Lx,0.5_prec*Lx)
    do i = 1,maxLevel+2
      call controller%Adapt(modelobj,adapted)
      if(.not. adapted) exit
      call modelobj%SphericalSoundWave(rhoprime,Lr,0.5_prec*Lx,0.5_prec*Lx,0.5_prec*Lx)
    enddo
  endif
  print*,"initial adaptation: nElem =",modelobj%mesh%nElem, &
    ", max level =",controller%forest%MaxLevel()
  if(controller%forest%MaxLevel() < 1) then
    print*,"ERROR: the indicator did not refine around the pulse."
    stop 1
  endif

  call modelobj%CalculateEntropy()
  e0 = modelobj%entropy

  ! Snapshot the initial condition (epoch snapshots follow from ForwardStep's file IO).
  call modelobj%WriteModel()
  call modelobj%IncrementIOCounter()

  ! ---- Time-step through adaptation epochs ----
  ! The epoch loop is timed in two buckets: time integration (ForwardStep) and adaptation (Adapt).
  !
  ! What the boundaries actually guarantee, on a GPU build. c1clk is clean: ForwardStep ends with
  ! an entropy reduction that reads back to the host, so the queue is drained there. c2clk is NOT
  ! guaranteed clean: on a GPU build Adapt ends by launching the device solution transfer and does
  ! not synchronize afterwards, so some of that kernel's time can land in the NEXT epoch's
  ! ForwardStep bucket instead of this epoch's Adapt bucket. No synchronize is added here rather
  ! than perturbing what is being measured; the effect is bounded by one transfer kernel per epoch,
  ! which is orders of magnitude below the adaptation cost it is being compared against, and
  ! tForwardStep is reported so any leakage into it is visible.
  call system_clock(count_rate=crate)
  tFwd = 0.0_real64
  tAdapt = 0.0_real64
  nStepsTotal = 0
  nAdaptEpochs = 0
  nElemTransferred = 0
  do epoch = 1,nEpochs
    dt = controller%RecommendedTimeStep(dtBase)
    call system_clock(c0clk)
    ! ioInterval = (tn - t) exactly, so ForwardStep's io count is exactly 1 per epoch regardless
    ! of accumulated floating-point drift in t.
    call modelobj%ForwardStep(tn=modelobj%t+epochLength,dt=dt, &
                              ioInterval=(modelobj%t+epochLength)-modelobj%t)
    call system_clock(c1clk)
    call controller%Adapt(modelobj,adapted)
    call system_clock(c2clk)
    tFwd = tFwd+real(c1clk-c0clk,real64)/real(crate,real64)
    tAdapt = tAdapt+real(c2clk-c1clk,real64)/real(crate,real64)
    nStepsTotal = nStepsTotal+int(epochLength/dt,int64)
    if(adapted) then
      nAdaptEpochs = nAdaptEpochs+1
      nElemTransferred = nElemTransferred+int(modelobj%mesh%nElem,int64)
    endif
    print*,"epoch",epoch,": t =",modelobj%t,", dt =",dt, &
      ", nElem =",modelobj%mesh%nElem,", adapted =",adapted
  enddo

  ! ---- Wall-clock summary (fixed BENCH_ keys for machine parsing) ----
  tSim = real(nEpochs,real64)*real(epochLength,real64)
  write(output_unit,'(A)') "BENCH_case = lineareuler3d_amr_sphericalsoundwave"
  write(output_unit,'(A,I0)') "BENCH_nEpochs = ",nEpochs
  write(output_unit,'(A,I0)') "BENCH_nSteps = ",nStepsTotal
  write(output_unit,'(A,I0)') "BENCH_nElemFinal = ",modelobj%mesh%nElem
  write(output_unit,'(A,ES16.7E3)') "BENCH_simTime_s = ",tSim
  write(output_unit,'(A,ES16.7E3)') "BENCH_tForwardStep_s = ",tFwd
  write(output_unit,'(A,ES16.7E3)') "BENCH_tAdapt_s = ",tAdapt
  write(output_unit,'(A,ES16.7E3)') "BENCH_tEpochLoop_s = ",tFwd+tAdapt
  write(output_unit,'(A,ES16.7E3)') "BENCH_wallPerStep_s = ",tFwd/real(nStepsTotal,real64)
  write(output_unit,'(A,ES16.7E3)') "BENCH_wallPerSimTime = ",tFwd/tSim
  write(output_unit,'(A,ES16.7E3)') "BENCH_wallPerSimTimeIncAMR = ",(tFwd+tAdapt)/tSim
  write(output_unit,'(A,ES16.7E3)') "BENCH_amrFraction = ",tAdapt/(tFwd+tAdapt)
  ! An epoch whose leaf set is unchanged costs only an indicator pass, so the total above mixes
  ! two very different quantities. The per-adaptation figure below is the one to compare between
  ! builds, with the transferred element count as its denominator.
  write(output_unit,'(A,I0)') "BENCH_nAdaptEpochs = ",nAdaptEpochs
  write(output_unit,'(A,I0)') "BENCH_elemTransferred = ",nElemTransferred
  if(nAdaptEpochs > 0) then
    write(output_unit,'(A,ES16.7E3)') "BENCH_tAdaptPerAdaptation_s = ", &
      tAdapt/real(nAdaptEpochs,real64)
  endif
  ! Migration volume (issue #167). These are per-RANK totals over the whole run, and they are the
  ! primary comparison between the point-to-point migration and SELF_AMR_MIGRATE_GATHER=1: the
  ! gather path receives every element the rank does not own on every adapting epoch, the
  ! point-to-point path only the old elements its new range actually reads from a peer. Both
  ! paths count, so the two runs are directly comparable. Zero on a single rank, where nothing
  ! migrates at all.
  write(output_unit,'(A,I0)') "BENCH_migrateBytesRecv = ",controller%nMigrateBytesRecv
  write(output_unit,'(A,I0)') "BENCH_migrateBytesSent = ",controller%nMigrateBytesSent
  write(output_unit,'(A,I0)') "BENCH_migrateElemRemote = ",controller%nMigrateElemRemote
  ! Geometry reuse (AMR Stage 6c): the share of elements whose geometry was carried forward from
  ! the previous epoch instead of being regenerated.
  write(output_unit,'(A,I0)') "BENCH_geomReused = ",controller%nGeomReused
  write(output_unit,'(A,I0)') "BENCH_geomGenerated = ",controller%nGeomGenerated
  if(controller%nGeomReused+controller%nGeomGenerated > 0) then
    write(output_unit,'(A,ES16.7E3)') "BENCH_geomReuseFraction = ", &
      real(controller%nGeomReused,real64)/ &
      real(controller%nGeomReused+controller%nGeomGenerated,real64)
  endif

  ! ---- Integration checks ----
  call modelobj%CalculateEntropy()
  ef = modelobj%entropy
  print*,"initial, final acoustic energy :",e0,ef
  if(ef /= ef) then
    print*,"ERROR: acoustic energy is NaN."
    stop 1
  endif
  if(ef-e0 > entropyTol*max(abs(e0),1.0_prec)) then
    print*,"ERROR: acoustic energy grew over the adaptive run."
    stop 1
  endif

  call modelobj%Free()
  call controller%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram LinearEuler3D_AMR_SphericalSoundWave
