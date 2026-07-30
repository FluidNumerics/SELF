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

program LinearEuler2D_AMR_Ultrasound_PointSource
!! An ultrasound-range point-source wavelet in water, propagating through a 1 m x 1 m domain
!! on a dynamically adapting (quadtree/mortar) mesh - the demonstration problem for 2-D AMR
!! with the linear Euler model.
!!
!! Physical configuration
!! ----------------------
!! Water at rest: sound speed c0 = 1500 m/s, background density rho0 = 1000 kg/m^3. A Gaussian
!! pressure pulse of half-width Lr = 3.75 mm (peak ~225 Pa) is released at the domain centre;
!! it radiates an annular wave packet whose width is ~4*Lr, i.e. a dominant wavelength of
!! ~15 mm and dominant frequency f0 ~ c0/(4*Lr) = 100 kHz - the low ultrasound band. All four
!! domain boundaries carry impedance-matched radiation (non-reflecting outflow) conditions.
!!
!! Discretization and adaptivity
!! -----------------------------
!! The base mesh is 16 x 16 elements (h0 = 62.5 mm) with degree N = 7 - deliberately far too
!! coarse for the wavelet. By default the AMR controller refines up to 3 levels (h = 7.8 mm at
!! level 3, ~15 points per wavelength), driven by the Legendre modal-decay indicator on the
!! pressure field, with a two-element refine-flag halo. The adaptation cadence is
!! stepsPerEpoch = 160 timesteps, which keeps the wavefront inside the refined band: the front
!! crosses one finest-level element in ~167 steps, and the halo buys another ~330. The time
!! step follows the level-based bound dtBase / 2**MaxLevel per epoch.
!!
!! Reaching real ultrasound frequencies is a matter of refinement depth rather than base-mesh
!! size (see SELF_AMR_ULTRASOUND_MAXLEVEL below); maxLevel 7 puts f0 at 1.6 MHz. Note the
!! cold-start behaviour this exposes: indicator-driven refinement can only chase structure the
!! *current* mesh resolves, and the base mesh's nodal spacing is 8.9 mm, so a target pulse
!! narrower than that would sample as zero on the level-0 mesh and refinement would never
!! start. For deeper-than-default targets the initial condition is therefore released as a
!! width ramp - wide enough for the base mesh to see, halved as each level appears - until the
!! target half-width is reached.
!!
!! Output and visualization
!! ------------------------
!! One HDF5 snapshot (solution + geometry, self-describing) is written per epoch, plus the
!! initial condition; because every snapshot carries its own geometry, the changing mesh needs
!! no special handling downstream. Render the pressure field with the mesh skeleton overlaid
!! using:
!!
!!     python examples/linear_euler2d_amr_plot.py <output directory>
!!
!! Environment overrides (all optional; the defaults are the CI-sized configuration):
!!
!!   SELF_AMR_ULTRASOUND_EPOCHS    number of adaptation epochs (default 2 = 10 us).
!!                                 Set 60 for the classic full-domain movie at maxLevel 3.
!!   SELF_AMR_ULTRASOUND_MAXLEVEL  refinement-level cap (default 3). Lr and epochLength track
!!                                 it automatically, so this alone selects the frequency band.
!!   SELF_AMR_ULTRASOUND_LR        source half-width in m, overriding the maxLevel-derived
!!                                 value; f0 ~ c0/(4*Lr).
!!   SELF_AMR_ULTRASOUND_EPOCHLEN  adaptation cadence in s, overriding stepsPerEpoch * dt.
!!
!! The run prints a CONFIG_* block (resolved resolution, f0, points per wavelength) and a
!! BENCH_* block (wall-clock split between time integration and adaptation) so a sweep over
!! maxLevel is interpretable and machine-parseable.
!!
!! The run doubles as an integration check: it verifies that refinement occurred, that the
!! acoustic energy stays finite and non-increasing (upwind flux + radiation boundaries are
!! dissipative; the solution transfer is conservative), and that the solution is NaN-free.

  use iso_fortran_env,only:int64,real64,output_unit
  use self_data
  use self_lineareuler2d
  use self_mesh_2d
  use SELF_Geometry_2D
  use SELF_AMRController_2D

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 7
  integer,parameter :: defaultMaxLevel = 3
  real(prec),parameter :: Lx = 1.0_prec ! domain extent (m)
  integer,parameter :: nElemX = 16 ! base-mesh elements per direction
  real(prec),parameter :: c0 = 1500.0_prec ! sound speed in water (m/s)
  real(prec),parameter :: rho0 = 1000.0_prec ! background density (kg/m^3)
  real(prec),parameter :: defaultLr = 3.75e-3_prec ! pulse half-width (m) at defaultMaxLevel
  real(prec),parameter :: rhoprime = 1.0e-4_prec ! density-anomaly amplitude (kg/m^3), ~225 Pa
  real(prec),parameter :: dtBase = 2.5e-7_prec ! stable on the level-0 mesh (s)
  ! Timesteps per adaptation epoch. The wave crosses one finest-level element in ~167 steps at
  ! this dt, and nHalo = 2 gives ~330 steps of margin, so the cadence must stay O(100) steps:
  ! epochLength is derived from this rather than fixed, so it tracks maxLevel automatically.
  integer,parameter :: stepsPerEpoch = 160
  ! CI-sized default (debug builds run every example ~5-10x slower than release; see issue
  ! #157). Two epochs still exercise the full loop: initial static refinement to the level
  ! cap, two dynamic adaptations, and the entropy/NaN checks.
  integer,parameter :: defaultEpochs = 2

  type(LinearEuler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  type(AMRController2D) :: controller
  integer :: bcids(1:4)
  integer :: nEpochs,epoch,i,envstat
  integer :: maxLevel
  logical :: adapted
  real(prec) :: e0,ef,dt
  real(prec) :: Lr,LrRamp,epochLength,hFinest,lambda,f0,ppw
  character(32) :: envstr
  ! Wall-clock instrumentation of the epoch loop. system_clock with an integer(int64) count
  ! is a monotonic wall clock; ForwardStep's own TIMER macro (src/SELF_Macros.h) expands to
  ! cpu_time in non-multithreaded builds and so cannot be used for GPU timing.
  integer(int64) :: c0clk,c1clk,c2clk,crate
  integer(int64) :: nStepsTotal
  real(real64) :: tFwd,tAdapt,tSim

  ! Number of adaptation epochs: CI-sized by default, override for movie production.
  nEpochs = defaultEpochs
  call get_environment_variable("SELF_AMR_ULTRASOUND_EPOCHS",envstr,status=envstat)
  if(envstat == 0) then
    read(envstr,*) nEpochs
  endif

  ! ---- Refinement depth and the quantities that must track it ----
  ! Reaching ultrasound length scales is a matter of refinement depth, not base-mesh size: the
  ! base mesh spacing is Lx/nElemX = 62.5 mm and each level halves it, so the finest spacing is
  ! 62.5 mm / 2**maxLevel.
  !
  ! Two other quantities must move with maxLevel, or extra depth buys nothing:
  !
  !   Lr  - the source half-width sets the dominant frequency, f0 ~ c0/(4 Lr). Holding Lr fixed
  !         while raising maxLevel just refines the same pulse until it is resolved and then
  !         stops, so Lr is halved per level. This holds points-per-wavelength constant at the
  !         default value (~15.4), and raises f0 geometrically with depth.
  !
  !   epochLength - dt = dtBase/2**maxLevel, so a fixed epochLength would put 2**k times more
  !         timesteps in an epoch and the wavefront would leave the refined patch mid-epoch.
  !         It is therefore derived as stepsPerEpoch * dt.
  !
  ! At maxLevel = defaultMaxLevel both reduce exactly to the original constants (Lr = 3.75 mm,
  ! epochLength = 5.0e-6 s), so default behaviour is unchanged.
  maxLevel = defaultMaxLevel
  call get_environment_variable("SELF_AMR_ULTRASOUND_MAXLEVEL",envstr,status=envstat)
  if(envstat == 0) then
    read(envstr,*) maxLevel
  endif

  Lr = defaultLr*2.0_prec**(defaultMaxLevel-maxLevel)
  call get_environment_variable("SELF_AMR_ULTRASOUND_LR",envstr,status=envstat)
  if(envstat == 0) then
    read(envstr,*) Lr
  endif

  epochLength = real(stepsPerEpoch,prec)*dtBase/2.0_prec**maxLevel
  call get_environment_variable("SELF_AMR_ULTRASOUND_EPOCHLEN",envstr,status=envstat)
  if(envstat == 0) then
    read(envstr,*) epochLength
  endif

  ! Resolution actually being requested, reported so a sweep over maxLevel is interpretable.
  hFinest = Lx/real(nElemX,prec)/2.0_prec**maxLevel
  lambda = 4.0_prec*Lr ! dominant wavelength of the wavelet
  f0 = c0/lambda
  ppw = real(controlDegree+1,prec)*lambda/hFinest
  write(output_unit,'(A,I0)') "CONFIG_maxLevel = ",maxLevel
  write(output_unit,'(A,ES16.7E3)') "CONFIG_hFinest_m = ",hFinest
  write(output_unit,'(A,ES16.7E3)') "CONFIG_Lr_m = ",Lr
  write(output_unit,'(A,ES16.7E3)') "CONFIG_lambda_m = ",lambda
  write(output_unit,'(A,ES16.7E3)') "CONFIG_f0_Hz = ",f0
  write(output_unit,'(A,ES16.7E3)') "CONFIG_pointsPerWavelength = ",ppw
  write(output_unit,'(A,ES16.7E3)') "CONFIG_epochLength_s = ",epochLength

  bcids(1:4) = [SELF_BC_RADIATION, & ! south
                SELF_BC_RADIATION, & ! east
                SELF_BC_RADIATION, & ! north
                SELF_BC_RADIATION] ! west

  call mesh%StructuredMesh(nElemX,nElemX,1,1,Lx/real(nElemX,prec),Lx/real(nElemX,prec),bcids)
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
  call modelobj%SetTimeIntegrator(integrator)

  call controller%Init(modelobj,refineThreshold=-3.0_prec,coarsenThreshold=-8.0_prec, &
                       ivar=3,maxLevel=maxLevel,nHalo=2)

  ! ---- Initial condition + static refinement to the pulse ----
  ! After each mesh change the pulse is re-evaluated analytically on the new mesh so the
  ! indicator sees the true field, not its coarse-mesh interpolant; the loop converges once
  ! the pulse is resolved (or the level cap is reached).
  !
  ! Cold start: indicator-driven refinement can only chase structure the current mesh can
  ! actually see. The base-mesh nodal spacing is Lx/(nElemX*controlDegree) = 8.9 mm, so a pulse
  ! narrower than that samples as ~0 on the level-0 mesh, the indicator finds nothing, and
  ! refinement never starts (observed: maxLevel >= 7 halts at "max level = 0"). The fix is a
  ! width ramp - begin with a pulse the base mesh resolves and halve its width as each level
  ! appears, so there is always resolvable structure to chase, until the target Lr is reached.
  ! Only engaged for deeper-than-default targets, so default behaviour is bit-for-bit unchanged.
  if(maxLevel > defaultMaxLevel) then
    LrRamp = defaultLr*2.0_prec**defaultMaxLevel ! widest: resolvable on the level-0 mesh
    call modelobj%SphericalSoundWave(rhoprime,LrRamp,0.5_prec*Lx,0.5_prec*Lx,c0)
    do i = 1,2*maxLevel+4
      call controller%Adapt(modelobj,adapted)
      LrRamp = max(Lr,defaultLr*2.0_prec**(defaultMaxLevel-controller%forest%MaxLevel()))
      call modelobj%SphericalSoundWave(rhoprime,LrRamp,0.5_prec*Lx,0.5_prec*Lx,c0)
      if(.not. adapted .and. LrRamp <= Lr) exit
    enddo
    print*,"cold-start ramp finished: Lr =",LrRamp,"(target",Lr,")"
  else
    call modelobj%SphericalSoundWave(rhoprime,Lr,0.5_prec*Lx,0.5_prec*Lx,c0)
    do i = 1,maxLevel+2
      call controller%Adapt(modelobj,adapted)
      if(.not. adapted) exit
      call modelobj%SphericalSoundWave(rhoprime,Lr,0.5_prec*Lx,0.5_prec*Lx,c0)
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
  ! The epoch loop is timed in two buckets: time integration (ForwardStep, which also writes
  ! one snapshot per epoch) and adaptation (Adapt). Adapt opens with solution%UpdateHost(), a
  ! blocking device-to-host copy, so the c1clk boundary lands after the GPU has drained.
  call system_clock(count_rate=crate)
  tFwd = 0.0_real64
  tAdapt = 0.0_real64
  nStepsTotal = 0
  do epoch = 1,nEpochs
    dt = controller%RecommendedTimeStep(dtBase)
    call system_clock(c0clk)
    ! ioInterval = (tn - t) exactly, so ForwardStep's io count is exactly 1 per epoch
    ! regardless of accumulated floating-point drift in t.
    call modelobj%ForwardStep(tn=modelobj%t+epochLength,dt=dt, &
                              ioInterval=(modelobj%t+epochLength)-modelobj%t)
    call system_clock(c1clk)
    call controller%Adapt(modelobj,adapted)
    call system_clock(c2clk)
    tFwd = tFwd+real(c1clk-c0clk,real64)/real(crate,real64)
    tAdapt = tAdapt+real(c2clk-c1clk,real64)/real(crate,real64)
    nStepsTotal = nStepsTotal+int(epochLength/dt,int64)
    print*,"epoch",epoch,": t =",modelobj%t,", dt =",dt, &
      ", nElem =",modelobj%mesh%nElem,", adapted =",adapted
  enddo

  ! ---- Wall-clock summary (fixed BENCH_ keys for machine parsing) ----
  tSim = real(nEpochs,real64)*real(epochLength,real64)
  write(output_unit,'(A)') "BENCH_case = ultrasound_pointsource"
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

  ! ---- Integration checks ----
  call modelobj%CalculateEntropy()
  ef = modelobj%entropy
  print*,"initial, final acoustic energy :",e0,ef
  if(ef /= ef) then
    print*,"ERROR: acoustic energy is NaN."
    stop 1
  endif
  if(ef > e0*(1.0_prec+10.0_prec*epsilon(1.0_prec))) then
    print*,"ERROR: acoustic energy grew over the adaptive run."
    stop 1
  endif

  call modelobj%Free()
  call controller%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram LinearEuler2D_AMR_Ultrasound_PointSource
