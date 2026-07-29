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
!! coarse for the wavelet. The AMR controller refines up to 3 levels (h = 7.8 mm at level 3,
!! ~15 points per wavelength), driven by the Legendre modal-decay indicator on the pressure
!! field, with a two-element refine-flag halo. The adaptation cadence (one epoch = 5 us) keeps
!! the wavefront inside the refined band: the front moves c0 * 5 us = 7.5 mm < 2 * h_fine per
!! epoch. The time step follows the level-based bound dtBase / 2**MaxLevel per epoch.
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
!! By default the run is CI-sized (2 epochs = 10 us). For a full-domain movie set
!! SELF_AMR_ULTRASOUND_EPOCHS=60 (300 us; the front approaches the radiation boundaries)
!! before running.
!!
!! The run doubles as an integration check: it verifies that refinement occurred, that the
!! acoustic energy stays finite and non-increasing (upwind flux + radiation boundaries are
!! dissipative; the solution transfer is conservative), and that the solution is NaN-free.

  use self_data
  use self_lineareuler2d
  use self_mesh_2d
  use SELF_Geometry_2D
  use SELF_AMRController_2D

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 7
  integer,parameter :: maxLevel = 3
  real(prec),parameter :: Lx = 1.0_prec ! domain extent (m)
  integer,parameter :: nElemX = 16 ! base-mesh elements per direction
  real(prec),parameter :: c0 = 1500.0_prec ! sound speed in water (m/s)
  real(prec),parameter :: rho0 = 1000.0_prec ! background density (kg/m^3)
  real(prec),parameter :: Lr = 3.75e-3_prec ! pulse half-width (m); f0 ~ c0/(4 Lr) = 100 kHz
  real(prec),parameter :: rhoprime = 1.0e-4_prec ! density-anomaly amplitude (kg/m^3), ~225 Pa
  real(prec),parameter :: dtBase = 2.5e-7_prec ! stable on the level-0 mesh (s)
  real(prec),parameter :: epochLength = 5.0e-6_prec ! adaptation cadence (s)
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
  logical :: adapted
  real(prec) :: e0,ef,dt
  character(32) :: envstr

  ! Number of adaptation epochs: CI-sized by default, override for movie production.
  nEpochs = defaultEpochs
  call get_environment_variable("SELF_AMR_ULTRASOUND_EPOCHS",envstr,status=envstat)
  if(envstat == 0) then
    read(envstr,*) nEpochs
  endif

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
  call modelobj%SphericalSoundWave(rhoprime,Lr,0.5_prec*Lx,0.5_prec*Lx,c0)
  do i = 1,maxLevel+2
    call controller%Adapt(modelobj,adapted)
    if(.not. adapted) exit
    call modelobj%SphericalSoundWave(rhoprime,Lr,0.5_prec*Lx,0.5_prec*Lx,c0)
  enddo
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
  do epoch = 1,nEpochs
    dt = controller%RecommendedTimeStep(dtBase)
    ! ioInterval = (tn - t) exactly, so ForwardStep's io count is exactly 1 per epoch
    ! regardless of accumulated floating-point drift in t.
    call modelobj%ForwardStep(tn=modelobj%t+epochLength,dt=dt, &
                              ioInterval=(modelobj%t+epochLength)-modelobj%t)
    call controller%Adapt(modelobj,adapted)
    print*,"epoch",epoch,": t =",modelobj%t,", dt =",dt, &
      ", nElem =",modelobj%mesh%nElem,", adapted =",adapted
  enddo

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
