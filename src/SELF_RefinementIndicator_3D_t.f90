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

module SELF_RefinementIndicator_3D_t
!! Modal (Legendre spectral-decay) refinement indicator for 3-D spectral element solutions.
!!
!! ## Purpose
!!
!! Provides the *trigger* mechanism for adaptive mesh refinement (AMR): a per-element scalar
!! that measures how much of a solution field's energy sits in its highest polynomial modes.
!! A well-resolved (smooth) field has energy that decays rapidly with mode number, so almost
!! all of its energy lives in the low modes; an under-resolved field (a steep front, a
!! discontinuity, or simply too-coarse an element) leaves significant energy in the top modes.
!! Comparing the top-mode energy fraction against thresholds flags each element for refinement,
!! coarsening, or no change.
!!
!! ## Method (Legendre modal-energy / spectral-decay indicator)
!!
!! Following Persson & Peraire (2006), "Sub-cell shock capturing for discontinuous Galerkin
!! methods", AIAA 2006-112, and the robustified variant of Hennemann, Rueda-Ramirez, Hindenlang
!! & Gassner (2021), J. Comput. Phys. 426, 109935 (the shock indicator used in Trixi.jl), the
!! nodal solution on an element is expanded in a tensor-product Legendre modal basis
!!
!!   u(xi,eta,zeta) = sum_{p=0}^{N} sum_{q=0}^{N} sum_{r=0}^{N}
!!                    uhat(p,q,r) Ltilde_p(xi) Ltilde_q(eta) Ltilde_r(zeta),
!!
!! where Ltilde_p are the L2-normalized Legendre polynomials on [-1,1]
!! (int_{-1}^{1} Ltilde_p Ltilde_q dxi = delta_pq). With that normalization the total modal
!! energy equals the exact L2 energy of the field on the reference element,
!!
!!   E_tot = sum_{p,q,r} uhat(p,q,r)^2 = ||u||_{L2([-1,1]^3)}^2 .
!!
!! Two "clipped" energies are formed by dropping the highest one and two modes in each direction,
!!
!!   E_clip1 = sum_{p,q,r <= N-1} uhat(p,q,r)^2 ,   E_clip2 = sum_{p,q,r <= N-2} uhat(p,q,r)^2 ,
!!
!! and the smoothness ratio is the larger of the top-shell and next-shell energy fractions,
!!
!!   S_e = max( (E_tot - E_clip1)/E_tot , (E_clip1 - E_clip2)/E_clip1 ) .
!!
!! The second term (present only for N >= 2) guards against the odd/even parity dropouts that a
!! single-mode measure can suffer for symmetric data - this is the Hennemann-Gassner refinement
!! of the original single-mode Persson-Peraire estimate. The indicator returned per element is
!! the base-10 logarithm sigma_e = log10(S_e); a near-machine-zero floor keeps it finite for a
!! perfectly resolved (or identically zero) field.
!!
!! ## Amplitude gate (relative energy floor)
!!
!! S_e is a ratio and therefore scale-free: it says nothing about how much energy an element
!! carries, only how that energy is distributed across modes. An element holding a negligible
!! share of the field's energy - low-amplitude grid-scale residue in the wake of a passing wave,
!! say - still shows a flat modal spectrum and so reads as under-resolved forever. Guarding the
!! ratio with an absolute floor at machine epsilon does not help: E_tot carries the SQUARE of the
!! field amplitude, so a wake at 1e-5 of the peak amplitude sits ~1e10 above epsilon and still
!! takes the ratio branch. Such elements are never flagged for coarsening and the mesh behind a
!! front is never released.
!!
!! The fix is a second, RELATIVE floor. A per-element gate energy
!!
!!   g_e = sum_v w_v E_tot,e,v
!!
!! is formed from the modal energies of the driving variables with non-negative weights w_v. Since
!! E_tot,e,v is the exact L2 energy of variable v on the reference element, g_e with w taken from
!! a quadratic entropy function is the discrete entropy integral over the element: a convex
!! function of the state, and therefore an amplitude measure that is meaningful across variables of
!! very different magnitude. Comparing it against a field-wide energy scale,
!!
!!   effective floor = max( epsilon, relativeEnergyFloor * energyScale )
!!   g_e <= effective floor   ->   S_e := 0   (element is quiescent, hence resolved)
!!
!! gives an amplitude gate that is independent of the modal shape. Because energy goes as
!! amplitude squared, relativeEnergyFloor is 10**(dB/10) in amplitude terms: the default 1e-12
!! treats amplitudes below 1e-6 (-120 dB) of the field scale as quiescent.
!!
!! Raising the floor is NOT free, and the useful range is narrower than it looks. The gate cannot
!! distinguish residue from the low-amplitude flank of a feature that is genuinely under-resolved,
!! so an aggressive floor also suppresses refinement DEPTH: on the ultrasound benchmark a floor of
!! 1e-8 stopped the source pulse from reaching the level cap, which doubled the level-based time
!! step and produced a mesh roughly twice as large by the end of the run. Measure before raising
!! it; see the note on SELF_AMR_DEFAULT_RELFLOOR.
!!
!! The gate is a single hard cut by default, which reintroduces on the energy axis exactly the
!! thrashing the two sigma thresholds exist to prevent: an element whose energy drifts across that
!! one value flips COARSEN <-> REFINE on successive epochs. Supplying significantEnergyFloor opens
!! a hysteresis band instead,
!!
!!   g_e <= quiescent floor            ->  COARSEN, whatever the spectrum says
!!   quiescent < g_e <= significant    ->  SELF_AMR_KEEP  (hold the mesh)
!!   g_e > significant floor           ->  the spectrum decides, as before
!!
!! the middle zone reading "too weak to be worth spending levels on, too strong to declare
!! resolved". Inside the band sigma_e still reports the true spectrum; only the flag is held.
!!
!! energyScale defaults to the largest gate energy over the elements (globally reduced when
!! Estimate is given an MPI communicator), which makes the gate track a decaying front. It can be
!! pinned with SetEnergyScale for callers that need the flag decision to be independent of a
!! floating-point reduction, and the gate itself can be supplied outright through the optional
!! gate argument to Estimate (for example an exactly integrated nodal entropy).
!!
!! Note that raising coarsenThreshold is NOT a substitute: in a low-amplitude wake sigma_e is near
!! 0, so any threshold high enough to coarsen the wake also stops a genuine front from refining.
!! The two decisions cannot be separated by thresholds alone.
!!
!! ## Trigger semantics
!!
!! With user-supplied thresholds sigma_refine > sigma_coarsen the per-element flag is
!!
!!   sigma_e > sigma_refine   -> SELF_AMR_REFINE  (+1)   top modes carry too much energy
!!   sigma_e < sigma_coarsen  -> SELF_AMR_COARSEN (-1)   field is over-resolved on this element
!!   otherwise                -> SELF_AMR_KEEP     (0)
!!
!! Larger (closer to zero) sigma_e means a less smooth / less resolved field. Recommended
!! starting values for double precision are sigma_refine ~ -3.0 and sigma_coarsen ~ -8.0; they
!! are problem dependent and are deliberately left as required arguments rather than hidden
!! defaults.
!!
!! ## Units and input ranges
!!
!!   - The indicator is dimensionless (an energy fraction); the driving field may carry any units.
!!   - Requires N >= 1 (a degree-0 element has no modal spectrum to measure).
!!   - The nodal->modal transform is the exact inverse of the Legendre Vandermonde built from the
!!     interpolant control points, so the indicator is independent of the control-node type
!!     (Legendre-Gauss or Legendre-Gauss-Lobatto) and is exact for polynomial data.
!!
!! This module contains the portable (CPU) implementation using do concurrent over elements.
!! The backend extension modules (SELF_RefinementIndicator_3D) add device storage and a GPU
!! kernel while preserving the mathematics bit-for-bit in structure.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Scalar_3D
  use iso_c_binding
  use mpi

  implicit none

  ! Per-element refinement flags returned in this%flag(:)
  integer,parameter :: SELF_AMR_COARSEN = -1
  integer,parameter :: SELF_AMR_KEEP = 0
  integer,parameter :: SELF_AMR_REFINE = 1

  ! Reduce the indicator over all variables when this sentinel is passed as the
  ! driving-variable index to Estimate.
  integer,parameter :: SELF_AMR_ALLVARS = 0

  ! Default relative energy floor of the amplitude gate (see the module header). 1e-12 in energy
  ! is 1e-6 (-120 dB) in amplitude relative to the field scale.
  !
  ! This default is deliberately conservative, and measured rather than argued. On the ultrasound
  ! point-source benchmark (examples/linear_euler2d_amr_ultrasound_pointsource, 30 epochs,
  ! maxLevel 2) it gives 1.40x fewer elements on average than no gate at all, with an initial
  ! adaptation - and hence a level-based time step - identical to the ungated run. Raising it to
  ! 1e-8 was measured to be actively HARMFUL there: the gate then also suppresses the source
  ! pulse's skirt, the forest never reaches the level cap, RecommendedTimeStep doubles dt, and the
  ! resulting time-integration error drives so much later refinement that the final mesh is
  ! roughly twice the ungated one (3304 elements against 1732). An aggressive amplitude gate
  ! trades against refinement DEPTH; that trade is what this value is set to avoid.
  !
  ! Note this sits below epsilon in real32 (1.2e-7), so in a single-precision build the absolute
  ! term of the effective floor dominates and the relative gate is inactive unless a caller raises
  ! it. That is the safe direction: real32 fields cannot represent -120 dB structure anyway.
  real(prec),parameter :: SELF_AMR_DEFAULT_RELFLOOR = 1.0e-12_prec

  type :: RefinementIndicator3D_t
    integer :: N = 0
      !! Polynomial degree of the interpolant the indicator is built for.
    integer :: nElem = 0
      !! Number of (rank-local) elements the indicator arrays are sized for.
    real(prec) :: refineThreshold = 0.0_prec
      !! Elements with sigma_e above this value are flagged SELF_AMR_REFINE.
    real(prec) :: coarsenThreshold = 0.0_prec
      !! Elements with sigma_e below this value are flagged SELF_AMR_COARSEN.
    real(prec) :: relativeEnergyFloor = SELF_AMR_DEFAULT_RELFLOOR
      !! Elements whose gate energy is at or below relativeEnergyFloor*energyScale are treated as
      !! quiescent (hence perfectly resolved) regardless of their modal shape. Energy goes as
      !! amplitude squared, so this is 10**(dB/10) in amplitude terms: 1e-8 gates amplitudes
      !! below 1e-4 (-80 dB) of the field scale. Set to 0 to recover the pure absolute
      !! (machine-epsilon) floor.
    real(prec) :: significantEnergyFloor = SELF_AMR_DEFAULT_RELFLOOR
      !! Upper edge of the hysteresis band on the energy axis. Elements between
      !! relativeEnergyFloor and this fraction of the energy scale are flagged SELF_AMR_KEEP
      !! whatever their spectrum: too weak to justify spending levels on, too strong to declare
      !! resolved. Equal to relativeEnergyFloor by default, which collapses the band to a single
      !! hard cut. See SetRelativeEnergyFloor.
    real(prec) :: energyScale = 0.0_prec
      !! Squared field scale the relative floor is measured against, in the units of the gate
      !! energy. Meaningful only when energyScaleIsSet is true; otherwise the scale is recomputed
      !! from the current field on every Estimate.
    logical :: energyScaleIsSet = .false.
      !! Whether energyScale was pinned by SetEnergyScale (true) or is computed automatically as
      !! the largest gate energy over the elements (false, the default).
    logical :: energyWeightsSet = .false.
      !! Whether energyWeight was supplied by SetEnergyWeights (true) or is regenerated from the
      !! driving-variable index on every Estimate (false, the default).
    integer :: nVarWeights = 0
      !! Allocated length of energyWeight (the solution variable count it was resolved for).
    real(prec),pointer,contiguous,dimension(:) :: energyWeight => null()
      !! Non-negative weights w_v of the gate energy g_e = sum_v w_v E_tot,e,v, indexed by
      !! solution variable. Resolved lazily, because the variable count is a property of the
      !! solution field and is not known at Init.
    real(prec),pointer,contiguous,dimension(:) :: gate => null()
      !! Per-element gate energy g_e from the most recent Estimate. Retained as a diagnostic: it
      !! is the quantity the amplitude gate actually compared against the effective floor.
    real(prec),pointer,contiguous,dimension(:,:) :: Pmodal => null()
      !! Nodal-to-modal transform. Pmodal(ii,p) is the (p,ii) entry of the inverse Legendre
      !! Vandermonde in the L2-normalized basis, so the 1-D modal coefficients are
      !! uhat(p) = sum_ii Pmodal(ii,p) * u(ii) (first index summed, SELF matrix convention).
    real(prec),pointer,contiguous,dimension(:) :: indicator => null()
      !! Per-element indicator value sigma_e = log10(S_e).
    integer,pointer,contiguous,dimension(:) :: flag => null()
      !! Per-element refinement flag: SELF_AMR_REFINE / SELF_AMR_KEEP / SELF_AMR_COARSEN.

  contains

    procedure,public :: Init => Init_RefinementIndicator3D_t
    procedure,public :: Free => Free_RefinementIndicator3D_t
    procedure,public :: SetThresholds => SetThresholds_RefinementIndicator3D_t
    procedure,public :: SetRelativeEnergyFloor => SetRelativeEnergyFloor_RefinementIndicator3D_t
    procedure,public :: SetEnergyScale => SetEnergyScale_RefinementIndicator3D_t
    procedure,public :: ClearEnergyScale => ClearEnergyScale_RefinementIndicator3D_t
    procedure,public :: SetEnergyWeights => SetEnergyWeights_RefinementIndicator3D_t
    procedure,public :: UpdateHost => UpdateHost_RefinementIndicator3D_t
    procedure,public :: UpdateDevice => UpdateDevice_RefinementIndicator3D_t
    procedure,public :: Estimate => Estimate_RefinementIndicator3D_t
    procedure,public :: CountFlagged => CountFlagged_RefinementIndicator3D_t

  endtype RefinementIndicator3D_t

contains

  subroutine Init_RefinementIndicator3D_t(this,interp,nElem,refineThreshold,coarsenThreshold)
    !! Allocate the indicator for an interpolant of degree interp%N and nElem elements and
    !! precompute the nodal->modal transform matrix from the interpolant control points.
    !!
    !! The amplitude gate is left at its defaults here (relativeEnergyFloor =
    !! SELF_AMR_DEFAULT_RELFLOOR, automatic energy scale, unit weights); a caller that has tuned
    !! those must re-apply them after any re-Init, since intent(out) resets them.
    implicit none
    class(RefinementIndicator3D_t),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nElem
    real(prec),intent(in) :: refineThreshold
    real(prec),intent(in) :: coarsenThreshold

    if(interp%N < 1) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : RefinementIndicator3D requires interpolant degree N >= 1, got ',interp%N
      stop 1
    endif

    if(refineThreshold <= coarsenThreshold) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : refineThreshold must be greater than coarsenThreshold.'
      stop 1
    endif

    this%N = interp%N
    this%nElem = nElem
    this%refineThreshold = refineThreshold
    this%coarsenThreshold = coarsenThreshold

    allocate(this%Pmodal(1:interp%N+1,1:interp%N+1))
    allocate(this%indicator(1:nElem))
    allocate(this%flag(1:nElem))
    allocate(this%gate(1:nElem))

    call BuildModalTransform(interp%controlPoints,interp%N,this%Pmodal)

    this%indicator = 0.0_prec
    this%flag = SELF_AMR_KEEP
    this%gate = 0.0_prec

  endsubroutine Init_RefinementIndicator3D_t

  subroutine Free_RefinementIndicator3D_t(this)
    implicit none
    class(RefinementIndicator3D_t),intent(inout) :: this

    this%N = 0
    this%nElem = 0
    if(associated(this%Pmodal)) deallocate(this%Pmodal)
    if(associated(this%indicator)) deallocate(this%indicator)
    if(associated(this%flag)) deallocate(this%flag)
    if(associated(this%gate)) deallocate(this%gate)
    if(associated(this%energyWeight)) deallocate(this%energyWeight)
    this%Pmodal => null()
    this%indicator => null()
    this%flag => null()
    this%gate => null()
    this%energyWeight => null()
    this%nVarWeights = 0
    this%energyWeightsSet = .false.
    this%energyScaleIsSet = .false.
    this%energyScale = 0.0_prec

  endsubroutine Free_RefinementIndicator3D_t

  subroutine SetThresholds_RefinementIndicator3D_t(this,refineThreshold,coarsenThreshold)
    !! Update the refine/coarsen thresholds without rebuilding the transform.
    implicit none
    class(RefinementIndicator3D_t),intent(inout) :: this
    real(prec),intent(in) :: refineThreshold
    real(prec),intent(in) :: coarsenThreshold

    if(refineThreshold <= coarsenThreshold) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : refineThreshold must be greater than coarsenThreshold.'
      stop 1
    endif
    this%refineThreshold = refineThreshold
    this%coarsenThreshold = coarsenThreshold

  endsubroutine SetThresholds_RefinementIndicator3D_t

  subroutine SetRelativeEnergyFloor_RefinementIndicator3D_t(this,relativeEnergyFloor, &
                                                            significantEnergyFloor)
    !! Set the relative energy floor of the amplitude gate. An element whose gate energy g
    !! satisfies
    !!
    !!   g <= max( epsilon(1.0_prec), relativeEnergyFloor*energyScale )
    !!
    !! carries no resolvable signal at the scale of the field and is reported perfectly resolved
    !! (sigma_e = log10(epsilon) -> SELF_AMR_COARSEN) whatever the shape of its modal spectrum.
    !!
    !! Dimensionless, and an ENERGY fraction, so it is the square of the corresponding amplitude
    !! fraction: 1e-12 gates amplitudes below 1e-6 (-120 dB) of the field scale, 1e-8 gates
    !! amplitudes below 1e-4 (-80 dB). Must lie in [0,1); 0 disables the relative floor and
    !! restores the pure absolute (machine-epsilon) guard.
    !!
    !! Raising this beyond the default trades against refinement depth - the gate cannot tell
    !! residue from the flank of an under-resolved feature - and has been measured to cost more
    !! than it saves on a propagating-wave problem. See SELF_AMR_DEFAULT_RELFLOOR.
    !!
    !! significantEnergyFloor, when supplied, is the UPPER edge of a hysteresis band on the energy
    !! axis, exactly analogous to the refine/coarsen band on sigma_e:
    !!
    !!   g <= quiescent floor                    -> COARSEN, whatever the spectrum says
    !!   quiescent floor < g <= significant      -> SELF_AMR_KEEP
    !!   g > significant floor                   -> the spectrum decides, as before
    !!
    !! Without it the amplitude gate is a single hard cut, so an element whose energy drifts across
    !! that one value flips COARSEN <-> REFINE on successive epochs - the thrashing the two sigma
    !! thresholds exist to prevent, reintroduced on the other axis. The middle zone says "too weak
    !! to be worth spending levels on, too strong to declare resolved", which is the honest answer
    !! there and is stable under a small change in amplitude.
    !!
    !! Must satisfy quiescent floor <= significant floor < 1. It defaults to the quiescent floor,
    !! which collapses the band to the single cut and reproduces the ungapped behaviour exactly.
    !!
    !! KEEP THE BAND NARROW. Its upper edge is functionally "do not spend levels below this
    !! energy", which is the same knob as the floor itself and carries the same cost in refinement
    !! depth. Measured on the ultrasound benchmark's initial adaptation (see
    !! SELF_AMR_DEFAULT_RELFLOOR), with a quiescent floor of 1e-12:
    !!
    !!   significant floor   1e-12   3e-12   1e-11   3e-11   1e-10
    !!   elements / level    328/2   328/2   268/1   268/1   268/1
    !!
    !! A band up to ~3x the floor leaves refinement depth untouched; at 10x it collapses a level,
    !! exactly as raising the floor to 1e-10 does. The band is a thrash damper, not a savings knob:
    !! widen it only as far as is needed to stop flags oscillating.
    implicit none
    class(RefinementIndicator3D_t),intent(inout) :: this
    real(prec),intent(in) :: relativeEnergyFloor
    real(prec),intent(in),optional :: significantEnergyFloor

    if(relativeEnergyFloor < 0.0_prec .or. relativeEnergyFloor >= 1.0_prec) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : relativeEnergyFloor must satisfy 0 <= floor < 1.'
      stop 1
    endif
    this%relativeEnergyFloor = relativeEnergyFloor
    this%significantEnergyFloor = relativeEnergyFloor
    if(present(significantEnergyFloor)) then
      if(significantEnergyFloor < relativeEnergyFloor .or. &
         significantEnergyFloor >= 1.0_prec) then
        print*,__FILE__,':',__LINE__, &
          ' : Error : significantEnergyFloor must satisfy relativeEnergyFloor <= floor < 1.'
        stop 1
      endif
      this%significantEnergyFloor = significantEnergyFloor
    endif

  endsubroutine SetRelativeEnergyFloor_RefinementIndicator3D_t

  subroutine SetEnergyScale_RefinementIndicator3D_t(this,energyScale)
    !! Pin the energy scale that normalizes the relative floor, in the units of the gate energy
    !! (squared field amplitude times the reference-element volume). Two reasons to use this rather
    !! than the automatic maximum over elements:
    !!
    !!   - determinism: the automatic scale is a floating-point reduction, and under MPI it is a
    !!     collective whose result can differ at round-off between rank counts, which propagates
    !!     into the flags;
    !!   - a poor automatic normalizer: a strong steady background would otherwise set the floor
    !!     for a weak transient that is the feature of interest.
    !!
    !! Must be >= 0; 0 makes the gate collapse onto the absolute floor.
    implicit none
    class(RefinementIndicator3D_t),intent(inout) :: this
    real(prec),intent(in) :: energyScale

    if(energyScale < 0.0_prec) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : energyScale must be non-negative.'
      stop 1
    endif
    this%energyScale = energyScale
    this%energyScaleIsSet = .true.

  endsubroutine SetEnergyScale_RefinementIndicator3D_t

  subroutine ClearEnergyScale_RefinementIndicator3D_t(this)
    !! Return to the automatic energy scale (the largest gate energy over the elements).
    implicit none
    class(RefinementIndicator3D_t),intent(inout) :: this

    this%energyScale = 0.0_prec
    this%energyScaleIsSet = .false.

  endsubroutine ClearEnergyScale_RefinementIndicator3D_t

  subroutine SetEnergyWeights_RefinementIndicator3D_t(this,w)
    !! Set the per-variable weights w_v >= 0 of the gate energy g_e = sum_v w_v E_tot,e,v, where
    !! E_tot,e,v is the exact L2 energy of variable v on the reference element.
    !!
    !! E_tot is a convex quadratic functional of the state, so g_e with these weights is a
    !! discrete quadratic entropy integral over the element: taking w from the diagonal of the
    !! entropy Hessian makes the gate an entropy (energy) measure rather than a raw sum of
    !! squared variables, which is what makes it meaningful for a system whose variables carry
    !! different units and magnitudes. For LinearEuler3D, whose entropy density is
    !! 0.5*rho0*(u^2 + v^2 + w^2) + 0.5*P^2/(rho0 c^2) and whose variables 5 and 6 are the
    !! time-constant background fields c and rho0,
    !!
    !!   w = [ 0.5*rho0, 0.5*rho0, 0.5*rho0, 0.5/(rho0*c0**2), 0.0, 0.0 ]
    !!
    !! is the entropy-weighted gate, and the zero weights keep the (large) background fields from
    !! setting the scale.
    !!
    !! size(w) must equal the solution's nVar at Estimate time. Every weight must be >= 0 and at
    !! least one must be > 0, else no element would ever clear the gate. Note that a weight on a
    !! variable outside the driving-variable selection makes Estimate transform that variable too
    !! (it is needed for the gate), which costs one extra modal transform per element per such
    !! variable.
    implicit none
    class(RefinementIndicator3D_t),intent(inout) :: this
    real(prec),intent(in) :: w(:)
    ! Local
    integer :: v
    logical :: anyPositive

    if(size(w) < 1) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : SetEnergyWeights requires at least one weight.'
      stop 1
    endif
    anyPositive = .false.
    do v = 1,size(w)
      if(w(v) < 0.0_prec) then
        print*,__FILE__,':',__LINE__, &
          ' : Error : energy weights must be non-negative.'
        stop 1
      endif
      if(w(v) > 0.0_prec) anyPositive = .true.
    enddo
    if(.not. anyPositive) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : at least one energy weight must be positive.'
      stop 1
    endif

    if(this%nVarWeights /= size(w)) then
      if(associated(this%energyWeight)) deallocate(this%energyWeight)
      allocate(this%energyWeight(1:size(w)))
      this%nVarWeights = size(w)
    endif
    do v = 1,size(w)
      this%energyWeight(v) = w(v)
    enddo
    this%energyWeightsSet = .true.

  endsubroutine SetEnergyWeights_RefinementIndicator3D_t

  subroutine ResolveEnergyWeights(this,nVar,ivar)
    !! Ensure energyWeight is associated and sized to nVar. Weights supplied by SetEnergyWeights
    !! must match nVar; otherwise the default weights are regenerated: unit weight on the driving
    !! variable, or on every variable when ivar is SELF_AMR_ALLVARS. Resolved here rather than at
    !! Init because the variable count belongs to the solution field, not to the indicator.
    implicit none
    class(RefinementIndicator3D_t),intent(inout) :: this
    integer,intent(in) :: nVar
    integer,intent(in) :: ivar
    ! Local
    integer :: v

    if(this%energyWeightsSet) then
      if(this%nVarWeights /= nVar) then
        print*,__FILE__,':',__LINE__, &
          ' : Error : energy weights were set for a different variable count.'
        stop 1
      endif
      return
    endif

    if(this%nVarWeights /= nVar) then
      if(associated(this%energyWeight)) deallocate(this%energyWeight)
      allocate(this%energyWeight(1:nVar))
      this%nVarWeights = nVar
    endif
    do v = 1,nVar
      if(ivar == SELF_AMR_ALLVARS .or. v == ivar) then
        this%energyWeight(v) = 1.0_prec
      else
        this%energyWeight(v) = 0.0_prec
      endif
    enddo

  endsubroutine ResolveEnergyWeights

  subroutine ResolveEnergyScale(this,energyScale,comm)
    !! Determine the energy scale the relative floor is measured against: the pinned value when
    !! SetEnergyScale was used, otherwise the largest gate energy over the (rank-local) elements.
    !! When comm is present the maximum is taken over the communicator, so every rank applies the
    !! same floor and the flags - hence the adapted mesh - do not depend on the decomposition.
    !! One small collective per indicator evaluation, i.e. per adaptation epoch, never inside the
    !! time-stepping loop.
    implicit none
    class(RefinementIndicator3D_t),intent(in) :: this
    real(prec),intent(out) :: energyScale
    integer,intent(in),optional :: comm
    ! Local
    integer :: iel,ierror,mpiPrec
    real(prec) :: gmax

    if(this%energyScaleIsSet) then
      energyScale = this%energyScale
      return
    endif

    gmax = 0.0_prec
    do iel = 1,this%nElem
      gmax = max(gmax,this%gate(iel))
    enddo

    if(present(comm)) then
      if(prec == real32) then
        mpiPrec = MPI_FLOAT
      else
        mpiPrec = MPI_DOUBLE
      endif
      call mpi_allreduce(gmax,energyScale,1,mpiPrec,MPI_MAX,comm,ierror)
    else
      energyScale = gmax
    endif

  endsubroutine ResolveEnergyScale

  subroutine FinalizeIndicator(this,energyScale)
    !! Second phase of an estimate, shared by every backend so that the flags are identical
    !! whichever computed the spectra: apply the amplitude gate, take the log10 and set the
    !! refine/keep/coarsen flags.
    !!
    !! On entry indicator(iel) holds the RAW smoothness ratio S_e from the first phase and
    !! gate(iel) the element's gate energy; on exit indicator(iel) = log10(max(S_e,epsilon)) with
    !! S_e forced to zero on quiescent elements. Nothing outside Estimate observes the intermediate
    !! state.
    !!
    !! Elements inside the energy hysteresis band keep their true spectral sigma_e - it stays a
    !! faithful diagnostic of the spectrum - but have their flag forced to SELF_AMR_KEEP. So for
    !! those elements the flag is deliberately NOT the value thresholding sigma_e would give; that
    !! is the whole point of the band.
    implicit none
    class(RefinementIndicator3D_t),intent(inout) :: this
    real(prec),intent(in) :: energyScale
    ! Local
    integer :: iel
    real(prec) :: effFloor,effSignificant,se

    ! The absolute term is the safety net for a field that is identically zero; the relative term
    ! is what releases low-amplitude (but far-above-epsilon) residue behind a passing front.
    effFloor = max(epsilon(1.0_prec),this%relativeEnergyFloor*energyScale)
    ! Upper edge of the band; never below the lower edge, so a degenerate band stays degenerate
    ! even when the absolute term is what sets the lower edge.
    effSignificant = max(effFloor,this%significantEnergyFloor*energyScale)

    do iel = 1,this%nElem
      if(this%gate(iel) <= effFloor) then
        se = 0.0_prec ! quiescent at the scale of the field: perfectly resolved
      else
        se = this%indicator(iel)
      endif
      this%indicator(iel) = log10(max(se,epsilon(1.0_prec)))
      if(this%gate(iel) > effFloor .and. this%gate(iel) <= effSignificant) then
        ! Inside the energy hysteresis band: too weak to justify spending levels on, too strong to
        ! declare resolved. Holding the mesh here is what keeps an element whose amplitude drifts
        ! across the gate from flipping REFINE <-> COARSEN on successive epochs.
        this%flag(iel) = SELF_AMR_KEEP
      elseif(this%indicator(iel) > this%refineThreshold) then
        this%flag(iel) = SELF_AMR_REFINE
      elseif(this%indicator(iel) < this%coarsenThreshold) then
        this%flag(iel) = SELF_AMR_COARSEN
      else
        this%flag(iel) = SELF_AMR_KEEP
      endif
    enddo

  endsubroutine FinalizeIndicator

  subroutine CheckEstimateArguments(this,solution,ivar,gate)
    !! Argument validation shared by the portable and backend Estimate implementations.
    implicit none
    class(RefinementIndicator3D_t),intent(in) :: this
    class(Scalar3D),intent(in) :: solution
    integer,intent(in) :: ivar
    real(prec),intent(in),optional :: gate(:)

    if(solution%N /= this%N) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : solution degree does not match indicator degree.'
      stop 1
    endif
    if(solution%nElem /= this%nElem) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : solution element count does not match indicator element count.'
      stop 1
    endif
    if(ivar < 0 .or. ivar > solution%nVar) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : driving-variable index out of range.'
      stop 1
    endif
    if(present(gate)) then
      if(size(gate) < this%nElem) then
        print*,__FILE__,':',__LINE__, &
          ' : Error : caller-supplied gate array is smaller than the element count.'
        stop 1
      endif
    endif

  endsubroutine CheckEstimateArguments

  subroutine UpdateHost_RefinementIndicator3D_t(this)
    implicit none
    class(RefinementIndicator3D_t),intent(inout) :: this
    if(.false.) this%N = this%N ! CPU stub; suppress unused-dummy-argument warning
  endsubroutine UpdateHost_RefinementIndicator3D_t

  subroutine UpdateDevice_RefinementIndicator3D_t(this)
    implicit none
    class(RefinementIndicator3D_t),intent(inout) :: this
    if(.false.) this%N = this%N ! CPU stub; suppress unused-dummy-argument warning
  endsubroutine UpdateDevice_RefinementIndicator3D_t

  subroutine Estimate_RefinementIndicator3D_t(this,solution,ivar,comm,gate)
    !! Compute the per-element modal-energy indicator sigma_e and refine/keep/coarsen flag from
    !! the nodal solution field. ivar selects the driving variable in [1,solution%nVar]; passing
    !! SELF_AMR_ALLVARS (=0) reduces the indicator over all variables by taking, per element, the
    !! largest (least smooth) smoothness ratio S_e before the log10.
    !!
    !! The estimate runs in two phases. The first is element-local and parallel: it forms each
    !! element's modal spectrum, its raw smoothness ratio S_e, and its gate energy
    !! g = sum_v w_v E_tot,v. The second applies the amplitude gate, the log10 and the thresholds,
    !! and is deferred because the gate needs a field-wide energy scale (a reduction, and under
    !! MPI a collective) that no element-local pass can supply. Between the phases indicator(:)
    !! transiently holds the raw ratio rather than its logarithm.
    !!
    !! comm, when present, is the MPI communicator over which the automatic energy scale is
    !! maximized, so that the flags do not depend on the domain decomposition. gate, when present,
    !! replaces the weighted-sum gate energy with a caller-computed per-element value - for
    !! example an exactly integrated nodal entropy - and is used both for the scale and for the
    !! gate comparison.
    implicit none
    class(RefinementIndicator3D_t),intent(inout) :: this
    class(Scalar3D),intent(in) :: solution
    integer,intent(in) :: ivar
    integer,intent(in),optional :: comm
    real(prec),intent(in),optional :: gate(:)
    ! Local
    integer :: iel
    real(prec) :: energyScale

    call CheckEstimateArguments(this,solution,ivar,gate)
    call ResolveEnergyWeights(this,solution%nVar,ivar)

    do concurrent(iel=1:this%nElem)
      block
        integer :: j,k,p,q,rr,ii,v,v0,v1,Np
        real(prec) :: tmp(1:this%N+1,1:this%N+1,1:this%N+1)
        real(prec) :: uhat(1:this%N+1,1:this%N+1,1:this%N+1)
        real(prec) :: acc,etot,eclip1,eclip2,r1,r2,se,semax,g,w
        logical :: needSe,needG
        ! Absolute (machine-epsilon) guard on the energy RATIOS below. This is not the amplitude
        ! gate: it only keeps 0/0 out of r1 and r2. The amplitude gate is relative and is applied
        ! in FinalizeIndicator, once the field-wide energy scale is known.
        real(prec),parameter :: energyFloor = epsilon(1.0_prec)

        Np = this%N+1
        if(ivar == SELF_AMR_ALLVARS) then
          v0 = 1
          v1 = solution%nVar
        else
          v0 = ivar
          v1 = ivar
        endif

        semax = 0.0_prec
        g = 0.0_prec
        do v = 1,solution%nVar
          ! A variable is transformed if it drives the indicator, or if it carries a non-zero gate
          ! weight (its energy is then needed for the gate). With the default weights these
          ! coincide, so the work is exactly what it was before the gate existed.
          w = this%energyWeight(v)
          needSe = (v >= v0 .and. v <= v1)
          needG = (w /= 0.0_prec)
          if(.not.(needSe .or. needG)) cycle

          ! Pass 1 (xi direction): tmp(p,j,k) = sum_i Pmodal(i,p) * u(i,j,k)
          do k = 1,Np
            do j = 1,Np
              do p = 1,Np
                acc = 0.0_prec
                do ii = 1,Np
                  acc = acc+this%Pmodal(ii,p)*solution%interior(ii,j,k,iel,v)
                enddo
                tmp(p,j,k) = acc
              enddo
            enddo
          enddo
          ! Pass 2 (eta direction): uhat(p,q,k) = sum_j Pmodal(j,q) * tmp(p,j,k)
          do k = 1,Np
            do q = 1,Np
              do p = 1,Np
                acc = 0.0_prec
                do ii = 1,Np
                  acc = acc+this%Pmodal(ii,q)*tmp(p,ii,k)
                enddo
                uhat(p,q,k) = acc
              enddo
            enddo
          enddo
          ! Pass 3 (zeta direction): tmp(p,q,r) = sum_k Pmodal(k,r) * uhat(p,q,k), so tmp holds
          ! the full 3-D modal coefficients uhat(p,q,r) after this pass. The two buffers alternate
          ! (tmp -> uhat -> tmp) so only two element-sized temporaries are needed, mirroring the
          ! shared-memory budget of the GPU kernel.
          do rr = 1,Np
            do q = 1,Np
              do p = 1,Np
                acc = 0.0_prec
                do ii = 1,Np
                  acc = acc+this%Pmodal(ii,rr)*uhat(p,q,ii)
                enddo
                tmp(p,q,rr) = acc
              enddo
            enddo
          enddo

          ! Modal energies. eclip1 drops the highest mode in each direction,
          ! eclip2 drops the highest two (only meaningful for N >= 2).
          etot = 0.0_prec
          eclip1 = 0.0_prec
          eclip2 = 0.0_prec
          do rr = 1,Np
            do q = 1,Np
              do p = 1,Np
                etot = etot+tmp(p,q,rr)*tmp(p,q,rr)
                if(p <= Np-1 .and. q <= Np-1 .and. rr <= Np-1) then
                  eclip1 = eclip1+tmp(p,q,rr)*tmp(p,q,rr)
                endif
                if(p <= Np-2 .and. q <= Np-2 .and. rr <= Np-2) then
                  eclip2 = eclip2+tmp(p,q,rr)*tmp(p,q,rr)
                endif
              enddo
            enddo
          enddo

          ! Gate energy: g = sum_v w_v E_tot,v, the discrete (quadratic) entropy integral over
          ! the element when the weights come from an entropy Hessian.
          if(needG) g = g+w*etot

          if(needSe) then
            if(etot <= energyFloor) then
              ! Field is (near) identically zero on this element: perfectly resolved.
              se = 0.0_prec
            else
              r1 = (etot-eclip1)/etot
              if(this%N >= 2 .and. eclip1 > energyFloor) then
                r2 = (eclip1-eclip2)/eclip1
              else
                r2 = 0.0_prec
              endif
              se = max(r1,r2)
            endif

            semax = max(semax,se)
          endif
        enddo

        ! Raw ratio; FinalizeIndicator gates it and takes the log10 in place.
        this%indicator(iel) = semax
        this%gate(iel) = g
      endblock
    enddo

    if(present(gate)) then
      do iel = 1,this%nElem
        this%gate(iel) = gate(iel)
      enddo
    endif

    call ResolveEnergyScale(this,energyScale,comm)
    call FinalizeIndicator(this,energyScale)

  endsubroutine Estimate_RefinementIndicator3D_t

  function CountFlagged_RefinementIndicator3D_t(this,flagValue) result(n)
    !! Count the rank-local elements currently carrying the requested flag value
    !! (SELF_AMR_REFINE / SELF_AMR_KEEP / SELF_AMR_COARSEN). Provided as a convenience for
    !! drivers and tests; a global count across MPI ranks is the caller's responsibility.
    implicit none
    class(RefinementIndicator3D_t),intent(in) :: this
    integer,intent(in) :: flagValue
    integer :: n
    ! Local
    integer :: iel

    n = 0
    do iel = 1,this%nElem
      if(this%flag(iel) == flagValue) n = n+1
    enddo

  endfunction CountFlagged_RefinementIndicator3D_t

  subroutine BuildModalTransform(controlPoints,N,Pmodal)
    !! Build the nodal->modal transform Pmodal for a 1-D degree-N interpolant whose nodes are
    !! controlPoints(1:N+1). The transform is the exact inverse of the L2-normalized Legendre
    !! Vandermonde V(i,p) = Ltilde_{p-1}(x_i); Pmodal(ii,p) stores V^{-1}(p,ii) so that
    !! uhat(p) = sum_ii Pmodal(ii,p) * u(ii). The Vandermonde inverse is formed in double
    !! precision for conditioning and cast back to the working precision.
    implicit none
    integer,intent(in) :: N
    real(prec),intent(in) :: controlPoints(1:N+1)
    real(prec),intent(out) :: Pmodal(1:N+1,1:N+1)
    ! Local
    integer :: i,p
    real(real64) :: V(1:N+1,1:N+1)
    real(real64) :: Vinv(1:N+1,1:N+1)

    ! Vandermonde in the normalized Legendre basis: column p (degree p-1) evaluated at node i.
    do i = 1,N+1
      do p = 1,N+1
        V(i,p) = NormalizedLegendre(p-1,real(controlPoints(i),real64))
      enddo
    enddo

    call InvertMatrix(V,Vinv,N+1)

    ! Store the transpose of the inverse so the summed (input node) index comes first,
    ! matching the SELF matrix-application convention.
    do p = 1,N+1
      do i = 1,N+1
        Pmodal(i,p) = real(Vinv(p,i),prec)
      enddo
    enddo

  endsubroutine BuildModalTransform

  pure function NormalizedLegendre(p,x) result(Lp)
    !! L2-normalized Legendre polynomial Ltilde_p(x) = L_p(x)*sqrt((2p+1)/2) on [-1,1],
    !! evaluated with the standard three-term recurrence. The normalization gives
    !! int_{-1}^{1} Ltilde_p Ltilde_q dx = delta_pq.
    implicit none
    integer,intent(in) :: p
    real(real64),intent(in) :: x
    real(real64) :: Lp
    ! Local
    integer :: k
    real(real64) :: lkm1,lk,lkp1

    if(p == 0) then
      Lp = 1.0_real64
    elseif(p == 1) then
      Lp = x
    else
      lkm1 = 1.0_real64
      lk = x
      do k = 1,p-1
        lkp1 = (real(2*k+1,real64)*x*lk-real(k,real64)*lkm1)/real(k+1,real64)
        lkm1 = lk
        lk = lkp1
      enddo
      Lp = lk
    endif

    Lp = Lp*sqrt((2.0_real64*real(p,real64)+1.0_real64)/2.0_real64)

  endfunction NormalizedLegendre

  subroutine InvertMatrix(A,Ainv,n)
    !! Invert the n-by-n matrix A by Gauss-Jordan elimination with partial pivoting, in double
    !! precision. Used once at initialization for the small (N+1) Legendre Vandermonde; A is a
    !! well-conditioned Vandermonde in an orthonormal basis, so a direct solve is appropriate.
    implicit none
    integer,intent(in) :: n
    real(real64),intent(in) :: A(1:n,1:n)
    real(real64),intent(out) :: Ainv(1:n,1:n)
    ! Local
    real(real64) :: M(1:n,1:n)
    integer :: i,j,k,piv
    real(real64) :: pmax,factor,tmp

    M = A
    Ainv = 0.0_real64
    do i = 1,n
      Ainv(i,i) = 1.0_real64
    enddo

    do k = 1,n
      ! Partial pivot: find the largest-magnitude entry in column k at or below the diagonal.
      piv = k
      pmax = abs(M(k,k))
      do i = k+1,n
        if(abs(M(i,k)) > pmax) then
          pmax = abs(M(i,k))
          piv = i
        endif
      enddo
      if(pmax <= 0.0_real64) then
        print*,__FILE__,':',__LINE__,' : Error : singular Vandermonde in modal transform.'
        stop 1
      endif
      if(piv /= k) then
        do j = 1,n
          tmp = M(k,j); M(k,j) = M(piv,j); M(piv,j) = tmp
          tmp = Ainv(k,j); Ainv(k,j) = Ainv(piv,j); Ainv(piv,j) = tmp
        enddo
      endif

      ! Normalize the pivot row.
      factor = M(k,k)
      do j = 1,n
        M(k,j) = M(k,j)/factor
        Ainv(k,j) = Ainv(k,j)/factor
      enddo

      ! Eliminate column k from every other row.
      do i = 1,n
        if(i /= k) then
          factor = M(i,k)
          do j = 1,n
            M(i,j) = M(i,j)-factor*M(k,j)
            Ainv(i,j) = Ainv(i,j)-factor*Ainv(k,j)
          enddo
        endif
      enddo
    enddo

  endsubroutine InvertMatrix

endmodule SELF_RefinementIndicator_3D_t
