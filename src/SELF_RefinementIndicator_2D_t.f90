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

module SELF_RefinementIndicator_2D_t
!! Modal (Legendre spectral-decay) refinement indicator for 2-D spectral element solutions.
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
!!   u(xi,eta) = sum_{p=0}^{N} sum_{q=0}^{N} uhat(p,q) Ltilde_p(xi) Ltilde_q(eta),
!!
!! where Ltilde_p are the L2-normalized Legendre polynomials on [-1,1]
!! (int_{-1}^{1} Ltilde_p Ltilde_q dxi = delta_pq). With that normalization the total modal
!! energy equals the exact L2 energy of the field on the reference element,
!!
!!   E_tot = sum_{p,q} uhat(p,q)^2 = ||u||_{L2([-1,1]^2)}^2 .
!!
!! Two "clipped" energies are formed by dropping the highest one and two modes in each direction,
!!
!!   E_clip1 = sum_{p,q <= N-1} uhat(p,q)^2 ,   E_clip2 = sum_{p,q <= N-2} uhat(p,q)^2 ,
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
!! The backend extension modules (SELF_RefinementIndicator_2D) add device storage and a GPU
!! kernel while preserving the mathematics bit-for-bit in structure.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Scalar_2D
  use iso_c_binding

  implicit none

  ! Per-element refinement flags returned in this%flag(:)
  integer,parameter :: SELF_AMR_COARSEN = -1
  integer,parameter :: SELF_AMR_KEEP = 0
  integer,parameter :: SELF_AMR_REFINE = 1

  ! Reduce the indicator over all variables when this sentinel is passed as the
  ! driving-variable index to Estimate.
  integer,parameter :: SELF_AMR_ALLVARS = 0

  type :: RefinementIndicator2D_t
    integer :: N = 0
      !! Polynomial degree of the interpolant the indicator is built for.
    integer :: nElem = 0
      !! Number of (rank-local) elements the indicator arrays are sized for.
    real(prec) :: refineThreshold = 0.0_prec
      !! Elements with sigma_e above this value are flagged SELF_AMR_REFINE.
    real(prec) :: coarsenThreshold = 0.0_prec
      !! Elements with sigma_e below this value are flagged SELF_AMR_COARSEN.
    real(prec),pointer,contiguous,dimension(:,:) :: Pmodal => null()
      !! Nodal-to-modal transform. Pmodal(ii,p) is the (p,ii) entry of the inverse Legendre
      !! Vandermonde in the L2-normalized basis, so the 1-D modal coefficients are
      !! uhat(p) = sum_ii Pmodal(ii,p) * u(ii) (first index summed, SELF matrix convention).
    real(prec),pointer,contiguous,dimension(:) :: indicator => null()
      !! Per-element indicator value sigma_e = log10(S_e).
    integer,pointer,contiguous,dimension(:) :: flag => null()
      !! Per-element refinement flag: SELF_AMR_REFINE / SELF_AMR_KEEP / SELF_AMR_COARSEN.

  contains

    procedure,public :: Init => Init_RefinementIndicator2D_t
    procedure,public :: Free => Free_RefinementIndicator2D_t
    procedure,public :: SetThresholds => SetThresholds_RefinementIndicator2D_t
    procedure,public :: UpdateHost => UpdateHost_RefinementIndicator2D_t
    procedure,public :: UpdateDevice => UpdateDevice_RefinementIndicator2D_t
    procedure,public :: Estimate => Estimate_RefinementIndicator2D_t
    procedure,public :: CountFlagged => CountFlagged_RefinementIndicator2D_t

  endtype RefinementIndicator2D_t

contains

  subroutine Init_RefinementIndicator2D_t(this,interp,nElem,refineThreshold,coarsenThreshold)
    !! Allocate the indicator for an interpolant of degree interp%N and nElem elements and
    !! precompute the nodal->modal transform matrix from the interpolant control points.
    implicit none
    class(RefinementIndicator2D_t),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nElem
    real(prec),intent(in) :: refineThreshold
    real(prec),intent(in) :: coarsenThreshold

    if(interp%N < 1) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : RefinementIndicator2D requires interpolant degree N >= 1, got ',interp%N
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

    call BuildModalTransform(interp%controlPoints,interp%N,this%Pmodal)

    this%indicator = 0.0_prec
    this%flag = SELF_AMR_KEEP

  endsubroutine Init_RefinementIndicator2D_t

  subroutine Free_RefinementIndicator2D_t(this)
    implicit none
    class(RefinementIndicator2D_t),intent(inout) :: this

    this%N = 0
    this%nElem = 0
    if(associated(this%Pmodal)) deallocate(this%Pmodal)
    if(associated(this%indicator)) deallocate(this%indicator)
    if(associated(this%flag)) deallocate(this%flag)
    this%Pmodal => null()
    this%indicator => null()
    this%flag => null()

  endsubroutine Free_RefinementIndicator2D_t

  subroutine SetThresholds_RefinementIndicator2D_t(this,refineThreshold,coarsenThreshold)
    !! Update the refine/coarsen thresholds without rebuilding the transform.
    implicit none
    class(RefinementIndicator2D_t),intent(inout) :: this
    real(prec),intent(in) :: refineThreshold
    real(prec),intent(in) :: coarsenThreshold

    if(refineThreshold <= coarsenThreshold) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : refineThreshold must be greater than coarsenThreshold.'
      stop 1
    endif
    this%refineThreshold = refineThreshold
    this%coarsenThreshold = coarsenThreshold

  endsubroutine SetThresholds_RefinementIndicator2D_t

  subroutine UpdateHost_RefinementIndicator2D_t(this)
    implicit none
    class(RefinementIndicator2D_t),intent(inout) :: this
    if(.false.) this%N = this%N ! CPU stub; suppress unused-dummy-argument warning
  endsubroutine UpdateHost_RefinementIndicator2D_t

  subroutine UpdateDevice_RefinementIndicator2D_t(this)
    implicit none
    class(RefinementIndicator2D_t),intent(inout) :: this
    if(.false.) this%N = this%N ! CPU stub; suppress unused-dummy-argument warning
  endsubroutine UpdateDevice_RefinementIndicator2D_t

  subroutine Estimate_RefinementIndicator2D_t(this,solution,ivar)
    !! Compute the per-element modal-energy indicator sigma_e and refine/keep/coarsen flag from
    !! the nodal solution field. ivar selects the driving variable in [1,solution%nVar]; passing
    !! SELF_AMR_ALLVARS (=0) reduces the indicator over all variables by taking, per element, the
    !! largest (least smooth) smoothness ratio S_e before the log10.
    implicit none
    class(RefinementIndicator2D_t),intent(inout) :: this
    class(Scalar2D),intent(in) :: solution
    integer,intent(in) :: ivar
    ! Local
    integer :: iel
    real(prec),parameter :: energyFloor = epsilon(1.0_prec)

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

    do concurrent(iel=1:this%nElem)
      block
        integer :: j,p,q,ii,v,v0,v1,Np
        real(prec) :: tmp(1:this%N+1,1:this%N+1)
        real(prec) :: uhat(1:this%N+1,1:this%N+1)
        real(prec) :: acc,etot,eclip1,eclip2,r1,r2,se,semax

        Np = this%N+1
        if(ivar == SELF_AMR_ALLVARS) then
          v0 = 1
          v1 = solution%nVar
        else
          v0 = ivar
          v1 = ivar
        endif

        semax = 0.0_prec
        do v = v0,v1
          ! Pass 1 (xi direction): tmp(p,j) = sum_i Pmodal(i,p) * u(i,j)
          do j = 1,Np
            do p = 1,Np
              acc = 0.0_prec
              do ii = 1,Np
                acc = acc+this%Pmodal(ii,p)*solution%interior(ii,j,iel,v)
              enddo
              tmp(p,j) = acc
            enddo
          enddo
          ! Pass 2 (eta direction): uhat(p,q) = sum_j Pmodal(j,q) * tmp(p,j)
          do q = 1,Np
            do p = 1,Np
              acc = 0.0_prec
              do ii = 1,Np
                acc = acc+this%Pmodal(ii,q)*tmp(p,ii)
              enddo
              uhat(p,q) = acc
            enddo
          enddo

          ! Modal energies. eclip1 drops the highest mode in each direction,
          ! eclip2 drops the highest two (only meaningful for N >= 2).
          etot = 0.0_prec
          eclip1 = 0.0_prec
          eclip2 = 0.0_prec
          do q = 1,Np
            do p = 1,Np
              etot = etot+uhat(p,q)*uhat(p,q)
              if(p <= Np-1 .and. q <= Np-1) eclip1 = eclip1+uhat(p,q)*uhat(p,q)
              if(p <= Np-2 .and. q <= Np-2) eclip2 = eclip2+uhat(p,q)*uhat(p,q)
            enddo
          enddo

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
        enddo

        this%indicator(iel) = log10(max(semax,energyFloor))
        if(this%indicator(iel) > this%refineThreshold) then
          this%flag(iel) = SELF_AMR_REFINE
        elseif(this%indicator(iel) < this%coarsenThreshold) then
          this%flag(iel) = SELF_AMR_COARSEN
        else
          this%flag(iel) = SELF_AMR_KEEP
        endif
      endblock
    enddo

  endsubroutine Estimate_RefinementIndicator2D_t

  function CountFlagged_RefinementIndicator2D_t(this,flagValue) result(n)
    !! Count the rank-local elements currently carrying the requested flag value
    !! (SELF_AMR_REFINE / SELF_AMR_KEEP / SELF_AMR_COARSEN). Provided as a convenience for
    !! drivers and tests; a global count across MPI ranks is the caller's responsibility.
    implicit none
    class(RefinementIndicator2D_t),intent(in) :: this
    integer,intent(in) :: flagValue
    integer :: n
    ! Local
    integer :: iel

    n = 0
    do iel = 1,this%nElem
      if(this%flag(iel) == flagValue) n = n+1
    enddo

  endfunction CountFlagged_RefinementIndicator2D_t

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

endmodule SELF_RefinementIndicator_2D_t
