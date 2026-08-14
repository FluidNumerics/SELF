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

program test

  implicit none
  integer :: exit_code

  exit_code = refinement_indicator_3d_relative_floor()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function refinement_indicator_3d_relative_floor() result(r)
    !! Validates the amplitude gate (relative energy floor) of the 3-D modal-decay refinement
    !! indicator, on both Legendre-Gauss and Legendre-Gauss-Lobatto nodes.
    !!
    !! The smoothness ratio S_e is a ratio and therefore scale-free: grid-scale ringing at 1e-6 of
    !! the field's amplitude has exactly the same S_e = 1 as a full-amplitude discontinuity. Under
    !! an absolute (machine-epsilon) guard alone, such low-amplitude residue reads as
    !! under-resolved forever and is never released, which is what pins the mesh behind a passing
    !! wave. The gate compares each element's energy against a fraction of the field's energy
    !! scale instead.
    !!
    !! Every element carries the SAME shape - the pure highest tensor mode
    !! Ltilde_N(xi)Ltilde_N(eta)Ltilde_N(zeta), for which S_e = 1 exactly - scaled to a different
    !! amplitude, so the flags can only differ through the amplitude gate. With the normalized
    !! basis the modal energy of an element of amplitude a is exactly a**2:
    !!
    !!   element 1 : a = 1      -> E = 1      (sets the automatic energy scale)
    !!   element 2 : a = 1e-2   -> E = 1e-4
    !!   element 3 : a = 0.5    -> E = 0.25
    !!   element 4 : a = 0      -> E = 0      (the absolute-floor path)
    !!   element 5 : a = 1e-7   -> E = 1e-14
    !!
    !! Element 5 sits two orders below the default floor and element 2 eight orders above it, so
    !! the default case has wide margins in real64. In real32 the default floor is below epsilon
    !! and the absolute term takes over, which gates element 5 (1e-14 < 1.2e-7) and leaves element
    !! 2 (1e-4) alone - the same verdicts, so the assertions hold in both precisions.
    !!
    !! Variable 2 is zero everywhere except on element 1, where it is a large constant. It stands
    !! in for a background field (LinearEuler3D carries the sound speed and background density as
    !! variables 5 and 6) whose magnitude would otherwise set the gate's energy scale for every
    !! other variable - which is what the per-variable gate weights exist to prevent.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Scalar_3D
    use SELF_RefinementIndicator_3D

    implicit none

    integer,parameter :: N = 7
    integer,parameter :: nel = 5
    integer,parameter :: nvar = 2
    real(prec),parameter :: refineThreshold = -3.0_prec
    real(prec),parameter :: coarsenThreshold = -8.0_prec
    ! Amplitudes of the identical top-mode shape on each element; the modal energy is the square.
    real(prec),parameter :: amp(1:nel) = [1.0_prec,1.0e-2_prec,0.5_prec,0.0_prec,1.0e-7_prec]
    ! A background-field amplitude on variable 2 of element 1, far above every element's
    ! variable-1 energy: 100**2 * 8 = 8e4 against a largest variable-1 energy of 1.
    real(prec),parameter :: background = 100.0_prec
    ! Test floor: 1e-2 of the energy scale (1), so elements 2, 4 and 5 are gated and 1 and 3 are
    ! not, with two orders of magnitude of margin either side in both real32 and real64.
    real(prec),parameter :: testFloor = 1.0e-2_prec
    type(Lagrange),target :: interp
    type(Scalar3D) :: sol
    type(RefinementIndicator3D) :: amr
    integer :: nodeType,i,j,k,iel
    real(prec) :: x,y,z
    real(prec) :: gateOverride(1:nel)

    r = 0

    do nodeType = 1,2 ! GAUSS = 1, GAUSS_LOBATTO = 2

      call interp%Init(N=N,controlNodeType=nodeType,M=N,targetNodeType=UNIFORM)
      call sol%Init(interp,nvar,nel)

      do iel = 1,nel
        do k = 1,N+1
          do j = 1,N+1
            do i = 1,N+1
              x = interp%controlPoints(i)
              y = interp%controlPoints(j)
              z = interp%controlPoints(k)
              sol%interior(i,j,k,iel,1) = amp(iel)*nlegendre(N,x)*nlegendre(N,y)*nlegendre(N,z)
              if(iel == 1) then
                sol%interior(i,j,k,iel,2) = background
              else
                sol%interior(i,j,k,iel,2) = 0.0_prec
              endif
            enddo
          enddo
        enddo
      enddo
      call sol%UpdateDevice() ! no-op on CPU; pushes interior to device on GPU

      call amr%Init(interp,nel,refineThreshold,coarsenThreshold)

      ! ---- 1. The shipped default floor ----
      ! Energy scale 1, so the effective floor is SELF_AMR_DEFAULT_RELFLOOR: element 5 (1e-14) is
      ! quiescent, element 2 (1e-4) is not. This is the reported bug in miniature - without the
      ! relative floor element 5 sits ~1e2 above machine epsilon in real64 and is flagged REFINE
      ! forever.
      call amr%Estimate(sol,ivar=1)
      call ReportFlags(amr,nel,nodeType,"default floor")
      if(amr%flag(5) /= SELF_AMR_COARSEN) then
        print*,"FAIL: negligible-amplitude element not gated by the default relative floor"
        r = 1
      endif
      if(amr%flag(1) /= SELF_AMR_REFINE .or. amr%flag(2) /= SELF_AMR_REFINE .or. &
         amr%flag(3) /= SELF_AMR_REFINE) then
        print*,"FAIL: the default relative floor suppressed a resolvable feature"
        r = 1
      endif
      if(amr%flag(4) /= SELF_AMR_COARSEN) then
        print*,"FAIL: identically zero element not flagged COARSEN"
        r = 1
      endif

      ! ---- 2. An explicit floor gates exactly the elements below it ----
      call amr%SetRelativeEnergyFloor(testFloor)
      call amr%Estimate(sol,ivar=1)
      call ReportFlags(amr,nel,nodeType,"floor = 1e-2")
      if(amr%flag(1) /= SELF_AMR_REFINE .or. amr%flag(3) /= SELF_AMR_REFINE) then
        print*,"FAIL: elements above the floor should still be flagged REFINE"
        r = 1
      endif
      if(amr%flag(2) /= SELF_AMR_COARSEN .or. amr%flag(5) /= SELF_AMR_COARSEN) then
        print*,"FAIL: elements below the floor should be flagged COARSEN"
        r = 1
      endif
      ! Flag tally: elements 1 and 3 refine, 2, 4 and 5 coarsen, none kept.
      if(amr%CountFlagged(SELF_AMR_REFINE) /= 2 .or. &
         amr%CountFlagged(SELF_AMR_COARSEN) /= 3 .or. &
         amr%CountFlagged(SELF_AMR_KEEP) /= 0) then
        print*,"FAIL: unexpected flag tally under the relative floor"
        r = 1
      endif

      ! ---- 3. Disabling the relative floor restores the previous behaviour ----
      ! Confirms the gate (and not some other change) is what releases element 2: on shape alone
      ! it is indistinguishable from element 1.
      call amr%SetRelativeEnergyFloor(0.0_prec)
      call amr%Estimate(sol,ivar=1)
      call ReportFlags(amr,nel,nodeType,"floor = 0 (absolute only)")
      if(amr%flag(2) /= SELF_AMR_REFINE) then
        print*,"FAIL: with the relative floor disabled, element 2 should be REFINE again"
        r = 1
      endif

      ! ---- 4. A pinned energy scale replaces the automatic one ----
      call amr%SetRelativeEnergyFloor(testFloor)
      call amr%SetEnergyScale(1.0e6_prec) ! effective floor 1e4, above every element's energy
      call amr%Estimate(sol,ivar=1)
      call ReportFlags(amr,nel,nodeType,"pinned scale 1e6")
      if(amr%CountFlagged(SELF_AMR_COARSEN) /= nel) then
        print*,"FAIL: a pinned scale above every element energy should gate every element"
        r = 1
      endif
      call amr%ClearEnergyScale()
      call amr%Estimate(sol,ivar=1)
      if(amr%flag(3) /= SELF_AMR_REFINE .or. amr%flag(2) /= SELF_AMR_COARSEN) then
        print*,"FAIL: ClearEnergyScale did not restore the automatic scale"
        r = 1
      endif

      ! ---- 5. A caller-supplied gate replaces the weighted-sum gate ----
      ! Uniform gate energy 1 against a scale of 1: nothing is quiescent.
      gateOverride = 1.0_prec
      call amr%Estimate(sol,ivar=1,gate=gateOverride)
      call ReportFlags(amr,nel,nodeType,"caller gate = 1")
      if(amr%flag(2) /= SELF_AMR_REFINE) then
        print*,"FAIL: caller-supplied gate above the floor should not gate element 2"
        r = 1
      endif
      ! A zero gate is quiescent everywhere, whatever the spectra say.
      gateOverride = 0.0_prec
      call amr%Estimate(sol,ivar=1,gate=gateOverride)
      call ReportFlags(amr,nel,nodeType,"caller gate = 0")
      if(amr%CountFlagged(SELF_AMR_COARSEN) /= nel) then
        print*,"FAIL: an identically zero caller gate should gate every element"
        r = 1
      endif

      ! ---- 6. Per-variable gate weights ----
      ! Reducing over all variables with unit weights sums energies of different magnitudes, so
      ! the background field on element 1 (8e4) sets the scale and drags the effective floor up to
      ! 8e2 - swamping element 3 (0.25), which is a genuine feature.
      call amr%Estimate(sol,ivar=SELF_AMR_ALLVARS)
      call ReportFlags(amr,nel,nodeType,"all variables, unit weights")
      if(amr%flag(3) /= SELF_AMR_COARSEN) then
        print*,"FAIL: expected the background variable to swamp the unweighted gate"
        r = 1
      endif
      ! Giving the background field zero weight excludes it from the gate (this is what an
      ! entropy-weighted gate does for the LinearEuler3D background fields) and element 3 is
      ! judged on its own energy again.
      call amr%SetEnergyWeights([1.0_prec,0.0_prec])
      call amr%Estimate(sol,ivar=SELF_AMR_ALLVARS)
      call ReportFlags(amr,nel,nodeType,"all variables, background weighted out")
      if(amr%flag(3) /= SELF_AMR_REFINE) then
        print*,"FAIL: zero-weighting the background field did not restore element 3"
        r = 1
      endif
      if(amr%flag(2) /= SELF_AMR_COARSEN .or. amr%flag(1) /= SELF_AMR_REFINE) then
        print*,"FAIL: weighted gate disagrees with the single-variable gate"
        r = 1
      endif

      ! ---- 7. Host/device sync entry points round-trip (CPU no-ops) ----
      call amr%UpdateDevice()
      call amr%Estimate(sol,ivar=1)
      call amr%UpdateHost()
      if(amr%flag(1) /= SELF_AMR_REFINE) then
        print*,"FAIL: flags did not survive an UpdateDevice/UpdateHost round trip"
        r = 1
      endif

      call amr%Free()
      call sol%Free()
      call interp%Free()

    enddo

    if(r == 0) print*,"REFINEMENT INDICATOR 3D RELATIVE FLOOR CHECKS PASSED"

  endfunction refinement_indicator_3d_relative_floor

  subroutine ReportFlags(amr,nel,nodeType,label)
    !! Print the per-element indicator, gate energy and flag so a failure is diagnosable from the
    !! ctest log without a rebuild.
    use SELF_Constants
    use SELF_RefinementIndicator_3D
    implicit none
    type(RefinementIndicator3D),intent(in) :: amr
    integer,intent(in) :: nel
    integer,intent(in) :: nodeType
    character(*),intent(in) :: label
    ! Local
    integer :: k

    print*,"node type ",nodeType," : ",label
    do k = 1,nel
      print*,"  element, sigma, gate, flag :",k,amr%indicator(k),amr%gate(k),amr%flag(k)
    enddo

  endsubroutine ReportFlags

  pure function nlegendre(p,x) result(Lp)
    !! L2-normalized Legendre polynomial on [-1,1], matching the basis used by the indicator.
    use SELF_Constants
    implicit none
    integer,intent(in) :: p
    real(prec),intent(in) :: x
    real(prec) :: Lp
    integer :: k
    real(prec) :: lkm1,lk,lkp1

    if(p == 0) then
      Lp = 1.0_prec
    elseif(p == 1) then
      Lp = x
    else
      lkm1 = 1.0_prec
      lk = x
      do k = 1,p-1
        lkp1 = (real(2*k+1,prec)*x*lk-real(k,prec)*lkm1)/real(k+1,prec)
        lkm1 = lk
        lk = lkp1
      enddo
      Lp = lk
    endif
    Lp = Lp*sqrt((2.0_prec*real(p,prec)+1.0_prec)/2.0_prec)

  endfunction nlegendre

endprogram test
