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

  exit_code = refinement_indicator_3d_energy_hysteresis()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function refinement_indicator_3d_energy_hysteresis() result(r)
    !! Validates the hysteresis band on the amplitude gate's energy axis.
    !!
    !! A single hard gate reintroduces, on the energy axis, exactly the thrashing that the two
    !! sigma thresholds exist to prevent: an element whose amplitude drifts across that one value
    !! flips COARSEN <-> REFINE on successive epochs, so the mesh churns without the solution
    !! changing meaningfully. Opening a band between a quiescent floor and a significant floor
    !! makes the middle zone report SELF_AMR_KEEP - too weak to justify spending levels on, too
    !! strong to declare resolved - which is stable under a small change in amplitude.
    !!
    !! Every element carries the same pure highest tensor mode (S_e = 1 exactly, so the spectrum
    !! always says REFINE) scaled to a different amplitude, so only the energy axis can change a
    !! flag. Element 1 has amplitude 1 and sets the energy scale; the others sample the three
    !! zones of a band spanning [1e-12, 1e-8] of that scale:
    !!
    !!   element 2 : E = 1e-14  below the quiescent floor      -> COARSEN
    !!   element 3 : E = 1e-10  inside the band                -> KEEP
    !!   element 4 : E = 1e-6   above the significant floor    -> REFINE
    !!
    !! The thrashing check then walks element 3's amplitude across the OLD single-cut boundary and
    !! requires the flag not to move.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Scalar_3D
    use SELF_RefinementIndicator_3D

    implicit none

    integer,parameter :: N = 7
    integer,parameter :: nel = 4
    integer,parameter :: nvar = 1
    real(prec),parameter :: refineThreshold = -3.0_prec
    real(prec),parameter :: coarsenThreshold = -8.0_prec
    real(prec),parameter :: quiescentFloor = 1.0e-12_prec
    real(prec),parameter :: significantFloor = 1.0e-8_prec
    ! Amplitudes: energies 1, 1e-14, 1e-10, 1e-6 (the square of each).
    real(prec),parameter :: amp(1:nel) = &
                            [1.0_prec,1.0e-7_prec,1.0e-5_prec,1.0e-3_prec]
    type(Lagrange),target :: interp
    type(Scalar3D) :: sol
    type(RefinementIndicator3D) :: amr
    integer :: nodeType,k,flagBefore
    real(prec) :: drift

    r = 0

    do nodeType = 1,2 ! GAUSS = 1, GAUSS_LOBATTO = 2

      call interp%Init(N=N,controlNodeType=nodeType,M=N,targetNodeType=UNIFORM)
      call sol%Init(interp,nvar,nel)
      call amr%Init(interp,nel,refineThreshold,coarsenThreshold)

      call SetAmplitudes(sol,interp,amp)
      call amr%SetRelativeEnergyFloor(quiescentFloor,significantEnergyFloor=significantFloor)
      call amr%Estimate(sol,ivar=1)

      print*,"node type ",nodeType," : band [1e-12,1e-8]"
      do k = 1,nel
        print*,"  element, sigma, gate, flag :",k,amr%indicator(k),amr%gate(k),amr%flag(k)
      enddo

      ! ---- 1. The three zones ----
      if(amr%flag(1) /= SELF_AMR_REFINE .or. amr%flag(4) /= SELF_AMR_REFINE) then
        print*,"FAIL: elements above the significant floor should be flagged REFINE"
        r = 1
      endif
      if(amr%flag(2) /= SELF_AMR_COARSEN) then
        print*,"FAIL: element below the quiescent floor should be flagged COARSEN"
        r = 1
      endif
      if(amr%flag(3) /= SELF_AMR_KEEP) then
        print*,"FAIL: element inside the energy band should be flagged KEEP"
        r = 1
      endif

      ! Inside the band the reported sigma must still be the TRUE spectral value (S_e = 1 here),
      ! not the quiescent floor value: the band holds the flag, it does not falsify the diagnostic.
      if(abs(amr%indicator(3)) > 1.0e-6_prec) then
        print*,"FAIL: sigma inside the band should still report the true spectrum, got", &
          amr%indicator(3)
        r = 1
      endif

      ! ---- 2. No thrashing across the old single-cut boundary ----
      ! With one cut at 1e-10 the element below would COARSEN and the element above would REFINE.
      ! Inside a band that brackets it, a factor-of-ten drift either way must not move the flag.
      flagBefore = amr%flag(3)
      do k = -1,1
        drift = 10.0_prec**real(k,prec)
        call SetAmplitudes(sol,interp,[amp(1),amp(2),amp(3)*sqrt(drift),amp(4)])
        call amr%Estimate(sol,ivar=1)
        if(amr%flag(3) /= flagBefore) then
          print*,"FAIL: flag thrashed as the in-band energy drifted by",drift, &
            " : ",flagBefore," -> ",amr%flag(3)
          r = 1
        endif
      enddo

      ! ---- 3. A degenerate band reproduces the single hard cut ----
      ! Band collapsed onto a single value of 1e-9, which brackets element 3's energy (1e-10) and
      ! its tenfold-amplitude version (1e-8) with an order of magnitude either side. The flag then
      ! flips COARSEN -> REFINE across that one cut: this is the thrashing the band exists to damp,
      ! asserted here rather than assumed, and it is the same drift that left the flag untouched in
      ! check 2 above.
      call SetAmplitudes(sol,interp,amp)
      call amr%SetRelativeEnergyFloor(1.0e-9_prec)
      call amr%Estimate(sol,ivar=1)
      if(amr%flag(3) /= SELF_AMR_COARSEN) then
        print*,"FAIL: at a degenerate band an element below the cut should be COARSEN"
        r = 1
      endif
      call SetAmplitudes(sol,interp,[amp(1),amp(2),amp(3)*10.0_prec,amp(4)])
      call amr%Estimate(sol,ivar=1)
      if(amr%flag(3) /= SELF_AMR_REFINE) then
        print*,"FAIL: at a degenerate band a tenfold rise should flip the element to REFINE"
        r = 1
      endif

      call amr%Free()
      call sol%Free()
      call interp%Free()

    enddo

    if(r == 0) print*,"REFINEMENT INDICATOR 3D ENERGY HYSTERESIS CHECKS PASSED"

  endfunction refinement_indicator_3d_energy_hysteresis

  subroutine SetAmplitudes(sol,interp,amp)
    !! Fill every element with the pure highest tensor Legendre mode at the given per-element
    !! amplitude, so the modal energy of element e is exactly amp(e)**2 and the smoothness ratio
    !! is exactly 1 on every element.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Scalar_3D
    implicit none
    type(Scalar3D),intent(inout) :: sol
    type(Lagrange),intent(in) :: interp
    real(prec),intent(in) :: amp(:)
    ! Local
    integer :: i,j,k,iel

    do iel = 1,size(amp)
      do k = 1,interp%N+1
        do j = 1,interp%N+1
          do i = 1,interp%N+1
            sol%interior(i,j,k,iel,1) = amp(iel)* &
                                        nlegendre(interp%N,interp%controlPoints(i))* &
                                        nlegendre(interp%N,interp%controlPoints(j))* &
                                        nlegendre(interp%N,interp%controlPoints(k))
          enddo
        enddo
      enddo
    enddo
    call sol%UpdateDevice() ! no-op on CPU; pushes interior to device on GPU

  endsubroutine SetAmplitudes

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
