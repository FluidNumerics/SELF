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

  exit_code = refinement_indicator_2d_spectraldecay()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function refinement_indicator_2d_spectraldecay() result(r)
    !! Validates the 2-D Legendre modal-decay refinement indicator against analytic fields whose
    !! modal spectra are known exactly, on both Legendre-Gauss and Legendre-Gauss-Lobatto nodes:
    !!
    !!   1. A constant field has energy only in the (0,0) mode, so the top-shell energy fraction
    !!      is zero: the indicator sits at its machine-zero floor and the element is flagged
    !!      COARSEN.
    !!   2. A pure highest tensor mode Ltilde_N(xi)*Ltilde_N(eta) has all of its energy in the top
    !!      shell, so S_e = 1 and the indicator is 0 (to roundoff): the element is flagged REFINE.
    !!      This simultaneously checks that the nodal->modal transform recovers a single mode
    !!      exactly.
    !!   3. A low-degree polynomial (degree <= N-2 per direction) has no energy in the top two
    !!      shells, so it too floors out and is flagged COARSEN.
    !!   4. A steep hyperbolic-tangent front is under-resolved at this degree and is flagged
    !!      REFINE.
    !!
    !! The indicator is built from the exact inverse Legendre Vandermonde, so cases 1-3 hold to
    !! roundoff independent of the control-node type.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Scalar_2D
    use SELF_RefinementIndicator_2D

    implicit none

    integer,parameter :: N = 7
    integer,parameter :: nel = 4
    integer,parameter :: nvar = 1
    real(prec),parameter :: refineThreshold = -3.0_prec
    real(prec),parameter :: coarsenThreshold = -8.0_prec
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 1.0e-9_prec
#else
    real(prec),parameter :: tolerance = 1.0e-3_prec
#endif
    type(Lagrange),target :: interp
    type(Scalar2D) :: sol
    type(RefinementIndicator2D) :: amr
    integer :: nodeType,i,j,iel
    real(prec) :: x,y

    r = 0

    do nodeType = 1,2 ! GAUSS = 1, GAUSS_LOBATTO = 2

      call interp%Init(N=N,controlNodeType=nodeType,M=N,targetNodeType=UNIFORM)
      call sol%Init(interp,nvar,nel)

      do iel = 1,nel
        do j = 1,N+1
          do i = 1,N+1
            x = interp%controlPoints(i)
            y = interp%controlPoints(j)
            select case(iel)
            case(1) ! constant -> resolved -> coarsen
              sol%interior(i,j,iel,1) = 2.5_prec
            case(2) ! pure highest tensor mode -> S_e = 1 -> sigma = 0 -> refine
              sol%interior(i,j,iel,1) = nlegendre(N,x)*nlegendre(N,y)
            case(3) ! low-degree polynomial (deg <= N-2) -> resolved -> coarsen
              sol%interior(i,j,iel,1) = 1.0_prec+0.5_prec*x+0.25_prec*x*y-0.3_prec*y*y
            case(4) ! steep front -> under-resolved -> refine
              sol%interior(i,j,iel,1) = tanh(20.0_prec*(x-0.1_prec))
            endselect
          enddo
        enddo
      enddo
      call sol%UpdateDevice() ! no-op on CPU; pushes interior to device on GPU

      call amr%Init(interp,nel,refineThreshold,coarsenThreshold)
      call amr%Estimate(sol,ivar=1)

      print*,"node type :",nodeType
      do iel = 1,nel
        print*,"  element, sigma, flag :",iel,amr%indicator(iel),amr%flag(iel)
      enddo

      ! Case 1: constant field -> coarsen, indicator below the coarsen threshold.
      if(amr%flag(1) /= SELF_AMR_COARSEN) then
        print*,"FAIL: constant field not flagged COARSEN"
        r = 1
      endif

      ! Case 2: pure top tensor mode -> refine, and sigma == log10(1) == 0 to roundoff.
      if(amr%flag(2) /= SELF_AMR_REFINE) then
        print*,"FAIL: pure top-mode field not flagged REFINE"
        r = 1
      endif
      if(abs(amr%indicator(2)) > tolerance) then
        print*,"FAIL: pure top-mode indicator should be 0, got",amr%indicator(2)
        r = 1
      endif

      ! Case 3: low-degree polynomial -> coarsen.
      if(amr%flag(3) /= SELF_AMR_COARSEN) then
        print*,"FAIL: low-degree polynomial not flagged COARSEN"
        r = 1
      endif

      ! Case 4: steep front -> refine.
      if(amr%flag(4) /= SELF_AMR_REFINE) then
        print*,"FAIL: steep-front field not flagged REFINE"
        r = 1
      endif

      ! Flag tally sanity: 2 refine, 2 coarsen, 0 keep.
      if(amr%CountFlagged(SELF_AMR_REFINE) /= 2 .or. &
         amr%CountFlagged(SELF_AMR_COARSEN) /= 2 .or. &
         amr%CountFlagged(SELF_AMR_KEEP) /= 0) then
        print*,"FAIL: unexpected flag tally"
        r = 1
      endif

      call amr%Free()
      call sol%Free()
      call interp%Free()

    enddo

  endfunction refinement_indicator_2d_spectraldecay

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
