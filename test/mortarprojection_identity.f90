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

  exit_code = mortarprojection_identity()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function mortarprojection_identity() result(r)
    !! Verifies the discrete identities of the 2:1 mortar operators on Lagrange_t :
    !!
    !!  1. Forward-backward consistency : sum_k P_k R_k = I, i.e. a polynomial trace
    !!     restricted to the two sub-edges and projected back is recovered exactly.
    !!  2. Discrete conservation : sum_i w_i (P_k g)_i = (1/2) sum_j w_j g_j for any
    !!     sub-edge trace g.
    !!
    !! Both identities are checked for Legendre-Gauss and Legendre-Gauss-Lobatto
    !! control points.
    use SELF_Constants
    use SELF_Lagrange

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 15
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(Lagrange),target :: interp
    integer :: nodeType,i,ii,jj,k
    real(prec) :: idmat,idErr,consErr,bigSum,smallSum
    real(prec) :: g(1:controlDegree+1)

    r = 0

    do nodeType = 1,2 ! GAUSS = 1, GAUSS_LOBATTO = 2

      call interp%Init(N=controlDegree, &
                       controlNodeType=nodeType, &
                       M=targetDegree, &
                       targetNodeType=UNIFORM)

      ! Identity check : (sum_k P_k R_k)(i,ii) = delta(i,ii)
      idErr = 0.0_prec
      do i = 1,controlDegree+1
        do ii = 1,controlDegree+1
          idmat = 0.0_prec
          do k = 1,2
            do jj = 1,controlDegree+1
              idmat = idmat+interp%mortarP(jj,i,k)*interp%mortarR(ii,jj,k)
            enddo
          enddo
          if(ii == i) idmat = idmat-1.0_prec
          idErr = max(idErr,abs(idmat))
        enddo
      enddo
      print*,"node type, max |sum_k P_k R_k - I| :",nodeType,idErr

      ! Conservation check with an arbitrary (non-polynomial-symmetric) trace
      do ii = 1,controlDegree+1
        g(ii) = sin(3.0_prec*real(ii,prec))+0.3_prec*real(ii,prec)**2
      enddo

      consErr = 0.0_prec
      do k = 1,2
        bigSum = 0.0_prec
        do i = 1,controlDegree+1
          do ii = 1,controlDegree+1
            bigSum = bigSum+interp%qWeights(i)*interp%mortarP(ii,i,k)*g(ii)
          enddo
        enddo
        smallSum = 0.0_prec
        do ii = 1,controlDegree+1
          smallSum = smallSum+0.5_prec*interp%qWeights(ii)*g(ii)
        enddo
        consErr = max(consErr,abs(bigSum-smallSum))
      enddo
      print*,"node type, max conservation defect :",nodeType,consErr

      if(idErr > tolerance .or. consErr > tolerance) then
        r = 1
      endif

      call interp%Free()

    enddo

  endfunction mortarprojection_identity
endprogram test
