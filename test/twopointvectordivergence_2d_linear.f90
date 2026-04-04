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
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
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

  exit_code = twopointvectordivergence_2d_linear()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function twopointvectordivergence_2d_linear() result(r)
    !! Verifies that the split-form divergence is exact for the field
    !!   V = xi1 * e_xi1 + xi2 * e_xi2
    !! on the reference element, for which the exact divergence is 2.
    !!
    !! The two-point arithmetic-mean flux for each computational direction is:
    !!   F^1(nn,i,j) = (xi1_i + xi1_nn) / 2    (xi^1 sum pairs)
    !!   F^2(nn,i,j) = (xi2_j + xi2_nn) / 2    (xi^2 sum pairs)
    !!
    !! Because the split-form sum with an arithmetic-mean two-point flux
    !! reproduces the standard spectral-element derivative exactly, this test
    !! verifies that the Divergence kernel correctly applies the derivative
    !! matrix with the factor-of-two scaling.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Scalar_2D
    use SELF_TwoPointVector_2D

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 3
    integer,parameter :: nelem = 100
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(TwoPointVector2D) :: f
    type(Scalar2D) :: df
    type(Lagrange),target :: interp
    integer :: i,j,nn,iel,ivar

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call f%Init(interp,nvar,nelem)
    call df%Init(interp,nvar,nelem)

    ! Set up two-point arithmetic-mean fluxes for V = (xi1, xi2).
    ! interior(nn,i,j,iel,ivar,1) = xi^1 two-point flux for pair (i,j)-(nn,j)
    ! interior(nn,i,j,iel,ivar,2) = xi^2 two-point flux for pair (i,j)-(i,nn)
    do ivar = 1,nvar
      do iel = 1,nelem
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            do nn = 1,controlDegree+1
              f%interior(nn,i,j,iel,ivar,1) = &
                0.5_prec*(interp%controlPoints(i)+interp%controlPoints(nn))
              f%interior(nn,i,j,iel,ivar,2) = &
                0.5_prec*(interp%controlPoints(j)+interp%controlPoints(nn))
            enddo
          enddo
        enddo
      enddo
    enddo
    call f%UpdateDevice()

#ifdef ENABLE_GPU
    call f%Divergence(df%interior_gpu)
#else
    call f%Divergence(df%interior)
#endif
    call df%UpdateHost()

    print*,"absmax (nabla.F) :",maxval(df%interior)
    df%interior = abs(df%interior-2.0_prec)
    print*,"absmax error     :",maxval(df%interior)

    if(maxval(df%interior) <= tolerance) then
      r = 0
    else
      r = 1
    endif

    call f%free()
    call df%free()
    call interp%free()

  endfunction twopointvectordivergence_2d_linear

endprogram test
