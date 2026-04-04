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

  exit_code = twopointvectordivergence_2d_rotational()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function twopointvectordivergence_2d_rotational() result(r)
    !! Verifies that the split-form divergence is identically zero for the
    !! purely rotational field V = -xi2 * e_xi1 + xi1 * e_xi2.
    !!
    !! The exact divergence is d(-xi2)/d(xi1) + d(xi1)/d(xi2) = 0.
    !!
    !! Two-point arithmetic-mean fluxes (linear components, no projection error):
    !!   F^1(nn,i,j) = (-xi2_j + (-xi2_j))/2 = -xi2_j    (independent of nn)
    !!   F^2(nn,i,j) = (xi1_i + xi1_i)/2     =  xi1_i    (independent of nn)
    !!
    !! Since neither flux depends on the two-point index nn, the split-form
    !! sum reduces to the standard D-matrix column sum, which is zero.
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

    ! V = (-xi2, xi1): two-point fluxes are independent of nn.
    ! F^1(nn,i,j) = -xi2_j  (xi1-component doesn't depend on xi1)
    ! F^2(nn,i,j) =  xi1_i  (xi2-component doesn't depend on xi2)
    do ivar = 1,nvar
      do iel = 1,nelem
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            do nn = 1,controlDegree+1
              f%interior(nn,i,j,iel,ivar,1) = -interp%controlPoints(j)
              f%interior(nn,i,j,iel,ivar,2) = interp%controlPoints(i)
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

    df%interior = abs(df%interior-0.0_prec)

    if(maxval(df%interior) <= tolerance) then
      r = 0
    else
      print*,"absmax error :",maxval(df%interior)
      r = 1
    endif

    call f%free()
    call df%free()
    call interp%free()

  endfunction twopointvectordivergence_2d_rotational

endprogram test
