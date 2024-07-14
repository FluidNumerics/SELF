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

  exit_code = scalargridinterp_1d_cpu_constant()
  stop exit_code

contains
  integer function scalargridinterp_1d_cpu_constant() result(r)
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Scalar_1D

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 1
    integer,parameter :: nelem = 100
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(Scalar1D) :: f
    type(Scalar1D) :: fTarget
    type(Lagrange),target :: interp
    type(Lagrange),target :: interpTarget
    real(prec) :: imat

    ! Create an interpolant
    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call interpTarget%Init(N=targetDegree, &
                           controlNodeType=UNIFORM, &
                           M=targetDegree, &
                           targetNodeType=UNIFORM)

    ! Initialize scalars
    call f%Init(interp,nvar,nelem)
    call fTarget%Init(interpTarget,nvar,nelem)

    ! ! Set the source scalar (on the control grid) to a non-zero constant
    f%interior(:,:,:) = 1.0_prec
    if(f%backend == 'cpu')then
      call f%GridInterp(fTarget%interior)
    else
      call f%UpdateDevice()
      call f%GridInterp(fTarget%interior_gpu)
      call fTarget%UpdateHost()
    endif
    ! ! Calculate diff from exact
    fTarget%interior = abs(fTarget%interior-1.0_prec)

    if(maxval(fTarget%interior) <= tolerance) then
      r = 0
    else
      print*,"max error : ",maxval(fTarget%interior)
      r = 1
    endif

    call f%free()
    call fTarget%free()
    call interp%free()
    call interpTarget%free()

  endfunction scalargridinterp_1d_cpu_constant
endprogram test
