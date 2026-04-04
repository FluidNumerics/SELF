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

  exit_code = twopointvectordivergence_3d_constant()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function twopointvectordivergence_3d_constant() result(r)
    !! Verifies that the 3-D split-form divergence of a constant two-point
    !! vector field is identically zero.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Scalar_3D
    use SELF_TwoPointVector_3D

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 3
    integer,parameter :: nelem = 64
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(TwoPointVector3D) :: f
    type(Scalar3D) :: df
    type(Lagrange),target :: interp

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS_LOBATTO, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call f%Init(interp,nvar,nelem)
    call df%Init(interp,nvar,nelem)

    f%interior = 1.0_prec
    call f%UpdateDevice()

#ifdef ENABLE_GPU
    call f%Divergence(df%interior_gpu)
#else
    call f%Divergence(df%interior)
#endif
    call df%UpdateHost()

    ! Add surface term for each of three coordinate directions
    block
      integer :: i,j,k,iEl,iVar
      do concurrent(i=1:controlDegree+1,j=1:controlDegree+1, &
                    k=1:controlDegree+1,iEl=1:nelem,iVar=1:nvar)

        df%interior(i,j,k,iEl,iVar) = df%interior(i,j,k,iEl,iVar)+ &
                                      (interp%bMatrix(i,2)*(+1.0_prec)+ &
                                       interp%bMatrix(i,1)*(-1.0_prec))/interp%qWeights(i)+ &
                                      (interp%bMatrix(j,2)*(+1.0_prec)+ &
                                       interp%bMatrix(j,1)*(-1.0_prec))/interp%qWeights(j)+ &
                                      (interp%bMatrix(k,2)*(+1.0_prec)+ &
                                       interp%bMatrix(k,1)*(-1.0_prec))/interp%qWeights(k)

      enddo
    endblock

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

  endfunction twopointvectordivergence_3d_constant

endprogram test
