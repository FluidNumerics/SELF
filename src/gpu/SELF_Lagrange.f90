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

module SELF_Lagrange

  use iso_fortran_env
  use iso_c_binding

  use SELF_Constants
  use SELF_Lagrange_t
  use SELF_GPU

  implicit none

  type,extends(Lagrange_t),public :: Lagrange
    character(3) :: backend = 'gpu'
    type(c_ptr) :: qWeights_gpu
    type(c_ptr) :: iMatrix_gpu
    type(c_ptr) :: dMatrix_gpu
    type(c_ptr) :: dgMatrix_gpu
    type(c_ptr) :: bMatrix_gpu

  contains

    procedure,public :: Init => Init_Lagrange
    procedure,public :: Free => Free_Lagrange

  endtype Lagrange

contains

  subroutine Init_Lagrange(this,N,controlNodeType,M,targetNodeType)
    !! Initialize an instance of the Lagrange class
    !! On output, all of the attributes for the Lagrange class are allocated and values are initialized according to the number of
    !! control points, number of target points, and the types for the control and target nodes.
    !! If a GPU is available, device pointers for the Lagrange attributes are allocated and initialized.
    implicit none
    class(Lagrange),intent(out) :: this
    !! Lagrange class instance
    integer,intent(in)          :: N
    !! The number of control points for interpolant
    integer,intent(in)          :: M
    !! The number of target points for the interpolant
    integer,intent(in)          :: controlNodeType
    !! The integer code specifying the type of control points. Parameters are defined in SELF_Constants.f90. One of GAUSS(=1),
    !! GAUSS_LOBATTO(=2), or UNIFORM(=3)
    integer,intent(in)          :: targetNodeType
    !! The integer code specifying the type of target points. Parameters are defined in SELF_Constants.f90. One of GAUSS(=1),
    !! GAUSS_LOBATTO(=2), or UNIFORM(=3)
    ! -------!
    ! Local
    real(prec) :: q(0:M)

    if(.not. GPUAvailable()) then
      print*,__FILE__,':',__LINE__,' : Error : Attempt to use GPU extension, but GPU is not available.'
      stop 1
    endif

    this%N = N
    this%M = M
    this%controlNodeType = controlNodeType
    this%targetNodeType = targetNodeType
    allocate(this%controlPoints(1:N+1), &
             this%targetPoints(1:M+1), &
             this%bWeights(1:N+1), &
             this%qWeights(1:N+1), &
             this%iMatrix(1:N+1,1:M+1), &
             this%dMatrix(1:N+1,1:N+1), &
             this%dgMatrix(1:N+1,1:N+1), &
             this%bMatrix(1:N+1,1:2))

    if(controlNodeType == GAUSS .or. controlNodeType == GAUSS_LOBATTO) then

      call LegendreQuadrature(N, &
                              this%controlPoints, &
                              this%qWeights, &
                              controlNodeType)

    elseif(controlNodeType == CHEBYSHEV_GAUSS .or. controlNodeType == CHEBYSHEV_GAUSS_LOBATTO) then

      call ChebyshevQuadrature(N, &
                               this%controlPoints, &
                               this%qWeights, &
                               controlNodeType)

    elseif(controlNodeType == UNIFORM) then

      this%controlPoints = UniformPoints(-1.0_prec,1.0_prec,0,N)
      this%qWeights = 2.0_prec/real(N,prec)

    endif

    ! Target Points
    if(targetNodeType == GAUSS .or. targetNodeType == GAUSS_LOBATTO) then

      call LegendreQuadrature(M, &
                              this%targetPoints, &
                              q, &
                              targetNodeType)

    elseif(targetNodeType == UNIFORM) then

      this%targetPoints = UniformPoints(-1.0_prec,1.0_prec,0,M)

    endif

    call this%CalculateBarycentricWeights()
    call this%CalculateInterpolationMatrix()
    call this%CalculateDerivativeMatrix()
    this%bMatrix(1:N+1,1) = this%CalculateLagrangePolynomials(-1.0_prec)
    this%bMatrix(1:N+1,2) = this%CalculateLagrangePolynomials(1.0_prec)

    call gpuCheck(hipMalloc(this%iMatrix_gpu,sizeof(this%iMatrix)))
    call gpuCheck(hipMalloc(this%dMatrix_gpu,sizeof(this%dMatrix)))
    call gpuCheck(hipMalloc(this%dgMatrix_gpu,sizeof(this%dgMatrix)))
    call gpuCheck(hipMalloc(this%bMatrix_gpu,sizeof(this%bMatrix)))
    call gpuCheck(hipMalloc(this%qWeights_gpu,sizeof(this%qWeights)))

    call gpuCheck(hipMemcpy(this%iMatrix_gpu, &
                            c_loc(this%iMatrix), &
                            sizeof(this%iMatrix), &
                            hipMemcpyHostToDevice))

    call gpuCheck(hipMemcpy(this%dMatrix_gpu, &
                            c_loc(this%dMatrix), &
                            sizeof(this%dMatrix), &
                            hipMemcpyHostToDevice))

    call gpuCheck(hipMemcpy(this%dgMatrix_gpu, &
                            c_loc(this%dgMatrix), &
                            sizeof(this%dgMatrix), &
                            hipMemcpyHostToDevice))

    call gpuCheck(hipMemcpy(this%bMatrix_gpu, &
                            c_loc(this%bMatrix), &
                            sizeof(this%bMatrix), &
                            hipMemcpyHostToDevice))

    call gpuCheck(hipMemcpy(this%qWeights_gpu, &
                            c_loc(this%qWeights), &
                            sizeof(this%qWeights), &
                            hipMemcpyHostToDevice))

  endsubroutine Init_Lagrange

  subroutine Free_Lagrange(this)
    !! Frees all memory (host and device) associated with an instance of the Lagrange class
    implicit none
    class(Lagrange),intent(inout) :: this
    !! Lagrange class instance

    deallocate(this%controlPoints)
    deallocate(this%targetPoints)
    deallocate(this%bWeights)
    deallocate(this%qWeights)
    deallocate(this%iMatrix)
    deallocate(this%dMatrix)
    deallocate(this%dgMatrix)
    deallocate(this%bMatrix)

    call gpuCheck(hipFree(this%iMatrix_gpu))
    call gpuCheck(hipFree(this%dMatrix_gpu))
    call gpuCheck(hipFree(this%dgMatrix_gpu))
    call gpuCheck(hipFree(this%bMatrix_gpu))
    call gpuCheck(hipFree(this%qWeights_gpu))

  endsubroutine Free_Lagrange

endmodule SELF_Lagrange
