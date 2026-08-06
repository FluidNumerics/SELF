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

module SELF_Tensor_2D

  use SELF_Constants
  use SELF_Tensor_2D_t
  use SELF_GPU
  use iso_c_binding
  use iso_fortran_env

  implicit none

  type,extends(Tensor2D_t),public :: Tensor2D
    character(3) :: backend = "gpu"
    type(c_ptr) :: interior_gpu
    type(c_ptr) :: boundary_gpu
    type(c_ptr) :: extBoundary_gpu
    !! Bytes currently allocated for each device buffer, so Resize can reuse them (Stage 6b/6c).
    integer(c_size_t) :: alloc_interior = 0
    integer(c_size_t) :: alloc_boundary = 0
    integer(c_size_t) :: alloc_extBoundary = 0

  contains

    procedure,public :: Init => Init_Tensor2D
    procedure,public :: Resize => Resize_Tensor2D
    procedure,public :: Free => Free_Tensor2D

    procedure,public :: UpdateHost => UpdateHost_Tensor2D
    procedure,public :: UpdateDevice => UpdateDevice_Tensor2D

  endtype Tensor2D

contains

  subroutine Init_Tensor2D(this,interp,nVar,nElem)
    implicit none
    class(Tensor2D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! local
    integer :: i

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    call this%MapArrays(interp%N+1,nVar,nElem)

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:4*nVar))

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%extBoundary = 0.0_prec

    ! Initialize equation parser
    ! This is done to prevent segmentation faults that arise
    ! when building with amdflang that are traced back to
    ! feqparse_functions.f90 : finalize routine
    ! When the equation parser is not initialized, the
    ! functions are not allocated, which I think are the
    ! source of the segfault - joe@fluidnumerics.com
    do i = 1,4*nvar
      this%eqn(i) = EquationParser('f=0',(/'x','y','z','t'/))
    enddo

    call EnsureDeviceBuffers_Tensor2D(this)

    call this%UpdateDevice()

  endsubroutine Init_Tensor2D

  subroutine EnsureDeviceBuffers_Tensor2D(this)
    !! Grow each device buffer to hold the current logical array, reusing the existing allocation
    !! when the bytes already fit (AMR Stage 6b/6c). A device pointer carries no shape, so only
    !! byte capacity has to be tracked.
    implicit none
    class(Tensor2D),intent(inout) :: this

    call EnsureDeviceBuffer(this%interior_gpu,this%alloc_interior,sizeof(this%interior))
    call EnsureDeviceBuffer(this%boundary_gpu,this%alloc_boundary,sizeof(this%boundary))
    call EnsureDeviceBuffer(this%extBoundary_gpu,this%alloc_extBoundary, &
                            sizeof(this%extBoundary))

  endsubroutine EnsureDeviceBuffers_Tensor2D

  subroutine Resize_Tensor2D(this,interp,nVar,nElem)
    !! Rebind to a new element count, reusing host pools and device buffers where they fit.
    !! Deliberately does NOT call UpdateDevice: Init uploads freshly zeroed arrays, which is pure
    !! waste in the adaptive loop because the arrays are rewritten immediately after.
    implicit none
    class(Tensor2D),intent(inout) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    call Resize_Tensor2D_t(this,interp,nVar,nElem)
    call EnsureDeviceBuffers_Tensor2D(this)

  endsubroutine Resize_Tensor2D

  subroutine Free_Tensor2D(this)
    implicit none
    class(Tensor2D),intent(inout) :: this

    this%interp => null()
    this%nVar = 0
    this%nElem = 0

    ! Host storage is owned by the pools in the parent type (Stage 6b/6c), so the parent Free
    ! releases it rather than deallocating these pointers.
    call Free_Tensor2D_t(this)

    call gpuCheck(hipFree(this%interior_gpu))
    call gpuCheck(hipFree(this%boundary_gpu))
    call gpuCheck(hipFree(this%extBoundary_gpu))
    this%interior_gpu = c_null_ptr
    this%boundary_gpu = c_null_ptr
    this%extBoundary_gpu = c_null_ptr
    this%alloc_interior = 0
    this%alloc_boundary = 0
    this%alloc_extBoundary = 0

  endsubroutine Free_Tensor2D

  subroutine UpdateHost_Tensor2D(this)
    implicit none
    class(Tensor2D),intent(inout) :: this

    call gpuCheck(hipMemcpy(c_loc(this%interior),this%interior_gpu,sizeof(this%interior),hipMemcpyDeviceToHost))
    call gpuCheck(hipMemcpy(c_loc(this%boundary),this%boundary_gpu,sizeof(this%boundary),hipMemcpyDeviceToHost))
    call gpuCheck(hipMemcpy(c_loc(this%extboundary),this%extboundary_gpu,sizeof(this%extboundary),hipMemcpyDeviceToHost))

  endsubroutine UpdateHost_Tensor2D

  subroutine UpdateDevice_Tensor2D(this)
    implicit none
    class(Tensor2D),intent(inout) :: this

    call gpuCheck(hipMemcpy(this%interior_gpu,c_loc(this%interior),sizeof(this%interior),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%boundary_gpu,c_loc(this%boundary),sizeof(this%boundary),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%extboundary_gpu,c_loc(this%extboundary),sizeof(this%extboundary),hipMemcpyHostToDevice))

  endsubroutine UpdateDevice_Tensor2D

endmodule SELF_Tensor_2D
