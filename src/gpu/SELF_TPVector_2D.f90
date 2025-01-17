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

module SELF_TPVector_2D

  use SELF_Constants
  use SELF_TPVector_2D_t
  use SELF_GPU
  use SELF_GPUBLAS
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(TPVector2D_t),public :: TPVector2D
    character(3) :: backend = "gpu"
    type(c_ptr) :: interior_gpu
    type(c_ptr) :: boundary_gpu
    type(c_ptr) :: riemann_gpu

  contains

    procedure,public :: Init => Init_TPVector2D
    procedure,public :: Free => Free_TPVector2D

    procedure,public :: UpdateHost => UpdateHost_TPVector2D
    procedure,public :: UpdateDevice => UpdateDevice_TPVector2D

  endtype TPVector2D

contains

  subroutine Init_TPVector2D(this,interp,nVar,nElem)
    implicit none
    class(TPVector2D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! local
    integer :: i
    integer(c_size_t) :: workSize

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar,1:2), &
             this%boundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
             this%riemann(1:interp%N+1,1:4,1:nelem,1:nvar))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:2*nVar))

    ! Initialize equation parser
    ! This is done to prevent segmentation faults that arise
    ! when building with amdflang that are traced back to
    ! feqparse_functions.f90 : finalize routine
    ! When the equation parser is not initialized, the
    ! functions are not allocated, which I think are the
    ! source of the segfault - joe@fluidnumerics.com
    do i = 1,2*nvar
      this%eqn(i) = EquationParser('f=0',(/'x','y','z','t'/))
    enddo

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%riemann = 0.0_prec

    call gpuCheck(hipMalloc(this%interior_gpu,sizeof(this%interior)))
    call gpuCheck(hipMalloc(this%boundary_gpu,sizeof(this%boundary)))
    call gpuCheck(hipMalloc(this%riemann_gpu,sizeof(this%riemann)))

    call this%UpdateDevice()

  endsubroutine Init_TPVector2D

  subroutine Free_TPVector2D(this)
    implicit none
    class(TPVector2D),intent(inout) :: this

    this%interp => null()
    this%nVar = 0
    this%nElem = 0

    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%riemann)

    deallocate(this%meta)
    deallocate(this%eqn)

    call gpuCheck(hipFree(this%interior_gpu))
    call gpuCheck(hipFree(this%boundary_gpu))
    call gpuCheck(hipFree(this%riemann_gpu))

  endsubroutine Free_TPVector2D

  subroutine UpdateHost_TPVector2D(this)
    implicit none
    class(TPVector2D),intent(inout) :: this

    call gpuCheck(hipMemcpy(c_loc(this%interior),this%interior_gpu,sizeof(this%interior),hipMemcpyDeviceToHost))
    call gpuCheck(hipMemcpy(c_loc(this%boundary),this%boundary_gpu,sizeof(this%boundary),hipMemcpyDeviceToHost))
    call gpuCheck(hipMemcpy(c_loc(this%riemann),this%riemann_gpu,sizeof(this%riemann),hipMemcpyDeviceToHost))

  endsubroutine UpdateHost_TPVector2D

  subroutine UpdateDevice_TPVector2D(this)
    implicit none
    class(TPVector2D),intent(inout) :: this

    call gpuCheck(hipMemcpy(this%interior_gpu,c_loc(this%interior),sizeof(this%interior),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%boundary_gpu,c_loc(this%boundary),sizeof(this%boundary),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%riemann_gpu,c_loc(this%riemann),sizeof(this%riemann),hipMemcpyHostToDevice))

  endsubroutine UpdateDevice_TPVector2D

endmodule SELF_TPVector_2D
