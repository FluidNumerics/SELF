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

module SELF_TwoPointVector_3D

  use SELF_Constants
  use SELF_TwoPointVector_3D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(TwoPointVector3D_t),public :: TwoPointVector3D
    character(3) :: backend = "gpu"
    type(c_ptr) :: interior_gpu

  contains

    procedure,public :: Init => Init_TwoPointVector3D
    procedure,public :: Free => Free_TwoPointVector3D

    procedure,public :: UpdateHost => UpdateHost_TwoPointVector3D
    procedure,public :: UpdateDevice => UpdateDevice_TwoPointVector3D

    generic,public :: Divergence => Divergence_TwoPointVector3D
    procedure,private :: Divergence_TwoPointVector3D

  endtype TwoPointVector3D

contains

  subroutine Init_TwoPointVector3D(this,interp,nVar,nElem)
    implicit none
    class(TwoPointVector3D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer :: i

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:interp%N+1, &
                           1:nElem,1:nVar,1:3))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:3*nVar))

    ! Initialize equation parser to prevent segmentation faults with amdflang
    do i = 1,3*nVar
      this%eqn(i) = EquationParser('f=0',(/'x','y','z','t'/))
    enddo

    this%interior = 0.0_prec

    call gpuCheck(hipMalloc(this%interior_gpu,sizeof(this%interior)))

    call this%UpdateDevice()

  endsubroutine Init_TwoPointVector3D

  subroutine Free_TwoPointVector3D(this)
    implicit none
    class(TwoPointVector3D),intent(inout) :: this

    this%interp => null()
    this%nVar = 0
    this%nElem = 0

    deallocate(this%interior)
    deallocate(this%meta)
    deallocate(this%eqn)

    call gpuCheck(hipFree(this%interior_gpu))

  endsubroutine Free_TwoPointVector3D

  subroutine UpdateHost_TwoPointVector3D(this)
    implicit none
    class(TwoPointVector3D),intent(inout) :: this

    call gpuCheck(hipMemcpy(c_loc(this%interior),this%interior_gpu, &
                            sizeof(this%interior),hipMemcpyDeviceToHost))

  endsubroutine UpdateHost_TwoPointVector3D

  subroutine UpdateDevice_TwoPointVector3D(this)
    implicit none
    class(TwoPointVector3D),intent(inout) :: this

    call gpuCheck(hipMemcpy(this%interior_gpu,c_loc(this%interior), &
                            sizeof(this%interior),hipMemcpyHostToDevice))

  endsubroutine UpdateDevice_TwoPointVector3D

  subroutine Divergence_TwoPointVector3D(this,df)
    !! GPU implementation of the reference-element split-form divergence.
    !! df must be a device pointer to a scalar field of size
    !! (N+1)^3 * nElem * nVar.
    implicit none
    class(TwoPointVector3D),intent(in) :: this
    type(c_ptr),intent(inout) :: df

    call TwoPointVectorDivergence_3D_gpu(this%interior_gpu,df, &
                                         this%interp%dSplitMatrix_gpu, &
                                         this%interp%N,this%nVar,this%nElem)

  endsubroutine Divergence_TwoPointVector3D

endmodule SELF_TwoPointVector_3D
