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

module SELF_RefinementIndicator_2D

  use SELF_Constants
  use SELF_RefinementIndicator_2D_t
  use SELF_Scalar_2D
  use SELF_GPU
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(RefinementIndicator2D_t),public :: RefinementIndicator2D
    character(3) :: backend = "gpu"
    type(c_ptr) :: Pmodal_gpu = c_null_ptr
    type(c_ptr) :: indicator_gpu = c_null_ptr
    type(c_ptr) :: flag_gpu = c_null_ptr

  contains

    procedure,public :: Init => Init_RefinementIndicator2D
    procedure,public :: Free => Free_RefinementIndicator2D
    procedure,public :: UpdateHost => UpdateHost_RefinementIndicator2D
    procedure,public :: UpdateDevice => UpdateDevice_RefinementIndicator2D
    procedure,public :: Estimate => Estimate_RefinementIndicator2D

  endtype RefinementIndicator2D

contains

  subroutine Init_RefinementIndicator2D(this,interp,nElem,refineThreshold,coarsenThreshold)
    implicit none
    class(RefinementIndicator2D),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nElem
    real(prec),intent(in) :: refineThreshold
    real(prec),intent(in) :: coarsenThreshold

    ! Reuse the portable host initialization (allocation + modal transform build).
    call Init_RefinementIndicator2D_t(this,interp,nElem,refineThreshold,coarsenThreshold)

    ! The device kernel uses per-thread scratch sized to (N+1) <= 16 (AMR2D_MAXNP in
    ! SELF_Refinement.cpp). Fail early rather than overrun that bound on the device.
    if(interp%N+1 > 16) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : GPU refinement indicator supports degree N <= 15 '// &
        '(raise AMR2D_MAXNP in SELF_Refinement.cpp to extend).'
      stop 1
    endif

    call gpuCheck(hipMalloc(this%Pmodal_gpu,sizeof(this%Pmodal)))
    call gpuCheck(hipMalloc(this%indicator_gpu,sizeof(this%indicator)))
    call gpuCheck(hipMalloc(this%flag_gpu,sizeof(this%flag)))

    call this%UpdateDevice()

  endsubroutine Init_RefinementIndicator2D

  subroutine Free_RefinementIndicator2D(this)
    implicit none
    class(RefinementIndicator2D),intent(inout) :: this

    if(c_associated(this%Pmodal_gpu)) call gpuCheck(hipFree(this%Pmodal_gpu))
    if(c_associated(this%indicator_gpu)) call gpuCheck(hipFree(this%indicator_gpu))
    if(c_associated(this%flag_gpu)) call gpuCheck(hipFree(this%flag_gpu))
    this%Pmodal_gpu = c_null_ptr
    this%indicator_gpu = c_null_ptr
    this%flag_gpu = c_null_ptr

    call Free_RefinementIndicator2D_t(this)

  endsubroutine Free_RefinementIndicator2D

  subroutine UpdateHost_RefinementIndicator2D(this)
    implicit none
    class(RefinementIndicator2D),intent(inout) :: this

    call gpuCheck(hipMemcpy(c_loc(this%indicator),this%indicator_gpu, &
                            sizeof(this%indicator),hipMemcpyDeviceToHost))
    call gpuCheck(hipMemcpy(c_loc(this%flag),this%flag_gpu, &
                            sizeof(this%flag),hipMemcpyDeviceToHost))

  endsubroutine UpdateHost_RefinementIndicator2D

  subroutine UpdateDevice_RefinementIndicator2D(this)
    implicit none
    class(RefinementIndicator2D),intent(inout) :: this

    call gpuCheck(hipMemcpy(this%Pmodal_gpu,c_loc(this%Pmodal), &
                            sizeof(this%Pmodal),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%indicator_gpu,c_loc(this%indicator), &
                            sizeof(this%indicator),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%flag_gpu,c_loc(this%flag), &
                            sizeof(this%flag),hipMemcpyHostToDevice))

  endsubroutine UpdateDevice_RefinementIndicator2D

  subroutine Estimate_RefinementIndicator2D(this,solution,ivar)
    !! GPU path: one device thread per element computes the modal transform and energy
    !! ratios, writes the indicator and flag on device, then copies both back to the host so
    !! the (host-side) mesh-adaptation logic can consume them directly.
    implicit none
    class(RefinementIndicator2D),intent(inout) :: this
    class(Scalar2D),intent(in) :: solution
    integer,intent(in) :: ivar

    if(solution%N /= this%N) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : solution degree does not match indicator degree.'
      stop 1
    endif
    if(solution%nElem /= this%nElem) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : solution element count does not match indicator element count.'
      stop 1
    endif
    if(ivar < 0 .or. ivar > solution%nVar) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : driving-variable index out of range.'
      stop 1
    endif

    call RefinementIndicator_2D_gpu(this%Pmodal_gpu,solution%interior_gpu, &
                                    this%indicator_gpu,this%flag_gpu, &
                                    this%refineThreshold,this%coarsenThreshold, &
                                    this%N,solution%nVar,ivar,this%nElem)

    call this%UpdateHost()

  endsubroutine Estimate_RefinementIndicator2D

endmodule SELF_RefinementIndicator_2D
