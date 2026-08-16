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

module SELF_RefinementIndicator_3D

  use SELF_Constants
  use SELF_RefinementIndicator_3D_t
  use SELF_Scalar_3D
  use SELF_GPU
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(RefinementIndicator3D_t),public :: RefinementIndicator3D
    character(3) :: backend = "gpu"
    type(c_ptr) :: Pmodal_gpu = c_null_ptr
    type(c_ptr) :: indicator_gpu = c_null_ptr
    type(c_ptr) :: flag_gpu = c_null_ptr
    type(c_ptr) :: gate_gpu = c_null_ptr
    type(c_ptr) :: energyWeight_gpu = c_null_ptr
    integer :: nVarWeights_gpu = 0
      !! Allocated length of energyWeight_gpu. Zero until the first Estimate fixes the solution's
      !! variable count, which is not known at Init.

  contains

    procedure,public :: Init => Init_RefinementIndicator3D
    procedure,public :: Free => Free_RefinementIndicator3D
    procedure,public :: UpdateHost => UpdateHost_RefinementIndicator3D
    procedure,public :: UpdateDevice => UpdateDevice_RefinementIndicator3D
    procedure,public :: Estimate => Estimate_RefinementIndicator3D

  endtype RefinementIndicator3D

contains

  subroutine Init_RefinementIndicator3D(this,interp,nElem,refineThreshold,coarsenThreshold)
    implicit none
    class(RefinementIndicator3D),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nElem
    real(prec),intent(in) :: refineThreshold
    real(prec),intent(in) :: coarsenThreshold

    ! Reuse the portable host initialization (allocation + modal transform build).
    call Init_RefinementIndicator3D_t(this,interp,nElem,refineThreshold,coarsenThreshold)

    ! The device kernel uses per-block shared scratch sized to (N+1) <= 12 (AMR3D_MAXNP in
    ! SELF_Refinement.cpp; two (N+1)^3 shared buffers must fit the 48 KB static shared-memory
    ! limit in double precision). Fail early rather than overrun that bound on the device.
    if(interp%N+1 > 12) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : GPU refinement indicator supports degree N <= 11 '// &
        '(raise AMR3D_MAXNP in SELF_Refinement.cpp to extend).'
      stop 1
    endif

    call gpuCheck(hipMalloc(this%Pmodal_gpu,sizeof(this%Pmodal)))
    call gpuCheck(hipMalloc(this%indicator_gpu,sizeof(this%indicator)))
    call gpuCheck(hipMalloc(this%flag_gpu,sizeof(this%flag)))
    call gpuCheck(hipMalloc(this%gate_gpu,sizeof(this%gate)))
    ! The gate weights are sized to the solution's variable count, which the first Estimate
    ! establishes; nothing is allocated for them here.
    this%energyWeight_gpu = c_null_ptr
    this%nVarWeights_gpu = 0

    call this%UpdateDevice()

  endsubroutine Init_RefinementIndicator3D

  subroutine Free_RefinementIndicator3D(this)
    implicit none
    class(RefinementIndicator3D),intent(inout) :: this

    if(c_associated(this%Pmodal_gpu)) call gpuCheck(hipFree(this%Pmodal_gpu))
    if(c_associated(this%indicator_gpu)) call gpuCheck(hipFree(this%indicator_gpu))
    if(c_associated(this%flag_gpu)) call gpuCheck(hipFree(this%flag_gpu))
    if(c_associated(this%gate_gpu)) call gpuCheck(hipFree(this%gate_gpu))
    if(c_associated(this%energyWeight_gpu)) call gpuCheck(hipFree(this%energyWeight_gpu))
    this%Pmodal_gpu = c_null_ptr
    this%indicator_gpu = c_null_ptr
    this%flag_gpu = c_null_ptr
    this%gate_gpu = c_null_ptr
    this%energyWeight_gpu = c_null_ptr
    this%nVarWeights_gpu = 0

    call Free_RefinementIndicator3D_t(this)

  endsubroutine Free_RefinementIndicator3D

  subroutine UpdateHost_RefinementIndicator3D(this)
    implicit none
    class(RefinementIndicator3D),intent(inout) :: this

    call gpuCheck(hipMemcpy(c_loc(this%indicator),this%indicator_gpu, &
                            sizeof(this%indicator),hipMemcpyDeviceToHost))
    call gpuCheck(hipMemcpy(c_loc(this%flag),this%flag_gpu, &
                            sizeof(this%flag),hipMemcpyDeviceToHost))
    call gpuCheck(hipMemcpy(c_loc(this%gate),this%gate_gpu, &
                            sizeof(this%gate),hipMemcpyDeviceToHost))

  endsubroutine UpdateHost_RefinementIndicator3D

  subroutine UpdateDevice_RefinementIndicator3D(this)
    implicit none
    class(RefinementIndicator3D),intent(inout) :: this

    call gpuCheck(hipMemcpy(this%Pmodal_gpu,c_loc(this%Pmodal), &
                            sizeof(this%Pmodal),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%indicator_gpu,c_loc(this%indicator), &
                            sizeof(this%indicator),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%flag_gpu,c_loc(this%flag), &
                            sizeof(this%flag),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%gate_gpu,c_loc(this%gate), &
                            sizeof(this%gate),hipMemcpyHostToDevice))

  endsubroutine UpdateDevice_RefinementIndicator3D

  subroutine Estimate_RefinementIndicator3D(this,solution,ivar,comm,gate)
    !! GPU path: one device block per element computes the modal transform, the raw smoothness
    !! ratio and the gate energy; those are copied back to the host, where the shared second phase
    !! applies the amplitude gate, the log10 and the thresholds.
    !!
    !! The second phase is host-side (rather than a second device kernel) because the amplitude
    !! gate needs a field-wide energy scale - a reduction, and an MPI collective when the mesh is
    !! decomposed - that an element-local device block cannot supply. Running it through the same
    !! FinalizeIndicator the portable backend uses also makes the flags identical on CPU and GPU
    !! by construction. The extra traffic is one per-element array each way, once per adaptation
    !! epoch; nothing is added to the time-stepping loop.
    implicit none
    class(RefinementIndicator3D),intent(inout) :: this
    class(Scalar3D),intent(in) :: solution
    integer,intent(in) :: ivar
    integer,intent(in),optional :: comm
    real(prec),intent(in),optional :: gate(:)
    ! Local
    integer :: iel
    real(prec) :: energyScale

    call CheckEstimateArguments(this,solution,ivar,gate)
    call ResolveEnergyWeights(this,solution%nVar,ivar)

    if(this%nVarWeights_gpu /= this%nVarWeights) then
      if(c_associated(this%energyWeight_gpu)) call gpuCheck(hipFree(this%energyWeight_gpu))
      call gpuCheck(hipMalloc(this%energyWeight_gpu,sizeof(this%energyWeight)))
      this%nVarWeights_gpu = this%nVarWeights
    endif
    call gpuCheck(hipMemcpy(this%energyWeight_gpu,c_loc(this%energyWeight), &
                            sizeof(this%energyWeight),hipMemcpyHostToDevice))

    ! First phase on the device. indicator_gpu transiently holds the RAW ratio S_e, not log10(S_e).
    call RefinementIndicator_3D_gpu(this%Pmodal_gpu,solution%interior_gpu, &
                                    this%energyWeight_gpu, &
                                    this%indicator_gpu,this%gate_gpu, &
                                    this%N,solution%nVar,ivar,this%nElem)

    call this%UpdateHost()

    if(present(gate)) then
      do iel = 1,this%nElem
        this%gate(iel) = gate(iel)
      enddo
    endif

    call ResolveEnergyScale(this,energyScale,comm)
    call FinalizeIndicator(this,energyScale)

    ! Re-establish the device mirrors of indicator/flag/gate, which the host phase has just
    ! rewritten, so the device copies are never left stale.
    call this%UpdateDevice()

  endsubroutine Estimate_RefinementIndicator3D

endmodule SELF_RefinementIndicator_3D
