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

module SELF_ECDGModel3D

  use SELF_ECDGModel3D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(ECDGModel3D_t),public :: ECDGModel3D

  contains

    procedure :: CalculateTendency => CalculateTendency_ECDGModel3D

  endtype ECDGModel3D

contains

  subroutine CalculateTendency_ECDGModel3D(this)
    implicit none
    class(ECDGModel3D),intent(inout) :: this
    ! Local
    integer :: ndof

    call this%solution%BoundaryInterp()
    call this%solution%SideExchange(this%mesh)

    call this%PreTendency()
    call this%SetBoundaryCondition()

    if(this%gradient_enabled) then
      call this%CalculateSolutionGradient()
      call this%SetGradientBoundaryCondition()
      call this%solutionGradient%AverageSides()
    endif

    call this%SourceMethod()
    call this%BoundaryFlux()

    call this%TwoPointFluxMethod()

    call this%twoPointFlux%MappedDivergence(this%fluxDivergence%interior_gpu)

    call ECDGSurfaceContribution_3D_gpu( &
      this%flux%boundarynormal_gpu, &
      this%geometry%J%interior_gpu, &
      this%solution%interp%bMatrix_gpu, &
      this%solution%interp%qWeights_gpu, &
      this%fluxDivergence%interior_gpu, &
      this%solution%interp%N,this%solution%nVar,this%mesh%nElem)

    ndof = this%solution%nVar* &
           this%mesh%nElem* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)

    call CalculateDSDt_gpu(this%fluxDivergence%interior_gpu, &
                           this%source%interior_gpu, &
                           this%dSdt%interior_gpu,ndof)

  endsubroutine CalculateTendency_ECDGModel3D

endmodule SELF_ECDGModel3D
