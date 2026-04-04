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

module SELF_ECAdvection2D

  use SELF_ECAdvection2D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(ECAdvection2D_t),public :: ECAdvection2D

  contains

    procedure :: SetBoundaryCondition => SetBoundaryCondition_ECAdvection2D
    procedure :: BoundaryFlux => BoundaryFlux_ECAdvection2D
    procedure :: TwoPointFluxMethod => TwoPointFluxMethod_ECAdvection2D
    procedure :: SourceMethod => SourceMethod_ECAdvection2D

  endtype ECAdvection2D

  interface
    subroutine setboundarycondition_ecadvection2d_gpu(extboundary,boundary,sideinfo,N,nel,nvar) &
      bind(c,name="setboundarycondition_ecadvection2d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,sideinfo
      integer(c_int),value :: N,nel,nvar
    endsubroutine setboundarycondition_ecadvection2d_gpu
  endinterface

  interface
    subroutine boundaryflux_ecadvection2d_gpu(fb,fextb,nhat,nscale,flux,u,v,N,nel,nvar) &
      bind(c,name="boundaryflux_ecadvection2d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fb,fextb,nhat,nscale,flux
      real(c_prec),value :: u,v
      integer(c_int),value :: N,nel,nvar
    endsubroutine boundaryflux_ecadvection2d_gpu
  endinterface

  interface
    subroutine twopointfluxmethod_ecadvection2d_gpu(f,s,dsdx,u,v,N,nvar,nel) &
      bind(c,name="twopointfluxmethod_ecadvection2d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: f,s,dsdx
      real(c_prec),value :: u,v
      integer(c_int),value :: N,nvar,nel
    endsubroutine twopointfluxmethod_ecadvection2d_gpu
  endinterface

contains

  subroutine SetBoundaryCondition_ECAdvection2D(this)
    !! Mirror BC on GPU: extBoundary = boundary at all domain faces.
    implicit none
    class(ECAdvection2D),intent(inout) :: this

    call setboundarycondition_ecadvection2d_gpu( &
      this%solution%extboundary_gpu, &
      this%solution%boundary_gpu, &
      this%mesh%sideinfo_gpu, &
      this%solution%interp%N, &
      this%solution%nelem, &
      this%solution%nvar)

  endsubroutine SetBoundaryCondition_ECAdvection2D

  subroutine BoundaryFlux_ECAdvection2D(this)
    !! LLF Riemann flux on GPU — fully device-resident.
    implicit none
    class(ECAdvection2D),intent(inout) :: this

    call boundaryflux_ecadvection2d_gpu( &
      this%solution%boundary_gpu, &
      this%solution%extboundary_gpu, &
      this%geometry%nhat%boundary_gpu, &
      this%geometry%nscale%boundary_gpu, &
      this%flux%boundarynormal_gpu, &
      this%u,this%v, &
      this%solution%interp%N, &
      this%solution%nelem, &
      this%solution%nvar)

  endsubroutine BoundaryFlux_ECAdvection2D

  subroutine TwoPointFluxMethod_ECAdvection2D(this)
    !! Contravariant EC two-point flux on GPU — fully device-resident.
    implicit none
    class(ECAdvection2D),intent(inout) :: this

    call twopointfluxmethod_ecadvection2d_gpu( &
      this%twoPointFlux%interior_gpu, &
      this%solution%interior_gpu, &
      this%geometry%dsdx%interior_gpu, &
      this%u,this%v, &
      this%solution%interp%N, &
      this%solution%nvar, &
      this%solution%nelem)

  endsubroutine TwoPointFluxMethod_ECAdvection2D

  subroutine SourceMethod_ECAdvection2D(this)
    !! No source term — zero the device array without touching the host.
    implicit none
    class(ECAdvection2D),intent(inout) :: this

    call gpuCheck(hipMemset(this%source%interior_gpu, &
                            0, &
                            sizeof(this%source%interior)))

  endsubroutine SourceMethod_ECAdvection2D

endmodule SELF_ECAdvection2D
