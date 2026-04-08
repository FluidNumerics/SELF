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

module SELF_ECAdvection3D

  use SELF_ECAdvection3D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use SELF_BoundaryConditions
  use iso_c_binding

  implicit none

  type,extends(ECAdvection3D_t),public :: ECAdvection3D

  contains

    procedure :: AdditionalInit => AdditionalInit_ECAdvection3D
    procedure :: BoundaryFlux => BoundaryFlux_ECAdvection3D
    procedure :: TwoPointFluxMethod => TwoPointFluxMethod_ECAdvection3D
    procedure :: SourceMethod => SourceMethod_ECAdvection3D

  endtype ECAdvection3D

  interface
    subroutine hbc3d_mirror_ecadvection3d_gpu(extboundary,boundary, &
                                              elements,sides,nBoundaries,N,nel,nvar) &
      bind(c,name="hbc3d_mirror_ecadvection3d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,elements,sides
      integer(c_int),value :: nBoundaries,N,nel,nvar
    endsubroutine hbc3d_mirror_ecadvection3d_gpu
  endinterface

  interface
    subroutine boundaryflux_ecadvection3d_gpu(fb,fextb,nhat,nscale,flux,u,v,w,N,nel,nvar) &
      bind(c,name="boundaryflux_ecadvection3d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fb,fextb,nhat,nscale,flux
      real(c_prec),value :: u,v,w
      integer(c_int),value :: N,nel,nvar
    endsubroutine boundaryflux_ecadvection3d_gpu
  endinterface

  interface
    subroutine twopointfluxmethod_ecadvection3d_gpu(f,s,dsdx,u,v,w,N,nvar,nel) &
      bind(c,name="twopointfluxmethod_ecadvection3d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: f,s,dsdx
      real(c_prec),value :: u,v,w
      integer(c_int),value :: N,nvar,nel
    endsubroutine twopointfluxmethod_ecadvection3d_gpu
  endinterface

contains

  subroutine AdditionalInit_ECAdvection3D(this)
    implicit none
    class(ECAdvection3D),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    ! Call parent _t AdditionalInit (registers CPU mirror BC)
    call AdditionalInit_ECAdvection3D_t(this)

    ! Re-register with GPU-accelerated version
    bcfunc => hbc3d_Mirror_ECAdvection3D_GPU_wrapper
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

  endsubroutine AdditionalInit_ECAdvection3D

  subroutine hbc3d_Mirror_ECAdvection3D_GPU_wrapper(bc,mymodel)
    !! GPU-accelerated mirror BC for 3D EC Advection.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(ECAdvection3D)
      if(bc%nBoundaries > 0) then
        call hbc3d_mirror_ecadvection3d_gpu( &
          m%solution%extBoundary_gpu, &
          m%solution%boundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N, &
          m%solution%nElem,m%solution%nvar)
      endif
    endselect

  endsubroutine hbc3d_Mirror_ECAdvection3D_GPU_wrapper

  subroutine BoundaryFlux_ECAdvection3D(this)
    !! LLF Riemann flux on GPU — fully device-resident.
    implicit none
    class(ECAdvection3D),intent(inout) :: this

    call boundaryflux_ecadvection3d_gpu( &
      this%solution%boundary_gpu, &
      this%solution%extboundary_gpu, &
      this%geometry%nhat%boundary_gpu, &
      this%geometry%nscale%boundary_gpu, &
      this%flux%boundarynormal_gpu, &
      this%u,this%v,this%w, &
      this%solution%interp%N, &
      this%solution%nelem, &
      this%solution%nvar)

  endsubroutine BoundaryFlux_ECAdvection3D

  subroutine TwoPointFluxMethod_ECAdvection3D(this)
    !! Contravariant EC two-point flux on GPU — fully device-resident.
    implicit none
    class(ECAdvection3D),intent(inout) :: this

    call twopointfluxmethod_ecadvection3d_gpu( &
      this%twoPointFlux%interior_gpu, &
      this%solution%interior_gpu, &
      this%geometry%dsdx%interior_gpu, &
      this%u,this%v,this%w, &
      this%solution%interp%N, &
      this%solution%nvar, &
      this%solution%nelem)

  endsubroutine TwoPointFluxMethod_ECAdvection3D

  subroutine SourceMethod_ECAdvection3D(this)
    !! No source term — upload the zero-initialised host array to device.
    implicit none
    class(ECAdvection3D),intent(inout) :: this

    call gpuCheck(hipMemcpy(this%source%interior_gpu, &
                            c_loc(this%source%interior), &
                            sizeof(this%source%interior), &
                            hipMemcpyHostToDevice))

  endsubroutine SourceMethod_ECAdvection3D

endmodule SELF_ECAdvection3D
