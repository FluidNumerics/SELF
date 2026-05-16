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
  use SELF_ECDGModel2D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use SELF_BoundaryConditions
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use iso_c_binding

  implicit none

  type,extends(ECAdvection2D_t),public :: ECAdvection2D

  contains

    procedure :: Init => Init_ECAdvection2D
    procedure :: Free => Free_ECAdvection2D
    procedure :: AdditionalInit => AdditionalInit_ECAdvection2D
    procedure :: BoundaryFlux => BoundaryFlux_ECAdvection2D
    procedure :: TwoPointFluxMethod => TwoPointFluxMethod_ECAdvection2D
    procedure :: SourceMethod => SourceMethod_ECAdvection2D

  endtype ECAdvection2D

  interface
    subroutine hbc2d_mirror_ecadvection2d_gpu(extboundary,boundary, &
                                              elements,sides,nBoundaries,N,nel,nvar) &
      bind(c,name="hbc2d_mirror_ecadvection2d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,elements,sides
      integer(c_int),value :: nBoundaries,N,nel,nvar
    endsubroutine hbc2d_mirror_ecadvection2d_gpu
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

  subroutine Init_ECAdvection2D(this,mesh,geometry)
    !! Initialize EC Advection 2D, then upload BC element/side arrays to GPU.
    implicit none
    class(ECAdvection2D),intent(out) :: this
    type(Mesh2D),intent(in),target :: mesh
    type(SEMQuad),intent(in),target :: geometry
    ! Local
    type(BoundaryCondition),pointer :: bc

    call Init_ECDGModel2D_t(this,mesh,geometry)

    ! Upload hyperbolic BC element/side arrays to device
    bc => this%hyperbolicBCs%head
    do while(associated(bc))
      if(bc%nBoundaries > 0) then
        call gpuCheck(hipMalloc(bc%elements_gpu,sizeof(bc%elements)))
        call gpuCheck(hipMemcpy(bc%elements_gpu,c_loc(bc%elements), &
                                sizeof(bc%elements),hipMemcpyHostToDevice))
        call gpuCheck(hipMalloc(bc%sides_gpu,sizeof(bc%sides)))
        call gpuCheck(hipMemcpy(bc%sides_gpu,c_loc(bc%sides), &
                                sizeof(bc%sides),hipMemcpyHostToDevice))
      endif
      bc => bc%next
    enddo

  endsubroutine Init_ECAdvection2D

  subroutine Free_ECAdvection2D(this)
    !! Free EC Advection 2D, including GPU BC arrays.
    implicit none
    class(ECAdvection2D),intent(inout) :: this
    ! Local
    type(BoundaryCondition),pointer :: bc

    bc => this%hyperbolicBCs%head
    do while(associated(bc))
      if(c_associated(bc%elements_gpu)) call gpuCheck(hipFree(bc%elements_gpu))
      if(c_associated(bc%sides_gpu)) call gpuCheck(hipFree(bc%sides_gpu))
      bc%elements_gpu = c_null_ptr
      bc%sides_gpu = c_null_ptr
      bc => bc%next
    enddo

    call Free_ECDGModel2D_t(this)

  endsubroutine Free_ECAdvection2D

  subroutine AdditionalInit_ECAdvection2D(this)
    implicit none
    class(ECAdvection2D),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    ! Call parent _t AdditionalInit (registers CPU mirror BC)
    call AdditionalInit_ECAdvection2D_t(this)

    ! Re-register with GPU-accelerated version
    bcfunc => hbc2d_Mirror_ECAdvection2D_GPU_wrapper
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

  endsubroutine AdditionalInit_ECAdvection2D

  subroutine hbc2d_Mirror_ECAdvection2D_GPU_wrapper(bc,mymodel)
    !! GPU-accelerated mirror BC for 2D EC Advection.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(ECAdvection2D)
      if(bc%nBoundaries > 0) then
        call hbc2d_mirror_ecadvection2d_gpu( &
          m%solution%extBoundary_gpu, &
          m%solution%boundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N, &
          m%solution%nElem,m%solution%nvar)
      endif
    endselect

  endsubroutine hbc2d_Mirror_ECAdvection2D_GPU_wrapper

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
    !! No source term — upload the zero-initialised host array to device.
    implicit none
    class(ECAdvection2D),intent(inout) :: this

    call gpuCheck(hipMemcpy(this%source%interior_gpu, &
                            c_loc(this%source%interior), &
                            sizeof(this%source%interior), &
                            hipMemcpyHostToDevice))

  endsubroutine SourceMethod_ECAdvection2D

endmodule SELF_ECAdvection2D
