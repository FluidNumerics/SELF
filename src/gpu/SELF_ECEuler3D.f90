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

module SELF_ECEuler3D

  use SELF_ECEuler3D_t
  use SELF_ECDGModel3D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use SELF_BoundaryConditions
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use iso_c_binding

  implicit none

  type,extends(ECEuler3D_t),public :: ECEuler3D

  contains

    procedure :: Init => Init_ECEuler3D
    procedure :: Free => Free_ECEuler3D
    procedure :: AdditionalInit => AdditionalInit_ECEuler3D
    procedure :: BoundaryFlux => BoundaryFlux_ECEuler3D
    procedure :: TwoPointFluxMethod => TwoPointFluxMethod_ECEuler3D
    procedure :: SourceMethod => SourceMethod_ECEuler3D

  endtype ECEuler3D

  interface
    subroutine hbc3d_nonormalflow_eceuler3d_gpu(extboundary,boundary, &
                                                nhat,elements,sides, &
                                                nBoundaries,N,nel,nvar) &
      bind(c,name="hbc3d_nonormalflow_eceuler3d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,nhat,elements,sides
      integer(c_int),value :: nBoundaries,N,nel,nvar
    endsubroutine hbc3d_nonormalflow_eceuler3d_gpu
  endinterface

  interface
    subroutine boundaryflux_eceuler3d_gpu(fb,fextb,nhat,nscale,flux, &
                                          p0,Rd,gamma,N,nel) &
      bind(c,name="boundaryflux_eceuler3d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fb,fextb,nhat,nscale,flux
      real(c_prec),value :: p0,Rd,gamma
      integer(c_int),value :: N,nel
    endsubroutine boundaryflux_eceuler3d_gpu
  endinterface

  interface
    subroutine twopointfluxmethod_eceuler3d_gpu(f,s,dsdx,p0,Rd,gamma, &
                                                N,nvar,nel) &
      bind(c,name="twopointfluxmethod_eceuler3d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: f,s,dsdx
      real(c_prec),value :: p0,Rd,gamma
      integer(c_int),value :: N,nvar,nel
    endsubroutine twopointfluxmethod_eceuler3d_gpu
  endinterface

  interface
    subroutine sourcemethod_eceuler3d_gpu(source,solution,g,N,nel) &
      bind(c,name="sourcemethod_eceuler3d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: source,solution
      real(c_prec),value :: g
      integer(c_int),value :: N,nel
    endsubroutine sourcemethod_eceuler3d_gpu
  endinterface

contains

  subroutine Init_ECEuler3D(this,mesh,geometry)
    implicit none
    class(ECEuler3D),intent(out) :: this
    type(Mesh3D),intent(in),target :: mesh
    type(SEMHex),intent(in),target :: geometry
    ! Local
    type(BoundaryCondition),pointer :: bc

    call Init_ECDGModel3D_t(this,mesh,geometry)

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

  endsubroutine Init_ECEuler3D

  subroutine Free_ECEuler3D(this)
    implicit none
    class(ECEuler3D),intent(inout) :: this
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

    call Free_ECDGModel3D_t(this)

  endsubroutine Free_ECEuler3D

  subroutine AdditionalInit_ECEuler3D(this)
    implicit none
    class(ECEuler3D),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    ! Call parent _t AdditionalInit (registers CPU BC)
    call AdditionalInit_ECEuler3D_t(this)

    ! Re-register with GPU-accelerated version
    bcfunc => hbc3d_NoNormalFlow_ECEuler3D_GPU_wrapper
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

  endsubroutine AdditionalInit_ECEuler3D

  subroutine hbc3d_NoNormalFlow_ECEuler3D_GPU_wrapper(bc,mymodel)
    !! GPU-accelerated no-normal-flow BC for 3D EC Euler.
    !! Reflects normal momentum, mirrors density and rho*theta.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(ECEuler3D)
      if(bc%nBoundaries > 0) then
        call hbc3d_nonormalflow_eceuler3d_gpu( &
          m%solution%extBoundary_gpu, &
          m%solution%boundary_gpu, &
          m%geometry%nhat%boundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N, &
          m%solution%nElem,m%solution%nvar)
      endif
    endselect

  endsubroutine hbc3d_NoNormalFlow_ECEuler3D_GPU_wrapper

  subroutine BoundaryFlux_ECEuler3D(this)
    !! LLF Riemann flux on GPU — fully device-resident.
    implicit none
    class(ECEuler3D),intent(inout) :: this

    call boundaryflux_eceuler3d_gpu( &
      this%solution%boundary_gpu, &
      this%solution%extboundary_gpu, &
      this%geometry%nhat%boundary_gpu, &
      this%geometry%nscale%boundary_gpu, &
      this%flux%boundarynormal_gpu, &
      this%p0,this%Rd,this%cp/this%cv, &
      this%solution%interp%N, &
      this%solution%nelem)

  endsubroutine BoundaryFlux_ECEuler3D

  subroutine TwoPointFluxMethod_ECEuler3D(this)
    !! Kennedy-Gruber two-point flux on GPU — fully device-resident.
    implicit none
    class(ECEuler3D),intent(inout) :: this

    call twopointfluxmethod_eceuler3d_gpu( &
      this%twoPointFlux%interior_gpu, &
      this%solution%interior_gpu, &
      this%geometry%dsdx%interior_gpu, &
      this%p0,this%Rd,this%cp/this%cv, &
      this%solution%interp%N, &
      this%solution%nvar, &
      this%solution%nelem)

  endsubroutine TwoPointFluxMethod_ECEuler3D

  subroutine SourceMethod_ECEuler3D(this)
    !! Gravitational source term on GPU — fully device-resident.
    implicit none
    class(ECEuler3D),intent(inout) :: this

    call sourcemethod_eceuler3d_gpu( &
      this%source%interior_gpu, &
      this%solution%interior_gpu, &
      this%g, &
      this%solution%interp%N, &
      this%solution%nelem)

  endsubroutine SourceMethod_ECEuler3D

endmodule SELF_ECEuler3D
