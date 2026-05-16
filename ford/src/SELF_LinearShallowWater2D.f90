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

module self_LinearShallowWater2D

  use self_LinearShallowWater2D_t
  use SELF_GPU
  use SELF_BoundaryConditions

  implicit none

  type,extends(LinearShallowWater2D_t) :: LinearShallowWater2D
  contains
    procedure :: AdditionalInit => AdditionalInit_LinearShallowWater2D
    procedure :: boundaryflux => boundaryflux_LinearShallowWater2D
    procedure :: fluxmethod => fluxmethod_LinearShallowWater2D
    procedure :: sourcemethod => sourcemethod_LinearShallowWater2D

  endtype LinearShallowWater2D

  interface
    subroutine hbc2d_nonormalflow_linearshallowwater2d_gpu(extboundary,boundary,nhat, &
                                                           elements,sides,nBoundaries,N,nel) &
      bind(c,name="hbc2d_nonormalflow_linearshallowwater2d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,nhat,elements,sides
      integer(c_int),value :: nBoundaries,N,nel
    endsubroutine hbc2d_nonormalflow_linearshallowwater2d_gpu
  endinterface

  interface
    subroutine hbc2d_radiation_linearshallowwater2d_gpu(extboundary, &
                                                        elements,sides,nBoundaries,N,nel) &
      bind(c,name="hbc2d_radiation_linearshallowwater2d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,elements,sides
      integer(c_int),value :: nBoundaries,N,nel
    endsubroutine hbc2d_radiation_linearshallowwater2d_gpu
  endinterface

  interface
    subroutine boundaryflux_LinearShallowWater2D_gpu(fb,fextb,nhat,nscale,flux,g,H,N,nel,nvar) &
      bind(c,name="boundaryflux_LinearShallowWater2D_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fb,fextb,flux,nhat,nscale
      real(c_prec),value :: g,H
      integer(c_int),value :: N,nel,nvar
    endsubroutine boundaryflux_LinearShallowWater2D_gpu
  endinterface

  interface
    subroutine fluxmethod_LinearShallowWater2D_gpu(solution,flux,g,H,N,nel,nvar) &
      bind(c,name="fluxmethod_LinearShallowWater2D_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution,flux
      real(c_prec),value :: g,H
      integer(c_int),value :: N,nel,nvar
    endsubroutine fluxmethod_LinearShallowWater2D_gpu
  endinterface

  interface
    subroutine sourcemethod_LinearShallowWater2D_gpu(solution,source,fCori,Cd,N,nel,nvar) &
      bind(c,name="sourcemethod_LinearShallowWater2D_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution,source,fCori
      real(c_prec),value :: Cd
      integer(c_int),value :: N,nel,nvar
    endsubroutine sourcemethod_LinearShallowWater2D_gpu
  endinterface

contains

  subroutine AdditionalInit_LinearShallowWater2D(this)
    implicit none
    class(LinearShallowWater2D),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    ! Call parent _t AdditionalInit (registers CPU BCs + initializes fCori)
    call AdditionalInit_LinearShallowWater2D_t(this)

    ! Re-register with GPU-accelerated versions
    bcfunc => hbc2d_NoNormalFlow_LinearShallowWater2D_GPU_wrapper
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

    bcfunc => hbc2d_Radiation_LinearShallowWater2D_GPU_wrapper
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_RADIATION,"radiation",bcfunc)

  endsubroutine AdditionalInit_LinearShallowWater2D

  subroutine hbc2d_NoNormalFlow_LinearShallowWater2D_GPU_wrapper(bc,mymodel)
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(LinearShallowWater2D)
      if(bc%nBoundaries > 0) then
        call hbc2d_nonormalflow_linearshallowwater2d_gpu( &
          m%solution%extBoundary_gpu, &
          m%solution%boundary_gpu, &
          m%geometry%nhat%boundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N,m%solution%nElem)
      endif
    endselect

  endsubroutine hbc2d_NoNormalFlow_LinearShallowWater2D_GPU_wrapper

  subroutine hbc2d_Radiation_LinearShallowWater2D_GPU_wrapper(bc,mymodel)
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(LinearShallowWater2D)
      if(bc%nBoundaries > 0) then
        call hbc2d_radiation_linearshallowwater2d_gpu( &
          m%solution%extBoundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N,m%solution%nElem)
      endif
    endselect

  endsubroutine hbc2d_Radiation_LinearShallowWater2D_GPU_wrapper

  subroutine boundaryflux_LinearShallowWater2D(this)
    implicit none
    class(LinearShallowWater2D),intent(inout) :: this

    call boundaryflux_LinearShallowWater2D_gpu(this%solution%boundary_gpu, &
                                               this%solution%extBoundary_gpu, &
                                               this%geometry%nhat%boundary_gpu, &
                                               this%geometry%nscale%boundary_gpu, &
                                               this%flux%boundaryNormal_gpu, &
                                               this%g,this%H, &
                                               this%solution%interp%N, &
                                               this%solution%nelem, &
                                               this%solution%nvar)

  endsubroutine boundaryflux_LinearShallowWater2D

  subroutine fluxmethod_LinearShallowWater2D(this)
    implicit none
    class(LinearShallowWater2D),intent(inout) :: this

    call fluxmethod_LinearShallowWater2D_gpu(this%solution%interior_gpu, &
                                             this%flux%interior_gpu, &
                                             this%g,this%H, &
                                             this%solution%interp%N, &
                                             this%solution%nelem, &
                                             this%solution%nvar)

  endsubroutine fluxmethod_LinearShallowWater2D

  subroutine sourcemethod_LinearShallowWater2D(this)
    implicit none
    class(LinearShallowWater2D),intent(inout) :: this

    call sourcemethod_LinearShallowWater2D_gpu(this%solution%interior_gpu, &
                                               this%source%interior_gpu, &
                                               this%fCori%interior_gpu, &
                                               this%Cd, &
                                               this%solution%interp%N, &
                                               this%solution%nelem, &
                                               this%solution%nvar)

  endsubroutine sourcemethod_LinearShallowWater2D

endmodule self_LinearShallowWater2D
