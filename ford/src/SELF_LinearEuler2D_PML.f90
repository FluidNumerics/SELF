! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2026 Fluid Numerics LLC
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

module self_LinearEuler2D_PML

  use self_LinearEuler2D_PML_t
  use SELF_GPU
  use SELF_BoundaryConditions

  implicit none

  type,extends(LinearEuler2D_PML_t) :: LinearEuler2D_PML
  contains
    procedure :: AdditionalInit => AdditionalInit_LinearEuler2D_PML
    procedure :: boundaryflux => boundaryflux_LinearEuler2D_PML
    procedure :: fluxmethod => fluxmethod_LinearEuler2D_PML
    procedure :: sourcemethod => sourcemethod_LinearEuler2D_PML

  endtype LinearEuler2D_PML

  interface
    subroutine hbc2d_nonormalflow_lineareuler2d_pml_gpu(extboundary,boundary,nhat, &
                                                        elements,sides,nBoundaries,N,nel) &
      bind(c,name="hbc2d_nonormalflow_lineareuler2d_pml_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,nhat,elements,sides
      integer(c_int),value :: nBoundaries,N,nel
    endsubroutine hbc2d_nonormalflow_lineareuler2d_pml_gpu
  endinterface

  interface
    subroutine hbc2d_radiation_lineareuler2d_pml_gpu(extboundary,boundary, &
                                                     elements,sides,nBoundaries,N,nel) &
      bind(c,name="hbc2d_radiation_lineareuler2d_pml_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,elements,sides
      integer(c_int),value :: nBoundaries,N,nel
    endsubroutine hbc2d_radiation_lineareuler2d_pml_gpu
  endinterface

  interface
    subroutine fluxmethod_LinearEuler2D_PML_gpu(solution,flux,rho0,N,nel,nvar) &
      bind(c,name="fluxmethod_LinearEuler2D_PML_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution,flux
      real(c_prec),value :: rho0
      integer(c_int),value :: N,nel,nvar
    endsubroutine fluxmethod_LinearEuler2D_PML_gpu
  endinterface

  interface
    subroutine boundaryflux_LinearEuler2D_PML_gpu(fb,fextb,nhat,nscale,flux,rho0,N,nel,nvar) &
      bind(c,name="boundaryflux_LinearEuler2D_PML_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fb,fextb,flux,nhat,nscale
      real(c_prec),value :: rho0
      integer(c_int),value :: N,nel,nvar
    endsubroutine boundaryflux_LinearEuler2D_PML_gpu
  endinterface

  interface
    subroutine sourcemethod_LinearEuler2D_PML_gpu(source,solution,sigma_x,sigma_y,N,nel,nvar) &
      bind(c,name="sourcemethod_LinearEuler2D_PML_gpu")
      use iso_c_binding
      type(c_ptr),value :: source,solution,sigma_x,sigma_y
      integer(c_int),value :: N,nel,nvar
    endsubroutine sourcemethod_LinearEuler2D_PML_gpu
  endinterface

contains

  subroutine AdditionalInit_LinearEuler2D_PML(this)
    implicit none
    class(LinearEuler2D_PML),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    ! Initialise the parent (allocates sigma_x, sigma_y; registers
    ! the CPU PML BC handlers). We immediately overwrite the BC
    ! method pointers with GPU-resident wrappers below.
    call AdditionalInit_LinearEuler2D_PML_t(this)

    bcfunc => hbc2d_NoNormalFlow_LinearEuler2D_PML_GPU_wrapper
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

    bcfunc => hbc2d_Radiation_LinearEuler2D_PML_GPU_wrapper
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_RADIATION,"radiation",bcfunc)

  endsubroutine AdditionalInit_LinearEuler2D_PML

  subroutine hbc2d_NoNormalFlow_LinearEuler2D_PML_GPU_wrapper(bc,mymodel)
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(LinearEuler2D_PML)
      if(bc%nBoundaries > 0) then
        call hbc2d_nonormalflow_lineareuler2d_pml_gpu( &
          m%solution%extBoundary_gpu, &
          m%solution%boundary_gpu, &
          m%geometry%nhat%boundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N,m%solution%nElem)
      endif
    endselect

  endsubroutine hbc2d_NoNormalFlow_LinearEuler2D_PML_GPU_wrapper

  subroutine hbc2d_Radiation_LinearEuler2D_PML_GPU_wrapper(bc,mymodel)
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(LinearEuler2D_PML)
      if(bc%nBoundaries > 0) then
        call hbc2d_radiation_lineareuler2d_pml_gpu( &
          m%solution%extBoundary_gpu, &
          m%solution%boundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N,m%solution%nElem)
      endif
    endselect

  endsubroutine hbc2d_Radiation_LinearEuler2D_PML_GPU_wrapper

  subroutine sourcemethod_LinearEuler2D_PML(this)
    implicit none
    class(LinearEuler2D_PML),intent(inout) :: this

    call sourcemethod_LinearEuler2D_PML_gpu( &
      this%source%interior_gpu, &
      this%solution%interior_gpu, &
      this%sigma_x%interior_gpu, &
      this%sigma_y%interior_gpu, &
      this%solution%interp%N,this%solution%nelem,this%solution%nvar)

  endsubroutine sourcemethod_LinearEuler2D_PML

  subroutine boundaryflux_LinearEuler2D_PML(this)
    implicit none
    class(LinearEuler2D_PML),intent(inout) :: this

    call boundaryflux_LinearEuler2D_PML_gpu(this%solution%boundary_gpu, &
                                            this%solution%extBoundary_gpu, &
                                            this%geometry%nhat%boundary_gpu, &
                                            this%geometry%nscale%boundary_gpu, &
                                            this%flux%boundarynormal_gpu, &
                                            this%rho0,this%solution%interp%N, &
                                            this%solution%nelem,this%solution%nvar)

  endsubroutine boundaryflux_LinearEuler2D_PML

  subroutine fluxmethod_LinearEuler2D_PML(this)
    implicit none
    class(LinearEuler2D_PML),intent(inout) :: this

    call fluxmethod_LinearEuler2D_PML_gpu(this%solution%interior_gpu, &
                                          this%flux%interior_gpu, &
                                          this%rho0,this%solution%interp%N,this%solution%nelem, &
                                          this%solution%nvar)

  endsubroutine fluxmethod_LinearEuler2D_PML

endmodule self_LinearEuler2D_PML
