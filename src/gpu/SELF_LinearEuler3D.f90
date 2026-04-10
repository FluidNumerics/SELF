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

module self_LinearEuler3D

  use self_LinearEuler3D_t
  use SELF_GPU
  use SELF_BoundaryConditions

  implicit none

  type,extends(LinearEuler3D_t) :: LinearEuler3D
  contains
    procedure :: AdditionalInit => AdditionalInit_LinearEuler3D
    procedure :: boundaryflux => boundaryflux_LinearEuler3D
    procedure :: fluxmethod => fluxmethod_LinearEuler3D

  endtype LinearEuler3D

  interface
    subroutine hbc3d_radiation_lineareuler3d_gpu(extboundary, &
                                                 elements,sides,nBoundaries,N,nel) &
      bind(c,name="hbc3d_radiation_lineareuler3d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,elements,sides
      integer(c_int),value :: nBoundaries,N,nel
    endsubroutine hbc3d_radiation_lineareuler3d_gpu
  endinterface

  interface
    subroutine fluxmethod_LinearEuler3D_gpu(solution,flux,rho0,c,N,nel,nvar) &
      bind(c,name="fluxmethod_LinearEuler3D_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution,flux
      real(c_prec),value :: rho0,c
      integer(c_int),value :: N,nel,nvar
    endsubroutine fluxmethod_LinearEuler3D_gpu
  endinterface

  interface
    subroutine boundaryflux_LinearEuler3D_gpu(fb,fextb,nhat,nscale,flux,rho0,c,N,nel) &
      bind(c,name="boundaryflux_LinearEuler3D_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fb,fextb,flux,nhat,nscale
      real(c_prec),value :: rho0,c
      integer(c_int),value :: N,nel
    endsubroutine boundaryflux_LinearEuler3D_gpu
  endinterface

contains

  subroutine AdditionalInit_LinearEuler3D(this)
    implicit none
    class(LinearEuler3D),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    bcfunc => hbc3d_Radiation_LinearEuler3D_GPU_wrapper
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_RADIATION,"radiation",bcfunc)

  endsubroutine AdditionalInit_LinearEuler3D

  subroutine hbc3d_Radiation_LinearEuler3D_GPU_wrapper(bc,mymodel)
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(LinearEuler3D)
      if(bc%nBoundaries > 0) then
        call hbc3d_radiation_lineareuler3d_gpu( &
          m%solution%extBoundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N,m%solution%nElem)
      endif
    endselect

  endsubroutine hbc3d_Radiation_LinearEuler3D_GPU_wrapper

  subroutine boundaryflux_LinearEuler3D(this)
    implicit none
    class(LinearEuler3D),intent(inout) :: this

    call boundaryflux_LinearEuler3D_gpu(this%solution%boundary_gpu, &
                                        this%solution%extBoundary_gpu, &
                                        this%geometry%nhat%boundary_gpu, &
                                        this%geometry%nscale%boundary_gpu, &
                                        this%flux%boundarynormal_gpu, &
                                        this%rho0,this%c,this%solution%interp%N, &
                                        this%solution%nelem)

  endsubroutine boundaryflux_LinearEuler3D

  subroutine fluxmethod_LinearEuler3D(this)
    implicit none
    class(LinearEuler3D),intent(inout) :: this

    call fluxmethod_LinearEuler3D_gpu(this%solution%interior_gpu, &
                                      this%flux%interior_gpu, &
                                      this%rho0,this%c,this%solution%interp%N,this%solution%nelem, &
                                      this%solution%nvar)

  endsubroutine fluxmethod_LinearEuler3D

endmodule self_LinearEuler3D
