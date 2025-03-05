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

module self_advection_diffusion_2d

  use self_advection_diffusion_2d_t

  implicit none

  type,extends(advection_diffusion_2d_t) :: advection_diffusion_2d

  contains
    procedure :: setboundarycondition => setboundarycondition_advection_diffusion_2d
    procedure :: setgradientboundarycondition => setgradientboundarycondition_advection_diffusion_2d
    procedure :: boundaryflux => boundaryflux_advection_diffusion_2d
    procedure :: fluxmethod => fluxmethod_advection_diffusion_2d

  endtype advection_diffusion_2d

  interface
    subroutine setboundarycondition_advection_diffusion_2d_gpu(extboundary,boundary,sideinfo,N,nel,nvar) &
      bind(c,name="setboundarycondition_advection_diffusion_2d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,sideinfo
      integer(c_int),value :: N,nel,nvar
    endsubroutine setboundarycondition_advection_diffusion_2d_gpu
  endinterface

  interface
    subroutine setgradientboundarycondition_advection_diffusion_2d_gpu(extboundary,boundary,sideinfo,N,nel,nvar) &
      bind(c,name="setgradientboundarycondition_advection_diffusion_2d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,sideinfo
      integer(c_int),value :: N,nel,nvar
    endsubroutine setgradientboundarycondition_advection_diffusion_2d_gpu
  endinterface

  interface
    subroutine fluxmethod_advection_diffusion_2d_gpu(solution,solutiongradient,flux,u,v,nu,N,nel,nvar) &
      bind(c,name="fluxmethod_advection_diffusion_2d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution,solutiongradient,flux
      real(c_prec),value :: u,v,nu
      integer(c_int),value :: N,nel,nvar
    endsubroutine fluxmethod_advection_diffusion_2d_gpu
  endinterface

  interface
    subroutine boundaryflux_advection_diffusion_2d_gpu(fb,fextb,dfavg,nhat,nscale,flux,u,v,nu,N,nel,nvar) &
      bind(c,name="boundaryflux_advection_diffusion_2d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fb,fextb,dfavg,flux,nhat,nscale
      real(c_prec),value :: u,v,nu
      integer(c_int),value :: N,nel,nvar
    endsubroutine boundaryflux_advection_diffusion_2d_gpu
  endinterface

contains

  subroutine setboundarycondition_advection_diffusion_2d(this)
    !! Boundary conditions are set to periodic boundary conditions
    implicit none
    class(advection_diffusion_2d),intent(inout) :: this

    call setboundarycondition_advection_diffusion_2d_gpu(this%solution%extboundary_gpu, &
                                                         this%solution%boundary_gpu,this%mesh%sideInfo_gpu,this%solution%interp%N, &
                                                         this%solution%nelem,this%solution%nvar)

  endsubroutine setboundarycondition_advection_diffusion_2d

  subroutine setgradientboundarycondition_advection_diffusion_2d(this)
    !! Gradient boundary conditions are set to periodic boundary conditions
    implicit none
    class(advection_diffusion_2d),intent(inout) :: this

    call setgradientboundarycondition_advection_diffusion_2d_gpu( &
      this%solutiongradient%extboundary_gpu, &
      this%solutiongradient%boundary_gpu,this%mesh%sideInfo_gpu, &
      this%solution%interp%N,this%solution%nelem,this%solution%nvar)

  endsubroutine setgradientboundarycondition_advection_diffusion_2d

  subroutine fluxmethod_advection_diffusion_2d(this)
    implicit none
    class(advection_diffusion_2d),intent(inout) :: this

    call fluxmethod_advection_diffusion_2d_gpu(this%solution%interior_gpu, &
                                               this%solutiongradient%interior_gpu,this%flux%interior_gpu, &
                                               this%u,this%v,this%nu,this%solution%interp%N,this%solution%nelem, &
                                               this%solution%nvar)

  endsubroutine fluxmethod_advection_diffusion_2d

  subroutine boundaryflux_advection_diffusion_2d(this)
    ! this method uses an linear upwind solver for the
    ! advective flux and the bassi-rebay method for the
    ! diffusive fluxes
    implicit none
    class(advection_diffusion_2d),intent(inout) :: this

    call boundaryflux_advection_diffusion_2d_gpu(this%solution%boundary_gpu, &
                                                 this%solution%extBoundary_gpu,this%solutionGradient%avgBoundary_gpu, &
                                                 this%geometry%nhat%boundary_gpu,this%geometry%nscale%boundary_gpu, &
                                                 this%flux%boundarynormal_gpu,this%u,this%v,this%nu,this%solution%interp%N, &
                                                 this%solution%nelem,this%solution%nvar)

  endsubroutine boundaryflux_advection_diffusion_2d

endmodule self_advection_diffusion_2d
