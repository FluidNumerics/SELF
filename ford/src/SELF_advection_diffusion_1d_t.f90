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
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUsLESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARIsLG IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module self_advection_diffusion_1d_t

  use self_model
  use self_dgmodel1d
  use self_mesh

  implicit none

  type,extends(dgmodel1d) :: advection_diffusion_1d_t
    real(prec) :: nu ! diffusion coefficient
    real(prec) :: u ! constant velocity

  contains
    procedure :: riemannflux1d => riemannflux1d_advection_diffusion_1d_t
    procedure :: flux1d => flux1d_advection_diffusion_1d_t
    procedure :: entropy_func => entropy_func_advection_diffusion_1d_t

  endtype advection_diffusion_1d_t

contains

  pure function entropy_func_advection_diffusion_1d_t(this,s) result(e)
    class(advection_diffusion_1d_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%solution%nvar)
    real(prec) :: e
    ! Local
    integer :: ivar

    e = 0.0_prec
    do ivar = 1,this%solution%nvar
      e = e+0.5_prec*s(ivar)*s(ivar)
    enddo

  endfunction entropy_func_advection_diffusion_1d_t

  pure function riemannflux1d_advection_diffusion_1d_t(this,sL,sR,dsdx,nhat) result(flux)
    class(advection_diffusion_1d_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%solution%nvar)
    real(prec),intent(in) :: sR(1:this%solution%nvar)
    real(prec),intent(in) :: dsdx(1:this%solution%nvar)
    real(prec),intent(in) :: nhat
    real(prec) :: flux(1:this%solution%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%solution%nvar
      flux(ivar) = 0.5_prec*(this%u*nhat*(sL(ivar)+sR(ivar))+ &
                             abs(this%u*nhat)*(sL(ivar)-sR(ivar)))- & ! advective flux
                   this%nu*dsdx(ivar)*nhat ! diffusive flux
    enddo

  endfunction riemannflux1d_advection_diffusion_1d_t

  pure function flux1d_advection_diffusion_1d_t(this,s,dsdx) result(flux)
    class(advection_diffusion_1d_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%solution%nvar)
    real(prec),intent(in) :: dsdx(1:this%solution%nvar)
    real(prec) :: flux(1:this%solution%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%solution%nvar
      flux(ivar) = this%u*s(ivar)-this%nu*dsdx(ivar) ! advective flux + diffusive flux
    enddo

  endfunction flux1d_advection_diffusion_1d_t

endmodule self_advection_diffusion_1d_t
