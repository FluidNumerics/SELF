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

module self_advection_diffusion_3d_t

  use self_model
  use self_dgmodel3d
  use self_mesh

  implicit none

  type,extends(dgmodel3d) :: advection_diffusion_3d_t
    real(prec) :: nu ! diffusion coefficient
    real(prec) :: u ! constant x-component of velocity
    real(prec) :: v ! constant y-component of velocity
    real(prec) :: w ! constant z-component of velocity

  contains

    procedure :: riemannflux3d => riemannflux3d_advection_diffusion_3d_t
    procedure :: flux3d => flux3d_advection_diffusion_3d_t
    procedure :: entropy_func => entropy_func_advection_diffusion_3d_t

  endtype advection_diffusion_3d_t

contains

  pure function entropy_func_advection_diffusion_3d_t(this,s) result(e)
    class(advection_diffusion_3d_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%solution%nvar)
    real(prec) :: e
! Local
    integer :: ivar

    e = 0.0_prec
    do ivar = 1,this%solution%nvar
      e = e+0.5_prec*s(ivar)*s(ivar)
    enddo

  endfunction entropy_func_advection_diffusion_3d_t

  pure function flux3d_advection_diffusion_3d_t(this,s,dsdx) result(flux)
    class(advection_diffusion_3d_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%solution%nvar)
    real(prec),intent(in) :: dsdx(1:this%solution%nvar,1:3)
    real(prec) :: flux(1:this%solution%nvar,1:3)
! Local
    integer :: ivar

    do ivar = 1,this%solution%nvar
      flux(ivar,1) = this%u*s(ivar)-this%nu*dsdx(ivar,1) ! advective flux + diffusive flux
      flux(ivar,2) = this%v*s(ivar)-this%nu*dsdx(ivar,2) ! advective flux + diffusive flux
      flux(ivar,3) = this%w*s(ivar)-this%nu*dsdx(ivar,3) ! advective flux + diffusive flux
    enddo

  endfunction flux3d_advection_diffusion_3d_t

  pure function riemannflux3d_advection_diffusion_3d_t(this,sL,sR,dsdx,nhat) result(flux)
    class(advection_diffusion_3d_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec),intent(in) :: nhat(1:3)
    real(prec) :: flux(1:this%nvar)
! Local
    integer :: ivar
    real(prec) :: un,dsdn

    un = this%u*nhat(1)+this%v*nhat(2)+this%w*nhat(3)

    do ivar = 1,this%nvar
      dsdn = dsdx(ivar,1)*nhat(1)+dsdx(ivar,2)*nhat(2)+dsdx(ivar,3)*nhat(3)
      flux(ivar) = 0.5_prec*( &
                   (sL(ivar)+sR(ivar))+abs(un)*(sL(ivar)-sR(ivar)))- & ! advective flux
                   this%nu*dsdn
    enddo

  endfunction riemannflux3d_advection_diffusion_3d_t

endmodule self_advection_diffusion_3d_t
