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

module SELF_ECAdvection3D_t
  !! Entropy-conserving linear scalar advection equation in 3-D.
  !!
  !! Solves:
  !!   du/dt + u_a * du/dx + v_a * du/dy + w_a * du/dz = 0
  !!
  !! Entropy function: eta(u) = u^2 / 2
  !!
  !! EC two-point flux (arithmetic mean):
  !!   F^EC(uL, uR) = (u_a, v_a, w_a) * (uL+uR)/2
  !!
  !! Surface Riemann flux (Local Lax-Friedrichs / Rusanov):
  !!   F_Riemann = 0.5 * (a.n) * (uL+uR) - 0.5 * |a| * (uR-uL)
  !! where |a| = sqrt(u^2 + v^2 + w^2) is the maximum wave speed.

  use SELF_ECDGModel3D_t
  use SELF_mesh

  implicit none

  type,extends(ECDGModel3D_t) :: ECAdvection3D_t

    real(prec) :: u ! x-component of advection velocity
    real(prec) :: v ! y-component of advection velocity
    real(prec) :: w ! z-component of advection velocity

  contains

    procedure :: entropy_func => entropy_func_ECAdvection3D_t
    procedure :: twopointflux3d => twopointflux3d_ECAdvection3D_t
    procedure :: riemannflux3d => riemannflux3d_ECAdvection3D_t
    procedure :: hbc3d_NoNormalFlow => hbc3d_NoNormalFlow_ECAdvection3D_t

  endtype ECAdvection3D_t

contains

  pure function entropy_func_ECAdvection3D_t(this,s) result(e)
    !! Quadratic entropy: eta(u) = u^2 / 2
    class(ECAdvection3D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: e

    e = 0.5_prec*s(1)*s(1)
    if(.false.) e = e+this%u ! suppress unused-dummy-argument warning

  endfunction entropy_func_ECAdvection3D_t

  pure function twopointflux3d_ECAdvection3D_t(this,sL,sR) result(flux)
    !! Arithmetic-mean two-point flux for linear advection.
    !! Entropy-conserving with respect to eta(u) = u^2/2.
    class(ECAdvection3D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec) :: flux(1:this%nvar,1:3)
    ! Local
    real(prec) :: savg

    savg = 0.5_prec*(sL(1)+sR(1))
    flux(1,1) = this%u*savg
    flux(1,2) = this%v*savg
    flux(1,3) = this%w*savg

  endfunction twopointflux3d_ECAdvection3D_t

  pure function riemannflux3d_ECAdvection3D_t(this,sL,sR,dsdx,nhat) result(flux)
    !! Local Lax-Friedrichs (Rusanov) Riemann flux for linear advection.
    !! lambda_max = sqrt(u^2 + v^2 + w^2).
    class(ECAdvection3D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec),intent(in) :: nhat(1:3)
    real(prec) :: flux(1:this%nvar)
    ! Local
    real(prec) :: un,lam

    un = this%u*nhat(1)+this%v*nhat(2)+this%w*nhat(3)
    lam = sqrt(this%u*this%u+this%v*this%v+this%w*this%w)
    flux(1) = 0.5_prec*(un*(sL(1)+sR(1))-lam*(sR(1)-sL(1)))
    if(.false.) flux(1) = flux(1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction riemannflux3d_ECAdvection3D_t

  pure function hbc3d_NoNormalFlow_ECAdvection3D_t(this,s,nhat) result(exts)
    !! Mirror boundary condition: sets extBoundary = interior state.
    class(ECAdvection3D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: nhat(1:3)
    real(prec) :: exts(1:this%nvar)

    exts = s
    if(.false.) exts(1) = exts(1)+nhat(1) ! suppress unused-dummy-argument warning

  endfunction hbc3d_NoNormalFlow_ECAdvection3D_t

endmodule SELF_ECAdvection3D_t
