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

module SELF_ECAdvection2D_t
  !! Entropy-conserving linear scalar advection equation in 2-D.
  !!
  !! Solves:
  !!   du/dt + u_a * du/dx + v_a * du/dy = 0
  !!
  !! Entropy function: eta(u) = u^2 / 2  (same as Trixi.jl)
  !!
  !! EC two-point flux (arithmetic mean, entropy-conserving for eta = u^2/2):
  !!   F^EC(uL, uR) = (u_a * (uL+uR)/2,  v_a * (uL+uR)/2)
  !!
  !! Surface Riemann flux (Godunov/upwind):
  !!   F_Riemann = 0.5 * (a.n) * (uL+uR) - 0.5 * |a.n| * (uR-uL)
  !! This is dissipative (entropy-stable) and reduces to the central flux
  !! when uL = uR (no dissipation at no-normal-flow or mirror boundaries).

  use SELF_ECDGModel2D_t
  use SELF_mesh

  implicit none

  type,extends(ECDGModel2D_t) :: ECAdvection2D_t

    real(prec) :: u ! x-component of advection velocity
    real(prec) :: v ! y-component of advection velocity

  contains

    procedure :: entropy_func => entropy_func_ECAdvection2D_t
    procedure :: twopointflux2d => twopointflux2d_ECAdvection2D_t
    procedure :: riemannflux2d => riemannflux2d_ECAdvection2D_t

    ! Mirror BC: sets extBoundary = interior boundary state so that
    ! sR = sL at all domain faces.  Used to eliminate boundary dissipation
    ! when testing entropy conservation.
    procedure :: hbc2d_NoNormalFlow => hbc2d_NoNormalFlow_ECAdvection2D_t

  endtype ECAdvection2D_t

contains

  pure function entropy_func_ECAdvection2D_t(this,s) result(e)
    !! Quadratic entropy: eta(u) = u^2 / 2
    class(ECAdvection2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: e

    e = 0.5_prec*s(1)*s(1)
    if(.false.) e = e+this%u ! suppress unused-dummy-argument warning

  endfunction entropy_func_ECAdvection2D_t

  pure function twopointflux2d_ECAdvection2D_t(this,sL,sR) result(flux)
    !! Arithmetic-mean two-point flux for linear advection.
    !! Entropy-conserving with respect to eta(u) = u^2/2.
    class(ECAdvection2D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec) :: flux(1:this%nvar,1:2)
    ! Local
    real(prec) :: savg

    savg = 0.5_prec*(sL(1)+sR(1))
    flux(1,1) = this%u*savg
    flux(1,2) = this%v*savg

  endfunction twopointflux2d_ECAdvection2D_t

  pure function riemannflux2d_ECAdvection2D_t(this,sL,sR,dsdx,nhat) result(flux)
    !! Godunov (upwind) Riemann flux for linear advection.
    !! Entropy-stable: dissipates entropy at outflow faces.
    class(ECAdvection2D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: flux(1:this%nvar)
    ! Local
    real(prec) :: un

    un = this%u*nhat(1)+this%v*nhat(2)
    flux(1) = 0.5_prec*(un*(sL(1)+sR(1))-abs(un)*(sR(1)-sL(1)))
    if(.false.) flux(1) = flux(1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction riemannflux2d_ECAdvection2D_t

  pure function hbc2d_NoNormalFlow_ECAdvection2D_t(this,s,nhat) result(exts)
    !! Mirror boundary condition: sets extBoundary = interior state.
    !! With the upwind Riemann flux, this gives sR = sL at the boundary,
    !! so the Riemann flux reduces to the central flux (a.n)*s — no upwind
    !! dissipation.  Use this BC when testing entropy conservation.
    class(ECAdvection2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: exts(1:this%nvar)

    exts = s
    if(.false.) exts(1) = exts(1)+nhat(1) ! suppress unused-dummy-argument warning

  endfunction hbc2d_NoNormalFlow_ECAdvection2D_t

endmodule SELF_ECAdvection2D_t
