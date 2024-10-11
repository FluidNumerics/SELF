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

module self_Burgers1D_t

  use self_model
  use self_dgmodel1d
  use self_mesh

  implicit none

  type,extends(dgmodel1d) :: Burgers1D_t
    ! Add any additional attributes here that are specific to your model
    real(prec) :: nu = 0.0_prec ! Diffusivity/viscosity

  contains
    procedure :: SetMetadata => SetMetadata_Burgers1D_t
    procedure :: entropy_func => entropy_func_Burgers1D_t
    procedure :: flux1d => flux1d_Burgers1D_t
    procedure :: riemannflux1d => riemannflux1d_Burgers1D_t

  endtype Burgers1D_t

contains
  subroutine SetMetadata_Burgers1D_t(this)
    implicit none
    class(Burgers1D_t),intent(inout) :: this

    call this%solution%SetName(1,"s")
    call this%solution%SetUnits(1,"[null]")

  endsubroutine SetMetadata_Burgers1D_t

  pure function entropy_func_Burgers1D_t(this,s) result(e)
    class(Burgers1D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%solution%nvar)
    real(prec) :: e

    e = 0.5_prec*s(1)*s(1)

  endfunction entropy_func_Burgers1D_t

  pure function flux1d_Burgers1D_t(this,s,dsdx) result(flux)
    class(Burgers1D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%solution%nvar)
    real(prec),intent(in) :: dsdx(1:this%solution%nvar)
    real(prec) :: flux(1:this%solution%nvar)

    flux(1) = 0.5_prec*s(1)*s(1)-this%nu*dsdx(1)

  endfunction flux1d_Burgers1D_t

  pure function riemannflux1d_Burgers1D_t(this,sL,sR,dsdx,nhat) result(flux)
    class(Burgers1D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%solution%nvar)
    real(prec),intent(in) :: sR(1:this%solution%nvar)
    real(prec),intent(in) :: dsdx(1:this%solution%nvar)
    real(prec),intent(in) :: nhat
    real(prec) :: flux(1:this%solution%nvar)
    ! Local
    real(prec) :: fL,fR,cmax

    ! Local Lax-Friedrich's flux
    fL = 0.5_prec*sL(1)*sL(1)*nhat
    fR = 0.5_prec*sR(1)*sR(1)*nhat
    cmax = max(abs(sL(1)*nhat),abs(sR(1)*nhat)) ! maximum wave speed

    flux(1) = 0.5_prec*(fL+fR)+cmax*(sL(1)-sR(1)) & ! advective flux
              -this%nu*dsdx(1)*nhat

  endfunction riemannflux1d_Burgers1D_t

endmodule self_Burgers1D_t
