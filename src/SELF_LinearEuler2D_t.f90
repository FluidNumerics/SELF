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

module self_LinearEuler2D_t
!! This module defines a class that can be used to solve the Linear Euler
!! equations in 2-D. The Linear Euler Equations, here, are the Euler equations
!! linearized about a motionless background state.
!!
!! The conserved variables are

!! \begin{equation}
!! \vec{s} = \begin{pmatrix}
!!     \rho \\
!!      u \\
!!      v \\
!!      p
!!  \end{pmatrix}
!! \end{equation}
!!
!! The conservative flux is
!!
!! \begin{equation}
!! \overleftrightarrow{f} = \begin{pmatrix}
!!     \rho_0 u \hat{x} + \rho_0 v \hat{y} \\
!!      \frac{p}{\rho_0} \hat{x} \\
!!      \frac{p}{\rho_0} \hat{y} \\
!!      c^2 \rho_0 ( u \hat{x} + v \hat{y} )
!!  \end{pmatrix}
!! \end{equation}
!!
!! and the source terms are null.
!!

  use self_model
  use self_dgmodel2d
  use self_mesh

  implicit none

  type,extends(dgmodel2d) :: LinearEuler2D_t
    ! Add any additional attributes here that are specific to your model
    real(prec) :: rho0 = 1.0_prec ! Reference density
    real(prec) :: c = 1.0_prec ! Sound speed
    real(prec) :: g = 0.0_prec ! gravitational acceleration (y-direction only)

  contains
    procedure :: SetNumberOfVariables => SetNumberOfVariables_LinearEuler2D_t
    procedure :: SetMetadata => SetMetadata_LinearEuler2D_t
    procedure :: entropy_func => entropy_func_LinearEuler2D_t
    procedure :: hbc2d_NoNormalFlow => hbc2d_NoNormalFlow_LinearEuler2D_t
    procedure :: flux2d => flux2d_LinearEuler2D_t
    procedure :: riemannflux2d => riemannflux2d_LinearEuler2D_t
    !procedure :: source2d => source2d_LinearEuler2D_t
    procedure :: SphericalSoundWave => SphericalSoundWave_LinearEuler2D_t

  endtype LinearEuler2D_t

contains

  subroutine SetNumberOfVariables_LinearEuler2D_t(this)
    implicit none
    class(LinearEuler2D_t),intent(inout) :: this

    this%nvar = 4

  endsubroutine SetNumberOfVariables_LinearEuler2D_t

  subroutine SetMetadata_LinearEuler2D_t(this)
    implicit none
    class(LinearEuler2D_t),intent(inout) :: this

    call this%solution%SetName(1,"rho") ! Density
    call this%solution%SetUnits(1,"kg⋅m⁻³")

    call this%solution%SetName(2,"u") ! x-velocity component
    call this%solution%SetUnits(2,"m⋅s⁻¹")

    call this%solution%SetName(3,"v") ! y-velocity component
    call this%solution%SetUnits(3,"m⋅s⁻¹")

    call this%solution%SetName(4,"P") ! Pressure
    call this%solution%SetUnits(4,"kg⋅m⁻¹⋅s⁻²")

  endsubroutine SetMetadata_LinearEuler2D_t

  pure function entropy_func_LinearEuler2D_t(this,s) result(e)
    !! The entropy function is the sum of kinetic and internal energy
    !! For the linear model, this is
    !!
    !! \begin{equation}
    !!   e = \frac{1}{2} \left( \rho_0*( u^2 + v^2 ) + \frac{P^2}{\rho_0 c^2} \right)
    class(LinearEuler2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: e

    e = 0.5_prec*this%rho0*(s(2)*s(2)+s(3)*(3))+ &
        0.5_prec*(s(4)*s(4)/(this%rho0*this%c*this%c))

  endfunction entropy_func_LinearEuler2D_t

  pure function hbc2d_NoNormalFlow_LinearEuler2D_t(this,s,nhat) result(exts)
    class(LinearEuler2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: exts(1:this%nvar)
    ! Local
    integer :: ivar

    exts(1) = s(1) ! density
    exts(2) = (nhat(2)**2-nhat(1)**2)*s(2)-2.0_prec*nhat(1)*nhat(2)*s(3) ! u
    exts(3) = (nhat(1)**2-nhat(2)**2)*s(3)-2.0_prec*nhat(1)*nhat(2)*s(2) ! v
    exts(4) = s(4) ! p

  endfunction hbc2d_NoNormalFlow_LinearEuler2D_t

  pure function flux2d_LinearEuler2D_t(this,s,dsdx) result(flux)
    class(LinearEuler2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec) :: flux(1:this%nvar,1:2)

    flux(1,1) = this%rho0*s(2) ! density, x flux ; rho0*u
    flux(1,2) = this%rho0*s(3) ! density, y flux ; rho0*v
    flux(2,1) = s(4)/this%rho0 ! x-velocity, x flux; p/rho0
    flux(2,2) = 0.0_prec ! x-velocity, y flux; 0
    flux(3,1) = 0.0_prec ! y-velocity, x flux; 0
    flux(3,2) = s(4)/this%rho0 ! y-velocity, y flux; p/rho0
    flux(4,1) = this%c*this%c*this%rho0*s(2) ! pressure, x flux : rho0*c^2*u
    flux(4,2) = this%c*this%c*this%rho0*s(3) ! pressure, y flux : rho0*c^2*v

  endfunction flux2d_LinearEuler2D_t

  pure function riemannflux2d_LinearEuler2D_t(this,sL,sR,dsdx,nhat) result(flux)
    !! Uses a local lax-friedrich's upwind flux
    !! The max eigenvalue is taken as the sound speed
    class(LinearEuler2D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: flux(1:this%nvar)
    ! Local
    real(prec) :: fL(1:this%nvar)
    real(prec) :: fR(1:this%nvar)
    real(prec) :: u,v,p,c,rho0

    u = sL(2)
    v = sL(3)
    p = sL(4)
    rho0 = this%rho0
    c = this%c
    fL(1) = rho0*(u*nhat(1)+v*nhat(2)) ! density
    fL(2) = p*nhat(1)/rho0 ! u
    fL(3) = p*nhat(2)/rho0 ! v
    fL(4) = rho0*c*c*(u*nhat(1)+v*nhat(2)) ! pressure

    u = sR(2)
    v = sR(3)
    p = sR(4)
    fR(1) = rho0*(u*nhat(1)+v*nhat(2)) ! density
    fR(2) = p*nhat(1)/rho0 ! u
    fR(3) = p*nhat(2)/rho0 ! v
    fR(4) = rho0*c*c*(u*nhat(1)+v*nhat(2)) ! pressure

    flux(1:4) = 0.5_prec*(fL(1:4)+fR(1:4))+c*(sL(1:4)-sR(1:4))

  endfunction riemannflux2d_LinearEuler2D_t

  subroutine SphericalSoundWave_LinearEuler2D_t(this,rhoprime,Lr,x0,y0)
    !! This subroutine sets the initial condition for a weak blast wave
    !! problem. The initial condition is given by
    !!
    !! \begin{equation}
    !! \begin{aligned}
    !! \rho &= \rho_0 + \rho' \exp\left( -\ln(2) \frac{(x-x_0)^2 + (y-y_0)^2}{L_r^2} \right)
    !! u &= 0 \\
    !! v &= 0 \\
    !! E &= \frac{P_0}{\gamma - 1} + E \exp\left( -\ln(2) \frac{(x-x_0)^2 + (y-y_0)^2}{L_e^2} \right)
    !! \end{aligned}
    !! \end{equation}
    !!
    implicit none
    class(LinearEuler2D_t),intent(inout) :: this
    real(prec),intent(in) ::  rhoprime,Lr,x0,y0
    ! Local
    integer :: i,j,iEl
    real(prec) :: x,y,rho,r,E

    print*,__FILE__," : Configuring weak blast wave initial condition. "
    print*,__FILE__," : rhoprime = ",rhoprime
    print*,__FILE__," : Lr = ",Lr
    print*,__FILE__," : x0 = ",x0
    print*,__FILE__," : y0 = ",y0

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem)
      x = this%geometry%x%interior(i,j,iEl,1,1)-x0
      y = this%geometry%x%interior(i,j,iEl,1,2)-y0
      r = sqrt(x**2+y**2)

      rho = (rhoprime)*exp(-log(2.0_prec)*r**2/Lr**2)

      this%solution%interior(i,j,iEl,1) = rho
      this%solution%interior(i,j,iEl,2) = 0.0_prec
      this%solution%interior(i,j,iEl,3) = 0.0_prec
      this%solution%interior(i,j,iEl,4) = rho*this%c*this%c

    enddo

    call this%ReportMetrics()

  endsubroutine SphericalSoundWave_LinearEuler2D_t

endmodule self_LinearEuler2D_t
