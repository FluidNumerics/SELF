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
!! The solution variables are
!!
!! \begin{equation}
!! \vec{s} = \begin{pmatrix}
!!     \rho \\
!!      u \\
!!      v \\
!!      p \\
!!      c
!!  \end{pmatrix}
!! \end{equation}
!!
!! The sound speed \(c\) is carried as a solution variable so that it
!! can vary in space. Its flux and source are identically zero, so the
!! sound speed is held fixed in time at each spatial location.
!!
!! The conservative flux is
!!
!! \begin{equation}
!! \overleftrightarrow{f} = \begin{pmatrix}
!!     \rho_0 u \hat{x} + \rho_0 v \hat{y} \\
!!      \frac{p}{\rho_0} \hat{x} \\
!!      \frac{p}{\rho_0} \hat{y} \\
!!      c^2 \rho_0 ( u \hat{x} + v \hat{y} ) \\
!!      \vec{0}
!!  \end{pmatrix}
!! \end{equation}
!!
!! and the source terms are null.
!!

  use self_model
  use self_dgmodel2d
  use self_mesh
  use SELF_BoundaryConditions

  implicit none

  type,extends(dgmodel2d) :: LinearEuler2D_t
    ! Add any additional attributes here that are specific to your model
    real(prec) :: rho0 = 1.0_prec ! Reference density
    real(prec) :: g = 0.0_prec ! gravitational acceleration (y-direction only)

  contains
    procedure :: SetNumberOfVariables => SetNumberOfVariables_LinearEuler2D_t
    procedure :: SetMetadata => SetMetadata_LinearEuler2D_t
    procedure :: AdditionalInit => AdditionalInit_LinearEuler2D_t
    procedure :: entropy_func => entropy_func_LinearEuler2D_t
    procedure :: flux2d => flux2d_LinearEuler2D_t
    procedure :: riemannflux2d => riemannflux2d_LinearEuler2D_t
    !procedure :: source2d => source2d_LinearEuler2D_t
    procedure :: SphericalSoundWave => SphericalSoundWave_LinearEuler2D_t

  endtype LinearEuler2D_t

contains

  subroutine SetNumberOfVariables_LinearEuler2D_t(this)
    implicit none
    class(LinearEuler2D_t),intent(inout) :: this

    this%nvar = 5

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

    call this%solution%SetName(5,"c") ! Sound speed
    call this%solution%SetUnits(5,"m⋅s⁻¹")

  endsubroutine SetMetadata_LinearEuler2D_t

  pure function entropy_func_LinearEuler2D_t(this,s) result(e)
    !! The entropy function is the sum of kinetic and internal energy
    !! For the linear model, this is
    !!
    !! \begin{equation}
    !!   e = \frac{1}{2} \left( \rho_0*( u^2 + v^2 ) + \frac{P^2}{\rho_0 c^2} \right)
    !!
    !! where the sound speed c is taken from s(5).
    class(LinearEuler2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: e

    e = 0.5_prec*this%rho0*(s(2)*s(2)+s(3)*(3))+ &
        0.5_prec*(s(4)*s(4)/(this%rho0*s(5)*s(5)))

  endfunction entropy_func_LinearEuler2D_t

  subroutine AdditionalInit_LinearEuler2D_t(this)
    implicit none
    class(LinearEuler2D_t),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    bcfunc => hbc2d_NoNormalFlow_LinearEuler2D
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

  endsubroutine AdditionalInit_LinearEuler2D_t

  subroutine hbc2d_NoNormalFlow_LinearEuler2D(bc,mymodel)
    !! No-normal-flow boundary condition for 2D linear Euler equations.
    !! Reflects the velocity vector about the boundary normal while
    !! preserving density, pressure, and sound speed.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel
    ! Local
    integer :: n,i,iEl,j
    real(prec) :: nhat(1:2),s(1:5)

    select type(m => mymodel)
    class is(LinearEuler2D_t)
      do n = 1,bc%nBoundaries
        iEl = bc%elements(n)
        j = bc%sides(n)
        do i = 1,m%solution%interp%N+1
          nhat = m%geometry%nhat%boundary(i,j,iEl,1,1:2)
          s = m%solution%boundary(i,j,iEl,1:5)
          m%solution%extBoundary(i,j,iEl,1) = s(1) ! density
          m%solution%extBoundary(i,j,iEl,2) = &
            (nhat(2)**2-nhat(1)**2)*s(2)-2.0_prec*nhat(1)*nhat(2)*s(3) ! u
          m%solution%extBoundary(i,j,iEl,3) = &
            (nhat(1)**2-nhat(2)**2)*s(3)-2.0_prec*nhat(1)*nhat(2)*s(2) ! v
          m%solution%extBoundary(i,j,iEl,4) = s(4) ! p
          m%solution%extBoundary(i,j,iEl,5) = s(5) ! c
        enddo
      enddo
    endselect

  endsubroutine hbc2d_NoNormalFlow_LinearEuler2D

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
    flux(4,1) = s(5)*s(5)*this%rho0*s(2) ! pressure, x flux : rho0*c^2*u
    flux(4,2) = s(5)*s(5)*this%rho0*s(3) ! pressure, y flux : rho0*c^2*v
    flux(5,1) = 0.0_prec ! sound speed, x flux; 0 (c is held fixed in time)
    flux(5,2) = 0.0_prec ! sound speed, y flux; 0 (c is held fixed in time)
    if(.false.) flux(1,1) = flux(1,1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction flux2d_LinearEuler2D_t

  pure function riemannflux2d_LinearEuler2D_t(this,sL,sR,dsdx,nhat) result(flux)
    !! Characteristic-decomposition (impedance-matched) interface flux for
    !! linear acoustics with possibly discontinuous sound speed.
    !!
    !! The normal-flux Jacobian has eigenstructure
    !!   +c : right-going acoustic mode, W_+ = rho0*c*u_n + p
    !!   -c : left-going  acoustic mode, W_- = -rho0*c*u_n + p
    !!    0 : entropy density mode      , W_0 = rho' - p/c^2
    !!    0 : tangential vorticity mode , u_t
    !! Upwinding each mode by its characteristic direction at the face,
    !!   W_+|* = W_+|_L   (transmitted from left,  at speed +c_L)
    !!   W_-|* = W_-|_R   (transmitted from right, at speed -c_R)
    !! yields the impedance-matched interface state
    !!   u_n* = (Z_L u_n,L + Z_R u_n,R + (p_L - p_R)) / (Z_L + Z_R)
    !!   p*   = (Z_R p_L + Z_L p_R + Z_L Z_R (u_n,L - u_n,R)) / (Z_L + Z_R)
    !! with Z = rho0*c. This is exact upwind / Godunov for the linearised
    !! acoustic system and reduces correctly to Fresnel reflection /
    !! transmission across an impedance jump (c_L .ne. c_R). LLF with
    !! cmax = max(c_L, c_R) over-dissipates tangential and entropy modes
    !! and at high polynomial order fails to stably handle the
    !! impedance mismatch (aliasing instability at material interfaces).
    !!
    !! The pressure flux uses an averaged c^2; this is a pragmatic
    !! treatment of the non-conservative product rho0*c^2*div(v) at a
    !! face where c jumps. A fully path-conservative (Castro-Pares)
    !! treatment would use each side's own c^2 in its surface integral.
    class(LinearEuler2D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: flux(1:this%nvar)
    ! Local
    real(prec) :: rho0,cL,cR,ZL,ZR,unL,unR,pL,pR,un_star,p_star,c2_avg

    rho0 = this%rho0
    cL = sL(5)
    cR = sR(5)
    ZL = rho0*cL
    ZR = rho0*cR

    unL = sL(2)*nhat(1)+sL(3)*nhat(2)
    unR = sR(2)*nhat(1)+sR(3)*nhat(2)
    pL = sL(4)
    pR = sR(4)

    un_star = (ZL*unL+ZR*unR+(pL-pR))/(ZL+ZR)
    p_star = (ZR*pL+ZL*pR+ZL*ZR*(unL-unR))/(ZL+ZR)
    c2_avg = 0.5_prec*(cL*cL+cR*cR)

    flux(1) = rho0*un_star
    flux(2) = p_star*nhat(1)/rho0
    flux(3) = p_star*nhat(2)/rho0
    flux(4) = rho0*c2_avg*un_star
    flux(5) = 0.0_prec
    if(.false.) flux(1) = flux(1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction riemannflux2d_LinearEuler2D_t

  subroutine SphericalSoundWave_LinearEuler2D_t(this,rhoprime,Lr,x0,y0,c)
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
    !! The sound speed `c` (passed as an argument since it is no longer a
    !! scalar model attribute) is set uniformly across the domain.
    implicit none
    class(LinearEuler2D_t),intent(inout) :: this
    real(prec),intent(in) ::  rhoprime,Lr,x0,y0,c
    ! Local
    integer :: i,j,iEl
    real(prec) :: x,y,rho,r

    print*,__FILE__," : Configuring weak blast wave initial condition. "
    print*,__FILE__," : rhoprime = ",rhoprime
    print*,__FILE__," : Lr = ",Lr
    print*,__FILE__," : x0 = ",x0
    print*,__FILE__," : y0 = ",y0
    print*,__FILE__," : c = ",c

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem)
      x = this%geometry%x%interior(i,j,iEl,1,1)-x0
      y = this%geometry%x%interior(i,j,iEl,1,2)-y0
      r = sqrt(x**2+y**2)

      rho = (rhoprime)*exp(-log(2.0_prec)*r**2/Lr**2)

      this%solution%interior(i,j,iEl,1) = rho
      this%solution%interior(i,j,iEl,2) = 0.0_prec
      this%solution%interior(i,j,iEl,3) = 0.0_prec
      this%solution%interior(i,j,iEl,4) = rho*c*c
      this%solution%interior(i,j,iEl,5) = c

    enddo

    call this%ReportMetrics()
    call this%solution%UpdateDevice()

  endsubroutine SphericalSoundWave_LinearEuler2D_t

endmodule self_LinearEuler2D_t
