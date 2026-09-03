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

module self_LinearEuler3D_t
!! This module defines a class that can be used to solve the Linear Euler
!! equations in 3-D. The Linear Euler Equations, here, are the Euler equations
!! linearized about a motionless background state. The sound speed, the
!! background density, and the relaxation rate are carried as solution
!! variables (static in time, possibly spatially varying) so that
!! heterogeneous media and sponge/damping layers are supported, mirroring the
!! 2-D model.
!!
!! The density anomaly is not carried as a solution variable: for a motionless
!! background state it is slaved to the pressure through \(\rho = p/c^2\) and
!! never feeds back into the velocity or pressure dynamics, so only the
!! velocity components and the pressure are forward-stepped. If the density
!! anomaly is needed as a diagnostic, it can be recovered pointwise as
!! \(\rho = p/c^2\).
!!
!! The solution variables are
!!
!! \begin{equation}
!! \vec{s} = \begin{pmatrix}
!!      u \\
!!      v \\
!!      w \\
!!      p \\
!!      c \\
!!      \rho_0 \\
!!      \sigma
!!  \end{pmatrix}
!! \end{equation}
!!
!! The conservative flux is
!!
!! \begin{equation}
!! \overleftrightarrow{f} = \begin{pmatrix}
!!      \frac{p}{\rho_0} \hat{x} \\
!!      \frac{p}{\rho_0} \hat{y} \\
!!      \frac{p}{\rho_0} \hat{z} \\
!!      c^2 \rho_0 ( u \hat{x} + v \hat{y} + w \hat{z}) \\
!!      0 \\
!!      0 \\
!!      0
!!  \end{pmatrix}
!! \end{equation}
!!
!! and the source term is the linear relaxation ("sponge"/damping) term
!!
!! \begin{equation}
!! \vec{q} = -\sigma ( u, v, w, p, 0, 0, 0 )^T
!! \end{equation}
!!
!! where \(\sigma \ge 0\) [s\(^{-1}\)] is the inverse of a local damping time
!! scale. With \(\sigma = 0\) everywhere (the default, since the solution array
!! is zero-initialized) the source vanishes identically and the model reduces
!! to the undamped linear Euler system. The sound speed, background density
!! and relaxation rate variables have zero flux and zero source; they are set
!! by the initial condition and preserved in time. This is entropy-stable for
!! piecewise-constant material regions aligned with element boundaries.
!! The Riemann flux is the impedance-matched (characteristic) flux, identical
!! in form to the 2-D model's; see riemannflux3d_LinearEuler3D_t.
!!

  use self_model
  use self_dgmodel3D
  use self_mesh
  use SELF_BoundaryConditions

  implicit none

  type,extends(dgmodel3D) :: LinearEuler3D_t
    ! Add any additional attributes here that are specific to your model
    real(prec) :: rho0 = 1.0_prec ! Reference density (used to fill variable 6 in initial conditions)
    real(prec) :: c = 1.0_prec ! Reference sound speed (used to fill variable 5 in initial conditions)
    real(prec) :: g = 0.0_prec ! gravitational acceleration (y-direction only)

  contains
    procedure :: SourceMethod => sourcemethod_LinearEuler3D_t
    procedure :: SetNumberOfVariables => SetNumberOfVariables_LinearEuler3D_t
    procedure :: SetMetadata => SetMetadata_LinearEuler3D_t
    procedure :: AdditionalInit => AdditionalInit_LinearEuler3D_t
    procedure :: entropy_func => entropy_func_LinearEuler3D_t
    procedure :: flux3D => flux3D_LinearEuler3D_t
    procedure :: riemannflux3D => riemannflux3D_LinearEuler3D_t
    procedure :: SphericalSoundWave => SphericalSoundWave_LinearEuler3D_t

  endtype LinearEuler3D_t

contains

  subroutine SetNumberOfVariables_LinearEuler3D_t(this)
    implicit none
    class(LinearEuler3D_t),intent(inout) :: this

    this%nvar = 7
    ! Only the first four variables (u, v, w, P) are advanced in time. The
    ! fifth, sixth and seventh variables - the sound speed c, the background
    ! density rho0, and the relaxation rate sigma - are spatially-varying but
    ! time-constant background fields: their flux and source are identically
    ! zero, so they are excluded from time integration rather than stepped to
    ! a no-op.
    this%nstepped = 4

  endsubroutine SetNumberOfVariables_LinearEuler3D_t

  subroutine SetMetadata_LinearEuler3D_t(this)
    implicit none
    class(LinearEuler3D_t),intent(inout) :: this

    call this%solution%SetName(1,"u") ! x-velocity component
    call this%solution%SetUnits(1,"m⋅s⁻¹")

    call this%solution%SetName(2,"v") ! y-velocity component
    call this%solution%SetUnits(2,"m⋅s⁻¹")

    call this%solution%SetName(3,"w") ! z-velocity component
    call this%solution%SetUnits(3,"m⋅s⁻¹")

    call this%solution%SetName(4,"P") ! Pressure
    call this%solution%SetUnits(4,"kg⋅m⁻¹⋅s⁻²")

    call this%solution%SetName(5,"c") ! Sound speed (static; possibly heterogeneous)
    call this%solution%SetUnits(5,"m⋅s⁻¹")

    call this%solution%SetName(6,"rho0") ! Background density (static; possibly heterogeneous)
    call this%solution%SetUnits(6,"kg⋅m⁻³")

    call this%solution%SetName(7,"sigma") ! Relaxation rate (inverse damping time scale)
    call this%solution%SetUnits(7,"s⁻¹")

  endsubroutine SetMetadata_LinearEuler3D_t

  subroutine AdditionalInit_LinearEuler3D_t(this)
    !! Register the (CPU) radiation boundary condition. GPU builds overwrite
    !! this registration with the device kernel in AdditionalInit_LinearEuler3D.
    implicit none
    class(LinearEuler3D_t),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    bcfunc => hbc3d_Radiation_LinearEuler3D
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_RADIATION,"radiation",bcfunc)

  endsubroutine AdditionalInit_LinearEuler3D_t

  subroutine hbc3d_Radiation_LinearEuler3D(bc,mymodel)
    !! Radiation BC: zero acoustic perturbation in the exterior state; the
    !! sound speed (variable 5), background density (variable 6) and relaxation
    !! rate (variable 7) are copied from the interior side so the Riemann solver
    !! sees a consistent c and rho0 and the exterior state stays a faithful
    !! mirror of the interior one for diagnostics. The relaxation rate carries
    !! no flux, so its exterior value never enters the tendency.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel
    ! Local
    integer :: n,i,j,iEl,s

    select type(m => mymodel)
    class is(LinearEuler3D_t)
      do n = 1,bc%nBoundaries
        iEl = bc%elements(n)
        s = bc%sides(n)
        do j = 1,m%solution%interp%N+1
          do i = 1,m%solution%interp%N+1
            m%solution%extBoundary(i,j,s,iEl,1:4) = 0.0_prec
            m%solution%extBoundary(i,j,s,iEl,5) = m%solution%boundary(i,j,s,iEl,5) ! c preserved
            m%solution%extBoundary(i,j,s,iEl,6) = m%solution%boundary(i,j,s,iEl,6) ! rho0 preserved
            m%solution%extBoundary(i,j,s,iEl,7) = m%solution%boundary(i,j,s,iEl,7) ! sigma preserved
          enddo
        enddo
      enddo
    endselect

  endsubroutine hbc3d_Radiation_LinearEuler3D

  pure function entropy_func_LinearEuler3D_t(this,s) result(e)
    !! The entropy function is the sum of kinetic and internal energy
    !! For the linear model, this is
    !!
    !! \begin{equation}
    !!   e = \frac{1}{2} \left( \rho_0*( u^2 + v^2 + w^2 ) + \frac{P^2}{\rho_0 c^2} \right)
    !!
    !! where the sound speed c is taken from s(5) and the background density
    !! rho0 from s(6).
    class(LinearEuler3D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: e

    e = 0.5_prec*s(6)*(s(1)*s(1)+s(2)*s(2)+s(3)*s(3))+ &
        0.5_prec*(s(4)*s(4)/(s(6)*s(5)*s(5)))

  endfunction entropy_func_LinearEuler3D_t

  subroutine sourcemethod_LinearEuler3D_t(this)
    !! Spatially-dependent linear relaxation (sponge / damping layer) source.
    !!
    !! The relaxation rate sigma = s(7), in units of inverse seconds, is the
    !! reciprocal of a local damping time scale. It contributes
    !!
    !!   q_i = -sigma * s_i ,   i = 1,2,3,4   (u, v, w, P)
    !!
    !! to the source term of the conservation law s_t + div(f) = q, relaxing
    !! the acoustic state toward the (motionless, unperturbed) background at
    !! rate sigma. Because a single rate is applied to the velocity and the
    !! pressure alike, the acoustic energy density decays as exp(-2*sigma*t)
    !! wherever sigma is constant, which is exactly the behaviour wanted from
    !! a sponge layer: no impedance mismatch is introduced between the damped
    !! and undamped regions by the source term itself.
    !!
    !! The static background variables (c, rho0, sigma) carry no source, so
    !! they are held fixed in time.
    !!
    !! Units: sigma [s^-1]; q_1, q_2, q_3 [m⋅s^-2]; q_4 [kg⋅m^-1⋅s^-3].
    !! Expected input range: sigma >= 0. A negative sigma amplifies the
    !! solution and is not a supported configuration.
    !!
    !! Note: sourcemethod is overridden rather than source3d so that the
    !! solution gradient - which this model does not use - is never read.
    implicit none
    class(LinearEuler3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,k,iel
    real(prec) :: sigma

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%solution%N+1,iel=1:this%mesh%nElem)

      sigma = this%solution%interior(i,j,k,iel,7)

      this%source%interior(i,j,k,iel,1) = -this%solution%interior(i,j,k,iel,1)*sigma
      this%source%interior(i,j,k,iel,2) = -this%solution%interior(i,j,k,iel,2)*sigma
      this%source%interior(i,j,k,iel,3) = -this%solution%interior(i,j,k,iel,3)*sigma
      this%source%interior(i,j,k,iel,4) = -this%solution%interior(i,j,k,iel,4)*sigma
      this%source%interior(i,j,k,iel,5) = 0.0_prec ! sound speed; held fixed in time
      this%source%interior(i,j,k,iel,6) = 0.0_prec ! background density; held fixed in time
      this%source%interior(i,j,k,iel,7) = 0.0_prec ! relaxation rate; held fixed in time

    enddo

  endsubroutine sourcemethod_LinearEuler3D_t

  pure function flux3D_LinearEuler3D_t(this,s,dsdx) result(flux)
    class(LinearEuler3D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec) :: flux(1:this%nvar,1:3)

    flux(1,1) = s(4)/s(6) ! x-velocity, x flux; p/rho0
    flux(1,2) = 0.0_prec ! x-velocity, y flux; 0
    flux(1,3) = 0.0_prec ! x-velocity, z flux; 0

    flux(2,1) = 0.0_prec ! y-velocity, x flux; 0
    flux(2,2) = s(4)/s(6) ! y-velocity, y flux; p/rho0
    flux(2,3) = 0.0_prec ! y-velocity, z flux; 0

    flux(3,1) = 0.0_prec ! z-velocity, x flux; 0
    flux(3,2) = 0.0_prec ! z-velocity, y flux; 0
    flux(3,3) = s(4)/s(6) ! z-velocity, z flux; p/rho0

    flux(4,1) = s(5)*s(5)*s(6)*s(1) ! pressure, x flux : rho0*c^2*u
    flux(4,2) = s(5)*s(5)*s(6)*s(2) ! pressure, y flux : rho0*c^2*v
    flux(4,3) = s(5)*s(5)*s(6)*s(3) ! pressure, z flux : rho0*c^2*w

    flux(5,1) = 0.0_prec ! sound speed, x flux; 0 (c held fixed in time)
    flux(5,2) = 0.0_prec ! sound speed, y flux; 0
    flux(5,3) = 0.0_prec ! sound speed, z flux; 0

    flux(6,1) = 0.0_prec ! background density, x flux; 0 (rho0 held fixed in time)
    flux(6,2) = 0.0_prec ! background density, y flux; 0
    flux(6,3) = 0.0_prec ! background density, z flux; 0

    flux(7,1) = 0.0_prec ! relaxation rate, x flux; 0 (sigma held fixed in time)
    flux(7,2) = 0.0_prec ! relaxation rate, y flux; 0
    flux(7,3) = 0.0_prec ! relaxation rate, z flux; 0
    if(.false.) flux(1,1) = flux(1,1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction flux3D_LinearEuler3D_t

  pure function riemannflux3D_LinearEuler3D_t(this,sL,sR,dsdx,nhat) result(flux)
    !! Impedance-matched (characteristic/Godunov) Riemann flux, identical in
    !! form to the 2-D model's. The interface states are resolved with the
    !! per-side acoustic impedances Z = rho0*c (each side using its own
    !! background density rho0 and sound speed c), so material interfaces in a
    !! heterogeneous field are handled with the physically correct
    !! transmission/reflection:
    !!
    !!   un* = (ZL*unL + ZR*unR + (pL - pR)) / (ZL + ZR)
    !!   p*  = (ZR*pL + ZL*pR + ZL*ZR*(unL - unR)) / (ZL + ZR)
    !!
    !! The reconstructed momentum/pressure fluxes use the arithmetic
    !! averages rho0_avg and c2_avg at the face (see the 2-D model for the
    !! rationale). The sound speed, background density and relaxation rate
    !! variables carry zero flux.
    class(LinearEuler3D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec),intent(in) :: nhat(1:3)
    real(prec) :: flux(1:this%nvar)
    ! Local
    real(prec) :: rho0L,rho0R,rho0_avg,cL,cR,ZL,ZR,unL,unR,pL,pR,un_star,p_star,c2_avg

    rho0L = sL(6)
    rho0R = sR(6)
    rho0_avg = 0.5_prec*(rho0L+rho0R)
    cL = sL(5)
    cR = sR(5)
    ZL = rho0L*cL
    ZR = rho0R*cR

    unL = sL(1)*nhat(1)+sL(2)*nhat(2)+sL(3)*nhat(3)
    unR = sR(1)*nhat(1)+sR(2)*nhat(2)+sR(3)*nhat(3)
    pL = sL(4)
    pR = sR(4)

    un_star = (ZL*unL+ZR*unR+(pL-pR))/(ZL+ZR)
    p_star = (ZR*pL+ZL*pR+ZL*ZR*(unL-unR))/(ZL+ZR)
    c2_avg = 0.5_prec*(cL*cL+cR*cR)

    flux(1) = p_star*nhat(1)/rho0_avg
    flux(2) = p_star*nhat(2)/rho0_avg
    flux(3) = p_star*nhat(3)/rho0_avg
    flux(4) = rho0_avg*c2_avg*un_star
    flux(5) = 0.0_prec
    flux(6) = 0.0_prec
    flux(7) = 0.0_prec
    if(.false.) flux(1) = flux(1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction riemannflux3D_LinearEuler3D_t

  subroutine SphericalSoundWave_LinearEuler3D_t(this,rhoprime,Lr,x0,y0,z0)
    !! This subroutine sets the initial condition for a weak blast wave
    !! problem. The initial condition is given by
    !!
    !! \begin{equation}
    !! \begin{aligned}
    !! u &= 0 \\
    !! v &= 0 \\
    !! w &= 0 \\
    !! p &= \rho' c^2 \exp\left( -\ln(2) \frac{(x-x_0)^2 + (y-y_0)^2 + (z-z_0)^2}{L_r^2} \right)
    !! \end{aligned}
    !! \end{equation}
    !!
    !! `rhoprime` is the amplitude of the (diagnostic) density anomaly that the
    !! pressure pulse corresponds to through the acoustic relation p = rho*c^2;
    !! the density anomaly itself is not a solution variable. The relaxation
    !! rate (variable 7) is left untouched, so a sponge profile may be
    !! installed either before or after this initializer is called; it is zero
    !! by default.
    !!
    implicit none
    class(LinearEuler3D_t),intent(inout) :: this
    real(prec),intent(in) ::  rhoprime,Lr,x0,y0,z0
    ! Local
    integer :: i,j,k,iEl
    real(prec) :: x,y,z,rho,r

    print*,__FILE__," : Configuring weak blast wave initial condition. "
    print*,__FILE__," : rhoprime = ",rhoprime
    print*,__FILE__," : Lr = ",Lr
    print*,__FILE__," : x0 = ",x0
    print*,__FILE__," : y0 = ",y0
    print*,__FILE__," : z0 = ",z0

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%solution%N+1,iel=1:this%mesh%nElem)
      x = this%geometry%x%interior(i,j,k,iEl,1,1)-x0
      y = this%geometry%x%interior(i,j,k,iEl,1,2)-y0
      z = this%geometry%x%interior(i,j,k,iEl,1,3)-z0
      r = sqrt(x**2+y**2+z**2)

      rho = (rhoprime)*exp(-log(2.0_prec)*r**2/Lr**2)

      this%solution%interior(i,j,k,iEl,1) = 0.0_prec
      this%solution%interior(i,j,k,iEl,2) = 0.0_prec
      this%solution%interior(i,j,k,iEl,3) = 0.0_prec
      this%solution%interior(i,j,k,iEl,4) = rho*this%c*this%c
      this%solution%interior(i,j,k,iEl,5) = this%c ! uniform background sound speed
      this%solution%interior(i,j,k,iEl,6) = this%rho0 ! uniform background density

    enddo

    call this%ReportMetrics()
    call this%solution%UpdateDevice()

  endsubroutine SphericalSoundWave_LinearEuler3D_t

endmodule self_LinearEuler3D_t
