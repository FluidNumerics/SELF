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

module SELF_ESAtmo3D_t
  !! Entropy-conserving compressible Euler equations in 3-D with potential
  !! temperature formulation (FastEddy / atmospheric LES equation set).
  !!
  !! Conservative variables:
  !!
  !!   U = [rho, rho*u, rho*v, rho*w, rho*theta]
  !!
  !! where theta is potential temperature.
  !!
  !! Governing equations:
  !!
  !!   d(rho)/dt     + div(rho * v)         = 0
  !!   d(rho*u)/dt   + div(rho*u * v) + dp/dx = 0
  !!   d(rho*v)/dt   + div(rho*v * v) + dp/dy = 0
  !!   d(rho*w)/dt   + div(rho*w * v) + dp/dz = -rho * g
  !!   d(rho*th)/dt  + div(rho*th * v)      = 0
  !!
  !! Equation of state (dry ideal gas):
  !!
  !!   p = p0 * (rho * Rd * theta / p0)^(cp/cv)
  !!
  !! Two-point volume flux: Kennedy-Gruber (kinetic energy preserving)
  !!
  !! Surface Riemann flux: Local Lax-Friedrichs (Rusanov)

  use SELF_ECDGModel3D
  use SELF_MappedScalar_3D
  use SELF_mesh
  use SELF_BoundaryConditions

  implicit none

  type,extends(ECDGModel3D) :: ESAtmo3D_t

    real(prec) :: p0 = 100000.0_prec ! Reference pressure [Pa]
    real(prec) :: Rd = 287.0_prec ! Gas constant for dry air [J/(kg*K)]
    real(prec) :: cp = 1004.0_prec ! Specific heat at constant pressure [J/(kg*K)]
    real(prec) :: cv = 717.0_prec ! Specific heat at constant volume [J/(kg*K)]
    real(prec) :: g = 9.81_prec ! Gravitational acceleration [m/s^2]

    !! Constant-coefficient Laplacian diffusion. The interface flux is
    !! BR1 central + Nitsche-style jump penalty (i.e. SIPG-style):
    !!
    !!   f_R^diff(iVar) = -coeff(iVar) * (avg_grad . n) * nmag
    !!                    + tau(iVar) * (uL - uR) * nmag
    !!
    !! with tau(iVar) = eta_penalty * coeff(iVar) * (N+1)^2 / length_scale.
    !! The penalty is what makes the discrete operator coercive on the
    !! checkerboard / odd-even modes that pure BR1 cannot damp.
    !!
    !!   nu             : kinematic diffusivity for momentum [m^2/s]
    !!   kappa          : thermal diffusivity for rho*theta  [m^2/s]
    !!   eta_penalty    : dimensionless SIPG penalty (default 4.0)
    !!   length_scale   : characteristic element length [m], filled by
    !!                    SetDiffusion from the geometry (mean J^(1/3)*2).
    !!
    !! All default to zero (inviscid). Setting nu or kappa > 0 via
    !! SetDiffusion enables the gradient pipeline so that the model
    !! can read solutionGradient inside the diffusive flux methods.
    real(prec) :: nu = 0.0_prec
    real(prec) :: kappa = 0.0_prec
    real(prec) :: eta_penalty = 4.0_prec
    real(prec) :: length_scale = 0.0_prec

    !! Diffusive-flux scratch buffer. Holds the constant-coefficient
    !! Laplacian flux F_diff(i,j,k,iel,iVar,d) = -coeff_iVar * d(s_iVar)/dx_d
    !! with coeff_iVar = 0 for rho, nu for rhou/rhov/rhow, and kappa for
    !! rhotheta. The boundaryNormal slot is filled with the BR1 central
    !! flux F_R^diff = -coeff * (avg_grad . n) * nmag (using
    !! solutionGradient%avgBoundary, populated by AverageSides).
    !!
    !! MappedDGDivergence on this buffer gives the parabolic divergence
    !! contribution that is then accumulated into fluxDivergence.
    type(MappedVector3D) :: diffFlux
    type(MappedScalar3D) :: diffDiv

  contains

    procedure :: SetNumberOfVariables => SetNumberOfVariables_ESAtmo3D_t
    procedure :: SetMetadata => SetMetadata_ESAtmo3D_t
    procedure :: entropy_func => entropy_func_ESAtmo3D_t
    procedure :: twopointflux3d => twopointflux3d_ESAtmo3D_t
    procedure :: TwoPointFluxMethod => TwoPointFluxMethod_ESAtmo3D_t
    procedure :: riemannflux3d => riemannflux3d_ESAtmo3D_t
    procedure :: BoundaryFlux => BoundaryFlux_ESAtmo3D_t
    procedure :: SourceMethod => SourceMethod_ESAtmo3D_t
    procedure :: AdditionalInit => AdditionalInit_ESAtmo3D_t
    procedure :: AdditionalFree => AdditionalFree_ESAtmo3D_t
    procedure :: SetHydrostaticBalance => SetHydrostaticBalance_ESAtmo3D_t
    procedure :: AddThermalBubble => AddThermalBubble_ESAtmo3D_t

    !! Constant-coefficient Laplacian / Bassi-Rebay diffusion hooks
    procedure :: SetDiffusion => SetDiffusion_ESAtmo3D_t
    procedure :: DiffusiveFluxMethod => DiffusiveFluxMethod_ESAtmo3D_t
    procedure :: DiffusiveBoundaryFlux => DiffusiveBoundaryFlux_ESAtmo3D_t
    procedure :: CalculateTendency => CalculateTendency_ESAtmo3D_t

  endtype ESAtmo3D_t

contains

  subroutine SetNumberOfVariables_ESAtmo3D_t(this)
    !! Six conserved variables: (rho, rho*u, rho*v, rho*w, rho*theta, Phi),
    !! where Phi = g*z is the geopotential. Phi has zero flux (volume and
    !! surface) so its tendency is identically zero; it is carried in the
    !! state vector solely so that the Souza et al. (2023) non-conservative
    !! gravity flux differencing in SourceMethod can read it node-locally.
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this

    this%nvar = 6

  endsubroutine SetNumberOfVariables_ESAtmo3D_t

  subroutine SetMetadata_ESAtmo3D_t(this)
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this

    call this%solution%SetName(1,"rho")
    call this%solution%SetUnits(1,"kg/m^3")

    call this%solution%SetName(2,"rhou")
    call this%solution%SetUnits(2,"kg/(m^2 s)")

    call this%solution%SetName(3,"rhov")
    call this%solution%SetUnits(3,"kg/(m^2 s)")

    call this%solution%SetName(4,"rhow")
    call this%solution%SetUnits(4,"kg/(m^2 s)")

    call this%solution%SetName(5,"rhotheta")
    call this%solution%SetUnits(5,"kg K/m^3")

    call this%solution%SetName(6,"phi")
    call this%solution%SetUnits(6,"m^2/s^2")

  endsubroutine SetMetadata_ESAtmo3D_t

  pure function entropy_func_ESAtmo3D_t(this,s) result(e)
    !! Mathematical entropy: total energy density (kinetic + internal).
    !!
    !!   e = 0.5*(rhou^2 + rhov^2 + rhow^2)/rho + p/(gamma - 1)
    !!
    !! where p = p0 * (rho * Rd * theta / p0)^gamma and gamma = cp/cv.
    class(ESAtmo3D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: e
    ! Local
    real(prec) :: rho,theta,p,gamma

    rho = s(1)
    theta = s(5)/rho
    gamma = this%cp/this%cv
    p = this%p0*(rho*this%Rd*theta/this%p0)**gamma

    e = 0.5_prec*(s(2)*s(2)+s(3)*s(3)+s(4)*s(4))/rho+ &
        p/(gamma-1.0_prec)

  endfunction entropy_func_ESAtmo3D_t

  pure function log_mean(a,b) result(am)
    !! Numerically stable logarithmic mean (Ismail-Roe 2009 / Ranocha 2018).
    !!
    !!   <a>_log = (a - b) / (ln(a) - ln(b))   for a /= b
    !!   <a>_log = (a + b) / 2                 for a == b
    !!
    !! Falls back to a Taylor series in u = f^2 with f = (a-b)/(a+b)
    !! when a and b are close, to avoid 0/0.
    real(prec),intent(in) :: a,b
    real(prec) :: am
    ! Local
    real(prec) :: zeta,f,u,F_F
    real(prec),parameter :: eps_lm = 1.0e-2_prec

    zeta = a/b
    f = (zeta-1.0_prec)/(zeta+1.0_prec)
    u = f*f
    if(u < eps_lm) then
      F_F = 1.0_prec+u/3.0_prec+u*u/5.0_prec+u*u*u/7.0_prec
    else
      F_F = log(zeta)/(2.0_prec*f)
    endif
    am = (a+b)/(2.0_prec*F_F)
  endfunction log_mean

  pure function twopointflux3d_ESAtmo3D_t(this,sL,sR) result(flux)
    !! Souza et al. (2023, JAMES) entropy-conservative two-point flux for
    !! compressible Euler in (rho, rho*v, rho*theta) variables with
    !! p = p0*(rho*Rd*theta/p0)^gamma.
    !!
    !!   f_d(rho)     = <rho>_log * <v_d>
    !!   f_d(rho*v_i) = <rho>_log * <v_i> * <v_d> + <p> * delta_{id}
    !!   f_d(rho*th)  = <rho*theta>_log * <v_d>
    !!   f_d(Phi)     = 0  (geopotential is carried, not advected)
    !!
    !! <a>_log is the logarithmic mean (see log_mean), <a> the arithmetic
    !! mean. Pressure here is the *total* pressure; gravity is handled by
    !! the Souza non-conservative term in SourceMethod, not by a flux split.
    class(ESAtmo3D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec) :: flux(1:this%nvar,1:3)
    ! Local
    real(prec) :: rhoL,uL,vL,wL,thetaL,pL,rthL
    real(prec) :: rhoR,uR,vR,wR,thetaR,pR,rthR
    real(prec) :: rho_log,rth_log,u_avg,v_avg,w_avg,p_avg
    real(prec) :: gamma

    gamma = this%cp/this%cv

    ! Left primitive variables
    rhoL = sL(1)
    rthL = sL(5)
    thetaL = rthL/rhoL
    uL = sL(2)/rhoL
    vL = sL(3)/rhoL
    wL = sL(4)/rhoL
    pL = this%p0*(rhoL*this%Rd*thetaL/this%p0)**gamma

    ! Right primitive variables
    rhoR = sR(1)
    rthR = sR(5)
    thetaR = rthR/rhoR
    uR = sR(2)/rhoR
    vR = sR(3)/rhoR
    wR = sR(4)/rhoR
    pR = this%p0*(rhoR*this%Rd*thetaR/this%p0)**gamma

    ! Logarithmic and arithmetic means
    rho_log = log_mean(rhoL,rhoR)
    rth_log = log_mean(rthL,rthR)
    u_avg = 0.5_prec*(uL+uR)
    v_avg = 0.5_prec*(vL+vR)
    w_avg = 0.5_prec*(wL+wR)
    p_avg = 0.5_prec*(pL+pR)

    ! x-direction flux (d=1)
    flux(1,1) = rho_log*u_avg
    flux(2,1) = rho_log*u_avg*u_avg+p_avg
    flux(3,1) = rho_log*v_avg*u_avg
    flux(4,1) = rho_log*w_avg*u_avg
    flux(5,1) = rth_log*u_avg
    flux(6,1) = 0.0_prec

    ! y-direction flux (d=2)
    flux(1,2) = rho_log*v_avg
    flux(2,2) = rho_log*u_avg*v_avg
    flux(3,2) = rho_log*v_avg*v_avg+p_avg
    flux(4,2) = rho_log*w_avg*v_avg
    flux(5,2) = rth_log*v_avg
    flux(6,2) = 0.0_prec

    ! z-direction flux (d=3)
    flux(1,3) = rho_log*w_avg
    flux(2,3) = rho_log*u_avg*w_avg
    flux(3,3) = rho_log*v_avg*w_avg
    flux(4,3) = rho_log*w_avg*w_avg+p_avg
    flux(5,3) = rth_log*w_avg
    flux(6,3) = 0.0_prec

  endfunction twopointflux3d_ESAtmo3D_t

  subroutine TwoPointFluxMethod_ESAtmo3D_t(this)
    !! Pre-projected scalar contravariant two-point Souza et al. (2023) EC
    !! flux. For each node pair (a, b) along reference direction r:
    !!
    !!   Fc^r_v_i = sum_d 0.5*(Ja^r_d(a) + Ja^r_d(b)) * f_d(v_i, sL=s(a), sR=s(b))
    !!
    !! Gravity is NOT split into the pressure flux here; it is handled by
    !! the Souza non-conservative gravity flux differencing in SourceMethod
    !! (which uses the geopotential carried as state variable index 6).
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this
    ! Local
    integer :: nn,i,j,k,d,iEl,iVar
    real(prec) :: sL(1:this%nvar),sR(1:this%nvar)
    real(prec) :: Fphys(1:this%nvar,1:3)
    real(prec) :: Fc

    do concurrent(nn=1:this%solution%N+1,i=1:this%solution%N+1, &
                  j=1:this%solution%N+1,k=1:this%solution%N+1, &
                  iEl=1:this%mesh%nElem)

      sL = this%solution%interior(i,j,k,iEl,1:this%nvar)

      ! ------------- xi^1: pair (i,j,k)-(nn,j,k) -------------
      sR = this%solution%interior(nn,j,k,iEl,1:this%nvar)
      Fphys = this%twopointflux3d(sL,sR)
      do iVar = 1,this%nvar
        Fc = 0.0_prec
        do d = 1,3
          Fc = Fc+0.5_prec*( &
               this%geometry%dsdx%interior(i,j,k,iEl,1,d,1)+ &
               this%geometry%dsdx%interior(nn,j,k,iEl,1,d,1))* &
               Fphys(iVar,d)
        enddo
        this%twoPointFlux%interior(nn,i,j,k,iEl,iVar,1) = Fc
      enddo

      ! ------------- xi^2: pair (i,j,k)-(i,nn,k) -------------
      sR = this%solution%interior(i,nn,k,iEl,1:this%nvar)
      Fphys = this%twopointflux3d(sL,sR)
      do iVar = 1,this%nvar
        Fc = 0.0_prec
        do d = 1,3
          Fc = Fc+0.5_prec*( &
               this%geometry%dsdx%interior(i,j,k,iEl,1,d,2)+ &
               this%geometry%dsdx%interior(i,nn,k,iEl,1,d,2))* &
               Fphys(iVar,d)
        enddo
        this%twoPointFlux%interior(nn,i,j,k,iEl,iVar,2) = Fc
      enddo

      ! ------------- xi^3: pair (i,j,k)-(i,j,nn) -------------
      sR = this%solution%interior(i,j,nn,iEl,1:this%nvar)
      Fphys = this%twopointflux3d(sL,sR)
      do iVar = 1,this%nvar
        Fc = 0.0_prec
        do d = 1,3
          Fc = Fc+0.5_prec*( &
               this%geometry%dsdx%interior(i,j,k,iEl,1,d,3)+ &
               this%geometry%dsdx%interior(i,j,nn,iEl,1,d,3))* &
               Fphys(iVar,d)
        enddo
        this%twoPointFlux%interior(nn,i,j,k,iEl,iVar,3) = Fc
      enddo

    enddo

  endsubroutine TwoPointFluxMethod_ESAtmo3D_t

  pure function riemannflux3d_ESAtmo3D_t(this,sL,sR,dsdx,nhat) result(flux)
    !! Local Lax-Friedrichs (Rusanov) Riemann flux.
    !!
    !!   F* = 0.5*(fL.n + fR.n) - 0.5*lambda_max*(sR - sL)
    !!
    !! where lambda_max = max(|vL.n| + cL, |vR.n| + cR)
    !! and c = sqrt(gamma * p / rho) is the sound speed.
    class(ESAtmo3D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec),intent(in) :: nhat(1:3)
    real(prec) :: flux(1:this%nvar)
    ! Local
    real(prec) :: rhoL,uL,vL,wL,thetaL,pL,unL,cL
    real(prec) :: rhoR,uR,vR,wR,thetaR,pR,unR,cR
    real(prec) :: fL(1:5),fR(1:5)
    real(prec) :: gamma,lam

    gamma = this%cp/this%cv

    ! Left state
    rhoL = sL(1)
    uL = sL(2)/rhoL
    vL = sL(3)/rhoL
    wL = sL(4)/rhoL
    thetaL = sL(5)/rhoL
    pL = this%p0*(rhoL*this%Rd*thetaL/this%p0)**gamma
    unL = uL*nhat(1)+vL*nhat(2)+wL*nhat(3)
    cL = sqrt(gamma*pL/rhoL)

    ! Right state
    rhoR = sR(1)
    uR = sR(2)/rhoR
    vR = sR(3)/rhoR
    wR = sR(4)/rhoR
    thetaR = sR(5)/rhoR
    pR = this%p0*(rhoR*this%Rd*thetaR/this%p0)**gamma
    unR = uR*nhat(1)+vR*nhat(2)+wR*nhat(3)
    cR = sqrt(gamma*pR/rhoR)

    ! Normal flux from left state
    fL(1) = rhoL*unL
    fL(2) = sL(2)*unL+pL*nhat(1)
    fL(3) = sL(3)*unL+pL*nhat(2)
    fL(4) = sL(4)*unL+pL*nhat(3)
    fL(5) = sL(5)*unL

    ! Normal flux from right state
    fR(1) = rhoR*unR
    fR(2) = sR(2)*unR+pR*nhat(1)
    fR(3) = sR(3)*unR+pR*nhat(2)
    fR(4) = sR(4)*unR+pR*nhat(3)
    fR(5) = sR(5)*unR

    ! Maximum wave speed
    lam = max(abs(unL)+cL,abs(unR)+cR)

    ! LLF flux
    flux(1:5) = 0.5_prec*(fL(1:5)+fR(1:5))- &
                0.5_prec*lam*(sR(1:5)-sL(1:5))
    ! Geopotential is carried in the state vector but has no flux.
    flux(6) = 0.0_prec

    if(.false.) flux(1) = flux(1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction riemannflux3d_ESAtmo3D_t

  subroutine BoundaryFlux_ESAtmo3D_t(this)
    !! LMARS (Low-Mach Approximate Riemann Solver, Chen et al. 2013)
    !! interface flux. No hydrostatic pressure split: gravity is folded
    !! into SourceMethod via the Souza non-conservative form using the
    !! geopotential carried in the state vector (variable index 6).
    !!
    !!   un* = 0.5*(unL + unR) - (pR - pL) / (2 * rho_bar * c_bar)
    !!   p*  = 0.5*(pL + pR)   - 0.5 * rho_bar * c_bar * (unR - unL)
    !!
    !! Upwind on sign of un* selects which side supplies the conserved
    !! state for the advective part. Pressure adds to normal momentum.
    !! Geopotential (var 6) has zero flux.
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,k,iel
    real(prec) :: nhat(1:3),nmag
    real(prec) :: rhoL,uL,vL,wL,thetaL,pL,unL,cL
    real(prec) :: rhoR,uR,vR,wR,thetaR,pR,unR,cR
    real(prec) :: rho_bar,c_bar,rc,un_star,p_star,gamma
    real(prec) :: sL(1:5),sR(1:5),s_up(1:5)

    gamma = this%cp/this%cv

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:6,iel=1:this%mesh%nElem)

      nhat = this%geometry%nHat%boundary(i,j,k,iEl,1,1:3)
      nmag = this%geometry%nScale%boundary(i,j,k,iEl,1)
      sL = this%solution%boundary(i,j,k,iel,1:5)
      sR = this%solution%extBoundary(i,j,k,iel,1:5)

      rhoL = sL(1)
      uL = sL(2)/rhoL
      vL = sL(3)/rhoL
      wL = sL(4)/rhoL
      thetaL = sL(5)/rhoL
      pL = this%p0*(rhoL*this%Rd*thetaL/this%p0)**gamma
      unL = uL*nhat(1)+vL*nhat(2)+wL*nhat(3)
      cL = sqrt(gamma*pL/rhoL)

      rhoR = sR(1)
      uR = sR(2)/rhoR
      vR = sR(3)/rhoR
      wR = sR(4)/rhoR
      thetaR = sR(5)/rhoR
      pR = this%p0*(rhoR*this%Rd*thetaR/this%p0)**gamma
      unR = uR*nhat(1)+vR*nhat(2)+wR*nhat(3)
      cR = sqrt(gamma*pR/rhoR)

      rho_bar = 0.5_prec*(rhoL+rhoR)
      c_bar = 0.5_prec*(cL+cR)
      rc = rho_bar*c_bar

      un_star = 0.5_prec*(unL+unR)-(pR-pL)/(2.0_prec*rc)
      p_star = 0.5_prec*(pL+pR)-0.5_prec*rc*(unR-unL)

      ! Upwind on un_star (advective part of LMARS)
      if(un_star >= 0.0_prec) then
        s_up = sL
      else
        s_up = sR
      endif

      this%flux%boundaryNormal(i,j,k,iEl,1) = (s_up(1)*un_star)*nmag
      this%flux%boundaryNormal(i,j,k,iEl,2) = (s_up(2)*un_star+p_star*nhat(1))*nmag
      this%flux%boundaryNormal(i,j,k,iEl,3) = (s_up(3)*un_star+p_star*nhat(2))*nmag
      this%flux%boundaryNormal(i,j,k,iEl,4) = (s_up(4)*un_star+p_star*nhat(3))*nmag
      this%flux%boundaryNormal(i,j,k,iEl,5) = (s_up(5)*un_star)*nmag
      this%flux%boundaryNormal(i,j,k,iEl,6) = 0.0_prec

    enddo

  endsubroutine BoundaryFlux_ESAtmo3D_t

  subroutine SourceMethod_ESAtmo3D_t(this)
    !! Souza et al. (2023) non-conservative gravity flux differencing.
    !!
    !! The rho*w equation carries the body force -rho * partial_z Phi where
    !! Phi = g*z is the geopotential, which we keep as state variable index
    !! 6 (no flux, no source of its own — its tendency is identically zero).
    !! On a curvilinear mesh, the SBP-EC strong-form flux differencing for
    !! the non-conservative product rho * partial_z Phi reads:
    !!
    !!   [rho * d_z Phi]_i = (1/J_i) sum_r sum_j D_split^r[i_r, j]
    !!                       * 0.5*(Ja^r_z(i) + Ja^r_z(j))
    !!                       * <rho>_log(s_i, s_j) * (Phi_j - Phi_i)
    !!
    !! and we set source(rho*w) = - that result. The other variables
    !! (rho, rho*u, rho*v, rho*theta, Phi) all have zero source.
    !!
    !! Why this form replaces the WB hydrostatic split: gravity now lives
    !! on the same node-pair carrier as the EC volume flux (every node
    !! pair, including those that span an element interface), so the
    !! gravity / pressure-flux topology mismatch that produced
    !! z-element-aligned banding is gone. Well-balancing is recovered
    !! automatically: in the hydrostatic state the EC volume flux gives
    !! d_z p_hyd and the source gives -rho_hyd*g, which cancel pointwise.
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,k,iEl,nn
    real(prec) :: rho_ijk,phi_ijk,rho_p,phi_p
    real(prec) :: rho_log,Ja_avg,acc,jac

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%solution%N+1,iEl=1:this%mesh%nElem)

      rho_ijk = this%solution%interior(i,j,k,iEl,1)
      phi_ijk = this%solution%interior(i,j,k,iEl,6)
      acc = 0.0_prec

      do nn = 1,this%solution%N+1
        ! xi^1 partner (nn, j, k)
        rho_p = this%solution%interior(nn,j,k,iEl,1)
        phi_p = this%solution%interior(nn,j,k,iEl,6)
        rho_log = log_mean(rho_ijk,rho_p)
        Ja_avg = 0.5_prec*(this%geometry%dsdx%interior(i,j,k,iEl,1,3,1)+ &
                           this%geometry%dsdx%interior(nn,j,k,iEl,1,3,1))
        acc = acc+this%solution%interp%dSplitMatrix(nn,i)* &
              Ja_avg*rho_log*(phi_p-phi_ijk)

        ! xi^2 partner (i, nn, k)
        rho_p = this%solution%interior(i,nn,k,iEl,1)
        phi_p = this%solution%interior(i,nn,k,iEl,6)
        rho_log = log_mean(rho_ijk,rho_p)
        Ja_avg = 0.5_prec*(this%geometry%dsdx%interior(i,j,k,iEl,1,3,2)+ &
                           this%geometry%dsdx%interior(i,nn,k,iEl,1,3,2))
        acc = acc+this%solution%interp%dSplitMatrix(nn,j)* &
              Ja_avg*rho_log*(phi_p-phi_ijk)

        ! xi^3 partner (i, j, nn)
        rho_p = this%solution%interior(i,j,nn,iEl,1)
        phi_p = this%solution%interior(i,j,nn,iEl,6)
        rho_log = log_mean(rho_ijk,rho_p)
        Ja_avg = 0.5_prec*(this%geometry%dsdx%interior(i,j,k,iEl,1,3,3)+ &
                           this%geometry%dsdx%interior(i,j,nn,iEl,1,3,3))
        acc = acc+this%solution%interp%dSplitMatrix(nn,k)* &
              Ja_avg*rho_log*(phi_p-phi_ijk)
      enddo

      jac = this%geometry%J%interior(i,j,k,iEl,1)
      this%source%interior(i,j,k,iEl,1) = 0.0_prec
      this%source%interior(i,j,k,iEl,2) = 0.0_prec
      this%source%interior(i,j,k,iEl,3) = 0.0_prec
      this%source%interior(i,j,k,iEl,4) = -acc/jac
      this%source%interior(i,j,k,iEl,5) = 0.0_prec
      this%source%interior(i,j,k,iEl,6) = 0.0_prec

    enddo

  endsubroutine SourceMethod_ESAtmo3D_t

  subroutine pbc3d_NoStress_ESAtmo3D(bc,mymodel)
    !! Parabolic boundary condition: zero diffusive flux normal to the wall
    !! (no-stress for momentum, no-heat-flux for rho*theta).
    !!
    !! Reflects the normal component of the interior solution gradient:
    !!
    !!   grad_ext = grad_int - 2 * (grad_int . n) * n
    !!
    !! After AverageSides() this gives avgGrad . n = 0 at every wall node,
    !! so the BR1 diffusive boundary flux f^diff = -coeff*(avgGrad.n)*nmag
    !! vanishes identically on no-normal-flow walls. Mirroring the gradient
    !! directly (grad_ext = grad_int) gives a non-zero flux through the
    !! wall and is NOT a no-stress condition.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel
    ! Local
    integer :: n,i,j,iEl,k,iVar
    real(prec) :: nhat(1:3),gn,g(1:3)

    select type(m => mymodel)
    class is(ESAtmo3D_t)
      do n = 1,bc%nBoundaries
        iEl = bc%elements(n)
        k = bc%sides(n)
        do j = 1,m%solution%interp%N+1
          do i = 1,m%solution%interp%N+1
            nhat = m%geometry%nHat%boundary(i,j,k,iEl,1,1:3)
            do iVar = 1,m%nvar
              g(1:3) = m%solutionGradient%boundary(i,j,k,iEl,iVar,1:3)
              gn = g(1)*nhat(1)+g(2)*nhat(2)+g(3)*nhat(3)
              m%solutionGradient%extBoundary(i,j,k,iEl,iVar,1) = g(1)-2.0_prec*gn*nhat(1)
              m%solutionGradient%extBoundary(i,j,k,iEl,iVar,2) = g(2)-2.0_prec*gn*nhat(2)
              m%solutionGradient%extBoundary(i,j,k,iEl,iVar,3) = g(3)-2.0_prec*gn*nhat(3)
            enddo
          enddo
        enddo
      enddo
    endselect

  endsubroutine pbc3d_NoStress_ESAtmo3D

  subroutine AdditionalInit_ESAtmo3D_t(this)
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    bcfunc => hbc3d_NoNormalFlow_ESAtmo3D
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

    ! Parabolic BC for the same wall tag: zero diffusive normal flux
    ! (no-stress / no-heat-flux). hyperbolicBCs and parabolicBCs are
    ! independent linked lists, so registering the same tag here does
    ! not clobber the hyperbolic registration above.
    bcfunc => pbc3d_NoStress_ESAtmo3D
    call this%parabolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

    ! Diffusive-flux scratch buffers. Always allocated (memory cost is
    ! modest); only used when nu>0 or kappa>0 (and gradient_enabled is
    ! therefore .true.).
    call this%diffFlux%Init(this%solution%interp,this%nvar,this%mesh%nElem)
    call this%diffFlux%AssociateGeometry(this%geometry)
    call this%diffDiv%Init(this%solution%interp,this%nvar,this%mesh%nElem)

  endsubroutine AdditionalInit_ESAtmo3D_t

  subroutine AdditionalFree_ESAtmo3D_t(this)
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this

    call this%diffFlux%Free()
    call this%diffDiv%Free()

  endsubroutine AdditionalFree_ESAtmo3D_t

  subroutine SetDiffusion_ESAtmo3D_t(this,nu,kappa,eta_penalty)
    !! Set the constant-coefficient Laplacian diffusion coefficients
    !! (kinematic momentum diffusivity and thermal diffusivity, both
    !! in m^2/s) and the dimensionless SIPG jump penalty. Setting nu
    !! or kappa > 0 enables the gradient pipeline so that the diffusive
    !! flux methods receive solutionGradient.
    !!
    !! length_scale is computed from the volume Jacobian: for a hex
    !! reference cell [-1,1]^3 with volume 8, the physical element
    !! volume is 8*<J> so the characteristic edge length is
    !! 2*<J>^(1/3). The smallest length over the mesh is conservative.
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this
    real(prec),intent(in) :: nu
    real(prec),intent(in) :: kappa
    real(prec),intent(in),optional :: eta_penalty
    ! Local
    real(prec) :: jmin

    this%nu = nu
    this%kappa = kappa
    if(present(eta_penalty)) this%eta_penalty = eta_penalty
    if(nu > 0.0_prec .or. kappa > 0.0_prec) then
      this%gradient_enabled = .true.
    endif

    jmin = minval(this%geometry%J%interior)
    this%length_scale = 2.0_prec*jmin**(1.0_prec/3.0_prec)

  endsubroutine SetDiffusion_ESAtmo3D_t

  subroutine DiffusiveFluxMethod_ESAtmo3D_t(this)
    !! Fill diffFlux%interior with the constant-coefficient Laplacian
    !! flux at every interior node:
    !!
    !!   F_d(rho)      = 0                              (no mass diffusion)
    !!   F_d(rho*v_i)  = -nu    * d(rho*v_i)/dx_d
    !!   F_d(rho*theta)= -kappa * d(rho*theta)/dx_d
    !!
    !! solutionGradient%interior(i,j,k,iel,iVar,d) holds d(s_iVar)/dx_d.
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,k,iel,d

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%solution%N+1,iel=1:this%mesh%nElem,d=1:3)
      this%diffFlux%interior(i,j,k,iel,1,d) = 0.0_prec
      this%diffFlux%interior(i,j,k,iel,2,d) = &
        -this%nu*this%solutionGradient%interior(i,j,k,iel,2,d)
      this%diffFlux%interior(i,j,k,iel,3,d) = &
        -this%nu*this%solutionGradient%interior(i,j,k,iel,3,d)
      this%diffFlux%interior(i,j,k,iel,4,d) = &
        -this%nu*this%solutionGradient%interior(i,j,k,iel,4,d)
      this%diffFlux%interior(i,j,k,iel,5,d) = &
        -this%kappa*this%solutionGradient%interior(i,j,k,iel,5,d)
      ! Geopotential carries no diffusion.
      this%diffFlux%interior(i,j,k,iel,6,d) = 0.0_prec
    enddo

  endsubroutine DiffusiveFluxMethod_ESAtmo3D_t

  subroutine DiffusiveBoundaryFlux_ESAtmo3D_t(this)
    !! Fill diffFlux%boundaryNormal with the SIPG-stabilised BR1 flux:
    !!
    !!   f_R^diff(iVar) = -coeff(iVar) * (avg_grad . n) * nmag
    !!                    + tau(iVar) * (uL - uR) * nmag
    !!
    !! tau(iVar) = eta_penalty * coeff(iVar) * (N+1)^2 / length_scale.
    !! avg_grad is solutionGradient%avgBoundary (populated by
    !! AverageSides()); uL, uR are solution%boundary, %extBoundary.
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,k,iel
    real(prec) :: nhat(1:3),nmag,gradn,coeff,tau,np2
    real(prec) :: uL,uR

    np2 = real((this%solution%interp%N+1)**2,prec)

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:6,iel=1:this%mesh%nElem)
      nhat = this%geometry%nHat%boundary(i,j,k,iel,1,1:3)
      nmag = this%geometry%nScale%boundary(i,j,k,iel,1)

      ! rho — no mass diffusion, no penalty
      this%diffFlux%boundaryNormal(i,j,k,iel,1) = 0.0_prec

      ! Momentum equations use nu
      coeff = this%nu
      tau = this%eta_penalty*coeff*np2/this%length_scale

      gradn = this%solutionGradient%avgBoundary(i,j,k,iel,2,1)*nhat(1)+ &
              this%solutionGradient%avgBoundary(i,j,k,iel,2,2)*nhat(2)+ &
              this%solutionGradient%avgBoundary(i,j,k,iel,2,3)*nhat(3)
      uL = this%solution%boundary(i,j,k,iel,2)
      uR = this%solution%extBoundary(i,j,k,iel,2)
      this%diffFlux%boundaryNormal(i,j,k,iel,2) = (-coeff*gradn+tau*(uL-uR))*nmag

      gradn = this%solutionGradient%avgBoundary(i,j,k,iel,3,1)*nhat(1)+ &
              this%solutionGradient%avgBoundary(i,j,k,iel,3,2)*nhat(2)+ &
              this%solutionGradient%avgBoundary(i,j,k,iel,3,3)*nhat(3)
      uL = this%solution%boundary(i,j,k,iel,3)
      uR = this%solution%extBoundary(i,j,k,iel,3)
      this%diffFlux%boundaryNormal(i,j,k,iel,3) = (-coeff*gradn+tau*(uL-uR))*nmag

      gradn = this%solutionGradient%avgBoundary(i,j,k,iel,4,1)*nhat(1)+ &
              this%solutionGradient%avgBoundary(i,j,k,iel,4,2)*nhat(2)+ &
              this%solutionGradient%avgBoundary(i,j,k,iel,4,3)*nhat(3)
      uL = this%solution%boundary(i,j,k,iel,4)
      uR = this%solution%extBoundary(i,j,k,iel,4)
      this%diffFlux%boundaryNormal(i,j,k,iel,4) = (-coeff*gradn+tau*(uL-uR))*nmag

      ! rho*theta uses kappa
      coeff = this%kappa
      tau = this%eta_penalty*coeff*np2/this%length_scale

      gradn = this%solutionGradient%avgBoundary(i,j,k,iel,5,1)*nhat(1)+ &
              this%solutionGradient%avgBoundary(i,j,k,iel,5,2)*nhat(2)+ &
              this%solutionGradient%avgBoundary(i,j,k,iel,5,3)*nhat(3)
      uL = this%solution%boundary(i,j,k,iel,5)
      uR = this%solution%extBoundary(i,j,k,iel,5)
      this%diffFlux%boundaryNormal(i,j,k,iel,5) = (-coeff*gradn+tau*(uL-uR))*nmag

      ! Geopotential — no diffusion, no penalty.
      this%diffFlux%boundaryNormal(i,j,k,iel,6) = 0.0_prec
    enddo

  endsubroutine DiffusiveBoundaryFlux_ESAtmo3D_t

  subroutine CalculateTendency_ESAtmo3D_t(this)
    !! ESAtmo3D tendency = EC inviscid pipeline (parent) + optional
    !! constant-coefficient Laplacian diffusion (BR1 weak-form DG).
    !!
    !! When nu>0 or kappa>0 the parent's CalculateTendency already
    !! computes solutionGradient (because gradient_enabled was set by
    !! SetDiffusion); we then fill diffFlux from solutionGradient,
    !! compute its DG divergence, and accumulate into fluxDivergence
    !! before forming dSdt.
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,k,iEl,iVar
    real(prec) :: bMi1,bMi2,bMj1,bMj2,bMk1,bMk2,qwi,qwj,qwk,jac

    call this%solution%BoundaryInterp()
    call this%solution%SideExchange(this%mesh)

    call this%PreTendencyHook()
    call this%SetBoundaryCondition()

    if(this%gradient_enabled) then
      call this%CalculateSolutionGradient()
      call this%SetGradientBoundaryCondition()
      call this%solutionGradient%AverageSides()
    endif

    call this%SourceMethod()
    call this%BoundaryFlux()

    call this%TwoPointFluxMethod()
    call this%twoPointFlux%UpdateDevice()
    call this%twoPointFlux%MappedDivergence(this%fluxDivergence%interior)

    ! EC surface contribution (1/J) M^{-1} B^T f_Riemann
    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%solution%N+1,iEl=1:this%mesh%nElem, &
                  iVar=1:this%solution%nVar)

      bMi1 = this%solution%interp%bMatrix(i,2)
      bMi2 = this%solution%interp%bMatrix(i,1)
      bMj1 = this%solution%interp%bMatrix(j,2)
      bMj2 = this%solution%interp%bMatrix(j,1)
      bMk1 = this%solution%interp%bMatrix(k,2)
      bMk2 = this%solution%interp%bMatrix(k,1)
      qwi = this%solution%interp%qWeights(i)
      qwj = this%solution%interp%qWeights(j)
      qwk = this%solution%interp%qWeights(k)
      jac = this%geometry%J%interior(i,j,k,iEl,1)

      this%fluxDivergence%interior(i,j,k,iEl,iVar) = &
        this%fluxDivergence%interior(i,j,k,iEl,iVar)+ &
        (bMk1*this%flux%boundaryNormal(i,j,6,iEl,iVar)+ &
         bMk2*this%flux%boundaryNormal(i,j,1,iEl,iVar))/(qwk*jac)+ &
        (bMi1*this%flux%boundaryNormal(j,k,3,iEl,iVar)+ &
         bMi2*this%flux%boundaryNormal(j,k,5,iEl,iVar))/(qwi*jac)+ &
        (bMj1*this%flux%boundaryNormal(i,k,4,iEl,iVar)+ &
         bMj2*this%flux%boundaryNormal(i,k,2,iEl,iVar))/(qwj*jac)
    enddo

    ! Add the diffusive contribution if any coefficient is nonzero.
    if(this%nu > 0.0_prec .or. this%kappa > 0.0_prec) then
      call this%DiffusiveFluxMethod()
      call this%DiffusiveBoundaryFlux()
      call this%diffFlux%MappedDGDivergence(this%diffDiv%interior)
      do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                    k=1:this%solution%N+1,iEl=1:this%mesh%nElem, &
                    iVar=1:this%solution%nVar)
        ! Subtract because MappedDGDivergence returns div(F); we want
        ! the diffusion to APPEAR in dSdt as +div(F_diff), so the
        ! contribution to fluxDivergence (which is later subtracted to
        ! form dSdt) is -div(F_diff). diffFlux already includes the
        ! minus sign of -nu*grad / -kappa*grad, so MappedDGDivergence
        ! returns -nu*Lap (the right sign for fluxDivergence).
        this%fluxDivergence%interior(i,j,k,iEl,iVar) = &
          this%fluxDivergence%interior(i,j,k,iEl,iVar)+ &
          this%diffDiv%interior(i,j,k,iEl,iVar)
      enddo
    endif

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%solution%N+1,iEl=1:this%mesh%nElem, &
                  iVar=1:this%solution%nVar)

      this%dSdt%interior(i,j,k,iEl,iVar) = &
        this%source%interior(i,j,k,iEl,iVar)- &
        this%fluxDivergence%interior(i,j,k,iEl,iVar)
    enddo

  endsubroutine CalculateTendency_ESAtmo3D_t

  subroutine hbc3d_NoNormalFlow_ESAtmo3D(bc,mymodel)
    !! No-normal-flow (wall) boundary condition.
    !!
    !! The normal component of momentum is negated while tangential
    !! components are mirrored (preserved). Density and potential
    !! temperature (tracers) are mirrored from the interior.
    !!
    !! Exterior momentum:
    !!   (rho*v)_ext = (rho*v)_int - 2 * ((rho*v)_int . nhat) * nhat
    !!
    !! This ensures v_ext . n = -(v_int . n) so the Riemann solver
    !! sees zero normal velocity at the wall.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel
    ! Local
    integer :: n,i,j,iEl,k
    real(prec) :: nhat(1:3),rhovn

    select type(m => mymodel)
    class is(ESAtmo3D_t)
      do n = 1,bc%nBoundaries
        iEl = bc%elements(n)
        k = bc%sides(n)
        do j = 1,m%solution%interp%N+1
          do i = 1,m%solution%interp%N+1

            nhat = m%geometry%nHat%boundary(i,j,k,iEl,1,1:3)

            ! Mirror density (tracer)
            m%solution%extBoundary(i,j,k,iEl,1) = &
              m%solution%boundary(i,j,k,iEl,1)

            ! Reflect momentum: negate normal component, preserve tangential
            ! (rho*v)_ext = (rho*v)_int - 2*((rho*v)_int . n)*n
            rhovn = m%solution%boundary(i,j,k,iEl,2)*nhat(1)+ &
                    m%solution%boundary(i,j,k,iEl,3)*nhat(2)+ &
                    m%solution%boundary(i,j,k,iEl,4)*nhat(3)

            m%solution%extBoundary(i,j,k,iEl,2) = &
              m%solution%boundary(i,j,k,iEl,2)-2.0_prec*rhovn*nhat(1)
            m%solution%extBoundary(i,j,k,iEl,3) = &
              m%solution%boundary(i,j,k,iEl,3)-2.0_prec*rhovn*nhat(2)
            m%solution%extBoundary(i,j,k,iEl,4) = &
              m%solution%boundary(i,j,k,iEl,4)-2.0_prec*rhovn*nhat(3)

            ! Mirror potential temperature (tracer)
            m%solution%extBoundary(i,j,k,iEl,5) = &
              m%solution%boundary(i,j,k,iEl,5)

            ! Mirror geopotential — Phi only depends on z, so the
            ! mirror across a no-normal-flow wall is the same value.
            m%solution%extBoundary(i,j,k,iEl,6) = &
              m%solution%boundary(i,j,k,iEl,6)

          enddo
        enddo
      enddo
    endselect

  endsubroutine hbc3d_NoNormalFlow_ESAtmo3D

  subroutine SetHydrostaticBalance_ESAtmo3D_t(this,theta0)
    !! Initialise a hydrostatically balanced atmosphere with uniform
    !! potential temperature theta0, zero velocity, and the geopotential
    !! Phi = g*z carried as state variable index 6.
    !!
    !! Exner function: pi(z) = 1 - g*z / (cp*theta0)
    !! rho(z)         = p0/(Rd*theta0) * pi(z)^(cv/Rd)
    !! rho*theta(z)   = p0/Rd * pi(z)^(cv/Rd)
    !!
    !! Boundary and extBoundary buffers of the solution are populated
    !! analytically from per-face z so that face values are bit-exact
    !! with what BoundaryInterp would produce (and extBoundary at walls
    !! mirrors). SideExchange will overwrite extBoundary at interior
    !! element interfaces with the neighbour's value — for a smooth
    !! profile this is a no-op.
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this
    real(prec),intent(in) :: theta0
    ! Local
    integer :: i,j,k,iEl
    integer :: side
    real(prec) :: z,exner,rho

    print*,__FILE__," : Setting hydrostatic balance with theta0 = ",theta0

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%solution%N+1,iEl=1:this%mesh%nElem)

      z = this%geometry%x%interior(i,j,k,iEl,1,3)
      exner = 1.0_prec-this%g*z/(this%cp*theta0)
      rho = this%p0/(this%Rd*theta0)*exner**(this%cv/this%Rd)

      this%solution%interior(i,j,k,iEl,1) = rho
      this%solution%interior(i,j,k,iEl,2) = 0.0_prec
      this%solution%interior(i,j,k,iEl,3) = 0.0_prec
      this%solution%interior(i,j,k,iEl,4) = 0.0_prec
      this%solution%interior(i,j,k,iEl,5) = rho*theta0
      this%solution%interior(i,j,k,iEl,6) = this%g*z

    enddo

    ! Populate boundary + extBoundary buffers analytically from the
    ! per-face z-coordinates. Geopotential at the boundary mirrors
    ! by construction (Phi = g*z is single-valued in z).
    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  side=1:6,iEl=1:this%mesh%nElem)

      z = this%geometry%x%boundary(i,j,side,iEl,1,3)
      exner = 1.0_prec-this%g*z/(this%cp*theta0)
      rho = this%p0/(this%Rd*theta0)*exner**(this%cv/this%Rd)

      this%solution%boundary(i,j,side,iEl,1) = rho
      this%solution%boundary(i,j,side,iEl,2) = 0.0_prec
      this%solution%boundary(i,j,side,iEl,3) = 0.0_prec
      this%solution%boundary(i,j,side,iEl,4) = 0.0_prec
      this%solution%boundary(i,j,side,iEl,5) = rho*theta0
      this%solution%boundary(i,j,side,iEl,6) = this%g*z

      this%solution%extBoundary(i,j,side,iEl,1) = rho
      this%solution%extBoundary(i,j,side,iEl,2) = 0.0_prec
      this%solution%extBoundary(i,j,side,iEl,3) = 0.0_prec
      this%solution%extBoundary(i,j,side,iEl,4) = 0.0_prec
      this%solution%extBoundary(i,j,side,iEl,5) = rho*theta0
      this%solution%extBoundary(i,j,side,iEl,6) = this%g*z

    enddo

    call this%ReportMetrics()
    call this%solution%UpdateDevice()

  endsubroutine SetHydrostaticBalance_ESAtmo3D_t

  subroutine AddThermalBubble_ESAtmo3D_t(this,dtheta,r0,x0,y0,z0)
    !! Adds a pressure-balanced warm bubble perturbation.
    !!
    !! The potential temperature perturbation has a cos^2 profile:
    !!
    !!   theta'(x,y,z) = dtheta * cos^2(pi*r / (2*r0))  for r <= r0
    !!   theta'(x,y,z) = 0                                for r >  r0
    !!
    !! where r = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2).
    !!
    !! To maintain pressure balance (p unchanged), the density is
    !! adjusted while rho*theta remains constant:
    !!
    !!   theta_new = theta_old + theta'
    !!   rho_new   = rho_old * theta_old / theta_new
    !!   (rho*theta)_new = rho_old * theta_old = (rho*theta)_old
    !!
    !! The buoyancy force arises from the density deficit in the
    !! gravitational source term -rho*g.
    implicit none
    class(ESAtmo3D_t),intent(inout) :: this
    real(prec),intent(in) :: dtheta ! Perturbation amplitude [K]
    real(prec),intent(in) :: r0 ! Bubble radius [m]
    real(prec),intent(in) :: x0,y0,z0 ! Bubble center [m]
    ! Local
    integer :: i,j,k,iEl
    real(prec) :: x,y,z,r,rho_old,theta_old,theta_new,thetap

    print*,__FILE__," : Adding thermal bubble perturbation"
    print*,__FILE__," : dtheta = ",dtheta
    print*,__FILE__," : r0     = ",r0
    print*,__FILE__," : center = ",x0,y0,z0

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%solution%N+1,iEl=1:this%mesh%nElem)

      x = this%geometry%x%interior(i,j,k,iEl,1,1)-x0
      y = this%geometry%x%interior(i,j,k,iEl,1,2)-y0
      z = this%geometry%x%interior(i,j,k,iEl,1,3)-z0
      r = sqrt(x*x+y*y+z*z)

      if(r <= r0) then
        rho_old = this%solution%interior(i,j,k,iEl,1)
        theta_old = this%solution%interior(i,j,k,iEl,5)/rho_old
        thetap = dtheta*cos(pi*r/(2.0_prec*r0))**2
        theta_new = theta_old+thetap

        ! Adjust density to maintain pressure balance: rho*theta = const
        this%solution%interior(i,j,k,iEl,1) = rho_old*theta_old/theta_new
        ! rho*theta is unchanged (pressure-balanced perturbation)
      endif

    enddo

    call this%ReportMetrics()
    call this%solution%UpdateDevice()

  endsubroutine AddThermalBubble_ESAtmo3D_t

endmodule SELF_ESAtmo3D_t
