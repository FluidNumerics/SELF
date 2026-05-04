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

module SELF_ESAtmo2D_t
  !! Entropy-conserving compressible Euler equations in 2-D with potential
  !! temperature formulation (FastEddy / atmospheric LES equation set).
  !!
  !! Conservative variables:
  !!
  !!   U = [rho, rho*u, rho*v, rho*theta, Phi]
  !!
  !! where theta is potential temperature and Phi = g*y is the geopotential
  !! (gravity acts in the +y direction; y is the vertical coordinate).
  !!
  !! Governing equations:
  !!
  !!   d(rho)/dt     + div(rho * v)         = 0
  !!   d(rho*u)/dt   + div(rho*u * v) + dp/dx = 0
  !!   d(rho*v)/dt   + div(rho*v * v) + dp/dy = -rho * g
  !!   d(rho*th)/dt  + div(rho*th * v)      = 0
  !!
  !! Equation of state (dry ideal gas):
  !!
  !!   p = p0 * (rho * Rd * theta / p0)^(cp/cv)
  !!
  !! Two-point volume flux: Souza et al. (2023) entropy-conservative.
  !!
  !! Surface Riemann flux: LMARS (Chen et al. 2013).

  use SELF_ECDGModel2D
  use SELF_MappedScalar_2D
  use SELF_mesh
  use SELF_BoundaryConditions

  implicit none

  type,extends(ECDGModel2D) :: ESAtmo2D_t

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
    !!
    !!   nu             : kinematic diffusivity for momentum [m^2/s]
    !!   kappa          : thermal diffusivity for rho*theta  [m^2/s]
    !!   eta_penalty    : dimensionless SIPG penalty (default 4.0)
    !!   length_scale   : characteristic element length [m], filled by
    !!                    SetDiffusion from the geometry (mean J^(1/2)*2).
    real(prec) :: nu = 0.0_prec
    real(prec) :: kappa = 0.0_prec
    real(prec) :: eta_penalty = 4.0_prec
    real(prec) :: length_scale = 0.0_prec

    !! Diffusive-flux scratch buffer. Holds the constant-coefficient
    !! Laplacian flux F_diff(i,j,iel,iVar,d) = -coeff_iVar * d(s_iVar)/dx_d.
    type(MappedVector2D) :: diffFlux
    type(MappedScalar2D) :: diffDiv

  contains

    procedure :: SetNumberOfVariables => SetNumberOfVariables_ESAtmo2D_t
    procedure :: SetMetadata => SetMetadata_ESAtmo2D_t
    procedure :: entropy_func => entropy_func_ESAtmo2D_t
    procedure :: twopointflux2d => twopointflux2d_ESAtmo2D_t
    procedure :: TwoPointFluxMethod => TwoPointFluxMethod_ESAtmo2D_t
    procedure :: riemannflux2d => riemannflux2d_ESAtmo2D_t
    procedure :: BoundaryFlux => BoundaryFlux_ESAtmo2D_t
    procedure :: SourceMethod => SourceMethod_ESAtmo2D_t
    procedure :: AdditionalInit => AdditionalInit_ESAtmo2D_t
    procedure :: AdditionalFree => AdditionalFree_ESAtmo2D_t
    procedure :: SetHydrostaticBalance => SetHydrostaticBalance_ESAtmo2D_t
    procedure :: AddThermalBubble => AddThermalBubble_ESAtmo2D_t

    !! Constant-coefficient Laplacian / Bassi-Rebay diffusion hooks
    procedure :: SetDiffusion => SetDiffusion_ESAtmo2D_t
    procedure :: DiffusiveFluxMethod => DiffusiveFluxMethod_ESAtmo2D_t
    procedure :: DiffusiveBoundaryFlux => DiffusiveBoundaryFlux_ESAtmo2D_t
    procedure :: CalculateTendency => CalculateTendency_ESAtmo2D_t

  endtype ESAtmo2D_t

contains

  subroutine SetNumberOfVariables_ESAtmo2D_t(this)
    !! Five conserved variables: (rho, rho*u, rho*v, rho*theta, Phi),
    !! where Phi = g*y is the geopotential. Phi has zero flux (volume
    !! and surface) so its tendency is identically zero; it is carried
    !! in the state vector solely so that the Souza et al. (2023)
    !! non-conservative gravity flux differencing in SourceMethod can
    !! read it node-locally.
    implicit none
    class(ESAtmo2D_t),intent(inout) :: this

    this%nvar = 5

  endsubroutine SetNumberOfVariables_ESAtmo2D_t

  subroutine SetMetadata_ESAtmo2D_t(this)
    implicit none
    class(ESAtmo2D_t),intent(inout) :: this

    call this%solution%SetName(1,"rho")
    call this%solution%SetUnits(1,"kg/m^2")

    call this%solution%SetName(2,"rhou")
    call this%solution%SetUnits(2,"kg/(m s)")

    call this%solution%SetName(3,"rhov")
    call this%solution%SetUnits(3,"kg/(m s)")

    call this%solution%SetName(4,"rhotheta")
    call this%solution%SetUnits(4,"kg K/m^2")

    call this%solution%SetName(5,"phi")
    call this%solution%SetUnits(5,"m^2/s^2")

  endsubroutine SetMetadata_ESAtmo2D_t

  pure function entropy_func_ESAtmo2D_t(this,s) result(e)
    !! Mathematical entropy: total energy density (kinetic + internal).
    !!
    !!   e = 0.5*(rhou^2 + rhov^2)/rho + p/(gamma - 1)
    !!
    !! where p = p0 * (rho * Rd * theta / p0)^gamma and gamma = cp/cv.
    class(ESAtmo2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: e
    ! Local
    real(prec) :: rho,theta,p,gamma

    rho = s(1)
    theta = s(4)/rho
    gamma = this%cp/this%cv
    p = this%p0*(rho*this%Rd*theta/this%p0)**gamma

    e = 0.5_prec*(s(2)*s(2)+s(3)*s(3))/rho+ &
        p/(gamma-1.0_prec)

  endfunction entropy_func_ESAtmo2D_t

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

  pure function twopointflux2d_ESAtmo2D_t(this,sL,sR) result(flux)
    !! Souza et al. (2023, JAMES) entropy-conservative two-point flux for
    !! 2-D compressible Euler in (rho, rho*v, rho*theta) variables with
    !! p = p0*(rho*Rd*theta/p0)^gamma.
    !!
    !!   f_d(rho)     = <rho>_log * <v_d>
    !!   f_d(rho*v_i) = <rho>_log * <v_i> * <v_d> + <p> * delta_{id}
    !!   f_d(rho*th)  = <rho*theta>_log * <v_d>
    !!   f_d(Phi)     = 0  (geopotential is carried, not advected)
    !!
    !! Pressure here is the *total* pressure; gravity is handled by the
    !! Souza non-conservative term in SourceMethod, not by a flux split.
    class(ESAtmo2D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec) :: flux(1:this%nvar,1:2)
    ! Local
    real(prec) :: rhoL,uL,vL,thetaL,pL,rthL
    real(prec) :: rhoR,uR,vR,thetaR,pR,rthR
    real(prec) :: rho_log,rth_log,u_avg,v_avg,p_avg
    real(prec) :: gamma

    gamma = this%cp/this%cv

    ! Left primitive variables
    rhoL = sL(1)
    rthL = sL(4)
    thetaL = rthL/rhoL
    uL = sL(2)/rhoL
    vL = sL(3)/rhoL
    pL = this%p0*(rhoL*this%Rd*thetaL/this%p0)**gamma

    ! Right primitive variables
    rhoR = sR(1)
    rthR = sR(4)
    thetaR = rthR/rhoR
    uR = sR(2)/rhoR
    vR = sR(3)/rhoR
    pR = this%p0*(rhoR*this%Rd*thetaR/this%p0)**gamma

    ! Logarithmic and arithmetic means
    rho_log = log_mean(rhoL,rhoR)
    rth_log = log_mean(rthL,rthR)
    u_avg = 0.5_prec*(uL+uR)
    v_avg = 0.5_prec*(vL+vR)
    p_avg = 0.5_prec*(pL+pR)

    ! x-direction flux (d=1)
    flux(1,1) = rho_log*u_avg
    flux(2,1) = rho_log*u_avg*u_avg+p_avg
    flux(3,1) = rho_log*v_avg*u_avg
    flux(4,1) = rth_log*u_avg
    flux(5,1) = 0.0_prec

    ! y-direction flux (d=2)
    flux(1,2) = rho_log*v_avg
    flux(2,2) = rho_log*u_avg*v_avg
    flux(3,2) = rho_log*v_avg*v_avg+p_avg
    flux(4,2) = rth_log*v_avg
    flux(5,2) = 0.0_prec

  endfunction twopointflux2d_ESAtmo2D_t

  subroutine TwoPointFluxMethod_ESAtmo2D_t(this)
    !! Pre-projected scalar contravariant two-point Souza et al. (2023) EC
    !! flux. For each node pair (a, b) along reference direction r:
    !!
    !!   Fc^r_v_i = sum_d 0.5*(Ja^r_d(a) + Ja^r_d(b)) * f_d(v_i, sL=s(a), sR=s(b))
    !!
    !! Gravity is NOT split into the pressure flux here; it is handled by
    !! the Souza non-conservative gravity flux differencing in SourceMethod
    !! (which uses the geopotential carried as state variable index 5).
    implicit none
    class(ESAtmo2D_t),intent(inout) :: this
    ! Local
    integer :: nn,i,j,d,iEl,iVar
    real(prec) :: sL(1:this%nvar),sR(1:this%nvar)
    real(prec) :: Fphys(1:this%nvar,1:2)
    real(prec) :: Fc

    do concurrent(nn=1:this%solution%N+1,i=1:this%solution%N+1, &
                  j=1:this%solution%N+1,iEl=1:this%mesh%nElem)

      sL = this%solution%interior(i,j,iEl,1:this%nvar)

      ! ------------- xi^1: pair (i,j)-(nn,j) -------------
      sR = this%solution%interior(nn,j,iEl,1:this%nvar)
      Fphys = this%twopointflux2d(sL,sR)
      do iVar = 1,this%nvar
        Fc = 0.0_prec
        do d = 1,2
          Fc = Fc+0.5_prec*( &
               this%geometry%dsdx%interior(i,j,iEl,1,d,1)+ &
               this%geometry%dsdx%interior(nn,j,iEl,1,d,1))* &
               Fphys(iVar,d)
        enddo
        this%twoPointFlux%interior(nn,i,j,iEl,iVar,1) = Fc
      enddo

      ! ------------- xi^2: pair (i,j)-(i,nn) -------------
      sR = this%solution%interior(i,nn,iEl,1:this%nvar)
      Fphys = this%twopointflux2d(sL,sR)
      do iVar = 1,this%nvar
        Fc = 0.0_prec
        do d = 1,2
          Fc = Fc+0.5_prec*( &
               this%geometry%dsdx%interior(i,j,iEl,1,d,2)+ &
               this%geometry%dsdx%interior(i,nn,iEl,1,d,2))* &
               Fphys(iVar,d)
        enddo
        this%twoPointFlux%interior(nn,i,j,iEl,iVar,2) = Fc
      enddo

    enddo

  endsubroutine TwoPointFluxMethod_ESAtmo2D_t

  pure function riemannflux2d_ESAtmo2D_t(this,sL,sR,dsdx,nhat) result(flux)
    !! Local Lax-Friedrichs (Rusanov) Riemann flux. Provided as a fallback;
    !! the model overrides BoundaryFlux directly with the LMARS solver.
    !!
    !!   F* = 0.5*(fL.n + fR.n) - 0.5*lambda_max*(sR - sL)
    !!
    !! where lambda_max = max(|vL.n| + cL, |vR.n| + cR)
    !! and c = sqrt(gamma * p / rho) is the sound speed.
    class(ESAtmo2D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: flux(1:this%nvar)
    ! Local
    real(prec) :: rhoL,uL,vL,thetaL,pL,unL,cL
    real(prec) :: rhoR,uR,vR,thetaR,pR,unR,cR
    real(prec) :: fL(1:4),fR(1:4)
    real(prec) :: gamma,lam

    gamma = this%cp/this%cv

    ! Left state
    rhoL = sL(1)
    uL = sL(2)/rhoL
    vL = sL(3)/rhoL
    thetaL = sL(4)/rhoL
    pL = this%p0*(rhoL*this%Rd*thetaL/this%p0)**gamma
    unL = uL*nhat(1)+vL*nhat(2)
    cL = sqrt(gamma*pL/rhoL)

    ! Right state
    rhoR = sR(1)
    uR = sR(2)/rhoR
    vR = sR(3)/rhoR
    thetaR = sR(4)/rhoR
    pR = this%p0*(rhoR*this%Rd*thetaR/this%p0)**gamma
    unR = uR*nhat(1)+vR*nhat(2)
    cR = sqrt(gamma*pR/rhoR)

    ! Normal flux from left state
    fL(1) = rhoL*unL
    fL(2) = sL(2)*unL+pL*nhat(1)
    fL(3) = sL(3)*unL+pL*nhat(2)
    fL(4) = sL(4)*unL

    ! Normal flux from right state
    fR(1) = rhoR*unR
    fR(2) = sR(2)*unR+pR*nhat(1)
    fR(3) = sR(3)*unR+pR*nhat(2)
    fR(4) = sR(4)*unR

    ! Maximum wave speed
    lam = max(abs(unL)+cL,abs(unR)+cR)

    ! LLF flux
    flux(1:4) = 0.5_prec*(fL(1:4)+fR(1:4))- &
                0.5_prec*lam*(sR(1:4)-sL(1:4))
    ! Geopotential is carried in the state vector but has no flux.
    flux(5) = 0.0_prec

    if(.false.) flux(1) = flux(1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction riemannflux2d_ESAtmo2D_t

  subroutine BoundaryFlux_ESAtmo2D_t(this)
    !! LMARS (Low-Mach Approximate Riemann Solver, Chen et al. 2013)
    !! interface flux. No hydrostatic pressure split: gravity is folded
    !! into SourceMethod via the Souza non-conservative form using the
    !! geopotential carried in the state vector (variable index 5).
    !!
    !!   un* = 0.5*(unL + unR) - (pR - pL) / (2 * rho_bar * c_bar)
    !!   p*  = 0.5*(pL + pR)   - 0.5 * rho_bar * c_bar * (unR - unL)
    !!
    !! Upwind on sign of un* selects which side supplies the conserved
    !! state for the advective part. Pressure adds to normal momentum.
    !! Geopotential (var 5) has zero flux.
    implicit none
    class(ESAtmo2D_t),intent(inout) :: this
    ! Local
    integer :: i,k,iel
    real(prec) :: nhat(1:2),nmag
    real(prec) :: rhoL,uL,vL,thetaL,pL,unL,cL
    real(prec) :: rhoR,uR,vR,thetaR,pR,unR,cR
    real(prec) :: rho_bar,c_bar,rc,un_star,p_star,gamma
    real(prec) :: sL(1:4),sR(1:4),s_up(1:4)

    gamma = this%cp/this%cv

    do concurrent(i=1:this%solution%N+1,k=1:4,iel=1:this%mesh%nElem)

      nhat = this%geometry%nHat%boundary(i,k,iEl,1,1:2)
      nmag = this%geometry%nScale%boundary(i,k,iEl,1)
      sL = this%solution%boundary(i,k,iel,1:4)
      sR = this%solution%extBoundary(i,k,iel,1:4)

      rhoL = sL(1)
      uL = sL(2)/rhoL
      vL = sL(3)/rhoL
      thetaL = sL(4)/rhoL
      pL = this%p0*(rhoL*this%Rd*thetaL/this%p0)**gamma
      unL = uL*nhat(1)+vL*nhat(2)
      cL = sqrt(gamma*pL/rhoL)

      rhoR = sR(1)
      uR = sR(2)/rhoR
      vR = sR(3)/rhoR
      thetaR = sR(4)/rhoR
      pR = this%p0*(rhoR*this%Rd*thetaR/this%p0)**gamma
      unR = uR*nhat(1)+vR*nhat(2)
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

      this%flux%boundaryNormal(i,k,iEl,1) = (s_up(1)*un_star)*nmag
      this%flux%boundaryNormal(i,k,iEl,2) = (s_up(2)*un_star+p_star*nhat(1))*nmag
      this%flux%boundaryNormal(i,k,iEl,3) = (s_up(3)*un_star+p_star*nhat(2))*nmag
      this%flux%boundaryNormal(i,k,iEl,4) = (s_up(4)*un_star)*nmag
      this%flux%boundaryNormal(i,k,iEl,5) = 0.0_prec

    enddo

  endsubroutine BoundaryFlux_ESAtmo2D_t

  subroutine SourceMethod_ESAtmo2D_t(this)
    !! Souza et al. (2023) non-conservative gravity flux differencing.
    !!
    !! The rho*v equation carries the body force -rho * partial_y Phi where
    !! Phi = g*y is the geopotential (state variable index 5; no flux, no
    !! source of its own). On a curvilinear mesh, the SBP-EC strong-form
    !! flux differencing for the non-conservative product rho * partial_y Phi
    !! reads:
    !!
    !!   [rho * d_y Phi]_i = (1/J_i) sum_r sum_j D_split^r[i_r, j]
    !!                       * 0.5*(Ja^r_y(i) + Ja^r_y(j))
    !!                       * <rho>_log(s_i, s_j) * (Phi_j - Phi_i)
    !!
    !! and we set source(rho*v) = - that result. The other variables
    !! (rho, rho*u, rho*theta, Phi) all have zero source.
    implicit none
    class(ESAtmo2D_t),intent(inout) :: this
    ! Local
    integer :: i,j,iEl,nn
    real(prec) :: rho_ijk,phi_ijk,rho_p,phi_p
    real(prec) :: rho_log,Ja_avg,acc,jac

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iEl=1:this%mesh%nElem)

      rho_ijk = this%solution%interior(i,j,iEl,1)
      phi_ijk = this%solution%interior(i,j,iEl,5)
      acc = 0.0_prec

      do nn = 1,this%solution%N+1
        ! xi^1 partner (nn, j) — uses Ja^1_y (d=2, r=1)
        rho_p = this%solution%interior(nn,j,iEl,1)
        phi_p = this%solution%interior(nn,j,iEl,5)
        rho_log = log_mean(rho_ijk,rho_p)
        Ja_avg = 0.5_prec*(this%geometry%dsdx%interior(i,j,iEl,1,2,1)+ &
                           this%geometry%dsdx%interior(nn,j,iEl,1,2,1))
        acc = acc+this%solution%interp%dSplitMatrix(nn,i)* &
              Ja_avg*rho_log*(phi_p-phi_ijk)

        ! xi^2 partner (i, nn) — uses Ja^2_y (d=2, r=2)
        rho_p = this%solution%interior(i,nn,iEl,1)
        phi_p = this%solution%interior(i,nn,iEl,5)
        rho_log = log_mean(rho_ijk,rho_p)
        Ja_avg = 0.5_prec*(this%geometry%dsdx%interior(i,j,iEl,1,2,2)+ &
                           this%geometry%dsdx%interior(i,nn,iEl,1,2,2))
        acc = acc+this%solution%interp%dSplitMatrix(nn,j)* &
              Ja_avg*rho_log*(phi_p-phi_ijk)
      enddo

      jac = this%geometry%J%interior(i,j,iEl,1)
      this%source%interior(i,j,iEl,1) = 0.0_prec
      this%source%interior(i,j,iEl,2) = 0.0_prec
      this%source%interior(i,j,iEl,3) = -acc/jac
      this%source%interior(i,j,iEl,4) = 0.0_prec
      this%source%interior(i,j,iEl,5) = 0.0_prec

    enddo

  endsubroutine SourceMethod_ESAtmo2D_t

  subroutine pbc2d_NoStress_ESAtmo2D(bc,mymodel)
    !! Parabolic boundary condition: zero diffusive flux normal to the wall
    !! (no-stress for momentum, no-heat-flux for rho*theta).
    !!
    !! Reflects the normal component of the interior solution gradient:
    !!
    !!   grad_ext = grad_int - 2 * (grad_int . n) * n
    !!
    !! After AverageSides() this gives avgGrad . n = 0 at every wall node,
    !! so the BR1 diffusive boundary flux f^diff = -coeff*(avgGrad.n)*nmag
    !! vanishes identically on no-normal-flow walls.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel
    ! Local
    integer :: n,i,iEl,k,iVar
    real(prec) :: nhat(1:2),gn,g(1:2)

    select type(m => mymodel)
    class is(ESAtmo2D_t)
      do n = 1,bc%nBoundaries
        iEl = bc%elements(n)
        k = bc%sides(n)
        do i = 1,m%solution%interp%N+1
          nhat = m%geometry%nHat%boundary(i,k,iEl,1,1:2)
          do iVar = 1,m%nvar
            g(1:2) = m%solutionGradient%boundary(i,k,iEl,iVar,1:2)
            gn = g(1)*nhat(1)+g(2)*nhat(2)
            m%solutionGradient%extBoundary(i,k,iEl,iVar,1) = g(1)-2.0_prec*gn*nhat(1)
            m%solutionGradient%extBoundary(i,k,iEl,iVar,2) = g(2)-2.0_prec*gn*nhat(2)
          enddo
        enddo
      enddo
    endselect

  endsubroutine pbc2d_NoStress_ESAtmo2D

  subroutine AdditionalInit_ESAtmo2D_t(this)
    implicit none
    class(ESAtmo2D_t),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    bcfunc => hbc2d_NoNormalFlow_ESAtmo2D
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

    ! Parabolic BC for the same wall tag: zero diffusive normal flux
    ! (no-stress / no-heat-flux). hyperbolicBCs and parabolicBCs are
    ! independent linked lists, so registering the same tag here does
    ! not clobber the hyperbolic registration above.
    bcfunc => pbc2d_NoStress_ESAtmo2D
    call this%parabolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

    ! Diffusive-flux scratch buffers. Always allocated (memory cost is
    ! modest); only used when nu>0 or kappa>0 (and gradient_enabled is
    ! therefore .true.).
    call this%diffFlux%Init(this%solution%interp,this%nvar,this%mesh%nElem)
    call this%diffFlux%AssociateGeometry(this%geometry)
    call this%diffDiv%Init(this%solution%interp,this%nvar,this%mesh%nElem)

  endsubroutine AdditionalInit_ESAtmo2D_t

  subroutine AdditionalFree_ESAtmo2D_t(this)
    implicit none
    class(ESAtmo2D_t),intent(inout) :: this

    call this%diffFlux%Free()
    call this%diffDiv%Free()

  endsubroutine AdditionalFree_ESAtmo2D_t

  subroutine SetDiffusion_ESAtmo2D_t(this,nu,kappa,eta_penalty)
    !! Set the constant-coefficient Laplacian diffusion coefficients
    !! (kinematic momentum diffusivity and thermal diffusivity, both
    !! in m^2/s) and the dimensionless SIPG jump penalty. Setting nu
    !! or kappa > 0 enables the gradient pipeline so that the diffusive
    !! flux methods receive solutionGradient.
    !!
    !! length_scale is computed from the area Jacobian: for a quad
    !! reference cell [-1,1]^2 with area 4, the physical element area
    !! is 4*<J> so the characteristic edge length is 2*<J>^(1/2).
    !! The smallest length over the mesh is conservative.
    implicit none
    class(ESAtmo2D_t),intent(inout) :: this
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
    this%length_scale = 2.0_prec*sqrt(jmin)

  endsubroutine SetDiffusion_ESAtmo2D_t

  subroutine DiffusiveFluxMethod_ESAtmo2D_t(this)
    !! Fill diffFlux%interior with the constant-coefficient Laplacian
    !! flux at every interior node:
    !!
    !!   F_d(rho)      = 0                              (no mass diffusion)
    !!   F_d(rho*v_i)  = -nu    * d(rho*v_i)/dx_d
    !!   F_d(rho*theta)= -kappa * d(rho*theta)/dx_d
    !!   F_d(Phi)      = 0
    !!
    !! solutionGradient%interior(i,j,iel,iVar,d) holds d(s_iVar)/dx_d.
    implicit none
    class(ESAtmo2D_t),intent(inout) :: this
    ! Local
    integer :: i,j,iel,d

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem,d=1:2)
      this%diffFlux%interior(i,j,iel,1,d) = 0.0_prec
      this%diffFlux%interior(i,j,iel,2,d) = &
        -this%nu*this%solutionGradient%interior(i,j,iel,2,d)
      this%diffFlux%interior(i,j,iel,3,d) = &
        -this%nu*this%solutionGradient%interior(i,j,iel,3,d)
      this%diffFlux%interior(i,j,iel,4,d) = &
        -this%kappa*this%solutionGradient%interior(i,j,iel,4,d)
      ! Geopotential carries no diffusion.
      this%diffFlux%interior(i,j,iel,5,d) = 0.0_prec
    enddo

  endsubroutine DiffusiveFluxMethod_ESAtmo2D_t

  subroutine DiffusiveBoundaryFlux_ESAtmo2D_t(this)
    !! Fill diffFlux%boundaryNormal with the SIPG-stabilised BR1 flux:
    !!
    !!   f_R^diff(iVar) = -coeff(iVar) * (avg_grad . n) * nmag
    !!                    + tau(iVar) * (uL - uR) * nmag
    !!
    !! tau(iVar) = eta_penalty * coeff(iVar) * (N+1)^2 / length_scale.
    !! avg_grad is solutionGradient%avgBoundary (populated by
    !! AverageSides()); uL, uR are solution%boundary, %extBoundary.
    implicit none
    class(ESAtmo2D_t),intent(inout) :: this
    ! Local
    integer :: i,k,iel
    real(prec) :: nhat(1:2),nmag,gradn,coeff,tau,np2
    real(prec) :: uL,uR

    np2 = real((this%solution%interp%N+1)**2,prec)

    do concurrent(i=1:this%solution%N+1,k=1:4,iel=1:this%mesh%nElem)
      nhat = this%geometry%nHat%boundary(i,k,iel,1,1:2)
      nmag = this%geometry%nScale%boundary(i,k,iel,1)

      ! rho — no mass diffusion, no penalty
      this%diffFlux%boundaryNormal(i,k,iel,1) = 0.0_prec

      ! Momentum equations use nu
      coeff = this%nu
      tau = this%eta_penalty*coeff*np2/this%length_scale

      gradn = this%solutionGradient%avgBoundary(i,k,iel,2,1)*nhat(1)+ &
              this%solutionGradient%avgBoundary(i,k,iel,2,2)*nhat(2)
      uL = this%solution%boundary(i,k,iel,2)
      uR = this%solution%extBoundary(i,k,iel,2)
      this%diffFlux%boundaryNormal(i,k,iel,2) = (-coeff*gradn+tau*(uL-uR))*nmag

      gradn = this%solutionGradient%avgBoundary(i,k,iel,3,1)*nhat(1)+ &
              this%solutionGradient%avgBoundary(i,k,iel,3,2)*nhat(2)
      uL = this%solution%boundary(i,k,iel,3)
      uR = this%solution%extBoundary(i,k,iel,3)
      this%diffFlux%boundaryNormal(i,k,iel,3) = (-coeff*gradn+tau*(uL-uR))*nmag

      ! rho*theta uses kappa
      coeff = this%kappa
      tau = this%eta_penalty*coeff*np2/this%length_scale

      gradn = this%solutionGradient%avgBoundary(i,k,iel,4,1)*nhat(1)+ &
              this%solutionGradient%avgBoundary(i,k,iel,4,2)*nhat(2)
      uL = this%solution%boundary(i,k,iel,4)
      uR = this%solution%extBoundary(i,k,iel,4)
      this%diffFlux%boundaryNormal(i,k,iel,4) = (-coeff*gradn+tau*(uL-uR))*nmag

      ! Geopotential — no diffusion, no penalty.
      this%diffFlux%boundaryNormal(i,k,iel,5) = 0.0_prec
    enddo

  endsubroutine DiffusiveBoundaryFlux_ESAtmo2D_t

  subroutine CalculateTendency_ESAtmo2D_t(this)
    !! ESAtmo2D tendency = EC inviscid pipeline (parent) + optional
    !! constant-coefficient Laplacian diffusion (BR1 weak-form DG).
    !!
    !! When nu>0 or kappa>0 the parent's CalculateTendency already
    !! computes solutionGradient (because gradient_enabled was set by
    !! SetDiffusion); we then fill diffFlux from solutionGradient,
    !! compute its DG divergence, and accumulate into fluxDivergence
    !! before forming dSdt.
    implicit none
    class(ESAtmo2D_t),intent(inout) :: this
    ! Local
    integer :: i,j,iEl,iVar
    real(prec) :: bMi1,bMi2,bMj1,bMj2,qwi,qwj,jac

    call this%solution%BoundaryInterp()
    call this%solution%SideExchange(this%mesh)

    call this%PreTendency()
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
    ! Side ordering: 1=South, 2=East, 3=North, 4=West
    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iEl=1:this%mesh%nElem,iVar=1:this%solution%nVar)

      bMi1 = this%solution%interp%bMatrix(i,2) ! East boundary value at node i
      bMi2 = this%solution%interp%bMatrix(i,1) ! West boundary value at node i
      bMj1 = this%solution%interp%bMatrix(j,2) ! North boundary value at node j
      bMj2 = this%solution%interp%bMatrix(j,1) ! South boundary value at node j
      qwi = this%solution%interp%qWeights(i)
      qwj = this%solution%interp%qWeights(j)
      jac = this%geometry%J%interior(i,j,iEl,1)

      this%fluxDivergence%interior(i,j,iEl,iVar) = &
        this%fluxDivergence%interior(i,j,iEl,iVar)+ &
        (bMi1*this%flux%boundaryNormal(j,2,iEl,iVar)+ &
         bMi2*this%flux%boundaryNormal(j,4,iEl,iVar))/(qwi*jac)+ &
        (bMj1*this%flux%boundaryNormal(i,3,iEl,iVar)+ &
         bMj2*this%flux%boundaryNormal(i,1,iEl,iVar))/(qwj*jac)
    enddo

    ! Add the diffusive contribution if any coefficient is nonzero.
    if(this%nu > 0.0_prec .or. this%kappa > 0.0_prec) then
      call this%DiffusiveFluxMethod()
      call this%DiffusiveBoundaryFlux()
      call this%diffFlux%MappedDGDivergence(this%diffDiv%interior)
      do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                    iEl=1:this%mesh%nElem,iVar=1:this%solution%nVar)
        this%fluxDivergence%interior(i,j,iEl,iVar) = &
          this%fluxDivergence%interior(i,j,iEl,iVar)+ &
          this%diffDiv%interior(i,j,iEl,iVar)
      enddo
    endif

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iEl=1:this%mesh%nElem,iVar=1:this%solution%nVar)

      this%dSdt%interior(i,j,iEl,iVar) = &
        this%source%interior(i,j,iEl,iVar)- &
        this%fluxDivergence%interior(i,j,iEl,iVar)
    enddo

  endsubroutine CalculateTendency_ESAtmo2D_t

  subroutine hbc2d_NoNormalFlow_ESAtmo2D(bc,mymodel)
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
    integer :: n,i,iEl,k
    real(prec) :: nhat(1:2),rhovn

    select type(m => mymodel)
    class is(ESAtmo2D_t)
      do n = 1,bc%nBoundaries
        iEl = bc%elements(n)
        k = bc%sides(n)
        do i = 1,m%solution%interp%N+1

          nhat = m%geometry%nHat%boundary(i,k,iEl,1,1:2)

          ! Mirror density (tracer)
          m%solution%extBoundary(i,k,iEl,1) = &
            m%solution%boundary(i,k,iEl,1)

          ! Reflect momentum: negate normal component, preserve tangential
          ! (rho*v)_ext = (rho*v)_int - 2*((rho*v)_int . n)*n
          rhovn = m%solution%boundary(i,k,iEl,2)*nhat(1)+ &
                  m%solution%boundary(i,k,iEl,3)*nhat(2)

          m%solution%extBoundary(i,k,iEl,2) = &
            m%solution%boundary(i,k,iEl,2)-2.0_prec*rhovn*nhat(1)
          m%solution%extBoundary(i,k,iEl,3) = &
            m%solution%boundary(i,k,iEl,3)-2.0_prec*rhovn*nhat(2)

          ! Mirror potential temperature (tracer)
          m%solution%extBoundary(i,k,iEl,4) = &
            m%solution%boundary(i,k,iEl,4)

          ! Mirror geopotential — Phi only depends on y, so the
          ! mirror across a no-normal-flow wall is the same value.
          m%solution%extBoundary(i,k,iEl,5) = &
            m%solution%boundary(i,k,iEl,5)

        enddo
      enddo
    endselect

  endsubroutine hbc2d_NoNormalFlow_ESAtmo2D

  subroutine SetHydrostaticBalance_ESAtmo2D_t(this,theta0)
    !! Initialise a hydrostatically balanced atmosphere with uniform
    !! potential temperature theta0, zero velocity, and the geopotential
    !! Phi = g*y carried as state variable index 5.
    !!
    !! Exner function: pi(y) = 1 - g*y / (cp*theta0)
    !! rho(y)         = p0/(Rd*theta0) * pi(y)^(cv/Rd)
    !! rho*theta(y)   = p0/Rd * pi(y)^(cv/Rd)
    !!
    !! Boundary and extBoundary buffers of the solution are populated
    !! analytically from per-face y so that face values are bit-exact
    !! with what BoundaryInterp would produce (and extBoundary at walls
    !! mirrors). SideExchange will overwrite extBoundary at interior
    !! element interfaces with the neighbour's value — for a smooth
    !! profile this is a no-op.
    implicit none
    class(ESAtmo2D_t),intent(inout) :: this
    real(prec),intent(in) :: theta0
    ! Local
    integer :: i,j,iEl
    integer :: side
    real(prec) :: y,exner,rho

    print*,__FILE__," : Setting hydrostatic balance with theta0 = ",theta0

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iEl=1:this%mesh%nElem)

      y = this%geometry%x%interior(i,j,iEl,1,2)
      exner = 1.0_prec-this%g*y/(this%cp*theta0)
      rho = this%p0/(this%Rd*theta0)*exner**(this%cv/this%Rd)

      this%solution%interior(i,j,iEl,1) = rho
      this%solution%interior(i,j,iEl,2) = 0.0_prec
      this%solution%interior(i,j,iEl,3) = 0.0_prec
      this%solution%interior(i,j,iEl,4) = rho*theta0
      this%solution%interior(i,j,iEl,5) = this%g*y

    enddo

    ! Populate boundary + extBoundary buffers analytically from the
    ! per-face y-coordinates. Geopotential at the boundary mirrors
    ! by construction (Phi = g*y is single-valued in y).
    do concurrent(i=1:this%solution%N+1,side=1:4,iEl=1:this%mesh%nElem)

      y = this%geometry%x%boundary(i,side,iEl,1,2)
      exner = 1.0_prec-this%g*y/(this%cp*theta0)
      rho = this%p0/(this%Rd*theta0)*exner**(this%cv/this%Rd)

      this%solution%boundary(i,side,iEl,1) = rho
      this%solution%boundary(i,side,iEl,2) = 0.0_prec
      this%solution%boundary(i,side,iEl,3) = 0.0_prec
      this%solution%boundary(i,side,iEl,4) = rho*theta0
      this%solution%boundary(i,side,iEl,5) = this%g*y

      this%solution%extBoundary(i,side,iEl,1) = rho
      this%solution%extBoundary(i,side,iEl,2) = 0.0_prec
      this%solution%extBoundary(i,side,iEl,3) = 0.0_prec
      this%solution%extBoundary(i,side,iEl,4) = rho*theta0
      this%solution%extBoundary(i,side,iEl,5) = this%g*y

    enddo

    call this%ReportMetrics()
    call this%solution%UpdateDevice()

  endsubroutine SetHydrostaticBalance_ESAtmo2D_t

  subroutine AddThermalBubble_ESAtmo2D_t(this,dtheta,r0,x0,y0)
    !! Adds a pressure-balanced warm bubble perturbation.
    !!
    !! The potential temperature perturbation has a cos^2 profile:
    !!
    !!   theta'(x,y) = dtheta * cos^2(pi*r / (2*r0))  for r <= r0
    !!   theta'(x,y) = 0                                for r >  r0
    !!
    !! where r = sqrt((x-x0)^2 + (y-y0)^2).
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
    class(ESAtmo2D_t),intent(inout) :: this
    real(prec),intent(in) :: dtheta ! Perturbation amplitude [K]
    real(prec),intent(in) :: r0 ! Bubble radius [m]
    real(prec),intent(in) :: x0,y0 ! Bubble center [m]
    ! Local
    integer :: i,j,iEl
    real(prec) :: x,y,r,rho_old,theta_old,theta_new,thetap

    print*,__FILE__," : Adding thermal bubble perturbation"
    print*,__FILE__," : dtheta = ",dtheta
    print*,__FILE__," : r0     = ",r0
    print*,__FILE__," : center = ",x0,y0

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iEl=1:this%mesh%nElem)

      x = this%geometry%x%interior(i,j,iEl,1,1)-x0
      y = this%geometry%x%interior(i,j,iEl,1,2)-y0
      r = sqrt(x*x+y*y)

      if(r <= r0) then
        rho_old = this%solution%interior(i,j,iEl,1)
        theta_old = this%solution%interior(i,j,iEl,4)/rho_old
        thetap = dtheta*cos(pi*r/(2.0_prec*r0))**2
        theta_new = theta_old+thetap

        ! Adjust density to maintain pressure balance: rho*theta = const
        this%solution%interior(i,j,iEl,1) = rho_old*theta_old/theta_new
        ! rho*theta is unchanged (pressure-balanced perturbation)
      endif

    enddo

    call this%ReportMetrics()
    call this%solution%UpdateDevice()

  endsubroutine AddThermalBubble_ESAtmo2D_t

endmodule SELF_ESAtmo2D_t
