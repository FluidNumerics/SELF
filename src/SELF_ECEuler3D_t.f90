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

module SELF_ECEuler3D_t
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
  use SELF_mesh
  use SELF_BoundaryConditions

  implicit none

  type,extends(ECDGModel3D) :: ECEuler3D_t

    real(prec) :: p0 = 100000.0_prec ! Reference pressure [Pa]
    real(prec) :: Rd = 287.0_prec ! Gas constant for dry air [J/(kg*K)]
    real(prec) :: cp = 1004.0_prec ! Specific heat at constant pressure [J/(kg*K)]
    real(prec) :: cv = 717.0_prec ! Specific heat at constant volume [J/(kg*K)]
    real(prec) :: g = 9.81_prec ! Gravitational acceleration [m/s^2]

  contains

    procedure :: SetNumberOfVariables => SetNumberOfVariables_ECEuler3D_t
    procedure :: SetMetadata => SetMetadata_ECEuler3D_t
    procedure :: entropy_func => entropy_func_ECEuler3D_t
    procedure :: twopointflux3d => twopointflux3d_ECEuler3D_t
    procedure :: riemannflux3d => riemannflux3d_ECEuler3D_t
    procedure :: SourceMethod => SourceMethod_ECEuler3D_t
    procedure :: AdditionalInit => AdditionalInit_ECEuler3D_t
    procedure :: SetHydrostaticBalance => SetHydrostaticBalance_ECEuler3D_t
    procedure :: AddThermalBubble => AddThermalBubble_ECEuler3D_t

  endtype ECEuler3D_t

contains

  subroutine SetNumberOfVariables_ECEuler3D_t(this)
    implicit none
    class(ECEuler3D_t),intent(inout) :: this

    this%nvar = 5

  endsubroutine SetNumberOfVariables_ECEuler3D_t

  subroutine SetMetadata_ECEuler3D_t(this)
    implicit none
    class(ECEuler3D_t),intent(inout) :: this

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

  endsubroutine SetMetadata_ECEuler3D_t

  pure function entropy_func_ECEuler3D_t(this,s) result(e)
    !! Mathematical entropy: total energy density (kinetic + internal).
    !!
    !!   e = 0.5*(rhou^2 + rhov^2 + rhow^2)/rho + p/(gamma - 1)
    !!
    !! where p = p0 * (rho * Rd * theta / p0)^gamma and gamma = cp/cv.
    class(ECEuler3D_t),intent(in) :: this
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

  endfunction entropy_func_ECEuler3D_t

  pure function twopointflux3d_ECEuler3D_t(this,sL,sR) result(flux)
    !! Kennedy-Gruber split-form two-point flux (kinetic energy preserving).
    !!
    !! Uses arithmetic means of primitive variables:
    !!
    !!   f_d(rho)    = {{rho}} * {{v_d}}
    !!   f_d(rho*vi) = {{rho}} * {{vi}} * {{v_d}} + {{p}} * delta_{id}
    !!   f_d(rho*th) = {{rho}} * {{theta}} * {{v_d}}
    !!
    !! where {{.}} = arithmetic mean.
    class(ECEuler3D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec) :: flux(1:this%nvar,1:3)
    ! Local
    real(prec) :: rhoL,uL,vL,wL,thetaL,pL
    real(prec) :: rhoR,uR,vR,wR,thetaR,pR
    real(prec) :: rho_avg,u_avg,v_avg,w_avg,theta_avg,p_avg
    real(prec) :: gamma

    gamma = this%cp/this%cv

    ! Left primitive variables
    rhoL = sL(1)
    uL = sL(2)/rhoL
    vL = sL(3)/rhoL
    wL = sL(4)/rhoL
    thetaL = sL(5)/rhoL
    pL = this%p0*(rhoL*this%Rd*thetaL/this%p0)**gamma

    ! Right primitive variables
    rhoR = sR(1)
    uR = sR(2)/rhoR
    vR = sR(3)/rhoR
    wR = sR(4)/rhoR
    thetaR = sR(5)/rhoR
    pR = this%p0*(rhoR*this%Rd*thetaR/this%p0)**gamma

    ! Arithmetic averages
    rho_avg = 0.5_prec*(rhoL+rhoR)
    u_avg = 0.5_prec*(uL+uR)
    v_avg = 0.5_prec*(vL+vR)
    w_avg = 0.5_prec*(wL+wR)
    theta_avg = 0.5_prec*(thetaL+thetaR)
    p_avg = 0.5_prec*(pL+pR)

    ! x-direction flux (d=1)
    flux(1,1) = rho_avg*u_avg
    flux(2,1) = rho_avg*u_avg*u_avg+p_avg
    flux(3,1) = rho_avg*v_avg*u_avg
    flux(4,1) = rho_avg*w_avg*u_avg
    flux(5,1) = rho_avg*theta_avg*u_avg

    ! y-direction flux (d=2)
    flux(1,2) = rho_avg*v_avg
    flux(2,2) = rho_avg*u_avg*v_avg
    flux(3,2) = rho_avg*v_avg*v_avg+p_avg
    flux(4,2) = rho_avg*w_avg*v_avg
    flux(5,2) = rho_avg*theta_avg*v_avg

    ! z-direction flux (d=3)
    flux(1,3) = rho_avg*w_avg
    flux(2,3) = rho_avg*u_avg*w_avg
    flux(3,3) = rho_avg*v_avg*w_avg
    flux(4,3) = rho_avg*w_avg*w_avg+p_avg
    flux(5,3) = rho_avg*theta_avg*w_avg

  endfunction twopointflux3d_ECEuler3D_t

  pure function riemannflux3d_ECEuler3D_t(this,sL,sR,dsdx,nhat) result(flux)
    !! Local Lax-Friedrichs (Rusanov) Riemann flux.
    !!
    !!   F* = 0.5*(fL.n + fR.n) - 0.5*lambda_max*(sR - sL)
    !!
    !! where lambda_max = max(|vL.n| + cL, |vR.n| + cR)
    !! and c = sqrt(gamma * p / rho) is the sound speed.
    class(ECEuler3D_t),intent(in) :: this
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

    if(.false.) flux(1) = flux(1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction riemannflux3d_ECEuler3D_t

  subroutine SourceMethod_ECEuler3D_t(this)
    !! Gravitational source term: S = [0, 0, 0, -rho*g, 0]
    implicit none
    class(ECEuler3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,k,iEl

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%solution%N+1,iEl=1:this%mesh%nElem)

      this%source%interior(i,j,k,iEl,1) = 0.0_prec
      this%source%interior(i,j,k,iEl,2) = 0.0_prec
      this%source%interior(i,j,k,iEl,3) = 0.0_prec
      this%source%interior(i,j,k,iEl,4) = &
        -this%solution%interior(i,j,k,iEl,1)*this%g
      this%source%interior(i,j,k,iEl,5) = 0.0_prec

    enddo

  endsubroutine SourceMethod_ECEuler3D_t

  subroutine AdditionalInit_ECEuler3D_t(this)
    implicit none
    class(ECEuler3D_t),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    bcfunc => hbc3d_NoNormalFlow_ECEuler3D
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

  endsubroutine AdditionalInit_ECEuler3D_t

  subroutine hbc3d_NoNormalFlow_ECEuler3D(bc,mymodel)
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
    class is(ECEuler3D_t)
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

          enddo
        enddo
      enddo
    endselect

  endsubroutine hbc3d_NoNormalFlow_ECEuler3D

  subroutine SetHydrostaticBalance_ECEuler3D_t(this,theta0)
    !! Sets the solution to a hydrostatically balanced atmosphere with
    !! uniform potential temperature theta0 and zero velocity.
    !!
    !! The Exner function is:
    !!   pi(z) = 1 - g*z / (cp*theta0)
    !!
    !! From which:
    !!   rho(z)       = p0/(Rd*theta0) * pi(z)^(cv/Rd)
    !!   rho*theta(z) = p0/Rd * pi(z)^(cv/Rd)
    !!   rho*u = rho*v = rho*w = 0
    implicit none
    class(ECEuler3D_t),intent(inout) :: this
    real(prec),intent(in) :: theta0
    ! Local
    integer :: i,j,k,iEl
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

    enddo

    call this%ReportMetrics()
    call this%solution%UpdateDevice()

  endsubroutine SetHydrostaticBalance_ECEuler3D_t

  subroutine AddThermalBubble_ECEuler3D_t(this,dtheta,r0,x0,y0,z0)
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
    class(ECEuler3D_t),intent(inout) :: this
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

  endsubroutine AddThermalBubble_ECEuler3D_t

endmodule SELF_ECEuler3D_t
