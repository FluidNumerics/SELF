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
!! linearized about a motionless background state. The sound speed and the
!! background density are carried as solution variables (static in time,
!! possibly spatially varying) so that heterogeneous media are supported,
!! mirroring the 2-D model.
!!
!! The solution variables are
!!
!! \begin{equation}
!! \vec{s} = \begin{pmatrix}
!!     \rho \\
!!      u \\
!!      v \\
!!      w \\
!!      p \\
!!      c \\
!!      \rho_0
!!  \end{pmatrix}
!! \end{equation}
!!
!! The conservative flux is
!!
!! \begin{equation}
!! \overleftrightarrow{f} = \begin{pmatrix}
!!     \rho_0 u \hat{x} + \rho_0 v \hat{y} + \rho_0 w \hat{z} \\
!!      \frac{p}{\rho_0} \hat{x} \\
!!      \frac{p}{\rho_0} \hat{y} \\
!!      \frac{p}{\rho_0} \hat{z} \\
!!      c^2 \rho_0 ( u \hat{x} + v \hat{y} + w \hat{z}) \\
!!      0 \\
!!      0
!!  \end{pmatrix}
!! \end{equation}
!!
!! and the source terms are null. The sound speed and background density
!! variables have zero flux and zero source; they are set by the initial
!! condition and preserved in time. This is entropy-stable for
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
    real(prec) :: rho0 = 1.0_prec ! Reference density (used to fill variable 7 in initial conditions)
    real(prec) :: c = 1.0_prec ! Reference sound speed (used to fill variable 6 in initial conditions)
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
    ! Only the first five variables (rho, u, v, w, P) are advanced in time. The
    ! sixth and seventh variables, the sound speed c and background density
    ! rho0, are spatially-varying but time-constant background fields: their
    ! flux and source are identically zero, so they are excluded from time
    ! integration rather than stepped to a no-op.
    this%nstepped = 5

  endsubroutine SetNumberOfVariables_LinearEuler3D_t

  subroutine SetMetadata_LinearEuler3D_t(this)
    implicit none
    class(LinearEuler3D_t),intent(inout) :: this

    call this%solution%SetName(1,"rho") ! Density
    call this%solution%SetUnits(1,"kg⋅m⁻³")

    call this%solution%SetName(2,"u") ! x-velocity component
    call this%solution%SetUnits(2,"m⋅s⁻¹")

    call this%solution%SetName(3,"v") ! y-velocity component
    call this%solution%SetUnits(3,"m⋅s⁻¹")

    call this%solution%SetName(4,"w") ! z-velocity component
    call this%solution%SetUnits(4,"m⋅s⁻¹")

    call this%solution%SetName(5,"P") ! Pressure
    call this%solution%SetUnits(5,"kg⋅m⁻¹⋅s⁻²")

    call this%solution%SetName(6,"c") ! Sound speed (static; possibly heterogeneous)
    call this%solution%SetUnits(6,"m⋅s⁻¹")

    call this%solution%SetName(7,"rho0") ! Background density (static; possibly heterogeneous)
    call this%solution%SetUnits(7,"kg⋅m⁻³")

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
    !! sound speed (variable 6) and background density (variable 7) are copied
    !! from the interior side so the Riemann solver sees a consistent c and
    !! rho0.
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
            m%solution%extBoundary(i,j,s,iEl,1:5) = 0.0_prec
            m%solution%extBoundary(i,j,s,iEl,6) = m%solution%boundary(i,j,s,iEl,6) ! c preserved
            m%solution%extBoundary(i,j,s,iEl,7) = m%solution%boundary(i,j,s,iEl,7) ! rho0 preserved
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
    !! where the sound speed c is taken from s(6) and the background density
    !! rho0 from s(7).
    class(LinearEuler3D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: e

    e = 0.5_prec*s(7)*(s(2)*s(2)+s(3)*s(3)+s(4)*s(4))+ &
        0.5_prec*(s(5)*s(5)/(s(7)*s(6)*s(6)))

  endfunction entropy_func_LinearEuler3D_t

  subroutine sourcemethod_LinearEuler3D_t(this)
    implicit none
    class(LinearEuler3D_t),intent(inout) :: this

    if(.false.) this%nvar = this%nvar ! suppress unused-dummy-argument warning

  endsubroutine sourcemethod_LinearEuler3D_t

  pure function flux3D_LinearEuler3D_t(this,s,dsdx) result(flux)
    class(LinearEuler3D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec) :: flux(1:this%nvar,1:3)

    flux(1,1) = s(7)*s(2) ! density, x flux ; rho0*u
    flux(1,2) = s(7)*s(3) ! density, y flux ; rho0*v
    flux(1,3) = s(7)*s(4) ! density, y flux ; rho0*w

    flux(2,1) = s(5)/s(7) ! x-velocity, x flux; p/rho0
    flux(2,2) = 0.0_prec ! x-velocity, y flux; 0
    flux(2,3) = 0.0_prec ! x-velocity, z flux; 0

    flux(3,1) = 0.0_prec ! y-velocity, x flux; 0
    flux(3,2) = s(5)/s(7) ! y-velocity, y flux; p/rho0
    flux(3,3) = 0.0_prec ! y-velocity, z flux; 0

    flux(4,1) = 0.0_prec ! z-velocity, x flux; 0
    flux(4,2) = 0.0_prec ! z-velocity, y flux; 0
    flux(4,3) = s(5)/s(7) ! z-velocity, z flux; p/rho0

    flux(5,1) = s(6)*s(6)*s(7)*s(2) ! pressure, x flux : rho0*c^2*u
    flux(5,2) = s(6)*s(6)*s(7)*s(3) ! pressure, y flux : rho0*c^2*v
    flux(5,3) = s(6)*s(6)*s(7)*s(4) ! pressure, z flux : rho0*c^2*w

    flux(6,1) = 0.0_prec ! sound speed, x flux; 0 (c held fixed in time)
    flux(6,2) = 0.0_prec ! sound speed, y flux; 0
    flux(6,3) = 0.0_prec ! sound speed, z flux; 0

    flux(7,1) = 0.0_prec ! background density, x flux; 0 (rho0 held fixed in time)
    flux(7,2) = 0.0_prec ! background density, y flux; 0
    flux(7,3) = 0.0_prec ! background density, z flux; 0
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
    !! The reconstructed density/momentum/pressure fluxes use the arithmetic
    !! averages rho0_avg and c2_avg at the face (see the 2-D model for the
    !! rationale). The sound speed and background density variables carry zero
    !! flux.
    class(LinearEuler3D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec),intent(in) :: nhat(1:3)
    real(prec) :: flux(1:this%nvar)
    ! Local
    real(prec) :: rho0L,rho0R,rho0_avg,cL,cR,ZL,ZR,unL,unR,pL,pR,un_star,p_star,c2_avg

    rho0L = sL(7)
    rho0R = sR(7)
    rho0_avg = 0.5_prec*(rho0L+rho0R)
    cL = sL(6)
    cR = sR(6)
    ZL = rho0L*cL
    ZR = rho0R*cR

    unL = sL(2)*nhat(1)+sL(3)*nhat(2)+sL(4)*nhat(3)
    unR = sR(2)*nhat(1)+sR(3)*nhat(2)+sR(4)*nhat(3)
    pL = sL(5)
    pR = sR(5)

    un_star = (ZL*unL+ZR*unR+(pL-pR))/(ZL+ZR)
    p_star = (ZR*pL+ZL*pR+ZL*ZR*(unL-unR))/(ZL+ZR)
    c2_avg = 0.5_prec*(cL*cL+cR*cR)

    flux(1) = rho0_avg*un_star
    flux(2) = p_star*nhat(1)/rho0_avg
    flux(3) = p_star*nhat(2)/rho0_avg
    flux(4) = p_star*nhat(3)/rho0_avg
    flux(5) = rho0_avg*c2_avg*un_star
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
    !! \rho &= \rho_0 + \rho' \exp\left( -\ln(2) \frac{(x-x_0)^2 + (y-y_0)^2}{L_r^2} \right)
    !! u &= 0 \\
    !! v &= 0 \\
    !! E &= \frac{P_0}{\gamma - 1} + E \exp\left( -\ln(2) \frac{(x-x_0)^2 + (y-y_0)^2}{L_e^2} \right)
    !! \end{aligned}
    !! \end{equation}
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

      this%solution%interior(i,j,k,iEl,1) = rho
      this%solution%interior(i,j,k,iEl,2) = 0.0_prec
      this%solution%interior(i,j,k,iEl,3) = 0.0_prec
      this%solution%interior(i,j,k,iEl,4) = 0.0_prec
      this%solution%interior(i,j,k,iEl,5) = rho*this%c*this%c
      this%solution%interior(i,j,k,iEl,6) = this%c ! uniform background sound speed
      this%solution%interior(i,j,k,iEl,7) = this%rho0 ! uniform background density

    enddo

    call this%ReportMetrics()
    call this%solution%UpdateDevice()

  endsubroutine SphericalSoundWave_LinearEuler3D_t

endmodule self_LinearEuler3D_t
