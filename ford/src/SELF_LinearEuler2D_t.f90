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
!!      p \\
!!      c \\
!!      \rho_0
!!  \end{pmatrix}
!! \end{equation}
!!
!! The sound speed \(c\) and the background density \(\rho_0\) are carried as
!! solution variables so that they can vary in space (heterogeneous media).
!! Their flux and source are identically zero, so they are held fixed in time
!! at each spatial location. This is entropy-stable for piecewise-constant
!! material regions aligned with element boundaries: interiors have
!! \(\nabla \rho_0 = \nabla c = 0\) (so the flux-divergence form is exact) and
!! the impedance-matched Riemann flux handles the jumps at faces.
!!
!! The conservative flux is
!!
!! \begin{equation}
!! \overleftrightarrow{f} = \begin{pmatrix}
!!      \frac{p}{\rho_0} \hat{x} \\
!!      \frac{p}{\rho_0} \hat{y} \\
!!      c^2 \rho_0 ( u \hat{x} + v \hat{y} ) \\
!!      \vec{0} \\
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
    real(prec) :: rho0 = 1.0_prec ! Reference density (used to fill variable 5 in initial conditions)
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
    ! Only the first three variables (u, v, P) are advanced in time. The
    ! fourth and fifth variables, the sound speed c and background density rho0,
    ! are spatially-varying but time-constant background fields: their flux and
    ! source are identically zero, so they are excluded from time integration
    ! rather than stepped to a no-op.
    this%nstepped = 3

  endsubroutine SetNumberOfVariables_LinearEuler2D_t

  subroutine SetMetadata_LinearEuler2D_t(this)
    implicit none
    class(LinearEuler2D_t),intent(inout) :: this

    call this%solution%SetName(1,"u") ! x-velocity component
    call this%solution%SetUnits(1,"m⋅s⁻¹")

    call this%solution%SetName(2,"v") ! y-velocity component
    call this%solution%SetUnits(2,"m⋅s⁻¹")

    call this%solution%SetName(3,"P") ! Pressure
    call this%solution%SetUnits(3,"kg⋅m⁻¹⋅s⁻²")

    call this%solution%SetName(4,"c") ! Sound speed
    call this%solution%SetUnits(4,"m⋅s⁻¹")

    call this%solution%SetName(5,"rho0") ! Background density (static; possibly heterogeneous)
    call this%solution%SetUnits(5,"kg⋅m⁻³")

  endsubroutine SetMetadata_LinearEuler2D_t

  pure function entropy_func_LinearEuler2D_t(this,s) result(e)
    !! The entropy function is the sum of kinetic and internal energy
    !! For the linear model, this is
    !!
    !! \begin{equation}
    !!   e = \frac{1}{2} \left( \rho_0*( u^2 + v^2 ) + \frac{P^2}{\rho_0 c^2} \right)
    !!
    !! where the sound speed c is taken from s(4) and the background density
    !! rho0 from s(5).
    class(LinearEuler2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: e

    e = 0.5_prec*s(5)*(s(1)*s(1)+s(2)*s(2))+ &
        0.5_prec*(s(3)*s(3)/(s(5)*s(4)*s(4)))

  endfunction entropy_func_LinearEuler2D_t

  subroutine AdditionalInit_LinearEuler2D_t(this)
    !! Register the (CPU) no-normal-flow and radiation boundary conditions.
    !! GPU builds call this parent and then overwrite both registrations with
    !! the device kernels in AdditionalInit_LinearEuler2D.
    implicit none
    class(LinearEuler2D_t),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    bcfunc => hbc2d_NoNormalFlow_LinearEuler2D
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

    bcfunc => hbc2d_Radiation_LinearEuler2D
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_RADIATION,"radiation",bcfunc)

  endsubroutine AdditionalInit_LinearEuler2D_t

  subroutine hbc2d_Radiation_LinearEuler2D(bc,mymodel)
    !! Radiation BC: zero acoustic perturbation in the exterior state; the
    !! sound speed (variable 4) and background density (variable 5) are copied
    !! from the interior side so the Riemann solver sees a consistent c and
    !! rho0 (impedance-matched, non-reflecting outflow).
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel
    ! Local
    integer :: n,i,iEl,j

    select type(m => mymodel)
    class is(LinearEuler2D_t)
      do n = 1,bc%nBoundaries
        iEl = bc%elements(n)
        j = bc%sides(n)
        do i = 1,m%solution%interp%N+1
          m%solution%extBoundary(i,j,iEl,1:3) = 0.0_prec
          m%solution%extBoundary(i,j,iEl,4) = m%solution%boundary(i,j,iEl,4) ! c preserved
          m%solution%extBoundary(i,j,iEl,5) = m%solution%boundary(i,j,iEl,5) ! rho0 preserved
        enddo
      enddo
    endselect

  endsubroutine hbc2d_Radiation_LinearEuler2D

  subroutine hbc2d_NoNormalFlow_LinearEuler2D(bc,mymodel)
    !! No-normal-flow boundary condition for 2D linear Euler equations.
    !! Reflects the velocity vector about the boundary normal while
    !! preserving pressure, sound speed, and background density.
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
          m%solution%extBoundary(i,j,iEl,1) = &
            (nhat(2)**2-nhat(1)**2)*s(1)-2.0_prec*nhat(1)*nhat(2)*s(2) ! u
          m%solution%extBoundary(i,j,iEl,2) = &
            (nhat(1)**2-nhat(2)**2)*s(2)-2.0_prec*nhat(1)*nhat(2)*s(1) ! v
          m%solution%extBoundary(i,j,iEl,3) = s(3) ! p
          m%solution%extBoundary(i,j,iEl,4) = s(4) ! c
          m%solution%extBoundary(i,j,iEl,5) = s(5) ! rho0
        enddo
      enddo
    endselect

  endsubroutine hbc2d_NoNormalFlow_LinearEuler2D

  pure function flux2d_LinearEuler2D_t(this,s,dsdx) result(flux)
    class(LinearEuler2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec) :: flux(1:this%nvar,1:2)

    flux(1,1) = s(3)/s(5) ! x-velocity, x flux; p/rho0
    flux(1,2) = 0.0_prec ! x-velocity, y flux; 0
    flux(2,1) = 0.0_prec ! y-velocity, x flux; 0
    flux(2,2) = s(3)/s(5) ! y-velocity, y flux; p/rho0
    flux(3,1) = s(4)*s(4)*s(5)*s(1) ! pressure, x flux : rho0*c^2*u
    flux(3,2) = s(4)*s(4)*s(5)*s(2) ! pressure, y flux : rho0*c^2*v
    flux(4,1) = 0.0_prec ! sound speed, x flux; 0 (c is held fixed in time)
    flux(4,2) = 0.0_prec ! sound speed, y flux; 0 (c is held fixed in time)
    flux(5,1) = 0.0_prec ! background density, x flux; 0 (rho0 is held fixed in time)
    flux(5,2) = 0.0_prec ! background density, y flux; 0 (rho0 is held fixed in time)
    if(.false.) flux(1,1) = flux(1,1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction flux2d_LinearEuler2D_t

  pure function riemannflux2d_LinearEuler2D_t(this,sL,sR,dsdx,nhat) result(flux)
    !! Characteristic-decomposition (impedance-matched) interface flux for
    !! linear acoustics with possibly discontinuous sound speed.
    !!
    !! The normal-flux Jacobian has eigenstructure
    !!   +c : right-going acoustic mode, W_+ = rho0*c*u_n + p
    !!   -c : left-going  acoustic mode, W_- = -rho0*c*u_n + p
    !!    0 : tangential vorticity mode , u_t
    !! Upwinding each mode by its characteristic direction at the face,
    !!   W_+|* = W_+|_L   (transmitted from left,  at speed +c_L)
    !!   W_-|* = W_-|_R   (transmitted from right, at speed -c_R)
    !! yields the impedance-matched interface state
    !!   u_n* = (Z_L u_n,L + Z_R u_n,R + (p_L - p_R)) / (Z_L + Z_R)
    !!   p*   = (Z_R p_L + Z_L p_R + Z_L Z_R (u_n,L - u_n,R)) / (Z_L + Z_R)
    !! with the per-side acoustic impedance Z = rho0*c (each side using its own
    !! background density rho0 and sound speed c). This is exact upwind /
    !! Godunov for the linearised acoustic system and reduces correctly to
    !! Fresnel reflection / transmission across an impedance jump (Z_L .ne. Z_R,
    !! whether from a density jump, a sound-speed jump, or both). LLF with
    !! cmax = max(c_L, c_R) over-dissipates the tangential mode
    !! and at high polynomial order fails to stably handle the
    !! impedance mismatch (aliasing instability at material interfaces).
    !!
    !! The reconstructed momentum/pressure fluxes need a single rho0
    !! (and c^2) at the face, but rho0 and c are two-valued across a material
    !! interface. We use the arithmetic averages rho0_avg and c2_avg; this is a
    !! pragmatic treatment of the non-conservative products p/rho0 and
    !! rho0*c^2*div(v) at a face where the coefficients jump. A fully
    !! path-conservative (Castro-Pares) treatment would use each side's own
    !! coefficients in its surface integral. For piecewise-constant material
    !! regions the interior is exactly entropy-conservative and the impedance
    !! solver above carries the (entropy-stable) interface dissipation.
    class(LinearEuler2D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: flux(1:this%nvar)
    ! Local
    real(prec) :: rho0L,rho0R,rho0_avg,cL,cR,ZL,ZR,unL,unR,pL,pR,un_star,p_star,c2_avg

    rho0L = sL(5)
    rho0R = sR(5)
    rho0_avg = 0.5_prec*(rho0L+rho0R)
    cL = sL(4)
    cR = sR(4)
    ZL = rho0L*cL
    ZR = rho0R*cR

    unL = sL(1)*nhat(1)+sL(2)*nhat(2)
    unR = sR(1)*nhat(1)+sR(2)*nhat(2)
    pL = sL(3)
    pR = sR(3)

    un_star = (ZL*unL+ZR*unR+(pL-pR))/(ZL+ZR)
    p_star = (ZR*pL+ZL*pR+ZL*ZR*(unL-unR))/(ZL+ZR)
    c2_avg = 0.5_prec*(cL*cL+cR*cR)

    flux(1) = p_star*nhat(1)/rho0_avg
    flux(2) = p_star*nhat(2)/rho0_avg
    flux(3) = rho0_avg*c2_avg*un_star
    flux(4) = 0.0_prec
    flux(5) = 0.0_prec
    if(.false.) flux(1) = flux(1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction riemannflux2d_LinearEuler2D_t

  subroutine SphericalSoundWave_LinearEuler2D_t(this,rhoprime,Lr,x0,y0,c)
    !! This subroutine sets the initial condition for a weak blast wave
    !! problem. The initial condition is given by
    !!
    !! \begin{equation}
    !! \begin{aligned}
    !! u &= 0 \\
    !! v &= 0 \\
    !! p &= \rho' c^2 \exp\left( -\ln(2) \frac{(x-x_0)^2 + (y-y_0)^2}{L_r^2} \right)
    !! \end{aligned}
    !! \end{equation}
    !!
    !! `rhoprime` is the amplitude of the (diagnostic) density anomaly that the
    !! pressure pulse corresponds to through the acoustic relation p = rho*c^2;
    !! the density anomaly itself is not a solution variable.
    !! The sound speed `c` (passed as an argument since it is not a
    !! scalar model attribute) and the background density `this%rho0` are set
    !! uniformly across the domain.
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

      this%solution%interior(i,j,iEl,1) = 0.0_prec
      this%solution%interior(i,j,iEl,2) = 0.0_prec
      this%solution%interior(i,j,iEl,3) = rho*c*c
      this%solution%interior(i,j,iEl,4) = c
      this%solution%interior(i,j,iEl,5) = this%rho0 ! uniform background density

    enddo

    call this%ReportMetrics()
    call this%solution%UpdateDevice()

  endsubroutine SphericalSoundWave_LinearEuler2D_t

endmodule self_LinearEuler2D_t
