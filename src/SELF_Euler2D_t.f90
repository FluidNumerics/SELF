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

module self_Euler2D_t
!! This module defines a class that can be used to solve the  Euler
!! equations in 2-D. The  Euler Equations, here, are the Euler equations
!! ized about a motionless background state.
!!
!! The conserved variables are

!! \begin{equation}
!! \vec{s} = \begin{pmatrix}
!!     \rho \\
!!     \rho u \\
!!     \rho v \\
!!     \rho E
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

  type,extends(dgmodel2d) :: Euler2D_t
    ! Add any additional attributes here that are specific to your model
    type(MappedScalar2D)   :: primitive
    type(MappedVector2D)   :: primitiveGradient
    type(MappedScalar2D)   :: diagnostics
    integer :: ndiagnostics = 2

    logical :: primitive_gradient_enabled = .false.
    real(prec) :: Cp = 1.4_prec ! Specific heat at constant pressure
    real(prec) :: Cv = 1.0_prec ! Specific heat at constant volume
    real(prec) :: gamma = 1.4_prec ! Ratio of specific heats (Cp/Cv)
    real(prec) :: nu = 0.0_prec ! Dynamic viscosity
    real(prec) :: kappa = 0.0_prec ! Thermal diffusivity
    real(prec) :: g = 0.0_prec ! gravitational acceleration (y-direction only)

  contains

    ! Setup / Book-keeping methods
    procedure :: AdditionalInit => AdditionalInit_Euler2D_t
    procedure :: AdditionalFree => AdditionalFree_Euler2D_t
    procedure :: SetNumberOfVariables => SetNumberOfVariables_Euler2D_t
    procedure :: SetMetadata => SetMetadata_Euler2D_t

    ! Pre-tendency methods
    procedure :: CalculateDiagnostics => CalculateDiagnostics_Euler2D_t
    procedure :: ConservativeToPrimitive => ConservativeToPrimitive_Euler2D_t
    
    ! Model definition methods
    procedure :: entropy_func => entropy_func_Euler2D_t
    
    procedure :: hbc2d_NoNormalFlow => hbc2d_NoNormalFlow_Euler2D_t
   ! procedure :: hbc2d_Radiation => hbc2d_Radiation_Euler2D_t

    procedure :: flux2d => flux2d_Euler2D_t
    procedure :: riemannflux2d => riemannflux2d_Euler2D_t
    !procedure :: source2d => source2d_Euler2D_t

    ! Additional support methods
    procedure :: PrimitiveToConservative => PrimitiveToConservative_Euler2D_t
    procedure,private :: pressure
    procedure,private :: speedofsound

    ! Example initial conditions
    procedure :: WeakBlastWave => WeakBlastWave_Euler2D_t

    
  endtype Euler2D_t

contains

 subroutine AdditionalInit_Euler2D_t(this)
    implicit none
    class(Euler2D_t),intent(inout) :: this

    call this%primitive%Init(this%geometry%x%interp,this%nvar,this%mesh%nElem)
    call this%primitiveGradient%Init(this%geometry%x%interp,this%nvar,this%mesh%nElem)
    call this%diagnostics%Init(this%geometry%x%interp,this%ndiagnostics,this%mesh%nElem)

  endsubroutine AdditionalInit_Euler2D_t

  subroutine AdditionalFree_Euler2D_t(this)
    implicit none
    class(Euler2D_t),intent(inout) :: this

    call this%primitive%Free()
    call this%primitiveGradient%Free()
    call this%diagnostics%Free()

  endsubroutine AdditionalFree_Euler2D_t

 subroutine SetNumberOfVariables_Euler2D_t(this)
    implicit none
    class(Euler2D_t),intent(inout) :: this

    this%nvar = 4

  endsubroutine SetNumberOfVariables_Euler2D_t

  subroutine SetMetadata_Euler2D_t(this)
    implicit none
    class(Euler2D_t),intent(inout) :: this

    call this%solution%SetName(1,"ρ") ! Density
    call this%solution%SetUnits(1,"kg⋅m⁻³")

    call this%solution%SetName(2,"ρu") ! z-momentum
    call this%solution%SetUnits(2,"(kg⋅m⁻³)(m⋅s⁻¹)")

    call this%solution%SetName(3,"ρv") ! y-momentum
    call this%solution%SetUnits(3,"(kg⋅m⁻³)(m⋅s⁻¹)")

    call this%solution%SetName(4,"ρE") ! Total energy
    call this%solution%SetUnits(4,"(kg⋅m⁻³)(m²⋅s⁻²)")

    call this%primitive%SetName(1,"ρ") ! Density
    call this%primitive%SetUnits(1,"kg⋅m⁻³")

    call this%primitive%SetName(2,"ρu") ! z-momentum
    call this%primitive%SetUnits(2,"(kg⋅m⁻³)(m⋅s⁻¹)")

    call this%primitive%SetName(3,"ρv") ! y-momentum
    call this%primitive%SetUnits(3,"(kg⋅m⁻³)(m⋅s⁻¹)")

    call this%primitive%SetName(4,"P") ! Pressure
    call this%primitive%SetUnits(4,"kg⋅m⁻¹⋅s⁻²")

    call this%diagnostics%SetName(1,"c") ! Speed of sound
    call this%diagnostics%SetUnits(1,"m⋅s⁻¹")

    call this%diagnostics%SetName(2,"T") ! in-situ temperature
    call this%diagnostics%SetUnits(2,"K")

  endsubroutine SetMetadata_Euler2D_t

  subroutine CalculateDiagnostics_Euler2D_t(this)
    implicit none
    class(Euler2D_t),intent(inout) :: this
    ! Local
    integer :: i,j,iEl
    real(prec) :: rho,u,v,p,c

    do concurrent(i=1:this%diagnostics%N+1,j=1:this%diagnostics%N+1, &
      iel=1:this%mesh%nElem)
      rho = this%primitive%interior(i,j,iEl,1)
      u = this%primitive%interior(i,j,iEl,2)
      v = this%primitive%interior(i,j,iEl,3)
      p = this%primitive%interior(i,j,iEl,4)
      c =  sqrt(this%gamma*p/rho) 
      this%diagnostics%interior(i,j,iEl,1) = c ! Speed of sound
      this%diagnostics%interior(i,j,iEl,2) = p/(rho*c*c) ! Temperature
    enddo

    call this%diagnostics%BoundaryInterp()

  endsubroutine CalculateDiagnostics_Euler2D_t

  subroutine ConservativeToPrimitive_Euler2D_t(this)
    implicit none
    class(Euler2D_t),intent(inout) :: this
    ! Local
    integer :: i,j,iEl

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
      iel=1:this%mesh%nElem)
      this%primitive%interior(i,j,iEl,1) = this%solution%interior(i,j,iEl,1) ! density
      this%primitive%interior(i,j,iEl,2) = this%solution%interior(i,j,iEl,2)/this%solution%interior(i,j,iEl,1) ! x-velocity
      this%primitive%interior(i,j,iEl,3) = this%solution%interior(i,j,iEl,3)/this%solution%interior(i,j,iEl,1) ! y-velocity
      this%primitive%interior(i,j,iEl,4) = (this%solution%interior(i,j,iEl,4)- &
       0.5_prec*(this%solution%interior(i,j,iEl,2)**2+this%solution%interior(i,j,iEl,3)**2)/this%solution%interior(i,j,iEl,1)) &
       *(this%gamma - 1.0_prec) ! pressure
    enddo
  
    call this%primitive%BoundaryInterp()
    call this%primitive%SideExchange(this%mesh)

  endsubroutine ConservativeToPrimitive_Euler2D_t

  subroutine PrimitiveToConservative_Euler2D_t(this)
    implicit none
    class(Euler2D_t),intent(inout) :: this
    ! Local
    integer :: i,j,iEl

    do concurrent(i=1:this%primitive%N+1,j=1:this%primitive%N+1, &
      iel=1:this%mesh%nElem)
      this%solution%interior(i,j,iEl,1) = this%primitive%interior(i,j,iEl,1) ! density
      this%solution%interior(i,j,iEl,2) = this%primitive%interior(i,j,iEl,2)*this%primitive%interior(i,j,iEl,1) ! x-momentum
      this%solution%interior(i,j,iEl,3) = this%primitive%interior(i,j,iEl,3)*this%primitive%interior(i,j,iEl,1) ! y-momentum
      this%solution%interior(i,j,iEl,4) = this%primitive%interior(i,j,iEl,4)/(this%gamma - 1.0_prec) + &
        0.5_prec*this%primitive%interior(i,j,iEl,1)*(this%primitive%interior(i,j,iEl,2)**2+ &
        this%primitive%interior(i,j,iEl,3)**2) ! total energy
    enddo

    call this%solution%BoundaryInterp()
    call this%solution%SideExchange(this%mesh)

  endsubroutine PrimitiveToConservative_Euler2D_t

  pure function entropy_func_Euler2D_t(this,s) result(e)
    !! The entropy function is proportional to the thermodynamic
    !! entropy. For the Euler equations, the entropy function is
    !! given by 
    !!
    !!    e = -(ln(p) - γ*ln(ρ))*ρ/(γ-1)
    !! 
    !! where p is the pressure, ρ is the density, and γ is the ratio
    !! of specific heats.
    !!
    class(Euler2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: e
    ! Local
    real(prec) :: p

    p = (s(4)-0.5_prec*(s(2)**2+s(3)**2)/s(1))*(this%gamma - 1.0_prec) ! pressure
    e = -(log(p) - this%gamma*log(s(1)))*s(1)/(this%gamma - 1.0_prec) ! mathematical entropy

  endfunction entropy_func_Euler2D_t

  pure function hbc2d_NoNormalFlow_Euler2D_t(this,s,nhat) result(exts)
    class(Euler2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: exts(1:this%nvar)
    ! Local
    integer :: ivar

    exts(1) = s(1) ! density
    exts(2) = (nhat(2)**2-nhat(1)**2)*s(2)-2.0_prec*nhat(1)*nhat(2)*s(3) ! \rho*u
    exts(3) = (nhat(1)**2-nhat(2)**2)*s(3)-2.0_prec*nhat(1)*nhat(2)*s(2) ! \rho*v
    exts(4) = s(4) ! energy

  endfunction hbc2d_NoNormalFlow_Euler2D_t

  pure function pressure(this,s) result(p)
    class(Euler2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: p

    p = (s(4)-0.5_prec*(s(2)**2+s(3)**2)/s(1))*(this%gamma - 1.0_prec)

  endfunction pressure

  pure function speedofsound(this,s) result(c)
    class(Euler2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: c

    c = sqrt(this%gamma*this%pressure(s)/s(1))
  
  endfunction speedofsound

  pure function flux2d_Euler2D_t(this,s,dsdx) result(flux)
    class(Euler2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec) :: flux(1:this%nvar,1:2)
    ! Local
    real(prec) :: p

    ! Computes the pressure for an ideal gas
    p = this%pressure(s)

    flux(1,1) = s(2) ! density, x flux ; rho*u
    flux(1,2) = s(3) ! density, y flux ; rho*v
    flux(2,1) = s(2)*s(2)/s(1) + p ! x-momentum, x flux; \rho*u*u + p
    flux(2,2) = s(2)*s(3)/s(1) ! x-momentum, y flux; \rho*u*v
    flux(3,1) = s(2)*s(3)/s(1) ! y-momentum, x flux; \rho*v*u
    flux(3,2) = s(3)*s(3)/s(1) + p ! y-momentum, y flux; \rho*v*v + p
    flux(4,1) = (s(4) + p)*s(2)/s(1) ! total energy, x flux : (\rho*E + p)*u
    flux(4,2) = (s(4) + p)*s(3)/s(1) ! total energy, y flux : (\rho*E + p)*v

  endfunction flux2d_Euler2D_t

  pure function riemannflux2d_Euler2D_t(this,sL,sR,dsdx,nhat) result(flux)
    !! Uses a local lax-friedrich's upwind flux
    !! The max eigenvalue is taken as 
    !!
    !!  max( |uL + cL|, |uL - cL|, |uR + cL|, |uR - cR| )
    !!
    class(Euler2D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: flux(1:this%nvar)
    ! Local
    real(prec) :: fL(1:this%nvar)
    real(prec) :: fR(1:this%nvar)
    real(prec) :: rhoun,v,p,c,rho0
    real(prec) :: unl, unr, cL,cR
    real(prec) :: lambda

    rhoun = sL(2)*nhat(1) + sL(3)*nhat(2)    
    p = this%pressure(sL)
    fL(1) = rhoun ! density
    unl = rhoun/sL(1)
    fL(2) = sL(2)*unl + p*nhat(1) ! x-momentum
    fL(3) = sL(2)*unl + p*nhat(2) ! y-momentum
    fL(4) = (sL(4) + p)*unl! total energy
    

    rhoun = sR(2)*nhat(1) + sR(3)*nhat(2)    
    p = this%pressure(sR)
    unr = rhoun/sR(1)
    fR(1) = rhoun ! density
    fR(2) = sR(2)*unr + p*nhat(1) ! x-momentum
    fR(3) = sR(2)*unr + p*nhat(2) ! y-momentum
    fR(4) = (sR(4) + p)*unr ! total energy

    cL = this%speedofsound(sL)
    cR = this%speedofsound(sR)
    lambda = max( abs(unr + cR), abs(unl - cL), abs(unl + cL), abs(unr - cR) )

    flux(1:4) = 0.5_prec*(fL(1:4)+fR(1:4))+lambda*(sL(1:4)-sR(1:4))

  endfunction riemannflux2d_Euler2D_t

  subroutine WeakBlastWave_Euler2D_t(this,rho0,rhoprime,Lr,P0,Eprime,Le,x0,y0)
    !! This subroutine sets the initial condition for a weak blast wave
    !! problem. The initial condition is given by
    !!
    !! \begin{equation}
    !! \begin{aligned}
    !! \rho &= \begin{cases}
    !!   \rho_0 & \text{if } r > 0.5 \\
    !!  \rho_0 + \frac{\rho'}{2\pi L_r^2} \exp\left( -\ln(2) \frac{(x-x_0)^2 + (y-y_0)^2}{L_r^2} \right) & \text{otherwise}
    !!  \end{cases} \\
    !! u &= 0 \\
    !! v &= 0 \\
    !! E &= \begin{cases}
    !!  \frac{P_0}{\gamma - 1} & \text{if } r > 0.5 \\
    !!  \frac{P_0}{\gamma - 1} + \frac{E'}{2\pi L_e^2} \exp\left( -\ln(2) \frac{(x-x_0)^2 + (y-y_0)^2}{L_e^2} \right) & \text{otherwise}
    !!  \end{cases}
    !! \end{aligned}
    !! \end{equation}
    !!
    implicit none
    class(Euler2D_t),intent(inout) :: this
    real(prec),intent(in) :: rho0,rhoprime,Lr,P0,Eprime,Le,x0,y0
    ! Local
    integer :: i,j,iEl
    real(prec) :: x,y,rho,r,E

    print*, __FILE__, " : Configuring weak blast wave initial condition. "
    print*, __FILE__, " : rho0 = ", rho0
    print*, __FILE__, " : rhoprime = ", rhoprime
    print*, __FILE__, " : Lr = ", Lr
    print*, __FILE__, " : P0 = ", P0
    print*, __FILE__, " : Eprime = ", Eprime
    print*, __FILE__, " : Le = ", Le
    print*, __FILE__, " : x0 = ", x0
    print*, __FILE__, " : y0 = ", y0

    do concurrent(i=1:this%primitive%N+1,j=1:this%primitive%N+1, &
      iel=1:this%mesh%nElem)
      x = this%geometry%x%interior(i,j,iEl,1,1) - x0
      y = this%geometry%x%interior(i,j,iEl,1,2) - y0
      r  = sqrt(x**2 + y**2)
      if (r > 0.5_prec) then
          rho = rho0
          E = P0/(this%gamma - 1.0_prec)
      else
          rho = rho0 + (rhoprime/(2.0_prec*pi*Lr**2))*exp(-log(2.0_prec)*(x**2 + y**2)/Lr**2)
          E = P0/(this%gamma - 1.0_prec) + (Eprime/(2.0_prec*pi*Le**2))*exp(-log(2.0_prec)*(x**2 + y**2)/Le**2)
      end if

      this%solution%interior(i,j,iEl,1) = rho
      this%solution%interior(i,j,iEl,2) = 0.0_prec
      this%solution%interior(i,j,iEl,3) = 0.0_prec
      this%solution%interior(i,j,iEl,4) = E

    enddo

    call this%ConservativeToPrimitive()
    call this%CalculateDiagnostics()

    print*, __FILE__, " : Max sound speed = ", maxval(this%diagnostics%interior(:,:,:,1))

  endsubroutine WeakBlastWave_Euler2D_t

endmodule self_Euler2D_t
