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

module self_GFDLES3D_t
!! This module defines a class that can be used to solve the filtered
!! compressible navier-stokes equations in 3-D
!!
!! The conserved variables are
!!
!! \begin{equation}
!! \vec{s} = \begin{pmatrix}
!!     \rho \\
!!     \rho u \\
!!     \rho v \\
!!     \rho w \\
!!     \rho \theta
!!  \end{pmatrix}
!! \end{equation}
!!
!! The conservative flux is
!!
!! \begin{equation}
!! \overleftrightarrow{f} = \begin{pmatrix}
!!     \rho u \hat{x} + \rho v \hat{y} + \rho w \hat{z} \\
!!      \vec{u} \rho u + \frac{p}{\rho_0} \hat{x} \\
!!      \vec{u} \rho v + \frac{p}{\rho_0} \hat{y} \\
!!      \vec{u} \rho w + \frac{p}{\rho_0} \hat{z} \\
!!      \vec{u} \rho \theta
!!  \end{pmatrix}
!! \end{equation}
!!
!! and the source terms include the graviational acceleration, where
!! gravity is assumed a constant acting in the z-direction
!!
!! \begin{equation}
!! \overleftrightarrow{f} = \begin{pmatrix}
!!      0 \\
!!      0 \\
!!      0 \\
!!      -\rho g \\
!!      0
!!  \end{pmatrix}
!! \end{equation}
!!
!!  ...Subgrid-scale closure...
!!
  use self_model
  use self_dgmodel3D
  use self_mesh

  implicit none

  type,extends(dgmodel3D) :: GFDLES3D_t
    type(MappedScalar3D)   :: primitive
    type(MappedScalar3D)   :: diagnostics
    integer :: ndiagnostics

    ! Model parameters
    real(prec) :: p0 = 10.0_prec**(5) ! Reference pressure for potential temperature ()
    real(prec) :: Cp = 1.005_prec*10.0_prec**(3) ! Specific heat at constant pressure (J/kg-K)
    real(prec) :: Cv = 0.718_prec*10.0_prec**(3) ! Specific heat at constant volume (J/kg-K)
    real(prec) :: R = 287.04_prec ! Gas constant (Cp - Cv)
    real(prec) :: gamma = 1.399721448_prec ! Ratio of specific heats (Cp/Cv)
    real(prec) :: nu = 0.0_prec ! Dynamic viscosity
    real(prec) :: kappa = 0.0_prec ! Thermal diffusivity
    real(prec) :: g = 9.81_prec ! gravitational acceleration (z-direction only)

    logical :: sgs_enabled = .false.
  contains

    ! Setup / Book-keeping methods
    procedure :: AdditionalInit => AdditionalInit_GFDLES3D_t
    procedure :: AdditionalFree => AdditionalFree_GFDLES3D_t
    procedure :: SetMetadata => SetMetadata_GFDLES3D_t
    procedure :: SetNumberOfVariables => SetNumberOfVariables_GFDLES3D_t

    ! File IO methods
    !procedure :: AdditionalOutput => AdditionalOutput_GFDLES3D_t

    ! Pre-tendency methods
    procedure :: CalculateDiagnostics => CalculateDiagnostics_GFDLES3D_t
    procedure :: ConservativeToPrimitive => ConservativeToPrimitive_GFDLES3D_t
    procedure :: SetPrimitiveBoundaryCondition => setprimitiveboundarycondition_GFDLES3D_t
    procedure :: PreTendency => PreTendency_GFDLES3D_t

    ! Model method overrides
    !procedure :: hbc2d_NoNormalFlow => hbc2d_NoNormalFlow_GFDLES3D_t
    !procedure :: pbc2d_NoNormalFlow => pbc2d_NoNormalFlow_GFDLES3D_t

    !procedure :: SourceMethod => sourcemethod_GFDLES3D_t
    procedure :: entropy_func => entropy_func_GFDLES3D_t
    procedure :: flux3D => flux3D_GFDLES3D_t
    procedure :: riemannflux3D => riemannflux3D_GFDLES3D_t

    ! Additional support methods
    procedure :: ReportUserMetrics => ReportUserMetrics_GFDLES3D_t
    procedure :: PrimitiveToConservative => PrimitiveToConservative_GFDLES3D_t
    procedure,private :: pressure
    !procedure,private :: temperature
    procedure,private :: speedofsound

    ! Example Initial Conditions
    !procedure :: SphericalSoundWave => SphericalSoundWave_GFDLES3D_t

  endtype GFDLES3D_t

contains

  subroutine AdditionalInit_GFDLES3D_t(this)
    implicit none
    class(GFDLES3D_t),intent(inout) :: this

    call this%primitive%Init(this%geometry%x%interp,this%nvar,this%mesh%nElem)
    call this%primitive%AssociateGeometry(this%geometry)
    call this%diagnostics%Init(this%geometry%x%interp,this%ndiagnostics,this%mesh%nElem)

  endsubroutine AdditionalInit_GFDLES3D_t

  subroutine AdditionalFree_GFDLES3D_t(this)
    implicit none
    class(GFDLES3D_t),intent(inout) :: this

    call this%primitive%Free()
    !call this%primitiveGradient%Free()
    call this%diagnostics%Free()

  endsubroutine AdditionalFree_GFDLES3D_t

  subroutine SetNumberOfVariables_GFDLES3D_t(this)
    implicit none
    class(GFDLES3D_t),intent(inout) :: this

    this%nvar = 5
    this%ndiagnostics = 4

  endsubroutine SetNumberOfVariables_GFDLES3D_t

  subroutine SetMetadata_GFDLES3D_t(this)
    implicit none
    class(GFDLES3D_t),intent(inout) :: this

    call this%solution%SetName(1,"ρ") ! Density
    call this%solution%SetDescription(1,"Density") ! Density
    call this%solution%SetUnits(1,"kg⋅m⁻³")

    call this%solution%SetName(2,"ρu") ! x-momentum
    call this%solution%SetDescription(2,"x-momentum")
    call this%solution%SetUnits(2,"(kg⋅m⁻³)(m⋅s⁻¹)")

    call this%solution%SetName(3,"ρv") ! y-momentum
    call this%solution%SetDescription(3,"y-momentum")
    call this%solution%SetUnits(3,"(kg⋅m⁻³)(m⋅s⁻¹)")

    call this%solution%SetName(4,"ρw") ! z-momentum
    call this%solution%SetDescription(4,"z-momentum")
    call this%solution%SetUnits(4,"(kg⋅m⁻³)(m⋅s⁻¹)")

    call this%solution%SetName(5,"ρθ") ! Density weighted potential temperature
    call this%solution%SetDescription(5,"Density weighted potential temperature")
    call this%solution%SetUnits(5,"(kg⋅m⁻³)(m²⋅s⁻²)")

    call this%primitive%SetName(1,"ρ") ! Density
    call this%primitive%SetDescription(1,"Density") ! Density
    call this%primitive%SetUnits(1,"kg⋅m⁻³")

    call this%primitive%SetName(2,"u") ! x-velocity
    call this%primitive%SetDescription(2,"x-velocity")
    call this%primitive%SetUnits(2,"(m⋅s⁻¹)")

    call this%primitive%SetName(3,"v") ! y-momentum
    call this%primitive%SetDescription(3,"y-velocity")
    call this%primitive%SetUnits(3,"(m⋅s⁻¹)")

    call this%primitive%SetName(4,"w") ! z-momentum
    call this%primitive%SetDescription(4,"z-velocity")
    call this%primitive%SetUnits(4,"(m⋅s⁻¹)")

    call this%primitive%SetName(5,"θ") ! in-situ temperature
    call this%primitive%SetDescription(5,"Potential temperature")
    call this%primitive%SetUnits(5,"K")

    call this%diagnostics%SetName(1,"c") ! Speed of sound
    call this%diagnostics%SetDescription(1,"Speed of sound")
    call this%diagnostics%SetUnits(1,"m⋅s⁻¹")

    call this%diagnostics%SetName(2,"P") ! Pressure
    call this%diagnostics%SetDescription(2,"Pressure")
    call this%diagnostics%SetUnits(2,"kg⋅m⁻¹⋅s⁻²")

    call this%diagnostics%SetName(3,"ρK") ! kinetic energy
    call this%diagnostics%SetDescription(3,"Kinetic energy")
    call this%diagnostics%SetUnits(3,"(kg⋅m⁻³)(m²⋅s⁻²)")

    call this%diagnostics%SetName(4,"CFL-J") ! kinetic energy
    call this%diagnostics%SetDescription(4,"CFL number using the |u|*dt/\sqrt{J}")
    call this%diagnostics%SetUnits(4,"-")

  endsubroutine SetMetadata_GFDLES3D_t

  subroutine ReportUserMetrics_GFDLES3D_t(this)
  !! Base method for reporting the entropy of a model
  !! to stdout. Only override this procedure if additional
  !! reporting is needed. Alternatively, if you think
  !! additional reporting would be valuable for all models,
  !! open a pull request with modifications to this base
  !! method.
    implicit none
    class(GFDLES3D_t),intent(inout) :: this
    ! Local
    character(len=20) :: modelTime
    character(len=20) :: minv,maxv
    character(len=:),allocatable :: str
    integer :: ivar

    call this%ConservativeToPrimitive()
    call this%CalculateDiagnostics()

    ! Copy the time and entropy to a string
    write(modelTime,"(ES16.7E3)") this%t

    do ivar = 1,this%nvar
      write(maxv,"(ES16.7E3)") maxval(this%primitive%interior(:,:,:,:,ivar))
      write(minv,"(ES16.7E3)") minval(this%primitive%interior(:,:,:,:,ivar))

      ! Write the output to STDOUT
      open(output_unit,ENCODING='utf-8')
      write(output_unit,'(1x, A," : ")',ADVANCE='no') __FILE__
      str = 'tᵢ ='//trim(modelTime)
      write(output_unit,'(A)',ADVANCE='no') str
 str = '  |  min('//trim(this%primitive%meta(ivar)%name)//'), max('//trim(this%primitive%meta(ivar)%name)//') = '//minv//" , "//maxv
      write(output_unit,'(A)',ADVANCE='yes') str
    enddo

    do ivar = 1,this%ndiagnostics
      write(maxv,"(ES16.7E3)") maxval(this%diagnostics%interior(:,:,:,:,ivar))
      write(minv,"(ES16.7E3)") minval(this%diagnostics%interior(:,:,:,:,ivar))

      ! Write the output to STDOUT
      open(output_unit,ENCODING='utf-8')
      write(output_unit,'(1x,A," : ")',ADVANCE='no') __FILE__
      str = 'tᵢ ='//trim(modelTime)
      write(output_unit,'(A)',ADVANCE='no') str
      str = '  |  min('//trim(this%diagnostics%meta(ivar)%name)//'), max('//trim(this%diagnostics%meta(ivar)%name)//') = '//minv//" , "//maxv
      write(output_unit,'(A)',ADVANCE='yes') str
    enddo

  endsubroutine ReportUserMetrics_GFDLES3D_t

  subroutine setprimitiveboundarycondition_GFDLES3D_t(this)
    !! Boundary conditions for the solution are set to
    !! 0 for the external state to provide radiation type
    !! boundary conditions.
    implicit none
    class(GFDLES3D_t),intent(inout) :: this
    ! local
    integer :: i,j,k,iEl,e2,bcid
    real(prec) :: nhat(1:3)

    do concurrent(k=1:6,iel=1:this%mesh%nElem)

      bcid = this%mesh%sideInfo(5,j,iEl) ! Boundary Condition ID
      e2 = this%mesh%sideInfo(3,j,iEl) ! Neighboring Element ID

      if(e2 == 0) then
        if(bcid == SELF_BC_PRESCRIBED) then
          ! To do : need to set different prescribed function for the primitive variables
          do j = 1,this%solution%interp%N+1 ! Loop over quadrature points
            do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
              nhat = this%geometry%nhat%boundary(i,j,k,iEl,1,1:3)

              this%primitive%extBoundary(i,j,k,iEl,1:this%nvar) = &
                this%hbc3d_Prescribed(this%primitive%boundary(i,j,k,iEl,1:this%nvar),this%t)
            enddo
          enddo

        elseif(bcid == SELF_BC_RADIATION) then
          ! To do : need to set different prescribed function for the primitive variables
          do j = 1,this%solution%interp%N+1 ! Loop over quadrature points
            do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
              nhat = this%geometry%nhat%boundary(i,j,k,iEl,1,1:3)

              this%primitive%extBoundary(i,j,k,iEl,1:this%nvar) = &
                this%hbc3d_Radiation(this%primitive%boundary(i,j,k,iEl,1:this%nvar),nhat)
            enddo
          enddo

        elseif(bcid == SELF_BC_NONORMALFLOW) then
          do j = 1,this%solution%interp%N+1 ! Loop over quadrature points
            do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
              nhat = this%geometry%nhat%boundary(i,j,k,iEl,1,1:3)

              this%primitive%extBoundary(i,j,k,iEl,1:this%nvar) = &
                this%hbc3d_NoNormalFlow(this%primitive%boundary(i,j,k,iEl,1:this%nvar),nhat)
            enddo
          enddo

        endif
      endif

    enddo

  endsubroutine setprimitiveboundarycondition_GFDLES3D_t
  subroutine PreTendency_GFDLES3D_t(this)
    implicit none
    class(GFDLES3D_t),intent(inout) :: this

    if(this%sgs_enabled) then
      call this%ConservativeToPrimitive()
      call this%primitive%BoundaryInterp()
      call this%primitive%SideExchange(this%mesh)

      call this%SetPrimitiveBoundaryCondition()

      call this%primitive%AverageSides()

      ! Compute the gradient of the primitive variables
      ! and store the result in the solutionGradient property.
      call this%primitive%MappedDGGradient(this%solutionGradient%interior)
      call this%solutionGradient%BoundaryInterp()
      call this%solutionGradient%SideExchange(this%mesh)
      call this%SetGradientBoundaryCondition()
      call this%solutionGradient%AverageSides()
    endif

  endsubroutine PreTendency_GFDLES3D_t

  subroutine CalculateDiagnostics_GFDLES3D_t(this)
    implicit none
    class(GFDLES3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,k,iEl
    real(prec) :: c,e,ke
    real(prec) :: s(1:this%nvar)

    do concurrent(i=1:this%diagnostics%N+1,j=1:this%diagnostics%N+1, &
                  k=1:this%diagnostics%N+1,iel=1:this%mesh%nElem)
      s(1:this%nvar) = this%solution%interior(i,j,k,iEl,1:this%nvar)
      c = this%speedofsound(s)
      ke = 0.5_prec*(s(2)**2+s(3)**2+s(4)**2)/s(1) ! kinetic energy (kg⋅m²⋅s⁻²)
      this%diagnostics%interior(i,j,k,iEl,1) = c ! Speed of sound
      this%diagnostics%interior(i,j,k,iEl,2) = this%pressure(s) ! Pressure (total)
      this%diagnostics%interior(i,j,k,iEl,3) = ke ! kinetic energy
      this%diagnostics%interior(i,j,k,iEl,4) = (sqrt(ke)+c)*this%dt/sqrt(this%geometry%J%interior(i,j,k,iEl,1)) ! CFL number
    enddo

  endsubroutine CalculateDiagnostics_GFDLES3D_t

  subroutine ConservativeToPrimitive_GFDLES3D_t(this)
    implicit none
    class(GFDLES3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,k,iEl
    real(prec) :: s(1:this%nvar)

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%diagnostics%N+1,iel=1:this%mesh%nElem)
      s(1:this%nvar) = this%solution%interior(i,j,k,iEl,1:this%nvar)
      this%primitive%interior(i,j,k,iEl,1) = s(1) ! density
      this%primitive%interior(i,j,k,iEl,2) = s(2)/s(1) ! x-velocity
      this%primitive%interior(i,j,k,iEl,3) = s(3)/s(1) ! y-velocity
      this%primitive%interior(i,j,k,iEl,4) = s(4)/s(1) ! z-velocity
      this%primitive%interior(i,j,k,iEl,5) = s(5)/s(1) ! Potential temperature
    enddo

  endsubroutine ConservativeToPrimitive_GFDLES3D_t

  subroutine PrimitiveToConservative_GFDLES3D_t(this)
    implicit none
    class(GFDLES3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,k,iEl
    real(prec) :: s(1:this%nvar)

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%diagnostics%N+1,iel=1:this%mesh%nElem)
      s(1:this%nvar) = this%primitive%interior(i,j,k,iEl,1:this%nvar)
      this%solution%interior(i,j,k,iEl,1) = s(1) ! density
      this%solution%interior(i,j,k,iEl,2) = s(2)*s(1) ! x-momentum
      this%solution%interior(i,j,k,iEl,3) = s(3)*s(1) ! y-momentum
      this%solution%interior(i,j,k,iEl,4) = s(4)*s(1) ! z-momentum
      this%solution%interior(i,j,k,iEl,5) = s(5)*s(1) ! Density weighted Potential temperature
    enddo

  endsubroutine PrimitiveToConservative_GFDLES3D_t

  pure function entropy_func_GFDLES3D_t(this,s) result(e)
    !! The entropy function is the sum of kinetic and internal energy
    !! For the linear model, this is
    !!
    !! \begin{equation}
    !!   e = \frac{1}{2} \left( \rho_0*( u^2 + v^2 ) + \frac{P^2}{\rho_0 c^2} \right)
    class(GFDLES3D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: e
    ! Local
    real(prec) :: ke,ie,pe

    ke = 0.5_prec*(s(2)*s(2)+s(3)*(3)+s(4)*s(4))/s(1) ! kinetic energy
    !pe = s(1)*this%g*z Potential energy
    ie = this%Cv*s(5)*(this%pressure(s)/this%p0)**(this%R/this%Cp) ! internal energy = rho*Cv*T

    e = ke+ie
  endfunction entropy_func_GFDLES3D_t

  ! pure function hbc3D_NoNormalFlow_GFDLES3D_t(this,s,nhat) result(exts)
  !   class(GFDLES3D_t),intent(in) :: this
  !   real(prec),intent(in) :: s(1:this%nvar)
  !   real(prec),intent(in) :: nhat(1:2)
  !   real(prec) :: exts(1:this%nvar)
  !   ! Local
  !   integer :: ivar

  !   exts(1) = s(1) ! density
  !   exts(2) = (nhat(2)**2-nhat(1)**2)*s(2)-2.0_prec*nhat(1)*nhat(2)*s(3) ! u
  !   exts(3) = (nhat(1)**2-nhat(2)**2)*s(3)-2.0_prec*nhat(1)*nhat(2)*s(2) ! v
  !   exts(4) = (nhat(1)**2-nhat(2)**2)*s(3)-2.0_prec*nhat(1)*nhat(2)*s(2) ! w
  !   exts(5) = s(4) ! p

  ! endfunction hbc3D_NoNormalFlow_GFDLES3D_t

  pure function pressure(this,s) result(p)
    class(GFDLES3D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: p

    p = (this%R*s(5)*(this%p0)**(-this%R/this%Cp))**(this%gamma)

  endfunction pressure

  !  pure function temperature(this, s) result(t)
  !     class(GFDLES3D_t), intent(in) :: this
  !     real(prec), intent(in) :: s(1:this%nvar)
  !     real(prec) :: t

  !     t = (s(4) - 0.5_prec*(s(2)**2 + s(3)**2)/s(1))/(s(1)*this%Cv) ! temperature = e/Cv

  !  end function temperature

  pure function speedofsound(this,s) result(c)
    class(GFDLES3D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: c

    c = sqrt(this%gamma*this%pressure(s)/s(1))

  endfunction speedofsound

  pure function flux3d_GFDLES3D_t(this,s,dsdx) result(flux)
    class(GFDLES3D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec) :: flux(1:this%nvar,1:3)
    ! Local
    real(prec) :: p,nu,kappa,u,v,w
    real(prec) :: tau_11,tau_12,tau_13
    real(prec) :: tau_22,tau_23
    real(prec) :: tau_33

    ! Computes the pressure for an ideal gas
    p = this%pressure(s)
    u = s(2)/s(1)
    v = s(3)/s(1)
    w = s(4)/s(1)
    ! LEFT OFF HERE !!
    flux(1,1) = s(2) ! density, x flux ; rho*u
    flux(1,2) = s(3) ! density, y flux ; rho*v
    flux(1,3) = s(3) ! density, z flux ; rho*w

    flux(2,1) = s(2)*u+p ! x-momentum, x flux; \rho*u*u + p
    flux(2,2) = s(2)*v ! x-momentum, y flux; \rho*u*v
    flux(3,1) = s(2)*u ! y-momentum, x flux; \rho*v*u
    flux(3,2) = s(3)*v+p ! y-momentum, y flux; \rho*v*v + p
    flux(4,1) = (s(4)+p)*s(2)/s(1) ! total energy, x flux : (\rho*E + p)*u
    flux(4,2) = (s(4)+p)*s(3)/s(1) ! total energy, y flux : (\rho*E + p)*v

    if(this%sgs_enabled) then
      ! Viscous and difussive terms
      ! Recall that the solutionGradient now contains
      ! the primitive variable gradients
      ! Calculate the stress tensor
      nu = this%nu
      kappa = this%kappa
      tau_11 = 4.0_prec*dsdx(2,1)/3.0_prec-2.0_prec*dsdx(3,2)/3.0_prec
      tau_12 = dsdx(2,2)+dsdx(3,1)
      !tau_21 = tau_12
      tau_22 = 4.0_prec*dsdx(3,2)/3.0_prec-2.0_prec*dsdx(2,1)/3.0_prec

      flux(2,1) = flux(2,1)-nu*tau_11 ! x-momentum, x flux
      flux(2,2) = flux(2,2)-nu*tau_12 ! x-momentum, y flux (-tau_21*nu = -tau_12*nu)
      flux(3,1) = flux(3,1)-nu*tau_12 ! y-momentum, x flux
      flux(3,2) = flux(3,2)-nu*tau_22 ! y-momentum, y flux
      flux(4,1) = flux(4,1)-(kappa*dsdx(4,1)+u*tau_11+v*tau_12) ! total energy, x flux = -(kappa*dTdx + u*tau_11 + v*tau_12)
      flux(4,2) = flux(4,2)-(kappa*dsdx(4,2)+u*tau_11+v*tau_12) ! total energy, y flux = -(kappa*dTdy + u*tau_12 + v*tau_22)
    endif

  endfunction flux3d_GFDLES3D_t

  pure function riemannflux3D_GFDLES3D_t(this,sL,sR,dsdx,nhat) result(flux)
    !! Uses a local lax-friedrich's upwind flux
    !! The max eigenvalue is taken as the sound speed
    class(GFDLES3D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec),intent(in) :: nhat(1:3)
    real(prec) :: flux(1:this%nvar)
    ! Local
    real(prec) :: fL(1:this%nvar)
    real(prec) :: fR(1:this%nvar)
    real(prec) :: u,v,w,p,c,rho0

    u = sL(2)
    v = sL(3)
    w = sL(4)
    p = sL(5)
    rho0 = 1.0 !this%rho0
    c = 1.0 !this%c
    fL(1) = rho0*(u*nhat(1)+v*nhat(2)+w*nhat(3)) ! density
    fL(2) = p*nhat(1)/rho0 ! u
    fL(3) = p*nhat(2)/rho0 ! v
    fL(4) = p*nhat(3)/rho0 ! w
    fL(5) = rho0*c*c*(u*nhat(1)+v*nhat(2)+w*nhat(3)) ! pressure

    u = sR(2)
    v = sR(3)
    w = sR(4)
    p = sR(5)
    fR(1) = rho0*(u*nhat(1)+v*nhat(2)+w*nhat(3)) ! density
    fR(2) = p*nhat(1)/rho0 ! u
    fR(3) = p*nhat(2)/rho0 ! v'
    fR(4) = p*nhat(3)/rho0 ! w
    fR(5) = rho0*c*c*(u*nhat(1)+v*nhat(2)+w*nhat(3)) ! pressure

    flux(1:5) = 0.5_prec*(fL(1:5)+fR(1:5))+c*(sL(1:5)-sR(1:5))

  endfunction riemannflux3D_GFDLES3D_t

  ! subroutine SphericalSoundWave_GFDLES3D_t(this,rhoprime,Lr,x0,y0,z0)
  ! !! This subroutine sets the initial condition for a weak blast wave
  ! !! problem. The initial condition is given by
  ! !!
  ! !! \begin{equation}
  ! !! \begin{aligned}
  ! !! \rho &= \rho_0 + \rho' \exp\left( -\ln(2) \frac{(x-x_0)^2 + (y-y_0)^2}{L_r^2} \right)
  ! !! u &= 0 \\
  ! !! v &= 0 \\
  ! !! E &= \frac{P_0}{\gamma - 1} + E \exp\left( -\ln(2) \frac{(x-x_0)^2 + (y-y_0)^2}{L_e^2} \right)
  ! !! \end{aligned}
  ! !! \end{equation}
  ! !!
  !   implicit none
  !   class(GFDLES3D_t),intent(inout) :: this
  !   real(prec),intent(in) ::  rhoprime,Lr,x0,y0,z0
  !   ! Local
  !   integer :: i,j,k,iEl
  !   real(prec) :: x,y,z,rho,r,E

  !   print*,__FILE__," : Configuring weak blast wave initial condition. "
  !   print*,__FILE__," : rhoprime = ",rhoprime
  !   print*,__FILE__," : Lr = ",Lr
  !   print*,__FILE__," : x0 = ",x0
  !   print*,__FILE__," : y0 = ",y0
  !   print*,__FILE__," : z0 = ",z0

  !   do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
  !                 k=1:this%solution%N+1,iel=1:this%mesh%nElem)
  !     x = this%geometry%x%interior(i,j,k,iEl,1,1)-x0
  !     y = this%geometry%x%interior(i,j,k,iEl,1,2)-y0
  !     z = this%geometry%x%interior(i,j,k,iEl,1,3)-z0
  !     r = sqrt(x**2+y**2+z**2)

  !     rho = (rhoprime)*exp(-log(2.0_prec)*r**2/Lr**2)

  !     this%solution%interior(i,j,k,iEl,1) = rho
  !     this%solution%interior(i,j,k,iEl,2) = 0.0_prec
  !     this%solution%interior(i,j,k,iEl,3) = 0.0_prec
  !     this%solution%interior(i,j,k,iEl,4) = 0.0_prec
  !     this%solution%interior(i,j,k,iEl,5) = rho*this%c*this%c

  !   enddo

  !   call this%ReportMetrics()
  !   call this%solution%UpdateDevice()

  ! endsubroutine SphericalSoundWave_GFDLES3D_t

endmodule self_GFDLES3D_t
