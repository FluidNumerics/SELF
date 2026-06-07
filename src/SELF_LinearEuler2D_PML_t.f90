! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2026 Fluid Numerics LLC
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

module self_LinearEuler2D_PML_t
!! Linear Euler 2D with a Perfectly Matched Layer (PML) absorbing region.
!!
!! Formulation: Hu (2001), "A Stable, Perfectly Matched Layer for Linearized
!! Euler Equations in Unsplit Physical Variables" (JCP 173, 455-480), in the
!! ADE (Auxiliary Differential Equation) form for the no-mean-flow case.
!!
!! In the PML region the dynamics for the acoustic state
!! q = (rho', u, v, p) become
!!
!!   dq/dt + A dq/dx + B dq/dy + (sigma_x + sigma_y) q + sigma_x*sigma_y phi = 0
!!   dphi/dt = q
!!
!! where phi = (phi_rho, phi_u, phi_v, phi_P) is a 4-component auxiliary
!! that accumulates the time integral of q. The sigma_x(x), sigma_y(y)
!! damping coefficients vanish in the interior and rise toward the outer
!! boundary, leaving the original linear Euler 2D dynamics unaltered where
!! sigma = 0.
!!
!! Solution layout (nvar = 9):
!!   s(1) = rho     (density perturbation)
!!   s(2) = u       (x-velocity)
!!   s(3) = v       (y-velocity)
!!   s(4) = P       (pressure perturbation)
!!   s(5) = c       (sound speed, static; per-node)
!!   s(6) = phi_rho (auxiliary; time-integral of rho)
!!   s(7) = phi_u   (auxiliary; time-integral of u)
!!   s(8) = phi_v   (auxiliary; time-integral of v)
!!   s(9) = phi_P   (auxiliary; time-integral of P)
!!
!! The auxiliary variables phi_* have identically zero flux; they evolve
!! only via the source term. Sound speed c is also static (zero flux,
!! zero source) and is preserved from the parent LinearEuler2D model.
!!
!! PML elements are identified by a material name beginning with the
!! configurable prefix (default "pml") in the mesh's material table.
!! Non-PML elements always carry sigma_x = sigma_y = 0 and therefore
!! reduce to the parent linear Euler 2D dynamics (modulo the inert
!! auxiliary variables, which remain identically zero there).

  use self_model
  use self_dgmodel2d
  use self_mesh
  use self_lineareuler2d_t
  use SELF_BoundaryConditions
  use SELF_MappedScalar_2D

  implicit none

  type,extends(LinearEuler2D_t) :: LinearEuler2D_PML_t

    type(MappedScalar2D) :: sigma_x ! sigma_x(x,y) damping coefficient, one variable per node
    type(MappedScalar2D) :: sigma_y ! sigma_y(x,y) damping coefficient, one variable per node

    ! PML region tagging and ramp parameters (set by SetPMLProfile)
    real(prec) :: pml_x_min = 0.0_prec
    real(prec) :: pml_x_max = 0.0_prec
    real(prec) :: pml_y_min = 0.0_prec
    real(prec) :: pml_y_max = 0.0_prec
    real(prec) :: pml_width = 1.0_prec
    real(prec) :: pml_sigma_max = 0.0_prec
    integer :: pml_ramp_exponent = 3
    character(LEN=SELF_MESH_MATNAME_LENGTH) :: pml_material_prefix = "pml"

  contains
    procedure :: SetNumberOfVariables => SetNumberOfVariables_LinearEuler2D_PML_t
    procedure :: SetMetadata => SetMetadata_LinearEuler2D_PML_t
    procedure :: AdditionalInit => AdditionalInit_LinearEuler2D_PML_t
    procedure :: AdditionalFree => AdditionalFree_LinearEuler2D_PML_t
    procedure :: flux2d => flux2d_LinearEuler2D_PML_t
    procedure :: riemannflux2d => riemannflux2d_LinearEuler2D_PML_t
    procedure :: sourcemethod => sourcemethod_LinearEuler2D_PML_t
    procedure :: SetPMLProfile => SetPMLProfile_LinearEuler2D_PML_t

  endtype LinearEuler2D_PML_t

contains

  subroutine SetNumberOfVariables_LinearEuler2D_PML_t(this)
    implicit none
    class(LinearEuler2D_PML_t),intent(inout) :: this

    this%nvar = 9

  endsubroutine SetNumberOfVariables_LinearEuler2D_PML_t

  subroutine SetMetadata_LinearEuler2D_PML_t(this)
    implicit none
    class(LinearEuler2D_PML_t),intent(inout) :: this

    ! Reuse parent metadata for the first five variables.
    call SetMetadata_LinearEuler2D_t(this)

    call this%solution%SetName(6,"phi_rho")
    call this%solution%SetUnits(6,"kg⋅m⁻³⋅s")

    call this%solution%SetName(7,"phi_u")
    call this%solution%SetUnits(7,"m")

    call this%solution%SetName(8,"phi_v")
    call this%solution%SetUnits(8,"m")

    call this%solution%SetName(9,"phi_P")
    call this%solution%SetUnits(9,"kg⋅m⁻¹⋅s⁻¹")

    call this%sigma_x%SetName(1,"sigma_x")
    call this%sigma_x%SetUnits(1,"s⁻¹")
    call this%sigma_y%SetName(1,"sigma_y")
    call this%sigma_y%SetUnits(1,"s⁻¹")

  endsubroutine SetMetadata_LinearEuler2D_PML_t

  subroutine AdditionalInit_LinearEuler2D_PML_t(this)
    implicit none
    class(LinearEuler2D_PML_t),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    ! Allocate the per-node PML damping fields (single variable each).
    call this%sigma_x%Init(this%geometry%x%interp,1,this%mesh%nElem)
    call this%sigma_y%Init(this%geometry%x%interp,1,this%mesh%nElem)

    ! Register PML-aware BC methods. The PML-aware versions reuse the
    ! parent acoustic-side behaviour for variables 1-5 and zero the
    ! auxiliary variables 6-9 in extBoundary so that interior-side
    ! and exterior-side states are consistent with the zero-flux
    ! treatment of phi_*.
    bcfunc => hbc2d_NoNormalFlow_LinearEuler2D_PML
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

    bcfunc => hbc2d_Radiation_LinearEuler2D_PML
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_RADIATION,"radiation",bcfunc)

  endsubroutine AdditionalInit_LinearEuler2D_PML_t

  subroutine AdditionalFree_LinearEuler2D_PML_t(this)
    implicit none
    class(LinearEuler2D_PML_t),intent(inout) :: this

    call this%sigma_x%Free()
    call this%sigma_y%Free()

  endsubroutine AdditionalFree_LinearEuler2D_PML_t

  subroutine SetPMLProfile_LinearEuler2D_PML_t(this,x_interior_min,x_interior_max, &
                                               y_interior_min,y_interior_max, &
                                               pml_width,sigma_max,ramp_exponent)
    !! Populate the per-node sigma_x, sigma_y fields. Only nodes inside
    !! elements whose material name starts with this%pml_material_prefix
    !! receive non-zero damping; all other nodes are forced to zero so
    !! the interior solution remains unmodified.
    implicit none
    class(LinearEuler2D_PML_t),intent(inout) :: this
    real(prec),intent(in) :: x_interior_min,x_interior_max
    real(prec),intent(in) :: y_interior_min,y_interior_max
    real(prec),intent(in) :: pml_width
    real(prec),intent(in) :: sigma_max
    integer,intent(in),optional :: ramp_exponent
    ! Local
    integer :: i,j,iel,matid,p,prefix_len
    real(prec) :: x,y,dx,dy,sx,sy
    character(LEN=SELF_MESH_MATNAME_LENGTH) :: matname
    logical :: is_pml

    this%pml_x_min = x_interior_min
    this%pml_x_max = x_interior_max
    this%pml_y_min = y_interior_min
    this%pml_y_max = y_interior_max
    this%pml_width = pml_width
    this%pml_sigma_max = sigma_max
    if(present(ramp_exponent)) then
      this%pml_ramp_exponent = ramp_exponent
    endif
    p = this%pml_ramp_exponent
    prefix_len = len_trim(this%pml_material_prefix)

    ! Default everything to zero, then fill the PML elements.
    this%sigma_x%interior(:,:,:,1) = 0.0_prec
    this%sigma_y%interior(:,:,:,1) = 0.0_prec

    do iel = 1,this%mesh%nElem
      matid = this%mesh%elemMaterial(iel)
      matname = this%mesh%materialNames(matid)
      is_pml = (len_trim(matname) >= prefix_len) .and. &
               (matname(1:prefix_len) == this%pml_material_prefix(1:prefix_len))
      if(.not. is_pml) cycle

      do j = 1,this%solution%N+1
        do i = 1,this%solution%N+1
          x = this%geometry%x%interior(i,j,iel,1,1)
          y = this%geometry%x%interior(i,j,iel,1,2)

          ! Distance into PML measured from the interior box edge.
          dx = max(0.0_prec,x_interior_min-x,x-x_interior_max)
          dy = max(0.0_prec,y_interior_min-y,y-y_interior_max)

          ! Polynomial ramp normalised by pml_width so sigma -> sigma_max
          ! at the far edge of a layer of thickness pml_width.
          sx = sigma_max*(dx/pml_width)**p
          sy = sigma_max*(dy/pml_width)**p

          this%sigma_x%interior(i,j,iel,1) = sx
          this%sigma_y%interior(i,j,iel,1) = sy
        enddo
      enddo
    enddo

    call this%sigma_x%UpdateDevice()
    call this%sigma_y%UpdateDevice()

  endsubroutine SetPMLProfile_LinearEuler2D_PML_t

  pure function flux2d_LinearEuler2D_PML_t(this,s,dsdx) result(flux)
    !! Interior flux. Variables 1-5 use the parent linear Euler 2D
    !! flux; auxiliary variables 6-9 carry zero flux in both directions
    !! and are evolved purely by the PML source term.
    class(LinearEuler2D_PML_t),intent(in) :: this
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
    flux(5,1) = 0.0_prec ! sound speed, x flux; 0 (c held fixed in time)
    flux(5,2) = 0.0_prec ! sound speed, y flux; 0
    flux(6:9,1) = 0.0_prec ! PML auxiliaries evolve only via source term
    flux(6:9,2) = 0.0_prec
    if(.false.) flux(1,1) = flux(1,1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction flux2d_LinearEuler2D_PML_t

  pure function riemannflux2d_LinearEuler2D_PML_t(this,sL,sR,dsdx,nhat) result(flux)
    !! Impedance-matched Riemann flux for acoustic variables 1-5, with
    !! zero flux returned for the auxiliary variables 6-9. The acoustic
    !! formula is identical to the parent LinearEuler2D model; see
    !! riemannflux2d_LinearEuler2D_t for the derivation.
    class(LinearEuler2D_PML_t),intent(in) :: this
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
    flux(6:9) = 0.0_prec
    if(.false.) flux(1) = flux(1)+dsdx(1,1) ! suppress unused-dummy-argument warning

  endfunction riemannflux2d_LinearEuler2D_PML_t

  subroutine sourcemethod_LinearEuler2D_PML_t(this)
    !! Hu (2001) unsplit PML source term, evaluated per node. In the
    !! interior (sigma_x = sigma_y = 0) this leaves the acoustic
    !! variables untouched and the auxiliaries integrate q in time but
    !! never re-enter the dynamics (since the coupling coefficient
    !! sigma_x*sigma_y is zero). Inside the PML, the (sigma_x + sigma_y)
    !! damping plus the sigma_x*sigma_y*phi term produce the correct
    !! perfectly-matched behaviour for the linear Euler 2D system.
    !!
    !! Note: we override sourcemethod rather than source2d because the
    !! per-node sigma_x(i,j,iel), sigma_y(i,j,iel) lookups are not
    !! reachable from the pure source2d(s,dsdx) signature.
    implicit none
    class(LinearEuler2D_PML_t),intent(inout) :: this
    ! Local
    integer :: i,j,iel
    real(prec) :: sx,sy
    real(prec) :: s(1:this%nvar)

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem)

      sx = this%sigma_x%interior(i,j,iel,1)
      sy = this%sigma_y%interior(i,j,iel,1)
      s = this%solution%interior(i,j,iel,1:this%nvar)

      this%source%interior(i,j,iel,1) = -(sx+sy)*s(1)-sx*sy*s(6)
      this%source%interior(i,j,iel,2) = -(sx+sy)*s(2)-sx*sy*s(7)
      this%source%interior(i,j,iel,3) = -(sx+sy)*s(3)-sx*sy*s(8)
      this%source%interior(i,j,iel,4) = -(sx+sy)*s(4)-sx*sy*s(9)
      this%source%interior(i,j,iel,5) = 0.0_prec
      this%source%interior(i,j,iel,6) = s(1)
      this%source%interior(i,j,iel,7) = s(2)
      this%source%interior(i,j,iel,8) = s(3)
      this%source%interior(i,j,iel,9) = s(4)

    enddo

  endsubroutine sourcemethod_LinearEuler2D_PML_t

  subroutine hbc2d_NoNormalFlow_LinearEuler2D_PML(bc,mymodel)
    !! No-normal-flow BC for the PML-augmented linear Euler model.
    !! Variables 1-5 are treated identically to the parent
    !! LinearEuler2D no-normal-flow BC. Auxiliary variables 6-9 are
    !! given a zero exterior state; since they carry zero Riemann flux
    !! the exterior value is mathematically irrelevant, but zeroing it
    !! keeps the boundary state clean for diagnostics.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel
    ! Local
    integer :: n,i,iEl,j
    real(prec) :: nhat(1:2),s(1:5)

    select type(m => mymodel)
    class is(LinearEuler2D_PML_t)
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
          m%solution%extBoundary(i,j,iEl,6:9) = 0.0_prec ! PML auxiliaries
        enddo
      enddo
    endselect

  endsubroutine hbc2d_NoNormalFlow_LinearEuler2D_PML

  subroutine hbc2d_Radiation_LinearEuler2D_PML(bc,mymodel)
    !! Radiation BC for the PML-augmented linear Euler model: zero
    !! acoustic perturbation in the exterior state, sound speed copied
    !! from interior so the Riemann solver sees a consistent c, and
    !! auxiliary variables set to zero.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel
    ! Local
    integer :: n,i,iEl,j

    select type(m => mymodel)
    class is(LinearEuler2D_PML_t)
      do n = 1,bc%nBoundaries
        iEl = bc%elements(n)
        j = bc%sides(n)
        do i = 1,m%solution%interp%N+1
          m%solution%extBoundary(i,j,iEl,1) = 0.0_prec
          m%solution%extBoundary(i,j,iEl,2) = 0.0_prec
          m%solution%extBoundary(i,j,iEl,3) = 0.0_prec
          m%solution%extBoundary(i,j,iEl,4) = 0.0_prec
          m%solution%extBoundary(i,j,iEl,5) = m%solution%boundary(i,j,iEl,5) ! c preserved
          m%solution%extBoundary(i,j,iEl,6:9) = 0.0_prec
        enddo
      enddo
    endselect

  endsubroutine hbc2d_Radiation_LinearEuler2D_PML

endmodule self_LinearEuler2D_PML_t
