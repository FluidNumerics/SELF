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

module SELF_ECDGModel2D_t
  !! Entropy-conserving DGSEM base class for 2-D conservation laws.
  !!
  !! Extends DGModel2D_t by replacing the standard volume flux divergence with
  !! the split-form (two-point) EC volume term following
  !!   Gassner, Winters, Kopriva (2016)
  !!   Winters, Kopriva, Gassner, Hindenlang (2021)
  !!
  !! Semidiscretization:
  !!   du/dt = source - (1/J) * [ 2 sum_n D_{n,i} F~^1(n,i,j)
  !!                             + 2 sum_n D_{n,j} F~^2(n,i,j)
  !!                             + M^{-1} B^T f_Riemann ]
  !!
  !! where F~^r is the Jacobian-weighted contravariant EC two-point flux in the
  !! r-th computational direction (computed by TwoPointFluxMethod via twopointflux2d),
  !! and f_Riemann is a standard Riemann solver flux at element faces (computed
  !! by BoundaryFlux via riemannflux2d).

  use SELF_DGModel2D_t
  use SELF_MappedTwoPointVector_2D

  implicit none

  type,extends(DGModel2D_t) :: ECDGModel2D_t

    type(MappedTwoPointVector2D) :: twoPointFlux

  contains

    procedure :: Init => Init_ECDGModel2D_t
    procedure :: Free => Free_ECDGModel2D_t

    procedure :: twopointflux2d => twopointflux2d_ECDGModel2D_t
    procedure :: TwoPointFluxMethod => TwoPointFluxMethod_ECDGModel2D_t
    procedure :: CalculateTendency => CalculateTendency_ECDGModel2D_t

  endtype ECDGModel2D_t

contains

  subroutine Init_ECDGModel2D_t(this,mesh,geometry)
    implicit none
    class(ECDGModel2D_t),intent(out) :: this
    type(Mesh2D),intent(in),target :: mesh
    type(SEMQuad),intent(in),target :: geometry

    ! Initialise all parent fields (solution, flux, source, ...)
    call Init_DGModel2D_t(this,mesh,geometry)

    ! Additional two-point flux field
    call this%twoPointFlux%Init(geometry%x%interp,this%nvar,mesh%nElem)
    call this%twoPointFlux%AssociateGeometry(geometry)

  endsubroutine Init_ECDGModel2D_t

  subroutine Free_ECDGModel2D_t(this)
    implicit none
    class(ECDGModel2D_t),intent(inout) :: this

    call this%twoPointFlux%DissociateGeometry()
    call this%twoPointFlux%Free()
    call Free_DGModel2D_t(this)

  endsubroutine Free_ECDGModel2D_t

  pure function twopointflux2d_ECDGModel2D_t(this,sL,sR) result(flux)
    !! Entropy-conserving two-point flux function.
    !!
    !! Returns the physical-space two-point EC flux for the state pair (sL, sR).
    !! flux(ivar, d) is the d-th physical component (d=1: x, d=2: y) for variable ivar.
    !!
    !! Override this procedure in a concrete model to supply an EC flux.
    !! The stub returns zero (no flux), which is safe but trivial.
    implicit none
    class(ECDGModel2D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec) :: flux(1:this%nvar,1:2)

    flux = 0.0_prec
    if(.false.) flux(1,1) = flux(1,1)+sL(1)+sR(1) ! suppress unused-argument warning

  endfunction twopointflux2d_ECDGModel2D_t

  subroutine TwoPointFluxMethod_ECDGModel2D_t(this)
    !! Computes pre-projected SCALAR contravariant two-point fluxes for all
    !! node pairs and stores them in twoPointFlux%interior(n,i,j,iel,ivar,r).
    !!
    !! Following Trixi.jl for curved meshes, each computational direction r
    !! uses the correct partner node AND the correct averaged metric Ja^r:
    !!
    !!   r=1: pair (i,j)-(nn,j), Ja^1_avg = 0.5*(Ja^1(i,j) + Ja^1(nn,j))
    !!        interior(nn,i,j,...,1) = sum_d Ja^1_avg(d) * F_EC_d(sL, sR)
    !!
    !!   r=2: pair (i,j)-(i,nn), Ja^2_avg = 0.5*(Ja^2(i,j) + Ja^2(i,nn))
    !!        interior(nn,i,j,...,2) = sum_d Ja^2_avg(d) * F_EC_d(sL, sR)
    !!
    !! The result is a SCALAR per (nn,i,j,iel,ivar,r) — no cross-contamination
    !! between directions.
    implicit none
    class(ECDGModel2D_t),intent(inout) :: this
    ! Local
    integer :: nn,i,j,d,iEl,iVar
    real(prec) :: sL(1:this%nvar),sR(1:this%nvar)
    real(prec) :: Fphys(1:this%nvar,1:2)
    real(prec) :: Fc

    do concurrent(nn=1:this%solution%N+1,i=1:this%solution%N+1, &
                  j=1:this%solution%N+1,iEl=1:this%mesh%nElem)

      sL = this%solution%interior(i,j,iEl,1:this%nvar)

      ! xi^1: pair (i,j)-(nn,j), project onto avg(Ja^1)
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

      ! xi^2: pair (i,j)-(i,nn), project onto avg(Ja^2)
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

  endsubroutine TwoPointFluxMethod_ECDGModel2D_t

  subroutine CalculateTendency_ECDGModel2D_t(this)
    !! Computes du/dt = source - EC-DG flux divergence.
    !!
    !! Volume term : MappedDivergence of twoPointFlux (EC split-form sum, /J)
    !! Surface term: (1/J) * M^{-1} B^T f_Riemann  (weak-form Riemann flux)
    implicit none
    class(ECDGModel2D_t),intent(inout) :: this
    ! Local
    integer :: i,j,iEl,iVar
    real(prec) :: bM1,bM2,qwi,qwj,jac

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
    call this%BoundaryFlux() ! fills flux%boundaryNormal with Riemann fluxes

    ! EC volume divergence (includes J weighting via MappedDivergence)
    call this%TwoPointFluxMethod()
    call this%twoPointFlux%UpdateDevice()
    call this%twoPointFlux%MappedDivergence(this%fluxDivergence%interior)

    ! Add weak-form DG surface term: (1/J) * M^{-1} B^T f_Riemann
    ! Boundary side ordering: 1=South, 2=East, 3=North, 4=West
    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iEl=1:this%mesh%nElem,iVar=1:this%solution%nVar)

      bM1 = this%solution%interp%bMatrix(i,2) ! right/north boundary value at node i
      bM2 = this%solution%interp%bMatrix(i,1) ! left/south boundary value at node i
      qwi = this%solution%interp%qWeights(i)
      qwj = this%solution%interp%qWeights(j)
      jac = this%geometry%J%interior(i,j,iEl,1)

      ! East (side 2) and West (side 4) contributions — vary in xi^1 (i-direction)
      this%fluxDivergence%interior(i,j,iEl,iVar) = &
        this%fluxDivergence%interior(i,j,iEl,iVar)+ &
        (bM1*this%flux%boundaryNormal(j,2,iEl,iVar)+ &
         bM2*this%flux%boundaryNormal(j,4,iEl,iVar))/(qwi*jac)+ &
        (this%solution%interp%bMatrix(j,2)*this%flux%boundaryNormal(i,3,iEl,iVar)+ &
         this%solution%interp%bMatrix(j,1)*this%flux%boundaryNormal(i,1,iEl,iVar))/(qwj*jac)

    enddo

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iEl=1:this%mesh%nElem,iVar=1:this%solution%nVar)

      this%dSdt%interior(i,j,iEl,iVar) = &
        this%source%interior(i,j,iEl,iVar)- &
        this%fluxDivergence%interior(i,j,iEl,iVar)

    enddo

  endsubroutine CalculateTendency_ECDGModel2D_t

endmodule SELF_ECDGModel2D_t
