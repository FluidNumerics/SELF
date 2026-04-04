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

module SELF_ECDGModel3D_t
  !! Entropy-conserving DGSEM base class for 3-D conservation laws.
  !! See SELF_ECDGModel2D_t for design documentation.

  use SELF_DGModel3D_t
  use SELF_MappedTwoPointVector_3D

  implicit none

  type,extends(DGModel3D_t) :: ECDGModel3D_t

    type(MappedTwoPointVector3D) :: twoPointFlux

  contains

    procedure :: Init => Init_ECDGModel3D_t
    procedure :: Free => Free_ECDGModel3D_t

    procedure :: twopointflux3d => twopointflux3d_ECDGModel3D_t
    procedure :: TwoPointFluxMethod => TwoPointFluxMethod_ECDGModel3D_t
    procedure :: CalculateTendency => CalculateTendency_ECDGModel3D_t

  endtype ECDGModel3D_t

contains

  subroutine Init_ECDGModel3D_t(this,mesh,geometry)
    implicit none
    class(ECDGModel3D_t),intent(out) :: this
    type(Mesh3D),intent(in),target :: mesh
    type(SEMHex),intent(in),target :: geometry

    call Init_DGModel3D_t(this,mesh,geometry)

    call this%twoPointFlux%Init(geometry%x%interp,this%nvar,mesh%nElem)
    call this%twoPointFlux%AssociateGeometry(geometry)

  endsubroutine Init_ECDGModel3D_t

  subroutine Free_ECDGModel3D_t(this)
    implicit none
    class(ECDGModel3D_t),intent(inout) :: this

    call this%twoPointFlux%DissociateGeometry()
    call this%twoPointFlux%Free()
    call Free_DGModel3D_t(this)

  endsubroutine Free_ECDGModel3D_t

  pure function twopointflux3d_ECDGModel3D_t(this,sL,sR) result(flux)
    !! Entropy-conserving two-point flux function (3-D).
    !! flux(ivar, d) is the d-th physical component (d=1:x, d=2:y, d=3:z).
    !! Override in concrete models. Stub returns zero.
    implicit none
    class(ECDGModel3D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec) :: flux(1:this%nvar,1:3)

    flux = 0.0_prec
    if(.false.) flux(1,1) = flux(1,1)+sL(1)+sR(1) ! suppress unused-argument warning

  endfunction twopointflux3d_ECDGModel3D_t

  subroutine TwoPointFluxMethod_ECDGModel3D_t(this)
    !! Fills twoPointFlux%interior(n,i,j,k,iel,ivar,d) with the physical-space
    !! EC two-point fluxes between nodes (i,j,k) and (n,j,k).
    implicit none
    class(ECDGModel3D_t),intent(inout) :: this
    ! Local
    integer :: nn,i,j,k,iEl
    real(prec) :: sL(1:this%nvar),sR(1:this%nvar)
    real(prec) :: F(1:this%nvar,1:3)

    do concurrent(nn=1:this%solution%N+1,i=1:this%solution%N+1, &
                  j=1:this%solution%N+1,k=1:this%solution%N+1, &
                  iEl=1:this%mesh%nElem)

      sL = this%solution%interior(i,j,k,iEl,1:this%nvar)
      sR = this%solution%interior(nn,j,k,iEl,1:this%nvar)
      F = this%twopointflux3d(sL,sR)
      this%twoPointFlux%interior(nn,i,j,k,iEl,1:this%nvar,1) = F(1:this%nvar,1)
      this%twoPointFlux%interior(nn,i,j,k,iEl,1:this%nvar,2) = F(1:this%nvar,2)
      this%twoPointFlux%interior(nn,i,j,k,iEl,1:this%nvar,3) = F(1:this%nvar,3)

    enddo

  endsubroutine TwoPointFluxMethod_ECDGModel3D_t

  subroutine CalculateTendency_ECDGModel3D_t(this)
    implicit none
    class(ECDGModel3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,k,iEl,iVar
    real(prec) :: bMi1,bMi2,bMj1,bMj2,bMk1,bMk2,qwi,qwj,qwk,jac

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

    ! Add (1/J) * M^{-1} B^T f_Riemann
    ! 3D boundary side ordering: 1=Bottom, 2=South, 3=East, 4=North, 5=West, 6=Top
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
        (bMk1*this%flux%boundaryNormal(i,j,6,iEl,iVar)+ & ! top
         bMk2*this%flux%boundaryNormal(i,j,1,iEl,iVar))/(qwk*jac)+ & ! bottom
        (bMi1*this%flux%boundaryNormal(j,k,3,iEl,iVar)+ & ! east
         bMi2*this%flux%boundaryNormal(j,k,5,iEl,iVar))/(qwi*jac)+ & ! west
        (bMj1*this%flux%boundaryNormal(i,k,4,iEl,iVar)+ & ! north
         bMj2*this%flux%boundaryNormal(i,k,2,iEl,iVar))/(qwj*jac) ! south

    enddo

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%solution%N+1,iEl=1:this%mesh%nElem, &
                  iVar=1:this%solution%nVar)

      this%dSdt%interior(i,j,k,iEl,iVar) = &
        this%source%interior(i,j,k,iEl,iVar)- &
        this%fluxDivergence%interior(i,j,k,iEl,iVar)

    enddo

  endsubroutine CalculateTendency_ECDGModel3D_t

endmodule SELF_ECDGModel3D_t
