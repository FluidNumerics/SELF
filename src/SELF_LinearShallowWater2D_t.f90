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
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module self_LinearShallowWater2D_t
  use self_model
  use self_dgmodel2d
  use self_mesh

  implicit none

  type,extends(dgmodel2d) :: LinearShallowWater2D_t
    real(prec) :: H = 0.0_prec ! uniform resting depth
    real(prec) :: g = 0.0_prec ! acceleration due to gravity
    type(MappedScalar2D) :: fCori ! The coriolis parameter

  contains
    procedure :: AdditionalInit => AdditionalInit_LinearShallowWater2D_t
    procedure :: AdditionalFree => AdditionalFree_LinearShallowWater2D_t
    procedure :: SetNumberOfVariables => SetNumberOfVariables_LinearShallowWater2D_t
    procedure :: SetMetadata => SetMetadata_LinearShallowWater2D_t
    procedure :: entropy_func => entropy_func_LinearShallowWater2D_t
    procedure :: flux2d => flux2d_LinearShallowWater2D_t
    procedure :: riemannflux2d => riemannflux2d_LinearShallowWater2D_t
    procedure :: hbc2d_NoNormalFlow => hbc2d_NoNormalFlow_LinearShallowWater2D_t
    procedure :: sourcemethod => sourcemethod_LinearShallowWater2D_t
    ! Custom methods
    generic,public :: SetCoriolis => SetCoriolis_fplane_LinearShallowWater2D_t, &
      SetCoriolis_betaplane_LinearShallowWater2D_t
    procedure,private :: SetCoriolis_fplane_LinearShallowWater2D_t
    procedure,private :: SetCoriolis_betaplane_LinearShallowWater2D_t

    procedure,public :: DiagnoseGeostrophicVelocity => DiagnoseGeostrophicVelocity_LinearShallowWater2D_t

  endtype LinearShallowWater2D_t

contains

  subroutine AdditionalInit_LinearShallowWater2D_t(this)
    implicit none
    class(LinearShallowWater2D_t),intent(inout) :: this

    call this%fCori%Init(this%geometry%x%interp, &
                         1,this%mesh%nElem)

  endsubroutine AdditionalInit_LinearShallowWater2D_t

  subroutine AdditionalFree_LinearShallowWater2D_t(this)
    implicit none
    class(LinearShallowWater2D_t),intent(inout) :: this

    call this%fCori%Free()

  endsubroutine AdditionalFree_LinearShallowWater2D_t

  subroutine SetNumberOfVariables_LinearShallowWater2D_t(this)
    implicit none
    class(LinearShallowWater2D_t),intent(inout) :: this

    this%nvar = 3

  endsubroutine SetNumberOfVariables_LinearShallowWater2D_t

  subroutine SetCoriolis_fplane_LinearShallowWater2D_t(this,f0)
    implicit none
    class(LinearShallowWater2D_t),intent(inout) :: this
    real(prec),intent(in) :: f0
    ! Local
    integer :: iel
    integer :: i
    integer :: j

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem)
      this%fCori%interior(i,j,iel,1) = f0
    enddo
    call this%fCori%UpdateDevice()

  endsubroutine SetCoriolis_fplane_LinearShallowWater2D_t

  subroutine SetCoriolis_betaplane_LinearShallowWater2D_t(this,f0,beta)
    implicit none
    class(LinearShallowWater2D_t),intent(inout) :: this
    real(prec),intent(in) :: f0
    real(prec),intent(in) :: beta
    ! Local
    integer :: iel
    integer :: i
    integer :: j
    real(prec) :: y

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem)
      y = this%geometry%x%interior(i,j,iel,1,2)
      this%fCori%interior(i,j,iel,1) = f0+beta*y
    enddo
    call this%fCori%UpdateDevice()

  endsubroutine SetCoriolis_betaplane_LinearShallowWater2D_t

  subroutine DiagnoseGeostrophicVelocity_LinearShallowWater2D_t(this)
    implicit none
    class(LinearShallowWater2D_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: i
    integer :: j
    real(prec) :: dpdx,dpdy,f

    ! We assume here that the velocity field is identically zero
    ! everywhere and the only field that is set is the free surface height
    ! with a non-zero coriolis parameter.
    ! In this case, we have that the tendency calculation will give
    ! the gradient in the free surface, consistent with the DG approximation
    this%solution%interior(:,:,:,1) = 0.0_prec ! Set u=0
    this%solution%interior(:,:,:,2) = 0.0_prec ! Set v=0
    call this%solution%UpdateDevice()
    call this%CalculateTendency()
    call this%dSdt%UpdateHost()

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem)

      dpdx = -this%dSdt%interior(i,j,iel,1)
      dpdy = -this%dSdt%interior(i,j,iel,2)
      f = this%fCori%interior(i,j,iel,1)
      this%solution%interior(i,j,iel,1) = -dpdy/f ! u
      this%solution%interior(i,j,iel,2) = dpdx/f ! v
    enddo

    call this%solution%UpdateDevice()

  endsubroutine DiagnoseGeostrophicVelocity_LinearShallowWater2D_t

  subroutine SetMetadata_LinearShallowWater2D_t(this)
    implicit none
    class(LinearShallowWater2D_t),intent(inout) :: this

    call this%solution%SetName(1,"u")
    call this%solution%SetUnits(1,"m/s")
    call this%solution%SetName(2,"v")
    call this%solution%SetUnits(2,"m/s")
    call this%solution%SetName(3,"eta")
    call this%solution%SetUnits(3,"m")
    call this%fCori%SetName(1,"f")
    call this%fCori%SetUnits(1,"1/s")

  endsubroutine SetMetadata_LinearShallowWater2D_t

  pure function entropy_func_LinearShallowWater2D_t(this,s) result(e)
    class(LinearShallowWater2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%solution%nvar)
    real(prec) :: e

    e = 0.5_prec*(this%H*s(1)*s(1)+ &
                  this%H*s(2)*s(2)+ &
                  this%g*s(3)*s(3))

  endfunction entropy_func_LinearShallowWater2D_t

  pure function flux2d_LinearShallowWater2D_t(this,s,dsdx) result(flux)
    class(LinearShallowWater2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%solution%nvar)
    real(prec),intent(in) :: dsdx(1:this%solution%nvar,1:2)
    real(prec) :: flux(1:this%solution%nvar,1:2)

    flux(1,1) = this%g*s(3)
    flux(1,2) = 0.0_prec
    flux(2,1) = 0.0_prec
    flux(2,2) = this%g*s(3)
    flux(3,1) = this%H*s(1)
    flux(3,2) = this%H*s(2)

  endfunction flux2d_LinearShallowWater2D_t

  pure function riemannflux2d_LinearShallowWater2D_t(this,sL,sR,dsdx,nhat) result(flux)
    class(LinearShallowWater2D_t),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%solution%nVar)
    real(prec),intent(in) :: sR(1:this%solution%nVar)
    real(prec),intent(in) :: dsdx(1:this%solution%nVar,1:2)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: flux(1:this%solution%nVar)
    ! Local
    real(prec) :: c
    real(prec) :: unL
    real(prec) :: unR

    c = sqrt(this%g*this%H)

    unL = sL(1)*nhat(1)+sL(2)*nhat(2)
    unR = sR(1)*nhat(1)+sR(2)*nhat(2)

    flux(1) = 0.5_prec*(this%g*(sL(3)+sR(3))+c*(unL-unR))*nhat(1)
    flux(2) = 0.5_prec*(this%g*(sL(3)+sR(3))+c*(unL-unR))*nhat(2)
    flux(3) = 0.5_prec*(this%H*(unL+unR)+c*(sL(3)-sR(3)))

  endfunction riemannflux2d_LinearShallowWater2D_t

  pure function hbc2d_NoNormalFlow_LinearShallowWater2D_t(this,s,nhat) result(exts)
    class(LinearShallowWater2D_t),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: exts(1:this%nvar)
    ! Local
    integer :: ivar

    exts(1) = (nhat(2)**2-nhat(1)**2)*s(1)-2.0_prec*nhat(1)*nhat(2)*s(2) ! u
    exts(2) = (nhat(1)**2-nhat(2)**2)*s(2)-2.0_prec*nhat(1)*nhat(2)*s(1) ! v
    exts(3) = s(3) ! eta

  endfunction hbc2d_NoNormalFlow_LinearShallowWater2D_t

  subroutine sourcemethod_LinearShallowWater2D_t(this)
    implicit none
    class(LinearShallowWater2D_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: i
    integer :: j
    real(prec) :: s(1:this%nvar)

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem)

      s = this%solution%interior(i,j,iel,1:this%nvar)

      this%source%interior(i,j,iel,1) = this%fCori%interior(i,j,iel,1)*s(2) ! du/dt = f*v
      this%source%interior(i,j,iel,2) = -this%fCori%interior(i,j,iel,1)*s(1) ! dv/dt = -f*u

    enddo

  endsubroutine sourcemethod_LinearShallowWater2D_t

endmodule self_LinearShallowWater2D_t
