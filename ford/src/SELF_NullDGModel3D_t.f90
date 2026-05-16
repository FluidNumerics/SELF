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

module self_NullDGModel3D_t

  use self_model
  use self_dgmodel3d
  use self_mesh

  implicit none

  type,extends(dgmodel3d) :: NullDGModel3D_t
    ! Add any additional attributes here that are specific to your model

  contains
    !   procedure :: SetNumberOfVariables => SetNumberOfVariables_NullDGModel3D_t
    !   procedure :: hbc3d_Prescribed => hbc3d_Generic_Model
    !   procedure :: hbc3d_Radiation => hbc3d_Generic_Model
    !   procedure :: hbc3d_NoNormalFlow => hbc3d_Generic_Model
    !   procedure :: pbc3d_Prescribed => pbc3d_Generic_Model
    !   procedure :: pbc3d_Radiation => pbc3d_Generic_Model
    !   procedure :: pbc3d_NoNormalFlow => pbc3d_Generic_Model
    !   procedure :: SetMetadata => SetMetadata_NullDGModel3D_t
    !   procedure :: pretendency => pretendency_NullDGModel3D_t
    !   procedure :: entropy_func => entropy_func_NullDGModel3D_t
    !   procedure :: flux3d => flux3d_NullDGModel3D_t
    !   procedure :: riemannflux3d => riemannflux3d_NullDGModel3D_t
    !   procedure :: source3d => source3d_NullDGModel3D_t

  endtype NullDGModel3D_t

contains

  ! subroutine SetNumberOfVariables_NullDGModel3D_t(this)
  !   implicit none
  !   class(Model),intent(inout) :: this
  !     this%nvar = 1
  ! endsubroutine SetNumberOfVariables_NullDGModel3D_t

  ! subroutine SetMetadata_NullDGModel3D_t(this)
  !   implicit none
  !   class(NullDGModel3D_t),intent(inout) :: this
  !   ! Local
  !   integer :: ivar
  !   character(LEN=3) :: ivarChar
  !   character(LEN=25) :: varname

  !   do ivar = 1,this%nvar
  !     write(ivarChar,'(I3.3)') ivar
  !     varname = "solution"//trim(ivarChar)
  !     call this%solution%SetName(ivar,varname)
  !     call this%solution%SetUnits(ivar,"[null]")
  !   enddo

  ! endsubroutine SetMetadata_NullDGModel3D_t

  !   pure function bcGeneric_NullDGModel3D_t(this,s) result(exts)
  !   class(NullDGModel3D_t),intent(in) :: this
  !   real(prec),intent(in) :: s(1:this%nvar)
  !   real(prec) :: exts(1:this%nvar)
  !   ! Local
  !   integer :: ivar

  !   do ivar = 1,this%nvar
  !     exts(ivar) = 0.0_prec
  !   enddo

  ! endfunction bcGeneric_NullDGModel3D_t

  ! pure function bcGrad3dGeneric_NullDGModel3D_t(this,dsdx) result(extDsdx)
  !   class(NullDGModel3D_t),intent(in) :: this
  !   real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
  !   real(prec) :: extDsdx(1:this%nvar,1:3)
  !   ! Local
  !   integer :: ivar

  !   do ivar = 1,this%nvar
  !     extDsdx(ivar,1:3) = dsdx(ivar,1:3)
  !   enddo

  ! endfunction bcGrad3dGeneric_NullDGModel3D_t
  ! pure function entropy_func_NullDGModel3D_t(this,s) result(e)
  !   class(NullDGModel3D_t),intent(in) :: this
  !   real(prec),intent(in) :: s(1:this%solution%nvar)
  !   real(prec) :: e

  !   e = 0.0_prec

  ! endfunction entropy_func_NullDGModel3D_t

  ! pure function flux3d_NullDGModel3D_t(this,s,dsdx) result(flux)
  !   class(NullDGModel3D_t),intent(in) :: this
  !   real(prec),intent(in) :: s(1:this%solution%nvar)
  !   real(prec),intent(in) :: dsdx(1:this%solution%nvar,1:3)
  !   real(prec) :: flux(1:this%solution%nvar,1:3)

  !   flux(1:this%nvar,1:3) = 0.0_prec

  ! endfunction flux3d_NullDGModel3D_t

  ! pure function riemannflux3d_NullDGModel3D_t(this,sL,sR,dsdx,nhat) result(flux)
  !   class(NullDGModel3D_t),intent(in) :: this
  !   real(prec),intent(in) :: sL(1:this%solution%nvar)
  !   real(prec),intent(in) :: sR(1:this%solution%nvar)
  !   real(prec),intent(in) :: dsdx(1:this%solution%nvar,1:3)
  !   real(prec),intent(in) :: nhat(1:3)
  !   real(prec) :: flux(1:this%solution%nvar)

  !   flux(1:this%nvar) = 0.0_prec

  ! endfunction riemannflux3d_NullDGModel3D_t

  ! subroutine PreTendency_NulDGModel3D_t(this)
  ! !! PreTendency is a template routine that is used to house any additional calculations
  ! !! that you want to execute at the beginning of the tendency calculation routine.
  ! !! This default PreTendency simply returns back to the caller without executing any instructions
  ! !!
  ! !! The intention is to provide a method that can be overridden through type-extension, to handle
  ! !! any steps that need to be executed before proceeding with the usual tendency calculation methods.
  ! !!
  !   implicit none
  !   class(Model),intent(inout) :: this

  !   return

  ! endsubroutine PreTendency_NulDGModel3D_t

  ! pure function source3d_NullDGModel3D_t(this,s,dsdx) result(source)
  !   class(NullDGModel3D_t),intent(in) :: this
  !   real(prec),intent(in) :: s(1:this%solution%nvar)
  !   real(prec),intent(in) :: dsdx(1:this%solution%nvar,1:3)
  !   real(prec) :: source(1:this%solution%nvar)

  !   source(1:this%nvar) = 0.0_prec

  ! endfunction source3d_NullDGModel3D_t

endmodule self_NullDGModel3D_t
