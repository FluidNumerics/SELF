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

module self_NullDGModel1D_t

  use self_model
  use self_dgmodel1d
  use self_mesh

  implicit none

  type,extends(dgmodel1d) :: NullDGModel1D_t
    ! Add any additional attributes here that are specific to your model

  contains
    !   procedure :: SetNumberOfVariables => SetNumberOfVariables_NullDGModel1D_t
    !   procedure :: hbc1d_Prescribed => hbc1d_Prescribed_Model
    !   procedure :: hbc1d_Radiation => hbc1d_Generic_Model
    !   procedure :: hbc1d_NoNormalFlow => hbc1d_Generic_Model
    !   procedure :: pbc1d_Prescribed => pbc1d_Prescribed_Model
    !   procedure :: pbc1d_Radiation => pbc1d_Generic_Model
    !   procedure :: pbc1d_NoNormalFlow => pbc1d_Generic_Model
    !   procedure :: SetMetadata => SetMetadata_NullDGModel1D_t
    !   procedure :: pretendency => pretendency_NullDGModel1D_t
    !   procedure :: entropy_func => entropy_func_NullDGModel1D_t
    !   procedure :: flux1d => flux1d_NullDGModel1D_t
    !   procedure :: riemannflux1d => riemannflux1d_NullDGModel1D_t
    !   procedure :: source1d => source1d_NullDGModel1D_t

  endtype NullDGModel1D_t

contains

  ! subroutine SetNumberOfVariables_NullDGModel1D_t(this)
  !   implicit none
  !   class(Model),intent(inout) :: this
  !     this%nvar = 1
  ! endsubroutine SetNumberOfVariables_NullDGModel1D_t

  ! subroutine SetMetadata_NullDGModel1D_t(this)
  !   implicit none
  !   class(NullDGModel1D_t),intent(inout) :: this
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

  ! endsubroutine SetMetadata_NullDGModel1D_t

  ! pure function entropy_func_NullDGModel1D_t(this,s) result(e)
  !   class(NullDGModel1D_t),intent(in) :: this
  !   real(prec),intent(in) :: s(1:this%solution%nvar)
  !   real(prec) :: e

  !   e = 0.0_prec

  ! endfunction entropy_func_NullDGModel1D_t

  ! pure function flux1d_NullDGModel1D_t(this,s,dsdx) result(flux)
  !   class(NullDGModel1D_t),intent(in) :: this
  !   real(prec),intent(in) :: s(1:this%solution%nvar)
  !   real(prec),intent(in) :: dsdx(1:this%solution%nvar)
  !   real(prec) :: flux(1:this%solution%nvar)

  !   flux(1:this%nvar) = 0.0_prec

  ! endfunction flux1d_NullDGModel1D_t

  ! pure function riemannflux1d_NullDGModel1D_t(this,sL,sR,dsdx,nhat) result(flux)
  !   class(NullDGModel1D_t),intent(in) :: this
  !   real(prec),intent(in) :: sL(1:this%solution%nvar)
  !   real(prec),intent(in) :: sR(1:this%solution%nvar)
  !   real(prec),intent(in) :: dsdx(1:this%solution%nvar)
  !   real(prec),intent(in) :: nhat
  !   real(prec) :: flux(1:this%solution%nvar)

  !   flux(1:this%nvar) = 0.0_prec

  ! endfunction riemannflux1d_NullDGModel1D_t

  ! subroutine PreTendency_NulDGModel1D_t(this)
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

  ! endsubroutine PreTendency_NulDGModel1D_t

  ! pure function source1d_NullDGModel1D_t(this,s,dsdx) result(source)
  !   class(NullDGModel1D_t),intent(in) :: this
  !   real(prec),intent(in) :: s(1:this%solution%nvar)
  !   real(prec),intent(in) :: dsdx(1:this%solution%nvar)
  !   real(prec) :: source(1:this%solution%nvar)

  !   source(1:this%nvar) = 0.0_prec

  ! endfunction source1d_NullDGModel1D_t

endmodule self_NullDGModel1D_t
