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

module SELF_MappedScalar_1D_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Scalar_1D
  use SELF_Mesh_1D
  use SELF_Geometry_1D
  use SELF_HDF5
  use HDF5

  use FEQParse

  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type,extends(Scalar1D),public :: MappedScalar1D_t
    logical :: geometry_associated = .false.
    type(Geometry1D),pointer :: geometry => null()

  contains
    procedure,public :: AssociateGeometry => AssociateGeometry_MappedScalar1D_t
    procedure,public :: DissociateGeometry => DissociateGeometry_MappedScalar1D_t
    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar1D_t

    procedure,public :: SideExchange => SideExchange_MappedScalar1D_t

    generic,public :: MappedDerivative => MappedDerivative_MappedScalar1D_t
    procedure,private :: MappedDerivative_MappedScalar1D_t

    generic,public :: MappedDGDerivative => MappedDGDerivative_MappedScalar1D_t
    procedure,private :: MappedDGDerivative_MappedScalar1D_t

  endtype MappedScalar1D_t

contains

! ---------------------- Scalars ---------------------- !

  subroutine AssociateGeometry_MappedScalar1D_t(this,geometry)
    implicit none
    class(MappedScalar1D_t),intent(inout) :: this
    type(Geometry1D),target,intent(in) :: geometry

    if(.not. associated(this%geometry)) then
      this%geometry => geometry
      this%geometry_associated = .true.
    endif

  endsubroutine AssociateGeometry_MappedScalar1D_t

  subroutine DissociateGeometry_MappedScalar1D_t(this)
    implicit none
    class(MappedScalar1D_t),intent(inout) :: this

    if(associated(this%geometry)) then
      this%geometry => null()
      this%geometry_associated = .false.
    endif

  endsubroutine DissociateGeometry_MappedScalar1D_t

  subroutine SetInteriorFromEquation_MappedScalar1D_t(this,time)
    !!  Sets the this % interior attribute using the eqn attribute,
    !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar1D_t),intent(inout) :: this
    real(prec),intent(in) :: time
    ! Local
    integer :: iVar

    do ivar = 1,this%nvar
      this%interior(:,:,ivar) = this%eqn(ivar)%evaluate(this%geometry%x%interior)
    enddo

  endsubroutine SetInteriorFromEquation_MappedScalar1D_t

  subroutine SideExchange_MappedScalar1D_t(this,mesh)
    implicit none
    class(MappedScalar1D_t),intent(inout) :: this
    type(Mesh1D),intent(inout) :: mesh
    ! Local
    integer :: e1,e2,s1,s2
    integer :: ivar

    do concurrent(e1=1:mesh%nElem,ivar=1:this%nvar)

      if(e1 == 1) then

        s1 = 2
        e2 = e1+1
        s2 = 1
        this%extBoundary(s1,e1,ivar) = this%boundary(s2,e2,ivar)

      elseif(e1 == mesh%nElem) then

        s1 = 1
        e2 = e1-1
        s2 = 2
        this%extBoundary(s1,e1,ivar) = this%boundary(s2,e2,ivar)

      else

        s1 = 1
        e2 = e1-1
        s2 = 2
        this%extBoundary(s1,e1,ivar) = this%boundary(s2,e2,ivar)

        s1 = 2
        e2 = e1+1
        s2 = 1
        this%extBoundary(s1,e1,ivar) = this%boundary(s2,e2,ivar)

      endif

    enddo

  endsubroutine SideExchange_MappedScalar1D_t

  subroutine MappedDerivative_MappedScalar1D_t(this,dF)
    implicit none
    class(MappedScalar1D_t),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,ii
    real(prec) :: dfloc

    do concurrent(i=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar)

      dfloc = 0.0_prec
      do ii = 1,this%N+1
        dfloc = dfloc+this%interp%dMatrix(ii,i)*this%interior(ii,iel,ivar)
      enddo
      df(i,iel,ivar) = dfloc/this%geometry%dxds%interior(i,iEl,1)

    enddo

  endsubroutine MappedDerivative_MappedScalar1D_t

  subroutine MappedDGDerivative_MappedScalar1D_t(this,dF)
    implicit none
    class(MappedScalar1D_t),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,ii
    real(prec) :: dfloc

    do concurrent(i=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar)

      dfloc = 0.0_prec
      do ii = 1,this%N+1
        dfloc = dfloc+this%interp%dgMatrix(ii,i)*this%interior(ii,iel,ivar)
      enddo

      dfloc = dfloc+(this%boundarynormal(2,iel,ivar)*this%interp%bMatrix(i,2)+ &
                     this%boundarynormal(1,iel,ivar)*this%interp%bMatrix(i,1))/ &
              this%interp%qWeights(i)

      df(i,iel,ivar) = dfloc/this%geometry%dxds%interior(i,iEl,1)

    enddo

  endsubroutine MappedDGDerivative_MappedScalar1D_t

endmodule SELF_MappedScalar_1D_t
