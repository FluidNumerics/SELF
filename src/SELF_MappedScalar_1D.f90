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

module SELF_MappedScalar_1D

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Scalar_1D
  use SELF_Mesh_1D
  use SELF_Geometry_1D
  use SELF_MPI
  use SELF_HDF5
  use HDF5

  use FEQParse

  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type,extends(Scalar1D),public :: MappedScalar1D

  contains
    procedure,public :: SideExchange => SideExchange_MappedScalar1D
    procedure,public :: BassiRebaySides => BassiRebaySides_MappedScalar1D

    generic,public :: Derivative => Derivative_MappedScalar1D
    procedure,private :: Derivative_MappedScalar1D
    generic,public :: DGDerivative => DGDerivative_MappedScalar1D
    procedure,private :: DGDerivative_MappedScalar1D
    generic,public :: BRDerivative => BRDerivative_MappedScalar1D
    procedure,private :: BRDerivative_MappedScalar1D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar1D

  endtype MappedScalar1D

contains

! ---------------------- Scalars ---------------------- !

  subroutine SetInteriorFromEquation_MappedScalar1D(this,geometry,time)
    !!  Sets the this % interior attribute using the eqn attribute,
    !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar1D),intent(inout) :: this
    type(Geometry1D),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,iEl,iVar

    do ivar = 1,this%nvar
      this%interior(:,:,ivar) = this%eqn(ivar)%evaluate(geometry%x%interior)
    enddo

  endsubroutine SetInteriorFromEquation_MappedScalar1D

  subroutine SideExchange_MappedScalar1D(this,mesh,decomp)
    implicit none
    class(MappedScalar1D),intent(inout) :: this
    type(Mesh1D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,ivar
    integer :: neighborRank
    integer :: rankId,offset

    !$omp target map(to:decomp % offsetElem, this % boundary) map(tofrom:this % extBoundary)
    !$omp teams distribute parallel do collapse(2) num_threads(256)
    do ivar = 1,this%nvar
      do e1 = 1,mesh%nElem

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
    enddo
    !$omp end target

  endsubroutine SideExchange_MappedScalar1D

  subroutine BassiRebaySides_MappedScalar1D(this)
    implicit none
    class(MappedScalar1D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i

    !$omp target map(to:this % boundary, this % extBoundary) map(from:this % avgBoundary)
    !$omp teams distribute parallel do collapse(2)
    do iel = 1,this%nElem
      do ivar = 1,this%nVar

        ! Left side - we account for the -\hat{x} normal
        this%avgBoundary(1,iel,ivar) = -0.5_prec*( &
                                       this%boundary(1,iel,ivar)+ &
                                       this%extBoundary(1,iel,ivar))

        ! Right side - we account for the +\hat{x} normal
        this%avgBoundary(2,iel,ivar) = 0.5_prec*( &
                                       this%boundary(2,iel,ivar)+ &
                                       this%extBoundary(2,iel,ivar))
      enddo
    enddo
    !$omp end do

  endsubroutine BassiRebaySides_MappedScalar1D

  subroutine Derivative_MappedScalar1D(this,geometry,dF)
    implicit none
    class(MappedScalar1D),intent(in) :: this
    type(Geometry1D),intent(in) :: geometry
    real(prec),intent(out) :: df(1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i

    call this%interp%Derivative_1D(this%interior, &
                                   df, &
                                   this%nVar, &
                                   this%nElem)
    !$omp target map(to:geometry % dxds % interior) map(tofrom:df)
    !$omp teams distribute parallel do collapse(3)
    do iEl = 1,this%nElem
      do iVar = 1,this%nVar
        do i = 1,this%interp%N+1
          df(i,iEl,iVar) = df(i,iEl,iVar)/ &
                           geometry%dxds%interior(i,iEl,1)
        enddo
      enddo
    enddo
    !$omp end target
  endsubroutine Derivative_MappedScalar1D

  subroutine DGDerivative_MappedScalar1D(this,geometry,dF)
    implicit none
    class(MappedScalar1D),intent(in) :: this
    type(Geometry1D),intent(in) :: geometry
    real(prec),intent(out) :: df(1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i

    call this%interp%DGDerivative_1D(this%interior, &
                                     this%boundary, &
                                     df, &
                                     this%nVar, &
                                     this%nElem)

    !$omp target map(to:geometry % dxds % interior) map(tofrom:df)
    !$omp teams distribute parallel do collapse(3)
    do iEl = 1,this%nElem
      do iVar = 1,this%nVar
        do i = 1,this%interp%N+1
          df(i,iEl,iVar) = df(i,iEl,iVar)/ &
                           geometry%dxds%interior(i,iEl,1)
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine DGDerivative_MappedScalar1D

  subroutine BRDerivative_MappedScalar1D(this,geometry,dF)
    implicit none
    class(MappedScalar1D),intent(in) :: this
    type(Geometry1D),intent(in) :: geometry
    real(prec),intent(out) :: df(1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i

    call this%interp%DGDerivative_1D(this%interior, &
                                     this%avgboundary, &
                                     df, &
                                     this%nVar, &
                                     this%nElem)

    !$omp target map(to:geometry % dxds % interior) map(tofrom:df)
    !$omp teams distribute parallel do collapse(3)
    do iEl = 1,this%nElem
      do iVar = 1,this%nVar
        do i = 1,this%interp%N+1
          df(i,iEl,iVar) = df(i,iEl,iVar)/ &
                           geometry%dxds%interior(i,iEl,1)
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine BRDerivative_MappedScalar1D

endmodule SELF_MappedScalar_1D
