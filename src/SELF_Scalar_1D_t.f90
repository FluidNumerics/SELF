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

module SELF_Scalar_1D_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Metadata
  use FEQParse
  use SELF_Data
  use SELF_HDF5

  use HDF5
  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

! ---------------------- Scalars ---------------------- !
  type,extends(SELF_DataObj),public :: Scalar1D_t

    real(prec),pointer,contiguous,dimension(:,:,:) :: interior
    real(prec),pointer,contiguous,dimension(:,:,:) :: boundary
    real(prec),pointer,contiguous,dimension(:,:,:) :: boundarynormal
    real(prec),pointer,contiguous,dimension(:,:,:) :: extBoundary
    real(prec),pointer,contiguous,dimension(:,:,:) :: avgBoundary

  contains

    procedure,public :: Init => Init_Scalar1D_t
    procedure,public :: Free => Free_Scalar1D_t

    procedure,public :: UpdateHost => UpdateHost_Scalar1D_t
    procedure,public :: UpdateDevice => UpdateDevice_Scalar1D_t

    procedure,public :: AverageSides => AverageSides_Scalar1D_t
    procedure,public :: BoundaryInterp => BoundaryInterp_Scalar1D_t
    generic,public :: GridInterp => GridInterp_Scalar1D_t
    procedure,private :: GridInterp_Scalar1D_t
    generic,public :: Derivative => Derivative_Scalar1D_t
    procedure,private :: Derivative_Scalar1D_t
    generic,public :: WriteHDF5 => WriteHDF5_Scalar1D_t
    procedure,private :: WriteHDF5_Scalar1D_t

  endtype Scalar1D_t

contains

! -- Scalar1D_t -- !

  subroutine Init_Scalar1D_t(this,interp,nVar,nElem)
    implicit none
    class(Scalar1D_t),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,1:nelem,1:nvar), &
             this%boundary(1:2,1:nelem,1:nvar), &
             this%boundarynormal(1:2,1:nelem,1:nvar), &
             this%extBoundary(1:2,1:nelem,1:nvar), &
             this%avgBoundary(1:2,1:nelem,1:nvar))

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%boundarynormal = 0.0_prec
    this%extBoundary = 0.0_prec
    this%avgBoundary = 0.0_prec

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:nVar))

  endsubroutine Init_Scalar1D_t

  subroutine Free_Scalar1D_t(this)
    implicit none
    class(Scalar1D_t),intent(inout) :: this

    this%interp => null()
    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%boundarynormal)
    deallocate(this%extBoundary)
    deallocate(this%avgBoundary)
    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_Scalar1D_t

  subroutine UpdateHost_Scalar1D_t(this)
    implicit none
    class(Scalar1D_t),intent(inout) :: this

  endsubroutine UpdateHost_Scalar1D_t

  subroutine UpdateDevice_Scalar1D_t(this)
    implicit none
    class(Scalar1D_t),intent(inout) :: this

  endsubroutine UpdateDevice_Scalar1D_t

  subroutine AverageSides_Scalar1D_t(this)
    implicit none
    class(Scalar1D_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: ivar

    do concurrent(iel=1:this%nElem,ivar=1:this%nVar)
      ! Left side - we account for the -\hat{x} normal
      this%avgboundary(1,iel,ivar) = 0.5_prec*( &
                                     this%boundary(1,iel,ivar)+ &
                                     this%extBoundary(1,iel,ivar))

      ! Right side - we account for the +\hat{x} normal
      this%avgboundary(2,iel,ivar) = 0.5_prec*( &
                                     this%boundary(2,iel,ivar)+ &
                                     this%extBoundary(2,iel,ivar))
    enddo

  endsubroutine AverageSides_Scalar1D_t

  subroutine BoundaryInterp_Scalar1D_t(this)
    implicit none
    class(Scalar1D_t),intent(inout) :: this
    ! Local
    integer :: ii,iel,ivar
    real(prec) :: fbl,fbr

    do concurrent(iel=1:this%nElem,ivar=1:this%nVar)
      fbl = 0.0_prec
      fbr = 0.0_prec
      do ii = 1,this%N+1
        fbl = fbl+this%interp%bMatrix(ii,1)*this%interior(ii,iel,ivar) ! West
        fbr = fbr+this%interp%bMatrix(ii,2)*this%interior(ii,iel,ivar) ! East
      enddo
      this%boundary(1,iel,ivar) = fbl
      this%boundary(2,iel,ivar) = fbr
    enddo

  endsubroutine BoundaryInterp_Scalar1D_t

  subroutine GridInterp_Scalar1D_t(this,f)
    implicit none
    class(Scalar1D_t),intent(in) :: this
    real(prec),intent(inout) :: f(1:this%M+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iel,ivar,i,ii
    real(prec) :: floc

    do concurrent(i=1:this%M+1,iel=1:this%nElem,ivar=1:this%nVar)
      floc = 0.0_prec
      do ii = 1,this%N+1
        floc = floc+this%interp%iMatrix(ii,i)*this%interior(ii,iel,ivar)
      enddo
      f(i,iel,ivar) = floc
    enddo

  endsubroutine GridInterp_Scalar1D_t

  subroutine Derivative_Scalar1D_t(this,df)
    implicit none
    class(Scalar1D_t),intent(in) :: this
    real(prec),intent(inout) :: df(1:this%N+1,1:this%nelem,1:this%nvar)

    ! Local
    integer :: i,ii,iel,ivar
    real(prec) :: dfloc

    do concurrent(i=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar)
      dfloc = 0.0_prec
      do ii = 1,this%N+1
        dfloc = dfloc+this%interp%dMatrix(ii,i)*this%interior(ii,iel,ivar)
      enddo
      df(i,iel,ivar) = dfloc
    enddo

  endsubroutine Derivative_Scalar1D_t

  subroutine WriteHDF5_Scalar1D_t(this,fileId,group)
    implicit none
    class(Scalar1D_t),intent(in) :: this
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: group
    ! Local
    integer :: ivar

    call CreateGroup_HDF5(fileId,trim(group))

    do ivar = 1,this%nVar
      call this%meta(ivar)%WriteHDF5(group,ivar,fileId)
    enddo

    call WriteArray_HDF5(fileId,trim(group)//"/interior", &
                         this%interior)

    call WriteArray_HDF5(fileId,trim(group)//"/boundary", &
                         this%boundary)

  endsubroutine WriteHDF5_Scalar1D_t

endmodule SELF_Scalar_1D_t
