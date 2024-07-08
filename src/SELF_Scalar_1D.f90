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

module SELF_Scalar_1D

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
  type,extends(SELF_DataObj),public :: Scalar1D

    real(prec),allocatable,dimension(:,:,:) :: interior
    real(prec),allocatable,dimension(:,:,:) :: boundary
    real(prec),allocatable,dimension(:,:,:) :: extBoundary
    real(prec),allocatable,dimension(:,:,:) :: avgBoundary
    real(prec),allocatable,dimension(:,:,:) :: jumpBoundary

  contains

    procedure,public :: Init => Init_Scalar1D
    procedure,public :: Free => Free_Scalar1D

    procedure,public :: BoundaryInterp => BoundaryInterp_Scalar1D
    procedure,public :: GridInterp => GridInterp_Scalar1D
    generic,public :: Derivative => Derivative_Scalar1D
    procedure,private :: Derivative_Scalar1D

    generic,public :: WriteHDF5 => WriteHDF5_Scalar1D,WriteHDF5_MPI_Scalar1D
    procedure,private :: WriteHDF5_Scalar1D
    procedure,private :: WriteHDF5_MPI_Scalar1D

  endtype Scalar1D

contains

! -- Scalar1D -- !

  subroutine Init_Scalar1D(this,interp,nVar,nElem)
    implicit none
    class(Scalar1D),intent(out) :: this
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
             this%extBoundary(1:2,1:nelem,1:nvar), &
             this%avgBoundary(2,1:nelem,1:nvar), &
             this%jumpBoundary(1:2,1:nelem,1:nvar))

    !$omp target enter data map(alloc: this)
    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)
    !$omp target enter data map(alloc: this % avgBoundary)
    !$omp target enter data map(alloc: this % jumpBoundary)
    !$omp target enter data map(to: this % interp)

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:nVar))

  endsubroutine Init_Scalar1D

  subroutine Free_Scalar1D(this)
    implicit none
    class(Scalar1D),intent(inout) :: this

    this%interp => null()
    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%extBoundary)
    deallocate(this%avgBoundary)
    deallocate(this%jumpBoundary)
    deallocate(this%meta)
    deallocate(this%eqn)

    !$omp target exit data map(delete: this % interior)
    !$omp target exit data map(delete: this % boundary)
    !$omp target exit data map(delete: this % extBoundary)
    !$omp target exit data map(delete: this % avgBoundary)
    !$omp target exit data map(delete: this % jumpBoundary)
    !$omp target exit data map(release: this % interp)
    !$omp target exit data map(delete: this)

  endsubroutine Free_Scalar1D

  subroutine BoundaryInterp_Scalar1D(this)
    implicit none
    class(Scalar1D),intent(inout) :: this
    ! Local
    integer :: ii,iel,ivar
    real(prec) :: fb(1:2)

    !$omp target map(to:this%interior,this%interp%bMatrix) map(from:this%boundary)
    !$omp teams loop collapse(2)
    do iel = 1,this%nelem
      do ivar = 1,this%nvar
        fb(1:2) = 0.0_prec
        do ii = 1,this%N+1
          fb(1) = fb(1)+this%interp%bMatrix(ii,1)*this%interior(ii,iel,ivar) ! West
          fb(2) = fb(2)+this%interp%bMatrix(ii,2)*this%interior(ii,iel,ivar) ! East
        enddo
        this%boundary(1:2,iel,ivar) = fb(1:2)
      enddo
    enddo
    !$omp end target
    !call self_hipblas_matrixop_1d(this % bMatrix,f,fTarget,2,this % N + 1,nvars*nelems,handle)

  endsubroutine BoundaryInterp_Scalar1D

  function GridInterp_Scalar1D(this) result(f)
    implicit none
    class(Scalar1D),intent(in) :: this
    real(prec) :: f(1:this%M+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iel,ivar,i,ii
    real(prec) :: floc

    !$omp target map(to:this % interior, this % interp % iMatrix) map(from:f)
    !$omp teams loop bind(teams) collapse(3)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do i = 1,this%M+1
          floc = 0.0_prec
          !$omp loop bind(thread)
          do ii = 1,this%N+1
            floc = floc+this%interp%iMatrix(ii,i)*this%interior(ii,iel,ivar)
          enddo
          f(i,iel,ivar) = floc
        enddo
      enddo
    enddo
    !$omp end target

  endfunction GridInterp_Scalar1D

  function Derivative_Scalar1D(this) result(df)
    implicit none
    class(Scalar1D),intent(in) :: this
    real(prec) :: df(1:this%N+1,1:this%nelem,1:this%nvar)

    ! Local
    integer :: i,ii,iel,ivar
    real(prec) :: dfloc

    !$omp target map(to:this%interior,this%interp%dMatrix) map(from:df)
    !$omp teams loop collapse(3)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do i = 1,this%N+1

          dfloc = 0.0_prec
          do ii = 1,this%N+1
            dfloc = dfloc+this%interp%dMatrix(ii,i)*this%interior(ii,iel,ivar)
          enddo
          df(i,iel,ivar) = dfloc

        enddo
      enddo
    enddo
    !$omp end target

    !call self_hipblas_matrixop_1d(this % dMatrix,f,df,this % N + 1,this % N + 1,nvars*nelems,handle)

  endfunction Derivative_Scalar1D

  subroutine WriteHDF5_MPI_Scalar1D(this,fileId,group,elemoffset,nglobalelem)
    implicit none
    class(Scalar1D),intent(in) :: this
    character(*),intent(in) :: group
    integer(HID_T),intent(in) :: fileId
    integer,intent(in) :: elemoffset
    integer,intent(in) :: nglobalelem
    ! Local
    integer(HID_T) :: offset(1:3)
    integer(HID_T) :: bOffset(1:3)
    integer(HID_T) :: globalDims(1:3)
    integer(HID_T) :: bGlobalDims(1:3)
    integer :: ivar

    offset(1:3) = (/0,0,elemoffset/)
    globalDims(1:3) = (/this%interp%N+1, &
                        this%nVar, &
                        nGlobalElem/)

    ! Offsets and dimensions for element boundary data
    bOffset(1:3) = (/0,0,elemoffset/)
    bGlobalDims(1:3) = (/this%nVar, &
                         2, &
                         nGlobalElem/)

    call CreateGroup_HDF5(fileId,trim(group))

    do ivar = 1,this%nVar
      call this%meta(ivar)%WriteHDF5(group,ivar,fileId)
    enddo

    call WriteArray_HDF5(fileId,trim(group)//"/interior", &
                         this%interior,offset,globalDims)

    call WriteArray_HDF5(fileId,trim(group)//"/boundary", &
                         this%boundary,bOffset,bGlobalDims)

  endsubroutine WriteHDF5_MPI_Scalar1D

  subroutine WriteHDF5_Scalar1D(this,fileId,group)
    implicit none
    class(Scalar1D),intent(in) :: this
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

  endsubroutine WriteHDF5_Scalar1D

endmodule SELF_Scalar_1D
