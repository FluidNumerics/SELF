! SELF_Data.F90
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Data_1D

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

    real(prec),pointer,dimension(:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:) :: extBoundary
    real(prec),pointer,dimension(:,:,:) :: avgBoundary
    real(prec),pointer,dimension(:,:,:) :: jumpBoundary

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

    allocate(this%interior(1:interp%N+1,1:nelem,1:nvar), &
             this%boundary(1:2,1:nelem,1:nvar), &
             this%extBoundary(1:2,1:nelem,1:nvar), &
             this%avgBoundary(2,1:nelem,1:nvar), &
             this%jumpBoundary(1:2,1:nelem,1:nvar))

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)
    !$omp target enter data map(alloc: this % avgBoundary)
    !$omp target enter data map(alloc: this % jumpBoundary)

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

  endsubroutine Free_Scalar1D

  subroutine BoundaryInterp_Scalar1D(this)
    implicit none
    class(Scalar1D),intent(inout) :: this

    call this%interp%ScalarBoundaryInterp_1D(this%interior, &
                                             this%boundary, &
                                             this%nVar, &
                                             this%nElem)

  endsubroutine BoundaryInterp_Scalar1D

  subroutine GridInterp_Scalar1D(this,that)
    implicit none
    class(Scalar1D),intent(in) :: this
    type(Scalar1D),intent(inout) :: that

    call this%interp%ScalarGridInterp_1D(this%interior, &
                                         that%interior, &
                                         this%nVar, &
                                         this%nElem)

  endsubroutine GridInterp_Scalar1D

  subroutine Derivative_Scalar1D(this,that)
    implicit none
    class(Scalar1D),intent(in) :: this
    type(Scalar1D),intent(inout) :: that

    call this%interp%Derivative_1D(this%interior, &
                                   that%interior, &
                                   this%nVar, &
                                   this%nElem)

  endsubroutine Derivative_Scalar1D

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

endmodule SELF_Data_1D
