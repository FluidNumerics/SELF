! SELF_Data.F90
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Scalar_2D

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Metadata
  use FEQParse
  use SELF_HDF5
  use SELF_Data

  use HDF5
  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type,extends(SELF_DataObj),public :: Scalar2D

    real(prec),pointer,dimension(:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:) :: interpWork
    real(prec),pointer,dimension(:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:) :: extBoundary
    real(prec),pointer,dimension(:,:,:,:) :: avgBoundary
    real(prec),pointer,dimension(:,:,:,:) :: jumpBoundary

  contains

    procedure,public :: Init => Init_Scalar2D
    procedure,public :: Free => Free_Scalar2D

    procedure,public :: BoundaryInterp => BoundaryInterp_Scalar2D
    procedure,public :: GridInterp => GridInterp_Scalar2D
    generic,public :: Gradient => Gradient_Scalar2D
    procedure,private :: Gradient_Scalar2D

    generic,public :: WriteHDF5 => WriteHDF5_MPI_Scalar2D,WriteHDF5_Scalar2D
    procedure,private :: WriteHDF5_MPI_Scalar2D
    procedure,private :: WriteHDF5_Scalar2D

  endtype Scalar2D

contains

  subroutine Init_Scalar2D(this,interp,nVar,nElem)
    implicit none
    class(Scalar2D),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,interp%N+1,nelem,nvar), &
             this%interpWork(1:interp%M+1,1:interp%N+1,1:nelem,1:nvar), &
             this%boundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
             this%extBoundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
             this%avgBoundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
             this%jumpBoundary(1:interp%N+1,1:4,1:nelem,1:nvar))

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % interpWork)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)
    !$omp target enter data map(alloc: this % avgBoundary)
    !$omp target enter data map(alloc: this % jumpBoundary)

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:nVar))

  endsubroutine Init_Scalar2D

  subroutine Free_Scalar2D(this)
    implicit none
    class(Scalar2D),intent(inout) :: this

    this%nVar = 0
    this%nElem = 0
    this%interp => null()
    deallocate(this%interior)
    deallocate(this%interpWork)
    deallocate(this%boundary)
    deallocate(this%extBoundary)
    deallocate(this%avgBoundary)
    deallocate(this%jumpBoundary)
    deallocate(this%meta)
    deallocate(this%eqn)

    !$omp target exit data map(delete: this % interior)
    !$omp target exit data map(delete: this % interpWork)
    !$omp target exit data map(delete: this % boundary)
    !$omp target exit data map(delete: this % extBoundary)
    !$omp target exit data map(delete: this % avgBoundary)
    !$omp target exit data map(delete: this % jumpBoundary)

  endsubroutine Free_Scalar2D

  subroutine BoundaryInterp_Scalar2D(this)
    implicit none
    class(Scalar2D),intent(inout) :: this

    call this%interp%ScalarBoundaryInterp_2D(this%interior, &
                                             this%boundary, &
                                             this%nVar, &
                                             this%nElem)

  endsubroutine BoundaryInterp_Scalar2D

  subroutine GridInterp_Scalar2D(this,that)
    implicit none
    class(Scalar2D),intent(in) :: this
    type(Scalar2D),intent(inout) :: that

    call this%interp%ScalarGridInterp_2D(this%interior, &
                                         that%interior, &
                                         this%nVar, &
                                         this%nElem)

  endsubroutine GridInterp_Scalar2D

  subroutine Gradient_Scalar2D(this,df)
    implicit none
    class(Scalar2D),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:2)

    call this%interp%ScalarGradient_2D(this%interior, &
                                       df, &
                                       this%nVar, &
                                       this%nElem)

  endsubroutine Gradient_Scalar2D

  subroutine WriteHDF5_MPI_Scalar2D(this,fileId,group,elemoffset,nglobalelem)
    implicit none
    class(Scalar2D),intent(in) :: this
    character(*),intent(in) :: group
    integer(HID_T),intent(in) :: fileId
    integer,intent(in) :: elemoffset
    integer,intent(in) :: nglobalelem
    ! Local
    integer(HID_T) :: offset(1:4)
    integer(HID_T) :: bOffset(1:4)
    integer(HID_T) :: globalDims(1:4)
    integer(HID_T) :: bGlobalDims(1:4)
    integer :: ivar

    offset(1:4) = (/0,0,0,elemoffset/)
    globalDims(1:4) = (/this%interp%N+1, &
                        this%interp%N+1, &
                        this%nVar, &
                        nglobalelem/)

    ! Offsets and dimensions for element boundary data
    bOffset(1:4) = (/0,0,0,elemoffset/)
    bGlobalDims(1:4) = (/this%interp%N+1, &
                         this%nVar, &
                         4, &
                         nglobalelem/)

    call CreateGroup_HDF5(fileId,trim(group))

    do ivar = 1,this%nVar
      call this%meta(ivar)%WriteHDF5(group,ivar,fileId)
    enddo

    call WriteArray_HDF5(fileId,trim(group)//"/interior", &
                         this%interior,offset,globalDims)

    call WriteArray_HDF5(fileId,trim(group)//"/boundary", &
                         this%boundary,bOffset,bGlobalDims)

  endsubroutine WriteHDF5_MPI_Scalar2D

  subroutine WriteHDF5_Scalar2D(this,fileId,group)
    implicit none
    class(Scalar2D),intent(in) :: this
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

  endsubroutine WriteHDF5_Scalar2D

endmodule SELF_Scalar_2D
