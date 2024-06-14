! SELF_Data.F90
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Vector_2D

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

  type,extends(SELF_DataObj),public :: Vector2D

    real(prec),pointer,dimension(:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:) :: extBoundary
    real(prec),pointer,dimension(:,:,:,:,:) :: avgBoundary
    real(prec),pointer,dimension(:,:,:,:) :: boundaryNormal

  contains

    procedure,public :: Init => Init_Vector2D
    procedure,public :: Free => Free_Vector2D

    procedure,public :: BoundaryInterp => BoundaryInterp_Vector2D
    procedure,public :: GridInterp => GridInterp_Vector2D
    procedure,public :: Gradient => Gradient_Vector2D
    generic,public :: Divergence => Divergence_Vector2D
    procedure,private :: Divergence_Vector2D

    generic,public :: SetEquation => SetEquation_Vector2D
    procedure,private :: SetEquation_Vector2D

    generic,public :: WriteHDF5 => WriteHDF5_MPI_Vector2D,WriteHDF5_Vector2D
    procedure,private :: WriteHDF5_MPI_Vector2D
    procedure,private :: WriteHDF5_Vector2D

  endtype Vector2D

contains

  subroutine Init_Vector2D(this,interp,nVar,nElem)
    implicit none
    class(Vector2D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:nelem,1:nvar,1:2), &
             this%boundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2), &
             this%extBoundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2), &
             this%avgBoundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2), &
             this%boundaryNormal(1:interp%N+1,1:4,1:nelem,1:nvar))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:2*nVar))

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)
    !$omp target enter data map(alloc: this % avgBoundary)
    !$omp target enter data map(alloc: this % boundaryNormal)

  endsubroutine Init_Vector2D

  subroutine Free_Vector2D(this)
    implicit none
    class(Vector2D),intent(inout) :: this

    this%interp => null()
    this%nVar = 0
    this%nElem = 0

    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%boundaryNormal)
    deallocate(this%extBoundary)
    deallocate(this%avgBoundary)

    !$omp target exit data map(delete: this % interior)
    !$omp target exit data map(delete: this % boundary)
    !$omp target exit data map(delete: this % boundaryNormal)
    !$omp target exit data map(delete: this % extBoundary)
    !$omp target exit data map(delete: this % avgBoundary)

    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_Vector2D

  subroutine SetEquation_Vector2D(this,idir,ivar,eqnChar)
    !! Sets the equation parser for the `idir` direction and `ivar-th` variable
    implicit none
    class(Vector2D),intent(inout) :: this
    integer,intent(in) :: idir,ivar
    character(*),intent(in) :: eqnChar

    this%eqn(idir+2*(ivar-1)) = EquationParser(trim(eqnChar), &
                                               (/'x','y','z','t'/))

  endsubroutine SetEquation_Vector2D

  subroutine GridInterp_Vector2D(this,that)
    implicit none
    class(Vector2D),intent(in) :: this
    type(Vector2D),intent(inout) :: that

    call this%interp%VectorGridInterp_2D(this%interior, &
                                         that%interior, &
                                         this%nVar, &
                                         this%nElem)

  endsubroutine GridInterp_Vector2D

  subroutine BoundaryInterp_Vector2D(this)
    implicit none
    class(Vector2D),intent(inout) :: this

    call this%interp%VectorBoundaryInterp_2D(this%interior, &
                                             this%boundary, &
                                             this%nVar, &
                                             this%nElem)

  endsubroutine BoundaryInterp_Vector2D

  subroutine Gradient_Vector2D(this,df)
    implicit none
    class(Vector2D),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:2,1:2)

    call this%interp%VectorGradient_2D(this%interior, &
                                       df, &
                                       this%nVar, &
                                       this%nElem)

  endsubroutine Gradient_Vector2D

  subroutine Divergence_Vector2D(this,df)
    implicit none
    class(Vector2D),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)

    call this%interp%VectorDivergence_2D(this%interior, &
                                         df, &
                                         this%nVar, &
                                         this%nElem)

  endsubroutine Divergence_Vector2D

  subroutine WriteHDF5_MPI_Vector2D(this,fileId,group,elemoffset,nglobalelem)
    implicit none
    class(Vector2D),intent(in) :: this
    character(*),intent(in) :: group
    integer(HID_T),intent(in) :: fileId
    integer,intent(in) :: elemoffset
    integer,intent(in) :: nglobalelem
    ! Local
    integer(HID_T) :: offset(1:5)
    integer(HID_T) :: bOffset(1:5)
    integer(HID_T) :: globalDims(1:5)
    integer(HID_T) :: bGlobalDims(1:5)
    integer :: ivar

    offset(1:5) = (/0,0,0,0,elemoffset/)
    globalDims(1:5) = (/2, &
                        this%interp%N+1, &
                        this%interp%N+1, &
                        this%nVar, &
                        nglobalelem/)

    ! Offsets and dimensions for element boundary data
    bOffset(1:5) = (/0,0,0,0,elemoffset/)
    bGlobalDims(1:5) = (/2, &
                         this%interp%N+1, &
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

  endsubroutine WriteHDF5_MPI_Vector2D

  subroutine WriteHDF5_Vector2D(this,fileId,group)
    implicit none
    class(Vector2D),intent(in) :: this
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

  endsubroutine WriteHDF5_Vector2D

endmodule SELF_Vector_2D
