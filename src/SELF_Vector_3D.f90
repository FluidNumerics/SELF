! SELF_Data.F90
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Vector_3D

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

  type,extends(SELF_DataObj),public :: Vector3D

    real(prec),pointer,dimension(:,:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:,:) :: extBoundary
    real(prec),pointer,dimension(:,:,:,:,:,:) :: avgBoundary
    real(prec),pointer,dimension(:,:,:,:,:) :: boundaryNormal

  contains

    procedure,public :: Init => Init_Vector3D
    procedure,public :: Free => Free_Vector3D

    procedure,public :: BoundaryInterp => BoundaryInterp_Vector3D
    procedure,public :: GridInterp => GridInterp_Vector3D
    procedure,public :: Gradient => Gradient_Vector3D
    generic,public :: Divergence => Divergence_Vector3D
    procedure,private :: Divergence_Vector3D

    generic,public :: SetEquation => SetEquation_Vector3D
    procedure,private :: SetEquation_Vector3D

    generic,public :: WriteHDF5 => WriteHDF5_MPI_Vector3D,WriteHDF5_Vector3D
    procedure,private :: WriteHDF5_MPI_Vector3D
    procedure,private :: WriteHDF5_Vector3D

  endtype Vector3D

contains

  subroutine Init_Vector3D(this,interp,nVar,nElem)
    implicit none
    class(Vector3D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar,1:3), &
             this%boundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar,1:3), &
             this%extBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar,1:3), &
             this%avgBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar,1:3), &
             this%boundaryNormal(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:3*nVar))

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)
    !$omp target enter data map(alloc: this % avgBoundary)
    !$omp target enter data map(alloc: this % boundaryNormal)

  endsubroutine Init_Vector3D

  subroutine Free_Vector3D(this)
    implicit none
    class(Vector3D),intent(inout) :: this

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
    !$omp target exit data map(delete: this % extBoundary)
    !$omp target exit data map(delete: this % avgBoundary)
    !$omp target exit data map(delete: this % boundaryNormal)

    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_Vector3D

  subroutine SetEquation_Vector3D(this,idir,ivar,eqnChar)
    !! Sets the equation parser for the `idir` direction and `ivar-th` variable
    implicit none
    class(Vector3D),intent(inout) :: this
    integer,intent(in) :: idir,ivar
    character(*),intent(in) :: eqnChar

    this%eqn(idir+3*(ivar-1)) = EquationParser(trim(eqnChar), &
                                               (/'x','y','z','t'/))

  endsubroutine SetEquation_Vector3D

  subroutine GridInterp_Vector3D(this,that)
    implicit none
    class(Vector3D),intent(in) :: this
    type(Vector3D),intent(inout) :: that

    call this%interp%VectorGridInterp_3D(this%interior, &
                                         that%interior, &
                                         this%nVar, &
                                         this%nElem)

  endsubroutine GridInterp_Vector3D

  subroutine BoundaryInterp_Vector3D(this)
    implicit none
    class(Vector3D),intent(inout) :: this

    call this%interp%VectorBoundaryInterp_3D(this%interior, &
                                             this%boundary, &
                                             this%nVar, &
                                             this%nElem)

  endsubroutine BoundaryInterp_Vector3D

  subroutine Gradient_Vector3D(this,df)
    implicit none
    class(Vector3D),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:3,1:3)

    call this%interp%VectorGradient_3D(this%interior, &
                                       df, &
                                       this%nVar, &
                                       this%nElem)

  endsubroutine Gradient_Vector3D

  subroutine Divergence_Vector3D(this,df)
    implicit none
    class(Vector3D),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)

    call this%interp%VectorDivergence_3D(this%interior, &
                                         df, &
                                         this%nVar, &
                                         this%nElem)

  endsubroutine Divergence_Vector3D

  subroutine WriteHDF5_MPI_Vector3D(this,fileId,group,elemoffset,nglobalelem)
    implicit none
    class(Vector3D),intent(in) :: this
    character(*),intent(in) :: group
    integer(HID_T),intent(in) :: fileId
    integer,intent(in) :: elemoffset
    integer,intent(in) :: nglobalelem
    ! Local
    integer(HID_T) :: offset(1:6)
    integer(HID_T) :: bOffset(1:6)
    integer(HID_T) :: globalDims(1:6)
    integer(HID_T) :: bGlobalDims(1:6)
    integer :: ivar

    offset(1:6) = (/0,0,0,0,0,elemoffset/)
    globalDims(1:6) = (/3, &
                        this%interp%N+1, &
                        this%interp%N+1, &
                        this%interp%N+1, &
                        this%nVar, &
                        nglobalelem/)

    ! Offsets and dimensions for element boundary data
    bOffset(1:6) = (/0,0,0,0,0,elemoffset/)
    bGlobalDims(1:6) = (/3, &
                         this%interp%N+1, &
                         this%interp%N+1, &
                         this%nVar, &
                         6, &
                         nglobalelem/)

    call CreateGroup_HDF5(fileId,trim(group))

    do ivar = 1,this%nVar
      call this%meta(ivar)%WriteHDF5(group,ivar,fileId)
    enddo

    call WriteArray_HDF5(fileId,trim(group)//"/interior", &
                         this%interior,offset,globalDims)

    call WriteArray_HDF5(fileId,trim(group)//"/boundary", &
                         this%boundary,bOffset,bGlobalDims)

  endsubroutine WriteHDF5_MPI_Vector3D

  subroutine WriteHDF5_Vector3D(this,fileId,group)
    implicit none
    class(Vector3D),intent(in) :: this
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

  endsubroutine WriteHDF5_Vector3D

endmodule SELF_Vector_3D
