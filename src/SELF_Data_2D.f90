! SELF_Data.F90
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Data_2D

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

  type,extends(SELF_DataObj),public :: Tensor2D

    real(prec),pointer,dimension(:,:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:,:) :: extBoundary

  contains

    procedure,public :: Init => Init_Tensor2D
    procedure,public :: Free => Free_Tensor2D

    procedure,public :: BoundaryInterp => BoundaryInterp_Tensor2D
    procedure,public :: Divergence => Divergence_Tensor2D
    procedure,public :: DGDivergence => DGDivergence_Tensor2D
    procedure,public :: Determinant => Determinant_Tensor2D

  endtype Tensor2D

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
    type(Vector2D),intent(inout) :: df

    call this%interp%ScalarGradient_2D(this%interior, &
                                       df%interior, &
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

  subroutine Init_Vector2D(this,interp,nVar,nElem)
    implicit none
    class(Vector2D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer :: N

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    N = interp%N

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
    type(Tensor2D),intent(inout) :: df

    call this%interp%VectorGradient_2D(this%interior, &
                                       df%interior, &
                                       this%nVar, &
                                       this%nElem)

  endsubroutine Gradient_Vector2D

  subroutine Divergence_Vector2D(this,that)
    implicit none
    class(Vector2D),intent(in) :: this
    type(Scalar2D),intent(inout) :: that

    call this%interp%VectorDivergence_2D(this%interior, &
                                         that%interior, &
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

  subroutine Init_Tensor2D(this,interp,nVar,nElem)
    implicit none
    class(Tensor2D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer :: N

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    N = interp%N

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:nelem,1:nvar,1:2,1:2), &
             this%boundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2,1:2), &
             this%extBoundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2,1:2))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:4*nVar))

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)

  endsubroutine Init_Tensor2D

  subroutine Free_Tensor2D(this)
    implicit none
    class(Tensor2D),intent(inout) :: this

    this%interp => null()
    this%nVar = 0
    this%nElem = 0

    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%extBoundary)

    !$omp target exit data map(delete: this % interior)
    !$omp target exit data map(delete: this % boundary)
    !$omp target exit data map(delete: this % extBoundary)

    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_Tensor2D

  subroutine BoundaryInterp_Tensor2D(this)
    implicit none
    class(Tensor2D),intent(inout) :: this

    call this%interp%TensorBoundaryInterp_2D(this%interior, &
                                             this%boundary, &
                                             this%nVar, &
                                             this%nElem)

  endsubroutine BoundaryInterp_Tensor2D

  subroutine Divergence_Tensor2D(this,that)
    implicit none
    class(Tensor2D),intent(in) :: this
    class(Vector2D),intent(inout) :: that

    call this%interp%TensorDivergence_2D(this%interior, &
                                         that%interior, &
                                         this%nVar, &
                                         this%nElem)

  endsubroutine Divergence_Tensor2D

  subroutine DGDivergence_Tensor2D(this,that)
    implicit none
    class(Tensor2D),intent(in) :: this
    class(Vector2D),intent(inout) :: that

    call this%interp%TensorDGDivergence_2D(this%interior, &
                                           this%boundary, &
                                           that%interior, &
                                           this%nVar, &
                                           this%nElem)

  endsubroutine DGDivergence_Tensor2D

  subroutine Determinant_Tensor2D(this,that)
    implicit none
    class(Tensor2D),intent(in) :: this
    type(Scalar2D),intent(inout) :: that
    ! Local
    integer :: iEl,iVar,i,j

    do iVar = 1,this%nVar
      do iEl = 1,this%nElem
        do j = 1,this%interp%N+1
          do i = 1,this%interp%N+1

            that%interior(i,j,iEl,iVar) = this%interior(i,j,iEl,iVar,1,1)* &
                                          this%interior(i,j,iEl,iVar,2,2)- &
                                          this%interior(i,j,iEl,iVar,1,2)* &
                                          this%interior(i,j,iEl,iVar,2,1)

          enddo
        enddo
      enddo
    enddo

  endsubroutine Determinant_Tensor2D

endmodule SELF_Data_2D
