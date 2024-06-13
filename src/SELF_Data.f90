! SELF_Data.F90
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Data

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Metadata
  use FEQParse
  use SELF_HDF5

  use HDF5
  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type,public :: SELF_DataObj
  !! The SELF_DataObj class is a base class for all data objects in SELF.
  !! A data object in SELF is a multidimensional array of data, represented
  !! on both host and device, that is associated with an interpolant, metadata,
  !! and (optionally) an equation string.
  !! Type extensions of the SELF_DataObj include scalars, vectors, and tensors
  !! in 1-D, 2-D, and 3-D using the storage patterns that are expected for
  !! derivative and interpolation operations defined in SELF_Lagrange.f90
  !! Additionally, each extended type has the necessary attributes to store
  !! information on element interiors and element boundaries, both of which
  !! are commonly used for spectral element solvers.

    integer :: nVar
    integer :: nElem
    type(Lagrange),pointer :: interp
    type(Metadata),allocatable :: meta(:)
    type(EquationParser),allocatable :: eqn(:)

  contains

    ! Procedures for setting metadata for
    procedure,public :: SetName => SetName_DataObj
    procedure,public :: SetDescription => SetDescription_DataObj
    procedure,public :: SetUnits => SetUnits_DataObj
    generic,public :: SetEquation => SetEquation_DataObj
    procedure,private :: SetEquation_DataObj

  endtype SELF_DataObj

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

  type,extends(SELF_DataObj),public :: Scalar3D

    real(prec),pointer,dimension(:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:) :: extBoundary
    real(prec),pointer,dimension(:,:,:,:,:) :: avgBoundary
    real(prec),pointer,dimension(:,:,:,:,:) :: jumpBoundary

    real(prec),pointer,dimension(:,:,:,:,:) :: interpWork1
    real(prec),pointer,dimension(:,:,:,:,:) :: interpWork2

  contains

    procedure,public :: Init => Init_Scalar3D
    procedure,public :: Free => Free_Scalar3D

    procedure,public :: BoundaryInterp => BoundaryInterp_Scalar3D
    procedure,public :: GridInterp => GridInterp_Scalar3D
    generic,public :: Gradient => Gradient_Scalar3D
    procedure,private :: Gradient_Scalar3D

    generic,public :: WriteHDF5 => WriteHDF5_MPI_Scalar3D,WriteHDF5_Scalar3D
    procedure,private :: WriteHDF5_MPI_Scalar3D
    procedure,private :: WriteHDF5_Scalar3D

  endtype Scalar3D

! ! ! ---------------------- Vectors ---------------------- !

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
! ! ---------------------- Tensors ---------------------- !

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

  type,extends(SELF_DataObj),public :: Tensor3D

    real(prec),pointer,dimension(:,:,:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:,:,:) :: extBoundary

  contains

    procedure,public :: Init => Init_Tensor3D
    procedure,public :: Free => Free_Tensor3D

    procedure,public :: BoundaryInterp => BoundaryInterp_Tensor3D
    procedure,public :: Divergence => Divergence_Tensor3D
    procedure,public :: DGDivergence => DGDivergence_Tensor3D
    procedure,public :: Determinant => Determinant_Tensor3D

  endtype Tensor3D

  integer,parameter :: selfStrongForm = 0
  integer,parameter :: selfWeakDGForm = 1
  integer,parameter :: selfWeakCGForm = 2
  integer,parameter :: selfWeakBRForm = 3

contains

! -- DataObj -- !

  subroutine SetName_DataObj(this,ivar,name)
    !! Set the name of the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: name

    call this%meta(ivar)%SetName(name)

  endsubroutine SetName_DataObj

  subroutine SetDescription_DataObj(this,ivar,description)
    !! Set the description of the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: description

    call this%meta(ivar)%SetDescription(description)

  endsubroutine SetDescription_DataObj

  subroutine SetUnits_DataObj(this,ivar,units)
    !! Set the units of the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: units

    call this%meta(ivar)%SetUnits(units)

  endsubroutine SetUnits_DataObj

  subroutine SetEquation_DataObj(this,ivar,eqnChar)
    !! Sets the equation parser for the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: eqnChar

    this%eqn(ivar) = EquationParser(trim(eqnChar), &
                                    (/'x','y','z','t'/))

  endsubroutine SetEquation_DataObj

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

! ! -- Scalar2D -- !

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

! ! -- Scalar3D -- !

  subroutine Init_Scalar3D(this,interp,nVar,nElem)
    implicit none
    class(Scalar3D),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar), &
             this%interpWork1(1:interp%M+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar), &
             this%interpWork2(1:interp%M+1,1:interp%M+1,1:interp%N+1,1:nelem,1:nvar), &
             this%boundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%extBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%avgBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%jumpBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar))

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % interpWork1)
    !$omp target enter data map(alloc: this % interpWork2)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)
    !$omp target enter data map(alloc: this % avgBoundary)
    !$omp target enter data map(alloc: this % jumpBoundary)

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:nVar))

  endsubroutine Init_Scalar3D

  subroutine Free_Scalar3D(this)
    implicit none
    class(Scalar3D),intent(inout) :: this

    this%nVar = 0
    this%nElem = 0
    this%interp => null()
    deallocate(this%interior)
    deallocate(this%interpWork1)
    deallocate(this%interpWork2)
    deallocate(this%boundary)
    deallocate(this%extBoundary)
    deallocate(this%avgBoundary)
    deallocate(this%jumpBoundary)
    deallocate(this%meta)
    deallocate(this%eqn)

    !$omp target exit data map(delete: this % interior)
    !$omp target exit data map(delete: this % interpWork1)
    !$omp target exit data map(delete: this % interpWork2)
    !$omp target exit data map(delete: this % boundary)
    !$omp target exit data map(delete: this % extBoundary)
    !$omp target exit data map(delete: this % avgBoundary)
    !$omp target exit data map(delete: this % jumpBoundary)

  endsubroutine Free_Scalar3D

  subroutine BoundaryInterp_Scalar3D(this)
    implicit none
    class(Scalar3D),intent(inout) :: this

    call this%interp%ScalarBoundaryInterp_3D(this%interior, &
                                             this%boundary, &
                                             this%nVar, &
                                             this%nElem)

  endsubroutine BoundaryInterp_Scalar3D

  subroutine GridInterp_Scalar3D(this,that)
    implicit none
    class(Scalar3D),intent(in) :: this
    type(Scalar3D),intent(inout) :: that

    call this%interp%ScalarGridInterp_3D(this%interior, &
                                         that%interior, &
                                         this%nVar, &
                                         this%nElem)

  endsubroutine GridInterp_Scalar3D

  subroutine Gradient_Scalar3D(this,df)
    implicit none
    class(Scalar3D),intent(in) :: this
    type(Vector3D),intent(inout) :: df

    call this%interp%ScalarGradient_3D(this%interior, &
                                       df%interior, &
                                       this%nVar, &
                                       this%nElem)

  endsubroutine Gradient_Scalar3D

  subroutine WriteHDF5_MPI_Scalar3D(this,fileId,group,elemoffset,nglobalelem)
    implicit none
    class(Scalar3D),intent(in) :: this
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
    globalDims(1:5) = (/this%interp%N+1, &
                        this%interp%N+1, &
                        this%interp%N+1, &
                        this%nVar, &
                        nglobalelem/)

    ! Offsets and dimensions for element boundary data
    bOffset(1:5) = (/0,0,0,0,elemoffset/)
    bGlobalDims(1:5) = (/this%interp%N+1, &
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

  endsubroutine WriteHDF5_MPI_Scalar3D

  subroutine WriteHDF5_Scalar3D(this,fileId,group)
    implicit none
    class(Scalar3D),intent(in) :: this
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

  endsubroutine WriteHDF5_Scalar3D

! -- Vector2D -- !

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

! ! -- Vector3D -- !
  subroutine Init_Vector3D(this,interp,nVar,nElem)
    implicit none
    class(Vector3D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer :: N

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    N = interp%N

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
    type(Tensor3D),intent(inout) :: df

    call this%interp%VectorGradient_3D(this%interior, &
                                       df%interior, &
                                       this%nVar, &
                                       this%nElem)

  endsubroutine Gradient_Vector3D

  subroutine Divergence_Vector3D(this,that)
    implicit none
    class(Vector3D),intent(in) :: this
    type(Scalar3D),intent(inout) :: that

    call this%interp%VectorDivergence_3D(this%interior, &
                                         that%interior, &
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

! ! -- Tensor2D -- !

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

! ! -- Tensor3D -- !

  subroutine Init_Tensor3D(this,interp,nVar,nElem)
    implicit none
    class(Tensor3D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer :: N

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    N = interp%N

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar,1:3,1:3), &
             this%boundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar,1:3,1:3), &
             this%extBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar,1:3,1:3))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:9*nVar))

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)

  endsubroutine Init_Tensor3D

  subroutine Free_Tensor3D(this)
    implicit none
    class(Tensor3D),intent(inout) :: this

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

  endsubroutine Free_Tensor3D

  subroutine BoundaryInterp_Tensor3D(this)
    implicit none
    class(Tensor3D),intent(inout) :: this

    call this%interp%TensorBoundaryInterp_3D(this%interior, &
                                             this%boundary, &
                                             this%nVar, &
                                             this%nElem)

  endsubroutine BoundaryInterp_Tensor3D

  subroutine Divergence_Tensor3D(this,that)
    implicit none
    class(Tensor3D),intent(in) :: this
    class(Vector3D),intent(inout) :: that

    call this%interp%TensorDivergence_3D(this%interior, &
                                         that%interior, &
                                         this%nVar, &
                                         this%nElem)

  endsubroutine Divergence_Tensor3D

  subroutine DGDivergence_Tensor3D(this,that)
    implicit none
    class(Tensor3D),intent(in) :: this
    class(Vector3D),intent(inout) :: that

    call this%interp%TensorDGDivergence_3D(this%interior, &
                                           this%boundary, &
                                           that%interior, &
                                           this%nVar, &
                                           this%nElem)

  endsubroutine DGDivergence_Tensor3D

  subroutine Determinant_Tensor3D(this,that)
    implicit none
    class(Tensor3D),intent(in) :: this
    type(Scalar3D),intent(inout) :: that
    ! Local
    integer :: iEl,iVar,i,j,k

    do iEl = 1,this%nElem
      do iVar = 1,this%nVar
        do k = 1,this%interp%N+1
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1

              that%interior(i,j,k,iEl,iVar) = &
                this%interior(i,j,k,iEl,iVar,1,1)* &
                (this%interior(i,j,k,iEl,iVar,2,2)* &
                 this%interior(i,j,k,iEl,iVar,3,3)- &
                 this%interior(i,j,k,iEl,iVar,2,3)* &
                 this%interior(i,j,k,iEl,iVar,3,2))- &
                this%interior(i,j,k,iEl,iVar,2,1)* &
                (this%interior(i,j,k,iEl,iVar,1,2)* &
                 this%interior(i,j,k,iEl,iVar,3,3)- &
                 this%interior(i,j,k,iEl,iVar,1,3)* &
                 this%interior(i,j,k,iEl,iVar,3,2))+ &
                this%interior(i,j,k,iEl,iVar,3,1)* &
                (this%interior(i,j,k,iEl,iVar,1,2)* &
                 this%interior(i,j,k,iEl,iVar,2,3)- &
                 this%interior(i,j,k,iEl,iVar,1,3)* &
                 this%interior(i,j,k,iEl,iVar,2,2))

            enddo
          enddo
        enddo
      enddo
    enddo

  endsubroutine Determinant_Tensor3D

endmodule SELF_Data
