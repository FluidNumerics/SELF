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

  end type SELF_DataObj

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

    procedure,public :: UpdateDevice => UpdateDevice_Scalar1D

    generic,public :: BoundaryInterp => BoundaryInterp_Scalar1D_cpu,BoundaryInterp_Scalar1D_gpu
    procedure,private :: BoundaryInterp_Scalar1D_cpu
    procedure,private :: BoundaryInterp_Scalar1D_gpu

    generic,public :: GridInterp => GridInterp_Scalar1D_cpu,GridInterp_Scalar1D_gpu
    procedure,private :: GridInterp_Scalar1D_cpu
    procedure,private :: GridInterp_Scalar1D_gpu

    generic,public :: Derivative => Derivative_Scalar1D_cpu,Derivative_Scalar1D_gpu
    procedure,private :: Derivative_Scalar1D_cpu
    procedure,private :: Derivative_Scalar1D_gpu

    ! GENERIC,PUBLIC :: WriteHDF5 => WriteHDF5_Scalar1D, WriteHDF5_MPI_Scalar1D
    ! PROCEDURE, PRIVATE :: WriteHDF5_Scalar1D
    ! PROCEDURE, PRIVATE :: WriteHDF5_MPI_Scalar1D

  end type Scalar1D

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
    procedure,public :: UpdateDevice => UpdateDevice_Scalar2D

    generic,public :: BoundaryInterp => BoundaryInterp_Scalar2D_cpu,BoundaryInterp_Scalar2D_gpu
    procedure,private :: BoundaryInterp_Scalar2D_cpu
    procedure,private :: BoundaryInterp_Scalar2D_gpu

    generic,public :: GridInterp => GridInterp_Scalar2D_cpu,GridInterp_Scalar2D_gpu
    procedure,private :: GridInterp_Scalar2D_cpu
    procedure,private :: GridInterp_Scalar2D_gpu

    generic,public :: Gradient => Gradient_Scalar2D_cpu,Gradient_Scalar2D_gpu
    procedure,private :: Gradient_Scalar2D_cpu
    procedure,private :: Gradient_Scalar2D_gpu

    ! GENERIC,PUBLIC :: WriteHDF5 => WriteHDF5_MPI_Scalar2D, WriteHDF5_Scalar2D
    ! PROCEDURE, PRIVATE :: WriteHDF5_MPI_Scalar2D
    ! PROCEDURE, PRIVATE :: WriteHDF5_Scalar2D

  end type Scalar2D

  type,extends(SELF_DataObj),public :: Scalar3D

    real(prec),pointer,dimension(:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:) :: extBoundary
    real(prec),pointer,dimension(:,:,:,:,:) :: avgBoundary
    real(prec),pointer,dimension(:,:,:,:,:) :: jumpBoundary

    real(prec),pointer,private,dimension(:,:,:,:,:) :: interpWork1
    real(prec),pointer,private,dimension(:,:,:,:,:) :: interpWork2

  contains

    procedure,public :: Init => Init_Scalar3D
    procedure,public :: Free => Free_Scalar3D
    procedure,public :: UpdateDevice => UpdateDevice_Scalar3D

    generic,public :: BoundaryInterp => BoundaryInterp_Scalar3D_cpu,BoundaryInterp_Scalar3D_gpu
    procedure,private :: BoundaryInterp_Scalar3D_cpu
    procedure,private :: BoundaryInterp_Scalar3D_gpu

    generic,public :: GridInterp => GridInterp_Scalar3D_cpu,GridInterp_Scalar3D_gpu
    procedure,private :: GridInterp_Scalar3D_cpu
    procedure,private :: GridInterp_Scalar3D_gpu

    generic,public :: Gradient => Gradient_Scalar3D_cpu,Gradient_Scalar3D_gpu
    procedure,private :: Gradient_Scalar3D_cpu
    procedure,private :: Gradient_Scalar3D_gpu

    ! GENERIC,PUBLIC :: WriteHDF5 => WriteHDF5_MPI_Scalar3D, WriteHDF5_Scalar3D
    ! PROCEDURE, PRIVATE :: WriteHDF5_MPI_Scalar3D
    ! PROCEDURE, PRIVATE :: WriteHDF5_Scalar3D

  end type Scalar3D

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
    procedure,public :: UpdateDevice => UpdateDevice_Vector2D

    generic,public :: BoundaryInterp => BoundaryInterp_Vector2D_cpu,BoundaryInterp_Vector2D_gpu
    procedure,private :: BoundaryInterp_Vector2D_cpu
    procedure,private :: BoundaryInterp_Vector2D_gpu

    generic,public :: GridInterp => GridInterp_Vector2D_cpu!,GridInterp_Vector2D_gpu
    procedure,private :: GridInterp_Vector2D_cpu
    !procedure,private :: GridInterp_Vector2D_gpu

    generic,public :: Gradient => Gradient_Vector2D_gpu,Gradient_Vector2D_cpu
    procedure,private :: Gradient_Vector2D_gpu
    procedure,private :: Gradient_Vector2D_cpu

    generic,public :: Divergence => Divergence_Vector2D_gpu,Divergence_Vector2D_cpu
    procedure,private :: Divergence_Vector2D_gpu
    procedure,private :: Divergence_Vector2D_cpu

    ! GENERIC,PUBLIC :: SetEquation => SetEquation_Vector2D
    ! PROCEDURE,PRIVATE :: SetEquation_Vector2D

    ! GENERIC,PUBLIC :: WriteHDF5 => WriteHDF5_MPI_Vector2D, WriteHDF5_Vector2D
    ! PROCEDURE, PRIVATE :: WriteHDF5_MPI_Vector2D
    ! PROCEDURE, PRIVATE :: WriteHDF5_Vector2D

  end type Vector2D

  type,extends(SELF_DataObj),public :: Vector3D

    real(prec),pointer,dimension(:,:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:,:) :: extBoundary
    real(prec),pointer,dimension(:,:,:,:,:,:) :: avgBoundary
    real(prec),pointer,dimension(:,:,:,:,:) :: boundaryNormal

  contains

    procedure,public :: Init => Init_Vector3D
    procedure,public :: Free => Free_Vector3D
    procedure,public :: UpdateDevice => UpdateDevice_Vector3D

    generic,public :: BoundaryInterp => BoundaryInterp_Vector3D_cpu,BoundaryInterp_Vector3D_gpu
    procedure,private :: BoundaryInterp_Vector3D_cpu
    procedure,private :: BoundaryInterp_Vector3D_gpu

    generic,public :: GridInterp => GridInterp_Vector3D_cpu!,GridInterp_Vector3D_gpu
    procedure,private :: GridInterp_Vector3D_cpu
    !procedure,private :: GridInterp_Vector3D_gpu

    generic,public :: Gradient => Gradient_Vector3D_gpu,Gradient_Vector3D_cpu
    procedure,private :: Gradient_Vector3D_gpu
    procedure,private :: Gradient_Vector3D_cpu

    generic,public :: Divergence => Divergence_Vector3D_gpu,Divergence_Vector3D_cpu
    procedure,private :: Divergence_Vector3D_gpu
    procedure,private :: Divergence_Vector3D_cpu

    ! GENERIC,PUBLIC :: SetEquation => SetEquation_Vector3D
    ! PROCEDURE,PRIVATE :: SetEquation_Vector3D

    ! GENERIC,PUBLIC :: WriteHDF5 => WriteHDF5_MPI_Vector3D, WriteHDF5_Vector3D
    ! PROCEDURE, PRIVATE :: WriteHDF5_MPI_Vector3D
    ! PROCEDURE, PRIVATE :: WriteHDF5_Vector3D

  end type Vector3D
! ! ---------------------- Tensors ---------------------- !

  type,extends(SELF_DataObj),public :: Tensor2D

    real(prec),pointer,dimension(:,:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:,:) :: extBoundary

  contains

    procedure,public :: Init => Init_Tensor2D
    procedure,public :: Free => Free_Tensor2D
    procedure,public :: UpdateDevice => UpdateDevice_Tensor2D

    generic,public :: BoundaryInterp => BoundaryInterp_Tensor2D_cpu,BoundaryInterp_Tensor2D_gpu
    procedure,private :: BoundaryInterp_Tensor2D_cpu
    procedure,private :: BoundaryInterp_Tensor2D_gpu

    generic,public :: Divergence => Divergence_Tensor2D_gpu,Divergence_Tensor2D_cpu
    procedure,private :: Divergence_Tensor2D_gpu
    procedure,private :: Divergence_Tensor2D_cpu

    procedure,public :: Determinant => Determinant_Tensor2D

  end type Tensor2D

  type,extends(SELF_DataObj),public :: Tensor3D

    real(prec),pointer,dimension(:,:,:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:,:,:) :: extBoundary

  contains

    procedure,public :: Init => Init_Tensor3D
    procedure,public :: Free => Free_Tensor3D
    procedure,public :: UpdateDevice => UpdateDevice_Tensor3D

    generic,public :: BoundaryInterp => BoundaryInterp_Tensor3D_cpu,BoundaryInterp_Tensor3D_gpu
    procedure,private :: BoundaryInterp_Tensor3D_cpu
    procedure,private :: BoundaryInterp_Tensor3D_gpu

    procedure,public :: Determinant => Determinant_Tensor3D

  end type Tensor3D

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

    call this % meta(ivar) % SetName(name)

  end subroutine SetName_DataObj

  subroutine SetDescription_DataObj(this,ivar,description)
    !! Set the description of the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: description

    call this % meta(ivar) % SetDescription(description)

  end subroutine SetDescription_DataObj

  subroutine SetUnits_DataObj(this,ivar,units)
    !! Set the units of the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: units

    call this % meta(ivar) % SetUnits(units)

  end subroutine SetUnits_DataObj

  subroutine SetEquation_DataObj(this,ivar,eqnChar)
    !! Sets the equation parser for the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: eqnChar

    this % eqn(ivar) = EquationParser(trim(eqnChar), &
                                      (/'x','y','z','t'/))

  end subroutine SetEquation_DataObj

! -- Scalar1D -- !

  subroutine Init_Scalar1D(this,interp,nVar,nElem)
    implicit none
    class(Scalar1D),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this % interp => interp
    this % nVar = nVar
    this % nElem = nElem

    call hipcheck(hipMallocManaged(this % interior,interp % N + 1,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % boundary,2,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % extBoundary,2,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % avgBoundary,2,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % jumpBoundary,2,nelem,nvar,hipMemAttachGlobal))

    allocate (this % meta(1:nVar))
    allocate (this % eqn(1:nVar))

  end subroutine Init_Scalar1D

  subroutine Free_Scalar1D(this)
    implicit none
    class(Scalar1D),intent(inout) :: this

    this % interp => null()
    call hipcheck(hipFree(this % interior))
    call hipcheck(hipFree(this % boundary))
    call hipcheck(hipFree(this % extBoundary))
    call hipcheck(hipFree(this % avgBoundary))
    call hipcheck(hipFree(this % jumpBoundary))
    deallocate (this % meta)
    deallocate (this % eqn)

  end subroutine Free_Scalar1D

  subroutine UpdateDevice_Scalar1D(this)
    implicit none
    class(Scalar1D),intent(inout) :: this

    call hipcheck(hipMemPrefetchAsync(c_loc(this % interior),sizeof(this % interior),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % boundary),sizeof(this % boundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % extBoundary),sizeof(this % extBoundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % avgBoundary),sizeof(this % avgBoundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % jumpBoundary),sizeof(this % jumpBoundary),0,c_null_ptr))

  end subroutine UpdateDevice_Scalar1D

  subroutine BoundaryInterp_Scalar1D_cpu(this)
    implicit none
    class(Scalar1D),intent(inout) :: this

    call this % interp % ScalarBoundaryInterp_1D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem)

  end subroutine BoundaryInterp_Scalar1D_cpu

  subroutine BoundaryInterp_Scalar1D_gpu(this,handle)
    implicit none
    class(Scalar1D),intent(inout) :: this
    type(c_ptr),intent(inout) :: handle

    call this % interp % ScalarBoundaryInterp_1D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem, &
                                                 handle)

  end subroutine BoundaryInterp_Scalar1D_gpu

  subroutine GridInterp_Scalar1D_cpu(this,that)
    implicit none
    class(Scalar1D),intent(in) :: this
    type(Scalar1D),intent(inout) :: that

    call this % interp % ScalarGridInterp_1D(this % interior, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem)

  end subroutine GridInterp_Scalar1D_cpu

  subroutine GridInterp_Scalar1D_gpu(this,that,handle)
    implicit none
    class(Scalar1D),intent(in) :: this
    type(Scalar1D),intent(inout) :: that
    type(c_ptr),intent(inout) :: handle

    call this % interp % ScalarGridInterp_1D(this % interior, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem, &
                                             handle)

  end subroutine GridInterp_Scalar1D_gpu

  subroutine Derivative_Scalar1D_cpu(this,that)
    implicit none
    class(Scalar1D),intent(in) :: this
    type(Scalar1D),intent(inout) :: that

    call this % interp % Derivative_1D(this % interior, &
                                       that % interior, &
                                       this % nVar, &
                                       this % nElem)

  end subroutine Derivative_Scalar1D_cpu

  subroutine Derivative_Scalar1D_gpu(this,that,handle)
    implicit none
    class(Scalar1D),intent(in) :: this
    type(Scalar1D),intent(inout) :: that
    type(c_ptr),intent(inout) :: handle

    call this % interp % Derivative_1D(this % interior, &
                                       that % interior, &
                                       this % nVar, &
                                       this % nElem, &
                                       handle)

  end subroutine Derivative_Scalar1D_gpu

  ! SUBROUTINE WriteHDF5_MPI_Scalar1D(this,fileId,group,elemoffset,nglobalelem)
  !   IMPLICIT NONE
  !   CLASS(Scalar1D), INTENT(in) :: this
  !   CHARACTER(*), INTENT(in) :: group
  !   INTEGER(HID_T), INTENT(in) :: fileId
  !   INTEGER, INTENT(in) :: elemoffset
  !   INTEGER, INTENT(in) :: nglobalelem
  !   ! Local
  !   INTEGER(HID_T) :: offset(1:3)
  !   INTEGER(HID_T) :: bOffset(1:3)
  !   INTEGER(HID_T) :: globalDims(1:3)
  !   INTEGER(HID_T) :: bGlobalDims(1:3)
  !   INTEGER :: ivar

  !     offset(1:3) = (/0,0,elemoffset/)
  !     globalDims(1:3) = (/this % interp % N + 1, &
  !                         this % nVar, &
  !                         nGlobalElem/)

  !     ! Offsets and dimensions for element boundary data
  !     bOffset(1:3) = (/0,0,elemoffset/)
  !     bGlobalDims(1:3) = (/this % nVar, &
  !                          2, &
  !                          nGlobalElem/)

  !     CALL CreateGroup_HDF5(fileId,TRIM(group))

  !     DO ivar = 1, this % nVar
  !       CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
  !     ENDDO

  !     CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
  !                          this % interior,offset,globalDims)

  !     CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
  !                          this % boundary,bOffset,bGlobalDims)

  ! END SUBROUTINE WriteHDF5_MPI_Scalar1D

  ! SUBROUTINE WriteHDF5_Scalar1D(this,fileId,group)
  !   IMPLICIT NONE
  !   CLASS(Scalar1D), INTENT(in) :: this
  !   INTEGER(HID_T), INTENT(in) :: fileId
  !   CHARACTER(*), INTENT(in) :: group
  !   ! Local
  !   INTEGER :: ivar

  !     CALL CreateGroup_HDF5(fileId,TRIM(group))

  !     DO ivar = 1, this % nVar
  !       CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
  !     ENDDO

  !     CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
  !                          this % interior)

  !     CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
  !                          this % boundary)

  ! END SUBROUTINE WriteHDF5_Scalar1D

! ! -- Scalar2D -- !

  subroutine Init_Scalar2D(this,interp,nVar,nElem)
    implicit none
    class(Scalar2D),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this % interp => interp
    this % nVar = nVar
    this % nElem = nElem

    call hipcheck(hipMallocManaged(this % interior,interp % N + 1,interp % N + 1,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % interpWork,interp % M + 1,interp % N + 1,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % boundary,interp % N + 1,4,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % extBoundary,interp % N + 1,4,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % avgBoundary,interp % N + 1,4,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % jumpBoundary,interp % N + 1,4,nelem,nvar,hipMemAttachGlobal))

    allocate (this % meta(1:nVar))
    allocate (this % eqn(1:nVar))

  end subroutine Init_Scalar2D

  subroutine Free_Scalar2D(this)
    implicit none
    class(Scalar2D),intent(inout) :: this

    this % nVar = 0
    this % nElem = 0
    this % interp => null()
    call hipcheck(hipFree(this % interior))
    call hipcheck(hipFree(this % interpWork))
    call hipcheck(hipFree(this % boundary))
    call hipcheck(hipFree(this % extBoundary))
    call hipcheck(hipFree(this % avgBoundary))
    call hipcheck(hipFree(this % jumpBoundary))
    deallocate (this % meta)
    deallocate (this % eqn)

  end subroutine Free_Scalar2D

  subroutine UpdateDevice_Scalar2D(this)
    implicit none
    class(Scalar2D),intent(inout) :: this

    call hipcheck(hipMemPrefetchAsync(c_loc(this % interior),sizeof(this % interior),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % boundary),sizeof(this % boundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % extBoundary),sizeof(this % extBoundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % avgBoundary),sizeof(this % avgBoundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % jumpBoundary),sizeof(this % jumpBoundary),0,c_null_ptr))

  end subroutine UpdateDevice_Scalar2D

  subroutine BoundaryInterp_Scalar2D_cpu(this)
    implicit none
    class(Scalar2D),intent(inout) :: this

    call this % interp % ScalarBoundaryInterp_2D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem)

  end subroutine BoundaryInterp_Scalar2D_cpu

  subroutine BoundaryInterp_Scalar2D_gpu(this,handle)
    implicit none
    class(Scalar2D),intent(inout) :: this
    type(c_ptr),intent(in) :: handle

    call this % interp % ScalarBoundaryInterp_2D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem, &
                                                 handle)

  end subroutine BoundaryInterp_Scalar2D_gpu

  subroutine GridInterp_Scalar2D_cpu(this,that)
    implicit none
    class(Scalar2D),intent(in) :: this
    type(Scalar2D),intent(inout) :: that

    call this % interp % ScalarGridInterp_2D(this % interior, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem)

  end subroutine GridInterp_Scalar2D_cpu

  subroutine GridInterp_Scalar2D_gpu(this,that,handle)
    implicit none
    class(Scalar2D),intent(inout) :: this
    type(Scalar2D),intent(inout) :: that
    type(c_ptr),intent(inout) :: handle

    call this % interp % ScalarGridInterp_2D(this % interior, &
                                             this % interpWork, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem, &
                                             handle)

  end subroutine GridInterp_Scalar2D_gpu

  subroutine Gradient_Scalar2D_cpu(this,df)
    implicit none
    class(Scalar2D),intent(in) :: this
    type(Vector2D),intent(inout) :: df

    call this % interp % ScalarGradient_2D(this % interior, &
                                           df % interior, &
                                           this % nVar, &
                                           this % nElem)

  end subroutine Gradient_Scalar2D_cpu

  subroutine Gradient_Scalar2D_gpu(this,df,handle)
    implicit none
    class(Scalar2D),intent(in) :: this
    type(Vector2D),intent(inout) :: df
    type(c_ptr),intent(inout) :: handle

    call this % interp % ScalarGradient_2D(this % interior, &
                                           df % interior, &
                                           this % nVar, &
                                           this % nElem, &
                                           handle)

  end subroutine Gradient_Scalar2D_gpu

!   SUBROUTINE WriteHDF5_MPI_Scalar2D(this,fileId,group,elemoffset,nglobalelem)
!     IMPLICIT NONE
!     CLASS(Scalar2D), INTENT(in) :: this
!     CHARACTER(*), INTENT(in) :: group
!     INTEGER(HID_T), INTENT(in) :: fileId
!     INTEGER, INTENT(in) :: elemoffset
!     INTEGER, INTENT(in) :: nglobalelem
!     ! Local
!     INTEGER(HID_T) :: offset(1:4)
!     INTEGER(HID_T) :: bOffset(1:4)
!     INTEGER(HID_T) :: globalDims(1:4)
!     INTEGER(HID_T) :: bGlobalDims(1:4)
!     INTEGER :: ivar

!       offset(1:4) = (/0,0,0,elemoffset/)
!       globalDims(1:4) = (/this % interp % N + 1, &
!                           this % interp % N + 1, &
!                           this % nVar, &
!                           nglobalelem/)

!       ! Offsets and dimensions for element boundary data
!       bOffset(1:4) = (/0,0,0,elemoffset/)
!       bGlobalDims(1:4) = (/this % interp % N + 1, &
!                            this % nVar, &
!                            4, &
!                            nglobalelem/)

!       CALL CreateGroup_HDF5(fileId,TRIM(group))

!       DO ivar = 1, this % nVar
!         CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
!       ENDDO

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
!                            this % interior,offset,globalDims)

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
!                            this % boundary,bOffset,bGlobalDims)

!   END SUBROUTINE WriteHDF5_MPI_Scalar2D

!   SUBROUTINE WriteHDF5_Scalar2D(this,fileId,group)
!     IMPLICIT NONE
!     CLASS(Scalar2D), INTENT(in) :: this
!     INTEGER(HID_T), INTENT(in) :: fileId
!     CHARACTER(*), INTENT(in) :: group
!     ! Local
!     INTEGER :: ivar

!       CALL CreateGroup_HDF5(fileId,TRIM(group))

!       DO ivar = 1, this % nVar
!         CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
!       ENDDO

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
!                            this % interior)

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
!                            this % boundary)

!   END SUBROUTINE WriteHDF5_Scalar2D

! ! -- Scalar3D -- !

  subroutine Init_Scalar3D(this,interp,nVar,nElem)
    implicit none
    class(Scalar3D),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this % interp => interp
    this % nVar = nVar
    this % nElem = nElem

    call hipcheck(hipMallocManaged(this % interior,interp % N + 1,interp % N + 1,interp % N + 1,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % interpWork1,interp % M + 1,interp % N + 1,interp % N + 1,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % interpWork2,interp % M + 1,interp % M + 1,interp % N + 1,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % boundary,interp % N + 1,interp % N + 1,6,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % extBoundary,interp % N + 1,interp % N + 1,6,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % avgBoundary,interp % N + 1,interp % N + 1,6,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % jumpBoundary,interp % N + 1,interp % N + 1,6,nelem,nvar,hipMemAttachGlobal))

    allocate (this % meta(1:nVar))
    allocate (this % eqn(1:nVar))

  end subroutine Init_Scalar3D

  subroutine Free_Scalar3D(this)
    implicit none
    class(Scalar3D),intent(inout) :: this

    this % nVar = 0
    this % nElem = 0
    this % interp => null()
    call hipcheck(hipFree(this % interior))
    call hipcheck(hipFree(this % interpWork1))
    call hipcheck(hipFree(this % interpWork2))
    call hipcheck(hipFree(this % boundary))
    call hipcheck(hipFree(this % extBoundary))
    call hipcheck(hipFree(this % avgBoundary))
    call hipcheck(hipFree(this % jumpBoundary))
    deallocate (this % meta)
    deallocate (this % eqn)

  end subroutine Free_Scalar3D

  subroutine UpdateDevice_Scalar3D(this)
    implicit none
    class(Scalar3D),intent(inout) :: this

    call hipcheck(hipMemPrefetchAsync(c_loc(this % interior),sizeof(this % interior),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % boundary),sizeof(this % boundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % extBoundary),sizeof(this % extBoundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % avgBoundary),sizeof(this % avgBoundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % jumpBoundary),sizeof(this % jumpBoundary),0,c_null_ptr))

  end subroutine UpdateDevice_Scalar3D

  subroutine BoundaryInterp_Scalar3D_cpu(this)
    implicit none
    class(Scalar3D),intent(inout) :: this

    call this % interp % ScalarBoundaryInterp_3D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem)

  end subroutine BoundaryInterp_Scalar3D_cpu

  subroutine BoundaryInterp_Scalar3D_gpu(this,handle)
    implicit none
    class(Scalar3D),intent(inout) :: this
    type(c_ptr),intent(in) :: handle

    call this % interp % ScalarBoundaryInterp_3D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem, &
                                                 handle)

  end subroutine BoundaryInterp_Scalar3D_gpu

  subroutine GridInterp_Scalar3D_cpu(this,that)
    implicit none
    class(Scalar3D),intent(in) :: this
    type(Scalar3D),intent(inout) :: that

    call this % interp % ScalarGridInterp_3D(this % interior, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem)

  end subroutine GridInterp_Scalar3D_cpu

  subroutine GridInterp_Scalar3D_gpu(this,that,handle)
    implicit none
    class(Scalar3D),intent(inout) :: this
    type(Scalar3D),intent(inout) :: that
    type(c_ptr),intent(inout) :: handle

    call this % interp % ScalarGridInterp_3D(this % interior, &
                                             this % interpWork1, &
                                             this % interpWork2, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem, &
                                             handle)

  end subroutine GridInterp_Scalar3D_gpu

  subroutine Gradient_Scalar3D_cpu(this,df)
    implicit none
    class(Scalar3D),intent(in) :: this
    type(Vector3D),intent(inout) :: df

    call this % interp % ScalarGradient_3D(this % interior, &
                                           df % interior, &
                                           this % nVar, &
                                           this % nElem)

  end subroutine Gradient_Scalar3D_cpu

  subroutine Gradient_Scalar3D_gpu(this,df,handle)
    implicit none
    class(Scalar3D),intent(in) :: this
    type(Vector3D),intent(inout) :: df
    type(c_ptr),intent(inout) :: handle

    call this % interp % ScalarGradient_3D(this % interior, &
                                           df % interior, &
                                           this % nVar, &
                                           this % nElem, &
                                           handle)

  end subroutine Gradient_Scalar3D_gpu

!   SUBROUTINE WriteHDF5_MPI_Scalar3D(this,fileId,group,elemoffset,nglobalelem)
!     IMPLICIT NONE
!     CLASS(Scalar3D), INTENT(in) :: this
!     CHARACTER(*), INTENT(in) :: group
!     INTEGER(HID_T), INTENT(in) :: fileId
!     INTEGER, INTENT(in) :: elemoffset
!     INTEGER, INTENT(in) :: nglobalelem
!     ! Local
!     INTEGER(HID_T) :: offset(1:5)
!     INTEGER(HID_T) :: bOffset(1:5)
!     INTEGER(HID_T) :: globalDims(1:5)
!     INTEGER(HID_T) :: bGlobalDims(1:5)
!     INTEGER :: ivar

!       offset(1:5) = (/0,0,0,0,elemoffset/)
!       globalDims(1:5) = (/this % interp % N + 1, &
!                           this % interp % N + 1, &
!                           this % interp % N + 1, &
!                           this % nVar, &
!                           nglobalelem/)

!       ! Offsets and dimensions for element boundary data
!       bOffset(1:5) = (/0,0,0,0,elemoffset/)
!       bGlobalDims(1:5) = (/this % interp % N + 1, &
!                            this % interp % N + 1, &
!                            this % nVar, &
!                            6, &
!                            nglobalelem/)

!       CALL CreateGroup_HDF5(fileId,TRIM(group))

!       DO ivar = 1, this % nVar
!         CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
!       ENDDO

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
!                            this % interior,offset,globalDims)

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
!                            this % boundary,bOffset,bGlobalDims)

!   END SUBROUTINE WriteHDF5_MPI_Scalar3D

!   SUBROUTINE WriteHDF5_Scalar3D(this,fileId,group)
!     IMPLICIT NONE
!     CLASS(Scalar3D), INTENT(in) :: this
!     INTEGER(HID_T), INTENT(in) :: fileId
!     CHARACTER(*), INTENT(in) :: group
!     ! Local
!     INTEGER :: ivar

!       CALL CreateGroup_HDF5(fileId,TRIM(group))

!       DO ivar = 1, this % nVar
!         CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
!       ENDDO

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
!                            this % interior)

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
!                            this % boundary)

!   END SUBROUTINE WriteHDF5_Scalar3D

! -- Vector2D -- !

  subroutine Init_Vector2D(this,interp,nVar,nElem)
    implicit none
    class(Vector2D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer :: N

    this % interp => interp
    this % nVar = nVar
    this % nElem = nElem
    N = interp % N

    call hipcheck(hipMallocManaged(this % interior,interp % N + 1,interp % N + 1,nelem,nvar,2,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % boundary,interp % N + 1,4,nelem,nvar,2,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % extBoundary,interp % N + 1,4,nelem,nvar,2,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % avgBoundary,interp % N + 1,4,nelem,nvar,2,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % boundaryNormal,interp % N + 1,4,nelem,nvar,hipMemAttachGlobal))

    allocate (this % meta(1:nVar))
    allocate (this % eqn(1:2*nVar))

  end subroutine Init_Vector2D

  subroutine Free_Vector2D(this)
    implicit none
    class(Vector2D),intent(inout) :: this

    this % interp => null()
    this % nVar = 0
    this % nElem = 0

    call hipcheck(hipFree(this % interior))
    call hipcheck(hipFree(this % boundary))
    call hipcheck(hipFree(this % boundaryNormal))
    call hipcheck(hipFree(this % extBoundary))
    call hipcheck(hipFree(this % avgBoundary))

    deallocate (this % meta)
    deallocate (this % eqn)

  end subroutine Free_Vector2D

!   SUBROUTINE SetEquation_Vector2D(this,idir,ivar,eqnChar)
!     !! Sets the equation parser for the `idir` direction and `ivar-th` variable
!     IMPLICIT NONE
!     CLASS(Vector2D),INTENT(inout) :: this
!     INTEGER,INTENT(in) :: idir,ivar
!     CHARACTER(*),INTENT(in) :: eqnChar

!     this % eqn(idir+2*(ivar-1)) = EquationParser( TRIM(eqnChar), &
!                                               (/'x','y','z','t'/) )

!   END SUBROUTINE SetEquation_Vector2D

  subroutine UpdateDevice_Vector2D(this)
    implicit none
    class(Vector2D),intent(inout) :: this

    call hipcheck(hipMemPrefetchAsync(c_loc(this % interior),sizeof(this % interior),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % boundary),sizeof(this % boundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % boundaryNormal),sizeof(this % boundaryNormal),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % extBoundary),sizeof(this % extBoundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % avgBoundary),sizeof(this % avgBoundary),0,c_null_ptr))

  end subroutine UpdateDevice_Vector2D

  subroutine GridInterp_Vector2D_cpu(this,that)
    implicit none
    class(Vector2D),intent(in) :: this
    type(Vector2D),intent(inout) :: that

    call this % interp % VectorGridInterp_2D(this % interior, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem)

  end subroutine GridInterp_Vector2D_cpu

  subroutine BoundaryInterp_Vector2D_cpu(this)
    implicit none
    class(Vector2D),intent(inout) :: this

    call this % interp % VectorBoundaryInterp_2D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem)

  end subroutine BoundaryInterp_Vector2D_cpu

  subroutine BoundaryInterp_Vector2D_gpu(this,handle)
    implicit none
    class(Vector2D),intent(inout) :: this
    type(c_ptr),intent(in) :: handle

    call this % interp % VectorBoundaryInterp_2D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem, &
                                                 handle)

  end subroutine BoundaryInterp_Vector2D_gpu

  subroutine Gradient_Vector2D_cpu(this,df)
    implicit none
    class(Vector2D),intent(in) :: this
    type(Tensor2D),intent(inout) :: df

    call this % interp % VectorGradient_2D(this % interior, &
                                           df % interior, &
                                           this % nVar, &
                                           this % nElem)

  end subroutine Gradient_Vector2D_cpu

  subroutine Gradient_Vector2D_gpu(this,df,handle)
    implicit none
    class(Vector2D),intent(in) :: this
    type(Tensor2D),intent(inout) :: df
    type(c_ptr),intent(inout) :: handle

    call this % interp % VectorGradient_2D(this % interior, &
                                           df % interior, &
                                           this % nVar, &
                                           this % nElem, &
                                           handle)

  end subroutine Gradient_Vector2D_gpu

  subroutine Divergence_Vector2D_cpu(this,that)
    implicit none
    class(Vector2D),intent(in) :: this
    type(Scalar2D),intent(inout) :: that

    call this % interp % VectorDivergence_2D(this % interior, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem)

  end subroutine Divergence_Vector2D_cpu

  subroutine Divergence_Vector2D_gpu(this,that,handle)
    implicit none
    class(Vector2D),intent(in) :: this
    type(Scalar2D),intent(inout) :: that
    type(c_ptr),intent(inout) :: handle

    call this % interp % VectorDivergence_2D(this % interior, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem, &
                                             handle)

  end subroutine Divergence_Vector2D_gpu

!   SUBROUTINE WriteHDF5_MPI_Vector2D(this,fileId,group,elemoffset,nglobalelem)
!     IMPLICIT NONE
!     CLASS(Vector2D), INTENT(in) :: this
!     CHARACTER(*), INTENT(in) :: group
!     INTEGER(HID_T), INTENT(in) :: fileId
!     INTEGER, INTENT(in) :: elemoffset
!     INTEGER, INTENT(in) :: nglobalelem
!     ! Local
!     INTEGER(HID_T) :: offset(1:5)
!     INTEGER(HID_T) :: bOffset(1:5)
!     INTEGER(HID_T) :: globalDims(1:5)
!     INTEGER(HID_T) :: bGlobalDims(1:5)
!     INTEGER :: ivar

!       offset(1:5) = (/0,0,0,0,elemoffset/)
!       globalDims(1:5) = (/2, &
!                           this % interp % N + 1, &
!                           this % interp % N + 1, &
!                           this % nVar, &
!                           nglobalelem/)

!       ! Offsets and dimensions for element boundary data
!       bOffset(1:5) = (/0,0,0,0,elemoffset/)
!       bGlobalDims(1:5) = (/2, &
!                            this % interp % N + 1, &
!                            this % nVar, &
!                            4, &
!                            nglobalelem/)

!       CALL CreateGroup_HDF5(fileId,TRIM(group))

!       DO ivar = 1, this % nVar
!         CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
!       ENDDO

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
!                            this % interior,offset,globalDims)

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
!                            this % boundary,bOffset,bGlobalDims)

!   END SUBROUTINE WriteHDF5_MPI_Vector2D

!   SUBROUTINE WriteHDF5_Vector2D(this,fileId,group)
!     IMPLICIT NONE
!     CLASS(Vector2D), INTENT(in) :: this
!     INTEGER(HID_T), INTENT(in) :: fileId
!     CHARACTER(*), INTENT(in) :: group
!     ! Local
!     INTEGER :: ivar

!       CALL CreateGroup_HDF5(fileId,TRIM(group))

!       DO ivar = 1, this % nVar
!         CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
!       ENDDO

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
!                            this % interior)

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
!                            this % boundary)

!   END SUBROUTINE WriteHDF5_Vector2D

! ! -- Vector3D -- !
  subroutine Init_Vector3D(this,interp,nVar,nElem)
    implicit none
    class(Vector3D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer :: N

    this % interp => interp
    this % nVar = nVar
    this % nElem = nElem
    N = interp % N

    call hipcheck(hipMallocManaged(this % interior,interp % N + 1,interp % N + 1,interp % N + 1,nelem,nvar,3,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % boundary,interp % N + 1,interp % N + 1,6,nelem,nvar,3,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % extBoundary,interp % N + 1,interp % N + 1,6,nelem,nvar,3,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % avgBoundary,interp % N + 1,interp % N + 1,6,nelem,nvar,3,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % boundaryNormal,interp % N + 1,interp % N + 1,6,nelem,nvar,hipMemAttachGlobal))

    allocate (this % meta(1:nVar))
    allocate (this % eqn(1:3*nVar))

  end subroutine Init_Vector3D

  subroutine Free_Vector3D(this)
    implicit none
    class(Vector3D),intent(inout) :: this

    this % interp => null()
    this % nVar = 0
    this % nElem = 0

    call hipcheck(hipFree(this % interior))
    call hipcheck(hipFree(this % boundary))
    call hipcheck(hipFree(this % boundaryNormal))
    call hipcheck(hipFree(this % extBoundary))
    call hipcheck(hipFree(this % avgBoundary))

    deallocate (this % meta)
    deallocate (this % eqn)

  end subroutine Free_Vector3D

!   SUBROUTINE SetEquation_Vector3D(this,idir,ivar,eqnChar)
!     !! Sets the equation parser for the `idir` direction and `ivar-th` variable
!     IMPLICIT NONE
!     CLASS(Vector3D),INTENT(inout) :: this
!     INTEGER,INTENT(in) :: idir,ivar
!     CHARACTER(*),INTENT(in) :: eqnChar

!     this % eqn(idir+2*(ivar-1)) = EquationParser( TRIM(eqnChar), &
!                                               (/'x','y','z','t'/) )

!   END SUBROUTINE SetEquation_Vector3D

  subroutine UpdateDevice_Vector3D(this)
    implicit none
    class(Vector3D),intent(inout) :: this

    call hipcheck(hipMemPrefetchAsync(c_loc(this % interior),sizeof(this % interior),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % boundary),sizeof(this % boundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % boundaryNormal),sizeof(this % boundaryNormal),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % extBoundary),sizeof(this % extBoundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % avgBoundary),sizeof(this % avgBoundary),0,c_null_ptr))
    

  end subroutine UpdateDevice_Vector3D

  subroutine GridInterp_Vector3D_cpu(this,that)
    implicit none
    class(Vector3D),intent(in) :: this
    type(Vector3D),intent(inout) :: that

    call this % interp % VectorGridInterp_3D(this % interior, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem)

  end subroutine GridInterp_Vector3D_cpu

  subroutine BoundaryInterp_Vector3D_cpu(this)
    implicit none
    class(Vector3D),intent(inout) :: this

    call this % interp % VectorBoundaryInterp_3D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem)

  end subroutine BoundaryInterp_Vector3D_cpu

  subroutine BoundaryInterp_Vector3D_gpu(this,handle)
    implicit none
    class(Vector3D),intent(inout) :: this
    type(c_ptr),intent(in) :: handle

    call this % interp % VectorBoundaryInterp_3D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem, &
                                                 handle)

  end subroutine BoundaryInterp_Vector3D_gpu

  subroutine Gradient_Vector3D_cpu(this,df)
    implicit none
    class(Vector3D),intent(in) :: this
    type(Tensor3D),intent(inout) :: df

    call this % interp % VectorGradient_3D(this % interior, &
                                           df % interior, &
                                           this % nVar, &
                                           this % nElem)

  end subroutine Gradient_Vector3D_cpu

  subroutine Gradient_Vector3D_gpu(this,df,handle)
    implicit none
    class(Vector3D),intent(in) :: this
    type(Tensor3D),intent(inout) :: df
    type(c_ptr),intent(inout) :: handle

    call this % interp % VectorGradient_3D(this % interior, &
                                           df % interior, &
                                           this % nVar, &
                                           this % nElem, &
                                           handle)

  end subroutine Gradient_Vector3D_gpu

  subroutine Divergence_Vector3D_cpu(this,that)
    implicit none
    class(Vector3D),intent(in) :: this
    type(Scalar3D),intent(inout) :: that

    call this % interp % VectorDivergence_3D(this % interior, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem)

  end subroutine Divergence_Vector3D_cpu

  subroutine Divergence_Vector3D_gpu(this,that,handle)
    implicit none
    class(Vector3D),intent(in) :: this
    type(Scalar3D),intent(inout) :: that
    type(c_ptr),intent(inout) :: handle

    call this % interp % VectorDivergence_3D(this % interior, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem, &
                                             handle)

  end subroutine Divergence_Vector3D_gpu

!   SUBROUTINE WriteHDF5_MPI_Vector3D(this,fileId,group,elemoffset,nglobalelem)
!     IMPLICIT NONE
!     CLASS(Vector3D), INTENT(in) :: this
!     CHARACTER(*), INTENT(in) :: group
!     INTEGER(HID_T), INTENT(in) :: fileId
!     INTEGER, INTENT(in) :: elemoffset
!     INTEGER, INTENT(in) :: nglobalelem
!     ! Local
!     INTEGER(HID_T) :: offset(1:6)
!     INTEGER(HID_T) :: bOffset(1:6)
!     INTEGER(HID_T) :: globalDims(1:6)
!     INTEGER(HID_T) :: bGlobalDims(1:6)
!     INTEGER :: ivar

!       offset(1:6) = (/0,0,0,0,0,elemoffset/)
!       globalDims(1:6) = (/3, &
!                           this % interp % N + 1, &
!                           this % interp % N + 1, &
!                           this % interp % N + 1, &
!                           this % nVar, &
!                           nglobalelem/)

!       ! Offsets and dimensions for element boundary data
!       bOffset(1:6) = (/0,0,0,0,0,elemoffset/)
!       bGlobalDims(1:6) = (/3, &
!                            this % interp % N + 1, &
!                            this % interp % N + 1, &
!                            this % nVar, &
!                            6, &
!                            nglobalelem/)

!       CALL CreateGroup_HDF5(fileId,TRIM(group))

!       DO ivar = 1, this % nVar
!         CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
!       ENDDO

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
!                            this % interior,offset,globalDims)

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
!                            this % boundary,bOffset,bGlobalDims)

!   END SUBROUTINE WriteHDF5_MPI_Vector3D

!   SUBROUTINE WriteHDF5_Vector3D(this,fileId,group)
!     IMPLICIT NONE
!     CLASS(Vector3D), INTENT(in) :: this
!     INTEGER(HID_T), INTENT(in) :: fileId
!     CHARACTER(*), INTENT(in) :: group
!     ! Local
!     INTEGER :: ivar

!       CALL CreateGroup_HDF5(fileId,TRIM(group))

!       DO ivar = 1, this % nVar
!         CALL this % meta(ivar) % WriteHDF5( group, ivar, fileId )
!       ENDDO

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/interior", &
!                            this % interior)

!       CALL WriteArray_HDF5(fileId,TRIM(group)//"/boundary", &
!                            this % boundary)

!   END SUBROUTINE WriteHDF5_Vector3D

! ! -- Tensor2D -- !

  subroutine Init_Tensor2D(this,interp,nVar,nElem)
    implicit none
    class(Tensor2D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer :: N

    this % interp => interp
    this % nVar = nVar
    this % nElem = nElem
    N = interp % N

    call hipcheck(hipMallocManaged(this % interior,interp % N + 1,interp % N + 1,nelem,nvar,2,2,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % boundary,interp % N + 1,4,nelem,nvar,2,2,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % extBoundary,interp % N + 1,4,nelem,nvar,2,2,hipMemAttachGlobal))

    allocate (this % meta(1:nVar))
    allocate (this % eqn(1:4*nVar))

  end subroutine Init_Tensor2D

  subroutine Free_Tensor2D(this)
    implicit none
    class(Tensor2D),intent(inout) :: this

    this % interp => null()
    this % nVar = 0
    this % nElem = 0

    call hipcheck(hipFree(this % interior))
    call hipcheck(hipFree(this % boundary))
    call hipcheck(hipFree(this % extBoundary))

    deallocate (this % meta)
    deallocate (this % eqn)

  end subroutine Free_Tensor2D

  subroutine UpdateDevice_Tensor2D(this)
    implicit none
    class(Tensor2D),intent(inout) :: this

    call hipcheck(hipMemPrefetchAsync(c_loc(this % interior),sizeof(this % interior),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % boundary),sizeof(this % boundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % extBoundary),sizeof(this % extBoundary),0,c_null_ptr))

  end subroutine UpdateDevice_Tensor2D

  subroutine BoundaryInterp_Tensor2D_cpu(this)
    implicit none
    class(Tensor2D),intent(inout) :: this

    call this % interp % TensorBoundaryInterp_2D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem)

  end subroutine BoundaryInterp_Tensor2D_cpu

  subroutine BoundaryInterp_Tensor2D_gpu(this,handle)
    implicit none
    class(Tensor2D),intent(inout) :: this
    type(c_ptr),intent(in) :: handle

    call this % interp % TensorBoundaryInterp_2D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem, &
                                                 handle)

  end subroutine BoundaryInterp_Tensor2D_gpu

  subroutine Divergence_Tensor2D_cpu(this,that)
    implicit none
    class(Tensor2D),intent(in) :: this
    class(Vector2D),intent(inout) :: that

    call this % interp % TensorDivergence_2D(this % interior, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem)

  end subroutine Divergence_Tensor2D_cpu

  subroutine Divergence_Tensor2D_gpu(this,that,handle)
    implicit none
    class(Tensor2D),intent(in) :: this
    class(Vector2D),intent(inout) :: that
    type(c_ptr),intent(inout) :: handle

    call this % interp % TensorDivergence_2D(this % interior, &
                                             that % interior, &
                                             this % nVar, &
                                             this % nElem, &
                                             handle)

  end subroutine Divergence_Tensor2D_gpu

  subroutine Determinant_Tensor2D(this,that)
    implicit none
    class(Tensor2D),intent(in) :: this
    type(Scalar2D),intent(inout) :: that
    ! Local
    integer :: iEl,iVar,i,j
    do iVar = 1,this % nVar
      do iEl = 1,this % nElem
        do j = 1,this % interp % N + 1
          do i = 1,this % interp % N + 1

            that % interior(i,j,iEl,iVar) = this % interior(i,j,iEl,iVar,1,1)* &
                                               this % interior(i,j,iEl,iVar,2,2) - &
                                               this % interior(i,j,iEl,iVar,1,2)* &
                                               this % interior(i,j,iEl,iVar,2,1)

          end do
        end do
      end do
    end do

  end subroutine Determinant_Tensor2D

! ! -- Tensor3D -- !

  subroutine Init_Tensor3D(this,interp,nVar,nElem)
    implicit none
    class(Tensor3D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer :: N

    this % interp => interp
    this % nVar = nVar
    this % nElem = nElem
    N = interp % N

    call hipcheck(hipMallocManaged(this % interior,interp % N + 1,interp % N + 1,interp % N +1,nelem,nvar,3,3,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % boundary,interp % N + 1,interp % N + 1,6,nelem,nvar,3,3,hipMemAttachGlobal))
   call hipcheck(hipMallocManaged(this % extBoundary,interp % N + 1,interp % N + 1,6,nelem,nvar,3,3,hipMemAttachGlobal))

    allocate (this % meta(1:nVar))
    allocate (this % eqn(1:9*nVar))

  end subroutine Init_Tensor3D

  subroutine Free_Tensor3D(this)
    implicit none
    class(Tensor3D),intent(inout) :: this

    this % interp => null()
    this % nVar = 0
    this % nElem = 0

    call hipcheck(hipFree(this % interior))
    call hipcheck(hipFree(this % boundary))
    call hipcheck(hipFree(this % extBoundary))

    deallocate (this % meta)
    deallocate (this % eqn)

  end subroutine Free_Tensor3D

  subroutine UpdateDevice_Tensor3D(this)
    implicit none
    class(Tensor3D),intent(inout) :: this

    call hipcheck(hipMemPrefetchAsync(c_loc(this % interior),sizeof(this % interior),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % boundary),sizeof(this % boundary),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % extBoundary),sizeof(this % extBoundary),0,c_null_ptr))

  end subroutine UpdateDevice_Tensor3D

  subroutine BoundaryInterp_Tensor3D_cpu(this)
    implicit none
    class(Tensor3D),intent(inout) :: this

    call this % interp % TensorBoundaryInterp_3D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem)

  end subroutine BoundaryInterp_Tensor3D_cpu

  subroutine BoundaryInterp_Tensor3D_gpu(this,handle)
    implicit none
    class(Tensor3D),intent(inout) :: this
    type(c_ptr),intent(in) :: handle

    call this % interp % TensorBoundaryInterp_3D(this % interior, &
                                                 this % boundary, &
                                                 this % nVar, &
                                                 this % nElem, &
                                                 handle)

  end subroutine BoundaryInterp_Tensor3D_gpu

  subroutine Determinant_Tensor3D(this,that)
    implicit none
    class(Tensor3D),intent(in) :: this
    type(Scalar3D),intent(inout) :: that
    ! Local
    integer :: iEl,iVar,i,j,k

    do iEl = 1,this % nElem
      do iVar = 1,this % nVar
        do k = 1,this % interp % N+1
          do j = 1,this % interp % N+1
            do i = 1,this % interp % N+1

              that % interior(i,j,k,iEl,iVar) = &
                this % interior(i,j,k,iEl,iVar,1,1)* &
                (this % interior(i,j,k,iEl,iVar,2,2)* &
                 this % interior(i,j,k,iEl,iVar,3,3) - &
                 this % interior(i,j,k,iEl,iVar,2,3)* &
                 this % interior(i,j,k,iEl,iVar,3,2)) - &
                this % interior(i,j,k,iEl,iVar,2,1)* &
                (this % interior(i,j,k,iEl,iVar,1,2)* &
                 this % interior(i,j,k,iEl,iVar,3,3) - &
                 this % interior(i,j,k,iEl,iVar,1,3)* &
                 this % interior(i,j,k,iEl,iVar,3,2)) + &
                this % interior(i,j,k,iEl,iVar,3,1)* &
                (this % interior(i,j,k,iEl,iVar,1,2)* &
                 this % interior(i,j,k,iEl,iVar,2,3) - &
                 this % interior(i,j,k,iEl,iVar,1,3)* &
                 this % interior(i,j,k,iEl,iVar,2,2))

            end do
          end do
        end do
      end do
    end do

  end subroutine Determinant_Tensor3D

end module SELF_Data
