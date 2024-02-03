! SELF_Lagrange.f90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Lagrange

  use iso_fortran_env
  use iso_c_binding
  use hipfort
  use hipfort_check
  use hipfort_hipmalloc
  use hipfort_hipblas

  use SELF_Constants
  use SELF_SupportRoutines
  use SELF_Quadrature
  !use SELF_HDF5
  !use HDF5

  use hipfort_hipblas

  use iso_c_binding
  implicit none

  type,public :: Lagrange
    !! A data structure for working with Lagrange Interpolating Polynomials in one, two, and three dimensions.
    !! The Lagrange data-structure stores the information necessary to interpolate between two
    !! sets of grid-points and to estimate the derivative of data at native grid points. Routines for
    !! multidimensional interpolation are based on the tensor product of 1-D interpolants. It is
    !! assumed that the polynomial degree (and the interpolation nodes) are the same in each direction.
    !! This assumption permits the storage of only one array of interpolation nodes and barycentric
    !! weights and is what allows this data structure to be flexible.

    integer :: N
      !! The number of control points.

    integer :: controlNodeType

    integer :: M
      !! The number of target points.

    integer :: targetNodeType

    real(prec),pointer,dimension(:) :: controlPoints
      !! The set of nodes in one dimension where data is known.
      !! To create higher dimension interpolation and differentiation operators, structured grids in two and three
      !! dimensions are created by tensor products of the controlPoints. This design decision implies that all
      !! spectral element methods supported by the Lagrange class have the same polynomial degree in each
      !! computational/spatial dimension. In practice, the controlPoints are the Legendre-Gauss, Legendre-Gauss-Lobatto,
      !! Legendre-Gauss-Radau, Chebyshev-Gauss, Chebyshev-Gauss-Lobatto, or Chebyshev-Gauss-Radau quadrature points over
      !! the domain [-1,1] (computational space). The Init routine for this class restricts controlPoints to one of
      !! these quadrature types or uniform points on [-1,1].

    real(prec),pointer,dimension(:) :: targetPoints
      !! The set of nodes in one dimension where data is to be interpolated to. To create higher dimension interpolation
      !! and differentiation operators, structured grids in two and three dimensions are created by tensor products of
      !! the targetPoints. In practice, the targetPoints are set to a uniformly distributed set of points between [-1,1]
      !! (computational space) to allow for interpolation from unevenly spaced quadrature points to a plotting grid.

    real(prec),pointer,dimension(:) :: bWeights
      !! The barycentric weights that are calculated from the controlPoints and used for interpolation.

    real(prec),pointer,dimension(:) :: qWeights
      !! The quadrature weights for discrete integration. The quadradture weights depend on the type of controlPoints
      !! provided; one of Legendre-Gauss, Legendre-Gauss-Lobatto, Legendre-Gauss-Radau, Chebyshev-Gauss,
      !! Chebyshev-Gauss-Lobatto, Chebyshev-Gauss Radau, or Uniform. If Uniform, the quadrature weights are constant
      !! $$dx = \frac{2.0}{N+1}$$.

    real(prec),pointer,dimension(:,:) :: iMatrix
      !! The interpolation matrix (transpose) for mapping data from the control grid to the target grid.

    real(prec),pointer,dimension(:,:) :: dMatrix
      !! The derivative matrix for mapping function nodal values to a nodal values of the derivative estimate. The
      !! dMatrix is based on a strong form of the derivative.

    real(prec),pointer,dimension(:,:) :: dgMatrix
      !! The derivative matrix for mapping function nodal values to a nodal values of the derivative estimate. The dgMatrix is based
      !! on a weak form of the derivative. It must be used with bMatrix to account for boundary contributions in the weak form.

    real(prec),pointer,dimension(:,:) :: bMatrix
      !! The boundary interpolation matrix that is used to map a grid of nodal values at the control points to the element boundaries.

  contains

    procedure,public :: Init => Init_Lagrange
    procedure,public :: Free => Free_Lagrange

    procedure,public :: UpdateDevice => UpdateDevice_Lagrange

    generic,public :: ScalarGridInterp_1D => ScalarGridInterp_1D_cpu,ScalarGridInterp_1D_gpu
    procedure,private :: ScalarGridInterp_1D_cpu,ScalarGridInterp_1D_gpu

    ! GENERIC,PUBLIC :: ScalarGridInterp_2D => ScalarGridInterp_2D_cpu,ScalarGridInterp_2D_gpu
    ! PROCEDURE,PRIVATE :: ScalarGridInterp_2D_cpu,ScalarGridInterp_2D_gpu

    ! GENERIC,PUBLIC :: VectorGridInterp_2D => VectorGridInterp_2D_cpu,VectorGridInterp_2D_gpu
    ! PROCEDURE,PRIVATE :: VectorGridInterp_2D_cpu,VectorGridInterp_2D_gpu

    ! GENERIC,PUBLIC :: ScalarGridInterp_3D => ScalarGridInterp_3D_cpu,ScalarGridInterp_3D_gpu
    ! PROCEDURE,PRIVATE :: ScalarGridInterp_3D_cpu,ScalarGridInterp_3D_gpu

    ! GENERIC,PUBLIC :: VectorGridInterp_3D => VectorGridInterp_3D_cpu,VectorGridInterp_3D_gpu
    ! PROCEDURE,PRIVATE :: VectorGridInterp_3D_cpu,VectorGridInterp_3D_gpu

    ! generic,public :: ScalarBoundaryInterp_1D => ScalarBoundaryInterp_1D_cpu,ScalarBoundaryInterp_1D_gpu
    ! procedure,private :: ScalarBoundaryInterp_1D_cpu,ScalarBoundaryInterp_1D_gpu

    ! GENERIC,PUBLIC :: ScalarBoundaryInterp_2D => ScalarBoundaryInterp_2D_cpu,ScalarBoundaryInterp_2D_gpu
    ! PROCEDURE,PRIVATE :: ScalarBoundaryInterp_2D_cpu,ScalarBoundaryInterp_2D_gpu

    ! GENERIC,PUBLIC :: VectorBoundaryInterp_2D => VectorBoundaryInterp_2D_cpu,VectorBoundaryInterp_2D_gpu
    ! PROCEDURE,PRIVATE :: VectorBoundaryInterp_2D_cpu,VectorBoundaryInterp_2D_gpu

    ! GENERIC,PUBLIC :: TensorBoundaryInterp_2D => TensorBoundaryInterp_2D_cpu,TensorBoundaryInterp_2D_gpu
    ! PROCEDURE,PRIVATE :: TensorBoundaryInterp_2D_cpu,TensorBoundaryInterp_2D_gpu

    ! GENERIC,PUBLIC :: ScalarBoundaryInterp_3D => ScalarBoundaryInterp_3D_cpu,ScalarBoundaryInterp_3D_gpu
    ! PROCEDURE,PRIVATE :: ScalarBoundaryInterp_3D_cpu,ScalarBoundaryInterp_3D_gpu

    ! GENERIC,PUBLIC :: VectorBoundaryInterp_3D => VectorBoundaryInterp_3D_cpu,VectorBoundaryInterp_3D_gpu
    ! PROCEDURE,PRIVATE :: VectorBoundaryInterp_3D_cpu,VectorBoundaryInterp_3D_gpu

    ! GENERIC,PUBLIC :: TensorBoundaryInterp_3D => TensorBoundaryInterp_3D_cpu,TensorBoundaryInterp_3D_gpu
    ! PROCEDURE,PRIVATE :: TensorBoundaryInterp_3D_cpu,TensorBoundaryInterp_3D_gpu

    generic,public :: Derivative_1D => Derivative_1D_cpu,Derivative_1D_gpu
    procedure,private :: Derivative_1D_cpu,Derivative_1D_gpu

    ! generic,public :: DGDerivative_1D => DGDerivative_1D_cpu,DGDerivative_1D_gpu
    ! procedure,private :: DGDerivative_1D_cpu,DGDerivative_1D_gpu

    ! GENERIC,PUBLIC :: ScalarGradient_2D => ScalarGradient_2D_cpu,ScalarGradient_2D_gpu
    ! PROCEDURE,PRIVATE :: ScalarGradient_2D_cpu,ScalarGradient_2D_gpu

    ! GENERIC,PUBLIC :: VectorGradient_2D => VectorGradient_2D_cpu,VectorGradient_2D_gpu
    ! PROCEDURE,PRIVATE :: VectorGradient_2D_cpu,VectorGradient_2D_gpu

    ! GENERIC,PUBLIC :: VectorDivergence_2D => VectorDivergence_2D_cpu,VectorDivergence_2D_gpu
    ! PROCEDURE,PRIVATE :: VectorDivergence_2D_cpu,VectorDivergence_2D_gpu

    ! GENERIC,PUBLIC :: VectorDGDivergence_2D => VectorDGDivergence_2D_cpu,VectorDGDivergence_2D_gpu
    ! PROCEDURE,PRIVATE :: VectorDGDivergence_2D_cpu,VectorDGDivergence_2D_gpu

    ! GENERIC,PUBLIC :: ScalarGradient_3D => ScalarGradient_3D_cpu,ScalarGradient_3D_gpu
    ! PROCEDURE,PRIVATE :: ScalarGradient_3D_cpu,ScalarGradient_3D_gpu

    ! GENERIC,PUBLIC :: VectorGradient_3D => VectorGradient_3D_cpu,VectorGradient_3D_gpu
    ! PROCEDURE,PRIVATE :: VectorGradient_3D_cpu,VectorGradient_3D_gpu

    ! GENERIC,PUBLIC :: VectorDivergence_3D => VectorDivergence_3D_cpu,VectorDivergence_3D_gpu
    ! PROCEDURE,PRIVATE :: VectorDivergence_3D_cpu,VectorDivergence_3D_gpu

    ! GENERIC,PUBLIC :: VectorDGDivergence_3D => VectorDGDivergence_3D_cpu,VectorDGDivergence_3D_gpu
    ! PROCEDURE,PRIVATE :: VectorDGDivergence_3D_cpu,VectorDGDivergence_3D_gpu

    !procedure,public :: WriteHDF5 => WriteHDF5_Lagrange
    procedure,private :: CalculateBarycentricWeights
    procedure,private :: CalculateInterpolationMatrix
    procedure,private :: CalculateDerivativeMatrix
    procedure,private :: CalculateLagrangePolynomials

  end type Lagrange


  ! /////////////// !
  ! Boundary Interpolation Routines

  ! interface
  !   subroutine ScalarBoundaryInterp_1D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
  !     bind(c,name="ScalarBoundaryInterp_1D_gpu_wrapper")
  !     use iso_c_binding
  !     implicit none
  !     type(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
  !     integer(c_int),value :: N,nVar,nEl
  !   end subroutine ScalarBoundaryInterp_1D_gpu_wrapper
  ! end interface

  ! INTERFACE
  !   SUBROUTINE ScalarBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
  !     bind(c,name="ScalarBoundaryInterp_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE ScalarBoundaryInterp_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
  !     bind(c,name="VectorBoundaryInterp_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorBoundaryInterp_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE TensorBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
  !     bind(c,name="TensorBoundaryInterp_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE TensorBoundaryInterp_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE ScalarBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
  !     bind(c,name="ScalarBoundaryInterp_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE ScalarBoundaryInterp_3D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
  !     bind(c,name="VectorBoundaryInterp_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorBoundaryInterp_3D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE TensorBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
  !     bind(c,name="TensorBoundaryInterp_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE TensorBoundaryInterp_3D_gpu_wrapper
  ! END INTERFACE

  ! /////////////// !

  ! interface
  !   subroutine Derivative_1D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="Derivative_1D_gpu_wrapper")
  !     use iso_c_binding
  !     implicit none
  !     type(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     integer(c_int),value :: N,nVar,nEl
  !   end subroutine Derivative_1D_gpu_wrapper
  ! end interface

  ! interface
  !   subroutine DGDerivative_1D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="DGDerivative_1D_gpu_wrapper")
  !     use iso_c_binding
  !     implicit none
  !     type(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
  !     integer(c_int),value :: N,nVar,nEl
  !   end subroutine DGDerivative_1D_gpu_wrapper
  ! end interface

  ! INTERFACE
  !   SUBROUTINE ScalarGradient_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="ScalarGradient_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE ScalarGradient_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorGradient_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorGradient_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorGradient_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorDGGradient_2D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorDGGradient_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorDGGradient_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorDivergence_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorDivergence_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorDivergence_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !  SUBROUTINE VectorDGDivergence_2D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorDGDivergence_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorDGDivergence_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE ScalarGradient_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="ScalarGradient_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE ScalarGradient_3D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorGradient_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorGradient_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorGradient_3D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorDivergence_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorDivergence_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorDivergence_3D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !  SUBROUTINE VectorDGDivergence_3D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorDGDivergence_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorDGDivergence_3D_gpu_wrapper
  ! END INTERFACE

contains

  subroutine Init_Lagrange(this,N,controlNodeType,M,targetNodeType)
    !! Initialize an instance of the Lagrange class
    !! On output, all of the attributes for the Lagrange class are allocated and values are initialized according to the number of
    !! control points, number of target points, and the types for the control and target nodes.
    !! If a GPU is available, device pointers for the Lagrange attributes are allocated and initialized.
    implicit none
    class(Lagrange),intent(out) :: this
    !! Lagrange class instance
    integer,intent(in)          :: N
    !! The number of control points for interpolant
    integer,intent(in)          :: M
    !! The number of target points for the interpolant
    integer,intent(in)          :: controlNodeType
    !! The integer code specifying the type of control points. Parameters are defined in SELF_Constants.f90. One of GAUSS(=1),
    !! GAUSS_LOBATTO(=2), or UNIFORM(=3)
    integer,intent(in)          :: targetNodeType
    !! The integer code specifying the type of target points. Parameters are defined in SELF_Constants.f90. One of GAUSS(=1),
    !! GAUSS_LOBATTO(=2), or UNIFORM(=3)
    ! -------!
    ! Local
    real(prec) :: q(0:M)

    this % N = N
    this % M = M
    this % controlNodeType = controlNodeType
    this % targetNodeType = targetNodeType

    call hipcheck(hipMallocManaged(this % controlPoints,N + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % targetPoints,M + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % bWeights,N + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % qWeights,N + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % iMatrix,N + 1,M + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % dMatrix,N + 1,N + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % dgMatrix,N + 1,N + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % bMatrix,N + 1,2,hipMemAttachGlobal))

    if (controlNodeType == GAUSS .or. controlNodeType == GAUSS_LOBATTO) then

      call LegendreQuadrature(N, &
                              this % controlPoints, &
                              this % qWeights, &
                              controlNodeType)

    elseif (controlNodeType == CHEBYSHEV_GAUSS .or. controlNodeType == CHEBYSHEV_GAUSS_LOBATTO) then

      call ChebyshevQuadrature(N, &
                               this % controlPoints, &
                               this % qWeights, &
                               controlNodeType)

    elseif (controlNodeType == UNIFORM) then

      this % controlPoints = UniformPoints(-1.0_prec,1.0_prec,0,N)
      this % qWeights = 2.0_prec/real(N,prec)

    end if

    ! Target Points
    if (targetNodeType == GAUSS .or. targetNodeType == GAUSS_LOBATTO) then

      call LegendreQuadrature(M, &
                              this % targetPoints, &
                              q, &
                              targetNodeType)

    elseif (targetNodeType == UNIFORM) then

      this % targetPoints  = UniformPoints(-1.0_prec,1.0_prec,0,M)

    end if

    call this % CalculateBarycentricWeights()
    call this % CalculateInterpolationMatrix()
    call this % CalculateDerivativeMatrix()
    this % bMatrix (1:N + 1,1) = this % CalculateLagrangePolynomials(-1.0_prec)
    this % bMatrix (1:N + 1,2) = this % CalculateLagrangePolynomials(1.0_prec)

    call this % UpdateDevice()

  end subroutine Init_Lagrange

  subroutine Free_Lagrange(this)
    !! Frees all memory (host and device) associated with an instance of the Lagrange class
    implicit none
    class(Lagrange),intent(inout) :: this
    !! Lagrange class instance

    call hipcheck(hipFree(this % controlPoints))
    call hipcheck(hipFree(this % targetPoints))
    call hipcheck(hipFree(this % bWeights))
    call hipcheck(hipFree(this % qWeights))
    call hipcheck(hipFree(this % iMatrix))
    call hipcheck(hipFree(this % dMatrix))
    call hipcheck(hipFree(this % dgMatrix))
    call hipcheck(hipFree(this % bMatrix))

  end subroutine Free_Lagrange

  subroutine UpdateDevice_Lagrange(this)
    !! Copy the Lagrange attributes from the host (CPU) to the device (GPU)
    implicit none
    class(Lagrange),intent(inout) :: this
    !! Lagrange class instance

    call hipcheck(hipMemPrefetchAsync(c_loc(this % controlPoints),sizeof(this % controlPoints),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % targetPoints),sizeof(this % targetPoints),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % bWeights),sizeof(this % bWeights),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % qWeights),sizeof(this % qWeights),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % iMatrix),sizeof(this % iMatrix),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % dMatrix),sizeof(this % dMatrix),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % dgMatrix),sizeof(this % dgMatrix),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % bMatrix),sizeof(this % bMatrix),0,c_null_ptr))

  end subroutine UpdateDevice_Lagrange

  subroutine ScalarGridInterp_1D_cpu(this,f,fInterp,nVariables,nElements)
    !! Host (CPU) implementation of the ScalarGridInterp_1D interface.
    !! In most cases, you should use the `ScalarGridInterp_1D` generic interface,
    !! rather than calling this routine directly.
    !! Interpolate a scalar-1D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{m,ivar,iel} = \sum_{i=0}^N f_{i,ivar,iel} I_{i,m} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in)     :: nVariables
    !! The number of variables/functions that are interpolated
    integer,intent(in)     :: nElements
    !! The number of spectral elements in the SEM grid
    real(prec), intent(in)  :: f(1:this % N+1,1:nVariables,1:nElements)
    !! (Input) Array of function values, defined on the control grid
    real(prec), intent(out) :: fInterp(1:this % M+1,1:nVariables,1:nElements)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: iVar,iEl,i,ii
    real(prec) :: floc

    do iEl = 1,nElements
      do iVar = 1,nVariables
        do i = 1,this % M+1
          floc = 0.0_prec
          do ii = 1,this % N+1
            floc = floc + this % iMatrix (ii,i)*f(ii,iVar,iEl)
          end do
          fInterp(i,iVar,iEl) = floc
        end do
      end do
    end do

  end subroutine ScalarGridInterp_1D_cpu

  subroutine ScalarGridInterp_1D_gpu(this,f,fInterp,nVariables,nElements,hipblas_handle)
    !! Device (GPU) implementation of the ScalarGridInterp_1D interface.
    !! In most cases, you should use the `ScalarGridInterp_1D` generic interface,
    !! rather than calling this routine directly.
    !! This routine calls hip/SELF_Lagrange.cpp:ScalarGridInterp_1D_gpu_wrapper
    !! Interpolate a scalar-1D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{m,ivar,iel} = \sum_{i=0}^N f_{i,ivar,iel} I_{i,m} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in) :: nVariables
    !! The number of variables/functions that are interpolated
    integer,intent(in) :: nElements
    !! The number of spectral elements in the SEM grid
    real(prec), pointer, intent(in)  :: f(:,:,:)
    !! (Input) Array of function values, defined on the control grid
    real(prec), pointer, intent(inout) :: fInterp(:,:,:)
    !! (Output) Array of function values, defined on the target grid
    type(c_ptr), intent(inout) :: hipblas_handle
    ! Local
    integer(c_int) :: m
    integer(c_int) :: n
    integer(c_int) :: k
    real(c_prec) :: alpha
    integer(c_int) :: lda
    integer(c_int) :: ldb
    integer(c_int) :: ldc
    real(c_prec) :: beta
 
    m = this % M+1 ! number of rows of A^T
    n = nvariables*nelements ! number of columns of B
    k = this % N+1! number of columns of A^T
    alpha = 1.0_c_prec
    lda = k ! leading dimension of A (interoplation matrix)
    ldb = k ! leading dimension of B (f)
    ldc = m ! leading dimension of C (fTarget)
    beta = 0.0_c_prec
 
#ifdef DOUBLE_PRECISION
    call hipblasCheck(hipblasDgemm(hipblas_handle,&
        HIPBLAS_OP_T, HIPBLAS_OP_N, &
        m, n, k, alpha, &
        c_loc(this % iMatrix), lda, &
        c_loc(f), ldb, &
        beta, &
        c_loc(fInterp), ldc))
#else
    call hipblasCheck(hipblasSgemm(hipblas_handle,&
        HIPBLAS_OP_T, HIPBLAS_OP_N, &
        m, n, k, alpha, &
        c_loc(this % iMatrix), lda, &
        c_loc(f), ldb, &
        beta, &
        c_loc(fInterp), ldc))
#endif

  end subroutine ScalarGridInterp_1D_gpu

!   subroutine ScalarGridInterp_2D_cpu(this,f,fNew,nVariables,nElements)
!     !! Host (CPU) implementation of the ScalarGridInterp_2D interface.
!     !! In most cases, you should use the `ScalarGridInterp_2D` generic interface,
!     !! rather than calling this routine directly.
!     !! Interpolate a scalar-2D (real) array from the control grid to the target grid.
!     !! The control and target grids are the ones associated with an initialized
!     !! Lagrange instance.
!     !!
!     !! Interpolation is applied using a series of matrix-vector multiplications, using
!     !! the Lagrange class's interpolation matrix
!     !!
!     !! $$ \tilde{f}_{m,n,ivar,iel} = \sum_{j=0}^N \sum_{i=0}^N f_{i,j,ivar,iel} I_{i,m} I_{j,n} $$
!     !!
!     implicit none
!     class(Lagrange),intent(in) :: this
!     !! Lagrange class instance
!     integer,intent(in)     :: nVariables
!     !! The number of variables/functions that are interpolated
!     integer,intent(in)     :: nElements
!     !! The number of spectral elements in the SEM grid
!     real(prec),intent(in)  :: f(0:this % N,0:this % N,1:nVariables,1:nElements)
!     !! (Input) Array of function values, defined on the control grid
!     real(prec),intent(out) :: fNew(0:this % M,0:this % M,1:nVariables,1:nElements)
!     !! (Output) Array of function values, defined on the target grid
!     ! Local
!     integer :: i,j,ii,jj,iEl,iVar
!     real(prec) :: fi,fij

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do j = 0,this % M
!           do i = 0,this % M

!             fij = 0.0_prec
!             do jj = 0,this % N
!               fi = 0.0_prec
!               do ii = 0,this % N
!                 fi = fi + f(ii,jj,iVar,iEl)*this % iMatrix (ii,i)
!               end do
!               fij = fij + fi*this % iMatrix (jj,j)
!             end do
!             fNew(i,j,iVar,iEl) = fij

!           end do
!         end do
!       end do
!     end do

!   end subroutine ScalarGridInterp_2D_cpu

!   subroutine ScalarGridInterp_2D_gpu(this,f_dev,fInterp_dev,nVariables,nElements)
!     !! Device (GPU) implementation of the ScalarGridInterp_2D interface.
!     !! In most cases, you should use the `ScalarGridInterp_2D` generic interface,
!     !! rather than calling this routine directly.
!     !! This routine calls hip/SELF_Lagrange.cpp:ScalarGridInterp_2D_gpu_wrapper
!     !! Interpolate a scalar-2D (real) array from the control grid to the target grid.
!     !! The control and target grids are the ones associated with an initialized
!     !! Lagrange instance.
!     !!
!     !! Interpolation is applied using a series of matrix-vector multiplications, using
!     !! the Lagrange class's interpolation matrix
!     !!
!     !! $$ \tilde{f}_{m,n,ivar,iel} = \sum_{j=0}^N \sum_{i=0}^N f_{i,j,ivar,iel} I_{i,m} I_{j,n} $$
!     !!
!     implicit none
!     class(Lagrange),intent(in) :: this
!     !! Lagrange class instance
!     integer,intent(in) :: nVariables
!     !! The number of variables/functions that are interpolated
!     integer,intent(in) :: nElements
!     !! The number of spectral elements in the SEM grid
!     type(c_ptr),intent(in)  :: f_dev
!     !! (Input) Array of function values, defined on the control grid
!     type(c_ptr),intent(out) :: fInterp_dev
!     !! (Output) Array of function values, defined on the target grid

!     call ScalarGridInterp_2D_gpu_wrapper(this % iMatrix % deviceData, &
!                                          f_dev,fInterp_dev, &
!                                          this % N,this % M, &
!                                          nVariables,nElements)

!   end subroutine ScalarGridInterp_2D_gpu

!   subroutine VectorGridInterp_2D_cpu(this,f,fNew,nVariables,nElements)
!     !! Host (CPU) implementation of the VectorGridInterp_2D interface.
!     !! In most cases, you should use the `VectorGridInterp_2D` generic interface,
!     !! rather than calling this routine directly.
!     !! Interpolate a vector-2D (real) array from the control grid to the target grid.
!     !! The control and target grids are the ones associated with an initialized
!     !! Lagrange instance.
!     !!
!     !! Interpolation is applied using a series of matrix-vector multiplications, using
!     !! the Lagrange class's interpolation matrix
!     !!
!     !! $$ \tilde{f}_{dir,m,n,ivar,iel} = \sum_{j=0}^N \sum_{i=0}^N f_{dir,i,j,ivar,iel} I_{i,m} I_{j,n} $$
!     !!
!     implicit none
!     class(Lagrange),intent(in) :: this
!     !! Lagrange class instance
!     integer,intent(in)     :: nVariables
!     !! The number of variables/functions that are interpolated
!     integer,intent(in)     :: nElements
!     !! The number of spectral elements in the SEM grid
!     real(prec),intent(in)  :: f(1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!     !! (Input) Array of function values, defined on the control grid
!     real(prec),intent(out) :: fNew(1:2,0:this % M,0:this % M,1:nVariables,1:nElements)
!     !! (Output) Array of function values, defined on the target grid
!     ! Local
!     integer :: i,j,ii,jj,iEl,iVar
!     real(prec) :: fi(1:2)

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do j = 0,this % M
!           do i = 0,this % M

!             fNew(1,i,j,iVar,iEl) = 0.0_prec
!             fNew(2,i,j,iVar,iEl) = 0.0_prec

!             do jj = 0,this % N

!               fi(1:2) = 0.0_prec
!               do ii = 0,this % N
!                 fi(1:2) = fi(1:2) + f(1:2,ii,jj,iVar,iEl)*this % iMatrix (ii,i)
!               end do

!               fNew(1:2,i,j,iVar,iEl) = fNew(1:2,i,j,iVar,iEl) + fi(1:2)*this % iMatrix (jj,j)

!             end do

!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorGridInterp_2D_cpu
! !
!   subroutine VectorGridInterp_2D_gpu(this,f_dev,fInterp_dev,nVariables,nElements)
!     !! Device (GPU) implementation of the VectorGridInterp_2D interface.
!     !! In most cases, you should use the `VectorGridInterp_2D` generic interface,
!     !! rather than calling this routine directly.
!     !! This routine calls hip/SELF_Lagrange.cpp:VectorGridInterp_2D_gpu_wrapper
!     !! Interpolate a vector-2D (real) array from the control grid to the target grid.
!     !! The control and target grids are the ones associated with an initialized
!     !! Lagrange instance.
!     !!
!     !! Interpolation is applied using a series of matrix-vector multiplications, using
!     !! the Lagrange class's interpolation matrix
!     !!
!     !! $$ \tilde{f}_{dir,m,n,ivar,iel} = \sum_{j=0}^N \sum_{i=0}^N f_{dir,i,j,ivar,iel} I_{i,m} I_{j,n} $$
!     !!
!     implicit none
!     class(Lagrange),intent(in) :: this
!     !! Lagrange class instance
!     integer,intent(in) :: nVariables
!     !! The number of variables/functions that are interpolated
!     integer,intent(in) :: nElements
!     !! The number of spectral elements in the SEM grid
!     type(c_ptr),intent(in)  :: f_dev
!     !! (Input) Array of function values, defined on the control grid
!     type(c_ptr),intent(out) :: fInterp_dev
!     !! (Output) Array of function values, defined on the target grid

!     call VectorGridInterp_2D_gpu_wrapper(this % iMatrix % deviceData, &
!                                          f_dev,fInterp_dev, &
!                                          this % N,this % M, &
!                                          nVariables,nElements)

!   end subroutine VectorGridInterp_2D_gpu

!   subroutine ScalarGridInterp_3D_cpu(this,f,fInterp,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nVariables,nElements
!     real(prec),intent(in)  :: f(0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out) :: fInterp(0:this % M,0:this % M,0:this % M,1:nVariables,1:nElements)
!     ! Local
!     integer :: iEl,iVar,i,j,k,ii,jj,kk
!     real(prec) :: fi,fij,fijk

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do k = 0,this % M
!           do j = 0,this % M
!             do i = 0,this % M

!               fijk = 0.0_prec
!               do kk = 0,this % N
!                 fij = 0.0_prec
!                 do jj = 0,this % N
!                   fi = 0.0_prec
!                   do ii = 0,this % N
!                     fi = fi + f(ii,jj,kk,iVar,iEl)*this % iMatrix (ii,i)
!                   end do
!                   fij = fij + fi*this % iMatrix (jj,j)
!                 end do
!                 fijk = fijk + fij*this % iMatrix (kk,k)
!               end do
!               fInterp(i,j,k,iVar,iEl) = fijk

!             end do
!           end do
!         end do
!       end do
!     end do

!   end subroutine ScalarGridInterp_3D_cpu
! !
!   subroutine ScalarGridInterp_3D_gpu(this,f_dev,fInterp_dev,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in) :: nVariables,nElements
!     type(c_ptr),intent(in)  :: f_dev
!     type(c_ptr),intent(out) :: fInterp_dev

!     call ScalarGridInterp_3D_gpu_wrapper(this % iMatrix % deviceData, &
!                                          f_dev,fInterp_dev, &
!                                          this % N,this % M, &
!                                          nVariables,nElements)

!   end subroutine ScalarGridInterp_3D_gpu

!   subroutine VectorGridInterp_3D_cpu(this,f,fInterp,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nVariables,nElements
!     real(prec),intent(in)  :: f(1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out) :: fInterp(1:3,0:this % M,0:this % M,0:this % M,1:nVariables,1:nElements)
!     ! Local
!     integer :: iEl,iVar,i,j,k,ii,jj,kk
!     real(prec) :: fi(1:3),fij(1:3)

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do k = 0,this % M
!           do j = 0,this % M
!             do i = 0,this % M

!               fInterp(1:3,i,j,k,iVar,iEl) = 0.0_prec
!               do kk = 0,this % N
!                 fij(1:3) = 0.0_prec
!                 do jj = 0,this % N
!                   fi(1:3) = 0.0_prec
!                   do ii = 0,this % N
!                     fi(1:3) = fi(1:3) + f(1:3,ii,jj,kk,iVar,iEl)*this % iMatrix (ii,i)
!                   end do
!                   fij(1:3) = fij(1:3) + fi(1:3)*this % iMatrix (jj,j)
!                 end do
!                 fInterp(1:3,i,j,k,iVar,iEl) = fInterp(1:3,i,j,k,iVar,iEl) + fij(1:3)*this % iMatrix (kk,k)
!               end do

!             end do
!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorGridInterp_3D_cpu
! !
!   subroutine VectorGridInterp_3D_gpu(this,f_dev,fInterp_dev,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in) :: nVariables,nElements
!     type(c_ptr),intent(in)  :: f_dev
!     type(c_ptr),intent(out) :: fInterp_dev

!     call VectorGridInterp_3D_gpu_wrapper(this % iMatrix % deviceData, &
!                                          f_dev,fInterp_dev, &
!                                          this % N,this % M, &
!                                          nVariables,nElements)

!   end subroutine VectorGridInterp_3D_gpu

! !   SUBROUTINE TensorGridInterp_3D_cpu(this,f,fInterp,nVariables,nElements)
! !     IMPLICIT NONE
! !     CLASS(Lagrange),INTENT(in) :: this
! !     INTEGER,INTENT(in)     :: nVariables,nElements
! !     REAL(prec),INTENT(in)  :: f(1:3,1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
! !     REAL(prec),INTENT(out) :: fInterp(1:3,1:3,0:this % M,0:this % M,0:this % M,1:nVariables,1:nElements)
! !     ! Local
! !     INTEGER :: iEl,iVar,i,j,k,ii,jj,kk
! !     REAL(prec) :: fi(1:3,1:3),fij(1:3,1:3)

! !     DO iEl = 1,nElements
! !       DO iVar = 1,nVariables
! !         DO k = 0,this % M
! !           DO j = 0,this % M
! !             DO i = 0,this % M

! !               fInterp(1:3,1:3,i,j,k,iVar,iEl) = 0.0_prec
! !               DO kk = 0,this % N
! !                 fij(1:3,1:3) = 0.0_prec
! !                 DO jj = 0,this % N
! !                   fi(1:3,1:3) = 0.0_prec
! !                   DO ii = 0,this % N
! !                     fi(1:3,1:3) = fi(1:3,1:3) + f(1:3,1:3,ii,jj,kk,iVar,iEl)*this % iMatrix (ii,i)
! !                   END DO
! !                   fij(1:3,1:3) = fij(1:3,1:3) + fi(1:3,1:3)*this % iMatrix (jj,j)
! !                 END DO
! !                 fInterp(1:3,1:3,i,j,k,iVar,iEl) = fInterp(1:3,1:3,i,j,k,iVar,iEl) + &
! !                                                   fij(1:3,1:3)*this % iMatrix (kk,k)
! !               END DO

! !             END DO
! !           END DO
! !         END DO
! !       END DO
! !     END DO

! !   END SUBROUTINE TensorGridInterp_3D_cpu
! ! !
! !   SUBROUTINE TensorGridInterp_3D_gpu(this,f_dev,fInterp_dev,nVariables,nElements)
! !     IMPLICIT NONE
! !     CLASS(Lagrange),INTENT(in) :: this
! !     INTEGER,INTENT(in) :: nVariables,nElements
! !     TYPE(c_ptr),INTENT(in)  :: f_dev
! !     TYPE(c_ptr),INTENT(out) :: fInterp_dev

! !     CALL TensorGridInterp_3D_gpu_wrapper(this % iMatrix % deviceData, &
! !                                          f_dev,fInterp_dev, &
! !                                          this % N,this % M, &
! !                                          nVariables,nElements)

! !   END SUBROUTINE TensorGridInterp_3D_gpu
! ================================================================================================ !
!
! Derivative_1D
!
!   Calculates the derivative of the Lagrange interpolant given a set of nodal function values at
!   the native interpolation nodes
!
!   Given a set of nodal values at the interpolation nodes, the derivative of a function through
!   the interpolation nodes can be estimated by
!
!                       f'_a = \sum_{i=0}^N f_{i} l'_i(\xi_a),   a=0,1,2,...,N
!
!   where l_i(\xi) are the Lagrange interpolating polynomials through the interpolation points.
!   The derivative matrix is D_{a,i} = l'_i(\xi_a) maps an array of nodal values at the interpolation
!   nodes to its estimated derivative. This routine serves as a wrapper to call either the CUDA
!   kernel (if CUDA is enabled) or the CPU version.
!
!   Usage :
!
!     TYPE(Lagrange) :: interp
!     INTEGER        :: nVariables, nElements
!     REAL(prec)     :: f(0:interp % N,1:nVariables,1:nElements)
!     REAL(prec)     :: derF(0:interp % N,1:nVariables,1:nElements)
!
!       CALL interp % Derivative_1D( f, derF, nVariables, nElements )
!
!     * If CUDA is enabled, the fnative and ftarget arrays must be CUDA device variables.
!
!   Parameters :
!
!     interp (in)
!       A previously constructed Lagrange data-structure.
!
!     f (in)
!       Array of function nodal values at the native interpolation nodes.
!
!     nVariables (in)
!
!     nElements (in)
!
!     derF (out)
!      Array of derivative values at the target interpolation nodes.
!
! ================================================================================================ !

  subroutine Derivative_1D_cpu(this,f,df,nVariables,nElements)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nVariables,nElements
    real(prec),intent(in)  :: f(1:this % N+1,1:nVariables,1:nElements)
    real(prec),intent(out) :: df(1:this % N+1,1:nVariables,1:nElements)
    ! Local
    integer :: i,ii,iVar,iEl
    real(prec) :: dfloc

    do iEl = 1,nElements
      do iVar = 1,nVariables
        do i = 1,this % N+1

          dfloc = 0.0_prec
          do ii = 1,this % N+1
           dfloc = dfloc + this % dMatrix (ii,i)*f(ii,iVar,iEl)
          end do
          df(i,iVar,iEl) = dfloc

        end do
      end do
    end do

  end subroutine Derivative_1D_cpu

  subroutine Derivative_1D_gpu(this,f,df,nVariables,nElements,hipblas_handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in) :: nVariables,nElements
    real(prec), pointer, intent(in)  :: f(:,:,:)
    real(prec), pointer, intent(out) :: df(:,:,:)
    type(c_ptr), intent(inout) :: hipblas_handle
    ! Local
    integer(c_int) :: m
    integer(c_int) :: n
    integer(c_int) :: k
    real(c_prec) :: alpha
    integer(c_int) :: lda
    integer(c_int) :: ldb
    integer(c_int) :: ldc
    real(c_prec) :: beta
 
    m = this % N+1 ! number of rows of A^T
    n = nvariables*nelements ! number of columns of B
    k = this % N+1! number of columns of A^T
    alpha = 1.0_c_prec
    lda = k ! leading dimension of A (interoplation matrix)
    ldb = k ! leading dimension of B (f)
    ldc = m ! leading dimension of C (fTarget)
    beta = 0.0_c_prec
 
#ifdef DOUBLE_PRECISION
    call hipblasCheck(hipblasDgemm(hipblas_handle,&
        HIPBLAS_OP_T, HIPBLAS_OP_N, &
        m, n, k, alpha, &
        c_loc(this % dMatrix), lda, &
        c_loc(f), ldb, &
        beta, &
        c_loc(df), ldc))
#else
    call hipblasCheck(hipblasSgemm(hipblas_handle,&
        HIPBLAS_OP_T, HIPBLAS_OP_N, &
        m, n, k, alpha, &
        c_loc(this % dMatrix), lda, &
        c_loc(f), ldb, &
        beta, &
        c_loc(df), ldc))
#endif

  end subroutine Derivative_1D_gpu

!   subroutine DGDerivative_1D_cpu(this,f,bf,df,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nVariables,nElements
!     real(prec),intent(in)  :: f(0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(in)  :: bf(1:nVariables,1:2,1:nElements)
!     real(prec),intent(out) :: df(0:this % N,1:nVariables,1:nElements)
!     ! Local
!     integer :: i,ii,iVar,iEl

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do i = 0,this % N

!           ! Interior Derivative Matrix Application
!           df(i,iVar,iEl) = 0.0_prec
!           do ii = 0,this % N
!             df(i,iVar,iEl) = df(i,iVar,iEl) + this % dgMatrix (ii,i)*f(ii,iVar,iEl)
!           end do

!           ! Boundary Contribution
!           df(i,iVar,iEl) = df(i,iVar,iEl) + (bf(iVar,2,iEl)*this % bMatrix (i,1) + &
!                                              bf(iVar,1,iEl)*this % bMatrix (i,0))/ &
!                            this % qWeights (i)

!         end do

!       end do
!     end do

!   end subroutine DGDerivative_1D_cpu

!   subroutine DGDerivative_1D_gpu(this,f_dev,bf_dev,df_dev,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in) :: nVariables,nElements
!     type(c_ptr),intent(in)  :: f_dev
!     type(c_ptr),intent(in)  :: bf_dev
!     type(c_ptr),intent(out) :: df_dev

!     call DGDerivative_1D_gpu_wrapper(this % dgMatrix % deviceData, &
!                                      this % bMatrix % deviceData, &
!                                      this % qWeights % deviceData, &
!                                      f_dev,bf_dev,df_dev, &
!                                      this % N, &
!                                      nVariables,nElements)

!   end subroutine DGDerivative_1D_gpu

! ! ================================================================================================ !
! !
! ! CalculateGradient_2D
! !
! !   Calculates the gradient of a 2-D function, represented by a 2-D array of nodal values.
! !
! !   Given a set of nodal values at the interpolation nodes, the gradient of a function through
! !   the interpolation nodes can be estimated by
! !
! !                       (df/dx)_{a,b} = \sum_{i=0}^N f_{i,b} l'_i(\xi_a),   a,b=0,1,2,...,N
! !                       (df/dy)_{a,b} = \sum_{j=0}^N f_{a,j} l'_j(\xi_b),   a,b=0,1,2,...,N
! !
! !   where l_i(\xi) are the Lagrange interpolating polynomials through the interpolation points.
! !   The derivative matrix is D_{a,i} = l'_i(\xi_a) maps an array of nodal values at the interpolation
! !   nodes to its estimated derivative. This routine serves as a wrapper to call either the CUDA
! !   kernel (if CUDA is enabled) or the CPU version.
! !
! !   Usage :
! !
! !     TYPE(Lagrange) :: interp
! !     INTEGER        :: nVariables, nElements
! !     REAL(prec)     :: f(0:interp % N,0:interp % N,1:nVariables,1:nElements)
! !     REAL(prec)     :: gradF(1:2,0:interp % N,0:interp % N,1:nVariables,1:nElements)
! !
! !       CALL interp % CalculateGradient_2D( f, gradF, nVariables, nElements )
! !
! !     * If CUDA is enabled, the fnative and ftarget arrays must be CUDA device variables.
! !
! !   Parameters :
! !
! !     interp (in)
! !       A previously constructed Lagrange data-structure.
! !
! !     f (in)
! !       Array of function nodal values at the native interpolation nodes.
! !
! !     nVariables (in)
! !
! !     nElements (in)
! !
! !     gradF (out)
! !      Array of derivative values at the target interpolation nodes.
! !
! ! ================================================================================================ !
! !
!   subroutine ScalarGradient_2D_cpu(this,f,gradF,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nVariables,nElements
!     real(prec),intent(in)  :: f(0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out) :: gradF(1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!     ! Local
!     integer    :: i,j,ii,iVar,iEl

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do j = 0,this % N
!           do i = 0,this % N

!             gradF(1,i,j,iVar,iEl) = 0.0_prec
!             gradF(2,i,j,iVar,iEl) = 0.0_prec
!             do ii = 0,this % N
!               gradF(1,i,j,iVar,iEl) = gradF(1,i,j,iVar,iEl) + this % dMatrix (ii,i)*f(ii,j,iVar,iEl)
!               gradF(2,i,j,iVar,iEl) = gradF(2,i,j,iVar,iEl) + this % dMatrix (ii,j)*f(i,ii,iVar,iEl)
!             end do

!           end do
!         end do
!       end do
!     end do

!   end subroutine ScalarGradient_2D_cpu

!   subroutine ScalarGradient_2D_gpu(this,f_dev,gradF_dev,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nVariables,nElements
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(out)    :: gradF_dev

!     call ScalarGradient_2D_gpu_wrapper(this % dMatrix % deviceData, &
!                                        f_dev,gradF_dev,this % N, &
!                                        nVariables,nElements)

!   end subroutine ScalarGradient_2D_gpu
! !
! !
!   ! SUBROUTINE ScalarDGGradient_2D_cpu(this,f,bf,gradF,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nVariables,nElements
!   !   REAL(prec),INTENT(in)  :: f(0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   REAL(prec),INTENT(in)  :: bf(0:this % N,1:nVariables,1:4,1:nElements)
!   !   REAL(prec),INTENT(out) :: gradF(1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   ! Local
!   !   INTEGER    :: i,j,ii,iVar,iEl

!   !   DO iEl = 1,nElements
!   !     DO iVar = 1,nVariables
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           gradF(1,i,j,iVar,iEl) = 0.0_prec
!   !           gradF(2,i,j,iVar,iEl) = 0.0_prec
!   !           DO ii = 0,this % N
!   !             gradF(1,i,j,iVar,iEl) = gradF(1,i,j,iVar,iEl) + this % dgMatrix (ii,i)*f(ii,j,iVar,iEl)
!   !             gradF(2,i,j,iVar,iEl) = gradF(2,i,j,iVar,iEl) + this % dgMatrix (ii,j)*f(i,ii,iVar,iEl)
!   !           END DO

!   !           ! Boundary Contribution
!   !           gradF(1,i,j,iVar,iEl) = gradF(1,i,j,iVar,iEl) + (bf(j,iVar,2,iEl)*this % bMatrix (i,1) + &
!   !                                                            bf(j,iVar,4,iEl)*this % bMatrix (i,0))/ &
!   !                                   this % qWeights (i)

!   !           gradF(2,i,j,iVar,iEl) = gradF(2,i,j,iVar,iEl) + (bf(i,iVar,3,iEl)*this % bMatrix (j,1) + &
!   !                                                            bf(i,iVar,1,iEl)*this % bMatrix (j,0))/ &
!   !                                   this % qWeights (j)

!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE ScalarDGGradient_2D_cpu

!   ! SUBROUTINE ScalarDGGradient_2D_gpu(this,f_dev,bf_dev,gradF_dev,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nVariables,nElements
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(in)     :: bf_dev
!   !   TYPE(c_ptr),INTENT(out)    :: gradF_dev

!   !   CALL ScalarDGGradient_2D_gpu_wrapper(this % dgMatrix % deviceData, &
!   !                                        this % bMatrix % deviceData, &
!   !                                        this % qWeights % deviceData, &
!   !                                        f_dev,bf_dev,gradF_dev,this % N, &
!   !                                        nVariables,nElements)

!   ! END SUBROUTINE ScalarDGGradient_2D_gpu

!   subroutine VectorGradient_2D_cpu(this,f,gradF,nVariables,nElements)
!     !
!     ! Input : Vector(1:2,...)
!     ! Output : Tensor(1:2,1:2,....)
!     !          > Tensor(1,1) = d/ds1( Vector(1,...) )
!     !          > Tensor(2,1) = d/ds1( Vector(2,...) )
!     !          > Tensor(1,2) = d/ds2( Vector(1,...) )
!     !          > Tensor(2,2) = d/ds2( Vector(2,...) )
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nVariables,nElements
!     real(prec),intent(in)  :: f(1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out) :: gradF(1:2,1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!     ! Local
!     integer    :: i,j,ii,iVar,iEl
!     real(prec) :: gf(1:2,1:2)

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do j = 0,this % N
!           do i = 0,this % N

!             gf(1,1) = 0.0_prec
!             gf(2,1) = 0.0_prec
!             gf(1,2) = 0.0_prec
!             gf(2,2) = 0.0_prec
!             do ii = 0,this % N
!               gf(1,1) = gf(1,1) + this % dMatrix (ii,i)*f(1,ii,j,iVar,iEl)
!               gf(2,1) = gf(2,1) + this % dMatrix (ii,i)*f(2,ii,j,iVar,iEl)
!               gf(1,2) = gf(1,2) + this % dMatrix (ii,j)*f(1,i,ii,iVar,iEl)
!               gf(2,2) = gf(2,2) + this % dMatrix (ii,j)*f(2,i,ii,iVar,iEl)
!             end do
!             gradF(1:2,1:2,i,j,iVar,iEl) = gf(1:2,1:2)

!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorGradient_2D_cpu

!   subroutine VectorGradient_2D_gpu(this,f_dev,gradF_dev,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nVariables,nElements
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(out)    :: gradF_dev

!     call VectorGradient_2D_gpu_wrapper(this % dMatrix % deviceData, &
!                                        f_dev,gradF_dev,this % N, &
!                                        nVariables,nElements)

!   end subroutine VectorGradient_2D_gpu

!   ! SUBROUTINE VectorDGGradient_2D_cpu(this,f,bf,gradF,nVariables,nElements)
!   !   !
!   !   ! Input : Vector(1:2,...)
!   !   ! Output : Tensor(1:2,1:2,....)
!   !   !          > Tensor(1,1) = d/ds1( Vector(1,...) )
!   !   !          > Tensor(2,1) = d/ds1( Vector(2,...) )
!   !   !          > Tensor(1,2) = d/ds2( Vector(1,...) )
!   !   !          > Tensor(2,2) = d/ds2( Vector(2,...) )
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nVariables,nElements
!   !   REAL(prec),INTENT(in)  :: f(1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   REAL(prec),INTENT(in)  :: bf(1:2,0:this % N,1:nVariables,1:4,1:nElements)
!   !   REAL(prec),INTENT(out) :: gradF(1:2,1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   ! Local
!   !   INTEGER    :: i,j,ii,iVar,iEl

!   !   DO iEl = 1,nElements
!   !     DO iVar = 1,nVariables
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           gradF(1,1,i,j,iVar,iEl) = 0.0_prec
!   !           gradF(2,1,i,j,iVar,iEl) = 0.0_prec
!   !           gradF(1,2,i,j,iVar,iEl) = 0.0_prec
!   !           gradF(2,2,i,j,iVar,iEl) = 0.0_prec
!   !           DO ii = 0,this % N
!   !             gradF(1,1,i,j,iVar,iEl) = gradF(1,1,i,j,iVar,iEl) + this % dgMatrix (ii,i)*f(1,ii,j,iVar,iEl)
!   !             gradF(2,1,i,j,iVar,iEl) = gradF(2,1,i,j,iVar,iEl) + this % dgMatrix (ii,i)*f(2,ii,j,iVar,iEl)
!   !             gradF(1,2,i,j,iVar,iEl) = gradF(1,2,i,j,iVar,iEl) + this % dgMatrix (ii,j)*f(1,i,ii,iVar,iEl)
!   !             gradF(2,2,i,j,iVar,iEl) = gradF(2,2,i,j,iVar,iEl) + this % dgMatrix (ii,j)*f(2,i,ii,iVar,iEl)
!   !           END DO
!   !           gradF(1,1,i,j,iVar,iEl) = gradF(1,1,i,j,iVar,iEl) + (this % bMatrix (i,1)*bf(1,j,iVar,2,iEl) + &
!   !                                                                this % bMatrix (i,0)*bf(1,j,iVar,4,iEl))/ &
!   !                                     this % qWeights (i)

!   !           gradF(2,1,i,j,iVar,iEl) = gradF(2,1,i,j,iVar,iEl) + (this % bMatrix (i,1)*bf(2,j,iVar,2,iEl) + &
!   !                                                                this % bMatrix (i,0)*bf(2,j,iVar,4,iEl))/ &
!   !                                     this % qWeights (i)

!   !           gradF(1,2,i,j,iVar,iEl) = gradF(1,2,i,j,iVar,iEl) + (this % bMatrix (j,1)*bf(1,i,iVar,3,iEl) + &
!   !                                                                this % bMatrix (j,0)*bf(1,i,iVar,1,iEl))/ &
!   !                                     this % qWeights (j)

!   !           gradF(2,2,i,j,iVar,iEl) = gradF(2,2,i,j,iVar,iEl) + (this % bMatrix (j,1)*bf(2,i,iVar,3,iEl) + &
!   !                                                                this % bMatrix (j,0)*bf(2,i,iVar,1,iEl))/ &
!   !                                     this % qWeights (j)

!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE VectorDGGradient_2D_cpu

!   ! SUBROUTINE VectorDGGradient_2D_gpu(this,f_dev,bf_dev,gradF_dev,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nVariables,nElements
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(in)     :: bf_dev
!   !   TYPE(c_ptr),INTENT(out)    :: gradF_dev

!   !   CALL VectorDGGradient_2D_gpu_wrapper(this % dMatrix % deviceData, &
!   !                                        this % bMatrix % deviceData, &
!   !                                        this % qWeights % deviceData, &
!   !                                        f_dev,bf_dev,gradF_dev,this % N, &
!   !                                        nVariables,nElements)

!   ! END SUBROUTINE VectorDGGradient_2D_gpu

!   subroutine VectorDivergence_2D_cpu(this,f,dF,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nVariables,nElements
!     real(prec),intent(in)  :: f(1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out) :: dF(0:this % N,0:this % N,1:nVariables,1:nElements)
!     ! Local
!     integer    :: i,j,ii,iVar,iEl

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do j = 0,this % N
!           do i = 0,this % N

!             dF(i,j,iVar,iEl) = 0.0_prec
!             do ii = 0,this % N
!               dF(i,j,iVar,iEl) = dF(i,j,iVar,iEl) + this % dMatrix (ii,i)*f(1,ii,j,iVar,iEl) + &
!                                  this % dMatrix (ii,j)*f(2,i,ii,iVar,iEl)
!             end do

!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorDivergence_2D_cpu

!   subroutine VectorDivergence_2D_gpu(this,f_dev,dF_dev,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nVariables,nElements
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(out)    :: dF_dev

!     call VectorDivergence_2D_gpu_wrapper(this % dMatrix % deviceData, &
!                                          f_dev,dF_dev,this % N, &
!                                          nVariables,nElements)

!   end subroutine VectorDivergence_2D_gpu

!   subroutine VectorDGDivergence_2D_cpu(this,f,bF,dF,nVariables,nElements)
!     ! Assumes bF is the vector component in the direction normal to the boundary
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nVariables,nElements
!     real(prec),intent(in)  :: f(1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(in)  :: bF(0:this % N,1:nVariables,1:4,1:nElements)
!     real(prec),intent(out) :: dF(0:this % N,0:this % N,1:nVariables,1:nElements)
!     ! Local
!     real(prec) :: dfLoc
!     integer    :: i,j,ii,iVar,iEl

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do j = 0,this % N
!           do i = 0,this % N

!             dfLoc = 0.0_prec
!             do ii = 0,this % N
!               dfLoc = dfLoc + this % dgMatrix (ii,i)*f(1,ii,j,iVar,iEl) + &
!                       this % dgMatrix (ii,j)*f(2,i,ii,iVar,iEl)
!             end do

!             dfLoc = dfLoc + (this % bMatrix (i,1)*bF(j,iVar,2,iEl) + &
!                              this % bMatrix (i,0)*bF(j,iVar,4,iEl))/ &
!                     this % qWeights (i) + &
!                     (this % bMatrix (j,1)*bF(i,iVar,3,iEl) + &
!                      this % bMatrix (j,0)*bF(i,iVar,1,iEl))/ &
!                     this % qWeights (j)
!             dF(i,j,iVar,iEl) = dFLoc

!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorDGDivergence_2D_cpu

!   subroutine VectorDGDivergence_2D_gpu(this,f_dev,bF_dev,dF_dev,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nVariables,nElements
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(in)     :: bF_dev
!     type(c_ptr),intent(out)    :: dF_dev

!     call VectorDGDivergence_2D_gpu_wrapper(this % dgMatrix % deviceData, &
!                                            this % bMatrix % deviceData, &
!                                            this % qWeights % deviceData, &
!                                            f_dev,bF_dev,dF_dev,this % N, &
!                                            nVariables,nElements)

!   end subroutine VectorDGDivergence_2D_gpu

!   ! SUBROUTINE VectorCurl_2D_cpu(this,f,dF,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nVariables,nElements
!   !   REAL(prec),INTENT(in)  :: f(1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   REAL(prec),INTENT(out) :: dF(0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   ! Local
!   !   INTEGER    :: i,j,ii,iVar,iEl

!   !   DO iEl = 1,nElements
!   !     DO iVar = 1,nVariables
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           dF(i,j,iVar,iEl) = 0.0_prec
!   !           DO ii = 0,this % N
!   !             dF(i,j,iVar,iEl) = dF(i,j,iVar,iEl) + this % dMatrix (ii,j)*f(1,i,ii,iVar,iEl) - &
!   !                                this % dMatrix (ii,i)*f(2,ii,j,iVar,iEl)
!   !           END DO

!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE VectorCurl_2D_cpu

!   ! SUBROUTINE VectorCurl_2D_gpu(this,f_dev,dF_dev,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nVariables,nElements
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL VectorCurl_2D_gpu_wrapper(this % dMatrix % deviceData, &
!   !                                  f_dev,dF_dev,this % N, &
!   !                                  nVariables,nElements)

!   ! END SUBROUTINE VectorCurl_2D_gpu

!   ! SUBROUTINE P2VectorDivergence_2D_cpu(this,f,dF,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nVariables,nElements
!   !   REAL(prec),INTENT(in)  :: f(1:2,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   REAL(prec),INTENT(out) :: dF(0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   ! Local
!   !   INTEGER    :: i,j,n,iVar,iEl
!   !   REAL(prec) :: dfloc

!   !   DO iEl = 1,nElements
!   !     DO iVar = 1,nVariables
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           dfloc = 0.0_prec
!   !           DO n = 0,this % N
!   !             dfloc = dfloc + this % dMatrix (n,i)*f(1,n,i,j,iVar,iEl) + &
!   !                             this % dMatrix (n,j)*f(2,n,i,j,iVar,iEl)
!   !           END DO

!   !           dF(i,j,iVar,iEl) = 2.0_prec*dfloc

!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE P2VectorDivergence_2D_cpu

!   ! SUBROUTINE P2VectorDivergence_2D_gpu(this,f_dev,dF_dev,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nVariables,nElements
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL P2VectorDivergence_2D_gpu_wrapper(this % dMatrix % deviceData, &
!   !                                        f_dev,dF_dev,this % N, &
!   !                                        nVariables,nElements)

!   ! END SUBROUTINE P2VectorDivergence_2D_gpu

!   ! SUBROUTINE P2VectorDGDivergence_2D_cpu(this,f,bF,dF,nVariables,nElements)
!   !   ! Assumes bF is the vector component in the direction normal to the boundary
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nVariables,nElements
!   !   REAL(prec),INTENT(in)  :: f(1:2,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   REAL(prec),INTENT(in)  :: bF(0:this % N,1:nVariables,1:4,1:nElements)
!   !   REAL(prec),INTENT(out) :: dF(0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   ! Local
!   !   REAL(prec) :: dfLoc
!   !   INTEGER    :: i,j,n,iVar,iEl

!   !   DO iEl = 1,nElements
!   !     DO iVar = 1,nVariables
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           dfLoc = 0.0_prec
!   !           DO n = 0,this % N
!   !             dfLoc = dfLoc + this % dgMatrix (n,i)*f(1,n,i,j,iVar,iEl) + &
!   !                             this % dgMatrix (n,j)*f(2,n,i,j,iVar,iEl)
!   !           END DO

!   !           dfLoc = dfLoc + (this % bMatrix (i,1)*bF(j,iVar,2,iEl) + &
!   !                            this % bMatrix (i,0)*bF(j,iVar,4,iEl))/ &
!   !                              this % qWeights (i) + &
!   !                           (this % bMatrix (j,1)*bF(i,iVar,3,iEl) + &
!   !                            this % bMatrix (j,0)*bF(i,iVar,1,iEl))/ &
!   !                              this % qWeights (j)
!   !           dF(i,j,iVar,iEl) = dFLoc

!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE P2VectorDGDivergence_2D_cpu

!   ! SUBROUTINE P2VectorDGDivergence_2D_gpu(this,f_dev,bF_dev,dF_dev,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nVariables,nElements
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(in)     :: bF_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL P2VectorDGDivergence_2D_gpu_wrapper(this % dgMatrix % deviceData, &
!   !                                          this % bMatrix % deviceData, &
!   !                                          this % qWeights % deviceData, &
!   !                                          f_dev,bF_dev,dF_dev,this % N, &
!   !                                          nVariables,nElements)

!   ! END SUBROUTINE P2VectorDGDivergence_2D_gpu

!   ! SUBROUTINE TensorDivergence_2D_cpu(this,f,dF,nVariables,nElements)
!   !   ! Note that the divergence is taken over the first dimension (row dimension) of the tensor matrix
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nVariables,nElements
!   !   REAL(prec),INTENT(in)  :: f(1:2,1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   REAL(prec),INTENT(out) :: dF(1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   ! Local
!   !   INTEGER    :: i,j,ii,iVar,iEl

!   !   DO iEl = 1,nElements
!   !     DO iVar = 1,nVariables
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           dF(1,i,j,iVar,iEl) = 0.0_prec
!   !           dF(2,i,j,iVar,iEl) = 0.0_prec
!   !           DO ii = 0,this % N
!   !             dF(1,i,j,iVar,iEl) = dF(1,i,j,iVar,iEl) + this % dMatrix (ii,i)*f(1,1,ii,j,iVar,iEl) + &
!   !                                  this % dMatrix (ii,j)*f(2,1,i,ii,iVar,iEl)
!   !             dF(2,i,j,iVar,iEl) = dF(2,i,j,iVar,iEl) + this % dMatrix (ii,i)*f(1,2,ii,j,iVar,iEl) + &
!   !                                  this % dMatrix (ii,j)*f(2,2,i,ii,iVar,iEl)
!   !           END DO

!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE TensorDivergence_2D_cpu

!   ! SUBROUTINE TensorDivergence_2D_gpu(this,f_dev,dF_dev,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nVariables,nElements
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL TensorDivergence_2D_gpu_wrapper(this % dMatrix % deviceData, &
!   !                                        f_dev,dF_dev,this % N, &
!   !                                        nVariables,nElements)

!   ! END SUBROUTINE TensorDivergence_2D_gpu

!   ! SUBROUTINE TensorDGDivergence_2D_cpu(this,f,bF,dF,nVariables,nElements)
!   !   ! Note that the divergence is taken over the first dimension (row dimension) of the tensor matrix
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nVariables,nElements
!   !   REAL(prec),INTENT(in)  :: f(1:2,1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   REAL(prec),INTENT(in)  :: bf(1:2,1:2,0:this % N,1:nVariables,1:4,1:nElements)
!   !   REAL(prec),INTENT(out) :: dF(1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   ! Local
!   !   INTEGER    :: i,j,ii,iVar,iEl

!   !   DO iEl = 1,nElements
!   !     DO iVar = 1,nVariables
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           dF(1,i,j,iVar,iEl) = 0.0_prec
!   !           dF(2,i,j,iVar,iEl) = 0.0_prec
!   !           DO ii = 0,this % N
!   !             dF(1,i,j,iVar,iEl) = dF(1,i,j,iVar,iEl) + this % dgMatrix (ii,i)*f(1,1,ii,j,iVar,iEl) + &
!   !                                  this % dgMatrix (ii,j)*f(2,1,i,ii,iVar,iEl)
!   !             dF(2,i,j,iVar,iEl) = dF(2,i,j,iVar,iEl) + this % dgMatrix (ii,i)*f(1,2,ii,j,iVar,iEl) + &
!   !                                  this % dgMatrix (ii,j)*f(2,2,i,ii,iVar,iEl)
!   !           END DO

!   !           dF(1,i,j,iVar,iEl) = dF(1,i,j,iVar,iEl) + (this % bMatrix (i,1)*bf(1,1,j,iVar,2,iEl) + &
!   !                                                      this % bMatrix (i,0)*bf(1,1,j,iVar,4,iEl))/ &
!   !                                this % qWeights (i) + &
!   !                                (this % bMatrix (j,1)*bf(2,1,i,iVar,3,iEl) + &
!   !                                 this % bMatrix (j,0)*bf(2,1,i,iVar,1,iEl))/ &
!   !                                this % qWeights (j)

!   !           dF(2,i,j,iVar,iEl) = dF(2,i,j,iVar,iEl) + (this % bMatrix (i,1)*bf(1,2,j,iVar,2,iEl) + &
!   !                                                      this % bMatrix (i,0)*bf(1,2,j,iVar,4,iEl))/ &
!   !                                this % qWeights (i) + &
!   !                                (this % bMatrix (j,1)*bf(2,2,i,iVar,3,iEl) + &
!   !                                 this % bMatrix (j,0)*bf(2,2,i,iVar,1,iEl))/ &
!   !                                this % qWeights (j)
!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE TensorDGDivergence_2D_cpu

!   ! SUBROUTINE TensorDGDivergence_2D_gpu(this,f_dev,bF_dev,dF_dev,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nVariables,nElements
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(in)     :: bf_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL TensorDGDivergence_2D_gpu_wrapper(this % dgMatrix % deviceData, &
!   !                                          this % bMatrix % deviceData, &
!   !                                          this % qWeights % deviceData, &
!   !                                          f_dev,bF_dev,dF_dev,this % N, &
!   !                                          nVariables,nElements)

!   ! END SUBROUTINE TensorDGDivergence_2D_gpu

!   subroutine ScalarGradient_3D_cpu(this,f,gradF,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nVariables,nElements
!     real(prec),intent(in)  :: f(0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out) :: gradF(1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     ! Local
!     integer    :: i,j,k,ii,iVar,iEl
!     real(prec) :: gf(1:3)

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do k = 0,this % N
!           do j = 0,this % N
!             do i = 0,this % N

!               gF(1) = 0.0_prec
!               gF(2) = 0.0_prec
!               gF(3) = 0.0_prec
!               do ii = 0,this % N
!                 gF(1) = gF(1) + this % dMatrix (ii,i)*f(ii,j,k,iVar,iEl)
!                 gF(2) = gF(2) + this % dMatrix (ii,j)*f(i,ii,k,iVar,iEl)
!                 gF(3) = gF(3) + this % dMatrix (ii,k)*f(i,j,ii,iVar,iEl)
!               end do

!               gradF(1,i,j,k,iVar,iEl) = gF(1)
!               gradF(2,i,j,k,iVar,iEl) = gF(2)
!               gradF(3,i,j,k,iVar,iEl) = gF(3)

!             end do
!           end do
!         end do
!       end do
!     end do

!   end subroutine ScalarGradient_3D_cpu

!   subroutine ScalarGradient_3D_gpu(this,f_dev,gradF_dev,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nVariables,nElements
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(out)    :: gradF_dev

!     call ScalarGradient_3D_gpu_wrapper(this % dMatrix % deviceData, &
!                                        f_dev,gradF_dev,this % N, &
!                                        nVariables,nElements)

!   end subroutine ScalarGradient_3D_gpu
! !
!   subroutine VectorGradient_3D_cpu(this,f,gradF,nVariables,nElements)
!     !
!     ! Input : Vector(1:3,...)
!     ! Output : Tensor(1:3,1:3,....)
!     !          > Tensor(1,1) = d/ds1( Vector(1,...) )
!     !          > Tensor(2,1) = d/ds1( Vector(2,...) )
!     !          > Tensor(3,1) = d/ds1( Vector(3,...) )
!     !          > Tensor(1,2) = d/ds2( Vector(1,...) )
!     !          > Tensor(2,2) = d/ds2( Vector(2,...) )
!     !          > Tensor(3,2) = d/ds2( Vector(3,...) )
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nVariables,nElements
!     real(prec),intent(in)  :: f(1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out) :: gradF(1:3,1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     ! Local
!     integer    :: i,j,k,ii,iVar,iEl
!     real(prec) :: gF(1:3,1:3)

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do k = 0,this % N
!           do j = 0,this % N
!             do i = 0,this % N

!               gF = 0.0_prec
!               do ii = 0,this % N
!                 gF(1,1) = gF(1,1) + this % dMatrix (ii,i)*f(1,ii,j,k,iVar,iEl)
!                 gF(2,1) = gF(2,1) + this % dMatrix (ii,i)*f(2,ii,j,k,iVar,iEl)
!                 gF(3,1) = gF(3,1) + this % dMatrix (ii,i)*f(3,ii,j,k,iVar,iEl)
!                 gF(1,2) = gF(1,2) + this % dMatrix (ii,j)*f(1,i,ii,k,iVar,iEl)
!                 gF(2,2) = gF(2,2) + this % dMatrix (ii,j)*f(2,i,ii,k,iVar,iEl)
!                 gF(3,2) = gF(3,2) + this % dMatrix (ii,j)*f(3,i,ii,k,iVar,iEl)
!                 gF(1,3) = gF(1,3) + this % dMatrix (ii,k)*f(1,i,j,ii,iVar,iEl)
!                 gF(2,3) = gF(2,3) + this % dMatrix (ii,k)*f(2,i,j,ii,iVar,iEl)
!                 gF(3,3) = gF(3,3) + this % dMatrix (ii,k)*f(3,i,j,ii,iVar,iEl)
!               end do

!               gradF(1,1,i,j,k,iVar,iEl) = gF(1,1)
!               gradF(2,1,i,j,k,iVar,iEl) = gF(2,1)
!               gradF(3,1,i,j,k,iVar,iEl) = gF(3,1)
!               gradF(1,2,i,j,k,iVar,iEl) = gF(1,2)
!               gradF(2,2,i,j,k,iVar,iEl) = gF(2,2)
!               gradF(3,2,i,j,k,iVar,iEl) = gF(3,2)
!               gradF(1,3,i,j,k,iVar,iEl) = gF(1,3)
!               gradF(2,3,i,j,k,iVar,iEl) = gF(2,3)
!               gradF(3,3,i,j,k,iVar,iEl) = gF(3,3)

!             end do
!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorGradient_3D_cpu

!   subroutine VectorGradient_3D_gpu(this,f_dev,gradF_dev,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nVariables,nElements
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(out)    :: gradF_dev

!     call VectorGradient_3D_gpu_wrapper(this % dMatrix % deviceData, &
!                                        f_dev,gradF_dev,this % N, &
!                                        nVariables,nElements)

!   end subroutine VectorGradient_3D_gpu

!   subroutine VectorDivergence_3D_cpu(this,f,dF,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nVariables,nElements
!     real(prec),intent(in)  :: f(1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out) :: dF(0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     ! Local
!     integer    :: i,j,k,ii,iVar,iEl

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do k = 0,this % N
!           do j = 0,this % N
!             do i = 0,this % N

!               dF(i,j,k,iVar,iEl) = 0.0_prec
!               do ii = 0,this % N
!                 dF(i,j,k,iVar,iEl) = dF(i,j,k,iVar,iEl) + this % dMatrix (ii,i)*f(1,ii,j,k,iVar,iEl) + &
!                                      this % dMatrix (ii,j)*f(2,i,ii,k,iVar,iEl) + &
!                                      this % dMatrix (ii,k)*f(3,i,j,ii,iVar,iEl)
!               end do

!             end do
!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorDivergence_3D_cpu

!   subroutine VectorDivergence_3D_gpu(this,f_dev,dF_dev,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nVariables,nElements
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(out)    :: dF_dev

!     call VectorDivergence_3D_gpu_wrapper(this % dMatrix % deviceData, &
!                                          f_dev,dF_dev,this % N, &
!                                          nVariables,nElements)

!   end subroutine VectorDivergence_3D_gpu

!   subroutine VectorDGDivergence_3D_cpu(this,f,bF,dF,nVariables,nElements)
!     ! Assumes bF is the vector component in the direction normal to the element boundaries
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nVariables,nElements
!     real(prec),intent(in)  :: f(1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(in)  :: bf(0:this % N,0:this % N,1:nVariables,1:6,1:nElements)
!     real(prec),intent(out) :: dF(0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     ! Local
!     integer    :: i,j,k,ii,iVar,iEl

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do k = 0,this % N
!           do j = 0,this % N
!             do i = 0,this % N

!               dF(i,j,k,iVar,iEl) = 0.0_prec
!               do ii = 0,this % N
!                 dF(i,j,k,iVar,iEl) = dF(i,j,k,iVar,iEl) + this % dgMatrix (ii,i)*f(1,ii,j,k,iVar,iEl) + &
!                                      this % dgMatrix (ii,j)*f(2,i,ii,k,iVar,iEl) + &
!                                      this % dgMatrix (ii,k)*f(3,i,j,ii,iVar,iEl)
!               end do

!               dF(i,j,k,iVar,iEl) = dF(i,j,k,iVar,iEl) + (this % bMatrix (i,1)*bF(j,k,iVar,3,iEl) + & ! east
!                                                          this % bMatrix (i,0)*bF(j,k,iVar,5,iEl))/ &  ! west
!                                    this % qWeights (i) + &
!                                    (this % bMatrix (j,1)*bF(i,k,iVar,4,iEl) + & ! north
!                                     this % bMatrix (j,0)*bF(i,k,iVar,2,iEl))/ &  ! south
!                                    this % qWeights (j) + &
!                                    (this % bMatrix (k,1)*bF(i,j,iVar,6,iEl) + & ! top
!                                     this % bMatrix (k,0)*bF(i,j,iVar,1,iEl))/ &  ! bottom
!                                    this % qWeights (k)

!             end do
!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorDGDivergence_3D_cpu

!   subroutine VectorDGDivergence_3D_gpu(this,f_dev,bF_dev,dF_dev,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nVariables,nElements
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(in)     :: bF_dev
!     type(c_ptr),intent(out)    :: dF_dev

!     call VectorDGDivergence_3D_gpu_wrapper(this % dgMatrix % deviceData, &
!                                            this % bMatrix % deviceData, &
!                                            this % qWeights % deviceData, &
!                                            f_dev,bF_dev,dF_dev,this % N, &
!                                            nVariables,nElements)

!   end subroutine VectorDGDivergence_3D_gpu

!   ! SUBROUTINE VectorCurl_3D_cpu(this,f,dF,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nVariables,nElements
!   !   REAL(prec),INTENT(in)  :: f(1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   REAL(prec),INTENT(out) :: dF(1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   ! Local
!   !   INTEGER    :: i,j,k,ii,iVar,iEl

!   !   DO iEl = 1,nElements
!   !     DO iVar = 1,nVariables
!   !       DO k = 0,this % N
!   !         DO j = 0,this % N
!   !           DO i = 0,this % N

!   !             dF(1,i,j,k,iVar,iEl) = 0.0_prec
!   !             dF(2,i,j,k,iVar,iEl) = 0.0_prec
!   !             dF(3,i,j,k,iVar,iEl) = 0.0_prec
!   !             DO ii = 0,this % N
!   !               dF(1,i,j,k,iVar,iEl) = dF(1,i,j,k,iVar,iEl) + this % dMatrix (ii,j)*f(3,i,ii,k,iVar,iEl) - &
!   !                                      this % dMatrix (ii,k)*f(2,i,j,ii,iVar,iEl)
!   !               dF(2,i,j,k,iVar,iEl) = dF(2,i,j,k,iVar,iEl) + this % dMatrix (ii,k)*f(1,i,j,ii,iVar,iEl) - &
!   !                                      this % dMatrix (ii,i)*f(3,ii,j,k,iVar,iEl)
!   !               dF(3,i,j,k,iVar,iEl) = dF(3,i,j,k,iVar,iEl) + this % dMatrix (ii,i)*f(2,ii,j,k,iVar,iEl) - &
!   !                                      this % dMatrix (ii,j)*f(1,i,ii,k,iVar,iEl)
!   !             END DO

!   !           END DO
!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE VectorCurl_3D_cpu

!   ! SUBROUTINE VectorCurl_3D_gpu(this,f_dev,dF_dev,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nVariables,nElements
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL VectorCurl_3D_gpu_wrapper(this % dMatrix % deviceData, &
!   !                                  f_dev,dF_dev,this % N, &
!   !                                  nVariables,nElements)

!   ! END SUBROUTINE VectorCurl_3D_gpu

!   ! SUBROUTINE TensorDivergence_3D_cpu(this,f,dF,nVariables,nElements)
!   !   ! Note that the divergence is taken over the first dimension (row dimension) of the tensor matrix
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nVariables,nElements
!   !   REAL(prec),INTENT(in)  :: f(1:3,1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   REAL(prec),INTENT(out) :: dF(1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   ! Local
!   !   INTEGER    :: i,j,k,ii,iVar,iEl

!   !   DO iEl = 1,nElements
!   !     DO iVar = 1,nVariables
!   !       DO k = 0,this % N
!   !         DO j = 0,this % N
!   !           DO i = 0,this % N

!   !             dF(1,i,j,k,iVar,iEl) = 0.0_prec
!   !             dF(2,i,j,k,iVar,iEl) = 0.0_prec
!   !             dF(3,i,j,k,iVar,iEl) = 0.0_prec
!   !             DO ii = 0,this % N
!   !               dF(1,i,j,k,iVar,iEl) = dF(1,i,j,k,iVar,iEl) + &
!   !                                      this % dMatrix (ii,i)*f(1,1,ii,j,k,iVar,iEl) + &
!   !                                      this % dMatrix (ii,j)*f(2,1,i,ii,k,iVar,iEl) + &
!   !                                      this % dMatrix (ii,k)*f(3,1,i,j,ii,iVar,iEl)

!   !               dF(2,i,j,k,iVar,iEl) = dF(2,i,j,k,iVar,iEl) + &
!   !                                      this % dMatrix (ii,i)*f(1,2,ii,j,k,iVar,iEl) + &
!   !                                      this % dMatrix (ii,j)*f(2,2,i,ii,k,iVar,iEl) + &
!   !                                      this % dMatrix (ii,k)*f(3,2,i,j,ii,iVar,iEl)

!   !               dF(3,i,j,k,iVar,iEl) = dF(3,i,j,k,iVar,iEl) + &
!   !                                      this % dMatrix (ii,i)*f(1,3,ii,j,k,iVar,iEl) + &
!   !                                      this % dMatrix (ii,j)*f(2,3,i,ii,k,iVar,iEl) + &
!   !                                      this % dMatrix (ii,k)*f(3,3,i,j,ii,iVar,iEl)
!   !             END DO

!   !           END DO
!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE TensorDivergence_3D_cpu

!   ! SUBROUTINE TensorDivergence_3D_gpu(this,f_dev,dF_dev,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nVariables,nElements
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL TensorDivergence_3D_gpu_wrapper(this % dMatrix % deviceData, &
!   !                                        f_dev,dF_dev,this % N, &
!   !                                        nVariables,nElements)

!   ! END SUBROUTINE TensorDivergence_3D_gpu

!   ! SUBROUTINE TensorDGDivergence_3D_cpu(this,f,bF,dF,nVariables,nElements)
!   !   ! Note that the divergence is taken over the first dimension (row dimension) of the tensor matrix
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nVariables,nElements
!   !   REAL(prec),INTENT(in)  :: f(1:3,1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   REAL(prec),INTENT(in)  :: bF(1:3,1:3,0:this % N,0:this % N,1:nVariables,1:6,1:nElements)
!   !   REAL(prec),INTENT(out) :: dF(1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!   !   ! Local
!   !   INTEGER    :: i,j,k,ii,iVar,iEl

!   !   DO iEl = 1,nElements
!   !     DO iVar = 1,nVariables
!   !       DO k = 0,this % N
!   !         DO j = 0,this % N
!   !           DO i = 0,this % N

!   !             dF(1,i,j,k,iVar,iEl) = 0.0_prec
!   !             dF(2,i,j,k,iVar,iEl) = 0.0_prec
!   !             dF(3,i,j,k,iVar,iEl) = 0.0_prec
!   !             DO ii = 0,this % N
!   !               dF(1,i,j,k,iVar,iEl) = dF(1,i,j,k,iVar,iEl) + &
!   !                                      this % dgMatrix (ii,i)*f(1,1,ii,j,k,iVar,iEl) + &
!   !                                      this % dgMatrix (ii,j)*f(2,1,i,ii,k,iVar,iEl) + &
!   !                                      this % dgMatrix (ii,k)*f(3,1,i,j,ii,iVar,iEl)

!   !               dF(2,i,j,k,iVar,iEl) = dF(2,i,j,k,iVar,iEl) + &
!   !                                      this % dgMatrix (ii,i)*f(1,2,ii,j,k,iVar,iEl) + &
!   !                                      this % dgMatrix (ii,j)*f(2,2,i,ii,k,iVar,iEl) + &
!   !                                      this % dgMatrix (ii,k)*f(3,2,i,j,ii,iVar,iEl)

!   !               dF(3,i,j,k,iVar,iEl) = dF(3,i,j,k,iVar,iEl) + &
!   !                                      this % dgMatrix (ii,i)*f(1,3,ii,j,k,iVar,iEl) + &
!   !                                      this % dgMatrix (ii,j)*f(2,3,i,ii,k,iVar,iEl) + &
!   !                                      this % dgMatrix (ii,k)*f(3,3,i,j,ii,iVar,iEl)
!   !             END DO

!   !             dF(1,i,j,k,iVar,iEl) = dF(1,i,j,k,iVar,iEl) + (this % bMatrix (i,1)*bF(1,1,j,k,iVar,3,iEl) + & ! east
!   !                                                    this % bMatrix (i,0)*bF(1,1,j,k,iVar,5,iEl))/ &  ! west
!   !                                    this % qWeights (i) + &
!   !                                    (this % bMatrix (j,1)*bF(2,1,i,k,iVar,4,iEl) + & ! north
!   !                                     this % bMatrix (j,0)*bF(2,1,i,k,iVar,2,iEl))/ &  ! south
!   !                                    this % qWeights (j) + &
!   !                                    (this % bMatrix (k,1)*bF(3,1,i,j,iVar,6,iEl) + & ! top
!   !                                     this % bMatrix (k,0)*bF(3,1,i,j,iVar,1,iEl))/ &  ! bottom
!   !                                    this % qWeights (k)

!   !             dF(2,i,j,k,iVar,iEl) = dF(2,i,j,k,iVar,iEl) + (this % bMatrix (i,1)*bF(1,2,j,k,iVar,3,iEl) + & ! east
!   !                                                    this % bMatrix (i,0)*bF(1,2,j,k,iVar,5,iEl))/ &  ! west
!   !                                    this % qWeights (i) + &
!   !                                    (this % bMatrix (j,1)*bF(2,2,i,k,iVar,4,iEl) + & ! north
!   !                                     this % bMatrix (j,0)*bF(2,2,i,k,iVar,2,iEl))/ &  ! south
!   !                                    this % qWeights (j) + &
!   !                                    (this % bMatrix (k,1)*bF(3,2,i,j,iVar,6,iEl) + & ! top
!   !                                     this % bMatrix (k,0)*bF(3,2,i,j,iVar,1,iEl))/ &  ! bottom
!   !                                    this % qWeights (k)

!   !             dF(3,i,j,k,iVar,iEl) = dF(3,i,j,k,iVar,iEl) + (this % bMatrix (i,1)*bF(1,3,j,k,iVar,3,iEl) + & ! east
!   !                                                    this % bMatrix (i,0)*bF(1,3,j,k,iVar,5,iEl))/ &  ! west
!   !                                    this % qWeights (i) + &
!   !                                    (this % bMatrix (j,1)*bF(2,3,i,k,iVar,4,iEl) + & ! north
!   !                                     this % bMatrix (j,0)*bF(2,3,i,k,iVar,2,iEl))/ &  ! south
!   !                                    this % qWeights (j) + &
!   !                                    (this % bMatrix (k,1)*bF(3,3,i,j,iVar,6,iEl) + & ! top
!   !                                     this % bMatrix (k,0)*bF(3,3,i,j,iVar,1,iEl))/ &  ! bottom
!   !                                    this % qWeights (k)

!   !           END DO
!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE TensorDGDivergence_3D_cpu

!   ! SUBROUTINE TensorDGDivergence_3D_gpu(this,f_dev,bF_dev,dF_dev,nVariables,nElements)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nVariables,nElements
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(in)     :: bF_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL TensorDGDivergence_3D_gpu_wrapper(this % dgMatrix % deviceData, &
!   !                                          this % bMatrix % deviceData, &
!   !                                          this % qWeights % deviceData, &
!   !                                          f_dev,bF_dev,dF_dev,this % N, &
!   !                                          nVariables,nElements)

!   ! END SUBROUTINE TensorDGDivergence_3D_gpu
!   ! /////////////////////////////// !
!   ! Boundary Interpolation Routines !

!   subroutine ScalarBoundaryInterp_1D_cpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nVariables,nElements
!     real(prec),intent(in)      :: f(0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out)     :: fBound(1:nVariables,1:2,1:nElements)
!     ! Local
!     integer :: ii,iVar,iEl
!     real(prec) :: fb(1:2)

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         fb(1:2) = 0.0_prec
!         do ii = 0,this % N
!           fb(1) = fb(1) + this % bMatrix (ii,0)*f(ii,iVar,iEl) ! West
!           fb(2) = fb(2) + this % bMatrix (ii,1)*f(ii,iVar,iEl) ! East
!         end do
!         fBound(iVar,1:2,iEl) = fb(1:2)
!       end do
!     end do

!   end subroutine ScalarBoundaryInterp_1D_cpu

!   subroutine ScalarBoundaryInterp_1D_gpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nVariables,nElements
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fBound

!     call ScalarBoundaryInterp_1D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fBound,this % N,nVariables,nElements)

!   end subroutine ScalarBoundaryInterp_1D_gpu

!   subroutine ScalarBoundaryInterp_2D_cpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nVariables,nElements
!     real(prec),intent(in)      :: f(0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out)     :: fBound(0:this % N,1:nVariables,1:4,1:nElements)
!     ! Local
!     integer :: i,ii,iVar,iEl
!     real(prec) :: fb(1:4)

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do i = 0,this % N

!           fb(1:4) = 0.0_prec

!           do ii = 0,this % N
!             fb(1) = fb(1) + this % bMatrix (ii,0)*f(i,ii,iVar,iEl) ! South
!             fb(2) = fb(2) + this % bMatrix (ii,1)*f(ii,i,iVar,iEl) ! East
!             fb(3) = fb(3) + this % bMatrix (ii,1)*f(i,ii,iVar,iEl) ! North
!             fb(4) = fb(4) + this % bMatrix (ii,0)*f(ii,i,iVar,iEl) ! West
!           end do

!           fBound(i,iVar,1:4,iEl) = fb(1:4)

!         end do
!       end do
!     end do

!   end subroutine ScalarBoundaryInterp_2D_cpu

!   subroutine ScalarBoundaryInterp_2D_gpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nVariables,nElements
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fBound

!     call ScalarBoundaryInterp_2D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fBound,this % N,nVariables,nElements)

!   end subroutine ScalarBoundaryInterp_2D_gpu

!   subroutine VectorBoundaryInterp_2D_cpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nVariables,nElements
!     real(prec),intent(in)  :: f(1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out)  :: fBound(1:2,0:this % N,1:nVariables,1:4,1:nElements)
!     ! Local
!     integer :: i,ii,idir,iVar,iEl
!     real(prec) :: fb(1:2,1:4)

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do i = 0,this % N

!           fb(1:2,1:4) = 0.0_prec
!           do ii = 0,this % N
!             do idir = 1,2
!               fb(idir,1) = fb(idir,1) + this % bMatrix (ii,0)*f(idir,i,ii,iVar,iEl) ! South
!               fb(idir,2) = fb(idir,2) + this % bMatrix (ii,1)*f(idir,ii,i,iVar,iEl) ! East
!               fb(idir,3) = fb(idir,3) + this % bMatrix (ii,1)*f(idir,i,ii,iVar,iEl) ! North
!               fb(idir,4) = fb(idir,4) + this % bMatrix (ii,0)*f(idir,ii,i,iVar,iEl) ! West
!             end do
!           end do

!           do idir = 1,2
!             fBound(idir,i,iVar,1:4,iEl) = fb(idir,1:4)
!           end do

!         end do
!       end do
!     end do

!   end subroutine VectorBoundaryInterp_2D_cpu

!   subroutine VectorBoundaryInterp_2D_gpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nVariables,nElements
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fBound

!     call VectorBoundaryInterp_2D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fBound,this % N,nVariables,nElements)

!   end subroutine VectorBoundaryInterp_2D_gpu

!   subroutine TensorBoundaryInterp_2D_cpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nVariables,nElements
!     real(prec),intent(in)  :: f(1:2,1:2,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out)  :: fBound(1:2,1:2,0:this % N,1:nVariables,1:4,1:nElements)
!     ! Local
!     integer :: i,ii,idir,jdir,iVar,iEl
!     real(prec) :: fb(1:2,1:2,1:4)

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do i = 0,this % N

!           fb(1:2,1:2,1:4) = 0.0_prec
!           do ii = 0,this % N
!             do jdir = 1,2
!               do idir = 1,2
!                 fb(idir,jdir,1) = fb(idir,jdir,1) + this % bMatrix (ii,0)*f(idir,jdir,i,ii,iVar,iEl) ! South
!                 fb(idir,jdir,2) = fb(idir,jdir,2) + this % bMatrix (ii,1)*f(idir,jdir,ii,i,iVar,iEl) ! East
!                 fb(idir,jdir,3) = fb(idir,jdir,3) + this % bMatrix (ii,1)*f(idir,jdir,i,ii,iVar,iEl) ! North
!                 fb(idir,jdir,4) = fb(idir,jdir,4) + this % bMatrix (ii,0)*f(idir,jdir,ii,i,iVar,iEl) ! West
!               end do
!             end do
!           end do

!           do jdir = 1,2
!             do idir = 1,2
!               fBound(idir,jdir,i,iVar,1:4,iEl) = fb(idir,jdir,1:4)
!             end do
!           end do

!         end do
!       end do
!     end do

!   end subroutine TensorBoundaryInterp_2D_cpu

!   subroutine TensorBoundaryInterp_2D_gpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nVariables,nElements
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fBound

!     call TensorBoundaryInterp_2D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fBound,this % N,nVariables,nElements)

!   end subroutine TensorBoundaryInterp_2D_gpu

!   subroutine ScalarBoundaryInterp_3D_cpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nVariables,nElements
!     real(prec),intent(in)      :: f(0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out)     :: fBound(0:this % N,0:this % N,1:nVariables,1:6,1:nElements)
!     ! Local
!     integer :: i,j,ii,iVar,iEl
!     real(prec) :: fb(1:6)

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do j = 0,this % N
!           do i = 0,this % N

!             fb(1:6) = 0.0_prec

!             do ii = 0,this % N
!               fb(1) = fb(1) + this % bMatrix (ii,0)*f(i,j,ii,iVar,iEl) ! Bottom
!               fb(2) = fb(2) + this % bMatrix (ii,0)*f(i,ii,j,iVar,iEl) ! South
!               fb(3) = fb(3) + this % bMatrix (ii,1)*f(ii,i,j,iVar,iEl) ! East
!               fb(4) = fb(4) + this % bMatrix (ii,1)*f(i,ii,j,iVar,iEl) ! North
!               fb(5) = fb(5) + this % bMatrix (ii,0)*f(ii,i,j,iVar,iEl) ! West
!               fb(6) = fb(6) + this % bMatrix (ii,1)*f(i,j,ii,iVar,iEl) ! Top
!             end do

!             fBound(i,j,iVar,1:6,iEl) = fb(1:6)

!           end do
!         end do
!       end do
!     end do

!   end subroutine ScalarBoundaryInterp_3D_cpu

!   subroutine ScalarBoundaryInterp_3D_gpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nVariables,nElements
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fBound

!     call ScalarBoundaryInterp_3D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fBound,this % N,nVariables,nElements)

!   end subroutine ScalarBoundaryInterp_3D_gpu

!   subroutine VectorBoundaryInterp_3D_cpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nVariables,nElements
!     real(prec),intent(in)  :: f(1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out)  :: fBound(1:3,0:this % N,0:this % N,1:nVariables,1:6,1:nElements)
!     ! Local
!     integer :: i,j,ii,idir,iVar,iEl
!     real(prec) :: fb(1:3,1:6)

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do j = 0,this % N
!           do i = 0,this % N

!             fb(1:3,1:6) = 0.0_prec
!             do ii = 0,this % N
!               do idir = 1,3
!                 fb(idir,1) = fb(idir,1) + this % bMatrix (ii,0)*f(idir,i,j,ii,iVar,iEl) ! Bottom
!                 fb(idir,2) = fb(idir,2) + this % bMatrix (ii,0)*f(idir,i,ii,j,iVar,iEl) ! South
!                 fb(idir,3) = fb(idir,3) + this % bMatrix (ii,1)*f(idir,ii,i,j,iVar,iEl) ! East
!                 fb(idir,4) = fb(idir,4) + this % bMatrix (ii,1)*f(idir,i,ii,j,iVar,iEl) ! North
!                 fb(idir,5) = fb(idir,5) + this % bMatrix (ii,0)*f(idir,ii,i,j,iVar,iEl) ! West
!                 fb(idir,6) = fb(idir,6) + this % bMatrix (ii,1)*f(idir,i,j,ii,iVar,iEl) ! Top
!               end do
!             end do

!             do idir = 1,3
!               fBound(idir,i,j,iVar,1:6,iEl) = fb(idir,1:6)
!             end do

!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorBoundaryInterp_3D_cpu

!   subroutine VectorBoundaryInterp_3D_gpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nVariables,nElements
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fBound

!     call VectorBoundaryInterp_3D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fBound,this % N,nVariables,nElements)

!   end subroutine VectorBoundaryInterp_3D_gpu

!   subroutine TensorBoundaryInterp_3D_cpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nVariables,nElements
!     real(prec),intent(in)  :: f(1:3,1:3,0:this % N,0:this % N,0:this % N,1:nVariables,1:nElements)
!     real(prec),intent(out)  :: fBound(1:3,1:3,0:this % N,0:this % N,1:nVariables,1:6,1:nElements)
!     ! Local
!     integer :: i,j,ii,idir,jdir,iVar,iEl
!     real(prec) :: fb(1:3,1:3,1:6)

!     do iEl = 1,nElements
!       do iVar = 1,nVariables
!         do j = 0,this % N
!           do i = 0,this % N

!             fb(1:3,1:3,1:6) = 0.0_prec
!             do ii = 0,this % N
!               do jdir = 1,3
!                 do idir = 1,3
!                   fb(idir,jdir,1) = fb(idir,jdir,1) + this % bMatrix (ii,0)*f(idir,jdir,i,j,ii,iVar,iEl) ! Bottom
!                   fb(idir,jdir,2) = fb(idir,jdir,2) + this % bMatrix (ii,0)*f(idir,jdir,i,ii,j,iVar,iEl) ! South
!                   fb(idir,jdir,3) = fb(idir,jdir,3) + this % bMatrix (ii,1)*f(idir,jdir,ii,i,j,iVar,iEl) ! East
!                   fb(idir,jdir,4) = fb(idir,jdir,4) + this % bMatrix (ii,1)*f(idir,jdir,i,ii,j,iVar,iEl) ! North
!                   fb(idir,jdir,5) = fb(idir,jdir,5) + this % bMatrix (ii,0)*f(idir,jdir,ii,i,j,iVar,iEl) ! West
!                   fb(idir,jdir,6) = fb(idir,jdir,6) + this % bMatrix (ii,1)*f(idir,jdir,i,j,ii,iVar,iEl) ! Top
!                 end do
!               end do
!             end do

!             do jdir = 1,3
!               do idir = 1,3
!                 fBound(idir,jdir,i,j,iVar,1:6,iEl) = fb(idir,jdir,1:6)
!               end do
!             end do

!           end do
!         end do
!       end do
!     end do

!   end subroutine TensorBoundaryInterp_3D_cpu

!   subroutine TensorBoundaryInterp_3D_gpu(this,f,fBound,nVariables,nElements)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nVariables,nElements
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fBound

!     call TensorBoundaryInterp_3D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fBound,this % N,nVariables,nElements)

!   end subroutine TensorBoundaryInterp_3D_gpu

! ================================================================================================ !
!
! CalculateBarycentricWeights (PRIVATE)
!
!   A PRIVATE routine that calculates and stores the barycentric weights for the Lagrange
!   data-structure.
!
!   This routine is from Alg. 30 on pg. 75 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  subroutine CalculateBarycentricWeights(this)
    implicit none
    class(Lagrange),intent(inout) :: this
    ! Local
    integer :: i,j
    real(real64) :: bWeights(0:this % N)
    real(real64) :: controlPoints(0:this % N)

    do i = 0,this % N
      bWeights(i) = 1.0_real64
      controlPoints(i) = real(this % controlPoints (i+1),real64)
    end do

    ! Computes the product w_k = w_k*(s_k - s_j), k /= j
    do j = 1,this % N
      do i = 0,j - 1

        bWeights(i) = bWeights(i)*(controlPoints(i) - controlPoints(j))
        bWeights(j) = bWeights(j)*(controlPoints(j) - controlPoints(i))

      end do
    end do

    do j = 0,this % N
      bWeights(j) = 1.0_prec/bWeights(j)
      this % bWeights (j + 1) = real(bWeights(j),prec)
    end do

  end subroutine CalculateBarycentricWeights

! ================================================================================================ !
!
! CalculateInterpolationMatrix (PRIVATE)
!
!   A PRIVATE routine that fills in the interpolation matrix for the Lagrange data structure.
!
!   This function is from Alg. 32 on pg. 76 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  subroutine CalculateInterpolationMatrix(this)
    implicit none
    class(Lagrange),intent(inout) :: this
    ! Local
    integer    :: row,col
    logical    :: rowHasMatch
    real(real64) :: temp1,temp2
    real(real64) :: iMatrix(0:this % M,0:this % N)
    real(real64) :: bWeights(0:this % N)
    real(real64) :: controlPoints(0:this % N)
    real(real64) :: targetPoints(0:this % M)

    do col = 0,this % N
      controlPoints(col) = real(this % controlPoints (col+1),real64)
      bWeights(col) = real(this % bWeights (col+1),real64)
    end do
    do row = 0,this % M
      targetPoints(row) = real(this % targetPoints (row+1),real64)
    end do

    do row = 0,this % M

      rowHasMatch = .false.

      do col = 0,this % N

        iMatrix(row,col) = 0.0_real64

        if (AlmostEqual(targetPoints(row),controlPoints(col))) then
          rowHasMatch = .true.
          iMatrix(row,col) = 1.0_real64
        end if

      end do

      if (.not. (rowHasMatch)) then

        temp1 = 0.0_real64

        do col = 0,this % N
          temp2 = bWeights(col)/ &
                  (targetPoints(row) - &
                   controlPoints(col))
          iMatrix(row,col) = temp2
          temp1 = temp1 + temp2
        end do

        do col = 0,this % N
          iMatrix(row,col) = iMatrix(row,col)/temp1
        end do

      end if

    end do

    do row = 0,this % M
      do col = 0,this % N
        this % iMatrix (col + 1,row + 1) = real(iMatrix(row,col),prec)
      end do
    end do

  end subroutine CalculateInterpolationMatrix

! ================================================================================================ !
!
! CalculateDerivativeMatrix (PRIVATE)
!
!   Calculates and stores the derivative matrix and its transpose.
!   Generates a matrix that can be used to approximate derivatives at the interpolation nodes.
!
!   This function is from Alg. 37 on pg. 82 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  subroutine CalculateDerivativeMatrix(this)
    implicit none
    class(Lagrange),intent(inout) :: this
    ! Local
    integer      :: row,col
    real(real64) :: dmat(0:this % N,0:this % N)
    real(real64) :: dgmat(0:this % N,0:this % N)
    real(real64) :: bWeights(0:this % N)
    real(real64) :: qWeights(0:this % N)
    real(real64) :: controlPoints(0:this % N)

    do row = 0,this % N
      bWeights(row) = real(this % bWeights (row+1),real64)
      qWeights(row) = real(this % qWeights (row+1),real64)
      controlPoints(row) = real(this % controlPoints (row+1),real64)
    end do

    do row = 0,this % N

      dmat(row,row) = 0.0_prec

      do col = 0,this % N

        if (.not. (col == row)) then

          dmat(row,col) = bWeights(col)/ &
                          (bWeights(row)* &
                           (controlPoints(row) - &
                            controlPoints(col)))

          dmat(row,row) = dmat(row,row) - dmat(row,col)

        end if

      end do

    end do

    do row = 0,this % N
      do col = 0,this % N
        dgmat(row,col) = -dmat(col,row)* &
                         qWeights(col)/ &
                         qWeights(row)
      end do
    end do

    do row = 0,this % N
      do col = 0,this % N
        this % dMatrix (row + 1,col + 1) = real(dmat(col,row),prec)
        this % dgMatrix (row + 1,col + 1) = real(dgmat(col,row),prec)
      end do
    end do

  end subroutine CalculateDerivativeMatrix

! ================================================================================================ !
!
! CalculateLagrangePolynomials
!
!   Evaluates each of the 1-D Lagrange interpolating polynomials at a specified point.
!
!   This function is from Alg. 34 on pg. 77 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  function CalculateLagrangePolynomials(this,sE) result(lAtS)
    implicit none
    class(Lagrange) :: this
    real(prec)      :: sE
    real(prec)      :: lAtS(0:this % N)
    ! Local
    integer    :: j
    logical    :: xMatchesNode
    real(real64) :: temp1,temp2
    real(real64) :: sELocal
    real(real64) :: controlPoints(0:this % N)
    real(real64) :: bWeights(0:this % N)
    real(real64) :: lS(0:this % N)

    sELocal = real(sE,real64)
    do j = 0,this % N
      controlPoints(j) = real(this % controlPoints (j+1),real64)
      bWeights(j) = real(this % bWeights (j+1),real64)
    end do

    xMatchesNode = .false.

    do j = 0,this % N

      lS(j) = 0.0_real64
      if (AlmostEqual(sELocal,controlPoints(j))) then
        lS(j) = 1.0_real64
        xMatchesNode = .true.
      end if

    end do

    if (xMatchesNode) then
      do j = 0,this % N
        lAtS(j) = real(lS(j),prec)
      end do
      return
    end if

    temp1 = 0.0_real64

    do j = 0,this % N
      temp2 = bWeights(j)/(sE - controlPoints(j))
      lS(j) = temp2
      temp1 = temp1 + temp2
    end do

    lS = lS/temp1

    do j = 0,this % N
      lAtS(j) = real(lS(j),prec)
    end do

  end function CalculateLagrangePolynomials

  ! subroutine WriteHDF5_Lagrange(this,fileId)
  !   implicit none
  !   class(Lagrange),intent(in) :: this
  !   integer(HID_T),intent(in) :: fileId

  !   call CreateGroup_HDF5(fileId,'/interp')

  !   call WriteArray_HDF5(fileId,'/interp/controlpoints', &
  !                        this % controlPoints)

  !   call WriteArray_HDF5(fileId,'/interp/qweights', &
  !                        this % qWeights)

  !   call WriteArray_HDF5(fileId,'/interp/dgmatrix', &
  !                        this % dgMatrix)

  !   call WriteArray_HDF5(fileId,'/interp/dmatrix', &
  !                        this % dMatrix)

  !   call WriteArray_HDF5(fileId,'/interp/bmatrix', &
  !                        this % bMatrix)

  !   call WriteArray_HDF5(fileId,'/interp/imatrix', &
  !                        this % iMatrix)

  ! end subroutine WriteHDF5_Lagrange

end module SELF_Lagrange
