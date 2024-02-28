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

    generic,public :: ScalarGridInterp_2D => ScalarGridInterp_2D_cpu,ScalarGridInterp_2D_gpu
    procedure,private :: ScalarGridInterp_2D_cpu,ScalarGridInterp_2D_gpu

    GENERIC,PUBLIC :: VectorGridInterp_2D => VectorGridInterp_2D_cpu
    PROCEDURE,PRIVATE :: VectorGridInterp_2D_cpu

    generic,public :: ScalarGridInterp_3D => ScalarGridInterp_3D_cpu,ScalarGridInterp_3D_gpu
    procedure,private :: ScalarGridInterp_3D_cpu,ScalarGridInterp_3D_gpu

    GENERIC,PUBLIC :: VectorGridInterp_3D => VectorGridInterp_3D_cpu
    PROCEDURE,PRIVATE :: VectorGridInterp_3D_cpu

    generic,public :: ScalarBoundaryInterp_1D => ScalarBoundaryInterp_1D_cpu,ScalarBoundaryInterp_1D_gpu
    procedure,private :: ScalarBoundaryInterp_1D_cpu,ScalarBoundaryInterp_1D_gpu

    generic,public :: ScalarBoundaryInterp_2D => ScalarBoundaryInterp_2D_cpu,ScalarBoundaryInterp_2D_gpu
    procedure,private :: ScalarBoundaryInterp_2D_cpu,ScalarBoundaryInterp_2D_gpu

    generic,public :: ScalarBoundaryInterp_3D => ScalarBoundaryInterp_3D_cpu,ScalarBoundaryInterp_3D_gpu
    procedure,private :: ScalarBoundaryInterp_3D_cpu,ScalarBoundaryInterp_3D_gpu

    generic,public :: VectorBoundaryInterp_2D => VectorBoundaryInterp_2D_cpu,VectorBoundaryInterp_2D_gpu
    procedure,private :: VectorBoundaryInterp_2D_cpu,VectorBoundaryInterp_2D_gpu

    generic,public :: VectorBoundaryInterp_3D => VectorBoundaryInterp_3D_cpu,VectorBoundaryInterp_3D_gpu
    procedure,private :: VectorBoundaryInterp_3D_cpu,VectorBoundaryInterp_3D_gpu

    generic,public :: TensorBoundaryInterp_2D => TensorBoundaryInterp_2D_cpu,TensorBoundaryInterp_2D_gpu
    procedure,private :: TensorBoundaryInterp_2D_cpu,TensorBoundaryInterp_2D_gpu

    generic,public :: TensorBoundaryInterp_3D => TensorBoundaryInterp_3D_cpu,TensorBoundaryInterp_3D_gpu
    procedure,private :: TensorBoundaryInterp_3D_cpu,TensorBoundaryInterp_3D_gpu

    generic,public :: Derivative_1D => Derivative_1D_cpu,Derivative_1D_gpu
    procedure,private :: Derivative_1D_cpu,Derivative_1D_gpu

    ! generic,public :: DGDerivative_1D => DGDerivative_1D_cpu,DGDerivative_1D_gpu
    ! procedure,private :: DGDerivative_1D_cpu,DGDerivative_1D_gpu

    generic,public :: ScalarGradient_2D => ScalarGradient_2D_cpu,ScalarGradient_2D_gpu
    procedure,private :: ScalarGradient_2D_cpu,ScalarGradient_2D_gpu

    generic,public :: ScalarGradient_3D => ScalarGradient_3D_cpu,ScalarGradient_3D_gpu
    procedure,private :: ScalarGradient_3D_cpu,ScalarGradient_3D_gpu

    generic,public :: VectorGradient_2D => VectorGradient_2D_cpu,VectorGradient_2D_gpu
    procedure,private :: VectorGradient_2D_cpu,VectorGradient_2D_gpu

    generic,public :: VectorDivergence_2D => VectorDivergence_2D_cpu,VectorDivergence_2D_gpu
    procedure,private :: VectorDivergence_2D_cpu,VectorDivergence_2D_gpu

    generic,public :: TensorDivergence_2D => TensorDivergence_2D_cpu,TensorDivergence_2D_gpu
    procedure,private :: TensorDivergence_2D_cpu,TensorDivergence_2D_gpu

    ! GENERIC,PUBLIC :: VectorDGDivergence_2D => VectorDGDivergence_2D_cpu,VectorDGDivergence_2D_gpu
    ! PROCEDURE,PRIVATE :: VectorDGDivergence_2D_cpu,VectorDGDivergence_2D_gpu

    generic,public :: VectorGradient_3D => VectorGradient_3D_cpu,VectorGradient_3D_gpu
    procedure,private :: VectorGradient_3D_cpu,VectorGradient_3D_gpu

    generic,public :: VectorDivergence_3D => VectorDivergence_3D_cpu,VectorDivergence_3D_gpu
    procedure,private :: VectorDivergence_3D_cpu,VectorDivergence_3D_gpu

    generic,public :: TensorDivergence_3D => TensorDivergence_3D_cpu,TensorDivergence_3D_gpu
    procedure,private :: TensorDivergence_3D_cpu,TensorDivergence_3D_gpu

    ! GENERIC,PUBLIC :: VectorDGDivergence_3D => VectorDGDivergence_3D_cpu,VectorDGDivergence_3D_gpu
    ! PROCEDURE,PRIVATE :: VectorDGDivergence_3D_cpu,VectorDGDivergence_3D_gpu

    !procedure,public :: WriteHDF5 => WriteHDF5_Lagrange
    procedure,private :: CalculateBarycentricWeights
    procedure,private :: CalculateInterpolationMatrix
    procedure,private :: CalculateDerivativeMatrix
    procedure,private :: CalculateLagrangePolynomials

  end type Lagrange

  interface
    subroutine ScalarBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="ScalarBoundaryInterp_2D_gpu_wrapper")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: bMatrix_dev,f_dev,fBound_dev
      integer(c_int),value :: N,nVar,nEl
    end subroutine ScalarBoundaryInterp_2D_gpu_wrapper
  end interface

  interface
    subroutine VectorBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="VectorBoundaryInterp_2D_gpu_wrapper")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: bMatrix_dev,f_dev,fBound_dev
      integer(c_int),value :: N,nVar,nEl
    end subroutine VectorBoundaryInterp_2D_gpu_wrapper
  end interface

  interface
    subroutine VectorBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="VectorBoundaryInterp_3D_gpu_wrapper")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: bMatrix_dev,f_dev,fBound_dev
      integer(c_int),value :: N,nVar,nEl
    end subroutine VectorBoundaryInterp_3D_gpu_wrapper
  end interface

  interface
    subroutine TensorBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="TensorBoundaryInterp_2D_gpu_wrapper")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: bMatrix_dev,f_dev,fBound_dev
      integer(c_int),value :: N,nVar,nEl
    end subroutine TensorBoundaryInterp_2D_gpu_wrapper
  end interface

  interface
    subroutine ScalarBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="ScalarBoundaryInterp_3D_gpu_wrapper")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: bMatrix_dev,f_dev,fBound_dev
      integer(c_int),value :: N,nVar,nEl
    end subroutine ScalarBoundaryInterp_3D_gpu_wrapper
  end interface

  interface
    subroutine TensorBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="TensorBoundaryInterp_3D_gpu_wrapper")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: bMatrix_dev,f_dev,fBound_dev
      integer(c_int),value :: N,nVar,nEl
    end subroutine TensorBoundaryInterp_3D_gpu_wrapper
  end interface

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

      this % targetPoints = UniformPoints(-1.0_prec,1.0_prec,0,M)

    end if

    call this % CalculateBarycentricWeights()
    call this % CalculateInterpolationMatrix()
    call this % CalculateDerivativeMatrix()
    this % bMatrix(1:N + 1,1) = this % CalculateLagrangePolynomials(-1.0_prec)
    this % bMatrix(1:N + 1,2) = this % CalculateLagrangePolynomials(1.0_prec)

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

  subroutine self_hipblas_matrixop_1d(A,f,Af,opArows,opAcols,bcols,handle)
    real(prec),pointer,intent(in) :: A(:,:)
    real(prec),pointer,intent(in) :: f(:,:,:)
    real(prec),pointer,intent(inout) :: Af(:,:,:)
    integer,intent(in) :: opArows,opAcols,bcols
    type(c_ptr),intent(inout) :: handle
    ! Local
    integer(c_int) :: m
    integer(c_int) :: n
    integer(c_int) :: k
    real(c_prec) :: alpha
    integer(c_int) :: lda
    integer(c_int) :: ldb
    integer(c_int) :: ldc
    real(c_prec) :: beta

    m = opArows ! number of rows of A^T
    n = bcols ! number of columns of B
    k = opAcols! number of columns of A^T
    alpha = 1.0_c_prec
    lda = k ! leading dimension of A (matrix)
    ldb = k ! leading dimension of B (f)
    ldc = m ! leading dimension of C (Af)
    beta = 0.0_c_prec
#ifdef DOUBLE_PRECISION
    call hipblasCheck(hipblasDgemm(handle, &
                                   HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                   m,n,k,alpha, &
                                   c_loc(A),lda, &
                                   c_loc(f),ldb, &
                                   beta, &
                                   c_loc(Af),ldc))
#else
    call hipblasCheck(hipblasSgemm(handle, &
                                   HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                   m,n,k,alpha, &
                                   c_loc(A),lda, &
                                   c_loc(f),ldb, &
                                   beta, &
                                   c_loc(Af),ldc))
#endif
  end subroutine self_hipblas_matrixop_1d

  subroutine self_hipblas_matrixop_dim1_2d(A,f,Af,controldegree,targetdegree,nvars,nelems,handle)
    real(prec),pointer,intent(in) :: A(:,:)
    real(prec),pointer,intent(in) :: f(:,:,:,:)
    real(prec),pointer,intent(inout) :: Af(:,:,:,:)
    integer,intent(in) :: controldegree,targetdegree,nvars,nelems
    type(c_ptr),intent(inout) :: handle
    ! Local
    integer(c_int) :: m
    integer(c_int) :: n
    integer(c_int) :: k
    real(c_prec) :: alpha
    integer(c_int) :: lda
    integer(c_int) :: ldb
    integer(c_int) :: ldc
    real(c_prec) :: beta

    m = targetdegree + 1 ! number of rows of A^T
    n = nvars*nelems*(controldegree + 1) ! number of columns of B
    k = controldegree + 1! number of columns of A^T
    alpha = 1.0_c_prec
    lda = k ! leading dimension of A (interoplation matrix)
    ldb = k ! leading dimension of B (f)
    ldc = m ! leading dimension of C (fTarget)
    beta = 0.0_c_prec

#ifdef DOUBLE_PRECISION
    ! First pass interpolates in the first quadrature dimension
    call hipblasCheck(hipblasDgemm(handle, &
                                   HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                   m,n,k,alpha, &
                                   c_loc(A),lda, &
                                   c_loc(f),ldb,beta, &
                                   c_loc(Af),ldc))
#else
    ! First pass interpolates in the first quadrature dimension
    call hipblasCheck(hipblasSgemm(handle, &
                                   HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                   m,n,k,alpha, &
                                   c_loc(A),lda, &
                                   c_loc(f),ldb,beta, &
                                   c_loc(Af),ldc))
#endif

  end subroutine self_hipblas_matrixop_dim1_2d

  subroutine self_hipblas_matrixop_dim2_2d(A,f,Af,beta,controldegree,targetdegree,nvars,nelems,handle)
    real(prec),pointer,intent(in) :: A(:,:)
    real(prec),pointer,intent(in) :: f(:,:,:,:)
    real(prec),pointer,intent(inout) :: Af(:,:,:,:)
    real(c_prec),intent(in) :: beta
    integer,intent(in) :: controldegree,targetdegree,nvars,nelems
    type(c_ptr),intent(inout) :: handle
    ! Local
    integer(c_int) :: m
    integer(c_int) :: n
    real(c_prec) :: alpha
    integer(c_int) :: lda

    integer :: i
    integer(c_int64_t) :: strideA
    integer(c_int) :: incx
    integer(c_int64_t) :: stridex
    integer(c_int) :: incy
    integer(c_int64_t) :: stridey
    integer(c_int) :: batchCount

    m = controldegree + 1 ! number of rows of A
    n = targetdegree + 1 ! number of columns of A
    alpha = 1.0_c_prec
    lda = m ! leading dimension of A
    strideA = 0 ! stride for the batches of A (no stride)
    incx = targetdegree + 1 !
    stridex = (controldegree + 1)*(targetdegree + 1)
    !beta = 0.0_c_prec
    incy = targetdegree + 1
    stridey = (targetdegree + 1)*(targetdegree + 1)
    batchCount = nvars*nelems
    do i = 0,targetdegree
#ifdef DOUBLE_PRECISION
      call hipblasCheck(hipblasDgemvStridedBatched(handle, &
                                                   HIPBLAS_OP_T, &
                                                   m,n,alpha, &
                                                   c_loc(A),lda,strideA, &
                                                   c_loc(f(1 + i,1,1,1)),incx,stridex,beta, &
                                                   c_loc(Af(1 + i,1,1,1)),incy,stridey,batchCount))
#else
      call hipblasCheck(hipblasSgemvStridedBatched(handle, &
                                                   HIPBLAS_OP_T, &
                                                   m,n,alpha, &
                                                   c_loc(A),lda,strideA, &
                                                   c_loc(f(1 + i,1,1,1)),incx,stridex,beta, &
                                                   c_loc(Af(1 + i,1,1,1)),incy,stridey,batchCount))
#endif
    end do

  end subroutine self_hipblas_matrixop_dim2_2d

  subroutine self_hipblas_matrixop_dim1_3d(A,f,Af,controldegree,targetdegree,nvars,nelems,handle)
    real(prec),pointer,intent(in) :: A(:,:)
    real(prec),pointer,intent(in) :: f(:,:,:,:,:)
    real(prec),pointer,intent(inout) :: Af(:,:,:,:,:)
    integer,intent(in) :: controldegree,targetdegree,nvars,nelems
    type(c_ptr),intent(inout) :: handle
    ! Local
    integer(c_int) :: m
    integer(c_int) :: n
    integer(c_int) :: k
    real(c_prec) :: alpha
    integer(c_int) :: lda
    integer(c_int) :: ldb
    integer(c_int) :: ldc
    real(c_prec) :: beta

    m = targetdegree + 1 ! number of rows of A^T
    n = nvars*nelems*(controldegree + 1)*(controldegree + 1) ! number of columns of B
    k = controldegree + 1! number of columns of A^T
    alpha = 1.0_c_prec
    lda = k ! leading dimension of A (interoplation matrix)
    ldb = k ! leading dimension of B (f)
    ldc = m ! leading dimension of C (fTarget)
    beta = 0.0_c_prec

#ifdef DOUBLE_PRECISION
    ! First pass interpolates in the first quadrature dimension
    call hipblasCheck(hipblasDgemm(handle, &
                                   HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                   m,n,k,alpha, &
                                   c_loc(A),lda, &
                                   c_loc(f),ldb,beta, &
                                   c_loc(Af),ldc))
#else
    ! First pass interpolates in the first quadrature dimension
    call hipblasCheck(hipblasSgemm(handle, &
                                   HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                   m,n,k,alpha, &
                                   c_loc(A),lda, &
                                   c_loc(f),ldb,beta, &
                                   c_loc(Af),ldc))
#endif

  end subroutine self_hipblas_matrixop_dim1_3d

  subroutine self_hipblas_matrixop_dim2_3d(A,f,Af,beta,controldegree,targetdegree,nvars,nelems,handle)
    real(prec),pointer,intent(in) :: A(:,:)
    real(prec),pointer,intent(in) :: f(:,:,:,:,:)
    real(prec),pointer,intent(inout) :: Af(:,:,:,:,:)
    real(c_prec),intent(in) :: beta
    integer,intent(in) :: controldegree,targetdegree,nvars,nelems
    type(c_ptr),intent(inout) :: handle
    ! Local
    integer(c_int) :: m
    integer(c_int) :: n
    real(c_prec) :: alpha
    integer(c_int) :: lda

    integer :: i
    integer(c_int64_t) :: strideA
    integer(c_int) :: incx
    integer(c_int64_t) :: stridex
    integer(c_int) :: incy
    integer(c_int64_t) :: stridey
    integer(c_int) :: batchCount

    m = controldegree + 1 ! number of rows of A
    n = targetdegree + 1 ! number of columns of A
    alpha = 1.0_c_prec
    lda = m ! leading dimension of A
    strideA = 0 ! stride for the batches of A (no stride)
    incx = targetdegree + 1 !
    stridex = (controldegree + 1)*(targetdegree + 1)
    !beta = 0.0_c_prec
    incy = targetdegree + 1
    stridey = (targetdegree + 1)*(targetdegree + 1)
    batchCount = (controldegree + 1)*nvars*nelems
    do i = 0,targetdegree
#ifdef DOUBLE_PRECISION
      call hipblasCheck(hipblasDgemvStridedBatched(handle, &
                                                   HIPBLAS_OP_T, &
                                                   m,n,alpha, &
                                                   c_loc(A),lda,strideA, &
                                                   c_loc(f(1 + i,1,1,1,1)),incx,stridex,beta, &
                                                   c_loc(Af(1 + i,1,1,1,1)),incy,stridey,batchCount))
#else
      call hipblasCheck(hipblasSgemvStridedBatched(handle, &
                                                   HIPBLAS_OP_T, &
                                                   m,n,alpha, &
                                                   c_loc(A),lda,strideA, &
                                                   c_loc(f(1 + i,1,1,1,1)),incx,stridex,beta, &
                                                   c_loc(Af(1 + i,1,1,1,1)),incy,stridey,batchCount))
#endif
    end do

  end subroutine self_hipblas_matrixop_dim2_3d

  subroutine self_hipblas_matrixop_dim3_3d(A,f,Af,beta,controldegree,targetdegree,nvars,nelems,handle)
    real(prec),pointer,intent(in) :: A(:,:)
    real(prec),pointer,intent(in) :: f(:,:,:,:,:)
    real(prec),pointer,intent(inout) :: Af(:,:,:,:,:)
    real(c_prec),intent(in) :: beta
    integer,intent(in) :: controldegree,targetdegree,nvars,nelems
    type(c_ptr),intent(inout) :: handle
    ! Local
    integer(c_int) :: m
    integer(c_int) :: n
    real(c_prec) :: alpha
    integer(c_int) :: lda

    integer :: i,j
    integer(c_int64_t) :: strideA
    integer(c_int) :: incx
    integer(c_int64_t) :: stridex
    integer(c_int) :: incy
    integer(c_int64_t) :: stridey
    integer(c_int) :: batchCount

    m = controldegree + 1 ! number of rows of A
    n = targetdegree + 1 ! number of columns of A
    alpha = 1.0_c_prec
    lda = m ! leading dimension of A
    strideA = 0 ! stride for the batches of A (no stride)
    incx = (targetdegree + 1)*(targetdegree + 1) !
    stridex = (controldegree + 1)*(targetdegree + 1)*(targetdegree + 1)
    !beta = 0.0_c_prec
    incy = (targetdegree + 1)*(targetdegree + 1)
    stridey = (targetdegree + 1)*(targetdegree + 1)*(targetdegree + 1)
    batchCount = nvars*nelems
    do j = 0,targetdegree
      do i = 0,targetdegree
#ifdef DOUBLE_PRECISION
        call hipblasCheck(hipblasDgemvStridedBatched(handle, &
                                                     HIPBLAS_OP_T, &
                                                     m,n,alpha, &
                                                     c_loc(A),lda,strideA, &
                                                     c_loc(f(1 + i,1 + j,1,1,1)),incx,stridex,beta, &
                                                     c_loc(Af(1 + i,1 + j,1,1,1)),incy,stridey,batchCount))
#else
        call hipblasCheck(hipblasSgemvStridedBatched(handle, &
                                                     HIPBLAS_OP_T, &
                                                     m,n,alpha, &
                                                     c_loc(A),lda,strideA, &
                                                     c_loc(f(1 + i,1 + j,1,1,1)),incx,stridex,beta, &
                                                     c_loc(Af(1 + i,1 + j,1,1,1)),incy,stridey,batchCount))
#endif
      end do
    end do

  end subroutine self_hipblas_matrixop_dim3_3d

  subroutine ScalarGridInterp_1D_cpu(this,f,fTarget,nvars,nelems)
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
    !! $$ \tilde{f}_{m,iel,ivar} = \sum_{i=0}^N f_{i,iel,ivar} I_{i,m} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in)     :: nvars
    !! The number of variables/functions that are interpolated
    integer,intent(in)     :: nelems
    !! The number of spectral elements in the SEM grid
    real(prec),intent(in)  :: f(1:this % N + 1,1:nelems,1:nvars)
    !! (Input) Array of function values, defined on the control grid
    real(prec),intent(out) :: fTarget(1:this % M + 1,1:nelems,1:nvars)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: iel,ivar,i,ii
    real(prec) :: floc

    do ivar = 1,nvars
      do iel = 1,nelems
        do i = 1,this % M + 1
          floc = 0.0_prec
          do ii = 1,this % N + 1
            floc = floc + this % iMatrix(ii,i)*f(ii,iel,ivar)
          end do
          fTarget(i,iel,ivar) = floc
        end do
      end do
    end do

  end subroutine ScalarGridInterp_1D_cpu

  subroutine ScalarGridInterp_1D_gpu(this,f,fTarget,nvars,nelems,handle)
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
    !! $$ \tilde{f}_{m,iel,ivar} = \sum_{i=0}^N f_{i,iel,ivar} I_{i,m} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in) :: nvars
    !! The number of variables/functions that are interpolated
    integer,intent(in) :: nelems
    !! The number of spectral elements in the SEM grid
    real(prec),pointer,intent(in)  :: f(:,:,:)
    !! (Input) Array of function values, defined on the control grid
    real(prec),pointer,intent(inout) :: fTarget(:,:,:)
    !! (Output) Array of function values, defined on the target grid
    type(c_ptr),intent(inout) :: handle

    call self_hipblas_matrixop_1d(this % iMatrix,f,fTarget,this % M + 1,this % N + 1,nvars*nelems,handle)

  end subroutine ScalarGridInterp_1D_gpu

  subroutine ScalarGridInterp_2D_cpu(this,f,fTarget,nvars,nelems)
    !! Host (CPU) implementation of the ScalarGridInterp_2D interface.
    !! In most cases, you should use the `ScalarGridInterp_2D` generic interface,
    !! rather than calling this routine directly.
    !! Interpolate a scalar-2D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{m,n,iel,ivar} = \sum_{j=0}^N \sum_{i=0}^N f_{i,j,iel,ivar} I_{i,m} I_{j,n} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in)     :: nvars
    !! The number of variables/functions that are interpolated
    integer,intent(in)     :: nelems
    !! The number of spectral elements in the SEM grid
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars)
    !! (Input) Array of function values, defined on the control grid
    real(prec),intent(inout) :: fTarget(1:this % M + 1,1:this % M + 1,1:nelems,1:nvars)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: i,j,ii,jj,iel,ivar
    real(prec) :: fi,fij

    do ivar = 1,nvars
      do iel = 1,nelems
        do j = 1,this % M + 1
          do i = 1,this % M + 1

            fij = 0.0_prec
            do jj = 1,this % N + 1
              fi = 0.0_prec
              do ii = 1,this % N + 1
                fi = fi + f(ii,jj,iel,ivar)*this % iMatrix(ii,i)
              end do
              fij = fij + fi*this % iMatrix(jj,j)
            end do
            fTarget(i,j,iel,ivar) = fij

          end do
        end do
      end do
    end do

  end subroutine ScalarGridInterp_2D_cpu

  subroutine ScalarGridInterp_2D_gpu(this,f,fInt,fTarget,nvars,nelems,handle)
    !! Device (GPU) implementation of the ScalarGridInterp_2D interface.
    !! In most cases, you should use the `ScalarGridInterp_2D` generic interface,
    !! rather than calling this routine directly.
    !! This routine calls hip/SELF_Lagrange.cpp:ScalarGridInterp_2D_gpu_wrapper
    !! Interpolate a scalar-2D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{m,n,iel,ivar} = \sum_{j=0}^N \sum_{i=0}^N f_{i,j,iel,ivar} I_{i,m} I_{j,n} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in) :: nvars
    !! The number of variables/functions that are interpolated
    integer,intent(in) :: nelems
    !! The number of spectral elements in the SEM grid
    real(prec),pointer,intent(in)  :: f(:,:,:,:)
    !! (Input) Array of function values, defined on the control grid
    real(prec),pointer,intent(inout) :: fInt(:,:,:,:)
    !! (Inout) workspace array for handling intermediate values interpolated in one direction
    real(prec),pointer,intent(inout) :: fTarget(:,:,:,:)
    !! (Output) Array of function values, defined on the target grid
    type(c_ptr),intent(inout) :: handle

    call self_hipblas_matrixop_dim1_2d(this % iMatrix,f,fInt,this % N,this % M,nvars,nelems,handle)
    call self_hipblas_matrixop_dim2_2d(this % iMatrix,fInt,fTarget,0.0_c_prec,this % N,this % M,nvars,nelems,handle)

  end subroutine ScalarGridInterp_2D_gpu

  subroutine ScalarGridInterp_3D_cpu(this,f,fTarget,nvars,nelems)
    !! Host (CPU) implementation of the ScalarGridInterp_3D interface.
    !! In most cases, you should use the `ScalarGridInterp_3D` generic interface,
    !! rather than calling this routine directly.
    !! Interpolate a scalar-3D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{m,n,iel,ivar} = \sum_{j=0}^N \sum_{i=0}^N f_{i,j,iel,ivar} I_{i,m} I_{j,n} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in)     :: nvars
    !! The number of variables/functions that are interpolated
    integer,intent(in)     :: nelems
    !! The number of spectral elements in the SEM grid
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars)
    !! (Input) Array of function values, defined on the control grid
    real(prec),intent(inout) :: fTarget(1:this % M + 1,1:this % M + 1,1:this % M + 1,1:nelems,1:nvars)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: i,j,k,ii,jj,kk,iel,ivar
    real(prec) :: fi,fij,fijk

    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this % M + 1
          do j = 1,this % M + 1
            do i = 1,this % M + 1

              fijk = 0.0_prec
              do kk = 1,this % N + 1
                fij = 0.0_prec
                do jj = 1,this % N + 1
                  fi = 0.0_prec
                  do ii = 1,this % N + 1
                    fi = fi + f(ii,jj,kk,iel,ivar)*this % iMatrix(ii,i)
                  end do
                  fij = fij + fi*this % iMatrix(jj,j)
                end do
                fijk = fijk + fij*this % iMatrix(kk,k)
              end do
              fTarget(i,j,k,iel,ivar) = fijk

            end do
          end do
        end do
      end do
    end do

  end subroutine ScalarGridInterp_3D_cpu

  subroutine ScalarGridInterp_3D_gpu(this,f,fInt1,fInt2,fTarget,nvars,nelems,handle)
    !! Device (GPU) implementation of the ScalarGridInterp_3D interface.
    !! In most cases, you should use the `ScalarGridInterp_3D` generic interface,
    !! rather than calling this routine directly.
    !! This routine uses hipblas calls to execute the grid interoplation
    !! Interpolate a scalar-3D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{m,n,iel,ivar} = \sum_{j=0}^N \sum_{i=0}^N f_{i,j,iel,ivar} I_{i,m} I_{j,n} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in) :: nvars
    !! The number of variables/functions that are interpolated
    integer,intent(in) :: nelems
    !! The number of spectral elements in the SEM grid
    real(prec),pointer,intent(in)  :: f(:,:,:,:,:)
    !! (Input) Array of function values, defined on the control grid
    real(prec),pointer,intent(inout) :: fInt1(:,:,:,:,:)
    ! (Input) Array of function values, defined on the control grid
    real(prec),pointer,intent(inout) :: fInt2(:,:,:,:,:)
    !! (Inout) workspace array for handling intermediate values interpolated in one direction
    real(prec),pointer,intent(inout) :: fTarget(:,:,:,:,:)
    !! (Output) Array of function values, defined on the target grid
    type(c_ptr),intent(inout) :: handle

    call self_hipblas_matrixop_dim1_3d(this % iMatrix,f,fInt1,this % N,this % M,nvars,nelems,handle)
    call self_hipblas_matrixop_dim2_3d(this % iMatrix,fInt1,fInt2,0.0_c_prec,this % N,this % M,nvars,nelems,handle)
    call self_hipblas_matrixop_dim3_3d(this % iMatrix,fInt2,fTarget,0.0_c_prec,this % N,this % M,nvars,nelems,handle)

  end subroutine ScalarGridInterp_3D_gpu

  subroutine VectorGridInterp_2D_cpu(this,f,fTarget,nvars,nelems)
    !! Host (CPU) implementation of the ScalarGridInterp_2D interface.
    !! In most cases, you should use the `ScalarGridInterp_2D` generic interface,
    !! rather than calling this routine directly.
    !! Interpolate a scalar-2D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{m,n,iel,ivar} = \sum_{j=0}^N \sum_{i=0}^N f_{i,j,iel,ivar} I_{i,m} I_{j,n} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in)     :: nvars
    !! The number of variables/functions that are interpolated
    integer,intent(in)     :: nelems
    !! The number of spectral elements in the SEM grid
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:2)
    !! (Input) Array of function values, defined on the control grid
    real(prec),intent(inout) :: fTarget(1:this % M + 1,1:this % M + 1,1:nelems,1:nvars,1:2)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: i,j,ii,jj,iel,ivar,idir
    real(prec) :: fi,fij

    do idir = 1,2
      do ivar = 1,nvars
        do iel = 1,nelems
          do j = 1,this % M + 1
            do i = 1,this % M + 1

              fij = 0.0_prec
              do jj = 1,this % N + 1
                fi = 0.0_prec
                do ii = 1,this % N + 1
                  fi = fi + f(ii,jj,iel,ivar,idir)*this % iMatrix(ii,i)
                end do
                fij = fij + fi*this % iMatrix(jj,j)
              end do
              fTarget(i,j,iel,ivar,idir) = fij

            end do
          end do
        end do
      end do
    end do

  end subroutine VectorGridInterp_2D_cpu

  subroutine VectorGridInterp_3D_cpu(this,f,fTarget,nvars,nelems)
    !! Host (CPU) implementation of the ScalarGridInterp_3D interface.
    !! In most cases, you should use the `ScalarGridInterp_3D` generic interface,
    !! rather than calling this routine directly.
    !! Interpolate a scalar-3D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{m,n,iel,ivar} = \sum_{j=0}^N \sum_{i=0}^N f_{i,j,iel,ivar} I_{i,m} I_{j,n} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in)     :: nvars
    !! The number of variables/functions that are interpolated
    integer,intent(in)     :: nelems
    !! The number of spectral elements in the SEM grid
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:3)
    !! (Input) Array of function values, defined on the control grid
    real(prec),intent(inout) :: fTarget(1:this % M + 1,1:this % M + 1,1:this % M + 1,1:nelems,1:nvars,1:3)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: i,j,k,ii,jj,kk,iel,ivar,idir
    real(prec) :: fi,fij,fijk

    do idir = 1,3
      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this % M + 1
            do j = 1,this % M + 1
              do i = 1,this % M + 1

                fijk = 0.0_prec
                do kk = 1,this % N + 1
                  fij = 0.0_prec
                  do jj = 1,this % N + 1
                    fi = 0.0_prec
                    do ii = 1,this % N + 1
                      fi = fi + f(ii,jj,kk,iel,ivar,idir)*this % iMatrix(ii,i)
                    end do
                    fij = fij + fi*this % iMatrix(jj,j)
                  end do
                  fijk = fijk + fij*this % iMatrix(kk,k)
                end do
                fTarget(i,j,k,iel,ivar,idir) = fijk

              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine VectorGridInterp_3D_cpu

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
!     INTEGER        :: nvars, nelems
!     REAL(prec)     :: f(0:interp % N,1:nelems,1:nvars)
!     REAL(prec)     :: derF(0:interp % N,1:nelems,1:nvars)
!
!       CALL interp % Derivative_1D( f, derF, nvars, nelems )
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
!     nvars (in)
!
!     nelems (in)
!
!     derF (out)
!      Array of derivative values at the target interpolation nodes.
!
! ================================================================================================ !

  subroutine Derivative_1D_cpu(this,f,df,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:nelems,1:nvars)
    real(prec),intent(out) :: df(1:this % N + 1,1:nelems,1:nvars)
    ! Local
    integer :: i,ii,iel,ivar
    real(prec) :: dfloc

    do iel = 1,nelems
      do ivar = 1,nvars
        do i = 1,this % N + 1

          dfloc = 0.0_prec
          do ii = 1,this % N + 1
            dfloc = dfloc + this % dMatrix(ii,i)*f(ii,iel,ivar)
          end do
          df(i,iel,ivar) = dfloc

        end do
      end do
    end do

  end subroutine Derivative_1D_cpu

  subroutine Derivative_1D_gpu(this,f,df,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in) :: nvars,nelems
    real(prec),pointer,intent(in)  :: f(:,:,:)
    real(prec),pointer,intent(out) :: df(:,:,:)
    type(c_ptr),intent(inout) :: handle

    call self_hipblas_matrixop_1d(this % dMatrix,f,df,this % N + 1,this % N + 1,nvars*nelems,handle)

  end subroutine Derivative_1D_gpu

  ! subroutine DGDerivative_1D_cpu(this,f,bf,df,nvars,nelems)
  !   implicit none
  !   class(Lagrange),intent(in) :: this
  !   integer,intent(in)     :: nvars,nelems
  !   real(prec),intent(in)  :: f(1:this % N+1,1:nelems,1:nvars)
  !   real(prec),intent(in)  :: bf(1:nvars,1:2,1:nelems)
  !   real(prec),intent(out) :: df(1:this % N+1,1:nelems,1:nvars)
  !   ! Local
  !   integer :: i,ii,iel,ivar

  !   do ivar = 1,nvars
  !     do iel = 1,nelems
  !       do i = 1,this % N+1

  !         ! Interior Derivative Matrix Application
  !         df(i,iel,ivar) = 0.0_prec
  !         do ii = 1,this % N+1
  !           df(i,iel,ivar) = df(i,iel,ivar) + this % dgMatrix (ii,i)*f(ii,iel,ivar)
  !         end do

  !         ! Boundary Contribution
  !         df(i,iel,ivar) = df(i,iel,ivar) + (bf(ivar,2,iel)*this % bMatrix (i,1) + &
  !                                            bf(ivar,1,iel)*this % bMatrix (i,0))/ &
  !                          this % qWeights (i)

  !       end do

  !     end do
  !   end do

  ! end subroutine DGDerivative_1D_cpu

  ! subroutine DGDerivative_1D_gpu(this,f_dev,bf_dev,df_dev,nvars,nelems)
  !   implicit none
  !   class(Lagrange),intent(in) :: this
  !   integer,intent(in) :: nvars,nelems
  !   type(c_ptr),intent(in)  :: f_dev
  !   type(c_ptr),intent(in)  :: bf_dev
  !   type(c_ptr),intent(out) :: df_dev

  !   call DGDerivative_1D_gpu_wrapper(this % dgMatrix, &
  !                                    this % bMatrix, &
  !                                    this % qWeights, &
  !                                    f_dev,bf_dev,df_dev, &
  !                                    this % N, &
  !                                    nvars,nelems)

  ! end subroutine DGDerivative_1D_gpu

! ================================================================================================ !
!
! CalculateGradient_2D
!
!   Calculates the gradient of a 2-D function, represented by a 2-D array of nodal values.
!
!   Given a set of nodal values at the interpolation nodes, the gradient of a function through
!   the interpolation nodes can be estimated by
!
!                       (df/dx)_{a,b} = \sum_{i=0}^N f_{i,b} l'_i(\xi_a),   a,b=0,1,2,...,N
!                       (df/dy)_{a,b} = \sum_{j=0}^N f_{a,j} l'_j(\xi_b),   a,b=0,1,2,...,N
!
!   where l_i(\xi) are the Lagrange interpolating polynomials through the interpolation points.
!   The derivative matrix is D_{a,i} = l'_i(\xi_a) maps an array of nodal values at the interpolation
!   nodes to its estimated derivative. This routine serves as a wrapper to call either the CUDA
!   kernel (if CUDA is enabled) or the CPU version.
!
!   Usage :
!
!     TYPE(Lagrange) :: interp
!     INTEGER        :: nvars, nelems
!     REAL(prec)     :: f(0:interp % N,0:interp % N,1:nelems,1:nvars)
!     REAL(prec)     :: gradF(1:2,0:interp % N,0:interp % N,1:nelems,1:nvars)
!
!       CALL interp % CalculateGradient_2D( f, gradF, nvars, nelems )
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
!     nvars (in)
!
!     nelems (in)
!
!     gradF (out)
!      Array of derivative values at the target interpolation nodes.
!
! ================================================================================================ !
!
  subroutine ScalarGradient_2D_cpu(this,f,gradF,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars)
    real(prec),intent(out) :: gradF(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:2)
    ! Local
    integer    :: i,j,ii,iel,ivar
    real(prec) :: df1,df2

    do ivar = 1,nvars
      do iel = 1,nelems
        do j = 1,this % N + 1
          do i = 1,this % N + 1

            df1 = 0.0_prec
            df2 = 0.0_prec
            do ii = 1,this % N + 1
              df1 = df1 + this % dMatrix(ii,i)*f(ii,j,iel,ivar)
              df2 = df2 + this % dMatrix(ii,j)*f(i,ii,iel,ivar)
            end do
            gradf(i,j,iel,ivar,1) = df1
            gradf(i,j,iel,ivar,2) = df2

          end do
        end do
      end do
    end do

  end subroutine ScalarGradient_2D_cpu

  subroutine ScalarGradient_2D_gpu(this,f,df,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),pointer,intent(in) :: f(:,:,:,:)
    real(prec),pointer,intent(inout) :: df(:,:,:,:,:)
    type(c_ptr),intent(inout) :: handle
    ! local
    real(prec),pointer :: dfloc(:,:,:,:)

    dfloc(1:,1:,1:,1:) => df(1:,1:,1:,1:,1)
    call self_hipblas_matrixop_dim1_2d(this % dMatrix,f,dfloc,this % N,this % N,nvars,nelems,handle)
    dfloc(1:,1:,1:,1:) => df(1:,1:,1:,1:,2)
    call self_hipblas_matrixop_dim2_2d(this % dMatrix,f,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)
    dfloc => null()

  end subroutine ScalarGradient_2D_gpu

  subroutine ScalarGradient_3D_cpu(this,f,gradF,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars)
    real(prec),intent(out) :: gradF(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:3)
    ! Local
    integer    :: i,j,k,ii,iel,ivar
    real(prec) :: df1,df2,df3

    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this % N + 1
          do j = 1,this % N + 1
            do i = 1,this % N + 1

              df1 = 0.0_prec
              df2 = 0.0_prec
              df3 = 0.0_prec
              do ii = 1,this % N + 1
                df1 = df1 + this % dMatrix(ii,i)*f(ii,j,k,iel,ivar)
                df2 = df2 + this % dMatrix(ii,j)*f(i,ii,k,iel,ivar)
                df3 = df3 + this % dMatrix(ii,k)*f(i,j,ii,iel,ivar)
              end do
              gradf(i,j,k,iel,ivar,1) = df1
              gradf(i,j,k,iel,ivar,2) = df2
              gradf(i,j,k,iel,ivar,3) = df3
            end do
          end do
        end do
      end do
    end do

  end subroutine ScalarGradient_3D_cpu

  subroutine ScalarGradient_3D_gpu(this,f,df,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),pointer,intent(in) :: f(:,:,:,:,:)
    real(prec),pointer,intent(inout) :: df(:,:,:,:,:,:)
    type(c_ptr),intent(inout) :: handle
    ! local
    real(prec),pointer :: dfloc(:,:,:,:,:)

    dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,1)
    call self_hipblas_matrixop_dim1_3d(this % dMatrix,f,dfloc,this % N,this % N,nvars,nelems,handle)
    dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,2)
    call self_hipblas_matrixop_dim2_3d(this % dMatrix,f,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)
    dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,3)
    call self_hipblas_matrixop_dim3_3d(this % dMatrix,f,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)
    dfloc => null()

  end subroutine ScalarGradient_3D_gpu

  subroutine VectorGradient_2D_cpu(this,f,gradF,nvars,nelems)
    ! Input : Vector(1:2,...)
    ! Output : Tensor(1:2,1:2,....)
    !          > Tensor(1,1) = d/ds1( Vector(1,...) )
    !          > Tensor(2,1) = d/ds1( Vector(2,...) )
    !          > Tensor(1,2) = d/ds2( Vector(1,...) )
    !          > Tensor(2,2) = d/ds2( Vector(2,...) )
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:2)
    real(prec),intent(out) :: gradF(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:2,1:2)
    ! Local
    integer    :: i,j,ii,idir,iel,ivar
    real(prec) :: dfds1,dfds2

    do idir = 1,2
      do ivar = 1,nvars
        do iel = 1,nelems
          do j = 1,this % N + 1
            do i = 1,this % N + 1

              dfds1 = 0.0_prec
              dfds2 = 0.0_prec
              do ii = 1,this % N + 1
                dfds1 = dfds1 + this % dMatrix(ii,i)*f(ii,j,iel,ivar,idir)
                dfds2 = dfds2 + this % dMatrix(ii,j)*f(i,ii,iel,ivar,idir)
              end do
              gradf(i,j,iel,ivar,idir,1) = dfds1
              gradf(i,j,iel,ivar,idir,2) = dfds2

            end do
          end do
        end do
      end do
    end do

  end subroutine VectorGradient_2D_cpu

  subroutine VectorGradient_2D_gpu(this,f,df,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),pointer,intent(in) :: f(:,:,:,:,:)
    real(prec),pointer,intent(inout) :: df(:,:,:,:,:,:)
    type(c_ptr),intent(inout) :: handle
    ! local
    real(prec),pointer :: dfloc(:,:,:,:)
    real(prec),pointer :: floc(:,:,:,:)
    integer :: idir

    do idir = 1,2
      floc(1:,1:,1:,1:) => f(1:,1:,1:,1:,idir)

      dfloc(1:,1:,1:,1:) => df(1:,1:,1:,1:,idir,1)
      call self_hipblas_matrixop_dim1_2d(this % dMatrix,floc,dfloc,this % N,this % N,nvars,nelems,handle)

      dfloc(1:,1:,1:,1:) => df(1:,1:,1:,1:,idir,2)
  call self_hipblas_matrixop_dim2_2d(this % dMatrix,floc,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)

    end do

    dfloc => null()

  end subroutine VectorGradient_2D_gpu

  subroutine VectorGradient_3D_cpu(this,f,gradF,nvars,nelems)
    ! Input : Vector(1:3,...)
    ! Output : Tensor(1:3,1:3,....)
    !          > Tensor(1,1) = d/ds1( Vector(1,...) )
    !          > Tensor(2,1) = d/ds1( Vector(2,...) )
    !          > Tensor(1,2) = d/ds2( Vector(1,...) )
    !          > Tensor(2,2) = d/ds2( Vector(2,...) )
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:3)
    real(prec),intent(out) :: gradF(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:3,1:3)
    ! Local
    integer    :: i,j,k,ii,idir,iel,ivar
    real(prec) :: dfds1,dfds2,dfds3

    do idir = 1,3
      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this % N + 1
            do j = 1,this % N + 1
              do i = 1,this % N + 1

                dfds1 = 0.0_prec
                dfds2 = 0.0_prec
                dfds3 = 0.0_prec
                do ii = 1,this % N + 1
                  dfds1 = dfds1 + this % dMatrix(ii,i)*f(ii,j,k,iel,ivar,idir)
                  dfds2 = dfds2 + this % dMatrix(ii,j)*f(i,ii,k,iel,ivar,idir)
                  dfds3 = dfds3 + this % dMatrix(ii,k)*f(i,j,ii,iel,ivar,idir)
                end do
                gradf(i,j,k,iel,ivar,idir,1) = dfds1
                gradf(i,j,k,iel,ivar,idir,2) = dfds2
                gradf(i,j,k,iel,ivar,idir,3) = dfds3

              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine VectorGradient_3D_cpu

  subroutine VectorGradient_3D_gpu(this,f,df,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),pointer,intent(in) :: f(:,:,:,:,:,:)
    real(prec),pointer,intent(inout) :: df(:,:,:,:,:,:,:)
    type(c_ptr),intent(inout) :: handle
    ! local
    real(prec),pointer :: dfloc(:,:,:,:,:)
    real(prec),pointer :: floc(:,:,:,:,:)
    integer :: idir

    do idir = 1,3
      floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,idir)

      dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,idir,1)
      call self_hipblas_matrixop_dim1_3d(this % dMatrix,floc,dfloc,this % N,this % N,nvars,nelems,handle)

      dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,idir,2)
  call self_hipblas_matrixop_dim2_3d(this % dMatrix,floc,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)

      dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,idir,3)
  call self_hipblas_matrixop_dim3_3d(this % dMatrix,floc,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)

    end do

    dfloc => null()

  end subroutine VectorGradient_3D_gpu

  subroutine VectorDivergence_2D_cpu(this,f,dF,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:2)
    real(prec),intent(out) :: dF(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars)
    ! Local
    integer    :: i,j,ii,iel,ivar
    real(prec) :: dfLoc

    do ivar = 1,nvars
      do iel = 1,nelems
        do j = 1,this % N + 1
          do i = 1,this % N + 1

            dfLoc = 0.0_prec
            do ii = 1,this % N + 1
              dfLoc = dfLoc + this % dMatrix(ii,i)*f(ii,j,iel,ivar,1)
            end do
            dF(i,j,iel,ivar) = dfLoc

          end do
        end do
      end do
    end do

    do ivar = 1,nvars
      do iel = 1,nelems
        do j = 1,this % N + 1
          do i = 1,this % N + 1

            dfLoc = 0.0_prec
            do ii = 1,this % N + 1
              dfLoc = dfLoc + this % dMatrix(ii,j)*f(i,ii,iel,ivar,2)
            end do
            dF(i,j,iel,ivar) = dF(i,j,iel,ivar) + dfLoc

          end do
        end do
      end do
    end do

  end subroutine VectorDivergence_2D_cpu

  subroutine VectorDivergence_2D_gpu(this,f,df,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),pointer,intent(in) :: f(:,:,:,:,:)
    real(prec),pointer,intent(inout) :: df(:,:,:,:)
    type(c_ptr),intent(inout) :: handle
    ! local
    real(prec),pointer :: floc(:,:,:,:)

    floc(1:,1:,1:,1:) => f(1:,1:,1:,1:,1)
    call self_hipblas_matrixop_dim1_2d(this % dMatrix,floc,df,this % N,this % N,nvars,nelems,handle)
    floc(1:,1:,1:,1:) => f(1:,1:,1:,1:,2)
    call self_hipblas_matrixop_dim2_2d(this % dMatrix,floc,df,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    floc => null()

  end subroutine VectorDivergence_2D_gpu

  subroutine TensorDivergence_2D_cpu(this,f,dF,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:2,1:2)
    real(prec),intent(out) :: dF(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:2)
    ! Local
    integer    :: i,j,ii,iel,ivar,idir
    real(prec) :: dfLoc

    do idir = 1, 2
      do ivar = 1,nvars
        do iel = 1,nelems
          do j = 1,this % N + 1
            do i = 1,this % N + 1

              dfLoc = 0.0_prec
              do ii = 1,this % N + 1
                dfLoc = dfLoc + this % dMatrix(ii,i)*f(ii,j,iel,ivar,idir,1)
              end do
              dF(i,j,iel,ivar,idir) = dfLoc

            end do
          end do
        end do
      end do

      do ivar = 1,nvars
        do iel = 1,nelems
          do j = 1,this % N + 1
            do i = 1,this % N + 1

              dfLoc = 0.0_prec
              do ii = 1,this % N + 1
                dfLoc = dfLoc + this % dMatrix(ii,j)*f(i,ii,iel,ivar,idir,2)
              end do
              dF(i,j,iel,ivar,idir) = dF(i,j,iel,ivar,idir) + dfLoc

            end do
          end do
        end do
      end do
    enddo

  end subroutine TensorDivergence_2D_cpu

  subroutine TensorDivergence_2D_gpu(this,f,df,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),pointer,intent(in) :: f(:,:,:,:,:,:)
    real(prec),pointer,intent(inout) :: df(:,:,:,:,:)
    type(c_ptr),intent(inout) :: handle
    ! local
    real(prec),pointer :: dfloc(:,:,:,:)
    real(prec),pointer :: floc(:,:,:,:)

    dfloc(1:,1:,1:,1:) => df(1:,1:,1:,1:,1)
    floc(1:,1:,1:,1:) => f(1:,1:,1:,1:,1,1)
    call self_hipblas_matrixop_dim1_2d(this % dMatrix,floc,dfloc,this % N,this % N,nvars,nelems,handle)
    floc(1:,1:,1:,1:) => f(1:,1:,1:,1:,1,2)
    call self_hipblas_matrixop_dim2_2d(this % dMatrix,floc,dfloc,1.0_c_prec,this % N,this % N,nvars,nelems,handle)

    dfloc(1:,1:,1:,1:) => df(1:,1:,1:,1:,2)
    floc(1:,1:,1:,1:) => f(1:,1:,1:,1:,2,1)
    call self_hipblas_matrixop_dim1_2d(this % dMatrix,floc,dfloc,this % N,this % N,nvars,nelems,handle)
    floc(1:,1:,1:,1:) => f(1:,1:,1:,1:,2,2)
    call self_hipblas_matrixop_dim2_2d(this % dMatrix,floc,dfloc,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    floc => null()

  end subroutine TensorDivergence_2D_gpu

  subroutine VectorDivergence_3D_cpu(this,f,dF,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:3)
    real(prec),intent(out) :: dF(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars)
    ! Local
    integer    :: i,j,k,ii,iel,ivar
    real(prec) :: dfLoc

    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this % N + 1
          do j = 1,this % N + 1
            do i = 1,this % N + 1

              dfLoc = 0.0_prec
              do ii = 1,this % N + 1
                dfLoc = dfLoc + this % dMatrix(ii,i)*f(ii,j,k,iel,ivar,1)
              end do
              dF(i,j,k,iel,ivar) = dfLoc

            end do
          end do
        end do
      end do
    end do

    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this % N + 1
          do j = 1,this % N + 1
            do i = 1,this % N + 1

              dfLoc = 0.0_prec
              do ii = 1,this % N + 1
                dfLoc = dfLoc + this % dMatrix(ii,j)*f(i,ii,k,iel,ivar,2)
              end do
              dF(i,j,k,iel,ivar) = dF(i,j,k,iel,ivar) + dfLoc

            end do
          end do
        end do
      end do
    end do

    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this % N + 1
          do j = 1,this % N + 1
            do i = 1,this % N + 1

              dfLoc = 0.0_prec
              do ii = 1,this % N + 1
                dfLoc = dfLoc + this % dMatrix(ii,k)*f(i,j,ii,iel,ivar,3)
              end do
              dF(i,j,k,iel,ivar) = dF(i,j,k,iel,ivar) + dfLoc

            end do
          end do
        end do
      end do
    end do

  end subroutine VectorDivergence_3D_cpu

  subroutine VectorDivergence_3D_gpu(this,f,df,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),pointer,intent(in) :: f(:,:,:,:,:,:)
    real(prec),pointer,intent(inout) :: df(:,:,:,:,:)
    type(c_ptr),intent(inout) :: handle
    ! local
    real(prec),pointer :: floc(:,:,:,:,:)

    floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,1)
    call self_hipblas_matrixop_dim1_3d(this % dMatrix,floc,df,this % N,this % N,nvars,nelems,handle)
    floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,2)
    call self_hipblas_matrixop_dim2_3d(this % dMatrix,floc,df,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,3)
    call self_hipblas_matrixop_dim2_3d(this % dMatrix,floc,df,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    floc => null()

  end subroutine VectorDivergence_3D_gpu

  subroutine TensorDivergence_3D_cpu(this,f,dF,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:3,1:3)
    real(prec),intent(out) :: dF(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:3)
    ! Local
    integer    :: i,j,k,ii,iel,ivar,idir
    real(prec) :: dfLoc

    do idir = 1,3
      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this % N + 1
            do j = 1,this % N + 1
              do i = 1,this % N + 1

                dfLoc = 0.0_prec
                do ii = 1,this % N + 1
                  dfLoc = dfLoc + this % dMatrix(ii,i)*f(ii,j,k,iel,ivar,idir,1)
                end do
                dF(i,j,k,iel,ivar,idir) = dfLoc

              end do
            end do
          end do
        end do
      end do

      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this % N + 1
            do j = 1,this % N + 1
              do i = 1,this % N + 1

                dfLoc = 0.0_prec
                do ii = 1,this % N + 1
                  dfLoc = dfLoc + this % dMatrix(ii,j)*f(i,ii,k,iel,ivar,idir,2)
                end do
                dF(i,j,k,iel,ivar,idir) = dF(i,j,k,iel,ivar,idir) + dfLoc

              end do
            end do
          end do
        end do
      end do

      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this % N + 1
            do j = 1,this % N + 1
              do i = 1,this % N + 1

                dfLoc = 0.0_prec
                do ii = 1,this % N + 1
                  dfLoc = dfLoc + this % dMatrix(ii,k)*f(i,j,ii,iel,ivar,idir,3)
                end do
                dF(i,j,k,iel,ivar,idir) = dF(i,j,k,iel,ivar,idir) + dfLoc

              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine TensorDivergence_3D_cpu

  subroutine TensorDivergence_3D_gpu(this,f,df,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),pointer,intent(in) :: f(:,:,:,:,:,:,:)
    real(prec),pointer,intent(inout) :: df(:,:,:,:,:,:)
    type(c_ptr),intent(inout) :: handle
    ! local
    real(prec),pointer :: floc(:,:,:,:,:)
    real(prec),pointer :: dfloc(:,:,:,:,:)
    integer :: idir

    do idir = 1,3
      dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,idir)
      floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,idir,1)
      call self_hipblas_matrixop_dim1_3d(this % dMatrix,floc,dfloc,this % N,this % N,nvars,nelems,handle)
      floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,idir,2)
      call self_hipblas_matrixop_dim2_3d(this % dMatrix,floc,dfloc,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
      floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,idir,3)
      call self_hipblas_matrixop_dim2_3d(this % dMatrix,floc,dfloc,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    enddo

    floc => null()

  end subroutine TensorDivergence_3D_gpu

!   ! /////////////////////////////// !
!   ! Boundary Interpolation Routines !

  subroutine ScalarBoundaryInterp_1D_cpu(this,f,fTarget,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),intent(in)      :: f(1:this % N + 1,1:nelems,1:nvars)
    real(prec),intent(inout)   :: fTarget(1:2,1:nelems,1:nvars)
    ! Local
    integer :: ii,iel,ivar
    real(prec) :: fb(1:2)

    do iel = 1,nelems
      do ivar = 1,nvars
        fb(1:2) = 0.0_prec
        do ii = 1,this % N + 1
          fb(1) = fb(1) + this % bMatrix(ii,1)*f(ii,iel,ivar) ! West
          fb(2) = fb(2) + this % bMatrix(ii,2)*f(ii,iel,ivar) ! East
        end do
        fTarget(1:2,iel,ivar) = fb(1:2)
      end do
    end do

  end subroutine ScalarBoundaryInterp_1D_cpu

  subroutine ScalarBoundaryInterp_1D_gpu(this,f,fTarget,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),pointer,intent(in)      :: f(:,:,:)
    real(prec),pointer,intent(inout)     :: fTarget(:,:,:)
    type(c_ptr),intent(inout) :: handle

    call self_hipblas_matrixop_1d(this % bMatrix,f,fTarget,2,this % N + 1,nvars*nelems,handle)

  end subroutine ScalarBoundaryInterp_1D_gpu

  subroutine ScalarBoundaryInterp_2D_cpu(this,f,fTarget,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),intent(in)      :: f(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars)
    real(prec),intent(out)     :: fTarget(1:this % N + 1,1:4,1:nelems,1:nvars)
    ! Local
    integer :: i,ii,iel,ivar
    real(prec) :: fb(1:4)

    do iel = 1,nelems
      do ivar = 1,nvars
        do i = 1,this % N + 1

          fb(1:4) = 0.0_prec

          do ii = 1,this % N + 1
            fb(1) = fb(1) + this % bMatrix(ii,1)*f(i,ii,iel,ivar) ! South
            fb(2) = fb(2) + this % bMatrix(ii,2)*f(ii,i,iel,ivar) ! East
            fb(3) = fb(3) + this % bMatrix(ii,2)*f(i,ii,iel,ivar) ! North
            fb(4) = fb(4) + this % bMatrix(ii,1)*f(ii,i,iel,ivar) ! West
          end do

          fTarget(i,1:4,iel,ivar) = fb(1:4)

        end do
      end do
    end do

  end subroutine ScalarBoundaryInterp_2D_cpu

  subroutine ScalarBoundaryInterp_2D_gpu(this,f,fTarget,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),pointer,intent(in)  :: f(:,:,:,:)
    real(prec),pointer,intent(inout)  :: fTarget(:,:,:,:)
    type(c_ptr),intent(in) :: handle

    call ScalarBoundaryInterp_2D_gpu_wrapper(c_loc(this % bMatrix), &
                                             c_loc(f),c_loc(fTarget), &
                                             this % N,nvars,nelems)

  end subroutine ScalarBoundaryInterp_2D_gpu

  subroutine ScalarBoundaryInterp_3D_cpu(this,f,fTarget,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),intent(in)      :: f(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars)
    real(prec),intent(out)     :: fTarget(1:this % N + 1,1:this % N + 1,1:6,1:nelems,1:nvars)
    ! Local
    integer :: i,j,ii,iel,ivar
    real(prec) :: fb(1:6)

    do iel = 1,nelems
      do ivar = 1,nvars
        do j = 1,this % N + 1
          do i = 1,this % N + 1

            fb(1:6) = 0.0_prec

            do ii = 1,this % N + 1
              fb(1) = fb(1) + this % bMatrix(ii,1)*f(i,j,ii,iel,ivar) ! Bottom
              fb(2) = fb(2) + this % bMatrix(ii,1)*f(i,ii,j,iel,ivar) ! South
              fb(3) = fb(3) + this % bMatrix(ii,2)*f(ii,i,j,iel,ivar) ! East
              fb(4) = fb(4) + this % bMatrix(ii,2)*f(i,ii,j,iel,ivar) ! North
              fb(5) = fb(5) + this % bMatrix(ii,1)*f(ii,i,j,iel,ivar) ! West
              fb(6) = fb(6) + this % bMatrix(ii,2)*f(i,j,ii,iel,ivar) ! Top
            end do

            fTarget(i,j,1:6,iel,ivar) = fb(1:6)

          end do
        end do
      end do
    end do

  end subroutine ScalarBoundaryInterp_3D_cpu

  subroutine ScalarBoundaryInterp_3D_gpu(this,f,fTarget,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),pointer,intent(in)  :: f(:,:,:,:,:)
    real(prec),pointer,intent(inout)  :: fTarget(:,:,:,:,:)
    type(c_ptr),intent(in) :: handle

    call ScalarBoundaryInterp_3D_gpu_wrapper(c_loc(this % bMatrix), &
                                             c_loc(f),c_loc(fTarget), &
                                             this % N,nvars,nelems)

  end subroutine ScalarBoundaryInterp_3D_gpu

  subroutine VectorBoundaryInterp_2D_cpu(this,f,fTarget,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:2)
    real(prec),intent(out)  :: fTarget(1:this % N + 1,1:4,1:nelems,1:nvars,1:2)
    ! Local
    integer :: i,ii,idir,iel,ivar
    real(prec) :: fb(1:4)

    do idir = 1,2
      do ivar = 1,nvars
        do iel = 1,nelems
          do i = 1,this % N + 1

            fb(1:4) = 0.0_prec
            do ii = 1,this % N + 1
              fb(1) = fb(1) + this % bMatrix(ii,1)*f(i,ii,iel,ivar,idir) ! South
              fb(2) = fb(2) + this % bMatrix(ii,2)*f(ii,i,iel,ivar,idir) ! East
              fb(3) = fb(3) + this % bMatrix(ii,2)*f(i,ii,iel,ivar,idir) ! North
              fb(4) = fb(4) + this % bMatrix(ii,1)*f(ii,i,iel,ivar,idir) ! West
            end do

            fTarget(i,1:4,iel,ivar,idir) = fb(1:4)

          end do
        end do
      end do
    end do

  end subroutine VectorBoundaryInterp_2D_cpu

  subroutine VectorBoundaryInterp_2D_gpu(this,f,fTarget,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),pointer,intent(in)  :: f(:,:,:,:,:)
    real(prec),pointer,intent(inout)  :: fTarget(:,:,:,:,:)
    type(c_ptr),intent(in) :: handle

    call VectorBoundaryInterp_2D_gpu_wrapper(c_loc(this % bMatrix), &
                                             c_loc(f),c_loc(fTarget), &
                                             this % N,nvars,nelems)

  end subroutine VectorBoundaryInterp_2D_gpu

  subroutine VectorBoundaryInterp_3D_cpu(this,f,fTarget,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:3)
    real(prec),intent(out)  :: fTarget(1:this % N + 1,1:this % N + 1,1:6,1:nelems,1:nvars,1:3)
    ! Local
    integer :: i,j,ii,idir,iel,ivar
    real(prec) :: fb(1:6)

    do idir = 1,3
      do ivar = 1,nvars
        do iel = 1,nelems
          do j = 1,this % N + 1
            do i = 1,this % N + 1

              fb(1:6) = 0.0_prec
              do ii = 1,this % N + 1
                fb(1) = fb(1) + this % bMatrix(ii,1)*f(i,j,ii,iel,ivar,idir) ! Bottom
                fb(2) = fb(2) + this % bMatrix(ii,1)*f(i,ii,j,iel,ivar,idir) ! South
                fb(3) = fb(3) + this % bMatrix(ii,2)*f(ii,i,j,iel,ivar,idir) ! East
                fb(4) = fb(4) + this % bMatrix(ii,2)*f(i,ii,j,iel,ivar,idir) ! North
                fb(5) = fb(5) + this % bMatrix(ii,1)*f(ii,i,j,iel,ivar,idir) ! West
                fb(6) = fb(6) + this % bMatrix(ii,2)*f(i,j,ii,iel,ivar,idir) ! Bottom
              end do

              fTarget(i,j,1:6,iel,ivar,idir) = fb(1:6)

            end do
          end do
        end do
      end do
    end do

  end subroutine VectorBoundaryInterp_3D_cpu

  subroutine VectorBoundaryInterp_3D_gpu(this,f,fTarget,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),pointer,intent(in)  :: f(:,:,:,:,:,:)
    real(prec),pointer,intent(inout)  :: fTarget(:,:,:,:,:,:)
    type(c_ptr),intent(in) :: handle

    call VectorBoundaryInterp_3D_gpu_wrapper(c_loc(this % bMatrix), &
                                             c_loc(f),c_loc(fTarget), &
                                             this % N,nvars,nelems)

  end subroutine VectorBoundaryInterp_3D_gpu

  subroutine TensorBoundaryInterp_2D_cpu(this,f,fTarget,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:2,1:2)
    real(prec),intent(out)  :: fTarget(1:this % N + 1,1:4,1:nelems,1:nvars,1:2,1:2)
    ! Local
    integer :: i,ii,idir,jdir,iel,ivar
    real(prec) :: fb(1:4)

    do jdir = 1,2
      do idir = 1,2
        do ivar = 1,nvars
          do iel = 1,nelems
            do i = 1,this % N + 1

              fb(1:4) = 0.0_prec
              do ii = 1,this % N + 1
                fb(1) = fb(1) + this % bMatrix(ii,1)*f(i,ii,iel,ivar,idir,jdir) ! South
                fb(2) = fb(2) + this % bMatrix(ii,2)*f(ii,i,iel,ivar,idir,jdir) ! East
                fb(3) = fb(3) + this % bMatrix(ii,2)*f(i,ii,iel,ivar,idir,jdir) ! North
                fb(4) = fb(4) + this % bMatrix(ii,1)*f(ii,i,iel,ivar,idir,jdir) ! West
              end do

              fTarget(i,1:4,iel,ivar,idir,jdir) = fb(1:4)

            end do
          end do
        end do
      end do
    end do

  end subroutine TensorBoundaryInterp_2D_cpu

  subroutine TensorBoundaryInterp_2D_gpu(this,f,fTarget,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),pointer,intent(in)  :: f(:,:,:,:,:,:)
    real(prec),pointer,intent(inout)  :: fTarget(:,:,:,:,:,:)
    type(c_ptr),intent(in) :: handle

    call TensorBoundaryInterp_2D_gpu_wrapper(c_loc(this % bMatrix), &
                                             c_loc(f),c_loc(fTarget), &
                                             this % N,nvars,nelems)

  end subroutine TensorBoundaryInterp_2D_gpu

  subroutine TensorBoundaryInterp_3D_cpu(this,f,fTarget,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:this % N + 1,1:this % N + 1,1:nelems,1:nvars,1:3,1:3)
    real(prec),intent(out)  :: fTarget(1:this % N + 1,1:this % N + 1,1:6,1:nelems,1:nvars,1:3,1:3)
    ! Local
    integer :: i,j,ii,idir,jdir,iel,ivar
    real(prec) :: fb(1:6)

    do jdir = 1,3
      do idir = 1,3
        do ivar = 1,nvars
          do iel = 1,nelems
            do j = 1,this % N + 1
              do i = 1,this % N + 1

                fb(1:6) = 0.0_prec
                do ii = 1,this % N + 1
                  fb(1) = fb(1) + this % bMatrix(ii,1)*f(i,j,ii,iel,ivar,idir,jdir) ! Bottom
                  fb(2) = fb(2) + this % bMatrix(ii,1)*f(i,ii,j,iel,ivar,idir,jdir) ! South
                  fb(3) = fb(3) + this % bMatrix(ii,2)*f(ii,i,j,iel,ivar,idir,jdir) ! East
                  fb(4) = fb(4) + this % bMatrix(ii,2)*f(i,ii,j,iel,ivar,idir,jdir) ! North
                  fb(5) = fb(5) + this % bMatrix(ii,1)*f(ii,i,j,iel,ivar,idir,jdir) ! West
                  fb(6) = fb(6) + this % bMatrix(ii,2)*f(i,j,ii,iel,ivar,idir,jdir) ! Bottom
                end do

                fTarget(i,j,1:6,iel,ivar,idir,jdir) = fb(1:6)

              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine TensorBoundaryInterp_3D_cpu

  subroutine TensorBoundaryInterp_3D_gpu(this,f,fTarget,nvars,nelems,handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),pointer,intent(in)  :: f(:,:,:,:,:,:,:)
    real(prec),pointer,intent(inout)  :: fTarget(:,:,:,:,:,:,:)
    type(c_ptr),intent(in) :: handle

    call TensorBoundaryInterp_3D_gpu_wrapper(c_loc(this % bMatrix), &
                                             c_loc(f),c_loc(fTarget), &
                                             this % N,nvars,nelems)

  end subroutine TensorBoundaryInterp_3D_gpu

  ! subroutine TensorBoundaryInterp_2D_cpu(this,f,fTarget,nvars,nelems)
  !   implicit none
  !   class(Lagrange),intent(in) :: this
  !   integer,intent(in)  :: nvars,nelems
  !   real(prec),intent(in)  :: f(1:this % N+1,1:this % N+1,1:nelems,1:nvars,1:2,1:2)
  !   real(prec),intent(out)  :: fTarget(1:this % N+1,1:4,1:nelems,1:nvars,1:2,1:2)
  !   ! Local
  !   integer :: i,ii,idir,jdir,iel,ivar
  !   real(prec) :: fb(1:4)
  !   do jdir = 1,2
  !     do idir = 1,2
  !       do ivar = 1,nvars
  !         do iel = 1,nelems
  !           do i = 1,this % N+1

  !             fb(1:4) = 0.0_prec
  !             do ii = 1,this % N+1

  !               fb(1) = fb(1) + this % bMatrix (ii,0)*f(i,ii,iel,ivar,idir,jdir) ! South
  !               fb(2) = fb(2) + this % bMatrix (ii,1)*f(ii,i,iel,ivar,idir,jdir) ! East
  !               fb(3) = fb(3) + this % bMatrix (ii,1)*f(i,ii,iel,ivar,idir,jdir) ! North
  !               fb(4) = fb(4) + this % bMatrix (ii,0)*f(ii,i,iel,ivar,idir,jdir) ! West

  !             end do

  !             fTarget(i,1:4,iel,ivar,idir,jdir) = fb(1:4)
  !           end do
  !         end do
  !       end do
  !     end do
  !   end do

  ! end subroutine TensorBoundaryInterp_2D_cpu

  ! subroutine TensorBoundaryInterp_2D_gpu(this,f,fTarget,nvars,nelems)
  !   implicit none
  !   class(Lagrange),intent(in) :: this
  !   integer,intent(in)  :: nvars,nelems
  !   real(prec), pointer, intent(in)  :: f(:,:,:,:,:,:)
  !   real(prec), pointer, intent(inout)  :: fTarget(:,:,:,:,:,:)

  !   call TensorBoundaryInterp_2D_gpu_wrapper(c_loc(this % bMatrix), &
  !                                            c_loc(f),c_loc(fTarget),&
  !                                            this % N,nvars,nelems)

  ! end subroutine TensorBoundaryInterp_2D_gpu

  ! subroutine VectorBoundaryInterp_3D_cpu(this,f,fTarget,nvars,nelems)
  !   implicit none
  !   class(Lagrange),intent(in) :: this
  !   integer,intent(in)  :: nvars,nelems
  !   real(prec),intent(in)  :: f(1:this % N+1,1:this % N+1,1:this % N+1,1:nelems,1:nvars,1:3)
  !   real(prec),intent(out)  :: fTarget(1:this % N+1,1:this % N+1,1:6,1:nelems,1:nvars,1:3)
  !   ! Local
  !   integer :: i,j,ii,idir,iel,ivar
  !   real(prec) :: fb(1:6)

  !   do idir = 1,3
  !     do ivar = 1,nvars
  !       do iel = 1,nelems
  !         do j = 1,this % N+1
  !           do i = 1,this % N+1

  !             fb(1:6) = 0.0_prec
  !             do ii = 1,this % N+1
  !               fb(1) = fb(1) + this % bMatrix (ii,0)*f(i,j,ii,iel,ivar,idir) ! Bottom
  !               fb(2) = fb(2) + this % bMatrix (ii,0)*f(i,ii,j,iel,ivar,idir) ! South
  !               fb(3) = fb(3) + this % bMatrix (ii,1)*f(ii,i,j,iel,ivar,idir) ! East
  !               fb(4) = fb(4) + this % bMatrix (ii,1)*f(i,ii,j,iel,ivar,idir) ! North
  !               fb(5) = fb(5) + this % bMatrix (ii,0)*f(ii,i,j,iel,ivar,idir) ! West
  !               fb(6) = fb(6) + this % bMatrix (ii,1)*f(i,j,ii,iel,ivar,idir) ! Top
  !             end do
  !             fTarget(i,j,1:6,iel,ivar,idir) = fb(1:6)

  !           end do
  !         end do
  !       end do
  !     end do
  !   end do

  ! end subroutine VectorBoundaryInterp_3D_cpu

  ! subroutine VectorBoundaryInterp_3D_gpu(this,f,fTarget,nvars,nelems)
  !   implicit none
  !   class(Lagrange),intent(in) :: this
  !   integer,intent(in)  :: nvars,nelems
  !   real(prec), pointer, intent(in)  :: f(:,:,:,:,:,:)
  !   real(prec), pointer, intent(inout)  :: fTarget(:,:,:,:,:,:)

  !   call VectorBoundaryInterp_3D_gpu_wrapper(c_loc(this % bMatrix), &
  !                                            c_loc(f),c_loc(fTarget),&
  !                                            this % N,nvars,nelems)

  ! end subroutine VectorBoundaryInterp_3D_gpu

  ! subroutine TensorBoundaryInterp_3D_cpu(this,f,fTarget,nvars,nelems)
  !   implicit none
  !   class(Lagrange),intent(in) :: this
  !   integer,intent(in)  :: nvars,nelems
  !   real(prec),intent(in)  :: f(1:this % N+1,1:this % N+1,1:this % N+1,1:nelems,1:nvars,1:3,1:3)
  !   real(prec),intent(out)  :: fTarget(1:this % N+1,1:this % N+1,1:6,1:nelems,1:nvars,1:3,1:3)
  !   ! Local
  !   integer :: i,j,ii,idir,jdir,iel,ivar
  !   real(prec) :: fb(1:6)

  !   do jdir = 1,3
  !     do idir = 1,3
  !   do iel = 1,nelems
  !     do ivar = 1,nvars
  !       do j = 1,this % N+1
  !         do i = 1,this % N+1

  !           fb(1:6) = 0.0_prec
  !           do ii = 1,this % N+1

  !                 fb(1) = fb(1) + this % bMatrix (ii,0)*f(i,j,ii,iel,ivar,idir,jdir) ! Bottom
  !                 fb(2) = fb(2) + this % bMatrix (ii,0)*f(i,ii,j,iel,ivar,idir,jdir) ! South
  !                 fb(3) = fb(3) + this % bMatrix (ii,1)*f(ii,i,j,iel,ivar,idir,jdir) ! East
  !                 fb(4) = fb(4) + this % bMatrix (ii,1)*f(i,ii,j,iel,ivar,idir,jdir) ! North
  !                 fb(5) = fb(5) + this % bMatrix (ii,0)*f(ii,i,j,iel,ivar,idir,jdir) ! West
  !                 fb(6) = fb(6) + this % bMatrix (ii,1)*f(i,j,ii,iel,ivar,idir,jdir) ! Top

  !           end do

  !           fTarget(i,j,1:6,iel,ivar,idir,jdir) = fb(1:6)
  !             end do
  !           end do

  !         end do
  !       end do
  !     end do
  !   end do

  ! end subroutine TensorBoundaryInterp_3D_cpu

  ! subroutine TensorBoundaryInterp_3D_gpu(this,f,fTarget,nvars,nelems)
  !   implicit none
  !   class(Lagrange),intent(in) :: this
  !   integer,intent(in)  :: nvars,nelems
  !   real(prec), pointer, intent(in)  :: f(:,:,:,:,:,:,:)
  !   real(prec), pointer, intent(inout)  :: fTarget(:,:,:,:,:,:,:)

  !   call TensorBoundaryInterp_3D_gpu_wrapper(c_loc(this % bMatrix), &
  !                                            c_loc(f),c_loc(fTarget), &
  !                                            this % N,nvars,nelems)

  ! end subroutine TensorBoundaryInterp_3D_gpu

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
      controlPoints(i) = real(this % controlPoints(i + 1),real64)
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
      this % bWeights(j + 1) = real(bWeights(j),prec)
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
      controlPoints(col) = real(this % controlPoints(col + 1),real64)
      bWeights(col) = real(this % bWeights(col + 1),real64)
    end do
    do row = 0,this % M
      targetPoints(row) = real(this % targetPoints(row + 1),real64)
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
        this % iMatrix(col + 1,row + 1) = real(iMatrix(row,col),prec)
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
      bWeights(row) = real(this % bWeights(row + 1),real64)
      qWeights(row) = real(this % qWeights(row + 1),real64)
      controlPoints(row) = real(this % controlPoints(row + 1),real64)
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
        this % dMatrix(row + 1,col + 1) = real(dmat(col,row),prec)
        this % dgMatrix(row + 1,col + 1) = real(dgmat(col,row),prec)
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
      controlPoints(j) = real(this % controlPoints(j + 1),real64)
      bWeights(j) = real(this % bWeights(j + 1),real64)
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
