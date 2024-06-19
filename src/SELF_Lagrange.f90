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

module SELF_Lagrange

  use iso_fortran_env
  use iso_c_binding
  ! use hipfort
  ! use hipfort_check
  ! use hipfort_hipmalloc
  ! use hipfort_hipblas

  use SELF_Constants
  use SELF_SupportRoutines
  use SELF_Quadrature
  use SELF_HDF5
  use HDF5

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

    type(c_ptr) :: blas_handle = c_null_ptr
      !! A handle for working with hipblas

    real(prec),allocatable,dimension(:) :: controlPoints
      !! The set of nodes in one dimension where data is known.
      !! To create higher dimension interpolation and differentiation operators, structured grids in two and three
      !! dimensions are created by tensor products of the controlPoints. This design decision implies that all
      !! spectral element methods supported by the Lagrange class have the same polynomial degree in each
      !! computational/spatial dimension. In practice, the controlPoints are the Legendre-Gauss, Legendre-Gauss-Lobatto,
      !! Legendre-Gauss-Radau, Chebyshev-Gauss, Chebyshev-Gauss-Lobatto, or Chebyshev-Gauss-Radau quadrature points over
      !! the domain [-1,1] (computational space). The Init routine for this class restricts controlPoints to one of
      !! these quadrature types or uniform points on [-1,1].

    real(prec),allocatable,dimension(:) :: targetPoints
      !! The set of nodes in one dimension where data is to be interpolated to. To create higher dimension interpolation
      !! and differentiation operators, structured grids in two and three dimensions are created by tensor products of
      !! the targetPoints. In practice, the targetPoints are set to a uniformly distributed set of points between [-1,1]
      !! (computational space) to allow for interpolation from unevenly spaced quadrature points to a plotting grid.

    real(prec),allocatable,dimension(:) :: bWeights
      !! The barycentric weights that are calculated from the controlPoints and used for interpolation.

    real(prec),allocatable,dimension(:) :: qWeights
      !! The quadrature weights for discrete integration. The quadradture weights depend on the type of controlPoints
      !! provided; one of Legendre-Gauss, Legendre-Gauss-Lobatto, Legendre-Gauss-Radau, Chebyshev-Gauss,
      !! Chebyshev-Gauss-Lobatto, Chebyshev-Gauss Radau, or Uniform. If Uniform, the quadrature weights are constant
      !! $$dx = \frac{2.0}{N+1}$$.

    real(prec),allocatable,dimension(:,:) :: iMatrix
      !! The interpolation matrix (transpose) for mapping data from the control grid to the target grid.

    real(prec),allocatable,dimension(:,:) :: dMatrix
      !! The derivative matrix for mapping function nodal values to a nodal values of the derivative estimate. The
      !! dMatrix is based on a strong form of the derivative.

    real(prec),allocatable,dimension(:,:) :: dgMatrix
      !! The derivative matrix for mapping function nodal values to a nodal values of the derivative estimate. The dgMatrix is based
      !! on a weak form of the derivative. It must be used with bMatrix to account for boundary contributions in the weak form.

    real(prec),allocatable,dimension(:,:) :: bMatrix
      !! The boundary interpolation matrix that is used to map a grid of nodal values at the control points to the element boundaries.

  contains

    procedure,public :: Init => Init_Lagrange
    procedure,public :: Free => Free_Lagrange

    procedure,public :: ScalarGridInterp_3D
    procedure,public :: VectorGridInterp_3D
    procedure,public :: ScalarBoundaryInterp_3D
    procedure,public :: VectorBoundaryInterp_3D
    procedure,public :: TensorBoundaryInterp_2D
    procedure,public :: TensorBoundaryInterp_3D
    procedure,public :: ScalarGradient_3D
    procedure,public :: VectorGradient_3D
    procedure,public :: VectorDivergence_3D
    procedure,public :: TensorDivergence_3D
    procedure,public :: TensorDGDivergence_3D
    procedure,public :: VectorDGDivergence_3D

    procedure,public :: WriteHDF5 => WriteHDF5_Lagrange

    procedure,private :: CalculateBarycentricWeights
    procedure,private :: CalculateInterpolationMatrix
    procedure,private :: CalculateDerivativeMatrix
    procedure,private :: CalculateLagrangePolynomials

  endtype Lagrange

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

    this%N = N
    this%M = M
    this%controlNodeType = controlNodeType
    this%targetNodeType = targetNodeType
    allocate(this%controlPoints(1:N+1), &
             this%targetPoints(1:M+1), &
             this%bWeights(1:N+1), &
             this%qWeights(1:N+1), &
             this%iMatrix(1:N+1,1:M+1), &
             this%dMatrix(1:N+1,1:N+1), &
             this%dgMatrix(1:N+1,1:N+1), &
             this%bMatrix(1:N+1,1:2))

    !$omp target enter data map(alloc: this)
    !$omp target enter data map(alloc: this % controlPoints)
    !$omp target enter data map(alloc: this % targetPoints)
    !$omp target enter data map(alloc: this % bWeights)
    !$omp target enter data map(alloc: this % qWeights)
    !$omp target enter data map(alloc: this % iMatrix)
    !$omp target enter data map(alloc: this % dMatrix)
    !$omp target enter data map(alloc: this % dgMatrix)
    !$omp target enter data map(alloc: this % bMatrix)

    if(controlNodeType == GAUSS .or. controlNodeType == GAUSS_LOBATTO) then

      call LegendreQuadrature(N, &
                              this%controlPoints, &
                              this%qWeights, &
                              controlNodeType)

    elseif(controlNodeType == CHEBYSHEV_GAUSS .or. controlNodeType == CHEBYSHEV_GAUSS_LOBATTO) then

      call ChebyshevQuadrature(N, &
                               this%controlPoints, &
                               this%qWeights, &
                               controlNodeType)

    elseif(controlNodeType == UNIFORM) then

      this%controlPoints = UniformPoints(-1.0_prec,1.0_prec,0,N)
      this%qWeights = 2.0_prec/real(N,prec)

    endif

    ! Target Points
    if(targetNodeType == GAUSS .or. targetNodeType == GAUSS_LOBATTO) then

      call LegendreQuadrature(M, &
                              this%targetPoints, &
                              q, &
                              targetNodeType)

    elseif(targetNodeType == UNIFORM) then

      this%targetPoints = UniformPoints(-1.0_prec,1.0_prec,0,M)

    endif

    call this%CalculateBarycentricWeights()
    call this%CalculateInterpolationMatrix()
    call this%CalculateDerivativeMatrix()
    this%bMatrix(1:N+1,1) = this%CalculateLagrangePolynomials(-1.0_prec)
    this%bMatrix(1:N+1,2) = this%CalculateLagrangePolynomials(1.0_prec)

    !$omp target update to(this % controlPoints)
    !$omp target update to(this % targetPoints)
    !$omp target update to(this % bWeights)
    !$omp target update to(this % qWeights)
    !$omp target update to(this % iMatrix)
    !$omp target update to(this % dMatrix)
    !$omp target update to(this % dgMatrix)
    !$omp target update to(this % bMatrix)

  endsubroutine Init_Lagrange

  subroutine Free_Lagrange(this)
    !! Frees all memory (host and device) associated with an instance of the Lagrange class
    implicit none
    class(Lagrange),intent(inout) :: this
    !! Lagrange class instance

    deallocate(this%controlPoints)
    deallocate(this%targetPoints)
    deallocate(this%bWeights)
    deallocate(this%qWeights)
    deallocate(this%iMatrix)
    deallocate(this%dMatrix)
    deallocate(this%dgMatrix)
    deallocate(this%bMatrix)

    !$omp target exit data map(delete: this % controlPoints)
    !$omp target exit data map(delete: this % targetPoints)
    !$omp target exit data map(delete: this % bWeights)
    !$omp target exit data map(delete: this % qWeights)
    !$omp target exit data map(delete: this % iMatrix)
    !$omp target exit data map(delete: this % dMatrix)
    !$omp target exit data map(delete: this % dgMatrix)
    !$omp target exit data map(delete: this % bMatrix)
    !$omp target exit data map(delete: this)

  endsubroutine Free_Lagrange

!   subroutine self_hipblas_matrixop_1d(A,f,Af,opArows,opAcols,bcols,handle)
!     real(prec),pointer,intent(in) :: A(:,:)
!     real(prec),pointer,intent(in) :: f(:,:,:)
!     real(prec),pointer,intent(inout) :: Af(:,:,:)
!     integer,intent(in) :: opArows,opAcols,bcols
!     type(c_ptr),intent(inout) :: handle
!     ! Local
!     integer(c_int) :: m
!     integer(c_int) :: n
!     integer(c_int) :: k
!     real(c_prec) :: alpha
!     integer(c_int) :: lda
!     integer(c_int) :: ldb
!     integer(c_int) :: ldc
!     real(c_prec) :: beta

!     m = opArows ! number of rows of A^T
!     n = bcols ! number of columns of B
!     k = opAcols! number of columns of A^T
!     alpha = 1.0_c_prec
!     lda = k ! leading dimension of A (matrix)
!     ldb = k ! leading dimension of B (f)
!     ldc = m ! leading dimension of C (Af)
!     beta = 0.0_c_prec
! #ifdef DOUBLE_PRECISION
!     call hipblasCheck(hipblasDgemm(handle, &
!                                    HIPBLAS_OP_T,HIPBLAS_OP_N, &
!                                    m,n,k,alpha, &
!                                    c_loc(A),lda, &
!                                    c_loc(f),ldb, &
!                                    beta, &
!                                    c_loc(Af),ldc))
! #else
!     call hipblasCheck(hipblasSgemm(handle, &
!                                    HIPBLAS_OP_T,HIPBLAS_OP_N, &
!                                    m,n,k,alpha, &
!                                    c_loc(A),lda, &
!                                    c_loc(f),ldb, &
!                                    beta, &
!                                    c_loc(Af),ldc))
! #endif
!   end subroutine self_hipblas_matrixop_1d

!   subroutine self_hipblas_matrixop_dim1_2d(A,f,Af,controldegree,targetdegree,nvars,nelems,handle)
!     real(prec),pointer,intent(in) :: A(:,:)
!     real(prec),pointer,intent(in) :: f(:,:,:,:)
!     real(prec),pointer,intent(inout) :: Af(:,:,:,:)
!     integer,intent(in) :: controldegree,targetdegree,nvars,nelems
!     type(c_ptr),intent(inout) :: handle
!     ! Local
!     integer(c_int) :: m
!     integer(c_int) :: n
!     integer(c_int) :: k
!     real(c_prec) :: alpha
!     integer(c_int) :: lda
!     integer(c_int) :: ldb
!     integer(c_int) :: ldc
!     real(c_prec) :: beta

!     m = targetdegree + 1 ! number of rows of A^T
!     n = nvars*nelems*(controldegree + 1) ! number of columns of B
!     k = controldegree + 1! number of columns of A^T
!     alpha = 1.0_c_prec
!     lda = k ! leading dimension of A (interpolation/derivative matrix)
!     ldb = k ! leading dimension of B (f)
!     ldc = m ! leading dimension of C (fTarget)
!     beta = 0.0_c_prec

! #ifdef DOUBLE_PRECISION
!     ! First pass interpolates in the first quadrature dimension
!     call hipblasCheck(hipblasDgemm(handle, &
!                                    HIPBLAS_OP_T,HIPBLAS_OP_N, &
!                                    m,n,k,alpha, &
!                                    c_loc(A),lda, &
!                                    c_loc(f),ldb,beta, &
!                                    c_loc(Af),ldc))
! #else
!     ! First pass interpolates in the first quadrature dimension
!     call hipblasCheck(hipblasSgemm(handle, &
!                                    HIPBLAS_OP_T,HIPBLAS_OP_N, &
!                                    m,n,k,alpha, &
!                                    c_loc(A),lda, &
!                                    c_loc(f),ldb,beta, &
!                                    c_loc(Af),ldc))
! #endif

!   end subroutine self_hipblas_matrixop_dim1_2d

!   subroutine self_hipblas_matrixop_dim2_2d(A,f,Af,beta,controldegree,targetdegree,nvars,nelems,handle)
!     real(prec),pointer,intent(in) :: A(:,:)
!     real(prec),pointer,intent(in) :: f(:,:,:,:)
!     real(prec),pointer,intent(inout) :: Af(:,:,:,:)
!     real(c_prec),intent(in) :: beta
!     integer,intent(in) :: controldegree,targetdegree,nvars,nelems
!     type(c_ptr),intent(inout) :: handle
!     ! Local
!     integer(c_int) :: m
!     integer(c_int) :: n
!     real(c_prec) :: alpha
!     integer(c_int) :: lda

!     integer :: i
!     integer(c_int64_t) :: strideA
!     integer(c_int) :: incx
!     integer(c_int64_t) :: stridex
!     integer(c_int) :: incy
!     integer(c_int64_t) :: stridey
!     integer(c_int) :: batchCount

!     m = controldegree + 1 ! number of rows of A
!     n = targetdegree + 1 ! number of columns of A
!     alpha = 1.0_c_prec
!     lda = m ! leading dimension of A
!     strideA = 0 ! stride for the batches of A (no stride)
!     incx = targetdegree + 1 !
!     stridex = (controldegree + 1)*(targetdegree + 1)
!     incy = targetdegree + 1
!     stridey = (targetdegree + 1)*(targetdegree + 1)
!     batchCount = nvars*nelems
!     do i = 0,targetdegree
! #ifdef DOUBLE_PRECISION
!       call hipblasCheck(hipblasDgemvStridedBatched(handle, &
!                                                    HIPBLAS_OP_T, &
!                                                    m,n,alpha, &
!                                                    c_loc(A),lda,strideA, &
!                                                    c_loc(f(1 + i,1,1,1)),incx,stridex,beta, &
!                                                    c_loc(Af(1 + i,1,1,1)),incy,stridey,batchCount))
! #else
!       call hipblasCheck(hipblasSgemvStridedBatched(handle, &
!                                                    HIPBLAS_OP_T, &
!                                                    m,n,alpha, &
!                                                    c_loc(A),lda,strideA, &
!                                                    c_loc(f(1 + i,1,1,1)),incx,stridex,beta, &
!                                                    c_loc(Af(1 + i,1,1,1)),incy,stridey,batchCount))
! #endif
!     end do

!   end subroutine self_hipblas_matrixop_dim2_2d

!   subroutine self_hipblas_matrixop_dim1_3d(A,f,Af,controldegree,targetdegree,nvars,nelems,handle)
!     real(prec),pointer,intent(in) :: A(:,:)
!     real(prec),pointer,intent(in) :: f(:,:,:,:,:)
!     real(prec),pointer,intent(inout) :: Af(:,:,:,:,:)
!     integer,intent(in) :: controldegree,targetdegree,nvars,nelems
!     type(c_ptr),intent(inout) :: handle
!     ! Local
!     integer(c_int) :: m
!     integer(c_int) :: n
!     integer(c_int) :: k
!     real(c_prec) :: alpha
!     integer(c_int) :: lda
!     integer(c_int) :: ldb
!     integer(c_int) :: ldc
!     real(c_prec) :: beta

!     m = targetdegree + 1 ! number of rows of A^T
!     n = nvars*nelems*(controldegree + 1)*(controldegree + 1) ! number of columns of B
!     k = controldegree + 1! number of columns of A^T
!     alpha = 1.0_c_prec
!     lda = k ! leading dimension of A (interoplation matrix)
!     ldb = k ! leading dimension of B (f)
!     ldc = m ! leading dimension of C (fTarget)
!     beta = 0.0_c_prec

! #ifdef DOUBLE_PRECISION
!     ! First pass interpolates in the first quadrature dimension
!     call hipblasCheck(hipblasDgemm(handle, &
!                                    HIPBLAS_OP_T,HIPBLAS_OP_N, &
!                                    m,n,k,alpha, &
!                                    c_loc(A),lda, &
!                                    c_loc(f),ldb,beta, &
!                                    c_loc(Af),ldc))
! #else
!     ! First pass interpolates in the first quadrature dimension
!     call hipblasCheck(hipblasSgemm(handle, &
!                                    HIPBLAS_OP_T,HIPBLAS_OP_N, &
!                                    m,n,k,alpha, &
!                                    c_loc(A),lda, &
!                                    c_loc(f),ldb,beta, &
!                                    c_loc(Af),ldc))
! #endif

!   end subroutine self_hipblas_matrixop_dim1_3d

!   subroutine self_hipblas_matrixop_dim2_3d(A,f,Af,beta,controldegree,targetdegree,nvars,nelems,handle)
!     real(prec),pointer,intent(in) :: A(:,:)
!     real(prec),pointer,intent(in) :: f(:,:,:,:,:)
!     real(prec),pointer,intent(inout) :: Af(:,:,:,:,:)
!     real(c_prec),intent(in) :: beta
!     integer,intent(in) :: controldegree,targetdegree,nvars,nelems
!     type(c_ptr),intent(inout) :: handle
!     ! Local
!     integer(c_int) :: m
!     integer(c_int) :: n
!     real(c_prec) :: alpha
!     integer(c_int) :: lda

!     integer :: i
!     integer(c_int64_t) :: strideA
!     integer(c_int) :: incx
!     integer(c_int64_t) :: stridex
!     integer(c_int) :: incy
!     integer(c_int64_t) :: stridey
!     integer(c_int) :: batchCount

!     m = controldegree + 1 ! number of rows of A
!     n = targetdegree + 1 ! number of columns of A
!     alpha = 1.0_c_prec
!     lda = m ! leading dimension of A
!     strideA = 0 ! stride for the batches of A (no stride)
!     incx = targetdegree + 1 !
!     stridex = (controldegree + 1)*(targetdegree + 1)
!     !beta = 0.0_c_prec
!     incy = targetdegree + 1
!     stridey = (targetdegree + 1)*(targetdegree + 1)
!     batchCount = (controldegree + 1)*nvars*nelems
!     do i = 0,targetdegree
! #ifdef DOUBLE_PRECISION
!       call hipblasCheck(hipblasDgemvStridedBatched(handle, &
!                                                    HIPBLAS_OP_T, &
!                                                    m,n,alpha, &
!                                                    c_loc(A),lda,strideA, &
!                                                    c_loc(f(1 + i,1,1,1,1)),incx,stridex,beta, &
!                                                    c_loc(Af(1 + i,1,1,1,1)),incy,stridey,batchCount))
! #else
!       call hipblasCheck(hipblasSgemvStridedBatched(handle, &
!                                                    HIPBLAS_OP_T, &
!                                                    m,n,alpha, &
!                                                    c_loc(A),lda,strideA, &
!                                                    c_loc(f(1 + i,1,1,1,1)),incx,stridex,beta, &
!                                                    c_loc(Af(1 + i,1,1,1,1)),incy,stridey,batchCount))
! #endif
!     end do

!   end subroutine self_hipblas_matrixop_dim2_3d

!   subroutine self_hipblas_matrixop_dim3_3d(A,f,Af,beta,controldegree,targetdegree,nvars,nelems,handle)
!     real(prec),pointer,intent(in) :: A(:,:)
!     real(prec),pointer,intent(in) :: f(:,:,:,:,:)
!     real(prec),pointer,intent(inout) :: Af(:,:,:,:,:)
!     real(c_prec),intent(in) :: beta
!     integer,intent(in) :: controldegree,targetdegree,nvars,nelems
!     type(c_ptr),intent(inout) :: handle
!     ! Local
!     integer(c_int) :: m
!     integer(c_int) :: n
!     real(c_prec) :: alpha
!     integer(c_int) :: lda

!     integer :: i,j
!     integer(c_int64_t) :: strideA
!     integer(c_int) :: incx
!     integer(c_int64_t) :: stridex
!     integer(c_int) :: incy
!     integer(c_int64_t) :: stridey
!     integer(c_int) :: batchCount

!     m = controldegree + 1 ! number of rows of A
!     n = targetdegree + 1 ! number of columns of A
!     alpha = 1.0_c_prec
!     lda = m ! leading dimension of A
!     strideA = 0 ! stride for the batches of A (no stride)
!     incx = (targetdegree + 1)*(targetdegree + 1) !
!     stridex = (controldegree + 1)*(targetdegree + 1)*(targetdegree + 1)
!     incy = (targetdegree + 1)*(targetdegree + 1)
!     stridey = (targetdegree + 1)*(targetdegree + 1)*(targetdegree + 1)
!     batchCount = nvars*nelems
!     do j = 0,targetdegree
!       do i = 0,targetdegree
! #ifdef DOUBLE_PRECISION
!         call hipblasCheck(hipblasDgemvStridedBatched(handle, &
!                                                      HIPBLAS_OP_T, &
!                                                      m,n,alpha, &
!                                                      c_loc(A),lda,strideA, &
!                                                      c_loc(f(1 + i,1 + j,1,1,1)),incx,stridex,beta, &
!                                                      c_loc(Af(1 + i,1 + j,1,1,1)),incy,stridey,batchCount))
! #else
!         call hipblasCheck(hipblasSgemvStridedBatched(handle, &
!                                                      HIPBLAS_OP_T, &
!                                                      m,n,alpha, &
!                                                      c_loc(A),lda,strideA, &
!                                                      c_loc(f(1 + i,1 + j,1,1,1)),incx,stridex,beta, &
!                                                      c_loc(Af(1 + i,1 + j,1,1,1)),incy,stridey,batchCount))
! #endif
!       end do
!     end do

!   end subroutine self_hipblas_matrixop_dim3_3d

  subroutine ScalarGridInterp_3D(this,f,fTarget,nvars,nelems)
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
    real(prec),intent(in)  :: f(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars)
    !! (Input) Array of function values, defined on the control grid
    real(prec),intent(out) :: fTarget(1:this%M+1,1:this%M+1,1:this%M+1,1:nelems,1:nvars)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: i,j,k,ii,jj,kk,iel,ivar
    real(prec) :: fi,fij,fijk

    !$omp target map(to:f,this % iMatrix) map(from:fTarget)
    !$omp teams distribute parallel do collapse(5)
    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this%M+1
          do j = 1,this%M+1
            do i = 1,this%M+1

              fijk = 0.0_prec
              do kk = 1,this%N+1
                fij = 0.0_prec
                do jj = 1,this%N+1
                  fi = 0.0_prec
                  do ii = 1,this%N+1
                    fi = fi+f(ii,jj,kk,iel,ivar)*this%iMatrix(ii,i)
                  enddo
                  fij = fij+fi*this%iMatrix(jj,j)
                enddo
                fijk = fijk+fij*this%iMatrix(kk,k)
              enddo
              fTarget(i,j,k,iel,ivar) = fijk

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

    ! call self_hipblas_matrixop_dim1_3d(this % iMatrix,f,fInt1,this % N,this % M,nvars,nelems,handle)
    ! call self_hipblas_matrixop_dim2_3d(this % iMatrix,fInt1,fInt2,0.0_c_prec,this % N,this % M,nvars,nelems,handle)
    ! call self_hipblas_matrixop_dim3_3d(this % iMatrix,fInt2,fTarget,0.0_c_prec,this % N,this % M,nvars,nelems,handle)

  endsubroutine ScalarGridInterp_3D

  subroutine VectorGridInterp_3D(this,f,fTarget,nvars,nelems)
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
    real(prec),intent(in)  :: f(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:3)
    !! (Input) Array of function values, defined on the control grid
    real(prec),intent(out) :: fTarget(1:this%M+1,1:this%M+1,1:this%M+1,1:nelems,1:nvars,1:3)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: i,j,k,ii,jj,kk,iel,ivar,idir
    real(prec) :: fi,fij,fijk

    !$omp target map(to:f,this % iMatrix) map(from:fTarget)
    !$omp teams distribute parallel do collapse(6)
    do idir = 1,3
      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this%M+1
            do j = 1,this%M+1
              do i = 1,this%M+1

                fijk = 0.0_prec
                do kk = 1,this%N+1
                  fij = 0.0_prec
                  do jj = 1,this%N+1
                    fi = 0.0_prec
                    do ii = 1,this%N+1
                      fi = fi+f(ii,jj,kk,iel,ivar,idir)*this%iMatrix(ii,i)
                    enddo
                    fij = fij+fi*this%iMatrix(jj,j)
                  enddo
                  fijk = fijk+fij*this%iMatrix(kk,k)
                enddo
                fTarget(i,j,k,iel,ivar,idir) = fijk

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine VectorGridInterp_3D

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
!     REAL(prec)     :: df(1:2,0:interp % N,0:interp % N,1:nelems,1:nvars)
!
!       CALL interp % CalculateGradient_2D( f, df, nvars, nelems )
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
!     df (out)
!      Array of derivative values at the target interpolation nodes.
!
! ================================================================================================ !
!

  subroutine ScalarGradient_3D(this,f,df,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars)
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:3)
    ! Local
    integer    :: i,j,k,ii,iel,ivar
    real(prec) :: df1,df2,df3

    !$omp target map(to:f,this % dMatrix) map(from:df)
    !$omp teams distribute parallel do collapse(5)
    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

              df1 = 0.0_prec
              df2 = 0.0_prec
              df3 = 0.0_prec
              do ii = 1,this%N+1
                df1 = df1+this%dMatrix(ii,i)*f(ii,j,k,iel,ivar)
                df2 = df2+this%dMatrix(ii,j)*f(i,ii,k,iel,ivar)
                df3 = df3+this%dMatrix(ii,k)*f(i,j,ii,iel,ivar)
              enddo
              df(i,j,k,iel,ivar,1) = df1
              df(i,j,k,iel,ivar,2) = df2
              df(i,j,k,iel,ivar,3) = df3
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

    ! dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,1)
    ! call self_hipblas_matrixop_dim1_3d(this % dMatrix,f,dfloc,this % N,this % N,nvars,nelems,handle)
    ! dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,2)
    ! call self_hipblas_matrixop_dim2_3d(this % dMatrix,f,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)
    ! dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,3)
    ! call self_hipblas_matrixop_dim3_3d(this % dMatrix,f,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)
    ! dfloc => null()

  endsubroutine ScalarGradient_3D

  subroutine VectorGradient_3D(this,f,df,nvars,nelems)
    ! Input : Vector(1:3,...)
    ! Output : Tensor(1:3,1:3,....)
    !          > Tensor(1,1) = d/ds1( Vector(1,...) )
    !          > Tensor(2,1) = d/ds1( Vector(2,...) )
    !          > Tensor(1,2) = d/ds2( Vector(1,...) )
    !          > Tensor(2,2) = d/ds2( Vector(2,...) )
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:3)
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:3,1:3)
    ! Local
    integer    :: i,j,k,ii,idir,iel,ivar
    real(prec) :: dfds1,dfds2,dfds3

    !$omp target map(to:f,this % dMatrix) map(from:df)
    !$omp teams distribute parallel do collapse(6)
    do idir = 1,3
      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this%N+1
            do j = 1,this%N+1
              do i = 1,this%N+1

                dfds1 = 0.0_prec
                dfds2 = 0.0_prec
                dfds3 = 0.0_prec
                do ii = 1,this%N+1
                  dfds1 = dfds1+this%dMatrix(ii,i)*f(ii,j,k,iel,ivar,idir)
                  dfds2 = dfds2+this%dMatrix(ii,j)*f(i,ii,k,iel,ivar,idir)
                  dfds3 = dfds3+this%dMatrix(ii,k)*f(i,j,ii,iel,ivar,idir)
                enddo
                df(i,j,k,iel,ivar,idir,1) = dfds1
                df(i,j,k,iel,ivar,idir,2) = dfds2
                df(i,j,k,iel,ivar,idir,3) = dfds3

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

    !   do idir = 1,3
    !     floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,idir)
    !     dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,idir,1)
    !     call self_hipblas_matrixop_dim1_3d(this % dMatrix,floc,dfloc,this % N,this % N,nvars,nelems,handle)
    !     dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,idir,2)
    !     call self_hipblas_matrixop_dim2_3d(this % dMatrix,floc,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)
    !     dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,idir,3)
    !     call self_hipblas_matrixop_dim3_3d(this % dMatrix,floc,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)
    !   end do
    !   dfloc => null()

  endsubroutine VectorGradient_3D

  subroutine VectorDivergence_3D(this,f,dF,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:3)
    real(prec),intent(out) :: dF(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars)
    ! Local
    integer    :: i,j,k,ii,iel,ivar
    real(prec) :: dfLoc

    !$omp target map(to:f,this % dMatrix) map(from:df)
    !$omp teams
    !$omp distribute parallel do collapse(5)
    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

              dfLoc = 0.0_prec
              do ii = 1,this%N+1
                dfLoc = dfLoc+this%dMatrix(ii,i)*f(ii,j,k,iel,ivar,1)
              enddo
              dF(i,j,k,iel,ivar) = dfLoc

            enddo
          enddo
        enddo
      enddo
    enddo

    !$omp distribute parallel do collapse(5)
    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

              dfLoc = 0.0_prec
              do ii = 1,this%N+1
                dfLoc = dfLoc+this%dMatrix(ii,j)*f(i,ii,k,iel,ivar,2)
              enddo
              dF(i,j,k,iel,ivar) = dF(i,j,k,iel,ivar)+dfLoc

            enddo
          enddo
        enddo
      enddo
    enddo

    !$omp distribute parallel do collapse(5)
    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

              dfLoc = 0.0_prec
              do ii = 1,this%N+1
                dfLoc = dfLoc+this%dMatrix(ii,k)*f(i,j,ii,iel,ivar,3)
              enddo
              dF(i,j,k,iel,ivar) = dF(i,j,k,iel,ivar)+dfLoc

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target

    ! ! local
    ! real(prec),pointer :: floc(:,:,:,:,:)

    ! floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,1)
    ! call self_hipblas_matrixop_dim1_3d(this % dMatrix,floc,df,this % N,this % N,nvars,nelems,handle)
    ! floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,2)
    ! call self_hipblas_matrixop_dim2_3d(this % dMatrix,floc,df,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    ! floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,3)
    ! call self_hipblas_matrixop_dim2_3d(this % dMatrix,floc,df,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    ! floc => null()

  endsubroutine VectorDivergence_3D

  subroutine VectorDGDivergence_3D(this,f,bf,dF,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:3)
    real(prec),intent(in)  :: bf(1:this%N+1,1:this%N+1,1:6,1:nelems,1:nvars)
    real(prec),intent(out) :: dF(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars)
    ! Local
    integer    :: i,j,k,ii,iel,ivar
    real(prec) :: dfLoc

    !$omp target map(to:f,bf,this % dgMatrix,this % bMatrix,this % qWeights) map(from:df)
    !$omp teams
    !$omp distribute parallel do collapse(5)
    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

              dfLoc = 0.0_prec
              do ii = 1,this%N+1
                dfLoc = dfLoc+this%dgMatrix(ii,i)*f(ii,j,k,iel,ivar,1)
              enddo
              dF(i,j,k,iel,ivar) = dfLoc+(this%bMatrix(i,2)*bf(j,k,3,iel,ivar)+ & ! east
                                          this%bMatrix(i,1)*bf(j,k,5,iel,ivar))/ & ! west
                                   this%qweights(i)

            enddo
          enddo
        enddo
      enddo
    enddo

    !$omp distribute parallel do collapse(5)
    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

              dfLoc = 0.0_prec
              do ii = 1,this%N+1
                dfLoc = dfLoc+this%dgMatrix(ii,j)*f(i,ii,k,iel,ivar,2)
              enddo
              dF(i,j,k,iel,ivar) = dF(i,j,k,iel,ivar)+dfLoc+ &
                                   (this%bMatrix(j,2)*bf(i,k,4,iel,ivar)+ & ! north
                                    this%bMatrix(j,1)*bf(i,k,2,iel,ivar))/ & ! south
                                   this%qweights(j)

            enddo
          enddo
        enddo
      enddo
    enddo

    !$omp distribute parallel do collapse(5)
    do ivar = 1,nvars
      do iel = 1,nelems
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

              dfLoc = 0.0_prec
              do ii = 1,this%N+1
                dfLoc = dfLoc+this%dgMatrix(ii,k)*f(i,j,ii,iel,ivar,3)
              enddo
              dF(i,j,k,iel,ivar) = dF(i,j,k,iel,ivar)+dfLoc+ &
                                   (this%bMatrix(k,2)*bf(i,j,6,iel,ivar)+ & ! top
                                    this%bMatrix(k,1)*bf(i,j,1,iel,ivar))/ & ! bottom
                                   this%qweights(k)

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target

    !     ! local
    ! real(prec),pointer :: floc(:,:,:,:,:)

    ! floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,1)
    ! call self_hipblas_matrixop_dim1_3d(this % dgMatrix,floc,df,this % N,this % N,nvars,nelems,handle)
    ! floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,2)
    ! call self_hipblas_matrixop_dim2_3d(this % dgMatrix,floc,df,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    ! floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,3)
    ! call self_hipblas_matrixop_dim2_3d(this % dgMatrix,floc,df,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    ! floc => null()

    ! ! Add the boundary contributions
    ! call VectorDGDivergence_BoundaryContribution_3D(c_loc(this % bMatrix),&
    !                                                 c_loc(this % qWeights),&
    !                                                 c_loc(bf), c_loc(df),&
    !                                                 this % N, nvars, nelems)
  endsubroutine VectorDGDivergence_3D

  subroutine TensorDivergence_3D(this,f,dF,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:3,1:3)
    real(prec),intent(out) :: dF(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:3)
    ! Local
    integer    :: i,j,k,ii,iel,ivar,idir
    real(prec) :: dfLoc

    !$omp target map(to:f,this % dMatrix) map(from:df)
    !$omp teams
    do idir = 1,3
      !$omp distribute parallel do collapse(5)
      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this%N+1
            do j = 1,this%N+1
              do i = 1,this%N+1

                dfLoc = 0.0_prec
                do ii = 1,this%N+1
                  dfLoc = dfLoc+this%dMatrix(ii,i)*f(ii,j,k,iel,ivar,idir,1)
                enddo
                dF(i,j,k,iel,ivar,idir) = dfLoc

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp distribute parallel do collapse(5)
      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this%N+1
            do j = 1,this%N+1
              do i = 1,this%N+1

                dfLoc = 0.0_prec
                do ii = 1,this%N+1
                  dfLoc = dfLoc+this%dMatrix(ii,j)*f(i,ii,k,iel,ivar,idir,2)
                enddo
                dF(i,j,k,iel,ivar,idir) = dF(i,j,k,iel,ivar,idir)+dfLoc

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp distribute parallel do collapse(5)
      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this%N+1
            do j = 1,this%N+1
              do i = 1,this%N+1

                dfLoc = 0.0_prec
                do ii = 1,this%N+1
                  dfLoc = dfLoc+this%dMatrix(ii,k)*f(i,j,ii,iel,ivar,idir,3)
                enddo
                dF(i,j,k,iel,ivar,idir) = dF(i,j,k,iel,ivar,idir)+dfLoc

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target
    ! ! local
    ! integer :: idir
    ! real(prec),pointer :: dfloc(:,:,:,:,:)
    ! real(prec),pointer :: floc(:,:,:,:,:)

    ! do idir = 1,3
    !   dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,idir)   ! Set the gradient direction pointer

    !   floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,idir,1)   ! Set the interior pointer
    !   call self_hipblas_matrixop_dim1_3d(this % dMatrix,floc,dfloc,this % N,this % N,nvars,nelems,handle)

    !   floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,idir,2)   ! Set the interior pointer
    !   call self_hipblas_matrixop_dim2_3d(this % dMatrix,floc,dfloc,1.0_c_prec,this % N,this % N,nvars,nelems,handle)

    !   floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,idir,3)   ! Set the interior pointer
    !   call self_hipblas_matrixop_dim3_3d(this % dMatrix,floc,dfloc,1.0_c_prec,this % N,this % N,nvars,nelems,handle)

    ! enddo

    ! nullify(floc)
    ! nullify(dfloc)
  endsubroutine TensorDivergence_3D

  subroutine TensorDGDivergence_3D(this,f,bf,dF,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:3,1:3)
    real(prec),intent(in)  :: bf(1:this%N+1,1:this%N+1,1:6,1:nelems,1:nvars,1:3,1:3)
    real(prec),intent(out) :: dF(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:3)
    ! Local
    integer    :: i,j,k,ii,iel,ivar,idir
    real(prec) :: dfLoc

    !$omp target map(to:f,bf,this % dgMatrix,this % bMatrix,this % qWeights) map(from:df)
    !$omp teams
    do idir = 1,3
      !$omp distribute parallel do collapse(5)
      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this%N+1
            do j = 1,this%N+1
              do i = 1,this%N+1

                dfLoc = 0.0_prec
                do ii = 1,this%N+1
                  dfLoc = dfLoc+this%dgMatrix(ii,i)*f(ii,j,k,iel,ivar,idir,1)
                enddo
                dF(i,j,k,iel,ivar,idir) = dfLoc+(this%bMatrix(i,2)*bf(j,k,3,iel,ivar,idir,1)+ & ! east
                                                 this%bMatrix(i,1)*bf(j,k,5,iel,ivar,idir,1))/ & ! west
                                          this%qweights(i)

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp distribute parallel do collapse(5)
      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this%N+1
            do j = 1,this%N+1
              do i = 1,this%N+1

                dfLoc = 0.0_prec
                do ii = 1,this%N+1
                  dfLoc = dfLoc+this%dgMatrix(ii,j)*f(i,ii,k,iel,ivar,idir,2)
                enddo
                dF(i,j,k,iel,ivar,idir) = dF(i,j,k,iel,ivar,idir)+dfLoc+ &
                                          (this%bMatrix(j,2)*bf(i,k,4,iel,ivar,idir,2)+ & ! north
                                           this%bMatrix(j,1)*bf(i,k,2,iel,ivar,idir,2))/ & ! south
                                          this%qweights(j)

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp distribute parallel do collapse(5)
      do ivar = 1,nvars
        do iel = 1,nelems
          do k = 1,this%N+1
            do j = 1,this%N+1
              do i = 1,this%N+1

                dfLoc = 0.0_prec
                do ii = 1,this%N+1
                  dfLoc = dfLoc+this%dgMatrix(ii,k)*f(i,j,ii,iel,ivar,idir,3)
                enddo
                dF(i,j,k,iel,ivar,idir) = dF(i,j,k,iel,ivar,idir)+dfLoc+ &
                                          (this%bMatrix(k,2)*bf(i,j,6,iel,ivar,idir,3)+ & ! top
                                           this%bMatrix(k,1)*bf(i,j,1,iel,ivar,idir,3))/ & ! bottom
                                          this%qweights(k)

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target
    ! local
    ! real(prec),pointer :: bfloc(:,:,:,:,:)
    ! real(prec),pointer :: dfloc(:,:,:,:,:)
    ! real(prec),pointer :: floc(:,:,:,:,:)

    ! do idir = 1,3
    !   dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,idir)   ! Set the gradient direction pointer

    !   floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,idir,1)   ! Set the interior pointer
    !   bfloc(1:,1:,1:,1:,1:) => bf(1:,1:,1:,1:,1:,idir,1) ! Set the boundary pointer
    !   call self_hipblas_matrixop_dim1_3d(this % dgMatrix,floc,dfloc,this % N,this % N,nvars,nelems,handle)
    !   call DGDivergence_BoundaryContribution_dim1_3D_gpu(c_loc(this % bMatrix),&
    !                                                     c_loc(this % qWeights),&
    !                                                     c_loc(bfloc), c_loc(dfloc),&
    !                                                     this % N, nvars, nelems)

    !   floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,idir,2)   ! Set the interior pointer
    !   bfloc(1:,1:,1:,1:,1:) => bf(1:,1:,1:,1:,1:,idir,2) ! Set the boundary pointer
    !   call self_hipblas_matrixop_dim2_3d(this % dgMatrix,floc,dfloc,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    !   call DGDivergence_BoundaryContribution_dim2_3D_gpu(c_loc(this % bMatrix),&
    !                                                     c_loc(this % qWeights),&
    !                                                     c_loc(bfloc), c_loc(dfloc),&
    !                                                     this % N, nvars, nelems)

    !   floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,idir,3)   ! Set the interior pointer
    !   bfloc(1:,1:,1:,1:,1:) => bf(1:,1:,1:,1:,1:,idir,3) ! Set the boundary pointer
    !   call self_hipblas_matrixop_dim3_3d(this % dgMatrix,floc,dfloc,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    !   call DGDivergence_BoundaryContribution_dim3_3D_gpu(c_loc(this % bMatrix),&
    !                                                     c_loc(this % qWeights),&
    !                                                     c_loc(bfloc), c_loc(dfloc),&
    !                                                     this % N, nvars, nelems)

    ! enddo

    ! nullify(floc)
    ! nullify(dfloc)
    ! nullify(bfloc)

  endsubroutine TensorDGDivergence_3D

!   ! /////////////////////////////// !
!   ! Boundary Interpolation Routines !


  subroutine ScalarBoundaryInterp_3D(this,f,fTarget,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),intent(in)      :: f(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars)
    real(prec),intent(out)     :: fTarget(1:this%N+1,1:this%N+1,1:6,1:nelems,1:nvars)
    ! Local
    integer :: i,j,ii,iel,ivar
    real(prec) :: fb(1:6)

    !$omp target map(to:f,this % bMatrix) map(from:fTarget)
    !$omp teams distribute parallel do collapse(4)
    do iel = 1,nelems
      do ivar = 1,nvars
        do j = 1,this%N+1
          do i = 1,this%N+1

            fb(1:6) = 0.0_prec

            do ii = 1,this%N+1
              fb(1) = fb(1)+this%bMatrix(ii,1)*f(i,j,ii,iel,ivar) ! Bottom
              fb(2) = fb(2)+this%bMatrix(ii,1)*f(i,ii,j,iel,ivar) ! South
              fb(3) = fb(3)+this%bMatrix(ii,2)*f(ii,i,j,iel,ivar) ! East
              fb(4) = fb(4)+this%bMatrix(ii,2)*f(i,ii,j,iel,ivar) ! North
              fb(5) = fb(5)+this%bMatrix(ii,1)*f(ii,i,j,iel,ivar) ! West
              fb(6) = fb(6)+this%bMatrix(ii,2)*f(i,j,ii,iel,ivar) ! Top
            enddo

            fTarget(i,j,1:6,iel,ivar) = fb(1:6)

          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine ScalarBoundaryInterp_3D

  subroutine VectorBoundaryInterp_3D(this,f,fTarget,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),intent(in)  :: f(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:3)
    real(prec),intent(out)  :: fTarget(1:this%N+1,1:this%N+1,1:6,1:nelems,1:nvars,1:3)
    ! Local
    integer :: i,j,ii,idir,iel,ivar
    real(prec) :: fb(1:6)

    !$omp target map(to:f,this % bMatrix) map(from:fTarget)
    !$omp teams distribute parallel do collapse(5)
    do idir = 1,3
      do ivar = 1,nvars
        do iel = 1,nelems
          do j = 1,this%N+1
            do i = 1,this%N+1

              fb(1:6) = 0.0_prec
              do ii = 1,this%N+1
                fb(1) = fb(1)+this%bMatrix(ii,1)*f(i,j,ii,iel,ivar,idir) ! Bottom
                fb(2) = fb(2)+this%bMatrix(ii,1)*f(i,ii,j,iel,ivar,idir) ! South
                fb(3) = fb(3)+this%bMatrix(ii,2)*f(ii,i,j,iel,ivar,idir) ! East
                fb(4) = fb(4)+this%bMatrix(ii,2)*f(i,ii,j,iel,ivar,idir) ! North
                fb(5) = fb(5)+this%bMatrix(ii,1)*f(ii,i,j,iel,ivar,idir) ! West
                fb(6) = fb(6)+this%bMatrix(ii,2)*f(i,j,ii,iel,ivar,idir) ! Bottom
              enddo

              fTarget(i,j,1:6,iel,ivar,idir) = fb(1:6)

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine VectorBoundaryInterp_3D

  subroutine TensorBoundaryInterp_2D(this,f,fTarget,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),intent(in)  :: f(1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:2,1:2)
    real(prec),intent(out)  :: fTarget(1:this%N+1,1:4,1:nelems,1:nvars,1:2,1:2)
    ! Local
    integer :: i,ii,idir,jdir,iel,ivar
    real(prec) :: fb(1:4)

    !$omp target map(to:f,this % bMatrix) map(from:fTarget)
    !$omp teams distribute parallel do collapse(5)
    do jdir = 1,2
      do idir = 1,2
        do ivar = 1,nvars
          do iel = 1,nelems
            do i = 1,this%N+1

              fb(1:4) = 0.0_prec
              do ii = 1,this%N+1
                fb(1) = fb(1)+this%bMatrix(ii,1)*f(i,ii,iel,ivar,idir,jdir) ! South
                fb(2) = fb(2)+this%bMatrix(ii,2)*f(ii,i,iel,ivar,idir,jdir) ! East
                fb(3) = fb(3)+this%bMatrix(ii,2)*f(i,ii,iel,ivar,idir,jdir) ! North
                fb(4) = fb(4)+this%bMatrix(ii,1)*f(ii,i,iel,ivar,idir,jdir) ! West
              enddo

              fTarget(i,1:4,iel,ivar,idir,jdir) = fb(1:4)

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine TensorBoundaryInterp_2D

  subroutine TensorBoundaryInterp_3D(this,f,fTarget,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),intent(in)  :: f(1:this%N+1,1:this%N+1,1:this%N+1,1:nelems,1:nvars,1:3,1:3)
    real(prec),intent(out)  :: fTarget(1:this%N+1,1:this%N+1,1:6,1:nelems,1:nvars,1:3,1:3)
    ! Local
    integer :: i,j,ii,idir,jdir,iel,ivar
    real(prec) :: fb(1:6)

    !$omp target map(to:f,this % bMatrix) map(from:fTarget)
    !$omp teams distribute parallel do collapse(6)
    do jdir = 1,3
      do idir = 1,3
        do ivar = 1,nvars
          do iel = 1,nelems
            do j = 1,this%N+1
              do i = 1,this%N+1

                fb(1:6) = 0.0_prec
                do ii = 1,this%N+1
                  fb(1) = fb(1)+this%bMatrix(ii,1)*f(i,j,ii,iel,ivar,idir,jdir) ! Bottom
                  fb(2) = fb(2)+this%bMatrix(ii,1)*f(i,ii,j,iel,ivar,idir,jdir) ! South
                  fb(3) = fb(3)+this%bMatrix(ii,2)*f(ii,i,j,iel,ivar,idir,jdir) ! East
                  fb(4) = fb(4)+this%bMatrix(ii,2)*f(i,ii,j,iel,ivar,idir,jdir) ! North
                  fb(5) = fb(5)+this%bMatrix(ii,1)*f(ii,i,j,iel,ivar,idir,jdir) ! West
                  fb(6) = fb(6)+this%bMatrix(ii,2)*f(i,j,ii,iel,ivar,idir,jdir) ! Bottom
                enddo

                fTarget(i,j,1:6,iel,ivar,idir,jdir) = fb(1:6)

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine TensorBoundaryInterp_3D

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
    real(real64) :: bWeights(0:this%N)
    real(real64) :: controlPoints(0:this%N)

    do i = 0,this%N
      bWeights(i) = 1.0_real64
      controlPoints(i) = real(this%controlPoints(i+1),real64)
    enddo

    ! Computes the product w_k = w_k*(s_k - s_j), k /= j
    do j = 1,this%N
      do i = 0,j-1

        bWeights(i) = bWeights(i)*(controlPoints(i)-controlPoints(j))
        bWeights(j) = bWeights(j)*(controlPoints(j)-controlPoints(i))

      enddo
    enddo

    do j = 0,this%N
      bWeights(j) = 1.0_prec/bWeights(j)
      this%bWeights(j+1) = real(bWeights(j),prec)
    enddo

  endsubroutine CalculateBarycentricWeights

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
    real(real64) :: iMatrix(0:this%M,0:this%N)
    real(real64) :: bWeights(0:this%N)
    real(real64) :: controlPoints(0:this%N)
    real(real64) :: targetPoints(0:this%M)

    do col = 0,this%N
      controlPoints(col) = real(this%controlPoints(col+1),real64)
      bWeights(col) = real(this%bWeights(col+1),real64)
    enddo
    do row = 0,this%M
      targetPoints(row) = real(this%targetPoints(row+1),real64)
    enddo

    do row = 0,this%M

      rowHasMatch = .false.

      do col = 0,this%N

        iMatrix(row,col) = 0.0_real64

        if(AlmostEqual(targetPoints(row),controlPoints(col))) then
          rowHasMatch = .true.
          iMatrix(row,col) = 1.0_real64
        endif

      enddo

      if(.not.(rowHasMatch)) then

        temp1 = 0.0_real64

        do col = 0,this%N
          temp2 = bWeights(col)/ &
                  (targetPoints(row)- &
                   controlPoints(col))
          iMatrix(row,col) = temp2
          temp1 = temp1+temp2
        enddo

        do col = 0,this%N
          iMatrix(row,col) = iMatrix(row,col)/temp1
        enddo

      endif

    enddo

    do row = 0,this%M
      do col = 0,this%N
        this%iMatrix(col+1,row+1) = real(iMatrix(row,col),prec)
      enddo
    enddo

  endsubroutine CalculateInterpolationMatrix

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
    real(real64) :: dmat(0:this%N,0:this%N)
    real(real64) :: dgmat(0:this%N,0:this%N)
    real(real64) :: bWeights(0:this%N)
    real(real64) :: qWeights(0:this%N)
    real(real64) :: controlPoints(0:this%N)

    do row = 0,this%N
      bWeights(row) = real(this%bWeights(row+1),real64)
      qWeights(row) = real(this%qWeights(row+1),real64)
      controlPoints(row) = real(this%controlPoints(row+1),real64)
    enddo

    do row = 0,this%N

      dmat(row,row) = 0.0_prec

      do col = 0,this%N

        if(.not.(col == row)) then

          dmat(row,col) = bWeights(col)/ &
                          (bWeights(row)* &
                           (controlPoints(row)- &
                            controlPoints(col)))

          dmat(row,row) = dmat(row,row)-dmat(row,col)

        endif

      enddo

    enddo

    do row = 0,this%N
      do col = 0,this%N
        dgmat(row,col) = -dmat(col,row)* &
                         qWeights(col)/ &
                         qWeights(row)
      enddo
    enddo

    do row = 0,this%N
      do col = 0,this%N
        this%dMatrix(row+1,col+1) = real(dmat(col,row),prec)
        this%dgMatrix(row+1,col+1) = real(dgmat(col,row),prec)
      enddo
    enddo

  endsubroutine CalculateDerivativeMatrix

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
    real(prec)      :: lAtS(0:this%N)
    ! Local
    integer    :: j
    logical    :: xMatchesNode
    real(real64) :: temp1,temp2
    real(real64) :: sELocal
    real(real64) :: controlPoints(0:this%N)
    real(real64) :: bWeights(0:this%N)
    real(real64) :: lS(0:this%N)

    sELocal = real(sE,real64)
    do j = 0,this%N
      controlPoints(j) = real(this%controlPoints(j+1),real64)
      bWeights(j) = real(this%bWeights(j+1),real64)
    enddo

    xMatchesNode = .false.

    do j = 0,this%N

      lS(j) = 0.0_real64
      if(AlmostEqual(sELocal,controlPoints(j))) then
        lS(j) = 1.0_real64
        xMatchesNode = .true.
      endif

    enddo

    if(xMatchesNode) then
      do j = 0,this%N
        lAtS(j) = real(lS(j),prec)
      enddo
      return
    endif

    temp1 = 0.0_real64

    do j = 0,this%N
      temp2 = bWeights(j)/(sE-controlPoints(j))
      lS(j) = temp2
      temp1 = temp1+temp2
    enddo

    lS = lS/temp1

    do j = 0,this%N
      lAtS(j) = real(lS(j),prec)
    enddo

  endfunction CalculateLagrangePolynomials

  subroutine WriteHDF5_Lagrange(this,fileId)
    implicit none
    class(Lagrange),intent(in) :: this
    integer(HID_T),intent(in) :: fileId

    call CreateGroup_HDF5(fileId,'/interp')

    call WriteArray_HDF5(fileId,'/interp/controlpoints', &
                         this%controlPoints)

    call WriteArray_HDF5(fileId,'/interp/qweights', &
                         this%qWeights)

    call WriteArray_HDF5(fileId,'/interp/dgmatrix', &
                         this%dgMatrix)

    call WriteArray_HDF5(fileId,'/interp/dmatrix', &
                         this%dMatrix)

    call WriteArray_HDF5(fileId,'/interp/bmatrix', &
                         this%bMatrix)

    call WriteArray_HDF5(fileId,'/interp/imatrix', &
                         this%iMatrix)

  endsubroutine WriteHDF5_Lagrange

endmodule SELF_Lagrange
