! SELF_Lagrange.f90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self-fluids@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE SELF_Lagrange

  USE SELF_Constants
  USE SELF_Memory
  USE SELF_SupportRoutines
  USE SELF_Quadrature

  USE hipfort
  USE ISO_C_BINDING

  IMPLICIT NONE

!INCLUDE 'SELF_Macros.h'

! =============================================================================================== !
!
!
! The Lagrange data-structure stores the information necessary to interpolate between two
! sets of grid-points and to estimate the derivative of data at native grid points. Routines for
! multidimensional interpolation are based on the tensor product of 1-D interpolants. It is
! assumed that the polynomial degree (and the interpolation nodes) are the same in each direction.
! This assumption permits the storage of only one array of interpolation nodes and barycentric
! weights and is what allows this data structure to be flexible.
!
! =============================================================================================== !

  TYPE,PUBLIC :: Lagrange
    !> A data structure for working with Lagrange Interpolating Polynomials in one, two, and three dimensions.
    !>
    !> attribute : controlPoints(0:N) : The set of nodes in one dimension where data is known. To create higher dimension interpolation and differentiation operators, structured grids in two and three dimensions are created by tensor products of the controlPoints. This design decision implies that all Spectral Element Methods supported by the Lagrange class have the same polynomial degree in each computational/spatial dimension. In practice, the controlPoints are the Legendre-Gauss, Legendre-Gauss-Lobatto, Legendre-Gauss-Radau, Chebyshev-Gauss, Chebyshev-Gauss-Lobatto, or Chebyshev-Gauss-Radau quadrature points over the domain [-1,1] (computational space). The Init routine for this class restricts controlPoints to one of these quadrature types or uniform points on [-1,1].
    !
    !> attribute : targetPoints(0:M) : The set of nodes in one dimension where data is to be interpolated to. To create higher dimension interpolation and differentiation operators, structured grids in two and three dimensions are created by tensor products of the targetPoints. In practice, the targetPoints are set to a uniformly distributed set of points between [-1,1] (computational space) to allow for interpolation from unevenly spaced quadrature points to a plotting grid.
    !
    !> attribute : bWeights(0:N) : The barycentric weights that are calculated from the controlPoints and used for interpolation.
    !
    !> attribute : qWeights(0:N) : The quadrature weights for discrete integration. The quadradture weights depend on the type of controlPoints provided; one of Legendre-Gauss, Legendre-Gauss-Lobatto, Legendre-Gauss-Radau, Chebyshev-Gauss, Chebyshev-Gauss-Lobatto, Chebyshev-Gauss Radau, or Uniform. If Uniform, the quadrature weights are constant = dx = 2.0/(N+1).
    !
    !> attribute : iMatrix(0:N,0:M) :
    !> attribute : dMatrix(0:N,0:N) : TO DO
    !> attribute : dgMatrix(0:N,0:N) : TO DO
    !> attribute : bMatrix(0:N,1:2) : TO DO

    INTEGER :: N
    INTEGER :: M
    TYPE(hfReal_r1) :: controlPoints
    TYPE(hfReal_r1) :: targetPoints
    TYPE(hfReal_r1) :: bWeights
    TYPE(hfReal_r1) :: qWeights
    TYPE(hfReal_r2) :: iMatrix
    TYPE(hfReal_r2) :: dMatrix
    TYPE(hfReal_r2) :: dgMatrix
    TYPE(hfReal_r2) :: bMatrix

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Lagrange
    PROCEDURE,PUBLIC :: Free => Free_Lagrange

#ifdef GPU
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Lagrange
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Lagrange
#endif

    GENERIC,PUBLIC :: ScalarGridInterp_1D => ScalarGridInterp_1D_cpu,ScalarGridInterp_1D_gpu
    PROCEDURE,PRIVATE :: ScalarGridInterp_1D_cpu,ScalarGridInterp_1D_gpu

    GENERIC,PUBLIC :: ScalarGridInterp_2D => ScalarGridInterp_2D_cpu,ScalarGridInterp_2D_gpu
    PROCEDURE,PRIVATE :: ScalarGridInterp_2D_cpu,ScalarGridInterp_2D_gpu

    GENERIC,PUBLIC :: VectorGridInterp_2D => VectorGridInterp_2D_cpu,VectorGridInterp_2D_gpu
    PROCEDURE,PRIVATE :: VectorGridInterp_2D_cpu,VectorGridInterp_2D_gpu

    GENERIC,PUBLIC :: TensorGridInterp_2D => TensorGridInterp_2D_cpu,TensorGridInterp_2D_gpu
    PROCEDURE,PRIVATE :: TensorGridInterp_2D_cpu,TensorGridInterp_2D_gpu

    GENERIC,PUBLIC :: ScalarGridInterp_3D => ScalarGridInterp_3D_cpu,ScalarGridInterp_3D_gpu
    PROCEDURE,PRIVATE :: ScalarGridInterp_3D_cpu,ScalarGridInterp_3D_gpu

    GENERIC,PUBLIC :: VectorGridInterp_3D => VectorGridInterp_3D_cpu,VectorGridInterp_3D_gpu
    PROCEDURE,PRIVATE :: VectorGridInterp_3D_cpu,VectorGridInterp_3D_gpu

    GENERIC,PUBLIC :: TensorGridInterp_3D => TensorGridInterp_3D_cpu,TensorGridInterp_3D_gpu
    PROCEDURE,PRIVATE :: TensorGridInterp_3D_cpu,TensorGridInterp_3D_gpu

    GENERIC,PUBLIC :: ScalarBoundaryInterp_1D => ScalarBoundaryInterp_1D_cpu,ScalarBoundaryInterp_1D_gpu
    PROCEDURE,PRIVATE :: ScalarBoundaryInterp_1D_cpu,ScalarBoundaryInterp_1D_gpu

    GENERIC,PUBLIC :: ScalarBoundaryInterp_2D => ScalarBoundaryInterp_2D_cpu,ScalarBoundaryInterp_2D_gpu
    PROCEDURE,PRIVATE :: ScalarBoundaryInterp_2D_cpu,ScalarBoundaryInterp_2D_gpu

    GENERIC,PUBLIC :: VectorBoundaryInterp_2D => VectorBoundaryInterp_2D_cpu,VectorBoundaryInterp_2D_gpu
    PROCEDURE,PRIVATE :: VectorBoundaryInterp_2D_cpu,VectorBoundaryInterp_2D_gpu

    GENERIC,PUBLIC :: TensorBoundaryInterp_2D => TensorBoundaryInterp_2D_cpu,TensorBoundaryInterp_2D_gpu
    PROCEDURE,PRIVATE :: TensorBoundaryInterp_2D_cpu,TensorBoundaryInterp_2D_gpu

    GENERIC,PUBLIC :: ScalarBoundaryInterp_3D => ScalarBoundaryInterp_3D_cpu,ScalarBoundaryInterp_3D_gpu
    PROCEDURE,PRIVATE :: ScalarBoundaryInterp_3D_cpu,ScalarBoundaryInterp_3D_gpu

    GENERIC,PUBLIC :: VectorBoundaryInterp_3D => VectorBoundaryInterp_3D_cpu,VectorBoundaryInterp_3D_gpu
    PROCEDURE,PRIVATE :: VectorBoundaryInterp_3D_cpu,VectorBoundaryInterp_3D_gpu

    GENERIC,PUBLIC :: TensorBoundaryInterp_3D => TensorBoundaryInterp_3D_cpu,TensorBoundaryInterp_3D_gpu
    PROCEDURE,PRIVATE :: TensorBoundaryInterp_3D_cpu,TensorBoundaryInterp_3D_gpu

    GENERIC,PUBLIC :: Derivative_1D => Derivative_1D_cpu,Derivative_1D_gpu
    PROCEDURE,PRIVATE :: Derivative_1D_cpu,Derivative_1D_gpu

    GENERIC,PUBLIC :: DGDerivative_1D => DGDerivative_1D_cpu,DGDerivative_1D_gpu
    PROCEDURE,PRIVATE :: DGDerivative_1D_cpu,DGDerivative_1D_gpu

    GENERIC,PUBLIC :: ScalarGradient_2D => ScalarGradient_2D_cpu,ScalarGradient_2D_gpu
    PROCEDURE,PRIVATE :: ScalarGradient_2D_cpu,ScalarGradient_2D_gpu

    GENERIC,PUBLIC :: VectorGradient_2D => VectorGradient_2D_cpu,VectorGradient_2D_gpu
    PROCEDURE,PRIVATE :: VectorGradient_2D_cpu,VectorGradient_2D_gpu

    GENERIC,PUBLIC :: VectorDivergence_2D => VectorDivergence_2D_cpu,VectorDivergence_2D_gpu
    PROCEDURE,PRIVATE :: VectorDivergence_2D_cpu,VectorDivergence_2D_gpu

    GENERIC,PUBLIC :: VectorCurl_2D => VectorCurl_2D_cpu,VectorCurl_2D_gpu
    PROCEDURE,PRIVATE :: VectorCurl_2D_cpu,VectorCurl_2D_gpu

    GENERIC,PUBLIC :: TensorDivergence_2D => TensorDivergence_2D_cpu,TensorDivergence_2D_gpu
    PROCEDURE,PRIVATE :: TensorDivergence_2D_cpu,TensorDivergence_2D_gpu

    GENERIC,PUBLIC :: ScalarGradient_3D => ScalarGradient_3D_cpu,ScalarGradient_3D_gpu
    PROCEDURE,PRIVATE :: ScalarGradient_3D_cpu,ScalarGradient_3D_gpu

    GENERIC,PUBLIC :: VectorGradient_3D => VectorGradient_3D_cpu,VectorGradient_3D_gpu
    PROCEDURE,PRIVATE :: VectorGradient_3D_cpu,VectorGradient_3D_gpu

    GENERIC,PUBLIC :: VectorDivergence_3D => VectorDivergence_3D_cpu,VectorDivergence_3D_gpu
    PROCEDURE,PRIVATE :: VectorDivergence_3D_cpu,VectorDivergence_3D_gpu

    GENERIC,PUBLIC :: VectorCurl_3D => VectorCurl_3D_cpu,VectorCurl_3D_gpu
    PROCEDURE,PRIVATE :: VectorCurl_3D_cpu,VectorCurl_3D_gpu

    GENERIC,PUBLIC :: TensorDivergence_3D => TensorDivergence_3D_cpu,TensorDivergence_3D_gpu
    PROCEDURE,PRIVATE :: TensorDivergence_3D_cpu,TensorDivergence_3D_gpu

    PROCEDURE,PRIVATE :: CalculateBarycentricWeights
    PROCEDURE,PRIVATE :: CalculateInterpolationMatrix
    PROCEDURE,PRIVATE :: CalculateDerivativeMatrix
    PROCEDURE,PRIVATE :: CalculateLagrangePolynomials

  END TYPE Lagrange

  INTERFACE
    SUBROUTINE ScalarGridInterp_1D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="ScalarGridInterp_1D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER,VALUE :: N,M,nVar,nEl
    END SUBROUTINE ScalarGridInterp_1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarGridInterp_2D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="ScalarGridInterp_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER,VALUE :: N,M,nVar,nEl
    END SUBROUTINE ScalarGridInterp_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorGridInterp_2D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="VectorGridInterp_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER,VALUE :: N,M,nVar,nEl
    END SUBROUTINE VectorGridInterp_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorGridInterp_2D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="TensorGridInterp_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER,VALUE :: N,M,nVar,nEl
    END SUBROUTINE TensorGridInterp_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarGridInterp_3D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="ScalarGridInterp_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER,VALUE :: N,M,nVar,nEl
    END SUBROUTINE ScalarGridInterp_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorGridInterp_3D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="VectorGridInterp_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER,VALUE :: N,M,nVar,nEl
    END SUBROUTINE VectorGridInterp_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorGridInterp_3D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="TensorGridInterp_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER,VALUE :: N,M,nVar,nEl
    END SUBROUTINE TensorGridInterp_3D_gpu_wrapper
  END INTERFACE

  ! /////////////// !
  ! Boundary Interpolation Routines

  INTERFACE
    SUBROUTINE ScalarBoundaryInterp_1D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="ScalarBoundaryInterp_1D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ScalarBoundaryInterp_1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="ScalarBoundaryInterp_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ScalarBoundaryInterp_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="VectorBoundaryInterp_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE VectorBoundaryInterp_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="TensorBoundaryInterp_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE TensorBoundaryInterp_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="ScalarBoundaryInterp_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ScalarBoundaryInterp_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="VectorBoundaryInterp_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE VectorBoundaryInterp_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="TensorBoundaryInterp_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE TensorBoundaryInterp_3D_gpu_wrapper
  END INTERFACE

  ! /////////////// !

  INTERFACE
    SUBROUTINE Derivative_1D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="Derivative_1D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE Derivative_1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE DGDerivative_1D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
      bind(c,name="DGDerivative_1D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE DGDerivative_1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarGradient_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="ScalarGradient_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ScalarGradient_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorGradient_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorGradient_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE VectorGradient_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorDivergence_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorDivergence_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE VectorDivergence_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorCurl_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorCurl_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE VectorCurl_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorDivergence_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="TensorDivergence_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE TensorDivergence_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarGradient_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="ScalarGradient_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE ScalarGradient_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorGradient_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorGradient_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE VectorGradient_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorDivergence_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorDivergence_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE VectorDivergence_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorCurl_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorCurl_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE VectorCurl_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorDivergence_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="TensorDivergence_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER,VALUE :: N,nVar,nEl
    END SUBROUTINE TensorDivergence_3D_gpu_wrapper
  END INTERFACE

CONTAINS

! ================================================================================================ !
!
! Init_Lagrange
!
!   A manual constructor for the Lagrange class that allocates memory and fills in data
!   for the attributes of the Lagrange class.
!
!   The Init subroutine allocates memory for the interpolation and target points, barycentric
!   weights, interpolation matrix, and derivative matrix.
!
!   Usage :
!
!     TYPE(Lagrange) :: interp
!     INTEGER        :: N, M
!     REAL(prec)     :: interpNodes(0:N), targetNodes(0:M+1)
!
!     CALL interp % Init( N, M, interpNodes, targetNodes )
!
!   Input/Output :
!
!     myPoly  (out)
!       The Lagrange data structure to be constructed
!
!     N (in)
!       The degree of the polynomial interpolant
!
!     M (in)
!       M+1 is the number of target grid points. The upper bound of the targetNodes array
!
!     interpNodes(0:N) (in)
!       The interpolation nodes.
!
!     targetNodes(0:M) (in)
!       The target nodes. The iMatrix will map the a function at the interpolation
!       nodes onto the target nodes.
!
! =============================================================================================== !

  SUBROUTINE Init_Lagrange(myPoly,N,controlNodeType,M,targetNodeType)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(out) :: myPoly
    INTEGER,INTENT(in)          :: N,M
    INTEGER,INTENT(in)          :: controlNodeType,targetNodeType
    ! Local
    REAL(prec) :: q(0:M)

    myPoly % N = N
    myPoly % M = M

    CALL myPoly % controlPoints % Alloc(loBound=0, &
                                        upBound=N)

    CALL myPoly % targetPoints % Alloc(loBound=0, &
                                       upBound=M)

    CALL myPoly % bWeights % Alloc(loBound=0, &
                                   upBound=N)

    CALL myPoly % qWeights % Alloc(loBound=0, &
                                   upBound=N)

    CALL myPoly % iMatrix % Alloc(loBound=(/0,0/), &
                                  upBound=(/N,M/))

    CALL myPoly % dMatrix % Alloc(loBound=(/0,0/), &
                                  upBound=(/N,N/))

    CALL myPoly % dgMatrix % Alloc(loBound=(/0,0/), &
                                   upBound=(/N,N/))

    CALL myPoly % bMatrix % Alloc(loBound=(/0,0/), &
                                  upBound=(/N,1/))

    IF (controlNodeType == GAUSS .OR. controlNodeType == GAUSS_LOBATTO) THEN

      CALL LegendreQuadrature(N, &
                              myPoly % controlPoints % hostData, &
                              myPoly % qWeights % hostData, &
                              controlNodeType)

    ELSEIF (controlNodeType == UNIFORM) THEN

      myPoly % controlPoints % hostData = UniformPoints(-1.0_prec,1.0_prec,0,N)
      myPoly % qWeights % hostData = 2.0_prec/REAL(N,prec)

    END IF

    ! Target Points
    IF (targetNodeType == GAUSS .OR. targetNodeType == GAUSS_LOBATTO) THEN

      CALL LegendreQuadrature(M, &
                              myPoly % targetPoints % hostData, &
                              q, &
                              targetNodeType)

    ELSEIF (targetNodeType == UNIFORM) THEN

      myPoly % targetPoints % hostData = UniformPoints(-1.0_prec,1.0_prec,0,M)

    END IF

    CALL myPoly % CalculateBarycentricWeights()
    CALL myPoly % CalculateInterpolationMatrix()
    CALL myPoly % CalculateDerivativeMatrix()
    myPoly % bMatrix % hostData(0:N,0) = myPoly % CalculateLagrangePolynomials(-1.0_prec)
    myPoly % bMatrix % hostData(0:N,1) = myPoly % CalculateLagrangePolynomials(1.0_prec)

#ifdef GPU
    CALL myPoly % UpdateDevice()
#endif

  END SUBROUTINE Init_Lagrange

! ================================================================================================ !
!
! Free_Lagrange
!
!   A manual destructor for the Lagrange class that deallocates the memory held by its attributes.
!
! ================================================================================================ !

  SUBROUTINE Free_Lagrange(myPoly)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(inout) :: myPoly

    CALL myPoly % controlPoints % Free()
    CALL myPoly % targetPoints % Free()
    CALL myPoly % bWeights % Free()
    CALL myPoly % qWeights % Free()
    CALL myPoly % iMatrix % Free()
    CALL myPoly % dMatrix % Free()
    CALL myPoly % dgMatrix % Free()
    CALL myPoly % bMatrix % Free()

  END SUBROUTINE Free_Lagrange

#ifdef GPU
  SUBROUTINE UpdateDevice_Lagrange(myPoly)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(inout) :: myPoly

    CALL myPoly % controlPoints % UpdateDevice()
    CALL myPoly % targetPoints % UpdateDevice()
    CALL myPoly % bWeights % UpdateDevice()
    CALL myPoly % qWeights % UpdateDevice()
    CALL myPoly % iMatrix % UpdateDevice()
    CALL myPoly % dMatrix % UpdateDevice()
    CALL myPoly % dgMatrix % UpdateDevice()
    CALL myPoly % bMatrix % UpdateDevice()

  END SUBROUTINE UpdateDevice_Lagrange

  SUBROUTINE UpdateHost_Lagrange(myPoly)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(inout) :: myPoly

    CALL myPoly % controlPoints % UpdateHost()
    CALL myPoly % targetPoints % UpdateHost()
    CALL myPoly % bWeights % UpdateHost()
    CALL myPoly % qWeights % UpdateHost()
    CALL myPoly % iMatrix % UpdateHost()
    CALL myPoly % dMatrix % UpdateHost()
    CALL myPoly % dgMatrix % UpdateHost()
    CALL myPoly % bMatrix % UpdateHost()

  END SUBROUTINE UpdateHost_Lagrange
#endif

! ================================================================================================ !
!
! ScalarGridInterp_1D
!
!   Interpolates an array of nodal values at the native nodes to nodal values at the target nodes.
!
!   We can write the operations of interpolating data from one set of points to another (in this
!   case from "controlPoints" to "targetPoints") as
!
!              fnew_a = \sum_{i=0}^N f_i l_i(\xi_a),   a=0,1,2,...,M
!
!   where l_i(\xi) are the Lagrange interpolating polynomials through the interpolation points.
!   The interpolation matrix is T_{a,i} = l_i(\xi_a) maps an array of nodal values from the native
!   interpolation nodes to the target nodes. This routine serves as a wrapper to call either the CUDA
!   kernel (if CUDA is enabled) or the CPU version.
!
!   Usage :
!
!     TYPE(Lagrange) :: interp
!     INTEGER        :: nVariables, nElements
!     REAL(prec)     :: f(0:interp % N,1:nVariables,1:nElements)
!     REAL(prec)     :: fNew(0:interp % M,1:nVariables,1:nElements)
!
!       CALL interp % GridInterpolation_1D( fnative, fNew, nVariables, nElements )
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
!     fNew (out)
!      Array of function nodal values at the target interpolation nodes.
!
! ================================================================================================ !

  SUBROUTINE ScalarGridInterp_1D_cpu(myPoly,f,fInterp,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: fInterp(0:myPoly % M,1:nVariables,1:nElements)
    ! Local
    INTEGER :: iVar,iEl,i,ii

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO i = 0,myPoly % M
          fInterp(i,iVar,iEl) = 0.0_prec
          DO ii = 0,myPoly % N
            fInterp(i,iVar,iEl) = fInterp(i,iVar,iEl) + myPoly % iMatrix % hostData(ii,i)*f(ii,iVar,iEl)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ScalarGridInterp_1D_cpu

  SUBROUTINE ScalarGridInterp_1D_gpu(myPoly,f_dev,fInterp_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in) :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f_dev
    TYPE(c_ptr),INTENT(out) :: fInterp_dev

    CALL ScalarGridInterp_1D_gpu_wrapper(myPoly % iMatrix % deviceData, &
                                         f_dev,fInterp_dev, &
                                         myPoly % N,myPoly % M, &
                                         nVariables,nElements)

  END SUBROUTINE ScalarGridInterp_1D_gpu
!
! ================================================================================================ !
!
! ScalarGridInterp_2D
!
!   Interpolates an array of nodal values at the native nodes to nodal values at the target nodes.
!
!   We can write the operations of interpolating data from one set of points to another (in this
!   case from "controlPoints" to "targetPoints") as
!
!              fnew_{a,b} = \sum_{i,j=0}^N f_{i,j} l_i(\xi_a) l_j(\xi_b),   a,b=0,1,2,...,M
!
!   where l_i(\xi) are the Lagrange interpolating polynomials through the interpolation points.
!   The interpolation matrix is T_{a,i} = l_i(\xi_a) maps an array of nodal values from the native
!   interpolation nodes to the target nodes. This routine serves as a wrapper to call either the CUDA
!   kernel (if CUDA is enabled) or the CPU version.
!
!   Usage :
!
!     TYPE(Lagrange) :: interp
!     INTEGER        :: nVariables, nElements
!     REAL(prec)     :: f(0:interp % N,0:interp % N,1:nVariables,1:nElements)
!     REAL(prec)     :: fNew(0:interp % M,0:interp % M,1:nVariables,1:nElements)
!
!       CALL interp % GridInterpolation_2D( fnative, fNew, nVariables, nElements )
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
!     fNew (out)
!      Array of function nodal values at the target interpolation nodes.
!
! ================================================================================================ !

  SUBROUTINE ScalarGridInterp_2D_cpu(myPoly,f,fNew,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: fNew(0:myPoly % M,0:myPoly % M,1:nVariables,1:nElements)
    ! Local
    INTEGER :: i,j,ii,jj,p,iEl,iVar
    REAL(prec) :: fi,fij

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % M
          DO i = 0,myPoly % M

            fij = 0.0_prec
            DO jj = 0,myPoly % N
              fi = 0.0_prec
              DO ii = 0,myPoly % N
                fi = fi + f(ii,jj,iVar,iEl)*myPoly % iMatrix % hostData(ii,i)
              END DO
              fij = fij + fi*myPoly % iMatrix % hostData(jj,j)
            END DO
            fNew(i,j,iVar,iEl) = fij

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ScalarGridInterp_2D_cpu
!
  SUBROUTINE ScalarGridInterp_2D_gpu(myPoly,f_dev,fInterp_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in) :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f_dev
    TYPE(c_ptr),INTENT(out) :: fInterp_dev

    CALL ScalarGridInterp_2D_gpu_wrapper(myPoly % iMatrix % deviceData, &
                                         f_dev,fInterp_dev, &
                                         myPoly % N,myPoly % M, &
                                         nVariables,nElements)

  END SUBROUTINE ScalarGridInterp_2D_gpu

  SUBROUTINE VectorGridInterp_2D_cpu(myPoly,f,fNew,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: fNew(1:2,0:myPoly % M,0:myPoly % M,1:nVariables,1:nElements)
    ! Local
    INTEGER :: i,j,ii,jj,p,iEl,iVar
    REAL(prec) :: fi(1:2)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % M
          DO i = 0,myPoly % M

            fNew(1,i,j,iVar,iEl) = 0.0_prec
            fNew(2,i,j,iVar,iEl) = 0.0_prec

            DO jj = 0,myPoly % N

              fi(1:2) = 0.0_prec
              DO ii = 0,myPoly % N
                fi(1:2) = fi(1:2) + f(1:2,ii,jj,iVar,iEl)*myPoly % iMatrix % hostData(ii,i)
              END DO

              fNew(1:2,i,j,iVar,iEl) = fNew(1:2,i,j,iVar,iEl) + fi(1:2)*myPoly % iMatrix % hostData(jj,j)

            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE VectorGridInterp_2D_cpu
!
  SUBROUTINE VectorGridInterp_2D_gpu(myPoly,f_dev,fInterp_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in) :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f_dev
    TYPE(c_ptr),INTENT(out) :: fInterp_dev

    CALL VectorGridInterp_2D_gpu_wrapper(myPoly % iMatrix % deviceData, &
                                         f_dev,fInterp_dev, &
                                         myPoly % N,myPoly % M, &
                                         nVariables,nElements)

  END SUBROUTINE VectorGridInterp_2D_gpu

  SUBROUTINE TensorGridInterp_2D_cpu(myPoly,f,fNew,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:2,1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: fNew(1:2,1:2,0:myPoly % M,0:myPoly % M,1:nVariables,1:nElements)
    ! Local
    INTEGER :: i,j,ii,jj,p,iEl,iVar
    REAL(prec) :: fi(1:2,1:2)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % M
          DO i = 0,myPoly % M

            fNew(1:2,1:2,i,j,iVar,iEl) = 0.0_prec

            DO jj = 0,myPoly % N

              fi(1:2,1:2) = 0.0_prec
              DO ii = 0,myPoly % N
                fi(1:2,1:2) = fi(1:2,1:2) + f(1:2,1:2,ii,jj,iVar,iEl)*myPoly % iMatrix % hostData(ii,i)
              END DO

              fNew(1:2,1:2,i,j,iVar,iEl) = fNew(1:2,1:2,i,j,iVar,iEl) + fi(1:2,1:2)*myPoly % iMatrix % hostData(jj,j)

            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE TensorGridInterp_2D_cpu
!
  SUBROUTINE TensorGridInterp_2D_gpu(myPoly,f_dev,fInterp_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in) :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f_dev
    TYPE(c_ptr),INTENT(out) :: fInterp_dev

    CALL TensorGridInterp_2D_gpu_wrapper(myPoly % iMatrix % deviceData, &
                                         f_dev,fInterp_dev, &
                                         myPoly % N,myPoly % M, &
                                         nVariables,nElements)

  END SUBROUTINE TensorGridInterp_2D_gpu

! ================================================================================================ !
!
! GridInterpolate_3D
!
!   Interpolates an array of nodal values at the native nodes to nodal values at the target nodes.
!
!   We can write the operations of interpolating data from one set of points to another (in this
!   case from "controlPoints" to "targetPoints") as
!
!   fnew_{a,b,c} = \sum_{i,j,k=0}^N f_{i,j,k} l_i(\xi_a) l_j(\xi_b) l_k(\xi_c),   a,b,c=0,1,2,...,M
!
!   where l_i(\xi) are the Lagrange interpolating polynomials through the interpolation points.
!   The interpolation matrix is T_{a,i} = l_i(\xi_a) maps an array of nodal values from the native
!   interpolation nodes to the target nodes. This routine serves as a wrapper to call either the CUDA
!   kernel (if CUDA is enabled) or the CPU version.
!
!   Usage :
!
!     TYPE(Lagrange) :: interp
!     INTEGER        :: nVariables, nElements
!     REAL(prec)     :: f(0:interp % N,0:interp % N,0:interp % N,1:nVariables,1:nElements)
!     REAL(prec)     :: fNew(0:interp % M,0:interp % M,0:interp % M,1:nVariables,1:nElements)
!
!       CALL interp % GridInterpolation_3D( fnative, fNew, nVariables, nElements )
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
!     fNew (out)
!      Array of function nodal values at the target interpolation nodes.
!
! ================================================================================================ !

  SUBROUTINE ScalarGridInterp_3D_cpu(myPoly,f,fInterp,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: fInterp(0:myPoly % M,0:myPoly % M,0:myPoly % M,1:nVariables,1:nElements)
    ! Local
    INTEGER :: iEl,iVar,i,j,k,ii,jj,kk
    REAL(prec) :: fi,fij,fijk

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO k = 0,myPoly % M
          DO j = 0,myPoly % M
            DO i = 0,myPoly % M

              fijk = 0.0_prec
              DO kk = 0,myPoly % N
                fij = 0.0_prec
                DO jj = 0,myPoly % N
                  fi = 0.0_prec
                  DO ii = 0,myPoly % N
                    fi = fi + f(ii,jj,kk,iVar,iEl)*myPoly % iMatrix % hostData(ii,i)
                  END DO
                  fij = fij + fi*myPoly % iMatrix % hostData(jj,j)
                END DO
                fijk = fijk + fij*myPoly % iMatrix % hostData(kk,k)
              END DO
              fInterp(i,j,k,iVar,iEl) = fijk

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ScalarGridInterp_3D_cpu
!
  SUBROUTINE ScalarGridInterp_3D_gpu(myPoly,f_dev,fInterp_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in) :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f_dev
    TYPE(c_ptr),INTENT(out) :: fInterp_dev

    CALL ScalarGridInterp_3D_gpu_wrapper(myPoly % iMatrix % deviceData, &
                                         f_dev,fInterp_dev, &
                                         myPoly % N,myPoly % M, &
                                         nVariables,nElements)

  END SUBROUTINE ScalarGridInterp_3D_gpu

  SUBROUTINE VectorGridInterp_3D_cpu(myPoly,f,fInterp,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: fInterp(1:3,0:myPoly % M,0:myPoly % M,0:myPoly % M,1:nVariables,1:nElements)
    ! Local
    INTEGER :: iEl,iVar,i,j,k,ii,jj,kk
    REAL(prec) :: fi(1:3),fij(1:3)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO k = 0,myPoly % M
          DO j = 0,myPoly % M
            DO i = 0,myPoly % M

              fInterp(1:3,i,j,k,iVar,iEl) = 0.0_prec
              DO kk = 0,myPoly % N
                fij(1:3) = 0.0_prec
                DO jj = 0,myPoly % N
                  fi(1:3) = 0.0_prec
                  DO ii = 0,myPoly % N
                    fi(1:3) = fi(1:3) + f(1:3,ii,jj,kk,iVar,iEl)*myPoly % iMatrix % hostData(ii,i)
                  END DO
                  fij(1:3) = fij(1:3) + fi(1:3)*myPoly % iMatrix % hostData(jj,j)
                END DO
                fInterp(1:3,i,j,k,iVar,iEl) = fInterp(1:3,i,j,k,iVar,iEl) + fij(1:3)*myPoly % iMatrix % hostData(kk,k)
              END DO

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE VectorGridInterp_3D_cpu
!
  SUBROUTINE VectorGridInterp_3D_gpu(myPoly,f_dev,fInterp_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in) :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f_dev
    TYPE(c_ptr),INTENT(out) :: fInterp_dev

    CALL VectorGridInterp_3D_gpu_wrapper(myPoly % iMatrix % deviceData, &
                                         f_dev,fInterp_dev, &
                                         myPoly % N,myPoly % M, &
                                         nVariables,nElements)

  END SUBROUTINE VectorGridInterp_3D_gpu
!
  SUBROUTINE TensorGridInterp_3D_cpu(myPoly,f,fInterp,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:3,1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: fInterp(1:3,1:3,0:myPoly % M,0:myPoly % M,0:myPoly % M,1:nVariables,1:nElements)
    ! Local
    INTEGER :: iEl,iVar,i,j,k,ii,jj,kk
    REAL(prec) :: fi(1:3,1:3),fij(1:3,1:3)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO k = 0,myPoly % M
          DO j = 0,myPoly % M
            DO i = 0,myPoly % M

              fInterp(1:3,1:3,i,j,k,iVar,iEl) = 0.0_prec
              DO kk = 0,myPoly % N
                fij(1:3,1:3) = 0.0_prec
                DO jj = 0,myPoly % N
                  fi(1:3,1:3) = 0.0_prec
                  DO ii = 0,myPoly % N
                    fi(1:3,1:3) = fi(1:3,1:3) + f(1:3,1:3,ii,jj,kk,iVar,iEl)*myPoly % iMatrix % hostData(ii,i)
                  END DO
                  fij(1:3,1:3) = fij(1:3,1:3) + fi(1:3,1:3)*myPoly % iMatrix % hostData(jj,j)
                END DO
                fInterp(1:3,1:3,i,j,k,iVar,iEl) = fInterp(1:3,1:3,i,j,k,iVar,iEl) + &
                                                  fij(1:3,1:3)*myPoly % iMatrix % hostData(kk,k)
              END DO

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE TensorGridInterp_3D_cpu
!
  SUBROUTINE TensorGridInterp_3D_gpu(myPoly,f_dev,fInterp_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in) :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f_dev
    TYPE(c_ptr),INTENT(out) :: fInterp_dev

    CALL TensorGridInterp_3D_gpu_wrapper(myPoly % iMatrix % deviceData, &
                                         f_dev,fInterp_dev, &
                                         myPoly % N,myPoly % M, &
                                         nVariables,nElements)

  END SUBROUTINE TensorGridInterp_3D_gpu

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

  SUBROUTINE Derivative_1D_cpu(myPoly,f,df,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: df(0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER :: i,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO i = 0,myPoly % N

          df(i,iVar,iEl) = 0.0_prec
          DO ii = 0,myPoly % N
            df(i,iVar,iEl) = df(i,iVar,iEl) + myPoly % dMatrix % hostData(ii,i)*f(ii,iVar,iEl)
          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE Derivative_1D_cpu

  SUBROUTINE Derivative_1D_gpu(myPoly,f_dev,df_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in) :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f_dev
    TYPE(c_ptr),INTENT(out) :: df_dev

    CALL Derivative_1D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                   f_dev,df_dev, &
                                   myPoly % N, &
                                   nVariables,nElements)

  END SUBROUTINE Derivative_1D_gpu

  SUBROUTINE DGDerivative_1D_cpu(myPoly,f,bf,df,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(in)  :: bf(1:nVariables,1:2,1:nElements)
    REAL(prec),INTENT(out) :: df(0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER :: i,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO i = 0,myPoly % N

          ! Interior Derivative Matrix Application
          df(i,iVar,iEl) = 0.0_prec
          DO ii = 0,myPoly % N
            df(i,iVar,iEl) = df(i,iVar,iEl) + myPoly % dgMatrix % hostData(ii,i)*f(ii,iVar,iEl)
          END DO

          ! Boundary Contribution
          df(i,iVar,iEl) = df(i,iVar,iEl) + (bf(iVar,2,iEl)*myPoly % bMatrix % hostData(i,1) - &
                                             bf(iVar,1,iEl)*myPoly % bMatrix % hostData(i,0))/&
                                             myPoly % qWeights % hostData(i)

        END DO

      END DO
    END DO

  END SUBROUTINE DGDerivative_1D_cpu

  SUBROUTINE DGDerivative_1D_gpu(myPoly,f_dev,bf_dev,df_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in) :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f_dev
    TYPE(c_ptr),INTENT(in)  :: bf_dev
    TYPE(c_ptr),INTENT(out) :: df_dev

    CALL DGDerivative_1D_gpu_wrapper(myPoly % dgMatrix % deviceData, &
                                     myPoly % bMatrix % deviceData, &
                                     myPoly % qWeights % deviceData, &
                                     f_dev,bf_dev,df_dev, &
                                     myPoly % N, &
                                     nVariables,nElements)

  END SUBROUTINE DGDerivative_1D_gpu
!
!! ================================================================================================ !
!!
!! CalculateGradient_2D
!!
!!   Calculates the gradient of a 2-D function, represented by a 2-D array of nodal values.
!!
!!   Given a set of nodal values at the interpolation nodes, the gradient of a function through
!!   the interpolation nodes can be estimated by
!!
!!                       (df/dx)_{a,b} = \sum_{i=0}^N f_{i,b} l'_i(\xi_a),   a,b=0,1,2,...,N
!!                       (df/dy)_{a,b} = \sum_{j=0}^N f_{a,j} l'_j(\xi_b),   a,b=0,1,2,...,N
!!
!!   where l_i(\xi) are the Lagrange interpolating polynomials through the interpolation points.
!!   The derivative matrix is D_{a,i} = l'_i(\xi_a) maps an array of nodal values at the interpolation
!!   nodes to its estimated derivative. This routine serves as a wrapper to call either the CUDA
!!   kernel (if CUDA is enabled) or the CPU version.
!!
!!   Usage :
!!
!!     TYPE(Lagrange) :: interp
!!     INTEGER        :: nVariables, nElements
!!     REAL(prec)     :: f(0:interp % N,0:interp % N,1:nVariables,1:nElements)
!!     REAL(prec)     :: gradF(1:2,0:interp % N,0:interp % N,1:nVariables,1:nElements)
!!
!!       CALL interp % CalculateGradient_2D( f, gradF, nVariables, nElements )
!!
!!     * If CUDA is enabled, the fnative and ftarget arrays must be CUDA device variables.
!!
!!   Parameters :
!!
!!     interp (in)
!!       A previously constructed Lagrange data-structure.
!!
!!     f (in)
!!       Array of function nodal values at the native interpolation nodes.
!!
!!     nVariables (in)
!!
!!     nElements (in)
!!
!!     gradF (out)
!!      Array of derivative values at the target interpolation nodes.
!!
!! ================================================================================================ !
!
  SUBROUTINE ScalarGradient_2D_cpu(myPoly,f,gradF,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: gradF(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % N
          DO i = 0,myPoly % N

            gradF(1,i,j,iVar,iEl) = 0.0_prec
            gradF(2,i,j,iVar,iEl) = 0.0_prec
            DO ii = 0,myPoly % N
              gradF(1,i,j,iVar,iEl) = gradF(1,i,j,iVar,iEl) + myPoly % dMatrix % hostData(ii,i)*f(ii,j,iVar,iEl)
              gradF(2,i,j,iVar,iEl) = gradF(2,i,j,iVar,iEl) + myPoly % dMatrix % hostData(ii,j)*f(i,ii,iVar,iEl)
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ScalarGradient_2D_cpu

  SUBROUTINE ScalarGradient_2D_gpu(myPoly,f_dev,gradF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(out)    :: gradF_dev
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    CALL ScalarGradient_2D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                       f_dev,gradF_dev,myPoly % N, &
                                       nVariables,nElements)

  END SUBROUTINE ScalarGradient_2D_gpu
!
  SUBROUTINE VectorGradient_2D_cpu(myPoly,f,gradF,nVariables,nElements)
    !
    ! Input : Vector(1:2,...)
    ! Output : Tensor(1:2,1:2,....)
    !          > Tensor(1,1) = d/ds1( Vector(1,...) )
    !          > Tensor(2,1) = d/ds1( Vector(2,...) )
    !          > Tensor(1,2) = d/ds2( Vector(1,...) )
    !          > Tensor(2,2) = d/ds2( Vector(2,...) )
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: gradF(1:2,1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % N
          DO i = 0,myPoly % N

            gradF(1,1,i,j,iVar,iEl) = 0.0_prec
            gradF(2,1,i,j,iVar,iEl) = 0.0_prec
            gradF(1,2,i,j,iVar,iEl) = 0.0_prec
            gradF(2,2,i,j,iVar,iEl) = 0.0_prec
            DO ii = 0,myPoly % N
              gradF(1,1,i,j,iVar,iEl) = gradF(1,1,i,j,iVar,iEl) + myPoly % dMatrix % hostData(ii,i)*f(1,ii,j,iVar,iEl)
              gradF(2,1,i,j,iVar,iEl) = gradF(2,1,i,j,iVar,iEl) + myPoly % dMatrix % hostData(ii,i)*f(2,ii,j,iVar,iEl)
              gradF(1,2,i,j,iVar,iEl) = gradF(1,2,i,j,iVar,iEl) + myPoly % dMatrix % hostData(ii,j)*f(1,i,ii,iVar,iEl)
              gradF(2,2,i,j,iVar,iEl) = gradF(2,2,i,j,iVar,iEl) + myPoly % dMatrix % hostData(ii,j)*f(2,i,ii,iVar,iEl)
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE VectorGradient_2D_cpu

  SUBROUTINE VectorGradient_2D_gpu(myPoly,f_dev,gradF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(out)    :: gradF_dev
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    CALL VectorGradient_2D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                       f_dev,gradF_dev,myPoly % N, &
                                       nVariables,nElements)

  END SUBROUTINE VectorGradient_2D_gpu

  SUBROUTINE VectorDivergence_2D_cpu(myPoly,f,dF,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: dF(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % N
          DO i = 0,myPoly % N

            dF(i,j,iVar,iEl) = 0.0_prec
            DO ii = 0,myPoly % N
              dF(i,j,iVar,iEl) = dF(i,j,iVar,iEl) + myPoly % dMatrix % hostData(ii,i)*f(1,ii,j,iVar,iEl) + &
                                 myPoly % dMatrix % hostData(ii,j)*f(2,i,ii,iVar,iEl)
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE VectorDivergence_2D_cpu

  SUBROUTINE VectorDivergence_2D_gpu(myPoly,f_dev,dF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(out)    :: dF_dev
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    CALL VectorDivergence_2D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                         f_dev,dF_dev,myPoly % N, &
                                         nVariables,nElements)

  END SUBROUTINE VectorDivergence_2D_gpu

  SUBROUTINE VectorCurl_2D_cpu(myPoly,f,dF,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: dF(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % N
          DO i = 0,myPoly % N

            dF(i,j,iVar,iEl) = 0.0_prec
            DO ii = 0,myPoly % N
              dF(i,j,iVar,iEl) = dF(i,j,iVar,iEl) + myPoly % dMatrix % hostData(ii,j)*f(1,i,ii,iVar,iEl) - &
                                 myPoly % dMatrix % hostData(ii,i)*f(2,ii,j,iVar,iEl)
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE VectorCurl_2D_cpu

  SUBROUTINE VectorCurl_2D_gpu(myPoly,f_dev,dF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(out)    :: dF_dev
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    CALL VectorCurl_2D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                   f_dev,dF_dev,myPoly % N, &
                                   nVariables,nElements)

  END SUBROUTINE VectorCurl_2D_gpu

  SUBROUTINE TensorDivergence_2D_cpu(myPoly,f,dF,nVariables,nElements)
    ! Note that the divergence is taken over the first dimension (row dimension) of the tensor matrix
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:2,1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: dF(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % N
          DO i = 0,myPoly % N

            dF(1,i,j,iVar,iEl) = 0.0_prec
            dF(2,i,j,iVar,iEl) = 0.0_prec
            DO ii = 0,myPoly % N
              dF(1,i,j,iVar,iEl) = dF(1,i,j,iVar,iEl) + myPoly % dMatrix % hostData(ii,i)*f(1,1,ii,j,iVar,iEl) + &
                                   myPoly % dMatrix % hostData(ii,j)*f(2,1,i,ii,iVar,iEl)
              dF(2,i,j,iVar,iEl) = dF(2,i,j,iVar,iEl) + myPoly % dMatrix % hostData(ii,i)*f(1,2,ii,j,iVar,iEl) + &
                                   myPoly % dMatrix % hostData(ii,j)*f(2,2,i,ii,iVar,iEl)
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE TensorDivergence_2D_cpu

  SUBROUTINE TensorDivergence_2D_gpu(myPoly,f_dev,dF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(out)    :: dF_dev
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    CALL TensorDivergence_2D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                         f_dev,dF_dev,myPoly % N, &
                                         nVariables,nElements)

  END SUBROUTINE TensorDivergence_2D_gpu

  SUBROUTINE ScalarGradient_3D_cpu(myPoly,f,gradF,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: gradF(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER    :: i,j,k,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO k = 0,myPoly % N
          DO j = 0,myPoly % N
            DO i = 0,myPoly % N

              gradF(1,i,j,k,iVar,iEl) = 0.0_prec
              gradF(2,i,j,k,iVar,iEl) = 0.0_prec
              gradF(3,i,j,k,iVar,iEl) = 0.0_prec
              DO ii = 0,myPoly % N
                gradF(1,i,j,k,iVar,iEl) = gradF(1,i,j,k,iVar,iEl) + myPoly % dMatrix % hostData(ii,i)*f(ii,j,k,iVar,iEl)
                gradF(2,i,j,k,iVar,iEl) = gradF(2,i,j,k,iVar,iEl) + myPoly % dMatrix % hostData(ii,j)*f(i,ii,k,iVar,iEl)
                gradF(3,i,j,k,iVar,iEl) = gradF(3,i,j,k,iVar,iEl) + myPoly % dMatrix % hostData(ii,k)*f(i,j,ii,iVar,iEl)
              END DO

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ScalarGradient_3D_cpu

  SUBROUTINE ScalarGradient_3D_gpu(myPoly,f_dev,gradF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(out)    :: gradF_dev
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    CALL ScalarGradient_3D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                       f_dev,gradF_dev,myPoly % N, &
                                       nVariables,nElements)

  END SUBROUTINE ScalarGradient_3D_gpu
!
  SUBROUTINE VectorGradient_3D_cpu(myPoly,f,gradF,nVariables,nElements)
    !
    ! Input : Vector(1:3,...)
    ! Output : Tensor(1:3,1:3,....)
    !          > Tensor(1,1) = d/ds1( Vector(1,...) )
    !          > Tensor(2,1) = d/ds1( Vector(2,...) )
    !          > Tensor(3,1) = d/ds1( Vector(3,...) )
    !          > Tensor(1,2) = d/ds2( Vector(1,...) )
    !          > Tensor(2,2) = d/ds2( Vector(2,...) )
    !          > Tensor(3,2) = d/ds2( Vector(3,...) )
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: gradF(1:3,1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER    :: i,j,k,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO k = 0,myPoly % N
          DO j = 0,myPoly % N
            DO i = 0,myPoly % N

              gradF(1,1,i,j,k,iVar,iEl) = 0.0_prec
              gradF(2,1,i,j,k,iVar,iEl) = 0.0_prec
              gradF(3,1,i,j,k,iVar,iEl) = 0.0_prec
              gradF(1,2,i,j,k,iVar,iEl) = 0.0_prec
              gradF(2,2,i,j,k,iVar,iEl) = 0.0_prec
              gradF(3,2,i,j,k,iVar,iEl) = 0.0_prec
              gradF(1,3,i,j,k,iVar,iEl) = 0.0_prec
              gradF(2,3,i,j,k,iVar,iEl) = 0.0_prec
              gradF(3,3,i,j,k,iVar,iEl) = 0.0_prec
              DO ii = 0,myPoly % N
                gradF(1,1,i,j,k,iVar,iEl) = gradF(1,1,i,j,k,iVar,iEl) + &
                                            myPoly % dMatrix % hostData(ii,i)*f(1,ii,j,k,iVar,iEl)

                gradF(2,1,i,j,k,iVar,iEl) = gradF(2,1,i,j,k,iVar,iEl) + &
                                            myPoly % dMatrix % hostData(ii,i)*f(2,ii,j,k,iVar,iEl)

                gradF(3,1,i,j,k,iVar,iEl) = gradF(3,1,i,j,k,iVar,iEl) + &
                                            myPoly % dMatrix % hostData(ii,i)*f(3,ii,j,k,iVar,iEl)

                gradF(1,2,i,j,k,iVar,iEl) = gradF(1,2,i,j,k,iVar,iEl) + &
                                            myPoly % dMatrix % hostData(ii,j)*f(1,i,ii,k,iVar,iEl)

                gradF(2,2,i,j,k,iVar,iEl) = gradF(2,2,i,j,k,iVar,iEl) + &
                                            myPoly % dMatrix % hostData(ii,j)*f(2,i,ii,k,iVar,iEl)

                gradF(3,2,i,j,k,iVar,iEl) = gradF(3,2,i,j,k,iVar,iEl) + &
                                            myPoly % dMatrix % hostData(ii,j)*f(3,i,ii,k,iVar,iEl)

                gradF(1,3,i,j,k,iVar,iEl) = gradF(1,3,i,j,k,iVar,iEl) + &
                                            myPoly % dMatrix % hostData(ii,k)*f(1,i,j,ii,iVar,iEl)

                gradF(2,3,i,j,k,iVar,iEl) = gradF(2,3,i,j,k,iVar,iEl) + &
                                            myPoly % dMatrix % hostData(ii,k)*f(2,i,j,ii,iVar,iEl)

                gradF(3,3,i,j,k,iVar,iEl) = gradF(3,3,i,j,k,iVar,iEl) + &
                                            myPoly % dMatrix % hostData(ii,k)*f(3,i,j,ii,iVar,iEl)

              END DO

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE VectorGradient_3D_cpu

  SUBROUTINE VectorGradient_3D_gpu(myPoly,f_dev,gradF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(out)    :: gradF_dev
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    CALL VectorGradient_3D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                       f_dev,gradF_dev,myPoly % N, &
                                       nVariables,nElements)

  END SUBROUTINE VectorGradient_3D_gpu

  SUBROUTINE VectorDivergence_3D_cpu(myPoly,f,dF,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: dF(0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER    :: i,j,k,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO k = 0,myPoly % N
          DO j = 0,myPoly % N
            DO i = 0,myPoly % N

              dF(i,j,k,iVar,iEl) = 0.0_prec
              DO ii = 0,myPoly % N
                dF(i,j,k,iVar,iEl) = dF(i,j,k,iVar,iEl) + myPoly % dMatrix % hostData(ii,i)*f(1,ii,j,k,iVar,iEl) + &
                                     myPoly % dMatrix % hostData(ii,j)*f(2,i,ii,k,iVar,iEl) + &
                                     myPoly % dMatrix % hostData(ii,k)*f(3,i,j,ii,iVar,iEl)
              END DO

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE VectorDivergence_3D_cpu

  SUBROUTINE VectorDivergence_3D_gpu(myPoly,f_dev,dF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(out)    :: dF_dev
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    CALL VectorDivergence_3D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                         f_dev,dF_dev,myPoly % N, &
                                         nVariables,nElements)

  END SUBROUTINE VectorDivergence_3D_gpu

  SUBROUTINE VectorCurl_3D_cpu(myPoly,f,dF,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: dF(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER    :: i,j,k,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO k = 0,myPoly % N
          DO j = 0,myPoly % N
            DO i = 0,myPoly % N

              dF(1,i,j,k,iVar,iEl) = 0.0_prec
              dF(2,i,j,k,iVar,iEl) = 0.0_prec
              dF(3,i,j,k,iVar,iEl) = 0.0_prec
              DO ii = 0,myPoly % N
                dF(1,i,j,k,iVar,iEl) = dF(1,i,j,k,iVar,iEl) + myPoly % dMatrix % hostData(ii,j)*f(3,i,ii,k,iVar,iEl) - &
                                       myPoly % dMatrix % hostData(ii,k)*f(2,i,j,ii,iVar,iEl)
                dF(2,i,j,k,iVar,iEl) = dF(2,i,j,k,iVar,iEl) + myPoly % dMatrix % hostData(ii,k)*f(1,i,j,ii,iVar,iEl) - &
                                       myPoly % dMatrix % hostData(ii,i)*f(3,ii,j,k,iVar,iEl)
                dF(3,i,j,k,iVar,iEl) = dF(3,i,j,k,iVar,iEl) + myPoly % dMatrix % hostData(ii,i)*f(2,ii,j,k,iVar,iEl) - &
                                       myPoly % dMatrix % hostData(ii,j)*f(1,i,ii,k,iVar,iEl)
              END DO

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE VectorCurl_3D_cpu

  SUBROUTINE VectorCurl_3D_gpu(myPoly,f_dev,dF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(out)    :: dF_dev
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    CALL VectorCurl_3D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                   f_dev,dF_dev,myPoly % N, &
                                   nVariables,nElements)

  END SUBROUTINE VectorCurl_3D_gpu

  SUBROUTINE TensorDivergence_3D_cpu(myPoly,f,dF,nVariables,nElements)
    ! Note that the divergence is taken over the first dimension (row dimension) of the tensor matrix
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:3,1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: dF(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER    :: i,j,k,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO k = 0,myPoly % N
          DO j = 0,myPoly % N
            DO i = 0,myPoly % N

              dF(1,i,j,k,iVar,iEl) = 0.0_prec
              dF(2,i,j,k,iVar,iEl) = 0.0_prec
              dF(3,i,j,k,iVar,iEl) = 0.0_prec
              DO ii = 0,myPoly % N
                dF(1,i,j,k,iVar,iEl) = dF(1,i,j,k,iVar,iEl) + &
                                       myPoly % dMatrix % hostData(ii,i)*f(1,1,ii,j,k,iVar,iEl) + &
                                       myPoly % dMatrix % hostData(ii,j)*f(2,1,i,ii,k,iVar,iEl) + &
                                       myPoly % dMatrix % hostData(ii,k)*f(3,1,i,j,ii,iVar,iEl)

                dF(2,i,j,k,iVar,iEl) = dF(2,i,j,k,iVar,iEl) + &
                                       myPoly % dMatrix % hostData(ii,i)*f(1,2,ii,j,k,iVar,iEl) + &
                                       myPoly % dMatrix % hostData(ii,j)*f(2,2,i,ii,k,iVar,iEl) + &
                                       myPoly % dMatrix % hostData(ii,k)*f(3,2,i,j,ii,iVar,iEl)

                dF(3,i,j,k,iVar,iEl) = dF(3,i,j,k,iVar,iEl) + &
                                       myPoly % dMatrix % hostData(ii,i)*f(1,3,ii,j,k,iVar,iEl) + &
                                       myPoly % dMatrix % hostData(ii,j)*f(2,3,i,ii,k,iVar,iEl) + &
                                       myPoly % dMatrix % hostData(ii,k)*f(3,3,i,j,ii,iVar,iEl)
              END DO

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE TensorDivergence_3D_cpu

  SUBROUTINE TensorDivergence_3D_gpu(myPoly,f_dev,dF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(out)    :: dF_dev
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl

    CALL TensorDivergence_3D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                         f_dev,dF_dev,myPoly % N, &
                                         nVariables,nElements)

  END SUBROUTINE TensorDivergence_3D_gpu

  ! /////////////////////////////// !
  ! Boundary Interpolation Routines !

  SUBROUTINE ScalarBoundaryInterp_1D_cpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    REAL(prec),INTENT(in)      :: f(0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out)     :: fBound(1:nVariables,1:2,1:nElements)
    ! Local
    INTEGER :: i,j,ii,iVar,iEl
    REAL(prec) :: fb(1:2)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        fb(1:2) = 0.0_prec
        DO ii = 0,myPoly % N
          fb(1) = fb(1) + myPoly % bMatrix % hostData(ii,0)*f(ii,iVar,iEl) ! West
          fb(2) = fb(2) + myPoly % bMatrix % hostData(ii,1)*f(ii,iVar,iEl) ! East
        END DO
        fBound(iVar,1:2,iEl) = fb(1:2)
      END DO
    END DO

  END SUBROUTINE ScalarBoundaryInterp_1D_cpu

  SUBROUTINE ScalarBoundaryInterp_1D_gpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)  :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f
    TYPE(c_ptr),INTENT(out)  :: fBound

    CALL ScalarBoundaryInterp_1D_gpu_wrapper(myPoly % bMatrix % deviceData, &
                                             f,fBound,myPoly % N,nVariables,nElements)

  END SUBROUTINE ScalarBoundaryInterp_1D_gpu

  SUBROUTINE ScalarBoundaryInterp_2D_cpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    REAL(prec),INTENT(in)      :: f(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out)     :: fBound(0:myPoly % N,1:nVariables,1:4,1:nElements)
    ! Local
    INTEGER :: i,j,ii,iVar,iEl
    REAL(prec) :: fb(1:4)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO i = 0,myPoly % N

          fb(1:4) = 0.0_prec

          DO ii = 0,myPoly % N
            fb(1) = fb(1) + myPoly % bMatrix % hostData(ii,0)*f(i,ii,iVar,iEl) ! South
            fb(2) = fb(2) + myPoly % bMatrix % hostData(ii,1)*f(ii,i,iVar,iEl) ! East
            fb(3) = fb(3) + myPoly % bMatrix % hostData(ii,1)*f(i,ii,iVar,iEl) ! North
            fb(4) = fb(4) + myPoly % bMatrix % hostData(ii,0)*f(ii,i,iVar,iEl) ! West
          END DO

          fBound(i,iVar,1:4,iEl) = fb(1:4)

        END DO
      END DO
    END DO

  END SUBROUTINE ScalarBoundaryInterp_2D_cpu

  SUBROUTINE ScalarBoundaryInterp_2D_gpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)  :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f
    TYPE(c_ptr),INTENT(out)  :: fBound

    CALL ScalarBoundaryInterp_2D_gpu_wrapper(myPoly % bMatrix % deviceData, &
                                             f,fBound,myPoly % N,nVariables,nElements)

  END SUBROUTINE ScalarBoundaryInterp_2D_gpu

  SUBROUTINE VectorBoundaryInterp_2D_cpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)  :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out)  :: fBound(1:2,0:myPoly % N,1:nVariables,1:4,1:nElements)
    ! Local
    INTEGER :: i,j,ii,idir,iVar,iEl
    REAL(prec) :: fb(1:2,1:4)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO i = 0,myPoly % N

          fb(1:2,1:4) = 0.0_prec
          DO ii = 0,myPoly % N
            DO idir = 1,2
              fb(idir,1) = fb(idir,1) + myPoly % bMatrix % hostData(ii,0)*f(idir,i,ii,iVar,iEl) ! South
              fb(idir,2) = fb(idir,2) + myPoly % bMatrix % hostData(ii,1)*f(idir,ii,i,iVar,iEl) ! East
              fb(idir,3) = fb(idir,3) + myPoly % bMatrix % hostData(ii,1)*f(idir,i,ii,iVar,iEl) ! North
              fb(idir,4) = fb(idir,4) + myPoly % bMatrix % hostData(ii,0)*f(idir,ii,i,iVar,iEl) ! West
            END DO
          END DO

          DO idir = 1,2
            fBound(idir,i,iVar,1:4,iEl) = fb(idir,1:4)
          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE VectorBoundaryInterp_2D_cpu

  SUBROUTINE VectorBoundaryInterp_2D_gpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)  :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f
    TYPE(c_ptr),INTENT(out)  :: fBound

    CALL VectorBoundaryInterp_2D_gpu_wrapper(myPoly % bMatrix % deviceData, &
                                             f,fBound,myPoly % N,nVariables,nElements)

  END SUBROUTINE VectorBoundaryInterp_2D_gpu

  SUBROUTINE TensorBoundaryInterp_2D_cpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)  :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:2,1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out)  :: fBound(1:2,1:2,0:myPoly % N,1:nVariables,1:4,1:nElements)
    ! Local
    INTEGER :: i,j,ii,idir,jdir,iVar,iEl
    REAL(prec) :: fb(1:2,1:2,1:4)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO i = 0,myPoly % N

          fb(1:2,1:2,1:4) = 0.0_prec
          DO ii = 0,myPoly % N
            DO jdir = 1,2
              DO idir = 1,2
                fb(idir,jdir,1) = fb(idir,jdir,1) + myPoly % bMatrix % hostData(ii,0)*f(idir,jdir,i,ii,iVar,iEl) ! South
                fb(idir,jdir,2) = fb(idir,jdir,2) + myPoly % bMatrix % hostData(ii,1)*f(idir,jdir,ii,i,iVar,iEl) ! East
                fb(idir,jdir,3) = fb(idir,jdir,3) + myPoly % bMatrix % hostData(ii,1)*f(idir,jdir,i,ii,iVar,iEl) ! North
                fb(idir,jdir,4) = fb(idir,jdir,4) + myPoly % bMatrix % hostData(ii,0)*f(idir,jdir,ii,i,iVar,iEl) ! West
              END DO
            END DO
          END DO

          DO jdir = 1,2
            DO idir = 1,2
              fBound(idir,jdir,i,iVar,1:4,iEl) = fb(idir,jdir,1:4)
            END DO
          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE TensorBoundaryInterp_2D_cpu

  SUBROUTINE TensorBoundaryInterp_2D_gpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)  :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f
    TYPE(c_ptr),INTENT(out)  :: fBound

    CALL TensorBoundaryInterp_2D_gpu_wrapper(myPoly % bMatrix % deviceData, &
                                             f,fBound,myPoly % N,nVariables,nElements)

  END SUBROUTINE TensorBoundaryInterp_2D_gpu

  SUBROUTINE ScalarBoundaryInterp_3D_cpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    REAL(prec),INTENT(in)      :: f(0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out)     :: fBound(0:myPoly % N,0:myPoly % N,1:nVariables,1:6,1:nElements)
    ! Local
    INTEGER :: i,j,ii,iVar,iEl
    REAL(prec) :: fb(1:6)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % N
          DO i = 0,myPoly % N

            fb(1:6) = 0.0_prec

            DO ii = 0,myPoly % N
              fb(1) = fb(1) + myPoly % bMatrix % hostData(ii,0)*f(i,ii,j,iVar,iEl) ! South
              fb(2) = fb(2) + myPoly % bMatrix % hostData(ii,1)*f(ii,i,j,iVar,iEl) ! East
              fb(3) = fb(3) + myPoly % bMatrix % hostData(ii,1)*f(i,ii,j,iVar,iEl) ! North
              fb(4) = fb(4) + myPoly % bMatrix % hostData(ii,0)*f(ii,i,j,iVar,iEl) ! West
              fb(5) = fb(5) + myPoly % bMatrix % hostData(ii,0)*f(i,j,ii,iVar,iEl) ! Bottom
              fb(6) = fb(6) + myPoly % bMatrix % hostData(ii,1)*f(i,j,ii,iVar,iEl) ! Top
            END DO

            fBound(i,j,iVar,1:6,iEl) = fb(1:6)

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ScalarBoundaryInterp_3D_cpu

  SUBROUTINE ScalarBoundaryInterp_3D_gpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)  :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f
    TYPE(c_ptr),INTENT(out)  :: fBound

    CALL ScalarBoundaryInterp_3D_gpu_wrapper(myPoly % bMatrix % deviceData, &
                                             f,fBound,myPoly % N,nVariables,nElements)

  END SUBROUTINE ScalarBoundaryInterp_3D_gpu

  SUBROUTINE VectorBoundaryInterp_3D_cpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)  :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out)  :: fBound(1:3,0:myPoly % N,0:myPoly % N,1:nVariables,1:6,1:nElements)
    ! Local
    INTEGER :: i,j,ii,idir,iVar,iEl
    REAL(prec) :: fb(1:3,1:6)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % N
          DO i = 0,myPoly % N

            fb(1:3,1:6) = 0.0_prec
            DO ii = 0,myPoly % N
              DO idir = 1,3
                fb(idir,1) = fb(idir,1) + myPoly % bMatrix % hostData(ii,0)*f(idir,i,ii,j,iVar,iEl) ! South
                fb(idir,2) = fb(idir,2) + myPoly % bMatrix % hostData(ii,1)*f(idir,ii,i,j,iVar,iEl) ! East
                fb(idir,3) = fb(idir,3) + myPoly % bMatrix % hostData(ii,1)*f(idir,i,ii,j,iVar,iEl) ! North
                fb(idir,4) = fb(idir,4) + myPoly % bMatrix % hostData(ii,0)*f(idir,ii,i,j,iVar,iEl) ! West
                fb(idir,5) = fb(idir,5) + myPoly % bMatrix % hostData(ii,0)*f(idir,i,j,ii,iVar,iEl) ! Bottom
                fb(idir,6) = fb(idir,6) + myPoly % bMatrix % hostData(ii,1)*f(idir,i,j,ii,iVar,iEl) ! Top
              END DO
            END DO

            DO idir = 1,3
              fBound(idir,i,j,iVar,1:6,iEl) = fb(idir,1:6)
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE VectorBoundaryInterp_3D_cpu

  SUBROUTINE VectorBoundaryInterp_3D_gpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)  :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f
    TYPE(c_ptr),INTENT(out)  :: fBound

    CALL VectorBoundaryInterp_3D_gpu_wrapper(myPoly % bMatrix % deviceData, &
                                             f,fBound,myPoly % N,nVariables,nElements)

  END SUBROUTINE VectorBoundaryInterp_3D_gpu

  SUBROUTINE TensorBoundaryInterp_3D_cpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)  :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:3,1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out)  :: fBound(1:3,1:3,0:myPoly % N,0:myPoly % N,1:nVariables,1:6,1:nElements)
    ! Local
    INTEGER :: i,j,ii,idir,jdir,iVar,iEl
    REAL(prec) :: fb(1:3,1:3,1:6)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % N
          DO i = 0,myPoly % N

            fb(1:3,1:3,1:6) = 0.0_prec
            DO ii = 0,myPoly % N
              DO jdir = 1,3
                DO idir = 1,3
                  fb(idir,jdir,1) = fb(idir,jdir,1) + myPoly % bMatrix % hostData(ii,0)*f(idir,jdir,i,ii,j,iVar,iEl) ! South
                  fb(idir,jdir,2) = fb(idir,jdir,2) + myPoly % bMatrix % hostData(ii,1)*f(idir,jdir,ii,i,j,iVar,iEl) ! East
                  fb(idir,jdir,3) = fb(idir,jdir,3) + myPoly % bMatrix % hostData(ii,1)*f(idir,jdir,i,ii,j,iVar,iEl) ! North
                  fb(idir,jdir,4) = fb(idir,jdir,4) + myPoly % bMatrix % hostData(ii,0)*f(idir,jdir,ii,i,j,iVar,iEl) ! West
                  fb(idir,jdir,5) = fb(idir,jdir,5) + myPoly % bMatrix % hostData(ii,0)*f(idir,jdir,i,j,ii,iVar,iEl) ! Bottom
                  fb(idir,jdir,6) = fb(idir,jdir,6) + myPoly % bMatrix % hostData(ii,1)*f(idir,jdir,i,j,ii,iVar,iEl) ! Top
                END DO
              END DO
            END DO

            DO jdir = 1,3
              DO idir = 1,3
                fBound(idir,jdir,i,j,iVar,1:6,iEl) = fb(idir,jdir,1:6)
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE TensorBoundaryInterp_3D_cpu

  SUBROUTINE TensorBoundaryInterp_3D_gpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)  :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)  :: f
    TYPE(c_ptr),INTENT(out)  :: fBound

    CALL TensorBoundaryInterp_3D_gpu_wrapper(myPoly % bMatrix % deviceData, &
                                             f,fBound,myPoly % N,nVariables,nElements)

  END SUBROUTINE TensorBoundaryInterp_3D_gpu
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

  SUBROUTINE CalculateBarycentricWeights(myPoly)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(inout) :: myPoly
    ! Local
    INTEGER :: i,j

    DO i = 0,myPoly % N
      myPoly % bWeights % hostData(i) = 1.0_prec
    END DO

    ! Computes the product w_k = w_k*(s_k - s_j), k /= j
    DO j = 1,myPoly % N
      DO i = 0,j - 1

        myPoly % bWeights % hostData(i) = myPoly % bWeights % hostData(i)* &
                                          (myPoly % controlPoints % hostData(i) - myPoly % controlPoints % hostData(j))
        myPoly % bWeights % hostData(j) = myPoly % bWeights % hostData(j)* &
                                          (myPoly % controlPoints % hostData(j) - myPoly % controlPoints % hostData(i))

      END DO
    END DO

    DO j = 0,myPoly % N
      myPoly % bWeights % hostData(j) = 1.0_prec/myPoly % bWeights % hostData(j)
    END DO

  END SUBROUTINE CalculateBarycentricWeights

! ================================================================================================ !
!
! CalculateInterpolationMatrix (PRIVATE)
!
!   A PRIVATE routine that fills in the interpolation matrix for the Lagrange data structure.
!
!   This function is from Alg. 32 on pg. 76 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  SUBROUTINE CalculateInterpolationMatrix(myPoly)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(inout) :: myPoly
    ! Local
    REAL(prec) :: temp1,temp2
    INTEGER    :: row,col
    LOGICAL    :: rowHasMatch
    REAL(prec) :: iMatrix(0:myPoly % M,0:myPoly % N)

    DO row = 0,myPoly % M

      rowHasMatch = .FALSE.

      DO col = 0,myPoly % N

        iMatrix(row,col) = 0.0_prec

        IF (AlmostEqual(myPoly % targetPoints % hostData(row),myPoly % controlPoints % hostData(col))) THEN
          rowHasMatch = .TRUE.
          iMatrix(row,col) = 1.0_prec
        END IF

      END DO

      IF (.NOT. (rowHasMatch)) THEN

        temp1 = 0.0_prec

        DO col = 0,myPoly % N
          temp2 = myPoly % bWeights % hostData(col)/ &
                  (myPoly % targetPoints % hostData(row) - &
                   myPoly % controlPoints % hostData(col))
          iMatrix(row,col) = temp2
          temp1 = temp1 + temp2
        END DO

        DO col = 0,myPoly % N
          iMatrix(row,col) = iMatrix(row,col)/temp1
        END DO

      END IF

    END DO

    myPoly % iMatrix % hostData = TRANSPOSE(iMatrix)

  END SUBROUTINE CalculateInterpolationMatrix

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

  SUBROUTINE CalculateDerivativeMatrix(myPoly)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(inout) :: myPoly
    ! Local
    INTEGER    :: row,col
    REAL(prec) :: dmat(0:myPoly % N,0:myPoly % N)
    REAL(prec) :: dgmat(0:myPoly % N,0:myPoly % N)

    DO row = 0,myPoly % N

      dmat(row,row) = 0.0_prec

      DO col = 0,myPoly % N

        IF (.NOT. (col == row)) THEN

          dmat(row,col) = myPoly % bWeights % hostData(col)/ &
                          (myPoly % bWeights % hostData(row)* &
                           (myPoly % controlPoints % hostData(row) - &
                            myPoly % controlPoints % hostData(col)))

          dmat(row,row) = dmat(row,row) - dmat(row,col)

        END IF

      END DO

    END DO

    DO row = 0,myPoly % N
      DO col = 0,myPoly % N
        dgmat(row,col) = -dmat(col,row)*&
                          myPoly % qWeights % hostData(col)/&
                          myPoly % qWeights % hostData(row)
      ENDDO
    ENDDO

    myPoly % dMatrix % hostData = TRANSPOSE(dmat)
    myPoly % dgMatrix % hostData = TRANSPOSE(dgmat)

  END SUBROUTINE CalculateDerivativeMatrix

! ================================================================================================ !
!
! CalculateLagrangePolynomials
!
!   Evaluates each of the 1-D Lagrange interpolating polynomials at a specified point.
!
!   This function is from Alg. 34 on pg. 77 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  FUNCTION CalculateLagrangePolynomials(myPoly,sE) RESULT(lAtS)
    IMPLICIT NONE
    CLASS(Lagrange) :: myPoly
    REAL(prec)      :: sE
    REAL(prec)      :: lAtS(0:myPoly % N)
    ! Local
    REAL(prec) :: temp1,temp2
    INTEGER    :: j
    LOGICAL    :: xMatchesNode

    xMatchesNode = .FALSE.

    DO j = 0,myPoly % N

      lAtS(j) = 0.0_prec

      IF (AlmostEqual(sE,myPoly % controlPoints % hostData(j))) THEN
        lAtS(j) = 1.0_prec
        xMatchesNode = .TRUE.
      END IF

    END DO

    IF (xMatchesNode) THEN
      RETURN
    END IF

    temp1 = 0.0_prec

    DO j = 0,myPoly % N
      temp2 = myPoly % bWeights % hostData(j)/(sE - myPoly % controlPoints % hostData(j))
      lAtS(j) = temp2
      temp1 = temp1 + temp2
    END DO

    lAtS = lAtS/temp1

  END FUNCTION CalculateLagrangePolynomials

END MODULE SELF_Lagrange
