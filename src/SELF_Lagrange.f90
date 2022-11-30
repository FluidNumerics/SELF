! SELF_Lagrange.f90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE SELF_Lagrange

  USE ISO_FORTRAN_ENV
  USE SELF_Constants
  USE SELF_Memory
  USE SELF_SupportRoutines
  USE SELF_Quadrature

  USE hipfort
  USE hipfort_check

  USE ISO_C_BINDING
  IMPLICIT NONE

  TYPE,PUBLIC :: Lagrange
    !! A data structure for working with Lagrange Interpolating Polynomials in one, two, and three dimensions.
    !! The Lagrange data-structure stores the information necessary to interpolate between two
    !! sets of grid-points and to estimate the derivative of data at native grid points. Routines for
    !! multidimensional interpolation are based on the tensor product of 1-D interpolants. It is
    !! assumed that the polynomial degree (and the interpolation nodes) are the same in each direction.
    !! This assumption permits the storage of only one array of interpolation nodes and barycentric
    !! weights and is what allows this data structure to be flexible.

    INTEGER :: N
      !! The number of control points.

    INTEGER :: controlNodeType

    INTEGER :: M
      !! The number of target points.

    INTEGER :: targetNodeType

    TYPE(hfReal_r1) :: controlPoints
      !! The set of nodes in one dimension where data is known.
      !! To create higher dimension interpolation and differentiation operators, structured grids in two and three
      !! dimensions are created by tensor products of the controlPoints. This design decision implies that all
      !! spectral element methods supported by the Lagrange class have the same polynomial degree in each
      !! computational/spatial dimension. In practice, the controlPoints are the Legendre-Gauss, Legendre-Gauss-Lobatto,
      !! Legendre-Gauss-Radau, Chebyshev-Gauss, Chebyshev-Gauss-Lobatto, or Chebyshev-Gauss-Radau quadrature points over
      !! the domain [-1,1] (computational space). The Init routine for this class restricts controlPoints to one of
      !! these quadrature types or uniform points on [-1,1].

    TYPE(hfReal_r1) :: targetPoints
      !! The set of nodes in one dimension where data is to be interpolated to. To create higher dimension interpolation
      !! and differentiation operators, structured grids in two and three dimensions are created by tensor products of
      !! the targetPoints. In practice, the targetPoints are set to a uniformly distributed set of points between [-1,1]
      !! (computational space) to allow for interpolation from unevenly spaced quadrature points to a plotting grid.

    TYPE(hfReal_r1) :: bWeights
      !! The barycentric weights that are calculated from the controlPoints and used for interpolation.

    TYPE(hfReal_r1) :: qWeights
      !! The quadrature weights for discrete integration. The quadradture weights depend on the type of controlPoints
      !! provided; one of Legendre-Gauss, Legendre-Gauss-Lobatto, Legendre-Gauss-Radau, Chebyshev-Gauss,
      !! Chebyshev-Gauss-Lobatto, Chebyshev-Gauss Radau, or Uniform. If Uniform, the quadrature weights are constant
      !! $$dx = \frac{2.0}{N+1}$$.

    TYPE(hfReal_r2) :: iMatrix
      !! The interpolation matrix (transpose) for mapping data from the control grid to the target grid.

    TYPE(hfReal_r2) :: dMatrix
      !! The derivative matrix for mapping function nodal values to a nodal values of the derivative estimate. The
      !! dMatrix is based on a strong form of the derivative.

    TYPE(hfReal_r2) :: dgMatrix
      !! The derivative matrix for mapping function nodal values to a nodal values of the derivative estimate. The dgMatrix is based
      !! on a weak form of the derivative. It must be used with bMatrix to account for boundary contributions in the weak form.

    TYPE(hfReal_r2) :: bMatrix
      !! The boundary interpolation matrix that is used to map a grid of nodal values at the control points to the element boundaries.

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Lagrange
    PROCEDURE,PUBLIC :: Free => Free_Lagrange

    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Lagrange
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Lagrange

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

    GENERIC,PUBLIC :: ScalarDGGradient_2D => ScalarDGGradient_2D_cpu,ScalarDGGradient_2D_gpu
    PROCEDURE,PRIVATE :: ScalarDGGradient_2D_cpu,ScalarDGGradient_2D_gpu

    GENERIC,PUBLIC :: VectorGradient_2D => VectorGradient_2D_cpu,VectorGradient_2D_gpu
    PROCEDURE,PRIVATE :: VectorGradient_2D_cpu,VectorGradient_2D_gpu

    GENERIC,PUBLIC :: VectorDGGradient_2D => VectorDGGradient_2D_cpu,VectorDGGradient_2D_gpu
    PROCEDURE,PRIVATE :: VectorDGGradient_2D_cpu,VectorDGGradient_2D_gpu

    GENERIC,PUBLIC :: VectorDivergence_2D => VectorDivergence_2D_cpu,VectorDivergence_2D_gpu
    PROCEDURE,PRIVATE :: VectorDivergence_2D_cpu,VectorDivergence_2D_gpu

    GENERIC,PUBLIC :: VectorDGDivergence_2D => VectorDGDivergence_2D_cpu,VectorDGDivergence_2D_gpu
    PROCEDURE,PRIVATE :: VectorDGDivergence_2D_cpu,VectorDGDivergence_2D_gpu

    GENERIC,PUBLIC :: VectorCurl_2D => VectorCurl_2D_cpu,VectorCurl_2D_gpu
    PROCEDURE,PRIVATE :: VectorCurl_2D_cpu,VectorCurl_2D_gpu

    GENERIC,PUBLIC :: TensorDivergence_2D => TensorDivergence_2D_cpu,TensorDivergence_2D_gpu
    PROCEDURE,PRIVATE :: TensorDivergence_2D_cpu,TensorDivergence_2D_gpu

    GENERIC,PUBLIC :: TensorDGDivergence_2D => TensorDGDivergence_2D_cpu,TensorDGDivergence_2D_gpu
    PROCEDURE,PRIVATE :: TensorDGDivergence_2D_cpu,TensorDGDivergence_2D_gpu

    GENERIC,PUBLIC :: ScalarGradient_3D => ScalarGradient_3D_cpu,ScalarGradient_3D_gpu
    PROCEDURE,PRIVATE :: ScalarGradient_3D_cpu,ScalarGradient_3D_gpu

    GENERIC,PUBLIC :: VectorGradient_3D => VectorGradient_3D_cpu,VectorGradient_3D_gpu
    PROCEDURE,PRIVATE :: VectorGradient_3D_cpu,VectorGradient_3D_gpu

    GENERIC,PUBLIC :: VectorDivergence_3D => VectorDivergence_3D_cpu,VectorDivergence_3D_gpu
    PROCEDURE,PRIVATE :: VectorDivergence_3D_cpu,VectorDivergence_3D_gpu

    GENERIC,PUBLIC :: VectorDGDivergence_3D => VectorDGDivergence_3D_cpu,VectorDGDivergence_3D_gpu
    PROCEDURE,PRIVATE :: VectorDGDivergence_3D_cpu,VectorDGDivergence_3D_gpu

    GENERIC,PUBLIC :: VectorCurl_3D => VectorCurl_3D_cpu,VectorCurl_3D_gpu
    PROCEDURE,PRIVATE :: VectorCurl_3D_cpu,VectorCurl_3D_gpu

    GENERIC,PUBLIC :: TensorDivergence_3D => TensorDivergence_3D_cpu,TensorDivergence_3D_gpu
    PROCEDURE,PRIVATE :: TensorDivergence_3D_cpu,TensorDivergence_3D_gpu

    GENERIC,PUBLIC :: TensorDGDivergence_3D => TensorDGDivergence_3D_cpu,TensorDGDivergence_3D_gpu
    PROCEDURE,PRIVATE :: TensorDGDivergence_3D_cpu,TensorDGDivergence_3D_gpu

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
      INTEGER(C_INT),VALUE :: N,M,nVar,nEl
    END SUBROUTINE ScalarGridInterp_1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarGridInterp_2D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="ScalarGridInterp_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER(C_INT),VALUE :: N,M,nVar,nEl
    END SUBROUTINE ScalarGridInterp_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorGridInterp_2D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="VectorGridInterp_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER(C_INT),VALUE :: N,M,nVar,nEl
    END SUBROUTINE VectorGridInterp_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorGridInterp_2D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="TensorGridInterp_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER(C_INT),VALUE :: N,M,nVar,nEl
    END SUBROUTINE TensorGridInterp_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarGridInterp_3D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="ScalarGridInterp_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER(C_INT),VALUE :: N,M,nVar,nEl
    END SUBROUTINE ScalarGridInterp_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorGridInterp_3D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="VectorGridInterp_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER(C_INT),VALUE :: N,M,nVar,nEl
    END SUBROUTINE VectorGridInterp_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorGridInterp_3D_gpu_wrapper(iMatrixT_dev,f_dev,fInterp_dev,N,M,nVar,nEl) &
      bind(c,name="TensorGridInterp_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: iMatrixT_dev,f_dev,fInterp_dev
      INTEGER(C_INT),VALUE :: N,M,nVar,nEl
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
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ScalarBoundaryInterp_1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="ScalarBoundaryInterp_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ScalarBoundaryInterp_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="VectorBoundaryInterp_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE VectorBoundaryInterp_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="TensorBoundaryInterp_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE TensorBoundaryInterp_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="ScalarBoundaryInterp_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ScalarBoundaryInterp_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="VectorBoundaryInterp_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE VectorBoundaryInterp_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fBound_dev,N,nVar,nEl) &
      bind(c,name="TensorBoundaryInterp_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: bMatrix_dev,f_dev,fBound_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE TensorBoundaryInterp_3D_gpu_wrapper
  END INTERFACE

  ! /////////////// !

  INTERFACE
    SUBROUTINE Derivative_1D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="Derivative_1D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE Derivative_1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE DGDerivative_1D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
      bind(c,name="DGDerivative_1D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE DGDerivative_1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarGradient_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="ScalarGradient_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ScalarGradient_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarDGGradient_2D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
      bind(c,name="ScalarDGGradient_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ScalarDGGradient_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorGradient_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorGradient_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE VectorGradient_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorDGGradient_2D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorDGGradient_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE VectorDGGradient_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorDivergence_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorDivergence_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE VectorDivergence_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
   SUBROUTINE VectorDGDivergence_2D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorDGDivergence_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE VectorDGDivergence_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorCurl_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorCurl_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE VectorCurl_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorDivergence_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="TensorDivergence_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE TensorDivergence_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorDGDivergence_2D_gpu_wrapper(dMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
      bind(c,name="TensorDGDivergence_2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE TensorDGDivergence_2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ScalarGradient_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="ScalarGradient_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ScalarGradient_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorGradient_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorGradient_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE VectorGradient_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorDivergence_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorDivergence_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE VectorDivergence_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
   SUBROUTINE VectorDGDivergence_3D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorDGDivergence_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE VectorDGDivergence_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE VectorCurl_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="VectorCurl_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE VectorCurl_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE TensorDivergence_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
      bind(c,name="TensorDivergence_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE TensorDivergence_3D_gpu_wrapper
  END INTERFACE

  INTERFACE
   SUBROUTINE TensorDGDivergence_3D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
      bind(c,name="TensorDGDivergence_3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE TensorDGDivergence_3D_gpu_wrapper
  END INTERFACE

CONTAINS

  SUBROUTINE Init_Lagrange(myPoly,N,controlNodeType,M,targetNodeType)
    !! Initialize an instance of the Lagrange class
    !! On output, all of the attributes for the Lagrange class are allocated and values are initialized according to the number of
    !! control points, number of target points, and the types for the control and target nodes.
    !! If a GPU is available, device pointers for the Lagrange attributes are allocated and initialized.
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(out) :: myPoly
    !! Lagrange class instance
    INTEGER,INTENT(in)          :: N
    !! The number of control points for interpolant
    INTEGER,INTENT(in)          :: M
    !! The number of target points for the interpolant
    INTEGER,INTENT(in)          :: controlNodeType
    !! The integer code specifying the type of control points. Parameters are defined in SELF_Constants.f90. One of GAUSS(=1),
    !! GAUSS_LOBATTO(=2), or UNIFORM(=3)
    INTEGER,INTENT(in)          :: targetNodeType
    !! The integer code specifying the type of target points. Parameters are defined in SELF_Constants.f90. One of GAUSS(=1),
    !! GAUSS_LOBATTO(=2), or UNIFORM(=3)
    ! -------!
    ! Local
    REAL(prec) :: q(0:M)

    myPoly % N = N
    myPoly % M = M
    myPoly % controlNodeType = controlNodeType
    myPoly % targetNodeType = targetNodeType

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

    ELSEIF (controlNodeType == CHEBYSHEV_GAUSS .OR. controlNodeType == CHEBYSHEV_GAUSS_LOBATTO) THEN

      CALL ChebyshevQuadrature(N, &
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

    CALL myPoly % UpdateDevice()

  END SUBROUTINE Init_Lagrange

  SUBROUTINE Free_Lagrange(myPoly)
    !! Frees all memory (host and device) associated with an instance of the Lagrange class
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(inout) :: myPoly
    !! Lagrange class instance

    CALL myPoly % controlPoints % Free()
    CALL myPoly % targetPoints % Free()
    CALL myPoly % bWeights % Free()
    CALL myPoly % qWeights % Free()
    CALL myPoly % iMatrix % Free()
    CALL myPoly % dMatrix % Free()
    CALL myPoly % dgMatrix % Free()
    CALL myPoly % bMatrix % Free()

  END SUBROUTINE Free_Lagrange

  SUBROUTINE UpdateDevice_Lagrange(myPoly)
    !! Copy the Lagrange attributes from the host (CPU) to the device (GPU)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(inout) :: myPoly
    !! Lagrange class instance

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
    !! Copy the Lagrange attributes from the device (GPU) to the host (CPU)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(inout) :: myPoly
    !! Lagrange class instance

    CALL myPoly % controlPoints % UpdateHost()
    CALL myPoly % targetPoints % UpdateHost()
    CALL myPoly % bWeights % UpdateHost()
    CALL myPoly % qWeights % UpdateHost()
    CALL myPoly % iMatrix % UpdateHost()
    CALL myPoly % dMatrix % UpdateHost()
    CALL myPoly % dgMatrix % UpdateHost()
    CALL myPoly % bMatrix % UpdateHost()

  END SUBROUTINE UpdateHost_Lagrange

  SUBROUTINE ScalarGridInterp_1D_cpu(myPoly,f,fInterp,nVariables,nElements)
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
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    !! Lagrange class instance
    INTEGER,INTENT(in)     :: nVariables
    !! The number of variables/functions that are interpolated
    INTEGER,INTENT(in)     :: nElements
    !! The number of spectral elements in the SEM grid
    REAL(prec),INTENT(in)  :: f(0:myPoly % N,1:nVariables,1:nElements)
    !! (Input) Array of function values, defined on the control grid
    REAL(prec),INTENT(out) :: fInterp(0:myPoly % M,1:nVariables,1:nElements)
    !! (Output) Array of function values, defined on the target grid
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
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    !! Lagrange class instance
    INTEGER,INTENT(in) :: nVariables
    !! The number of variables/functions that are interpolated
    INTEGER,INTENT(in) :: nElements
    !! The number of spectral elements in the SEM grid
    TYPE(c_ptr),INTENT(in)  :: f_dev
    !! (Input) Array of function values, defined on the control grid
    TYPE(c_ptr),INTENT(out) :: fInterp_dev
    !! (Output) Array of function values, defined on the target grid

    CALL ScalarGridInterp_1D_gpu_wrapper(myPoly % iMatrix % deviceData, &
                                         f_dev,fInterp_dev, &
                                         myPoly % N,myPoly % M, &
                                         nVariables,nElements)

  END SUBROUTINE ScalarGridInterp_1D_gpu

  SUBROUTINE ScalarGridInterp_2D_cpu(myPoly,f,fNew,nVariables,nElements)
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
    !! $$ \tilde{f}_{m,n,ivar,iel} = \sum_{j=0}^N \sum_{i=0}^N f_{i,j,ivar,iel} I_{i,m} I_{j,n} $$
    !! 
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    !! Lagrange class instance
    INTEGER,INTENT(in)     :: nVariables
    !! The number of variables/functions that are interpolated
    INTEGER,INTENT(in)     :: nElements
    !! The number of spectral elements in the SEM grid
    REAL(prec),INTENT(in)  :: f(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    !! (Input) Array of function values, defined on the control grid
    REAL(prec),INTENT(out) :: fNew(0:myPoly % M,0:myPoly % M,1:nVariables,1:nElements)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    INTEGER :: i,j,ii,jj,iEl,iVar
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

  SUBROUTINE ScalarGridInterp_2D_gpu(myPoly,f_dev,fInterp_dev,nVariables,nElements)
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
    !! $$ \tilde{f}_{m,n,ivar,iel} = \sum_{j=0}^N \sum_{i=0}^N f_{i,j,ivar,iel} I_{i,m} I_{j,n} $$
    !! 
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    !! Lagrange class instance
    INTEGER,INTENT(in) :: nVariables
    !! The number of variables/functions that are interpolated
    INTEGER,INTENT(in) :: nElements
    !! The number of spectral elements in the SEM grid
    TYPE(c_ptr),INTENT(in)  :: f_dev
    !! (Input) Array of function values, defined on the control grid
    TYPE(c_ptr),INTENT(out) :: fInterp_dev
    !! (Output) Array of function values, defined on the target grid

    CALL ScalarGridInterp_2D_gpu_wrapper(myPoly % iMatrix % deviceData, &
                                         f_dev,fInterp_dev, &
                                         myPoly % N,myPoly % M, &
                                         nVariables,nElements)

  END SUBROUTINE ScalarGridInterp_2D_gpu

  SUBROUTINE VectorGridInterp_2D_cpu(myPoly,f,fNew,nVariables,nElements)
    !! Host (CPU) implementation of the VectorGridInterp_2D interface.
    !! In most cases, you should use the `VectorGridInterp_2D` generic interface,
    !! rather than calling this routine directly.
    !! Interpolate a vector-2D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized 
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{dir,m,n,ivar,iel} = \sum_{j=0}^N \sum_{i=0}^N f_{dir,i,j,ivar,iel} I_{i,m} I_{j,n} $$
    !! 
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    !! Lagrange class instance
    INTEGER,INTENT(in)     :: nVariables
    !! The number of variables/functions that are interpolated
    INTEGER,INTENT(in)     :: nElements
    !! The number of spectral elements in the SEM grid
    REAL(prec),INTENT(in)  :: f(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    !! (Input) Array of function values, defined on the control grid
    REAL(prec),INTENT(out) :: fNew(1:2,0:myPoly % M,0:myPoly % M,1:nVariables,1:nElements)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    INTEGER :: i,j,ii,jj,iEl,iVar
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
    !! Device (GPU) implementation of the VectorGridInterp_2D interface.
    !! In most cases, you should use the `VectorGridInterp_2D` generic interface,
    !! rather than calling this routine directly.
    !! This routine calls hip/SELF_Lagrange.cpp:VectorGridInterp_2D_gpu_wrapper
    !! Interpolate a vector-2D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized 
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{dir,m,n,ivar,iel} = \sum_{j=0}^N \sum_{i=0}^N f_{dir,i,j,ivar,iel} I_{i,m} I_{j,n} $$
    !! 
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    !! Lagrange class instance
    INTEGER,INTENT(in) :: nVariables
    !! The number of variables/functions that are interpolated
    INTEGER,INTENT(in) :: nElements
    !! The number of spectral elements in the SEM grid
    TYPE(c_ptr),INTENT(in)  :: f_dev
    !! (Input) Array of function values, defined on the control grid
    TYPE(c_ptr),INTENT(out) :: fInterp_dev
    !! (Output) Array of function values, defined on the target grid

    CALL VectorGridInterp_2D_gpu_wrapper(myPoly % iMatrix % deviceData, &
                                         f_dev,fInterp_dev, &
                                         myPoly % N,myPoly % M, &
                                         nVariables,nElements)

  END SUBROUTINE VectorGridInterp_2D_gpu

  SUBROUTINE TensorGridInterp_2D_cpu(myPoly,f,fNew,nVariables,nElements)
    !! Host (CPU) implementation of the TensorGridInterp_2D interface.
    !! In most cases, you should use the `TensorGridInterp_2D` generic interface,
    !! rather than calling this routine directly.
    !! Interpolate a tensor-2D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized 
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-tensor multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{row,col,m,n,ivar,iel} = \sum_{j=0}^N \sum_{i=0}^N f_{row,col,i,j,ivar,iel} I_{i,m} I_{j,n} $$
    !! 
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    !! Lagrange class instance
    INTEGER,INTENT(in)     :: nVariables
    !! The number of variables/functions that are interpolated
    INTEGER,INTENT(in)     :: nElements
    !! The number of spectral elements in the SEM grid
    REAL(prec),INTENT(in)  :: f(1:2,1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    !! (Input) Array of function values, defined on the control grid
    REAL(prec),INTENT(out) :: fNew(1:2,1:2,0:myPoly % M,0:myPoly % M,1:nVariables,1:nElements)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    INTEGER :: i,j,ii,jj,iEl,iVar
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
    !! Device (GPU) implementation of the TensorGridInterp_2D interface.
    !! In most cases, you should use the `TensorGridInterp_2D` generic interface,
    !! rather than calling this routine directly.
    !! This routine calls hip/SELF_Lagrange.cpp:TensorGridInterp_2D_gpu_wrapper
    !! Interpolate a tensor-2D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized 
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{row,col,m,n,ivar,iel} = \sum_{j=0}^N \sum_{i=0}^N f_{row,col,i,j,ivar,iel} I_{i,m} I_{j,n} $$
    !! 
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    !! Lagrange class instance
    INTEGER,INTENT(in) :: nVariables
    !! The number of variables/functions that are interpolated
    INTEGER,INTENT(in) :: nElements
    !! The number of spectral elements in the SEM grid
    TYPE(c_ptr),INTENT(in)  :: f_dev
    !! (Input) Array of function values, defined on the control grid
    TYPE(c_ptr),INTENT(out) :: fInterp_dev
    !! (Output) Array of function values, defined on the target grid

    CALL TensorGridInterp_2D_gpu_wrapper(myPoly % iMatrix % deviceData, &
                                         f_dev,fInterp_dev, &
                                         myPoly % N,myPoly % M, &
                                         nVariables,nElements)

  END SUBROUTINE TensorGridInterp_2D_gpu

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
          df(i,iVar,iEl) = df(i,iVar,iEl) + (bf(iVar,2,iEl)*myPoly % bMatrix % hostData(i,1) + &
                                             bf(iVar,1,iEl)*myPoly % bMatrix % hostData(i,0))/ &
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
!     INTEGER        :: nVariables, nElements
!     REAL(prec)     :: f(0:interp % N,0:interp % N,1:nVariables,1:nElements)
!     REAL(prec)     :: gradF(1:2,0:interp % N,0:interp % N,1:nVariables,1:nElements)
!
!       CALL interp % CalculateGradient_2D( f, gradF, nVariables, nElements )
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
!     gradF (out)
!      Array of derivative values at the target interpolation nodes.
!
! ================================================================================================ !
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

    CALL ScalarGradient_2D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                       f_dev,gradF_dev,myPoly % N, &
                                       nVariables,nElements)

  END SUBROUTINE ScalarGradient_2D_gpu
!
!
  SUBROUTINE ScalarDGGradient_2D_cpu(myPoly,f,bf,gradF,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(in)  :: bf(0:myPoly % N,1:nVariables,1:4,1:nElements)
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
              gradF(1,i,j,iVar,iEl) = gradF(1,i,j,iVar,iEl) + myPoly % dgMatrix % hostData(ii,i)*f(ii,j,iVar,iEl)
              gradF(2,i,j,iVar,iEl) = gradF(2,i,j,iVar,iEl) + myPoly % dgMatrix % hostData(ii,j)*f(i,ii,iVar,iEl)
            END DO

            ! Boundary Contribution
            gradF(1,i,j,iVar,iEl) = gradF(1,i,j,iVar,iEl) + (bf(j,iVar,2,iEl)*myPoly % bMatrix % hostData(i,1) + &
                                                             bf(j,iVar,4,iEl)*myPoly % bMatrix % hostData(i,0))/ &
                                    myPoly % qWeights % hostData(i)

            gradF(2,i,j,iVar,iEl) = gradF(2,i,j,iVar,iEl) + (bf(i,iVar,3,iEl)*myPoly % bMatrix % hostData(j,1) + &
                                                             bf(i,iVar,1,iEl)*myPoly % bMatrix % hostData(j,0))/ &
                                    myPoly % qWeights % hostData(j)

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ScalarDGGradient_2D_cpu

  SUBROUTINE ScalarDGGradient_2D_gpu(myPoly,f_dev,bf_dev,gradF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(in)     :: bf_dev
    TYPE(c_ptr),INTENT(out)    :: gradF_dev

    CALL ScalarDGGradient_2D_gpu_wrapper(myPoly % dgMatrix % deviceData, &
                                         myPoly % bMatrix % deviceData, &
                                         myPoly % qWeights % deviceData, &
                                         f_dev,bf_dev,gradF_dev,myPoly % N, &
                                         nVariables,nElements)

  END SUBROUTINE ScalarDGGradient_2D_gpu

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
    REAL(prec) :: gf(1:2,1:2)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % N
          DO i = 0,myPoly % N

            gf(1,1) = 0.0_prec
            gf(2,1) = 0.0_prec
            gf(1,2) = 0.0_prec
            gf(2,2) = 0.0_prec
            DO ii = 0,myPoly % N
              gf(1,1) = gf(1,1) + myPoly % dMatrix % hostData(ii,i)*f(1,ii,j,iVar,iEl)
              gf(2,1) = gf(2,1) + myPoly % dMatrix % hostData(ii,i)*f(2,ii,j,iVar,iEl)
              gf(1,2) = gf(1,2) + myPoly % dMatrix % hostData(ii,j)*f(1,i,ii,iVar,iEl)
              gf(2,2) = gf(2,2) + myPoly % dMatrix % hostData(ii,j)*f(2,i,ii,iVar,iEl)
            END DO
            gradF(1:2,1:2,i,j,iVar,iEl) = gf(1:2,1:2)

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

    CALL VectorGradient_2D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                       f_dev,gradF_dev,myPoly % N, &
                                       nVariables,nElements)

  END SUBROUTINE VectorGradient_2D_gpu

  SUBROUTINE VectorDGGradient_2D_cpu(myPoly,f,bf,gradF,nVariables,nElements)
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
    REAL(prec),INTENT(in)  :: bf(1:2,0:myPoly % N,1:nVariables,1:4,1:nElements)
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
              gradF(1,1,i,j,iVar,iEl) = gradF(1,1,i,j,iVar,iEl) + myPoly % dgMatrix % hostData(ii,i)*f(1,ii,j,iVar,iEl)
              gradF(2,1,i,j,iVar,iEl) = gradF(2,1,i,j,iVar,iEl) + myPoly % dgMatrix % hostData(ii,i)*f(2,ii,j,iVar,iEl)
              gradF(1,2,i,j,iVar,iEl) = gradF(1,2,i,j,iVar,iEl) + myPoly % dgMatrix % hostData(ii,j)*f(1,i,ii,iVar,iEl)
              gradF(2,2,i,j,iVar,iEl) = gradF(2,2,i,j,iVar,iEl) + myPoly % dgMatrix % hostData(ii,j)*f(2,i,ii,iVar,iEl)
            END DO
            gradF(1,1,i,j,iVar,iEl) = gradF(1,1,i,j,iVar,iEl) + (myPoly % bMatrix % hostData(i,1)*bf(1,j,iVar,2,iEl) + &
                                                                 myPoly % bMatrix % hostData(i,0)*bf(1,j,iVar,4,iEl))/ &
                                      myPoly % qWeights % hostData(i)

            gradF(2,1,i,j,iVar,iEl) = gradF(2,1,i,j,iVar,iEl) + (myPoly % bMatrix % hostData(i,1)*bf(2,j,iVar,2,iEl) + &
                                                                 myPoly % bMatrix % hostData(i,0)*bf(2,j,iVar,4,iEl))/ &
                                      myPoly % qWeights % hostData(i)

            gradF(1,2,i,j,iVar,iEl) = gradF(1,2,i,j,iVar,iEl) + (myPoly % bMatrix % hostData(j,1)*bf(1,i,iVar,3,iEl) + &
                                                                 myPoly % bMatrix % hostData(j,0)*bf(1,i,iVar,1,iEl))/ &
                                      myPoly % qWeights % hostData(j)

            gradF(2,2,i,j,iVar,iEl) = gradF(2,2,i,j,iVar,iEl) + (myPoly % bMatrix % hostData(j,1)*bf(2,i,iVar,3,iEl) + &
                                                                 myPoly % bMatrix % hostData(j,0)*bf(2,i,iVar,1,iEl))/ &
                                      myPoly % qWeights % hostData(j)

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE VectorDGGradient_2D_cpu

  SUBROUTINE VectorDGGradient_2D_gpu(myPoly,f_dev,bf_dev,gradF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(in)     :: bf_dev
    TYPE(c_ptr),INTENT(out)    :: gradF_dev

    CALL VectorDGGradient_2D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                         myPoly % bMatrix % deviceData, &
                                         myPoly % qWeights % deviceData, &
                                         f_dev,bf_dev,gradF_dev,myPoly % N, &
                                         nVariables,nElements)

  END SUBROUTINE VectorDGGradient_2D_gpu

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

    CALL VectorDivergence_2D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                         f_dev,dF_dev,myPoly % N, &
                                         nVariables,nElements)

  END SUBROUTINE VectorDivergence_2D_gpu

  SUBROUTINE VectorDGDivergence_2D_cpu(myPoly,f,bF,dF,nVariables,nElements)
    ! Assumes bF is the vector component in the direction normal to the boundary
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(in)  :: bF(0:myPoly % N,1:nVariables,1:4,1:nElements)
    REAL(prec),INTENT(out) :: dF(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    REAL(prec) :: dfLoc
    INTEGER    :: i,j,ii,iVar,iEl

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO j = 0,myPoly % N
          DO i = 0,myPoly % N

            dfLoc = 0.0_prec
            DO ii = 0,myPoly % N
              dfLoc = dfLoc + myPoly % dgMatrix % hostData(ii,i)*f(1,ii,j,iVar,iEl) + &
                              myPoly % dgMatrix % hostData(ii,j)*f(2,i,ii,iVar,iEl)
            END DO

            dfLoc = dfLoc + (myPoly % bMatrix % hostData(i,1)*bF(j,iVar,2,iEl) + &
                             myPoly % bMatrix % hostData(i,0)*bF(j,iVar,4,iEl))/ &
                               myPoly % qWeights % hostData(i) + &
                            (myPoly % bMatrix % hostData(j,1)*bF(i,iVar,3,iEl) + &
                             myPoly % bMatrix % hostData(j,0)*bF(i,iVar,1,iEl))/ &
                               myPoly % qWeights % hostData(j)
            dF(i,j,iVar,iEl) = dFLoc

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE VectorDGDivergence_2D_cpu

  SUBROUTINE VectorDGDivergence_2D_gpu(myPoly,f_dev,bF_dev,dF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(in)     :: bF_dev
    TYPE(c_ptr),INTENT(out)    :: dF_dev

    CALL VectorDGDivergence_2D_gpu_wrapper(myPoly % dgMatrix % deviceData, &
                                           myPoly % bMatrix % deviceData, &
                                           myPoly % qWeights % deviceData, &
                                           f_dev,bF_dev,dF_dev,myPoly % N, &
                                           nVariables,nElements)

  END SUBROUTINE VectorDGDivergence_2D_gpu

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

    CALL TensorDivergence_2D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                         f_dev,dF_dev,myPoly % N, &
                                         nVariables,nElements)

  END SUBROUTINE TensorDivergence_2D_gpu

  SUBROUTINE TensorDGDivergence_2D_cpu(myPoly,f,bF,dF,nVariables,nElements)
    ! Note that the divergence is taken over the first dimension (row dimension) of the tensor matrix
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:2,1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(in)  :: bf(1:2,1:2,0:myPoly % N,1:nVariables,1:4,1:nElements)
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
              dF(1,i,j,iVar,iEl) = dF(1,i,j,iVar,iEl) + myPoly % dgMatrix % hostData(ii,i)*f(1,1,ii,j,iVar,iEl) + &
                                   myPoly % dgMatrix % hostData(ii,j)*f(2,1,i,ii,iVar,iEl)
              dF(2,i,j,iVar,iEl) = dF(2,i,j,iVar,iEl) + myPoly % dgMatrix % hostData(ii,i)*f(1,2,ii,j,iVar,iEl) + &
                                   myPoly % dgMatrix % hostData(ii,j)*f(2,2,i,ii,iVar,iEl)
            END DO

            dF(1,i,j,iVar,iEl) = dF(1,i,j,iVar,iEl) + (myPoly % bMatrix % hostData(i,1)*bf(1,1,j,iVar,2,iEl) + &
                                                       myPoly % bMatrix % hostData(i,0)*bf(1,1,j,iVar,4,iEl))/ &
                                 myPoly % qWeights % hostData(i) + &
                                 (myPoly % bMatrix % hostData(j,1)*bf(2,1,i,iVar,3,iEl) + &
                                  myPoly % bMatrix % hostData(j,0)*bf(2,1,i,iVar,1,iEl))/ &
                                 myPoly % qWeights % hostData(j)

            dF(2,i,j,iVar,iEl) = dF(2,i,j,iVar,iEl) + (myPoly % bMatrix % hostData(i,1)*bf(1,2,j,iVar,2,iEl) + &
                                                       myPoly % bMatrix % hostData(i,0)*bf(1,2,j,iVar,4,iEl))/ &
                                 myPoly % qWeights % hostData(i) + &
                                 (myPoly % bMatrix % hostData(j,1)*bf(2,2,i,iVar,3,iEl) + &
                                  myPoly % bMatrix % hostData(j,0)*bf(2,2,i,iVar,1,iEl))/ &
                                 myPoly % qWeights % hostData(j)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE TensorDGDivergence_2D_cpu

  SUBROUTINE TensorDGDivergence_2D_gpu(myPoly,f_dev,bF_dev,dF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(in)     :: bf_dev
    TYPE(c_ptr),INTENT(out)    :: dF_dev

    CALL TensorDGDivergence_2D_gpu_wrapper(myPoly % dgMatrix % deviceData, &
                                           myPoly % bMatrix % deviceData, &
                                           myPoly % qWeights % deviceData, &
                                           f_dev,bF_dev,dF_dev,myPoly % N, &
                                           nVariables,nElements)

  END SUBROUTINE TensorDGDivergence_2D_gpu

  SUBROUTINE ScalarGradient_3D_cpu(myPoly,f,gradF,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out) :: gradF(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER    :: i,j,k,ii,iVar,iEl
    REAL(prec) :: gf(1:3)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO k = 0,myPoly % N
          DO j = 0,myPoly % N
            DO i = 0,myPoly % N

              gF(1) = 0.0_prec
              gF(2) = 0.0_prec
              gF(3) = 0.0_prec
              DO ii = 0,myPoly % N
                gF(1) = gF(1) + myPoly % dMatrix % hostData(ii,i)*f(ii,j,k,iVar,iEl)
                gF(2) = gF(2) + myPoly % dMatrix % hostData(ii,j)*f(i,ii,k,iVar,iEl)
                gF(3) = gF(3) + myPoly % dMatrix % hostData(ii,k)*f(i,j,ii,iVar,iEl)
              END DO

              gradF(1,i,j,k,iVar,iEl) = gF(1)
              gradF(2,i,j,k,iVar,iEl) = gF(2)
              gradF(3,i,j,k,iVar,iEl) = gF(3)

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
    REAL(prec) :: gF(1:3,1:3)

    DO iEl = 1,nElements
      DO iVar = 1,nVariables
        DO k = 0,myPoly % N
          DO j = 0,myPoly % N
            DO i = 0,myPoly % N

              gF = 0.0_prec
              DO ii = 0,myPoly % N
                gF(1,1) = gF(1,1) + myPoly % dMatrix % hostData(ii,i)*f(1,ii,j,k,iVar,iEl)
                gF(2,1) = gF(2,1) + myPoly % dMatrix % hostData(ii,i)*f(2,ii,j,k,iVar,iEl)
                gF(3,1) = gF(3,1) + myPoly % dMatrix % hostData(ii,i)*f(3,ii,j,k,iVar,iEl)
                gF(1,2) = gF(1,2) + myPoly % dMatrix % hostData(ii,j)*f(1,i,ii,k,iVar,iEl)
                gF(2,2) = gF(2,2) + myPoly % dMatrix % hostData(ii,j)*f(2,i,ii,k,iVar,iEl)
                gF(3,2) = gF(3,2) + myPoly % dMatrix % hostData(ii,j)*f(3,i,ii,k,iVar,iEl)
                gF(1,3) = gF(1,3) + myPoly % dMatrix % hostData(ii,k)*f(1,i,j,ii,iVar,iEl)
                gF(2,3) = gF(2,3) + myPoly % dMatrix % hostData(ii,k)*f(2,i,j,ii,iVar,iEl)
                gF(3,3) = gF(3,3) + myPoly % dMatrix % hostData(ii,k)*f(3,i,j,ii,iVar,iEl)
              END DO

              gradF(1,1,i,j,k,iVar,iEl) = gF(1,1) 
              gradF(2,1,i,j,k,iVar,iEl) = gF(2,1) 
              gradF(3,1,i,j,k,iVar,iEl) = gF(3,1) 
              gradF(1,2,i,j,k,iVar,iEl) = gF(1,2) 
              gradF(2,2,i,j,k,iVar,iEl) = gF(2,2) 
              gradF(3,2,i,j,k,iVar,iEl) = gF(3,2) 
              gradF(1,3,i,j,k,iVar,iEl) = gF(1,3) 
              gradF(2,3,i,j,k,iVar,iEl) = gF(2,3) 
              gradF(3,3,i,j,k,iVar,iEl) = gF(3,3) 

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

    CALL VectorDivergence_3D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                         f_dev,dF_dev,myPoly % N, &
                                         nVariables,nElements)

  END SUBROUTINE VectorDivergence_3D_gpu

  SUBROUTINE VectorDGDivergence_3D_cpu(myPoly,f,bF,dF,nVariables,nElements)
    ! Assumes bF is the vector component in the direction normal to the element boundaries
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(in)  :: bf(0:myPoly % N,0:myPoly % N,1:nVariables,1:6,1:nElements)
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
                dF(i,j,k,iVar,iEl) = dF(i,j,k,iVar,iEl) + myPoly % dgMatrix % hostData(ii,i)*f(1,ii,j,k,iVar,iEl) + &
                                     myPoly % dgMatrix % hostData(ii,j)*f(2,i,ii,k,iVar,iEl) + &
                                     myPoly % dgMatrix % hostData(ii,k)*f(3,i,j,ii,iVar,iEl)
              END DO

              dF(i,j,k,iVar,iEl) = dF(i,j,k,iVar,iEl) + (myPoly % bMatrix % hostData(i,1)*bF(j,k,iVar,3,iEl) + & ! east
                                                         myPoly % bMatrix % hostData(i,0)*bF(j,k,iVar,5,iEl))/ &  ! west
                                   myPoly % qWeights % hostData(i) + &
                                   (myPoly % bMatrix % hostData(j,1)*bF(i,k,iVar,4,iEl) + & ! north
                                    myPoly % bMatrix % hostData(j,0)*bF(i,k,iVar,2,iEl))/ &  ! south
                                   myPoly % qWeights % hostData(j) + &
                                   (myPoly % bMatrix % hostData(k,1)*bF(i,j,iVar,6,iEl) + & ! top
                                    myPoly % bMatrix % hostData(k,0)*bF(i,j,iVar,1,iEl))/ &  ! bottom
                                   myPoly % qWeights % hostData(k)

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE VectorDGDivergence_3D_cpu

  SUBROUTINE VectorDGDivergence_3D_gpu(myPoly,f_dev,bF_dev,dF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(in)     :: bF_dev
    TYPE(c_ptr),INTENT(out)    :: dF_dev

    CALL VectorDGDivergence_3D_gpu_wrapper(myPoly % dgMatrix % deviceData, &
                                           myPoly % bMatrix % deviceData, &
                                           myPoly % qWeights % deviceData, &
                                           f_dev,bF_dev,dF_dev,myPoly % N, &
                                           nVariables,nElements)

  END SUBROUTINE VectorDGDivergence_3D_gpu

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

    CALL TensorDivergence_3D_gpu_wrapper(myPoly % dMatrix % deviceData, &
                                         f_dev,dF_dev,myPoly % N, &
                                         nVariables,nElements)

  END SUBROUTINE TensorDivergence_3D_gpu

  SUBROUTINE TensorDGDivergence_3D_cpu(myPoly,f,bF,dF,nVariables,nElements)
    ! Note that the divergence is taken over the first dimension (row dimension) of the tensor matrix
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)     :: nVariables,nElements
    REAL(prec),INTENT(in)  :: f(1:3,1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(in)  :: bF(1:3,1:3,0:myPoly % N,0:myPoly % N,1:nVariables,1:6,1:nElements)
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
                                       myPoly % dgMatrix % hostData(ii,i)*f(1,1,ii,j,k,iVar,iEl) + &
                                       myPoly % dgMatrix % hostData(ii,j)*f(2,1,i,ii,k,iVar,iEl) + &
                                       myPoly % dgMatrix % hostData(ii,k)*f(3,1,i,j,ii,iVar,iEl)

                dF(2,i,j,k,iVar,iEl) = dF(2,i,j,k,iVar,iEl) + &
                                       myPoly % dgMatrix % hostData(ii,i)*f(1,2,ii,j,k,iVar,iEl) + &
                                       myPoly % dgMatrix % hostData(ii,j)*f(2,2,i,ii,k,iVar,iEl) + &
                                       myPoly % dgMatrix % hostData(ii,k)*f(3,2,i,j,ii,iVar,iEl)

                dF(3,i,j,k,iVar,iEl) = dF(3,i,j,k,iVar,iEl) + &
                                       myPoly % dgMatrix % hostData(ii,i)*f(1,3,ii,j,k,iVar,iEl) + &
                                       myPoly % dgMatrix % hostData(ii,j)*f(2,3,i,ii,k,iVar,iEl) + &
                                       myPoly % dgMatrix % hostData(ii,k)*f(3,3,i,j,ii,iVar,iEl)
              END DO

              dF(1,i,j,k,iVar,iEl) = dF(1,i,j,k,iVar,iEl) + (myPoly % bMatrix % hostData(i,1)*bF(1,1,j,k,iVar,3,iEl) + & ! east
                                                     myPoly % bMatrix % hostData(i,0)*bF(1,1,j,k,iVar,5,iEl))/ &  ! west
                                     myPoly % qWeights % hostData(i) + &
                                     (myPoly % bMatrix % hostData(j,1)*bF(2,1,i,k,iVar,4,iEl) + & ! north
                                      myPoly % bMatrix % hostData(j,0)*bF(2,1,i,k,iVar,2,iEl))/ &  ! south
                                     myPoly % qWeights % hostData(j) + &
                                     (myPoly % bMatrix % hostData(k,1)*bF(3,1,i,j,iVar,6,iEl) + & ! top
                                      myPoly % bMatrix % hostData(k,0)*bF(3,1,i,j,iVar,1,iEl))/ &  ! bottom
                                     myPoly % qWeights % hostData(k)

              dF(2,i,j,k,iVar,iEl) = dF(2,i,j,k,iVar,iEl) + (myPoly % bMatrix % hostData(i,1)*bF(1,2,j,k,iVar,3,iEl) + & ! east
                                                     myPoly % bMatrix % hostData(i,0)*bF(1,2,j,k,iVar,5,iEl))/ &  ! west
                                     myPoly % qWeights % hostData(i) + &
                                     (myPoly % bMatrix % hostData(j,1)*bF(2,2,i,k,iVar,4,iEl) + & ! north
                                      myPoly % bMatrix % hostData(j,0)*bF(2,2,i,k,iVar,2,iEl))/ &  ! south
                                     myPoly % qWeights % hostData(j) + &
                                     (myPoly % bMatrix % hostData(k,1)*bF(3,2,i,j,iVar,6,iEl) + & ! top
                                      myPoly % bMatrix % hostData(k,0)*bF(3,2,i,j,iVar,1,iEl))/ &  ! bottom
                                     myPoly % qWeights % hostData(k)

              dF(3,i,j,k,iVar,iEl) = dF(3,i,j,k,iVar,iEl) + (myPoly % bMatrix % hostData(i,1)*bF(1,3,j,k,iVar,3,iEl) + & ! east
                                                     myPoly % bMatrix % hostData(i,0)*bF(1,3,j,k,iVar,5,iEl))/ &  ! west
                                     myPoly % qWeights % hostData(i) + &
                                     (myPoly % bMatrix % hostData(j,1)*bF(2,3,i,k,iVar,4,iEl) + & ! north
                                      myPoly % bMatrix % hostData(j,0)*bF(2,3,i,k,iVar,2,iEl))/ &  ! south
                                     myPoly % qWeights % hostData(j) + &
                                     (myPoly % bMatrix % hostData(k,1)*bF(3,3,i,j,iVar,6,iEl) + & ! top
                                      myPoly % bMatrix % hostData(k,0)*bF(3,3,i,j,iVar,1,iEl))/ &  ! bottom
                                     myPoly % qWeights % hostData(k)

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE TensorDGDivergence_3D_cpu

  SUBROUTINE TensorDGDivergence_3D_gpu(myPoly,f_dev,bF_dev,dF_dev,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    TYPE(c_ptr),INTENT(in)     :: f_dev
    TYPE(c_ptr),INTENT(in)     :: bF_dev
    TYPE(c_ptr),INTENT(out)    :: dF_dev

    CALL TensorDGDivergence_3D_gpu_wrapper(myPoly % dgMatrix % deviceData, &
                                           myPoly % bMatrix % deviceData, &
                                           myPoly % qWeights % deviceData, &
                                           f_dev,bF_dev,dF_dev,myPoly % N, &
                                           nVariables,nElements)

  END SUBROUTINE TensorDGDivergence_3D_gpu
  ! /////////////////////////////// !
  ! Boundary Interpolation Routines !

  SUBROUTINE ScalarBoundaryInterp_1D_cpu(myPoly,f,fBound,nVariables,nElements)
    IMPLICIT NONE
    CLASS(Lagrange),INTENT(in) :: myPoly
    INTEGER,INTENT(in)         :: nVariables,nElements
    REAL(prec),INTENT(in)      :: f(0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec),INTENT(out)     :: fBound(1:nVariables,1:2,1:nElements)
    ! Local
    INTEGER :: ii,iVar,iEl
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
    INTEGER :: i,ii,iVar,iEl
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
    INTEGER :: i,ii,idir,iVar,iEl
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
    INTEGER :: i,ii,idir,jdir,iVar,iEl
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
              fb(1) = fb(1) + myPoly % bMatrix % hostData(ii,0)*f(i,j,ii,iVar,iEl) ! Bottom
              fb(2) = fb(2) + myPoly % bMatrix % hostData(ii,0)*f(i,ii,j,iVar,iEl) ! South
              fb(3) = fb(3) + myPoly % bMatrix % hostData(ii,1)*f(ii,i,j,iVar,iEl) ! East
              fb(4) = fb(4) + myPoly % bMatrix % hostData(ii,1)*f(i,ii,j,iVar,iEl) ! North
              fb(5) = fb(5) + myPoly % bMatrix % hostData(ii,0)*f(ii,i,j,iVar,iEl) ! West
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
                fb(idir,1) = fb(idir,1) + myPoly % bMatrix % hostData(ii,0)*f(idir,i,j,ii,iVar,iEl) ! Bottom
                fb(idir,2) = fb(idir,2) + myPoly % bMatrix % hostData(ii,0)*f(idir,i,ii,j,iVar,iEl) ! South
                fb(idir,3) = fb(idir,3) + myPoly % bMatrix % hostData(ii,1)*f(idir,ii,i,j,iVar,iEl) ! East
                fb(idir,4) = fb(idir,4) + myPoly % bMatrix % hostData(ii,1)*f(idir,i,ii,j,iVar,iEl) ! North
                fb(idir,5) = fb(idir,5) + myPoly % bMatrix % hostData(ii,0)*f(idir,ii,i,j,iVar,iEl) ! West
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
                  fb(idir,jdir,1) = fb(idir,jdir,1) + myPoly % bMatrix % hostData(ii,0)*f(idir,jdir,i,j,ii,iVar,iEl) ! Bottom
                  fb(idir,jdir,2) = fb(idir,jdir,2) + myPoly % bMatrix % hostData(ii,0)*f(idir,jdir,i,ii,j,iVar,iEl) ! South
                  fb(idir,jdir,3) = fb(idir,jdir,3) + myPoly % bMatrix % hostData(ii,1)*f(idir,jdir,ii,i,j,iVar,iEl) ! East
                  fb(idir,jdir,4) = fb(idir,jdir,4) + myPoly % bMatrix % hostData(ii,1)*f(idir,jdir,i,ii,j,iVar,iEl) ! North
                  fb(idir,jdir,5) = fb(idir,jdir,5) + myPoly % bMatrix % hostData(ii,0)*f(idir,jdir,ii,i,j,iVar,iEl) ! West
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
    REAL(real64) :: bWeights(0:myPoly % N)
    REAL(real64) :: controlPoints(0:myPoly % N)

    DO i = 0,myPoly % N
      bWeights(i) = 1.0_real64
      controlPoints(i) = REAL(myPoly % controlPoints % hostData(i),real64)
    END DO

    ! Computes the product w_k = w_k*(s_k - s_j), k /= j
    DO j = 1,myPoly % N
      DO i = 0,j - 1

        bWeights(i) = bWeights(i)*(controlPoints(i) - controlPoints(j))
        bWeights(j) = bWeights(j)*(controlPoints(j) - controlPoints(i))

      END DO
    END DO

    DO j = 0,myPoly % N
      bWeights(j) = 1.0_prec/bWeights(j)
      myPoly % bWeights % hostData(j) = REAL(bWeights(j),prec)
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
    INTEGER    :: row,col
    LOGICAL    :: rowHasMatch
    REAL(real64) :: temp1,temp2
    REAL(real64) :: iMatrix(0:myPoly % M,0:myPoly % N)
    REAL(real64) :: bWeights(0:myPoly % N)
    REAL(real64) :: controlPoints(0:myPoly % N)
    REAL(real64) :: targetPoints(0:myPoly % M)

    DO col = 0,myPoly % N
      controlPoints(col) = REAL(myPoly % controlPoints % hostData(col),real64)
      bWeights(col) = REAL(myPoly % bWeights % hostData(col),real64)
    END DO
    DO row = 0,myPoly % M
      targetPoints(row) = REAL(myPoly % targetPoints % hostData(row),real64)
    END DO

    DO row = 0,myPoly % M

      rowHasMatch = .FALSE.

      DO col = 0,myPoly % N

        iMatrix(row,col) = 0.0_real64

        IF (AlmostEqual(targetPoints(row),controlPoints(col))) THEN
          rowHasMatch = .TRUE.
          iMatrix(row,col) = 1.0_real64
        END IF

      END DO

      IF (.NOT. (rowHasMatch)) THEN

        temp1 = 0.0_real64

        DO col = 0,myPoly % N
          temp2 = bWeights(col)/ &
                  (targetPoints(row) - &
                   controlPoints(col))
          iMatrix(row,col) = temp2
          temp1 = temp1 + temp2
        END DO

        DO col = 0,myPoly % N
          iMatrix(row,col) = iMatrix(row,col)/temp1
        END DO

      END IF

    END DO

    DO row = 0,myPoly % M
      DO col = 0,myPoly % N
        myPoly % iMatrix % hostData(col,row) = REAL(iMatrix(row,col),prec)
      END DO
    END DO

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
    INTEGER      :: row,col
    REAL(real64) :: dmat(0:myPoly % N,0:myPoly % N)
    REAL(real64) :: dgmat(0:myPoly % N,0:myPoly % N)
    REAL(real64) :: bWeights(0:myPoly % N)
    REAL(real64) :: qWeights(0:myPoly % N)
    REAL(real64) :: controlPoints(0:myPoly % N)

    DO row = 0,myPoly % N
      bWeights(row) = REAL(myPoly % bWeights % hostData(row),real64)
      qWeights(row) = REAL(myPoly % qWeights % hostData(row),real64)
      controlPoints(row) = REAL(myPoly % controlPoints % hostData(row),real64)
    END DO

    DO row = 0,myPoly % N

      dmat(row,row) = 0.0_prec

      DO col = 0,myPoly % N

        IF (.NOT. (col == row)) THEN

          dmat(row,col) = bWeights(col)/ &
                          (bWeights(row)* &
                           (controlPoints(row) - &
                            controlPoints(col)))

          dmat(row,row) = dmat(row,row) - dmat(row,col)

        END IF

      END DO

    END DO

    DO row = 0,myPoly % N
      DO col = 0,myPoly % N
        dgmat(row,col) = -dmat(col,row)* &
                         qWeights(col)/ &
                         qWeights(row)
      END DO
    END DO

    DO row = 0,myPoly % N
      DO col = 0,myPoly % N
        myPoly % dMatrix % hostData(row,col) = REAL(dmat(col,row),prec)
        myPoly % dgMatrix % hostData(row,col) = REAL(dgmat(col,row),prec)
      END DO
    END DO

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
    INTEGER    :: j
    LOGICAL    :: xMatchesNode
    REAL(real64) :: temp1,temp2
    REAL(real64) :: sELocal
    REAL(real64) :: controlPoints(0:myPoly % N)
    REAL(real64) :: bWeights(0:myPoly % N)
    REAL(real64) :: lS(0:myPoly % N)

    sELocal = REAL(sE,real64)
    DO j = 0,myPoly % N
      controlPoints(j) = REAL(myPoly % controlPoints % hostData(j),real64)
      bWeights(j) = REAL(myPoly % bWeights % hostData(j),real64)
    END DO

    xMatchesNode = .FALSE.

    DO j = 0,myPoly % N

      lS(j) = 0.0_real64
      IF (AlmostEqual(sELocal,controlPoints(j))) THEN
        lS(j) = 1.0_real64
        xMatchesNode = .TRUE.
      END IF

    END DO

    IF (xMatchesNode) THEN
      DO j = 0,myPoly % N
        lAtS(j) = REAL(lS(j),prec)
      END DO
      RETURN
    END IF

    temp1 = 0.0_real64

    DO j = 0,myPoly % N
      temp2 = bWeights(j)/(sE - controlPoints(j))
      lS(j) = temp2
      temp1 = temp1 + temp2
    END DO

    lS = lS/temp1

    DO j = 0,myPoly % N
      lAtS(j) = REAL(lS(j),prec)
    END DO

  END FUNCTION CalculateLagrangePolynomials

END MODULE SELF_Lagrange
