! Lagrange_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE Lagrange_Class

!src/self
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
USE Quadrature

USE hip
USE iso_c_binding

IMPLICIT NONE


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

  TYPE, PUBLIC :: Lagrange
  !> A data structure for working with Lagrange Interpolating Polynomials in one, two, and three dimensions.
  !>
  !> attribute : controlPoints(0:N) : The set of nodes in one dimension where data is known. To create higher dimension interpolation and differentiation operators, structured grids in two and three dimensions are created by tensor products of the controlPoints. This design decision implies that all Spectral Element Methods supported by the Lagrange class have the same polynomial degree in each computational/spatial dimension. In practice, the controlPoints are the Legendre-Gauss, Legendre-Gauss-Lobatto, Legendre-Gauss-Radau, Chebyshev-Gauss, Chebyshev-Gauss-Lobatto, or Chebyshev-Gauss-Radau quadrature points over the domain [-1,1] (computational space). The Build routine for this class restricts controlPoints to one of these quadrature types or uniform points on [-1,1].
  !
  !> attribute : targetPoints(0:M) : The set of nodes in one dimension where data is to be interpolated to. To create higher dimension interpolation and differentiation operators, structured grids in two and three dimensions are created by tensor products of the targetPoints. In practice, the targetPoints are set to a uniformly distributed set of points between [-1,1] (computational space) to allow for interpolation from unevenly spaced quadrature points to a plotting grid.
  !
  !> attribute : bWeights(0:N) : The barycentric weights that are calculated from the controlPoints and used for interpolation.
  !
  !> attribute : qWeights(0:N) : The quadrature weights for discrete integration. The quadradture weights depend on the type of controlPoints provided; one of Legendre-Gauss, Legendre-Gauss-Lobatto, Legendre-Gauss-Radau, Chebyshev-Gauss, Chebyshev-Gauss-Lobatto, Chebyshev-Gauss Radau, or Uniform. If Uniform, the quadrature weights are constant = dx = 2.0/(N+1).
  ! 
  !> attribute : iMatrix(0:N,0:M) :
  !> attribute : dMatrix(0:N,0:N) : TO DO
  !> attribute : bMatrix(0:N,1:2) : TO DO

    INTEGER :: N     
    INTEGER :: M 
    REAL(prec), POINTER :: controlPoints(:)
    REAL(prec), POINTER :: targetPoints(:)
    REAL(prec), POINTER :: bWeights(:)
    REAL(prec), POINTER :: qWeights(:)
    REAL(prec), POINTER :: iMatrix(:,:)
    REAL(prec), POINTER :: dMatrix(:,:)  
    REAL(prec), POINTER :: bMatrix(:,:)

    TYPE(c_ptr) :: bWeights_dev
    TYPE(c_ptr) :: qWeights_dev
    TYPE(c_ptr) :: iMatrix_dev
    TYPE(c_ptr) :: dMatrix_dev
    TYPE(c_ptr) :: bMatrix_dev

    CONTAINS
      
      PROCEDURE, PUBLIC :: Build => Build_Lagrange
      PROCEDURE, PUBLIC :: Trash => Trash_Lagrange

      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_Lagrange
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_Lagrange
      
      GENERIC, PUBLIC :: ScalarGridInterp_1D => ScalarGridInterp_1D_cpu, ScalarGridInterp_1D_gpu
      PROCEDURE, PRIVATE :: ScalarGridInterp_1D_cpu, ScalarGridInterp_1D_gpu

      GENERIC, PUBLIC :: ScalarGridInterp_2D => ScalarGridInterp_2D_cpu, ScalarGridInterp_2D_gpu
      PROCEDURE, PRIVATE :: ScalarGridInterp_2D_cpu, ScalarGridInterp_2D_gpu

      GENERIC, PUBLIC :: VectorGridInterp_2D => VectorGridInterp_2D_cpu, VectorGridInterp_2D_gpu
      PROCEDURE, PRIVATE :: VectorGridInterp_2D_cpu, VectorGridInterp_2D_gpu

      GENERIC, PUBLIC :: TensorGridInterp_2D => TensorGridInterp_2D_cpu, TensorGridInterp_2D_gpu
      PROCEDURE, PRIVATE :: TensorGridInterp_2D_cpu, TensorGridInterp_2D_gpu

      GENERIC, PUBLIC :: ScalarGridInterp_3D => ScalarGridInterp_3D_cpu, ScalarGridInterp_3D_gpu
      PROCEDURE, PRIVATE :: ScalarGridInterp_3D_cpu, ScalarGridInterp_3D_gpu

      GENERIC, PUBLIC :: VectorGridInterp_3D => VectorGridInterp_3D_cpu, VectorGridInterp_3D_gpu
      PROCEDURE, PRIVATE :: VectorGridInterp_3D_cpu, VectorGridInterp_3D_gpu

      GENERIC, PUBLIC :: TensorGridInterp_3D => TensorGridInterp_3D_cpu, TensorGridInterp_3D_gpu
      PROCEDURE, PRIVATE :: TensorGridInterp_3D_cpu, TensorGridInterp_3D_gpu

      GENERIC, PUBLIC :: ScalarBoundaryInterp_1D => ScalarBoundaryInterp_1D_cpu, ScalarBoundaryInterp_1D_gpu
      PROCEDURE, PRIVATE :: ScalarBoundaryInterp_1D_cpu, ScalarBoundaryInterp_1D_gpu

      GENERIC, PUBLIC :: ScalarBoundaryInterp_2D => ScalarBoundaryInterp_2D_cpu, ScalarBoundaryInterp_2D_gpu
      PROCEDURE, PRIVATE :: ScalarBoundaryInterp_2D_cpu, ScalarBoundaryInterp_2D_gpu

      GENERIC, PUBLIC :: VectorBoundaryInterp_2D => VectorBoundaryInterp_2D_cpu, VectorBoundaryInterp_2D_gpu
      PROCEDURE, PRIVATE :: VectorBoundaryInterp_2D_cpu, VectorBoundaryInterp_2D_gpu

      GENERIC, PUBLIC :: TensorBoundaryInterp_2D => TensorBoundaryInterp_2D_cpu, TensorBoundaryInterp_2D_gpu
      PROCEDURE, PRIVATE :: TensorBoundaryInterp_2D_cpu, TensorBoundaryInterp_2D_gpu

      GENERIC, PUBLIC :: ScalarBoundaryInterp_3D => ScalarBoundaryInterp_3D_cpu, ScalarBoundaryInterp_3D_gpu
      PROCEDURE, PRIVATE :: ScalarBoundaryInterp_3D_cpu, ScalarBoundaryInterp_3D_gpu

      GENERIC, PUBLIC :: VectorBoundaryInterp_3D => VectorBoundaryInterp_3D_cpu, VectorBoundaryInterp_3D_gpu
      PROCEDURE, PRIVATE :: VectorBoundaryInterp_3D_cpu, VectorBoundaryInterp_3D_gpu

      GENERIC, PUBLIC :: TensorBoundaryInterp_3D => TensorBoundaryInterp_3D_cpu, TensorBoundaryInterp_3D_gpu
      PROCEDURE, PRIVATE :: TensorBoundaryInterp_3D_cpu, TensorBoundaryInterp_3D_gpu

      GENERIC, PUBLIC :: Derivative_1D => Derivative_1D_cpu, Derivative_1D_gpu
      PROCEDURE, PRIVATE :: Derivative_1D_cpu, Derivative_1D_gpu

      GENERIC, PUBLIC :: ScalarGradient_2D => ScalarGradient_2D_cpu, ScalarGradient_2D_gpu
      PROCEDURE, PRIVATE :: ScalarGradient_2D_cpu, ScalarGradient_2D_gpu

      GENERIC, PUBLIC :: VectorGradient_2D => VectorGradient_2D_cpu, VectorGradient_2D_gpu
      PROCEDURE, PRIVATE :: VectorGradient_2D_cpu, VectorGradient_2D_gpu

      GENERIC, PUBLIC :: VectorDivergence_2D => VectorDivergence_2D_cpu, VectorDivergence_2D_gpu
      PROCEDURE, PRIVATE :: VectorDivergence_2D_cpu, VectorDivergence_2D_gpu

      GENERIC, PUBLIC :: VectorCurl_2D => VectorCurl_2D_cpu, VectorCurl_2D_gpu
      PROCEDURE, PRIVATE :: VectorCurl_2D_cpu, VectorCurl_2D_gpu

      GENERIC, PUBLIC :: ScalarGradient_3D => ScalarGradient_3D_cpu, ScalarGradient_3D_gpu
      PROCEDURE, PRIVATE :: ScalarGradient_3D_cpu, ScalarGradient_3D_gpu

      GENERIC, PUBLIC :: VectorGradient_3D => VectorGradient_3D_cpu, VectorGradient_3D_gpu
      PROCEDURE, PRIVATE :: VectorGradient_3D_cpu, VectorGradient_3D_gpu

      GENERIC, PUBLIC :: VectorDivergence_3D => VectorDivergence_3D_cpu, VectorDivergence_3D_gpu
      PROCEDURE, PRIVATE :: VectorDivergence_3D_cpu, VectorDivergence_3D_gpu

      GENERIC, PUBLIC :: VectorCurl_3D => VectorCurl_3D_cpu, VectorCurl_3D_gpu
      PROCEDURE, PRIVATE :: VectorCurl_3D_cpu, VectorCurl_3D_gpu
   
      PROCEDURE, PRIVATE :: CalculateBarycentricWeights 
      PROCEDURE, PRIVATE :: CalculateInterpolationMatrix
      PROCEDURE, PRIVATE :: CalculateDerivativeMatrix
      PROCEDURE, PRIVATE :: CalculateLagrangePolynomials
      
    END TYPE Lagrange

    INTERFACE
      SUBROUTINE ScalarGridInterp_1D_gpu_wrapper(iMatrixT_dev, f_dev, fInterp_dev, N, M, nVar, nEl) bind(c,name="ScalarGridInterp_1D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: iMatrixT_dev, f_dev, fInterp_dev
        INTEGER, VALUE :: N, M, nVar, nEl
      END SUBROUTINE ScalarGridInterp_1D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE ScalarGridInterp_2D_gpu_wrapper(iMatrixT_dev, f_dev, fInterp_dev, N, M, nVar, nEl) bind(c,name="ScalarGridInterp_2D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: iMatrixT_dev, f_dev, fInterp_dev
        INTEGER, VALUE :: N, M, nVar, nEl
      END SUBROUTINE ScalarGridInterp_2D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE VectorGridInterp_2D_gpu_wrapper(iMatrixT_dev, f_dev, fInterp_dev, N, M, nVar, nEl) bind(c,name="VectorGridInterp_2D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: iMatrixT_dev, f_dev, fInterp_dev
        INTEGER, VALUE :: N, M, nVar, nEl
      END SUBROUTINE VectorGridInterp_2D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE TensorGridInterp_2D_gpu_wrapper(iMatrixT_dev, f_dev, fInterp_dev, N, M, nVar, nEl) bind(c,name="TensorGridInterp_2D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: iMatrixT_dev, f_dev, fInterp_dev
        INTEGER, VALUE :: N, M, nVar, nEl
      END SUBROUTINE TensorGridInterp_2D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE ScalarGridInterp_3D_gpu_wrapper(iMatrixT_dev, f_dev, fInterp_dev, N, M, nVar, nEl) bind(c,name="ScalarGridInterp_3D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: iMatrixT_dev, f_dev, fInterp_dev
        INTEGER, VALUE :: N, M, nVar, nEl
      END SUBROUTINE ScalarGridInterp_3D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE VectorGridInterp_3D_gpu_wrapper(iMatrixT_dev, f_dev, fInterp_dev, N, M, nVar, nEl) bind(c,name="VectorGridInterp_3D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: iMatrixT_dev, f_dev, fInterp_dev
        INTEGER, VALUE :: N, M, nVar, nEl
      END SUBROUTINE VectorGridInterp_3D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE TensorGridInterp_3D_gpu_wrapper(iMatrixT_dev, f_dev, fInterp_dev, N, M, nVar, nEl) bind(c,name="TensorGridInterp_3D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: iMatrixT_dev, f_dev, fInterp_dev
        INTEGER, VALUE :: N, M, nVar, nEl
      END SUBROUTINE TensorGridInterp_3D_gpu_wrapper
    END INTERFACE

    ! /////////////// !
    ! Boundary Interpolation Routines

    INTERFACE
      SUBROUTINE ScalarBoundaryInterp_1D_gpu_wrapper(bMatrix_dev, f_dev, fBound_dev, N, nVar, nEl) bind(c,name="ScalarBoundaryInterp_1D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: bMatrix_dev, f_dev, fBound_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE ScalarBoundaryInterp_1D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE ScalarBoundaryInterp_2D_gpu_wrapper(bMatrix_dev, f_dev, fBound_dev, N, nVar, nEl) bind(c,name="ScalarBoundaryInterp_2D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: bMatrix_dev, f_dev, fBound_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE ScalarBoundaryInterp_2D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE VectorBoundaryInterp_2D_gpu_wrapper(bMatrix_dev, f_dev, fBound_dev, N, nVar, nEl) bind(c,name="VectorBoundaryInterp_2D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: bMatrix_dev, f_dev, fBound_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE VectorBoundaryInterp_2D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE TensorBoundaryInterp_2D_gpu_wrapper(bMatrix_dev, f_dev, fBound_dev, N, nVar, nEl) bind(c,name="TensorBoundaryInterp_2D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: bMatrix_dev, f_dev, fBound_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE TensorBoundaryInterp_2D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE ScalarBoundaryInterp_3D_gpu_wrapper(bMatrix_dev, f_dev, fBound_dev, N, nVar, nEl) bind(c,name="ScalarBoundaryInterp_3D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: bMatrix_dev, f_dev, fBound_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE ScalarBoundaryInterp_3D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE VectorBoundaryInterp_3D_gpu_wrapper(bMatrix_dev, f_dev, fBound_dev, N, nVar, nEl) bind(c,name="VectorBoundaryInterp_3D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: bMatrix_dev, f_dev, fBound_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE VectorBoundaryInterp_3D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE TensorBoundaryInterp_3D_gpu_wrapper(bMatrix_dev, f_dev, fBound_dev, N, nVar, nEl) bind(c,name="TensorBoundaryInterp_3D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: bMatrix_dev, f_dev, fBound_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE TensorBoundaryInterp_3D_gpu_wrapper
    END INTERFACE

    ! /////////////// !

    INTERFACE
      SUBROUTINE Derivative_1D_gpu_wrapper(dMatrixT_dev, f_dev, df_dev, N, nVar, nEl) bind(c,name="Derivative_1D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: dMatrixT_dev, f_dev, df_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE Derivative_1D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE ScalarGradient_2D_gpu_wrapper(dMatrixT_dev, f_dev, df_dev, N, nVar, nEl) bind(c,name="ScalarGradient_2D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: dMatrixT_dev, f_dev, df_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE ScalarGradient_2D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE VectorGradient_2D_gpu_wrapper(dMatrixT_dev, f_dev, df_dev, N, nVar, nEl) bind(c,name="VectorGradient_2D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: dMatrixT_dev, f_dev, df_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE VectorGradient_2D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE VectorDivergence_2D_gpu_wrapper(dMatrixT_dev, f_dev, df_dev, N, nVar, nEl) bind(c,name="VectorDivergence_2D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: dMatrixT_dev, f_dev, df_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE VectorDivergence_2D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE VectorCurl_2D_gpu_wrapper(dMatrixT_dev, f_dev, df_dev, N, nVar, nEl) bind(c,name="VectorCurl_2D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: dMatrixT_dev, f_dev, df_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE VectorCurl_2D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE ScalarGradient_3D_gpu_wrapper(dMatrixT_dev, f_dev, df_dev, N, nVar, nEl) bind(c,name="ScalarGradient_3D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: dMatrixT_dev, f_dev, df_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE ScalarGradient_3D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE VectorGradient_3D_gpu_wrapper(dMatrixT_dev, f_dev, df_dev, N, nVar, nEl) bind(c,name="VectorGradient_3D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: dMatrixT_dev, f_dev, df_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE VectorGradient_3D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE VectorDivergence_3D_gpu_wrapper(dMatrixT_dev, f_dev, df_dev, N, nVar, nEl) bind(c,name="VectorDivergence_3D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: dMatrixT_dev, f_dev, df_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE VectorDivergence_3D_gpu_wrapper
    END INTERFACE

    INTERFACE
      SUBROUTINE VectorCurl_3D_gpu_wrapper(dMatrixT_dev, f_dev, df_dev, N, nVar, nEl) bind(c,name="VectorCurl_3D_gpu_wrapper")
        USE iso_c_binding
        IMPLICIT NONE
        TYPE(c_ptr) :: dMatrixT_dev, f_dev, df_dev
        INTEGER, VALUE :: N, nVar, nEl
      END SUBROUTINE VectorCurl_3D_gpu_wrapper
    END INTERFACE

 CONTAINS

! ================================================================================================ !
!
! Build_Lagrange
!
!   A manual constructor for the Lagrange class that allocates memory and fills in data 
!   for the attributes of the Lagrange class.
!  
!   The Build subroutine allocates memory for the interpolation and target points, barycentric
!   weights, interpolation matrix, and derivative matrix.
!
!   Usage :
!
!     TYPE(Lagrange) :: interp 
!     INTEGER        :: N, M
!     REAL(prec)     :: interpNodes(0:N), targetNodes(0:M+1) 
!
!     CALL interp % Build( N, M, interpNodes, targetNodes )
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

  SUBROUTINE Build_Lagrange( myPoly, N, controlNodeType, M, targetNodeType )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(out) :: myPoly
    INTEGER, INTENT(in)          :: N, M
    INTEGER, INTENT(in)          :: controlNodeType, targetNodeType
    ! Local
    REAL(prec) :: q(0:M)
    
      myPoly % N  = N
      myPoly % M  = M
     
      ALLOCATE( myPoly % controlPoints(0:N), &
                myPoly % bWeights(0:N), &
                myPoly % qWeights(0:N), &
                myPoly % targetPoints(0:M), &
                myPoly % iMatrix(0:N,0:M), &
                myPoly % dMatrix(0:N,0:N), &
                myPoly % bMatrix(0:N,0:1))

      myPoly % controlPoints = 0.0_prec
      myPoly % bWeights = 0.0_prec
      myPoly % qWeights = 0.0_prec
      myPoly % targetPoints  = 0.0_prec
      myPoly % iMatrix = 0.0_prec
      myPoly % dMatrix = 0.0_prec
      myPoly % bMatrix = 0.0_prec
      
     
      IF(controlNodeType == GAUSS .OR. controlNodeType == GAUSS_LOBATTO)THEN

        CALL LegendreQuadrature(N, & 
                                myPoly % controlPoints, &
                                myPoly % qWeights, &
                                controlNodeType )

      ELSEIF(controlNodeType == UNIFORM)THEN

        myPoly % controlPoints = UniformPoints(-1.0_prec,1.0_prec,N)
        myPoly % qWeights = 2.0_prec/REAL(N,prec)

      ENDIF

      ! Target Points
      IF(targetNodeType == GAUSS .OR. targetNodeType == GAUSS_LOBATTO)THEN

        CALL LegendreQuadrature(N, & 
                                myPoly % targetPoints, &
                                q, &
                                targetNodeType )

      ELSEIF(targetNodeType == UNIFORM)THEN

        myPoly % targetPoints = UniformPoints(-1.0_prec,1.0_prec,M)

      ENDIF


      CALL myPoly % CalculateBarycentricWeights( )
      CALL myPoly % CalculateInterpolationMatrix( )
      CALL myPoly % CalculateDerivativeMatrix( )
      myPoly % bMatrix(0:N,0) = myPoly % CalculateLagrangePolynomials( -1.0_prec)
      myPoly % bMatrix(0:N,1) = myPoly % CalculateLagrangePolynomials( 1.0_prec )


#ifdef GPU
!      IF( N > 7 )THEN
!         CALL Logging( WARN, 'Number of control points > 7 not fully supported for 3-D SEM operations.' )
!      ENDIF
      CALL hipCheck(hipMalloc(myPoly % bWeights_dev, SIZEOF(myPoly % bWeights)))
      CALL hipCheck(hipMalloc(myPoly % qWeights_dev, SIZEOF(myPoly % qWeights)))
      CALL hipCheck(hipMalloc(myPoly % iMatrix_dev, SIZEOF(myPoly % iMatrix)))
      CALL hipCheck(hipMalloc(myPoly % dMatrix_dev, SIZEOF(myPoly % dMatrix)))
      CALL hipCheck(hipMalloc(myPoly % bMatrix_dev, SIZEOF(myPoly % bMatrix)))

      CALL myPoly % UpdateDevice()
#endif
 
 END SUBROUTINE Build_Lagrange

! ================================================================================================ !
!
! Trash_Lagrange
!
!   A manual destructor for the Lagrange class that deallocates the memory held by its attributes. 
!
! ================================================================================================ ! 

  SUBROUTINE Trash_Lagrange( myPoly )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly

      DEALLOCATE( myPoly % controlPoints, &
                  myPoly % bWeights, &
                  myPoly % targetPoints, &
                  myPoly % iMatrix, &
                  myPoly % dMatrix, &
                  myPoly % bMatrix )
#ifdef GPU
      CALL hipCheck(hipFree(myPoly % bWeights_dev))
      CALL hipCheck(hipFree(myPoly % qWeights_dev))
      CALL hipCheck(hipFree(myPoly % iMatrix_dev))
      CALL hipCheck(hipFree(myPoly % dMatrix_dev))
      CALL hipCheck(hipFree(myPoly % bMatrix_dev))
#endif

  END SUBROUTINE Trash_Lagrange

  SUBROUTINE UpdateDevice_Lagrange(myPoly)
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly

      CALL hipCheck(hipMemcpy(myPoly % bWeights_dev, &
                                c_loc(myPoly % bWeights), &
                                SIZEOF(myPoly % bWeights), &
                                hipMemcpyHostToDevice))

      CALL hipCheck(hipMemcpy(myPoly % qWeights_dev, &
                                c_loc(myPoly % qWeights), &
                                SIZEOF(myPoly % qWeights), &
                                hipMemcpyHostToDevice))
                        
      CALL hipCheck(hipMemcpy(myPoly % iMatrix_dev, &
                                c_loc(myPoly % iMatrix), &
                                SIZEOF(myPoly % iMatrix), &
                                hipMemcpyHostToDevice))


      CALL hipCheck(hipMemcpy(myPoly % dMatrix_dev, &
                                c_loc(myPoly % dMatrix), &
                                SIZEOF(myPoly % dMatrix), &
                                hipMemcpyHostToDevice))

      CALL hipCheck(hipMemcpy(myPoly % bMatrix_dev, &
                                c_loc(myPoly % bMatrix), &
                                SIZEOF(myPoly % bMatrix), &
                                hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_Lagrange

  SUBROUTINE UpdateHost_Lagrange(myPoly)
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly

      CALL hipCheck(hipMemcpy(c_loc(myPoly % bWeights), &
                                myPoly % bWeights_dev, &
                                SIZEOF(myPoly % bWeights), &
                                hipMemcpyDeviceToHost))

      CALL hipCheck(hipMemcpy(c_loc(myPoly % qWeights), &
                                myPoly % qWeights_dev, &
                                SIZEOF(myPoly % qWeights), &
                                hipMemcpyDeviceToHost))
                        
      CALL hipCheck(hipMemcpy(c_loc(myPoly % iMatrix), &
                                myPoly % iMatrix_dev, &
                                SIZEOF(myPoly % iMatrix), &
                                hipMemcpyDeviceToHost))

      CALL hipCheck(hipMemcpy(c_loc(myPoly % dMatrix), &
                                myPoly % dMatrix_dev, &
                                SIZEOF(myPoly % dMatrix), &
                                hipMemcpyDeviceToHost))

      CALL hipCheck(hipMemcpy(c_loc(myPoly % bMatrix), &
                                myPoly % bMatrix_dev, &
                                SIZEOF(myPoly % bMatrix), &
                                hipMemcpyDeviceToHost))

  END SUBROUTINE UpdateHost_Lagrange

!! ================================================================================================ !
!!
!! Interpolate_1D_Lagrange
!!
!!   Interpolates an array of data onto a desired location using Lagrange interpolation.
!!
!!   Usage :
!!
!!     TYPE(Lagrange) :: interp
!!     REAL(prec)     :: f(0:interp % N) 
!!     REAL(prec)     :: sE 
!!     REAL(prec)     :: fAtSE 
!!         
!!       fAtSE = interp % Interpolate_1D( f, sE ) 
!! 
!!   Parameters :
!!  
!!     myPoly (in)
!!       A previously constructed Lagrange data-structure.
!!
!!     f (in)  
!!       An array of function nodal values located at the native interpolation nodes.
!!
!!     sE (in) 
!!       The location where you want to interpolate to.
!!
!!     fAtSE (out)
!!       The interpolant evaluated at sE.  
!!
!! ================================================================================================ ! 
!
!  FUNCTION Interpolate_1D_Lagrange( myPoly, f, sE ) RESULT( interpF )  
!    IMPLICIT NONE
!    CLASS(Lagrange) :: myPoly
!    REAL(prec)      :: sE
!    REAL(prec)      :: f(0:myPoly % N)
!    REAL(prec)      :: interpF
!    ! Local
!    REAL(prec) :: lAtS(0:myPoly % N)
!    INTEGER    :: i
!   
!      lAtS = myPoly % CalculateLagrangePolynomials( sE )
!
!      interpF = 0.0_prec
!      DO i = 0, myPoly % N
!        interpF = interpF + lAtS(i)*f(i)
!      ENDDO
!    
!  END FUNCTION Interpolate_1D_Lagrange
!
!! ================================================================================================ !
!!
!! Interpolate_2D_Lagrange
!!
!!   Interpolates an array of data onto a desired location using Lagrange interpolation.
!!
!!   Usage :
!!
!!     TYPE(Lagrange) :: interp
!!     REAL(prec)     :: f(0:interp % N,0:interp % N) 
!!     REAL(prec)     :: sE(1:2) 
!!     REAL(prec)     :: fAtSE 
!!         
!!       fAtSE = interp % Interpolate_2D( f, sE ) 
!! 
!!   Parameters :
!!  
!!     myPoly (in)
!!       A previously constructed Lagrange data-structure.
!!
!!     f (in)  
!!       An array of function nodal values located at the native interpolation nodes.
!!
!!     sE (in) 
!!       The location where you want to interpolate to.
!!
!!     fAtSE (out)
!!       The interpolant evaluated at sE.  
!!
!! ================================================================================================ ! 
!
!  FUNCTION Interpolate_2D_Lagrange( myPoly, f, sE ) RESULT( interpF )  
!    IMPLICIT NONE
!    CLASS(Lagrange) :: myPoly
!    REAL(prec)      :: sE(1:2)
!    REAL(prec)      :: f(0:myPoly % N, 0:myPoly % N)
!    REAL(prec)      :: interpF
!    ! Local
!    REAL(prec) :: fj
!    REAL(prec) :: ls(0:myPoly % N)
!    REAL(prec) :: lp(0:myPoly % N)
!    INTEGER    :: i, j
!
!      ls = myPoly % CalculateLagrangePolynomials( sE(1) ) 
!      lp = myPoly % CalculateLagrangePolynomials( sE(2) )
!      
!      interpF = 0.0_prec
!      DO j = 0, myPoly % N
!     
!        fj = 0.0_prec
!        DO i = 0, myPoly % N
!          fj = fj + f(i,j)*ls(i)
!        ENDDO
!            
!        interpF = interpF + fj*lp(j)
!
!      ENDDO
!      
! END FUNCTION Interpolate_2D_Lagrange
!
!! ================================================================================================ !
!!
!! Interpolate_3D_Lagrange
!!
!!   Interpolates an array of data onto a desired location using Lagrange interpolation.
!!
!!   Usage :
!!
!!     TYPE(Lagrange) :: interp
!!     REAL(prec)     :: f(0:interp % N,0:interp % N,0:interp % N) 
!!     REAL(prec)     :: sE(1:3) 
!!     REAL(prec)     :: fAtSE 
!!         
!!       fAtSE = interp % Interpolate_3D( f, sE ) 
!! 
!!   Parameters :
!!  
!!     myPoly (in)
!!       A previously constructed Lagrange data-structure.
!!
!!     f (in)  
!!       An array of function nodal values located at the native interpolation nodes.
!!
!!     sE (in) 
!!       The location where you want to interpolate to.
!!
!!     fAtSE (out)
!!       The interpolant evaluated at sE.  
!!
!! ================================================================================================ ! 
!
!  FUNCTION Interpolate_3D_Lagrange( myPoly, f, sE ) RESULT( interpF )  
!    IMPLICIT NONE
!    CLASS(Lagrange) :: myPoly
!    REAL(prec)      :: sE(1:3)
!    REAL(prec)      :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N)
!    REAL(prec)      :: interpF
!    ! Local
!    REAL(prec) :: fjk, fk
!    REAL(prec) :: ls(0:myPoly % N)
!    REAL(prec) :: lp(0:myPoly % N)
!    REAL(prec) :: lq(0:myPoly % N)
!    INTEGER    ::  i, j, k
!
!      ls = myPoly % CalculateLagrangePolynomials( sE(1) ) 
!      lp = myPoly % CalculateLagrangePolynomials( sE(2) )
!      lq = myPoly % CalculateLagrangePolynomials( sE(3) )
!      
!      interpF = 0.0_prec
!      DO k = 0, myPoly % N
!      
!         fk = 0.0_prec
!         DO j = 0, myPoly % N
!         
!            fjk = 0.0_prec
!            DO i = 0, myPoly % N
!               fjk = fjk + f(i,j,k)*ls(i)
!            ENDDO
!            
!            fk = fk + fjk*lp(j)
!         ENDDO
!         
!         interpF = interpF + fk*lq(k)
!      ENDDO
!      
!  END FUNCTION Interpolate_3D_Lagrange
!
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

  SUBROUTINE ScalarGridInterp_1D_cpu( myPoly, f, fInterp, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: fInterp(0:myPoly % M, 1:nVariables, 1:nElements)
    ! Local
    INTEGER :: iVar, iEl, i, ii

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO i = 0, myPoly % M
            fInterp(i,iVar,iEl) = 0.0_prec
            DO ii = 0, myPoly % N
              fInterp(i,iVar,iEl) = fInterp(i,iVar,iEl) + myPoly % iMatrix(ii,i)*f(ii,iVar,iEl)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE ScalarGridInterp_1D_cpu

  SUBROUTINE ScalarGridInterp_1D_gpu( myPoly, f_dev, fInterp_dev, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in) :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)  :: f_dev
    TYPE(c_ptr), INTENT(out) :: fInterp_dev

      CALL ScalarGridInterp_1D_gpu_wrapper(myPoly % iMatrix_dev, &
                                                f_dev, fInterp_dev, &
                                                myPoly % N, myPoly % M, &
                                                nVariables, nElements)

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

  SUBROUTINE ScalarGridInterp_2D_cpu( myPoly, f, fNew, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: fNew(0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
    ! Local
    INTEGER :: i, j, ii, jj, p, iEl, iVar
    REAL(prec) :: fi, fij
   
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % M
            DO i = 0, myPoly % M
             
              fij = 0.0_prec
              DO jj = 0, myPoly % N
                fi = 0.0_prec
                DO ii = 0, myPoly % N
                  fi = fi + f(ii,jj,iVar,iEl)*myPoly % iMatrix(ii,i)
                ENDDO
                fij = fij + fi*myPoly % iMatrix(jj,j)
              ENDDO
              fNew(i,j,iVar,iEl) = fij
                   
            ENDDO
          ENDDO
        ENDDO   
      ENDDO


  END SUBROUTINE ScalarGridInterp_2D_cpu
!
  SUBROUTINE ScalarGridInterp_2D_gpu( myPoly, f_dev, fInterp_dev, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in) :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)  :: f_dev
    TYPE(c_ptr), INTENT(out) :: fInterp_dev

      CALL ScalarGridInterp_2D_gpu_wrapper(myPoly % iMatrix_dev, &
                                           f_dev, fInterp_dev, &
                                           myPoly % N, myPoly % M, &
                                           nVariables, nElements)

  END SUBROUTINE ScalarGridInterp_2D_gpu

  SUBROUTINE VectorGridInterp_2D_cpu( myPoly, f, fNew, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:2,0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: fNew(1:2,0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
    ! Local
    INTEGER :: i, j, ii, jj, p, iEl, iVar
    REAL(prec) :: fi(1:2)
   
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % M
            DO i = 0, myPoly % M
             
              fNew(1,i,j,iVar,iEl) = 0.0_prec
              fNew(2,i,j,iVar,iEl) = 0.0_prec

              DO jj = 0, myPoly % N
                   
                fi(1:2) = 0.0_prec
                DO ii = 0, myPoly % N
                  fi(1:2)= fi(1:2) + f(1:2,ii,jj,iVar,iEl)*myPoly % iMatrix(ii,i)
                ENDDO
                      
                fNew(1:2,i,j,iVar,iEl) = fNew(1:2,i,j,iVar,iEl) + fi(1:2)*myPoly % iMatrix(jj,j)

              ENDDO
                   
            ENDDO
          ENDDO
        ENDDO   
      ENDDO


  END SUBROUTINE VectorGridInterp_2D_cpu
!
  SUBROUTINE VectorGridInterp_2D_gpu( myPoly, f_dev, fInterp_dev, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in) :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)  :: f_dev
    TYPE(c_ptr), INTENT(out) :: fInterp_dev

      CALL VectorGridInterp_2D_gpu_wrapper(myPoly % iMatrix_dev, &
                                                f_dev, fInterp_dev, &
                                                myPoly % N, myPoly % M, &
                                                nVariables, nElements)

  END SUBROUTINE VectorGridInterp_2D_gpu

  SUBROUTINE TensorGridInterp_2D_cpu( myPoly, f, fNew, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:2,1:2,0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: fNew(1:2,1:2,0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
    ! Local
    INTEGER :: i, j, ii, jj, p, iEl, iVar
    REAL(prec) :: fi(1:2,1:2)
   
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % M
            DO i = 0, myPoly % M
             
              fNew(1:2,1:2,i,j,iVar,iEl) = 0.0_prec

              DO jj = 0, myPoly % N
                   
                fi(1:2,1:2) = 0.0_prec
                DO ii = 0, myPoly % N
                  fi(1:2,1:2)= fi(1:2,1:2) + f(1:2,1:2,ii,jj,iVar,iEl)*myPoly % iMatrix(ii,i)
                ENDDO
                      
                fNew(1:2,1:2,i,j,iVar,iEl) = fNew(1:2,1:2,i,j,iVar,iEl) + fi(1:2,1:2)*myPoly % iMatrix(jj,j)

              ENDDO
                   
            ENDDO
          ENDDO
        ENDDO   
      ENDDO


  END SUBROUTINE TensorGridInterp_2D_cpu
!
  SUBROUTINE TensorGridInterp_2D_gpu( myPoly, f_dev, fInterp_dev, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in) :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)  :: f_dev
    TYPE(c_ptr), INTENT(out) :: fInterp_dev

      CALL TensorGridInterp_2D_gpu_wrapper(myPoly % iMatrix_dev, &
                                                f_dev, fInterp_dev, &
                                                myPoly % N, myPoly % M, &
                                                nVariables, nElements)

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

  SUBROUTINE ScalarGridInterp_3D_cpu( myPoly, f, fInterp, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: fInterp(0:myPoly % M, 0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
    ! Local
    INTEGER :: iEl, iVar, i, j, k, ii, jj, kk
    REAL(prec) :: fi, fij
   
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, myPoly % M
            DO j = 0, myPoly % M
              DO i = 0, myPoly % M
               
                fInterp(i,j,k,iVar,iEl) = 0.0_prec
                DO kk = 0, myPoly % N
                  fij = 0.0_prec
                  DO jj = 0, myPoly % N
                    fi = 0.0_prec
                    DO ii = 0, myPoly % N
                      fi = fi + f(ii,jj,kk,iVar,iEl)*myPoly % iMatrix(ii,i)
                    ENDDO
                    fij = fij + fi*myPoly % iMatrix(jj,j)
                  ENDDO
                  fInterp(i,j,k,iVar,iEl) = fInterp(i,j,k,iVar,iEl) + fij*myPoly % iMatrix(kk,k)
                ENDDO

              ENDDO
            ENDDO
          ENDDO
        ENDDO   
      ENDDO

  END SUBROUTINE ScalarGridInterp_3D_cpu
!
  SUBROUTINE ScalarGridInterp_3D_gpu( myPoly, f_dev, fInterp_dev, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in) :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)  :: f_dev
    TYPE(c_ptr), INTENT(out) :: fInterp_dev

      CALL ScalarGridInterp_3D_gpu_wrapper(myPoly % iMatrix_dev, &
                                                f_dev, fInterp_dev, &
                                                myPoly % N, myPoly % M, &
                                                nVariables, nElements)

  END SUBROUTINE ScalarGridInterp_3D_gpu

  SUBROUTINE VectorGridInterp_3D_cpu( myPoly, f, fInterp, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:3,0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: fInterp(1:3,0:myPoly % M, 0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
    ! Local
    INTEGER :: iEl, iVar, i, j, k, ii, jj, kk
    REAL(prec) :: fi(1:3), fij(1:3)
   
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, myPoly % M
            DO j = 0, myPoly % M
              DO i = 0, myPoly % M
               
                fInterp(1:3,i,j,k,iVar,iEl) = 0.0_prec
                DO kk = 0, myPoly % N
                  fij(1:3) = 0.0_prec
                  DO jj = 0, myPoly % N
                    fi(1:3) = 0.0_prec
                    DO ii = 0, myPoly % N
                      fi(1:3) = fi(1:3) + f(1:3,ii,jj,kk,iVar,iEl)*myPoly % iMatrix(ii,i)
                    ENDDO
                    fij(1:3) = fij(1:3) + fi(1:3)*myPoly % iMatrix(jj,j)
                  ENDDO
                  fInterp(1:3,i,j,k,iVar,iEl) = fInterp(1:3,i,j,k,iVar,iEl) + fij(1:3)*myPoly % iMatrix(kk,k)
                ENDDO

              ENDDO
            ENDDO
          ENDDO
        ENDDO   
      ENDDO

  END SUBROUTINE VectorGridInterp_3D_cpu
!
  SUBROUTINE VectorGridInterp_3D_gpu( myPoly, f_dev, fInterp_dev, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in) :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)  :: f_dev
    TYPE(c_ptr), INTENT(out) :: fInterp_dev

      CALL VectorGridInterp_3D_gpu_wrapper(myPoly % iMatrix_dev, &
                                                f_dev, fInterp_dev, &
                                                myPoly % N, myPoly % M, &
                                                nVariables, nElements)

  END SUBROUTINE VectorGridInterp_3D_gpu
!
  SUBROUTINE TensorGridInterp_3D_cpu( myPoly, f, fInterp, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:3,1:3,0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: fInterp(1:3,1:3,0:myPoly % M, 0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
    ! Local
    INTEGER :: iEl, iVar, i, j, k, ii, jj, kk
    REAL(prec) :: fi(1:3,1:3), fij(1:3,1:3)
   
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, myPoly % M
            DO j = 0, myPoly % M
              DO i = 0, myPoly % M
               
                fInterp(1:3,1:3,i,j,k,iVar,iEl) = 0.0_prec
                DO kk = 0, myPoly % N
                  fij(1:3,1:3) = 0.0_prec
                  DO jj = 0, myPoly % N
                    fi(1:3,1:3) = 0.0_prec
                    DO ii = 0, myPoly % N
                      fi(1:3,1:3) = fi(1:3,1:3) + f(1:3,1:3,ii,jj,kk,iVar,iEl)*myPoly % iMatrix(ii,i)
                    ENDDO
                    fij(1:3,1:3) = fij(1:3,1:3) + fi(1:3,1:3)*myPoly % iMatrix(jj,j)
                  ENDDO
                  fInterp(1:3,1:3,i,j,k,iVar,iEl) = fInterp(1:3,1:3,i,j,k,iVar,iEl) + fij(1:3,1:3)*myPoly % iMatrix(kk,k)
                ENDDO

              ENDDO
            ENDDO
          ENDDO
        ENDDO   
      ENDDO

  END SUBROUTINE TensorGridInterp_3D_cpu
!
  SUBROUTINE TensorGridInterp_3D_gpu( myPoly, f_dev, fInterp_dev, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in) :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)  :: f_dev
    TYPE(c_ptr), INTENT(out) :: fInterp_dev

      CALL TensorGridInterp_3D_gpu_wrapper(myPoly % iMatrix_dev, &
                                                f_dev, fInterp_dev, &
                                                myPoly % N, myPoly % M, &
                                                nVariables, nElements)

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

  SUBROUTINE Derivative_1D_cpu( myPoly, f, df, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: df(0:myPoly % N, 1:nVariables, 1:nElements)
    ! Local
    INTEGER :: i, ii, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO i = 0, myPoly % N
    
            df(i,iVar,iEl) = 0.0_prec 
            DO ii = 0, myPoly % N
              df(i,iVar,iEl) = df(i,iVar,iEl) + myPoly % dMatrix(ii,i)*f(ii,iVar,iEl)
            ENDDO
    
          ENDDO
        ENDDO
      ENDDO


  END SUBROUTINE Derivative_1D_cpu

  SUBROUTINE Derivative_1D_gpu( myPoly, f_dev, df_dev, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in) :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)  :: f_dev
    TYPE(c_ptr), INTENT(out) :: df_dev

      CALL Derivative_1D_gpu_wrapper(myPoly % dMatrix_dev, &
                                         f_dev, df_dev, &
                                         myPoly % N, &
                                         nVariables, nElements)

  END SUBROUTINE Derivative_1D_gpu
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
  SUBROUTINE ScalarGradient_2D_cpu( myPoly, f, gradF, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: gradF(1:2, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % N
            DO i = 0, myPoly % N
    
              gradF(1,i,j,iVar,iEl) = 0.0_prec 
              gradF(2,i,j,iVar,iEl) = 0.0_prec 
              DO ii = 0, myPoly % N
                gradF(1,i,j,iVar,iEl) = gradF(1,i,j,iVar,iEl) + myPoly % dMatrix(ii,i)*f(ii,j,iVar,iEl)
                gradF(2,i,j,iVar,iEl) = gradF(2,i,j,iVar,iEl) + myPoly % dMatrix(ii,j)*f(i,ii,iVar,iEl)
              ENDDO
    
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE ScalarGradient_2D_cpu

  SUBROUTINE ScalarGradient_2D_gpu( myPoly, f_dev, gradF_dev, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)         :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)     :: f_dev
    TYPE(c_ptr), INTENT(out)    :: gradF_dev
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      CALL ScalarGradient_2D_gpu_wrapper(myPoly % dMatrix_dev, &
                                         f_dev, gradF_dev, myPoly % N, &
                                         nVariables, nElements)

  END SUBROUTINE ScalarGradient_2D_gpu
!
  SUBROUTINE VectorGradient_2D_cpu( myPoly, f, gradF, nVariables, nElements )
  !
  ! Input : Vector(1:2,...)
  ! Output : Tensor(1:2,1:2,....)
  !          > Tensor(1,1) = d/ds1( Vector(1,...) )
  !          > Tensor(2,1) = d/ds1( Vector(2,...) )
  !          > Tensor(1,2) = d/ds2( Vector(1,...) )
  !          > Tensor(2,2) = d/ds2( Vector(2,...) )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:2,0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: gradF(1:2,1:2,0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % N
            DO i = 0, myPoly % N
    
              gradF(1,1,i,j,iVar,iEl) = 0.0_prec 
              gradF(2,1,i,j,iVar,iEl) = 0.0_prec 
              gradF(1,2,i,j,iVar,iEl) = 0.0_prec 
              gradF(2,2,i,j,iVar,iEl) = 0.0_prec 
              DO ii = 0, myPoly % N
                gradF(1,1,i,j,iVar,iEl) = gradF(1,1,i,j,iVar,iEl) + myPoly % dMatrix(ii,i)*f(1,ii,j,iVar,iEl)
                gradF(2,1,i,j,iVar,iEl) = gradF(2,1,i,j,iVar,iEl) + myPoly % dMatrix(ii,i)*f(2,ii,j,iVar,iEl)
                gradF(1,2,i,j,iVar,iEl) = gradF(1,2,i,j,iVar,iEl) + myPoly % dMatrix(ii,j)*f(1,i,ii,iVar,iEl)
                gradF(2,2,i,j,iVar,iEl) = gradF(2,2,i,j,iVar,iEl) + myPoly % dMatrix(ii,j)*f(2,i,ii,iVar,iEl)
              ENDDO
    
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE VectorGradient_2D_cpu

  SUBROUTINE VectorGradient_2D_gpu( myPoly, f_dev, gradF_dev, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)         :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)     :: f_dev
    TYPE(c_ptr), INTENT(out)    :: gradF_dev
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      CALL VectorGradient_2D_gpu_wrapper(myPoly % dMatrix_dev, &
                                         f_dev, gradF_dev, myPoly % N, &
                                         nVariables, nElements)

  END SUBROUTINE VectorGradient_2D_gpu

  SUBROUTINE VectorDivergence_2D_cpu( myPoly, f, dF, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:2,0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: dF(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % N
            DO i = 0, myPoly % N
    
              dF(i,j,iVar,iEl) = 0.0_prec 
              DO ii = 0, myPoly % N
                dF(i,j,iVar,iEl) = dF(i,j,iVar,iEl) + myPoly % dMatrix(ii,i)*f(1,ii,j,iVar,iEl) +&
                                                      myPoly % dMatrix(ii,j)*f(2,i,ii,iVar,iEl)
              ENDDO
    
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE VectorDivergence_2D_cpu

  SUBROUTINE VectorDivergence_2D_gpu( myPoly, f_dev, dF_dev, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)         :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)     :: f_dev
    TYPE(c_ptr), INTENT(out)    :: dF_dev
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      CALL VectorDivergence_2D_gpu_wrapper(myPoly % dMatrix_dev, &
                                         f_dev, dF_dev, myPoly % N, &
                                         nVariables, nElements)

  END SUBROUTINE VectorDivergence_2D_gpu

  SUBROUTINE VectorCurl_2D_cpu( myPoly, f, dF, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:2,0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: dF(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % N
            DO i = 0, myPoly % N
    
              dF(i,j,iVar,iEl) = 0.0_prec 
              DO ii = 0, myPoly % N
                dF(i,j,iVar,iEl) = dF(i,j,iVar,iEl) + myPoly % dMatrix(ii,j)*f(1,i,ii,iVar,iEl) -&
                                                      myPoly % dMatrix(ii,i)*f(2,ii,j,iVar,iEl)
              ENDDO
    
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE VectorCurl_2D_cpu

  SUBROUTINE VectorCurl_2D_gpu( myPoly, f_dev, dF_dev, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)         :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)     :: f_dev
    TYPE(c_ptr), INTENT(out)    :: dF_dev
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      CALL VectorCurl_2D_gpu_wrapper(myPoly % dMatrix_dev, &
                                         f_dev, dF_dev, myPoly % N, &
                                         nVariables, nElements)

  END SUBROUTINE VectorCurl_2D_gpu

  SUBROUTINE ScalarGradient_3D_cpu( myPoly, f, gradF, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: gradF(1:3, 0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    ! Local
    INTEGER    :: i, j, k, ii, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, myPoly % N
            DO j = 0, myPoly % N
              DO i = 0, myPoly % N
    
                gradF(1,i,j,k,iVar,iEl) = 0.0_prec 
                gradF(2,i,j,k,iVar,iEl) = 0.0_prec 
                gradF(3,i,j,k,iVar,iEl) = 0.0_prec 
                DO ii = 0, myPoly % N
                  gradF(1,i,j,k,iVar,iEl) = gradF(1,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,i)*f(ii,j,k,iVar,iEl)
                  gradF(2,i,j,k,iVar,iEl) = gradF(2,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,j)*f(i,ii,k,iVar,iEl)
                  gradF(3,i,j,k,iVar,iEl) = gradF(3,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,k)*f(i,j,ii,iVar,iEl)
                ENDDO
    
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE ScalarGradient_3D_cpu

  SUBROUTINE ScalarGradient_3D_gpu( myPoly, f_dev, gradF_dev, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)         :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)     :: f_dev
    TYPE(c_ptr), INTENT(out)    :: gradF_dev
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      CALL ScalarGradient_3D_gpu_wrapper(myPoly % dMatrix_dev, &
                                         f_dev, gradF_dev, myPoly % N, &
                                         nVariables, nElements)

  END SUBROUTINE ScalarGradient_3D_gpu
!
  SUBROUTINE VectorGradient_3D_cpu( myPoly, f, gradF, nVariables, nElements )
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
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:3,0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: gradF(1:3,1:3,0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    ! Local
    INTEGER    :: i, j, k, ii, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, myPoly % N
            DO j = 0, myPoly % N
              DO i = 0, myPoly % N
    
                gradF(1,1,i,j,k,iVar,iEl) = 0.0_prec 
                gradF(2,1,i,j,k,iVar,iEl) = 0.0_prec 
                gradF(3,1,i,j,k,iVar,iEl) = 0.0_prec 
                gradF(1,2,i,j,k,iVar,iEl) = 0.0_prec 
                gradF(2,2,i,j,k,iVar,iEl) = 0.0_prec 
                gradF(3,2,i,j,k,iVar,iEl) = 0.0_prec 
                gradF(1,3,i,j,k,iVar,iEl) = 0.0_prec 
                gradF(2,3,i,j,k,iVar,iEl) = 0.0_prec 
                gradF(3,3,i,j,k,iVar,iEl) = 0.0_prec 
                DO ii = 0, myPoly % N
                  gradF(1,1,i,j,k,iVar,iEl) = gradF(1,1,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,i)*f(1,ii,j,k,iVar,iEl)
                  gradF(2,1,i,j,k,iVar,iEl) = gradF(2,1,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,i)*f(2,ii,j,k,iVar,iEl)
                  gradF(3,1,i,j,k,iVar,iEl) = gradF(3,1,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,i)*f(3,ii,j,k,iVar,iEl)
                  gradF(1,2,i,j,k,iVar,iEl) = gradF(1,2,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,j)*f(1,i,ii,k,iVar,iEl)
                  gradF(2,2,i,j,k,iVar,iEl) = gradF(2,2,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,j)*f(2,i,ii,k,iVar,iEl)
                  gradF(3,2,i,j,k,iVar,iEl) = gradF(3,2,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,j)*f(3,i,ii,k,iVar,iEl)
                  gradF(1,3,i,j,k,iVar,iEl) = gradF(1,3,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,k)*f(1,i,j,ii,iVar,iEl)
                  gradF(2,3,i,j,k,iVar,iEl) = gradF(2,3,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,k)*f(2,i,j,ii,iVar,iEl)
                  gradF(3,3,i,j,k,iVar,iEl) = gradF(3,3,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,k)*f(3,i,j,ii,iVar,iEl)
                ENDDO
    
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE VectorGradient_3D_cpu

  SUBROUTINE VectorGradient_3D_gpu( myPoly, f_dev, gradF_dev, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)         :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)     :: f_dev
    TYPE(c_ptr), INTENT(out)    :: gradF_dev
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      CALL VectorGradient_3D_gpu_wrapper(myPoly % dMatrix_dev, &
                                         f_dev, gradF_dev, myPoly % N, &
                                         nVariables, nElements)

  END SUBROUTINE VectorGradient_3D_gpu

  SUBROUTINE VectorDivergence_3D_cpu( myPoly, f, dF, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:3,0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: dF(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    ! Local
    INTEGER    :: i, j, k, ii, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, myPoly % N
            DO j = 0, myPoly % N
              DO i = 0, myPoly % N
    
                dF(i,j,k,iVar,iEl) = 0.0_prec 
                DO ii = 0, myPoly % N
                  dF(i,j,k,iVar,iEl) = dF(i,j,k,iVar,iEl) + myPoly % dMatrix(ii,i)*f(1,ii,j,k,iVar,iEl) +&
                                                            myPoly % dMatrix(ii,j)*f(2,i,ii,k,iVar,iEl) +&
                                                            myPoly % dMatrix(ii,k)*f(3,i,j,ii,iVar,iEl)
                ENDDO
    
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE VectorDivergence_3D_cpu

  SUBROUTINE VectorDivergence_3D_gpu( myPoly, f_dev, dF_dev, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)         :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)     :: f_dev
    TYPE(c_ptr), INTENT(out)    :: dF_dev
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      CALL VectorDivergence_3D_gpu_wrapper(myPoly % dMatrix_dev, &
                                         f_dev, dF_dev, myPoly % N, &
                                         nVariables, nElements)

  END SUBROUTINE VectorDivergence_3D_gpu

  SUBROUTINE VectorCurl_3D_cpu( myPoly, f, dF, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:3,0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: dF(1:3,0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    ! Local
    INTEGER    :: i, j, k, ii, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, myPoly % N
            DO j = 0, myPoly % N
              DO i = 0, myPoly % N
    
                dF(1,i,j,k,iVar,iEl) = 0.0_prec 
                dF(2,i,j,k,iVar,iEl) = 0.0_prec 
                dF(3,i,j,k,iVar,iEl) = 0.0_prec 
                DO ii = 0, myPoly % N
                  dF(1,i,j,k,iVar,iEl) = dF(1,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,j)*f(3,i,ii,k,iVar,iEl) -&
                                                                myPoly % dMatrix(ii,k)*f(2,i,j,ii,iVar,iEl)
                  dF(2,i,j,k,iVar,iEl) = dF(2,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,k)*f(1,i,j,ii,iVar,iEl) -&
                                                                myPoly % dMatrix(ii,i)*f(3,ii,j,k,iVar,iEl)
                  dF(3,i,j,k,iVar,iEl) = dF(3,i,j,k,iVar,iEl) + myPoly % dMatrix(ii,i)*f(2,ii,j,k,iVar,iEl) -&
                                                                myPoly % dMatrix(ii,j)*f(1,i,ii,k,iVar,iEl)
                ENDDO  
    
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE VectorCurl_3D_cpu

  SUBROUTINE VectorCurl_3D_gpu( myPoly, f_dev, dF_dev, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)         :: nVariables, nElements
    TYPE(c_ptr), INTENT(in)     :: f_dev
    TYPE(c_ptr), INTENT(out)    :: dF_dev
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      CALL VectorCurl_3D_gpu_wrapper(myPoly % dMatrix_dev, &
                                         f_dev, dF_dev, myPoly % N, &
                                         nVariables, nElements)

  END SUBROUTINE VectorCurl_3D_gpu

  ! /////////////////////////////// !
  ! Boundary Interpolation Routines !

  SUBROUTINE ScalarBoundaryInterp_1D_cpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)         :: nVariables, nElements
    REAL(prec), INTENT(in)      :: f(0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out)     :: fBound(1:nVariables, 1:2, 1:nElements)
    ! Local
    INTEGER :: i, j, ii, iVar, iEl
    REAL(prec) :: fb(1:2)

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          fb(1:2) = 0.0_prec
          
          DO ii = 0, myPoly % N
            fb(1) = fb(1) + myPoly % bMatrix(ii,1)*f(ii,iVar,iEl) ! East
            fb(2) = fb(2) + myPoly % bMatrix(ii,0)*f(ii,iVar,iEl) ! West
          ENDDO

          fBound(iVar,1:2,iEl) = fb(1:2)
          
        ENDDO
      ENDDO

  END SUBROUTINE ScalarBoundaryInterp_1D_cpu

  SUBROUTINE ScalarBoundaryInterp_1D_gpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER       , INTENT(in)  :: nVariables, nElements
    TYPE(c_ptr)    , INTENT(in)  :: f
    TYPE(c_ptr)    , INTENT(out)  :: fBound

      CALL ScalarBoundaryInterp_1D_gpu_wrapper(myPoly % bMatrix_dev, &
                                               f, fBound, myPoly % N, nVariables, nElements)

  END SUBROUTINE ScalarBoundaryInterp_1D_gpu

  SUBROUTINE ScalarBoundaryInterp_2D_cpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)         :: nVariables, nElements
    REAL(prec), INTENT(in)      :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out)     :: fBound(0:myPoly % N, 1:nVariables, 1:4, 1:nElements)
    ! Local
    INTEGER :: i, j, ii, iVar, iEl
    REAL(prec) :: fb(1:4)

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO i = 0, myPoly % N
          
            fb(1:4) = 0.0_prec
            
            DO ii = 0, myPoly % N
              fb(1) = fb(1) + myPoly % bMatrix(ii,0)*f(i,ii,iVar,iEl) ! South
              fb(2) = fb(2) + myPoly % bMatrix(ii,1)*f(ii,i,iVar,iEl) ! East
              fb(3) = fb(3) + myPoly % bMatrix(ii,1)*f(i,ii,iVar,iEl) ! North
              fb(4) = fb(4) + myPoly % bMatrix(ii,0)*f(ii,i,iVar,iEl) ! West
            ENDDO

            fBound(i,iVar,1:4,iEl) = fb(1:4)
            
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE ScalarBoundaryInterp_2D_cpu

  SUBROUTINE ScalarBoundaryInterp_2D_gpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER       , INTENT(in)  :: nVariables, nElements
    TYPE(c_ptr)    , INTENT(in)  :: f
    TYPE(c_ptr)    , INTENT(out)  :: fBound

      CALL ScalarBoundaryInterp_2D_gpu_wrapper(myPoly % bMatrix_dev, &
                                               f, fBound, myPoly % N, nVariables, nElements)

  END SUBROUTINE ScalarBoundaryInterp_2D_gpu

  SUBROUTINE VectorBoundaryInterp_2D_cpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER       , INTENT(in)  :: nVariables, nElements
    REAL(prec)    , INTENT(in)  :: f(1:2,0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec)    , INTENT(out)  :: fBound(1:2,0:myPoly % N, 1:nVariables, 1:4, 1:nElements)
    ! Local
    INTEGER :: i, j, ii, idir, iVar, iEl
    REAL(prec) :: fb(1:2,1:4)

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO i = 0, myPoly % N

            fb(1:2,1:4) = 0.0_prec
            DO ii = 0, myPoly % N
              DO idir = 1, 2
                fb(idir,1) = fb(idir,1) + myPoly % bMatrix(ii,0)*f(idir,i,ii,iVar,iEl) ! South
                fb(idir,2) = fb(idir,2) + myPoly % bMatrix(ii,1)*f(idir,ii,i,iVar,iEl) ! East
                fb(idir,3) = fb(idir,3) + myPoly % bMatrix(ii,1)*f(idir,i,ii,iVar,iEl) ! North
                fb(idir,4) = fb(idir,4) + myPoly % bMatrix(ii,0)*f(idir,ii,i,iVar,iEl) ! West
              ENDDO
            ENDDO

            DO idir = 1, 2
              fBound(idir,i,iVar,1:4,iEl) = fb(idir,1:4)
            ENDDO
            
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE VectorBoundaryInterp_2D_cpu

  SUBROUTINE VectorBoundaryInterp_2D_gpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER       , INTENT(in)  :: nVariables, nElements
    TYPE(c_ptr)    , INTENT(in)  :: f
    TYPE(c_ptr)    , INTENT(out)  :: fBound

      CALL VectorBoundaryInterp_2D_gpu_wrapper(myPoly % bMatrix_dev, &
                                               f, fBound, myPoly % N, nVariables, nElements)

  END SUBROUTINE VectorBoundaryInterp_2D_gpu

  SUBROUTINE TensorBoundaryInterp_2D_cpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER       , INTENT(in)  :: nVariables, nElements
    REAL(prec)    , INTENT(in)  :: f(1:2,1:2,0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec)    , INTENT(out)  :: fBound(1:2,1:2,0:myPoly % N, 1:nVariables, 1:4, 1:nElements)
    ! Local
    INTEGER :: i, j, ii, idir, jdir, iVar, iEl
    REAL(prec) :: fb(1:2,1:2,1:4)

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO i = 0, myPoly % N

            fb(1:2,1:2,1:4) = 0.0_prec
            DO ii = 0, myPoly % N
              DO jdir = 1, 2
                DO idir = 1, 2
                  fb(idir,jdir,1) = fb(idir,jdir,1) + myPoly % bMatrix(ii,0)*f(idir,jdir,i,ii,iVar,iEl) ! South
                  fb(idir,jdir,2) = fb(idir,jdir,2) + myPoly % bMatrix(ii,1)*f(idir,jdir,ii,i,iVar,iEl) ! East
                  fb(idir,jdir,3) = fb(idir,jdir,3) + myPoly % bMatrix(ii,1)*f(idir,jdir,i,ii,iVar,iEl) ! North
                  fb(idir,jdir,4) = fb(idir,jdir,4) + myPoly % bMatrix(ii,0)*f(idir,jdir,ii,i,iVar,iEl) ! West
                ENDDO
              ENDDO
            ENDDO

            DO jdir = 1, 2
              DO idir = 1, 2
                fBound(idir,jdir,i,iVar,1:4,iEl) = fb(idir,jdir,1:4)
              ENDDO
            ENDDO
            
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE TensorBoundaryInterp_2D_cpu

  SUBROUTINE TensorBoundaryInterp_2D_gpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)  :: nVariables, nElements
    TYPE(c_ptr)    , INTENT(in)  :: f
    TYPE(c_ptr)    , INTENT(out)  :: fBound

      CALL TensorBoundaryInterp_2D_gpu_wrapper(myPoly % bMatrix_dev, &
                                               f, fBound, myPoly % N, nVariables, nElements)

  END SUBROUTINE TensorBoundaryInterp_2D_gpu

  SUBROUTINE ScalarBoundaryInterp_3D_cpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER, INTENT(in)         :: nVariables, nElements
    REAL(prec), INTENT(in)      :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out)     :: fBound(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:6, 1:nElements)
    ! Local
    INTEGER :: i, j, ii, iVar, iEl
    REAL(prec) :: fb(1:6)

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % N
            DO i = 0, myPoly % N
            
              fb(1:6) = 0.0_prec
              
              DO ii = 0, myPoly % N
                fb(1) = fb(1) + myPoly % bMatrix(ii,0)*f(i,ii,j,iVar,iEl) ! South
                fb(2) = fb(2) + myPoly % bMatrix(ii,1)*f(ii,i,j,iVar,iEl) ! East
                fb(3) = fb(3) + myPoly % bMatrix(ii,1)*f(i,ii,j,iVar,iEl) ! North
                fb(4) = fb(4) + myPoly % bMatrix(ii,0)*f(ii,i,j,iVar,iEl) ! West
                fb(5) = fb(5) + myPoly % bMatrix(ii,0)*f(i,j,ii,iVar,iEl) ! Bottom
                fb(6) = fb(6) + myPoly % bMatrix(ii,1)*f(i,j,ii,iVar,iEl) ! Top
              ENDDO

              fBound(i,j,iVar,1:6,iEl) = fb(1:6)
              
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE ScalarBoundaryInterp_3D_cpu

  SUBROUTINE ScalarBoundaryInterp_3D_gpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER       , INTENT(in)  :: nVariables, nElements
    TYPE(c_ptr)    , INTENT(in)  :: f
    TYPE(c_ptr)    , INTENT(out)  :: fBound

      CALL ScalarBoundaryInterp_3D_gpu_wrapper(myPoly % bMatrix_dev, &
                                               f, fBound, myPoly % N, nVariables, nElements)

  END SUBROUTINE ScalarBoundaryInterp_3D_gpu

  SUBROUTINE VectorBoundaryInterp_3D_cpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER       , INTENT(in)  :: nVariables, nElements
    REAL(prec)    , INTENT(in)  :: f(1:3,0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec)    , INTENT(out)  :: fBound(1:3,0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:6, 1:nElements)
    ! Local
    INTEGER :: i, j, ii, idir, iVar, iEl
    REAL(prec) :: fb(1:3,1:6)

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % N
            DO i = 0, myPoly % N

              fb(1:3,1:6) = 0.0_prec
              DO ii = 0, myPoly % N
                DO idir = 1, 3
                  fb(idir,1) = fb(idir,1) + myPoly % bMatrix(ii,0)*f(idir,i,ii,j,iVar,iEl) ! South
                  fb(idir,2) = fb(idir,2) + myPoly % bMatrix(ii,1)*f(idir,ii,i,j,iVar,iEl) ! East
                  fb(idir,3) = fb(idir,3) + myPoly % bMatrix(ii,1)*f(idir,i,ii,j,iVar,iEl) ! North
                  fb(idir,4) = fb(idir,4) + myPoly % bMatrix(ii,0)*f(idir,ii,i,j,iVar,iEl) ! West
                  fb(idir,5) = fb(idir,5) + myPoly % bMatrix(ii,0)*f(idir,i,j,ii,iVar,iEl) ! Bottom
                  fb(idir,6) = fb(idir,6) + myPoly % bMatrix(ii,1)*f(idir,i,j,ii,iVar,iEl) ! Top
                ENDDO
              ENDDO

              DO idir = 1, 3
                fBound(idir,i,j,iVar,1:6,iEl) = fb(idir,1:6)
              ENDDO
              
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE VectorBoundaryInterp_3D_cpu

  SUBROUTINE VectorBoundaryInterp_3D_gpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER       , INTENT(in)  :: nVariables, nElements
    TYPE(c_ptr)    , INTENT(in)  :: f
    TYPE(c_ptr)    , INTENT(out)  :: fBound

      CALL VectorBoundaryInterp_3D_gpu_wrapper(myPoly % bMatrix_dev, &
                                               f, fBound, myPoly % N, nVariables, nElements)

  END SUBROUTINE VectorBoundaryInterp_3D_gpu

  SUBROUTINE TensorBoundaryInterp_3D_cpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER       , INTENT(in)  :: nVariables, nElements
    REAL(prec)    , INTENT(in)  :: f(1:3,1:3,0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec)    , INTENT(out)  :: fBound(1:3,1:3,0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:6, 1:nElements)
    ! Local
    INTEGER :: i, j, ii, idir, jdir, iVar, iEl
    REAL(prec) :: fb(1:3,1:3,1:6)

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % N
            DO i = 0, myPoly % N

              fb(1:3,1:3,1:6) = 0.0_prec
              DO ii = 0, myPoly % N
                DO jdir = 1, 3
                  DO idir = 1, 3
                    fb(idir,jdir,1) = fb(idir,jdir,1) + myPoly % bMatrix(ii,0)*f(idir,jdir,i,ii,j,iVar,iEl) ! South
                    fb(idir,jdir,2) = fb(idir,jdir,2) + myPoly % bMatrix(ii,1)*f(idir,jdir,ii,i,j,iVar,iEl) ! East
                    fb(idir,jdir,3) = fb(idir,jdir,3) + myPoly % bMatrix(ii,1)*f(idir,jdir,i,ii,j,iVar,iEl) ! North
                    fb(idir,jdir,4) = fb(idir,jdir,4) + myPoly % bMatrix(ii,0)*f(idir,jdir,ii,i,j,iVar,iEl) ! West
                    fb(idir,jdir,5) = fb(idir,jdir,5) + myPoly % bMatrix(ii,0)*f(idir,jdir,i,j,ii,iVar,iEl) ! Bottom
                    fb(idir,jdir,6) = fb(idir,jdir,6) + myPoly % bMatrix(ii,1)*f(idir,jdir,i,j,ii,iVar,iEl) ! Top
                  ENDDO
                ENDDO
              ENDDO

              DO jdir = 1, 3
                DO idir = 1, 3
                  fBound(idir,jdir,i,j,iVar,1:6,iEl) = fb(idir,jdir,1:6)
                ENDDO
              ENDDO
              
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE TensorBoundaryInterp_3D_cpu

  SUBROUTINE TensorBoundaryInterp_3D_gpu( myPoly, f, fBound, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly
    INTEGER       , INTENT(in)  :: nVariables, nElements
    TYPE(c_ptr)    , INTENT(in)  :: f
    TYPE(c_ptr)    , INTENT(out)  :: fBound

      CALL TensorBoundaryInterp_3D_gpu_wrapper(myPoly % bMatrix_dev, &
                                               f, fBound, myPoly % N, nVariables, nElements)

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

  SUBROUTINE CalculateBarycentricWeights( myPoly )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly
    ! Local
    INTEGER :: i, j
   
      DO i = 0, myPoly % N
        myPoly % bWeights(i) = 1.0_prec
      ENDDO

      ! Computes the product w_k = w_k*(s_k - s_j), k /= j
      DO j = 1, myPoly % N
        DO i = 0, j-1

          myPoly % bWeights(i) = myPoly % bWeights(i)*&
                                           ( myPoly % controlPoints(i) - myPoly % controlPoints(j) )
          myPoly % bWeights(j) = myPoly % bWeights(j)*&
                                           ( myPoly % controlPoints(j) - myPoly % controlPoints(i) )

         ENDDO 
      ENDDO 
 
      DO j = 0, myPoly % N
        myPoly % bWeights(j) = 1.0_prec/myPoly % bWeights(j)
      ENDDO 

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

  SUBROUTINE CalculateInterpolationMatrix( myPoly )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly
    ! Local
    REAL(prec) :: temp1, temp2
    INTEGER    :: row, col
    LOGICAL    :: rowHasMatch 
    REAL(prec) :: iMatrix(0:myPoly % M, 0:myPoly % N)

      DO row = 0, myPoly % M

         rowHasMatch = .FALSE.
       
         DO col = 0, myPoly % N

            iMatrix(row,col) = 0.0_prec
           
            IF( AlmostEqual( myPoly % targetPoints(row), myPoly % controlPoints(col) ) )THEN
               rowHasMatch = .TRUE.
               iMatrix(row,col) = 1.0_prec
            ENDIF

         ENDDO 

         IF( .NOT.(rowHasMatch) )THEN 

            temp1 = 0.0_prec

            DO col = 0, myPoly % N        
               temp2 = myPoly % bWeights(col)/( myPoly % targetPoints(row) - myPoly % controlPoints(col) )
               iMatrix(row,col) = temp2
               temp1 = temp1 + temp2
            ENDDO 

            DO col = 0, myPoly % N 
               iMatrix(row,col) = iMatrix(row,col)/temp1
            ENDDO

         ENDIF 

      ENDDO

      myPoly % iMatrix = TRANSPOSE(iMatrix)

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

  SUBROUTINE CalculateDerivativeMatrix( myPoly )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly
    ! Local
    INTEGER    :: row, col
    REAL(prec) :: dmat(0:myPoly % N, 0:myPoly % N)

      DO row = 0, myPoly % N
         
        dmat(row,row) = 0.0_prec

        DO col = 0, myPoly % N
           
          IF( .NOT. (col == row) )THEN

            dmat(row,col) = myPoly % bWeights(col)/&
                                                 ( myPoly % bWeights(row)*&
                                                   ( myPoly % controlPoints(row) - &
                                                     myPoly % controlPoints(col) ) )

            dmat(row,row) = dmat(row,row) - dmat(row,col)

          ENDIF
        
        ENDDO 

      ENDDO 
      
      myPoly % dMatrix = TRANSPOSE( dmat )

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

  FUNCTION CalculateLagrangePolynomials( myPoly, sE ) RESULT( lAtS )  
    IMPLICIT NONE
    CLASS(Lagrange) :: myPoly
    REAL(prec)      :: sE
    REAL(prec)      :: lAtS(0:myPoly % N)
    ! Local
    REAL(prec) :: temp1, temp2
    INTEGER    :: j
    LOGICAL    :: xMatchesNode

      xMatchesNode = .FALSE.

      DO j = 0, myPoly % N
        
         lAtS(j) = 0.0_prec

         IF( AlmostEqual(sE, myPoly % controlPoints(j)) ) THEN
            lAtS(j) = 1.0_prec
            xMatchesNode = .TRUE.
         ENDIF 

      ENDDO

      IF( xMatchesNode )THEN 
         RETURN
      ENDIF

      temp1 = 0.0_prec
     
      DO j = 0, myPoly % N 
         temp2 = myPoly % bWeights(j)/(sE - myPoly % controlPoints(j))
         lAtS(j) = temp2
         temp1 = temp1 + temp2
      ENDDO 
  
      lAtS = lAtS/temp1 

  END FUNCTION CalculateLagrangePolynomials

END MODULE Lagrange_Class
