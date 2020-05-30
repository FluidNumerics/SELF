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

USE hip_fortran
USE iso_c_binding

IMPLICIT NONE


! =============================================================================================== !
!
! A data-structure for handling Lagrange interpolation in one, two, or three dimensions
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

    INTEGER :: N     
    INTEGER :: M 
    REAL(prec), POINTER :: controlPoints(:)
    REAL(prec), POINTER :: targetPoints(:)
    REAL(prec), POINTER :: bWeights(:)
    REAL(prec), POINTER :: qWeights(:)
    REAL(prec), POINTER :: iMatrix(:,:)
    REAL(prec), POINTER :: iMatrixTranspose(:,:)
    REAL(prec), POINTER :: dMatrix(:,:)  
    REAL(prec), POINTER :: dMatrixTranspose(:,:)  

    TYPE(c_ptr) :: bWeights_dev
    TYPE(c_ptr) :: qWeights_dev
    TYPE(c_ptr) :: iMatrix_dev
    TYPE(c_ptr) :: iMatrixTranspose_dev
    TYPE(c_ptr) :: dMatrix_dev
    TYPE(c_ptr) :: dMatrixTranspose_dev 

    CONTAINS
      
      PROCEDURE :: Build => Build_Lagrange
      PROCEDURE :: Trash => Trash_Lagrange

      PROCEDURE :: UpdateDevice => UpdateDevice_Lagrange
      PROCEDURE :: UpdateHost => UpdateHost_Lagrange
      
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

      GENERIC, PUBLIC :: Derivative_1D => Derivative_1D_cpu, Derivative_1D_gpu
      PROCEDURE, PRIVATE :: Derivative_1D_cpu, Derivative_1D_gpu

!      Takes in Scalar and Returns Vector
      GENERIC, PUBLIC :: ScalarGradient_2D => ScalarGradient_2D_cpu, ScalarGradient_2D_gpu
      PROCEDURE, PRIVATE :: ScalarGradient_2D_cpu, ScalarGradient_2D_gpu

!      Takes in Vector and Returns Tensor
      GENERIC, PUBLIC :: VectorGradient_2D => VectorGradient_2D_cpu, VectorGradient_2D_gpu
      PROCEDURE, PRIVATE :: VectorGradient_2D_cpu, VectorGradient_2D_gpu

!      Takes in Vector and Returns Scalar
      GENERIC, PUBLIC :: VectorDivergence_2D => VectorDivergence_2D_cpu, VectorDivergence_2D_gpu
      PROCEDURE, PRIVATE :: VectorDivergence_2D_cpu, VectorDivergence_2D_gpu

!      Takes in Vector and Returns Vector
      GENERIC, PUBLIC :: VectorCurl_2D => VectorCurl_2D_cpu, VectorCurl_2D_gpu
      PROCEDURE, PRIVATE :: VectorCurl_2D_cpu, VectorCurl_2D_gpu

!      Takes in Scalar and Returns Vector
      GENERIC, PUBLIC :: ScalarGradient_3D => ScalarGradient_3D_cpu, ScalarGradient_3D_gpu
      PROCEDURE, PRIVATE :: ScalarGradient_3D_cpu, ScalarGradient_3D_gpu

!      Takes in Vector and Returns Tensor
      GENERIC, PUBLIC :: VectorGradient_3D => VectorGradient_3D_cpu, VectorGradient_3D_gpu
      PROCEDURE, PRIVATE :: VectorGradient_3D_cpu, VectorGradient_3D_gpu

!      Takes in Vector and Returns Scalar
      GENERIC, PUBLIC :: VectorDivergence_3D => VectorDivergence_3D_cpu, VectorDivergence_3D_gpu
      PROCEDURE, PRIVATE :: VectorDivergence_3D_cpu, VectorDivergence_3D_gpu

!      Takes in Vector and Returns Vector
      GENERIC, PUBLIC :: VectorCurl_3D => VectorCurl_3D_cpu, VectorCurl_3D_gpu
      PROCEDURE, PRIVATE :: VectorCurl_3D_cpu, VectorCurl_3D_gpu

      PROCEDURE, PRIVATE :: CalculateBarycentricWeights  => CalculateBarycentricWeights_Lagrange
      PROCEDURE, PRIVATE :: CalculateInterpolationMatrix => CalculateInterpolationMatrix_Lagrange
      PROCEDURE, PRIVATE :: CalculateDerivativeMatrix    => CalculateDerivativeMatrix_Lagrange
      PROCEDURE          :: CalculateLagrangePolynomials  
      
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
                myPoly % iMatrix(0:M,0:N), &
                myPoly % iMatrixTranspose(0:N,0:M), &
                myPoly % dMatrix(0:N,0:N), &
                myPoly % dMatrixTranspose(0:N,0:N) )

      myPoly % controlPoints = 0.0_prec
      myPoly % bWeights = 0.0_prec
      myPoly % qWeights = 0.0_prec
      myPoly % targetPoints  = 0.0_prec
      myPoly % iMatrix = 0.0_prec
      myPoly % iMatrixTranspose = 0.0_prec
      myPoly % dMatrix = 0.0_prec
      myPoly % dMatrixTranspose  = 0.0_prec
      
     
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

#ifdef GPU
!      IF( N > 7 )THEN
!         CALL Logging( WARN, 'Number of control points > 7 not fully supported for 3-D SEM operations.' )
!      ENDIF
      CALL hfMalloc(myPoly % bWeights_dev, SIZEOF(myPoly % bWeights))
      CALL hfMalloc(myPoly % qWeights_dev, SIZEOF(myPoly % qWeights))
      CALL hfMalloc(myPoly % iMatrix_dev, SIZEOF(myPoly % iMatrix))
      CALL hfMalloc(myPoly % iMatrixTranspose_dev, SIZEOF(myPoly % iMatrixTranspose))
      CALL hfMalloc(myPoly % dMatrix_dev, SIZEOF(myPoly % dMatrix))
      CALL hfMalloc(myPoly % dMatrixTranspose_dev, SIZEOF(myPoly % dMatrixTranspose))

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
                  myPoly % iMatrixTranspose, &
                  myPoly % dMatrix, &
                  myPoly % dMatrixTranspose )
#ifdef GPU
      CALL hfFree(myPoly % bWeights_dev)
      CALL hfFree(myPoly % qWeights_dev)
      CALL hfFree(myPoly % iMatrix_dev)
      CALL hfFree(myPoly % iMatrixTranspose_dev)
      CALL hfFree(myPoly % dMatrix_dev)
      CALL hfFree(myPoly % dMatrixTranspose_dev)
#endif

  END SUBROUTINE Trash_Lagrange

  SUBROUTINE UpdateDevice_Lagrange(myPoly)
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly

      CALL hipFortran(hipMemcpy(myPoly % bWeights_dev, &
                                c_loc(myPoly % bWeights), &
                                SIZEOF(myPoly % bWeights), &
                                hipMemcpyHostToDevice))

      CALL hipFortran(hipMemcpy(myPoly % qWeights_dev, &
                                c_loc(myPoly % qWeights), &
                                SIZEOF(myPoly % qWeights), &
                                hipMemcpyHostToDevice))
                        
      CALL hipFortran(hipMemcpy(myPoly % iMatrix_dev, &
                                c_loc(myPoly % iMatrix), &
                                SIZEOF(myPoly % iMatrix), &
                                hipMemcpyHostToDevice))

      CALL hipFortran(hipMemcpy(myPoly % iMatrixTranspose_dev, &
                                c_loc(myPoly % iMatrixTranspose), &
                                SIZEOF(myPoly % iMatrixTranspose), &
                                hipMemcpyHostToDevice))

      CALL hipFortran(hipMemcpy(myPoly % dMatrix_dev, &
                                c_loc(myPoly % dMatrix), &
                                SIZEOF(myPoly % dMatrix), &
                                hipMemcpyHostToDevice))

      CALL hipFortran(hipMemcpy(myPoly % dMatrixTranspose_dev, &
                                c_loc(myPoly % dMatrixTranspose), &
                                SIZEOF(myPoly % dMatrixTranspose), &
                                hipMemcpyHostToDevice))

  END SUBROUTINE UpdateDevice_Lagrange

  SUBROUTINE UpdateHost_Lagrange(myPoly)
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly

      CALL hipFortran(hipMemcpy(c_loc(myPoly % bWeights), &
                                myPoly % bWeights_dev, &
                                SIZEOF(myPoly % bWeights), &
                                hipMemcpyDeviceToHost))

      CALL hipFortran(hipMemcpy(c_loc(myPoly % qWeights), &
                                myPoly % qWeights_dev, &
                                SIZEOF(myPoly % qWeights), &
                                hipMemcpyDeviceToHost))
                        
      CALL hipFortran(hipMemcpy(c_loc(myPoly % iMatrix), &
                                myPoly % iMatrix_dev, &
                                SIZEOF(myPoly % iMatrix), &
                                hipMemcpyDeviceToHost))

      CALL hipFortran(hipMemcpy(c_loc(myPoly % iMatrixTranspose), &
                                myPoly % iMatrixTranspose_dev, &
                                SIZEOF(myPoly % iMatrixTranspose), &
                                hipMemcpyDeviceToHost))

      CALL hipFortran(hipMemcpy(c_loc(myPoly % dMatrix), &
                                myPoly % dMatrix_dev, &
                                SIZEOF(myPoly % dMatrix), &
                                hipMemcpyDeviceToHost))

      CALL hipFortran(hipMemcpy(c_loc(myPoly % dMatrixTranspose), &
                                myPoly % dMatrixTranspose_dev, &
                                SIZEOF(myPoly % dMatrixTranspose), &
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
              fInterp(i,iVar,iEl) = fInterp(i,iVar,iEl) + myPoly % iMatrixTranspose(ii,i)*f(ii,iVar,iEl)
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

      CALL ScalarGridInterp_1D_gpu_wrapper(myPoly % iMatrixTranspose_dev, &
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
    REAL(prec) :: fa
   
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % M
            DO i = 0, myPoly % M
             
              fNew(i,j,iVar,iEl) = 0.0_prec

              DO jj = 0, myPoly % N
                   
                fa = 0.0_prec
                DO ii = 0, myPoly % N
                  fa = fa + f(ii,jj,iVar,iEl)*myPoly % iMatrixTranspose(ii,i)
                ENDDO
                      
                fNew(i,j,iVar,iEl) = fNew(i,j,iVar,iEl) + fa*myPoly % iMatrixTranspose(jj,j)

              ENDDO
                   
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

      CALL ScalarGridInterp_2D_gpu_wrapper(myPoly % iMatrixTranspose_dev, &
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
                  fi(1:2)= fi(1:2) + f(1:2,ii,jj,iVar,iEl)*myPoly % iMatrixTranspose(ii,i)
                ENDDO
                      
                fNew(1:2,i,j,iVar,iEl) = fNew(1:2,i,j,iVar,iEl) + fi(1:2)*myPoly % iMatrixTranspose(jj,j)

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

      CALL VectorGridInterp_2D_gpu_wrapper(myPoly % iMatrixTranspose_dev, &
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
                  fi(1:2,1:2)= fi(1:2,1:2) + f(1:2,1:2,ii,jj,iVar,iEl)*myPoly % iMatrixTranspose(ii,i)
                ENDDO
                      
                fNew(1:2,1:2,i,j,iVar,iEl) = fNew(1:2,1:2,i,j,iVar,iEl) + fi(1:2,1:2)*myPoly % iMatrixTranspose(jj,j)

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

      CALL TensorGridInterp_2D_gpu_wrapper(myPoly % iMatrixTranspose_dev, &
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
                      fi = fi + f(ii,jj,kk,iVar,iEl)*myPoly % iMatrixTranspose(ii,i)
                    ENDDO
                    fij = fij + fi*myPoly % iMatrixTranspose(jj,j)
                  ENDDO
                  fInterp(i,j,k,iVar,iEl) = fInterp(i,j,k,iVar,iEl) + fij*myPoly % iMatrixTranspose(kk,k)
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

      CALL ScalarGridInterp_3D_gpu_wrapper(myPoly % iMatrixTranspose_dev, &
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
                      fi(1:3) = fi(1:3) + f(1:3,ii,jj,kk,iVar,iEl)*myPoly % iMatrixTranspose(ii,i)
                    ENDDO
                    fij(1:3) = fij(1:3) + fi(1:3)*myPoly % iMatrixTranspose(jj,j)
                  ENDDO
                  fInterp(1:3,i,j,k,iVar,iEl) = fInterp(1:3,i,j,k,iVar,iEl) + fij(1:3)*myPoly % iMatrixTranspose(kk,k)
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

      CALL VectorGridInterp_3D_gpu_wrapper(myPoly % iMatrixTranspose_dev, &
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
                      fi(1:3,1:3) = fi(1:3,1:3) + f(1:3,1:3,ii,jj,kk,iVar,iEl)*myPoly % iMatrixTranspose(ii,i)
                    ENDDO
                    fij(1:3,1:3) = fij(1:3,1:3) + fi(1:3,1:3)*myPoly % iMatrixTranspose(jj,j)
                  ENDDO
                  fInterp(1:3,1:3,i,j,k,iVar,iEl) = fInterp(1:3,1:3,i,j,k,iVar,iEl) + fij(1:3,1:3)*myPoly % iMatrixTranspose(kk,k)
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

      CALL TensorGridInterp_3D_gpu_wrapper(myPoly % iMatrixTranspose_dev, &
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
              df(i,iVar,iEl) = df(i,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(ii,iVar,iEl)
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

      CALL Derivative_1D_gpu_wrapper(myPoly % dMatrixTranspose_dev, &
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
                gradF(1,i,j,iVar,iEl) = gradF(1,i,j,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(ii,j,iVar,iEl)
                gradF(2,i,j,iVar,iEl) = gradF(2,i,j,iVar,iEl) + myPoly % dMatrixTranspose(ii,j)*f(i,ii,iVar,iEl)
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

      CALL ScalarGradient_2D_gpu_wrapper(myPoly % dMatrixTranspose_dev, &
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
                gradF(1,1,i,j,iVar,iEl) = gradF(1,1,i,j,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(1,ii,j,iVar,iEl)
                gradF(2,1,i,j,iVar,iEl) = gradF(2,1,i,j,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(2,ii,j,iVar,iEl)
                gradF(1,2,i,j,iVar,iEl) = gradF(1,2,i,j,iVar,iEl) + myPoly % dMatrixTranspose(ii,j)*f(1,i,ii,iVar,iEl)
                gradF(2,2,i,j,iVar,iEl) = gradF(2,2,i,j,iVar,iEl) + myPoly % dMatrixTranspose(ii,j)*f(2,i,ii,iVar,iEl)
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

      CALL VectorGradient_2D_gpu_wrapper(myPoly % dMatrixTranspose_dev, &
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
                dF(i,j,iVar,iEl) = dF(i,j,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(1,ii,j,iVar,iEl) +&
                                                      myPoly % dMatrixTranspose(ii,j)*f(2,i,ii,iVar,iEl)
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

      CALL VectorDivergence_2D_gpu_wrapper(myPoly % dMatrixTranspose_dev, &
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
                dF(i,j,iVar,iEl) = dF(i,j,iVar,iEl) + myPoly % dMatrixTranspose(ii,j)*f(1,i,ii,iVar,iEl) -&
                                                      myPoly % dMatrixTranspose(ii,i)*f(2,ii,j,iVar,iEl)
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

      CALL VectorCurl_2D_gpu_wrapper(myPoly % dMatrixTranspose_dev, &
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
                  gradF(1,i,j,k,iVar,iEl) = gradF(1,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(ii,j,k,iVar,iEl)
                  gradF(2,i,j,k,iVar,iEl) = gradF(2,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,j)*f(i,ii,k,iVar,iEl)
                  gradF(3,i,j,k,iVar,iEl) = gradF(3,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,k)*f(i,j,ii,iVar,iEl)
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

      CALL ScalarGradient_3D_gpu_wrapper(myPoly % dMatrixTranspose_dev, &
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
                  gradF(1,1,i,j,k,iVar,iEl) = gradF(1,1,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(1,ii,j,k,iVar,iEl)
                  gradF(2,1,i,j,k,iVar,iEl) = gradF(2,1,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(2,ii,j,k,iVar,iEl)
                  gradF(3,1,i,j,k,iVar,iEl) = gradF(3,1,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(3,ii,j,k,iVar,iEl)
                  gradF(1,2,i,j,k,iVar,iEl) = gradF(1,2,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,j)*f(1,i,ii,k,iVar,iEl)
                  gradF(2,2,i,j,k,iVar,iEl) = gradF(2,2,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,j)*f(2,i,ii,k,iVar,iEl)
                  gradF(3,2,i,j,k,iVar,iEl) = gradF(3,2,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,j)*f(3,i,ii,k,iVar,iEl)
                  gradF(1,3,i,j,k,iVar,iEl) = gradF(1,3,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,k)*f(1,i,j,ii,iVar,iEl)
                  gradF(2,3,i,j,k,iVar,iEl) = gradF(2,3,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,k)*f(2,i,j,ii,iVar,iEl)
                  gradF(3,3,i,j,k,iVar,iEl) = gradF(3,3,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,k)*f(3,i,j,ii,iVar,iEl)
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

      CALL VectorGradient_3D_gpu_wrapper(myPoly % dMatrixTranspose_dev, &
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
                  dF(i,j,k,iVar,iEl) = dF(i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(1,ii,j,k,iVar,iEl) +&
                                                            myPoly % dMatrixTranspose(ii,j)*f(2,i,ii,k,iVar,iEl) +&
                                                            myPoly % dMatrixTranspose(ii,k)*f(3,i,j,ii,iVar,iEl)
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

      CALL VectorDivergence_3D_gpu_wrapper(myPoly % dMatrixTranspose_dev, &
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
                  dF(1,i,j,k,iVar,iEl) = dF(1,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,j)*f(3,i,ii,k,iVar,iEl) -&
                                                                myPoly % dMatrixTranspose(ii,k)*f(2,i,j,ii,iVar,iEl)
                  dF(2,i,j,k,iVar,iEl) = dF(2,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,k)*f(1,i,j,ii,iVar,iEl) -&
                                                                myPoly % dMatrixTranspose(ii,i)*f(3,ii,j,k,iVar,iEl)
                  dF(3,i,j,k,iVar,iEl) = dF(3,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(2,ii,j,k,iVar,iEl) -&
                                                                myPoly % dMatrixTranspose(ii,j)*f(1,i,ii,k,iVar,iEl)
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

      CALL VectorCurl_3D_gpu_wrapper(myPoly % dMatrixTranspose_dev, &
                                         f_dev, dF_dev, myPoly % N, &
                                         nVariables, nElements)

  END SUBROUTINE VectorCurl_3D_gpu
!! ================================================================================================ !
!!
!! CalculateDivergence_2D 
!!
!!   Calculates the gradient of a 2-D function, represented by a 2-D array of nodal values.
!!  
!!   Let \vec{F} = f \hat{x} + g \hat{y}. 
!!
!!   Given a set of nodal values at the interpolation nodes, the divergence of a vector function
!!   through the interpolation nodes can be estimated by  
!!
!!                      div( \vec{F} )_{a,b} = (df/dx)_{a,b} + (dg/dy)_{a,b}
!!   where
!!                       (df/dx)_{a,b} = \sum_{i=0}^N f_{i,b} l'_i(\xi_a),   a,b=0,1,2,...,N
!!                       (dg/dy)_{a,b} = \sum_{j=0}^N g_{a,j} l'_j(\xi_b),   a,b=0,1,2,...,N
!!             
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
!!     REAL(prec)     :: f(1:2,0:interp % N,0:interp % N,1:nVariables,1:nElements) 
!!     REAL(prec)     :: divF(0:interp % N,0:interp % N,1:nVariables,1:nElements) 
!!     
!!       CALL interp % CalculateDivergence_2D( f, divF, nVariables, nElements ) 
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
!!     divF (out) 
!!      Array of derivative values at the target interpolation nodes.
!!   
!! ================================================================================================ ! 
!
!  SUBROUTINE CalculateDivergence_2D( myPoly, f, divF, nVariables, nElements )
!    IMPLICIT NONE
!    CLASS(Lagrange), INTENT(in) :: myPoly
!
!#ifdef GPU
!    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
!    REAL(prec), DEVICE, INTENT(in)  :: f(1:2, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), DEVICE, INTENT(out) :: divF(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    TYPE(dim3) :: grid, tBlock
!    INTEGER    :: threadCount
!  
!      threadCount = MIN( 4*(ceiling( REAL(myPoly % N+1)/4 ) ), 8 )
!
!      tBlock = dim3( threadCount, threadCount, 1 )
!      grid   = dim3( nVariables, nElements, 1)
!
!      CALL CalculateDivergence_2D_CUDAKernel<<<grid, tBlock>>>( myPoly % dMatrixTranspose_dev, &
!                                                                f, divF, myPoly % N_dev, nVariables, nElements )
!
!#else
!    INTEGER, INTENT(in)     :: nVariables, nElements
!    REAL(prec), INTENT(in)  :: f(1:2, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), INTENT(out) :: divF(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!
!      divF = CalculateDivergence_2D_Lagrange( myPoly, f, nVariables, nElements )
!
!#endif
!
!  END SUBROUTINE CalculateDivergence_2D
!
!! ================================================================================================ !
!!
!! CalculateGradient_3D 
!!
!!   Calculates the gradient of a 3-D function, represented by a 3-D array of nodal values.
!!   
!!   Given a set of nodal values at the interpolation nodes, the gradient of a function through 
!!   the interpolation nodes can be estimated by  
!!
!!               (df/dx)_{a,b,c} = \sum_{i=0}^N f_{i,b,c} l'_i(\xi_a),   a,b,c=0,1,2,...,N
!!               (df/dy)_{a,b,c} = \sum_{j=0}^N f_{a,j,c} l'_j(\xi_b),   a,b,c=0,1,2,...,N
!!               (df/dz)_{a,b,c} = \sum_{k=0}^N f_{a,b,k} l'_k(\xi_c),   a,b,c=0,1,2,...,N
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
!!     REAL(prec)     :: f(0:interp % N,0:interp % N,0:interp % N,1:nVariables,1:nElements) 
!!     REAL(prec)     :: gradF(1:3,0:interp % N,0:interp % N,0:interp % N,1:nVariables,1:nElements) 
!!     
!!       CALL interp % CalculateGradient_3D( f, gradF, nVariables, nElements ) 
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
!  SUBROUTINE CalculateGradient_3D( myPoly, f, gradF, nVariables, nElements )
!    IMPLICIT NONE
!    CLASS(Lagrange), INTENT(in) :: myPoly
!
!#ifdef GPU
!    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
!    REAL(prec), DEVICE, INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), DEVICE, INTENT(out) :: gradF(1:3,0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    TYPE(dim3) :: grid, tBlock
!    INTEGER    :: threadCount
!  
!      threadCount = MIN( 4*(ceiling( REAL(myPoly % N+1)/4 ) ), 8 )
!
!      tBlock = dim3( threadCount, threadCount, threadCount )
!      grid   = dim3( nVariables, nElements, 1)
!
!      CALL CalculateGradient_3D_CUDAKernel<<<grid, tBlock>>>( myPoly % dMatrixTranspose_dev, &
!                                                              f, gradF, myPoly % N_dev, nVariables, nElements )
!
!#else
!    INTEGER, INTENT(in)     :: nVariables, nElements
!    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), INTENT(out) :: gradF(1:3, 0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!
!      gradF = CalculateGradient_3D_Lagrange( myPoly, f, nVariables, nElements )
!
!#endif
!
!  END SUBROUTINE CalculateGradient_3D
!
!! ================================================================================================ !
!!
!! CalculateDivergence_3D 
!!
!!   Calculates the gradient of a 3-D function, represented by a 3-D array of nodal values.
!!  
!!   Let \vec{F} = f \hat{x} + g \hat{y} + h \hat{z}. 
!!
!!   Given a set of nodal values at the interpolation nodes, the divergence of a vector function
!!   through the interpolation nodes can be estimated by  
!!
!!            div( \vec{F} )_{a,b,c} = (df/dx)_{a,b,c} + (dg/dy)_{a,b,c} + (dh/dy)_{a,b,c}
!!   where
!!             (df/dx)_{a,b,c} = \sum_{i=0}^N f_{i,b,c} l'_i(\xi_a),   a,b,c=0,1,2,...,N
!!             (dg/dy)_{a,b,c} = \sum_{j=0}^N g_{a,j,c} l'_j(\xi_b),   a,b,c=0,1,2,...,N
!!             (dg/dy)_{a,b,c} = \sum_{k=0}^N g_{a,b,k} l'_k(\xi_c),   a,b,c=0,1,2,...,N
!!             
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
!!     REAL(prec)     :: f(1:3,0:interp % N,0:interp % N,0:interp % N,1:nVariables,1:nElements) 
!!     REAL(prec)     :: divF(0:interp % N,0:interp % N,0:interp %N,1:nVariables,1:nElements) 
!!     
!!       CALL interp % CalculateDivergence_3D( f, divF, nVariables, nElements ) 
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
!!     divF (out) 
!!      Array of derivative values at the target interpolation nodes.
!!   
!! ================================================================================================ ! 
!
!  SUBROUTINE CalculateDivergence_3D( myPoly, f, divF, nVariables, nElements )
!    IMPLICIT NONE
!    CLASS(Lagrange), INTENT(in) :: myPoly
!
!#ifdef GPU
!    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
!    REAL(prec), DEVICE, INTENT(in)  :: f(1:3, 0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), DEVICE, INTENT(out) :: divF(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    TYPE(dim3) :: grid, tBlock
!    INTEGER    :: threadCount
!  
!      threadCount = MIN( 4*(ceiling( REAL(myPoly % N+1)/4 ) ), 8 )
!
!      tBlock = dim3( threadCount, threadCount, threadCount )
!      grid   = dim3( nVariables, nElements, 1)
!
!      CALL CalculateDivergence_3D_CUDAKernel<<<grid, tBlock>>>( myPoly % dMatrixTranspose_dev, &
!                                                                f, divF, myPoly % N_dev, nVariables, nElements )
!
!#else
!    INTEGER, INTENT(in)     :: nVariables, nElements
!    REAL(prec), INTENT(in)  :: f(1:3, 0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), INTENT(out) :: divF(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!
!      divF = CalculateDivergence_3D_Lagrange( myPoly, f, nVariables, nElements )
!
!#endif
!
!  END SUBROUTINE CalculateDivergence_3D
! 
! ================================================================================================ !
! ------------------------------------- PRIVATE ROUTINES ----------------------------------------- !
! ================================================================================================ !


! ================================================================================================ !
!
! CalculateBarycentricWeights_Lagrange (PRIVATE)
!
!   A PRIVATE routine that calculates and stores the barycentric weights for the Lagrange 
!   data-structure.
! 
!   This routine is from Alg. 30 on pg. 75 of D.A. Kopriva, 2009.
! 
! ================================================================================================ ! 

  SUBROUTINE CalculateBarycentricWeights_Lagrange( myPoly )
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

  END SUBROUTINE CalculateBarycentricWeights_Lagrange

! ================================================================================================ !
!
! CalculateInterpolationMatrix_Lagrange (PRIVATE) 
!
!   A PRIVATE routine that fills in the interpolation matrix for the Lagrange data structure.
!
!   This function is from Alg. 32 on pg. 76 of D.A. Kopriva, 2009.
!
! ================================================================================================ ! 

  SUBROUTINE CalculateInterpolationMatrix_Lagrange( myPoly )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly
    ! Local
    REAL(prec) :: temp1, temp2
    INTEGER    :: row, col
    LOGICAL    :: rowHasMatch 

      DO row = 0, myPoly % M

         rowHasMatch = .FALSE.
       
         DO col = 0, myPoly % N

            myPoly % iMatrix(row,col) = 0.0_prec
           
            IF( AlmostEqual( myPoly % targetPoints(row), myPoly % controlPoints(col) ) )THEN
               rowHasMatch = .TRUE.
               myPoly % iMatrix(row,col) = 1.0_prec
            ENDIF

         ENDDO 

         IF( .NOT.(rowHasMatch) )THEN 

            temp1 = 0.0_prec

            DO col = 0, myPoly % N        
               temp2 = myPoly % bWeights(col)/( myPoly % targetPoints(row) - myPoly % controlPoints(col) )
               myPoly % iMatrix(row,col) = temp2
               temp1 = temp1 + temp2
            ENDDO 

            DO col = 0, myPoly % N 
               myPoly % iMatrix(row,col) = myPoly % iMatrix(row,col)/temp1
            ENDDO

         ENDIF 

      ENDDO

      myPoly % iMatrixTranspose = TRANSPOSE( myPoly % iMatrix )

 END SUBROUTINE CalculateInterpolationMatrix_Lagrange

! ================================================================================================ !
!
! CalculateDerivativeMatrix_Lagrange (PRIVATE) 
!
!   Calculates and stores the derivative matrix and its transpose. 
!   Generates a matrix that can be used to approximate derivatives at the interpolation nodes.
!
!   This function is from Alg. 37 on pg. 82 of D.A. Kopriva, 2009.
!
! ================================================================================================ ! 

  SUBROUTINE CalculateDerivativeMatrix_Lagrange( myPoly )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly
    ! Local
    INTEGER    :: row, col

      DO row = 0, myPoly % N
         
        myPoly % dMatrix(row,row) = 0.0_prec

        DO col = 0, myPoly % N
           
          IF( .NOT. (col == row) )THEN

            myPoly % dMatrix(row,col) = myPoly % bWeights(col)/&
                                                 ( myPoly % bWeights(row)*&
                                                   ( myPoly % controlPoints(row) - &
                                                     myPoly % controlPoints(col) ) )

            myPoly % dMatrix(row,row) = myPoly % dMatrix(row,row) - myPoly % dMatrix(row,col)

          ENDIF
        
        ENDDO 

      ENDDO 
      
      myPoly % dMatrixTranspose = TRANSPOSE( myPoly % dMatrix )

  END SUBROUTINE CalculateDerivativeMatrix_Lagrange

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
!
!  FUNCTION GridInterpolation_1D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( fNew )  
!    IMPLICIT NONE
!    TYPE(Lagrange) :: myPoly
!    INTEGER         :: nVariables, nElements
!    REAL(prec)      :: f(0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec)      :: fNew(0:myPoly % M, 1:nVariables, 1:nElements)
!
!  END FUNCTION GridInterpolation_1D_Lagrange
!
!  FUNCTION GridInterpolation_2D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( fNew )  
!    IMPLICIT NONE
!    TYPE(Lagrange) :: myPoly
!    INTEGER         :: nVariables, nElements
!    REAL(prec)      :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec)      :: fNew(0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
!
!  END FUNCTION GridInterpolation_2D_Lagrange
!
!  FUNCTION GridInterpolation_3D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( fNew )  
!    IMPLICIT NONE
!    TYPE(Lagrange) :: myPoly
!    INTEGER         :: nElements, nVariables
!    REAL(prec)      :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec)      :: fNew(0:myPoly % M, 0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
!
!  END FUNCTION GridInterpolation_3D_Lagrange
!
!  FUNCTION Derivative_1D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( derF )  
!    IMPLICIT NONE
!    TYPE(Lagrange) :: myPoly
!    INTEGER         :: nVariables, nElements
!    REAL(prec)      :: f(0:myPoly % N,1:nVariables,1:nElements)
!    REAL(prec)      :: derF(0:myPoly % N,1:nVariables,1:nElements)
!    ! Local
!    INTEGER :: i, ii, iVar, iEl
!
!      DO iEl = 1, nElements
!        DO iVar = 1, nVariables
! 
!          DO i = 0, myPoly % N
!    
!            derF(i,iVar,iEl) = 0.0_prec 
!            DO ii = 0, myPoly % N
!              derF(i,iVar,iEl) = derF(i,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(ii,iVar,iEl)
!            ENDDO
!    
!          ENDDO
!
!        ENDDO
!      ENDDO
!
!  END FUNCTION Derivative_1D_Lagrange
!!
!  FUNCTION CalculateGradient_2D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( gradF )  
!
!    IMPLICIT NONE
!    TYPE(Lagrange) :: myPoly
!    INTEGER         :: nVariables, nElements
!    REAL(prec)      :: f(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
!    REAL(prec)      :: gradF(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
!
!  END FUNCTION CalculateGradient_2D_Lagrange
!!
!  FUNCTION CalculateDivergence_2D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( divF )  
!    IMPLICIT NONE
!    TYPE( Lagrange ) :: myPoly
!    INTEGER          :: nVariables, nElements
!    REAL(prec)       :: f(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
!    REAL(prec)       :: divF(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
!    ! Local
!    INTEGER :: i, j, iEl, iVar, ii
!  
!      
!      DO iEl = 1, nElements
!        DO iVar = 1, nVariables
!          DO j = 0, myPoly % N
!            DO i = 0, myPoly % N
!               
!              divf(i,j,iVar,iEl) = 0.0_prec
!
!              DO ii = 0, myPoly % N  
!                divf(i,j,iVar,iEl) = divf(i,j,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(1,ii,j,iVar,iEl) + &
!                                                            myPoly % dMatrixTranspose(ii,j)*f(2,i,ii,iVar,iEl)
!              ENDDO
!
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO
!      
!  END FUNCTION CalculateDivergence_2D_Lagrange
!!
!  FUNCTION CalculateGradient_3D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( gradF )  
!    IMPLICIT NONE
!    TYPE( Lagrange ) :: myPoly
!    INTEGER          :: nVariables, nElements
!    REAL(prec)       :: f(0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
!    REAL(prec)       :: gradF(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
!    ! Local
!    INTEGER :: i, j, k, iEl, iVar, ii
!  
!      
!      DO iEl = 1, nElements
!        DO iVar = 1, nVariables
!          DO k = 0, myPoly % N
!            DO j = 0, myPoly % N
!              DO i = 0, myPoly % N
!               
!                gradf(1,i,j,k,iVar,iEl) = 0.0_prec
!                gradf(2,i,j,k,iVar,iEl) = 0.0_prec
!                gradf(3,i,j,k,iVar,iEl) = 0.0_prec
!
!                DO ii = 0, myPoly % N  
!                  gradf(1,i,j,k,iVar,iEl) = gradf(1,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(ii,j,k,iVar,iEl)
!                  gradf(2,i,j,k,iVar,iEl) = gradf(2,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,j)*f(i,ii,k,iVar,iEl)
!                  gradf(3,i,j,k,iVar,iEl) = gradf(3,i,j,k,iVar,iEl) + myPoly % dMatrixTranspose(ii,k)*f(i,j,ii,iVar,iEl)
!                ENDDO
!
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO
!      
!  END FUNCTION CalculateGradient_3D_Lagrange
!!
!  FUNCTION CalculateDivergence_3D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( divF )  
!    IMPLICIT NONE
!    TYPE( Lagrange ) :: myPoly
!    INTEGER          :: nVariables, nElements
!    REAL(prec)       :: f(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
!    REAL(prec)       :: divF(0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
!    ! Local
!    INTEGER :: i, j, k, iEl, iVar, ii
!  
!      DO iEl = 1, nElements
!        DO iVar = 1, nVariables
!          DO k = 0, myPoly % N
!            DO j = 0, myPoly % N
!              DO i = 0, myPoly % N
!               
!                divf(i,j,k,iVar,iEl) = 0.0_prec
!
!                DO ii = 0, myPoly % N  
!                  divf(i,j,k,iVar,iEl) = divf(i,j,k,iVar,iEl) + &
!                                         myPoly % dMatrixTranspose(ii,i)*f(1,ii,j,k,iVar,iEl) + &
!                                         myPoly % dMatrixTranspose(ii,j)*f(2,i,ii,k,iVar,iEl) + &
!                                         myPoly % dMatrixTranspose(ii,k)*f(3,i,j,ii,iVar,iEl)
!                ENDDO
!
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO
!      
!  END FUNCTION CalculateDivergence_3D_Lagrange
!
!!#ifdef GPU
!!
!!  ATTRIBUTES(Global) SUBROUTINE GridInterpolation_1D_CUDAKernel( IntMatT, f, fNew, N, M, nVariables, nElems  )
!!    ! Currently only works for M,N <= 7
!!    IMPLICIT NONE
!!    INTEGER, INTENT(in)             :: N, M, nVariables, nElems
!!    REAL(prec), DEVICE, INTENT(in)  :: IntMatT(0:N,0:M)
!!    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,1:nVariables,1:nElems)
!!    REAL(prec), DEVICE, INTENT(out) :: fnew(0:M,1:nVariables,1:nElems)
!!    ! Local
!!    INTEGER            :: i, a, iEl, iVar
!!    REAL(prec), SHARED :: floc(0:7)
!!    REAL(prec)         :: fm
!!   
!!      iVar = blockIdx % x
!!      iEl  = blockIdx % y
!!      a = threadIdx % x-1 + (blockIDx % y-1)*blockDim % x
!!      
!!      IF( a <= M )THEN
!!
!!        ! Pre-fetch data
!!        IF( a <= N )THEN
!!          floc(a) = f(a,iVar,iEl)
!!        ENDIF
!!        
!!        CALL syncthreads( )
!!      
!!        fm = 0.0_prec
!!        DO i = 0, N
!!          fm = fm + floc(i)*IntMatT(i,a)
!!        ENDDO
!!               
!!        fnew(a,iVar,iEl) = fm
!!         
!!      ENDIF
!!
!!  END SUBROUTINE GridInterpolation_1D_CUDAKernel
!!!
!!  ATTRIBUTES(Global) SUBROUTINE GridInterpolation_2D_CUDAKernel( IntMatT, f, fNew, N, M, nVariables, nElems  )
!!    ! Currently only works for M,N <= 7
!!    IMPLICIT NONE
!!    INTEGER, INTENT(in)             :: N, M, nVariables, nElems
!!    REAL(prec), DEVICE, INTENT(in)  :: IntMatT(0:N,0:M)
!!    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,0:N,1:nVariables,1:nElems)
!!    REAL(prec), DEVICE, INTENT(out) :: fnew(0:M,0:M,1:nVariables,1:nElems)
!!    ! Local
!!    INTEGER            :: i, j, a, b, iEl, iVar
!!    REAL(prec), SHARED :: floc(0:7,0:7)
!!    REAL(prec)         :: fm, fmn
!!   
!!      iVar = blockIdx % x
!!      iEl  = blockIdx % y
!!      a = threadIdx % x-1 + (blockIDx % y-1)*blockDim % x
!!      b = threadIdx % y-1 + (blockIDx % y-1)*blockDim % y
!!      
!!      IF( a <= M .AND. b <= M )THEN
!!
!!        ! Pre-fetch data
!!        IF( a <= N .AND. b <= N )THEN
!!          floc(a,b) = f(a,b,iVar,iEl)
!!        ENDIF
!!        
!!        CALL syncthreads( )
!!               
!!        fmn = 0.0_prec
!!        DO j = 0, N
!!            
!!          fm = 0.0_prec
!!          DO i = 0, N
!!            fm = fm + floc(i,j)*IntMatT(i,a)
!!          ENDDO
!!               
!!          fmn = fmn + fm*IntMatT(j,b)
!!
!!        ENDDO
!!            
!!        fnew(a,b,iVar,iEl) = fmn
!!         
!!      ENDIF
!!
!!  END SUBROUTINE GridInterpolation_2D_CUDAKernel
!!
!!  ATTRIBUTES(Global) SUBROUTINE GridInterpolation_3D_CUDAKernel( IntMatT, f, fNew, N, M, nVariables, nElems  )
!!    ! Currently only works for M,N <= 7
!!    IMPLICIT NONE
!!    INTEGER, INTENT(in)             :: N, M, nVariables, nElems
!!    REAL(prec), DEVICE, INTENT(in)  :: IntMatT(0:N,0:M)
!!    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,0:N,0:N,1:nVariables,1:nElems)
!!    REAL(prec), DEVICE, INTENT(out) :: fnew(0:M,0:M,0:M,1:nVariables,1:nElems)
!!    ! Local
!!    INTEGER            :: i, j, k, a, b, c, iEl, iVar
!!    REAL(prec), SHARED :: floc(0:7,0:7,0:7)
!!    REAL(prec)         :: fm, fmn, fmnp
!!   
!!   
!!      iVar = blockIdx % x
!!      iEl  = blockIdx % y
!!      a = threadIdx % x-1 + (blockIDx % y-1)*blockDim % x
!!      b = threadIdx % y-1 + (blockIDx % y-1)*blockDim % y
!!      c = threadIdx % z-1 + (blockIDx % y-1)*blockDim % z
!!      
!!      IF( a <= M .AND. b <= M .AND. c <= M )THEN
!!
!!        ! Pre-fetch data
!!        IF( a <= N .AND. b <= N .AND. c <= N )THEN
!!          floc(a,b,c) = f(a,b,c,iVar,iEl)
!!        ENDIF
!!        
!!        CALL syncthreads( )
!!      
!!              
!!        fmnp = 0.0_prec
!!        DO k = 0, N
!!         
!!          fmn = 0.0_prec
!!          DO j = 0, N
!!          
!!            fm = 0.0_prec
!!            DO i = 0, N
!!              fm = fm + floc(i,j,k)*IntMatT(i,a)
!!            ENDDO
!!             
!!            fmn = fmn + fm*IntMatT(j,b)
!!          ENDDO
!!           
!!          fmnp = fmnp + fmn*IntMatT(k,c)
!!
!!        ENDDO
!!         
!!         fnew(a,b,c,iVar,iEl) = fmnp
!!         
!!      ENDIF
!!
!!  END SUBROUTINE GridInterpolation_3D_CUDAKernel
!!!
!!  ATTRIBUTES(Global) SUBROUTINE Derivative_1D_CUDAKernel( DMatT, f, derf, N, nVariables, nElems )
!!    ! Currently only works for M,N <= 7
!!    IMPLICIT NONE
!!    INTEGER, INTENT(in)             :: N, nVariables, nElems
!!    REAL(prec), DEVICE, INTENT(in)  :: DMatT(0:N,0:N)
!!    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,1:nVariables,1:nElems)
!!    REAL(prec), DEVICE, INTENT(out) :: derf(0:N,1:nVariables,1:nElems)
!!    ! Local
!!    INTEGER            :: i, j, iEl, iVar, ii
!!    REAL(prec), SHARED :: floc(0:7)
!!    REAL(prec)         :: df
!!  
!!  
!!      iVar = blockIdx % x
!!      iEl  = blockIdx % y
!!      i    = threadIdx % x-1
!!      
!!      IF( i <= N )THEN
!!      
!!         ! Pre-fetch into shared memory for this block
!!         floc(i) = f(i,iVar,iEl)
!!
!!         CALL syncthreads( )
!!
!!         df = 0.0_prec
!!
!!         DO ii = 0, N
!!            df = df + DMatT(ii,i)*floc(ii)
!!         ENDDO
!!                  
!!         derf(i,iVar,iEl) = df
!!
!!      ENDIF
!!    
!!  END SUBROUTINE Derivative_1D_CUDAKernel
!!!
!!  ATTRIBUTES(Global) SUBROUTINE CalculateGradient_2D_CUDAKernel( DMatT, f, gradf, N, nVariables, nElems )
!!    ! Currently only works for M,N <= 7
!!    IMPLICIT NONE
!!    INTEGER, INTENT(in)             :: N, nVariables, nElems
!!    REAL(prec), DEVICE, INTENT(in)  :: DMatT(0:N,0:N)
!!    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,0:N,1:nVariables,1:nElems)
!!    REAL(prec), DEVICE, INTENT(out) :: gradf(1:2,0:N,0:N,1:nVariables,1:nElems)
!!    ! Local
!!    INTEGER            :: i, j, iEl, iVar, ii
!!    REAL(prec), SHARED :: floc(0:7,0:7)
!!    REAL(prec)         :: df(1:2)
!!  
!!  
!!      iVar = blockIdx % x
!!      iEl  = blockIdx % y
!!      i    = threadIdx % x-1
!!      j    = threadIdx % y-1
!!      
!!      IF( i <= N .AND. j <= N )THEN
!!      
!!         ! Pre-fetch into shared memory for this block
!!         floc(i,j) = f(i,j,iVar,iEl)
!!
!!         CALL syncthreads( )
!!
!!         df(1) = 0.0_prec
!!         df(2) = 0.0_prec
!!
!!         DO ii = 0, N
!!            df(1) = df(1) + DMatT(ii,i)*floc(ii,j)
!!            df(2) = df(2) + DMatT(ii,j)*floc(i,ii)
!!         ENDDO
!!                  
!!         gradf(1,i,j,iVar,iEl) = df(1)
!!         gradf(2,i,j,iVar,iEl) = df(2)
!!
!!      ENDIF
!!    
!!  END SUBROUTINE CalculateGradient_2D_CUDAKernel
!!!
!!  ATTRIBUTES(Global) SUBROUTINE CalculateDivergence_2D_CUDAKernel( DMatT, f, divf, N, nVariables, nElems )
!!    ! Currently only works for M,N <= 7
!!    IMPLICIT NONE
!!    INTEGER, INTENT(in)             :: N, nVariables, nElems
!!    REAL(prec), DEVICE, INTENT(in)  :: DMatT(0:N,0:N)
!!    REAL(prec), DEVICE, INTENT(in)  :: f(1:2,0:N,0:N,1:nVariables,1:nElems)
!!    REAL(prec), DEVICE, INTENT(out) :: divf(0:N,0:N,1:nVariables,1:nElems)
!!    ! Local
!!    INTEGER            :: i, j, k, iEl, iVar, ii
!!    REAL(prec), SHARED :: floc(1:2,0:7,0:7)
!!    REAL(prec) :: df
!!  
!!  
!!      iVar = blockIdx % x
!!      iEl  = blockIdx % y
!!      i    = threadIdx % x-1
!!      j    = threadIdx % y-1
!!      
!!      IF( i <= N .AND. j <= N )THEN
!!      
!!         ! Pre-fetch into shared memory for this block
!!         floc(1,i,j) = f(1,i,j,iVar,iEl)
!!         floc(2,i,j) = f(2,i,j,iVar,iEl)
!!
!!         CALL syncthreads( )
!!
!!         df = 0.0_prec
!!         DO ii = 0, N
!!            df = df + DMatT(ii,i)*floc(1,ii,j) + &
!!                      DMatT(ii,j)*floc(2,i,ii)
!!         ENDDO
!!                  
!!         divf(i,j,iVar,iEl) = df
!!
!!      ENDIF
!!    
!!  END SUBROUTINE CalculateDivergence_2D_CUDAKernel 
!!!
!!  ATTRIBUTES(Global) SUBROUTINE CalculateGradient_3D_CUDAKernel( DMatT, f, gradf, N, nVariables, nElems )
!!    ! Currently only works for M,N <= 7
!!    IMPLICIT NONE
!!    INTEGER, INTENT(in)             :: N, nVariables, nElems
!!    REAL(prec), DEVICE, INTENT(in)  :: DMatT(0:N,0:N)
!!    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,0:N,0:N,1:nVariables,1:nElems)
!!    REAL(prec), DEVICE, INTENT(out) :: gradf(1:3,0:N,0:N,0:N,1:nVariables,1:nElems)
!!    ! Local
!!    INTEGER            :: i, j, k, iEl, iVar, ii
!!    REAL(prec), SHARED :: floc(0:7,0:7,0:7)
!!    REAL(prec)         :: df(1:3)
!!  
!!  
!!      iVar = blockIdx % x
!!      iEl  = blockIdx % y
!!      i    = threadIdx % x-1
!!      j    = threadIdx % y-1
!!      k    = threadIdx % z-1
!!      
!!      IF( i <= N .AND. j <= N .AND. k <= N )THEN
!!      
!!         ! Pre-fetch into shared memory for this block
!!         floc(i,j,k) = f(i,j,k,iVar,iEl)
!!
!!         CALL syncthreads( )
!!
!!         df(1) = 0.0_prec
!!         df(2) = 0.0_prec
!!         df(3) = 0.0_prec
!!
!!         DO ii = 0, N
!!            df(1) = df(1) + DMatT(ii,i)*floc(ii,j,k)
!!            df(2) = df(2) + DMatT(ii,j)*floc(i,ii,k)
!!            df(3) = df(3) + DMatT(ii,k)*floc(i,j,ii)
!!         ENDDO
!!                  
!!         gradf(1,i,j,k,iVar,iEl) = df(1)
!!         gradf(2,i,j,k,iVar,iEl) = df(2)
!!         gradf(3,i,j,k,iVar,iEl) = df(3)
!!
!!      ENDIF
!!    
!!  END SUBROUTINE CalculateGradient_3D_CUDAKernel
!!!
!!  ATTRIBUTES(Global) SUBROUTINE CalculateDivergence_3D_CUDAKernel( DMatT, f, divf, N, nVariables, nElems )
!!    ! Currently only works for M,N <= 7
!!    IMPLICIT NONE
!!    INTEGER, INTENT(in)             :: N, nVariables, nElems
!!    REAL(prec), DEVICE, INTENT(in)  :: DMatT(0:N,0:N)
!!    REAL(prec), DEVICE, INTENT(in)  :: f(1:3,0:N,0:N,0:N,1:nVariables,1:nElems)
!!    REAL(prec), DEVICE, INTENT(out) :: divf(0:N,0:N,0:N,1:nVariables,1:nElems)
!!    ! Local
!!    INTEGER            :: i, j, k, iEl, iVar, ii
!!    REAL(prec), SHARED :: floc(1:3,0:7,0:7,0:7)
!!    REAL(prec) :: df
!!  
!!  
!!      iVar = blockIdx % x
!!      iEl  = blockIdx % y
!!      i    = threadIdx % x-1
!!      j    = threadIdx % y-1
!!      k    = threadIdx % z-1
!!      
!!      IF( i <= N .AND. j <= N .AND. k <= N )THEN
!!      
!!         ! Pre-fetch into shared memory for this block
!!         floc(1,i,j,k) = f(1,i,j,k,iVar,iEl)
!!         floc(2,i,j,k) = f(2,i,j,k,iVar,iEl)
!!         floc(3,i,j,k) = f(3,i,j,k,iVar,iEl)
!!      ENDIF
!!      
!!         CALL syncthreads( )
!! 
!!      IF( i <= N .AND. j <= N .AND. k <= N )THEN
!!      
!!         df = 0.0_prec
!!         DO ii = 0, N
!!            df = df + DMatT(ii,i)*floc(1,ii,j,k) + &
!!                      DMatT(ii,j)*floc(2,i,ii,k) + &
!!                      DMatT(ii,k)*floc(3,i,j,ii)
!!         ENDDO
!!                  
!!         divf(i,j,k,iVar,iEl) = df
!!
!!      ENDIF
!!    
!! END SUBROUTINE CalculateDivergence_3D_CUDAKernel
!! 
!!#endif

END MODULE Lagrange_Class
