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
    REAL(prec), POINTER :: iMatrix(:,:)
    REAL(prec), POINTER :: iMatrixTranspose(:,:)
    REAL(prec), POINTER :: dMatrix(:,:)  
    REAL(prec), POINTER :: dMatrixTranspose(:,:)  

    TYPE(c_ptr) :: bWeights_dev
    TYPE(c_ptr) :: iMatrix_dev
    TYPE(c_ptr) :: iMatrixTranspose_dev
    TYPE(c_ptr) :: dMatrix_dev
    TYPE(c_ptr) :: dMatrixTranspose_dev 

    CONTAINS
      
      PROCEDURE :: Build => Build_Lagrange
      PROCEDURE :: Trash => Trash_Lagrange

      PROCEDURE :: UpdateDevice => UpdateDevice_Lagrange
      PROCEDURE :: UpdateHost => UpdateHost_Lagrange
      
!      PROCEDURE :: Interpolate_1D => Interpolate_1D_Lagrange
!      PROCEDURE :: Interpolate_2D => Interpolate_2D_Lagrange
!      PROCEDURE :: Interpolate_3D => Interpolate_3D_Lagrange

!      PROCEDURE :: GridInterpolation_1D
!      PROCEDURE :: GridInterpolation_2D
!      PROCEDURE :: GridInterpolation_3D

!      PROCEDURE :: GridDerivative_1D
!
!      PROCEDURE :: CalculateGradient_2D
!      PROCEDURE :: CalculateDivergence_2D
!!      PROCEDURE :: CalculateCurl_2D
!
!      PROCEDURE :: CalculateGradient_3D
!      PROCEDURE :: CalculateDivergence_3D
!!      PROCEDURE :: CalculateCurl_3D

      PROCEDURE, PRIVATE :: CalculateBarycentricWeights  => CalculateBarycentricWeights_Lagrange
      PROCEDURE, PRIVATE :: CalculateInterpolationMatrix => CalculateInterpolationMatrix_Lagrange
      PROCEDURE, PRIVATE :: CalculateDerivativeMatrix    => CalculateDerivativeMatrix_Lagrange
      PROCEDURE          :: CalculateLagrangePolynomials  
      
    END TYPE Lagrange

 
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

  SUBROUTINE Build_Lagrange( myPoly, N, M, interpNodes, targetNodes )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(out)   :: myPoly
    INTEGER, INTENT(in)            :: N, M
    REAL(prec), INTENT(in)         :: interpNodes(0:N), targetNodes(0:M)
    
      myPoly % N  = N
      myPoly % M  = M
     
      ALLOCATE( myPoly % controlPoints(0:N), &
                myPoly % bWeights(0:N), &
                myPoly % targetPoints(0:M), &
                myPoly % iMatrix(0:M,0:N), &
                myPoly % iMatrixTranspose(0:N,0:M), &
                myPoly % dMatrix(0:N,0:N), &
                myPoly % dMatrixTranspose(0:N,0:N) )

      myPoly % controlPoints = 0.0_prec
      myPoly % bWeights = 0.0_prec
      myPoly % targetPoints  = 0.0_prec
      myPoly % iMatrix = 0.0_prec
      myPoly % iMatrixTranspose = 0.0_prec
      myPoly % dMatrix = 0.0_prec
      myPoly % dMatrixTranspose  = 0.0_prec
      
      myPoly % controlPoints(0:N) = interpNodes(0:N)
      myPoly % targetPoints(0:M) = targetNodes(0:M)

      CALL myPoly % CalculateBarycentricWeights( )
      CALL myPoly % CalculateInterpolationMatrix( )
      CALL myPoly % CalculateDerivativeMatrix( )

#ifdef GPU
      CALL hfMalloc(myPoly % bWeights_dev, SIZEOF(myPoly % bWeights))
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
!! ================================================================================================ !
!!
!! GridInterpolation_1D 
!!
!!   Interpolates an array of nodal values at the native nodes to nodal values at the target nodes.
!!   
!!   We can write the operations of interpolating data from one set of points to another (in this 
!!   case from "controlPoints" to "targetPoints") as
!!
!!              fnew_a = \sum_{i=0}^N f_i l_i(\xi_a),   a=0,1,2,...,M
!!             
!!   where l_i(\xi) are the Lagrange interpolating polynomials through the interpolation points.
!!   The interpolation matrix is T_{a,i} = l_i(\xi_a) maps an array of nodal values from the native
!!   interpolation nodes to the target nodes. This routine serves as a wrapper to call either the CUDA
!!   kernel (if CUDA is enabled) or the CPU version. 
!! 
!!   Usage : 
!!
!!     TYPE(Lagrange) :: interp
!!     INTEGER        :: nVariables, nElements
!!     REAL(prec)     :: f(0:interp % N,1:nVariables,1:nElements) 
!!     REAL(prec)     :: fNew(0:interp % M,1:nVariables,1:nElements) 
!!     
!!       CALL interp % GridInterpolation_1D( fnative, fNew, nVariables, nElements ) 
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
!!     fNew (out) 
!!      Array of function nodal values at the target interpolation nodes.
!!   
!! ================================================================================================ ! 
!
!  SUBROUTINE GridInterpolation_1D( myPoly, f, fNew, nVariables, nElements )  
!    IMPLICIT NONE
!    CLASS(Lagrange), INTENT(in) :: myPoly
!
!#ifdef GPU
!    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
!    REAL(prec), DEVICE, INTENT(in)  :: f(0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), DEVICE, INTENT(out) :: fNew(0:myPoly % M, 1:nVariables, 1:nElements)
!    TYPE(dim3) :: grid, tBlock
!    INTEGER    :: threadCount
!  
!      threadCount = MIN( 4*(ceiling( REAL(myPoly % M+1)/4 ) ), 8 )
!
!      tBlock = dim3( threadCount, 1, 1 )
!      grid   = dim3( nVariables, nElements, 1)
!     
!      CALL GridInterpolation_1D_CUDAKernel<<<grid, tBlock>>>( myPoly % iMatrixTranspose_dev, &
!                                                                     f, fNew, myPoly % N_dev, myPoly % M_dev, &
!                                                                     nVariables, nElements )
!
!#else
!    INTEGER, INTENT(in)     :: nVariables, nElements
!    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), INTENT(out) :: fNew(0:myPoly % M, 1:nVariables, 1:nElements)
!
!      fNew = GridInterpolation_1D_Lagrange( myPoly, f, nVariables, nElements )
!
!#endif
!
!  END SUBROUTINE GridInterpolation_1D
!
!! ================================================================================================ !
!!
!! GridInterpolation_2D 
!!
!!   Interpolates an array of nodal values at the native nodes to nodal values at the target nodes.
!!   
!!   We can write the operations of interpolating data from one set of points to another (in this 
!!   case from "controlPoints" to "targetPoints") as
!!
!!              fnew_{a,b} = \sum_{i,j=0}^N f_{i,j} l_i(\xi_a) l_j(\xi_b),   a,b=0,1,2,...,M
!!             
!!   where l_i(\xi) are the Lagrange interpolating polynomials through the interpolation points.
!!   The interpolation matrix is T_{a,i} = l_i(\xi_a) maps an array of nodal values from the native
!!   interpolation nodes to the target nodes. This routine serves as a wrapper to call either the CUDA
!!   kernel (if CUDA is enabled) or the CPU version. 
!! 
!!   Usage : 
!!
!!     TYPE(Lagrange) :: interp
!!     INTEGER        :: nVariables, nElements
!!     REAL(prec)     :: f(0:interp % N,0:interp % N,1:nVariables,1:nElements) 
!!     REAL(prec)     :: fNew(0:interp % M,0:interp % M,1:nVariables,1:nElements) 
!!     
!!       CALL interp % GridInterpolation_2D( fnative, fNew, nVariables, nElements ) 
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
!!     fNew (out) 
!!      Array of function nodal values at the target interpolation nodes.
!!   
!! ================================================================================================ ! 
!
!  SUBROUTINE GridInterpolation_2D( myPoly, f, fNew, nVariables, nElements )  
!    IMPLICIT NONE
!    CLASS(Lagrange), INTENT(in) :: myPoly
!
!#ifdef GPU
!    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
!    REAL(prec), DEVICE, INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), DEVICE, INTENT(out) :: fNew(0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
!    TYPE(dim3) :: grid, tBlock
!    INTEGER    :: threadCount
!  
!      threadCount = MIN( 4*(ceiling( REAL(myPoly % M+1)/4 ) ), 8 )
!
!      tBlock = dim3( threadCount, threadCount, 1 )
!      grid   = dim3( nVariables, nElements, 1)
!     
!      CALL GridInterpolation_2D_CUDAKernel<<<grid, tBlock>>>( myPoly % iMatrixTranspose_dev, &
!                                                                     f, fNew, myPoly % N_dev, myPoly % M_dev, &
!                                                                     nVariables, nElements )
!
!#else
!    INTEGER, INTENT(in)     :: nVariables, nElements
!    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), INTENT(out) :: fNew(0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
!
!      fNew = GridInterpolation_2D_Lagrange( myPoly, f, nVariables, nElements )
!
!#endif
!
!  END SUBROUTINE GridInterpolation_2D
!
!! ================================================================================================ !
!!
!! GridInterpolation_3D 
!!
!!   Interpolates an array of nodal values at the native nodes to nodal values at the target nodes.
!!   
!!   We can write the operations of interpolating data from one set of points to another (in this 
!!   case from "controlPoints" to "targetPoints") as
!!
!!   fnew_{a,b,c} = \sum_{i,j,k=0}^N f_{i,j,k} l_i(\xi_a) l_j(\xi_b) l_k(\xi_c),   a,b,c=0,1,2,...,M
!!             
!!   where l_i(\xi) are the Lagrange interpolating polynomials through the interpolation points.
!!   The interpolation matrix is T_{a,i} = l_i(\xi_a) maps an array of nodal values from the native
!!   interpolation nodes to the target nodes. This routine serves as a wrapper to call either the CUDA
!!   kernel (if CUDA is enabled) or the CPU version. 
!! 
!!   Usage : 
!!
!!     TYPE(Lagrange) :: interp
!!     INTEGER        :: nVariables, nElements
!!     REAL(prec)     :: f(0:interp % N,0:interp % N,0:interp % N,1:nVariables,1:nElements) 
!!     REAL(prec)     :: fNew(0:interp % M,0:interp % M,0:interp % M,1:nVariables,1:nElements) 
!!     
!!       CALL interp % GridInterpolation_3D( fnative, fNew, nVariables, nElements ) 
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
!!     fNew (out) 
!!      Array of function nodal values at the target interpolation nodes.
!!   
!! ================================================================================================ ! 
!
!  SUBROUTINE GridInterpolation_3D( myPoly, f, fNew, nVariables, nElements )  
!    IMPLICIT NONE
!    CLASS(Lagrange), INTENT(in) :: myPoly
!
!#ifdef GPU
!    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
!    REAL(prec), DEVICE, INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), DEVICE, INTENT(out) :: fNew(0:myPoly % M, 0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
!    TYPE(dim3) :: grid, tBlock
!    INTEGER    :: threadCount
!  
!      threadCount = MIN( 4*(ceiling( REAL(myPoly % M+1)/4 ) ), 8 )
!
!      tBlock = dim3( threadCount, threadCount, threadCount )
!      grid   = dim3( nVariables, nElements, 1)
!      CALL GridInterpolation_3D_CUDAKernel<<<grid, tBlock>>>( myPoly % iMatrixTranspose_dev, &
!                                                                     f, fNew, myPoly % N_dev, myPoly % M_dev,&
!                                                                     nVariables, nElements )
!
!#else
!    INTEGER, INTENT(in)     :: nVariables, nElements
!    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), INTENT(out) :: fNew(0:myPoly % M, 0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
!
!      fNew = GridInterpolation_3D_Lagrange( myPoly, f, nVariables, nElements )
!
!#endif
!
!  END SUBROUTINE GridInterpolation_3D
!
!! ================================================================================================ !
!!
!! GridDerivative_1D 
!!
!!   Calculates the derivative of the Lagrange interpolant given a set of nodal function values at
!!   the native interpolation nodes
!!   
!!   Given a set of nodal values at the interpolation nodes, the derivative of a function through 
!!   the interpolation nodes can be estimated by  
!!
!!                       f'_a = \sum_{i=0}^N f_{i} l'_i(\xi_a),   a=0,1,2,...,N
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
!!     REAL(prec)     :: f(0:interp % N,1:nVariables,1:nElements) 
!!     REAL(prec)     :: derF(0:interp % N,1:nVariables,1:nElements) 
!!     
!!       CALL interp % GridDerivative_1D( f, derF, nVariables, nElements ) 
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
!!     derF (out) 
!!      Array of derivative values at the target interpolation nodes.
!!   
!! ================================================================================================ ! 
!
!  SUBROUTINE GridDerivative_1D( myPoly, f, derF, nVariables, nElements )
!    IMPLICIT NONE
!    CLASS(Lagrange), INTENT(in) :: myPoly
!
!#ifdef GPU
!    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
!    REAL(prec), DEVICE, INTENT(in)  :: f(0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), DEVICE, INTENT(out) :: derF(0:myPoly % N, 1:nVariables, 1:nElements)
!    TYPE(dim3) :: grid, tBlock
!    INTEGER    :: threadCount
!  
!      threadCount = MIN( 4*(ceiling( REAL(myPoly % N+1)/4 ) ), 8 )
!
!      tBlock = dim3( threadCount, 1, 1 )
!      grid   = dim3( nVariables, nElements, 1)
!     
!      CALL GridDerivative_1D_CUDAKernel<<<grid, tBlock>>>( myPoly % dMatrixTranspose_dev, &
!                                                                f, derF, myPoly % N_dev , nVariables, nElements )
!
!#else
!    INTEGER, INTENT(in)     :: nVariables, nElements
!    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), INTENT(out) :: derF(0:myPoly % N, 1:nVariables, 1:nElements)
!
!      derF = GridDerivative_1D_Lagrange( myPoly, f, nVariables, nElements )
!
!#endif
!
!  END SUBROUTINE GridDerivative_1D
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
!  SUBROUTINE CalculateGradient_2D( myPoly, f, gradF, nVariables, nElements )
!    IMPLICIT NONE
!    CLASS(Lagrange), INTENT(in) :: myPoly
!
!#ifdef GPU
!    INTEGER, DEVICE, INTENT(in)     :: nVariables, nElements
!    REAL(prec), DEVICE, INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), DEVICE, INTENT(out) :: gradF(1:2,0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    TYPE(dim3) :: grid, tBlock
!    INTEGER    :: threadCount
!  
!      threadCount = MIN( 4*(ceiling( REAL(myPoly % N+1)/4 ) ), 8 )
!
!      tBlock = dim3( threadCount, threadCount, 1 )
!      grid   = dim3( nVariables, nElements, 1)
!
!      CALL CalculateGradient_2D_CUDAKernel<<<grid, tBlock>>>( myPoly % dMatrixTranspose_dev, &
!                                                                f, gradF, myPoly % N_dev, nVariables, nElements )
!
!#else
!    INTEGER, INTENT(in)     :: nVariables, nElements
!    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec), INTENT(out) :: gradF(1:2, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!
!      gradF = CalculateGradient_2D_Lagrange( myPoly, f, nVariables, nElements )
!
!#endif
!
!  END SUBROUTINE CalculateGradient_2D
!
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
!    ! Local
!    INTEGER :: iVar, iEl, i, a 
!
!      DO iEl = 1, nElements
!        DO iVar = 1, nVariables
!
!          DO a = 0, myPoly % M
!
!            fNew(a,iVar,iEl) = 0.0_prec
!
!            DO i = 0, myPoly % N
!
!              fNew(a,iVar,iEl) = fNew(a,iVar,iEl) + myPoly % iMatrixTranspose(i,a)*f(i,iVar,iEl)
!
!            ENDDO
!
!          ENDDO
!
!        ENDDO
!      ENDDO
!
!  END FUNCTION GridInterpolation_1D_Lagrange
!
!  FUNCTION GridInterpolation_2D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( fNew )  
!    IMPLICIT NONE
!    TYPE(Lagrange) :: myPoly
!    INTEGER         :: nVariables, nElements
!    REAL(prec)      :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec)      :: fNew(0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
!    ! Local
!    INTEGER :: i, j, a, b, p, iEl, iVar
!    REAL(prec) :: fa
!   
!      DO iEl = 1, nElements
!        DO iVar = 1, nVariables
!         
!          DO b = 0, myPoly % M
!            DO a = 0, myPoly % M
!             
!              fNew(a,b,iVar,iEl) = 0.0_prec
!
!              DO j = 0, myPoly % N
!                   
!                fa = 0.0_prec
!                DO i = 0, myPoly % N
!                  fa = fa + f(i,j,iVar,iEl)*myPoly % iMatrixTranspose(i,a)
!                ENDDO
!                      
!                fNew(a,b,iVar,iEl) = fNew(a,b,iVar,iEl) + fa*myPoly % iMatrixTranspose(j,b)
!
!              ENDDO
!                   
!            ENDDO
!          ENDDO
!
!        ENDDO   
!      ENDDO
!
!  END FUNCTION GridInterpolation_2D_Lagrange
!
!  FUNCTION GridInterpolation_3D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( fNew )  
!    IMPLICIT NONE
!    TYPE(Lagrange) :: myPoly
!    INTEGER         :: nElements, nVariables
!    REAL(prec)      :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
!    REAL(prec)      :: fNew(0:myPoly % M, 0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
!    ! Local
!    INTEGER :: i, j, k, a, b, c, iEl, iVar
!    REAL(prec) :: fa, fab
!   
!      DO iEl = 1, nElements
!        DO iVar = 1, nVariables
!         
!         
!          DO c = 0, myPoly % M
!            DO b = 0, myPoly % M
!              DO a = 0, myPoly % M
!               
!                fNew(a,b,c,iVar,iEl) = 0.0_prec
!
!                DO k = 0, myPoly % N
!                  
!                  fab = 0.0_prec
!                  DO j = 0, myPoly % N
!                     
!                    fa = 0.0_prec
!                    DO i = 0, myPoly % N
!                      fa = fa + f(i,j,k,iVar,iEl)*myPoly % iMatrixTranspose(i,a)
!                    ENDDO
!                        
!                    fab = fab + fa*myPoly % iMatrixTranspose(j,b)
!
!                  ENDDO
!                     
!                  fNew(a,b,c,iVar,iEl) = fNew(a,b,c,iVar,iEl) + fab*myPoly % iMatrixTranspose(k,c)
!
!                ENDDO
!                  
!
!              ENDDO
!            ENDDO
!          ENDDO
!
!        ENDDO   
!      ENDDO
!
!  END FUNCTION GridInterpolation_3D_Lagrange
!
!  FUNCTION GridDerivative_1D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( derF )  
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
!  END FUNCTION GridDerivative_1D_Lagrange
!!
!  FUNCTION CalculateGradient_2D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( gradF )  
!
!    IMPLICIT NONE
!    TYPE(Lagrange) :: myPoly
!    INTEGER         :: nVariables, nElements
!    REAL(prec)      :: f(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
!    REAL(prec)      :: gradF(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
!    ! Local
!    INTEGER    :: i, j, ii, iVar, iEl
!
!      DO iEl = 1, nElements
!        DO iVar = 1, nVariables
!          DO j = 0, myPoly % N
!            DO i = 0, myPoly % N
!    
!              gradF(1,i,j,iVar,iEl) = 0.0_prec 
!              gradF(2,i,j,iVar,iEl) = 0.0_prec 
!              DO ii = 0, myPoly % N
!                gradF(1,i,j,iVar,iEl) = gradF(1,i,j,iVar,iEl) + myPoly % dMatrixTranspose(ii,i)*f(ii,j,iVar,iEl)
!                gradF(2,i,j,iVar,iEl) = gradF(2,i,j,iVar,iEl) + myPoly % dMatrixTranspose(ii,j)*f(i,ii,iVar,iEl)
!              ENDDO
!    
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO
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
!!  ATTRIBUTES(Global) SUBROUTINE GridDerivative_1D_CUDAKernel( DMatT, f, derf, N, nVariables, nElems )
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
!!  END SUBROUTINE GridDerivative_1D_CUDAKernel
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
