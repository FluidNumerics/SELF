! Lagrange_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE Lagrange_Class

!src/common
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines

#ifdef HAVE_CUDA
USE cudafor
#endif

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

    INTEGER                 :: N     
    INTEGER                 :: M 
    REAL(prec), ALLOCATABLE :: interpolationPoints(:)
    REAL(prec), ALLOCATABLE :: targetPoints(:)
    REAL(prec), ALLOCATABLE :: barycentricWeights(:)
    REAL(prec), ALLOCATABLE :: interpolationMatrix(:,:)
    REAL(prec), ALLOCATABLE :: interpolationMatrixTranspose(:,:)
    REAL(prec), ALLOCATABLE :: derivativeMatrix(:,:)  
    REAL(prec), ALLOCATABLE :: derivativeMatrixTranspose(:,:)  

#ifdef HAVE_CUDA
    INTEGER, ALLOCATABLE, DEVICE    :: N_dev, M_dev
    REAL(prec), ALLOCATABLE, DEVICE :: barycentricWeights_dev(:)  
    REAL(prec), ALLOCATABLE, DEVICE :: interpolationMatrix_dev(:,:)
    REAL(prec), ALLOCATABLE, DEVICE :: interpolationMatrixTranspose_dev(:,:)
    REAL(prec), ALLOCATABLE, DEVICE :: derivativeMatrix_dev(:,:)  
    REAL(prec), ALLOCATABLE, DEVICE :: derivativeMatrixTranspose_dev(:,:)  
#endif

    CONTAINS
      
      PROCEDURE :: Build => Build_Lagrange
      PROCEDURE :: Trash => Trash_Lagrange
      
      PROCEDURE :: Interpolate_1D => Interpolate_1D_Lagrange
      PROCEDURE :: Interpolate_2D => Interpolate_2D_Lagrange
      PROCEDURE :: Interpolate_3D => Interpolate_3D_Lagrange

      PROCEDURE :: ApplyInterpolationMatrix_1D
      PROCEDURE :: ApplyInterpolationMatrix_2D
      PROCEDURE :: ApplyInterpolationMatrix_3D

      PROCEDURE :: CalculateDerivative_1D

      PROCEDURE :: CalculateGradient_2D
      PROCEDURE :: CalculateDivergence_2D
!      PROCEDURE :: CalculateCurl_2D

      PROCEDURE :: CalculateGradient_3D
      PROCEDURE :: CalculateDivergence_3D
!      PROCEDURE :: CalculateCurl_3D

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
!       The target nodes. The interpolationMatrix will map the a function at the interpolation
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
     
      ALLOCATE( myPoly % interpolationPoints(0:N), &
                myPoly % barycentricWeights(0:N), &
                myPoly % targetPoints(0:M), &
                myPoly % interpolationMatrix(0:M,0:N), &
                myPoly % interpolationMatrixTranspose(0:N,0:M), &
                myPoly % derivativeMatrix(0:N,0:N), &
                myPoly % derivativeMatrixTranspose(0:N,0:N) )

      myPoly % interpolationPoints          = 0.0_prec
      myPoly % barycentricWeights           = 0.0_prec
      myPoly % targetPoints                 = 0.0_prec
      myPoly % interpolationMatrix          = 0.0_prec
      myPoly % interpolationMatrixTranspose = 0.0_prec
      myPoly % derivativeMatrix             = 0.0_prec
      myPoly % derivativeMatrixTranspose    = 0.0_prec
      
      myPoly % interpolationPoints(0:N) = interpNodes(0:N)
      myPoly % targetPoints(0:M)        = targetNodes(0:M)

      CALL myPoly % CalculateBarycentricWeights( )
      CALL myPoly % CalculateInterpolationMatrix( )
      CALL myPoly % CalculateDerivativeMatrix( )

#ifdef HAVE_CUDA
      ALLOCATE( myPoly % N_dev, myPoly % M_dev )
      ALLOCATE( myPoly % barycentricWeights_dev(0:N), &
                myPoly % interpolationMatrix_dev(0:M,0:N), &
                myPoly % interpolationMatrixTranspose_dev(0:N,0:M), &
                myPoly % derivativeMatrix_dev(0:N,0:N), &
                myPoly % derivativeMatrixTranspose_dev(0:N,0:N) )
 
      myPoly % N_dev = N
      myPoly % M_dev = M
      
      myPoly % barycentricWeights_dev           = myPoly % barycentricWeights
      myPoly % interpolationMatrix_dev          = myPoly % interpolationMatrix
      myPoly % interpolationMatrixTranspose_dev = myPoly % interpolationMatrixTranspose
      myPoly % derivativeMatrix_dev             = myPoly % derivativeMatrix
      myPoly % derivativeMatrixTranspose_dev    = myPoly % derivativeMatrixTranspose

#endif
 
 END SUBROUTINE Build_Lagrange

! ================================================================================================ !
!
! Trash_Lagrange
!
!   A manual destructor for the Lagrange class that deallocates the memory held by its attributes. 
! 
!   Usage :
!
!     TYPE(Lagrange) :: interp
!     
!       CALL interp % Trash( ) 
!  
!   Input/Output :
! 
!     interp (in/out)  
!       On input, a previously constructed Lagrange data structure. 
!       On output, the memory associated with this data-structure is freed.
!   
! ================================================================================================ ! 

  SUBROUTINE Trash_Lagrange( myPoly )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly

      DEALLOCATE( myPoly % interpolationPoints, &
                  myPoly % barycentricWeights, &
                  myPoly % targetPoints, &
                  myPoly % interpolationMatrix, &
                  myPoly % interpolationMatrixTranspose, &
                  myPoly % derivativeMatrix, &
                  myPoly % derivativeMatrixTranspose )
#ifdef HAVE_CUDA
      DEALLOCATE( myPoly % N_dev, myPoly % M_dev )
      DEALLOCATE( myPoly % barycentricWeights_dev, &
                  myPoly % interpolationMatrix_dev, &
                  myPoly % interpolationMatrixTranspose_dev, &
                  myPoly % derivativeMatrix_dev, &
                  myPoly % derivativeMatrixTranspose_dev )
#endif

  END SUBROUTINE Trash_Lagrange

! ================================================================================================ !
!
! Interpolate_1D_Lagrange
!
!   Interpolates an array of data onto a desired location using Lagrange interpolation.
!
!   Usage :
!
!     TYPE(Lagrange) :: interp
!     REAL(prec)     :: f(0:interp % N) 
!     REAL(prec)     :: sE 
!     REAL(prec)     :: fAtSE 
!         
!       fAtSE = interp % Interpolate_1D( f, sE ) 
! 
!   Parameters :
!  
!     myPoly (in)
!       A previously constructed Lagrange data-structure.
!
!     f (in)  
!       An array of function nodal values located at the native interpolation nodes.
!
!     sE (in) 
!       The location where you want to interpolate to.
!
!     fAtSE (out)
!       The interpolant evaluated at sE.  
!
! ================================================================================================ ! 

  FUNCTION Interpolate_1D_Lagrange( myPoly, f, sE ) RESULT( interpF )  
    IMPLICIT NONE
    CLASS(Lagrange) :: myPoly
    REAL(prec)      :: sE
    REAL(prec)      :: f(0:myPoly % N)
    REAL(prec)      :: interpF
    ! Local
    REAL(prec) :: lAtS(0:myPoly % N)
    INTEGER    :: i
   
      lAtS = myPoly % CalculateLagrangePolynomials( sE )

      interpF = 0.0_prec
      DO i = 0, myPoly % N
        interpF = interpF + lAtS(i)*f(i)
      ENDDO
    
  END FUNCTION Interpolate_1D_Lagrange

! ================================================================================================ !
!
! Interpolate_2D_Lagrange
!
!   Interpolates an array of data onto a desired location using Lagrange interpolation.
!
!   Usage :
!
!     TYPE(Lagrange) :: interp
!     REAL(prec)     :: f(0:interp % N,0:interp % N) 
!     REAL(prec)     :: sE(1:2) 
!     REAL(prec)     :: fAtSE 
!         
!       fAtSE = interp % Interpolate_2D( f, sE ) 
! 
!   Parameters :
!  
!     myPoly (in)
!       A previously constructed Lagrange data-structure.
!
!     f (in)  
!       An array of function nodal values located at the native interpolation nodes.
!
!     sE (in) 
!       The location where you want to interpolate to.
!
!     fAtSE (out)
!       The interpolant evaluated at sE.  
!
! ================================================================================================ ! 

  FUNCTION Interpolate_2D_Lagrange( myPoly, f, sE ) RESULT( interpF )  
    IMPLICIT NONE
    CLASS(Lagrange) :: myPoly
    REAL(prec)      :: sE(1:2)
    REAL(prec)      :: f(0:myPoly % N, 0:myPoly % N)
    REAL(prec)      :: interpF
    ! Local
    REAL(prec) :: fj
    REAL(prec) :: ls(0:myPoly % N)
    REAL(prec) :: lp(0:myPoly % N)
    INTEGER    :: i, j

      ls = myPoly % CalculateLagrangePolynomials( sE(1) ) 
      lp = myPoly % CalculateLagrangePolynomials( sE(2) )
      
      interpF = 0.0_prec
      DO j = 0, myPoly % N
     
        fj = 0.0_prec
        DO i = 0, myPoly % N
          fj = fj + f(i,j)*ls(i)
        ENDDO
            
        interpF = interpF + fj*lp(j)

      ENDDO
      
 END FUNCTION Interpolate_2D_Lagrange

! ================================================================================================ !
!
! Interpolate_3D_Lagrange
!
!   Interpolates an array of data onto a desired location using Lagrange interpolation.
!
!   Usage :
!
!     TYPE(Lagrange) :: interp
!     REAL(prec)     :: f(0:interp % N,0:interp % N,0:interp % N) 
!     REAL(prec)     :: sE(1:3) 
!     REAL(prec)     :: fAtSE 
!         
!       fAtSE = interp % Interpolate_3D( f, sE ) 
! 
!   Parameters :
!  
!     myPoly (in)
!       A previously constructed Lagrange data-structure.
!
!     f (in)  
!       An array of function nodal values located at the native interpolation nodes.
!
!     sE (in) 
!       The location where you want to interpolate to.
!
!     fAtSE (out)
!       The interpolant evaluated at sE.  
!
! ================================================================================================ ! 

  FUNCTION Interpolate_3D_Lagrange( myPoly, f, sE ) RESULT( interpF )  
    IMPLICIT NONE
    CLASS(Lagrange) :: myPoly
    REAL(prec)      :: sE(1:3)
    REAL(prec)      :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N)
    REAL(prec)      :: interpF
    ! Local
    REAL(prec) :: fjk, fk
    REAL(prec) :: ls(0:myPoly % N)
    REAL(prec) :: lp(0:myPoly % N)
    REAL(prec) :: lq(0:myPoly % N)
    INTEGER    ::  i, j, k

      ls = myPoly % CalculateLagrangePolynomials( sE(1) ) 
      lp = myPoly % CalculateLagrangePolynomials( sE(2) )
      lq = myPoly % CalculateLagrangePolynomials( sE(3) )
      
      interpF = 0.0_prec
      DO k = 0, myPoly % N
      
         fk = 0.0_prec
         DO j = 0, myPoly % N
         
            fjk = 0.0_prec
            DO i = 0, myPoly % N
               fjk = fjk + f(i,j,k)*ls(i)
            ENDDO
            
            fk = fk + fjk*lp(j)
         ENDDO
         
         interpF = interpF + fk*lq(k)
      ENDDO
      
  END FUNCTION Interpolate_3D_Lagrange

!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R ApplyInterpolationMatrix_1D 
! 
!> \fn ApplyInterpolationMatrix_1D_Lagrange 
!! Interpolates an array of nodal values at the native nodes to nodal values at the target nodes.
!! 
!! As described in calculateinterpolationmatrix_lagrange_1d, 
!! we can write the operations of interpolating data from one set of points to another (in this 
!! case from "s" to "so") as
!!            \f[ f_j = \sum_{i=0}^N f_i l_i(\xi_j), \hspace{2mm} j=0,1,2,...,M
!!            \f]
!! where \f$ l_i(\xi) \f$ are the Lagrange interpolating polynomials at the 
!! \f$ \lbrace \xi_i \rbrace_{i=0}^N \f$ nodes. This routine performs the matrix-multiply that
!! maps an array of nodal values from the native interpolation nodes to the target nodes 
!! 
!!  Usage : </H2> 
!! TYPE</B>(Lagrange) :: this 
!! REAL</B>(prec)     :: fnative(0:this % N) 
!! REAL</B>(prec)     :: ftarget(0:this % M) 
!!         .... 
!!     ftarget = this % ApplyInterpolationMatrix_1D( fnative ) 
!! 
!!   Parameters : </H2>
!!  <table> 
!!     in  myPoly  TYPE(Lagrange)  
!!                     A previously constructed Lagrange data-structure.
!!     in  f(0:myPoly % N)  REAL(prec) 
!!                     Array of function nodal values at the native interpolation nodes.
!!     out  fNew(0:myPoly % M)  REAL(prec)  
!!                     Array of function nodal values at the target interpolation nodes.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}

  SUBROUTINE ApplyInterpolationMatrix_1D( myPoly, f, fNew, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly

#ifdef HAVE_CUDA
    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: fNew(0:myPoly % M, 1:nVariables, 1:nElements)
    TYPE(dim3) :: grid, tBlock
    INTEGER    :: threadCount
  
      threadCount = MIN( 4*(ceiling( REAL(myPoly % M+1)/4 ) ), 8 )

      tBlock = dim3( threadCount, 1, 1 )
      grid   = dim3( nVariables, nElements, 1)
     
      CALL ApplyInterpolationMatrix_1D_CUDAKernel<<<grid, tBlock>>>( myPoly % interpolationMatrixTranspose_dev, &
                                                                     f, fNew, myPoly % N_dev, myPoly % M_dev, &
                                                                     nVariables, nElements )

#else
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: fNew(0:myPoly % M, 1:nVariables, 1:nElements)

      fNew = ApplyInterpolationMatrix_1D_Lagrange( myPoly, f, nVariables, nElements )

#endif

  END SUBROUTINE ApplyInterpolationMatrix_1D

!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R ApplyInterpolationMatrix_2D 
! 
!> \fn ApplyInterpolationMatrix_2D_Lagrange  
!! Interpolates an array of nodal values at the native nodes to nodal values at the target nodes.
!! 
!! As described in \ref calculateinterpolationmatrix_Lagrange, 
!! we can write the operations of interpolating data from one set of points to another (in this 
!! case from "s" to "so") as a sequence of matrix-matrix multiplications
!! \f[ 
!!       \tilde{f} = Tf 
!! \f]
!! \f[ 
!!       f_{target} = \tilde{f} T^T 
!! \f]
!!
!! This routine performs the matrix-multiplications that map an array of nodal values from the 
!! native interpolation nodes to the target nodes. 
!! 
!!  Usage : </H2> 
!! TYPE</B>(Lagrange) :: this 
!! REAL</B>(prec)     :: fnative(0:this % N,0:this % N) 
!! REAL</B>(prec)     :: ftarget(0:this % M,0:this % M) 
!!         .... 
!!     ftarget = this % ApplyInterpolationMatrix_2D( fnative ) 
!! 
!!   Parameters : </H2>
!!  <table> 
!!     in  myPoly  TYPE(Lagrange)  
!!                     A previously constructed Lagrange data-structure.
!!     in  f(0:myPoly % N,0:this % N)  REAL(prec) 
!!                     2-D Array of function nodal values at the native interpolation nodes.
!!     out  fNew(0:myPoly % M,0:this % M)  REAL(prec)  
!!                     2-D Array of function nodal values at the target interpolation nodes.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}

  SUBROUTINE ApplyInterpolationMatrix_2D( myPoly, f, fNew, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly

#ifdef HAVE_CUDA
    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: fNew(0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
    TYPE(dim3) :: grid, tBlock
    INTEGER    :: threadCount
  
      threadCount = MIN( 4*(ceiling( REAL(myPoly % M+1)/4 ) ), 8 )

      tBlock = dim3( threadCount, threadCount, 1 )
      grid   = dim3( nVariables, nElements, 1)
     
      CALL ApplyInterpolationMatrix_2D_CUDAKernel<<<grid, tBlock>>>( myPoly % interpolationMatrixTranspose_dev, &
                                                                     f, fNew, myPoly % N_dev, myPoly % M_dev, &
                                                                     nVariables, nElements )

#else
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: fNew(0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)

      fNew = ApplyInterpolationMatrix_2D_Lagrange( myPoly, f, nVariables, nElements )

#endif

  END SUBROUTINE ApplyInterpolationMatrix_2D

!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R ApplyInterpolationMatrix_3D 
! 
!> \fn ApplyInterpolationMatrix_3D_Lagrange  
!! Interpolates an array of nodal values at the native nodes to nodal values at the target nodes.
!! 
!! As described in \ref calculateinterpolationmatrix_lagrange_3D, 
!! we can write the operations of interpolating data from one set of points to another (in this 
!! case from "s" to "so") as a sequence of matrix-matrix multiplications. First a sequence of 2-D
!! interpolations are performed,
!! \f[ 
!!       \tilde{f}_k = Tf_kT ^T
!! \f]
!! Then, by collapsing the first two dimensions of the data, a single matrix-matrix multiplication
!! followed by an unpacking of the first two dimensions results in the data interpolated onto the 
!! target nodes
!! \f[ 
!!       f_{target} = UNPACK( PACK( \tilde{f} T^T ) ) 
!! \f]
!!
!! 
!!  Usage : </H2> 
!! TYPE</B>(Lagrange) :: this 
!! REAL</B>(prec)     :: fnative(0:this % N,0:this % N,0:this % N) 
!! REAL</B>(prec)     :: ftarget(0:this % M,0:this % M,0:this % M) 
!!         .... 
!!     ftarget = this % ApplyInterpolationMatrix_3D( fnative ) 
!! 
!!   Parameters : </H2>
!!  <table> 
!!     in  myPoly  TYPE(Lagrange)  
!!                     A previously constructed Lagrange data-structure.
!!     in  f(0:myPoly % N,0:this % N,0:this % N)  REAL(prec) 
!!                     3-D Array of function nodal values at the native interpolation nodes.
!!     out  fNew(0:myPoly % M,0:this % M,0:this % M)  REAL(prec)  
!!                     3-D Array of function nodal values at the target interpolation nodes.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}

  SUBROUTINE ApplyInterpolationMatrix_3D( myPoly, f, fNew, nVariables, nElements )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly

#ifdef HAVE_CUDA
    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: fNew(0:myPoly % M, 0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
    TYPE(dim3) :: grid, tBlock
    INTEGER    :: threadCount
  
      threadCount = MIN( 4*(ceiling( REAL(myPoly % M+1)/4 ) ), 8 )

      tBlock = dim3( threadCount, threadCount, threadCount )
      grid   = dim3( nVariables, nElements, 1)
      CALL ApplyInterpolationMatrix_3D_CUDAKernel<<<grid, tBlock>>>( myPoly % interpolationMatrixTranspose_dev, &
                                                                     f, fNew, myPoly % N_dev, myPoly % M_dev,&
                                                                     nVariables, nElements )

#else
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: fNew(0:myPoly % M, 0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)

      fNew = ApplyInterpolationMatrix_3D_Lagrange( myPoly, f, nVariables, nElements )

#endif

  END SUBROUTINE ApplyInterpolationMatrix_3D

!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R CalculateDerivative_1D 
! 
!> \fn CalculateDerivative_1D_Lagrange 
!! Calculates the derivative of the Lagrange interpolant given a set of nodal function values at
!! the native interpolation nodes
!! 
!! As described in calculatederivativematrix_lagrange_1d, 
!! given nodal values of an interpolant, the derivative can be estimated at the interpolation 
!! nodes using the summation
!!      \f[ \lbrace f' \rbrace_{j=0}^N = \sum_{i}^N( f_i l'_i(\xi_j) )
!!      \f]
!! Where \f$ l'_i(\xi_j) \f$ is the derivative of the \f$ i^{th} \f$ Lagrange interpolating 
!! polynomial at the \f$ j^{th} \f$ native interpolation node. This routine performs the 
!! matrix-multiply that results in an estimate of the derivative of a function whose nodal values
!! are the  \f$ f_i \f$ at the native interpolation nodes.
!! 
!!  Usage : </H2> 
!! TYPE</B>(Lagrange) :: this 
!! REAL</B>(prec)     :: fnative(0:this % N) 
!! REAL</B>(prec)     :: dfds(0:this % N) 
!!         .... 
!!     dfds = this % ApplyDerivativeMatrix_1D( fnative ) 
!! 
!!   Parameters : </H2>
!!  <table> 
!!     in  myPoly  TYPE(Lagrange)  
!!                     A previously constructed Lagrange data-structure.
!!     in  f(0:myPoly % N)  REAL(prec) 
!!                     Array of function nodal values at the native interpolation nodes.
!!     out  derF(0:myPoly % M)  REAL(prec)  
!!                     Array of estimated derivative values at the native interpolation nodes.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}

  SUBROUTINE CalculateDerivative_1D( myPoly, f, derF, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly

#ifdef HAVE_CUDA
    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: derF(0:myPoly % N, 1:nVariables, 1:nElements)
    TYPE(dim3) :: grid, tBlock
    INTEGER    :: threadCount
  
      threadCount = MIN( 4*(ceiling( REAL(myPoly % N+1)/4 ) ), 8 )

      tBlock = dim3( threadCount, 1, 1 )
      grid   = dim3( nVariables, nElements, 1)
     
      CALL CalculateDerivative_1D_CUDAKernel<<<grid, tBlock>>>( myPoly % derivativeMatrixTranspose_dev, &
                                                                f, derF, myPoly % N_dev , nVariables, nElements )

#else
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: derF(0:myPoly % N, 1:nVariables, 1:nElements)

      derF = CalculateDerivative_1D_Lagrange( myPoly, f, nVariables, nElements )

#endif

  END SUBROUTINE CalculateDerivative_1D

  SUBROUTINE CalculateGradient_2D( myPoly, f, gradF, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly

#ifdef HAVE_CUDA
    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: gradF(1:2,0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    TYPE(dim3) :: grid, tBlock
    INTEGER    :: threadCount
  
      threadCount = MIN( 4*(ceiling( REAL(myPoly % N+1)/4 ) ), 8 )

      tBlock = dim3( threadCount, threadCount, 1 )
      grid   = dim3( nVariables, nElements, 1)

      CALL CalculateGradient_2D_CUDAKernel<<<grid, tBlock>>>( myPoly % derivativeMatrixTranspose_dev, &
                                                                f, gradF, myPoly % N_dev, nVariables, nElements )

#else
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: gradF(1:2, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)

      gradF = CalculateGradient_2D_Lagrange( myPoly, f, nVariables, nElements )

#endif

  END SUBROUTINE CalculateGradient_2D


!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R CalculateDivergence_2D 
! 
!> \fn CalculateDivergence_2D_Lagrange  
!! Calculates the derivative of the Lagrange interpolant, in each computational direction, given a 
!! set of nodal function values at the native interpolation nodes.
!!
!! As described in \ref calculatederivativematrix_Lagrange, 
!! we can write the derivative calculations as a set of two matrix-matrix products
!! \f[ 
!!       \frac{\partial f}{\partial s} = D f  
!! \f]
!! \f[ 
!!       \frac{\partial f}{\partial p} = f D^T 
!! \f]
!! This routine performs the matrix-multiplications that result in the derivative of the interpolant
!! at the native interpolation nodes. This serves as an estimate of the derivative of the underlying
!! function.
!! 
!!  Usage : </H2> 
!! TYPE</B>(Lagrange) :: this 
!! REAL</B>(prec)     :: fnative(0:this % N,0:this % N) 
!! REAL</B>(prec)     :: derF(0:this % N,0:this % N,1:2) 
!!         .... 
!!     derF = this % ApplyDerivativeMatrix_2D( fnative ) 
!! 
!!   Parameters : </H2>
!!  <table> 
!!     in  myPoly  TYPE(Lagrange)  
!!                     A previously constructed Lagrange data-structure.
!!     in  f(0:myPoly % N,0:this % N)  REAL(prec) 
!!                     2-D Array of function nodal values at the native Derivative nodes.
!!     out  derF(0:myPoly % N,0:this % N,1:2)  REAL(prec)  
!!                     3-D Array containing the derivative (in each direction) of the interpolant at
!!                     the native interpolation nodes.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}

  SUBROUTINE CalculateDivergence_2D( myPoly, f, divF, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly

#ifdef HAVE_CUDA
    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(1:2, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: divF(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    TYPE(dim3) :: grid, tBlock
    INTEGER    :: threadCount
  
      threadCount = MIN( 4*(ceiling( REAL(myPoly % N+1)/4 ) ), 8 )

      tBlock = dim3( threadCount, threadCount, 1 )
      grid   = dim3( nVariables, nElements, 1)

      CALL CalculateDivergence_2D_CUDAKernel<<<grid, tBlock>>>( myPoly % derivativeMatrixTranspose_dev, &
                                                                f, divF, myPoly % N_dev, nVariables, nElements )

#else
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:2, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: divF(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)

      divF = CalculateDivergence_2D_Lagrange( myPoly, f, nVariables, nElements )

#endif

  END SUBROUTINE CalculateDivergence_2D

  SUBROUTINE CalculateGradient_3D( myPoly, f, gradF, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly

#ifdef HAVE_CUDA
    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: gradF(1:3,0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    TYPE(dim3) :: grid, tBlock
    INTEGER    :: threadCount
  
      threadCount = MIN( 4*(ceiling( REAL(myPoly % N+1)/4 ) ), 8 )

      tBlock = dim3( threadCount, threadCount, threadCount )
      grid   = dim3( nVariables, nElements, 1)

      CALL CalculateGradient_3D_CUDAKernel<<<grid, tBlock>>>( myPoly % derivativeMatrixTranspose_dev, &
                                                              f, gradF, myPoly % N_dev, nVariables, nElements )

#else
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: gradF(1:3, 0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)

      gradF = CalculateGradient_3D_Lagrange( myPoly, f, nVariables, nElements )

#endif

  END SUBROUTINE CalculateGradient_3D


!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R CalculateDivergence_3D 
! 
!> \fn CalculateDivergence_3D_Lagrange  
!! Calculates the derivative of the Lagrange interpolant, in each computational direction, given a 
!! set of nodal function values at the native interpolation nodes.
!!
!! As described in \ref calculatederivativematrix_lagrange_3D, to compute the derivative in a single
!! computational direction, the other array dimensions can be collapsed enabling a matrix-matrix 
!! product. After computing the product, the result can be unpacked into a 3-D array.
!!
!! This routine performs the matrix-multiplications that result in the derivative of the interpolant
!! at the native interpolation nodes. This serves as an estimate of the derivative of the underlying
!! function.
!! 
!!  Usage : </H2> 
!! TYPE</B>(Lagrange) :: this 
!! REAL</B>(prec)     :: fnative(0:this % N,0:this % N,0:this % N) 
!! REAL</B>(prec)     :: derF(0:this % N,0:this % N,0:this % N,1:3) 
!!         .... 
!!     derF = this % ApplyDerivativeMatrix_3D( fnative ) 
!! 
!!   Parameters : </H2>
!!  <table> 
!!     in  myPoly  TYPE(Lagrange)  
!!                     A previously constructed Lagrange data-structure.
!!     in  f(0:myPoly % N,0:this % N,0:this % N)  REAL(prec) 
!!                     3-D Array of function nodal values at the native Derivative nodes.
!!     out  derF(0:myPoly % N,0:this % N,1:2)  REAL(prec)  
!!                     4-D Array containing the derivative of the interpolant (in each direction) at
!!                     the nativeinterpolation nodes.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}

  SUBROUTINE CalculateDivergence_3D( myPoly, f, divF, nVariables, nElements )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(in) :: myPoly

#ifdef HAVE_CUDA
    INTEGER, DEVICE, INTENT(in)    :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(1:3, 0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: divF(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    TYPE(dim3) :: grid, tBlock
    INTEGER    :: threadCount
  
      threadCount = MIN( 4*(ceiling( REAL(myPoly % N+1)/4 ) ), 8 )

      tBlock = dim3( threadCount, threadCount, threadCount )
      grid   = dim3( nVariables, nElements, 1)

      CALL CalculateDivergence_3D_CUDAKernel<<<grid, tBlock>>>( myPoly % derivativeMatrixTranspose_dev, &
                                                                f, divF, myPoly % N_dev, nVariables, nElements )

#else
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:3, 0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: divF(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)

      divF = CalculateDivergence_3D_Lagrange( myPoly, f, nVariables, nElements )

#endif

  END SUBROUTINE CalculateDivergence_3D
 
! ================================================================================================ !
! ------------------------------------- PRIVATE ROUTINES ----------------------------------------- !
! ================================================================================================ !


!> \addtogroup Lagrange_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateBarycentricWeights (PRIVATE)
! 
!> \fn CalculateBarycentricWeights_Lagrange
!!  A PRIVATE routine that calculates and stores the barycentric weights for the Lagrange 
!!  data-structure.
!! 
!!  Calculates the barycentric weights from the interpolation nodes and stores them in the "barycentricWeights"
!!  attribute. This routine should be called after the native interpolation nodes have been 
!!  assigned.
!!
!! From a set of interpolation nodes \f$ \lbrace \xi_j \rbrace_{j=0}^N \f$, the Lagrange 
!! interpolating polynomials are given as
!!     \f[
!!           l_j = \prod_{i=0,\neq j}^N \frac{\xi-\xi_i}{\xi_j-\xi_i}
!!     \f]
!! For efficient interpolation with favorable round-off error, the "barycentric weights" are usually 
!! stored when performing Lagrange interpolation. The barycentric weights are
!!     \f[
!!           w_j = \prod_{i=0,\neq j}^N \frac{1}{\xi_j-\xi_i}
!!     \f]
!! 
!!   This routine is from Alg. 30 on pg. 75 of D.A. Kopriva, 2009.
!! 
!!  Usage : </H2> 
!! TYPE</B>(Lagrange) :: this 
!!         .... 
!!     CALL</B> this % CalculateBarycentricWeights( ) 
!! 
!!   Parameters : </H2>
!!  <table> 
!!     in/out  myPoly  TYPE(Lagrange) 
!!           On input</B>, myPoly is the Lagrange data structure is sent in with the 
!!           native interpolation nodes already filled in. 
!!           On output</B>, myPoly has the barycentric weights filled in.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}

  SUBROUTINE CalculateBarycentricWeights_Lagrange( myPoly )
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly
    ! Local
    INTEGER :: i, j
   
      DO i = 0, myPoly % N
        myPoly % barycentricWeights(i) = 1.0_prec
      ENDDO

      ! Computes the product w_k = w_k*(s_k - s_j), k /= j
      DO j = 1, myPoly % N
        DO i = 0, j-1

          myPoly % barycentricWeights(i) = myPoly % barycentricWeights(i)*&
                                           ( myPoly % interpolationPoints(i) - myPoly % interpolationPoints(j) )
          myPoly % barycentricWeights(j) = myPoly % barycentricWeights(j)*&
                                           ( myPoly % interpolationPoints(j) - myPoly % interpolationPoints(i) )

         ENDDO 
      ENDDO 
 
      DO j = 0, myPoly % N
        myPoly % barycentricWeights(j) = 1.0_prec/myPoly % barycentricWeights(j)
      ENDDO 

  END SUBROUTINE CalculateBarycentricWeights_Lagrange

!> \addtogroup Lagrange_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateInterpolationMatrix (PRIVATE) 
! 
!> \fn CalculateInterpolationMatrix_Lagrange 
!! A PRIVATE routine that fills in the interpolation matrix for the Lagrange data structure.
!! 
!! Given the interpolation formula
!!  \f[
!!       I_N(f) = \sum_{i=0}^N f_i l_i(\xi)
!!  \f] 
!! we can map \f$ \lbrace f_i \rbrace_{i=0}^N\f$ to \f$ \lbrace \tilde{f}_j \rbrace_{j=0}^M\f$ by
!! computing
!!  \f[
!!       \tilde{f}_j = \sum_{i=0}^N f_i l_i(\xi_j)
!!  \f]
!! Row j, column i of the "interpolation matrix" is 
!!  \f[
!!     T_{j,i} = l_i(\xi_j)
!!  \f]
!!
!!   This function is from Alg. 32 on pg. 76 of D.A. Kopriva, 2009.
!!
!! We can write the operation of interpolating data from the interpolation points to the target points as
!!
!!            \f[ f_{m,n} = \sum_{i,j=0}^N f_{i,j} l_i(\xi^1_m)l_j(\xi^2_n) , \hspace{2mm} m,n=0,1,2,...,M
!!            \f]
!!
!! where \f$ l_i(\xi^1) \f$ and \f$ l_j(\xi^2) \f$ are the Lagrange interpolating polynomials at the 
!! \f$ \lbrace ( \xi^1_i, \xi^2_j ) \rbrace_{i,j=0}^N \f$ nodes. Evaluation of the Lagrange 
!! interpolating polynomials at each of the new points 
!! \f$ \lbrace (\xi^1_m, \xi^2_n )\rbrace_{m,n=0}^M \f$ in each computational directions can be
!! written as two matrices, where 
!!       \f[ T^{(1)}_{m,i} = l_i(\xi^1_m).
!!        \f]
!! and 
!!       \f[ T^{(2)}_{n,j} = l_j(\xi^2_n).
!!        \f]
!! For simplicity, the native interpolation nodes are assumed identical in each computational
!! direction. Similarly, the target nodes are assumed identical in each direction. This assumption
!! implies that \f$ T^{(1)} \f$ and \f$ T^{(2)} \f$ are identical.
!!
!! In the SELF, interpolation onto the target grid (in 2-D) is executed via two matrix-matrix 
!! multiplications. The native data, \f$ f_{i,j} \f$ can be viewed as an \f$ N+1 \times N+1 \f$ 
!! matrix. Matrix multiplication (on the left) by \f$ T \f$ maps the native data to the target nodes
!! in the \f$ \xi^1 \f$ direction,
!!
!! \f[
!!      \tilde{f} = T f 
!! \f]
!! where \f$ \tilde{f} \f$ is now viewed as an \f$ M+1 \times N+1 \f$ matrix. Multiplication on the
!! right by the transpose of the interpolation matrix completes the 2-D interpolation onto the 
!! target nodes
!! \f[
!!      f_{target} = \tilde{f} T^{T} 
!! \f]
!!
!! Because of this implementation, this routine fills in the "interpolationMatrix" attribute with the interpolation
!! matrix, and the "interpolationMatrixTranspose" attribute with the transpose of the interpolation matrix.
!!
!! 
!!  Usage : </H2> 
!! TYPE</B>(Lagrange) :: this 
!!         .... 
!!     CALL</B> this % CalculateInterpolationMatrix(  ) 
!! 
!!   Parameters : </H2>
!!  <table> 
!!     in/out  myPoly  TYPE  
!!             A previously constructed Lagrange data-structure. 
!!             On input</B>, the native interpolation nodes, target interpolation nodes,
!!             and the barycentric weights must be defined. 
!!             On output</B>, the interpolation matrix and its transpose are filled in.
!!  </table>
!!
! ================================================================================================ ! 
!>@}

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

            myPoly % interpolationMatrix(row,col) = 0.0_prec
           
            IF( AlmostEqual( myPoly % targetPoints(row), myPoly % interpolationPoints(col) ) )THEN
               rowHasMatch = .TRUE.
               myPoly % interpolationMatrix(row,col) = 1.0_prec
            ENDIF

         ENDDO 

         IF( .NOT.(rowHasMatch) )THEN 

            temp1 = 0.0_prec

            DO col = 0, myPoly % N        
               temp2 = myPoly % barycentricWeights(col)/( myPoly % targetPoints(row) - myPoly % interpolationPoints(col) )
               myPoly % interpolationMatrix(row,col) = temp2
               temp1 = temp1 + temp2
            ENDDO 

            DO col = 0, myPoly % N 
               myPoly % InterpolationMatrix(row,col) = myPoly % InterpolationMatrix(row,col)/temp1
            ENDDO

         ENDIF 

      ENDDO

      myPoly % interpolationMatrixTranspose = TRANSPOSE( myPoly % interpolationMatrix )

 END SUBROUTINE CalculateInterpolationMatrix_Lagrange

!> \addtogroup Lagrange_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateDerivativeMatrix 
! 
!> \fn CalculateDerivativeMatrix_Lagrange  
!! Calculates and stores the derivative matrix and its transpose to estimate the derivative of a 
!! function ( in both computational directions ) at the native interpolation nodes.
!! 
!! Generates a matrix that can be used to approximate derivatives at the interpolation nodes.
!!
!! Differentiation of the Lagrange interpolating polynomial
!!     \f[
!!            I_N(f) = \frac{\sum_{i=0}^N f_i\frac{w_i}{\xi-\xi_i}}{\sum_{i=0}^N \frac{w_i}{\xi-\xi_i}}
!!     \f]
!! can be rearranged to give (Eq. 3.46 of Kopriva (2009), pg. 80  ) 
!!     \f[
!!            I'_N(f)|_{\xi_j} = \frac{\sum_{i=0}^N f_i\frac{w_i}{\xi-\xi_i}\frac{f_j-f_i}{\xi-\xi_i}}
!!                           {\sum_{i=0}^N \frac{w_i}{\xi-\xi_i}}
!!     \f]
!! when evaluated at each interpolation node
!!
!!   This function is from Alg. 37 on pg. 82 of D.A. Kopriva, 2009.
!!
!! Given nodal values of an interpolant, the derivative can be estimated at the interpolation 
!! nodes using the summations
!!      \f[ \lbrace \frac{\partial f}{\partial \xi^1} \rbrace_{m,n=0}^N = \sum_{i=0}^N( f_{i,n} l'_i(\xi^1_m) )
!!      \f]
!! and
!!      \f[ \lbrace \frac{\partial f}{\partial \xi^2} \rbrace_{m,n=0}^N = \sum_{j=0}^N( f_{m,j} l'_j(\xi^2_n) )
!!      \f]
!!
!! The native interpolation nodes are assumed identical in each direction so that the derivative
!! matrices in each direction are identical. We can write the derivatives as matrix-matrix products
!!      \f[ \lbrace \frac{\partial f}{\partial \xi^1} \rbrace_{m,n=0}^N = D f
!!      \f]
!! and
!!      \f[ \lbrace \frac{\partial f}{\partial \xi^2} \rbrace_{m,n=0}^N = f D^T
!!      \f] 
!! where 
!!      \f[ D_{j,i} = l'_i(\xi_j) 
!!      \f]
!! and \f$ f \f$ is the 2-D array of the native data that is viewed as a matrix.
!!
!! This subroutine calculates the derivative matrix and its transpose and stores them in the
!!  data-structure attributes "D" and "derivativeMatrixTranspose" respectively.
!! 
!! This subroutine depends on 
!!   Module \ref InterpolatioNupportRoutines.f90 : \ref derivativematrix
!!
!!  Usage : </H2> 
!! TYPE</B>(Lagrange) :: this 
!!         .... 
!!     CALL</B> this % CaclulateDerivativeMatrix(  ) 
!! 
!!   Parameters : </H2>
!!  <table> 
!!     in/out  myPoly  TYPE  
!!             A previously constructed Lagrange data-structure. 
!!             On input</B>, the native interpolation nodes
!!             and the barycentric weights must be defined. 
!!             On output</B>, the derivative matrix is filled in.
!!  </table>
! ================================================================================================ ! 
!>@}

  SUBROUTINE CalculateDerivativeMatrix_Lagrange( myPoly )  
    IMPLICIT NONE
    CLASS(Lagrange), INTENT(inout) :: myPoly
    ! Local
    INTEGER    :: row, col

      DO row = 0, myPoly % N
         
        myPoly % derivativeMatrix(row,row) = 0.0_prec

        DO col = 0, myPoly % N
           
          IF( .NOT. (col == row) )THEN

            myPoly % derivativeMatrix(row,col) = myPoly % barycentricWeights(col)/&
                                                 ( myPoly % barycentricWeights(row)*&
                                                   ( myPoly % interpolationPoints(row) - &
                                                     myPoly % interpolationPoints(col) ) )

            myPoly % derivativeMatrix(row,row) = myPoly % derivativeMatrix(row,row) - myPoly % derivativeMatrix(row,col)

          ENDIF
        
        ENDDO 

      ENDDO 
      
      myPoly % derivativeMatrixTranspose = TRANSPOSE( myPoly % derivativeMatrix )

  END SUBROUTINE CalculateDerivativeMatrix_Lagrange

!> \addtogroup Lagrange_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateLagrangePolynomials 
! 
!> \fn CalculateLagrangePolynomials  
!! Evaluates each of the 1-D Lagrange interpolating polynomials at a specified point. 
!! 
!! This function returns the value of each Lagrange interpolating polynomial associated with a 
!! set of interpolation nodes and barycentric weights at a specified point.
!! 
!! From a set of interpolation nodes \f$ \lbrace \xi_j \rbrace_{j=0}^N \f$, the Lagrange 
!! interpolating polynomials are given as
!!     \f[
!!           l_j = \prod_{i=0,\neq j}^N \frac{\xi-\xi_i}{\xi_j-\xi_i}
!!     \f] 
!! This function evaluates the Lagrange interpolating polynomials using the "Barycentric Formulation"
!! (Eq. 3.36 of Kopriva (2009), pg. 74 (with \f$ f_j = \delta_{i,j} \f$))
!!
!!     \f[
!!            l_j = \frac{\frac{w_j}{\xi-\xi_i}}{\sum_{i=0}^N \frac{w_i}{\xi-\xi_i}}
!!     \f] 
!!
!!   This function is from Alg. 34 on pg. 77 of D.A. Kopriva, 2009.
!!
!! The Lagrange interpolating polynomials are given by 
!!   \f[ l_j(\xi) = \prod_{i=0,\neq j}^N \left( \frac{\xi-\xi_i}{\xi_j-\xi_i} \right)  
!!   \f] 
!! where the \f$ \xi_i \f$ are the native nodes (the "s" attribute). Given an input value (\f$\xi\f$),
!! this function returns each of the Lagrange interpolating polynomials \f$ l_j(\xi) \f$ evaluated
!! at this point. This is useful if you have multiple arrays of data that are given at the same
!! native nodes and require interpolation onto a single point.
!!
!!  Usage : </H2> 
!! TYPE</B>(Lagrange) :: this 
!! REAL</B>(prec)        :: lAtS(0:this % N,0:this % N) 
!! REAL</B>(prec)        :: sE, pE 
!!         .... 
!!     lAtS = this % CalculateLagrangePolynomials( sE ) 
!! 
!!   Parameters : </H2>
!!  <table> 
!!     in  myPoly  TYPE(Lagrange)  
!!                     A previously constructed Lagrange structure. 
!!                     The interpolation nodes and barycentric weights are required to produce 
!!                     sensible output.
!!     in  sE  REAL(prec)  
!!                     Location to evaluate the Lagrange interpolating polyomials
!!     out  lAtS(0:myPoly % N)  REAL(prec) 
!!                      Array containing the value of each Lagrange interpolating polynomial 
!!                      at sE. 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}

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

         IF( AlmostEqual(sE, myPoly % interpolationPoints(j)) ) THEN
            lAtS(j) = 1.0_prec
            xMatchesNode = .TRUE.
         ENDIF 

      ENDDO

      IF( xMatchesNode )THEN 
         RETURN
      ENDIF

      temp1 = 0.0_prec
     
      DO j = 0, myPoly % N 
         temp2 = myPoly % barycentricWeights(j)/(sE - myPoly % interpolationPoints(j))
         lAtS(j) = temp2
         temp1 = temp1 + temp2
      ENDDO 
  
      lAtS = lAtS/temp1 

  END FUNCTION CalculateLagrangePolynomials

  FUNCTION ApplyInterpolationMatrix_1D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( fNew )  
    IMPLICIT NONE
    TYPE(Lagrange) :: myPoly
    INTEGER         :: nVariables, nElements
    REAL(prec)      :: f(0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec)      :: fNew(0:myPoly % M, 1:nVariables, 1:nElements)
    ! Local
    INTEGER :: iVar, iEl, i, a 

      DO iEl = 1, nElements
        DO iVar = 1, nVariables

          DO a = 0, myPoly % M

            fNew(a,iVar,iEl) = 0.0_prec

            DO i = 0, myPoly % N

              fNew(a,iVar,iEl) = fNew(a,iVar,iEl) + myPoly % interpolationMatrixTranspose(i,a)*f(i,iVar,iEl)

            ENDDO

          ENDDO

        ENDDO
      ENDDO

  END FUNCTION ApplyInterpolationMatrix_1D_Lagrange

  FUNCTION ApplyInterpolationMatrix_2D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( fNew )  
    IMPLICIT NONE
    TYPE(Lagrange) :: myPoly
    INTEGER         :: nVariables, nElements
    REAL(prec)      :: f(0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec)      :: fNew(0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
    ! Local
    INTEGER :: i, j, a, b, p, iEl, iVar
    REAL(prec) :: fa
   
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
         
          DO b = 0, myPoly % M
            DO a = 0, myPoly % M
             
              fNew(a,b,iVar,iEl) = 0.0_prec

              DO j = 0, myPoly % N
                   
                fa = 0.0_prec
                DO i = 0, myPoly % N
                  fa = fa + f(i,j,iVar,iEl)*myPoly % interpolationMatrixTranspose(i,a)
                ENDDO
                      
                fNew(a,b,iVar,iEl) = fNew(a,b,iVar,iEl) + fa*myPoly % interpolationMatrixTranspose(j,b)

              ENDDO
                   
            ENDDO
          ENDDO

        ENDDO   
      ENDDO

  END FUNCTION ApplyInterpolationMatrix_2D_Lagrange

  FUNCTION ApplyInterpolationMatrix_3D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( fNew )  
    IMPLICIT NONE
    TYPE(Lagrange) :: myPoly
    INTEGER         :: nElements, nVariables
    REAL(prec)      :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N, 1:nVariables, 1:nElements)
    REAL(prec)      :: fNew(0:myPoly % M, 0:myPoly % M, 0:myPoly % M, 1:nVariables, 1:nElements)
    ! Local
    INTEGER :: i, j, k, a, b, c, iEl, iVar
    REAL(prec) :: fa, fab
   
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
         
         
          DO c = 0, myPoly % M
            DO b = 0, myPoly % M
              DO a = 0, myPoly % M
               
                fNew(a,b,c,iVar,iEl) = 0.0_prec

                DO k = 0, myPoly % N
                  
                  fab = 0.0_prec
                  DO j = 0, myPoly % N
                     
                    fa = 0.0_prec
                    DO i = 0, myPoly % N
                      fa = fa + f(i,j,k,iVar,iEl)*myPoly % interpolationMatrixTranspose(i,a)
                    ENDDO
                        
                    fab = fab + fa*myPoly % interpolationMatrixTranspose(j,b)

                  ENDDO
                     
                  fNew(a,b,c,iVar,iEl) = fNew(a,b,c,iVar,iEl) + fab*myPoly % interpolationMatrixTranspose(k,c)

                ENDDO
                  

              ENDDO
            ENDDO
          ENDDO

        ENDDO   
      ENDDO

  END FUNCTION ApplyInterpolationMatrix_3D_Lagrange

  FUNCTION CalculateDerivative_1D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( derF )  
    IMPLICIT NONE
    TYPE(Lagrange) :: myPoly
    INTEGER         :: nVariables, nElements
    REAL(prec)      :: f(0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec)      :: derF(0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER :: i, ii, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
 
          DO i = 0, myPoly % N
    
            derF(i,iVar,iEl) = 0.0_prec 
            DO ii = 0, myPoly % N
              derF(i,iVar,iEl) = derF(i,iVar,iEl) + myPoly % derivativeMatrixTranspose(ii,i)*f(ii,iVar,iEl)
            ENDDO
    
          ENDDO

        ENDDO
      ENDDO

  END FUNCTION CalculateDerivative_1D_Lagrange
!
  FUNCTION CalculateGradient_2D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( gradF )  

    IMPLICIT NONE
    TYPE(Lagrange) :: myPoly
    INTEGER         :: nVariables, nElements
    REAL(prec)      :: f(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec)      :: gradF(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER    :: i, j, ii, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % N
            DO i = 0, myPoly % N
    
              gradF(1,i,j,iVar,iEl) = 0.0_prec 
              gradF(2,i,j,iVar,iEl) = 0.0_prec 
              DO ii = 0, myPoly % N
                gradF(1,i,j,iVar,iEl) = gradF(1,i,j,iVar,iEl) + myPoly % derivativeMatrixTranspose(ii,i)*f(ii,j,iVar,iEl)
                gradF(2,i,j,iVar,iEl) = gradF(2,i,j,iVar,iEl) + myPoly % derivativeMatrixTranspose(ii,j)*f(i,ii,iVar,iEl)
              ENDDO
    
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END FUNCTION CalculateGradient_2D_Lagrange
!
  FUNCTION CalculateDivergence_2D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( divF )  
    IMPLICIT NONE
    TYPE( Lagrange ) :: myPoly
    INTEGER          :: nVariables, nElements
    REAL(prec)       :: f(1:2,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec)       :: divF(0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER :: i, j, iEl, iVar, ii
  
      
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myPoly % N
            DO i = 0, myPoly % N
               
              divf(i,j,iVar,iEl) = 0.0_prec

              DO ii = 0, myPoly % N  
                divf(i,j,iVar,iEl) = divf(i,j,iVar,iEl) + myPoly % derivativeMatrixTranspose(ii,i)*f(1,ii,j,iVar,iEl) + &
                                                            myPoly % derivativeMatrixTranspose(ii,j)*f(2,i,ii,iVar,iEl)
              ENDDO

            ENDDO
          ENDDO
        ENDDO
      ENDDO
      
  END FUNCTION CalculateDivergence_2D_Lagrange
!
  FUNCTION CalculateGradient_3D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( gradF )  
    IMPLICIT NONE
    TYPE( Lagrange ) :: myPoly
    INTEGER          :: nVariables, nElements
    REAL(prec)       :: f(0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec)       :: gradF(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER :: i, j, k, iEl, iVar, ii
  
      
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, myPoly % N
            DO j = 0, myPoly % N
              DO i = 0, myPoly % N
               
                gradf(1,i,j,k,iVar,iEl) = 0.0_prec
                gradf(2,i,j,k,iVar,iEl) = 0.0_prec
                gradf(3,i,j,k,iVar,iEl) = 0.0_prec

                DO ii = 0, myPoly % N  
                  gradf(1,i,j,k,iVar,iEl) = gradf(1,i,j,k,iVar,iEl) + myPoly % derivativeMatrixTranspose(ii,i)*f(ii,j,k,iVar,iEl)
                  gradf(2,i,j,k,iVar,iEl) = gradf(2,i,j,k,iVar,iEl) + myPoly % derivativeMatrixTranspose(ii,j)*f(i,ii,k,iVar,iEl)
                  gradf(3,i,j,k,iVar,iEl) = gradf(3,i,j,k,iVar,iEl) + myPoly % derivativeMatrixTranspose(ii,k)*f(i,j,ii,iVar,iEl)
                ENDDO

              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      
  END FUNCTION CalculateGradient_3D_Lagrange
!
  FUNCTION CalculateDivergence_3D_Lagrange( myPoly, f, nVariables, nElements ) RESULT( divF )  
    IMPLICIT NONE
    TYPE( Lagrange ) :: myPoly
    INTEGER          :: nVariables, nElements
    REAL(prec)       :: f(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    REAL(prec)       :: divF(0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nVariables,1:nElements)
    ! Local
    INTEGER :: i, j, k, iEl, iVar, ii
  
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, myPoly % N
            DO j = 0, myPoly % N
              DO i = 0, myPoly % N
               
                divf(i,j,k,iVar,iEl) = 0.0_prec

                DO ii = 0, myPoly % N  
                  divf(i,j,k,iVar,iEl) = divf(i,j,k,iVar,iEl) + &
                                         myPoly % derivativeMatrixTranspose(ii,i)*f(1,ii,j,k,iVar,iEl) + &
                                         myPoly % derivativeMatrixTranspose(ii,j)*f(2,i,ii,k,iVar,iEl) + &
                                         myPoly % derivativeMatrixTranspose(ii,k)*f(3,i,j,ii,iVar,iEl)
                ENDDO

              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      
  END FUNCTION CalculateDivergence_3D_Lagrange

#ifdef HAVE_CUDA

  ATTRIBUTES(Global) SUBROUTINE ApplyInterpolationMatrix_1D_CUDAKernel( IntMatT, f, fNew, N, M, nVariables, nElems  )
    ! Currently only works for M,N <= 7
    IMPLICIT NONE
    INTEGER, INTENT(in)             :: N, M, nVariables, nElems
    REAL(prec), DEVICE, INTENT(in)  :: IntMatT(0:N,0:M)
    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,1:nVariables,1:nElems)
    REAL(prec), DEVICE, INTENT(out) :: fnew(0:M,1:nVariables,1:nElems)
    ! Local
    INTEGER            :: i, a, iEl, iVar
    REAL(prec), SHARED :: floc(0:7)
    REAL(prec)         :: fm
   
      iVar = blockIdx % x
      iEl  = blockIdx % y
      a = threadIdx % x-1 + (blockIDx % y-1)*blockDim % x
      
      IF( a <= M )THEN

        ! Pre-fetch data
        IF( a <= N )THEN
          floc(a) = f(a,iVar,iEl)
        ENDIF
        
        CALL syncthreads( )
      
        fm = 0.0_prec
        DO i = 0, N
          fm = fm + floc(i)*IntMatT(i,a)
        ENDDO
               
        fnew(a,iVar,iEl) = fm
         
      ENDIF

  END SUBROUTINE ApplyInterpolationMatrix_1D_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE ApplyInterpolationMatrix_2D_CUDAKernel( IntMatT, f, fNew, N, M, nVariables, nElems  )
    ! Currently only works for M,N <= 7
    IMPLICIT NONE
    INTEGER, INTENT(in)             :: N, M, nVariables, nElems
    REAL(prec), DEVICE, INTENT(in)  :: IntMatT(0:N,0:M)
    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,0:N,1:nVariables,1:nElems)
    REAL(prec), DEVICE, INTENT(out) :: fnew(0:M,0:M,1:nVariables,1:nElems)
    ! Local
    INTEGER            :: i, j, a, b, iEl, iVar
    REAL(prec), SHARED :: floc(0:7,0:7)
    REAL(prec)         :: fm, fmn
   
      iVar = blockIdx % x
      iEl  = blockIdx % y
      a = threadIdx % x-1 + (blockIDx % y-1)*blockDim % x
      b = threadIdx % y-1 + (blockIDx % y-1)*blockDim % y
      
      IF( a <= M .AND. b <= M )THEN

        ! Pre-fetch data
        IF( a <= N .AND. b <= N )THEN
          floc(a,b) = f(a,b,iVar,iEl)
        ENDIF
        
        CALL syncthreads( )
               
        fmn = 0.0_prec
        DO j = 0, N
            
          fm = 0.0_prec
          DO i = 0, N
            fm = fm + floc(i,j)*IntMatT(i,a)
          ENDDO
               
          fmn = fmn + fm*IntMatT(j,b)

        ENDDO
            
        fnew(a,b,iVar,iEl) = fmn
         
      ENDIF

  END SUBROUTINE ApplyInterpolationMatrix_2D_CUDAKernel

  ATTRIBUTES(Global) SUBROUTINE ApplyInterpolationMatrix_3D_CUDAKernel( IntMatT, f, fNew, N, M, nVariables, nElems  )
    ! Currently only works for M,N <= 7
    IMPLICIT NONE
    INTEGER, INTENT(in)             :: N, M, nVariables, nElems
    REAL(prec), DEVICE, INTENT(in)  :: IntMatT(0:N,0:M)
    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,0:N,0:N,1:nVariables,1:nElems)
    REAL(prec), DEVICE, INTENT(out) :: fnew(0:M,0:M,0:M,1:nVariables,1:nElems)
    ! Local
    INTEGER            :: i, j, k, a, b, c, iEl, iVar
    REAL(prec), SHARED :: floc(0:7,0:7,0:7)
    REAL(prec)         :: fm, fmn, fmnp
   
   
      iVar = blockIdx % x
      iEl  = blockIdx % y
      a = threadIdx % x-1 + (blockIDx % y-1)*blockDim % x
      b = threadIdx % y-1 + (blockIDx % y-1)*blockDim % y
      c = threadIdx % z-1 + (blockIDx % y-1)*blockDim % z
      
      IF( a <= M .AND. b <= M .AND. c <= M )THEN

        ! Pre-fetch data
        IF( a <= N .AND. b <= N .AND. c <= N )THEN
          floc(a,b,c) = f(a,b,c,iVar,iEl)
        ENDIF
        
        CALL syncthreads( )
      
              
        fmnp = 0.0_prec
        DO k = 0, N
         
          fmn = 0.0_prec
          DO j = 0, N
          
            fm = 0.0_prec
            DO i = 0, N
              fm = fm + floc(i,j,k)*IntMatT(i,a)
            ENDDO
             
            fmn = fmn + fm*IntMatT(j,b)
          ENDDO
           
          fmnp = fmnp + fmn*IntMatT(k,c)

        ENDDO
         
         fnew(a,b,c,iVar,iEl) = fmnp
         
      ENDIF

  END SUBROUTINE ApplyInterpolationMatrix_3D_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE CalculateDerivative_1D_CUDAKernel( DMatT, f, derf, N, nVariables, nElems )
    ! Currently only works for M,N <= 7
    IMPLICIT NONE
    INTEGER, INTENT(in)             :: N, nVariables, nElems
    REAL(prec), DEVICE, INTENT(in)  :: DMatT(0:N,0:N)
    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,1:nVariables,1:nElems)
    REAL(prec), DEVICE, INTENT(out) :: derf(0:N,1:nVariables,1:nElems)
    ! Local
    INTEGER            :: i, j, iEl, iVar, ii
    REAL(prec), SHARED :: floc(0:7)
    REAL(prec)         :: df
  
  
      iVar = blockIdx % x
      iEl  = blockIdx % y
      i    = threadIdx % x-1
      
      IF( i <= N )THEN
      
         ! Pre-fetch into shared memory for this block
         floc(i) = f(i,iVar,iEl)

         CALL syncthreads( )

         df = 0.0_prec

         DO ii = 0, N
            df = df + DMatT(ii,i)*floc(ii)
         ENDDO
                  
         derf(i,iVar,iEl) = df

      ENDIF
    
  END SUBROUTINE CalculateDerivative_1D_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE CalculateGradient_2D_CUDAKernel( DMatT, f, gradf, N, nVariables, nElems )
    ! Currently only works for M,N <= 7
    IMPLICIT NONE
    INTEGER, INTENT(in)             :: N, nVariables, nElems
    REAL(prec), DEVICE, INTENT(in)  :: DMatT(0:N,0:N)
    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,0:N,1:nVariables,1:nElems)
    REAL(prec), DEVICE, INTENT(out) :: gradf(1:2,0:N,0:N,1:nVariables,1:nElems)
    ! Local
    INTEGER            :: i, j, iEl, iVar, ii
    REAL(prec), SHARED :: floc(0:7,0:7)
    REAL(prec)         :: df(1:2)
  
  
      iVar = blockIdx % x
      iEl  = blockIdx % y
      i    = threadIdx % x-1
      j    = threadIdx % y-1
      
      IF( i <= N .AND. j <= N )THEN
      
         ! Pre-fetch into shared memory for this block
         floc(i,j) = f(i,j,iVar,iEl)

         CALL syncthreads( )

         df(1) = 0.0_prec
         df(2) = 0.0_prec

         DO ii = 0, N
            df(1) = df(1) + DMatT(ii,i)*floc(ii,j)
            df(2) = df(2) + DMatT(ii,j)*floc(i,ii)
         ENDDO
                  
         gradf(1,i,j,iVar,iEl) = df(1)
         gradf(2,i,j,iVar,iEl) = df(2)

      ENDIF
    
  END SUBROUTINE CalculateGradient_2D_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE CalculateDivergence_2D_CUDAKernel( DMatT, f, divf, N, nVariables, nElems )
    ! Currently only works for M,N <= 7
    IMPLICIT NONE
    INTEGER, INTENT(in)             :: N, nVariables, nElems
    REAL(prec), DEVICE, INTENT(in)  :: DMatT(0:N,0:N)
    REAL(prec), DEVICE, INTENT(in)  :: f(1:2,0:N,0:N,1:nVariables,1:nElems)
    REAL(prec), DEVICE, INTENT(out) :: divf(0:N,0:N,1:nVariables,1:nElems)
    ! Local
    INTEGER            :: i, j, k, iEl, iVar, ii
    REAL(prec), SHARED :: floc(1:2,0:7,0:7)
    REAL(prec) :: df
  
  
      iVar = blockIdx % x
      iEl  = blockIdx % y
      i    = threadIdx % x-1
      j    = threadIdx % y-1
      
      IF( i <= N .AND. j <= N )THEN
      
         ! Pre-fetch into shared memory for this block
         floc(1,i,j) = f(1,i,j,iVar,iEl)
         floc(2,i,j) = f(2,i,j,iVar,iEl)

         CALL syncthreads( )

         df = 0.0_prec
         DO ii = 0, N
            df = df + DMatT(ii,i)*floc(1,ii,j) + &
                      DMatT(ii,j)*floc(2,i,ii)
         ENDDO
                  
         divf(i,j,iVar,iEl) = df

      ENDIF
    
  END SUBROUTINE CalculateDivergence_2D_CUDAKernel 
!
  ATTRIBUTES(Global) SUBROUTINE CalculateGradient_3D_CUDAKernel( DMatT, f, gradf, N, nVariables, nElems )
    ! Currently only works for M,N <= 7
    IMPLICIT NONE
    INTEGER, INTENT(in)             :: N, nVariables, nElems
    REAL(prec), DEVICE, INTENT(in)  :: DMatT(0:N,0:N)
    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,0:N,0:N,1:nVariables,1:nElems)
    REAL(prec), DEVICE, INTENT(out) :: gradf(1:3,0:N,0:N,0:N,1:nVariables,1:nElems)
    ! Local
    INTEGER            :: i, j, k, iEl, iVar, ii
    REAL(prec), SHARED :: floc(0:7,0:7,0:7)
    REAL(prec)         :: df(1:3)
  
  
      iVar = blockIdx % x
      iEl  = blockIdx % y
      i    = threadIdx % x-1
      j    = threadIdx % y-1
      k    = threadIdx % z-1
      
      IF( i <= N .AND. j <= N .AND. k <= N )THEN
      
         ! Pre-fetch into shared memory for this block
         floc(i,j,k) = f(i,j,k,iVar,iEl)

         CALL syncthreads( )

         df(1) = 0.0_prec
         df(2) = 0.0_prec
         df(3) = 0.0_prec

         DO ii = 0, N
            df(1) = df(1) + DMatT(ii,i)*floc(ii,j,k)
            df(2) = df(2) + DMatT(ii,j)*floc(i,ii,k)
            df(3) = df(3) + DMatT(ii,k)*floc(i,j,ii)
         ENDDO
                  
         gradf(1,i,j,k,iVar,iEl) = df(1)
         gradf(2,i,j,k,iVar,iEl) = df(2)
         gradf(3,i,j,k,iVar,iEl) = df(3)

      ENDIF
    
  END SUBROUTINE CalculateGradient_3D_CUDAKernel
!
  ATTRIBUTES(Global) SUBROUTINE CalculateDivergence_3D_CUDAKernel( DMatT, f, divf, N, nVariables, nElems )
    ! Currently only works for M,N <= 7
    IMPLICIT NONE
    INTEGER, INTENT(in)             :: N, nVariables, nElems
    REAL(prec), DEVICE, INTENT(in)  :: DMatT(0:N,0:N)
    REAL(prec), DEVICE, INTENT(in)  :: f(1:3,0:N,0:N,0:N,1:nVariables,1:nElems)
    REAL(prec), DEVICE, INTENT(out) :: divf(0:N,0:N,0:N,1:nVariables,1:nElems)
    ! Local
    INTEGER            :: i, j, k, iEl, iVar, ii
    REAL(prec), SHARED :: floc(1:3,0:7,0:7,0:7)
    REAL(prec) :: df
  
  
      iVar = blockIdx % x
      iEl  = blockIdx % y
      i    = threadIdx % x-1
      j    = threadIdx % y-1
      k    = threadIdx % z-1
      
      IF( i <= N .AND. j <= N .AND. k <= N )THEN
      
         ! Pre-fetch into shared memory for this block
         floc(1,i,j,k) = f(1,i,j,k,iVar,iEl)
         floc(2,i,j,k) = f(2,i,j,k,iVar,iEl)
         floc(3,i,j,k) = f(3,i,j,k,iVar,iEl)
      ENDIF
      
         CALL syncthreads( )
 
      IF( i <= N .AND. j <= N .AND. k <= N )THEN
      
         df = 0.0_prec
         DO ii = 0, N
            df = df + DMatT(ii,i)*floc(1,ii,j,k) + &
                      DMatT(ii,j)*floc(2,i,ii,k) + &
                      DMatT(ii,k)*floc(3,i,j,ii)
         ENDDO
                  
         divf(i,j,k,iVar,iEl) = df

      ENDIF
    
 END SUBROUTINE CalculateDivergence_3D_CUDAKernel
 
#endif



END MODULE Lagrange_Class
