! NodalDG_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE NodalDG_Class

USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
USE Quadrature
USE Lagrange_Class

IMPLICIT NONE

! ================================================================================================ !
!
!  The NodalDG class contains attributes and procedure needed for performing interpolation, 
!  differentiation, and integration for spectral element methods
!  
!  An interpolant is formed that handles mapping between a computational "quadrature" mesh and a
!  uniform plotting mesh. Quadrature (integration) weights are stored for use with Galerkin type
!  methods. Galerkin derivative matrices are stored to facilitate the computation of weak derivatives.
!  Finally, an interpolation matrix and accompanying subroutine is provided to interpolate 3-D data, 
!  defined on the quadrature mesh, to the element boundaries.
!
! ================================================================================================ !

  TYPE NodalDG

    INTEGER                 :: N     
    INTEGER                 :: nTargetPoints
    TYPE(Lagrange)          :: interp 
    REAL(prec), ALLOCATABLE :: quadratureWeights(:)
    REAL(prec), ALLOCATABLE :: dgDerivativeMatrix(:,:)
    REAL(prec), ALLOCATABLE :: dgDerivativeMatrixTranspose(:,:)
    REAL(prec), ALLOCATABLE :: boundaryInterpolationMatrix(:,:)
      
#ifdef HAVE_CUDA
    INTEGER, ALLOCATABLE, DEVICE    :: N_dev  
    INTEGER, ALLOCATABLE, DEVICE    :: nTargetPoints_dev
    REAL(prec), DEVICE, ALLOCATABLE :: quadratureWeights_dev(:)
    REAL(prec), DEVICE, ALLOCATABLE :: dgDerivativeMatrix_dev(:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: dgDerivativeMatrixTranspose_dev(:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: boundaryInterpolationMatrix_dev(:,:)
#endif
    CONTAINS

      PROCEDURE :: Build => Build_NodalDG
      PROCEDURE :: Trash => Trash_NodalDG

      PROCEDURE :: CalculateFunctionsAtBoundaries_3D
      PROCEDURE :: DG_Divergence_3D
      PROCEDURE :: DG_Gradient_3D
      
    END TYPE NodalDG
    
 CONTAINS

! ================================================================================================ !
!
! Build_NodalDG 
!
!   Allocates space fills values for the NodalDG attributes using to the specified quadrature and 
!   approximation form.
! 
!   Usage :
!
!     TYPE(NodalDG) :: dgStorage
!     INTEGER       :: N, nTargetPoints, quadrature
!     REAL(prec)    :: targetPoints(0:nTargetPoints)
!
!       CALL this % Build( targetPoints N, nTargetPoints, quadrature ) 
! 
!   Input/Output :
!   
!      dgStorage (out) 
!        On output, the memory has been allocated for the attributes of the NodalDG structure, and
!        the attributes are set according to the targetPoints and quadrature specified. If CUDA
!        is enabled, device memory is allocated for the device attributes and device attributes are
!        copied from the host-side arrays. 
!       
!      targetPoints(0:nTargetPoints) (in)
!        The target interpolation points to aid in mapping from the quadrature mesh to a uniform
!        plotting mesh.
!
!      N  (in)
!        Polynomial degree of the DG Spectral Element Method
!
!      nTargetPoints (in)
!        The upper bound of the targetPoints array.
!
!      quadrature (in)
!        A flag for specifying the desired type of quadrature. Can be set to either GAUSS or
!        GAUSS_LOBATTO.
!   
! ================================================================================================ ! 

  SUBROUTINE Build_NodalDG( myNodal, targetPoints, N, nTargetPoints, quadrature  )
    IMPLICIT NONE
    CLASS(NodalDG), INTENT(out) :: myNodal
    INTEGER, INTENT(in)         :: N, nTargetPoints
    REAL(prec), INTENT(in)      :: targetPoints(0:nTargetPoints)
    INTEGER, INTENT(in)         :: quadrature
    ! Local
    INTEGER   :: i, j
    REAL(prec):: quadraturePoints(0:N), quadratureWeights(0:N)

      myNodal % N             = N
      myNodal % nTargetPoints = nTargetPoints
      
      ! Allocate space
      ALLOCATE( myNodal % dgDerivativeMatrix(0:N,0:N), &
                myNodal % dgDerivativeMatrixTranspose(0:N,0:N), &
                myNodal % boundaryInterpolationMatrix(0:N,0:1), &
                myNodal % quadratureWeights(0:N) )
                
      myNodal % dgDerivativeMatrix          = 0.0_prec
      myNodal % dgDerivativeMatrixTranspose = 0.0_prec
      myNodal % boundaryInterpolationMatrix = 0.0_prec
      myNodal % quadratureWeights           = 0.0_prec

      CALL LegendreQuadrature( N, quadraturePoints, quadratureWeights, quadrature )
      myNodal % quadratureWeights = quadratureWeights    
       
      CALL myNodal % interp % Build( N, nTargetPoints, quadraturePoints, targetPoints )
   
      myNodal % boundaryInterpolationMatrix(0:N,0) = myNodal % interp % CalculateLagrangePolynomials( -1.0_prec)
      myNodal % boundaryInterpolationMatrix(0:N,1) = myNodal % interp % CalculateLagrangePolynomials( 1.0_prec )

      
      DO j = 0, N 
        DO i = 0, N

          myNodal % dgDerivativeMatrix(i,j) = -myNodal % interp % derivativeMatrix(j,i)*&
                                               myNodal % quadratureWeights(j)/&
                                               myNodal % quadratureWeights(i)

        ENDDO
      ENDDO
         
      DO j = 0, N
        DO i = 0, N

          myNodal % dgDerivativeMatrixTranspose(j,i) = myNodal % dgDerivativeMatrix(i,j)

        ENDDO
      ENDDO
 
#ifdef HAVE_CUDA
      ALLOCATE( myNodal % N_dev, myNodal % nTargetPoints_dev )
      myNodal % N_dev             = N
      myNodal % nTargetPoints_dev = nTargetPoints
      
      
      ALLOCATE( myNodal % dgDerivativeMatrix_dev(0:N,0:N), &
                myNodal % dgDerivativeMatrixTranspose_dev(0:N,0:N), &
                myNodal % boundaryInterpolationMatrix_dev(0:N,0:1), &
                myNodal % quadratureWeights_dev(0:N) )
                
      myNodal % dgDerivativeMatrix_dev          = myNodal % dgDerivativeMatrix
      myNodal % dgDerivativeMatrixTranspose_dev = myNodal % dgDerivativeMatrixTranspose
      myNodal % boundaryInterpolationMatrix_dev = myNodal % boundaryInterpolationMatrix
      myNodal % quadratureWeights_dev           = myNodal % quadratureWeights
#endif

  END SUBROUTINE Build_NodalDG
  
! ================================================================================================ ! 
!
! Trash_NodalDG  
!   Frees memory held by the attributes of the NodalDG class. 
! 
!   Usage :
!
!     TYPE(NodalDG) :: dgStorage
!         
!       CALL dgStorage % Trash( ) 
! 
!    Input/Output:
!   
!      dgStorage (in/out)
!        On input, a NodalDG class that has previously been constructed. 
!        On output, the memory held by the attributes of this data-structure is deallocated. If CUDA
!        is enabled, device memory associated with device attributes is deallocated
!                                                           
! ================================================================================================ ! 

  SUBROUTINE Trash_NodalDG( myNodal)
    IMPLICIT NONE
    CLASS( NodalDG ), INTENT(inout) :: myNodal

      CALL myNodal % interp % Trash( )
      DEALLOCATE( myNodal % dgDerivativeMatrix, &
                  myNodal % dgDerivativeMatrixTranspose, &
                  myNodal % boundaryInterpolationMatrix, &
                  myNodal % quadratureWeights )
                  
#ifdef HAVE_CUDA
      DEALLOCATE( myNodal % N_dev, myNodal % nTargetPoints_dev )
      DEALLOCATE( myNodal % dgDerivativeMatrix_dev, &
                  myNodal % dgDerivativeMatrixTranspose_dev, &
                  myNodal % boundaryInterpolationMatrix_dev, &
                  myNodal % quadratureWeights_dev )
#endif


  END SUBROUTINE Trash_NodalDG

! ================================================================================================ ! 
! 
! CalculateFunctionsAtBoundaries_3D
!   Interpolates a set of functions, given on each element's quadrature mesh, to the element faces. 
! 
!   Usage :
!
!     TYPE(NodalDG) :: dgStorage
!     REAL(prec)    :: f(0:N,0:N,0:N,1:nVariables,1:nElements)
!     REAL(prec)    :: fAtBoundaries(0:N,0:N,1:nVariables,1:6,1:nElements)
!     INTEGER       :: nVariables, nElements
!     
!         
!       CALL dgStorage % CalculateFunctionsAtBoundaries( f, fAtBoundaries, nVariables, nElements ) 
! 
!    * Note that N refers to polynomial degree under the NodalDG class ( dgStorage % N )
!
!   Input/Output:
!   
!     If CUDA is enabled, all of the input and output must be device scalars and arrays
!
!     dgStorage (in)
!     
!     f(0:N,0:N,0:N,1:nVariables,1:nElements) (in)
!       Array of nodal function values. The first three dimensions are over the quadrature points
!       in 3-D, the fourth dimension is over the number of variables, and the last dimension
!       cycles over the elements in a spectral element mesh.
!
!     fAtBoundaries(0:N,0:N,1:nVariables,1:6,1:nElements) (out)
!       Array of nodal function values on the element faces. The first two dimensions are over 
!       the quadrature points on the element faces, the third dimension is over the number of
!       variables, the fourth dimension is over the faces of the element (south, east, north, west,
!       bottom, top), and the last dimension cycles over the elements in a spectral element mesh.
!      
!     nVariables (in)
!       The number of variables
!
!     nElements (in)
!       The number of elements in the spectral element mesh.
!       
! ================================================================================================ ! 

  SUBROUTINE CalculateFunctionsAtBoundaries_3D( myNodal, f, fAtBoundaries, nVariables, nElements )
    IMPLICIT NONE
    CLASS( NodalDG ), INTENT(in) :: myNodal 
#ifdef HAVE_CUDA
    INTEGER, DEVICE, INTENT(in)     :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(0:myNodal % N,0:myNodal % N,0:myNodal % N,1:nVariables,nElements)
    REAL(prec), DEVICE, INTENT(out) :: fAtBoundaries(0:myNodal % N,0:myNodal % N,1:nVariables,1:6,1:nElements)
    ! Local
    TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3(4*(ceiling( REAL(myNodal % N+1)/4 ) ), &
                    4*(ceiling( REAL(myNodal % N+1)/4 ) ) , &
                    nVariables )
      grid = dim3(nElements, 1, 1)  
      
      CALL CalculateFunctionsAtBoundaries_3D_CUDAKernel<<<grid, tBlock>>>( f, fAtBoundaries, &
                                                                           myNodal % boundaryInterpolationMatrix_dev, &
                                                                           myNodal % N_dev, nVariables, nElements )
      
#else
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myNodal % N,0:myNodal % N,0:myNodal % N,1:nVariables,nElements)
    REAL(prec), INTENT(out) :: fAtBoundaries(0:myNodal % N,0:myNodal % N,1:nVariables,1:6,1:nElements)
     
      fAtBoundaries = CalculateFunctionsAtBoundaries_3D_NodalDG( myNodal, f, nVariables, nElements )

#endif    

    
    
  END SUBROUTINE CalculateFunctionsAtBoundaries_3D
  
! ================================================================================================ ! 
! 
! DG_Divergence_3D
!   Calculates the weak form of the divergence of a vector function.
! 
!   Usage :
!
!     TYPE(NodalDG) :: dgStorage
!     REAL(prec)    :: f(1:3,0:N,0:N,0:N,1:nVariables,1:nElements)
!     REAL(prec)    :: fnAtBoundaries(0:N,0:N,1:nVariables,1:6,1:nElements)
!     REAL(prec)    :: divF(0:N,0:N,0:N,1:nVariables,1:nElements)
!     INTEGER       :: nVariables, nElements
!     
!         
!       CALL dgStorage % DG_Divergence_3D( f, fnAtBoundaries, divF, nVariables, nElements ) 
! 
!    * Note that N refers to polynomial degree under the NodalDG class ( dgStorage % N )
!
!   Input/Output :
!
!     If CUDA is enabled, all of the input and output must be device scalars and arrays
!
!     dgStorage (in)
!     
!     f(1:3,0:N,0:N,0:N,1:nVariables,1:nElements) (in)
!       Array of nodal (3-D)vector-function values. The first dimension is over the three spatial 
!       directions of the vector, the next three dimensions are over the quadrature points
!       in 3-D, the fifth dimension is over the number of variables, and the last dimension
!       cycles over the elements in a spectral element mesh.
!
!     fnAtBoundaries(0:N,0:N,1:nVariables,1:6,1:nElements) (in)
!       Array of normal flux values on the element faces. The first two dimensions are over 
!       the quadrature points on the element faces, the third dimension is over the number of
!       variables, the fourth dimension is over the faces of the element (south, east, north, west,
!       bottom, top), and the last dimension cycles over the elements in a spectral element mesh.
!      
!     divF(0:N,0:N,0:N,1:nVariables,1:nElements) (out)
!       The divergence of the vector functions over all of the elements
!
!     nVariables (in)
!       The number of variables
!
!     nElements (in)
!       The number of elements in the spectral element mesh.
!
! ================================================================================================ ! 

  SUBROUTINE DG_Divergence_3D( myNodal, f, fnAtBoundaries, divF, nVariables, nElements )
    IMPLICIT NONE
    CLASS( NodalDG ) :: myNodal
#ifdef HAVE_CUDA
    INTEGER, DEVICE, INTENT(in)     :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(1:3,0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: fnAtBoundaries(0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:6, 1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: divF(0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    ! Local
    TYPE(dim3) :: grid, tBlock
    INTEGER    :: threadCount
    
      threadCount = MIN( 4*(ceiling( REAL(myNodal % N+1)/4 ) ), 8 )

      tBlock = dim3( threadCount, threadCount, threadCount )
      grid = dim3( nVariables, nElements, 1) 
      
      CALL DG_Divergence_3D_CUDAKernel<<<grid, tBlock>>>( f, fnAtBoundaries, divF, &
                                                          myNodal % boundaryInterpolationMatrix_dev, &
                                                          myNodal % dgDerivativeMatrixTranspose_dev, &
                                                          myNodal % quadratureWeights_dev, &
                                                          myNodal % N_dev, nVariables, nElements )

#else
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(1:3,0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(in)  :: fnAtBoundaries(0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:6, 1:nElements)
    REAL(prec), INTENT(out) :: divf(0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)


      divF = DG_Divergence_3D_NodalDG( myNodal, f, fnAtBoundaries, nVariables, nElements )
    
#endif

  END SUBROUTINE DG_Divergence_3D

! ================================================================================================ ! 
! 
! DG_Gradient_3D
!   Calculates the weak form of the gradient of function.
! 
!   Usage :
!
!     TYPE(NodalDG) :: dgStorage
!     REAL(prec)    :: f(0:N,0:N,0:N,1:nVariables,1:nElements)
!     REAL(prec)    :: fAtBoundaries(0:N,0:N,1:nVariables,1:6,1:nElements)
!     REAL(prec)    :: gradF(1:3,0:N,0:N,0:N,1:nVariables,1:nElements)
!     INTEGER       :: nVariables, nElements
!     
!         
!       CALL dgStorage % DG_Gradient_3D( f, fAtBoundaries, gradF, nVariables, nElements ) 
! 
!    * Note that N refers to polynomial degree under the NodalDG class ( dgStorage % N )
!
!   Input/Output :
!
!     If CUDA is enabled, all of the input and output must be device scalars and arrays
!
!     dgStorage (in)
!     
!     f(0:N,0:N,0:N,1:nVariables,1:nElements) (in)
!       Array of nodalfunction values. The first three dimensions are over the quadrature points
!       in 3-D, the fifth dimension is over the number of variables, and the last dimension
!       cycles over the elements in a spectral element mesh.
!
!     fAtBoundaries(0:N,0:N,1:nVariables,1:6,1:nElements) (in)
!       Array of function values on the element faces (weighted by the face normal). The first two 
!       dimensions are over the quadrature points on the element faces, the third dimension is over 
!       the number of variables, the fourth dimension is over the faces of the element (south, east,
!       north, west, bottom, top), and the last dimension cycles over the elements in a spectral
!        element mesh.
!      
!     gradF(1:3,0:N,0:N,0:N,1:nVariables,1:nElements) (out)
!       The divergence of the vector functions over all of the elements
!
!     nVariables (in)
!       The number of variables
!
!     nElements (in)
!       The number of elements in the spectral element mesh.
!
! ================================================================================================ ! 

  SUBROUTINE DG_Gradient_3D( myNodal, f, fAtBoundaries, gradF, nVariables, nElements )
    IMPLICIT NONE
    CLASS( NodalDG ) :: myNodal
#ifdef HAVE_CUDA
    INTEGER, DEVICE, INTENT(in)     :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: fAtBoundaries(0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:6, 1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: gradF(1:3,0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    ! Local
    TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3( 4*(ceiling( REAL(myNodal % N+1)/4 ) ), &
                     4*(ceiling( REAL(myNodal % N+1)/4 ) ) , &
                     4*(ceiling( REAL(myNodal % N+1)/4 ) ) )
      grid = dim3( nVariables, nElements, 1)  

      CALL DG_Gradient_3D_CUDAKernel<<<grid, tBlock>>>( f, fAtBoundaries, gradF, &
                                                        myNodal % boundaryInterpolationMatrix_dev, &
                                                        myNodal % dgDerivativeMatrixTranspose_dev, &
                                                        myNodal % quadratureWeights_dev, &
                                                        myNodal % N_dev, nVariables, nElements )


#else
    INTEGER         :: nVariables, nElements
    REAL(prec)      :: f(0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    REAL(prec)      :: fAtBoundaries(0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:6, 1:nElements)
    REAL(prec)      :: gradf(1:3,0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)

      gradF = DG_Gradient_3D_NodalDG( myNodal, f, fAtBoundaries, nVariables, nElements )
    
#endif

  END SUBROUTINE DG_Gradient_3D

! ================================================================================================ !
! ------------------------------------- PRIVATE ROUTINES ----------------------------------------- !
! ================================================================================================ !


  FUNCTION CalculateFunctionsAtBoundaries_3D_NodalDG( myNodal, f, nVariables, nElements ) RESULT( fAtBoundaries )
    IMPLICIT NONE
    TYPE( NodalDG ) :: myNodal
    INTEGER         :: nVariables, nElements
    REAL(prec)      :: f(0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    REAL(prec)      :: fAtBoundaries(0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:6, 1:nElements)
    ! Local
    INTEGER :: i, j, k, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO j = 0, myNodal % N
            DO i = 0, myNodal % N
            
              fAtBoundaries(i,j,iVar,1:6,iEl) = 0.0_prec
              
              DO k = 0, myNodal % N
               
                fAtBoundaries(i,j,iVar,1,iEl) = fAtBoundaries(i,j,iVar,1,iEl) + myNodal % boundaryInterpolationMatrix(k,0)*f(k,i,j,iVar,iEl) ! South
                fAtBoundaries(i,j,iVar,2,iEl) = fAtBoundaries(i,j,iVar,2,iEl) + myNodal % boundaryInterpolationMatrix(k,1)*f(i,k,j,iVar,iEl) ! East
                fAtBoundaries(i,j,iVar,3,iEl) = fAtBoundaries(i,j,iVar,3,iEl) + myNodal % boundaryInterpolationMatrix(k,1)*f(k,i,j,iVar,iEl) ! North
                fAtBoundaries(i,j,iVar,4,iEl) = fAtBoundaries(i,j,iVar,4,iEl) + myNodal % boundaryInterpolationMatrix(k,0)*f(i,k,j,iVar,iEl) ! West
                fAtBoundaries(i,j,iVar,5,iEl) = fAtBoundaries(i,j,iVar,5,iEl) + myNodal % boundaryInterpolationMatrix(k,0)*f(i,j,k,iVar,iEl) ! Bottom
                fAtBoundaries(i,j,iVar,6,iEl) = fAtBoundaries(i,j,iVar,6,iEl) + myNodal % boundaryInterpolationMatrix(k,1)*f(i,j,k,iVar,iEl) ! Top
               
              ENDDO
              
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END FUNCTION CalculateFunctionsAtBoundaries_3D_NodalDG
 
  FUNCTION DG_Divergence_3D_NodalDG( myNodal, f, fnAtBoundaries, nVariables, nElements ) RESULT( divF )
    IMPLICIT NONE
    TYPE( NodalDG ) :: myNodal
    INTEGER         :: nVariables, nElements
    REAL(prec)      :: f(1:3,0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    REAL(prec)      :: fnAtBoundaries(0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:6, 1:nElements)
    REAL(prec)      :: divf(0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    ! Local
    INTEGER    :: ii, i, j, k, iVar, iEl
    REAL(prec) :: df

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, myNodal % N
            DO j = 0, myNodal % N
              DO i = 0, myNodal % N   
 
                df = 0.0_prec
                DO ii = 0, myNodal % N
                  df = df + myNodal % dgDerivativeMatrixTranspose(ii,i)*f(1,ii,j,k,iVar,iEl) + &
                            myNodal % dgDerivativeMatrixTranspose(ii,j)*f(2,i,ii,k,iVar,iEl) + &
                            myNodal % dgDerivativeMatrixTranspose(ii,k)*f(3,i,j,ii,iVar,iEl)
                ENDDO
                 
                divF(i,j,k,iVar,iEl) =  df+ ( fnAtBoundaries(i,k,iVar,1,iEl)*myNodal % boundaryInterpolationMatrix(j,0) + &
                                                fnAtBoundaries(i,k,iVar,3,iEl)*myNodal % boundaryInterpolationMatrix(j,1) )/&
                                              myNodal % quadratureWeights(j) + &
                                              ( fnAtBoundaries(j,k,iVar,4,iEl)*myNodal % boundaryInterpolationMatrix(i,0) + &
                                                fnAtBoundaries(j,k,iVar,2,iEl)*myNodal % boundaryInterpolationMatrix(i,1) )/&
                                              myNodal % quadratureWeights(i) + &
                                              ( fnAtBoundaries(i,j,iVar,5,iEl)*myNodal % boundaryInterpolationMatrix(k,0) + &
                                                fnAtBoundaries(i,j,iVar,6,iEl)*myNodal % boundaryInterpolationMatrix(k,1) )/&
                                              myNodal % quadratureWeights(k) 

              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
                            

  END FUNCTION DG_Divergence_3D_NodalDG

  FUNCTION DG_Gradient_3D_NodalDG( myNodal, f, fAtBoundaries, nVariables, nElements ) RESULT( gradF )
    IMPLICIT NONE
    TYPE( NodalDG ) :: myNodal
    INTEGER         :: nVariables, nElements
    REAL(prec)      :: f(0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    REAL(prec)      :: fAtBoundaries(0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:6, 1:nElements)
    REAL(prec)      :: gradf(1:3,0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    ! Local
    INTEGER    :: ii, i, j, k, iVar, iEl

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, myNodal % N
            DO j = 0, myNodal % N
              DO i = 0, myNodal % N   
 
                gradF(1:3,i,j,k,iVar,iEl) = 0.0_prec

                DO ii = 0, myNodal % N
                  gradF(1,i,j,k,iVar,iEl) = gradF(1,i,j,k,iVar,iEl) + myNodal % dgDerivativeMatrixTranspose(ii,i)*f(ii,j,k,iVar,iEl)
                  gradF(2,i,j,k,iVar,iEl) = gradF(2,i,j,k,iVar,iEl) + myNodal % dgDerivativeMatrixTranspose(ii,j)*f(i,ii,k,iVar,iEl)
                  gradF(3,i,j,k,iVar,iEl) = gradF(3,i,j,k,iVar,iEl) + myNodal % dgDerivativeMatrixTranspose(ii,k)*f(i,j,ii,iVar,iEl)
                ENDDO
                 
                gradF(1,i,j,k,iVar,iEl) = -( gradF(1,i,j,k,iVar,iEl) + ( fAtBoundaries(j,k,iVar,4,iEl)*myNodal % boundaryInterpolationMatrix(i,0) + &
                                                                         fAtBoundaries(j,k,iVar,2,iEl)*myNodal % boundaryInterpolationMatrix(i,1) )/&
                                                                       myNodal % quadratureWeights(i)  )

                gradF(2,i,j,k,iVar,iEl) = -( gradF(2,i,j,k,iVar,iEl) + ( fAtBoundaries(i,k,iVar,1,iEl)*myNodal % boundaryInterpolationMatrix(j,0) + &
                                                                         fAtBoundaries(i,k,iVar,3,iEl)*myNodal % boundaryInterpolationMatrix(j,1) )/&
                                                                       myNodal % quadratureWeights(j) )

                gradF(3,i,j,k,iVar,iEl) = -( gradF(3,i,j,k,iVar,iEl) + ( fAtBoundaries(i,j,iVar,5,iEl)*myNodal % boundaryInterpolationMatrix(k,0) + &
                                                                         fAtBoundaries(i,j,iVar,6,iEl)*myNodal % boundaryInterpolationMatrix(k,1) )/&
                                                                       myNodal % quadratureWeights(k) )

              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
                            

  END FUNCTION DG_Gradient_3D_NodalDG

#ifdef HAVE_CUDA
  ATTRIBUTES(Global) SUBROUTINE CalculateFunctionsAtBoundaries_3D_CUDAKernel( f, fAtBoundaries, boundaryMatrix, N, nVariables, nElements ) 
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: N, nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,0:N,0:N,1:nVariables,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryMatrix(0:N,0:1)
    REAL(prec), DEVICE, INTENT(out) :: fAtBoundaries(0:N,0:N,1:nVariables,1:6,1:nElements)
    ! Local
    INTEGER    :: iVar, iEl, i, j, k
    REAL(prec) :: bSol(1:6)

      iEl = blockIdx % x
      
      iVar = threadIdx % z
      k   = threadIdx % y-1
      j   = threadIdx % x-1
      
      IF( j <= N .AND. k <= N )THEN
      bSol(1:6) = 0.0_prec

      DO i = 0, N

        bSol(1) = bSol(1) + boundaryMatrix(i,0)*f(j,i,k,iVar,iEl) ! south
        bSol(2) = bSol(2) + boundaryMatrix(i,1)*f(i,j,k,iVar,iEl) ! east
        bSol(3) = bSol(3) + boundaryMatrix(i,1)*f(j,i,k,iVar,iEl) ! north
        bSol(4) = bSol(4) + boundaryMatrix(i,0)*f(i,j,k,iVar,iEl) ! west
        bSol(5) = bSol(5) + boundaryMatrix(i,0)*f(j,k,i,iVar,iEl) ! botom
        bSol(6) = bSol(6) + boundaryMatrix(i,1)*f(j,k,i,iVar,iEl) ! top

      ENDDO
               
      DO i = 1, 6
        fAtBoundaries(j,k,iVar,i,iEl) = bSol(i)
      ENDDO
      
      ENDIF
      
  END SUBROUTINE CalculateFunctionsAtBoundaries_3D_CUDAKernel

  ATTRIBUTES(Global) SUBROUTINE DG_Divergence_3D_CUDAKernel( f, fnAtBoundaries, divF, boundaryMatrix, dgDerivativeMatrixTranspose, quadratureWeights, N, nVariables, nElements )
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: N, nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(1:3,0:N,0:N,0:N,1:nVariables,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: fnAtBoundaries(0:N,0:N,1:nVariables,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryMatrix(0:N,0:1)
    REAL(prec), DEVICE, INTENT(in)  :: dgDerivativeMatrixTranspose(0:N,0:N)
    REAL(prec), DEVICE, INTENT(in)  :: quadratureWeights(0:N)
    REAL(prec), DEVICE, INTENT(out) :: divF(0:N,0:N,0:N,1:nVariables,1:nElements)
    ! Local
    INTEGER            :: i, j, k, iVar, iEl, ii
    REAL(prec)         :: df
    REAL(prec), SHARED :: fLocal(1:3,0:7,0:7,0:7)
    
    
      iVar = blockIDx % x
      iEl  = blockIDx % y
      
      i = threadIdx % x - 1
      j = threadIdx % y - 1
      k = threadIdx % z - 1
    
    
      IF( i <= N .AND. j <= N .AND. k <= N )THEN
      
        fLocal(1:3,i,j,k) = f(1:3,i,j,k,iVar,iEl)
        
      ENDIF
      
      CALL syncthreads( )
      
      IF( i <= N .AND. j <= N .AND. k <= N )THEN
      
        df = 0.0_prec
        DO ii = 0, N
          df = df + dgDerivativeMatrixTranspose(ii,i)*fLocal(1,ii,j,k) + &
                    dgDerivativeMatrixTranspose(ii,j)*fLocal(2,i,ii,k) + &
                    dgDerivativeMatrixTranspose(ii,k)*fLocal(3,i,j,ii)
        ENDDO
       
        divF(i,j,k,iVar,iEl) = ( df+ ( fnAtBoundaries(i,k,iVar,1,iEl)*boundaryMatrix(j,0) + &
                                        fnAtBoundaries(i,k,iVar,3,iEl)*boundaryMatrix(j,1) )/&
                                      quadratureWeights(j) + &
                                      ( fnAtBoundaries(j,k,iVar,4,iEl)*boundaryMatrix(i,0) + &
                                        fnAtBoundaries(j,k,iVar,2,iEl)*boundaryMatrix(i,1) )/&
                                      quadratureWeights(i) + &
                                      ( fnAtBoundaries(i,j,iVar,5,iEl)*boundaryMatrix(k,0) + &
                                        fnAtBoundaries(i,j,iVar,6,iEl)*boundaryMatrix(k,1) )/&
                                      quadratureWeights(k) )
                                      
      ENDIF          
                  
  END SUBROUTINE DG_Divergence_3D_CUDAKernel

  ATTRIBUTES(Global) SUBROUTINE DG_Gradient_3D_CUDAKernel( f, fAtBoundaries, gradF, boundaryMatrix, dgDerivativeMatrixTranspose, quadratureWeights, N, nVariables, nElements )
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)    :: N, nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,0:N,0:N,1:nVariables,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: fAtBoundaries(0:N,0:N,1:nVariables,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryMatrix(0:N,0:1)
    REAL(prec), DEVICE, INTENT(in)  :: dgDerivativeMatrixTranspose(0:N,0:N)
    REAL(prec), DEVICE, INTENT(in)  :: quadratureWeights(0:N)
    REAL(prec), DEVICE, INTENT(out) :: gradF(1:3,0:N,0:N,1:nVariables,1:6,1:nElements)
    ! Local
    INTEGER            :: i, j, k, iVar, iEl, ii
    REAL(prec)         :: df(1:3)
    REAL(prec), SHARED :: fLocal(0:7,0:7,0:7)
    
    
      iVar = blockIDx % x
      iEl  = blockIDx % y
      
      i = threadIdx % x - 1
      j = threadIdx % y - 1
      k = threadIdx % z - 1
    
      IF( i <= N .AND. j <= N .AND. k <= N )THEN
      fLocal(i,j,k) = f(i,j,k,iVar,iEl)
    
      CALL syncthreads( )
      
      df(1) = 0.0_prec
      df(2) = 0.0_prec
      df(3) = 0.0_prec

      DO ii = 0, N
        df(1) = df(1) + dgDerivativeMatrixTranspose(ii,i)*fLocal(ii,j,k)
        df(2) = df(2) + dgDerivativeMatrixTranspose(ii,j)*fLocal(i,ii,k)
        df(3) = df(3) + dgDerivativeMatrixTranspose(ii,k)*fLocal(i,j,ii)
      ENDDO
       
      gradF(1,i,j,k,iVar,iEl) = -( df(1) + ( fAtBoundaries(j,k,iVar,4,iEl)*boundaryMatrix(i,0) + &
                                             fAtBoundaries(j,k,iVar,2,iEl)*boundaryMatrix(i,1) )/&
                                           quadratureWeights(i) )

      gradF(2,i,j,k,iVar,iEl) = -( df(2) + ( fAtBoundaries(i,k,iVar,1,iEl)*boundaryMatrix(j,0) + &
                                             fAtBoundaries(i,k,iVar,3,iEl)*boundaryMatrix(j,1) )/&
                                           quadratureWeights(j) )

      gradF(3,i,j,k,iVar,iEl) = -( df(3) + ( fAtBoundaries(i,j,iVar,5,iEl)*boundaryMatrix(k,0) + &
                                             fAtBoundaries(i,j,iVar,6,iEl)*boundaryMatrix(k,1) )/&
                                           quadratureWeights(k) )
        
      ENDIF
                  
  END SUBROUTINE DG_Gradient_3D_CUDAKernel

#endif

END MODULE NodalDG_Class
