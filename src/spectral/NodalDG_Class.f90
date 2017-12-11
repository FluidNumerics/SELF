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
USE Lagrange_Class

IMPLICIT NONE

!> \addtogroup NodalDG_Class 
!! @{

!> \struct NodalDG
!!  The NodalDG class contains attributes needed for implementing spectral element methods
!!  in 3-D.
!!  
!!  An interpolant is formed that handles mapping between a computational "quadrature" mesh and a
!!  uniform plotting mesh. Quadrature (integration) weights are stored for use with Galerkin type
!!  methods. Galerkin derivative matrices (collocation derivative matrices weighted by appropriate
!!  ratios of the quadrature weights) are stored to facilitate the computation of weak derivatives.
!!  Finally, an interpolation matrix and accompanying subroutine is provided to interpolate 3-D data, 
!!  defined on the quadrature mesh, to the element boundaries.
!!
!! <H2> NodalDG </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> N <td> INTEGER  <td> Polynomial degree of the spectral element method
!!       <tr> <th> nPlot <td> INTEGER <td> Number uniform plotting points
!!       <tr> <th> interp <td> Lagrange <td> Lagrange interpolant
!!       <tr> <th> qWeight(0:N) <td> REAL(prec) <td> Quadrature integration weights
!!       <tr> <th> dMatS(0:N,0:N) <td> REAL(prec) <td> Either the DG or CG derivative matrix
!!       <tr> <th> dMatP(0:N,0:N) <td> REAL(prec) <td> Either the DG or CG derivative matrix transpose
!!       <tr> <th> bMat(0:1,0:N) <td> REAL(prec) <td> Matrix for interpolating data to element boundaries
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref NodalDG_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Build <td> Build_NodalDG
!!       <tr> <th> Trash <td> Trash_NodalDG
!!       <tr> <th> CalculateAtBoundaries_1D <td> CalculateAtBoundaries_1D_NodalDG
!!       <tr> <th> CalculateAtBoundaries_2D <td> CalculateAtBoundaries_2D_NodalDG
!!       <tr> <th> CalculateAtBoundaries_3D <td> CalculateAtBoundaries_3D_NodalDG
!!    </table>
!!

!>@}

  TYPE NodalDG
#ifdef HAVE_CUDA
    INTEGER, MANAGED        :: N    
    INTEGER, MANAGED        :: nTargetPoints
#else
    INTEGER                 :: N     
    INTEGER                 :: nTargetPoints
#endif
    TYPE(Lagrange)          :: interp 
    REAL(prec), ALLOCATABLE :: quadratureWeights(:)
    REAL(prec), ALLOCATABLE :: dgDerivativeMatrix(:,:)
    REAL(prec), ALLOCATABLE :: dgDerivativeMatrixTranspose(:,:)
    REAL(prec), ALLOCATABLE :: boundaryInterpolationMatrix(:,:)
      
#ifdef HAVE_CUDA
    REAL(prec), DEVICE, ALLOCATABLE :: quadratureWeights_dev(:)
    REAL(prec), DEVICE, ALLOCATABLE :: dgDerivativeMatrix_dev(:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: dgDerivativeMatrixTranspose_dev(:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: boundaryInterpolationMatrix_dev(:,:)
#endif
    CONTAINS

      ! Manual Constructors/Destructors
      PROCEDURE :: Build => Build_NodalDG
      PROCEDURE :: Trash => Trash_NodalDG

      ! Type-Specific
      PROCEDURE :: CalculateFunctionsAtBoundaries_3D
      PROCEDURE :: DG_Divergence_3D
      PROCEDURE :: DG_Gradient_3D
      
    END TYPE NodalDG
    
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup NodalDG_Class 
!! @{ 
! ================================================================================================ !
! S/R Build 
! 
!> \fn Build_NodalDG 
!!  Allocates space fills values for the NodalDG attributes using to the specified 
!!  quadrature and approximation form.
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodalDG) :: this <BR>
!! <B>INTEGER</B>               :: N, nPlot
!!         .... <BR>
!!     ! To build a  structure for Continuous Galerkin with Gauss-Lobatto quadrature <BR>
!!     <B>CALL</B> this % Build( N, nPlot, GAUSS_LOBATTO, CG ) <BR>
!!
!!     ! To build a  structure for Discontinuous Galerkin with Gauss quadrature <BR>
!!     <B>CALL</B> this % Build( N, nPlot, GAUSS, DG ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myNodal <td> NodalDG <td> On output, the attributes of the
!!                                                        NodalDG data structure are filled
!!                                                        in.
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the method.
!!   <tr> <td> in <th> nPlot <td> INTEGER <td> The number of uniform plotting points in each 
!!                                             computational direction.
!!   <tr> <td> in <th> quadrature <td> INTEGER <td> A flag for specifying the desired type of 
!!                                                  quadrature. Can be set to either GAUSS or
!!                                                  GAUSS_LOBATTO. See \ref ModelFlags.f90 for
!!                                                  flag definitions.
!!   <tr> <td> in <th> approxForm <td> INTEGER <td> A flag for specifying the type of method that 
!!                                                  you are using. Can be set to either CG or DG.
!!                                                  See \ref ModelFlags.f90 for flag definitions.
!! 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}

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
       
      CALL myNodal % interp % Build( N, nTargetPoints, quadratureWeights, targetPoints )
   
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
  
!> \addtogroup NodalDG_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash
! 
!> \fn Trash_NodalDG  
!! Frees memory held by the attributes of the NodalDG class. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodalDG) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myNodal <td> NodalDG <td>
!!                         On <B>input</B>, a NodalDG class that has previously been 
!!                         constructed. <BR>
!!                         On <B>output</B>, the memory held by the attributes of this 
!!                         data-structure have been freed.
!!                                                           
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}

  SUBROUTINE Trash_NodalDG( myNodal)
    IMPLICIT NONE
    CLASS( NodalDG ), INTENT(inout) :: myNodal

      CALL myNodal % interp % Trash( )
      DEALLOCATE( myNodal % dgDerivativeMatrix, &
                  myNodal % dgDerivativeMatrixTranspose, &
                  myNodal % boundaryInterpolationMatrix, &
                  myNodal % quadratureWeights )
                  
#ifdef HAVE_CUDA
      DEALLOCATE( myNodal % dgDerivativeMatrix_dev, &
                  myNodal % dgDerivativeMatrixTranspose_dev, &
                  myNodal % boundaryInterpolationMatrix_dev, &
                  myNodal % quadratureWeights_dev )
#endif


  END SUBROUTINE Trash_NodalDG


  SUBROUTINE CalculateFunctionsAtBoundaries_3D( myNodal, f, fAtBoundaries, nVariables, nElements )
    IMPLICIT NONE
    CLASS( NodalDG ), INTENT(in) :: myNodal 
#ifdef HAVE_CUDA
    INTEGER, MANAGED, INTENT(in)    :: nVariables, nElements
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
                                                                           myNodal % N, nVariables, nElements )
      
#else
    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: f(0:myNodal % N,0:myNodal % N,0:myNodal % N,1:nVariables,nElements)
    REAL(prec), INTENT(out) :: fAtBoundaries(0:myNodal % N,0:myNodal % N,1:nVariables,1:6,1:nElements)
     
      fAtBoundaries = CalculateFunctionsAtBoundaries_3D_NodalDG( myNodal, f, nVariables, nElements )

#endif    

    
    
  END SUBROUTINE CalculateFunctionsAtBoundaries_3D
  
  SUBROUTINE DG_Divergence_3D( myNodal, f, fnAtBoundaries, divF, nVariables, nElements )
    IMPLICIT NONE
    CLASS( NodalDG ) :: myNodal
#ifdef HAVE_CUDA
    INTEGER, MANAGED, INTENT(in)    :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(1:3,0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: fnAtBoundaries(0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:6, 1:nElements)
    REAL(prec), DEVICE, INTENT(out) :: divF(0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    ! Local
    TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3( 4*(ceiling( REAL(myNodal % N+1)/4 ) ), &
                     4*(ceiling( REAL(myNodal % N+1)/4 ) ) , &
                     4*(ceiling( REAL(myNodal % N+1)/4 ) ) )
      grid = dim3( nVariables, nElements, 1)  

      CALL DG_Divergence_3D_CUDAKernel<<<grid, tBlock>>>( f, fnAtBoundaries, divF, &
                                                          myNodal % boundaryInterpolationMatrix_dev, &
                                                          myNodal % dgDerivativeMatrixTranspose_dev, &
                                                          myNodal % quadratureWeights_dev, &
                                                          myNodal % N, nVariables, nElements )


#else
    INTEGER         :: nVariables, nElements
    REAL(prec)      :: f(1:3,0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)
    REAL(prec)      :: fnAtBoundaries(0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:6, 1:nElements)
    REAL(prec)      :: divf(0:myNodal % N, 0:myNodal % N, 0:myNodal % N, 1:nVariables, 1:nElements)


      divF = DG_Divergence_3D_NodalDG( myNodal, f, fnAtBoundaries, nVariables, nElements )
    
#endif

  END SUBROUTINE DG_Divergence_3D

  SUBROUTINE DG_Gradient_3D( myNodal, f, fAtBoundaries, gradF, nVariables, nElements )
    IMPLICIT NONE
    CLASS( NodalDG ) :: myNodal
#ifdef HAVE_CUDA
    INTEGER, MANAGED, INTENT(in)    :: nVariables, nElements
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
                                                        myNodal % N, nVariables, nElements )


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
                 
                divF(i,j,k,iVar,iEl) = -( df+ ( fnAtBoundaries(i,k,iVar,1,iEl)*myNodal % boundaryInterpolationMatrix(j,0) + &
                                                fnAtBoundaries(i,k,iVar,3,iEl)*myNodal % boundaryInterpolationMatrix(j,1) )/&
                                              myNodal % quadratureWeights(j) + &
                                              ( fnAtBoundaries(j,k,iVar,4,iEl)*myNodal % boundaryInterpolationMatrix(i,0) + &
                                                fnAtBoundaries(j,k,iVar,2,iEl)*myNodal % boundaryInterpolationMatrix(i,1) )/&
                                              myNodal % quadratureWeights(i) + &
                                              ( fnAtBoundaries(i,j,iVar,5,iEl)*myNodal % boundaryInterpolationMatrix(k,0) + &
                                                fnAtBoundaries(i,j,iVar,6,iEl)*myNodal % boundaryInterpolationMatrix(k,1) )/&
                                              myNodal % quadratureWeights(k) )

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
    INTEGER, MANAGED, INTENT(in)    :: N, nVariables, nElements
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
      
  END SUBROUTINE CalculateFunctionsAtBoundaries_3D_CUDAKernel

  ATTRIBUTES(Global) SUBROUTINE DG_Divergence_3D_CUDAKernel( f, fnAtBoundaries, divF, boundaryMatrix, dgDerivativeMatrixTranspose, quadratureWeights, N, nVariables, nElements )
    IMPLICIT NONE
    INTEGER, MANAGED, INTENT(in)    :: N, nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(1:3,0:N,0:N,0:N,1:nVariables,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: fnAtBoundaries(0:N,0:N,1:nVariables,1:6,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: boundaryMatrix(0:N,0:1)
    REAL(prec), DEVICE, INTENT(in)  :: dgDerivativeMatrixTranspose(0:N,0:N)
    REAL(prec), DEVICE, INTENT(in)  :: quadratureWeights(0:N)
    REAL(prec), DEVICE, INTENT(out) :: divF(0:N,0:N,1:nVariables,1:6,1:nElements)
    ! Local
    INTEGER            :: i, j, k, iVar, iEl, ii
    REAL(prec) :: df
    REAL(prec), SHARED :: fLocal(1:3,0:7,0:7,0:7)
    
    
      iVar = blockIDx % x
      iEl  = blockIDx % y
      
      i = threadIdx % x - 1
      j = threadIdx % y - 1
      k = threadIdx % z - 1
    
      fLocal(1:3,i,j,k) = f(1:3,i,j,k,iVar,iEl)
    
      CALL syncthreads( )
      
      df = 0.0_prec
      DO ii = 0, N
        df = df + dgDerivativeMatrixTranspose(ii,i)*fLocal(1,ii,j,k) + &
                  dgDerivativeMatrixTranspose(ii,j)*fLocal(2,i,ii,k) + &
                  dgDerivativeMatrixTranspose(ii,k)*fLocal(3,i,j,ii)
      ENDDO
       
      divF(i,j,k,iVar,iEl) = -( df+ ( fnAtBoundaries(i,k,iVar,1,iEl)*boundaryMatrix(j,0) + &
                                      fnAtBoundaries(i,k,iVar,3,iEl)*boundaryMatrix(j,1) )/&
                                    quadratureWeights(j) + &
                                    ( fnAtBoundaries(j,k,iVar,4,iEl)*boundaryMatrix(i,0) + &
                                      fnAtBoundaries(j,k,iVar,2,iEl)*boundaryMatrix(i,1) )/&
                                    quadratureWeights(i) + &
                                    ( fnAtBoundaries(i,j,iVar,5,iEl)*boundaryMatrix(k,0) + &
                                      fnAtBoundaries(i,j,iVar,6,iEl)*boundaryMatrix(k,1) )/&
                                    quadratureWeights(k) )
                  
                  
  END SUBROUTINE DG_Divergence_3D_CUDAKernel

  ATTRIBUTES(Global) SUBROUTINE DG_Gradient_3D_CUDAKernel( f, fAtBoundaries, gradF, boundaryMatrix, dgDerivativeMatrixTranspose, quadratureWeights, N, nVariables, nElements )
    IMPLICIT NONE
    INTEGER, MANAGED, INTENT(in)    :: N, nVariables, nElements
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
                  
                  
  END SUBROUTINE DG_Gradient_3D_CUDAKernel

#endif

END MODULE NodalDG_Class
