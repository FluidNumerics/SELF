! DGITIntegrator_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE DGITIntegrator_Class

! src/common/
 USE ModelPrecision
 USE ConstantsDictionary
 USE CommonRoutines
! src/interp/
 USE Quadrature
 USE Lagrange_Class

IMPLICIT NONE

   TYPE DGITIntegrator
      INTEGER                 :: N
      REAL(prec)              :: stepSize 
      TYPE(Lagrange)          :: interp            ! Lagrange interpolant
      REAL(prec), ALLOCATABLE :: rVector(:)        ! Weights for initial condition
      REAL(prec), ALLOCATABLE :: gMat(:,:)         ! Derivative matrix
      REAL(prec), ALLOCATABLE :: prolongMat(:)     ! Matrix for interpolating functions to the "right" side of an element 
      REAL(prec), ALLOCATABLE :: integratorMat(:,:)! Matrix for forward stepping DGIT solver
      CONTAINS

      ! Manual Constructors/Destructors
      PROCEDURE :: Build => Build_DGITIntegrator
      PROCEDURE :: Trash => Trash_DGITIntegrator

      
    END TYPE DGITIntegrator
    
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_DGITIntegrator( myNodal, N, stepSize, quadrature  )

   IMPLICIT NONE
   CLASS(DGITIntegrator), INTENT(out) :: myNodal
   INTEGER, INTENT(in)                :: N
   REAL(prec), INTENT(in)             :: stepSize
   INTEGER, INTENT(in)                :: quadrature
   !LOCAL
   INTEGER                 :: i, j, INFO
   REAL(prec), ALLOCATABLE :: tempS(:), tempQ(:), tempL(:), tempR(:), tempMat(:,:), WORK(:)
   INTEGER, ALLOCATABLE    :: IPIV(:)

      myNodal % N        = N
      myNodal % stepSize = stepSize
      
      ! Allocate space
      ALLOCATE( myNodal % rVector(0:N), &
                myNodal % gMat(0:N,0:N), &
                myNodal % prolongMat(0:N), &
                myNodal % integratorMat(0:N,0:N) )
                
      myNodal % rVector       = 0.0_prec
      myNodal % gMat          = 0.0_prec
      myNodal % prolongMat    = 0.0_prec
      myNodal % integratorMat = 0.0_prec

      ALLOCATE( tempS(0:N), tempQ(0:N), tempL(0:N), tempR(0:N), tempMat(0:N,0:N) )
      tempS   = 0.0_prec
      tempQ   = 0.0_prec
      tempL   = 0.0_prec
      tempR   = 0.0_prec
      tempMat = 0.0_prec

      ALLOCATE( IPIV(1:N+1), WORK(1:N+1) )
      ! Generate the quadrature
      CALL LegendreQuadrature( N, tempS, tempQ, quadrature )
       

      CALL myNodal % interp % Build( N, N, tempS, tempS )
   
      ! Calculate and store the interpolants evaluated at the endpoints
      tempR = myNodal % interp % CalculateLagrangePolynomials( 1.0_prec )
      tempL = myNodal % interp % CalculateLagrangePolynomials( -1.0_prec )
      
       ! For Discontinuous Galerkin, the matrix is transposed and multiplied by a ratio of quadrature
       ! weights.
      DO j = 0, N ! loop over the matrix rows
         DO i = 0, N ! loop over the matrix columns

            myNodal % gMat(i,j) = tempR(i)*tempR(j)/tempQ(i) - &
                                  myNodal % interp % D(j,i)*tempQ(j)/tempQ(i)
 
            tempMat(i,j) =stepsize*0.5_prec*( tempR(i)*tempR(j)/tempQ(i) - &
                                              myNodal % interp % D(j,i)*tempQ(j)/tempQ(i) )

         ENDDO
         myNodal % rVector(j)    = tempL(j)/tempQ(j)
         myNodal % prolongMat(j) = tempR(j)
         tempMat(j,j)            = tempMat(j,j) + 1.0_prec
      ENDDO
      
      

      ! Invert the integrator matrix
     ! CALL DGETRF( N+1, N+1, tempMat, N+1, IPIV, INFO )
     ! CALL DGETRI( N+1, tempMat, N+1, IPIV, WORK, N+1, INFO )

      myNodal % integratorMat = InvertSpectralOpMatrix( tempMat, N )
      myNodal % integratorMat = tempMat


      DEALLOCATE( tempS, tempQ, tempL, tempR, tempMat )
      DEALLOCATE( IPIV, WORK )

 END SUBROUTINE Build_DGITIntegrator
!
 SUBROUTINE Trash_DGITIntegrator( myNodal)

   IMPLICIT NONE
   CLASS(DGITIntegrator), INTENT(inout) :: myNodal

   CALL myNodal % interp % Trash( )
   DEALLOCATE( myNodal % rVector, myNodal % gMat, myNodal % prolongMat, myNodal % integratorMat )


 END SUBROUTINE Trash_DGITIntegrator
!
!
!==================================================================================================!
!--------------------------- Type Specific Routines -----------------------------------------------!
!==================================================================================================!
!
!
END MODULE DGITIntegrator_Class
