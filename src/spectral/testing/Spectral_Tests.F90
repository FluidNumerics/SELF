! Spectral_Tests.F90
!
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
!
!
!  This program exercises the Lagrange and NodalDG routines for performing interpolation and differentiation 
!  operations. A reference interpolant of degree 50 is used to estimate errors in lower order interpolation.
!
!  For 3-D cases, interpolation and differentation tests are performed on the following analytical
!  functions : 
!
!    (1)  f(x,y,z) = x*y*z (scalar function, exactness)
!
!         In this case, the interpolation should be exact for all polynomials of degree 1 or higher.
!         Interpolation errors and differentiation errors should be to machine precision accuracy.
!
!    (2)  f(x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z) (scalar function, exponential error decay)
!   
!         In this case, the interpolation errors should decay exponentially with polynomial degree N.
!
!    (3)  \vec{f}(x,y,z) = y*z \hat{x} + x*z \hat{y} + x*y \hat{z} ( exactness, divergence free )
!
!    (4)  \vec{f}(x,y,z) = sin(pi*y*z) \hat{x} + cos(pi*x*z) \hat{y} + sin(pi*x)*cos(pi*y) \hat{z} ( exponential error decay, divergence free )

PROGRAM Spectral_Tests


USE ModelPrecision
USE CommonRoutines
USE Lagrange_Class
USE NodalDG_Class
USE NodalDGSolution_3D_Class


IMPLICIT NONE

  INTEGER, PARAMETER :: polyLow = 3
  INTEGER, PARAMETER :: polyHigh = 20

  TYPE( NodalDG ) :: referenceInterpolant
  TYPE( NodalDG ) :: trialInterpolant
  
  TYPE( NodalDGSolution_3D ) :: referenceFunctions
  TYPE( NodalDGSolution_3D ) :: interpolatedFunctions
  TYPE( NodalDGSolution_3D ) :: trialFunctions
  
  REAL(prec) :: F1_interpolation_error(polyLow:polyHigh)
  REAL(prec) :: F1_divergence_strong_error(polyLow:polyHigh)
!  REAL(prec) :: F1_divergence_dgweak_error(polyLow:polyHigh)
  

  
    CALL referenceInterpolant % Build( targetPoints  = UniformPoints( -1.0_prec, 1.0_prec, 100 ), &
                                       N             = 50, &
                                       nTargetPoints = 100, &
                                       quadrature    = GAUSS  )
                                       
    CALL referenceFunctions % Build( N          = 50, &
                                     nEquations = 1, &
                                     nElements  = 1 )
                 
    CALL interpolatedFunctions % Build( N          = 50, &
                                        nEquations = 1, &
                                        nElements  = 1 )
                                     

    CALL F1_Case( )


    CALL interpolatedFunctions % Trash( )  
    CALL referenceFunctions % Trash( )
    CALL referenceInterpolant % Trash( )
  
  
CONTAINS

  ! Case 1 - Scalar Function, exactness
  SUBROUTINE F1_Case( )
    IMPLICIT NONE
    INTEGER :: i
  
      PRINT*, '  Case f_1(x,y,z) = x*y*z  '
      PRINT*, ' ----------------------------------------------------------------------------- '
      
      CALL F1( referenceInterpolant % interp, &
               referenceFunctions % solution, & 
               referenceFunctions % flux, 1, 1 )

#ifdef HAVE_CUDA
      CALL referenceFunctions % UpdateDevice( )
#endif
  
      ! The main loop loops over polynomial degree. For each polynomial degree, max error in interpolation and 
      ! differentiation operations are estimated and stored in the error arrays
      DO i = polyLow, polyHigh

        ! Allocate space for the NodalDG and DGSolution_3D structures (trial function)
        CALL trialInterpolant % Build( targetPoints  = referenceInterpolant % interp % interpolationPoints, &
                                       N             = i, &
                                       nTargetPoints = referenceInterpolant % N, &
                                       quadrature    = GAUSS  )
                                       
        CALL trialFunctions % Build( N          = i, &
                                     nEquations = 1, &
                                     nElements  = 1 )  
                    
        ! The example function for case 1 is filled in the "solution" attribute of the trialFunctions data structure.
        CALL F1( trialInterpolant % interp, &
                 trialFunctions % solution, & 
                 trialFunctions % flux, 1, 1 )

#ifdef HAVE_CUDA
        CALL trialFunctions % UpdateDevice( )
#endif
  
        ! Now that we have a lower order interpolant of the example function, it can be interpolated onto the 
        ! reference  ( O(50) )quadrature points so that a max error can be estimated
#ifdef HAVE_CUDA
        CALL trialInterpolant % interp % ApplyInterpolationMatrix_3D( f          = trialFunctions % solution_dev, &
                                                                      fNew       = interpolatedFunctions % solution_dev, & 
                                                                      nVariables = 1, &
                                                                      nElements  = 1) 
#else
        CALL trialInterpolant % interp % ApplyInterpolationMatrix_3D( f          = trialFunctions % solution, &
                                                                      fNew       = interpolatedFunctions % solution, & 
                                                                      nVariables = 1, &
                                                                      nElements  = 1) 
#endif

        F1_interpolation_error(i) = MAXVAL( ABS( interpolatedFunctions % solution - referenceFunctions % solution ) )
  
  
        ! Divergence test (strong form)
        CALL trialInterpolant % interp % CalculateDivergence_3D( f          = trialFunctions % flux, &
                                                                 divF       = trialFunctions % tendency, &
                                                                 nVariables = 1, &
                                                                 nElements  = 1 )
                                                                 
        F1_divergence_strong_error(i) = MAXVAL( ABS(trialFunctions % tendency) )
        
        PRINT*, '  ', i, F1_interpolation_error(i), F1_divergence_strong_error(i)
  
        CALL trialFunctions % Trash( )
        CALL trialInterpolant % Trash( )
  
      ENDDO
    
  END SUBROUTINE F1_Case
  
  SUBROUTINE F1( interp, f1f, f1Gradient, nVariables, nElements )
    IMPLICIT NONE
    TYPE( Lagrange ), INTENT(in) :: interp 
    INTEGER, INTENT(in)          :: nVariables, nElements
    REAL(prec), INTENT(out)      :: f1f(0:interp % N, 0:interp % N, 0:interp % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out)      :: f1Gradient(1:3,0:interp % N, 0:interp % N, 0:interp % N, 1:nVariables, 1:nElements)
    ! Local
    INTEGER :: i, j, k, iVar, iEl
    
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, interp % N
            DO j = 0, interp % N
              DO i = 0, interp % N
          
                f1f(i,j,k,iVar,iEl) = interp % interpolationPoints(i)*&
                                      interp % interpolationPoints(j)*&
                                      interp % interpolationPoints(k)
                        
                f1Gradient(1,i,j,k,iVar,iEl) = interp % interpolationPoints(j)*&
                                               interp % interpolationPoints(k)
                                  
                f1Gradient(2,i,j,k,iVar,iEl) = interp % interpolationPoints(i)*&
                                               interp % interpolationPoints(k)
                                  
                f1Gradient(3,i,j,k,iVar,iEl) = interp % interpolationPoints(i)*&
                                               interp % interpolationPoints(j)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
                                  
    
  END SUBROUTINE F1
  
END PROGRAM Spectral_Tests
