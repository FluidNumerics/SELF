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
!    (1) Exactness
!        (Scalar) f(x,y,z) = x*y*z 
!        (Vector) \vec{g}(x,y,z) = y*z \hat{x} + x*z \hat{y} + x*y \hat{z}
!                 
!                  ===> divergence( \vec{g} ) = 0.0
!
!         In this case, the interpolation should be exact for all polynomials of degree 1 or higher.
!         Interpolation errors and differentiation errors should be to machine precision accuracy.
!
!    (2) Exponential Error Decay
!        (Scalar)  f(x,y,z)       = sin(pi*x)*sin(pi*y)*sin(pi*z)
!        (Vector)  \vec{g}(x,y,z) = -pi*cos(pi*x)*sin(pi*y)*sin(pi*z) \hat{x} +
!                                   -pi*sin(pi*x)*cos(pi*y)*sin(pi*z) \hat{y} +
!                                   -pi*sin(pi*x)*sin(pi*y)*cos(pi*z) \hat{z}
!

PROGRAM Spectral_Tests


USE ModelPrecision
USE CommonRoutines
USE Lagrange_Class
USE NodalDG_Class
USE NodalDGSolution_3D_Class


IMPLICIT NONE

  INTEGER, PARAMETER :: polyLow = 3
#ifdef HAVE_CUDA
  INTEGER, PARAMETER :: nRef = 7
  INTEGER, PARAMETER :: polyHigh = 7
#else
  INTEGER, PARAMETER :: nRef = 50
  INTEGER, PARAMETER :: polyHigh = 20
#endif

#ifdef HAVE_CUDA
  INTEGER, ALLOCATABLE, DEVICE :: nVars, nElems
#endif

  TYPE( NodalDG ) :: referenceInterpolant
  TYPE( NodalDG ) :: trialInterpolant
  
  TYPE( NodalDGSolution_3D ) :: referenceFunctions
  TYPE( NodalDGSolution_3D ) :: interpolatedFunctions
  TYPE( NodalDGSolution_3D ) :: trialFunctions
  
  REAL(prec) :: F1_interpolation_error(polyLow:polyHigh)
  REAL(prec) :: F1_divergence_strong_error(polyLow:polyHigh)
  REAL(prec) :: F1_divergence_dgweak_error(polyLow:polyHigh)
  
#ifdef HAVE_CUDA
    ALLOCATE( nVars, nElems )
    nVars  = 1
    nElems = 1
#endif

    CALL referenceInterpolant % Build( targetPoints  = UniformPoints( -1.0_prec, 1.0_prec, 100 ), &
                                       N             = nRef, &
                                       nTargetPoints = 100, &
                                       quadrature    = GAUSS  )

    CALL referenceFunctions % Build( N          = nRef, &
                                     nEquations = 1, &
                                     nElements  = 1 )

    CALL interpolatedFunctions % Build( N          = nRef, &
                                        nEquations = 1, &
                                        nElements  = 1 )
                                     

    CALL F1_Case( )


    CALL interpolatedFunctions % Trash( )  
    CALL referenceFunctions % Trash( )
    CALL referenceInterpolant % Trash( )
#ifdef HAVE_CUDA
    DEALLOCATE( nVars, nElems )
#endif
  
CONTAINS

  ! Case 1 - Scalar Function, exactness
  SUBROUTINE F1_Case( )
    IMPLICIT NONE
    INTEGER :: i
  
      WRITE(*,*) '  Case f_1(x,y,z) = x*y*z  '
      WRITE(*,*) ' ------------------------------------------------------------------------------------------- '
      WRITE(*,*) ' Polynomial degree   Interpolation Error  Divergence Error (strong)  Divergence Error (weak)'
      WRITE(*,*) ' ------------------------------------------------------------------------------------------- '
      
      CALL F1( referenceInterpolant % interp, &
               referenceFunctions % solution, & 
               referenceFunctions % flux, &
               referenceFunctions % boundaryFlux, 1, 1 )

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
                 trialFunctions % flux, &
                 trialFunctions % boundaryFlux, &
                 1, 1 )
                 
        CALL RunTests 

        ! Report max errors !
        F1_interpolation_error(i)     = MAXVAL( ABS( interpolatedFunctions % solution - referenceFunctions % solution ) )                                                         
        F1_divergence_strong_error(i) = MAXVAL( ABS( trialFunctions % tendency ) )
        F1_divergence_dgweak_error(i) = MAXVAL( ABS( trialFunctions % source ) )

        WRITE(*,'(8x,I3,5x,2(4x,E17.5),8x,E17.5)') i, F1_interpolation_error(i), F1_divergence_strong_error(i), F1_divergence_dgweak_error(i)
        
  
        CALL trialFunctions % Trash( )
        CALL trialInterpolant % Trash( )
  
      ENDDO
    
  END SUBROUTINE F1_Case
  
  SUBROUTINE F1( interp, f1f, f1Gradient, boundaryFlux, nVariables, nElements )
    IMPLICIT NONE
    TYPE( Lagrange ), INTENT(in) :: interp 
    INTEGER, INTENT(in)          :: nVariables, nElements
    REAL(prec), INTENT(out)      :: f1f(0:interp % N, 0:interp % N, 0:interp % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out)      :: f1Gradient(1:3,0:interp % N, 0:interp % N, 0:interp % N, 1:nVariables, 1:nElements)
    REAL(prec), INTENT(out)      :: boundaryFlux(0:interp % N, 0:interp % N, 1:nVariables, 1:6, 1:nElements)
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
          
          DO k = 0, interp % N
            DO j = 0, interp % N
          
                ! south
                boundaryFlux(j,k,iVar,1,iEl) = -interp % interpolationPoints(j)*&
                                                interp % interpolationPoints(k)
                ! east
                boundaryFlux(j,k,iVar,2,iEl) = interp % interpolationPoints(j)*&
                                                interp % interpolationPoints(k)
                ! north
                boundaryFlux(j,k,iVar,3,iEl) = interp % interpolationPoints(j)*&
                                                interp % interpolationPoints(k)
                ! west
                boundaryFlux(j,k,iVar,4,iEl) = -interp % interpolationPoints(j)*&
                                                interp % interpolationPoints(k)
                ! bottom
                boundaryFlux(j,k,iVar,5,iEl) = -interp % interpolationPoints(j)*&
                                                interp % interpolationPoints(k)
                ! top
                boundaryFlux(j,k,iVar,6,iEl) = interp % interpolationPoints(j)*&
                                                interp % interpolationPoints(k)                                                                                                
            ENDDO
          ENDDO
          
        ENDDO
      ENDDO
                                  
    
  END SUBROUTINE F1
  
  SUBROUTINE RunTests( ) 
  
#ifdef HAVE_CUDA
        
        CALL trialFunctions % UpdateDevice( )
        
!        ! Interpolation test
        CALL trialInterpolant % interp % ApplyInterpolationMatrix_3D( f          = trialFunctions % solution_dev, &
                                                                      fNew       = interpolatedFunctions % solution_dev, & 
                                                                      nVariables = nVars, &
                                                                      nElements  = nElems)
                                                                      
        ! Divergence test (strong form)                                                              
        CALL trialInterpolant % interp % CalculateDivergence_3D( f          = trialFunctions % flux_dev, &
                                                                 divF       = trialFunctions % tendency_dev, &
                                                                 nVariables = nVars, &
                                                                 nElements  = nElems )
        ! Divergence test (weak form)                                                         
        CALL trialInterpolant % DG_Divergence_3D( f              = trialFunctions % flux_dev, &
                                                  fnAtBoundaries = trialFunctions % boundaryFlux_dev, & 
                                                  divF           = trialFunctions % source_dev, &
                                                  nVariables     = nVars, &
                                                  nElements      = nElems )
                                                  
        ! Gradient test (strong form)
        CALL trialInterpolant % interp % CalculateGradient_3D( f          = trialFunctions % solution_dev, &
                                                               gradF      = trialFunctions % flux_dev, &
                                                               nVariables = nVars, &
                                                               nElements  = nElems )
                                                                      
        CALL interpolatedFunctions % UpdateHost( )
        CALL trialFunctions % UpdateHost( )
        
#else
        ! Interpolation test
        CALL trialInterpolant % interp % ApplyInterpolationMatrix_3D( f          = trialFunctions % solution, &
                                                                      fNew       = interpolatedFunctions % solution, & 
                                                                      nVariables = 1, &
                                                                      nElements  = 1) 
        ! Divergence test (strong form) 
        CALL trialInterpolant % interp % CalculateDivergence_3D( f          = trialFunctions % flux, &
                                                                 divF       = trialFunctions % tendency, &
                                                                 nVariables = 1, &
                                                                 nElements  = 1 )
                                                                 
        ! Divergence test (weak form)                                                         
        CALL trialInterpolant % DG_Divergence_3D( f              = trialFunctions % flux, &
                                                  fnAtBoundaries = trialFunctions % boundaryFlux, & 
                                                  divF           = trialFunctions % source, &
                                                  nVariables     = 1, &
                                                  nElements      = 1 )
                                                  
        ! Gradient test (strong form)
        CALL trialInterpolant % interp % CalculateGradient_3D( f          = trialFunctions % solution, &
                                                               gradF      = trialFunctions % flux, &
                                                               nVariables = 1, &
                                                               nElements  = 1 )
#endif
  END SUBROUTINE RunTests
  
END PROGRAM Spectral_Tests
