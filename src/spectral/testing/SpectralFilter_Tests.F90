! SpectralFilter_Tests.F90
!
! Copyright 2018 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
!
!
! This program exercises the 3-D filtering routines using a cutoff filter and a
! rolloff filter. The filters also come equipped with a routine for diagnosing
! modal coefficients. Here, it is used to verify internal consistency of the
! filtering process. The filtering is demonstrated on a test function that is
! the sum of a “smooth” function with low wave number variability and a function
! with high wave number variability.
!
! —> How to verify these procedures ?? Use an identity. Set test function equal
! to the 7th Lagrange interpolating polynomial. Filter that function with a
! modal cutoff filter with the cutoff set to 6. The result should be zero to
! machine precision.
!
! ========================================================================================================= !

USE Quadrature
USE NodalDG_Class
USE NodalDGSolution_3D_Class
USE SpectralFilter_Class


PROGRAM SpectralFilter_Tests

  INTEGER, PARAMETER :: polyDegree = 7
  INTEGER, PARAMETER :: nCutoff = 6
 
  TYPE( SpectralFilter )     :: modalFilter
  TYPE( NodalDG )            :: dgStorage
  TYPE( NodalDGSolution_3D ) :: testFunction
  TYPE( NodalDGSolution_3D ) :: filteredFunction


    CALL dgStorage % Build( targetPoints  = UniformPoints( -1.0_prec, 1.0_prec, 10 ), &
                            N             = polyDegree, &
                            nTargetPoints = 10, &
                            quadrature    = GAUSS  )

    CALL testFunction % Build( N          = polyDegree, &
                               nEquations = 1, &
                               nElements  = 1 )

    CALL filteredFunction % Build( N          = polyDegree, &
                                   nEquations = 1, &
                                   nElements  = 1 )

    CALL modalFilter % Build( dgStorage % interp % interpolationPoints, &
                              dgStorage % interp % interpolationWeights, &
                              polyDegree, nCutoff, ModalCutoff )

 

END PROGRAM SpectralFilter_Tests
