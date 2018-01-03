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

PROGRAM SpectralFilter_Tests

END PROGRAM SpectralFilter_Tests
