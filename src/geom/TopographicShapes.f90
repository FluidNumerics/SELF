
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.

MODULE TopographicShapes

! src/COMMON/
  USE ModelPrecision

  IMPLICIT NONE

CONTAINS

  FUNCTION GaussianHill( x, y )
    IMPLICIT NONE
    REAL(prec)            :: GaussianHill
    REAL(prec),INTENT(in) :: x, y

    GaussianHill = 0.5_prec*exp( -( (x-0.5_prec)**2/(2.0_prec*0.1_prec**2) + (y-0.5_prec)**2/(2.0_prec*0.1_prec**2)  ) )

  END FUNCTION GaussianHill

END MODULE TopographicShapes
