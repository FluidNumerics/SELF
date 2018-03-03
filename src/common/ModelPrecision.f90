! ModelPrecision.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! ModelPrecision.f90 is part of the Spectral Element Libraries in Fortran for Fluids (SELF-Fluids).
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

 
MODULE ModelPrecision

INTEGER, PARAMETER :: sp   = SELECTED_REAL_KIND(6, 37)     ! 32-bit
INTEGER, PARAMETER :: dp   = SELECTED_REAL_KIND(15, 307)   ! 64-bit
INTEGER, PARAMETER :: qp   = SELECTED_REAL_KIND(33, 4931)  ! 128-bit

#ifdef DOUBLE_PRECISION
INTEGER, PARAMETER :: prec = dp
#else
INTEGER, PARAMETER :: prec = sp
#endif

END MODULE ModelPrecision
