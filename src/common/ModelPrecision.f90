! ModelPrecision.f90
! 
! Copyright 2017 Joseph Schoonover < schoonover.numerics@gmail.com > 
!
!
! ModelPrecision.f90 is part of the Spectral Element Libraries in Fortran for Fluids (SELF-Fluids).
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

 
MODULE ModelPrecision

INTEGER, PARAMETER :: sp   = SELECTED_REAL_KIND(6, 37)     ! 32-bit
INTEGER, PARAMETER :: dp   = SELECTED_REAL_KIND(15, 307)   ! 64-bit
INTEGER, PARAMETER :: qp   = SELECTED_REAL_KIND(33, 4931)  ! 128-bit
INTEGER, PARAMETER :: prec = dp                            ! Specify the precision here

END MODULE ModelPrecision
