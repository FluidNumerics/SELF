! ModelPrecision.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>, The Florida State University
! Copyright 2016 Joseph Schoonover <jschoonover@lanl.gov>, Los Alamos National Laboratory
!
! The SELF and accompanying documentation were produced in part under the 
! support of Florida State University and the National Science Foundation 
! through Grant OCE-1049131 during 2015 and in part  the support of the 
! Center for Nonlinear Studies and the Department of Energy through the 
! LANL/LDRD program in 2016.
!
! ModelPrecision.f90 is part of the Spectral Element Libraries in Fortran (SELF).
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file ModelPrecision.f90
!! Contains the \ref ModelPrecision module

!> \defgroup ModelPrecision ModelPrecision 
!! This module is used to set the precision of floating point numbers throughout the rest of the code.
!!
!! To set the precision to <B>single precision</B>, set <B>prec=sp</B> <BR>
!! To set the precision to <B>double precision</B>, set <B>prec=dp</B> <BR>
!! The precision must be set prior to compilation.

 
MODULE ModelPrecision

INTEGER, PARAMETER :: sp   = SELECTED_REAL_KIND(6, 37)     ! 32-bit
INTEGER, PARAMETER :: dp   = SELECTED_REAL_KIND(15, 307)   ! 64-bit
INTEGER, PARAMETER :: qp   = SELECTED_REAL_KIND(33, 4931)  ! 128-bit
INTEGER, PARAMETER :: prec = sp                            ! Specify the precision here

END MODULE ModelPrecision
