! SELF_Constants.f90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE SELF_Constants

  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  INCLUDE 'mpif.h'

#ifdef DOUBLE_PRECISION
  INTEGER,PARAMETER :: prec = real64
#else
  INTEGER,PARAMETER :: prec = real32
#endif

!*************************************************************!
! ------------------ MATHEMATICAL CONSTANTS ------------------!
! ************************************************************!
!                                                             !
! ------------------------------------------------------------!
  REAL(prec),PARAMETER :: pi = 4.0_prec*atan(1.0_prec)
  REAL(prec),PARAMETER :: TOL = epsilon(1.0_prec)

  REAL(prec),PARAMETER :: fillValue = -9999.99_prec
  INTEGER,PARAMETER    :: fillValueInt = -99999

!*************************************************************!
! ----------------- ROOT FINDER CONSTANTS --------------------!
! ************************************************************!
!                                                             !
! ------------------------------------------------------------!
  REAL(prec),PARAMETER :: tolerance = 10.0**(-10)
  INTEGER,PARAMETER    :: maxInverseIters = 1000
  REAL(prec),PARAMETER :: newtonTolerance = 10.0**(-8)
  INTEGER,PARAMETER    :: newtonMax = 500

!*************************************************************!
! ----------------- TIME STEPPING CONSTANTS ------------------!
! ************************************************************!
!                                                             !
! ------------------------------------------------------------!
! Runge-Kutta 3rd Order, low storage constants
  REAL(prec),PARAMETER :: rk3_a(1:3) = (/0.0_prec,-5.0_prec/9.0_prec,-153.0_prec/128.0_prec/)
  REAL(prec),PARAMETER :: rk3_b(1:3) = (/0.0_prec,1.0_prec/3.0_prec,3.0_prec/4.0_prec/)
  REAL(prec),PARAMETER :: rk3_g(1:3) = (/1.0_prec/3.0_prec,15.0_prec/16.0_prec,8.0_prec/15.0_prec/)

!*************************************************************!
! ------------------- PHYSICAL CONSTANTS ---------------------!
! ************************************************************!
!                                                             !
! ------------------------------------------------------------!
! Time conversion factors
  REAL(prec),PARAMETER   :: secondsToMinutes = 1.0_prec/60.0_prec                   ! conversion for seconds to minutes
  REAL(prec),PARAMETER   :: minutesToHours = 1.0_prec/60.0_prec                   ! conversion for minutes to hours
  REAL(prec),PARAMETER   :: hoursToDays = 1.0_prec/24.0_prec                   ! conversion for hours to days
  REAL(prec),PARAMETER   :: daysToMonths = 12.0_prec/365.25_prec                ! conversion for days to months
  REAL(prec),PARAMETER   :: monthsToYears = 1.0_prec/12.0_prec                   ! conversion for months to years
  REAL(prec),PARAMETER   :: daysToSeconds = 86400.0_prec

!==============================================!
! --------------- Quadrature------------------ !
!==============================================!
  INTEGER,PARAMETER :: GAUSS = 1
  INTEGER,PARAMETER :: GAUSS_LOBATTO = 2
  INTEGER,PARAMETER :: UNIFORM = 3
  INTEGER,PARAMETER :: DG = 2000
  INTEGER,PARAMETER :: CG = 2001

! Misc. INTEGER and CHARACTER flag definitions
  CHARACTER(1),PARAMETER :: nada = ' '
  CHARACTER(6),PARAMETER :: MsgFmt = '(2x,A)'
  INTEGER,PARAMETER :: self_FileNameLength = 500

END MODULE SELF_Constants
