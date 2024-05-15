! SELF_Constants.f90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE SELF_Constants

  USE ISO_C_BINDING
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE
  INCLUDE "mpif.h"

#include "SELF_Macros.h"

#ifdef DOUBLE_PRECISION
  INTEGER,PARAMETER :: prec = c_Double
  INTEGER,PARAMETER :: c_prec = c_Double
#else
  INTEGER,PARAMETER :: prec = c_float
  INTEGER,PARAMETER :: c_prec = C_float
#endif

!*************************************************************!
! ------------------ CHARACTER LENGTHS----- ------------------!
! ************************************************************!
!                                                             !
INTEGER, PARAMETER :: SELF_EQN_DEFAULT_LENGTH=100
INTEGER, PARAMETER :: SELF_FILE_DEFAULT_LENGTH=500

! ------------------------------------------------------------!
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
  INTEGER,PARAMETER    :: maxInverseIters = 1000
  REAL(prec),PARAMETER :: newtonTolerance = 10.0**(-8)
  INTEGER,PARAMETER    :: newtonMax = 500

!*************************************************************!
! ----------------- TIME STEPPING CONSTANTS ------------------!
! ************************************************************!
!                                                             !
! ------------------------------------------------------------!

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
  INTEGER,PARAMETER :: CHEBYSHEV_GAUSS = 3
  INTEGER,PARAMETER :: CHEBYSHEV_GAUSS_LOBATTO = 4
  INTEGER,PARAMETER :: UNIFORM = 5
  INTEGER,PARAMETER :: DG = 2000
  INTEGER,PARAMETER :: CG = 2001

! Misc. INTEGER and CHARACTER flag definitions
  CHARACTER(1),PARAMETER :: nada = ' '
  CHARACTER(6),PARAMETER :: MsgFmt = '(2x,A)'
  INTEGER,PARAMETER :: self_FileNameLength = 500
  INTEGER,PARAMETER :: self_TecplotHeaderLength = 500
  INTEGER,PARAMETER :: self_EquationLength = 210
  INTEGER,PARAMETER :: self_FormatLength=30
  INTEGER,PARAMETER :: self_QuadratureTypeCharLength = 50
  INTEGER,PARAMETER :: self_IntegratorTypeCharLength = 50

END MODULE SELF_Constants
