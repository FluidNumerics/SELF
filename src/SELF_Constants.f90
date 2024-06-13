! SELF_Constants.f90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Constants

  use iso_c_binding
  use iso_fortran_env

  implicit none
  include "mpif.h"

#include "SELF_Macros.h"

#ifdef DOUBLE_PRECISION
  integer,parameter :: prec = c_double
  integer,parameter :: c_prec = c_double
#else
  integer,parameter :: prec = c_float
  integer,parameter :: c_prec = c_float
#endif

!*************************************************************!
! ------------------ CHARACTER LENGTHS----- ------------------!
! ************************************************************!
!                                                             !
  integer,parameter :: SELF_EQN_DEFAULT_LENGTH = 100
  integer,parameter :: SELF_FILE_DEFAULT_LENGTH = 500

! ------------------------------------------------------------!
!*************************************************************!
! ------------------ MATHEMATICAL CONSTANTS ------------------!
! ************************************************************!
!                                                             !
! ------------------------------------------------------------!
  real(prec),parameter :: pi = 4.0_prec*atan(1.0_prec)
  real(prec),parameter :: TOL = epsilon(1.0_prec)

  real(prec),parameter :: fillValue = -9999.99_prec
  integer,parameter    :: fillValueInt = -99999

!*************************************************************!
! ----------------- ROOT FINDER CONSTANTS --------------------!
! ************************************************************!
!                                                             !
! ------------------------------------------------------------!
  integer,parameter    :: maxInverseIters = 1000
  real(prec),parameter :: newtonTolerance = 10.0**(-8)
  integer,parameter    :: newtonMax = 500

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
  real(prec),parameter   :: secondsToMinutes = 1.0_prec/60.0_prec ! conversion for seconds to minutes
  real(prec),parameter   :: minutesToHours = 1.0_prec/60.0_prec ! conversion for minutes to hours
  real(prec),parameter   :: hoursToDays = 1.0_prec/24.0_prec ! conversion for hours to days
  real(prec),parameter   :: daysToMonths = 12.0_prec/365.25_prec ! conversion for days to months
  real(prec),parameter   :: monthsToYears = 1.0_prec/12.0_prec ! conversion for months to years
  real(prec),parameter   :: daysToSeconds = 86400.0_prec

!==============================================!
! --------------- Quadrature------------------ !
!==============================================!
  integer,parameter :: GAUSS = 1
  integer,parameter :: GAUSS_LOBATTO = 2
  integer,parameter :: CHEBYSHEV_GAUSS = 3
  integer,parameter :: CHEBYSHEV_GAUSS_LOBATTO = 4
  integer,parameter :: UNIFORM = 5
  integer,parameter :: DG = 2000
  integer,parameter :: CG = 2001

! Misc. INTEGER and CHARACTER flag definitions
  character(1),parameter :: nada = ' '
  character(6),parameter :: MsgFmt = '(2x,A)'
  integer,parameter :: self_FileNameLength = 500
  integer,parameter :: self_TecplotHeaderLength = 500
  integer,parameter :: self_EquationLength = 210
  integer,parameter :: self_FormatLength = 30
  integer,parameter :: self_QuadratureTypeCharLength = 50
  integer,parameter :: self_IntegratorTypeCharLength = 50

endmodule SELF_Constants
