! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2024 Fluid Numerics LLC
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
!    the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Constants

  use iso_c_binding
  use iso_fortran_env

  implicit none
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
