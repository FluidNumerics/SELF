! SELF_SupportRoutines.f90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file SELF_SupportRoutines.f90
!! Contains the \ref SELF_SupportRoutines module

!> \defgroup SELF_SupportRoutines SELF_SupportRoutines
!! This module defines a set of general purpose routines.

MODULE SELF_SupportRoutines

  USE ISO_FORTRAN_ENV
  USE SELF_Constants

  IMPLICIT NONE

  INTERFACE AlmostEqual
    MODULE PROCEDURE AlmostEqual_r64
  END INTERFACE AlmostEqual


  REAL(prec),PRIVATE,PARAMETER :: tolerance = 10.0**(-10)

CONTAINS

!> \addtogroup SELF_SupportRoutines
!! @{
! ================================================================================================ !
! Function AlmostEqual
!
!> \fn AlmostEqual
!! Compares two floating point numbers and determines if they are equal (to machine precision).
!!
!!   This function is from Alg. 139 on pg. 359 of D.A. Kopriva, 2009, "Implementing Spectral Element
!!    Methods for Scientists and Engineers"
!!
!! <H2> Usage : </H2>
!! <B>Logical</B> :: AisB <BR>
!! <B>REAL</B>(prec) :: a, b <BR>
!!         .... <BR>
!!     AisB = AlmostEqual( a, b ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in <th> a <td> REAL(prec) <td> scalar
!!   <tr> <td> in <th> b <td> REAL(prec) <td> scalar
!!   <tr> <td> in <th> AisB <td> Logical <td>
!!                     <B>.TRUE.</B> IF a=b to machine precision <BR>
!!                     <B>.FALSE.</B> otherwise
!!  </table>
!!
! ================================================================================================ !
!>@}

  FUNCTION AlmostEqual_r64(a,b) RESULT(AisB)

    IMPLICIT NONE
    REAL(real64) :: a,b
    LOGICAL :: AisB

    IF (a == 0.0_real64 .OR. b == 0.0_real64) THEN
      IF (ABS(a - b) <= EPSILON(1.0_real64)) THEN
        AisB = .TRUE.
      ELSE
        AisB = .FALSE.
      END IF
    ELSE
      IF ((abs(a - b) <= EPSILON(1.0_real64)*abs(a)) .OR. (abs(a - b) <= EPSILON(1.0_real64)*abs(b))) THEN
        AisB = .TRUE.
      ELSE
        AisB = .FALSE.
      END IF
    END IF

  END FUNCTION AlmostEqual_r64

!> \addtogroup SELF_SupportRoutines
!! @{
! ================================================================================================ !
! S/R ForwardShift
!
!> \fn ForwardShift
!! Shift an array integers by one index forward, moving the last index to the first.
!!
!! Shifts the array entries as follows : <BR>
!!  myArray(1) <-- myArray(N) <BR>
!!  myArray(2) <-- myArray(1) <BR>
!!  myArray(3) <-- myArray(2) <BR>
!!
!! <H2> Usage : </H2>
!! <B>INTEGER</B> :: N
!! <B>INTEGER</B> :: myArray(1:N) <BR>
!!         .... <BR>
!!     <B>CALL</B> ForwardShift( myArray, N ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myArray(1:N) <td> INTEGER <td>
!!                         On <B>output</B>, the input array with elements shifted forward by
!!                         one index.
!!   <tr> <td> in <th> N <td> INTEGER <td>
!!                     The number of elements in the array
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE ForwardShift(myArray,N)

    IMPLICIT NONE
    INTEGER,INTENT(in)    :: N
    INTEGER,INTENT(inout) :: myArray(1:N)
    ! LOCAL
    INTEGER :: temp(1:N)

    temp = myArray
    myArray(1) = temp(N)
    myArray(2:N) = temp(1:N - 1)

  END SUBROUTINE ForwardShift
!
!> \addtogroup SELF_SupportRoutines
!! @{
! ================================================================================================ !
! S/R CompareArray
!
!> \fn CompareArray
!! Compares to INTEGER arrays and determines if they are identical.
!!
!! A logical is returned that specifies whether or not two arrays are identical. To determine
!! if the two arrays are identical, the sum of the difference between each element in the input
!! array is calculated. If the arrays are identical, each contribution to the sum is zero and hence
!! the sum is zero. If the sum is non-zero, the arrays are distinct.
!!
!! This routine is used in the \ref HexMeshClass module. A face of an element in an unstructured
!! mesh is identified by its four corner nodes. When identifying unique faces in an unstructured
!! mesh, we need to determine if two elements share a face. This can be accomplished by comparing
!! the four corner nodes (from each element) that define each face.
!!
!! <H2> Usage : </H2>
!! <B>INTEGER</B> :: N <BR>
!! <B>INTEGER</B> :: arrayOne(1:N) <BR>
!! <B>INTEGER</B> :: arrayTwo(1:N) <BR>
!! <B>LOGICAL</B> :: arraysMatch <BR>
!!         .... <BR>
!!     arraysMatch = CompareArray( arrayOne, arrayTwo, N ) <BR>
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in <th> arrayOne(1:N) <td> INTEGER <td>
!!   <tr> <td> in <th> arrayTwo(1:N) <td> INTEGER <td>
!!   <tr> <td> in <th> N <td> INTEGER <td>
!!   <tr> <td> out <th> arraysMatch <td> INTEGER <td>
!!
!!  </table>
!!
! ================================================================================================ !
!>@}

  FUNCTION CompareArray(arrayOne,arrayTwo,N) RESULT(arraysMatch)

    IMPLICIT NONE
    INTEGER :: N
    INTEGER :: arrayOne(1:N),arrayTwo(1:N)
    LOGICAL :: arraysMatch
    ! LOCAL
    INTEGER :: i,theSumOfDiffs

    theSumOfDiffs = 0

    DO i = 1,N
      theSumOfDiffs = theSumOfDiffs + ABS(arrayOne(i) - arrayTwo(i))
    END DO

    IF (theSumOfDiffs == 0) THEN
      arraysMatch = .TRUE.
    ELSE
      arraysMatch = .FALSE.
    END IF

  END FUNCTION CompareArray
!
!> \addtogroup SELF_SupportRoutines
!! @{
! ================================================================================================ !
! S/R UniformPoints
!
!> \fn UniformPoints
!! Generates a REAL(prec) array of N points evenly spaced between two points.
!!
!!
!! <H2> Usage : </H2>
!! <B>REAL</B>(prec) :: a <BR>
!! <B>REAL</B>(prec) :: b <BR>
!! <B>REAL</B>(prec) :: xU(0:N) <BR>
!! <B>INTEGER</B> :: N <BR>
!!         .... <BR>
!!     xU = UniformPoints( a, b, N ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in <th> a <td> REAL(prec) <td> Starting point of the interval
!!   <tr> <td> in <th> b <td> REAL(prec) <td> Ending point of the interval
!!   <tr> <td> in <th> N <td> INTEGER <td> The number of points in the interval \f$[a,b]\f$
!!   <tr> <td> in <th> xU(0:N) <td> REAL(prec) <td>
!!                     Array of evenly spaced points in the interval \f$[a,b]\f$
!!  </table>
!!
! ================================================================================================ !
!>@}

  FUNCTION UniformPoints(a,b,firstInd,lastInd) RESULT(xU)

    IMPLICIT NONE
    REAL(prec) :: a,b
    INTEGER    :: firstInd,lastInd
    REAL(prec) :: xU(firstInd:lastInd)
    ! LOCAL
    REAL(prec)    :: dx
    INTEGER :: i

    dx = (b - a)/REAL((lastInd - firstInd),prec)

    DO i = firstInd,lastInd

      xU(i) = a + dx*REAL(i - firstInd,prec)

    END DO

  END FUNCTION UniformPoints
!
  FUNCTION TimeStamp(time,units) RESULT(timeStampString)
    IMPLICIT NONE
    REAL(prec)    :: time
    CHARACTER(1)  :: units
    CHARACTER(13) :: timeStampString
    ! Local
    INTEGER      :: day,minute,hour,second,millisecond
    CHARACTER(4) :: dayStamp
    CHARACTER(2) :: hourStamp,minuteStamp,secondStamp
    CHARACTER(3) :: milliSecondStamp
    REAL(real64) :: time_real64

    time_real64 = REAL(time,real64)
    ! Units in "seconds"
    IF (units(1:1) == 's') THEN

      ! Obtain the day
      day = INT(time_real64/86400.0_real64)
      hour = INT((time_real64 &
                  - 86400.0_real64*day)/3600.0_real64)
      minute = INT((time_real64 &
                    - 3600.0_real64*hour &
                    - 86400.0_real64*day)/60.0_real64)
      second = INT((time_real64 &
                    - 60.0_real64*minute &
                    - 3600.0_real64*hour &
                    - 86400.0_real64*day))
      milliSecond = NINT(((time_real64 &
                           - 60.0_real64*minute &
                           - 3600.0_real64*hour &
                           - 86400.0_real64*day) &
                          - REAL(second,real64))*1000.0_real64)

      IF( milliSecond >= 1000 )THEN
        milliSecond = milliSecond - 1000
        second = second + 1
      ENDIF

      IF( second >= 60 )THEN
        second = second - 60
        minute = minute + 1
      ENDIF

      IF( minute >= 60 )THEN
        minute = minute - 60
        hour = hour + 1
      ENDIF

      IF( hour >= 24 )THEN
        hour = hour - 24
        day = day + 1
      ENDIF

      WRITE (dayStamp,'(I4.4)') day
      WRITE (hourStamp,'(I2.2)') hour
      WRITE (minuteStamp,'(I2.2)') minute
      WRITE (secondStamp,'(I2.2)') second
      WRITE (milliSecondStamp,'(I3.3)') millisecond
      timeStampString = dayStamp//hourStamp//minuteStamp//secondStamp//milliSecondStamp

      ! minutes
    ELSEIF (units(1:1) == 'm') THEN

      ! hours
    ELSEIF (units(1:1) == 'h') THEN

    END IF

  END FUNCTION TimeStamp

  integer function newunit(unit)
  !  https://fortranwiki.org/fortran/show/newunit
  integer, intent(out), optional :: unit
! local
  integer, parameter :: LUN_MIN=10, LUN_MAX=1000
  logical :: opened
  integer :: lun
! begin
  newunit=-1
  do lun=LUN_MIN,LUN_MAX
    inquire(unit=lun,opened=opened)
    if (.not. opened) then
      newunit=lun
      exit
    end if
  end do
  if (present(unit)) unit=newunit
  end function newunit

  FUNCTION UpperCase(str) RESULT(upper)

    Implicit None
    CHARACTER(*),INTENT(In) :: str
    CHARACTER(LEN(str))      :: Upper

    INTEGER :: ic,i

    CHARACTER(27),PARAMETER :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ '
    CHARACTER(27),PARAMETER :: low = 'abcdefghijklmnopqrstuvwxyz '

    DO i = 1,LEN(str)
      ic = INDEX(low,str(i:i))
      IF (ic > 0) THEN
        Upper(i:i) = cap(ic:ic)
      ELSE
        Upper(i:i) = str(i:i)
      END IF
    END DO

  END FUNCTION UpperCase

END MODULE SELF_SupportRoutines
