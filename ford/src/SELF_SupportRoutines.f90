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

!> \file SELF_SupportRoutines.f90
!! Contains the \ref SELF_SupportRoutines module

!> \defgroup SELF_SupportRoutines SELF_SupportRoutines
!! This module defines a set of general purpose routines.

module SELF_SupportRoutines

  use iso_fortran_env
  use SELF_Constants

  implicit none

  interface AlmostEqual
    module procedure AlmostEqual_r64
  endinterface AlmostEqual

  real(prec),private,parameter :: tolerance = 10.0**(-10)

contains

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

  function AlmostEqual_r64(a,b) result(AisB)

    implicit none
    real(real64) :: a,b
    logical :: AisB

    if(a == 0.0_real64 .or. b == 0.0_real64) then
      if(abs(a-b) <= epsilon(1.0_real64)) then
        AisB = .true.
      else
        AisB = .false.
      endif
    else
      if((abs(a-b) <= epsilon(1.0_real64)*abs(a)) .or. (abs(a-b) <= epsilon(1.0_real64)*abs(b))) then
        AisB = .true.
      else
        AisB = .false.
      endif
    endif

  endfunction AlmostEqual_r64

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
  subroutine ForwardShift(myArray,N)

    implicit none
    integer,intent(in)    :: N
    integer,intent(inout) :: myArray(1:N)
    ! LOCAL
    integer :: temp(1:N)

    temp = myArray
    myArray(1) = temp(N)
    myArray(2:N) = temp(1:N-1)

  endsubroutine ForwardShift
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

  function CompareArray(arrayOne,arrayTwo,N) result(arraysMatch)

    implicit none
    integer :: N
    integer :: arrayOne(1:N),arrayTwo(1:N)
    logical :: arraysMatch
    ! LOCAL
    integer :: i,theSumOfDiffs

    theSumOfDiffs = 0

    do i = 1,N
      theSumOfDiffs = theSumOfDiffs+abs(arrayOne(i)-arrayTwo(i))
    enddo

    if(theSumOfDiffs == 0) then
      arraysMatch = .true.
    else
      arraysMatch = .false.
    endif

  endfunction CompareArray
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

  function UniformPoints(a,b,firstInd,lastInd) result(xU)

    implicit none
    real(prec) :: a,b
    integer    :: firstInd,lastInd
    real(prec) :: xU(firstInd:lastInd)
    ! LOCAL
    real(prec)    :: dx
    integer :: i

    dx = (b-a)/real((lastInd-firstInd),prec)

    do i = firstInd,lastInd

      xU(i) = a+dx*real(i-firstInd,prec)

    enddo

  endfunction UniformPoints

  integer function newunit(unit)
    !  https://fortranwiki.org/fortran/show/newunit
    integer,intent(out),optional :: unit
! local
    integer,parameter :: LUN_MIN = 10,LUN_MAX = 1000
    logical :: opened
    integer :: lun
! begin
    newunit = -1
    do lun = LUN_MIN,LUN_MAX
      inquire(unit=lun,opened=opened)
      if(.not. opened) then
        newunit = lun
        exit
      endif
    enddo
    if(present(unit)) unit = newunit
  endfunction newunit

  function UpperCase(str) result(upper)

    implicit none
    character(*),intent(In) :: str
    character(len(str))      :: Upper

    integer :: ic,i

    character(27),parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ '
    character(27),parameter :: low = 'abcdefghijklmnopqrstuvwxyz '

    do i = 1,len(str)
      ic = index(low,str(i:i))
      if(ic > 0) then
        Upper(i:i) = cap(ic:ic)
      else
        Upper(i:i) = str(i:i)
      endif
    enddo

  endfunction UpperCase

endmodule SELF_SupportRoutines
