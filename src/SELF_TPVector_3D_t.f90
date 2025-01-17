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

module SELF_TPVector_3D_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Metadata
  use FEQParse
  use SELF_HDF5
  use SELF_Data

  use HDF5
  use iso_c_binding

  implicit none

  type,extends(SELF_DataObj),public :: TPVector3D_t

    real(prec),pointer,contiguous,dimension(:,:,:,:,:,:,:) :: interior
    real(prec),pointer,contiguous,dimension(:,:,:,:,:) :: boundary
    real(prec),pointer,contiguous,dimension(:,:,:,:,:) :: riemann

  contains

    procedure,public :: Init => Init_TPVector3D_t
    procedure,public :: Free => Free_TPVector3D_t

    procedure,public :: UpdateHost => UpdateHost_TPVector3D_t
    procedure,public :: UpdateDevice => UpdateDevice_TPVector3D_t

  endtype TPVector3D_t

contains

  subroutine Init_TPVector3D_t(this,interp,nVar,nElem)
    implicit none
    class(TPVector3D_t),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! local
    integer :: i

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar,1:3), &
             this%boundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%riemann(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:3*nVar))

    ! Initialize equation parser
    ! This is done to prevent segmentation faults that arise
    ! when building with amdflang that are traced back to
    ! feqparse_functions.f90 : finalize routine
    ! When the equation parser is not initialized, the
    ! functions are not allocated, which I think are the
    ! source of the segfault - joe@fluidnumerics.com
    do i = 1,3*nvar
      this%eqn(i) = EquationParser('f=0',(/'x','y','z','t'/))
    enddo

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%riemann = 0.0_prec

  endsubroutine Init_TPVector3D_t

  subroutine Free_TPVector3D_t(this)
    implicit none
    class(TPVector3D_t),intent(inout) :: this

    this%interp => null()
    this%nVar = 0
    this%nElem = 0

    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%riemann)

    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_TPVector3D_t

  subroutine UpdateHost_TPVector3D_t(this)
    implicit none
    class(TPVector3D_t),intent(inout) :: this

  endsubroutine UpdateHost_TPVector3D_t

  subroutine UpdateDevice_TPVector3D_t(this)
    implicit none
    class(TPVector3D_t),intent(inout) :: this

  endsubroutine UpdateDevice_TPVector3D_t

endmodule SELF_TPVector_3D_t
