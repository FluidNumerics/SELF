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

module SELF_Tensor_3D_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Metadata
  use FEQParse
  use SELF_HDF5
  use SELF_Data

  use HDF5
  use iso_c_binding

  implicit none

  type,extends(SELF_DataObj),public :: Tensor3D_t

    real(prec),pointer,contiguous,dimension(:,:,:,:,:,:,:) :: interior
    real(prec),pointer,contiguous,dimension(:,:,:,:,:,:,:) :: boundary
    real(prec),pointer,contiguous,dimension(:,:,:,:,:,:,:) :: extBoundary

  contains

    procedure,public :: Init => Init_Tensor3D_t
    procedure,public :: Free => Free_Tensor3D_t

    procedure,public :: BoundaryInterp => BoundaryInterp_Tensor3D_t

    procedure,public :: UpdateHost => UpdateHost_Tensor3D_t
    procedure,public :: UpdateDevice => UpdateDevice_Tensor3D_t

    generic,public :: Determinant => Determinant_Tensor3D_t
    procedure,private :: Determinant_Tensor3D_t

  endtype Tensor3D_t

contains

  subroutine Init_Tensor3D_t(this,interp,nVar,nElem)
    implicit none
    class(Tensor3D_t),intent(out) :: this
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

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar,1:3,1:3), &
             this%boundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar,1:3,1:3), &
             this%extBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar,1:3,1:3))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:9*nVar))

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%extBoundary = 0.0_prec

    ! Initialize equation parser
    ! This is done to prevent segmentation faults that arise
    ! when building with amdflang that are traced back to
    ! feqparse_functions.f90 : finalize routine
    ! When the equation parser is not initialized, the
    ! functions are not allocated, which I think are the
    ! source of the segfault - joe@fluidnumerics.com
    do i = 1,9*nvar
      this%eqn(i) = EquationParser('f=0',(/'x','y','z','t'/))
    enddo

  endsubroutine Init_Tensor3D_t

  subroutine Free_Tensor3D_t(this)
    implicit none
    class(Tensor3D_t),intent(inout) :: this

    this%interp => null()
    this%nVar = 0
    this%nElem = 0

    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%extBoundary)

    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_Tensor3D_t

  subroutine UpdateHost_Tensor3D_t(this)
    implicit none
    class(Tensor3D_t),intent(inout) :: this

  endsubroutine UpdateHost_Tensor3D_t

  subroutine UpdateDevice_Tensor3D_t(this)
    implicit none
    class(Tensor3D_t),intent(inout) :: this

  endsubroutine UpdateDevice_Tensor3D_t

  subroutine BoundaryInterp_Tensor3D_t(this)
    implicit none
    class(Tensor3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,ii,idir,jdir,iel,ivar
    real(prec) :: fbb,fbs,fbe,fbn,fbw,fbt

    do concurrent(i=1:this%N+1,j=1:this%N+1, &
                  iel=1:this%nelem,ivar=1:this%nvar, &
                  idir=1:3,jdir=1:3)

      fbb = 0.0_prec
      fbs = 0.0_prec
      fbe = 0.0_prec
      fbn = 0.0_prec
      fbw = 0.0_prec
      fbt = 0.0_prec

      do ii = 1,this%N+1
        fbb = fbb+this%interp%bMatrix(ii,1)*this%interior(i,j,ii,iel,ivar,idir,jdir) ! Bottom
        fbs = fbs+this%interp%bMatrix(ii,1)*this%interior(i,ii,j,iel,ivar,idir,jdir) ! South
        fbe = fbe+this%interp%bMatrix(ii,2)*this%interior(ii,i,j,iel,ivar,idir,jdir) ! East
        fbn = fbn+this%interp%bMatrix(ii,2)*this%interior(i,ii,j,iel,ivar,idir,jdir) ! North
        fbw = fbw+this%interp%bMatrix(ii,1)*this%interior(ii,i,j,iel,ivar,idir,jdir) ! West
        fbt = fbt+this%interp%bMatrix(ii,2)*this%interior(i,j,ii,iel,ivar,idir,jdir) ! Top
      enddo

      this%boundary(i,j,1,iel,ivar,idir,jdir) = fbb
      this%boundary(i,j,2,iel,ivar,idir,jdir) = fbs
      this%boundary(i,j,3,iel,ivar,idir,jdir) = fbe
      this%boundary(i,j,4,iel,ivar,idir,jdir) = fbn
      this%boundary(i,j,5,iel,ivar,idir,jdir) = fbw
      this%boundary(i,j,6,iel,ivar,idir,jdir) = fbt

    enddo

  endsubroutine BoundaryInterp_Tensor3D_t

  subroutine Determinant_Tensor3D_t(this,det)
    implicit none
    class(Tensor3D_t),intent(in) :: this
    real(prec),intent(out) :: det(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,j,k

    do concurrent(i=1:this%N+1,j=1:this%N+1, &
                  k=1:this%N+1,iel=1:this%nelem,ivar=1:this%nvar)

      det(i,j,k,iEl,iVar) = &
        this%interior(i,j,k,iEl,iVar,1,1)* &
        (this%interior(i,j,k,iEl,iVar,2,2)* &
         this%interior(i,j,k,iEl,iVar,3,3)- &
         this%interior(i,j,k,iEl,iVar,2,3)* &
         this%interior(i,j,k,iEl,iVar,3,2))- &
        this%interior(i,j,k,iEl,iVar,2,1)* &
        (this%interior(i,j,k,iEl,iVar,1,2)* &
         this%interior(i,j,k,iEl,iVar,3,3)- &
         this%interior(i,j,k,iEl,iVar,1,3)* &
         this%interior(i,j,k,iEl,iVar,3,2))+ &
        this%interior(i,j,k,iEl,iVar,3,1)* &
        (this%interior(i,j,k,iEl,iVar,1,2)* &
         this%interior(i,j,k,iEl,iVar,2,3)- &
         this%interior(i,j,k,iEl,iVar,1,3)* &
         this%interior(i,j,k,iEl,iVar,2,2))

    enddo

  endsubroutine Determinant_Tensor3D_t

endmodule SELF_Tensor_3D_t
