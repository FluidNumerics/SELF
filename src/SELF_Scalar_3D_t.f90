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

module SELF_Scalar_3D_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Metadata
  use FEQParse
  use SELF_HDF5
  use SELF_Data

  use HDF5
  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type,extends(SELF_DataObj),public :: Scalar3D_t

    real(prec),pointer,contiguous,dimension(:,:,:,:,:) :: interior
    real(prec),pointer,contiguous,dimension(:,:,:,:,:) :: boundary
    real(prec),pointer,contiguous,dimension(:,:,:,:,:) :: extBoundary
    real(prec),pointer,contiguous,dimension(:,:,:,:,:) :: avgBoundary
    real(prec),pointer,contiguous,dimension(:,:,:,:,:) :: boundarynormal

  contains

    procedure,public :: Init => Init_Scalar3D_t
    procedure,public :: Free => Free_Scalar3D_t

    procedure,public :: UpdateHost => UpdateHost_Scalar3D_t
    procedure,public :: UpdateDevice => UpdateDevice_Scalar3D_t

    procedure,public :: BoundaryInterp => BoundaryInterp_Scalar3D_t
    procedure,public :: AverageSides => AverageSides_Scalar3D_t
    generic,public :: GridInterp => GridInterp_Scalar3D_t
    procedure,private :: GridInterp_Scalar3D_t
    generic,public :: Gradient => Gradient_Scalar3D_t
    procedure,private :: Gradient_Scalar3D_t

    generic,public :: WriteHDF5 => WriteHDF5_MPI_Scalar3D_t,WriteHDF5_Scalar3D_t
    procedure,private :: WriteHDF5_MPI_Scalar3D_t
    procedure,private :: WriteHDF5_Scalar3D_t

  endtype Scalar3D_t

contains

  subroutine Init_Scalar3D_t(this,interp,nVar,nElem)
    implicit none
    class(Scalar3D_t),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar), &
             this%boundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%extBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%avgBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%boundarynormal(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:3*nvar))

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%extBoundary = 0.0_prec
    this%avgBoundary = 0.0_prec
    this%boundarynormal = 0.0_prec

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:nVar))

  endsubroutine Init_Scalar3D_t

  subroutine Free_Scalar3D_t(this)
    implicit none
    class(Scalar3D_t),intent(inout) :: this

    this%nVar = 0
    this%nElem = 0
    this%interp => null()
    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%extBoundary)
    deallocate(this%avgBoundary)
    deallocate(this%boundarynormal)
    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_Scalar3D_t

  subroutine UpdateHost_Scalar3D_t(this)
    implicit none
    class(Scalar3D_t),intent(inout) :: this

  endsubroutine UpdateHost_Scalar3D_t

  subroutine UpdateDevice_Scalar3D_t(this)
    implicit none
    class(Scalar3D_t),intent(inout) :: this

  endsubroutine UpdateDevice_Scalar3D_t

  subroutine BoundaryInterp_Scalar3D_t(this)
    implicit none
    class(Scalar3D_t),intent(inout) :: this
    ! Local
    integer :: i,j,ii,iel,ivar
    real(prec) :: fbb,fbs,fbe,fbn,fbw,fbt

    do concurrent(i=1:this%N+1,j=1:this%N+1, &
                  iel=1:this%nelem,ivar=1:this%nvar)

      fbb = 0.0_prec
      fbs = 0.0_prec
      fbe = 0.0_prec
      fbn = 0.0_prec
      fbw = 0.0_prec
      fbt = 0.0_prec

      do ii = 1,this%N+1
        fbb = fbb+this%interp%bMatrix(ii,1)*this%interior(i,j,ii,iel,ivar) ! Bottom
        fbs = fbs+this%interp%bMatrix(ii,1)*this%interior(i,ii,j,iel,ivar) ! South
        fbe = fbe+this%interp%bMatrix(ii,2)*this%interior(ii,i,j,iel,ivar) ! East
        fbn = fbn+this%interp%bMatrix(ii,2)*this%interior(i,ii,j,iel,ivar) ! North
        fbw = fbw+this%interp%bMatrix(ii,1)*this%interior(ii,i,j,iel,ivar) ! West
        fbt = fbt+this%interp%bMatrix(ii,2)*this%interior(i,j,ii,iel,ivar) ! Top
      enddo

      this%boundary(i,j,1,iel,ivar) = fbb
      this%boundary(i,j,2,iel,ivar) = fbs
      this%boundary(i,j,3,iel,ivar) = fbe
      this%boundary(i,j,4,iel,ivar) = fbn
      this%boundary(i,j,5,iel,ivar) = fbw
      this%boundary(i,j,6,iel,ivar) = fbt

    enddo

  endsubroutine BoundaryInterp_Scalar3D_t

  subroutine AverageSides_Scalar3D_t(this)
    implicit none
    class(Scalar3D_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i,j

    do concurrent(i=1:this%N+1,j=1:this%N+1, &
                  iside=1:6,iel=1:this%nelem,ivar=1:this%nvar)
      this%avgboundary(i,j,iside,iel,ivar) = 0.5_prec*( &
                                             this%boundary(i,j,iside,iel,ivar)+ &
                                             this%extBoundary(i,j,iside,iel,ivar))
    enddo

  endsubroutine AverageSides_Scalar3D_t

  subroutine GridInterp_Scalar3D_t(this,f)
    implicit none
    class(Scalar3D_t),intent(in) :: this
    real(prec),intent(out) :: f(1:this%M+1,1:this%M+1,1:this%M+1,1:this%nelem,1:this%nvar)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: i,j,k,ii,jj,kk,iel,ivar
    real(prec) :: fi,fij,fijk

    do concurrent(i=1:this%M+1,j=1:this%M+1, &
                  k=1:this%M+1,iel=1:this%nelem,ivar=1:this%nvar)

      fijk = 0.0_prec
      do kk = 1,this%N+1
        fij = 0.0_prec
        do jj = 1,this%N+1
          fi = 0.0_prec
          do ii = 1,this%N+1
            fi = fi+this%interior(ii,jj,kk,iel,ivar)*this%interp%iMatrix(ii,i)
          enddo
          fij = fij+fi*this%interp%iMatrix(jj,j)
        enddo
        fijk = fijk+fij*this%interp%iMatrix(kk,k)
      enddo
      f(i,j,k,iel,ivar) = fijk

    enddo

  endsubroutine GridInterp_Scalar3D_t

  subroutine Gradient_Scalar3D_t(this,df)
    implicit none
    class(Scalar3D_t),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:3)
    ! Local
    integer    :: i,j,k,ii,iel,ivar
    real(prec) :: df1,df2,df3

    do concurrent(i=1:this%N+1,j=1:this%N+1, &
                  k=1:this%N+1,iel=1:this%nelem,ivar=1:this%nvar)

      df1 = 0.0_prec
      df2 = 0.0_prec
      df3 = 0.0_prec
      do ii = 1,this%N+1
        df1 = df1+this%interp%dMatrix(ii,i)*this%interior(ii,j,k,iel,ivar)
        df2 = df2+this%interp%dMatrix(ii,j)*this%interior(i,ii,k,iel,ivar)
        df3 = df3+this%interp%dMatrix(ii,k)*this%interior(i,j,ii,iel,ivar)
      enddo
      df(i,j,k,iel,ivar,1) = df1
      df(i,j,k,iel,ivar,2) = df2
      df(i,j,k,iel,ivar,3) = df3

    enddo

  endsubroutine Gradient_Scalar3D_t

  subroutine WriteHDF5_MPI_Scalar3D_t(this,fileId,group,elemoffset,nglobalelem)
    implicit none
    class(Scalar3D_t),intent(in) :: this
    character(*),intent(in) :: group
    integer(HID_T),intent(in) :: fileId
    integer,intent(in) :: elemoffset
    integer,intent(in) :: nglobalelem
    ! Local
    integer(HID_T) :: offset(1:4)
    integer(HID_T) :: globalDims(1:4)
    integer :: ivar

    offset(1:4) = (/0,0,0,elemoffset/)
    globalDims(1:4) = (/this%interp%N+1, &
                        this%interp%N+1, &
                        this%interp%N+1, &
                        nglobalelem/)

    call CreateGroup_HDF5(fileId,trim(group))

    do ivar = 1,this%nVar
      !call this%meta(ivar)%WriteHDF5(group,ivar,fileId)
      call WriteArray_HDF5(fileId, &
                           trim(group)//"/"//trim(this%meta(ivar)%name), &
                           this%interior(:,:,:,:,ivar),offset,globalDims)
    enddo

  endsubroutine WriteHDF5_MPI_Scalar3D_t

  subroutine WriteHDF5_Scalar3D_t(this,fileId,group)
    implicit none
    class(Scalar3D_t),intent(in) :: this
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: group
    ! Local
    integer :: ivar

    call CreateGroup_HDF5(fileId,trim(group))

    do ivar = 1,this%nVar
      call this%meta(ivar)%WriteHDF5(group,ivar,fileId)
      call WriteArray_HDF5(fileId, &
                           trim(group)//trim(this%meta(ivar)%name), &
                           this%interior(:,:,:,:,ivar))
    enddo

  endsubroutine WriteHDF5_Scalar3D_t

endmodule SELF_Scalar_3D_t
