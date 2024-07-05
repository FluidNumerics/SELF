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

module SELF_Scalar_3D

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

  type,extends(SELF_DataObj),public :: Scalar3D

    real(prec),pointer,dimension(:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:) :: extBoundary
    real(prec),pointer,dimension(:,:,:,:,:) :: avgBoundary
    real(prec),pointer,dimension(:,:,:,:,:) :: jumpBoundary

    real(prec),pointer,dimension(:,:,:,:,:) :: interpWork1
    real(prec),pointer,dimension(:,:,:,:,:) :: interpWork2

  contains

    procedure,public :: Init => Init_Scalar3D
    procedure,public :: Free => Free_Scalar3D

    procedure,public :: BoundaryInterp => BoundaryInterp_Scalar3D
    procedure,public :: GridInterp => GridInterp_Scalar3D
    generic,public :: Gradient => Gradient_Scalar3D
    procedure,private :: Gradient_Scalar3D

    generic,public :: WriteHDF5 => WriteHDF5_MPI_Scalar3D,WriteHDF5_Scalar3D
    procedure,private :: WriteHDF5_MPI_Scalar3D
    procedure,private :: WriteHDF5_Scalar3D

  endtype Scalar3D

contains

  subroutine Init_Scalar3D(this,interp,nVar,nElem)
    implicit none
    class(Scalar3D),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar), &
             this%interpWork1(1:interp%M+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar), &
             this%interpWork2(1:interp%M+1,1:interp%M+1,1:interp%N+1,1:nelem,1:nvar), &
             this%boundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%extBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%avgBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%jumpBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar))

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % interpWork1)
    !$omp target enter data map(alloc: this % interpWork2)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)
    !$omp target enter data map(alloc: this % avgBoundary)
    !$omp target enter data map(alloc: this % jumpBoundary)

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:nVar))

  endsubroutine Init_Scalar3D

  subroutine Free_Scalar3D(this)
    implicit none
    class(Scalar3D),intent(inout) :: this

    this%nVar = 0
    this%nElem = 0
    this%interp => null()
    deallocate(this%interior)
    deallocate(this%interpWork1)
    deallocate(this%interpWork2)
    deallocate(this%boundary)
    deallocate(this%extBoundary)
    deallocate(this%avgBoundary)
    deallocate(this%jumpBoundary)
    deallocate(this%meta)
    deallocate(this%eqn)

    !$omp target exit data map(delete: this % interior)
    !$omp target exit data map(delete: this % interpWork1)
    !$omp target exit data map(delete: this % interpWork2)
    !$omp target exit data map(delete: this % boundary)
    !$omp target exit data map(delete: this % extBoundary)
    !$omp target exit data map(delete: this % avgBoundary)
    !$omp target exit data map(delete: this % jumpBoundary)

  endsubroutine Free_Scalar3D

  subroutine BoundaryInterp_Scalar3D(this)
    implicit none
    class(Scalar3D),intent(inout) :: this
    ! Local
    integer :: i,j,ii,iel,ivar
    real(prec) :: fb(1:6)

    !$omp target map(to:this%interior,this%interp%bMatrix) map(from:this%boundary)
    !$omp teams distribute parallel do collapse(4)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do j = 1,this%N+1
          do i = 1,this%N+1

            fb(1:6) = 0.0_prec

            do ii = 1,this%N+1
              fb(1) = fb(1)+this%interp%bMatrix(ii,1)*this%interior(i,j,ii,iel,ivar) ! Bottom
              fb(2) = fb(2)+this%interp%bMatrix(ii,1)*this%interior(i,ii,j,iel,ivar) ! South
              fb(3) = fb(3)+this%interp%bMatrix(ii,2)*this%interior(ii,i,j,iel,ivar) ! East
              fb(4) = fb(4)+this%interp%bMatrix(ii,2)*this%interior(i,ii,j,iel,ivar) ! North
              fb(5) = fb(5)+this%interp%bMatrix(ii,1)*this%interior(ii,i,j,iel,ivar) ! West
              fb(6) = fb(6)+this%interp%bMatrix(ii,2)*this%interior(i,j,ii,iel,ivar) ! Top
            enddo

            this%boundary(i,j,1:6,iel,ivar) = fb(1:6)

          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine BoundaryInterp_Scalar3D

  function GridInterp_Scalar3D(this) result(f)
    implicit none
    class(Scalar3D),intent(in) :: this
    real(prec) :: f(1:this%M+1,1:this%M+1,1:this%M+1,1:this%nelem,1:this%nvar)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: i,j,k,ii,jj,kk,iel,ivar
    real(prec) :: fi,fij,fijk

    !$omp target map(to:this%interior,this%interp%iMatrix) map(from:f)
    !$omp teams distribute parallel do collapse(5)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do k = 1,this%M+1
          do j = 1,this%M+1
            do i = 1,this%M+1

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
          enddo
        enddo
      enddo
    enddo
    !$omp end target

    ! call self_hipblas_matrixop_dim1_3d(this % iMatrix,f,fInt1,this % N,this % M,nvars,nelems,handle)
    ! call self_hipblas_matrixop_dim2_3d(this % iMatrix,fInt1,fInt2,0.0_c_prec,this % N,this % M,nvars,nelems,handle)
    ! call self_hipblas_matrixop_dim3_3d(this % iMatrix,fInt2,fTarget,0.0_c_prec,this % N,this % M,nvars,nelems,handle)


  endfunction GridInterp_Scalar3D

  function Gradient_Scalar3D(this) result(df)
    implicit none
    class(Scalar3D),intent(in) :: this
    real(prec) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:3)
    ! Local
    integer    :: i,j,k,ii,iel,ivar
    real(prec) :: df1,df2,df3

    !$omp target map(to:this%interior,this%interp%dMatrix) map(from:df)
    !$omp teams distribute parallel do collapse(5)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

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
          enddo
        enddo
      enddo
    enddo
    !$omp end target

    ! dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,1)
    ! call self_hipblas_matrixop_dim1_3d(this % dMatrix,f,dfloc,this % N,this % N,nvars,nelems,handle)
    ! dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,2)
    ! call self_hipblas_matrixop_dim2_3d(this % dMatrix,f,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)
    ! dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,3)
    ! call self_hipblas_matrixop_dim3_3d(this % dMatrix,f,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)
    ! dfloc => null()

  endfunction Gradient_Scalar3D

  subroutine WriteHDF5_MPI_Scalar3D(this,fileId,group,elemoffset,nglobalelem)
    implicit none
    class(Scalar3D),intent(in) :: this
    character(*),intent(in) :: group
    integer(HID_T),intent(in) :: fileId
    integer,intent(in) :: elemoffset
    integer,intent(in) :: nglobalelem
    ! Local
    integer(HID_T) :: offset(1:5)
    integer(HID_T) :: bOffset(1:5)
    integer(HID_T) :: globalDims(1:5)
    integer(HID_T) :: bGlobalDims(1:5)
    integer :: ivar

    offset(1:5) = (/0,0,0,0,elemoffset/)
    globalDims(1:5) = (/this%interp%N+1, &
                        this%interp%N+1, &
                        this%interp%N+1, &
                        this%nVar, &
                        nglobalelem/)

    ! Offsets and dimensions for element boundary data
    bOffset(1:5) = (/0,0,0,0,elemoffset/)
    bGlobalDims(1:5) = (/this%interp%N+1, &
                         this%interp%N+1, &
                         this%nVar, &
                         6, &
                         nglobalelem/)

    call CreateGroup_HDF5(fileId,trim(group))

    do ivar = 1,this%nVar
      call this%meta(ivar)%WriteHDF5(group,ivar,fileId)
    enddo

    call WriteArray_HDF5(fileId,trim(group)//"/interior", &
                         this%interior,offset,globalDims)

    call WriteArray_HDF5(fileId,trim(group)//"/boundary", &
                         this%boundary,bOffset,bGlobalDims)

  endsubroutine WriteHDF5_MPI_Scalar3D

  subroutine WriteHDF5_Scalar3D(this,fileId,group)
    implicit none
    class(Scalar3D),intent(in) :: this
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: group
    ! Local
    integer :: ivar

    call CreateGroup_HDF5(fileId,trim(group))

    do ivar = 1,this%nVar
      call this%meta(ivar)%WriteHDF5(group,ivar,fileId)
    enddo

    call WriteArray_HDF5(fileId,trim(group)//"/interior", &
                         this%interior)

    call WriteArray_HDF5(fileId,trim(group)//"/boundary", &
                         this%boundary)

  endsubroutine WriteHDF5_Scalar3D

endmodule SELF_Scalar_3D
