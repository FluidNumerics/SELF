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

module SELF_Vector_3D

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

  type,extends(SELF_DataObj),public :: Vector3D

    real(prec),pointer,dimension(:,:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:,:) :: extBoundary
    real(prec),pointer,dimension(:,:,:,:,:,:) :: avgBoundary
    real(prec),pointer,dimension(:,:,:,:,:) :: boundaryNormal

  contains

    procedure,public :: Init => Init_Vector3D
    procedure,public :: Free => Free_Vector3D

    procedure,public :: BoundaryInterp => BoundaryInterp_Vector3D
    procedure,public :: GridInterp => GridInterp_Vector3D
    procedure,public :: Gradient => Gradient_Vector3D
    procedure,public :: Curl => Curl_Vector3D
    generic,public :: Divergence => Divergence_Vector3D
    procedure,private :: Divergence_Vector3D

    generic,public :: SetEquation => SetEquation_Vector3D
    procedure,private :: SetEquation_Vector3D

    generic,public :: WriteHDF5 => WriteHDF5_MPI_Vector3D,WriteHDF5_Vector3D
    procedure,private :: WriteHDF5_MPI_Vector3D
    procedure,private :: WriteHDF5_Vector3D

  endtype Vector3D

contains

  subroutine Init_Vector3D(this,interp,nVar,nElem)
    implicit none
    class(Vector3D),intent(out) :: this
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

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar,1:3), &
             this%boundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar,1:3), &
             this%extBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar,1:3), &
             this%avgBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar,1:3), &
             this%boundaryNormal(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar))

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

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)
    !$omp target enter data map(alloc: this % avgBoundary)
    !$omp target enter data map(alloc: this % boundaryNormal)

  endsubroutine Init_Vector3D

  subroutine Free_Vector3D(this)
    implicit none
    class(Vector3D),intent(inout) :: this

    this%interp => null()
    this%nVar = 0
    this%nElem = 0

    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%boundaryNormal)
    deallocate(this%extBoundary)
    deallocate(this%avgBoundary)

    !$omp target exit data map(delete: this % interior)
    !$omp target exit data map(delete: this % boundary)
    !$omp target exit data map(delete: this % extBoundary)
    !$omp target exit data map(delete: this % avgBoundary)
    !$omp target exit data map(delete: this % boundaryNormal)

    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_Vector3D

  subroutine SetEquation_Vector3D(this,idir,ivar,eqnChar)
    !! Sets the equation parser for the `idir` direction and `ivar-th` variable
    implicit none
    class(Vector3D),intent(inout) :: this
    integer,intent(in) :: idir,ivar
    character(*),intent(in) :: eqnChar

    this%eqn(idir+3*(ivar-1)) = EquationParser(trim(eqnChar), &
                                               (/'x','y','z','t'/))

  endsubroutine SetEquation_Vector3D

  function GridInterp_Vector3D(this) result(f)
    implicit none
    class(Vector3D),intent(in) :: this
    real(prec) :: f(1:this%M+1,1:this%M+1,1:this%M+1,1:this%nelem,1:this%nvar,1:3)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: i,j,k,ii,jj,kk,iel,ivar,idir
    real(prec) :: fi,fij,fijk

    !$omp target map(to:this%interior,this%interp%iMatrix) map(from:f)
    !$omp teams loop bind(teams) collapse(6)
    do idir = 1,3
      do ivar = 1,this%nvar
        do iel = 1,this%nelem
          do k = 1,this%M+1
            do j = 1,this%M+1
              do i = 1,this%M+1

                fijk = 0.0_prec
                !$omp loop bind(thread)
                do kk = 1,this%N+1
                  fij = 0.0_prec
                  !$omp loop bind(thread)
                  do jj = 1,this%N+1
                    fi = 0.0_prec
                    !$omp loop bind(thread)
                    do ii = 1,this%N+1
                      fi = fi+this%interior(ii,jj,kk,iel,ivar,idir)*this%interp%iMatrix(ii,i)
                    enddo
                    fij = fij+fi*this%interp%iMatrix(jj,j)
                  enddo
                  fijk = fijk+fij*this%interp%iMatrix(kk,k)
                enddo
                f(i,j,k,iel,ivar,idir) = fijk

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endfunction GridInterp_Vector3D

  subroutine BoundaryInterp_Vector3D(this)
    implicit none
    class(Vector3D),intent(inout) :: this
    ! Local
    integer :: i,j,ii,idir,iel,ivar
    real(prec) :: fbb,fbs,fbe,fbn,fbw,fbt

    !$omp target map(to:this%interior,this%interp%bMatrix) map(from:this%boundary)
    !$omp teams loop bind(teams) collapse(5)
    do idir = 1,3
      do ivar = 1,this%nvar
        do iel = 1,this%nelem
          do j = 1,this%N+1
            do i = 1,this%N+1

              fbb = 0.0_prec
              fbs = 0.0_prec
              fbe = 0.0_prec
              fbn = 0.0_prec
              fbw = 0.0_prec
              fbt = 0.0_prec
              !$omp loop bind(thread)
              do ii = 1,this%N+1
                fbb = fbb+this%interp%bMatrix(ii,1)*this%interior(i,j,ii,iel,ivar,idir) ! Bottom
                fbs = fbs+this%interp%bMatrix(ii,1)*this%interior(i,ii,j,iel,ivar,idir) ! South
                fbe = fbe+this%interp%bMatrix(ii,2)*this%interior(ii,i,j,iel,ivar,idir) ! East
                fbn = fbn+this%interp%bMatrix(ii,2)*this%interior(i,ii,j,iel,ivar,idir) ! North
                fbw = fbw+this%interp%bMatrix(ii,1)*this%interior(ii,i,j,iel,ivar,idir) ! West
                fbt = fbt+this%interp%bMatrix(ii,2)*this%interior(i,j,ii,iel,ivar,idir) ! Top
              enddo

              this%boundary(i,j,1,iel,ivar,idir) = fbb
              this%boundary(i,j,2,iel,ivar,idir) = fbs
              this%boundary(i,j,3,iel,ivar,idir) = fbe
              this%boundary(i,j,4,iel,ivar,idir) = fbn
              this%boundary(i,j,5,iel,ivar,idir) = fbw
              this%boundary(i,j,6,iel,ivar,idir) = fbt

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine BoundaryInterp_Vector3D

  function Gradient_Vector3D(this) result(df)
    implicit none
    class(Vector3D),intent(in) :: this
    real(prec) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:3,1:3)
    ! Local
    integer    :: i,j,k,ii,idir,iel,ivar
    real(prec) :: dfds1,dfds2,dfds3

    !$omp target map(to:this%interior,this%interp%dMatrix) map(from:df)
    !$omp teams loop bind(teams) collapse(6)
    do idir = 1,3
      do ivar = 1,this%nvar
        do iel = 1,this%nelem
          do k = 1,this%N+1
            do j = 1,this%N+1
              do i = 1,this%N+1

                dfds1 = 0.0_prec
                dfds2 = 0.0_prec
                dfds3 = 0.0_prec
                !$omp loop bind(thread)
                do ii = 1,this%N+1
                  dfds1 = dfds1+this%interp%dMatrix(ii,i)*this%interior(ii,j,k,iel,ivar,idir)
                  dfds2 = dfds2+this%interp%dMatrix(ii,j)*this%interior(i,ii,k,iel,ivar,idir)
                  dfds3 = dfds3+this%interp%dMatrix(ii,k)*this%interior(i,j,ii,iel,ivar,idir)
                enddo
                df(i,j,k,iel,ivar,idir,1) = dfds1
                df(i,j,k,iel,ivar,idir,2) = dfds2
                df(i,j,k,iel,ivar,idir,3) = dfds3

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

    !   do idir = 1,3
    !     floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,idir)
    !     dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,idir,1)
    !     call self_hipblas_matrixop_dim1_3d(this % dMatrix,floc,dfloc,this % N,this % N,nvars,nelems,handle)
    !     dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,idir,2)
    !     call self_hipblas_matrixop_dim2_3d(this % dMatrix,floc,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)
    !     dfloc(1:,1:,1:,1:,1:) => df(1:,1:,1:,1:,1:,idir,3)
    !     call self_hipblas_matrixop_dim3_3d(this % dMatrix,floc,dfloc,0.0_c_prec,this % N,this % N,nvars,nelems,handle)
    !   end do
    !   dfloc => null()

  endfunction Gradient_Vector3D

  function Curl_Vector3D(this) result(curlf)
    implicit none
    class(Vector3D),intent(in) :: this
    real(prec) :: curlf(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:3)
    ! Local
    integer    :: i,j,k,ii,idir,iel,ivar
    real(prec) :: dfds1,dfds2,dfds3
    real(prec) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:3,1:3)

    !$omp target map(to:this%interior,this%interp%dMatrix) map(from:df)
    !$omp teams
    !$omp loop bind(teams) collapse(6)
    do idir = 1,3
      do ivar = 1,this%nvar
        do iel = 1,this%nelem
          do k = 1,this%N+1
            do j = 1,this%N+1
              do i = 1,this%N+1

                dfds1 = 0.0_prec
                dfds2 = 0.0_prec
                dfds3 = 0.0_prec
                !$omp loop bind(thread)
                do ii = 1,this%N+1
                  dfds1 = dfds1+this%interp%dMatrix(ii,i)*this%interior(ii,j,k,iel,ivar,idir)
                  dfds2 = dfds2+this%interp%dMatrix(ii,j)*this%interior(i,ii,k,iel,ivar,idir)
                  dfds3 = dfds3+this%interp%dMatrix(ii,k)*this%interior(i,j,ii,iel,ivar,idir)
                enddo
                df(i,j,k,iel,ivar,idir,1) = dfds1
                df(i,j,k,iel,ivar,idir,2) = dfds2
                df(i,j,k,iel,ivar,idir,3) = dfds3

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    !$omp loop collapse(5)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

              curlf(i,j,k,iel,ivar,1) = (df(i,j,k,iel,ivar,3,2)-df(i,j,k,iel,ivar,2,3))
              curlf(i,j,k,iel,ivar,2) = (df(i,j,k,iel,ivar,1,3)-df(i,j,k,iel,ivar,3,1))
              curlf(i,j,k,iel,ivar,3) = (df(i,j,k,iel,ivar,2,1)-df(i,j,k,iel,ivar,1,2))

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target

  endfunction Curl_Vector3D

  function Divergence_Vector3D(this) result(df)
    implicit none
    class(Vector3D),intent(in) :: this
    real(prec) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer    :: i,j,k,ii,iel,ivar
    real(prec) :: dfLoc

    !$omp target map(to:this%interior,this%interp%dMatrix) map(from:df)
    !$omp teams
    !$omp loop bind(teams) collapse(5)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

              dfLoc = 0.0_prec
              !$omp loop bind(thread)
              do ii = 1,this%N+1
                dfLoc = dfLoc+this%interp%dMatrix(ii,i)*this%interior(ii,j,k,iel,ivar,1)
              enddo
              dF(i,j,k,iel,ivar) = dfLoc

            enddo
          enddo
        enddo
      enddo
    enddo

    !$omp loop bind(teams) collapse(5)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

              dfLoc = 0.0_prec
              !$omp loop bind(thread)
              do ii = 1,this%N+1
                dfLoc = dfLoc+this%interp%dMatrix(ii,j)*this%interior(i,ii,k,iel,ivar,2)
              enddo
              dF(i,j,k,iel,ivar) = dF(i,j,k,iel,ivar)+dfLoc

            enddo
          enddo
        enddo
      enddo
    enddo

    !$omp loop bind(teams) collapse(5)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

              dfLoc = 0.0_prec
              !$omp loop bind(thread)
              do ii = 1,this%N+1
                dfLoc = dfLoc+this%interp%dMatrix(ii,k)*this%interior(i,j,ii,iel,ivar,3)
              enddo
              dF(i,j,k,iel,ivar) = dF(i,j,k,iel,ivar)+dfLoc

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target

    ! ! local
    ! real(prec),pointer :: floc(:,:,:,:,:)

    ! floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,1)
    ! call self_hipblas_matrixop_dim1_3d(this % dMatrix,floc,df,this % N,this % N,nvars,nelems,handle)
    ! floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,2)
    ! call self_hipblas_matrixop_dim2_3d(this % dMatrix,floc,df,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    ! floc(1:,1:,1:,1:,1:) => f(1:,1:,1:,1:,1:,3)
    ! call self_hipblas_matrixop_dim2_3d(this % dMatrix,floc,df,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    ! floc => null()

  endfunction Divergence_Vector3D

  subroutine WriteHDF5_MPI_Vector3D(this,fileId,group,elemoffset,nglobalelem)
    implicit none
    class(Vector3D),intent(in) :: this
    character(*),intent(in) :: group
    integer(HID_T),intent(in) :: fileId
    integer,intent(in) :: elemoffset
    integer,intent(in) :: nglobalelem
    ! Local
    integer(HID_T) :: offset(1:6)
    integer(HID_T) :: bOffset(1:6)
    integer(HID_T) :: globalDims(1:6)
    integer(HID_T) :: bGlobalDims(1:6)
    integer :: ivar

    offset(1:6) = (/0,0,0,0,0,elemoffset/)
    globalDims(1:6) = (/3, &
                        this%interp%N+1, &
                        this%interp%N+1, &
                        this%interp%N+1, &
                        this%nVar, &
                        nglobalelem/)

    ! Offsets and dimensions for element boundary data
    bOffset(1:6) = (/0,0,0,0,0,elemoffset/)
    bGlobalDims(1:6) = (/3, &
                         this%interp%N+1, &
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

  endsubroutine WriteHDF5_MPI_Vector3D

  subroutine WriteHDF5_Vector3D(this,fileId,group)
    implicit none
    class(Vector3D),intent(in) :: this
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

  endsubroutine WriteHDF5_Vector3D

endmodule SELF_Vector_3D
