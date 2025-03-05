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

module SELF_Vector_2D_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Metadata
  use FEQParse
  use SELF_HDF5
  use SELF_Data

  use HDF5
  use iso_c_binding

  implicit none

  type,extends(SELF_DataObj),public :: Vector2D_t

    real(prec),pointer,contiguous,dimension(:,:,:,:,:) :: interior
    real(prec),pointer,contiguous,dimension(:,:,:,:,:) :: boundary
    real(prec),pointer,contiguous,dimension(:,:,:,:,:) :: extBoundary
    real(prec),pointer,contiguous,dimension(:,:,:,:,:) :: avgBoundary
    real(prec),pointer,contiguous,dimension(:,:,:,:) :: boundaryNormal

  contains

    procedure,public :: Init => Init_Vector2D_t
    procedure,public :: Free => Free_Vector2D_t

    procedure,public :: UpdateHost => UpdateHost_Vector2D_t
    procedure,public :: UpdateDevice => UpdateDevice_Vector2D_t

    procedure,public :: BoundaryInterp => BoundaryInterp_Vector2D_t
    procedure,public :: AverageSides => AverageSides_Vector2D_t

    generic,public :: GridInterp => GridInterp_Vector2D_t
    procedure,private :: GridInterp_Vector2D_t

    generic,public :: Gradient => Gradient_Vector2D_t
    procedure,private :: Gradient_Vector2D_t

    generic,public :: Divergence => Divergence_Vector2D_t
    procedure,private :: Divergence_Vector2D_t

    generic,public :: SetEquation => SetEquation_Vector2D_t
    procedure,private :: SetEquation_Vector2D_t

    generic,public :: WriteHDF5 => WriteHDF5_MPI_Vector2D_t,WriteHDF5_Vector2D_t
    procedure,private :: WriteHDF5_MPI_Vector2D_t
    procedure,private :: WriteHDF5_Vector2D_t

  endtype Vector2D_t

contains

  subroutine Init_Vector2D_t(this,interp,nVar,nElem)
    implicit none
    class(Vector2D_t),intent(out) :: this
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

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:nelem,1:nvar,1:2), &
             this%boundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2), &
             this%extBoundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2), &
             this%avgBoundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2), &
             this%boundaryNormal(1:interp%N+1,1:4,1:nelem,1:nvar))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:2*nVar))

    ! Initialize equation parser
    ! This is done to prevent segmentation faults that arise
    ! when building with amdflang that are traced back to
    ! feqparse_functions.f90 : finalize routine
    ! When the equation parser is not initialized, the
    ! functions are not allocated, which I think are the
    ! source of the segfault - joe@fluidnumerics.com
    do i = 1,2*nvar
      this%eqn(i) = EquationParser('f=0',(/'x','y','z','t'/))
    enddo

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%boundarynormal = 0.0_prec
    this%extBoundary = 0.0_prec
    this%avgBoundary = 0.0_prec

  endsubroutine Init_Vector2D_t

  subroutine Free_Vector2D_t(this)
    implicit none
    class(Vector2D_t),intent(inout) :: this

    this%interp => null()
    this%nVar = 0
    this%nElem = 0

    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%boundaryNormal)
    deallocate(this%extBoundary)
    deallocate(this%avgBoundary)

    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_Vector2D_t

  subroutine UpdateHost_Vector2D_t(this)
    implicit none
    class(Vector2D_t),intent(inout) :: this

  endsubroutine UpdateHost_Vector2D_t

  subroutine UpdateDevice_Vector2D_t(this)
    implicit none
    class(Vector2D_t),intent(inout) :: this

  endsubroutine UpdateDevice_Vector2D_t

  subroutine SetEquation_Vector2D_t(this,idir,ivar,eqnChar)
    !! Sets the equation parser for the `idir` direction and `ivar-th` variable
    implicit none
    class(Vector2D_t),intent(inout) :: this
    integer,intent(in) :: idir,ivar
    character(*),intent(in) :: eqnChar

    this%eqn(idir+2*(ivar-1)) = EquationParser(trim(eqnChar), &
                                               (/'x','y','z','t'/))

  endsubroutine SetEquation_Vector2D_t

  subroutine GridInterp_Vector2D_t(this,f)
    implicit none
    class(Vector2D_t),intent(in) :: this
    real(prec),intent(out) :: f(1:this%M+1,1:this%M+1,1:this%nelem,1:this%nvar,1:2)
    ! Local
    integer :: i,j,ii,jj,iel,ivar,idir
    real(prec) :: fi,fij

    do concurrent(i=1:this%M+1,j=1:this%M+1,iel=1:this%nElem, &
                  ivar=1:this%nVar,idir=1:2)

      fij = 0.0_prec
      do jj = 1,this%N+1
        fi = 0.0_prec
        do ii = 1,this%N+1
          fi = fi+this%interior(ii,jj,iel,ivar,idir)*this%interp%iMatrix(ii,i)
        enddo
        fij = fij+fi*this%interp%iMatrix(jj,j)
      enddo
      f(i,j,iel,ivar,idir) = fij

    enddo

  endsubroutine GridInterp_Vector2D_t

  subroutine AverageSides_Vector2D_t(this)
    implicit none
    class(Vector2D_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i
    integer :: idir

    do concurrent(i=1:this%interp%N+1,iside=1:4,iel=1:this%nElem, &
                  ivar=1:this%nVar,idir=1:2)
      this%avgboundary(i,iside,iel,ivar,idir) = 0.5_prec*( &
                                                this%boundary(i,iside,iel,ivar,idir)+ &
                                                this%extBoundary(i,iside,iel,ivar,idir))
    enddo

  endsubroutine AverageSides_Vector2D_t

  subroutine BoundaryInterp_Vector2D_t(this)
    implicit none
    class(Vector2D_t),intent(inout) :: this
! Local
    integer :: i,ii,idir,iel,ivar
    real(prec) :: fbs,fbe,fbn,fbw

    do concurrent(i=1:this%N+1,iel=1:this%nelem, &
                  ivar=1:this%nvar,idir=1:2)

      fbs = 0.0_prec
      fbe = 0.0_prec
      fbn = 0.0_prec
      fbw = 0.0_prec
      do ii = 1,this%N+1
        fbs = fbs+this%interp%bMatrix(ii,1)*this%interior(i,ii,iel,ivar,idir) ! South
        fbe = fbe+this%interp%bMatrix(ii,2)*this%interior(ii,i,iel,ivar,idir) ! East
        fbn = fbn+this%interp%bMatrix(ii,2)*this%interior(i,ii,iel,ivar,idir) ! North
        fbw = fbw+this%interp%bMatrix(ii,1)*this%interior(ii,i,iel,ivar,idir) ! West
      enddo
      this%boundary(i,1,iel,ivar,idir) = fbs
      this%boundary(i,2,iel,ivar,idir) = fbe
      this%boundary(i,3,iel,ivar,idir) = fbn
      this%boundary(i,4,iel,ivar,idir) = fbw

    enddo

  endsubroutine BoundaryInterp_Vector2D_t

  subroutine Gradient_Vector2D_t(this,df)
    implicit none
    class(Vector2D_t),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:2,1:2)
    ! Local
    integer :: i,j,ii,iEl,iVar,idir
    real(prec) :: dfds1,dfds2

    do concurrent(i=1:this%N+1,j=1:this%N+1,iel=1:this%nElem, &
                  ivar=1:this%nVar,idir=1:2)

      dfds1 = 0.0_prec
      dfds2 = 0.0_prec
      do ii = 1,this%N+1
        dfds1 = dfds1+this%interp%dMatrix(ii,i)*this%interior(ii,j,iel,ivar,idir)
        dfds2 = dfds2+this%interp%dMatrix(ii,j)*this%interior(i,ii,iel,ivar,idir)
      enddo
      df(i,j,iel,ivar,idir,1) = dfds1
      df(i,j,iel,ivar,idir,2) = dfds2

    enddo

  endsubroutine Gradient_Vector2D_t

  subroutine Divergence_Vector2D_t(this,df)
    implicit none
    class(Vector2D_t),intent(in) :: this
    real(prec) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer    :: i,j,ii,iel,ivar
    real(prec) :: dfLoc

    do concurrent(i=1:this%N+1,j=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar)

      dfLoc = 0.0_prec
      do ii = 1,this%N+1
        dfLoc = dfLoc+this%interp%dMatrix(ii,i)*this%interior(ii,j,iel,ivar,1)+ &
                this%interp%dMatrix(ii,j)*this%interior(i,ii,iel,ivar,2)
      enddo
      dF(i,j,iel,ivar) = dfLoc

    enddo

  endsubroutine Divergence_Vector2D_t

  subroutine WriteHDF5_MPI_Vector2D_t(this,fileId,group,elemoffset,nglobalelem)
    implicit none
    class(Vector2D_t),intent(in) :: this
    character(*),intent(in) :: group
    integer(HID_T),intent(in) :: fileId
    integer,intent(in) :: elemoffset
    integer,intent(in) :: nglobalelem
    ! Local
    integer(HID_T) :: offset(1:3)
    integer(HID_T) :: globalDims(1:3)
    integer :: ivar,idir
    character(4) :: dimvar

    offset(1:3) = (/0,0,elemoffset/)
    globalDims(1:3) = (/this%interp%N+1, &
                        this%interp%N+1, &
                        nglobalelem/)

    call CreateGroup_HDF5(fileId,trim(group))

    do idir = 1,2
      write(dimvar,'(I1)') idir
      dimvar = "dim"//trim(dimvar)
      do ivar = 1,this%nVar
        !call this%meta(ivar)%WriteHDF5(group,ivar,fileId)
        call WriteArray_HDF5(fileId, &
                             trim(group)//"/"//trim(this%meta(ivar)%name)//"_"//dimvar, &
                             this%interior(:,:,:,ivar,idir),offset,globalDims)
      enddo
    enddo

  endsubroutine WriteHDF5_MPI_Vector2D_t

  subroutine WriteHDF5_Vector2D_t(this,fileId,group)
    implicit none
    class(Vector2D_t),intent(in) :: this
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: group
    ! Local
    integer :: ivar,idir
    character(4) :: dimvar

    call CreateGroup_HDF5(fileId,trim(group))

    do ivar = 1,this%nVar
      call this%meta(ivar)%WriteHDF5(group,ivar,fileId)
    enddo
    do idir = 1,2
      write(dimvar,'(I1)') idir
      dimvar = "dim"//trim(dimvar)
      do ivar = 1,this%nVar
        call WriteArray_HDF5(fileId, &
                             trim(group)//"/"//trim(this%meta(ivar)%name)//"_"//dimvar, &
                             this%interior(:,:,:,ivar,idir))
      enddo
    enddo

  endsubroutine WriteHDF5_Vector2D_t

endmodule SELF_Vector_2D_t
