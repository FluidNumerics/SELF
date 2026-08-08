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

module SELF_Scalar_2D_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Metadata
  use FEQParse
  use SELF_HDF5
  use SELF_Data
  use SELF_DataPool

  use HDF5
  use iso_c_binding

  implicit none

  type,extends(SELF_DataObj),public :: Scalar2D_t

    real(prec),pointer,contiguous,dimension(:,:,:,:) :: interior
    real(prec),pointer,contiguous,dimension(:,:,:,:) :: boundary
    real(prec),pointer,contiguous,dimension(:,:,:,:) :: extBoundary
    real(prec),pointer,contiguous,dimension(:,:,:,:) :: avgBoundary
    real(prec),pointer,contiguous,dimension(:,:,:,:) :: boundarynormal

    !! High-water-mark backing store for the arrays above (AMR Stage 6b; see SELF_DataPool).
    !! The public members are pointers remapped onto the leading part of these pools at the
    !! exact logical shape, so every extent, stride and shape() stays exactly as before while
    !! Resize can reuse the storage across an adaptation epoch. Nothing should index the pools.
    real(prec),pointer,contiguous :: pool_interior(:) => null()
    real(prec),pointer,contiguous :: pool_boundary(:) => null()
    real(prec),pointer,contiguous :: pool_extBoundary(:) => null()
    real(prec),pointer,contiguous :: pool_avgBoundary(:) => null()
    real(prec),pointer,contiguous :: pool_boundarynormal(:) => null()

  contains

    procedure,public :: Init => Init_Scalar2D_t
    procedure,public :: Resize => Resize_Scalar2D_t
    procedure,public :: MapArrays => MapArrays_Scalar2D_t
    procedure,public :: Free => Free_Scalar2D_t

    procedure,public :: UpdateHost => UpdateHost_Scalar2D_t
    procedure,public :: UpdateDevice => UpdateDevice_Scalar2D_t

    procedure,public :: BoundaryInterp => BoundaryInterp_Scalar2D_t
    procedure,public :: AverageSides => AverageSides_Scalar2D_t
    generic,public :: GridInterp => GridInterp_Scalar2D_t
    procedure,private :: GridInterp_Scalar2D_t
    generic,public :: Gradient => Gradient_Scalar2D_t
    procedure,private :: Gradient_Scalar2D_t

    generic,public :: WriteHDF5 => WriteHDF5_MPI_Scalar2D_t,WriteHDF5_Scalar2D_t
    procedure,private :: WriteHDF5_MPI_Scalar2D_t
    procedure,private :: WriteHDF5_Scalar2D_t

  endtype Scalar2D_t

contains

  subroutine Init_Scalar2D_t(this,interp,nVar,nElem)
    implicit none
    class(Scalar2D_t),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    call this%MapArrays(interp%N+1,nVar,nElem)

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%extBoundary = 0.0_prec
    this%avgBoundary = 0.0_prec
    this%boundarynormal = 0.0_prec

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:nVar))

  endsubroutine Init_Scalar2D_t

  subroutine MapArrays_Scalar2D_t(this,Np,nVar,nElem)
    !! Size the backing pools for (Np,nVar,nElem) and remap the public arrays onto them at the
    !! exact logical shape. Bounds-remapping a rank-1 contiguous target to a rank-4 pointer is
    !! what keeps every extent and stride identical to a plain allocate while allowing the pool
    !! underneath to be larger and reused. See SELF_DataPool.
    implicit none
    class(Scalar2D_t),intent(inout) :: this
    integer,intent(in) :: Np
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    call EnsurePool(this%pool_interior,Np*Np*nElem*nVar)
    call EnsurePool(this%pool_boundary,Np*4*nElem*nVar)
    call EnsurePool(this%pool_extBoundary,Np*4*nElem*nVar)
    call EnsurePool(this%pool_avgBoundary,Np*4*nElem*nVar)
    call EnsurePool(this%pool_boundarynormal,Np*4*nElem*2*nVar)

    this%interior(1:Np,1:Np,1:nElem,1:nVar) => this%pool_interior(1:Np*Np*nElem*nVar)
    this%boundary(1:Np,1:4,1:nElem,1:nVar) => this%pool_boundary(1:Np*4*nElem*nVar)
    this%extBoundary(1:Np,1:4,1:nElem,1:nVar) => this%pool_extBoundary(1:Np*4*nElem*nVar)
    this%avgBoundary(1:Np,1:4,1:nElem,1:nVar) => this%pool_avgBoundary(1:Np*4*nElem*nVar)
    this%boundarynormal(1:Np,1:4,1:nElem,1:2*nVar) => &
      this%pool_boundarynormal(1:Np*4*nElem*2*nVar)

    call EnsureIndexList(this%allElem,nElem)

  endsubroutine MapArrays_Scalar2D_t

  subroutine Resize_Scalar2D_t(this,interp,nVar,nElem)
    !! Rebind a live object to a new element count, reusing the existing storage when it fits
    !! (AMR Stage 6b). Unlike Init - which is intent(out) and therefore resets the whole object,
    !! reallocates, zeroes and, on GPU builds, uploads those zeros - this preserves metadata and
    !! equation parsers (neither depends on nElem) and touches storage only when it must grow.
    !!
    !! All arrays are zeroed, exactly as Init leaves them. This is not optional: a pool that has
    !! been reused (or freshly allocated at a larger capacity) holds indeterminate values, and
    !! the boundary arrays are not all fully rewritten before they are read - leaving them alone
    !! produced NaN entropy on the first step after an adaptation. The saving over Free + Init is
    !! the allocation itself, the device upload of zeros, and the metadata and equation-parser
    !! reconstruction; the zeroing is required for correctness and is kept.
    implicit none
    class(Scalar2D_t),intent(inout) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    if(nVar /= this%nVar) then
      ! Metadata and equation parsers are sized by nVar; a change means this is not a regrid.
      print*,__FILE__,':',__LINE__, &
        ' : Error : Resize cannot change nVar. Use Free followed by Init.'
      stop 1
    endif

    this%interp => interp
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    call this%MapArrays(interp%N+1,nVar,nElem)

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%extBoundary = 0.0_prec
    this%avgBoundary = 0.0_prec
    this%boundarynormal = 0.0_prec

  endsubroutine Resize_Scalar2D_t

  subroutine Free_Scalar2D_t(this)
    implicit none
    class(Scalar2D_t),intent(inout) :: this

    this%nVar = 0
    this%nElem = 0
    this%interp => null()
    ! The public arrays point into the pools rather than owning their storage, so they are
    ! nullified and the pools are released (AMR Stage 6b).
    this%interior => null()
    this%boundary => null()
    this%extBoundary => null()
    this%avgBoundary => null()
    this%boundarynormal => null()
    if(associated(this%pool_interior)) deallocate(this%pool_interior)
    if(associated(this%pool_boundary)) deallocate(this%pool_boundary)
    if(associated(this%pool_extBoundary)) deallocate(this%pool_extBoundary)
    if(associated(this%pool_avgBoundary)) deallocate(this%pool_avgBoundary)
    if(associated(this%pool_boundarynormal)) deallocate(this%pool_boundarynormal)
    if(associated(this%allElem)) deallocate(this%allElem)
    this%allElem => null()
    if(associated(this%allMortar)) deallocate(this%allMortar)
    this%allMortar => null()
    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_Scalar2D_t

  subroutine UpdateHost_Scalar2D_t(this)
    implicit none
    class(Scalar2D_t),intent(inout) :: this
    if(.false.) this%N = this%N ! CPU stub; suppress unused-dummy-argument warning
  endsubroutine UpdateHost_Scalar2D_t

  subroutine UpdateDevice_Scalar2D_t(this)
    implicit none
    class(Scalar2D_t),intent(inout) :: this
    if(.false.) this%N = this%N ! CPU stub; suppress unused-dummy-argument warning
  endsubroutine UpdateDevice_Scalar2D_t

  subroutine BoundaryInterp_Scalar2D_t(this,elems)
    !! Interpolate the element interior onto the four element edges.
    !!
    !! elems, when present, restricts the operation to those rank-local elements and leaves
    !! every other element's boundary trace untouched; local time stepping uses it to
    !! interpolate only the refinement level it is about to advance. Omitting it traverses
    !! 1:nElem exactly as before (see ResolveIndexList in SELF_Data).
    implicit none
    class(Scalar2D_t),intent(inout) :: this
    integer,pointer,contiguous,intent(in),optional :: elems(:)
    ! Local
    integer :: i,ii,ie,iel,ivar,ne
    integer,pointer,contiguous :: eidx(:)
    real(prec) :: fbs,fbe,fbn,fbw

    call ResolveIndexList(this%allElem,this%nElem,elems,eidx,ne)

    do concurrent(i=1:this%N+1,ie=1:ne,ivar=1:this%nvar)
      iel = eidx(ie)
      fbs = 0.0_prec
      fbe = 0.0_prec
      fbn = 0.0_prec
      fbw = 0.0_prec
      do ii = 1,this%N+1
        fbs = fbs+this%interp%bMatrix(ii,1)*this%interior(i,ii,iel,ivar) ! South
        fbe = fbe+this%interp%bMatrix(ii,2)*this%interior(ii,i,iel,ivar) ! East
        fbn = fbn+this%interp%bMatrix(ii,2)*this%interior(i,ii,iel,ivar) ! North
        fbw = fbw+this%interp%bMatrix(ii,1)*this%interior(ii,i,iel,ivar) ! West
      enddo

      this%boundary(i,1,iel,ivar) = fbs
      this%boundary(i,2,iel,ivar) = fbe
      this%boundary(i,3,iel,ivar) = fbn
      this%boundary(i,4,iel,ivar) = fbw

    enddo

  endsubroutine BoundaryInterp_Scalar2D_t

  subroutine AverageSides_Scalar2D_t(this,elems)
    !! Average the interior and exterior edge traces. elems restricts the operation to a
    !! subset of rank-local elements; omitting it covers 1:nElem as before.
    implicit none
    class(Scalar2D_t),intent(inout) :: this
    integer,pointer,contiguous,intent(in),optional :: elems(:)
    ! Local
    integer :: ie
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i
    integer :: ne
    integer,pointer,contiguous :: eidx(:)

    call ResolveIndexList(this%allElem,this%nElem,elems,eidx,ne)

    do concurrent(i=1:this%interp%N+1,iside=1:4,ie=1:ne,ivar=1:this%nVar)
      iel = eidx(ie)
      this%avgBoundary(i,iside,iel,ivar) = 0.5_prec*( &
                                           this%boundary(i,iside,iel,ivar)+ &
                                           this%extBoundary(i,iside,iel,ivar))
    enddo

  endsubroutine AverageSides_Scalar2D_t

  subroutine GridInterp_Scalar2D_t(this,f)
    implicit none
    class(Scalar2D_t),intent(in) :: this
    real(prec),intent(inout) :: f(1:this%M+1,1:this%M+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: i,j,ii,jj,iel,ivar
    real(prec) :: fi,fij

    do concurrent(i=1:this%M+1,j=1:this%M+1,iel=1:this%nElem,ivar=1:this%nVar)

      fij = 0.0_prec
      do jj = 1,this%N+1
        fi = 0.0_prec
        do ii = 1,this%N+1
          fi = fi+this%interior(ii,jj,iel,ivar)*this%interp%iMatrix(ii,i)
        enddo
        fij = fij+fi*this%interp%iMatrix(jj,j)
      enddo
      f(i,j,iel,ivar) = fij

    enddo

  endsubroutine GridInterp_Scalar2D_t

  subroutine Gradient_Scalar2D_t(this,df)
    implicit none
    class(Scalar2D_t),intent(in) :: this
    real(prec),intent(inout) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:2)
    ! Local
    integer    :: i,j,ii,iel,ivar
    real(prec) :: df1,df2

    do concurrent(i=1:this%N+1,j=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar)
      df1 = 0.0_prec
      df2 = 0.0_prec
      do ii = 1,this%N+1
        df1 = df1+this%interp%dMatrix(ii,i)*this%interior(ii,j,iel,ivar)
        df2 = df2+this%interp%dMatrix(ii,j)*this%interior(i,ii,iel,ivar)
      enddo
      df(i,j,iel,ivar,1) = df1
      df(i,j,iel,ivar,2) = df2
    enddo

  endsubroutine Gradient_Scalar2D_t

  subroutine WriteHDF5_MPI_Scalar2D_t(this,fileId,group,elemoffset,nglobalelem)
    implicit none
    class(Scalar2D_t),intent(in) :: this
    character(*),intent(in) :: group
    integer(HID_T),intent(in) :: fileId
    integer,intent(in) :: elemoffset
    integer,intent(in) :: nglobalelem
    ! Local
    integer(HID_T) :: offset(1:3)
    integer(HID_T) :: globalDims(1:3)
    integer :: ivar

    offset(1:3) = (/0,0,elemoffset/)
    globalDims(1:3) = (/this%interp%N+1, &
                        this%interp%N+1, &
                        nglobalelem/)

    call CreateGroup_HDF5(fileId,trim(group))

    do ivar = 1,this%nVar
      !call this%meta(ivar)%WriteHDF5(group,ivar,fileId)
      call WriteArray_HDF5(fileId, &
                           trim(group)//"/"//trim(this%meta(ivar)%name), &
                           this%interior(:,:,:,ivar),offset,globalDims)
    enddo

  endsubroutine WriteHDF5_MPI_Scalar2D_t

  subroutine WriteHDF5_Scalar2D_t(this,fileId,group)
    implicit none
    class(Scalar2D_t),intent(in) :: this
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: group
    ! Local
    integer :: ivar

    call CreateGroup_HDF5(fileId,trim(group))

    do ivar = 1,this%nVar
      call this%meta(ivar)%WriteHDF5(group,ivar,fileId)
      call WriteArray_HDF5(fileId, &
                           trim(group)//"/"//trim(this%meta(ivar)%name), &
                           this%interior(:,:,:,ivar))
    enddo

  endsubroutine WriteHDF5_Scalar2D_t

endmodule SELF_Scalar_2D_t
