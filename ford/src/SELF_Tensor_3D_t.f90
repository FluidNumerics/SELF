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
  use SELF_DataPool

  use HDF5
  use iso_c_binding

  implicit none

  type,extends(SELF_DataObj),public :: Tensor3D_t

    real(prec),pointer,contiguous,dimension(:,:,:,:,:,:,:) :: interior
    real(prec),pointer,contiguous,dimension(:,:,:,:,:,:,:) :: boundary
    real(prec),pointer,contiguous,dimension(:,:,:,:,:,:,:) :: extBoundary

    !! High-water-mark backing store for the arrays above (AMR Stage 6b/6c; see SELF_DataPool).
    !! The public members are pointers remapped onto the leading part of these pools at the exact
    !! logical shape, so every extent, stride and shape() is unchanged while Resize can reuse the
    !! storage across an adaptation epoch. Nothing should index the pools directly.
    real(prec),pointer,contiguous :: pool_interior(:) => null()
    real(prec),pointer,contiguous :: pool_boundary(:) => null()
    real(prec),pointer,contiguous :: pool_extBoundary(:) => null()

  contains

    procedure,public :: Init => Init_Tensor3D_t
    procedure,public :: Resize => Resize_Tensor3D_t
    procedure,public :: MapArrays => MapArrays_Tensor3D_t
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

    call this%MapArrays(interp%N+1,nVar,nElem)

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

    ! The public arrays point into the pools rather than owning storage (Stage 6b/6c), so they
    ! are nullified and the pools released.
    this%interior => null()
    this%boundary => null()
    this%extBoundary => null()
    if(associated(this%pool_interior)) deallocate(this%pool_interior)
    if(associated(this%pool_boundary)) deallocate(this%pool_boundary)
    if(associated(this%pool_extBoundary)) deallocate(this%pool_extBoundary)

    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_Tensor3D_t

  subroutine MapArrays_Tensor3D_t(this,Np,nVar,nElem)
    !! Size the backing pools for (Np,nVar,nElem) and remap the public arrays onto them at the
    !! exact logical shape. See SELF_DataPool for why this is done with pointer remapping rather
    !! than by over-allocating the element dimension.
    implicit none
    class(Tensor3D_t),intent(inout) :: this
    integer,intent(in) :: Np
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer :: nInt,nBnd

    nInt = Np*Np*Np*nElem*nVar*3*3
    nBnd = Np*Np*6*nElem*nVar*3*3

    call EnsurePool(this%pool_interior,nInt)
    call EnsurePool(this%pool_boundary,nBnd)
    call EnsurePool(this%pool_extBoundary,nBnd)

    this%interior(1:Np,1:Np,1:Np,1:nElem,1:nVar,1:3,1:3) => this%pool_interior(1:nInt)
    this%boundary(1:Np,1:Np,1:6,1:nElem,1:nVar,1:3,1:3) => this%pool_boundary(1:nBnd)
    this%extBoundary(1:Np,1:Np,1:6,1:nElem,1:nVar,1:3,1:3) => this%pool_extBoundary(1:nBnd)

  endsubroutine MapArrays_Tensor3D_t

  subroutine Resize_Tensor3D_t(this,interp,nVar,nElem)
    !! Rebind a live object to a new element count, reusing existing storage when it fits (AMR
    !! Stage 6b/6c). Unlike Init - which is intent(out), so it resets the object, reallocates,
    !! zeroes and reconstructs the 9*nVar equation parsers - this preserves metadata and parsers
    !! (neither depends on nElem) and touches storage only when it must grow.
    !!
    !! All arrays are zeroed, exactly as Init leaves them. See Resize_Scalar3D_t for why that is
    !! required rather than merely tidy: a reused pool holds indeterminate values and the boundary
    !! arrays are not all rewritten before they are read.
    implicit none
    class(Tensor3D_t),intent(inout) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    if(nVar /= this%nVar) then
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

  endsubroutine Resize_Tensor3D_t

  subroutine UpdateHost_Tensor3D_t(this)
    implicit none
    class(Tensor3D_t),intent(inout) :: this
    if(.false.) this%N = this%N ! CPU stub; suppress unused-dummy-argument warning
  endsubroutine UpdateHost_Tensor3D_t

  subroutine UpdateDevice_Tensor3D_t(this)
    implicit none
    class(Tensor3D_t),intent(inout) :: this
    if(.false.) this%N = this%N ! CPU stub; suppress unused-dummy-argument warning
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
