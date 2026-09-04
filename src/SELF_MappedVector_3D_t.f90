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

module SELF_MappedVector_3D_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_DomainDecomposition
  use FEQParse
  use iso_c_binding

  implicit none

  type,extends(Vector3D),public :: MappedVector3D_t

    ! Mortar exchange work array, allocated on first use for meshes with 2:1
    ! nonconforming interfaces; same slot layout as MappedScalar3D_t%mortarBuff with a
    ! trailing direction index. MortarFluxCollect stages the boundaryNormal integrand
    ! in direction slot 1.
    real(prec),allocatable,dimension(:,:,:,:,:,:) :: mortarBuff

    logical :: geometry_associated = .false.
    type(SEMHex),pointer :: geometry => null()
  contains

    procedure,public :: AssociateGeometry => AssociateGeometry_MappedVector3D_t
    procedure,public :: DissociateGeometry => DissociateGeometry_MappedVector3D_t

    procedure,public :: SideExchange => SideExchange_MappedVector3D_t
    procedure,public :: MortarExchange => MortarExchange_MappedVector3D_t
    procedure,private :: MPIMortarExchangeAsync => MPIMortarExchangeAsync_MappedVector3D_t
    procedure,public :: MortarFluxCollect => MortarFluxCollect_MappedVector3D_t
    procedure,private :: MPIMortarFluxAsync => MPIMortarFluxAsync_MappedVector3D_t

    generic,public :: MappedDivergence => MappedDivergence_MappedVector3D_t
    procedure,private :: MappedDivergence_MappedVector3D_t

    generic,public :: MappedDGDivergence => MappedDGDivergence_MappedVector3D_t
    procedure,private :: MappedDGDivergence_MappedVector3D_t

    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedVector3D_t
    procedure,private :: ApplyFlip => ApplyFlip_MappedVector3D_t

    procedure,public :: Resize => Resize_MappedVector3D_t

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedVector3D_t

    !procedure,public :: WriteTecplot => WriteTecplot_MappedVector3D_t

  endtype MappedVector3D_t

contains

  subroutine Resize_MappedVector3D_t(this,interp,nVar,nElem)
    !! Rebind to a new element count, reusing storage where it fits (AMR Stage 6b). See
    !! Resize_MappedScalar3D_t for why the mortar staging buffer is invalidated rather than resized.
    implicit none
    class(MappedVector3D_t),intent(inout) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    call Resize_Vector3D_t(this,interp,nVar,nElem)
    if(allocated(this%mortarBuff)) deallocate(this%mortarBuff)

    ! The geometry binding refers to the PREVIOUS mesh's geometry, which the caller is
    ! about to destroy. It must be dropped here so the AssociateGeometry that follows a
    ! regrid actually rebinds: AssociateGeometry is a no-op when a geometry is already
    ! associated, so leaving the stale one in place silently kept every mapped operation
    ! reading freed metric terms - observed as NaN entropy on the first step after an
    ! adaptation. Free + Init used to hide this by nulling the pointer.
    call this%DissociateGeometry()

  endsubroutine Resize_MappedVector3D_t

  subroutine AssociateGeometry_MappedVector3D_t(this,geometry)
    implicit none
    class(MappedVector3D_t),intent(inout) :: this
    type(SEMHex),target,intent(in) :: geometry

    if(.not. associated(this%geometry)) then
      this%geometry => geometry
      this%geometry_associated = .true.
    endif

  endsubroutine AssociateGeometry_MappedVector3D_t

  subroutine DissociateGeometry_MappedVector3D_t(this)
    implicit none
    class(MappedVector3D_t),intent(inout) :: this

    if(associated(this%geometry)) then
      this%geometry => null()
      this%geometry_associated = .false.
    endif

  endsubroutine DissociateGeometry_MappedVector3D_t

  subroutine SetInteriorFromEquation_MappedVector3D_t(this,geometry,time)
    !!  Sets the this % interior attribute using the eqn attribute,
    !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedVector3D_t),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,j,k,iEl,iVar
    real(prec) :: x
    real(prec) :: y
    real(prec) :: z

    do iVar = 1,this%nVar
      do iEl = 1,this%nElem
        do k = 1,this%interp%N+1
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1

              ! Get the mesh positions
              x = geometry%x%interior(i,j,k,iEl,1,1)
              y = geometry%x%interior(i,j,k,iEl,1,2)
              z = geometry%x%interior(i,j,k,iEl,1,3)

              this%interior(i,j,k,iEl,iVar,1) = &
                this%eqn(1+3*(iVar-1))%Evaluate((/x,y,z,time/))

              this%interior(i,j,k,iEl,iVar,2) = &
                this%eqn(2+3*(iVar-1))%Evaluate((/x,y,z,time/))

              this%interior(i,j,k,iEl,iVar,3) = &
                this%eqn(3+3*(iVar-1))%Evaluate((/x,y,z,time/))

            enddo
          enddo
        enddo
      enddo
    enddo

  endsubroutine SetInteriorFromEquation_MappedVector3D_t

  subroutine MPIExchangeAsync_MappedVector3D_t(this,mesh)
    implicit none
    class(MappedVector3D_t),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: e1,s1,e2,s2,ivar,idir
    integer :: globalSideId,r2,tag
    integer :: iError
    integer :: msgCount

    ! A send and a receive per rank-remote face, per variable and direction.
    call mesh%decomp%ReserveMessages(2*3*this%nvar*mesh%decomp%CountRemoteSides(mesh%sideInfo,mesh%nElem,6))

    msgCount = 0

    do idir = 1,3
      do ivar = 1,this%nvar
        do e1 = 1,this%nElem
          do s1 = 1,6

            e2 = mesh%sideInfo(3,s1,e1) ! Neighbor Element
            if(e2 > 0) then
              r2 = mesh%decomp%elemToRank(e2) ! Neighbor Rank

              if(r2 /= mesh%decomp%rankId) then

                s2 = mesh%sideInfo(4,s1,e1)/10
                globalSideId = abs(mesh%sideInfo(2,s1,e1))
                tag = globalsideid+mesh%nUniqueSides*(ivar-1+this%nvar*(idir-1))

                msgCount = msgCount+1
                call MPI_IRECV(this%extBoundary(:,:,s1,e1,ivar,idir), &
                               (this%interp%N+1)*(this%interp%N+1), &
                               mesh%decomp%mpiPrec, &
                               r2,tag, &
                               mesh%decomp%mpiComm, &
                               mesh%decomp%requests(msgCount),iError)

                msgCount = msgCount+1
                call MPI_ISEND(this%boundary(:,:,s1,e1,ivar,idir), &
                               (this%interp%N+1)*(this%interp%N+1), &
                               mesh%decomp%mpiPrec, &
                               r2,tag, &
                               mesh%decomp%mpiComm, &
                               mesh%decomp%requests(msgCount),iError)
              endif
            endif

          enddo
        enddo
      enddo
    enddo

    mesh%decomp%msgCount = msgCount

  endsubroutine MPIExchangeAsync_MappedVector3D_t

  subroutine ApplyFlip_MappedVector3D_t(this,mesh)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedVector3D_t),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: e1,s1,e2,s2,idir
    integer :: i,i2,j,j2
    integer :: r2,flip,ivar
    integer :: bcid
    real(prec) :: extBuff(1:this%interp%N+1,1:this%interp%N+1)

    do idir = 1,3
      do ivar = 1,this%nvar
        do e1 = 1,this%nElem
          do s1 = 1,6

            e2 = mesh%sideInfo(3,s1,e1) ! Neighbor Element
            s2 = mesh%sideInfo(4,s1,e1)/10
            bcid = mesh%sideInfo(5,s1,e1)
            if(e2 > 0) then ! Interior Element
              r2 = mesh%decomp%elemToRank(e2) ! Neighbor Rank

              if(r2 /= mesh%decomp%rankId) then

                flip = mesh%sideInfo(4,s1,e1)-s2*10

                ! Need to update extBoundary with flip applied
                if(flip == 0) then

                  do j = 1,this%interp%N+1
                    do i = 1,this%interp%N+1
                      extBuff(i,j) = this%extBoundary(i,j,s1,e1,ivar,idir)
                    enddo
                  enddo

                else if(flip == 1) then

                  do j = 1,this%interp%N+1
                    do i = 1,this%interp%N+1
                      i2 = this%interp%N+2-i
                      j2 = j
                      extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar,idir)
                    enddo
                  enddo

                else if(flip == 2) then

                  do j = 1,this%interp%N+1
                    do i = 1,this%interp%N+1
                      i2 = this%interp%N+2-i
                      j2 = this%interp%N+2-j
                      extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar,idir)
                    enddo
                  enddo

                else if(flip == 3) then

                  do j = 1,this%interp%N+1
                    do i = 1,this%interp%N+1
                      i2 = i
                      j2 = this%interp%N+2-j
                      extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar,idir)
                    enddo
                  enddo

                else if(flip == 4) then

                  do j = 1,this%interp%N+1
                    do i = 1,this%interp%N+1
                      extBuff(i,j) = this%extBoundary(j,i,s1,e1,ivar,idir)
                    enddo
                  enddo

                else if(flip == 5) then

                  do j = 1,this%interp%N+1
                    do i = 1,this%interp%N+1
                      i2 = this%interp%N+2-j
                      j2 = i
                      extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar,idir)
                    enddo
                  enddo

                else if(flip == 6) then

                  do j = 1,this%interp%N+1
                    do i = 1,this%interp%N+1
                      i2 = this%interp%N+2-j
                      j2 = this%interp%N+2-i
                      extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar,idir)
                    enddo
                  enddo

                else if(flip == 7) then

                  do j = 1,this%interp%N+1
                    do i = 1,this%interp%N+1
                      i2 = j
                      j2 = this%interp%N+2-i
                      extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar,idir)
                    enddo
                  enddo

                endif

                do j = 1,this%interp%N+1
                  do i = 1,this%interp%N+1
                    this%extBoundary(i,j,s1,e1,ivar,idir) = extBuff(i,j)
                  enddo
                enddo

              endif

            endif

          enddo
        enddo
      enddo
    enddo

  endsubroutine ApplyFlip_MappedVector3D_t

  subroutine SideExchange_MappedVector3D_t(this,mesh)
    implicit none
    class(MappedVector3D_t),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip
    integer :: i,i2,j,j2,ivar
    integer :: r2
    integer :: rankId,offset
    integer :: idir
    integer,pointer :: elemtorank(:)

    elemtorank => mesh%decomp%elemToRank(:)
    rankId = mesh%decomp%rankId
    offset = mesh%decomp%offsetElem(rankId+1)

    if(mesh%decomp%mpiEnabled) then
      call this%MPIExchangeAsync(mesh)
    endif

    do concurrent(s1=1:6,e1=1:mesh%nElem,ivar=1:this%nvar,idir=1:3)

      e2Global = mesh%sideInfo(3,s1,e1)
      s2 = mesh%sideInfo(4,s1,e1)/10
      flip = mesh%sideInfo(4,s1,e1)-s2*10

      if(e2Global > 0) then

        r2 = elemToRank(e2Global)

        if(r2 == rankId) then

          e2 = e2Global-offset

          if(flip == 0) then

            do j = 1,this%interp%N+1
              do i = 1,this%interp%N+1
                this%extBoundary(i,j,s1,e1,ivar,idir) = &
                  this%boundary(i,j,s2,e2,ivar,idir)
              enddo
            enddo

          else if(flip == 1) then

            do j = 1,this%interp%N+1
              do i = 1,this%interp%N+1
                i2 = this%interp%N+2-i
                j2 = j
                this%extBoundary(i,j,s1,e1,ivar,idir) = &
                  this%boundary(i2,j2,s2,e2,ivar,idir)
              enddo
            enddo

          else if(flip == 2) then

            do j = 1,this%interp%N+1
              do i = 1,this%interp%N+1
                i2 = this%interp%N+2-i
                j2 = this%interp%N+2-j
                this%extBoundary(i,j,s1,e1,ivar,idir) = &
                  this%boundary(i2,j2,s2,e2,ivar,idir)
              enddo
            enddo

          else if(flip == 3) then

            do j = 1,this%interp%N+1
              do i = 1,this%interp%N+1
                i2 = i
                j2 = this%interp%N+2-j
                this%extBoundary(i,j,s1,e1,ivar,idir) = &
                  this%boundary(i2,j2,s2,e2,ivar,idir)
              enddo
            enddo

          else if(flip == 4) then

            do j = 1,this%interp%N+1
              do i = 1,this%interp%N+1
                this%extBoundary(i,j,s1,e1,ivar,idir) = &
                  this%boundary(j,i,s2,e2,ivar,idir)
              enddo
            enddo

          else if(flip == 5) then

            do j = 1,this%interp%N+1
              do i = 1,this%interp%N+1
                i2 = this%interp%N+2-j
                j2 = i
                this%extBoundary(i,j,s1,e1,ivar,idir) = this%boundary(i2,j2,s2,e2,ivar,idir)
              enddo
            enddo

          else if(flip == 6) then

            do j = 1,this%interp%N+1
              do i = 1,this%interp%N+1
                i2 = this%interp%N+2-j
                j2 = this%interp%N+2-i
                this%extBoundary(i,j,s1,e1,ivar,idir) = this%boundary(i2,j2,s2,e2,ivar,idir)
              enddo
            enddo

          else if(flip == 7) then

            do j = 1,this%interp%N+1
              do i = 1,this%interp%N+1
                i2 = j
                j2 = this%interp%N+2-i
                this%extBoundary(i,j,s1,e1,ivar,idir) = this%boundary(i2,j2,s2,e2,ivar,idir)
              enddo
            enddo

          endif
        endif
      endif

    enddo

    if(mesh%decomp%mpiEnabled) then
      call mesh%decomp%FinalizeMPIExchangeAsync()
      ! Apply side flips for data exchanged with MPI
      call this%ApplyFlip(mesh)
    endif

  endsubroutine SideExchange_MappedVector3D_t

  subroutine MPIMortarExchangeAsync_MappedVector3D_t(this,mesh)
    !! Posts the point-to-point messages required for mortar interfaces whose big and
    !! small elements reside on different ranks; vector analogue of the scalar
    !! MPIMortarExchangeAsync with one message per direction component.
    implicit none
    class(MappedVector3D_t),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: m,q,ivar,idir
    integer :: eB,sB,rB,eS,sS,rS
    integer :: globalSideId,tag
    integer :: offset
    integer :: iError
    integer :: msgCount

    ! A send and a receive per rank-crossing sub-face, per variable and direction.
    call mesh%decomp%ReserveMessages(2*3*this%nvar*mesh%decomp%CountRemoteMortarSides(mesh%mortarInfo,mesh%nMortars,4))

    msgCount = 0
    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)

    do idir = 1,3
      do ivar = 1,this%nvar
        do m = 1,mesh%nMortars

          eB = mesh%mortarInfo(1,m)
          sB = mesh%mortarInfo(2,m)
          rB = mesh%decomp%elemToRank(eB)

          do q = 1,4

            eS = mesh%mortarInfo(2*q+1,m)
            sS = mesh%mortarInfo(2*q+2,m)/10
            rS = mesh%decomp%elemToRank(eS)
            globalSideId = mesh%mortarInfo(10+q,m)
            tag = globalSideId+mesh%nUniqueSides*(ivar-1+this%nvar*(idir-1))

            if(rB == mesh%decomp%rankId .and. rS /= mesh%decomp%rankId) then

              msgCount = msgCount+1
              call MPI_IRECV(this%mortarBuff(:,:,4+q,m,ivar,idir), &
                             (this%interp%N+1)*(this%interp%N+1), &
                             mesh%decomp%mpiPrec, &
                             rS,tag, &
                             mesh%decomp%mpiComm, &
                             mesh%decomp%requests(msgCount),iError)

              msgCount = msgCount+1
              call MPI_ISEND(this%boundary(:,:,sB,eB-offset,ivar,idir), &
                             (this%interp%N+1)*(this%interp%N+1), &
                             mesh%decomp%mpiPrec, &
                             rS,tag, &
                             mesh%decomp%mpiComm, &
                             mesh%decomp%requests(msgCount),iError)

            elseif(rS == mesh%decomp%rankId .and. rB /= mesh%decomp%rankId) then

              msgCount = msgCount+1
              call MPI_IRECV(this%mortarBuff(:,:,q,m,ivar,idir), &
                             (this%interp%N+1)*(this%interp%N+1), &
                             mesh%decomp%mpiPrec, &
                             rB,tag, &
                             mesh%decomp%mpiComm, &
                             mesh%decomp%requests(msgCount),iError)

              msgCount = msgCount+1
              call MPI_ISEND(this%boundary(:,:,sS,eS-offset,ivar,idir), &
                             (this%interp%N+1)*(this%interp%N+1), &
                             mesh%decomp%mpiPrec, &
                             rB,tag, &
                             mesh%decomp%mpiComm, &
                             mesh%decomp%requests(msgCount),iError)

            endif

          enddo
        enddo
      enddo
    enddo

    mesh%decomp%msgCount = msgCount

  endsubroutine MPIMortarExchangeAsync_MappedVector3D_t

  subroutine MortarExchange_MappedVector3D_t(this,mesh)
    !! Fills the extBoundary attribute on all faces participating in a 2:1
    !! nonconforming (mortar) interface; vector analogue of the scalar MortarExchange
    !! (see MappedScalar3D_t for the algorithm description).
    implicit none
    class(MappedVector3D_t),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: m,q,ivar,idir,i,j,i2,j2
    integer :: eB,eS,sS,flip
    integer :: rankId,offset,N
    integer,pointer :: elemtorank(:)
    real(prec) :: extBuff(1:this%interp%N+1,1:this%interp%N+1)

    ! See https://github.com/FluidNumerics/SELF/issues/54 for the reason behind
    ! this pointer alias
    elemtorank => mesh%decomp%elemToRank(:)
    rankId = mesh%decomp%rankId
    offset = mesh%decomp%offsetElem(rankId+1)
    N = this%interp%N

    if(.not. allocated(this%mortarBuff)) then
      allocate(this%mortarBuff(1:N+1,1:N+1,1:8,1:mesh%nMortars,1:this%nvar,1:3))
      this%mortarBuff = 0.0_prec
    endif

    if(mesh%decomp%mpiEnabled) then
      call this%MPIMortarExchangeAsync(mesh)
    endif

    ! Stage rank-local traces in the big face's coordinates
    do concurrent(m=1:mesh%nMortars,ivar=1:this%nvar,idir=1:3)
      block
        integer :: i,j,i2,j2,q
        integer :: eB,sB,eS,sS,flip

        eB = mesh%mortarInfo(1,m)
        if(elemtorank(eB) == rankId) then
          sB = mesh%mortarInfo(2,m)
          do q = 1,4
            do j = 1,N+1
              do i = 1,N+1
                this%mortarBuff(i,j,q,m,ivar,idir) = &
                  this%boundary(i,j,sB,eB-offset,ivar,idir)
              enddo
            enddo
          enddo
        endif

        do q = 1,4
          eS = mesh%mortarInfo(2*q+1,m)
          if(elemtorank(eS) == rankId) then
            sS = mesh%mortarInfo(2*q+2,m)/10
            flip = mesh%mortarInfo(2*q+2,m)-10*sS
            do j = 1,N+1
              do i = 1,N+1
                call MortarFaceMap(i,j,N,flip,i2,j2)
                this%mortarBuff(i,j,4+q,m,ivar,idir) = &
                  this%boundary(i2,j2,sS,eS-offset,ivar,idir)
              enddo
            enddo
          endif
        enddo
      endblock
    enddo

    if(mesh%decomp%mpiEnabled) then
      call mesh%decomp%FinalizeMPIExchangeAsync()

      ! Reorient small-face traces received over MPI into the big face's coordinates
      do idir = 1,3
        do ivar = 1,this%nvar
          do m = 1,mesh%nMortars
            eB = mesh%mortarInfo(1,m)
            if(elemtorank(eB) == rankId) then
              do q = 1,4
                eS = mesh%mortarInfo(2*q+1,m)
                sS = mesh%mortarInfo(2*q+2,m)/10
                flip = mesh%mortarInfo(2*q+2,m)-10*sS
                if(elemtorank(eS) /= rankId .and. flip /= 0) then
                  do j = 1,N+1
                    do i = 1,N+1
                      call MortarFaceMap(i,j,N,flip,i2,j2)
                      extBuff(i,j) = this%mortarBuff(i2,j2,4+q,m,ivar,idir)
                    enddo
                  enddo
                  do j = 1,N+1
                    do i = 1,N+1
                      this%mortarBuff(i,j,4+q,m,ivar,idir) = extBuff(i,j)
                    enddo
                  enddo
                endif
              enddo
            endif
          enddo
        enddo
      enddo
    endif

    ! Compute external states :
    !  small faces get the restricted big-face trace (exact),
    !  the big face gets the L2 projection of the small-face traces
    do concurrent(m=1:mesh%nMortars,ivar=1:this%nvar,idir=1:3)
      block
        integer :: i,j,ii,jj,i2,j2,kx,ky,q
        integer :: eB,sB,eS,sS,flip
        real(prec) :: fm
        real(prec) :: tmp(1:N+1,1:N+1)
        real(prec) :: acc(1:N+1,1:N+1)

        do q = 1,4
          eS = mesh%mortarInfo(2*q+1,m)
          if(elemtorank(eS) == rankId) then
            sS = mesh%mortarInfo(2*q+2,m)/10
            flip = mesh%mortarInfo(2*q+2,m)-10*sS
            kx = mortarQuadKx(q)
            ky = mortarQuadKy(q)
            do jj = 1,N+1
              do i = 1,N+1
                fm = 0.0_prec
                do ii = 1,N+1
                  fm = fm+this%interp%mortarR(ii,i,kx)* &
                       this%mortarBuff(ii,jj,q,m,ivar,idir)
                enddo
                tmp(i,jj) = fm
              enddo
            enddo
            do j = 1,N+1
              do i = 1,N+1
                fm = 0.0_prec
                do jj = 1,N+1
                  fm = fm+this%interp%mortarR(jj,j,ky)*tmp(i,jj)
                enddo
                call MortarFaceMap(i,j,N,flip,i2,j2)
                this%extBoundary(i2,j2,sS,eS-offset,ivar,idir) = fm
              enddo
            enddo
          endif
        enddo

        eB = mesh%mortarInfo(1,m)
        if(elemtorank(eB) == rankId) then
          sB = mesh%mortarInfo(2,m)
          acc = 0.0_prec
          do q = 1,4
            kx = mortarQuadKx(q)
            ky = mortarQuadKy(q)
            do jj = 1,N+1
              do i = 1,N+1
                fm = 0.0_prec
                do ii = 1,N+1
                  fm = fm+this%interp%mortarP(ii,i,kx)* &
                       this%mortarBuff(ii,jj,4+q,m,ivar,idir)
                enddo
                tmp(i,jj) = fm
              enddo
            enddo
            do j = 1,N+1
              do i = 1,N+1
                fm = 0.0_prec
                do jj = 1,N+1
                  fm = fm+this%interp%mortarP(jj,j,ky)*tmp(i,jj)
                enddo
                acc(i,j) = acc(i,j)+fm
              enddo
            enddo
          enddo
          do j = 1,N+1
            do i = 1,N+1
              this%extBoundary(i,j,sB,eB-offset,ivar,idir) = acc(i,j)
            enddo
          enddo
        endif
      endblock
    enddo

  endsubroutine MortarExchange_MappedVector3D_t

  subroutine MPIMortarFluxAsync_MappedVector3D_t(this,mesh)
    !! Posts the one-directional messages for MortarFluxCollect : each remote small
    !! face sends its boundaryNormal trace to the big face's rank.
    implicit none
    class(MappedVector3D_t),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: m,q,ivar
    integer :: eB,rB,eS,sS,rS
    integer :: globalSideId,tag
    integer :: offset
    integer :: iError
    integer :: msgCount

    ! One message per rank-crossing sub-face, per variable.
    call mesh%decomp%ReserveMessages(this%nvar*mesh%decomp%CountRemoteMortarSides(mesh%mortarInfo,mesh%nMortars,4))

    msgCount = 0
    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)

    do ivar = 1,this%nvar
      do m = 1,mesh%nMortars

        eB = mesh%mortarInfo(1,m)
        rB = mesh%decomp%elemToRank(eB)

        do q = 1,4

          eS = mesh%mortarInfo(2*q+1,m)
          sS = mesh%mortarInfo(2*q+2,m)/10
          rS = mesh%decomp%elemToRank(eS)
          globalSideId = mesh%mortarInfo(10+q,m)
          tag = globalSideId+mesh%nUniqueSides*(ivar-1)

          if(rB == mesh%decomp%rankId .and. rS /= mesh%decomp%rankId) then

            msgCount = msgCount+1
            call MPI_IRECV(this%mortarBuff(:,:,4+q,m,ivar,1), &
                           (this%interp%N+1)*(this%interp%N+1), &
                           mesh%decomp%mpiPrec, &
                           rS,tag, &
                           mesh%decomp%mpiComm, &
                           mesh%decomp%requests(msgCount),iError)

          elseif(rS == mesh%decomp%rankId .and. rB /= mesh%decomp%rankId) then

            msgCount = msgCount+1
            call MPI_ISEND(this%boundaryNormal(:,:,sS,eS-offset,ivar), &
                           (this%interp%N+1)*(this%interp%N+1), &
                           mesh%decomp%mpiPrec, &
                           rB,tag, &
                           mesh%decomp%mpiComm, &
                           mesh%decomp%requests(msgCount),iError)

          endif

        enddo
      enddo
    enddo

    mesh%decomp%msgCount = msgCount

  endsubroutine MPIMortarFluxAsync_MappedVector3D_t

  subroutine MortarFluxCollect_MappedVector3D_t(this,mesh)
    !! Replaces the big-face boundaryNormal trace on each mortar interface with the L2
    !! projection of the four small faces' boundaryNormal traces.
    !!
    !! boundaryNormal holds the Riemann-solved surface-flux integrand f* . nHat * nScale
    !! (see BoundaryFlux in the DG models). Because the small faces' nScale is one
    !! quarter of the big face's and the sub-face coordinate Jacobian is 1/4, the
    !! projected big-face integrand is -4 * sum_q (P_kx x P_ky) g_q, where g_q are the
    !! small-face integrands and the sign accounts for the opposing outward normals.
    !! With this choice, the discrete surface integral of the big face equals minus the
    !! sum of the small faces' discrete surface integrals to roundoff, so the mortar
    !! interface is discretely conservative. Must be called after the model's
    !! BoundaryFlux and before the flux divergence is computed.
    implicit none
    class(MappedVector3D_t),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: m,q,ivar,i,j,i2,j2
    integer :: eB,eS,sS,flip
    integer :: rankId,offset,N
    integer,pointer :: elemtorank(:)
    real(prec) :: extBuff(1:this%interp%N+1,1:this%interp%N+1)

    ! See https://github.com/FluidNumerics/SELF/issues/54 for the reason behind
    ! this pointer alias
    elemtorank => mesh%decomp%elemToRank(:)
    rankId = mesh%decomp%rankId
    offset = mesh%decomp%offsetElem(rankId+1)
    N = this%interp%N

    if(.not. allocated(this%mortarBuff)) then
      allocate(this%mortarBuff(1:N+1,1:N+1,1:8,1:mesh%nMortars,1:this%nvar,1:3))
      this%mortarBuff = 0.0_prec
    endif

    if(mesh%decomp%mpiEnabled) then
      call this%MPIMortarFluxAsync(mesh)
    endif

    ! Stage rank-local small-face integrands in the big face's coordinates
    do concurrent(m=1:mesh%nMortars,ivar=1:this%nvar)
      block
        integer :: i,j,i2,j2,q
        integer :: eS,sS,flip

        do q = 1,4
          eS = mesh%mortarInfo(2*q+1,m)
          if(elemtorank(eS) == rankId) then
            sS = mesh%mortarInfo(2*q+2,m)/10
            flip = mesh%mortarInfo(2*q+2,m)-10*sS
            do j = 1,N+1
              do i = 1,N+1
                call MortarFaceMap(i,j,N,flip,i2,j2)
                this%mortarBuff(i,j,4+q,m,ivar,1) = &
                  this%boundaryNormal(i2,j2,sS,eS-offset,ivar)
              enddo
            enddo
          endif
        enddo
      endblock
    enddo

    if(mesh%decomp%mpiEnabled) then
      call mesh%decomp%FinalizeMPIExchangeAsync()

      ! Reorient small-face integrands received over MPI into the big face's coordinates
      do ivar = 1,this%nvar
        do m = 1,mesh%nMortars
          eB = mesh%mortarInfo(1,m)
          if(elemtorank(eB) == rankId) then
            do q = 1,4
              eS = mesh%mortarInfo(2*q+1,m)
              sS = mesh%mortarInfo(2*q+2,m)/10
              flip = mesh%mortarInfo(2*q+2,m)-10*sS
              if(elemtorank(eS) /= rankId .and. flip /= 0) then
                do j = 1,N+1
                  do i = 1,N+1
                    call MortarFaceMap(i,j,N,flip,i2,j2)
                    extBuff(i,j) = this%mortarBuff(i2,j2,4+q,m,ivar,1)
                  enddo
                enddo
                do j = 1,N+1
                  do i = 1,N+1
                    this%mortarBuff(i,j,4+q,m,ivar,1) = extBuff(i,j)
                  enddo
                enddo
              endif
            enddo
          endif
        enddo
      enddo
    endif

    ! Project the small-face integrands onto the big face's trace space. The factor of
    ! four converts the solution-space projection (the tensor-product mortarP carries
    ! the 1/4 sub-face Jacobian) into the integrand-space projection; the sign accounts
    ! for the opposing outward normals.
    do concurrent(m=1:mesh%nMortars,ivar=1:this%nvar)
      block
        integer :: i,j,ii,jj,kx,ky,q
        integer :: eB,sB
        real(prec) :: fm
        real(prec) :: tmp(1:N+1,1:N+1)
        real(prec) :: acc(1:N+1,1:N+1)

        eB = mesh%mortarInfo(1,m)
        if(elemtorank(eB) == rankId) then
          sB = mesh%mortarInfo(2,m)
          acc = 0.0_prec
          do q = 1,4
            kx = mortarQuadKx(q)
            ky = mortarQuadKy(q)
            do jj = 1,N+1
              do i = 1,N+1
                fm = 0.0_prec
                do ii = 1,N+1
                  fm = fm+this%interp%mortarP(ii,i,kx)* &
                       this%mortarBuff(ii,jj,4+q,m,ivar,1)
                enddo
                tmp(i,jj) = fm
              enddo
            enddo
            do j = 1,N+1
              do i = 1,N+1
                fm = 0.0_prec
                do jj = 1,N+1
                  fm = fm+this%interp%mortarP(jj,j,ky)*tmp(i,jj)
                enddo
                acc(i,j) = acc(i,j)+fm
              enddo
            enddo
          enddo
          do j = 1,N+1
            do i = 1,N+1
              this%boundaryNormal(i,j,sB,eB-offset,ivar) = -4.0_prec*acc(i,j)
            enddo
          enddo
        endif
      endblock
    enddo

  endsubroutine MortarFluxCollect_MappedVector3D_t

  subroutine MappedDivergence_MappedVector3D_t(this,df)
    ! Strong Form Operator
    !    !
    implicit none
    class(MappedVector3D_t),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,j,k,ii
    real(prec) :: dfLoc,Fx,Fy,Fz,Fc

    do concurrent(i=1:this%N+1,j=1:this%N+1, &
                  k=1:this%N+1,iel=1:this%nelem,ivar=1:this%nvar)

      dfLoc = 0.0_prec
      do ii = 1,this%N+1
        ! Convert from physical to computational space
        Fx = this%interior(ii,j,k,iEl,iVar,1)
        Fy = this%interior(ii,j,k,iEl,iVar,2)
        Fz = this%interior(ii,j,k,iEl,iVar,3)
        Fc = this%geometry%dsdx%interior(ii,j,k,iEl,1,1,1)*Fx+ &
             this%geometry%dsdx%interior(ii,j,k,iEl,1,2,1)*Fy+ &
             this%geometry%dsdx%interior(ii,j,k,iEl,1,3,1)*Fz
        dfLoc = dfLoc+this%interp%dMatrix(ii,i)*Fc
      enddo
      dF(i,j,k,iel,ivar) = dfLoc

    enddo

    do concurrent(i=1:this%N+1,j=1:this%N+1, &
                  k=1:this%N+1,iel=1:this%nelem,ivar=1:this%nvar)

      dfLoc = 0.0_prec
      do ii = 1,this%N+1
        ! Convert from physical to computational space
        Fx = this%interior(i,ii,k,iEl,iVar,1)
        Fy = this%interior(i,ii,k,iEl,iVar,2)
        Fz = this%interior(i,ii,k,iEl,iVar,3)
        Fc = this%geometry%dsdx%interior(i,ii,k,iEl,1,1,2)*Fx+ &
             this%geometry%dsdx%interior(i,ii,k,iEl,1,2,2)*Fy+ &
             this%geometry%dsdx%interior(i,ii,k,iEl,1,3,2)*Fz
        dfLoc = dfLoc+this%interp%dMatrix(ii,j)*Fc
      enddo
      dF(i,j,k,iel,ivar) = (dF(i,j,k,iel,ivar)+dfLoc)

    enddo

    do concurrent(i=1:this%N+1,j=1:this%N+1, &
                  k=1:this%N+1,iel=1:this%nelem,ivar=1:this%nvar)

      dfLoc = 0.0_prec
      do ii = 1,this%N+1
        ! Convert from physical to computational space
        Fx = this%interior(i,j,ii,iEl,iVar,1)
        Fy = this%interior(i,j,ii,iEl,iVar,2)
        Fz = this%interior(i,j,ii,iEl,iVar,3)
        Fc = this%geometry%dsdx%interior(i,j,ii,iEl,1,1,3)*Fx+ &
             this%geometry%dsdx%interior(i,j,ii,iEl,1,2,3)*Fy+ &
             this%geometry%dsdx%interior(i,j,ii,iEl,1,3,3)*Fz
        dfLoc = dfLoc+this%interp%dMatrix(ii,k)*Fc
      enddo
      dF(i,j,k,iel,ivar) = (dF(i,j,k,iel,ivar)+dfLoc)/this%geometry%J%interior(i,j,k,iEl,1)

    enddo

  endsubroutine MappedDivergence_MappedVector3D_t

  subroutine MappedDGDivergence_MappedVector3D_t(this,df)
      !! Computes the divergence of a 3-D vector using the weak form
      !! On input, the  attribute of the vector
      !! is assigned and the  attribute is set to the physical
      !! directions of the vector. This method will project the vector
      !! onto the contravariant basis vectors.
    implicit none
    class(MappedVector3D_t),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,j,k,ii
    real(prec) :: dfLoc,Fx,Fy,Fz,Fc

    do concurrent(i=1:this%N+1,j=1:this%N+1, &
                  k=1:this%N+1,iel=1:this%nelem,ivar=1:this%nvar)

      dfLoc = 0.0_prec
      do ii = 1,this%N+1
        ! Convert from physical to computational space
        Fx = this%interior(ii,j,k,iEl,iVar,1)
        Fy = this%interior(ii,j,k,iEl,iVar,2)
        Fz = this%interior(ii,j,k,iEl,iVar,3)
        Fc = this%geometry%dsdx%interior(ii,j,k,iEl,1,1,1)*Fx+ &
             this%geometry%dsdx%interior(ii,j,k,iEl,1,2,1)*Fy+ &
             this%geometry%dsdx%interior(ii,j,k,iEl,1,3,1)*Fz
        dfLoc = dfLoc+this%interp%dgMatrix(ii,i)*Fc
      enddo
      dfLoc = dfLoc+ &
              (this%interp%bMatrix(i,2)*this%boundaryNormal(j,k,3,iel,ivar)+ & ! east
               this%interp%bMatrix(i,1)*this%boundaryNormal(j,k,5,iel,ivar))/ & ! west
              this%interp%qweights(i)
      dF(i,j,k,iel,ivar) = dfLoc

    enddo

    do concurrent(i=1:this%N+1,j=1:this%N+1, &
                  k=1:this%N+1,iel=1:this%nelem,ivar=1:this%nvar)

      dfLoc = 0.0_prec
      do ii = 1,this%N+1
        ! Convert from physical to computational space
        Fx = this%interior(i,ii,k,iEl,iVar,1)
        Fy = this%interior(i,ii,k,iEl,iVar,2)
        Fz = this%interior(i,ii,k,iEl,iVar,3)
        Fc = this%geometry%dsdx%interior(i,ii,k,iEl,1,1,2)*Fx+ &
             this%geometry%dsdx%interior(i,ii,k,iEl,1,2,2)*Fy+ &
             this%geometry%dsdx%interior(i,ii,k,iEl,1,3,2)*Fz
        dfLoc = dfLoc+this%interp%dgMatrix(ii,j)*Fc
      enddo
      dfLoc = +dfLoc+ &
              (this%interp%bMatrix(j,2)*this%boundaryNormal(i,k,4,iel,ivar)+ & ! north
               this%interp%bMatrix(j,1)*this%boundaryNormal(i,k,2,iel,ivar))/ & ! south
              this%interp%qweights(j)
      dF(i,j,k,iel,ivar) = (dF(i,j,k,iel,ivar)+dfLoc)

    enddo

    do concurrent(i=1:this%N+1,j=1:this%N+1, &
                  k=1:this%N+1,iel=1:this%nelem,ivar=1:this%nvar)

      dfLoc = 0.0_prec
      do ii = 1,this%N+1
        ! Convert from physical to computational space
        Fx = this%interior(i,j,ii,iEl,iVar,1)
        Fy = this%interior(i,j,ii,iEl,iVar,2)
        Fz = this%interior(i,j,ii,iEl,iVar,3)
        Fc = this%geometry%dsdx%interior(i,j,ii,iEl,1,1,3)*Fx+ &
             this%geometry%dsdx%interior(i,j,ii,iEl,1,2,3)*Fy+ &
             this%geometry%dsdx%interior(i,j,ii,iEl,1,3,3)*Fz
        dfLoc = dfLoc+this%interp%dgMatrix(ii,k)*Fc
      enddo
      dfLoc = dfLoc+ &
              (this%interp%bMatrix(k,2)*this%boundaryNormal(i,j,6,iel,ivar)+ & ! top
               this%interp%bMatrix(k,1)*this%boundaryNormal(i,j,1,iel,ivar))/ & ! bottom
              this%interp%qweights(k)
      dF(i,j,k,iel,ivar) = (dF(i,j,k,iel,ivar)+dfLoc)/this%geometry%J%interior(i,j,k,iEl,1)

    enddo

  endsubroutine MappedDGDivergence_MappedVector3D_t

  ! subroutine WriteTecplot_MappedVector3D_t(this,geometry,filename)
  !   implicit none
  !   class(MappedVector3D_t),intent(inout) :: this
  !   type(SEMHex),intent(in) :: geometry
  !   character(*),intent(in),optional :: filename
  !   ! Local
  !   character(8) :: zoneID
  !   integer :: fUnit
  !   integer :: iEl,i,j,k,iVar
  !   character(LEN=self_FileNameLength) :: tecFile
  !   character(LEN=self_TecplotHeaderLength) :: tecHeader
  !   character(LEN=self_FormatLength) :: fmat
  !   character(13) :: timeStampString
  !   character(5) :: rankString
  !   real(prec) :: f(1:this%M+1,1:this%M+1,1:this%M+1,1:this%nelem,1:this%nvar,1:3)
  !   real(prec) :: x(1:this%M+1,1:this%M+1,1:this%M+1,1:this%nelem,1:this%nvar,1:3)

  !   if(present(filename)) then
  !     tecFile = filename
  !   else
  !     tecFile = "mappedvector.tec"
  !   endif

  !   ! Map the mesh positions to the target grid
  !   x = geometry%x%GridInterp()

  !   ! Map the solution to the target grid
  !   f = this%GridInterp()

  !   open(UNIT=NEWUNIT(fUnit), &
  !        FILE=trim(tecFile), &
  !        FORM='formatted', &
  !        STATUS='replace')

  !   tecHeader = 'VARIABLES = "X", "Y", "Z"'
  !   do iVar = 1,this%nVar
  !     tecHeader = trim(tecHeader)//', "'//trim(this%meta(iVar)%name)//'_x"'
  !   enddo
  !   do iVar = 1,this%nVar
  !     tecHeader = trim(tecHeader)//', "'//trim(this%meta(iVar)%name)//'_y"'
  !   enddo
  !   do iVar = 1,this%nVar
  !     tecHeader = trim(tecHeader)//', "'//trim(this%meta(iVar)%name)//'_z"'
  !   enddo

  !   write(fUnit,*) trim(tecHeader)

  !   ! Create format statement
  !   write(fmat,*) 3*this%nvar+3
  !   fmat = '('//trim(fmat)//'(ES16.7E3,1x))'

  !   do iEl = 1,this%nElem

  !     write(zoneID,'(I8.8)') iEl
  !     write(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',this%interp%M+1, &
  !       ', J=',this%interp%M+1,', K=',this%interp%M+1

  !     do k = 1,this%interp%M+1
  !       do j = 1,this%interp%M+1
  !         do i = 1,this%interp%M+1

  !           write(fUnit,fmat) x(i,j,k,iEl,1,1), &
  !             x(i,j,k,iEl,1,2), &
  !             x(i,j,k,iEl,1,3), &
  !             f(i,j,k,iEl,1:this%nvar,1), &
  !             f(i,j,k,iEl,1:this%nvar,2), &
  !             f(i,j,k,iEl,1:this%nvar,3)

  !         enddo
  !       enddo
  !     enddo

  !   enddo

  !   close(UNIT=fUnit)

  ! endsubroutine WriteTecplot_MappedVector3D_t

endmodule SELF_MappedVector_3D_t
