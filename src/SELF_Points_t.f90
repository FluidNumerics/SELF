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

module SELF_Points_t
!! A point-cloud container with type-bound procedures to (1) locate each
!! physical point inside a SEMQuad/SEMHex element-of-elements and recover its
!! reference (computational) coordinates via a Newton inverse-map, and (2)
!! sample MappedScalar2D / MappedScalar3D fields at the located points.
!!
!! The point-location uses a uniform-grid spatial hash built from per-element
!! axis-aligned bounding boxes of the discrete control-grid (geometry%x). For
!! each candidate element, a Newton iteration on
!!
!!   X(xi) = sum_{i,j[,k]} x_node(i,j[,k]) * l_i(s) * l_j(t) [* l_k(u)]
!!
!! is driven to convergence using the high-order Jacobian dX/dxi (interpolated
!! covariant-basis tensor dxds). Points that fall outside every element this
!! rank owns are marked with the sentinel elements(p) = 0 and their evaluated
!! field values are returned as zero.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Geometry_2D
  use SELF_Geometry_3D
  use SELF_MappedScalar_2D
  use SELF_MappedScalar_3D

  implicit none

  type,public :: Points_t
    integer :: nPoints = 0
    integer :: nDim = 0
    integer :: nCached = 0
      !! Polynomial degree at which the per-point Lagrange basis cache is
      !! valid. Zero means no cache. The cache is filled by LocatePoints and
      !! consumed by EvaluateScalar when the field's interpolant degree
      !! matches.
    real(prec),allocatable :: x(:,:)
      !! Physical coordinates, shape (1:nPoints, 1:nDim). User input.
    integer,pointer :: elements(:) => null()
      !! Element id (rank-local) containing each point; 0 = not found.
    real(prec),pointer :: coordinates(:,:) => null()
      !! Reference coordinates (s,t[,u]) in [-1,1]^nDim, shape (1:nPoints, 1:nDim).
      !! Defined only where elements(p) > 0.
    real(prec),pointer :: lS_cache(:,:) => null()
      !! Lagrange basis at coordinates(p,1), shape (0:nCached, 1:nPoints).
      !! Filled by LocatePoints for points with elements(p) > 0.
    real(prec),pointer :: lT_cache(:,:) => null()
      !! Lagrange basis at coordinates(p,2), shape (0:nCached, 1:nPoints).
    real(prec),pointer :: lU_cache(:,:) => null()
      !! Lagrange basis at coordinates(p,3), shape (0:nCached, 1:nPoints).
      !! Allocated only for nDim = 3.

  contains

    procedure,public :: Init => Init_Points_t
    procedure,public :: Free => Free_Points_t
    procedure,public :: SetPoints => SetPoints_Points_t

    generic,public :: LocatePoints => LocatePoints_2D_Points_t,LocatePoints_3D_Points_t
    procedure,public :: LocatePoints_2D_Points_t
    procedure,public :: LocatePoints_3D_Points_t

    generic,public :: EvaluateScalar => EvalScalar_2D_Points_t,EvalScalar_3D_Points_t
    procedure,public :: EvalScalar_2D_Points_t
    procedure,public :: EvalScalar_3D_Points_t

    generic,public :: DiracDelta => DiracDelta_2D_Points_t,DiracDelta_3D_Points_t
    procedure,public :: DiracDelta_2D_Points_t
    procedure,public :: DiracDelta_3D_Points_t

  endtype Points_t

contains

  subroutine Init_Points_t(this,nPoints,nDim)
    !! Allocate storage for a point cloud of size nPoints in nDim dimensions.
    !! nDim must be 2 or 3.
    implicit none
    class(Points_t),intent(out) :: this
    integer,intent(in) :: nPoints
    integer,intent(in) :: nDim

    if(nDim /= 2 .and. nDim /= 3) then
      print*,"SELF_Points_t::Init: nDim must be 2 or 3, got ",nDim
      stop 1
    endif

    this%nPoints = nPoints
    this%nDim = nDim
    allocate(this%x(1:nPoints,1:nDim))
    allocate(this%elements(1:nPoints))
    allocate(this%coordinates(1:nPoints,1:nDim))
    this%x = 0.0_prec
    this%elements = 0
    this%coordinates = 0.0_prec

  endsubroutine Init_Points_t

  subroutine Free_Points_t(this)
    implicit none
    class(Points_t),intent(inout) :: this

    if(allocated(this%x)) deallocate(this%x)
    if(associated(this%elements)) deallocate(this%elements)
    if(associated(this%coordinates)) deallocate(this%coordinates)
    if(associated(this%lS_cache)) deallocate(this%lS_cache)
    if(associated(this%lT_cache)) deallocate(this%lT_cache)
    if(associated(this%lU_cache)) deallocate(this%lU_cache)
    this%elements => null()
    this%coordinates => null()
    this%lS_cache => null()
    this%lT_cache => null()
    this%lU_cache => null()
    this%nPoints = 0
    this%nDim = 0
    this%nCached = 0

  endsubroutine Free_Points_t

  subroutine SetPoints_Points_t(this,xIn)
    !! Copy user-supplied physical coordinates into the cloud.
    !! xIn must be shape (nPoints, nDim) matching the init dimensions.
    implicit none
    class(Points_t),intent(inout) :: this
    real(prec),intent(in) :: xIn(:,:)
    ! Local
    integer :: p,d

    if(size(xIn,1) /= this%nPoints .or. size(xIn,2) /= this%nDim) then
      print*,"SELF_Points_t::SetPoints: shape mismatch"
      stop 1
    endif

    do p = 1,this%nPoints
      do d = 1,this%nDim
        this%x(p,d) = xIn(p,d)
      enddo
    enddo
    this%elements = 0
    this%nCached = 0

  endsubroutine SetPoints_Points_t

  subroutine LocatePoints_2D_Points_t(this,geometry)
    !! Locate each stored physical point inside the 2D SEMQuad geometry. On
    !! exit, elements(p) holds the (rank-local) element id and coordinates(p,:)
    !! holds the (s,t) reference coordinates, both for points that resolve.
    !! Points that resolve to no element are marked with elements(p) = 0.
    implicit none
    class(Points_t),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    ! Local
    integer :: p,iEl,iCand,nCand
    integer :: nElem,N
    integer :: ix,iy,nCellsX,nCellsY,nCellsTotal,cellId
    real(prec) :: gMin(2),gMax(2),cellSize(2),slack
    real(prec) :: bbDiag,padX,padY
    real(prec) :: xTarget(2),xi(2)
    real(prec),allocatable :: bbMin(:,:),bbMax(:,:)
    integer,allocatable :: cellStart(:),cellElems(:)
    integer :: totalEntries,ixLo,ixHi,iyLo,iyHi,k
    logical :: converged

    if(this%nDim /= 2) then
      print*,"SELF_Points_t::LocatePoints (2D): nDim must be 2"
      stop 1
    endif

    nElem = geometry%nElem
    N = geometry%x%interp%N

    this%elements = 0
    this%coordinates = 0.0_prec

    if(nElem == 0 .or. this%nPoints == 0) return

    ! --- Step 1: per-element axis-aligned bounding boxes --------------------
    allocate(bbMin(1:2,1:nElem),bbMax(1:2,1:nElem))
    call BuildElementBBoxes_2D(geometry,N,nElem,bbMin,bbMax)

    ! --- Step 2: global bbox + spatial-hash grid sizing ---------------------
    gMin(1) = minval(bbMin(1,:))
    gMin(2) = minval(bbMin(2,:))
    gMax(1) = maxval(bbMax(1,:))
    gMax(2) = maxval(bbMax(2,:))
    bbDiag = sqrt((gMax(1)-gMin(1))**2+(gMax(2)-gMin(2))**2)
    slack = max(bbDiag*1.0e-8_prec,tiny(1.0_prec)*1.0e6_prec)
    gMin = gMin-slack
    gMax = gMax+slack

    ! Aim for ~1 element per cell along each axis. Inflate by 25% to keep
    ! buckets short. Floor at 1.
    nCellsX = max(1,ceiling(sqrt(real(nElem,prec)*1.25_prec)))
    nCellsY = nCellsX
    nCellsTotal = nCellsX*nCellsY
    cellSize(1) = (gMax(1)-gMin(1))/real(nCellsX,prec)
    cellSize(2) = (gMax(2)-gMin(2))/real(nCellsY,prec)
    if(cellSize(1) <= 0.0_prec) cellSize(1) = 1.0_prec
    if(cellSize(2) <= 0.0_prec) cellSize(2) = 1.0_prec

    ! --- Step 3: build CSR-style cell->elements hash ------------------------
    allocate(cellStart(0:nCellsTotal))
    cellStart = 0

    ! Pass 1: count entries per cell.
    do iEl = 1,nElem
      padX = max((bbMax(1,iEl)-bbMin(1,iEl))*1.0e-8_prec,slack)
      padY = max((bbMax(2,iEl)-bbMin(2,iEl))*1.0e-8_prec,slack)
      ixLo = ClampCell((bbMin(1,iEl)-padX-gMin(1))/cellSize(1),nCellsX)
      ixHi = ClampCell((bbMax(1,iEl)+padX-gMin(1))/cellSize(1),nCellsX)
      iyLo = ClampCell((bbMin(2,iEl)-padY-gMin(2))/cellSize(2),nCellsY)
      iyHi = ClampCell((bbMax(2,iEl)+padY-gMin(2))/cellSize(2),nCellsY)
      do iy = iyLo,iyHi
        do ix = ixLo,ixHi
          cellId = ix+iy*nCellsX
          cellStart(cellId+1) = cellStart(cellId+1)+1
        enddo
      enddo
    enddo

    ! Prefix sum -> cellStart now indexes into cellElems.
    do k = 1,nCellsTotal
      cellStart(k) = cellStart(k)+cellStart(k-1)
    enddo
    totalEntries = cellStart(nCellsTotal)
    allocate(cellElems(1:totalEntries))

    ! Pass 2: place entries.
    do iEl = 1,nElem
      padX = max((bbMax(1,iEl)-bbMin(1,iEl))*1.0e-8_prec,slack)
      padY = max((bbMax(2,iEl)-bbMin(2,iEl))*1.0e-8_prec,slack)
      ixLo = ClampCell((bbMin(1,iEl)-padX-gMin(1))/cellSize(1),nCellsX)
      ixHi = ClampCell((bbMax(1,iEl)+padX-gMin(1))/cellSize(1),nCellsX)
      iyLo = ClampCell((bbMin(2,iEl)-padY-gMin(2))/cellSize(2),nCellsY)
      iyHi = ClampCell((bbMax(2,iEl)+padY-gMin(2))/cellSize(2),nCellsY)
      do iy = iyLo,iyHi
        do ix = ixLo,ixHi
          cellId = ix+iy*nCellsX
          cellStart(cellId) = cellStart(cellId)+1
          cellElems(cellStart(cellId)) = iEl
        enddo
      enddo
    enddo

    ! Restore cellStart to start-offsets (currently holds end-offsets).
    do k = nCellsTotal,1,-1
      cellStart(k) = cellStart(k-1)
    enddo
    cellStart(0) = 0

    ! --- Step 4: per-point location -----------------------------------------
    do p = 1,this%nPoints
      xTarget(1) = this%x(p,1)
      xTarget(2) = this%x(p,2)

      ! Quick reject: outside global bbox.
      if(xTarget(1) < gMin(1) .or. xTarget(1) > gMax(1) .or. &
         xTarget(2) < gMin(2) .or. xTarget(2) > gMax(2)) cycle

      ix = ClampCell((xTarget(1)-gMin(1))/cellSize(1),nCellsX)
      iy = ClampCell((xTarget(2)-gMin(2))/cellSize(2),nCellsY)
      cellId = ix+iy*nCellsX
      nCand = cellStart(cellId+1)-cellStart(cellId)

      do iCand = 1,nCand
        iEl = cellElems(cellStart(cellId)+iCand)

        ! Cheap AABB reject (with slack).
        padX = max((bbMax(1,iEl)-bbMin(1,iEl))*1.0e-8_prec,slack)
        padY = max((bbMax(2,iEl)-bbMin(2,iEl))*1.0e-8_prec,slack)
        if(xTarget(1) < bbMin(1,iEl)-padX .or. xTarget(1) > bbMax(1,iEl)+padX) cycle
        if(xTarget(2) < bbMin(2,iEl)-padY .or. xTarget(2) > bbMax(2,iEl)+padY) cycle

        call NewtonInverse_2D(geometry,iEl,N,xTarget,xi,converged)
        if(.not. converged) cycle

        ! Accept if reference coords lie in the (slightly inflated) bi-unit cube.
        if(abs(xi(1)) <= 1.0_prec+1.0e-6_prec .and. &
           abs(xi(2)) <= 1.0_prec+1.0e-6_prec) then
          this%elements(p) = iEl
          this%coordinates(p,1) = min(1.0_prec,max(-1.0_prec,xi(1)))
          this%coordinates(p,2) = min(1.0_prec,max(-1.0_prec,xi(2)))
          exit
        endif
      enddo
    enddo

    ! Fill per-point Lagrange basis cache: the basis values are a function of
    ! the reference coordinates only, so callers that sample many times reuse
    ! these without recomputing CalculateLagrangePolynomials per call.
    if(associated(this%lS_cache)) deallocate(this%lS_cache)
    if(associated(this%lT_cache)) deallocate(this%lT_cache)
    if(associated(this%lU_cache)) deallocate(this%lU_cache)
    this%lS_cache => null()
    this%lT_cache => null()
    this%lU_cache => null()
    allocate(this%lS_cache(0:N,1:this%nPoints))
    allocate(this%lT_cache(0:N,1:this%nPoints))
    this%lS_cache = 0.0_prec
    this%lT_cache = 0.0_prec
    do p = 1,this%nPoints
      if(this%elements(p) <= 0) cycle
      this%lS_cache(:,p) = geometry%x%interp%CalculateLagrangePolynomials(this%coordinates(p,1))
      this%lT_cache(:,p) = geometry%x%interp%CalculateLagrangePolynomials(this%coordinates(p,2))
    enddo
    this%nCached = N

    deallocate(bbMin,bbMax,cellStart,cellElems)

  endsubroutine LocatePoints_2D_Points_t

  subroutine LocatePoints_3D_Points_t(this,geometry)
    !! Locate each stored physical point inside the 3D SEMHex geometry. On
    !! exit, elements(p) and coordinates(p,1:3) hold the located element id
    !! and reference (s,t,u) for points that resolve; elements(p) = 0
    !! otherwise.
    implicit none
    class(Points_t),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    ! Local
    integer :: p,iEl,iCand,nCand
    integer :: nElem,N
    integer :: ix,iy,iz,nCellsX,nCellsY,nCellsZ,nCellsTotal,cellId
    real(prec) :: gMin(3),gMax(3),cellSize(3),slack,bbDiag
    real(prec) :: padX,padY,padZ
    real(prec) :: xTarget(3),xi(3)
    real(prec),allocatable :: bbMin(:,:),bbMax(:,:)
    integer,allocatable :: cellStart(:),cellElems(:)
    integer :: totalEntries,ixLo,ixHi,iyLo,iyHi,izLo,izHi,k
    real(prec) :: nC
    logical :: converged

    if(this%nDim /= 3) then
      print*,"SELF_Points_t::LocatePoints (3D): nDim must be 3"
      stop 1
    endif

    nElem = geometry%nElem
    N = geometry%x%interp%N

    this%elements = 0
    this%coordinates = 0.0_prec

    if(nElem == 0 .or. this%nPoints == 0) return

    ! --- Step 1: per-element AABBs ------------------------------------------
    allocate(bbMin(1:3,1:nElem),bbMax(1:3,1:nElem))
    call BuildElementBBoxes_3D(geometry,N,nElem,bbMin,bbMax)

    ! --- Step 2: global bbox + grid sizing ----------------------------------
    gMin(1) = minval(bbMin(1,:))
    gMin(2) = minval(bbMin(2,:))
    gMin(3) = minval(bbMin(3,:))
    gMax(1) = maxval(bbMax(1,:))
    gMax(2) = maxval(bbMax(2,:))
    gMax(3) = maxval(bbMax(3,:))
    bbDiag = sqrt((gMax(1)-gMin(1))**2+(gMax(2)-gMin(2))**2+(gMax(3)-gMin(3))**2)
    slack = max(bbDiag*1.0e-8_prec,tiny(1.0_prec)*1.0e6_prec)
    gMin = gMin-slack
    gMax = gMax+slack

    nC = real(nElem,prec)*1.25_prec
    nCellsX = max(1,ceiling(nC**(1.0_prec/3.0_prec)))
    nCellsY = nCellsX
    nCellsZ = nCellsX
    nCellsTotal = nCellsX*nCellsY*nCellsZ
    cellSize(1) = (gMax(1)-gMin(1))/real(nCellsX,prec)
    cellSize(2) = (gMax(2)-gMin(2))/real(nCellsY,prec)
    cellSize(3) = (gMax(3)-gMin(3))/real(nCellsZ,prec)
    if(cellSize(1) <= 0.0_prec) cellSize(1) = 1.0_prec
    if(cellSize(2) <= 0.0_prec) cellSize(2) = 1.0_prec
    if(cellSize(3) <= 0.0_prec) cellSize(3) = 1.0_prec

    ! --- Step 3: CSR-style cell->elements hash ------------------------------
    allocate(cellStart(0:nCellsTotal))
    cellStart = 0

    do iEl = 1,nElem
      padX = max((bbMax(1,iEl)-bbMin(1,iEl))*1.0e-8_prec,slack)
      padY = max((bbMax(2,iEl)-bbMin(2,iEl))*1.0e-8_prec,slack)
      padZ = max((bbMax(3,iEl)-bbMin(3,iEl))*1.0e-8_prec,slack)
      ixLo = ClampCell((bbMin(1,iEl)-padX-gMin(1))/cellSize(1),nCellsX)
      ixHi = ClampCell((bbMax(1,iEl)+padX-gMin(1))/cellSize(1),nCellsX)
      iyLo = ClampCell((bbMin(2,iEl)-padY-gMin(2))/cellSize(2),nCellsY)
      iyHi = ClampCell((bbMax(2,iEl)+padY-gMin(2))/cellSize(2),nCellsY)
      izLo = ClampCell((bbMin(3,iEl)-padZ-gMin(3))/cellSize(3),nCellsZ)
      izHi = ClampCell((bbMax(3,iEl)+padZ-gMin(3))/cellSize(3),nCellsZ)
      do iz = izLo,izHi
        do iy = iyLo,iyHi
          do ix = ixLo,ixHi
            cellId = ix+nCellsX*(iy+nCellsY*iz)
            cellStart(cellId+1) = cellStart(cellId+1)+1
          enddo
        enddo
      enddo
    enddo

    do k = 1,nCellsTotal
      cellStart(k) = cellStart(k)+cellStart(k-1)
    enddo
    totalEntries = cellStart(nCellsTotal)
    allocate(cellElems(1:totalEntries))

    do iEl = 1,nElem
      padX = max((bbMax(1,iEl)-bbMin(1,iEl))*1.0e-8_prec,slack)
      padY = max((bbMax(2,iEl)-bbMin(2,iEl))*1.0e-8_prec,slack)
      padZ = max((bbMax(3,iEl)-bbMin(3,iEl))*1.0e-8_prec,slack)
      ixLo = ClampCell((bbMin(1,iEl)-padX-gMin(1))/cellSize(1),nCellsX)
      ixHi = ClampCell((bbMax(1,iEl)+padX-gMin(1))/cellSize(1),nCellsX)
      iyLo = ClampCell((bbMin(2,iEl)-padY-gMin(2))/cellSize(2),nCellsY)
      iyHi = ClampCell((bbMax(2,iEl)+padY-gMin(2))/cellSize(2),nCellsY)
      izLo = ClampCell((bbMin(3,iEl)-padZ-gMin(3))/cellSize(3),nCellsZ)
      izHi = ClampCell((bbMax(3,iEl)+padZ-gMin(3))/cellSize(3),nCellsZ)
      do iz = izLo,izHi
        do iy = iyLo,iyHi
          do ix = ixLo,ixHi
            cellId = ix+nCellsX*(iy+nCellsY*iz)
            cellStart(cellId) = cellStart(cellId)+1
            cellElems(cellStart(cellId)) = iEl
          enddo
        enddo
      enddo
    enddo

    do k = nCellsTotal,1,-1
      cellStart(k) = cellStart(k-1)
    enddo
    cellStart(0) = 0

    ! --- Step 4: per-point location -----------------------------------------
    do p = 1,this%nPoints
      xTarget(1) = this%x(p,1)
      xTarget(2) = this%x(p,2)
      xTarget(3) = this%x(p,3)

      if(xTarget(1) < gMin(1) .or. xTarget(1) > gMax(1) .or. &
         xTarget(2) < gMin(2) .or. xTarget(2) > gMax(2) .or. &
         xTarget(3) < gMin(3) .or. xTarget(3) > gMax(3)) cycle

      ix = ClampCell((xTarget(1)-gMin(1))/cellSize(1),nCellsX)
      iy = ClampCell((xTarget(2)-gMin(2))/cellSize(2),nCellsY)
      iz = ClampCell((xTarget(3)-gMin(3))/cellSize(3),nCellsZ)
      cellId = ix+nCellsX*(iy+nCellsY*iz)
      nCand = cellStart(cellId+1)-cellStart(cellId)

      do iCand = 1,nCand
        iEl = cellElems(cellStart(cellId)+iCand)

        padX = max((bbMax(1,iEl)-bbMin(1,iEl))*1.0e-8_prec,slack)
        padY = max((bbMax(2,iEl)-bbMin(2,iEl))*1.0e-8_prec,slack)
        padZ = max((bbMax(3,iEl)-bbMin(3,iEl))*1.0e-8_prec,slack)
        if(xTarget(1) < bbMin(1,iEl)-padX .or. xTarget(1) > bbMax(1,iEl)+padX) cycle
        if(xTarget(2) < bbMin(2,iEl)-padY .or. xTarget(2) > bbMax(2,iEl)+padY) cycle
        if(xTarget(3) < bbMin(3,iEl)-padZ .or. xTarget(3) > bbMax(3,iEl)+padZ) cycle

        call NewtonInverse_3D(geometry,iEl,N,xTarget,xi,converged)
        if(.not. converged) cycle

        if(abs(xi(1)) <= 1.0_prec+1.0e-6_prec .and. &
           abs(xi(2)) <= 1.0_prec+1.0e-6_prec .and. &
           abs(xi(3)) <= 1.0_prec+1.0e-6_prec) then
          this%elements(p) = iEl
          this%coordinates(p,1) = min(1.0_prec,max(-1.0_prec,xi(1)))
          this%coordinates(p,2) = min(1.0_prec,max(-1.0_prec,xi(2)))
          this%coordinates(p,3) = min(1.0_prec,max(-1.0_prec,xi(3)))
          exit
        endif
      enddo
    enddo

    ! Fill per-point basis cache (see 2D variant for rationale).
    if(associated(this%lS_cache)) deallocate(this%lS_cache)
    if(associated(this%lT_cache)) deallocate(this%lT_cache)
    if(associated(this%lU_cache)) deallocate(this%lU_cache)
    this%lS_cache => null()
    this%lT_cache => null()
    this%lU_cache => null()
    allocate(this%lS_cache(0:N,1:this%nPoints))
    allocate(this%lT_cache(0:N,1:this%nPoints))
    allocate(this%lU_cache(0:N,1:this%nPoints))
    this%lS_cache = 0.0_prec
    this%lT_cache = 0.0_prec
    this%lU_cache = 0.0_prec
    do p = 1,this%nPoints
      if(this%elements(p) <= 0) cycle
      this%lS_cache(:,p) = geometry%x%interp%CalculateLagrangePolynomials(this%coordinates(p,1))
      this%lT_cache(:,p) = geometry%x%interp%CalculateLagrangePolynomials(this%coordinates(p,2))
      this%lU_cache(:,p) = geometry%x%interp%CalculateLagrangePolynomials(this%coordinates(p,3))
    enddo
    this%nCached = N

    deallocate(bbMin,bbMax,cellStart,cellElems)

  endsubroutine LocatePoints_3D_Points_t

  subroutine EvalScalar_2D_Points_t(this,scalar,values)
    !! Evaluate a 2D MappedScalar at all located points by tensor-product
    !! Lagrange interpolation at the stored reference coordinates. Points with
    !! elements(p) == 0 receive a value of zero. When LocatePoints has cached
    !! the per-point basis at the matching polynomial degree, this routine
    !! reuses it and skips Lagrange-polynomial evaluation altogether.
    implicit none
    class(Points_t),intent(in) :: this
    class(MappedScalar2D_t),intent(in) :: scalar
    real(prec),intent(out) :: values(1:this%nPoints,1:scalar%nVar)
    ! Local
    integer :: p,iEl,iVar,i,j,N
    real(prec) :: fij,fi
    real(prec),allocatable :: lS(:),lT(:)
    logical :: useCache

    if(this%nDim /= 2) then
      print*,"SELF_Points_t::EvaluateScalar (2D): nDim must be 2"
      stop 1
    endif

    N = scalar%interp%N
    useCache = (this%nCached == N) .and. associated(this%lS_cache) .and. associated(this%lT_cache)

    allocate(lS(0:N),lT(0:N))
    values = 0.0_prec

    do p = 1,this%nPoints
      iEl = this%elements(p)
      if(iEl <= 0) cycle

      if(useCache) then
        lS = this%lS_cache(:,p)
        lT = this%lT_cache(:,p)
      else
        lS = scalar%interp%CalculateLagrangePolynomials(this%coordinates(p,1))
        lT = scalar%interp%CalculateLagrangePolynomials(this%coordinates(p,2))
      endif

      do iVar = 1,scalar%nVar
        fij = 0.0_prec
        do j = 1,N+1
          fi = 0.0_prec
          do i = 1,N+1
            fi = fi+lS(i-1)*scalar%interior(i,j,iEl,iVar)
          enddo
          fij = fij+lT(j-1)*fi
        enddo
        values(p,iVar) = fij
      enddo
    enddo

    deallocate(lS,lT)

  endsubroutine EvalScalar_2D_Points_t

  subroutine EvalScalar_3D_Points_t(this,scalar,values)
    !! Evaluate a 3D MappedScalar at all located points by tensor-product
    !! Lagrange interpolation. Points with elements(p) == 0 receive zero.
    !! When LocatePoints has cached the basis at the matching degree, this
    !! routine reuses it.
    implicit none
    class(Points_t),intent(in) :: this
    class(MappedScalar3D_t),intent(in) :: scalar
    real(prec),intent(out) :: values(1:this%nPoints,1:scalar%nVar)
    ! Local
    integer :: p,iEl,iVar,i,j,k,N
    real(prec) :: fijk,fij,fi
    real(prec),allocatable :: lS(:),lT(:),lU(:)
    logical :: useCache

    if(this%nDim /= 3) then
      print*,"SELF_Points_t::EvaluateScalar (3D): nDim must be 3"
      stop 1
    endif

    N = scalar%interp%N
    useCache = (this%nCached == N) .and. associated(this%lS_cache) .and. &
               associated(this%lT_cache) .and. associated(this%lU_cache)

    allocate(lS(0:N),lT(0:N),lU(0:N))
    values = 0.0_prec

    do p = 1,this%nPoints
      iEl = this%elements(p)
      if(iEl <= 0) cycle

      if(useCache) then
        lS = this%lS_cache(:,p)
        lT = this%lT_cache(:,p)
        lU = this%lU_cache(:,p)
      else
        lS = scalar%interp%CalculateLagrangePolynomials(this%coordinates(p,1))
        lT = scalar%interp%CalculateLagrangePolynomials(this%coordinates(p,2))
        lU = scalar%interp%CalculateLagrangePolynomials(this%coordinates(p,3))
      endif

      do iVar = 1,scalar%nVar
        fijk = 0.0_prec
        do k = 1,N+1
          fij = 0.0_prec
          do j = 1,N+1
            fi = 0.0_prec
            do i = 1,N+1
              fi = fi+lS(i-1)*scalar%interior(i,j,k,iEl,iVar)
            enddo
            fij = fij+lT(j-1)*fi
          enddo
          fijk = fijk+lU(k-1)*fij
        enddo
        values(p,iVar) = fijk
      enddo
    enddo

    deallocate(lS,lT,lU)

  endsubroutine EvalScalar_3D_Points_t

  subroutine DiracDelta_2D_Points_t(this,geometry,scalar)
    !! Scatter a discrete Dirac delta of unit strength (S = 1) onto a 2D
    !! MappedScalar, one variable per stored point. Variable p of scalar
    !! receives the delta associated with point p.
    !!
    !! After the call, for the element iEl = elements(p) > 0:
    !!
    !!   scalar%interior(i,j,iEl,p) = l_i(xi_p) * l_j(eta_p)
    !!                              / ( w_i * w_j * J_0(p) )
    !!
    !! and scalar%interior(i,j,iEl',p) = 0 for all iEl' /= iEl. The J_0(p)
    !! denominator is the polynomial interpolation of geometry%J at the
    !! source's reference coordinates, J_0 = sum_{m,n} l_m(xi_p) l_n(eta_p)
    !! * J_{m,n}^{iEl}. For points with elements(p) == 0 (not located on this
    !! rank), variable p is filled with zeros.
    !!
    !! This is the post-mass-matrix form: it can be added directly to dSdt /
    !! the source term in a DGSEM model without further mass-matrix inversion,
    !! since SELF's MappedDGDivergence already divides by w_i w_j J_{ij}.
    !!
    !! Conservation property (exact for affine elements):
    !!   sum_{i,j} scalar%interior(i,j,iEl,p) * w_i * w_j * J_{ij}^{iEl} = 1
    !! For curved elements the equality holds up to LGL aliasing error.
    !!
    !! Requirements:
    !!   - scalar%nVar == this%nPoints
    !!   - LocatePoints must have been called against the same geometry. The
    !!     cached basis is reused when this%nCached == scalar%interp%N;
    !!     otherwise the per-point Lagrange basis is recomputed on the fly.
    !!
    !! Points that lie on a shared face/edge receive their contribution in the
    !! single element selected by LocatePoints — no S/2 face split is
    !! performed.
    implicit none
    class(Points_t),intent(in) :: this
    type(SEMQuad),intent(in) :: geometry
    class(MappedScalar2D_t),intent(inout) :: scalar
    ! Local
    integer :: p,iEl,i,j,m,n,N
    real(prec) :: J0,wi,wj
    real(prec),allocatable :: lS(:),lT(:)
    logical :: useCache

    if(this%nDim /= 2) then
      print*,"SELF_Points_t::DiracDelta (2D): nDim must be 2"
      stop 1
    endif
    if(scalar%nVar /= this%nPoints) then
      print*,"SELF_Points_t::DiracDelta (2D): scalar%nVar (",scalar%nVar, &
        ") must equal nPoints (",this%nPoints,")"
      stop 1
    endif

    N = scalar%interp%N
    useCache = (this%nCached == N) .and. associated(this%lS_cache) .and. &
               associated(this%lT_cache)

    ! Each variable column is owned exclusively by a single point. Zero the
    ! field and only fill the containing element for points we located.
    scalar%interior = 0.0_prec

    if(this%nPoints == 0) return

    allocate(lS(0:N),lT(0:N))

    do p = 1,this%nPoints
      iEl = this%elements(p)
      if(iEl <= 0) cycle

      if(useCache) then
        lS = this%lS_cache(:,p)
        lT = this%lT_cache(:,p)
      else
        lS = scalar%interp%CalculateLagrangePolynomials(this%coordinates(p,1))
        lT = scalar%interp%CalculateLagrangePolynomials(this%coordinates(p,2))
      endif

      ! Polynomial interpolation of the nodal Jacobian determinant at the
      ! source's reference coordinates.
      J0 = 0.0_prec
      do n = 1,N+1
        do m = 1,N+1
          J0 = J0+lS(m-1)*lT(n-1)*geometry%J%interior(m,n,iEl,1)
        enddo
      enddo
      if(abs(J0) < tiny(1.0_prec)*1.0e3_prec) then
        print*,"SELF_Points_t::DiracDelta (2D): vanishing Jacobian at point ",p
        stop 1
      endif

      ! Rank-1 scatter into the (i,j) tensor product, post-mass-matrix form.
      do j = 1,N+1
        wj = scalar%interp%qWeights(j)
        do i = 1,N+1
          wi = scalar%interp%qWeights(i)
          scalar%interior(i,j,iEl,p) = lS(i-1)*lT(j-1)/(wi*wj*J0)
        enddo
      enddo
    enddo

    deallocate(lS,lT)

  endsubroutine DiracDelta_2D_Points_t

  subroutine DiracDelta_3D_Points_t(this,geometry,scalar)
    !! Scatter a discrete Dirac delta of unit strength (S = 1) onto a 3D
    !! MappedScalar, one variable per stored point. See DiracDelta_2D_Points_t
    !! for the full specification; this is the direct 3D analogue.
    !!
    !! For the element iEl = elements(p) > 0:
    !!
    !!   scalar%interior(i,j,k,iEl,p) = l_i(xi_p) * l_j(eta_p) * l_k(zeta_p)
    !!                                / ( w_i * w_j * w_k * J_0(p) )
    !!
    !! and zero in all other elements. J_0(p) is the polynomial interpolation
    !! of geometry%J at the source point. Conservation:
    !!   sum_{i,j,k} scalar%interior(i,j,k,iEl,p) * w_i*w_j*w_k*J_{ijk} = 1
    !! (exact on affine elements).
    implicit none
    class(Points_t),intent(in) :: this
    type(SEMHex),intent(in) :: geometry
    class(MappedScalar3D_t),intent(inout) :: scalar
    ! Local
    integer :: p,iEl,i,j,k,m,n,l,N
    real(prec) :: J0,wi,wj,wk
    real(prec),allocatable :: lS(:),lT(:),lU(:)
    logical :: useCache

    if(this%nDim /= 3) then
      print*,"SELF_Points_t::DiracDelta (3D): nDim must be 3"
      stop 1
    endif
    if(scalar%nVar /= this%nPoints) then
      print*,"SELF_Points_t::DiracDelta (3D): scalar%nVar (",scalar%nVar, &
        ") must equal nPoints (",this%nPoints,")"
      stop 1
    endif

    N = scalar%interp%N
    useCache = (this%nCached == N) .and. associated(this%lS_cache) .and. &
               associated(this%lT_cache) .and. associated(this%lU_cache)

    scalar%interior = 0.0_prec

    if(this%nPoints == 0) return

    allocate(lS(0:N),lT(0:N),lU(0:N))

    do p = 1,this%nPoints
      iEl = this%elements(p)
      if(iEl <= 0) cycle

      if(useCache) then
        lS = this%lS_cache(:,p)
        lT = this%lT_cache(:,p)
        lU = this%lU_cache(:,p)
      else
        lS = scalar%interp%CalculateLagrangePolynomials(this%coordinates(p,1))
        lT = scalar%interp%CalculateLagrangePolynomials(this%coordinates(p,2))
        lU = scalar%interp%CalculateLagrangePolynomials(this%coordinates(p,3))
      endif

      ! Polynomial interpolation of the nodal Jacobian determinant.
      J0 = 0.0_prec
      do l = 1,N+1
        do n = 1,N+1
          do m = 1,N+1
            J0 = J0+lS(m-1)*lT(n-1)*lU(l-1)*geometry%J%interior(m,n,l,iEl,1)
          enddo
        enddo
      enddo
      if(abs(J0) < tiny(1.0_prec)*1.0e3_prec) then
        print*,"SELF_Points_t::DiracDelta (3D): vanishing Jacobian at point ",p
        stop 1
      endif

      do k = 1,N+1
        wk = scalar%interp%qWeights(k)
        do j = 1,N+1
          wj = scalar%interp%qWeights(j)
          do i = 1,N+1
            wi = scalar%interp%qWeights(i)
            scalar%interior(i,j,k,iEl,p) = lS(i-1)*lT(j-1)*lU(k-1)/ &
                                           (wi*wj*wk*J0)
          enddo
        enddo
      enddo
    enddo

    deallocate(lS,lT,lU)

  endsubroutine DiracDelta_3D_Points_t

  ! ===== Internal helpers =====================================================

  pure function ClampCell(rIdx,nCells) result(c)
    !! Convert a (signed, possibly out-of-range) real cell index to an integer
    !! cell index clamped to [0, nCells-1].
    implicit none
    real(prec),intent(in) :: rIdx
    integer,intent(in) :: nCells
    integer :: c

    c = int(floor(rIdx))
    if(c < 0) c = 0
    if(c > nCells-1) c = nCells-1

  endfunction ClampCell

  subroutine BuildElementBBoxes_2D(geometry,N,nElem,bbMin,bbMax)
    !! Axis-aligned bounding box of geometry%x%interior nodes for each element.
    implicit none
    type(SEMQuad),intent(in) :: geometry
    integer,intent(in) :: N,nElem
    real(prec),intent(out) :: bbMin(1:2,1:nElem),bbMax(1:2,1:nElem)
    ! Local
    integer :: iEl,i,j
    real(prec) :: xv,yv

    do iEl = 1,nElem
      bbMin(1,iEl) = huge(1.0_prec)
      bbMin(2,iEl) = huge(1.0_prec)
      bbMax(1,iEl) = -huge(1.0_prec)
      bbMax(2,iEl) = -huge(1.0_prec)
      do j = 1,N+1
        do i = 1,N+1
          xv = geometry%x%interior(i,j,iEl,1,1)
          yv = geometry%x%interior(i,j,iEl,1,2)
          if(xv < bbMin(1,iEl)) bbMin(1,iEl) = xv
          if(yv < bbMin(2,iEl)) bbMin(2,iEl) = yv
          if(xv > bbMax(1,iEl)) bbMax(1,iEl) = xv
          if(yv > bbMax(2,iEl)) bbMax(2,iEl) = yv
        enddo
      enddo
    enddo

  endsubroutine BuildElementBBoxes_2D

  subroutine BuildElementBBoxes_3D(geometry,N,nElem,bbMin,bbMax)
    implicit none
    type(SEMHex),intent(in) :: geometry
    integer,intent(in) :: N,nElem
    real(prec),intent(out) :: bbMin(1:3,1:nElem),bbMax(1:3,1:nElem)
    ! Local
    integer :: iEl,i,j,k
    real(prec) :: xv,yv,zv

    do iEl = 1,nElem
      bbMin(1,iEl) = huge(1.0_prec)
      bbMin(2,iEl) = huge(1.0_prec)
      bbMin(3,iEl) = huge(1.0_prec)
      bbMax(1,iEl) = -huge(1.0_prec)
      bbMax(2,iEl) = -huge(1.0_prec)
      bbMax(3,iEl) = -huge(1.0_prec)
      do k = 1,N+1
        do j = 1,N+1
          do i = 1,N+1
            xv = geometry%x%interior(i,j,k,iEl,1,1)
            yv = geometry%x%interior(i,j,k,iEl,1,2)
            zv = geometry%x%interior(i,j,k,iEl,1,3)
            if(xv < bbMin(1,iEl)) bbMin(1,iEl) = xv
            if(yv < bbMin(2,iEl)) bbMin(2,iEl) = yv
            if(zv < bbMin(3,iEl)) bbMin(3,iEl) = zv
            if(xv > bbMax(1,iEl)) bbMax(1,iEl) = xv
            if(yv > bbMax(2,iEl)) bbMax(2,iEl) = yv
            if(zv > bbMax(3,iEl)) bbMax(3,iEl) = zv
          enddo
        enddo
      enddo
    enddo

  endsubroutine BuildElementBBoxes_3D

  subroutine NewtonInverse_2D(geometry,iEl,N,xTarget,xi,converged)
    !! Newton inverse-map for the 2D SEM element iEl: solve X(xi) = xTarget for
    !! the reference coordinate xi in R^2, where X is the high-order
    !! interpolant. The Jacobian dX/dxi is taken from geometry%dxds (the
    !! covariant basis tensor), interpolated to xi via Lagrange basis.
    !!
    !! Index convention: dxds%interior(i,j,iEl,1,r,c) = d(x_r)/d(s_c).
    implicit none
    type(SEMQuad),intent(in) :: geometry
    integer,intent(in) :: iEl,N
    real(prec),intent(in) :: xTarget(2)
    real(prec),intent(out) :: xi(2)
    logical,intent(out) :: converged
    ! Local
    integer :: iter,i,j,r,c
    real(prec) :: lS(0:N),lT(0:N)
    real(prec) :: xCur(2),jac(2,2),res(2),delta(2)
    real(prec) :: w,det

    xi(1) = 0.0_prec
    xi(2) = 0.0_prec
    converged = .false.

    do iter = 1,newtonMax
      lS = geometry%x%interp%CalculateLagrangePolynomials(xi(1))
      lT = geometry%x%interp%CalculateLagrangePolynomials(xi(2))

      xCur(1) = 0.0_prec
      xCur(2) = 0.0_prec
      jac(1,1) = 0.0_prec
      jac(1,2) = 0.0_prec
      jac(2,1) = 0.0_prec
      jac(2,2) = 0.0_prec
      do j = 1,N+1
        do i = 1,N+1
          w = lS(i-1)*lT(j-1)
          xCur(1) = xCur(1)+w*geometry%x%interior(i,j,iEl,1,1)
          xCur(2) = xCur(2)+w*geometry%x%interior(i,j,iEl,1,2)
          do c = 1,2
            do r = 1,2
              jac(r,c) = jac(r,c)+w*geometry%dxds%interior(i,j,iEl,1,r,c)
            enddo
          enddo
        enddo
      enddo

      res(1) = xTarget(1)-xCur(1)
      res(2) = xTarget(2)-xCur(2)

      det = jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
      if(abs(det) < tiny(1.0_prec)*1.0e3_prec) return ! singular; abandon

      delta(1) = (jac(2,2)*res(1)-jac(1,2)*res(2))/det
      delta(2) = (-jac(2,1)*res(1)+jac(1,1)*res(2))/det
      xi(1) = xi(1)+delta(1)
      xi(2) = xi(2)+delta(2)

      if(max(abs(delta(1)),abs(delta(2))) < newtonTolerance) then
        converged = .true.
        return
      endif

      ! Guard against runaway iterates outside a generous box around [-1,1]^2.
      if(abs(xi(1)) > 5.0_prec .or. abs(xi(2)) > 5.0_prec) return
    enddo

  endsubroutine NewtonInverse_2D

  subroutine NewtonInverse_3D(geometry,iEl,N,xTarget,xi,converged)
    !! Newton inverse-map for the 3D SEM element iEl. See NewtonInverse_2D.
    implicit none
    type(SEMHex),intent(in) :: geometry
    integer,intent(in) :: iEl,N
    real(prec),intent(in) :: xTarget(3)
    real(prec),intent(out) :: xi(3)
    logical,intent(out) :: converged
    ! Local
    integer :: iter,i,j,k,r,c
    real(prec) :: lS(0:N),lT(0:N),lU(0:N)
    real(prec) :: xCur(3),jac(3,3),res(3),delta(3),inv(3,3)
    real(prec) :: w,det

    xi(1) = 0.0_prec
    xi(2) = 0.0_prec
    xi(3) = 0.0_prec
    converged = .false.

    do iter = 1,newtonMax
      lS = geometry%x%interp%CalculateLagrangePolynomials(xi(1))
      lT = geometry%x%interp%CalculateLagrangePolynomials(xi(2))
      lU = geometry%x%interp%CalculateLagrangePolynomials(xi(3))

      xCur(1) = 0.0_prec
      xCur(2) = 0.0_prec
      xCur(3) = 0.0_prec
      do c = 1,3
        do r = 1,3
          jac(r,c) = 0.0_prec
        enddo
      enddo

      do k = 1,N+1
        do j = 1,N+1
          do i = 1,N+1
            w = lS(i-1)*lT(j-1)*lU(k-1)
            xCur(1) = xCur(1)+w*geometry%x%interior(i,j,k,iEl,1,1)
            xCur(2) = xCur(2)+w*geometry%x%interior(i,j,k,iEl,1,2)
            xCur(3) = xCur(3)+w*geometry%x%interior(i,j,k,iEl,1,3)
            do c = 1,3
              do r = 1,3
                jac(r,c) = jac(r,c)+w*geometry%dxds%interior(i,j,k,iEl,1,r,c)
              enddo
            enddo
          enddo
        enddo
      enddo

      res(1) = xTarget(1)-xCur(1)
      res(2) = xTarget(2)-xCur(2)
      res(3) = xTarget(3)-xCur(3)

      ! Cofactor expansion for the 3x3 inverse.
      inv(1,1) = jac(2,2)*jac(3,3)-jac(2,3)*jac(3,2)
      inv(1,2) = jac(1,3)*jac(3,2)-jac(1,2)*jac(3,3)
      inv(1,3) = jac(1,2)*jac(2,3)-jac(1,3)*jac(2,2)
      inv(2,1) = jac(2,3)*jac(3,1)-jac(2,1)*jac(3,3)
      inv(2,2) = jac(1,1)*jac(3,3)-jac(1,3)*jac(3,1)
      inv(2,3) = jac(1,3)*jac(2,1)-jac(1,1)*jac(2,3)
      inv(3,1) = jac(2,1)*jac(3,2)-jac(2,2)*jac(3,1)
      inv(3,2) = jac(1,2)*jac(3,1)-jac(1,1)*jac(3,2)
      inv(3,3) = jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
      det = jac(1,1)*inv(1,1)+jac(1,2)*inv(2,1)+jac(1,3)*inv(3,1)
      if(abs(det) < tiny(1.0_prec)*1.0e3_prec) return

      delta(1) = (inv(1,1)*res(1)+inv(1,2)*res(2)+inv(1,3)*res(3))/det
      delta(2) = (inv(2,1)*res(1)+inv(2,2)*res(2)+inv(2,3)*res(3))/det
      delta(3) = (inv(3,1)*res(1)+inv(3,2)*res(2)+inv(3,3)*res(3))/det
      xi(1) = xi(1)+delta(1)
      xi(2) = xi(2)+delta(2)
      xi(3) = xi(3)+delta(3)

      if(max(abs(delta(1)),abs(delta(2)),abs(delta(3))) < newtonTolerance) then
        converged = .true.
        return
      endif

      if(abs(xi(1)) > 5.0_prec .or. abs(xi(2)) > 5.0_prec .or. abs(xi(3)) > 5.0_prec) return
    enddo

  endsubroutine NewtonInverse_3D

endmodule SELF_Points_t
