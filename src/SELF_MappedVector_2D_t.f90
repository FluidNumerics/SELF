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

module SELF_MappedVector_2D_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Vector_2D
  use SELF_Tensor_2D
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use SELF_DomainDecomposition
  use FEQParse
  use iso_c_binding

  implicit none

  type,extends(Vector2D),public :: MappedVector2D_t
    logical :: geometry_associated = .false.
    type(SEMQuad),pointer :: geometry => null()

    ! Mortar exchange work array, allocated on first use for meshes with 2:1
    ! nonconforming interfaces; same slot layout as MappedScalar2D_t%mortarBuff with a
    ! trailing physical-direction index. MortarFluxCollect reuses slots 3 and 4 at
    ! idir=1 to stage the small sides' boundaryNormal traces.
    real(prec),allocatable,dimension(:,:,:,:,:) :: mortarBuff

  contains

    procedure,public :: AssociateGeometry => AssociateGeometry_MappedVector2D_t
    procedure,public :: DissociateGeometry => DissociateGeometry_MappedVector2D_t

    procedure,public :: SideExchange => SideExchange_MappedVector2D_t
    procedure,public :: MortarExchange => MortarExchange_MappedVector2D_t
    procedure,public :: MortarFluxCollect => MortarFluxCollect_MappedVector2D_t
    procedure,private :: MPIMortarExchangeAsync => MPIMortarExchangeAsync_MappedVector2D_t
    procedure,private :: MPIMortarFluxAsync => MPIMortarFluxAsync_MappedVector2D_t

    generic,public :: MappedDivergence => MappedDivergence_MappedVector2D_t
    procedure,private :: MappedDivergence_MappedVector2D_t

    generic,public :: MappedDGDivergence => MappedDGDivergence_MappedVector2D_t
    procedure,private :: MappedDGDivergence_MappedVector2D_t

    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedVector2D_t
    procedure,private :: ApplyFlip => ApplyFlip_MappedVector2D_t

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedVector2D_t

  endtype MappedVector2D_t

contains

  subroutine AssociateGeometry_MappedVector2D_t(this,geometry)
    implicit none
    class(MappedVector2D_t),intent(inout) :: this
    type(SEMQuad),target,intent(in) :: geometry

    if(.not. associated(this%geometry)) then
      this%geometry => geometry
      this%geometry_associated = .true.
    endif

  endsubroutine AssociateGeometry_MappedVector2D_t

  subroutine DissociateGeometry_MappedVector2D_t(this)
    implicit none
    class(MappedVector2D_t),intent(inout) :: this

    if(associated(this%geometry)) then
      this%geometry => null()
      this%geometry_associated = .false.
    endif

  endsubroutine DissociateGeometry_MappedVector2D_t

  subroutine SetInteriorFromEquation_MappedVector2D_t(this,geometry,time)
  !!  Sets the this % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedVector2D_t),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,j,iEl,iVar
    real(prec) :: x
    real(prec) :: y

    do iVar = 1,this%nVar
      do iEl = 1,this%nElem
        do j = 1,this%interp%N+1
          do i = 1,this%interp%N+1

            ! Get the mesh positions
            x = geometry%x%interior(i,j,iEl,1,1)
            y = geometry%x%interior(i,j,iEl,1,2)

            this%interior(i,j,iEl,iVar,1) = &
              this%eqn(1+2*(iVar-1))%Evaluate((/x,y,0.0_prec,time/))

            this%interior(i,j,iEl,iVar,2) = &
              this%eqn(2+2*(iVar-1))%Evaluate((/x,y,0.0_prec,time/))

          enddo
        enddo
      enddo
    enddo

  endsubroutine SetInteriorFromEquation_MappedVector2D_t

  subroutine MPIExchangeAsync_MappedVector2D_t(this,mesh)
    implicit none
    class(MappedVector2D_t),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: e1,s1,e2,s2,ivar,idir
    integer :: globalSideId,r2,tag
    integer :: iError
    integer :: msgCount

    msgCount = 0

    do idir = 1,2
      do ivar = 1,this%nvar
        do e1 = 1,this%nElem
          do s1 = 1,4

            e2 = mesh%sideInfo(3,s1,e1) ! Neighbor Element
            if(e2 > 0) then
              r2 = mesh%decomp%elemToRank(e2) ! Neighbor Rank

              if(r2 /= mesh%decomp%rankId) then

                s2 = mesh%sideInfo(4,s1,e1)/10
                globalSideId = abs(mesh%sideInfo(2,s1,e1))
                ! create unique tag for each side and each variable
                tag = globalsideid+mesh%nUniqueSides*(ivar-1+this%nvar*(idir-1))

                msgCount = msgCount+1
                call MPI_IRECV(this%extBoundary(:,s1,e1,ivar,idir), &
                               (this%interp%N+1), &
                               mesh%decomp%mpiPrec, &
                               r2,tag, &
                               mesh%decomp%mpiComm, &
                               mesh%decomp%requests(msgCount),iError)

                msgCount = msgCount+1
                call MPI_ISEND(this%boundary(:,s1,e1,ivar,idir), &
                               (this%interp%N+1), &
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

  endsubroutine MPIExchangeAsync_MappedVector2D_t

  subroutine ApplyFlip_MappedVector2D_t(this,mesh)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedVector2D_t),intent(inout) :: this
    type(Mesh2D),intent(in) :: mesh
    ! Local
    integer :: e1,s1,e2,s2
    integer :: i,i2
    integer :: r2,flip,ivar,idir
    real(prec) :: extBuff(1:this%interp%N+1)

    do idir = 1,2
      do ivar = 1,this%nvar
        do e1 = 1,this%nElem
          do s1 = 1,4

            e2 = mesh%sideInfo(3,s1,e1) ! Neighbor Element (global id)

            if(e2 > 0) then ! Interior Element
              r2 = mesh%decomp%elemToRank(e2) ! Neighbor Rank

              if(r2 /= mesh%decomp%rankId) then

                s2 = mesh%sideInfo(4,s1,e1)/10
                flip = mesh%sideInfo(4,s1,e1)-s2*10

                ! Need to update extBoundary with flip applied
                if(flip == 1) then

                  do i = 1,this%interp%N+1
                    i2 = this%interp%N+2-i
                    extBuff(i) = this%extBoundary(i2,s1,e1,ivar,idir)
                  enddo
                  do i = 1,this%interp%N+1
                    this%extBoundary(i,s1,e1,ivar,idir) = extBuff(i)
                  enddo

                endif
              endif

            endif

          enddo
        enddo
      enddo
    enddo

  endsubroutine ApplyFlip_MappedVector2D_t

  subroutine SideExchange_MappedVector2D_t(this,mesh)
    implicit none
    class(MappedVector2D_t),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,ivar,idir
    integer :: r2
    integer :: rankId,offset
    integer,pointer :: elemtorank(:)

    ! This mapping is needed to resolve a build error with
    ! amdflang that appears to be caused by referencing
    ! the elemToRank attribute within the do concurrent
    ! https://github.com/FluidNumerics/SELF/issues/54
    elemtorank => mesh%decomp%elemToRank(:)

    rankId = mesh%decomp%rankId
    offset = mesh%decomp%offsetElem(rankId+1)

    if(mesh%decomp%mpiEnabled) then
      call this%MPIExchangeAsync(mesh)
    endif

    do concurrent(s1=1:4,e1=1:mesh%nElem,ivar=1:this%nvar,idir=1:2)

      e2Global = mesh%sideInfo(3,s1,e1)
      e2 = e2Global-offset
      s2 = mesh%sideInfo(4,s1,e1)/10
      flip = mesh%sideInfo(4,s1,e1)-s2*10
      bcid = mesh%sideInfo(5,s1,e1)

      if(e2Global > 0) then

        r2 = elemToRank(e2Global) ! Neighbor rank
        if(r2 == mesh%decomp%rankId) then

          if(flip == 0) then

            do i1 = 1,this%interp%N+1
              this%extBoundary(i1,s1,e1,ivar,idir) = &
                this%boundary(i1,s2,e2,ivar,idir)
            enddo

          elseif(flip == 1) then

            do i1 = 1,this%interp%N+1
              i2 = this%interp%N+2-i1
              this%extBoundary(i1,s1,e1,ivar,idir) = &
                this%boundary(i2,s2,e2,ivar,idir)
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

  endsubroutine SideExchange_MappedVector2D_t

  subroutine MPIMortarExchangeAsync_MappedVector2D_t(this,mesh)
    !! Vector analogue of the scalar mortar exchange message posting; each physical
    !! direction of each variable is exchanged as its own message.
    implicit none
    class(MappedVector2D_t),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: m,k,ivar,idir
    integer :: eB,sB,rB,eS,sS,rS
    integer :: globalSideId,tag
    integer :: offset
    integer :: iError
    integer :: msgCount

    msgCount = 0
    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)

    do idir = 1,2
      do ivar = 1,this%nvar
        do m = 1,mesh%nMortars

          eB = mesh%mortarInfo(1,m)
          sB = mesh%mortarInfo(2,m)
          rB = mesh%decomp%elemToRank(eB)

          do k = 1,2

            eS = mesh%mortarInfo(2*k+1,m)
            sS = mesh%mortarInfo(2*k+2,m)/10
            rS = mesh%decomp%elemToRank(eS)
            globalSideId = mesh%mortarInfo(6+k,m)
            tag = globalSideId+mesh%nUniqueSides*(ivar-1+this%nvar*(idir-1))

            if(rB == mesh%decomp%rankId .and. rS /= mesh%decomp%rankId) then

              msgCount = msgCount+1
              call MPI_IRECV(this%mortarBuff(:,2+k,m,ivar,idir), &
                             (this%interp%N+1), &
                             mesh%decomp%mpiPrec, &
                             rS,tag, &
                             mesh%decomp%mpiComm, &
                             mesh%decomp%requests(msgCount),iError)

              msgCount = msgCount+1
              call MPI_ISEND(this%boundary(:,sB,eB-offset,ivar,idir), &
                             (this%interp%N+1), &
                             mesh%decomp%mpiPrec, &
                             rS,tag, &
                             mesh%decomp%mpiComm, &
                             mesh%decomp%requests(msgCount),iError)

            elseif(rS == mesh%decomp%rankId .and. rB /= mesh%decomp%rankId) then

              msgCount = msgCount+1
              call MPI_IRECV(this%mortarBuff(:,k,m,ivar,idir), &
                             (this%interp%N+1), &
                             mesh%decomp%mpiPrec, &
                             rB,tag, &
                             mesh%decomp%mpiComm, &
                             mesh%decomp%requests(msgCount),iError)

              msgCount = msgCount+1
              call MPI_ISEND(this%boundary(:,sS,eS-offset,ivar,idir), &
                             (this%interp%N+1), &
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

  endsubroutine MPIMortarExchangeAsync_MappedVector2D_t

  subroutine MortarExchange_MappedVector2D_t(this,mesh)
    !! Fills the extBoundary attribute on all sides participating in a 2:1
    !! nonconforming (mortar) interface; vector analogue of the scalar MortarExchange
    !! (see MappedScalar2D_t for the algorithm description).
    implicit none
    class(MappedVector2D_t),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: m,k,ivar,idir,i,ii
    integer :: eB,sB,eS,sS,flip
    integer :: rankId,offset,N
    integer,pointer :: elemtorank(:)
    real(prec) :: fm
    real(prec) :: extBuff(1:this%interp%N+1)

    ! See https://github.com/FluidNumerics/SELF/issues/54 for the reason behind
    ! this pointer alias
    elemtorank => mesh%decomp%elemToRank(:)
    rankId = mesh%decomp%rankId
    offset = mesh%decomp%offsetElem(rankId+1)
    N = this%interp%N

    if(.not. allocated(this%mortarBuff)) then
      allocate(this%mortarBuff(1:N+1,1:4,1:mesh%nMortars,1:this%nvar,1:2))
      this%mortarBuff = 0.0_prec
    endif

    if(mesh%decomp%mpiEnabled) then
      call this%MPIMortarExchangeAsync(mesh)
    endif

    ! Stage rank-local traces in the big side's edge orientation
    do concurrent(m=1:mesh%nMortars,ivar=1:this%nvar,idir=1:2)

      eB = mesh%mortarInfo(1,m)
      if(elemtorank(eB) == rankId) then
        sB = mesh%mortarInfo(2,m)
        do i = 1,N+1
          this%mortarBuff(i,1,m,ivar,idir) = this%boundary(i,sB,eB-offset,ivar,idir)
          this%mortarBuff(i,2,m,ivar,idir) = this%boundary(i,sB,eB-offset,ivar,idir)
        enddo
      endif

      do k = 1,2
        eS = mesh%mortarInfo(2*k+1,m)
        if(elemtorank(eS) == rankId) then
          sS = mesh%mortarInfo(2*k+2,m)/10
          flip = mesh%mortarInfo(2*k+2,m)-10*sS
          if(flip == 0) then
            do i = 1,N+1
              this%mortarBuff(i,2+k,m,ivar,idir) = this%boundary(i,sS,eS-offset,ivar,idir)
            enddo
          else
            do i = 1,N+1
              this%mortarBuff(i,2+k,m,ivar,idir) = this%boundary(N+2-i,sS,eS-offset,ivar,idir)
            enddo
          endif
        endif
      enddo

    enddo

    if(mesh%decomp%mpiEnabled) then
      call mesh%decomp%FinalizeMPIExchangeAsync()

      ! Reorient small-side traces received over MPI into the big side's orientation
      do idir = 1,2
        do ivar = 1,this%nvar
          do m = 1,mesh%nMortars
            eB = mesh%mortarInfo(1,m)
            if(elemtorank(eB) == rankId) then
              do k = 1,2
                eS = mesh%mortarInfo(2*k+1,m)
                sS = mesh%mortarInfo(2*k+2,m)/10
                flip = mesh%mortarInfo(2*k+2,m)-10*sS
                if(elemtorank(eS) /= rankId .and. flip == 1) then
                  do i = 1,N+1
                    extBuff(i) = this%mortarBuff(N+2-i,2+k,m,ivar,idir)
                  enddo
                  do i = 1,N+1
                    this%mortarBuff(i,2+k,m,ivar,idir) = extBuff(i)
                  enddo
                endif
              enddo
            endif
          enddo
        enddo
      enddo
    endif

    ! Compute external states :
    !  small sides get the restricted big-side trace (exact),
    !  the big side gets the L2 projection of the small-side traces
    do concurrent(m=1:mesh%nMortars,ivar=1:this%nvar,idir=1:2)

      do k = 1,2
        eS = mesh%mortarInfo(2*k+1,m)
        if(elemtorank(eS) == rankId) then
          sS = mesh%mortarInfo(2*k+2,m)/10
          flip = mesh%mortarInfo(2*k+2,m)-10*sS
          do i = 1,N+1
            fm = 0.0_prec
            do ii = 1,N+1
              fm = fm+this%interp%mortarR(ii,i,k)*this%mortarBuff(ii,k,m,ivar,idir)
            enddo
            if(flip == 0) then
              this%extBoundary(i,sS,eS-offset,ivar,idir) = fm
            else
              this%extBoundary(N+2-i,sS,eS-offset,ivar,idir) = fm
            endif
          enddo
        endif
      enddo

      eB = mesh%mortarInfo(1,m)
      if(elemtorank(eB) == rankId) then
        sB = mesh%mortarInfo(2,m)
        do i = 1,N+1
          fm = 0.0_prec
          do k = 1,2
            do ii = 1,N+1
              fm = fm+this%interp%mortarP(ii,i,k)*this%mortarBuff(ii,2+k,m,ivar,idir)
            enddo
          enddo
          this%extBoundary(i,sB,eB-offset,ivar,idir) = fm
        enddo
      endif

    enddo

  endsubroutine MortarExchange_MappedVector2D_t

  subroutine MPIMortarFluxAsync_MappedVector2D_t(this,mesh)
    !! Posts the one-directional messages for MortarFluxCollect : each remote small
    !! side sends its boundaryNormal trace to the big side's rank.
    implicit none
    class(MappedVector2D_t),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: m,k,ivar
    integer :: eB,rB,eS,sS,rS
    integer :: globalSideId,tag
    integer :: offset
    integer :: iError
    integer :: msgCount

    msgCount = 0
    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)

    do ivar = 1,this%nvar
      do m = 1,mesh%nMortars

        eB = mesh%mortarInfo(1,m)
        rB = mesh%decomp%elemToRank(eB)

        do k = 1,2

          eS = mesh%mortarInfo(2*k+1,m)
          sS = mesh%mortarInfo(2*k+2,m)/10
          rS = mesh%decomp%elemToRank(eS)
          globalSideId = mesh%mortarInfo(6+k,m)
          tag = globalSideId+mesh%nUniqueSides*(ivar-1)

          if(rB == mesh%decomp%rankId .and. rS /= mesh%decomp%rankId) then

            msgCount = msgCount+1
            call MPI_IRECV(this%mortarBuff(:,2+k,m,ivar,1), &
                           (this%interp%N+1), &
                           mesh%decomp%mpiPrec, &
                           rS,tag, &
                           mesh%decomp%mpiComm, &
                           mesh%decomp%requests(msgCount),iError)

          elseif(rS == mesh%decomp%rankId .and. rB /= mesh%decomp%rankId) then

            msgCount = msgCount+1
            call MPI_ISEND(this%boundaryNormal(:,sS,eS-offset,ivar), &
                           (this%interp%N+1), &
                           mesh%decomp%mpiPrec, &
                           rB,tag, &
                           mesh%decomp%mpiComm, &
                           mesh%decomp%requests(msgCount),iError)

          endif

        enddo
      enddo
    enddo

    mesh%decomp%msgCount = msgCount

  endsubroutine MPIMortarFluxAsync_MappedVector2D_t

  subroutine MortarFluxCollect_MappedVector2D_t(this,mesh)
    !! Replaces the big-side boundaryNormal trace on each mortar interface with the L2
    !! projection of the two small sides' boundaryNormal traces.
    !!
    !! boundaryNormal holds the Riemann-solved surface-flux integrand f* . nHat * nScale
    !! (see BoundaryFlux in the DG models). Because the small sides' nScale is half the
    !! big side's and the sub-edge coordinate Jacobian is 1/2, the projected big-side
    !! integrand is -2 * sum_k P_k g_k, where g_k are the small-side integrands and the
    !! sign accounts for the opposing outward normals. With this choice, the discrete
    !! surface integral of the big side equals minus the sum of the small sides'
    !! discrete surface integrals to roundoff, so the mortar interface is discretely
    !! conservative. Must be called after the model's BoundaryFlux and before the flux
    !! divergence is computed.
    implicit none
    class(MappedVector2D_t),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: m,k,ivar,i,ii
    integer :: eB,sB,eS,sS,flip
    integer :: rankId,offset,N
    integer,pointer :: elemtorank(:)
    real(prec) :: fm
    real(prec) :: extBuff(1:this%interp%N+1)

    ! See https://github.com/FluidNumerics/SELF/issues/54 for the reason behind
    ! this pointer alias
    elemtorank => mesh%decomp%elemToRank(:)
    rankId = mesh%decomp%rankId
    offset = mesh%decomp%offsetElem(rankId+1)
    N = this%interp%N

    if(.not. allocated(this%mortarBuff)) then
      allocate(this%mortarBuff(1:N+1,1:4,1:mesh%nMortars,1:this%nvar,1:2))
      this%mortarBuff = 0.0_prec
    endif

    if(mesh%decomp%mpiEnabled) then
      call this%MPIMortarFluxAsync(mesh)
    endif

    ! Stage rank-local small-side integrands in the big side's edge orientation
    do concurrent(m=1:mesh%nMortars,ivar=1:this%nvar)

      do k = 1,2
        eS = mesh%mortarInfo(2*k+1,m)
        if(elemtorank(eS) == rankId) then
          sS = mesh%mortarInfo(2*k+2,m)/10
          flip = mesh%mortarInfo(2*k+2,m)-10*sS
          if(flip == 0) then
            do i = 1,N+1
              this%mortarBuff(i,2+k,m,ivar,1) = this%boundaryNormal(i,sS,eS-offset,ivar)
            enddo
          else
            do i = 1,N+1
              this%mortarBuff(i,2+k,m,ivar,1) = this%boundaryNormal(N+2-i,sS,eS-offset,ivar)
            enddo
          endif
        endif
      enddo

    enddo

    if(mesh%decomp%mpiEnabled) then
      call mesh%decomp%FinalizeMPIExchangeAsync()

      ! Reorient small-side integrands received over MPI into the big side's orientation
      do ivar = 1,this%nvar
        do m = 1,mesh%nMortars
          eB = mesh%mortarInfo(1,m)
          if(elemtorank(eB) == rankId) then
            do k = 1,2
              eS = mesh%mortarInfo(2*k+1,m)
              sS = mesh%mortarInfo(2*k+2,m)/10
              flip = mesh%mortarInfo(2*k+2,m)-10*sS
              if(elemtorank(eS) /= rankId .and. flip == 1) then
                do i = 1,N+1
                  extBuff(i) = this%mortarBuff(N+2-i,2+k,m,ivar,1)
                enddo
                do i = 1,N+1
                  this%mortarBuff(i,2+k,m,ivar,1) = extBuff(i)
                enddo
              endif
            enddo
          endif
        enddo
      enddo
    endif

    ! Project the small-side integrands onto the big side's trace space. The factor of
    ! two converts the solution-space projection (mortarP carries the 1/2 sub-edge
    ! Jacobian) into the integrand-space projection; the sign accounts for the opposing
    ! outward normals.
    do concurrent(m=1:mesh%nMortars,ivar=1:this%nvar)

      eB = mesh%mortarInfo(1,m)
      if(elemtorank(eB) == rankId) then
        sB = mesh%mortarInfo(2,m)
        do i = 1,N+1
          fm = 0.0_prec
          do k = 1,2
            do ii = 1,N+1
              fm = fm+this%interp%mortarP(ii,i,k)*this%mortarBuff(ii,2+k,m,ivar,1)
            enddo
          enddo
          this%boundaryNormal(i,sB,eB-offset,ivar) = -2.0_prec*fm
        enddo
      endif

    enddo

  endsubroutine MortarFluxCollect_MappedVector2D_t

  subroutine MappedDivergence_MappedVector2D_t(this,df)
    ! Strong Form Operator
    !    !
    implicit none
    class(MappedVector2D_t),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,j,ii
    real(prec) :: dfLoc,Fx,Fy,Fc

    do concurrent(i=1:this%N+1,j=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar)

      dfLoc = 0.0_prec
      do ii = 1,this%N+1
        ! Convert from physical to computational space
        Fx = this%interior(ii,j,iEl,iVar,1)
        Fy = this%interior(ii,j,iEl,iVar,2)
        Fc = this%geometry%dsdx%interior(ii,j,iEl,1,1,1)*Fx+ &
             this%geometry%dsdx%interior(ii,j,iEl,1,2,1)*Fy
        dfLoc = dfLoc+this%interp%dMatrix(ii,i)*Fc
      enddo
      dF(i,j,iel,ivar) = dfLoc

    enddo

    do concurrent(i=1:this%N+1,j=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar)

      dfLoc = 0.0_prec
      do ii = 1,this%N+1
        ! Convert from physical to computational space
        Fx = this%interior(i,ii,iEl,iVar,1)
        Fy = this%interior(i,ii,iEl,iVar,2)
        Fc = this%geometry%dsdx%interior(i,ii,iEl,1,1,2)*Fx+ &
             this%geometry%dsdx%interior(i,ii,iEl,1,2,2)*Fy
        dfLoc = dfLoc+this%interp%dMatrix(ii,j)*Fc
      enddo
      dF(i,j,iel,ivar) = (dF(i,j,iel,ivar)+dfLoc)/this%geometry%J%interior(i,j,iEl,1)

    enddo

  endsubroutine MappedDivergence_MappedVector2D_t

  subroutine MappedDGDivergence_MappedVector2D_t(this,df)
    !! Computes the divergence of a 2-D vector using the weak form
    !! On input, the  attribute of the vector
    !! is assigned and the  attribute is set to the physical
    !! directions of the vector. This method will project the vector
    !! onto the contravariant basis vectors.
    implicit none
    class(MappedVector2D_t),intent(in) :: this
    real(prec) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,j,ii
    real(prec) :: dfLoc,Fx,Fy,Fc

    do concurrent(i=1:this%N+1,j=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar)

      dfLoc = 0.0_prec
      do ii = 1,this%N+1
        ! Convert from physical to computational space
        Fx = this%interior(ii,j,iEl,iVar,1)
        Fy = this%interior(ii,j,iEl,iVar,2)
        Fc = this%geometry%dsdx%interior(ii,j,iEl,1,1,1)*Fx+ &
             this%geometry%dsdx%interior(ii,j,iEl,1,2,1)*Fy
        dfLoc = dfLoc+this%interp%dgMatrix(ii,i)*Fc
      enddo
      dF(i,j,iel,ivar) = dfLoc+ &
                         (this%interp%bMatrix(i,2)*this%boundaryNormal(j,2,iel,ivar)+ &
                          this%interp%bMatrix(i,1)*this%boundaryNormal(j,4,iel,ivar))/ &
                         this%interp%qweights(i)

    enddo

    do concurrent(i=1:this%N+1,j=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar)

      dfLoc = 0.0_prec
      do ii = 1,this%N+1
        ! Convert from physical to computational space
        Fx = this%interior(i,ii,iEl,iVar,1)
        Fy = this%interior(i,ii,iEl,iVar,2)
        Fc = this%geometry%dsdx%interior(i,ii,iEl,1,1,2)*Fx+ &
             this%geometry%dsdx%interior(i,ii,iEl,1,2,2)*Fy
        dfLoc = dfLoc+this%interp%dgMatrix(ii,j)*Fc
      enddo
      dfLoc = dfLoc+ &
              (this%interp%bMatrix(j,2)*this%boundaryNormal(i,3,iel,ivar)+ &
               this%interp%bMatrix(j,1)*this%boundaryNormal(i,1,iel,ivar))/ &
              this%interp%qweights(j)

      dF(i,j,iel,ivar) = (dF(i,j,iel,ivar)+dfLoc)/this%geometry%J%interior(i,j,iEl,1)

    enddo

  endsubroutine MappedDGDivergence_MappedVector2D_t

endmodule SELF_MappedVector_2D_t
