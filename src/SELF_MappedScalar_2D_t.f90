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

module SELF_MappedScalar_2D_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Scalar_2D
  use SELF_Tensor_2D
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use SELF_DomainDecomposition
  use FEQParse
  use iso_c_binding

  implicit none

  type,extends(Scalar2D),public :: MappedScalar2D_t
    logical :: geometry_associated = .false.
    type(SEMQuad),pointer :: geometry => null()

  contains

    procedure,public :: AssociateGeometry => AssociateGeometry_MappedScalar2D_t
    procedure,public :: DissociateGeometry => DissociateGeometry_MappedScalar2D_t

    procedure,public :: SideExchange => SideExchange_MappedScalar2D_t

    generic,public :: MappedGradient => MappedGradient_MappedScalar2D_t
    procedure,private :: MappedGradient_MappedScalar2D_t

    generic,public :: MappedDGGradient => MappedDGGradient_MappedScalar2D_t
    procedure,private :: MappedDGGradient_MappedScalar2D_t

    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedScalar2D_t
    procedure,private :: ApplyFlip => ApplyFlip_MappedScalar2D_t

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar2D_t

  endtype MappedScalar2D_t

contains

  subroutine AssociateGeometry_MappedScalar2D_t(this,geometry)
    implicit none
    class(MappedScalar2D_t),intent(inout) :: this
    type(SEMQuad),target,intent(in) :: geometry

    if(.not. associated(this%geometry)) then
      this%geometry => geometry
      this%geometry_associated = .true.
    endif

  endsubroutine AssociateGeometry_MappedScalar2D_t

  subroutine DissociateGeometry_MappedScalar2D_t(this)
    implicit none
    class(MappedScalar2D_t),intent(inout) :: this

    if(associated(this%geometry)) then
      this%geometry => null()
      this%geometry_associated = .false.
    endif

  endsubroutine DissociateGeometry_MappedScalar2D_t

  subroutine SetInteriorFromEquation_MappedScalar2D_t(this,geometry,time)
  !!  Sets the this % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar2D_t),intent(inout) :: this
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

            this%interior(i,j,iEl,iVar) = &
              this%eqn(iVar)%Evaluate((/x,y,0.0_prec,time/))

          enddo
        enddo
      enddo
    enddo

  endsubroutine SetInteriorFromEquation_MappedScalar2D_t

  subroutine MPIExchangeAsync_MappedScalar2D_t(this,mesh)
    implicit none
    class(MappedScalar2D_t),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: e1,s1,e2,s2,ivar
    integer :: globalSideId,r2,tag
    integer :: iError
    integer :: msgCount

    msgCount = 0

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
              tag = globalsideid+mesh%nUniqueSides*(ivar-1)

              msgCount = msgCount+1
              call MPI_IRECV(this%extBoundary(:,s1,e1,ivar), &
                             (this%interp%N+1), &
                             mesh%decomp%mpiPrec, &
                             r2,tag, &
                             mesh%decomp%mpiComm, &
                             mesh%decomp%requests(msgCount),iError)

              msgCount = msgCount+1
              call MPI_ISEND(this%boundary(:,s1,e1,ivar), &
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

    mesh%decomp%msgCount = msgCount

  endsubroutine MPIExchangeAsync_MappedScalar2D_t

  subroutine ApplyFlip_MappedScalar2D_t(this,mesh)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedScalar2D_t),intent(inout) :: this
    type(Mesh2D),intent(in) :: mesh
    ! Local
    integer :: e1,s1,e2,s2
    integer :: i,i2
    integer :: r2,flip,ivar
    integer :: globalSideId
    real(prec) :: extBuff(1:this%interp%N+1)

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
                  extBuff(i) = this%extBoundary(i2,s1,e1,ivar)
                enddo
                do i = 1,this%interp%N+1
                  this%extBoundary(i,s1,e1,ivar) = extBuff(i)
                enddo

              endif
            endif

          endif

        enddo
      enddo
    enddo

  endsubroutine ApplyFlip_MappedScalar2D_t

  subroutine SideExchange_MappedScalar2D_t(this,mesh)
    implicit none
    class(MappedScalar2D_t),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip
    integer :: i1,i2,ivar
    integer :: r2
    integer :: rankId,offset,N
    integer,pointer :: elemtorank(:)

    ! This mapping is needed to resolve a build error with
    ! amdflang that appears to be caused by referencing
    ! the elemToRank attribute within the do concurrent
    ! https://github.com/FluidNumerics/SELF/issues/54
    elemtorank => mesh%decomp%elemToRank(:)
    rankId = mesh%decomp%rankId
    offset = mesh%decomp%offsetElem(rankId+1)
    N = this%interp%N

    if(mesh%decomp%mpiEnabled) then
      call this%MPIExchangeAsync(mesh)
    endif

    do concurrent(s1=1:4,e1=1:mesh%nElem,ivar=1:this%nvar)

      e2Global = mesh%sideInfo(3,s1,e1)
      e2 = e2Global-offset
      s2 = mesh%sideInfo(4,s1,e1)/10
      flip = mesh%sideInfo(4,s1,e1)-s2*10

      if(e2Global > 0) then

        r2 = elemToRank(e2Global)
        if(r2 == rankId) then

          if(flip == 0) then
            do i1 = 1,N+1
              this%extBoundary(i1,s1,e1,ivar) = &
                this%boundary(i1,s2,e2,ivar)
            enddo

          elseif(flip == 1) then
            do i1 = 1,N+1
              i2 = N+2-i1
              this%extBoundary(i1,s1,e1,ivar) = &
                this%boundary(i2,s2,e2,ivar)
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

  endsubroutine SideExchange_MappedScalar2D_t

  subroutine MappedGradient_MappedScalar2D_t(this,df)
  !! Calculates the gradient of a function using the strong form of the gradient
  !! in mapped coordinates.
    implicit none
    class(MappedScalar2D_t),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:2)
    ! Local
    integer :: iEl,iVar,i,j,ii,idir
    real(prec) :: dfdx,ja

    do concurrent(i=1:this%N+1,j=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar,idir=1:2)

      dfdx = 0.0_prec
      do ii = 1,this%N+1
        ! dsdx(j,i) is contravariant vector i, component j
        ja = this%geometry%dsdx%interior(ii,j,iel,1,idir,1)
        dfdx = dfdx+this%interp%dMatrix(ii,i)*this%interior(ii,j,iel,ivar)*ja

      enddo

      df(i,j,iel,ivar,idir) = dfdx

    enddo

    do concurrent(i=1:this%N+1,j=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar,idir=1:2)

      dfdx = 0.0_prec
      do ii = 1,this%N+1
        ja = this%geometry%dsdx%interior(i,ii,iel,1,idir,2)
        dfdx = dfdx+this%interp%dMatrix(ii,j)*this%interior(i,ii,iel,ivar)*ja
      enddo

      df(i,j,iel,ivar,idir) = (df(i,j,iel,ivar,idir)+dfdx)/this%geometry%J%interior(i,j,iEl,1)

    enddo

  endsubroutine MappedGradient_MappedScalar2D_t

  subroutine MappedDGGradient_MappedScalar2D_t(this,df)
    !!
    implicit none
    class(MappedScalar2D_t),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:2)
    ! Local
    integer :: iEl,iVar,i,j,ii,idir
    real(prec) :: dfdx,dfdxb,ja,bfl,bfr

    do concurrent(i=1:this%N+1,j=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar,idir=1:2)

      dfdx = 0.0_prec
      do ii = 1,this%N+1
        ja = this%geometry%dsdx%interior(ii,j,iel,1,idir,1)
        dfdx = dfdx+this%interp%dgMatrix(ii,i)*this%interior(ii,j,iel,ivar)*ja
      enddo
      bfl = this%avgboundary(j,4,iel,ivar)*this%geometry%dsdx%boundary(j,4,iel,1,idir,1) ! west
      bfr = this%avgboundary(j,2,iel,ivar)*this%geometry%dsdx%boundary(j,2,iel,1,idir,1) ! east
      dfdxb = (this%interp%bMatrix(i,1)*bfl+this%interp%bMatrix(i,2)*bfr)/this%interp%qweights(i)
      df(i,j,iel,ivar,idir) = dfdx+dfdxb

    enddo

    do concurrent(i=1:this%N+1,j=1:this%N+1,iel=1:this%nElem,ivar=1:this%nVar,idir=1:2)

      dfdx = 0.0_prec
      do ii = 1,this%N+1
        ja = this%geometry%dsdx%interior(i,ii,iel,1,idir,2)
        dfdx = dfdx+this%interp%dgMatrix(ii,j)*this%interior(i,ii,iel,ivar)*ja
      enddo

      bfl = this%avgboundary(i,1,iel,ivar)*this%geometry%dsdx%boundary(i,1,iel,1,idir,2) ! south
      bfr = this%avgboundary(i,3,iel,ivar)*this%geometry%dsdx%boundary(i,3,iel,1,idir,2) ! north
      dfdxb = (this%interp%bMatrix(j,1)*bfl+this%interp%bMatrix(j,2)*bfr)/this%interp%qweights(j)

      df(i,j,iel,ivar,idir) = (df(i,j,iel,ivar,idir)+dfdx+dfdxb)/this%geometry%J%interior(i,j,iEl,1)

    enddo

  endsubroutine MappedDGGradient_MappedScalar2D_t

endmodule SELF_MappedScalar_2D_t
