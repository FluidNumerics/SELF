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

module SELF_MappedScalar_3D

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Scalar_3D
  use SELF_Tensor_3D
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_MPI
  use FEQParse
  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type,extends(Scalar3D),public :: MappedScalar3D

    type(Tensor3D) :: JaScalar ! contravariant weighted scalar
  contains

    procedure,public :: Init => Init_MappedScalar3D
    procedure,public :: Free => Free_MappedScalar3D
    procedure,public :: SideExchange => SideExchange_MappedScalar3D
    procedure,public :: BassiRebaySides => BassiRebaySides_MappedScalar3D

    procedure,public :: ContravariantWeightInterior => ContravariantWeightInterior_MappedScalar3D
    procedure,public :: ContravariantWeightAvgBoundary => ContravariantWeightAvgBoundary_MappedScalar3D

    generic,public :: Gradient => Gradient_MappedScalar3D
    procedure,private :: Gradient_MappedScalar3D

    generic,public :: BRGradient => BRGradient_MappedScalar3D
    procedure,private :: BRGradient_MappedScalar3D

    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedScalar3D
    procedure,private :: ApplyFlip => ApplyFlip_MappedScalar3D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar3D

  endtype MappedScalar3D

contains

  subroutine Init_MappedScalar3D(this,interp,nVar,nElem)
    implicit none
    class(MappedScalar3D),intent(out) :: this
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

    call this%JaScalar%Init(interp,nVar,nElem)

  endsubroutine Init_MappedScalar3D

  subroutine Free_MappedScalar3D(this)
    implicit none
    class(MappedScalar3D),intent(inout) :: this

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

    call this%JaScalar%Free()

  endsubroutine Free_MappedScalar3D

  subroutine SetInteriorFromEquation_MappedScalar3D(this,geometry,time)
  !!  Sets the this % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
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

              this%interior(i,j,k,iEl,iVar) = &
                this%eqn(iVar)%Evaluate((/x,y,z,time/))

            enddo
          enddo
        enddo
      enddo
    enddo

  endsubroutine SetInteriorFromEquation_MappedScalar3D

  subroutine MPIExchangeAsync_MappedScalar3D(this,decomp,mesh,resetCount)
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh3D),intent(in) :: mesh
    logical,intent(in) :: resetCount
    ! Local
    integer :: e1,s1,e2,s2,ivar
    integer :: globalSideId,r2
    integer :: iError
    integer :: msgCount

    if(decomp%mpiEnabled) then
      if(resetCount) then
        msgCount = 0
      else
        msgCount = decomp%msgCount
      endif

      do ivar = 1,this%nvar
        do e1 = 1,this%nElem
          do s1 = 1,6

            e2 = mesh%sideInfo(3,s1,e1) ! Neighbor Element
            if(e2 > 0) then
              r2 = decomp%elemToRank(e2) ! Neighbor Rank

              if(r2 /= decomp%rankId) then

                ! to do : create unique tag for each side and each variable
                ! tag = globalsideid + nglobalsides*ivar
                s2 = mesh%sideInfo(4,s1,e1)/10
                globalSideId = abs(mesh%sideInfo(2,s1,e1))

                msgCount = msgCount+1
                call MPI_IRECV(this%extBoundary(:,:,s1,e1,ivar), &
                               (this%interp%N+1)*(this%interp%N+1), &
                               decomp%mpiPrec, &
                               r2,globalSideId, &
                               decomp%mpiComm, &
                               decomp%requests(msgCount),iError)

                msgCount = msgCount+1
                call MPI_ISEND(this%boundary(:,:,s1,e1,ivar), &
                               (this%interp%N+1)*(this%interp%N+1), &
                               decomp%mpiPrec, &
                               r2,globalSideId, &
                               decomp%mpiComm, &
                               decomp%requests(msgCount),iError)
              endif
            endif

          enddo
        enddo
      enddo

      decomp%msgCount = msgCount
    endif

  endsubroutine MPIExchangeAsync_MappedScalar3D

  subroutine ApplyFlip_MappedScalar3D(this,decomp,mesh)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh3D),intent(in) :: mesh
    ! Local
    integer :: e1,s1,e2,s2
    integer :: i,i2,j,j2
    integer :: r2,flip,ivar
    integer :: globalSideId
    integer :: bcid
    real(prec) :: extBuff(1:this%interp%N+1,1:this%interp%N+1)

    if(decomp%mpiEnabled) then
      !$omp target map(to:mesh % sideInfo, decomp % elemToRank) map(tofrom:this % extBoundary)
      !$omp teams distribute parallel do collapse(3)
      do ivar = 1,this%nvar
        do e1 = 1,this%nElem
          do s1 = 1,6

            e2 = mesh%sideInfo(3,s1,e1) ! Neighbor Element
            s2 = mesh%sideInfo(4,s1,e1)/10
            bcid = mesh%sideInfo(5,s1,e1)
            if(s2 > 0 .or. bcid == 0) then ! Interior Element
              r2 = decomp%elemToRank(e2) ! Neighbor Rank

              if(r2 /= decomp%rankId) then

                flip = mesh%sideInfo(4,s1,e1)-s2*10
                globalSideId = mesh%sideInfo(2,s1,e1)

                ! Need to update extBoundary with flip applied
                if(flip == 1) then

                  do j = 1,this%interp%N+1
                    do i = 1,this%interp%N+1
                      i2 = j
                      j2 = this%interp%N+2-i
                      extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar)
                    enddo
                  enddo

                else if(flip == 2) then

                  do j = 1,this%interp%N+1
                    do i = 1,this%interp%N+1
                      i2 = this%interp%N+2-i
                      j2 = this%interp%N+2-j
                      extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar)
                    enddo
                  enddo

                else if(flip == 3) then

                  do j = 1,this%interp%N+1
                    do i = 1,this%interp%N+1
                      i2 = this%interp%N+2-j
                      j2 = i
                      extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar)
                    enddo
                  enddo

                else if(flip == 4) then

                  do j = 1,this%interp%N+1
                    do i = 1,this%interp%N+1
                      i2 = j
                      j2 = i
                      extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar)
                    enddo
                  enddo

                endif

                do j = 1,this%interp%N+1
                  do i = 1,this%interp%N+1
                    this%extBoundary(i,j,s1,e1,ivar) = extBuff(i,j)
                  enddo
                enddo

              endif

            endif

          enddo
        enddo
      enddo
      !$omp end target

    endif

  endsubroutine ApplyFlip_MappedScalar3D

  subroutine SideExchange_MappedScalar3D(this,mesh,decomp)
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(Mesh3D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,j1,j2,ivar
    integer :: neighborRank
    integer :: rankId,offset

    rankId = decomp%rankId
    offset = decomp%offsetElem(rankId+1)

    call this%MPIExchangeAsync(decomp,mesh,resetCount=.true.)

    !$omp target map(to: mesh % sideInfo, decomp % elemToRank) map(from: this % boundary) map(tofrom: this % extBoundary)
    !$omp teams distribute parallel do collapse(3)
    do ivar = 1,this%nvar
      do e1 = 1,mesh%nElem
        do s1 = 1,6
          e2Global = mesh%sideInfo(3,s1,e1)
          e2 = e2Global-offset
          s2 = mesh%sideInfo(4,s1,e1)/10
          flip = mesh%sideInfo(4,s1,e1)-s2*10
          bcid = mesh%sideInfo(5,s1,e1)

          ! If either s2 or e2 are equal to zero, then this is an exterior boundary and
          ! the extBoundary attribute is assigned by a boundary condition
          if(s2 > 0 .or. bcid == 0) then

            neighborRank = decomp%elemToRank(e2Global)

            if(neighborRank == decomp%rankId) then

              if(flip == 0) then

                do j1 = 1,this%interp%N+1
                  do i1 = 1,this%interp%N+1
                    this%extBoundary(i1,j1,s1,e1,ivar) = &
                      this%boundary(i1,j1,s2,e2,ivar)
                  enddo
                enddo

              else if(flip == 1) then

                do j1 = 1,this%interp%N+1
                  do i1 = 1,this%interp%N+1

                    i2 = j1
                    j2 = this%interp%N+2-i1
                    this%extBoundary(i1,j1,s1,e1,ivar) = &
                      this%boundary(i2,j2,s2,e2,ivar)

                  enddo
                enddo

              else if(flip == 2) then

                do j1 = 1,this%interp%N+1
                  do i1 = 1,this%interp%N+1
                    i2 = this%interp%N+2-i1
                    j2 = this%interp%N+2-j1
                    this%extBoundary(i1,j1,s1,e1,ivar) = &
                      this%boundary(i2,j2,s2,e2,ivar)
                  enddo
                enddo

              else if(flip == 3) then

                do j1 = 1,this%interp%N+1
                  do i1 = 1,this%interp%N+1
                    i2 = this%interp%N+2-j1
                    j2 = i1
                    this%extBoundary(i1,j1,s1,e1,ivar) = &
                      this%boundary(i2,j2,s2,e2,ivar)
                  enddo
                enddo

              else if(flip == 4) then

                do j1 = 1,this%interp%N+1
                  do i1 = 1,this%interp%N+1
                    i2 = j1
                    j2 = i1
                    this%extBoundary(i1,j1,s1,e1,ivar) = &
                      this%boundary(i2,j2,s2,e2,ivar)
                  enddo
                enddo

              endif

            endif

          endif

        enddo
      enddo
    enddo
    !$omp end target

    call decomp%FinalizeMPIExchangeAsync()

    ! Apply side flips for data exchanged with MPI
    call this%ApplyFlip(decomp,mesh)

  endsubroutine SideExchange_MappedScalar3D

  subroutine ContravariantWeightInterior_MappedScalar3D(this,geometry)
    !! Computes the scalar multiplied by the contravariant basis vectors
    !!
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    ! Local
    integer    :: i,j,k,iEl,iVar,row,col

    ! Interior
    !$omp target map(to: geometry % dsdx % interior, this % interior) map(from: this % JaScaalar % interior)
    !$omp teams distribute parallel do collapse(7) num_threads(256)
    do col = 1,3
      do row = 1,3
        do iVar = 1,this%nVar
          do iEl = 1,this%nElem
            do k = 1,this%interp%N+1
              do j = 1,this%interp%N+1
                do i = 1,this%interp%N+1

                  this%JaScalar%interior(i,j,k,iel,ivar,row,col) = geometry%dsdx%interior(i,j,k,iel,1,row,col)* &
                                                                   this%interior(i,j,k,iel,ivar)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine ContravariantWeightInterior_MappedScalar3D

  subroutine ContravariantWeightAvgBoundary_MappedScalar3D(this,geometry)
    !! Computes the scalar multiplied by the contravariant basis vectors
    !!
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    ! Local
    integer    :: i,j,k,iEl,iVar,row,col

    ! Interior
    !$omp target map(to:geometry % dsdx % boundary, this % avgBoundary) map(from: this % JaScalar % boundary)
    !$omp teams distribute parallel do collapse(7) num_threads(256)
    do col = 1,3
      do row = 1,3
        do iVar = 1,this%nVar
          do iEl = 1,this%nElem
            do k = 1,6
              do j = 1,this%interp%N+1
                do i = 1,this%interp%N+1

                  this%JaScalar%boundary(i,j,k,iel,ivar,row,col) = geometry%dsdx%boundary(i,j,k,iel,1,row,col)* &
                                                                   this%avgBoundary(i,j,k,iel,ivar)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine ContravariantWeightAvgBoundary_MappedScalar3D

  subroutine BassiRebaySides_MappedScalar3D(this)
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i,j

    !$omp target map(to:this % boundary, this % extBoundary) map(from:this % avgBoundary)
    !$omp teams distribute parallel do collapse(5) num_threads(256)
    do ivar = 1,this%nVar
      do iel = 1,this%nElem
        do iside = 1,6
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1
              this%avgBoundary(i,j,iside,iel,ivar) = 0.5_prec*( &
                                                     this%boundary(i,j,iside,iel,ivar)+ &
                                                     this%extBoundary(i,j,iside,iel,ivar))
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine BassiRebaySides_MappedScalar3D

  subroutine Gradient_MappedScalar3D(this,geometry,df)
  !! Calculates the gradient of a function using the strong form of the gradient
  !! in mapped coordinates.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:3)
    ! Local
    integer :: iEl,iVar,i,j,k,idir

    call this%ContravariantWeightInterior(geometry)
    call this%JaScalar%Divergence(df)

    !$omp target map(to: geometry % J % interior) map(tofrom: df)
    !$omp teams distribute parallel do collapse(6) num_threads(256)
    do idir = 1,3
      do iEl = 1,this%nElem
        do iVar = 1,this%nVar
          do k = 1,this%interp%N+1
            do j = 1,this%interp%N+1
              do i = 1,this%interp%N+1
                df(i,j,k,iEl,iVar,idir) = df(i,j,k,iEl,iVar,idir)/ &
                                          geometry%J%interior(i,j,k,iEl,1)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine Gradient_MappedScalar3D

  subroutine BRGradient_MappedScalar3D(this,geometry,df)
    !! Calculates the gradient of a function using the weak form of the gradient
    !! and the average boundary state.
    !! This method will compute the average boundary state from the
    !! and  attributes of
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:3)
    ! Local
    integer :: iEl,iVar,i,j,k,idir

    call this%BassiRebaySides()
    call this%ContravariantWeightInterior(geometry)
    call this%ContravariantWeightAvgBoundary(geometry)
    call this%JaScalar%DGDivergence(df)
    !$omp target map(to: geometry % J % interior) map(tofrom: df)
    !$omp teams distribute parallel do collapse(6) num_threads(256)
    do idir = 1,3
      do iEl = 1,this%nElem
        do iVar = 1,this%nVar
          do k = 1,this%interp%N+1
            do j = 1,this%interp%N+1
              do i = 1,this%interp%N+1
                df(i,j,k,iEl,iVar,idir) = df(i,j,k,iEl,iVar,idir)/ &
                                          geometry%J%interior(i,j,k,iEl,1)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target
  endsubroutine BRGradient_MappedScalar3D

endmodule SELF_MappedScalar_3D
