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

module SELF_MappedScalar_2D

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Scalar_2D
  use SELF_Tensor_2D
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use SELF_MPI
  use FEQParse
  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type,extends(Scalar2D),public :: MappedScalar2D

    !type(Tensor2D) :: JaScalar ! contravariant weighted scalar
  contains

    !procedure,public :: Init => Init_MappedScalar2D
    !procedure,public :: Free => Free_MappedScalar2D
    procedure,public :: SideExchange => SideExchange_MappedScalar2D
    procedure,public :: AverageSides => AverageSides_MappedScalar2D

    generic,public :: Gradient => Gradient_MappedScalar2D
    procedure,private :: Gradient_MappedScalar2D

    generic,public :: DGGradient => DGGradient_MappedScalar2D
    procedure,private :: DGGradient_MappedScalar2D

    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedScalar2D
    procedure,private :: ApplyFlip => ApplyFlip_MappedScalar2D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar2D

  endtype MappedScalar2D

contains

  ! subroutine Init_MappedScalar2D(this,interp,nVar,nElem)
  !   implicit none
  !   class(MappedScalar2D),intent(out) :: this
  !   type(Lagrange),intent(in),target :: interp
  !   integer,intent(in) :: nVar
  !   integer,intent(in) :: nElem

  !   this%interp => interp
  !   this%nVar = nVar
  !   this%nElem = nElem
  !   this%N = interp%N
  !   this%M = interp%M

  !   allocate(this%interior(1:interp%N+1,interp%N+1,nelem,nvar), &
  !            this%interpWork(1:interp%M+1,1:interp%N+1,1:nelem,1:nvar), &
  !            this%boundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
  !            this%extBoundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
  !            this%avgBoundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
  !            this%jumpBoundary(1:interp%N+1,1:4,1:nelem,1:nvar))

  !   !$omp target enter data map(alloc: this % interior)
  !   !$omp target enter data map(alloc: this % interpWork)
  !   !$omp target enter data map(alloc: this % boundary)
  !   !$omp target enter data map(alloc: this % extBoundary)
  !   !$omp target enter data map(alloc: this % avgBoundary)
  !   !$omp target enter data map(alloc: this % jumpBoundary)

  !   allocate(this%meta(1:nVar))
  !   allocate(this%eqn(1:nVar))

  !   call this%JaScalar%Init(interp,nVar,nElem)

  ! endsubroutine Init_MappedScalar2D

  ! subroutine Free_MappedScalar2D(this)
  !   implicit none
  !   class(MappedScalar2D),intent(inout) :: this

  !   this%nVar = 0
  !   this%nElem = 0
  !   this%interp => null()
  !   deallocate(this%interior)
  !   deallocate(this%interpWork)
  !   deallocate(this%boundary)
  !   deallocate(this%extBoundary)
  !   deallocate(this%avgBoundary)
  !   deallocate(this%jumpBoundary)
  !   deallocate(this%meta)
  !   deallocate(this%eqn)

  !   !$omp target exit data map(delete: this % interior)
  !   !$omp target exit data map(delete: this % interpWork)
  !   !$omp target exit data map(delete: this % boundary)
  !   !$omp target exit data map(delete: this % extBoundary)
  !   !$omp target exit data map(delete: this % avgBoundary)
  !   !$omp target exit data map(delete: this % jumpBoundary)

  !   call this%JaScalar%Free()

  ! endsubroutine Free_MappedScalar2D

  subroutine SetInteriorFromEquation_MappedScalar2D(this,geometry,time)
  !!  Sets the this % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar2D),intent(inout) :: this
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

  endsubroutine SetInteriorFromEquation_MappedScalar2D

  subroutine MPIExchangeAsync_MappedScalar2D(this,decomp,mesh,resetCount)
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh2D),intent(in) :: mesh
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
          do s1 = 1,4

            e2 = mesh%sideInfo(3,s1,e1) ! Neighbor Element
            if(e2 > 0) then
              r2 = decomp%elemToRank(e2) ! Neighbor Rank

              if(r2 /= decomp%rankId) then

                ! to do : create unique tag for each side and each variable
                ! tag = globalsideid + nglobalsides*ivar
                s2 = mesh%sideInfo(4,s1,e1)/10
                globalSideId = abs(mesh%sideInfo(2,s1,e1))

                msgCount = msgCount+1
                call MPI_IRECV(this%extBoundary(:,s1,e1,ivar), &
                               (this%interp%N+1), &
                               decomp%mpiPrec, &
                               r2,globalSideId, &
                               decomp%mpiComm, &
                               decomp%requests(msgCount),iError)

                msgCount = msgCount+1
                call MPI_ISEND(this%boundary(:,s1,e1,ivar), &
                               (this%interp%N+1), &
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

  endsubroutine MPIExchangeAsync_MappedScalar2D

  subroutine ApplyFlip_MappedScalar2D(this,decomp,mesh)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh2D),intent(in) :: mesh
    ! Local
    integer :: e1,s1,e2,s2
    integer :: i,i2
    integer :: r2,flip,ivar
    integer :: globalSideId
    integer :: bcid
    real(prec) :: extBuff(1:this%interp%N+1)

    if(decomp%mpiEnabled) then
      !$omp target map(to:mesh % sideInfo, decomp % elemToRank) map(tofrom:this % extBoundary)
      !$omp teams loop collapse(3)
      do ivar = 1,this%nvar
        do e1 = 1,this%nElem
          do s1 = 1,4

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
      !$omp end target

    endif

  endsubroutine ApplyFlip_MappedScalar2D

  subroutine SideExchange_MappedScalar2D(this,mesh,decomp)
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(Mesh2D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,ivar
    integer :: neighborRank
    integer :: rankId,offset

    rankId = decomp%rankId
    offset = decomp%offsetElem(rankId+1)

    call this%MPIExchangeAsync(decomp,mesh,resetCount=.true.)

    !$omp target map(to: mesh % sideInfo, decomp % elemToRank) map(from: this % boundary) map(tofrom: this % extBoundary)
    !$omp teams loop collapse(3)
    do ivar = 1,this%nvar
      do e1 = 1,mesh%nElem
        do s1 = 1,4
          e2Global = mesh%sideInfo(3,s1,e1)
          e2 = e2Global-offset
          s2 = mesh%sideInfo(4,s1,e1)/10
          flip = mesh%sideInfo(4,s1,e1)-s2*10
          bcid = mesh%sideInfo(5,s1,e1)

          if(s2 /= 0) then
            neighborRank = decomp%elemToRank(e2Global)

            if(neighborRank == decomp%rankId) then

              if(flip == 0) then

                do i1 = 1,this%interp%N+1
                  this%extBoundary(i1,s1,e1,ivar) = &
                    this%boundary(i1,s2,e2,ivar)
                enddo

              elseif(flip == 1) then

                do i1 = 1,this%interp%N+1
                  i2 = this%interp%N+2-i1
                  this%extBoundary(i1,s1,e1,ivar) = &
                    this%boundary(i2,s2,e2,ivar)
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

  endsubroutine SideExchange_MappedScalar2D

  subroutine AverageSides_MappedScalar2D(this)
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i

    !$omp target map(to:this % boundary, this % extBoundary) map(from:this % avgBoundary)
    !$omp teams loop collapse(4)
    do ivar = 1,this%nVar
      do iel = 1,this%nElem
        do iside = 1,4
          do i = 1,this%interp%N+1
            this%boundary(i,iside,iel,ivar) = 0.5_prec*( &
                                              this%boundary(i,iside,iel,ivar)+ &
                                              this%extBoundary(i,iside,iel,ivar))
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine AverageSides_MappedScalar2D

  function Gradient_MappedScalar2D(this,geometry) result(df)
  !! Calculates the gradient of a function using the strong form of the gradient
  !! in mapped coordinates.
    implicit none
    class(MappedScalar2D),intent(in) :: this
    type(SEMQuad),intent(in) :: geometry
    real(prec) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:2)
    ! Local
    integer :: iEl,iVar,i,j,ii,idir
    real(prec) :: dfdx,ja

    !$omp target map(to:geometry%J%interior,geometry%dsdx%interior,this%interior,this%interp%dMatrix) map(from:df)
    !$omp teams
    !$omp loop collapse(5)
    do idir = 1,2
      do iVar = 1,this%nVar
        do iEl = 1,this%nElem
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1

              dfdx = 0.0_prec
              do ii = 1,this%N+1
                ! dsdx(j,i) is contravariant vector i, component j
                ja = geometry%dsdx%interior(ii,j,iel,1,idir,1)
                dfdx = dfdx+this%interp%dMatrix(ii,i)*this%interior(ii,j,iel,ivar)*ja

              enddo

              df(i,j,iel,ivar,idir) = dfdx

            enddo
          enddo
        enddo
      enddo
    enddo

    !$omp loop collapse(5)
    do idir = 1,2
      do iVar = 1,this%nVar
        do iEl = 1,this%nElem
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1

              dfdx = 0.0_prec
              do ii = 1,this%N+1
                ja = geometry%dsdx%interior(i,ii,iel,1,idir,2)
                dfdx = dfdx+this%interp%dMatrix(ii,j)*this%interior(i,ii,iel,ivar)*ja
              enddo

              df(i,j,iel,ivar,idir) = (df(i,j,iel,ivar,idir)+dfdx)/geometry%J%interior(i,j,iEl,1)

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target

  endfunction Gradient_MappedScalar2D

  function DGGradient_MappedScalar2D(this,geometry) result(df)
    !!
    implicit none
    class(MappedScalar2D),intent(in) :: this
    type(SEMQuad),intent(in) :: geometry
    real(prec) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:2)
    ! Local
    integer :: iEl,iVar,i,j,ii,idir
    real(prec) :: dfdx,ja,bfl,bfr

    !$omp target map(to:geometry%J%interior,geometry%dsdx%interior,this%interior,this%interp%dgMatrix,this%interp%bmatrix,this%interp%qweights) map(from:df)
    !$omp teams
    !$omp loop collapse(5)
    do idir = 1,2
      do iVar = 1,this%nVar
        do iEl = 1,this%nElem
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1

              dfdx = 0.0_prec
              do ii = 1,this%N+1
                ! dsdx(j,i) is contravariant vector i, component j
                ja = geometry%dsdx%interior(ii,j,iel,1,idir,1)
                dfdx = dfdx+this%interp%dgMatrix(ii,i)*this%interior(ii,j,iel,ivar)*ja
              enddo

              bfl = this%boundary(j,4,iel,ivar)*geometry%dsdx%boundary(j,4,iel,1,idir,1) ! west
              bfr = this%boundary(j,2,iel,ivar)*geometry%dsdx%boundary(j,2,iel,1,idir,1) ! east
              df(i,j,iel,ivar,idir) = dfdx+(this%interp%bMatrix(i,1)*bfl+ &
                                            this%interp%bMatrix(i,2)*bfr)/this%interp%qweights(i)
            enddo
          enddo
        enddo
      enddo
    enddo

    !$omp loop collapse(5)
    do idir = 1,2
      do iVar = 1,this%nVar
        do iEl = 1,this%nElem
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1

              dfdx = 0.0_prec
              do ii = 1,this%N+1
                ! dsdx(j,i) is contravariant vector i, component j
                ja = geometry%dsdx%interior(i,ii,iel,1,idir,2)
                dfdx = dfdx+this%interp%dgMatrix(ii,j)*this%interior(i,ii,iel,ivar)*ja
              enddo

              bfl = this%boundary(i,1,iel,ivar)*geometry%dsdx%boundary(i,1,iel,1,idir,2) ! south
              bfr = this%boundary(i,3,iel,ivar)*geometry%dsdx%boundary(i,3,iel,1,idir,2) ! north
              dfdx = dfdx+(this%interp%bMatrix(j,1)*bfl+ &
                           this%interp%bMatrix(j,2)*bfr)/this%interp%qweights(j)

              df(i,j,iel,ivar,idir) = (df(i,j,iel,ivar,idir)+dfdx)/geometry%J%interior(i,j,iEl,1)
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target

  endfunction DGGradient_MappedScalar2D

endmodule SELF_MappedScalar_2D
