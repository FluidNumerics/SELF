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

module SELF_MappedVector_2D

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Scalar_2D
  use SELF_Vector_2D
  use SELF_Tensor_2D
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use SELF_MPI
  use FEQParse
  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type,extends(Vector2D),public :: MappedVector2D

  contains

    procedure,public :: SideExchange => SideExchange_MappedVector2D
    procedure,public :: AverageSides => AverageSides_MappedVector2D

    !procedure,public :: ContravariantProjection => ContravariantProjection_MappedVector2D

    generic,public :: Divergence => Divergence_MappedVector2D
    procedure,private :: Divergence_MappedVector2D

    generic,public :: DGDivergence => DGDivergence_MappedVector2D
    procedure,private :: DGDivergence_MappedVector2D

    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedVector2D
    procedure,private :: ApplyFlip => ApplyFlip_MappedVector2D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedVector2D

  endtype MappedVector2D

contains

  subroutine SetInteriorFromEquation_MappedVector2D(this,geometry,time)
  !!  Sets the this % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedVector2D),intent(inout) :: this
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

  endsubroutine SetInteriorFromEquation_MappedVector2D

  subroutine MPIExchangeAsync_MappedVector2D(this,decomp,mesh,resetCount)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh2D),intent(in) :: mesh
    logical,intent(in) :: resetCount
    ! Local
    integer :: e1,s1,e2,s2,ivar,idir
    integer :: globalSideId,r2
    integer :: iError
    integer :: msgCount

    if(decomp%mpiEnabled) then
      if(resetCount) then
        msgCount = 0
      else
        msgCount = decomp%msgCount
      endif

      do idir = 1,2
        do ivar = 1,this%nvar
          do e1 = 1,this%nElem
            do s1 = 1,4

              e2 = mesh%sideInfo(3,s1,e1) ! Neighbor Element
              if(e2 > 0) then
                r2 = decomp%elemToRank(e2) ! Neighbor Rank

                if(r2 /= decomp%rankId) then

                  ! to do : create unique tag for each side and each variable
                  ! tag = globalsideid + nglobalsides*(ivar + nvar*idir)
                  s2 = mesh%sideInfo(4,s1,e1)/10
                  globalSideId = abs(mesh%sideInfo(2,s1,e1))

                  msgCount = msgCount+1
                  call MPI_IRECV(this%extBoundary(:,s1,e1,ivar,idir), &
                                 (this%interp%N+1), &
                                 decomp%mpiPrec, &
                                 r2,globalSideId, &
                                 decomp%mpiComm, &
                                 decomp%requests(msgCount),iError)

                  msgCount = msgCount+1
                  call MPI_ISEND(this%boundary(:,s1,e1,ivar,idir), &
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
      enddo

      decomp%msgCount = msgCount
    endif

  endsubroutine MPIExchangeAsync_MappedVector2D

  subroutine ApplyFlip_MappedVector2D(this,decomp,mesh)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh2D),intent(in) :: mesh
    ! Local
    integer :: e1,s1,e2,s2
    integer :: i,i2
    integer :: r2,flip,ivar,idir
    integer :: globalSideId
    integer :: bcid
    real(prec) :: extBuff(1:this%interp%N+1)

    if(decomp%mpiEnabled) then
      !$omp target map(to:mesh % sideInfo, decomp % elemToRank) map(tofrom:this % extBoundary)
      !$omp teams loop collapse(4)
      do idir = 1,2
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
      !$omp end target

    endif

  endsubroutine ApplyFlip_MappedVector2D

  subroutine SideExchange_MappedVector2D(this,mesh,decomp)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(Mesh2D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,ivar,idir
    integer :: neighborRank
    integer :: rankId,offset

    rankId = decomp%rankId
    offset = decomp%offsetElem(rankId+1)

    call this%MPIExchangeAsync(decomp,mesh,resetCount=.true.)
    !$omp target map(to: mesh % sideInfo, decomp % elemToRank) map(from: this % boundary) map(tofrom: this % extBoundary)
    !$omp teams loop collapse(4)
    do idir = 1,2
      do ivar = 1,this%nvar
        do e1 = 1,mesh%nElem
          do s1 = 1,4
            e2Global = mesh%sideInfo(3,s1,e1)
            e2 = e2Global-offset
            s2 = mesh%sideInfo(4,s1,e1)/10
            flip = mesh%sideInfo(4,s1,e1)-s2*10
            bcid = mesh%sideInfo(5,s1,e1)

            if(s2 > 0 .or. bcid == 0) then

              neighborRank = decomp%elemToRank(e2Global)

              if(neighborRank == decomp%rankId) then

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
        enddo
      enddo
    enddo
    !$omp end target

    call decomp%FinalizeMPIExchangeAsync()

    ! Apply side flips for data exchanged with MPI
    call this%ApplyFlip(decomp,mesh)

  endsubroutine SideExchange_MappedVector2D

  subroutine AverageSides_MappedVector2D(this)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i
    integer :: idir

    !$omp target map(to:this % boundary, this % extBoundary) map(from:this % avgBoundary)
    !$omp teams loop collapse(5)
    do idir = 1,2
      do ivar = 1,this%nVar
        do iel = 1,this%nElem
          do iside = 1,4
            do i = 1,this%interp%N+1
              this%boundary(i,iside,iel,ivar,idir) = 0.5_prec*( &
                                                     this%boundary(i,iside,iel,ivar,idir)+ &
                                                     this%extBoundary(i,iside,iel,ivar,idir))
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine AverageSides_MappedVector2D

  function Divergence_MappedVector2D(this,geometry) result(df)
    ! Strong Form Operator
    !    !
    implicit none
    class(MappedVector2D),intent(in) :: this
    type(SEMQuad),intent(in) :: geometry
    real(prec) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,j,ii
    real(prec) :: dfLoc,Fx,Fy,Fc

    !$omp target map(to:geometry%dsdx%interior,this%interior,this%interp%dMatrix) map(from:df)
    !$omp teams
    !$omp loop collapse(4)
    do ivar = 1,this%nVar
      do iel = 1,this%nElem
        do j = 1,this%interp%N+1
          do i = 1,this%interp%N+1

            dfLoc = 0.0_prec
            do ii = 1,this%N+1
              ! Convert from physical to computational space
              Fx = this%interior(ii,j,iEl,iVar,1)
              Fy = this%interior(ii,j,iEl,iVar,2)
              Fc = geometry%dsdx%interior(ii,j,iEl,1,1,1)*Fx+ &
                   geometry%dsdx%interior(ii,j,iEl,1,2,1)*Fy
              dfLoc = dfLoc+this%interp%dMatrix(ii,i)*Fc
            enddo
            dF(i,j,iel,ivar) = dfLoc

          enddo
        enddo
      enddo
    enddo

    !$omp loop collapse(4)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do j = 1,this%N+1
          do i = 1,this%N+1

            dfLoc = 0.0_prec
            do ii = 1,this%N+1
              ! Convert from physical to computational space
              Fx = this%interior(i,ii,iEl,iVar,1)
              Fy = this%interior(i,ii,iEl,iVar,2)
              Fc = geometry%dsdx%interior(i,ii,iEl,1,1,2)*Fx+ &
                   geometry%dsdx%interior(i,ii,iEl,1,2,2)*Fy
              dfLoc = dfLoc+this%interp%dMatrix(ii,j)*Fc
            enddo
            dF(i,j,iel,ivar) = (dF(i,j,iel,ivar)+dfLoc)/geometry%J%interior(i,j,iEl,1)

          enddo
        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target

  endfunction Divergence_MappedVector2D

  function DGDivergence_MappedVector2D(this,geometry) result(df)
    !! Computes the divergence of a 2-D vector using the weak form
    !! On input, the  attribute of the vector
    !! is assigned and the  attribute is set to the physical
    !! directions of the vector. This method will project the vector
    !! onto the contravariant basis vectors.
    implicit none
    class(MappedVector2D),intent(in) :: this
    type(SEMQuad),intent(in) :: geometry
    real(prec) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,j,ii
    real(prec) :: dfLoc,Fx,Fy,Fc

    !$omp target map(to:geometry%dsdx%interior,this%interior,this%boundaryNormal,this%interp%dgMatrix,this%interp%bMatrix,this%interp%qweights) map(from:df)
    !$omp teams
    !$omp loop collapse(4)
    do ivar = 1,this%nVar
      do iel = 1,this%nElem
        do j = 1,this%interp%N+1
          do i = 1,this%interp%N+1

            dfLoc = 0.0_prec
            do ii = 1,this%N+1
              ! Convert from physical to computational space
              Fx = this%interior(ii,j,iEl,iVar,1)
              Fy = this%interior(ii,j,iEl,iVar,2)
              Fc = geometry%dsdx%interior(ii,j,iEl,1,1,1)*Fx+ &
                   geometry%dsdx%interior(ii,j,iEl,1,2,1)*Fy
              dfLoc = dfLoc+this%interp%dgMatrix(ii,i)*Fc
            enddo
            dF(i,j,iel,ivar) = dfLoc+ &
                               (this%interp%bMatrix(i,2)*this%boundaryNormal(j,2,iel,ivar)+ &
                                this%interp%bMatrix(i,1)*this%boundaryNormal(j,4,iel,ivar))/ &
                               this%interp%qweights(i)

          enddo
        enddo
      enddo
    enddo
    !!$omp end teams

    !$omp loop collapse(4)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do j = 1,this%N+1
          do i = 1,this%N+1

            dfLoc = 0.0_prec
            do ii = 1,this%N+1
              ! Convert from physical to computational space
              Fx = this%interior(i,ii,iEl,iVar,1)
              Fy = this%interior(i,ii,iEl,iVar,2)
              Fc = geometry%dsdx%interior(i,ii,iEl,1,1,2)*Fx+ &
                   geometry%dsdx%interior(i,ii,iEl,1,2,2)*Fy
              dfLoc = dfLoc+this%interp%dgMatrix(ii,j)*Fc
            enddo
            dfLoc = dfLoc+ &
                    (this%interp%bMatrix(j,2)*this%boundaryNormal(i,3,iel,ivar)+ &
                     this%interp%bMatrix(j,1)*this%boundaryNormal(i,1,iel,ivar))/ &
                    this%interp%qweights(j)

            dF(i,j,iel,ivar) = (dF(i,j,iel,ivar)+dfLoc)/geometry%J%interior(i,j,iEl,1)

          enddo
        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target

    ! ! Interior components of the vector divergence
    ! floc(1:,1:,1:,1:) => f(1:,1:,1:,1:,1)
    ! call self_hipblas_matrixop_dim1_2d(this % dgMatrix,floc,df,this % N,this % N,nvars,nelems,handle)
    ! floc(1:,1:,1:,1:) => f(1:,1:,1:,1:,2)
    ! call self_hipblas_matrixop_dim2_2d(this % dgMatrix,floc,df,1.0_c_prec,this % N,this % N,nvars,nelems,handle)
    ! floc => null()

    ! ! Add the boundary contributions
    ! call VectorDGDivergence_BoundaryContribution_2D(c_loc(this % bMatrix),&
    !                                                 c_loc(this % qWeights),&
    !                                                 c_loc(bf), c_loc(df),&
    !                                                 this % N, nvars, nelems)

  endfunction DGDivergence_MappedVector2D

endmodule SELF_MappedVector_2D
