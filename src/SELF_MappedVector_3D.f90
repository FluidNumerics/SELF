! SELF_MappedVector.F90
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_MappedVector_3D

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

  type,extends(Vector3D),public :: MappedVector3D

  contains

    procedure,public :: SideExchange => SideExchange_MappedVector3D
    procedure,public :: BassiRebaySides => BassiRebaySides_MappedVector3D

    procedure,public :: ContravariantProjection => ContravariantProjection_MappedVector3D

    generic,public :: Divergence => Divergence_MappedVector3D
    procedure,private :: Divergence_MappedVector3D

    generic,public :: DGDivergence => DGDivergence_MappedVector3D
    procedure,private :: DGDivergence_MappedVector3D

    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedVector3D
    procedure,private :: ApplyFlip => ApplyFlip_MappedVector3D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedVector3D

  endtype MappedVector3D

contains

  subroutine SetInteriorFromEquation_MappedVector3D(this,geometry,time)
    !!  Sets the this % interior attribute using the eqn attribute,
    !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedVector3D),intent(inout) :: this
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

  endsubroutine SetInteriorFromEquation_MappedVector3D

  subroutine MPIExchangeAsync_MappedVector3D(this,decomp,mesh,resetCount)
    implicit none
    class(MappedVector3D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh3D),intent(in) :: mesh
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

      do idir = 1,3
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
                  call MPI_IRECV(this%extBoundary(:,:,s1,e1,ivar,idir), &
                                 (this%interp%N+1)*(this%interp%N+1), &
                                 decomp%mpiPrec, &
                                 r2,globalSideId, &
                                 decomp%mpiComm, &
                                 decomp%requests(msgCount),iError)

                  msgCount = msgCount+1
                  call MPI_ISEND(this%boundary(:,:,s1,e1,ivar,idir), &
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
      enddo

      decomp%msgCount = msgCount
    endif

  endsubroutine MPIExchangeAsync_MappedVector3D

  subroutine ApplyFlip_MappedVector3D(this,decomp,mesh)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedVector3D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh3D),intent(in) :: mesh
    ! Local
    integer :: e1,s1,e2,s2
    integer :: i,j,i2,j2
    integer :: r2,flip,ivar
    integer :: globalSideId
    integer :: bcid,idir
    real(prec) :: extBuff(1:this%interp%N+1,1:this%interp%N+1)

    if(decomp%mpiEnabled) then

      !$omp target map(to:mesh % sideInfo, decomp % elemToRank) map(tofrom:this % extBoundary)
      !$omp teams distribute parallel do collapse(4)
      do idir = 1,3
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
                        i2 = this%interp%N+2-j
                        j2 = i
                        extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar,idir)
                      enddo
                    enddo

                  else if(flip == 4) then

                    do j = 1,this%interp%N+1
                      do i = 1,this%interp%N+1
                        i2 = j
                        j2 = i
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
      !$omp end target

    endif

  endsubroutine ApplyFlip_MappedVector3D

  subroutine SideExchange_MappedVector3D(this,mesh,decomp)
    implicit none
    class(MappedVector3D),intent(inout) :: this
    type(Mesh3D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,j1,j2,ivar
    integer :: neighborRank
    integer :: rankId,offset
    integer :: idir

    rankId = decomp%rankId
    offset = decomp%offsetElem(rankId+1)

    call this%MPIExchangeAsync(decomp,mesh,resetCount=.true.)

    !$omp target map(to: mesh % sideInfo, decomp % elemToRank) map(from: this % boundary) map(tofrom: this % extBoundary)
    !$omp teams distribute parallel do collapse(4)
    do idir = 1,3
      do ivar = 1,this%nvar
        do e1 = 1,mesh%nElem
          do s1 = 1,6
            e2Global = mesh%sideInfo(3,s1,e1)
            e2 = e2Global-offset
            s2 = mesh%sideInfo(4,s1,e1)/10
            flip = mesh%sideInfo(4,s1,e1)-s2*10
            bcid = mesh%sideInfo(5,s1,e1)

            if(s2 > 0 .or. bcid == 0) then

              neighborRank = decomp%elemToRank(e2Global)

              if(neighborRank == decomp%rankId) then

                if(flip == 0) then

                  do j1 = 1,this%interp%N+1
                    do i1 = 1,this%interp%N+1
                      this%extBoundary(i1,j1,s1,e1,ivar,idir) = &
                        this%boundary(i1,j1,s2,e2,ivar,idir)
                    enddo
                  enddo

                else if(flip == 1) then

                  do j1 = 1,this%interp%N+1
                    do i1 = 1,this%interp%N+1

                      i2 = j1
                      j2 = this%interp%N+2-i1
                      this%extBoundary(i1,j1,s1,e1,ivar,idir) = &
                        this%boundary(i2,j2,s2,e2,ivar,idir)

                    enddo
                  enddo

                else if(flip == 2) then

                  do j1 = 1,this%interp%N+1
                    do i1 = 1,this%interp%N+1
                      i2 = this%interp%N+2-i1
                      j2 = this%interp%N+2-j1
                      this%extBoundary(i1,j1,s1,e1,ivar,idir) = &
                        this%boundary(i2,j2,s2,e2,ivar,idir)
                    enddo
                  enddo

                else if(flip == 3) then

                  do j1 = 1,this%interp%N+1
                    do i1 = 1,this%interp%N+1
                      i2 = this%interp%N+2-j1
                      j2 = i1
                      this%extBoundary(i1,j1,s1,e1,ivar,idir) = &
                        this%boundary(i2,j2,s2,e2,ivar,idir)
                    enddo
                  enddo

                else if(flip == 4) then

                  do j1 = 1,this%interp%N+1
                    do i1 = 1,this%interp%N+1
                      i2 = j1
                      j2 = i1
                      this%extBoundary(i1,j1,s1,e1,ivar,idir) = &
                        this%boundary(i2,j2,s2,e2,ivar,idir)
                    enddo
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

  endsubroutine SideExchange_MappedVector3D

  subroutine BassiRebaySides_MappedVector3D(this)
    implicit none
    class(MappedVector3D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i,j
    integer :: idir

    !$omp target map(to:this % boundary, this % extBoundary) map(from:this % avgBoundary)
    !$omp teams distribute parallel do collapse(6) num_threads(256)
    do idir = 1,3
      do ivar = 1,this%nVar
        do iel = 1,this%nElem
          do iside = 1,6
            do j = 1,this%interp%N+1
              do i = 1,this%interp%N+1
                this%avgBoundary(i,j,iside,iel,ivar,idir) = 0.5_prec*( &
                                                            this%boundary(i,j,iside,iel,ivar,idir)+ &
                                                            this%extBoundary(i,j,iside,iel,ivar,idir))
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine BassiRebaySides_MappedVector3D

  subroutine ContravariantProjection_MappedVector3D(this,geometry)
    ! Takes a vector that has physical space coordinate directions (x,y,z) and projects the vector
    ! into the the contravariant basis vector directions. Keep in mind that the contravariant basis
    ! vectors are really the Jacobian weighted contravariant basis vectors
    implicit none
    class(MappedVector3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    ! Local
    integer :: i,j,k,iEl,iVar
    real(prec) :: Fx,Fy,Fz

    ! Assume that tensor(j,i) is vector i, component j
    ! => dot product is done along first dimension
    ! to project onto computational space
    !$omp target map(to:geometry % dsdx % interior) map(tofrom:this % interior)
    !$omp teams distribute parallel do collapse(5) num_threads(256)
    do iel = 1,this%nElem
      do ivar = 1,this%nVar
        do k = 1,this%interp%N+1
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1

              Fx = this%interior(i,j,k,iEl,iVar,1)
              Fy = this%interior(i,j,k,iEl,iVar,2)
              Fz = this%interior(i,j,k,iEl,iVar,3)

              this%interior(i,j,k,iEl,iVar,1) = &
                geometry%dsdx%interior(i,j,k,iEl,1,1,1)*Fx+ &
                geometry%dsdx%interior(i,j,k,iEl,1,2,1)*Fy+ &
                geometry%dsdx%interior(i,j,k,iEl,1,3,1)*Fz

              this%interior(i,j,k,iEl,iVar,2) = &
                geometry%dsdx%interior(i,j,k,iEl,1,1,2)*Fx+ &
                geometry%dsdx%interior(i,j,k,iEl,1,2,2)*Fy+ &
                geometry%dsdx%interior(i,j,k,iEl,1,3,2)*Fz

              this%interior(i,j,k,iEl,iVar,3) = &
                geometry%dsdx%interior(i,j,k,iEl,1,1,3)*Fx+ &
                geometry%dsdx%interior(i,j,k,iEl,1,2,3)*Fy+ &
                geometry%dsdx%interior(i,j,k,iEl,1,3,3)*Fz

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine ContravariantProjection_MappedVector3D

  subroutine Divergence_MappedVector3D(this,geometry,df)
    ! Strong Form Operator
    !    !
    implicit none
    class(MappedVector3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,j,k

    ! Convert from physical to computational space
    call this%ContravariantProjection(geometry)

    ! Compute the divergence
    call this%interp%VectorDivergence_3D(this%interior, &
                                         df, &
                                         this%nvar, &
                                         this%nelem)

    !$omp target map(to: geometry % J % interior) map(tofrom: df)
    !$omp teams distribute parallel do collapse(5) num_threads(256)
    do iVar = 1,this%nVar
      do iEl = 1,this%nElem
        do k = 1,this%interp%N+1
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1
              df(i,j,k,iEl,iVar) = df(i,j,k,iEl,iVar)/ &
                                   geometry%J%interior(i,j,k,iEl,1)
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine Divergence_MappedVector3D

  subroutine DGDivergence_MappedVector3D(this,geometry,df)
      !! Computes the divergence of a 3-D vector using the weak form
      !! On input, the  attribute of the vector
      !! is assigned and the  attribute is set to the physical
      !! directions of the vector. This method will project the vector
      !! onto the contravariant basis vectors.
    implicit none
    class(MappedVector3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,j,k

    ! Convert from physical to computational space
    call this%ContravariantProjection(geometry)

    ! Compute the divergence
    call this%interp%VectorDGDivergence_3D(this%interior, &
                                           this%boundaryNormal, &
                                           df, &
                                           this%nvar, &
                                           this%nelem)
    !$omp target map(to: geometry % J % interior) map(tofrom: df)
    !$omp teams distribute parallel do collapse(5) num_threads(256)
    do iVar = 1,this%nVar
      do iEl = 1,this%nElem
        do k = 1,this%interp%N+1
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1
              df(i,j,k,iEl,iVar) = df(i,j,k,iEl,iVar)/ &
                                   geometry%J%interior(i,j,k,iEl,1)
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine DGDivergence_MappedVector3D

endmodule SELF_MappedVector_3D