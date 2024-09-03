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
  use SELF_MPI
  use FEQParse
  use iso_c_binding

  implicit none

  type,extends(Vector3D),public :: MappedVector3D_t
    logical :: geometry_associated = .false.
    type(SEMHex),pointer :: geometry => null()
  contains

    procedure,public :: AssociateGeometry => AssociateGeometry_MappedVector3D_t
    procedure,public :: DissociateGeometry => DissociateGeometry_MappedVector3D_t

    procedure,public :: SideExchange => SideExchange_MappedVector3D_t

    generic,public :: MappedDivergence => MappedDivergence_MappedVector3D_t
    procedure,private :: MappedDivergence_MappedVector3D_t

    generic,public :: MappedDGDivergence => MappedDGDivergence_MappedVector3D_t
    procedure,private :: MappedDGDivergence_MappedVector3D_t

    !procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedVector3D_t
    !procedure,private :: ApplyFlip => ApplyFlip_MappedVector3D_t

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedVector3D_t

    !procedure,public :: WriteTecplot => WriteTecplot_MappedVector3D_t

  endtype MappedVector3D_t

contains

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

  ! subroutine MPIExchangeAsync_MappedVector3D_t(this,decomp,mesh,resetCount)
  !   implicit none
  !   class(MappedVector3D_t),intent(inout) :: this
  !   type(MPILayer),intent(inout) :: decomp
  !   type(Mesh3D),intent(in) :: mesh
  !   logical,intent(in) :: resetCount
  !   ! Local
  !   integer :: e1,s1,e2,s2,ivar,idir
  !   integer :: globalSideId,r2
  !   integer :: iError
  !   integer :: msgCount

  !   if(decomp%mpiEnabled) then
  !     if(resetCount) then
  !       msgCount = 0
  !     else
  !       msgCount = decomp%msgCount
  !     endif

  !     do idir = 1,3
  !       do ivar = 1,this%nvar
  !         do e1 = 1,this%nElem
  !           do s1 = 1,6

  !             e2 = mesh%sideInfo(3,s1,e1) ! Neighbor Element
  !             if(e2 > 0) then
  !               r2 = decomp%elemToRank(e2) ! Neighbor Rank

  !               if(r2 /= decomp%rankId) then

  !                 ! to do : create unique tag for each side and each variable
  !                 ! tag = globalsideid + nglobalsides*ivar
  !                 s2 = mesh%sideInfo(4,s1,e1)/10
  !                 globalSideId = abs(mesh%sideInfo(2,s1,e1))

  !                 msgCount = msgCount+1
  !                 call MPI_IRECV(this%extBoundary(:,:,s1,e1,ivar,idir), &
  !                                (this%interp%N+1)*(this%interp%N+1), &
  !                                decomp%mpiPrec, &
  !                                r2,globalSideId, &
  !                                decomp%mpiComm, &
  !                                decomp%requests(msgCount),iError)

  !                 msgCount = msgCount+1
  !                 call MPI_ISEND(this%boundary(:,:,s1,e1,ivar,idir), &
  !                                (this%interp%N+1)*(this%interp%N+1), &
  !                                decomp%mpiPrec, &
  !                                r2,globalSideId, &
  !                                decomp%mpiComm, &
  !                                decomp%requests(msgCount),iError)
  !               endif
  !             endif

  !           enddo
  !         enddo
  !       enddo
  !     enddo

  !     decomp%msgCount = msgCount
  !   endif

  ! endsubroutine MPIExchangeAsync_MappedVector3D_t

  ! subroutine ApplyFlip_MappedVector3D_t(this,decomp,mesh)
  !   ! Apply side flips to sides where MPI exchanges took place.
  !   implicit none
  !   class(MappedVector3D_t),intent(inout) :: this
  !   type(MPILayer),intent(inout) :: decomp
  !   type(Mesh3D),intent(in) :: mesh
  !   ! Local
  !   integer :: e1,s1,e2,s2
  !   integer :: i,j,i2,j2
  !   integer :: r2,flip,ivar
  !   integer :: globalSideId
  !   integer :: bcid,idir
  !   real(prec) :: extBuff(1:this%interp%N+1,1:this%interp%N+1)

  !   if(decomp%mpiEnabled) then

  !     !$omp target
  !     !$omp teams loop collapse(4)
  !     do idir = 1,3
  !       do ivar = 1,this%nvar
  !         do e1 = 1,this%nElem
  !           do s1 = 1,6

  !             e2 = mesh%sideInfo(3,s1,e1) ! Neighbor Element
  !             s2 = mesh%sideInfo(4,s1,e1)/10
  !             bcid = mesh%sideInfo(5,s1,e1)
  !             if(s2 > 0 .or. bcid == 0) then ! Interior Element
  !               r2 = decomp%elemToRank(e2) ! Neighbor Rank

  !               if(r2 /= decomp%rankId) then

  !                 flip = mesh%sideInfo(4,s1,e1)-s2*10
  !                 globalSideId = mesh%sideInfo(2,s1,e1)

  !                 ! Need to update extBoundary with flip applied
  !                 if(flip == 1) then

  !                   do j = 1,this%interp%N+1
  !                     do i = 1,this%interp%N+1
  !                       i2 = j
  !                       j2 = this%interp%N+2-i
  !                       extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar,idir)
  !                     enddo
  !                   enddo

  !                 else if(flip == 2) then

  !                   do j = 1,this%interp%N+1
  !                     do i = 1,this%interp%N+1
  !                       i2 = this%interp%N+2-i
  !                       j2 = this%interp%N+2-j
  !                       extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar,idir)
  !                     enddo
  !                   enddo

  !                 else if(flip == 3) then

  !                   do j = 1,this%interp%N+1
  !                     do i = 1,this%interp%N+1
  !                       i2 = this%interp%N+2-j
  !                       j2 = i
  !                       extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar,idir)
  !                     enddo
  !                   enddo

  !                 else if(flip == 4) then

  !                   do j = 1,this%interp%N+1
  !                     do i = 1,this%interp%N+1
  !                       i2 = j
  !                       j2 = i
  !                       extBuff(i,j) = this%extBoundary(i2,j2,s1,e1,ivar,idir)
  !                     enddo
  !                   enddo

  !                 endif

  !                 do j = 1,this%interp%N+1
  !                   do i = 1,this%interp%N+1
  !                     this%extBoundary(i,j,s1,e1,ivar,idir) = extBuff(i,j)
  !                   enddo
  !                 enddo

  !               endif

  !             endif

  !           enddo
  !         enddo
  !       enddo
  !     enddo
  !     !$omp end target

  !   endif

  ! endsubroutine ApplyFlip_MappedVector3D_t

  subroutine SideExchange_MappedVector3D_t(this,mesh,decomp)
    implicit none
    class(MappedVector3D_t),intent(inout) :: this
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

    !call this%MPIExchangeAsync(decomp,mesh,resetCount=.true.)

    !$omp target
    !$omp teams loop collapse(4)
    do idir = 1,3
      do ivar = 1,this%nvar
        do e1 = 1,mesh%nElem
          do s1 = 1,6
            e2Global = mesh%sideInfo(3,s1,e1)
            e2 = e2Global-offset
            s2 = mesh%sideInfo(4,s1,e1)/10
            flip = mesh%sideInfo(4,s1,e1)-s2*10
            bcid = mesh%sideInfo(5,s1,e1)

            if(e2Global /= 0) then

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

    !call decomp%FinalizeMPIExchangeAsync()

    ! Apply side flips for data exchanged with MPI
    !call this%ApplyFlip(decomp,mesh)

  endsubroutine SideExchange_MappedVector3D_t

  subroutine MappedDivergence_MappedVector3D_t(this,df)
    ! Strong Form Operator
    !    !
    implicit none
    class(MappedVector3D_t),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,j,k,ii
    real(prec) :: dfLoc,Fx,Fy,Fz,Fc

    !$omp target
    !$omp teams
    !$omp loop collapse(5)
    do ivar = 1,this%nVar
      do iel = 1,this%nElem
        do k = 1,this%interp%N+1
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1

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
          enddo
        enddo
      enddo
    enddo

    !$omp loop collapse(5)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

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
          enddo
        enddo
      enddo
    enddo

    !$omp loop collapse(5)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

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
          enddo
        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target

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

    !$omp target
    !$omp teams
    !$omp loop collapse(5)
    do ivar = 1,this%nVar
      do iel = 1,this%nElem
        do k = 1,this%interp%N+1
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1

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
          enddo
        enddo
      enddo
    enddo

    !$omp loop collapse(5)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

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
          enddo
        enddo
      enddo
    enddo

    !$omp loop collapse(5)
    do ivar = 1,this%nvar
      do iel = 1,this%nelem
        do k = 1,this%N+1
          do j = 1,this%N+1
            do i = 1,this%N+1

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
          enddo
        enddo
      enddo
    enddo
    !$omp end teams
    !$omp end target

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
