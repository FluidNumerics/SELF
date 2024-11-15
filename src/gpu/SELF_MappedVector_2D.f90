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

  use SELF_MappedVector_2D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use SELF_GPUBLAS
  use iso_c_binding

  implicit none

  type,extends(MappedVector2D_t),public :: MappedVector2D

  contains
    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedVector2D

    procedure,public :: SideExchange => SideExchange_MappedVector2D
    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedVector2D

    generic,public :: MappedDivergence => MappedDivergence_MappedVector2D
    procedure,private :: MappedDivergence_MappedVector2D

    generic,public :: MappedDGDivergence => MappedDGDivergence_MappedVector2D
    procedure,private :: MappedDGDivergence_MappedVector2D

  endtype MappedVector2D

  interface
    subroutine ContravariantProjection_2D_gpu(f,dsdx,N,nvar,nel) &
      bind(c,name="ContravariantProjection_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,dsdx
      integer(c_int),value :: N,nvar,nel
    endsubroutine ContravariantProjection_2D_gpu
  endinterface

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

    call gpuCheck(hipMemcpy(this%interior_gpu,c_loc(this%interior),sizeof(this%interior),hipMemcpyHostToDevice))

  endsubroutine SetInteriorFromEquation_MappedVector2D

  subroutine MPIExchangeAsync_MappedVector2D(this,mesh)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: e1,s1,e2,s2,ivar,idir
    integer :: globalSideId,r2,tag
    integer :: iError
    integer :: msgCount
    real(prec),pointer :: boundary(:,:,:,:,:)
    real(prec),pointer :: extboundary(:,:,:,:,:)

    msgCount = 0
    call c_f_pointer(this%boundary_gpu,boundary,[this%interp%N+1,4,this%nelem,this%nvar,2])
    call c_f_pointer(this%extboundary_gpu,extboundary,[this%interp%N+1,4,this%nelem,this%nvar,2])

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
                call MPI_IRECV(extBoundary(:,s1,e1,ivar,idir), &
                               (this%interp%N+1), &
                               mesh%decomp%mpiPrec, &
                               r2,tag, &
                               mesh%decomp%mpiComm, &
                               mesh%decomp%requests(msgCount),iError)

                msgCount = msgCount+1
                call MPI_ISEND(boundary(:,s1,e1,ivar,idir), &
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

  endsubroutine MPIExchangeAsync_MappedVector2D

  subroutine SideExchange_MappedVector2D(this,mesh)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,ivar,idir
    integer :: neighborRank
    integer :: offset

    offset = mesh%decomp%offsetElem(mesh%decomp%rankid+1)

    if(mesh%decomp%mpiEnabled) then
      call this%MPIExchangeAsync(mesh)
    endif

    ! Do the side exchange internal to this mpi process
    call SideExchange_2D_gpu(this%extboundary_gpu, &
                             this%boundary_gpu,mesh%sideinfo_gpu,mesh%decomp%elemToRank_gpu, &
                             mesh%decomp%rankid,offset,this%interp%N,2*this%nvar,this%nelem)

    if(mesh%decomp%mpiEnabled) then
      call mesh%decomp%FinalizeMPIExchangeAsync()
      ! Apply side flips for data exchanged with MPI
      call ApplyFlip_2D_gpu(this%extboundary_gpu,mesh%sideInfo_gpu, &
                            mesh%decomp%elemToRank_gpu,mesh%decomp%rankId, &
                            offset,this%interp%N,2*this%nVar,this%nElem)
    endif

  endsubroutine SideExchange_MappedVector2D

  subroutine MappedDivergence_MappedVector2D(this,df)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(c_ptr),intent(out) :: df
    ! Local
    real(prec),pointer :: f_p(:,:,:,:,:)
    type(c_ptr) :: fc

    ! Contravariant projection
    call ContravariantProjection_2D_gpu(this%interior_gpu, &
                                        this%geometry%dsdx%interior_gpu,this%interp%N,this%nvar,this%nelem)

    call Divergence_2D_gpu(this%interior_gpu,df,this%interp%dMatrix_gpu, &
                           this%interp%N,this%nvar,this%nelem)

    call JacobianWeight_2D_gpu(df,this%geometry%J%interior_gpu,this%interp%N,this%nVar,this%nelem)

  endsubroutine MappedDivergence_MappedVector2D

  subroutine MappedDGDivergence_MappedVector2D(this,df)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(c_ptr),intent(out) :: df

    ! Contravariant projection
    call ContravariantProjection_2D_gpu(this%interior_gpu, &
                                        this%geometry%dsdx%interior_gpu,this%interp%N,this%nvar,this%nelem)

    call Divergence_2D_gpu(this%interior_gpu,df,this%interp%dgMatrix_gpu, &
                           this%interp%N,this%nvar,this%nelem)

    ! Boundary terms --> TO DO : problem here when nvar > 1
    call DG_BoundaryContribution_2D_gpu(this%interp%bmatrix_gpu,this%interp%qweights_gpu, &
                                        this%boundarynormal_gpu,df,this%interp%N,this%nvar,this%nelem)

    call JacobianWeight_2D_gpu(df,this%geometry%J%interior_gpu,this%interp%N,this%nVar,this%nelem)

  endsubroutine MappedDGDivergence_MappedVector2D

endmodule SELF_MappedVector_2D
