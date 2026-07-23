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

    ! Packed device buffers for the aggregated MPI halo exchange. The side
    ! tables are shared across fields and live on mesh%decomp; the buffers
    ! are per-field (sized for 2*nvar variables, since both vector components
    ! are exchanged) and allocated lazily on the first exchange.
    type(c_ptr) :: halo_sendbuf_gpu = c_null_ptr ! packed device send buffer
    type(c_ptr) :: halo_recvbuf_gpu = c_null_ptr ! packed device receive buffer

  contains
    procedure,public :: Free => Free_MappedVector2D
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

  subroutine Free_MappedVector2D(this)
    implicit none
    class(MappedVector2D),intent(inout) :: this

    call Free_Vector2D(this)

    if(c_associated(this%halo_sendbuf_gpu)) call gpuCheck(hipFree(this%halo_sendbuf_gpu))
    if(c_associated(this%halo_recvbuf_gpu)) call gpuCheck(hipFree(this%halo_recvbuf_gpu))
    this%halo_sendbuf_gpu = c_null_ptr
    this%halo_recvbuf_gpu = c_null_ptr

  endsubroutine Free_MappedVector2D

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
  !! Post the aggregated halo exchange: one MPI_Irecv/MPI_Isend pair per
  !! neighboring rank, carrying every (side,variable,component) boundary
  !! trace shared with that rank in a single packed device buffer. The
  !! boundary array is laid out with the component index outermost, so the
  !! pack/unpack kernels treat the vector as 2*nvar scalar variables. Packed
  !! buffers are allocated on first use; the shared side tables are built by
  !! SideExchange before this is called.
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: n,npts,cnt,disp
    integer :: iError
    integer :: msgCount
    integer(c_size_t) :: worksize
    real(prec),pointer :: sendbuf(:)
    real(prec),pointer :: recvbuf(:)

    npts = (this%interp%N+1)*2*this%nvar

    if(.not. c_associated(this%halo_sendbuf_gpu)) then
      worksize = int(mesh%decomp%halo_nsides,c_size_t)* &
                 int(npts,c_size_t)*prec
      call gpuCheck(hipMalloc(this%halo_sendbuf_gpu,worksize))
      call gpuCheck(hipMalloc(this%halo_recvbuf_gpu,worksize))
    endif

    call HaloPack_2D_gpu(this%boundary_gpu,this%halo_sendbuf_gpu, &
                         mesh%decomp%halo_sides_gpu,this%interp%N,2*this%nvar, &
                         this%nelem,mesh%decomp%halo_nsides)

    call c_f_pointer(this%halo_sendbuf_gpu,sendbuf,[mesh%decomp%halo_nsides*npts])
    call c_f_pointer(this%halo_recvbuf_gpu,recvbuf,[mesh%decomp%halo_nsides*npts])

    msgCount = 0
    do n = 1,mesh%decomp%halo_nnbr

      cnt = (mesh%decomp%halo_offset(n+1)-mesh%decomp%halo_offset(n))*npts
      disp = mesh%decomp%halo_offset(n)*npts

      msgCount = msgCount+1
      call MPI_IRECV(recvbuf(disp+1),cnt, &
                     mesh%decomp%mpiPrec, &
                     mesh%decomp%halo_rank(n),0, &
                     mesh%decomp%mpiComm, &
                     mesh%decomp%requests(msgCount),iError)

      msgCount = msgCount+1
      call MPI_ISEND(sendbuf(disp+1),cnt, &
                     mesh%decomp%mpiPrec, &
                     mesh%decomp%halo_rank(n),0, &
                     mesh%decomp%mpiComm, &
                     mesh%decomp%requests(msgCount),iError)

    enddo

    mesh%decomp%msgCount = msgCount

  endsubroutine MPIExchangeAsync_MappedVector2D

  subroutine SideExchange_MappedVector2D(this,mesh)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: offset

    offset = mesh%decomp%offsetElem(mesh%decomp%rankid+1)

    if(mesh%decomp%mpiEnabled) then
      if(.not. mesh%decomp%halo_built) then
        call mesh%decomp%BuildHaloExchange(mesh%sideInfo,mesh%nElem,4)
      endif
      if(mesh%decomp%halo_nsides > 0) then
        call this%MPIExchangeAsync(mesh)
      endif
    endif

    ! The local (same-rank) side exchange runs on the device while the
    ! aggregated MPI messages are in flight.
    call SideExchange_2D_gpu(this%extboundary_gpu, &
                             this%boundary_gpu,mesh%sideinfo_gpu,mesh%decomp%elemToRank_gpu, &
                             mesh%decomp%rankid,offset,this%interp%N,2*this%nvar,this%nelem)

    if(mesh%decomp%mpiEnabled) then
      if(mesh%decomp%halo_nsides > 0) then
        call mesh%decomp%FinalizeMPIExchangeAsync()
        ! Unpack the received traces into extBoundary, applying side flips
        call HaloUnpack_2D_gpu(this%halo_recvbuf_gpu,this%extboundary_gpu, &
                               mesh%decomp%halo_sides_gpu,this%interp%N,2*this%nvar, &
                               this%nelem,mesh%decomp%halo_nsides)
      endif
    endif

  endsubroutine SideExchange_MappedVector2D

  subroutine MappedDivergence_MappedVector2D(this,df)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(c_ptr),intent(out) :: df
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

    ! Boundary terms
    call DG_BoundaryContribution_2D_gpu(this%interp%bmatrix_gpu,this%interp%qweights_gpu, &
                                        this%boundarynormal_gpu,df,this%interp%N,this%nvar,this%nelem)

    call JacobianWeight_2D_gpu(df,this%geometry%J%interior_gpu,this%interp%N,this%nVar,this%nelem)

  endsubroutine MappedDGDivergence_MappedVector2D

endmodule SELF_MappedVector_2D
