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

module SELF_MappedVector_3D

  use SELF_MappedVector_3D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(MappedVector3D_t),public :: MappedVector3D

    ! Packed device buffers for the aggregated MPI halo exchange. The side
    ! tables are shared across fields and live on mesh%decomp; the buffers
    ! are per-field (sized for 3*nvar variables, since all three vector
    ! components are exchanged) and allocated lazily on the first exchange.
    type(c_ptr) :: halo_sendbuf_gpu = c_null_ptr ! packed device send buffer
    type(c_ptr) :: halo_recvbuf_gpu = c_null_ptr ! packed device receive buffer

  contains

    procedure,public :: Free => Free_MappedVector3D
    procedure,public :: SideExchange => SideExchange_MappedVector3D
    procedure,public :: MPIExchangeAsync => MPIExchangeAsync_MappedVector3D

    generic,public :: MappedDivergence => MappedDivergence_MappedVector3D
    procedure,private :: MappedDivergence_MappedVector3D

    generic,public :: MappedDGDivergence => MappedDGDivergence_MappedVector3D
    procedure,private :: MappedDGDivergence_MappedVector3D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedVector3D

  endtype MappedVector3D

  interface
    subroutine ContravariantProjection_3D_gpu(f,dsdx,N,nvar,nel) &
      bind(c,name="ContravariantProjection_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,dsdx
      integer(c_int),value :: N,nvar,nel
    endsubroutine ContravariantProjection_3D_gpu
  endinterface

contains

  subroutine Free_MappedVector3D(this)
    implicit none
    class(MappedVector3D),intent(inout) :: this

    call Free_Vector3D(this)

    if(c_associated(this%halo_sendbuf_gpu)) call gpuCheck(hipFree(this%halo_sendbuf_gpu))
    if(c_associated(this%halo_recvbuf_gpu)) call gpuCheck(hipFree(this%halo_recvbuf_gpu))
    this%halo_sendbuf_gpu = c_null_ptr
    this%halo_recvbuf_gpu = c_null_ptr

  endsubroutine Free_MappedVector3D

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

    call gpuCheck(hipMemcpy(this%interior_gpu,c_loc(this%interior),sizeof(this%interior),hipMemcpyHostToDevice))

  endsubroutine SetInteriorFromEquation_MappedVector3D

  subroutine MPIExchangeAsync_MappedVector3D(this,mesh)
  !! Post the aggregated halo exchange: one MPI_Irecv/MPI_Isend pair per
  !! neighboring rank, carrying every (side,variable,component) boundary
  !! trace shared with that rank in a single packed device buffer. The
  !! boundary array is laid out with the component index outermost, so the
  !! pack/unpack kernels treat the vector as 3*nvar scalar variables. Packed
  !! buffers are allocated on first use; the shared side tables are built by
  !! SideExchange before this is called.
    implicit none
    class(MappedVector3D),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: n,npts,cnt,disp
    integer :: iError
    integer :: msgCount
    integer(c_size_t) :: worksize
    real(prec),pointer :: sendbuf(:)
    real(prec),pointer :: recvbuf(:)

    npts = (this%interp%N+1)*(this%interp%N+1)*3*this%nvar

    if(.not. c_associated(this%halo_sendbuf_gpu)) then
      worksize = int(mesh%decomp%halo_nsides,c_size_t)* &
                 int(npts,c_size_t)*prec
      call gpuCheck(hipMalloc(this%halo_sendbuf_gpu,worksize))
      call gpuCheck(hipMalloc(this%halo_recvbuf_gpu,worksize))
    endif

    call HaloPack_3D_gpu(this%boundary_gpu,this%halo_sendbuf_gpu, &
                         mesh%decomp%halo_sides_gpu,this%interp%N,3*this%nvar, &
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

  endsubroutine MPIExchangeAsync_MappedVector3D

  subroutine SideExchange_MappedVector3D(this,mesh)
    implicit none
    class(MappedVector3D),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: offset

    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)

    if(mesh%decomp%mpiEnabled) then
      if(.not. mesh%decomp%halo_built) then
        call mesh%decomp%BuildHaloExchange(mesh%sideInfo,mesh%nElem,6)
      endif
      if(mesh%decomp%halo_nsides > 0) then
        call this%MPIExchangeAsync(mesh)
      endif
    endif

    ! The local (same-rank) side exchange runs on the device while the
    ! aggregated MPI messages are in flight.
    call SideExchange_3D_gpu(this%extboundary_gpu, &
                             this%boundary_gpu,mesh%sideinfo_gpu,mesh%decomp%elemToRank_gpu, &
                             mesh%decomp%rankid,offset,this%interp%N,3*this%nvar,this%nelem)

    if(mesh%decomp%mpiEnabled) then
      if(mesh%decomp%halo_nsides > 0) then
        call mesh%decomp%FinalizeMPIExchangeAsync()
        ! Unpack the received traces into extBoundary, applying side flips
        call HaloUnpack_3D_gpu(this%halo_recvbuf_gpu,this%extboundary_gpu, &
                               mesh%decomp%halo_sides_gpu,this%interp%N,3*this%nvar, &
                               this%nelem,mesh%decomp%halo_nsides)
      endif
    endif

  endsubroutine SideExchange_MappedVector3D

  subroutine MappedDivergence_MappedVector3D(this,df)
    ! Strong Form Operator
    !    !
    implicit none
    class(MappedVector3D),intent(in) :: this
    type(c_ptr),intent(inout) :: df

    ! Contravariant projection
    call ContravariantProjection_3D_gpu(this%interior_gpu, &
                                        this%geometry%dsdx%interior_gpu,this%interp%N,this%nvar,this%nelem)

    call VectorDivergence_3D_gpu(this%interp%dMatrix_gpu,this%interior_gpu,df, &
                                 this%interp%N,this%nvar,this%nelem)

    call JacobianWeight_3D_gpu(df,this%geometry%J%interior_gpu,this%interp%N,this%nVar,this%nelem)

  endsubroutine MappedDivergence_MappedVector3D

  subroutine MappedDGDivergence_MappedVector3D(this,df)
      !! Computes the divergence of a 3-D vector using the weak form
      !! On input, the  attribute of the vector
      !! is assigned and the  attribute is set to the physical
      !! directions of the vector. This method will project the vector
      !! onto the contravariant basis vectors.
    implicit none
    class(MappedVector3D),intent(in) :: this
    type(c_ptr),intent(inout) :: df

    ! Contravariant projection
    call ContravariantProjection_3D_gpu(this%interior_gpu, &
                                        this%geometry%dsdx%interior_gpu,this%interp%N,this%nvar,this%nelem)

    call VectorDivergence_3D_gpu(this%interp%dgMatrix_gpu,this%interior_gpu,df, &
                                 this%interp%N,this%nvar,this%nelem)

    ! Boundary terms
    call DG_BoundaryContribution_3D_gpu(this%interp%bmatrix_gpu,this%interp%qweights_gpu, &
                                        this%boundarynormal_gpu,df,this%interp%N,this%nvar,this%nelem)

    call JacobianWeight_3D_gpu(df,this%geometry%J%interior_gpu,this%interp%N,this%nVar,this%nelem)

  endsubroutine MappedDGDivergence_MappedVector3D

endmodule SELF_MappedVector_3D
