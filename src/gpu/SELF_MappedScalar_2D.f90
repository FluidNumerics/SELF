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

  use SELF_MappedScalar_2D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(MappedScalar2D_t),public :: MappedScalar2D

    type(c_ptr) :: jas_gpu ! jacobian weighted scalar for gradient calculation

    ! Packed device buffers for the aggregated MPI halo exchange. The side
    ! tables are shared across fields and live on mesh%decomp; the buffers
    ! are per-field (sized by this field's variable count) and allocated
    ! lazily on the first exchange.
    type(c_ptr) :: halo_sendbuf_gpu = c_null_ptr ! packed device send buffer
    type(c_ptr) :: halo_recvbuf_gpu = c_null_ptr ! packed device receive buffer

  contains
    procedure,public :: Init => Init_MappedScalar2D
    procedure,public :: Free => Free_MappedScalar2D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar2D

    procedure,public :: SideExchange => SideExchange_MappedScalar2D
    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedScalar2D

    generic,public :: MappedGradient => MappedGradient_MappedScalar2D
    procedure,private :: MappedGradient_MappedScalar2D

    generic,public :: MappedDGGradient => MappedDGGradient_MappedScalar2D
    procedure,private :: MappedDGGradient_MappedScalar2D

  endtype MappedScalar2D

  interface
    subroutine ContravariantWeight_2D_gpu(f,dsdx,jaf,N,nvar,nel) &
      bind(c,name="ContravariantWeight_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,dsdx,jaf
      integer(c_int),value :: N,nvar,nel
    endsubroutine ContravariantWeight_2D_gpu
  endinterface

  interface
    subroutine NormalWeight_2D_gpu(fb,nhat,nscale,fbn,N,nvar,nel) &
      bind(c,name="NormalWeight_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: fb,nhat,nscale,fbn
      integer(c_int),value :: N,nvar,nel
    endsubroutine NormalWeight_2D_gpu
  endinterface

contains

  subroutine Init_MappedScalar2D(this,interp,nVar,nElem)
    implicit none
    class(MappedScalar2D),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer(c_size_t) :: workSize

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,interp%N+1,nelem,nvar), &
             this%boundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
             this%extBoundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
             this%avgBoundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
             this%boundarynormal(1:interp%N+1,1:4,1:nelem,1:2*nvar))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:nVar))

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%extBoundary = 0.0_prec
    this%avgBoundary = 0.0_prec
    this%boundarynormal = 0.0_prec

    call gpuCheck(hipMalloc(this%interior_gpu,sizeof(this%interior)))
    call gpuCheck(hipMalloc(this%boundary_gpu,sizeof(this%boundary)))
    call gpuCheck(hipMalloc(this%extBoundary_gpu,sizeof(this%extBoundary)))
    call gpuCheck(hipMalloc(this%avgBoundary_gpu,sizeof(this%avgBoundary)))
    call gpuCheck(hipMalloc(this%boundarynormal_gpu,sizeof(this%boundarynormal)))
    workSize = int(interp%N+1,c_size_t)*(interp%M+1)*nelem*nvar*prec
    call gpuCheck(hipMalloc(this%interpWork,workSize))
    workSize = int(interp%N+1,c_size_t)*(interp%N+1)*nelem*nvar*4*prec
    call gpuCheck(hipMalloc(this%jas_gpu,workSize))

    call hipblasCheck(hipblasCreate(this%blas_handle))

    call this%UpdateDevice()

  endsubroutine Init_MappedScalar2D

  subroutine Free_MappedScalar2D(this)
    implicit none
    class(MappedScalar2D),intent(inout) :: this

    this%nVar = 0
    this%nElem = 0
    this%interp => null()
    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%extBoundary)
    deallocate(this%avgBoundary)
    deallocate(this%boundarynormal)
    deallocate(this%meta)
    deallocate(this%eqn)

    call gpuCheck(hipFree(this%interior_gpu))
    call gpuCheck(hipFree(this%boundary_gpu))
    call gpuCheck(hipFree(this%extBoundary_gpu))
    call gpuCheck(hipFree(this%avgBoundary_gpu))
    call gpuCheck(hipFree(this%interpWork))
    call gpuCheck(hipFree(this%jas_gpu))
    call hipblasCheck(hipblasDestroy(this%blas_handle))

    if(c_associated(this%halo_sendbuf_gpu)) call gpuCheck(hipFree(this%halo_sendbuf_gpu))
    if(c_associated(this%halo_recvbuf_gpu)) call gpuCheck(hipFree(this%halo_recvbuf_gpu))
    this%halo_sendbuf_gpu = c_null_ptr
    this%halo_recvbuf_gpu = c_null_ptr

  endsubroutine Free_MappedScalar2D

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

    call gpuCheck(hipMemcpy(this%interior_gpu,c_loc(this%interior),sizeof(this%interior),hipMemcpyHostToDevice))

  endsubroutine SetInteriorFromEquation_MappedScalar2D

  subroutine MPIExchangeAsync_MappedScalar2D(this,mesh)
  !! Post the aggregated halo exchange: one MPI_Irecv/MPI_Isend pair per
  !! neighboring rank, carrying every (side,variable) boundary trace shared
  !! with that rank in a single packed device buffer. The pack kernel
  !! synchronizes the device before returning so the send buffer is complete
  !! before MPI_Isend is posted. Packed buffers are allocated on first use;
  !! the shared side tables are built by SideExchange before this is called.
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: n,npts,cnt,disp
    integer :: iError
    integer :: msgCount
    integer(c_size_t) :: worksize
    real(prec),pointer :: sendbuf(:)
    real(prec),pointer :: recvbuf(:)

    npts = (this%interp%N+1)*this%nvar

    if(.not. c_associated(this%halo_sendbuf_gpu)) then
      worksize = int(mesh%decomp%halo_nsides,c_size_t)* &
                 int(npts,c_size_t)*prec
      call gpuCheck(hipMalloc(this%halo_sendbuf_gpu,worksize))
      call gpuCheck(hipMalloc(this%halo_recvbuf_gpu,worksize))
    endif

    call HaloPack_2D_gpu(this%boundary_gpu,this%halo_sendbuf_gpu, &
                         mesh%decomp%halo_sides_gpu,this%interp%N,this%nvar, &
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

  endsubroutine MPIExchangeAsync_MappedScalar2D

  subroutine SideExchange_MappedScalar2D(this,mesh)
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: offset

    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)

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
                             mesh%decomp%rankid,offset,this%interp%N,this%nvar,this%nelem)

    if(mesh%decomp%mpiEnabled) then
      if(mesh%decomp%halo_nsides > 0) then
        call mesh%decomp%FinalizeMPIExchangeAsync()
        ! Unpack the received traces into extBoundary, applying side flips
        call HaloUnpack_2D_gpu(this%halo_recvbuf_gpu,this%extboundary_gpu, &
                               mesh%decomp%halo_sides_gpu,this%interp%N,this%nvar, &
                               this%nelem,mesh%decomp%halo_nsides)
      endif
    endif

  endsubroutine SideExchange_MappedScalar2D

  subroutine MappedGradient_MappedScalar2D(this,df)
  !! Calculates the gradient of a function using the strong form of the gradient
  !! in mapped coordinates.
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(c_ptr),intent(out) :: df
    ! Local
    real(prec),pointer :: f_p(:,:,:,:,:)
    type(c_ptr) :: fc

    call ContravariantWeight_2D_gpu(this%interior_gpu, &
                                    this%geometry%dsdx%interior_gpu,this%jas_gpu, &
                                    this%interp%N,this%nvar,this%nelem)

    ! From Vector divergence
    call c_f_pointer(this%jas_gpu,f_p,[this%interp%N+1,this%interp%N+1,this%nelem,2*this%nvar,2])

    fc = c_loc(f_p(1,1,1,1,1))
    call self_blas_matrixop_dim1_2d(this%interp%dMatrix_gpu,fc,df, &
                                    this%interp%N,this%interp%N,2*this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,2))
    call self_blas_matrixop_dim2_2d(this%interp%dMatrix_gpu,fc,df, &
                                    1.0_c_prec,this%interp%N,this%interp%N,2*this%nvar,this%nelem,this%blas_handle)

    f_p => null()

    call JacobianWeight_2D_gpu(df,this%geometry%J%interior_gpu,this%N,2*this%nVar,this%nelem)

  endsubroutine MappedGradient_MappedScalar2D

  subroutine MappedDGGradient_MappedScalar2D(this,df)
    !!
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(c_ptr),intent(inout) :: df
    ! Local
    real(prec),pointer :: f_p(:,:,:,:,:)
    type(c_ptr) :: fc

    call ContravariantWeight_2D_gpu(this%interior_gpu, &
                                    this%geometry%dsdx%interior_gpu,this%jas_gpu, &
                                    this%interp%N,this%nvar,this%nelem)

    ! From Vector divergence
    call c_f_pointer(this%jas_gpu,f_p,[this%interp%N+1,this%interp%N+1,this%nelem,2*this%nvar,2])

    fc = c_loc(f_p(1,1,1,1,1))
    call self_blas_matrixop_dim1_2d(this%interp%dgMatrix_gpu,fc,df, &
                                    this%interp%N,this%interp%N,2*this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,2))
    call self_blas_matrixop_dim2_2d(this%interp%dgMatrix_gpu,fc,df, &
                                    1.0_c_prec,this%interp%N,this%interp%N,2*this%nvar,this%nelem,this%blas_handle)

    f_p => null()

    ! Do the boundary terms
    call NormalWeight_2D_gpu(this%avgBoundary_gpu, &
                             this%geometry%nhat%boundary_gpu,this%geometry%nscale%boundary_gpu, &
                             this%boundarynormal_gpu, &
                             this%interp%N,this%nvar,this%nelem)

    call DG_BoundaryContribution_2D_gpu(this%interp%bmatrix_gpu,this%interp%qweights_gpu, &
                                        this%boundarynormal_gpu,df,this%interp%N,2*this%nvar,this%nelem)

    call JacobianWeight_2D_gpu(df,this%geometry%J%interior_gpu,this%N,2*this%nVar,this%nelem)

  endsubroutine MappedDGGradient_MappedScalar2D

endmodule SELF_MappedScalar_2D
