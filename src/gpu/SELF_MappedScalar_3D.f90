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

  use SELF_MappedScalar_3D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(MappedScalar3D_t),public :: MappedScalar3D

    type(c_ptr) :: jas_gpu ! jacobian weighted scalar for gradient calculation

    ! Packed device buffers and persistent MPI requests for the aggregated
    ! halo exchange. The side tables are shared across fields and live on
    ! mesh%decomp; the buffers and requests are per-field (sized by this
    ! field's variable count) and created lazily on the first exchange.
    type(c_ptr) :: halo_sendbuf_gpu = c_null_ptr ! packed device send buffer
    type(c_ptr) :: halo_recvbuf_gpu = c_null_ptr ! packed device receive buffer
    integer,allocatable :: halo_reqs(:) ! persistent requests; receives in 1:nnbr, sends in nnbr+1:2*nnbr
    integer :: halo_nactive = 0 ! variable count baked into halo_reqs
    integer :: halo_inflight = 0 ! variable count of the exchange in flight (0 = none)
    logical :: halo_static_done = .false. ! all variables have been exchanged at least once

  contains
    procedure,public :: Init => Init_MappedScalar3D
    procedure,public :: Free => Free_MappedScalar3D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar3D

    procedure,public :: SideExchange => SideExchange_MappedScalar3D
    procedure,public :: SideExchangeStart => SideExchangeStart_MappedScalar3D
    procedure,public :: SideExchangeFinish => SideExchangeFinish_MappedScalar3D

    generic,public :: MappedGradient => MappedGradient_MappedScalar3D
    procedure,private :: MappedGradient_MappedScalar3D

    generic,public :: MappedDGGradient => MappedDGGradient_MappedScalar3D
    procedure,private :: MappedDGGradient_MappedScalar3D

  endtype MappedScalar3D

  interface
    subroutine ContravariantWeight_3D_gpu(f,dsdx,jaf,N,nvar,nel) &
      bind(c,name="ContravariantWeight_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,dsdx,jaf
      integer(c_int),value :: N,nvar,nel
    endsubroutine ContravariantWeight_3D_gpu
  endinterface

  interface
    subroutine NormalWeight_3D_gpu(fb,nhat,nscale,fbn,N,nvar,nel) &
      bind(c,name="NormalWeight_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: fb,nhat,nscale,fbn
      integer(c_int),value :: N,nvar,nel
    endsubroutine NormalWeight_3D_gpu
  endinterface

contains

  subroutine Init_MappedScalar3D(this,interp,nVar,nElem)
    implicit none
    class(MappedScalar3D),intent(out) :: this
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

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar), &
             this%boundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%extBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%avgBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%boundarynormal(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:3*nvar))

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%extBoundary = 0.0_prec
    this%avgBoundary = 0.0_prec
    this%boundarynormal = 0.0_prec

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:nVar))

    call gpuCheck(hipMalloc(this%interior_gpu,sizeof(this%interior)))
    call gpuCheck(hipMalloc(this%boundary_gpu,sizeof(this%boundary)))
    call gpuCheck(hipMalloc(this%extBoundary_gpu,sizeof(this%extBoundary)))
    call gpuCheck(hipMalloc(this%avgBoundary_gpu,sizeof(this%avgBoundary)))
    call gpuCheck(hipMalloc(this%boundarynormal_gpu,sizeof(this%boundarynormal)))
    workSize = int(interp%N+1,c_size_t)*(interp%N+1)*(interp%N+1)*nelem*nvar*9*prec
    call gpuCheck(hipMalloc(this%jas_gpu,workSize))

    call this%UpdateDevice()

  endsubroutine Init_MappedScalar3D

  subroutine Free_MappedScalar3D(this)
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    ! Local
    integer :: n,iError
    logical :: mpiIsFinalized

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
    call gpuCheck(hipFree(this%boundarynormal_gpu))
    call gpuCheck(hipFree(this%jas_gpu))

    if(c_associated(this%halo_sendbuf_gpu)) call gpuCheck(hipFree(this%halo_sendbuf_gpu))
    if(c_associated(this%halo_recvbuf_gpu)) call gpuCheck(hipFree(this%halo_recvbuf_gpu))
    this%halo_sendbuf_gpu = c_null_ptr
    this%halo_recvbuf_gpu = c_null_ptr

    if(allocated(this%halo_reqs)) then
      ! Persistent requests can only be released while MPI is still
      ! initialized; if the mesh (and its MPI finalization) was freed first,
      ! MPI has reclaimed them already.
      call MPI_FINALIZED(mpiIsFinalized,iError)
      if(.not. mpiIsFinalized) then
        do n = 1,size(this%halo_reqs)
          call MPI_REQUEST_FREE(this%halo_reqs(n),iError)
        enddo
      endif
      deallocate(this%halo_reqs)
    endif
    this%halo_nactive = 0
    this%halo_inflight = 0
    this%halo_static_done = .false.

    if(c_associated(this%halo_sendbuf_gpu)) call gpuCheck(hipFree(this%halo_sendbuf_gpu))
    if(c_associated(this%halo_recvbuf_gpu)) call gpuCheck(hipFree(this%halo_recvbuf_gpu))
    this%halo_sendbuf_gpu = c_null_ptr
    this%halo_recvbuf_gpu = c_null_ptr

    if(allocated(this%halo_reqs)) then
      ! Persistent requests can only be released while MPI is still
      ! initialized; if the mesh (and its MPI finalization) was freed first,
      ! MPI has reclaimed them already.
      call MPI_FINALIZED(mpiIsFinalized,iError)
      if(.not. mpiIsFinalized) then
        do n = 1,size(this%halo_reqs)
          call MPI_REQUEST_FREE(this%halo_reqs(n),iError)
        enddo
      endif
      deallocate(this%halo_reqs)
    endif
    this%halo_nactive = 0
    this%halo_inflight = 0
    this%halo_static_done = .false.

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

    call gpuCheck(hipMemcpy(this%interior_gpu,c_loc(this%interior),sizeof(this%interior),hipMemcpyHostToDevice))

  endsubroutine SetInteriorFromEquation_MappedScalar3D

  subroutine SideExchangeStart_MappedScalar3D(this,mesh,nactive)
  !! Begin the aggregated halo exchange and launch the local (same-rank)
  !! side exchange kernel.
  !!
  !! One persistent MPI_Recv_init/MPI_Send_init pair per neighboring rank is
  !! (re)created the first time this is called (and again if the number of
  !! exchanged variables changes), then re-armed with MPI_Startall on every
  !! subsequent exchange; receives are always started before sends.
  !!
  !! When `nactive` is provided, only the leading `nactive` variables are
  !! exchanged - the variable index is the outermost dimension of the
  !! boundary arrays, so the leading variables form a contiguous prefix.
  !! The first exchange always carries all variables so that the boundary
  !! traces of static (non-stepped) variables are valid; static traces do
  !! not change thereafter, so later exchanges may carry the prefix only.
  !!
  !! The matching SideExchangeFinish must be called before extBoundary is
  !! read on interior faces. Between Start and Finish, callers may perform
  !! any work that does not read interior-face extBoundary entries; the MPI
  !! messages and the local exchange kernel launched here proceed
  !! concurrently with that work.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    integer,intent(in),optional :: nactive
    ! Local
    integer :: n,nnbr,npts,cnt,disp,nact
    integer :: iError
    integer :: offset
    integer(c_size_t) :: worksize
    real(prec),pointer :: sendbuf(:)
    real(prec),pointer :: recvbuf(:)

    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)

    if(mesh%decomp%mpiEnabled) then

      if(.not. mesh%decomp%halo_built) then
        call mesh%decomp%BuildHaloExchange(mesh%sideInfo,mesh%nElem,6)
      endif

      if(mesh%decomp%halo_nsides > 0) then

        nact = this%nvar
        if(present(nactive)) then
          if(this%halo_static_done .and. nactive >= 1 .and. nactive <= this%nvar) then
            nact = nactive
          endif
        endif

        if(.not. c_associated(this%halo_sendbuf_gpu)) then
          ! Buffers are sized for a full-variable exchange; prefix exchanges
          ! use the leading portion.
          worksize = int(mesh%decomp%halo_nsides,c_size_t)* &
                     int((this%interp%N+1)*(this%interp%N+1)*this%nvar,c_size_t)*prec
          call gpuCheck(hipMalloc(this%halo_sendbuf_gpu,worksize))
          call gpuCheck(hipMalloc(this%halo_recvbuf_gpu,worksize))
        endif

        ! Pack the send buffer; the pack kernel synchronizes the device
        ! before returning so the buffer is complete before the sends start.
        call HaloPack_3D_gpu(this%boundary_gpu,this%halo_sendbuf_gpu, &
                             mesh%decomp%halo_sides_gpu,this%interp%N,nact, &
                             this%nelem,mesh%decomp%halo_nsides)

        nnbr = mesh%decomp%halo_nnbr

        if(nact /= this%halo_nactive) then
          ! (Re)create the persistent requests for this variable count
          if(allocated(this%halo_reqs)) then
            do n = 1,size(this%halo_reqs)
              call MPI_REQUEST_FREE(this%halo_reqs(n),iError)
            enddo
            deallocate(this%halo_reqs)
          endif
          allocate(this%halo_reqs(1:2*nnbr))
          npts = (this%interp%N+1)*(this%interp%N+1)*nact
          call c_f_pointer(this%halo_sendbuf_gpu,sendbuf,[mesh%decomp%halo_nsides*npts])
          call c_f_pointer(this%halo_recvbuf_gpu,recvbuf,[mesh%decomp%halo_nsides*npts])
          do n = 1,nnbr
            cnt = (mesh%decomp%halo_offset(n+1)-mesh%decomp%halo_offset(n))*npts
            disp = mesh%decomp%halo_offset(n)*npts
            call MPI_RECV_INIT(recvbuf(disp+1),cnt,mesh%decomp%mpiPrec, &
                               mesh%decomp%halo_rank(n),0,mesh%decomp%mpiComm, &
                               this%halo_reqs(n),iError)
            call MPI_SEND_INIT(sendbuf(disp+1),cnt,mesh%decomp%mpiPrec, &
                               mesh%decomp%halo_rank(n),0,mesh%decomp%mpiComm, &
                               this%halo_reqs(nnbr+n),iError)
          enddo
          this%halo_nactive = nact
        endif

        ! Arm the receives before the sends
        call MPI_STARTALL(nnbr,this%halo_reqs(1:nnbr),iError)
        call MPI_STARTALL(nnbr,this%halo_reqs(nnbr+1:2*nnbr),iError)

        this%halo_inflight = nact

      endif
    endif

    ! The local (same-rank) side exchange runs on the device while the
    ! aggregated MPI messages are in flight.
    call SideExchange_3D_gpu(this%extboundary_gpu, &
                             this%boundary_gpu,mesh%sideinfo_gpu,mesh%decomp%elemToRank_gpu, &
                             mesh%decomp%rankid,offset,this%interp%N,this%nvar,this%nelem)

  endsubroutine SideExchangeStart_MappedScalar3D

  subroutine SideExchangeFinish_MappedScalar3D(this,mesh)
  !! Complete the aggregated halo exchange started by SideExchangeStart:
  !! wait on the persistent requests and unpack the received traces into
  !! extBoundary, applying side flips. The host blocks in MPI_Waitall (which
  !! also drives MPI progress) while previously launched device kernels
  !! continue to execute. No-op if no exchange is in flight.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: iError

    if(this%halo_inflight > 0) then

      call MPI_WAITALL(2*mesh%decomp%halo_nnbr,this%halo_reqs, &
                       MPI_STATUSES_IGNORE,iError)

      ! Unpack the received traces into extBoundary, applying side flips
      call HaloUnpack_3D_gpu(this%halo_recvbuf_gpu,this%extboundary_gpu, &
                             mesh%decomp%halo_sides_gpu,this%interp%N, &
                             this%halo_inflight,this%nelem,mesh%decomp%halo_nsides)

      if(this%halo_inflight == this%nvar) this%halo_static_done = .true.
      this%halo_inflight = 0

    endif

  endsubroutine SideExchangeFinish_MappedScalar3D

  subroutine SideExchange_MappedScalar3D(this,mesh)
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh

    call this%SideExchangeStart(mesh)
    call this%SideExchangeFinish(mesh)

  endsubroutine SideExchange_MappedScalar3D

  subroutine MappedGradient_MappedScalar3D(this,df)
  !! Calculates the gradient of a function using the strong form of the gradient
  !! in mapped coordinates.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(c_ptr),intent(out) :: df

    call ContravariantWeight_3D_gpu(this%interior_gpu, &
                                    this%geometry%dsdx%interior_gpu,this%jas_gpu, &
                                    this%interp%N,this%nvar,this%nelem)

    ! Strong-form divergence of the contravariant-weighted field (jas)
    call VectorDivergence_3D_gpu(this%interp%dMatrix_gpu,this%jas_gpu,df, &
                                 this%interp%N,3*this%nvar,this%nelem)

    call JacobianWeight_3D_gpu(df,this%geometry%J%interior_gpu,this%N,3*this%nVar,this%nelem)

  endsubroutine MappedGradient_MappedScalar3D

  subroutine MappedDGGradient_MappedScalar3D(this,df)
    !! Calculates the gradient of a function using the weak form of the gradient
    !! and the average boundary state.
    !! This method will compute the average boundary state from the
    !! and  attributes of
    implicit none
    class(MappedScalar3D),intent(in) :: this
    type(c_ptr),intent(inout) :: df

    call ContravariantWeight_3D_gpu(this%interior_gpu, &
                                    this%geometry%dsdx%interior_gpu,this%jas_gpu, &
                                    this%interp%N,this%nvar,this%nelem)

    ! Weak-form (DG) divergence of the contravariant-weighted field (jas)
    call VectorDivergence_3D_gpu(this%interp%dgMatrix_gpu,this%jas_gpu,df, &
                                 this%interp%N,3*this%nvar,this%nelem)

    ! Do the boundary terms
    call NormalWeight_3D_gpu(this%avgBoundary_gpu, &
                             this%geometry%nhat%boundary_gpu,this%geometry%nscale%boundary_gpu, &
                             this%boundarynormal_gpu, &
                             this%interp%N,this%nvar,this%nelem)

    call DG_BoundaryContribution_3D_gpu(this%interp%bmatrix_gpu,this%interp%qweights_gpu, &
                                        this%boundarynormal_gpu,df,this%interp%N,3*this%nvar,this%nelem)

    call JacobianWeight_3D_gpu(df,this%geometry%J%interior_gpu,this%N,3*this%nVar,this%nelem)

  endsubroutine MappedDGGradient_MappedScalar3D

endmodule SELF_MappedScalar_3D
