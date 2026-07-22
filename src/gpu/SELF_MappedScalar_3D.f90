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
  use iso_fortran_env,only:int64

  implicit none

  type,extends(MappedScalar3D_t),public :: MappedScalar3D

    type(c_ptr) :: jas_gpu ! jacobian weighted scalar for gradient calculation

    ! Aggregated MPI halo exchange tables. All locally-owned sides with a
    ! neighbor element on another rank are enumerated once, grouped by
    ! neighbor rank and sorted by global side id within each group, so that
    ! a single message per neighbor rank carries every (side,variable) trace.
    logical :: halo_built = .false.
    integer :: halo_nnbr = 0 ! Number of neighboring ranks
    integer :: halo_nsides = 0 ! Total number of sides exchanged with other ranks
    integer,allocatable :: halo_rank(:) ! (1:halo_nnbr) neighbor rank ids
    integer,allocatable :: halo_offset(:) ! (1:halo_nnbr+1) side-list offsets per neighbor
    type(c_ptr) :: halo_sides_gpu = c_null_ptr ! (element,side,flip) triplets, device copy
    type(c_ptr) :: halo_sendbuf_gpu = c_null_ptr ! packed device send buffer
    type(c_ptr) :: halo_recvbuf_gpu = c_null_ptr ! packed device receive buffer

  contains
    procedure,public :: Init => Init_MappedScalar3D
    procedure,public :: Free => Free_MappedScalar3D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar3D

    procedure,public :: SideExchange => SideExchange_MappedScalar3D
    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedScalar3D
    procedure,private :: SetupHaloExchange => SetupHaloExchange_MappedScalar3D

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
    workSize = int(interp%N+1,c_size_t)*(interp%N+1)*(interp%M+1)*nelem*nvar*prec
    call gpuCheck(hipMalloc(this%interpWork1,workSize))
    workSize = int(interp%N+1,c_size_t)*(interp%M+1)*(interp%M+1)*nelem*nvar*prec
    call gpuCheck(hipMalloc(this%interpWork2,workSize))
    workSize = int(interp%N+1,c_size_t)*(interp%N+1)*(interp%N+1)*nelem*nvar*9*prec
    call gpuCheck(hipMalloc(this%jas_gpu,workSize))

    call this%UpdateDevice()

    call hipblasCheck(hipblasCreate(this%blas_handle))

  endsubroutine Init_MappedScalar3D

  subroutine Free_MappedScalar3D(this)
    implicit none
    class(MappedScalar3D),intent(inout) :: this

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
    call gpuCheck(hipFree(this%interpWork1))
    call gpuCheck(hipFree(this%interpWork2))
    call gpuCheck(hipFree(this%jas_gpu))
    call hipblasCheck(hipblasDestroy(this%blas_handle))

    if(allocated(this%halo_rank)) deallocate(this%halo_rank)
    if(allocated(this%halo_offset)) deallocate(this%halo_offset)
    if(c_associated(this%halo_sides_gpu)) call gpuCheck(hipFree(this%halo_sides_gpu))
    if(c_associated(this%halo_sendbuf_gpu)) call gpuCheck(hipFree(this%halo_sendbuf_gpu))
    if(c_associated(this%halo_recvbuf_gpu)) call gpuCheck(hipFree(this%halo_recvbuf_gpu))
    this%halo_sides_gpu = c_null_ptr
    this%halo_sendbuf_gpu = c_null_ptr
    this%halo_recvbuf_gpu = c_null_ptr
    this%halo_built = .false.

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

  subroutine SetupHaloExchange_MappedScalar3D(this,mesh)
  !! Build the aggregated halo exchange tables for this field.
  !!
  !! Enumerates every locally-owned side whose neighbor element resides on
  !! another rank, groups the sides by neighbor rank, and sorts each group by
  !! global side id. Both ranks sharing an interface enumerate the interface
  !! sides in the same order (the global side id is shared), so the packed
  !! send buffer on one rank and the packed receive buffer on the other line
  !! up without any per-message metadata. The (element,side,flip) triplets
  !! are copied to the device for the pack/unpack kernels, and the packed
  !! send/receive device buffers are allocated here.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: e1,s1,e2,r2,n,nhalo
    integer :: rankId
    integer,allocatable :: elems(:),sides(:),flips(:),ranks(:)
    integer,allocatable,target :: halosides(:)
    integer(int64),allocatable :: keys(:)
    integer(c_size_t) :: worksize

    rankId = mesh%decomp%rankId

    nhalo = 0
    do e1 = 1,this%nElem
      do s1 = 1,6
        e2 = mesh%sideInfo(3,s1,e1)
        if(e2 > 0) then
          if(mesh%decomp%elemToRank(e2) /= rankId) then
            nhalo = nhalo+1
          endif
        endif
      enddo
    enddo

    this%halo_nsides = nhalo
    this%halo_nnbr = 0
    if(nhalo == 0) then
      this%halo_built = .true.
      return
    endif

    allocate(elems(1:nhalo),sides(1:nhalo),flips(1:nhalo), &
             ranks(1:nhalo),keys(1:nhalo))

    n = 0
    do e1 = 1,this%nElem
      do s1 = 1,6
        e2 = mesh%sideInfo(3,s1,e1)
        if(e2 > 0) then
          r2 = mesh%decomp%elemToRank(e2)
          if(r2 /= rankId) then
            n = n+1
            elems(n) = e1
            sides(n) = s1
            flips(n) = mesh%sideInfo(4,s1,e1)-10*(mesh%sideInfo(4,s1,e1)/10)
            ranks(n) = r2
            ! Composite sort key: neighbor rank (major), global side id (minor).
            keys(n) = int(r2,int64)*2147483648_int64+ &
                      int(abs(mesh%sideInfo(2,s1,e1)),int64)
          endif
        endif
      enddo
    enddo

    call HaloKeySort(keys,elems,sides,flips,ranks,nhalo)

    ! Count the neighboring ranks and record the segment offsets
    this%halo_nnbr = 1
    do n = 2,nhalo
      if(ranks(n) /= ranks(n-1)) this%halo_nnbr = this%halo_nnbr+1
    enddo

    allocate(this%halo_rank(1:this%halo_nnbr))
    allocate(this%halo_offset(1:this%halo_nnbr+1))

    this%halo_nnbr = 1
    this%halo_rank(1) = ranks(1)
    this%halo_offset(1) = 0
    do n = 2,nhalo
      if(ranks(n) /= ranks(n-1)) then
        this%halo_nnbr = this%halo_nnbr+1
        this%halo_rank(this%halo_nnbr) = ranks(n)
        this%halo_offset(this%halo_nnbr) = n-1
      endif
    enddo
    this%halo_offset(this%halo_nnbr+1) = nhalo

    ! Flatten the (element,side,flip) triplets with 0-based element and side
    ! indices for the device pack/unpack kernels
    allocate(halosides(1:3*nhalo))
    do n = 1,nhalo
      halosides(3*n-2) = elems(n)-1
      halosides(3*n-1) = sides(n)-1
      halosides(3*n) = flips(n)
    enddo

    call gpuCheck(hipMalloc(this%halo_sides_gpu,sizeof(halosides)))
    call gpuCheck(hipMemcpy(this%halo_sides_gpu,c_loc(halosides), &
                            sizeof(halosides),hipMemcpyHostToDevice))

    worksize = int(this%interp%N+1,c_size_t)*(this%interp%N+1)* &
               int(nhalo,c_size_t)*int(this%nvar,c_size_t)*prec
    call gpuCheck(hipMalloc(this%halo_sendbuf_gpu,worksize))
    call gpuCheck(hipMalloc(this%halo_recvbuf_gpu,worksize))

    deallocate(elems,sides,flips,ranks,keys,halosides)

    this%halo_built = .true.

  endsubroutine SetupHaloExchange_MappedScalar3D

  subroutine HaloKeySort(keys,elems,sides,flips,ranks,n)
  !! Heap sort of the halo side list by ascending 64-bit key, permuting the
  !! companion arrays alongside the keys. O(n log n), in place; the relative
  !! order of equal keys is irrelevant because global side ids are unique.
    implicit none
    integer,intent(in) :: n
    integer(int64),intent(inout) :: keys(1:n)
    integer,intent(inout) :: elems(1:n),sides(1:n),flips(1:n),ranks(1:n)
    ! Local
    integer :: i,last

    do i = n/2,1,-1
      call HaloSiftDown(keys,elems,sides,flips,ranks,i,n)
    enddo
    do last = n,2,-1
      call HaloSwap(keys,elems,sides,flips,ranks,1,last)
      call HaloSiftDown(keys,elems,sides,flips,ranks,1,last-1)
    enddo

  endsubroutine HaloKeySort

  subroutine HaloSiftDown(keys,elems,sides,flips,ranks,start,last)
    implicit none
    integer(int64),intent(inout) :: keys(:)
    integer,intent(inout) :: elems(:),sides(:),flips(:),ranks(:)
    integer,intent(in) :: start,last
    ! Local
    integer :: root,child

    root = start
    do while(2*root <= last)
      child = 2*root
      if(child < last) then
        if(keys(child) < keys(child+1)) child = child+1
      endif
      if(keys(root) < keys(child)) then
        call HaloSwap(keys,elems,sides,flips,ranks,root,child)
        root = child
      else
        return
      endif
    enddo

  endsubroutine HaloSiftDown

  subroutine HaloSwap(keys,elems,sides,flips,ranks,i,j)
    implicit none
    integer(int64),intent(inout) :: keys(:)
    integer,intent(inout) :: elems(:),sides(:),flips(:),ranks(:)
    integer,intent(in) :: i,j
    ! Local
    integer(int64) :: k
    integer :: t

    k = keys(i); keys(i) = keys(j); keys(j) = k
    t = elems(i); elems(i) = elems(j); elems(j) = t
    t = sides(i); sides(i) = sides(j); sides(j) = t
    t = flips(i); flips(i) = flips(j); flips(j) = t
    t = ranks(i); ranks(i) = ranks(j); ranks(j) = t

  endsubroutine HaloSwap

  subroutine MPIExchangeAsync_MappedScalar3D(this,mesh)
  !! Post the aggregated halo exchange: one MPI_Irecv/MPI_Isend pair per
  !! neighboring rank, carrying every (side,variable) boundary trace shared
  !! with that rank in a single packed device buffer. The pack kernel
  !! synchronizes the device before returning so the send buffer is complete
  !! before MPI_Isend is posted.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: n,npts,cnt,disp
    integer :: iError
    integer :: msgCount
    real(prec),pointer :: sendbuf(:)
    real(prec),pointer :: recvbuf(:)

    call HaloPack_3D_gpu(this%boundary_gpu,this%halo_sendbuf_gpu, &
                         this%halo_sides_gpu,this%interp%N,this%nvar, &
                         this%nelem,this%halo_nsides)

    npts = (this%interp%N+1)*(this%interp%N+1)*this%nvar
    call c_f_pointer(this%halo_sendbuf_gpu,sendbuf,[this%halo_nsides*npts])
    call c_f_pointer(this%halo_recvbuf_gpu,recvbuf,[this%halo_nsides*npts])

    msgCount = 0
    do n = 1,this%halo_nnbr

      cnt = (this%halo_offset(n+1)-this%halo_offset(n))*npts
      disp = this%halo_offset(n)*npts

      msgCount = msgCount+1
      call MPI_IRECV(recvbuf(disp+1),cnt, &
                     mesh%decomp%mpiPrec, &
                     this%halo_rank(n),0, &
                     mesh%decomp%mpiComm, &
                     mesh%decomp%requests(msgCount),iError)

      msgCount = msgCount+1
      call MPI_ISEND(sendbuf(disp+1),cnt, &
                     mesh%decomp%mpiPrec, &
                     this%halo_rank(n),0, &
                     mesh%decomp%mpiComm, &
                     mesh%decomp%requests(msgCount),iError)

    enddo

    mesh%decomp%msgCount = msgCount

  endsubroutine MPIExchangeAsync_MappedScalar3D

  subroutine SideExchange_MappedScalar3D(this,mesh)
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: offset

    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)

    if(mesh%decomp%mpiEnabled) then
      if(.not. this%halo_built) then
        call this%SetupHaloExchange(mesh)
      endif
      if(this%halo_nsides > 0) then
        call this%MPIExchangeAsync(mesh)
      endif
    endif

    ! The local (same-rank) side exchange runs on the device while the
    ! aggregated MPI messages are in flight.
    call SideExchange_3D_gpu(this%extboundary_gpu, &
                             this%boundary_gpu,mesh%sideinfo_gpu,mesh%decomp%elemToRank_gpu, &
                             mesh%decomp%rankid,offset,this%interp%N,this%nvar,this%nelem)

    if(mesh%decomp%mpiEnabled) then
      if(this%halo_nsides > 0) then
        call mesh%decomp%FinalizeMPIExchangeAsync()
        ! Unpack the received traces into extBoundary, applying side flips
        call HaloUnpack_3D_gpu(this%halo_recvbuf_gpu,this%extboundary_gpu, &
                               this%halo_sides_gpu,this%interp%N,this%nvar, &
                               this%nelem,this%halo_nsides)
      endif
    endif

  endsubroutine SideExchange_MappedScalar3D

  subroutine MappedGradient_MappedScalar3D(this,df)
  !! Calculates the gradient of a function using the strong form of the gradient
  !! in mapped coordinates.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(c_ptr),intent(out) :: df
    ! Local
    real(prec),pointer :: f_p(:,:,:,:,:,:)
    type(c_ptr) :: fc

    call ContravariantWeight_3D_gpu(this%interior_gpu, &
                                    this%geometry%dsdx%interior_gpu,this%jas_gpu, &
                                    this%interp%N,this%nvar,this%nelem)

    ! From Vector divergence
    call c_f_pointer(this%jas_gpu,f_p, &
                     [this%interp%N+1,this%interp%N+1,this%interp%N+1,this%nelem,3*this%nvar,3])

    fc = c_loc(f_p(1,1,1,1,1,1))
    call self_blas_matrixop_dim1_3d(this%interp%dMatrix_gpu,fc,df, &
                                    this%interp%N,this%interp%N,3*this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,1,2))
    call self_blas_matrixop_dim2_3d(this%interp%dMatrix_gpu,fc,df, &
                                    1.0_c_prec,this%interp%N,this%interp%N,3*this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,1,3))
    call self_blas_matrixop_dim3_3d(this%interp%dMatrix_gpu,fc,df, &
                                    1.0_c_prec,this%interp%N,this%interp%N,3*this%nvar,this%nelem,this%blas_handle)

    f_p => null()

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
    ! Local
    real(prec),pointer :: f_p(:,:,:,:,:,:)
    type(c_ptr) :: fc

    call ContravariantWeight_3D_gpu(this%interior_gpu, &
                                    this%geometry%dsdx%interior_gpu,this%jas_gpu, &
                                    this%interp%N,this%nvar,this%nelem)

    ! From Vector divergence
    call c_f_pointer(this%jas_gpu,f_p, &
                     [this%interp%N+1,this%interp%N+1,this%interp%N+1,this%nelem,3*this%nvar,3])

    fc = c_loc(f_p(1,1,1,1,1,1))
    call self_blas_matrixop_dim1_3d(this%interp%dgMatrix_gpu,fc,df, &
                                    this%interp%N,this%interp%N,3*this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,1,2))
    call self_blas_matrixop_dim2_3d(this%interp%dgMatrix_gpu,fc,df, &
                                    1.0_c_prec,this%interp%N,this%interp%N,3*this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,1,3))
    call self_blas_matrixop_dim3_3d(this%interp%dgMatrix_gpu,fc,df, &
                                    1.0_c_prec,this%interp%N,this%interp%N,3*this%nvar,this%nelem,this%blas_handle)

    f_p => null()

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
