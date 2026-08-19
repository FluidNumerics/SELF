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

module SELF_DGModel3D

  use SELF_DGModel3D_t
  use iso_fortran_env,only:int64
  use SELF_GPU
  use SELF_GPUInterfaces
  use SELF_BoundaryConditions
  use SELF_Geometry_3D
  use SELF_Mesh_3D
  use SELF_TransferPlan_3D

  implicit none

  type,extends(DGModel3D_t) :: DGModel3D
    !! Device-resident staging for the AMR solution transfer (the 3-D analogue of 2-D Stage 6a).
    !! These buffers persist across adaptation epochs and grow monotonically, so a settled run
    !! performs no allocation here at all; xferAllocBytes / planAllocElem record the current
    !! capacity.
    type(c_ptr) :: xferOld_gpu = c_null_ptr !! pre-regrid solution, staged device-side
    integer(c_size_t) :: xferAllocBytes = 0
    type(c_ptr) :: xferKind_gpu = c_null_ptr !! plan%sourceKind
    type(c_ptr) :: xferElem_gpu = c_null_ptr !! plan%sourceElem
    type(c_ptr) :: xferFamily_gpu = c_null_ptr !! plan%family
    type(c_ptr) :: xferDepth_gpu = c_null_ptr !! plan%depth
    type(c_ptr) :: xferPath_gpu = c_null_ptr !! plan%path
    integer :: planAllocElem = 0 !! new-element count the plan buffers are sized for
    integer :: planAllocStride = 0 !! plan%path leading dimension they are sized for
    integer :: xferNOld = 0 !! element count of the staged field (its device stride)
    !! Migrated old-element window, device-resident: the contiguous run of OLD elements this
    !! rank's new element range references, assembled by MigrateOldWindow and consumed by
    !! ApplyTransferPlan. Its own capacity and lifecycle, separate from xferOld_gpu, because the
    !! two have different lifetimes and different sizes - a window is the hull of a rank's new
    !! range in the OLD numbering and can be larger than the local old field - and because keeping
    !! them apart is what lets an unpaired stage/migrate be caught rather than silently served.
    type(c_ptr) :: xferWin_gpu = c_null_ptr
    integer(c_size_t) :: xferWinBytes = 0 !! current capacity, grown but never shrunk
    integer :: xferNWin = 0 !! window element count (its device stride); 0 means no window
    integer :: xferWinFirst = 0 !! global old index of the window's first element

  contains

    procedure :: Init => Init_DGModel3D
    procedure :: Free => Free_DGModel3D
    procedure :: Regrid => Regrid_DGModel3D

    procedure :: StageSolutionForTransfer => StageSolutionForTransfer_DGModel3D
    procedure :: ApplyTransferPlan => ApplyTransferPlan_DGModel3D
    procedure :: MigrateOldWindow => MigrateOldWindow_DGModel3D
    procedure :: DownloadOldWindow => DownloadOldWindow_DGModel3D

    procedure :: UpdateSolution => UpdateSolution_DGModel3D

    procedure :: CalculateEntropy => CalculateEntropy_DGModel3D
    procedure :: BoundaryFlux => BoundaryFlux_DGModel3D
    procedure :: FluxMethod => fluxmethod_DGModel3D
    procedure :: SourceMethod => sourcemethod_DGModel3D

    procedure :: UpdateGRK2 => UpdateGRK2_DGModel3D
    procedure :: UpdateGRK3 => UpdateGRK3_DGModel3D
    procedure :: UpdateGRK4 => UpdateGRK4_DGModel3D

    procedure :: CalculateSolutionGradient => CalculateSolutionGradient_DGModel3D
    procedure :: CalculateTendency => CalculateTendency_DGModel3D

  endtype DGModel3D

contains

  subroutine UpdateSolution_DGModel3D(this,dt)
    !! Computes a solution update as , where dt is either provided through the interface
    !! or taken as the Model's stored time step size (model % dt)
    implicit none
    class(DGModel3D),intent(inout) :: this
    real(prec),optional,intent(in) :: dt
    ! Local
    real(prec) :: dtLoc
    integer :: ndof

    if(present(dt)) then
      dtLoc = dt
    else
      dtLoc = this%dt
    endif
    ! The stepped variables occupy the leading nstepped entries of the
    ! variable dimension (the slowest-varying index of the interior array),
    ! so restricting ndof to nstepped updates exactly variables 1:nstepped.
    ndof = this%nstepped* &
           this%solution%nelem* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)

    ! Fused tendency + Euler update (dSdt = source - fluxDivergence formed in
    ! registers, no separate CalculateDSDt pass).
    call UpdateSolution_CalculateDSDt_gpu(this%solution%interior_gpu,this%fluxDivergence%interior_gpu, &
                                          this%source%interior_gpu,dtLoc,ndof)

  endsubroutine UpdateSolution_DGModel3D

  subroutine UpdateGRK2_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: ndof

    ! The stepped variables occupy the leading nstepped entries of the
    ! variable dimension (the slowest-varying index of the interior array),
    ! so restricting ndof to nstepped updates exactly variables 1:nstepped.
    ndof = this%nstepped* &
           this%solution%nelem* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)

    ! Fused tendency + RK stage update (no separate CalculateDSDt pass).
    call UpdateGRK_CalculateDSDt_gpu(this%worksol%interior_gpu,this%solution%interior_gpu, &
                                     this%fluxDivergence%interior_gpu,this%source%interior_gpu, &
                                     rk2_a(m),rk2_g(m),this%dt,ndof)

  endsubroutine UpdateGRK2_DGModel3D

  subroutine UpdateGRK3_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: ndof

    ! The stepped variables occupy the leading nstepped entries of the
    ! variable dimension (the slowest-varying index of the interior array),
    ! so restricting ndof to nstepped updates exactly variables 1:nstepped.
    ndof = this%nstepped* &
           this%solution%nelem* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)

    ! Fused tendency + RK stage update (no separate CalculateDSDt pass).
    call UpdateGRK_CalculateDSDt_gpu(this%worksol%interior_gpu,this%solution%interior_gpu, &
                                     this%fluxDivergence%interior_gpu,this%source%interior_gpu, &
                                     rk3_a(m),rk3_g(m),this%dt,ndof)

  endsubroutine UpdateGRK3_DGModel3D

  subroutine UpdateGRK4_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: ndof

    ! The stepped variables occupy the leading nstepped entries of the
    ! variable dimension (the slowest-varying index of the interior array),
    ! so restricting ndof to nstepped updates exactly variables 1:nstepped.
    ndof = this%nstepped* &
           this%solution%nelem* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)

    ! Fused tendency + RK stage update (no separate CalculateDSDt pass).
    call UpdateGRK_CalculateDSDt_gpu(this%worksol%interior_gpu,this%solution%interior_gpu, &
                                     this%fluxDivergence%interior_gpu,this%source%interior_gpu, &
                                     rk4_a(m),rk4_g(m),this%dt,ndof)

  endsubroutine UpdateGRK4_DGModel3D

  subroutine CalculateSolutionGradient_DGModel3D(this)
    implicit none
    class(DGModel3D),intent(inout) :: this

    call this%solution%AverageSides()

    call this%solution%MappedDGGradient(this%solutionGradient%interior_gpu)

    ! interpolate the solutiongradient to the element boundaries
    call this%solutionGradient%BoundaryInterp()

    ! perform the side exchange to populate the
    ! solutionGradient % extBoundary attribute
    call this%solutionGradient%SideExchange(this%mesh)

    ! populate the solutionGradient % extBoundary attribute on
    ! nonconforming (mortar) interfaces
    if(this%mesh%nMortars > 0) then
      call this%solutionGradient%MortarExchange(this%mesh)
    endif

  endsubroutine CalculateSolutionGradient_DGModel3D

  subroutine CalculateEntropy_DGModel3D(this)
    implicit none
    class(DGModel3D),intent(inout) :: this
    ! Local
    integer :: iel,i,j,k,ierror
    real(prec) :: e,jac
    real(prec) :: s(1:this%nvar)

    call gpuCheck(hipMemcpy(c_loc(this%solution%interior), &
                            this%solution%interior_gpu,sizeof(this%solution%interior), &
                            hipMemcpyDeviceToHost))

    e = 0.0_prec
    do iel = 1,this%geometry%nelem
      do k = 1,this%solution%interp%N+1
        do j = 1,this%solution%interp%N+1
          do i = 1,this%solution%interp%N+1
            jac = abs(this%geometry%J%interior(i,j,k,iel,1))
            s = this%solution%interior(i,j,k,iel,1:this%nvar)
            e = e+this%entropy_func(s)*jac* &
                this%solution%interp%qWeights(i)* &
                this%solution%interp%qWeights(j)* &
                this%solution%interp%qWeights(k)
          enddo
        enddo
      enddo
    enddo

    if(this%mesh%decomp%mpiEnabled) then
      call mpi_allreduce(e, &
                         this%entropy, &
                         1, &
                         this%mesh%decomp%mpiPrec, &
                         MPI_SUM, &
                         this%mesh%decomp%mpiComm, &
                         iError)
    else
      this%entropy = e
    endif

  endsubroutine CalculateEntropy_DGModel3D

  subroutine fluxmethod_DGModel3D(this)
    implicit none
    class(DGModel3D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: i,j,k
    real(prec) :: s(1:this%nvar),dsdx(1:this%nvar,1:3)

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%solution%N+1,iel=1:this%mesh%nElem)

      s = this%solution%interior(i,j,k,iel,1:this%nvar)
      dsdx = this%solutionGradient%interior(i,j,k,iel,1:this%nvar,1:3)
      this%flux%interior(i,j,k,iel,1:this%nvar,1:3) = this%flux3d(s,dsdx)

    enddo

    call gpuCheck(hipMemcpy(this%flux%interior_gpu, &
                            c_loc(this%flux%interior), &
                            sizeof(this%flux%interior), &
                            hipMemcpyHostToDevice))

  endsubroutine fluxmethod_DGModel3D

  subroutine BoundaryFlux_DGModel3D(this)
    ! this method uses an linear upwind solver for the
    ! advective flux and the bassi-rebay method for the
    ! diffusive fluxes
    implicit none
    class(DGModel3D),intent(inout) :: this
    ! Local
    integer :: i,j,k,iel
    real(prec) :: sL(1:this%nvar),sR(1:this%nvar)
    real(prec) :: dsdx(1:this%nvar,1:3)
    real(prec) :: nhat(1:3),nmag

    call gpuCheck(hipMemcpy(c_loc(this%solution%boundary), &
                            this%solution%boundary_gpu,sizeof(this%solution%boundary), &
                            hipMemcpyDeviceToHost))

    call gpuCheck(hipMemcpy(c_loc(this%solution%extboundary), &
                            this%solution%extboundary_gpu,sizeof(this%solution%extboundary), &
                            hipMemcpyDeviceToHost))

    call gpuCheck(hipMemcpy(c_loc(this%solutiongradient%avgboundary), &
                            this%solutiongradient%avgboundary_gpu,sizeof(this%solutiongradient%avgboundary), &
                            hipMemcpyDeviceToHost))

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:6,iel=1:this%mesh%nElem)
      ! Get the boundary normals on cell edges from the mesh geometry
      nhat = this%geometry%nHat%boundary(i,j,k,iEl,1,1:3)
      sL = this%solution%boundary(i,j,k,iel,1:this%nvar) ! interior solution
      sR = this%solution%extboundary(i,j,k,iel,1:this%nvar) ! exterior solution
      dsdx = this%solutiongradient%avgboundary(i,j,k,iel,1:this%nvar,1:3)
      nmag = this%geometry%nScale%boundary(i,j,k,iEl,1)

      this%flux%boundaryNormal(i,j,k,iEl,1:this%nvar) = this%riemannflux3d(sL,sR,dsdx,nhat)*nmag

    enddo

    call gpuCheck(hipMemcpy(this%flux%boundarynormal_gpu, &
                            c_loc(this%flux%boundarynormal), &
                            sizeof(this%flux%boundarynormal), &
                            hipMemcpyHostToDevice))

  endsubroutine BoundaryFlux_DGModel3D

  subroutine sourcemethod_DGModel3D(this)
    implicit none
    class(DGModel3D),intent(inout) :: this
    ! Local
    integer :: i,j,k,iel
    real(prec) :: s(1:this%nvar),dsdx(1:this%nvar,1:3)

    call gpuCheck(hipMemcpy(c_loc(this%solution%interior), &
                            this%solution%interior_gpu,sizeof(this%solution%interior), &
                            hipMemcpyDeviceToHost))

    call gpuCheck(hipMemcpy(c_loc(this%solutiongradient%interior), &
                            this%solutiongradient%interior_gpu,sizeof(this%solutiongradient%interior), &
                            hipMemcpyDeviceToHost))

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  k=1:this%solution%N+1,iel=1:this%mesh%nElem)

      s = this%solution%interior(i,j,k,iel,1:this%nvar)
      dsdx = this%solutionGradient%interior(i,j,k,iel,1:this%nvar,1:3)
      this%source%interior(i,j,k,iel,1:this%nvar) = this%source3d(s,dsdx)

    enddo

    call gpuCheck(hipMemcpy(this%source%interior_gpu, &
                            c_loc(this%source%interior), &
                            sizeof(this%source%interior), &
                            hipMemcpyHostToDevice))

  endsubroutine sourcemethod_DGModel3D

  subroutine Init_DGModel3D(this,mesh,geometry)
    !! Initialize the 3D DG model, then upload BC element/side arrays to GPU.
    implicit none
    class(DGModel3D),intent(out) :: this
    type(Mesh3D),intent(in),target :: mesh
    type(SEMHex),intent(in),target :: geometry
    ! Local
    type(BoundaryCondition),pointer :: bc

    call Init_DGModel3D_t(this,mesh,geometry)

    ! Upload hyperbolic BC element/side arrays to device
    bc => this%hyperbolicBCs%head
    do while(associated(bc))
      if(bc%nBoundaries > 0) then
        call gpuCheck(hipMalloc(bc%elements_gpu,sizeof(bc%elements)))
        call gpuCheck(hipMemcpy(bc%elements_gpu,c_loc(bc%elements), &
                                sizeof(bc%elements),hipMemcpyHostToDevice))
        call gpuCheck(hipMalloc(bc%sides_gpu,sizeof(bc%sides)))
        call gpuCheck(hipMemcpy(bc%sides_gpu,c_loc(bc%sides), &
                                sizeof(bc%sides),hipMemcpyHostToDevice))
      endif
      bc => bc%next
    enddo

    ! Upload parabolic BC element/side arrays to device
    bc => this%parabolicBCs%head
    do while(associated(bc))
      if(bc%nBoundaries > 0) then
        call gpuCheck(hipMalloc(bc%elements_gpu,sizeof(bc%elements)))
        call gpuCheck(hipMemcpy(bc%elements_gpu,c_loc(bc%elements), &
                                sizeof(bc%elements),hipMemcpyHostToDevice))
        call gpuCheck(hipMalloc(bc%sides_gpu,sizeof(bc%sides)))
        call gpuCheck(hipMemcpy(bc%sides_gpu,c_loc(bc%sides), &
                                sizeof(bc%sides),hipMemcpyHostToDevice))
      endif
      bc => bc%next
    enddo

  endsubroutine Init_DGModel3D

  subroutine Free_DGModel3D(this)
    !! Free the 3D DG model, including GPU BC arrays.
    implicit none
    class(DGModel3D),intent(inout) :: this
    ! Local
    type(BoundaryCondition),pointer :: bc

    ! Free hyperbolic BC device arrays
    bc => this%hyperbolicBCs%head
    do while(associated(bc))
      if(c_associated(bc%elements_gpu)) call gpuCheck(hipFree(bc%elements_gpu))
      if(c_associated(bc%sides_gpu)) call gpuCheck(hipFree(bc%sides_gpu))
      bc%elements_gpu = c_null_ptr
      bc%sides_gpu = c_null_ptr
      bc => bc%next
    enddo

    ! Free parabolic BC device arrays
    bc => this%parabolicBCs%head
    do while(associated(bc))
      if(c_associated(bc%elements_gpu)) call gpuCheck(hipFree(bc%elements_gpu))
      if(c_associated(bc%sides_gpu)) call gpuCheck(hipFree(bc%sides_gpu))
      bc%elements_gpu = c_null_ptr
      bc%sides_gpu = c_null_ptr
      bc => bc%next
    enddo

    ! Release the AMR transfer staging buffers. These are lazily created by
    ! StageSolutionForTransfer / ApplyTransferPlan and survive across adaptation epochs (and
    ! across Regrid, which runs between the two), so they are owned by the model and released
    ! only here.
    if(c_associated(this%xferOld_gpu)) call gpuCheck(hipFree(this%xferOld_gpu))
    if(c_associated(this%xferWin_gpu)) call gpuCheck(hipFree(this%xferWin_gpu))
    if(c_associated(this%xferKind_gpu)) call gpuCheck(hipFree(this%xferKind_gpu))
    if(c_associated(this%xferElem_gpu)) call gpuCheck(hipFree(this%xferElem_gpu))
    if(c_associated(this%xferFamily_gpu)) call gpuCheck(hipFree(this%xferFamily_gpu))
    if(c_associated(this%xferDepth_gpu)) call gpuCheck(hipFree(this%xferDepth_gpu))
    if(c_associated(this%xferPath_gpu)) call gpuCheck(hipFree(this%xferPath_gpu))
    this%xferOld_gpu = c_null_ptr
    this%xferWin_gpu = c_null_ptr
    this%xferKind_gpu = c_null_ptr
    this%xferElem_gpu = c_null_ptr
    this%xferFamily_gpu = c_null_ptr
    this%xferDepth_gpu = c_null_ptr
    this%xferPath_gpu = c_null_ptr
    this%xferAllocBytes = 0
    this%planAllocElem = 0
    this%planAllocStride = 0
    this%xferNOld = 0
    this%xferWinBytes = 0
    this%xferNWin = 0
    this%xferWinFirst = 0

    call Free_DGModel3D_t(this)

  endsubroutine Free_DGModel3D

  subroutine StageSolutionForTransfer_DGModel3D(this)
    !! Stage the pre-regrid solution on the DEVICE, replacing the base implementation's
    !! device-to-host copy of the whole field. Pair with ApplyTransferPlan, which consumes the
    !! staged copy:
    !!
    !!     call model%StageSolutionForTransfer()
    !!     call model%Regrid(newMesh,newGeom)
    !!     call model%ApplyTransferPlan(plan,interp,eFirst,eLast)
    !!
    !! The staging buffer is model-owned and grows monotonically, so a settled adapting run
    !! performs no allocation in this path; it is released in Free.
    implicit none
    class(DGModel3D),intent(inout) :: this
    ! Local
    integer(c_size_t) :: nbytes

    nbytes = int(this%solution%interp%N+1,c_size_t)*(this%solution%interp%N+1)* &
             (this%solution%interp%N+1)*this%solution%nElem*this%nvar*prec

    if(nbytes > this%xferAllocBytes) then
      if(c_associated(this%xferOld_gpu)) call gpuCheck(hipFree(this%xferOld_gpu))
      call gpuCheck(hipMalloc(this%xferOld_gpu,nbytes))
      this%xferAllocBytes = nbytes
    endif

    call gpuCheck(hipMemcpy(this%xferOld_gpu,this%solution%interior_gpu,nbytes, &
                            hipMemcpyDeviceToDevice))

    ! Remember the staged field's element count: it is the stride the transfer kernel must use
    ! to index the staged buffer, and it is no longer recoverable from the model after Regrid.
    this%xferNOld = this%solution%nElem

  endsubroutine StageSolutionForTransfer_DGModel3D

  subroutine ApplyTransferPlan_DGModel3D(this,plan,interp,eFirst,eLast,uGlobal,oldFirst)
    !! Apply the transfer plan on the device, writing solution%interior_gpu directly and moving
    !! no solution data across the PCIe/xGMI link - on any number of ranks.
    !!
    !! The old field is read from whichever device buffer the pairing routine left it in: the
    !! whole staged local field on a single rank (StageSolutionForTransfer), or this rank's
    !! migrated window on several (MigrateOldWindow), in which case the kernel is given the
    !! window's element stride and the global old index of its first element so it can rebase the
    !! plan's source indices. The two are mutually exclusive and one of them is required.
    !!
    !! uGlobal remains for a caller that supplies host-side old-field data - the Stage-5 v1
    !! allgather migration, and external callers - and takes the portable host path, forwarding
    !! oldFirst so the window is indexed correctly.
    implicit none
    ! target: matches the base declaration, which remaps a pointer onto its window buffer.
    class(DGModel3D),intent(inout),target :: this
    type(TransferPlan3D),intent(in),target :: plan
    type(Lagrange),intent(in) :: interp
    integer,intent(in) :: eFirst
    integer,intent(in) :: eLast
    real(prec),intent(in),optional,contiguous :: uGlobal(:,:,:,:,:)
    integer,intent(in),optional :: oldFirst
    ! Local
    integer :: nLocal,pathStride,li,c,src,o1,o2
    integer :: nOldStride,oldFirst0
    integer(c_size_t) :: nb
    type(c_ptr) :: uOld_gpu

    if(present(uGlobal)) then
      ! Caller-supplied host old field: the Stage-5 v1 allgather migration, and any external
      ! caller. Never combined with a migrated device window - that would be two sources for one
      ! apply, so it is a caller bug rather than something to resolve silently.
      if(this%xferNWin > 0) then
        print*,__FILE__,':',__LINE__, &
          ' : Error : ApplyTransferPlan given both uGlobal and a migrated window.'
        stop 1
      endif
      call ApplyTransferPlan_DGModel3D_t(this,plan,interp,eFirst,eLast,uGlobal,oldFirst)
      return
    endif

    ! Which device buffer holds the old field, and how the kernel indexes it.
    !
    ! The guards are on the element counts, not on c_associated(). Both buffers are persistent
    ! capacity buffers that stay allocated after use, so testing a pointer would only catch a
    ! missing stage before the FIRST epoch; from the second on, an unpaired call would silently
    ! transfer the previous epoch's field. xferNOld and xferNWin are set by
    ! StageSolutionForTransfer and MigrateOldWindow respectively and both cleared below, so they
    ! track the pairing the way the base implementation's allocatable transferStage does.
    if(this%xferNWin > 0 .and. this%xferNOld > 0) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : ApplyTransferPlan has both a staged local field and a migrated window.'
      stop 1
    elseif(this%xferNWin > 0) then
      ! Multi-rank: a migrated window of the old field, first element global old xferWinFirst.
      uOld_gpu = this%xferWin_gpu
      nOldStride = this%xferNWin
      oldFirst0 = this%xferWinFirst-1
    elseif(this%xferNOld > 0) then
      ! Single rank: the whole local old field, so the window base is zero and the kernel's
      ! indexing is exactly what it was before the window form existed.
      uOld_gpu = this%xferOld_gpu
      nOldStride = this%xferNOld
      oldFirst0 = 0
    else
      print*,__FILE__,':',__LINE__, &
        ' : Error : ApplyTransferPlan called without a staged solution or a migrated window.'
      stop 1
    endif

    ! The device kernel holds three (N+1)^3 working buffers in shared memory, sized to a compile
    ! time bound; guard the degree here rather than overrunning them.
    if(interp%N+1 > 12) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : the device solution transfer supports N+1 <= AMR3D_MAXNP (12).'
      stop 1
    endif

    nLocal = eLast-eFirst+1
    pathStride = size(plan%path,1)

    ! Window coverage. On the host path ApplyTransferPlanWindow checks every element's sources
    ! against the window and stops if one falls outside it. The kernel cannot: it has no way to
    ! report, and an out-of-window index there is a silent out-of-bounds device read rather than a
    ! wrong answer you can see. The same conditions are checked here instead - integer comparisons
    ! over this rank's own element range, once per adapting epoch.
    if(this%xferNWin > 0) then
      o1 = this%xferWinFirst
      o2 = o1+this%xferNWin-1
      do li = eFirst,eLast
        if(plan%sourceKind(li) == SELF_TRANSFER_RESTRICT) then
          do c = 1,8
            src = plan%family(c,li)
            if(src < o1 .or. src > o2) then
              print*,__FILE__,':',__LINE__, &
                ' : Error : a restricted child lies outside the migrated window.'
              stop 1
            endif
          enddo
        else
          src = plan%sourceElem(li)
          if(src < o1 .or. src > o2) then
            print*,__FILE__,':',__LINE__, &
              ' : Error : a transfer source lies outside the migrated window.'
            stop 1
          endif
        endif
      enddo
    endif

    ! (Re)size the device copies of the plan's integer arrays, reusing them across epochs.
    if(plan%nNew > this%planAllocElem .or. pathStride > this%planAllocStride) then
      if(c_associated(this%xferKind_gpu)) call gpuCheck(hipFree(this%xferKind_gpu))
      if(c_associated(this%xferElem_gpu)) call gpuCheck(hipFree(this%xferElem_gpu))
      if(c_associated(this%xferFamily_gpu)) call gpuCheck(hipFree(this%xferFamily_gpu))
      if(c_associated(this%xferDepth_gpu)) call gpuCheck(hipFree(this%xferDepth_gpu))
      if(c_associated(this%xferPath_gpu)) call gpuCheck(hipFree(this%xferPath_gpu))

      nb = int(plan%nNew,c_size_t)*c_int
      call gpuCheck(hipMalloc(this%xferKind_gpu,nb))
      call gpuCheck(hipMalloc(this%xferElem_gpu,nb))
      call gpuCheck(hipMalloc(this%xferDepth_gpu,nb))
      call gpuCheck(hipMalloc(this%xferFamily_gpu,8_c_size_t*nb))
      call gpuCheck(hipMalloc(this%xferPath_gpu,int(pathStride,c_size_t)*nb))

      ! All five are freed and re-allocated together, so both capacities describe the same set
      ! of buffers. They track the last plan rather than a high-water mark (as in the 2-D
      ! implementation this mirrors): a later epoch that grows either dimension re-allocates,
      ! and a settled mesh - the case that matters - trips neither condition.
      this%planAllocElem = plan%nNew
      this%planAllocStride = pathStride
    endif

    nb = int(plan%nNew,c_size_t)*c_int
    call gpuCheck(hipMemcpy(this%xferKind_gpu,c_loc(plan%sourceKind), &
                            nb,hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%xferElem_gpu,c_loc(plan%sourceElem), &
                            nb,hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%xferDepth_gpu,c_loc(plan%depth), &
                            nb,hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%xferFamily_gpu,c_loc(plan%family), &
                            8_c_size_t*nb,hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%xferPath_gpu,c_loc(plan%path), &
                            int(pathStride,c_size_t)*nb,hipMemcpyHostToDevice))

    call TransferSolution_3D_gpu(uOld_gpu,this%solution%interior_gpu, &
                                 this%xferKind_gpu,this%xferElem_gpu,this%xferFamily_gpu, &
                                 this%xferDepth_gpu,this%xferPath_gpu, &
                                 interp%mortarR_gpu,interp%mortarP_gpu, &
                                 pathStride,eFirst-1,oldFirst0,interp%N,this%nvar, &
                                 nOldStride,this%solution%nElem,nLocal)

    ! The staged field has been consumed. The buffer itself is kept (it is the capacity that
    ! makes a settled run allocation-free); only the pairing marker is cleared.
    this%xferNOld = 0
    this%xferNWin = 0

  endsubroutine ApplyTransferPlan_DGModel3D

  subroutine MigrateOldWindow_DGModel3D(this,winFirst,winLast,wFirst,wLast, &
                                        nBytesRecv,nBytesSent,nElemRemote)
    !! Assemble this rank's window of the pre-regrid old field IN DEVICE MEMORY, replacing the
    !! base implementation's host-memory migration. The rank's own run is copied device-to-device
    !! and the peers' runs are received directly into device memory, so no solution data crosses
    !! the host link and the windowed ApplyTransferPlan that follows runs as a device kernel.
    !!
    !! This hands device pointers to MPI, which requires a GPU-aware MPI. That is not a new
    !! requirement: the per-step aggregated halo exchange has always posted its sends and receives
    !! on device allocations (SideExchangeStart in SELF_MappedScalar_3D), with no host-staged
    !! fallback anywhere, and SELF_REQUIRE_GPU_AWARE_MPI makes the configure-time check fatal by
    !! default. A multi-rank GPU run that lacked it could never have taken a time step.
    !!
    !! Call it BEFORE Regrid, which releases the storage the sends read.
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: winFirst(:)
    integer,intent(in) :: winLast(:)
    integer,intent(in) :: wFirst
    integer,intent(in) :: wLast
    integer(int64),intent(inout) :: nBytesRecv
    integer(int64),intent(inout) :: nBytesSent
    integer(int64),intent(inout) :: nElemRemote
    ! Local
    integer :: perElem,nWinElem,nLocalOld,myFirst,a,b,msgCount
    integer(c_size_t) :: nbytes
    integer,allocatable :: requests(:)
    real(prec),pointer,contiguous :: win(:),uloc(:)

    perElem = (this%solution%interp%N+1)*(this%solution%interp%N+1)*(this%solution%interp%N+1)
    nWinElem = max(wLast-wFirst+1,0)
    nLocalOld = this%solution%nElem
    myFirst = this%mesh%decomp%offsetElem(this%mesh%decomp%rankId+1)+1

    ! Grow-only capacity, with a one-element floor so the descriptor below always has a target
    ! even for an empty window (a rank that owns no new elements still takes part, because its
    ! peers may need old elements it owns).
    nbytes = int(perElem,c_size_t)*int(max(nWinElem,1),c_size_t)*int(this%nvar,c_size_t)*prec
    if(nbytes > this%xferWinBytes) then
      if(c_associated(this%xferWin_gpu)) call gpuCheck(hipFree(this%xferWin_gpu))
      call gpuCheck(hipMalloc(this%xferWin_gpu,nbytes))
      this%xferWinBytes = nbytes
    endif

    ! Device allocations viewed as flat Fortran arrays, so the migration schedule can index them
    ! and MPI can be handed their addresses: the idiom SideExchangeStart already uses for the
    ! packed halo buffers.
    call c_f_pointer(this%xferWin_gpu,win, &
                     [perElem*max(nWinElem,1)*this%nvar])
    if(nLocalOld > 0) then
      call c_f_pointer(this%solution%interior_gpu,uloc,[perElem*nLocalOld*this%nvar])
    else
      ! This rank owned no old elements, so interior_gpu need not be a valid allocation. It posts
      ! no sends either (its old range intersects nobody's window), so uloc is never dereferenced;
      ! aim it somewhere valid rather than calling c_f_pointer on a possibly null pointer.
      uloc => win
    endif

    ! The run of my window that I already own, device-to-device. This ALSO synchronizes the
    ! device unconditionally - see MigrateWindowLocal_gpu - which is what makes it safe to post
    ! MPI against solution%interior_gpu on the next line. Unlike the host path, which overlaps
    ! its local copy with the messages, the copy happens first: the barrier has to precede the
    ! sends, and once per adapting epoch the lost overlap is not worth a second synchronization.
    call OwnedRun(this%mesh%decomp%offsetElem,this%mesh%decomp%rankId+1,wFirst,wLast,a,b)
    call MigrateWindowLocal_gpu(this%solution%interior_gpu,this%xferWin_gpu, &
                                perElem,this%nvar,max(nLocalOld,1),max(nWinElem,1), &
                                a-wFirst,a-myFirst,max(b-a+1,0))

    allocate(requests(1:2*this%nvar*this%mesh%decomp%nRanks))
    msgCount = 0
    call PostOldWindowExchange(this%mesh%decomp,perElem,this%nvar,nLocalOld,uloc, &
                               winFirst,winLast,wFirst,wLast,win,requests,msgCount, &
                               nBytesRecv,nBytesSent,nElemRemote)
    call FinishOldWindowExchange(requests,msgCount)
    deallocate(requests)

    win => null()
    uloc => null()
    this%xferWinFirst = wFirst
    this%xferNWin = nWinElem

  endsubroutine MigrateOldWindow_DGModel3D

  subroutine DownloadOldWindow_DGModel3D(this,wFirst,wLast,uWin)
    !! Copy the device-resident migrated window to the host, for the SELF_AMR_MIGRATE_VERIFY
    !! diagnostic only. This is exactly the host traffic the device path exists to avoid, so it
    !! must never appear on the default path.
    !!
    !! The comparison it feeds stays BITWISE. Everything that touched these values is data
    !! movement - a device-to-device copy, MPI byte transfers, and this download - so no
    !! arithmetic has been applied and exactness is available. The transfer APPLY is a different
    !! matter: it agrees between host and device only to round-off, because the device compiler
    !! contracts its multiply-accumulates into FMAs. The two must not be conflated.
    implicit none
    class(DGModel3D),intent(in) :: this
    integer,intent(in) :: wFirst
    integer,intent(in) :: wLast
    ! target: c_loc below requires it, and contiguous keeps the download a single memcpy.
    real(prec),intent(out),target,contiguous :: uWin(:,:,:,:,:)
    ! Local
    integer :: perElem,nWinElem
    integer(c_size_t) :: nbytes

    nWinElem = max(wLast-wFirst+1,0)
    if(nWinElem == 0) return

    perElem = (this%solution%interp%N+1)*(this%solution%interp%N+1)*(this%solution%interp%N+1)
    if(size(uWin,4) /= nWinElem .or. size(uWin,5) /= this%nvar) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : DownloadOldWindow given a buffer that does not match the window.'
      stop 1
    endif
    if(this%xferNWin /= nWinElem .or. this%xferWinFirst /= wFirst) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : DownloadOldWindow called without a matching migrated window.'
      stop 1
    endif

    nbytes = int(perElem,c_size_t)*int(nWinElem,c_size_t)*int(this%nvar,c_size_t)*prec
    call gpuCheck(hipMemcpy(c_loc(uWin),this%xferWin_gpu,nbytes,hipMemcpyDeviceToHost))

  endsubroutine DownloadOldWindow_DGModel3D

  subroutine Regrid_DGModel3D(this,mesh,geometry)
    !! GPU regrid: release the device copies of the old boundary-condition element/side
    !! arrays, rebuild the model storage and BC maps on the new mesh (base Regrid), then
    !! upload the new BC arrays to the device - the same device bookkeeping Init/Free do
    !! around the base Init/Free.
    implicit none
    class(DGModel3D),intent(inout) :: this
    type(Mesh3D),intent(in),target :: mesh
    type(SEMHex),intent(in),target :: geometry
    ! Local
    type(BoundaryCondition),pointer :: bc

    ! Free hyperbolic BC device arrays (old mesh)
    bc => this%hyperbolicBCs%head
    do while(associated(bc))
      if(c_associated(bc%elements_gpu)) call gpuCheck(hipFree(bc%elements_gpu))
      if(c_associated(bc%sides_gpu)) call gpuCheck(hipFree(bc%sides_gpu))
      bc%elements_gpu = c_null_ptr
      bc%sides_gpu = c_null_ptr
      bc => bc%next
    enddo

    ! Free parabolic BC device arrays (old mesh)
    bc => this%parabolicBCs%head
    do while(associated(bc))
      if(c_associated(bc%elements_gpu)) call gpuCheck(hipFree(bc%elements_gpu))
      if(c_associated(bc%sides_gpu)) call gpuCheck(hipFree(bc%sides_gpu))
      bc%elements_gpu = c_null_ptr
      bc%sides_gpu = c_null_ptr
      bc => bc%next
    enddo

    call Regrid_DGModel3D_t(this,mesh,geometry)

    ! Upload hyperbolic BC element/side arrays to device (new mesh)
    bc => this%hyperbolicBCs%head
    do while(associated(bc))
      if(bc%nBoundaries > 0) then
        call gpuCheck(hipMalloc(bc%elements_gpu,sizeof(bc%elements)))
        call gpuCheck(hipMemcpy(bc%elements_gpu,c_loc(bc%elements), &
                                sizeof(bc%elements),hipMemcpyHostToDevice))
        call gpuCheck(hipMalloc(bc%sides_gpu,sizeof(bc%sides)))
        call gpuCheck(hipMemcpy(bc%sides_gpu,c_loc(bc%sides), &
                                sizeof(bc%sides),hipMemcpyHostToDevice))
      endif
      bc => bc%next
    enddo

    ! Upload parabolic BC element/side arrays to device (new mesh)
    bc => this%parabolicBCs%head
    do while(associated(bc))
      if(bc%nBoundaries > 0) then
        call gpuCheck(hipMalloc(bc%elements_gpu,sizeof(bc%elements)))
        call gpuCheck(hipMemcpy(bc%elements_gpu,c_loc(bc%elements), &
                                sizeof(bc%elements),hipMemcpyHostToDevice))
        call gpuCheck(hipMalloc(bc%sides_gpu,sizeof(bc%sides)))
        call gpuCheck(hipMemcpy(bc%sides_gpu,c_loc(bc%sides), &
                                sizeof(bc%sides),hipMemcpyHostToDevice))
      endif
      bc => bc%next
    enddo

  endsubroutine Regrid_DGModel3D

  subroutine CalculateTendency_DGModel3D(this)
    implicit none
    class(DGModel3D),intent(inout) :: this

    call this%solution%BoundaryInterp()

    ! Post the halo exchange for the prognostic variables; static variables
    ! (nstepped+1:nvar) are carried in full by the first exchange only, and
    ! their boundary traces do not change thereafter. The MPI messages are
    ! in flight while the hooks, boundary conditions, source, and interior
    ! flux evaluations below execute.
    call this%solution%SideExchangeStart(this%mesh,this%nstepped)

    ! The hooks and methods between SideExchangeStart and SideExchangeFinish
    ! must not read extBoundary on interior (rank-shared) faces.
    ! SetBoundaryCondition writes extBoundary only on physical-boundary
    ! faces, which are disjoint from the faces the exchange fills.
    call this%PreTendencyHook() ! User-supplied
    call this%SetBoundaryCondition() ! User-supplied

    if(this%gradient_enabled) then
      ! The BR gradient consumes extBoundary (through the side averages), so
      ! the exchange must complete before the gradient is computed. The mortar
      ! exchange runs after SideExchangeFinish : it posts its own messages on
      ! mesh%decomp%requests and fills extBoundary on nonconforming faces,
      ! which the side averages also consume.
      call this%solution%SideExchangeFinish(this%mesh)
      if(this%mesh%nMortars > 0) then
        call this%solution%MortarExchange(this%mesh)
      endif
      call this%CalculateSolutionGradient()
      call this%SetGradientBoundaryCondition() ! User-supplied
      call this%solutionGradient%AverageSides()
      call this%SourceMethod() ! User supplied
      call this%BoundaryFlux() ! User supplied
      call this%FluxMethod() ! User supplied
    else
      ! Interior work overlaps with the halo exchange; BoundaryFlux is the
      ! first consumer of extBoundary and runs after the exchange completes.
      ! FluxMethod and BoundaryFlux write disjoint outputs (flux interior
      ! vs. boundarynormal), so this reordering does not change any
      ! floating-point results. The mortar exchange runs after
      ! SideExchangeFinish for the same reason.
      call this%SourceMethod() ! User supplied
      call this%FluxMethod() ! User supplied
      call this%solution%SideExchangeFinish(this%mesh)
      if(this%mesh%nMortars > 0) then
        call this%solution%MortarExchange(this%mesh)
      endif
      call this%BoundaryFlux() ! User supplied
    endif

    ! On mortar interfaces, replace the big face's surface-flux integrand with the
    ! projection of the small faces' integrands so that the interface is conservative
    if(this%mesh%nMortars > 0) then
      call this%flux%MortarFluxCollect(this%mesh)
    endif

    call this%flux%MappedDGDivergence(this%fluxDivergence%interior_gpu)

    ! The tendency (source - fluxDivergence) is formed inside the fused RK/Euler
    ! update kernels (UpdateGRK_CalculateDSDt_gpu / UpdateSolution_CalculateDSDt_gpu),
    ! so there is no separate CalculateDSDt pass here.

  endsubroutine CalculateTendency_DGModel3D

endmodule SELF_DGModel3D
