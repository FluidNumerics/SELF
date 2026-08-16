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

  contains

    procedure :: Init => Init_DGModel3D
    procedure :: Free => Free_DGModel3D
    procedure :: Regrid => Regrid_DGModel3D

    procedure :: StageSolutionForTransfer => StageSolutionForTransfer_DGModel3D
    procedure :: ApplyTransferPlan => ApplyTransferPlan_DGModel3D

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
    if(c_associated(this%xferKind_gpu)) call gpuCheck(hipFree(this%xferKind_gpu))
    if(c_associated(this%xferElem_gpu)) call gpuCheck(hipFree(this%xferElem_gpu))
    if(c_associated(this%xferFamily_gpu)) call gpuCheck(hipFree(this%xferFamily_gpu))
    if(c_associated(this%xferDepth_gpu)) call gpuCheck(hipFree(this%xferDepth_gpu))
    if(c_associated(this%xferPath_gpu)) call gpuCheck(hipFree(this%xferPath_gpu))
    this%xferOld_gpu = c_null_ptr
    this%xferKind_gpu = c_null_ptr
    this%xferElem_gpu = c_null_ptr
    this%xferFamily_gpu = c_null_ptr
    this%xferDepth_gpu = c_null_ptr
    this%xferPath_gpu = c_null_ptr
    this%xferAllocBytes = 0
    this%planAllocElem = 0
    this%planAllocStride = 0
    this%xferNOld = 0

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

  subroutine ApplyTransferPlan_DGModel3D(this,plan,interp,eFirst,eLast,uGlobal)
    !! Apply the transfer plan on the device, writing solution%interior_gpu directly and moving
    !! no solution data across the PCIe/xGMI link.
    !!
    !! uGlobal is the multi-rank escape hatch: when the caller has assembled a global old field
    !! on the host (the Stage-5 v1 allgather migration), this falls back to the portable host
    !! path. A device transfer only pays off on several ranks together with a point-to-point
    !! (Stage-5 v2) migration; on a single rank it stands alone, which is the case optimized
    !! here.
    implicit none
    class(DGModel3D),intent(inout) :: this
    type(TransferPlan3D),intent(in),target :: plan
    type(Lagrange),intent(in) :: interp
    integer,intent(in) :: eFirst
    integer,intent(in) :: eLast
    real(prec),intent(in),optional :: uGlobal(:,:,:,:,:)
    ! Local
    integer :: nLocal,pathStride
    integer(c_size_t) :: nb

    if(present(uGlobal)) then
      call ApplyTransferPlan_DGModel3D_t(this,plan,interp,eFirst,eLast,uGlobal)
      return
    endif

    ! The guard is on xferNOld, not on c_associated(xferOld_gpu). The staging buffer is a
    ! persistent capacity buffer and stays allocated after use, so testing the pointer would only
    ! catch a missing stage before the FIRST epoch; from the second on, an unpaired call would
    ! silently transfer the previous epoch's field. xferNOld is set by StageSolutionForTransfer
    ! and cleared below, so it tracks the stage/apply pairing the way the base implementation's
    ! allocatable transferStage does.
    if(this%xferNOld <= 0) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : ApplyTransferPlan called without a staged solution.'
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

    call TransferSolution_3D_gpu(this%xferOld_gpu,this%solution%interior_gpu, &
                                 this%xferKind_gpu,this%xferElem_gpu,this%xferFamily_gpu, &
                                 this%xferDepth_gpu,this%xferPath_gpu, &
                                 interp%mortarR_gpu,interp%mortarP_gpu, &
                                 pathStride,eFirst-1,interp%N,this%nvar, &
                                 this%xferNOld,this%solution%nElem,nLocal)

    ! The staged field has been consumed. The buffer itself is kept (it is the capacity that
    ! makes a settled run allocation-free); only the pairing marker is cleared.
    this%xferNOld = 0

  endsubroutine ApplyTransferPlan_DGModel3D

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
