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
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_ESAtmo2D

  use SELF_ESAtmo2D_t
  use SELF_ECDGModel2D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use SELF_BoundaryConditions
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use iso_c_binding

  implicit none

  type,extends(ESAtmo2D_t),public :: ESAtmo2D

  contains

    procedure :: Init => Init_ESAtmo2D
    procedure :: Free => Free_ESAtmo2D
    procedure :: AdditionalInit => AdditionalInit_ESAtmo2D
    procedure :: BoundaryFlux => BoundaryFlux_ESAtmo2D
    procedure :: TwoPointFluxMethod => TwoPointFluxMethod_ESAtmo2D
    procedure :: SourceMethod => SourceMethod_ESAtmo2D
    procedure :: DiffusiveFluxMethod => DiffusiveFluxMethod_ESAtmo2D
    procedure :: DiffusiveBoundaryFlux => DiffusiveBoundaryFlux_ESAtmo2D
    procedure :: CalculateTendency => CalculateTendency_ESAtmo2D

  endtype ESAtmo2D

  interface
    subroutine hbc2d_nonormalflow_esatmo2d_gpu(extboundary,boundary, &
                                               nhat,elements,sides, &
                                               nBoundaries,N,nel,nvar) &
      bind(c,name="hbc2d_nonormalflow_esatmo2d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,nhat,elements,sides
      integer(c_int),value :: nBoundaries,N,nel,nvar
    endsubroutine hbc2d_nonormalflow_esatmo2d_gpu
  endinterface

  interface
    subroutine pbc2d_nostress_esatmo2d_gpu(extgrad,grad, &
                                           nhat,elements,sides, &
                                           nBoundaries,N,nel,nvar) &
      bind(c,name="pbc2d_nostress_esatmo2d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extgrad,grad,nhat,elements,sides
      integer(c_int),value :: nBoundaries,N,nel,nvar
    endsubroutine pbc2d_nostress_esatmo2d_gpu
  endinterface

  interface
    subroutine boundaryflux_esatmo2d_gpu(fb,fextb,nhat,nscale,flux, &
                                         p0,Rd,gamma,N,nel) &
      bind(c,name="boundaryflux_esatmo2d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fb,fextb,nhat,nscale,flux
      real(c_prec),value :: p0,Rd,gamma
      integer(c_int),value :: N,nel
    endsubroutine boundaryflux_esatmo2d_gpu
  endinterface

  interface
    subroutine twopointfluxmethod_esatmo2d_gpu(f,s,dsdx, &
                                               p0,Rd,gamma, &
                                               N,nvar,nel) &
      bind(c,name="twopointfluxmethod_esatmo2d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: f,s,dsdx
      real(c_prec),value :: p0,Rd,gamma
      integer(c_int),value :: N,nvar,nel
    endsubroutine twopointfluxmethod_esatmo2d_gpu
  endinterface

  interface
    subroutine diffusiveflux_esatmo2d_gpu(diffFlux,grad,nu,kappa, &
                                          N,nvar,nel) &
      bind(c,name="diffusiveflux_esatmo2d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: diffFlux,grad
      real(c_prec),value :: nu,kappa
      integer(c_int),value :: N,nvar,nel
    endsubroutine diffusiveflux_esatmo2d_gpu
  endinterface

  interface
    subroutine diffusiveboundaryflux_esatmo2d_gpu(fluxN,avgGrad,uBnd,uExt, &
                                                  nhat,nscale, &
                                                  nu,kappa,tau_nu,tau_kappa, &
                                                  N,nvar,nel) &
      bind(c,name="diffusiveboundaryflux_esatmo2d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fluxN,avgGrad,uBnd,uExt,nhat,nscale
      real(c_prec),value :: nu,kappa,tau_nu,tau_kappa
      integer(c_int),value :: N,nvar,nel
    endsubroutine diffusiveboundaryflux_esatmo2d_gpu
  endinterface

  interface
    subroutine sourcemethod_esatmo2d_gpu(source,solution,dsdx,J,dSplit, &
                                         N,nvar,nel) &
      bind(c,name="sourcemethod_esatmo2d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: source,solution,dsdx,J,dSplit
      integer(c_int),value :: N,nvar,nel
    endsubroutine sourcemethod_esatmo2d_gpu
  endinterface

contains

  subroutine Init_ESAtmo2D(this,mesh,geometry)
    implicit none
    class(ESAtmo2D),intent(out) :: this
    type(Mesh2D),intent(in),target :: mesh
    type(SEMQuad),intent(in),target :: geometry
    ! Local
    type(BoundaryCondition),pointer :: bc

    call Init_ECDGModel2D_t(this,mesh,geometry)

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

    ! Upload parabolic BC element/side arrays to device. Required for the
    ! GPU parabolic BC kernels (e.g. pbc2d_NoStress_ESAtmo2D_GPU_wrapper).
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

  endsubroutine Init_ESAtmo2D

  subroutine Free_ESAtmo2D(this)
    implicit none
    class(ESAtmo2D),intent(inout) :: this
    ! Local
    type(BoundaryCondition),pointer :: bc

    bc => this%hyperbolicBCs%head
    do while(associated(bc))
      if(c_associated(bc%elements_gpu)) call gpuCheck(hipFree(bc%elements_gpu))
      if(c_associated(bc%sides_gpu)) call gpuCheck(hipFree(bc%sides_gpu))
      bc%elements_gpu = c_null_ptr
      bc%sides_gpu = c_null_ptr
      bc => bc%next
    enddo

    bc => this%parabolicBCs%head
    do while(associated(bc))
      if(c_associated(bc%elements_gpu)) call gpuCheck(hipFree(bc%elements_gpu))
      if(c_associated(bc%sides_gpu)) call gpuCheck(hipFree(bc%sides_gpu))
      bc%elements_gpu = c_null_ptr
      bc%sides_gpu = c_null_ptr
      bc => bc%next
    enddo

    call Free_ECDGModel2D_t(this)

  endsubroutine Free_ESAtmo2D

  subroutine AdditionalInit_ESAtmo2D(this)
    implicit none
    class(ESAtmo2D),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    ! Call parent _t AdditionalInit (registers CPU BC, allocates diff buffers)
    call AdditionalInit_ESAtmo2D_t(this)

    ! Re-register with GPU-accelerated versions for both lists.
    bcfunc => hbc2d_NoNormalFlow_ESAtmo2D_GPU_wrapper
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

    bcfunc => pbc2d_NoStress_ESAtmo2D_GPU_wrapper
    call this%parabolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

  endsubroutine AdditionalInit_ESAtmo2D

  subroutine hbc2d_NoNormalFlow_ESAtmo2D_GPU_wrapper(bc,mymodel)
    !! GPU-accelerated no-normal-flow BC for 2D Entropy-Stable Atmosphere.
    !! Reflects normal momentum, mirrors density, rho*theta, and Phi.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(ESAtmo2D)
      if(bc%nBoundaries > 0) then
        call hbc2d_nonormalflow_esatmo2d_gpu( &
          m%solution%extBoundary_gpu, &
          m%solution%boundary_gpu, &
          m%geometry%nhat%boundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N, &
          m%solution%nElem,m%solution%nvar)
      endif
    endselect

  endsubroutine hbc2d_NoNormalFlow_ESAtmo2D_GPU_wrapper

  subroutine pbc2d_NoStress_ESAtmo2D_GPU_wrapper(bc,mymodel)
    !! GPU-accelerated parabolic no-stress / no-heat-flux BC for 2D Entropy-Stable Atmosphere.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(ESAtmo2D)
      if(bc%nBoundaries > 0) then
        call pbc2d_nostress_esatmo2d_gpu( &
          m%solutionGradient%extBoundary_gpu, &
          m%solutionGradient%boundary_gpu, &
          m%geometry%nhat%boundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N, &
          m%solution%nElem,m%solution%nvar)
      endif
    endselect

  endsubroutine pbc2d_NoStress_ESAtmo2D_GPU_wrapper

  subroutine BoundaryFlux_ESAtmo2D(this)
    !! LMARS interface flux on GPU.
    implicit none
    class(ESAtmo2D),intent(inout) :: this

    call boundaryflux_esatmo2d_gpu( &
      this%solution%boundary_gpu, &
      this%solution%extboundary_gpu, &
      this%geometry%nhat%boundary_gpu, &
      this%geometry%nscale%boundary_gpu, &
      this%flux%boundarynormal_gpu, &
      this%p0,this%Rd,this%cp/this%cv, &
      this%solution%interp%N, &
      this%solution%nelem)

  endsubroutine BoundaryFlux_ESAtmo2D

  subroutine TwoPointFluxMethod_ESAtmo2D(this)
    !! Souza et al. (2023) entropy-conservative two-point flux on GPU.
    implicit none
    class(ESAtmo2D),intent(inout) :: this

    call twopointfluxmethod_esatmo2d_gpu( &
      this%twoPointFlux%interior_gpu, &
      this%solution%interior_gpu, &
      this%geometry%dsdx%interior_gpu, &
      this%p0,this%Rd,this%cp/this%cv, &
      this%solution%interp%N, &
      this%solution%nvar, &
      this%solution%nelem)

  endsubroutine TwoPointFluxMethod_ESAtmo2D

  subroutine SourceMethod_ESAtmo2D(this)
    !! Souza et al. (2023) non-conservative gravity flux differencing on
    !! GPU. The geopotential lives at solution(:,:,:,5); the source for
    !! rho*v is computed via the SBP-EC two-point form using log-mean
    !! density and the contravariant metric.
    implicit none
    class(ESAtmo2D),intent(inout) :: this

    call sourcemethod_esatmo2d_gpu( &
      this%source%interior_gpu, &
      this%solution%interior_gpu, &
      this%geometry%dsdx%interior_gpu, &
      this%geometry%J%interior_gpu, &
      this%solution%interp%dSplitMatrix_gpu, &
      this%solution%interp%N, &
      this%solution%nvar, &
      this%solution%nelem)

  endsubroutine SourceMethod_ESAtmo2D

  subroutine DiffusiveFluxMethod_ESAtmo2D(this)
    !! GPU-resident fill of diffFlux%interior with the constant-coefficient
    !! Laplacian flux F_d(iVar) = -coeff(iVar) * d(s_iVar)/dx_d.
    implicit none
    class(ESAtmo2D),intent(inout) :: this

    call diffusiveflux_esatmo2d_gpu( &
      this%diffFlux%interior_gpu, &
      this%solutionGradient%interior_gpu, &
      this%nu,this%kappa, &
      this%solution%interp%N, &
      this%solution%nvar, &
      this%solution%nelem)

  endsubroutine DiffusiveFluxMethod_ESAtmo2D

  subroutine DiffusiveBoundaryFlux_ESAtmo2D(this)
    !! GPU-resident fill of diffFlux%boundaryNormal with the SIPG-stabilised
    !! BR1 diffusive flux:
    !!   f = -coeff*(avg_grad . n)*nmag + tau*(uL - uR)*nmag
    implicit none
    class(ESAtmo2D),intent(inout) :: this
    ! Local
    real(prec) :: np2,tau_nu,tau_kappa

    np2 = real((this%solution%interp%N+1)**2,prec)
    tau_nu = this%eta_penalty*this%nu*np2/this%length_scale
    tau_kappa = this%eta_penalty*this%kappa*np2/this%length_scale

    call diffusiveboundaryflux_esatmo2d_gpu( &
      this%diffFlux%boundarynormal_gpu, &
      this%solutionGradient%avgBoundary_gpu, &
      this%solution%boundary_gpu, &
      this%solution%extBoundary_gpu, &
      this%geometry%nhat%boundary_gpu, &
      this%geometry%nscale%boundary_gpu, &
      this%nu,this%kappa,tau_nu,tau_kappa, &
      this%solution%interp%N, &
      this%solution%nvar, &
      this%solution%nelem)

  endsubroutine DiffusiveBoundaryFlux_ESAtmo2D

  subroutine CalculateTendency_ESAtmo2D(this)
    !! GPU-resident tendency for ESAtmo2D.
    implicit none
    class(ESAtmo2D),intent(inout) :: this
    ! Local
    integer :: ndof,ndof_diff

    call this%solution%BoundaryInterp()

    ! Post the halo exchange for the prognostic variables; the MPI messages
    ! are in flight while the hooks, boundary conditions, source, and the
    ! two-point volume flux and its divergence below execute.
    call this%solution%SideExchangeStart(this%mesh,this%nstepped)

    call this%PreTendencyHook()
    call this%SetBoundaryCondition()

    if(this%gradient_enabled) then
      ! The BR gradient consumes extBoundary (through the side averages), so
      ! the exchange must complete before the gradient is computed.
      call this%solution%SideExchangeFinish(this%mesh)
      call this%CalculateSolutionGradient()
      call this%SetGradientBoundaryCondition()
      call this%solutionGradient%AverageSides()
    endif

    call this%SourceMethod()
    call this%TwoPointFluxMethod()

    call this%twoPointFlux%MappedDivergence(this%fluxDivergence%interior_gpu)

    ! BoundaryFlux is the first consumer of extBoundary; the two-point volume
    ! flux and its divergence above overlap with the halo exchange (a no-op
    ! wait when gradients already finished it). BoundaryFlux writes only the
    ! Riemann trace (flux boundarynormal), disjoint from the volume-term
    ! outputs, so this reordering does not change any floating-point results.
    call this%solution%SideExchangeFinish(this%mesh)
    call this%BoundaryFlux()

    call ECDGSurfaceContribution_2D_gpu( &
      this%flux%boundarynormal_gpu, &
      this%geometry%J%interior_gpu, &
      this%solution%interp%bMatrix_gpu, &
      this%solution%interp%qWeights_gpu, &
      this%fluxDivergence%interior_gpu, &
      this%solution%interp%N,this%solution%nVar,this%mesh%nElem)

    if(this%nu > 0.0_prec .or. this%kappa > 0.0_prec) then
      call this%DiffusiveFluxMethod()
      call this%DiffusiveBoundaryFlux()
      call this%diffFlux%MappedDGDivergence(this%diffDiv%interior_gpu)
      ndof_diff = this%solution%nVar* &
                  this%mesh%nElem* &
                  (this%solution%interp%N+1)* &
                  (this%solution%interp%N+1)
      call AccumulateField_gpu(this%fluxDivergence%interior_gpu, &
                               this%diffDiv%interior_gpu, &
                               ndof_diff)
    endif

    ndof = this%solution%nVar* &
           this%mesh%nElem* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)

    call CalculateDSDt_gpu(this%fluxDivergence%interior_gpu, &
                           this%source%interior_gpu, &
                           this%dSdt%interior_gpu,ndof)

  endsubroutine CalculateTendency_ESAtmo2D

endmodule SELF_ESAtmo2D
