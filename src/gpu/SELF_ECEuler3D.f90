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

module SELF_ECEuler3D

  use SELF_ECEuler3D_t
  use SELF_ECDGModel3D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use SELF_BoundaryConditions
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use iso_c_binding

  implicit none

  type,extends(ECEuler3D_t),public :: ECEuler3D

  contains

    procedure :: Init => Init_ECEuler3D
    procedure :: Free => Free_ECEuler3D
    procedure :: AdditionalInit => AdditionalInit_ECEuler3D
    procedure :: BoundaryFlux => BoundaryFlux_ECEuler3D
    procedure :: TwoPointFluxMethod => TwoPointFluxMethod_ECEuler3D
    procedure :: SourceMethod => SourceMethod_ECEuler3D
    procedure :: DiffusiveFluxMethod => DiffusiveFluxMethod_ECEuler3D
    procedure :: DiffusiveBoundaryFlux => DiffusiveBoundaryFlux_ECEuler3D
    procedure :: CalculateTendency => CalculateTendency_ECEuler3D

  endtype ECEuler3D

  interface
    subroutine hbc3d_nonormalflow_eceuler3d_gpu(extboundary,boundary, &
                                                nhat,elements,sides, &
                                                nBoundaries,N,nel,nvar) &
      bind(c,name="hbc3d_nonormalflow_eceuler3d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,nhat,elements,sides
      integer(c_int),value :: nBoundaries,N,nel,nvar
    endsubroutine hbc3d_nonormalflow_eceuler3d_gpu
  endinterface

  interface
    subroutine pbc3d_nostress_eceuler3d_gpu(extgrad,grad, &
                                            nhat,elements,sides, &
                                            nBoundaries,N,nel,nvar) &
      bind(c,name="pbc3d_nostress_eceuler3d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extgrad,grad,nhat,elements,sides
      integer(c_int),value :: nBoundaries,N,nel,nvar
    endsubroutine pbc3d_nostress_eceuler3d_gpu
  endinterface

  interface
    subroutine boundaryflux_eceuler3d_gpu(fb,fextb,nhat,nscale,flux, &
                                          p0,Rd,gamma,N,nel) &
      bind(c,name="boundaryflux_eceuler3d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fb,fextb,nhat,nscale,flux
      real(c_prec),value :: p0,Rd,gamma
      integer(c_int),value :: N,nel
    endsubroutine boundaryflux_eceuler3d_gpu
  endinterface

  interface
    subroutine twopointfluxmethod_eceuler3d_gpu(f,s,dsdx, &
                                                p0,Rd,gamma, &
                                                N,nvar,nel) &
      bind(c,name="twopointfluxmethod_eceuler3d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: f,s,dsdx
      real(c_prec),value :: p0,Rd,gamma
      integer(c_int),value :: N,nvar,nel
    endsubroutine twopointfluxmethod_eceuler3d_gpu
  endinterface

  interface
    subroutine diffusiveflux_eceuler3d_gpu(diffFlux,grad,nu,kappa, &
                                           N,nvar,nel) &
      bind(c,name="diffusiveflux_eceuler3d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: diffFlux,grad
      real(c_prec),value :: nu,kappa
      integer(c_int),value :: N,nvar,nel
    endsubroutine diffusiveflux_eceuler3d_gpu
  endinterface

  interface
    subroutine diffusiveboundaryflux_eceuler3d_gpu(fluxN,avgGrad,uBnd,uExt, &
                                                   nhat,nscale, &
                                                   nu,kappa,tau_nu,tau_kappa, &
                                                   N,nvar,nel) &
      bind(c,name="diffusiveboundaryflux_eceuler3d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fluxN,avgGrad,uBnd,uExt,nhat,nscale
      real(c_prec),value :: nu,kappa,tau_nu,tau_kappa
      integer(c_int),value :: N,nvar,nel
    endsubroutine diffusiveboundaryflux_eceuler3d_gpu
  endinterface

  interface
    subroutine sourcemethod_eceuler3d_gpu(source,solution,dsdx,J,dSplit, &
                                          N,nvar,nel) &
      bind(c,name="sourcemethod_eceuler3d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: source,solution,dsdx,J,dSplit
      integer(c_int),value :: N,nvar,nel
    endsubroutine sourcemethod_eceuler3d_gpu
  endinterface

contains

  subroutine Init_ECEuler3D(this,mesh,geometry)
    implicit none
    class(ECEuler3D),intent(out) :: this
    type(Mesh3D),intent(in),target :: mesh
    type(SEMHex),intent(in),target :: geometry
    ! Local
    type(BoundaryCondition),pointer :: bc

    call Init_ECDGModel3D_t(this,mesh,geometry)

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
    ! GPU parabolic BC kernels (e.g. pbc3d_NoStress_ECEuler3D_GPU_wrapper).
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

  endsubroutine Init_ECEuler3D

  subroutine Free_ECEuler3D(this)
    implicit none
    class(ECEuler3D),intent(inout) :: this
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

    call Free_ECDGModel3D_t(this)

  endsubroutine Free_ECEuler3D

  subroutine AdditionalInit_ECEuler3D(this)
    implicit none
    class(ECEuler3D),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    ! Call parent _t AdditionalInit (registers CPU BC)
    call AdditionalInit_ECEuler3D_t(this)

    ! Re-register with GPU-accelerated versions for both lists.
    ! hyperbolicBCs and parabolicBCs are independent linked lists, so
    ! the same SELF_BC_NONORMALFLOW tag applies in both contexts.
    bcfunc => hbc3d_NoNormalFlow_ECEuler3D_GPU_wrapper
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

    bcfunc => pbc3d_NoStress_ECEuler3D_GPU_wrapper
    call this%parabolicBCs%RegisterBoundaryCondition( &
      SELF_BC_NONORMALFLOW,"no_normal_flow",bcfunc)

  endsubroutine AdditionalInit_ECEuler3D

  subroutine hbc3d_NoNormalFlow_ECEuler3D_GPU_wrapper(bc,mymodel)
    !! GPU-accelerated no-normal-flow BC for 3D EC Euler.
    !! Reflects normal momentum, mirrors density and rho*theta.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(ECEuler3D)
      if(bc%nBoundaries > 0) then
        call hbc3d_nonormalflow_eceuler3d_gpu( &
          m%solution%extBoundary_gpu, &
          m%solution%boundary_gpu, &
          m%geometry%nhat%boundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N, &
          m%solution%nElem,m%solution%nvar)
      endif
    endselect

  endsubroutine hbc3d_NoNormalFlow_ECEuler3D_GPU_wrapper

  subroutine pbc3d_NoStress_ECEuler3D_GPU_wrapper(bc,mymodel)
    !! GPU-accelerated parabolic no-stress / no-heat-flux BC for 3D EC Euler.
    !! Reflects the normal component of the solution gradient at every wall
    !! node so that BR1 averaging gives avgGrad . n = 0 (zero diffusive flux
    !! through the wall) for every variable.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(ECEuler3D)
      if(bc%nBoundaries > 0) then
        call pbc3d_nostress_eceuler3d_gpu( &
          m%solutionGradient%extBoundary_gpu, &
          m%solutionGradient%boundary_gpu, &
          m%geometry%nhat%boundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N, &
          m%solution%nElem,m%solution%nvar)
      endif
    endselect

  endsubroutine pbc3d_NoStress_ECEuler3D_GPU_wrapper

  subroutine BoundaryFlux_ECEuler3D(this)
    !! LMARS interface flux on GPU. No hydrostatic pressure split:
    !! gravity is handled by the Souza non-conservative source term
    !! using the geopotential carried as state variable index 6.
    implicit none
    class(ECEuler3D),intent(inout) :: this

    call boundaryflux_eceuler3d_gpu( &
      this%solution%boundary_gpu, &
      this%solution%extboundary_gpu, &
      this%geometry%nhat%boundary_gpu, &
      this%geometry%nscale%boundary_gpu, &
      this%flux%boundarynormal_gpu, &
      this%p0,this%Rd,this%cp/this%cv, &
      this%solution%interp%N, &
      this%solution%nelem)

  endsubroutine BoundaryFlux_ECEuler3D

  subroutine TwoPointFluxMethod_ECEuler3D(this)
    !! Souza et al. (2023) entropy-conservative two-point flux on GPU.
    !! Fully device-resident.
    implicit none
    class(ECEuler3D),intent(inout) :: this

    call twopointfluxmethod_eceuler3d_gpu( &
      this%twoPointFlux%interior_gpu, &
      this%solution%interior_gpu, &
      this%geometry%dsdx%interior_gpu, &
      this%p0,this%Rd,this%cp/this%cv, &
      this%solution%interp%N, &
      this%solution%nvar, &
      this%solution%nelem)

  endsubroutine TwoPointFluxMethod_ECEuler3D

  subroutine SourceMethod_ECEuler3D(this)
    !! Souza et al. (2023) non-conservative gravity flux differencing on
    !! GPU. The geopotential lives at solution(:,:,:,:,6); the source for
    !! rho*w is computed via the SBP-EC two-point form using log-mean
    !! density and the contravariant metric. Fully device-resident.
    implicit none
    class(ECEuler3D),intent(inout) :: this

    call sourcemethod_eceuler3d_gpu( &
      this%source%interior_gpu, &
      this%solution%interior_gpu, &
      this%geometry%dsdx%interior_gpu, &
      this%geometry%J%interior_gpu, &
      this%solution%interp%dSplitMatrix_gpu, &
      this%solution%interp%N, &
      this%solution%nvar, &
      this%solution%nelem)

  endsubroutine SourceMethod_ECEuler3D

  subroutine DiffusiveFluxMethod_ECEuler3D(this)
    !! GPU-resident fill of diffFlux%interior with the constant-coefficient
    !! Laplacian flux F_d(iVar) = -coeff(iVar) * d(s_iVar)/dx_d, where
    !! coeff is 0 for rho, nu for momentum, kappa for rho*theta.
    implicit none
    class(ECEuler3D),intent(inout) :: this

    call diffusiveflux_eceuler3d_gpu( &
      this%diffFlux%interior_gpu, &
      this%solutionGradient%interior_gpu, &
      this%nu,this%kappa, &
      this%solution%interp%N, &
      this%solution%nvar, &
      this%solution%nelem)

  endsubroutine DiffusiveFluxMethod_ECEuler3D

  subroutine DiffusiveBoundaryFlux_ECEuler3D(this)
    !! GPU-resident fill of diffFlux%boundaryNormal with the SIPG-stabilised
    !! BR1 diffusive flux:
    !!   f = -coeff*(avg_grad . n)*nmag + tau*(uL - uR)*nmag
    !! tau = eta_penalty*coeff*(N+1)^2/length_scale, computed on the host.
    implicit none
    class(ECEuler3D),intent(inout) :: this
    ! Local
    real(prec) :: np2,tau_nu,tau_kappa

    np2 = real((this%solution%interp%N+1)**2,prec)
    tau_nu = this%eta_penalty*this%nu*np2/this%length_scale
    tau_kappa = this%eta_penalty*this%kappa*np2/this%length_scale

    call diffusiveboundaryflux_eceuler3d_gpu( &
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

  endsubroutine DiffusiveBoundaryFlux_ECEuler3D

  subroutine CalculateTendency_ECEuler3D(this)
    !! GPU-resident tendency for ECEuler3D. The inviscid pipeline is
    !! identical to ECDGModel3D's GPU CalculateTendency; if either nu
    !! or kappa is positive, the constant-coefficient Laplacian
    !! divergence (BR1 weak-form) is then accumulated into
    !! fluxDivergence before forming dSdt.
    implicit none
    class(ECEuler3D),intent(inout) :: this
    ! Local
    integer :: ndof,ndof_diff

    call this%solution%BoundaryInterp()
    call this%solution%SideExchange(this%mesh)

    call this%PreTendency()
    call this%SetBoundaryCondition()

    if(this%gradient_enabled) then
      call this%CalculateSolutionGradient()
      call this%SetGradientBoundaryCondition()
      call this%solutionGradient%AverageSides()
    endif

    call this%SourceMethod()
    call this%BoundaryFlux()

    call this%TwoPointFluxMethod()

    call this%twoPointFlux%MappedDivergence(this%fluxDivergence%interior_gpu)

    call ECDGSurfaceContribution_3D_gpu( &
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
                  (this%solution%interp%N+1)* &
                  (this%solution%interp%N+1)
      call AccumulateField_gpu(this%fluxDivergence%interior_gpu, &
                               this%diffDiv%interior_gpu, &
                               ndof_diff)
    endif

    ndof = this%solution%nVar* &
           this%mesh%nElem* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)

    call CalculateDSDt_gpu(this%fluxDivergence%interior_gpu, &
                           this%source%interior_gpu, &
                           this%dSdt%interior_gpu,ndof)

  endsubroutine CalculateTendency_ECEuler3D

endmodule SELF_ECEuler3D
