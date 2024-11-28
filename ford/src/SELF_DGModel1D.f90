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

module SELF_DGModel1D

  use SELF_SupportRoutines
  use SELF_Metadata
  use SELF_Mesh_1D
  use SELF_MappedScalar_1D
  use SELF_HDF5
  use HDF5
  use FEQParse
  use SELF_Model
  use SELF_DGModel1D_t
  use SELF_GPU
  use SELF_GPUInterfaces

  implicit none

  type,extends(DGModel1D_t) :: DGModel1D

  contains

    procedure :: UpdateSolution => UpdateSolution_DGModel1D

    procedure :: CalculateEntropy => CalculateEntropy_DGModel1D
    procedure :: BoundaryFlux => BoundaryFlux_DGModel1D
    procedure :: FluxMethod => fluxmethod_DGModel1D
    procedure :: SourceMethod => sourcemethod_DGModel1D
    procedure :: SetBoundaryCondition => setboundarycondition_DGModel1D
    procedure :: SetGradientBoundaryCondition => setgradientboundarycondition_DGModel1D

    procedure :: UpdateGRK2 => UpdateGRK2_DGModel1D
    procedure :: UpdateGRK3 => UpdateGRK3_DGModel1D
    procedure :: UpdateGRK4 => UpdateGRK4_DGModel1D

    procedure :: CalculateSolutionGradient => CalculateSolutionGradient_DGModel1D
    procedure :: CalculateTendency => CalculateTendency_DGModel1D

  endtype DGModel1D

contains

  subroutine UpdateSolution_DGModel1D(this,dt)
    !! Computes a solution update as , where dt is either provided through the interface
    !! or taken as the Model's stored time step size (model % dt)
    implicit none
    class(DGModel1D),intent(inout) :: this
    real(prec),optional,intent(in) :: dt
    ! Local
    real(prec) :: dtLoc
    integer :: ndof

    if(present(dt)) then
      dtLoc = dt
    else
      dtLoc = this%dt
    endif
    ndof = this%solution%nvar*this%solution%nelem*(this%solution%interp%N+1)

    call UpdateSolution_gpu(this%solution%interior_gpu,this%dsdt%interior_gpu,dtLoc,ndof)

  endsubroutine UpdateSolution_DGModel1D

  subroutine UpdateGRK2_DGModel1D(this,m)
    implicit none
    class(DGModel1D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: ndof

    ndof = this%solution%nvar*this%solution%nelem*(this%solution%interp%N+1)
    call UpdateGRK_gpu(this%worksol%interior_gpu,this%solution%interior_gpu,this%dsdt%interior_gpu, &
                       rk2_a(m),rk2_g(m),this%dt,ndof)

  endsubroutine UpdateGRK2_DGModel1D

  subroutine UpdateGRK3_DGModel1D(this,m)
    implicit none
    class(DGModel1D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: ndof

    ndof = this%solution%nvar*this%solution%nelem*(this%solution%interp%N+1)
    call UpdateGRK_gpu(this%worksol%interior_gpu,this%solution%interior_gpu,this%dsdt%interior_gpu, &
                       rk3_a(m),rk3_g(m),this%dt,ndof)

  endsubroutine UpdateGRK3_DGModel1D

  subroutine UpdateGRK4_DGModel1D(this,m)
    implicit none
    class(DGModel1D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: ndof

    ndof = this%solution%nvar*this%solution%nelem*(this%solution%interp%N+1)
    call UpdateGRK_gpu(this%worksol%interior_gpu,this%solution%interior_gpu,this%dsdt%interior_gpu, &
                       rk4_a(m),rk4_g(m),this%dt,ndof)

  endsubroutine UpdateGRK4_DGModel1D

  subroutine CalculateSolutionGradient_DGModel1D(this)
    implicit none
    class(DGModel1D),intent(inout) :: this
    ! Local
    integer :: ndof

    call this%solution%AverageSides()

    ! Account for the outward pointing normal before computing dg derivative
    ndof = this%solution%nvar*this%solution%nelem*2
    call GradientNormal_1D_gpu(this%solution%boundarynormal_gpu, &
                               this%solution%avgBoundary_gpu,ndof)

    call this%solution%MappedDGDerivative(this%solutionGradient%interior_gpu)

    ! interpolate the solutiongradient to the element boundaries
    call this%solutionGradient%BoundaryInterp()

    ! perform the side exchange to populate the
    ! solutionGradient % extBoundary attribute
    call this%solutionGradient%SideExchange(this%mesh)

  endsubroutine CalculateSolutionGradient_DGModel1D

  subroutine CalculateEntropy_DGModel1D(this)
    implicit none
    class(DGModel1D),intent(inout) :: this
    ! Local
    integer :: iel,i,ivar
    real(prec) :: e,s(1:this%solution%nvar),J

    call gpuCheck(hipMemcpy(c_loc(this%solution%interior), &
                            this%solution%interior_gpu,sizeof(this%solution%interior), &
                            hipMemcpyDeviceToHost))

    e = 0.0_prec
    do iel = 1,this%geometry%nelem
      do i = 1,this%solution%interp%N+1
        J = this%geometry%dxds%interior(i,iel,1)
        s(1:this%solution%nvar) = this%solution%interior(i,iel,1:this%solution%nvar)
        e = e+this%entropy_func(s)*J
      enddo
    enddo

    this%entropy = e

  endsubroutine CalculateEntropy_DGModel1D

  subroutine setboundarycondition_DGModel1D(this)
    ! Here, we use the pre-tendency method to calculate the
    ! derivative of the solution using a bassi-rebay method
    ! We then do a boundary interpolation and side exchange
    ! on the gradient field
    implicit none
    class(DGModel1D),intent(inout) :: this
    ! local
    integer :: ivar
    integer :: N,nelem
    real(prec) :: x

    call gpuCheck(hipMemcpy(c_loc(this%solution%boundary), &
                            this%solution%boundary_gpu,sizeof(this%solution%boundary), &
                            hipMemcpyDeviceToHost))

    nelem = this%geometry%nelem ! number of elements in the mesh
    N = this%solution%interp%N ! polynomial degree
    ! left-most boundary
    if(this%mesh%bcid(1) == SELF_BC_PRESCRIBED) then

      x = this%geometry%x%boundary(1,1,1)
      this%solution%extBoundary(1,1,1:this%nvar) = &
        this%hbc1d_Prescribed(x,this%t)

    elseif(this%mesh%bcid(1) == SELF_BC_RADIATION) then

      this%solution%extBoundary(1,1,1:this%nvar) = &
        this%hbc1d_Radiation(this%solution%boundary(1,1,1:this%nvar),-1.0_prec)

    elseif(this%mesh%bcid(1) == SELF_BC_NONORMALFLOW) then

      this%solution%extBoundary(1,1,1:this%nvar) = &
        this%hbc1d_NoNormalFlow(this%solution%boundary(1,1,1:this%nvar),-1.0_prec)

    else ! Periodic

      this%solution%extBoundary(1,1,1:this%nvar) = this%solution%boundary(2,nelem,1:this%nvar)

    endif

    ! right-most boundary
    if(this%mesh%bcid(1) == SELF_BC_PRESCRIBED) then

      x = this%geometry%x%boundary(2,nelem,1)
      this%solution%extBoundary(2,nelem,1:this%nvar) = &
        this%hbc1d_Prescribed(x,this%t)

    elseif(this%mesh%bcid(1) == SELF_BC_RADIATION) then

      this%solution%extBoundary(2,nelem,1:this%nvar) = &
        this%hbc1d_Radiation(this%solution%boundary(2,nelem,1:this%nvar),-1.0_prec)

    elseif(this%mesh%bcid(1) == SELF_BC_NONORMALFLOW) then

      this%solution%extBoundary(2,nelem,1:this%nvar) = &
        this%hbc1d_NoNormalFlow(this%solution%boundary(2,nelem,1:this%nvar),-1.0_prec)

    else ! Periodic

      this%solution%extBoundary(2,nelem,1:this%nvar) = this%solution%boundary(1,1,1:this%nvar)

    endif

    call gpuCheck(hipMemcpy(this%solution%extBoundary_gpu, &
                            c_loc(this%solution%extBoundary), &
                            sizeof(this%solution%extBoundary), &
                            hipMemcpyHostToDevice))

  endsubroutine setboundarycondition_DGModel1D

  subroutine setgradientboundarycondition_DGModel1D(this)
    ! Here, we set the boundary conditions for the
    ! solution and the solution gradient at the left
    ! and right most boundaries.
    !
    ! Here, we use periodic boundary conditions
    implicit none
    class(DGModel1D),intent(inout) :: this
    ! local
    integer :: ivar
    integer :: nelem
    real(prec) :: x

    call gpuCheck(hipMemcpy(c_loc(this%solutiongradient%boundary), &
                            this%solutiongradient%boundary_gpu,sizeof(this%solutiongradient%boundary), &
                            hipMemcpyDeviceToHost))

    nelem = this%geometry%nelem ! number of elements in the mesh

    ! left-most boundary
    if(this%mesh%bcid(1) == SELF_BC_PRESCRIBED) then

      x = this%geometry%x%boundary(1,1,1)
      this%solutionGradient%extBoundary(1,1,1:this%nvar) = &
        this%pbc1d_Prescribed(x,this%t)

    elseif(this%mesh%bcid(1) == SELF_BC_RADIATION) then

      this%solutionGradient%extBoundary(1,1,1:this%nvar) = &
        this%pbc1d_Radiation(this%solutionGradient%boundary(1,1,1:this%nvar),-1.0_prec)

    elseif(this%mesh%bcid(1) == SELF_BC_NONORMALFLOW) then

      this%solutionGradient%extBoundary(1,1,1:this%nvar) = &
        this%pbc1d_NoNormalFlow(this%solutionGradient%boundary(1,1,1:this%nvar),-1.0_prec)

    else ! Periodic

      this%solutionGradient%extBoundary(1,1,1:this%nvar) = this%solutionGradient%boundary(2,nelem,1:this%nvar)

    endif

    ! right-most boundary
    if(this%mesh%bcid(1) == SELF_BC_PRESCRIBED) then

      x = this%geometry%x%boundary(2,nelem,1)
      this%solutionGradient%extBoundary(2,nelem,1:this%nvar) = &
        this%pbc1d_Prescribed(x,this%t)

    elseif(this%mesh%bcid(1) == SELF_BC_RADIATION) then

      this%solutionGradient%extBoundary(2,nelem,1:this%nvar) = &
        this%pbc1d_Radiation(this%solutionGradient%boundary(2,nelem,1:this%nvar),-1.0_prec)

    elseif(this%mesh%bcid(1) == SELF_BC_NONORMALFLOW) then

      this%solutionGradient%extBoundary(2,nelem,1:this%nvar) = &
        this%pbc1d_NoNormalFlow(this%solutionGradient%boundary(2,nelem,1:this%nvar),-1.0_prec)

    else ! Periodic

      this%solutionGradient%extBoundary(2,nelem,1:this%nvar) = this%solutionGradient%boundary(1,1,1:this%nvar)

    endif

    call gpuCheck(hipMemcpy(this%solutiongradient%extBoundary_gpu, &
                            c_loc(this%solutiongradient%extBoundary), &
                            sizeof(this%solutiongradient%extBoundary), &
                            hipMemcpyHostToDevice))

  endsubroutine setgradientboundarycondition_DGModel1D

  subroutine BoundaryFlux_DGModel1D(this)
    ! this method uses an linear upwind solver for the
    ! advective flux and the bassi-rebay method for the
    ! diffusive fluxes
    implicit none
    class(DGModel1D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    real(prec) :: fin(1:this%solution%nvar)
    real(prec) :: fout(1:this%solution%nvar)
    real(prec) :: dfdx(1:this%solution%nvar),nhat

    call gpuCheck(hipMemcpy(c_loc(this%solution%boundary), &
                            this%solution%boundary_gpu,sizeof(this%solution%boundary), &
                            hipMemcpyDeviceToHost))

    call gpuCheck(hipMemcpy(c_loc(this%solution%extboundary), &
                            this%solution%extboundary_gpu,sizeof(this%solution%extboundary), &
                            hipMemcpyDeviceToHost))

    call gpuCheck(hipMemcpy(c_loc(this%solutiongradient%avgboundary), &
                            this%solutiongradient%avgboundary_gpu,sizeof(this%solutiongradient%avgboundary), &
                            hipMemcpyDeviceToHost))

    do concurrent(iside=1:2,iel=1:this%mesh%nElem)

      ! set the normal velocity
      if(iside == 1) then
        nhat = -1.0_prec
      else
        nhat = 1.0_prec
      endif

      fin = this%solution%boundary(iside,iel,1:this%solution%nvar) ! interior solution
      fout = this%solution%extboundary(iside,iel,1:this%solution%nvar) ! exterior solution
      dfdx = this%solutionGradient%avgboundary(iside,iel,1:this%solution%nvar) ! average solution gradient (with direction taken into account)
      this%flux%boundarynormal(iside,iel,1:this%solution%nvar) = &
        this%riemannflux1d(fin,fout,dfdx,nhat)

    enddo

    call gpuCheck(hipMemcpy(this%flux%boundarynormal_gpu, &
                            c_loc(this%flux%boundarynormal), &
                            sizeof(this%flux%boundarynormal), &
                            hipMemcpyHostToDevice))

  endsubroutine BoundaryFlux_DGModel1D

  subroutine fluxmethod_DGModel1D(this)
    implicit none
    class(DGModel1D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: i
    real(prec) :: f(1:this%solution%nvar),dfdx(1:this%solution%nvar)

    do concurrent(i=1:this%solution%N+1,iel=1:this%mesh%nElem)

      f = this%solution%interior(i,iel,1:this%solution%nvar)
      dfdx = this%solutionGradient%interior(i,iel,1:this%solution%nvar)

      this%flux%interior(i,iel,1:this%solution%nvar) = &
        this%flux1d(f,dfdx)

    enddo

    call gpuCheck(hipMemcpy(this%flux%interior_gpu, &
                            c_loc(this%flux%interior), &
                            sizeof(this%flux%interior), &
                            hipMemcpyHostToDevice))

  endsubroutine fluxmethod_DGModel1D

  subroutine sourcemethod_DGModel1D(this)
    implicit none
    class(DGModel1D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: i
    real(prec) :: f(1:this%solution%nvar),dfdx(1:this%solution%nvar)

    call gpuCheck(hipMemcpy(c_loc(this%solution%interior), &
                            this%solution%interior_gpu,sizeof(this%solution%interior), &
                            hipMemcpyDeviceToHost))

    call gpuCheck(hipMemcpy(c_loc(this%solutiongradient%interior), &
                            this%solutiongradient%interior_gpu,sizeof(this%solutiongradient%interior), &
                            hipMemcpyDeviceToHost))

    do concurrent(i=1:this%solution%N+1,iel=1:this%mesh%nElem)

      f = this%solution%interior(i,iel,1:this%solution%nvar)
      dfdx = this%solutionGradient%interior(i,iel,1:this%solution%nvar)

      this%source%interior(i,iel,1:this%solution%nvar) = &
        this%source1d(f,dfdx)

    enddo

    call gpuCheck(hipMemcpy(this%source%interior_gpu, &
                            c_loc(this%source%interior), &
                            sizeof(this%source%interior), &
                            hipMemcpyHostToDevice))

  endsubroutine sourcemethod_DGModel1D

  subroutine CalculateTendency_DGModel1D(this)
    implicit none
    class(DGModel1D),intent(inout) :: this
    ! Local
    integer :: ndof

    call this%solution%BoundaryInterp()
    call this%solution%SideExchange(this%mesh)

    call this%PreTendency() ! User-supplied
    call this%SetBoundaryCondition() ! User-supplied

    if(this%gradient_enabled) then
      call this%solution%AverageSides()
      call this%CalculateSolutionGradient()
      call this%SetGradientBoundaryCondition() ! User-supplied
      call this%solutionGradient%AverageSides()
    endif

    call this%SourceMethod() ! User supplied
    call this%BoundaryFlux() ! User supplied
    call this%FluxMethod() ! User supplied

    call this%flux%MappedDGDerivative(this%fluxDivergence%interior_gpu)

    ndof = this%solution%nvar*this%solution%nelem*(this%solution%interp%N+1)
    call CalculateDSDt_gpu(this%fluxDivergence%interior_gpu,this%source%interior_gpu, &
                           this%dsdt%interior_gpu,ndof)

  endsubroutine CalculateTendency_DGModel1D

endmodule SELF_DGModel1D
