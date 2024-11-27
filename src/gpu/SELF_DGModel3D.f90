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

  implicit none

  type,extends(DGModel3D_t) :: DGModel3D

  contains

    procedure :: UpdateSolution => UpdateSolution_DGModel3D

    procedure :: CalculateEntropy => CalculateEntropy_DGModel3D
    procedure :: BoundaryFlux => BoundaryFlux_DGModel3D
    procedure :: FluxMethod => fluxmethod_DGModel3D
    procedure :: SourceMethod => sourcemethod_DGModel3D
    procedure :: SetBoundaryCondition => setboundarycondition_DGModel3D
    procedure :: SetGradientBoundaryCondition => setgradientboundarycondition_DGModel3D

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
    ndof = this%solution%nvar* &
           this%solution%nelem* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)

    call UpdateSolution_gpu(this%solution%interior_gpu,this%dsdt%interior_gpu,dtLoc,ndof)

  endsubroutine UpdateSolution_DGModel3D

  subroutine UpdateGRK2_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: ndof

    ndof = this%solution%nvar* &
           this%solution%nelem* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)

    call UpdateGRK_gpu(this%worksol%interior_gpu,this%solution%interior_gpu,this%dsdt%interior_gpu, &
                       rk2_a(m),rk2_g(m),this%dt,ndof)

  endsubroutine UpdateGRK2_DGModel3D

  subroutine UpdateGRK3_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: ndof

    ndof = this%solution%nvar* &
           this%solution%nelem* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)

    call UpdateGRK_gpu(this%worksol%interior_gpu,this%solution%interior_gpu,this%dsdt%interior_gpu, &
                       rk3_a(m),rk3_g(m),this%dt,ndof)

  endsubroutine UpdateGRK3_DGModel3D

  subroutine UpdateGRK4_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: ndof

    ndof = this%solution%nvar* &
           this%solution%nelem* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)

    call UpdateGRK_gpu(this%worksol%interior_gpu,this%solution%interior_gpu,this%dsdt%interior_gpu, &
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
            e = e+this%entropy_func(s)*jac
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

  subroutine setboundarycondition_DGModel3D(this)
    !! Boundary conditions for the solution are set to
    !! 0 for the external state to provide radiation type
    !! boundary conditions.
    implicit none
    class(DGModel3D),intent(inout) :: this
    ! local
    integer :: i,iEl,j,k,e2,bcid
    real(prec) :: nhat(1:3),x(1:3)

    call gpuCheck(hipMemcpy(c_loc(this%solution%boundary), &
                            this%solution%boundary_gpu,sizeof(this%solution%boundary), &
                            hipMemcpyDeviceToHost))

    call gpuCheck(hipMemcpy(c_loc(this%solution%extboundary), &
                            this%solution%extboundary_gpu,sizeof(this%solution%extboundary), &
                            hipMemcpyDeviceToHost))

    do concurrent(k=1:6,iel=1:this%mesh%nElem)

      bcid = this%mesh%sideInfo(5,k,iEl) ! Boundary Condition ID
      e2 = this%mesh%sideInfo(3,k,iEl) ! Neighboring Element ID

      if(e2 == 0) then
        if(bcid == SELF_BC_PRESCRIBED) then

          do j = 1,this%solution%interp%N+1 ! Loop over quadrature points
            do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
              x = this%geometry%x%boundary(i,j,k,iEl,1,1:3)

              this%solution%extBoundary(i,j,k,iEl,1:this%nvar) = &
                this%hbc3d_Prescribed(x,this%t)
            enddo
          enddo

        elseif(bcid == SELF_BC_RADIATION) then

          do j = 1,this%solution%interp%N+1 ! Loop over quadrature points
            do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
              nhat = this%geometry%nhat%boundary(i,j,k,iEl,1,1:3)

              this%solution%extBoundary(i,j,k,iEl,1:this%nvar) = &
                this%hbc3d_Radiation(this%solution%boundary(i,j,k,iEl,1:this%nvar),nhat)
            enddo
          enddo

        elseif(bcid == SELF_BC_NONORMALFLOW) then

          do j = 1,this%solution%interp%N+1 ! Loop over quadrature points
            do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
              nhat = this%geometry%nhat%boundary(i,j,k,iEl,1,1:3)

              this%solution%extBoundary(i,j,k,iEl,1:this%nvar) = &
                this%hbc3d_NoNormalFlow(this%solution%boundary(i,j,k,iEl,1:this%nvar),nhat)
            enddo
          enddo

        endif
      endif

    enddo

    call gpuCheck(hipMemcpy(this%solution%extBoundary_gpu, &
                            c_loc(this%solution%extBoundary), &
                            sizeof(this%solution%extBoundary), &
                            hipMemcpyHostToDevice))

  endsubroutine setboundarycondition_DGModel3D

  subroutine setgradientboundarycondition_DGModel3D(this)
    !! Boundary conditions for the solution are set to
    !! 0 for the external state to provide radiation type
    !! boundary conditions.
    implicit none
    class(DGModel3D),intent(inout) :: this
    ! local
    integer :: i,iEl,j,k,e2,bcid
    real(prec) :: dsdx(1:this%nvar,1:3)
    real(prec) :: nhat(1:3),x(1:3)

    call gpuCheck(hipMemcpy(c_loc(this%solutiongradient%boundary), &
                            this%solutiongradient%boundary_gpu,sizeof(this%solutiongradient%boundary), &
                            hipMemcpyDeviceToHost))

    call gpuCheck(hipMemcpy(c_loc(this%solutiongradient%extboundary), &
                            this%solutiongradient%extboundary_gpu,sizeof(this%solutiongradient%extboundary), &
                            hipMemcpyDeviceToHost))

    do concurrent(k=1:6,iel=1:this%mesh%nElem)

      bcid = this%mesh%sideInfo(5,k,iEl) ! Boundary Condition ID
      e2 = this%mesh%sideInfo(3,k,iEl) ! Neighboring Element ID

      if(e2 == 0) then
        if(bcid == SELF_BC_PRESCRIBED) then

          do j = 1,this%solutiongradient%interp%N+1 ! Loop over quadrature points
            do i = 1,this%solutiongradient%interp%N+1 ! Loop over quadrature points
              x = this%geometry%nhat%boundary(i,j,k,iEl,1,1:3)

              this%solutiongradient%extBoundary(i,j,k,iEl,1:this%nvar,1:3) = &
                this%pbc3d_Prescribed(x,this%t)
            enddo
          enddo

        elseif(bcid == SELF_BC_RADIATION) then

          do j = 1,this%solutiongradient%interp%N+1 ! Loop over quadrature points
            do i = 1,this%solutiongradient%interp%N+1 ! Loop over quadrature points
              nhat = this%geometry%nhat%boundary(i,j,k,iEl,1,1:3)

              dsdx = this%solutiongradient%boundary(i,j,k,iEl,1:this%nvar,1:3)

              this%solutiongradient%extBoundary(i,j,k,iEl,1:this%nvar,1:3) = &
                this%pbc3d_Radiation(dsdx,nhat)
            enddo
          enddo

        elseif(bcid == SELF_BC_NONORMALFLOW) then

          do j = 1,this%solutiongradient%interp%N+1 ! Loop over quadrature points
            do i = 1,this%solutiongradient%interp%N+1 ! Loop over quadrature points
              nhat = this%geometry%nhat%boundary(i,j,k,iEl,1,1:3)

              dsdx = this%solutiongradient%boundary(i,j,k,iEl,1:this%nvar,1:3)

              this%solutiongradient%extBoundary(i,j,k,iEl,1:this%nvar,1:3) = &
                this%pbc3d_NoNormalFlow(dsdx,nhat)
            enddo
          enddo

        endif
      endif

    enddo

    call gpuCheck(hipMemcpy(this%solutiongradient%extBoundary_gpu, &
                            c_loc(this%solutiongradient%extBoundary), &
                            sizeof(this%solutiongradient%extBoundary), &
                            hipMemcpyHostToDevice))

  endsubroutine setgradientboundarycondition_DGModel3D

  subroutine CalculateTendency_DGModel3D(this)
    implicit none
    class(DGModel3D),intent(inout) :: this
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

    call this%flux%MappedDGDivergence(this%fluxDivergence%interior_gpu)

    ndof = this%solution%nvar* &
           this%solution%nelem* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)* &
           (this%solution%interp%N+1)

    call CalculateDSDt_gpu(this%fluxDivergence%interior_gpu,this%source%interior_gpu, &
                           this%dsdt%interior_gpu,ndof)

  endsubroutine CalculateTendency_DGModel3D

endmodule SELF_DGModel3D
