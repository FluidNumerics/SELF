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
    call this%RiemannSolver() ! User supplied
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
