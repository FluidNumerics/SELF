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
    call this%solutionGradient%SideExchange(this%mesh, &
                                            this%decomp)

  endsubroutine CalculateSolutionGradient_DGModel1D

  subroutine CalculateTendency_DGModel1D(this)
    implicit none
    class(DGModel1D),intent(inout) :: this
    ! Local
    integer :: ndof

    call this%solution%BoundaryInterp()
    call this%solution%SideExchange(this%mesh,this%decomp)

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

    call this%flux%MappedDGDerivative(this%fluxDivergence%interior_gpu)

    ndof = this%solution%nvar*this%solution%nelem*(this%solution%interp%N+1)
    call CalculateDSDt_gpu(this%fluxDivergence%interior_gpu,this%source%interior_gpu, &
                           this%dsdt%interior_gpu,ndof)

  endsubroutine CalculateTendency_DGModel1D

endmodule SELF_DGModel1D
