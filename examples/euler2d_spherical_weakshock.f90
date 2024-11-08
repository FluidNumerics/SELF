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

program Euler_Example

  use self_data
  use self_Euler2D

  implicit none

! Adapted from https://github.com/trixi-framework/Trixi.jl/blob/main/src/equations/compressible_euler_2d.jl#L157
!   A weak blast wave taken from
! - Sebastian Hennemann, Gregor J. Gassner (2020)
!   A provably entropy stable subcell shock capturing approach for high order split form DG
!   [arXiv: 2008.12044](https://arxiv.org/abs/2008.12044)
  real(prec),parameter :: rho0 = 1.0_prec
  real(prec),parameter :: rhoprime = 0.5_prec
  real(prec),parameter :: Lr = 3.0_prec*10.0_prec**(-2)
  real(prec),parameter :: P0 = 10.0_prec**(-5)
  real(prec),parameter :: Eprime = 1.0_prec
  real(prec),parameter :: Le = 2.0_prec*10.0_prec**(-2)
  real(prec),parameter :: x0 = 0.5_prec
  real(prec),parameter :: y0 = 0.5_prec
  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 15
  real(prec),parameter :: dt = 2.0_prec*10.0_prec**(-4) ! time-step size
  real(prec),parameter :: endtime = 1.0_prec
  real(prec),parameter :: iointerval = 0.1_prec
  real(prec) :: e0,ef ! Initial and final entropy
  type(Euler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)

  ! Create a structured mesh
  bcids(1:4) = [SELF_BC_NONORMALFLOW, & ! South
                SELF_BC_NONORMALFLOW, & ! East
                SELF_BC_NONORMALFLOW, & ! North
                SELF_BC_NONORMALFLOW] ! West

  call mesh%StructuredMesh(10,10,2,2,0.05_prec,0.05_prec,bcids)

  ! Create an interpolant
  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  ! Generate geometry (metric terms) from the mesh elements
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  ! Initialize the model
  call modelobj%Init(mesh,geometry)
  modelobj%prescribed_bcs_enabled = .false. ! Disables prescribed boundary condition block for gpu accelerated implementations

  call modelobj%WeakBlastWave(rho0,rhoprime,Lr,P0,Eprime,Le,x0,y0)

  call modelobj%WriteModel()
  call modelobj%IncrementIOCounter()

  call modelobj%CalculateEntropy()
  call modelobj%ReportEntropy()
  e0 = modelobj%entropy
  ! Set the model's time integration method
  call modelobj%SetTimeIntegrator(integrator)

  ! forward step the model to `endtime` using a time step
  ! of `dt` and outputing model data every `iointerval`
  call modelobj%ForwardStep(endtime,dt,iointerval)

  ef = modelobj%entropy

  if(ef > e0) then
    print*,"Error: Final absmax greater than initial absmax! ",e0,ef
    stop 1
  endif
  ! Clean up
  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram Euler_Example
