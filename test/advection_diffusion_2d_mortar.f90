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

program advection_diffusion_2d_mortar
!! Regression test for parabolic (Bassi-Rebay) terms across 2:1 nonconforming
!! (mortar) interfaces. An advection-diffusion tracer blob is released on the big
!! element of a mortar interface and advected across it while diffusing, exercising
!! the full gradient path on a nonconforming mesh : the solution mortar exchange,
!! the solution-gradient mortar exchange inside CalculateSolutionGradient, the
!! project-then-average interface gradient, and the mortar flux collect for the
!! combined advective-diffusive flux. The DoubleMortarMesh includes a flip = 1
!! sub-edge, so trace reorientation is exercised on every path.
!!
!! Assertions:
!!  (a) entropy (0.5*s^2 integrated over the domain) does not increase : advection
!!      with no-normal-flow walls is entropy-neutral and diffusion is strictly
!!      dissipative, so any growth indicates an interface error;
!!  (b) the solution remains free of NaNs.

  use self_data
  use self_advection_diffusion_2d
  use self_mesh_2d

  implicit none
  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  real(prec),parameter :: dx = 0.1_prec
  real(prec),parameter :: u = 0.25_prec ! velocity
  real(prec),parameter :: v = 0.1_prec
  real(prec),parameter :: nu = 0.005_prec ! diffusivity
  real(prec),parameter :: dt = 2.0_prec*10.0_prec**(-4) ! time-step size
  real(prec),parameter :: endtime = 0.05_prec
  real(prec),parameter :: iointerval = 0.05_prec
  real(prec) :: e0,ef ! Initial and final entropy
  type(advection_diffusion_2d) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)

  ! The mesh is constructed before the interpolant : the mesh's domain
  ! decomposition assigns each MPI rank its GPU device, and the interpolant's
  ! device arrays must be allocated on that device.
  bcids(1:4) = [SELF_BC_NONORMALFLOW, & ! south
                SELF_BC_NONORMALFLOW, & ! east
                SELF_BC_NONORMALFLOW, & ! north
                SELF_BC_NONORMALFLOW] ! west
  call mesh%DoubleMortarMesh(dx,bcids)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%gradient_enabled = .true.
  modelobj%tecplot_enabled = .false.

  modelobj%u = u
  modelobj%v = v
  modelobj%nu = nu

  ! Tracer blob straddling both mortar interfaces at x = 2*dx
  call modelobj%solution%SetEquation(1,'f = exp( -( (x-0.2)^2 + (y-0.2)^2 )/0.005 )')
  call modelobj%solution%SetInteriorFromEquation(geometry,0.0_prec)

  call modelobj%CalculateEntropy()
  call modelobj%ReportEntropy()
  e0 = modelobj%entropy

  call modelobj%SetTimeIntegrator(integrator)
  call modelobj%ForwardStep(endtime,dt,iointerval)

  ef = modelobj%entropy

  if(ef /= ef) then
    print*,"Error: entropy is NaN after crossing the mortar interface."
    stop 1
  endif
  if(ef > e0) then
    print*,"Error: entropy grew across the mortar interface :",e0,ef
    stop 1
  endif
  print*,"initial, final entropy :",e0,ef

  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram advection_diffusion_2d_mortar
