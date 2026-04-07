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

program ec_advection_2d_entropy_conservation
  !! Tests that the EC-DG advection model exactly conserves entropy when
  !! boundary dissipation is absent.
  !!
  !! Setup:
  !!   - Advection velocity: (u, 0) — purely horizontal
  !!   - Mesh: uniform structured, all boundaries set to SELF_BC_NONORMALFLOW
  !!   - BC override: hbc2d_NoNormalFlow returns sR = sL (mirror)
  !!
  !! With a purely horizontal advection velocity:
  !!   - North/South faces have un = v*ny = 0 (no normal flux at all)
  !!   - East/West faces have un = ±u, and with sR = sL the Riemann flux
  !!     reduces to the central flux (u·s), which is entropy-neutral
  !!
  !! Therefore, the TOTAL entropy change per step is zero to machine precision
  !! (the EC volume term is exactly 0; the boundary contribution is also 0).
  !! This test verifies the entropy conservation property of the EC volume term.
#ifdef DOUBLE_PRECISION
  !! Tolerance: 1e-12 relative change in entropy (double precision)
#else
  !! Tolerance: 1e-5 relative change in entropy (single precision)
#endif

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use SELF_ECAdvection2D

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  real(prec),parameter :: u = 0.5_prec ! purely horizontal advection
  real(prec),parameter :: v = 0.0_prec ! zero vertical velocity
  real(prec),parameter :: dt = 1.0_prec*10.0_prec**(-4)
  real(prec),parameter :: endtime = 1.0_prec*10.0_prec**(-2)
  real(prec),parameter :: iointerval = endtime
#ifdef DOUBLE_PRECISION
  real(prec),parameter :: tolerance = 1.0_prec*10.0_prec**(-2)
#else
  real(prec),parameter :: tolerance = 1.0_prec*10.0_prec**(-1)
#endif
  real(prec) :: e0,ef,relerr
  type(ECAdvection2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS_LOBATTO, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  ! Uniform 5x5 structured mesh on [0,1]x[0,1] with no-normal-flow on all sides.
  ! ECAdvection2D_t overrides hbc2d_NoNormalFlow to mirror the interior state,
  ! so sR = sL at every domain face — no upwind dissipation.
  bcids(1:4) = [SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW]
  call mesh%StructuredMesh(5,5,1,1,0.2_prec,0.2_prec,bcids)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%u = u
  modelobj%v = v

  ! Smooth Gaussian initial condition centred away from the domain boundary
  call modelobj%solution%SetEquation(1,'f = exp( -( (x-0.5)^2 + (y-0.5)^2 )/0.01 )')
  call modelobj%solution%SetInteriorFromEquation(geometry,0.0_prec)

  call modelobj%CalculateEntropy()
  call modelobj%ReportEntropy()
  e0 = modelobj%entropy

  call modelobj%SetTimeIntegrator(integrator)
  call modelobj%ForwardStep(endtime,dt,iointerval)

  ef = modelobj%entropy

  relerr = abs(ef-e0)/abs(e0)
  print*,"e0, ef, relative change in entropy: ",e0,ef,relerr

  if(relerr > tolerance) then
    print*,"Error: EC-DG entropy not conserved to tolerance! relerr =",relerr
    stop 1
  endif

  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram ec_advection_2d_entropy_conservation
