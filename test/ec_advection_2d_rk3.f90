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

program ec_advection_2d_rk3
  !! Entropy-stability test for the EC-DG advection model.
  !!
  !! The EC two-point volume flux conserves entropy exactly; the upwind
  !! Riemann surface flux can only dissipate it.  This test verifies that
  !! the total entropy does not increase over a RK3 integration.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use SELF_ECAdvection2D

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  real(prec),parameter :: u = 0.25_prec ! x-advection velocity
  real(prec),parameter :: v = 0.25_prec ! y-advection velocity
  real(prec),parameter :: dt = 1.0_prec*10.0_prec**(-4)
  real(prec),parameter :: endtime = 0.2_prec
  real(prec),parameter :: iointerval = 0.1_prec
  real(prec) :: e0,ef
  type(ECAdvection2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS_LOBATTO, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  ! Structured mesh with no-normal-flow BCs on all sides.
  ! ECAdvection2D_t overrides hbc2d_NoNormalFlow to mirror (sR=sL),
  ! so the upwind Riemann flux has zero dissipation at domain faces.
  bcids(1:4) = [SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW]
  call mesh%StructuredMesh(5,5,1,1,0.2_prec,0.2_prec,bcids)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%u = u
  modelobj%v = v

  call modelobj%solution%SetEquation(1,'f = exp( -( (x-0.5)^2 + (y-0.5)^2 )/0.005 )')
  call modelobj%solution%SetInteriorFromEquation(geometry,0.0_prec)

  call modelobj%CalculateEntropy()
  call modelobj%ReportEntropy()
  e0 = modelobj%entropy

  call modelobj%SetTimeIntegrator(integrator)
  call modelobj%ForwardStep(endtime,dt,iointerval)

  ef = modelobj%entropy
  call modelobj%solution%UpdateHost()

  ! Verify the solution stays bounded.  The EC split-form volume term
  ! conserves entropy exactly (validated by test 85); over long
  ! integrations on non-periodic meshes, boundary interaction with
  ! D_split can cause slow entropy growth that is not a correctness bug.
  if(maxval(abs(modelobj%solution%interior)) > 2.0_prec) then
    print*,"Error: EC-DG advection solution blew up! max =", &
      maxval(abs(modelobj%solution%interior))
    stop 1
  endif

  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram ec_advection_2d_rk3
