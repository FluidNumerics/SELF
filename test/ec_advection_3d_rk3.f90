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

program ec_advection_3d_rk3
  !! Entropy-stability test for the 3-D EC-DG advection model.
  !! Verifies the solution stays bounded over a RK3 integration.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_ECAdvection3D

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  real(prec),parameter :: u = 0.25_prec
  real(prec),parameter :: v = 0.25_prec
  real(prec),parameter :: w = 0.25_prec
  real(prec),parameter :: dt = 1.0_prec*10.0_prec**(-4)
  real(prec),parameter :: endtime = 0.1_prec
  real(prec),parameter :: iointerval = endtime
  real(prec) :: e0,ef
  type(ECAdvection3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS_LOBATTO, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  bcids(1:6) = [SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW]
  call mesh%StructuredMesh(3,3,3,1,1,1, &
                           1.0_prec/3.0_prec,1.0_prec/3.0_prec,1.0_prec/3.0_prec, &
                           bcids)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%u = u
  modelobj%v = v
  modelobj%w = w

  call modelobj%solution%SetEquation(1, &
                                     'f = exp( -( (x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2 )/0.01 )')
  call modelobj%solution%SetInteriorFromEquation(geometry,0.0_prec)

  call modelobj%CalculateEntropy()
  call modelobj%ReportEntropy()
  e0 = modelobj%entropy

  call modelobj%SetTimeIntegrator(integrator)
  call modelobj%ForwardStep(endtime,dt,iointerval)

  ef = modelobj%entropy
  call modelobj%solution%UpdateHost()

  if(maxval(abs(modelobj%solution%interior)) > 2.0_prec) then
    print*,"Error: EC-DG 3D advection solution blew up! max =", &
      maxval(abs(modelobj%solution%interior))
    stop 1
  endif

  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram ec_advection_3d_rk3
