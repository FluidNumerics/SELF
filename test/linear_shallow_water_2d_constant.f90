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

program LinearShallowWater2D_constant
  use self_data
  use self_LinearShallowWater2D

  implicit none
  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'

  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  real(prec),parameter :: H = 1.0_prec ! uniform resting depth
  real(prec),parameter :: g = 9.8_prec ! acceleration due to gravity
  real(prec),parameter :: dt = 1.0_prec*10.0_prec**(-4) ! time-step size
  real(prec),parameter :: endtime = 0.2_prec
  real(prec),parameter :: iointerval = 0.1_prec
  real(prec) :: e0,ef ! Initial and final entropy
  type(LinearShallowWater2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  integer :: bcids(1:4)
  type(SEMQuad),target :: geometry
  character(LEN=255) :: WORKSPACE

  ! Set boundary conditions
  bcids(1:4) = [SELF_BC_PRESCRIBED, & ! South
                SELF_BC_PRESCRIBED, & ! East
                SELF_BC_PRESCRIBED, & ! North
                SELF_BC_PRESCRIBED] ! West

  ! Create a uniform block mesh
  call mesh%StructuredMesh(5,5,2,2,0.1_prec,0.1_prec,bcids)

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
  modelobj%gradient_enabled = .true.

  ! Set the resting surface height and gravity
  modelobj%H = H
  modelobj%g = g

  ! Set the initial condition
  call modelobj%solution%SetEquation(1,'f = 1.0')
  call modelobj%solution%SetEquation(2,'f = 1.0')
  call modelobj%solution%SetEquation(3,'f = 1.0')
  call modelobj%solution%SetInteriorFromEquation(geometry,0.0_prec)

  call modelobj%CalculateTendency()
  print*,"min, max (interior)", &
    minval(modelobj%solution%interior), &
    maxval(modelobj%solution%interior)

  call modelobj%CalculateEntropy()
  e0 = modelobj%entropy

  !Write the initial condition
  call modelobj%WriteModel()
  call modelobj%WriteTecplot()
  call modelobj%IncrementIOCounter()

  ! Set the model's time integration method
  call modelobj%SetTimeIntegrator(integrator)

  ! forward step the model to `endtime` using a time step
  ! of `dt` and outputing model data every `iointerval`
  call modelobj%ForwardStep(endtime,dt,iointerval)

  print*,"min, max (interior)", &
    minval(modelobj%solution%interior), &
    maxval(modelobj%solution%interior)

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

endprogram LinearShallowWater2D_constant