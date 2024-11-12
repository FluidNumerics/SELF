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

program LinearShallowWater2D_nonormalflow_model
  use self_data
  use self_LinearShallowWater2D
  use self_mesh_2d

  implicit none
  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3' ! Which integrator method
  integer,parameter :: controlDegree = 7                            ! Degree of control polynomial
  integer,parameter :: targetDegree = 16                            ! Degree of target polynomial
  real(prec),parameter :: dt = 0.5_prec*10.0_prec**(-4)             ! Time-step size
  real(prec),parameter :: endtime = 1.0_prec                        ! Final time
  real(prec),parameter :: iointerval = 0.05_prec                    ! How often to write .tec files
  
  real(prec) :: e0,ef                                               ! Initial and final entropy
  type(LinearShallowWater2D) :: modelobj                                  ! Shallow water model
  type(Lagrange),target :: interp                                   ! Interpolant
  integer :: bcids(1:4)                                             ! Boundary conditions for structured mesh
  type(Mesh2D),target :: mesh                                       ! Mesh class
  type(SEMQuad),target :: geometry                                  ! Geometry class
  character(LEN=255) :: WORKSPACE                                   ! Used for file I/O

  real(prec),parameter :: g = 1.0_prec                              ! Acceleration due to gravity
  real(prec),parameter :: H = 1.0_prec                              ! Uniform resting depth

  ! Set no normal flow boundary conditions
  bcids(1:4) = [SELF_BC_NONORMALFLOW,& ! South
                SELF_BC_NONORMALFLOW,& ! East
                SELF_BC_NONORMALFLOW,& ! North
                SELF_BC_NONORMALFLOW]  ! West

  ! Create a uniform block mesh
  call mesh % StructuredMesh(10,10,2,2,0.05_prec,0.05_prec,bcids)

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

  ! Set the resting surface height and gravity
  modelobj%H = H
  modelobj%g = g

  ! Set the initial conditions
  call modelobj%solution%SetEquation(1,'f = 0')
  call modelobj%solution%SetEquation(2,'f = 0')
  call modelobj%solution%SetEquation(3,'f = 0.001*exp( -( (x-0.5)^2 + (y-0.5)^2 )/0.01 )')
  call modelobj%solution%SetInteriorFromEquation(geometry,0.0_prec)

  call modelobj%CalculateEntropy()
  e0 = modelobj%entropy

  ! Set the model's time integration method
  call modelobj%SetTimeIntegrator(integrator)

  ! forward step the model to `endtime` using a time step
  ! of `dt` and outputing model data every `iointerval`
  call modelobj%ForwardStep(endtime,dt,iointerval)

  ef = modelobj%entropy

  if(ef > e0) then
    print*,"Error: Final entropy greater than initial entropy! ",e0,ef
    stop 1
  endif
  
  ! Clean up
  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram LinearShallowWater2D_nonormalflow_model