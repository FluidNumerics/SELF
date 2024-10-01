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

program advection_diffusion_3d_rk4

  use self_data
  use self_advection_diffusion_3d

  implicit none
  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk4'
  integer,parameter :: nvar = 1
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  real(prec),parameter :: u = 0.25_prec ! velocity
  real(prec),parameter :: v = 0.25_prec
  real(prec),parameter :: w = 0.25_prec
  real(prec),parameter :: nu = 0.001_prec ! diffusivity
  real(prec),parameter :: dt = 1.0_prec*10.0_prec**(-4) ! time-step size
  real(prec),parameter :: endtime = 0.01_prec
  real(prec),parameter :: iointerval = 0.01_prec
  type(advection_diffusion_3d) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  character(LEN=255) :: WORKSPACE

  ! Create a uniform block mesh
  call get_environment_variable("WORKSPACE",WORKSPACE)
  call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block3D/Block3D_mesh.h5")

  ! Create an interpolant
  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  ! Generate geometry (metric terms) from the mesh elements
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  ! Initialize the model
  call modelobj%Init(nvar,mesh,geometry)
  modelobj%gradient_enabled = .true.

  ! Set the velocity
  modelobj%u = u
  modelobj%v = v
  modelobj%w = w
  !Set the diffusivity
  modelobj%nu = nu

  ! Set the initial condition
  call modelobj%solution%SetEquation(1,'f = exp( -( (x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2 )/0.005 )')
  call modelobj%solution%SetInteriorFromEquation(geometry,0.0_prec)

  print*,"min, max (interior)", &
    minval(modelobj%solution%interior), &
    maxval(modelobj%solution%interior)

  ! Set the model's time integration method
  call modelobj%SetTimeIntegrator(integrator)

  ! forward step the model to `endtime` using a time step
  ! of `dt` and outputing model data every `iointerval`
  call modelobj%ForwardStep(endtime,dt,iointerval)

  print*,"min, max (interior)", &
    minval(modelobj%solution%interior), &
    maxval(modelobj%solution%interior)

  ! Clean up
  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram advection_diffusion_3d_rk4
