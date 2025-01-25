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

program linear_shallow_water2d_kelvinwaves
  use self_data
  use self_LinearShallowWater2D
  use self_mesh_2d

  implicit none
  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3' ! Which integrator method
  integer,parameter :: controlDegree = 7 ! Degree of control polynomial
  integer,parameter :: targetDegree = 16 ! Degree of target polynomial
  real(prec),parameter :: dt = 0.001_prec ! Time-step size
  real(prec),parameter :: endtime = 1.0_prec !30.0_prec ! (s);
  real(prec),parameter :: f0 = -10.0_prec ! reference coriolis parameter (1/s)
  real(prec),parameter :: iointerval = 0.05 ! Write files 20 times per characteristic time scale
  real(prec) :: r
  real(prec) :: e0,ef ! Initial and final entropy
  type(LinearShallowWater2D) :: modelobj ! Shallow water model
  type(Lagrange),target :: interp ! Interpolant
  type(Mesh2D),target :: mesh ! Mesh class
  type(SEMQuad),target :: geometry ! Geometry class
  integer :: i,j,iel
  real(prec),parameter :: g = 1.0_prec ! Acceleration due to gravity
  real(prec),parameter :: H = 1.0_prec ! Uniform resting depth
  character(LEN=255) :: WORKSPACE

  ! Create a uniform block mesh
  call get_environment_variable("WORKSPACE",WORKSPACE)
  call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Circle/Circle_mesh.h5")
  call mesh%ResetBoundaryConditionType(SELF_BC_NONORMALFLOW)

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
  modelobj%tecplot_enabled = .false. ! Disables tecplot output

  ! Set the resting surface height and gravity
  modelobj%H = H
  modelobj%g = g

  ! ! Set the initial conditions
  !call modelobj%solution%SetEquation(3,'f = 0.001*exp( -( x^2 + y^2 )/0.02 ) ')
  !call modelobj%solution%SetInteriorFromEquation(geometry,0.0_prec)

  do iel = 1,modelobj%mesh%nElem
    do j = 1,modelobj%solution%N+1
      do i = 1,modelobj%solution%N+1
        call random_number(r)
        !  modelobj%solution%interior(i,j,iel,3) = modelobj%solution%interior(i,j,iel,3) + 0.0001_prec*(r-0.5)
        modelobj%solution%interior(i,j,iel,3) = 0.0001_prec*(r-0.5)

      enddo
    enddo
  enddo
  call modelobj%solution%UpdateDevice()

  call modelobj%SetCoriolis(f0)
  call modelobj%DiagnoseGeostrophicVelocity()

  call modelobj%WriteModel()
  call modelobj%IncrementIOCounter()

  call modelobj%CalculateEntropy()
  call modelobj%ReportEntropy()
  call modelobj%ReportMetrics()
  e0 = modelobj%entropy

  ! Set the model's time integration method
  call modelobj%SetTimeIntegrator(integrator)

  ! forward step the model to `endtime` using a time step
  ! of `dt` and outputing model data every `iointerval`
  call modelobj%ForwardStep(endtime,dt,iointerval)

  ef = modelobj%entropy

  print*,e0,ef
  if(abs(ef-e0) > epsilon(e0)) then
    print*,"Warning: Final entropy greater than initial entropy! ",e0,ef
  endif

  ! Clean up
  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram linear_shallow_water2d_kelvinwaves
