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

program ec_euler3d_motionless
  !! Tests that the EC-DG Euler 3-D model maintains a motionless fluid
  !! in a closed domain (no gravity).
  !!
  !! Setup:
  !!   - Uniform state: rho=1.0, u=v=w=0, theta=300 K
  !!   - Mesh: 3x3x3 structured, all boundaries SELF_BC_NONORMALFLOW
  !!   - Gravity: g = 0 (no external forcing)
  !!   - Time integration: RK3, 100 time steps
  !!
  !! Pass criteria:
  !!   - Maximum absolute change in any solution variable < tolerance

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_ECEuler3D

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  real(prec),parameter :: dt = 1.0_prec*10.0_prec**(-4)
  real(prec),parameter :: endtime = 1.0_prec*10.0_prec**(-2)
  real(prec),parameter :: iointerval = endtime
  real(prec),parameter :: rho0 = 1.0_prec
  real(prec),parameter :: theta0 = 300.0_prec
#ifdef DOUBLE_PRECISION
  real(prec),parameter :: tolerance = 1.0_prec*10.0_prec**(-10)
#else
  real(prec),parameter :: tolerance = 1.0_prec*10.0_prec**(-4)
#endif
  type(ECEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)
  real(prec) :: maxerr
  integer :: ivar

  ! Create interpolant (Gauss-Lobatto for EC split form)
  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS_LOBATTO, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  ! Create structured mesh with no-normal-flow on all faces
  bcids(1:6) = [SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW]
  call mesh%StructuredMesh(3,3,3,1,1,1, &
                           1.0_prec/3.0_prec,1.0_prec/3.0_prec, &
                           1.0_prec/3.0_prec,bcids)

  ! Generate geometry
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  ! Initialize model
  call modelobj%Init(mesh,geometry)
  modelobj%g = 0.0_prec ! No gravity for this test

  ! Set uniform initial condition
  ! rho = 1.0
  call modelobj%solution%SetEquation(1,'f = 1.0')
  ! rhou = 0.0
  call modelobj%solution%SetEquation(2,'f = 0.0')
  ! rhov = 0.0
  call modelobj%solution%SetEquation(3,'f = 0.0')
  ! rhow = 0.0
  call modelobj%solution%SetEquation(4,'f = 0.0')
  ! rhotheta = rho0 * theta0 = 300.0
  call modelobj%solution%SetEquation(5,'f = 300.0')
  call modelobj%solution%SetInteriorFromEquation(geometry,0.0_prec)

  print*,"Initial state:"
  do ivar = 1,5
    print*,"  var ",ivar,": min/max = ", &
      minval(modelobj%solution%interior(:,:,:,:,ivar)), &
      maxval(modelobj%solution%interior(:,:,:,:,ivar))
  enddo

  ! Time integrate
  call modelobj%SetTimeIntegrator(integrator)
  modelobj%tecplot_enabled = .false.
  call modelobj%ForwardStep(endtime,dt,iointerval)

  ! Check that the solution has not changed
  maxerr = 0.0_prec
  maxerr = max(maxerr,maxval(abs( &
                             modelobj%solution%interior(:,:,:,:,1)-rho0)))
  maxerr = max(maxerr,maxval(abs( &
                             modelobj%solution%interior(:,:,:,:,2)-0.0_prec)))
  maxerr = max(maxerr,maxval(abs( &
                             modelobj%solution%interior(:,:,:,:,3)-0.0_prec)))
  maxerr = max(maxerr,maxval(abs( &
                             modelobj%solution%interior(:,:,:,:,4)-0.0_prec)))
  maxerr = max(maxerr,maxval(abs( &
                             modelobj%solution%interior(:,:,:,:,5)-rho0*theta0)))

  print*,"Final state:"
  do ivar = 1,5
    print*,"  var ",ivar,": min/max = ", &
      minval(modelobj%solution%interior(:,:,:,:,ivar)), &
      maxval(modelobj%solution%interior(:,:,:,:,ivar))
  enddo
  print*,"Maximum absolute error: ",maxerr

  if(maxerr > tolerance) then
    print*,"FAIL: motionless fluid not maintained! maxerr =",maxerr
    stop 1
  endif

  print*,"PASS: motionless fluid maintained to tolerance"

  ! Clean up
  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram ec_euler3d_motionless
