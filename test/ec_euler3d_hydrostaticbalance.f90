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

program ec_euler3d_hydrostaticbalance
  !! Tests that the EC-DG Euler 3-D model maintains a hydrostatically
  !! balanced atmosphere with uniform potential temperature.
  !!
  !! Setup:
  !!   - theta0 = 300 K (uniform potential temperature)
  !!   - Domain: [0, 1000m] x [0, 1000m] x [0, 1500m]
  !!   - Mesh: 2x2x3 elements (dx=dy=dz=500m)
  !!   - Gravity: g = 9.81 m/s^2
  !!   - All boundaries: no-normal-flow
  !!   - Time integration: RK3, dt=0.01s, 100 steps
  !!
  !! Hydrostatic profile:
  !!   pi(z) = 1 - g*z/(cp*theta0)
  !!   rho(z) = p0/(Rd*theta0) * pi(z)^(cv/Rd)
  !!
  !! Pass criteria:
  !!   - Maximum velocity magnitude stays below tolerance
  !!   - Relative change in density and rho*theta stay below tolerance

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_ECEuler3D

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  real(prec),parameter :: theta0 = 300.0_prec
  real(prec),parameter :: dt = 1.0_prec*10.0_prec**(-2)
  real(prec),parameter :: endtime = 1.0_prec*10.0_prec**(-1)
  real(prec),parameter :: iointerval = endtime
#ifdef DOUBLE_PRECISION
  real(prec),parameter :: tolerance = 1.0_prec*10.0_prec**(-4)
#else
  real(prec),parameter :: tolerance = 1.0_prec*10.0_prec**(0)
#endif
  type(ECEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)
  real(prec) :: maxvel,maxreldiff_rho,maxreldiff_rhotheta
  real(prec),allocatable :: rho0(:,:,:,:),rhotheta0(:,:,:,:)
  integer :: nN,nEl

  ! Create interpolant (Gauss-Lobatto for EC split form)
  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS_LOBATTO, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  ! Create structured mesh: 2x2x3 elements, dx=dy=dz=500m
  ! Domain: [0, 1000] x [0, 1000] x [0, 1500]
  bcids(1:6) = [SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW]
  call mesh%StructuredMesh(2,2,3,1,1,1, &
                           500.0_prec,500.0_prec,500.0_prec,bcids)

  ! Generate geometry
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  ! Initialize model
  call modelobj%Init(mesh,geometry)

  ! Set hydrostatic balance
  call modelobj%SetHydrostaticBalance(theta0)

  ! Store initial state for comparison
  nN = controlDegree+1
  nEl = mesh%nElem
  allocate(rho0(nN,nN,nN,nEl))
  allocate(rhotheta0(nN,nN,nN,nEl))
  rho0 = modelobj%solution%interior(:,:,:,:,1)
  rhotheta0 = modelobj%solution%interior(:,:,:,:,5)

  print*,"Initial state:"
  print*,"  rho min/max      = ", &
    minval(modelobj%solution%interior(:,:,:,:,1)), &
    maxval(modelobj%solution%interior(:,:,:,:,1))
  print*,"  rhotheta min/max = ", &
    minval(modelobj%solution%interior(:,:,:,:,5)), &
    maxval(modelobj%solution%interior(:,:,:,:,5))

  ! Time integrate
  call modelobj%SetTimeIntegrator(integrator)
  modelobj%tecplot_enabled = .false.
  call modelobj%ForwardStep(endtime,dt,iointerval)

  ! Check that the velocity remains near zero
  maxvel = max( &
           maxval(abs(modelobj%solution%interior(:,:,:,:,2))), &
           maxval(abs(modelobj%solution%interior(:,:,:,:,3))), &
           maxval(abs(modelobj%solution%interior(:,:,:,:,4))))

  ! Check relative change in density
  maxreldiff_rho = maxval(abs( &
                          modelobj%solution%interior(:,:,:,:,1)-rho0)/abs(rho0))

  ! Check relative change in rho*theta
  maxreldiff_rhotheta = maxval(abs( &
                               modelobj%solution%interior(:,:,:,:,5)-rhotheta0)/abs(rhotheta0))

  print*,"Final state:"
  print*,"  rho min/max      = ", &
    minval(modelobj%solution%interior(:,:,:,:,1)), &
    maxval(modelobj%solution%interior(:,:,:,:,1))
  print*,"  rhotheta min/max = ", &
    minval(modelobj%solution%interior(:,:,:,:,5)), &
    maxval(modelobj%solution%interior(:,:,:,:,5))
  print*,"  max |momentum|   = ",maxvel
  print*,"  max rel diff rho = ",maxreldiff_rho
  print*,"  max rel diff rth = ",maxreldiff_rhotheta

  if(maxvel > tolerance .or. &
     maxreldiff_rho > tolerance .or. &
     maxreldiff_rhotheta > tolerance) then
    print*,"FAIL: hydrostatic balance not maintained!"
    print*,"  maxvel            = ",maxvel," (tol = ",tolerance,")"
    print*,"  maxreldiff_rho    = ",maxreldiff_rho
    print*,"  maxreldiff_rhoth  = ",maxreldiff_rhotheta
    stop 1
  endif

  print*,"PASS: hydrostatic balance maintained to tolerance"

  ! Clean up
  deallocate(rho0,rhotheta0)
  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram ec_euler3d_hydrostaticbalance
