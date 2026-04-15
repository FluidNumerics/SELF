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

program ec_euler3d_thermal_bubble
  !! Rising warm bubble benchmark for the EC-DG compressible Euler model
  !! with potential temperature formulation.
  !!
  !! A pressure-balanced warm bubble (cos^2 profile) is placed in a
  !! hydrostatically balanced atmosphere with uniform theta0 = 300 K.
  !! The density deficit drives buoyant rising motion.
  !!
  !! Domain: [0, 1000m] x [0, 1000m] x [0, 1500m]
  !! Bubble center: (500, 500, 350) m
  !! Bubble radius: 250 m
  !! Perturbation: dtheta = 2.0 K
  !!
  !! Reference: Robert (1993), J. Comput. Phys.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_ECEuler3D

  implicit none

  ! Physical parameters
  real(prec),parameter :: theta0 = 300.0_prec ! Background pot. temp. [K]
  real(prec),parameter :: dtheta = 2.0_prec ! Bubble perturbation [K]
  real(prec),parameter :: bubble_r0 = 250.0_prec ! Bubble radius [m]
  real(prec),parameter :: bubble_x0 = 500.0_prec ! Bubble center x [m]
  real(prec),parameter :: bubble_y0 = 500.0_prec ! Bubble center y [m]
  real(prec),parameter :: bubble_z0 = 350.0_prec ! Bubble center z [m]

  ! Grid parameters
  integer,parameter :: nxPerTile = 4 ! Elements in x per tile
  integer,parameter :: nyPerTile = 4 ! Elements in y per tile
  integer,parameter :: nzPerTile = 6 ! Elements in z per tile
  real(prec),parameter :: dx = 250.0_prec ! Element size x [m]
  real(prec),parameter :: dy = 250.0_prec ! Element size y [m]
  real(prec),parameter :: dz = 250.0_prec ! Element size z [m]

  ! Time integration parameters
  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  real(prec),parameter :: dt = 1.0_prec*10.0_prec**(-2)
  real(prec),parameter :: endtime = 1.0_prec*10.0_prec**(-1)
  real(prec),parameter :: iointerval = endtime

  type(ECEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)
  real(prec) :: e0,ef

  ! Create interpolant (Gauss-Lobatto for EC split form)
  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS_LOBATTO, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  ! Create structured mesh
  ! Domain: [0, 1000] x [0, 1000] x [0, 1500]
  bcids(1:6) = [SELF_BC_NONORMALFLOW, & ! Bottom
                SELF_BC_NONORMALFLOW, & ! South
                SELF_BC_NONORMALFLOW, & ! East
                SELF_BC_NONORMALFLOW, & ! North
                SELF_BC_NONORMALFLOW, & ! West
                SELF_BC_NONORMALFLOW] ! Top
  call mesh%StructuredMesh(nxPerTile,nyPerTile,nzPerTile, &
                           1,1,1,dx,dy,dz,bcids)

  ! Generate geometry
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  ! Initialize model
  call modelobj%Init(mesh,geometry)
  modelobj%tecplot_enabled = .false.

  ! Set up hydrostatically balanced background
  call modelobj%SetHydrostaticBalance(theta0)

  ! Add warm bubble perturbation (pressure-balanced)
  call modelobj%AddThermalBubble(dtheta,bubble_r0, &
                                 bubble_x0,bubble_y0,bubble_z0)

  ! Write initial condition
  call modelobj%WriteModel()
  call modelobj%IncrementIOCounter()

  ! Report initial entropy
  call modelobj%CalculateEntropy()
  call modelobj%ReportEntropy()
  e0 = modelobj%entropy

  ! Time integrate
  call modelobj%SetTimeIntegrator(integrator)
  call modelobj%ForwardStep(endtime,dt,iointerval)

  ! Check final entropy
  ef = modelobj%entropy
  if(ef /= ef) then
    print*,"Error: Final entropy is inf or nan",ef
    stop 1
  endif

  print*,"Initial entropy: ",e0
  print*,"Final entropy:   ",ef

  ! Clean up
  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram ec_euler3d_thermal_bubble
