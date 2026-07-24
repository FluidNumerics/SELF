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

program LinearEuler2D_Mortar_SoundWave
!! Regression test for 2:1 nonconforming (mortar) interfaces in the 2D linear Euler
!! model. A spherical acoustic pulse is released inside the big element of the
!! three-element SimpleMortarMesh, centered on the mortar interface, and propagates
!! across it into the two refined elements.
!!
!! Assertions:
!!  (a) the model remains numerically stable across the nonconforming interface :
!!      the acoustic energy (entropy) stays finite and does not exceed its initial
!!      value (the interior upwind flux and radiation boundaries are dissipative);
!!  (b) the solution remains free of NaNs.

  use self_data
  use self_lineareuler2d
  use self_mesh_2d

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 5
  integer,parameter :: targetDegree = 10
  real(prec),parameter :: dx = 0.1_prec
  real(prec),parameter :: dt = 5.0e-4_prec
  real(prec),parameter :: endtime = 0.1_prec
  real(prec),parameter :: iointerval = 0.05_prec
  real(prec),parameter :: c0 = 1.0_prec
  real(prec),parameter :: rho0 = 1.0_prec
  real(prec),parameter :: amp = 1.0e-4_prec
  real(prec),parameter :: Lr = 0.02_prec

  type(LinearEuler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)
  real(prec) :: e0,ef

  ! Radiation (outflow) conditions on all domain boundaries
  bcids(1:4) = [SELF_BC_RADIATION, & ! south
                SELF_BC_RADIATION, & ! east
                SELF_BC_RADIATION, & ! north
                SELF_BC_RADIATION] ! west

  ! Three-element mesh with a single 2:1 mortar interface at x = 2*dx
  call mesh%SimpleMortarMesh(dx,bcids)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%prescribed_bcs_enabled = .false.
  modelobj%tecplot_enabled = .false.
  modelobj%rho0 = rho0

  ! Acoustic pulse centered on the mortar interface
  call modelobj%SphericalSoundWave(amp,Lr,2.0_prec*dx,dx,c0)

  call modelobj%CalculateEntropy()
  e0 = modelobj%entropy

  call modelobj%SetTimeIntegrator(integrator)
  call modelobj%ForwardStep(tn=endtime,dt=dt,ioInterval=iointerval)

  ef = modelobj%entropy

  if(ef /= ef) then
    print*,"Error: entropy is NaN after crossing the mortar interface."
    stop 1
  endif
  if(ef > e0) then
    print*,"Error: entropy grew across the mortar interface :",e0,ef
    stop 1
  endif
  print*,"initial, final entropy :",e0,ef

  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram LinearEuler2D_Mortar_SoundWave
