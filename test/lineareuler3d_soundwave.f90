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

program LinearEuler3D_SoundWave
!! Conforming-mesh regression baseline for the 3D linear Euler model. A spherical
!! acoustic pulse is released at the center of a 2x2x2-element structured mesh with
!! radiation (outflow) boundaries on all sides.
!!
!! Assertions:
!!  (a) the model remains numerically stable : the acoustic energy (entropy) stays
!!      finite and does not exceed its initial value (the interior upwind flux and
!!      radiation boundaries are dissipative);
!!  (b) the solution remains free of NaNs.
!!
!! This test pins the conforming-mesh behavior that the mortar and AMR regression
!! tests (lineareuler3d_mortar_soundwave, and the 3D AMR suite built on it) treat as
!! their baseline.

  use self_data
  use self_lineareuler3d
  use self_mesh_3d

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 5
  integer,parameter :: targetDegree = 10
  real(prec),parameter :: dx = 0.1_prec
  real(prec),parameter :: dt = 5.0e-4_prec
  real(prec),parameter :: endtime = 0.05_prec
  real(prec),parameter :: iointerval = 0.05_prec
  real(prec),parameter :: rho0 = 1.0_prec
  real(prec),parameter :: amp = 1.0e-4_prec
  real(prec),parameter :: Lr = 0.02_prec

  type(LinearEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)
  real(prec) :: e0,ef

  ! Radiation (outflow) conditions on all domain boundaries
  bcids(1:6) = SELF_BC_RADIATION

  call mesh%StructuredMesh(2,2,2,1,1,1,dx,dx,dx,bcids)

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

  ! Acoustic pulse at the domain center
  call modelobj%SphericalSoundWave(amp,Lr,dx,dx,dx)

  call modelobj%CalculateEntropy()
  e0 = modelobj%entropy

  call modelobj%SetTimeIntegrator(integrator)
  call modelobj%ForwardStep(tn=endtime,dt=dt,ioInterval=iointerval)

  ef = modelobj%entropy

  if(ef /= ef) then
    print*,"Error: entropy is NaN."
    stop 1
  endif
  if(ef > e0) then
    print*,"Error: entropy grew on the conforming mesh :",e0,ef
    stop 1
  endif
  print*,"initial, final entropy :",e0,ef

  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram LinearEuler3D_SoundWave
