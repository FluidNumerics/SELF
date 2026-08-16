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

program LinearEuler3D_Mortar_SoundWave
!! Regression test for 2:1 nonconforming (mortar) interfaces in the 3D linear Euler
!! model. A spherical acoustic pulse is released inside the big element of the
!! five-element SimpleMortarMesh, centered on the mortar interface, and propagates
!! across it into the four refined elements. The run is repeated on a mesh whose
!! small elements are rotated (mortar flips 1,4,6,3), exercising the flip-aware
!! staging in both the solution exchange and the conservative flux collection under
!! live dynamics.
!!
!! Assertions, for both meshes:
!!  (a) the model remains numerically stable across the nonconforming interface :
!!      the acoustic energy (entropy) stays finite and does not exceed its initial
!!      value (the interior upwind flux and radiation boundaries are dissipative);
!!  (b) the solution remains free of NaNs.

  use self_data
  use self_lineareuler3d
  use self_mesh_3d

  implicit none

  integer :: flipsA(1:4)
  integer :: keepAliveBcids(1:6)
  type(Mesh3D),target :: keepAlive
  integer :: r

  ! SELF tears down MPI when the number of live meshes drops to zero (see
  ! test/domaindecomposition_two_meshes.f90). This test builds two meshes in
  ! sequence, so an extra mesh is held open across both of them.
  keepAliveBcids(1:6) = SELF_BC_RADIATION
  call keepAlive%StructuredMesh(1,1,1,1,1,1,0.5_prec,0.5_prec,0.5_prec, &
                                keepAliveBcids)

  r = 0

  ! Axis-aligned : all mortar flips are 0
  r = r+RunSoundWave()

  ! Rotated small elements : mortar flips 1,4,6,3
  flipsA = [1,4,6,3]
  r = r+RunSoundWave(flipsA)

  call keepAlive%Free()

  if(r /= 0) then
    stop 1
  endif

contains

  integer function RunSoundWave(flips) result(r)
    implicit none
    integer,intent(in),optional :: flips(1:4)

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

    r = 0

    ! Radiation (outflow) conditions on all domain boundaries
    bcids(1:6) = SELF_BC_RADIATION

    ! Five-element mesh with a single 2:1 mortar interface at x = 2*dx
    if(present(flips)) then
      call mesh%SimpleMortarMesh(dx,bcids,flips=flips)
    else
      call mesh%SimpleMortarMesh(dx,bcids)
    endif

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
    call modelobj%SphericalSoundWave(amp,Lr,2.0_prec*dx,dx,dx)

    call modelobj%CalculateEntropy()
    e0 = modelobj%entropy

    call modelobj%SetTimeIntegrator(integrator)
    call modelobj%ForwardStep(tn=endtime,dt=dt,ioInterval=iointerval)

    ef = modelobj%entropy

    if(ef /= ef) then
      print*,"Error: entropy is NaN after crossing the mortar interface."
      r = 1
    elseif(ef > e0) then
      print*,"Error: entropy grew across the mortar interface :",e0,ef
      r = 1
    endif
    print*,"initial, final entropy :",e0,ef

    call modelobj%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()

  endfunction RunSoundWave

endprogram LinearEuler3D_Mortar_SoundWave
