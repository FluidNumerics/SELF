! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2026 Fluid Numerics LLC
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

program LinearEuler2D_PML_Absorption
!! Regression test for the Hu (2001) unsplit PML on the 2D linear
!! Euler model. A right-going Gaussian acoustic pulse (uniform in y)
!! is released near the western interior; it propagates as a planar
!! wave into an eastern PML region and is absorbed. The uniform-in-y
!! initial condition keeps the pulse planar so the only mechanism
!! reducing the interior amplitude is the PML absorption (not 2D
!! cylindrical spreading).
!!
!! Assertions:
!!  (a) the model remains numerically stable: entropy stays finite
!!      and does not exceed the initial value;
!!  (b) by the final time, the acoustic perturbation in the interior
!!      region (x <= x_interior_max) is below 5% of its initial peak,
!!      i.e. the PML absorbed >= 95% of the energy that crossed
!!      into it.

  use self_data
  use self_lineareuler2d_pml
  use self_mesh_2d

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 5
  integer,parameter :: targetDegree = 10
  real(prec),parameter :: dt = 5.0e-3_prec
  real(prec),parameter :: endtime = 5.0_prec
  real(prec),parameter :: iointerval = 1.0_prec
  real(prec),parameter :: c0 = 1.0_prec
  real(prec),parameter :: amp = 1.0e-3_prec
  real(prec),parameter :: pulse_width = 0.5_prec ! Gaussian halfwidth in x
  real(prec),parameter :: pulse_x0 = 1.5_prec
  real(prec),parameter :: x_interior_max = 4.0_prec
  real(prec),parameter :: pml_width = 2.0_prec
  real(prec),parameter :: sigma_max = 20.0_prec
  real(prec),parameter :: absorption_tol = 5.0e-2_prec

  type(LinearEuler2D_PML) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)
  integer :: i,j,iel
  real(prec) :: x,r2,p_loc,xmid
  real(prec) :: e0,ef,p_max_initial,p_max_interior

  ! Thin (in y) domain: 12 elements in x by 2 elements in y. The
  ! initial condition is uniform in y so the no-normal-flow walls
  ! on north/south sides preserve a planar wave.
  bcids(1:4) = [SELF_BC_NONORMALFLOW, & ! south
                SELF_BC_NONORMALFLOW, & ! east (behind PML)
                SELF_BC_NONORMALFLOW, & ! north
                SELF_BC_NONORMALFLOW] ! west
  call mesh%StructuredMesh(12,2,1,1,0.5_prec,0.5_prec,bcids)

  ! Tag the eastern strip (x > x_interior_max) as PML material.
  if(allocated(mesh%materialNames)) deallocate(mesh%materialNames)
  mesh%nMaterials = 2
  allocate(mesh%materialNames(1:2))
  mesh%materialNames(1) = "interior"
  mesh%materialNames(2) = "pml"
  do iel = 1,mesh%nElem
    xmid = sum(mesh%nodeCoords(1,:,:,iel))/real(size(mesh%nodeCoords,2)*size(mesh%nodeCoords,3),prec)
    if(xmid > x_interior_max) then
      mesh%elemMaterial(iel) = 2
    else
      mesh%elemMaterial(iel) = 1
    endif
  enddo

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%prescribed_bcs_enabled = .false.
  modelobj%tecplot_enabled = .false.
  modelobj%rho0 = 1.0_prec

  ! Cubic ramp from 0 to sigma_max across the PML thickness.
  call modelobj%SetPMLProfile(x_interior_min=0.0_prec, &
                              x_interior_max=x_interior_max, &
                              y_interior_min=0.0_prec, &
                              y_interior_max=1.0_prec, &
                              pml_width=pml_width, &
                              sigma_max=sigma_max, &
                              ramp_exponent=3)

  ! Initial planar (uniform in y) right-going acoustic pulse:
  !   p(x) = amp * exp(-((x-pulse_x0)/pulse_width)^2)
  !   u = p/(rho0*c),  v = 0
  ! which is the right eigenvector of A for the linear Euler system.
  do iel = 1,mesh%nElem
    do j = 1,modelobj%solution%N+1
      do i = 1,modelobj%solution%N+1
        x = modelobj%geometry%x%interior(i,j,iel,1,1)
        r2 = (x-pulse_x0)**2
        p_loc = amp*exp(-r2/(pulse_width*pulse_width))
        modelobj%solution%interior(i,j,iel,1) = p_loc/(modelobj%rho0*c0) ! u
        modelobj%solution%interior(i,j,iel,2) = 0.0_prec ! v
        modelobj%solution%interior(i,j,iel,3) = p_loc ! p
        modelobj%solution%interior(i,j,iel,4) = c0 ! c
        modelobj%solution%interior(i,j,iel,5:7) = 0.0_prec ! phi_*
      enddo
    enddo
  enddo
  call modelobj%solution%UpdateDevice()

  p_max_initial = maxval(modelobj%solution%interior(:,:,:,3))
  print*,"Initial max |p|: ",p_max_initial
  call modelobj%CalculateEntropy()
  e0 = modelobj%entropy
  print*,"Initial entropy: ",e0

  call modelobj%SetTimeIntegrator(integrator)
  call modelobj%ForwardStep(endtime,dt,iointerval)
  ef = modelobj%entropy
  print*,"Final entropy:   ",ef

  call modelobj%solution%UpdateHost()
  p_max_interior = 0.0_prec
  do iel = 1,mesh%nElem
    xmid = sum(mesh%nodeCoords(1,:,:,iel))/real(size(mesh%nodeCoords,2)*size(mesh%nodeCoords,3),prec)
    if(xmid > x_interior_max) cycle
    do j = 1,modelobj%solution%N+1
      do i = 1,modelobj%solution%N+1
        p_max_interior = max(p_max_interior,abs(modelobj%solution%interior(i,j,iel,3)))
      enddo
    enddo
  enddo
  print*,"Max |p| in interior at final time: ",p_max_interior
  print*,"Residual ratio:                    ",p_max_interior/p_max_initial

  if(ef /= ef) then
    print*,"Error: final entropy is NaN ",ef
    stop 1
  endif
  if(ef > e0) then
    print*,"Error: final entropy exceeds initial entropy (PML should dissipate) ",e0,ef
    stop 1
  endif
  if(p_max_interior > absorption_tol*p_max_initial) then
    print*,"Error: PML failed to absorb the pulse. Residual = ", &
      p_max_interior/p_max_initial," > tol = ",absorption_tol
    stop 1
  endif

  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram LinearEuler2D_PML_Absorption
