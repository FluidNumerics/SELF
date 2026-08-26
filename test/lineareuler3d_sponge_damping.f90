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

program LinearEuler3D_SpongeDamping
!! Regression test for the spatially-dependent relaxation (sponge / damping
!! layer) source term of the 3D linear Euler model.
!!
!! A Gaussian acoustic pulse is released at the center of a 4x4x4-element
!! structured mesh with radiation (outflow) boundaries on all six sides. The
!! outermost shell of elements carries a spatially varying relaxation rate
!! sigma (solution variable 7), ramped quadratically from sigma = 0 at the
!! inner edge of the shell (closer to the interior) up to sigma = sigma_max at
!! the domain boundary. The 2x2x2 block of interior elements has sigma = 0
!! identically, so the interior dynamics are the undamped linear Euler
!! dynamics.
!!
!! The final time is chosen so that the acoustic front has just reached the
!! outer boundary: in the undamped run essentially none of the pulse has left
!! the domain yet, so the comparison measures the sponge absorption rather
!! than radiation outflow.
!!
!! The same problem is integrated twice: case 1 with sigma = 0 everywhere (the
!! sigma = 0 regression path, which must reproduce the undamped model exactly)
!! and case 2 with the sponge profile.
!!
!! Assertions:
!!  (a) both runs remain numerically stable: the acoustic energy (entropy) is
!!      finite and does not exceed its initial value;
!!  (b) the relaxation rate is not itself modified by the time integration
!!      (variable 7 has zero flux and zero source, and is excluded from the
!!      stepped variables), so it is bitwise unchanged at the final time;
!!  (c) the sponge profile does not change the initial entropy - sigma is not
!!      part of the acoustic energy;
!!  (d) the sponge removes substantially more acoustic energy than radiation
!!      outflow alone, and the residual pressure in the outer quarter of the
!!      damped shell is a small fraction of the undamped run's.

  use self_data
  use self_lineareuler3d
  use self_mesh_3d

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 3
  integer,parameter :: targetDegree = 6
  integer,parameter :: nElemPerDim = 4
  real(prec),parameter :: dx = 0.1_prec
  real(prec),parameter :: dt = 1.0e-3_prec
  real(prec),parameter :: endtime = 0.19_prec
  real(prec),parameter :: iointerval = 0.19_prec
  real(prec),parameter :: rho0 = 1.0_prec
  real(prec),parameter :: c0 = 1.0_prec
  real(prec),parameter :: amp = 1.0e-4_prec
  real(prec),parameter :: Lr = 0.05_prec
  !! Sponge layer is exactly one element thick, so only the outermost shell of
  !! elements is damped.
  real(prec),parameter :: layer = dx
  real(prec),parameter :: sigma_max = 300.0_prec
  !! The sponge run must retain less than this fraction of the undamped run's
  !! final acoustic energy, and of its residual pressure in the outer quarter of
  !! the shell. The energy bound is the looser of the two: a pressure-only
  !! initial condition splits into an outgoing and an ingoing wave, and by the
  !! final time the ingoing half is still in the undamped interior, where the
  !! sponge cannot act on it.
  real(prec),parameter :: energy_ratio_tol = 0.5_prec
  real(prec),parameter :: pressure_ratio_tol = 0.1_prec

  type(LinearEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)
  integer :: i,j,k,iel,icase
  real(prec) :: x,y,z,r,d,sigma_amp,xmax
  real(prec) :: e0(1:2),ef(1:2),pshell(1:2),sigma_drift(1:2)
  real(prec),allocatable :: sigma0(:,:,:,:)

  ! Radiation (outflow) conditions on all domain boundaries
  bcids(1:6) = SELF_BC_RADIATION

  call mesh%StructuredMesh(nElemPerDim,nElemPerDim,nElemPerDim,1,1,1,dx,dx,dx,bcids)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  xmax = real(nElemPerDim,prec)*dx

  allocate(sigma0(1:controlDegree+1,1:controlDegree+1,1:controlDegree+1,1:mesh%nElem))

  do icase = 1,2

    if(icase == 1) then
      sigma_amp = 0.0_prec ! undamped reference
    else
      sigma_amp = sigma_max ! sponge layer
    endif

    call modelobj%Init(mesh,geometry)
    modelobj%prescribed_bcs_enabled = .false.
    modelobj%tecplot_enabled = .false.
    modelobj%rho0 = rho0
    modelobj%c = c0

    ! Gaussian pressure pulse at the domain center, at rest, in a uniform
    ! medium; the relaxation rate is ramped up inside the outermost shell of
    ! elements only.
    do iel = 1,mesh%nElem
      do k = 1,modelobj%solution%N+1
        do j = 1,modelobj%solution%N+1
          do i = 1,modelobj%solution%N+1
            x = modelobj%geometry%x%interior(i,j,k,iel,1,1)
            y = modelobj%geometry%x%interior(i,j,k,iel,1,2)
            z = modelobj%geometry%x%interior(i,j,k,iel,1,3)
            r = sqrt((x-0.5_prec*xmax)**2+(y-0.5_prec*xmax)**2+(z-0.5_prec*xmax)**2)

            ! Distance to the nearest domain boundary
            d = min(x,xmax-x,y,xmax-y,z,xmax-z)

            modelobj%solution%interior(i,j,k,iel,1) = 0.0_prec ! u
            modelobj%solution%interior(i,j,k,iel,2) = 0.0_prec ! v
            modelobj%solution%interior(i,j,k,iel,3) = 0.0_prec ! w
            modelobj%solution%interior(i,j,k,iel,4) = &
              amp*c0*c0*exp(-log(2.0_prec)*r*r/(Lr*Lr)) ! P
            modelobj%solution%interior(i,j,k,iel,5) = c0 ! c
            modelobj%solution%interior(i,j,k,iel,6) = rho0 ! rho0
            if(d < layer) then
              ! Quadratic ramp: zero at the inner edge of the shell, sigma_amp
              ! at the domain boundary.
              modelobj%solution%interior(i,j,k,iel,7) = &
                sigma_amp*((layer-d)/layer)**2
            else
              modelobj%solution%interior(i,j,k,iel,7) = 0.0_prec
            endif
          enddo
        enddo
      enddo
    enddo
    call modelobj%solution%UpdateDevice()

    sigma0 = modelobj%solution%interior(:,:,:,:,7)

    call modelobj%CalculateEntropy()
    e0(icase) = modelobj%entropy

    call modelobj%SetTimeIntegrator(integrator)
    call modelobj%ForwardStep(tn=endtime,dt=dt,ioInterval=iointerval)
    ef(icase) = modelobj%entropy

    call modelobj%solution%UpdateHost()

    ! Largest deviation of the relaxation rate from its initial value. The
    ! relaxation rate carries zero flux and zero source and is not one of the
    ! stepped variables, so this must be identically zero.
    sigma_drift(icase) = maxval(abs(modelobj%solution%interior(:,:,:,:,7)-sigma0))

    ! Residual pressure in the outer quarter of the (would-be) damped shell,
    ! nearest the domain boundary, where the relaxation rate exceeds 0.56*sigma_max.
    pshell(icase) = 0.0_prec
    do iel = 1,mesh%nElem
      do k = 1,modelobj%solution%N+1
        do j = 1,modelobj%solution%N+1
          do i = 1,modelobj%solution%N+1
            x = modelobj%geometry%x%interior(i,j,k,iel,1,1)
            y = modelobj%geometry%x%interior(i,j,k,iel,1,2)
            z = modelobj%geometry%x%interior(i,j,k,iel,1,3)
            d = min(x,xmax-x,y,xmax-y,z,xmax-z)
            if(d < 0.25_prec*layer) then
              pshell(icase) = max(pshell(icase), &
                                  abs(modelobj%solution%interior(i,j,k,iel,4)))
            endif
          enddo
        enddo
      enddo
    enddo

    print*,"case (1=undamped, 2=sponge) :",icase
    print*,"  initial, final entropy    :",e0(icase),ef(icase)
    print*,"  max |p| in the shell      :",pshell(icase)
    print*,"  max sigma drift           :",sigma_drift(icase)

    call modelobj%Free()

  enddo

  deallocate(sigma0)

  do icase = 1,2
    if(ef(icase) /= ef(icase)) then
      print*,"Error: entropy is NaN for case ",icase
      stop 1
    endif
    if(ef(icase) > e0(icase)) then
      print*,"Error: entropy grew for case ",icase,e0(icase),ef(icase)
      stop 1
    endif
    if(sigma_drift(icase) /= 0.0_prec) then
      print*,"Error: the relaxation rate was modified by time stepping, case ", &
        icase,sigma_drift(icase)
      stop 1
    endif
  enddo

  ! The sponge profile is not part of the acoustic energy, so both cases must
  ! start from the same entropy.
  if(e0(1) /= e0(2)) then
    print*,"Error: the sponge profile changed the initial entropy ",e0(1),e0(2)
    stop 1
  endif

  if(ef(2) > energy_ratio_tol*ef(1)) then
    print*,"Error: sponge layer failed to damp the acoustic energy. Ratio = ", &
      ef(2)/ef(1)," > tol = ",energy_ratio_tol
    stop 1
  endif

  if(pshell(2) > pressure_ratio_tol*pshell(1)) then
    print*,"Error: sponge layer failed to damp the pressure in the shell. Ratio = ", &
      pshell(2)/pshell(1)," > tol = ",pressure_ratio_tol
    stop 1
  endif

  print*,"energy ratio (sponge/undamped)   :",ef(2)/ef(1)
  print*,"pressure ratio (sponge/undamped) :",pshell(2)/pshell(1)

  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram LinearEuler3D_SpongeDamping
