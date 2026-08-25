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

program LinearEuler2D_SpongeDamping
!! Regression test for the spatially-dependent relaxation (sponge / damping
!! layer) source term of the 2D linear Euler model.
!!
!! A Gaussian acoustic pulse is released at the center of a 6x6-element
!! structured mesh with radiation (outflow) boundaries on all four sides. The
!! outermost two rings of elements carry a spatially varying relaxation rate
!! sigma (solution variable 6), ramped quadratically from sigma = 0 at the
!! inner edge of the layer (closer to the interior) up to sigma = sigma_max at
!! the domain boundary. The interior 2x2 block of elements has sigma = 0
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
!!      (variable 6 has zero flux and zero source, and is excluded from the
!!      stepped variables), so it is bitwise unchanged at the final time;
!!  (c) the sponge profile does not change the initial entropy - sigma is not
!!      part of the acoustic energy;
!!  (d) the sponge removes substantially more acoustic energy than radiation
!!      outflow alone, and the residual pressure in the outer quarter of the
!!      damped ring is a small fraction of the undamped run's.

  use self_data
  use self_lineareuler2d
  use self_mesh_2d

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 5
  integer,parameter :: targetDegree = 10
  integer,parameter :: nElemPerDim = 6
  real(prec),parameter :: dx = 0.1_prec
  real(prec),parameter :: dt = 5.0e-4_prec
  real(prec),parameter :: endtime = 0.28_prec
  real(prec),parameter :: iointerval = 0.28_prec
  real(prec),parameter :: rho0 = 1.0_prec
  real(prec),parameter :: c0 = 1.0_prec
  real(prec),parameter :: amp = 1.0e-4_prec
  real(prec),parameter :: Lr = 0.05_prec
  !! Sponge layer is exactly two elements thick, so the outermost two rings of
  !! elements are damped and the interior 2x2 block is not.
  real(prec),parameter :: layer = 2.0_prec*dx
  real(prec),parameter :: sigma_max = 100.0_prec
  !! The sponge run must retain less than this fraction of the undamped run's
  !! final acoustic energy, and of its residual pressure in the outer quarter of
  !! the ring.
  real(prec),parameter :: energy_ratio_tol = 0.25_prec
  real(prec),parameter :: pressure_ratio_tol = 0.1_prec

  type(LinearEuler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)
  integer :: i,j,iel,icase
  real(prec) :: x,y,r,d,sigma_amp,xmax,ymax
  real(prec) :: e0(1:2),ef(1:2),pring(1:2),sigma_drift(1:2)
  real(prec),allocatable :: sigma0(:,:,:)

  ! Radiation (outflow) conditions on all domain boundaries
  bcids(1:4) = SELF_BC_RADIATION

  call mesh%StructuredMesh(nElemPerDim,nElemPerDim,1,1,dx,dx,bcids)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  xmax = real(nElemPerDim,prec)*dx
  ymax = xmax

  allocate(sigma0(1:controlDegree+1,1:controlDegree+1,1:mesh%nElem))

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

    ! Gaussian pressure pulse at the domain center, at rest, in a uniform
    ! medium; the relaxation rate is ramped up inside the outermost ring of
    ! elements only.
    do iel = 1,mesh%nElem
      do j = 1,modelobj%solution%N+1
        do i = 1,modelobj%solution%N+1
          x = modelobj%geometry%x%interior(i,j,iel,1,1)
          y = modelobj%geometry%x%interior(i,j,iel,1,2)
          r = sqrt((x-0.5_prec*xmax)**2+(y-0.5_prec*ymax)**2)

          ! Distance to the nearest domain boundary
          d = min(x,xmax-x,y,ymax-y)

          modelobj%solution%interior(i,j,iel,1) = 0.0_prec ! u
          modelobj%solution%interior(i,j,iel,2) = 0.0_prec ! v
          modelobj%solution%interior(i,j,iel,3) = &
            amp*c0*c0*exp(-log(2.0_prec)*r*r/(Lr*Lr)) ! P
          modelobj%solution%interior(i,j,iel,4) = c0 ! c
          modelobj%solution%interior(i,j,iel,5) = rho0 ! rho0
          if(d < layer) then
            ! Quadratic ramp: zero at the inner edge of the ring, sigma_amp at
            ! the domain boundary.
            modelobj%solution%interior(i,j,iel,6) = &
              sigma_amp*((layer-d)/layer)**2
          else
            modelobj%solution%interior(i,j,iel,6) = 0.0_prec
          endif
        enddo
      enddo
    enddo
    call modelobj%solution%UpdateDevice()

    sigma0 = modelobj%solution%interior(:,:,:,6)

    call modelobj%CalculateEntropy()
    e0(icase) = modelobj%entropy

    call modelobj%SetTimeIntegrator(integrator)
    call modelobj%ForwardStep(tn=endtime,dt=dt,ioInterval=iointerval)
    ef(icase) = modelobj%entropy

    call modelobj%solution%UpdateHost()

    ! Largest deviation of the relaxation rate from its initial value. The
    ! relaxation rate carries zero flux and zero source and is not one of the
    ! stepped variables, so this must be identically zero.
    sigma_drift(icase) = maxval(abs(modelobj%solution%interior(:,:,:,6)-sigma0))

    ! Residual pressure in the outer quarter of the (would-be) damped ring,
    ! nearest the domain boundary, where the relaxation rate exceeds 0.56*sigma_max.
    pring(icase) = 0.0_prec
    do iel = 1,mesh%nElem
      do j = 1,modelobj%solution%N+1
        do i = 1,modelobj%solution%N+1
          x = modelobj%geometry%x%interior(i,j,iel,1,1)
          y = modelobj%geometry%x%interior(i,j,iel,1,2)
          d = min(x,xmax-x,y,ymax-y)
          if(d < 0.25_prec*layer) then
            pring(icase) = max(pring(icase),abs(modelobj%solution%interior(i,j,iel,3)))
          endif
        enddo
      enddo
    enddo

    print*,"case (1=undamped, 2=sponge) :",icase
    print*,"  initial, final entropy    :",e0(icase),ef(icase)
    print*,"  max |p| in the ring       :",pring(icase)
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

  if(pring(2) > pressure_ratio_tol*pring(1)) then
    print*,"Error: sponge layer failed to damp the pressure in the ring. Ratio = ", &
      pring(2)/pring(1)," > tol = ",pressure_ratio_tol
    stop 1
  endif

  print*,"energy ratio (sponge/undamped)   :",ef(2)/ef(1)
  print*,"pressure ratio (sponge/undamped) :",pring(2)/pring(1)

  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram LinearEuler2D_SpongeDamping
