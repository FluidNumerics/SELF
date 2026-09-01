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
module lineareuler2d_planewave45_model
!! Sinusoidal plane wave propagating at 45 degrees through a 2-D domain, used as a
!! temporal order-of-accuracy benchmark for the low-storage Runge-Kutta integrators.
!!
!! The linearised Euler equations for a uniform background (sound speed c, density rho0)
!!
!!   du/dt = -(1/rho0) dP/dx
!!   dv/dt = -(1/rho0) dP/dy
!!   dP/dt = -rho0 c^2 ( du/dx + dv/dy )
!!
!! admit the exact plane wave solution, for any wave vector k = (kx,ky) with
!! omega = c*|k|,
!!
!!   phase = kx*x + ky*y - omega*t
!!   P(x,y,t) = A sin(phase)
!!   u(x,y,t) = A/(rho0*c) * (kx/|k|) * sin(phase)
!!   v(x,y,t) = A/(rho0*c) * (ky/|k|) * sin(phase)
!!
!! Units: x,y [m], t [s], c [m s^-1], rho0 [kg m^-3], A [kg m^-1 s^-2], k [m^-1],
!! omega [s^-1]. Taking kx = ky sends the wave along the (1,1)/sqrt(2) diagonal, from
!! the lower-left towards the upper-right corner.
!!
!! Every domain boundary is tagged SELF_BC_PRESCRIBED and is handed the exact solution
!! evaluated at the model time. The exterior state is therefore a function of t, which is
!! what makes the Runge-Kutta stage times observable: an integrator that evaluates its
!! stages at the wrong times still converges, but only at first order.

  use self_lineareuler2d
  use SELF_BoundaryConditions
#ifdef ENABLE_GPU
  use iso_c_binding
  use SELF_Constants
#endif

  implicit none

  type,extends(lineareuler2d) :: lineareuler2d_planewave45
    real(prec) :: kx = 0.0_prec !! x wave number [m^-1]
    real(prec) :: ky = 0.0_prec !! y wave number [m^-1]
    real(prec) :: omega = 0.0_prec !! angular frequency, = c*sqrt(kx^2+ky^2) [s^-1]
    real(prec) :: amp = 0.0_prec !! peak pressure amplitude [kg m^-1 s^-2]
    real(prec) :: cref = 1.0_prec !! uniform background sound speed [m s^-1]

  contains

    procedure :: AdditionalInit => AdditionalInit_planewave45
    procedure :: SetExactSolution
    procedure :: MaxError

  endtype lineareuler2d_planewave45

#ifdef ENABLE_GPU
  interface
    subroutine hbc2d_planewave45_lineareuler2d_gpu(extboundary,xboundary, &
                                                   elements,sides,nBoundaries,N,nel, &
                                                   kx,ky,omega,amp,rho0,c,t) &
      bind(c,name="hbc2d_planewave45_lineareuler2d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: extboundary,xboundary,elements,sides
      integer(c_int),value :: nBoundaries,N,nel
      real(c_prec),value :: kx,ky,omega,amp,rho0,c,t
    endsubroutine hbc2d_planewave45_lineareuler2d_gpu
  endinterface
#endif

contains

  subroutine AdditionalInit_planewave45(this)
    !! Register the time-dependent prescribed boundary condition.
    implicit none
    class(lineareuler2d_planewave45),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    ! Register the parent class boundary conditions first. Going through the parent
    ! component (rather than calling AdditionalInit_LinearEuler2D_t directly) keeps the
    ! device-resident BC kernels that the GPU backend registers.
    call this%lineareuler2d%AdditionalInit()

#ifdef ENABLE_GPU
    bcfunc => hbc2d_Prescribed_planewave45_gpu
#else
    bcfunc => hbc2d_Prescribed_planewave45
#endif
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_PRESCRIBED,"prescribed",bcfunc)

  endsubroutine AdditionalInit_planewave45

  pure subroutine planewave_state(kx,ky,omega,amp,rho0,c,x,y,t,s)
    !! Exact plane wave state at (x,y,t), in the model's variable ordering
    !! s = (u, v, P, c, rho0, sigma).
    implicit none
    real(prec),intent(in) :: kx,ky,omega,amp,rho0,c,x,y,t
    real(prec),intent(out) :: s(1:6)
    ! Local
    real(prec) :: phase,sphase,kmag

    kmag = sqrt(kx*kx+ky*ky)
    phase = kx*x+ky*y-omega*t
    sphase = sin(phase)

    s(1) = amp*(kx/kmag)*sphase/(rho0*c) ! u
    s(2) = amp*(ky/kmag)*sphase/(rho0*c) ! v
    s(3) = amp*sphase ! P
    s(4) = c ! sound speed (uniform background)
    s(5) = rho0 ! background density (uniform)
    s(6) = 0.0_prec ! relaxation rate (no sponge)

  endsubroutine planewave_state

  subroutine SetExactSolution(this,t)
    !! Fill the interior solution with the exact plane wave at time t.
    implicit none
    class(lineareuler2d_planewave45),intent(inout) :: this
    real(prec),intent(in) :: t
    ! Local
    integer :: i,j,iel
    real(prec) :: s(1:6)

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem)

      call planewave_state(this%kx,this%ky,this%omega,this%amp,this%rho0,this%cref, &
                           this%geometry%x%interior(i,j,iel,1,1), &
                           this%geometry%x%interior(i,j,iel,1,2),t,s)

      this%solution%interior(i,j,iel,1:6) = s(1:6)

    enddo

    call this%solution%UpdateDevice()

  endsubroutine SetExactSolution

  function MaxError(this,t) result(err)
    !! Max-norm of the difference between the stepped solution and the exact plane wave
    !! at time t, over the three prognostic variables (u,v,P) and all interior nodes,
    !! normalised by the pressure amplitude. Assumes the host copy of the solution is
    !! current (call solution%UpdateHost() first on GPU builds).
    implicit none
    class(lineareuler2d_planewave45),intent(inout) :: this
    real(prec),intent(in) :: t
    real(prec) :: err
    ! Local
    integer :: i,j,iel,ivar
    real(prec) :: s(1:6)

    err = 0.0_prec
    do iel = 1,this%mesh%nElem
      do j = 1,this%solution%N+1
        do i = 1,this%solution%N+1

          call planewave_state(this%kx,this%ky,this%omega,this%amp,this%rho0,this%cref, &
                               this%geometry%x%interior(i,j,iel,1,1), &
                               this%geometry%x%interior(i,j,iel,1,2),t,s)

          do ivar = 1,3
            err = max(err,abs(this%solution%interior(i,j,iel,ivar)-s(ivar)))
          enddo

        enddo
      enddo
    enddo

    err = err/this%amp

  endfunction MaxError

  subroutine hbc2d_Prescribed_planewave45(bc,mymodel)
    !! Host implementation of the time-dependent prescribed boundary condition.
    implicit none
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel
    ! Local
    integer :: n,i,iEl,j
    real(prec) :: s(1:6)

    select type(m => mymodel)
    class is(lineareuler2d_planewave45)

      do n = 1,bc%nBoundaries
        iEl = bc%elements(n)
        j = bc%sides(n)
        do i = 1,m%solution%interp%N+1

          call planewave_state(m%kx,m%ky,m%omega,m%amp,m%rho0,m%cref, &
                               m%geometry%x%boundary(i,j,iEl,1,1), &
                               m%geometry%x%boundary(i,j,iEl,1,2),m%t,s)

          m%solution%extBoundary(i,j,iEl,1:6) = s(1:6)

        enddo
      enddo

    endselect

  endsubroutine hbc2d_Prescribed_planewave45

#ifdef ENABLE_GPU
  subroutine hbc2d_Prescribed_planewave45_gpu(bc,mymodel)
    !! Device implementation of the time-dependent prescribed boundary condition.
    !! The exterior state lives in solution%extBoundary_gpu; writing the host array
    !! instead would leave the device state untouched.
    implicit none
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    select type(m => mymodel)
    class is(lineareuler2d_planewave45)
      if(bc%nBoundaries > 0) then
        call hbc2d_planewave45_lineareuler2d_gpu( &
          m%solution%extBoundary_gpu, &
          m%geometry%x%boundary_gpu, &
          bc%elements_gpu,bc%sides_gpu, &
          bc%nBoundaries,m%solution%interp%N,m%solution%nElem, &
          m%kx,m%ky,m%omega,m%amp,m%rho0,m%cref,m%t)
      endif
    endselect

  endsubroutine hbc2d_Prescribed_planewave45_gpu
#endif

endmodule lineareuler2d_planewave45_model

program lineareuler2d_planewave45_dtconvergence
!! Temporal order-of-accuracy test for the low-storage Runge-Kutta integrators.
!!
!! A 45-degree sinusoidal plane wave is stepped across a structured 8x8 element mesh of
!! the unit square with SELF_BC_PRESCRIBED on all four sides, the exact solution being
!! supplied at the boundaries at every Runge-Kutta stage. Starting from a marginally
!! stable time step the run is repeated with dt successively halved, and the max-norm
!! error against the exact solution at the final time is recorded for each dt.
!!
!! Because the boundary data is time dependent, the error only falls at the scheme's
!! design order if the integrator evaluates its stages at the correct times. An integrator
!! that applies the stage-time offsets one stage late is still consistent, but violates the
!! first order condition sum(b_i c_i) = 1/2 and converges at first order: all three schemes
!! measure order 1.0 here if the stage time is assigned after the tendency instead of
!! before it.
!!
!! The observed order is measured by a least-squares fit of log(error) against log(dt),
!! restricted to the leading refinement levels on which halving dt still buys the design
!! order. That restriction is what stops the fit from being polluted by the levels at which
!! the (spectrally small) spatial discretisation error takes over; with rk4 the floor is
!! reached inside the sweep, and the affected level is dropped automatically.
!!
!! Which integrator is exercised is selected with the SELF_RK_CONV_INTEGRATOR environment
!! variable ("rk2", "rk3" or "rk4"); the default is "rk3".

  use self_data
  use self_mesh_2d
  use lineareuler2d_planewave45_model

  implicit none

  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 15
  integer,parameter :: nTile = 2 ! tiles in x and y
  real(prec),parameter :: c0 = 1.0_prec ! background sound speed
  real(prec),parameter :: rho0 = 1.0_prec ! background density
  real(prec),parameter :: amp = 1.0_prec*10.0_prec**(-2) ! peak pressure
  real(prec),parameter :: kx = 2.0_prec*pi ! x wave number
  real(prec),parameter :: ky = 2.0_prec*pi ! y wave number
  !! The coarsest level must already be resolved: a run that is unstable, or so far from
  !! the asymptotic regime that the "order" between levels is meaningless, is rejected
  !! outright rather than being allowed to produce a large apparent convergence rate.
  real(prec),parameter :: max_coarse_error = 0.25_prec
  integer,parameter :: minfit = 3 ! minimum levels required for a meaningful fit
  integer,parameter :: maxlev = 8

  character(SELF_INTEGRATOR_LENGTH) :: integrator
  character(len=255) :: envval
  integer :: nsteps0,nlev,k,nfit,nsteps
  integer :: nxPerTile ! elements per tile in x and y
  integer :: bcids(1:4)
  real(prec) :: dx ! element size; the domain is always the unit square
  real(prec) :: endtime
  real(prec) :: omega,target_order,order_tol,slope,ordr,refine_factor
  real(prec) :: err(1:maxlev),dtv(1:maxlev)
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry

  omega = c0*sqrt(kx*kx+ky*ky)

  call get_environment_variable("SELF_RK_CONV_INTEGRATOR",envval)
  if(len_trim(envval) == 0) then
    integrator = 'rk3'
  else
    integrator = trim(adjustl(envval))
  endif

  ! Each scheme gets the cheapest configuration that still leaves four refinement
  ! levels clear of the spatial error floor, because this test runs in the debug and
  ! coverage matrices too, where it costs roughly twenty times what it does here.
  !
  ! nsteps0 sets the coarsest time step, chosen just inside the stability limit for
  ! the mesh in use. Measured on the unit square at N=7 with c=1: at h=0.125 rk2 and
  ! rk3 lose stability between dt = 1.95e-3 and 2.6e-3 and rk4 between 3.9e-3 and
  ! 5.2e-3; the limits scale with h.
  !
  ! rk2 and rk3 run on 4x4 elements, whose floor near 3.6e-8 they stay well above.
  ! Fourth order eats four decades in four halvings, so rk4 needs the lower floor of
  ! an 8x8 mesh; it buys back the cost by integrating over half the time interval,
  ! which scales every error down together and leaves the observed order unchanged.
  select case(trim(integrator))
  case('rk2')
    nxPerTile = 2
    dx = 0.25_prec
    endtime = 0.25_prec
    nsteps0 = 64
    nlev = 4
    target_order = 2.0_prec
    order_tol = 0.4_prec
  case('rk3')
    nxPerTile = 2
    dx = 0.25_prec
    endtime = 0.25_prec
    nsteps0 = 64
    nlev = 4
    target_order = 3.0_prec
    order_tol = 0.4_prec
  case('rk4')
    nxPerTile = 4
    dx = 0.125_prec
    endtime = 0.0625_prec
    nsteps0 = 20
    nlev = 4
    target_order = 4.0_prec
    order_tol = 0.7_prec
  case default
    print*,"Error: unsupported integrator ",trim(integrator)
    stop 1
  endselect

  print*,"45-degree plane wave temporal convergence, integrator = ",trim(integrator)
  print*,"endtime, omega, elements per side, degree :",endtime,omega,nxPerTile*nTile,controlDegree

  ! The mesh, interpolant and geometry are identical at every refinement level, so they are
  ! built once. Freeing the last live mesh also finalizes MPI, which would end the run.
  bcids(1:4) = [SELF_BC_PRESCRIBED, & ! South
                SELF_BC_PRESCRIBED, & ! East
                SELF_BC_PRESCRIBED, & ! North
                SELF_BC_PRESCRIBED] ! West

  call mesh%StructuredMesh(nxPerTile,nxPerTile,nTile,nTile,dx,dx,bcids)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  do k = 1,nlev
    nsteps = nsteps0*2**(k-1)
    dtv(k) = endtime/real(nsteps,prec)
    call run_case(integrator,nsteps,mesh,geometry,err(k))
  enddo

  call geometry%free()
  call interp%free()
  call mesh%free()

  print*,""
  print*,"  level     nsteps               dt            error           order"
  do k = 1,nlev
    if(k == 1) then
      print*,k,nsteps0*2**(k-1),dtv(k),err(k)
    else
      ordr = log(err(k-1)/err(k))/log(2.0_prec)
      print*,k,nsteps0*2**(k-1),dtv(k),err(k),ordr
    endif
  enddo
  print*,""

  if(err(1) > max_coarse_error) then
    print*,"Error: the coarsest time step is not in the asymptotic regime. error =",err(1)
    print*,"       (limit ",max_coarse_error,"). The run is unstable or badly underresolved."
    stop 1
  endif

  do k = 1,nlev
    if(err(k) /= err(k)) then
      print*,"Error: error norm is not a number at level ",k
      stop 1
    endif
  enddo

  ! Keep the levels on which the temporal error still dominates. A level is retained only
  ! while halving dt still buys the order this integrator is supposed to deliver (to within
  ! the same tolerance the final assertion uses), so the fit stops of its own accord at the
  ! dt where the spatial discretisation error takes over. No floor value is hard-coded.
  refine_factor = 2.0_prec**(-(target_order-order_tol))
  print*,"Retaining levels while the error falls by at least a factor of",refine_factor

  nfit = 1
  do k = 2,nlev
    if(err(k) <= refine_factor*err(k-1)) then
      nfit = k
    else
      exit
    endif
  enddo

  print*,"Levels above the spatial error floor :",nfit

  ! Below the floor the error stops responding to dt and wanders by a few parts in a
  ! thousand, so monotonicity is only meaningful on the levels being fitted - where the
  ! retention test above already enforces it. What the remaining levels must show is
  ! that refining further settles onto the floor rather than diverging from it.
  do k = nfit+1,nlev
    if(err(k) > 2.0_prec*err(nfit)) then
      print*,"Error: refining past the spatial error floor made the error grow at level ", &
        k,err(nfit),err(k)
      stop 1
    endif
  enddo

  if(nfit < minfit) then
    print*,"Error: only ",nfit," consecutive refinement level(s) reached order ", &
      target_order-order_tol,"; ",minfit," are required."
    print*,"       Slope over the whole sweep :",fit_slope(dtv(1:nlev),err(1:nlev),nlev)
    print*,"       Either the integrator is not converging at its design rate (check that"
    print*,"       the Runge-Kutta stage times are consistent with the tableau), or the"
    print*,"       spatial error floor was reached before the order could be measured."
    stop 1
  endif

  slope = fit_slope(dtv(1:nfit),err(1:nfit),nfit)
  print*,"Observed temporal order (least-squares) :",slope
  print*,"Expected                                :",target_order," (tolerance ",order_tol,")"

  if(slope < target_order-order_tol) then
    print*,"Error: observed temporal order ",slope," is below ",target_order-order_tol
    print*,"       Runge-Kutta stage times are inconsistent with the tableau."
    stop 1
  endif

  print*,"PASS"

contains

  subroutine run_case(integrator,nsteps,mesh,geometry,err)
    !! Step the plane wave to endtime with dt = endtime/nsteps and return the
    !! amplitude-normalised max-norm error against the exact solution.
    implicit none
    character(SELF_INTEGRATOR_LENGTH),intent(in) :: integrator
    integer,intent(in) :: nsteps
    type(Mesh2D),target,intent(inout) :: mesh
    type(SEMQuad),target,intent(inout) :: geometry
    real(prec),intent(out) :: err
    ! Local
    type(lineareuler2d_planewave45) :: modelobj
    real(prec) :: dt

    dt = endtime/real(nsteps,prec)

    call modelobj%Init(mesh,geometry)
    modelobj%tecplot_enabled = .false.
    modelobj%rho0 = rho0
    modelobj%cref = c0
    modelobj%kx = kx
    modelobj%ky = ky
    modelobj%omega = omega
    modelobj%amp = amp

    call modelobj%SetExactSolution(0.0_prec)
    call modelobj%SetTimeIntegrator(integrator)

    ! A single IO interval keeps the run to one file write.
    call modelobj%ForwardStep(endtime,dt,endtime)

    call modelobj%solution%UpdateHost()
    err = modelobj%MaxError(endtime)

    call modelobj%free()

  endsubroutine run_case

  function fit_slope(dtv,err,n) result(slope)
    !! Least-squares slope of log(err) against log(dt), i.e. the observed order.
    implicit none
    integer,intent(in) :: n
    real(prec),intent(in) :: dtv(1:n)
    real(prec),intent(in) :: err(1:n)
    real(prec) :: slope
    ! Local
    integer :: k
    real(prec) :: x,y,sx,sy,sxx,sxy,rn

    sx = 0.0_prec
    sy = 0.0_prec
    sxx = 0.0_prec
    sxy = 0.0_prec
    do k = 1,n
      x = log(dtv(k))
      y = log(err(k))
      sx = sx+x
      sy = sy+y
      sxx = sxx+x*x
      sxy = sxy+x*y
    enddo
    rn = real(n,prec)
    slope = (rn*sxy-sx*sy)/(rn*sxx-sx*sx)

  endfunction fit_slope

endprogram lineareuler2d_planewave45_dtconvergence
