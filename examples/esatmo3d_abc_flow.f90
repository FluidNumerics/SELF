program esatmo3d_abc_flow
  !! Arnold-Beltrami-Childress (ABC) flow benchmark for the entropy-conserving
  !! DG compressible Euler / potential-temperature model (ESAtmo3D) on a triply
  !! periodic box T^3 = [0, 2*pi]^3.
  !!
  !! Background / motivation
  !! -----------------------
  !! The ABC velocity field
  !!
  !!   u = A sin(z) + C cos(y)
  !!   v = B sin(x) + A cos(z)
  !!   w = C sin(y) + B cos(x)
  !!
  !! with the canonical, symmetry-breaking parameters A = 1, B = sqrt(2/3),
  !! C = sqrt(1/3) is a Beltrami field (vorticity = velocity, omega = u) and
  !! hence an exact steady solution of the *incompressible* Euler equations: the
  !! nonlinear term (u.grad)u = grad(|u|^2/2) - u x omega reduces to a pure
  !! gradient because u x omega = 0 identically.
  !!
  !! Here we initialise that velocity field in the *compressible* potential-
  !! temperature model with, as requested:
  !!   - constant density            rho   = rho0
  !!   - constant potential temp.    theta = theta0  (so rho*theta = rho0*theta0)
  !!   - constant geopotential       Phi   = phi0
  !!   - no gravity                  g     = 0
  !!
  !! With constant rho and theta the equation of state gives a constant pressure,
  !! so there is no pressure gradient to balance grad(|u|^2/2). The ABC field is
  !! therefore only an *approximate* steady state of the compressible system. The
  !! approximation is excellent in the low-Mach regime realised here: the sound
  !! speed c = sqrt(gamma * Rd * theta0) ~ 347 m/s while |u| = O(1), giving a Mach
  !! number ~ 0.005. Departures from the initial state are dominated by numerical
  !! error and weak acoustic adjustment, exactly the quantities the diagnostics
  !! below are designed to expose.
  !!
  !! On the triply periodic box the dominant large-scale ABC instability band
  !! (|k| < 1) is absent from the discrete spectrum (the smallest representable
  !! wavenumber is |k| = 1), so a well-behaved solver should preserve the ABC
  !! initial condition for many advective turnover times ("geometric protection").
  !!
  !! Diagnostics (computed on the collocation grid, single-rank)
  !! -----------------------------------------------------------
  !!   E_err(t) = ||u(.,t) - u_ABC||_2 / ||u_ABC||_2          (pointwise error)
  !!   dE/E0    = (E(t)-E(0))/E(0),   E = (1/2) integral |u|^2 dV   (energy drift)
  !!   dH/H0    = (H(t)-H(0))/H(0),   H = integral u.omega dV       (helicity drift)
  !!   R(t)     = || u x omega ||_inf                          (Beltrami residual)
  !!
  !! For a Beltrami field R == 0 identically, so R(t) is the sharpest fingerprint
  !! of aliasing / loss of the Beltrami structure. Energy and helicity are ideal
  !! invariants of the Euler equations and should drift only slowly.
  !!
  !! Reference: "ABC Flow as a Benchmark for Spectral-Element Euler Solvers"
  !! (internal notes, 2026); Arnold (1965); Galloway & Frisch (1986).

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_ESAtmo3D

  implicit none

  ! ABC parameters (canonical, symmetry-breaking choice)
  real(prec),parameter :: A_abc = 1.0_prec
  real(prec),parameter :: B_abc = sqrt(2.0_prec/3.0_prec)
  real(prec),parameter :: C_abc = sqrt(1.0_prec/3.0_prec)

  ! Constant thermodynamic background
  real(prec),parameter :: rho0 = 1.0_prec ! constant density [kg/m^3]
  real(prec),parameter :: theta0 = 300.0_prec ! constant potential temperature [K]
  real(prec),parameter :: phi0 = 0.0_prec ! constant geopotential [m^2/s^2]

  ! Grid: triply periodic box [0, 2*pi]^3
  integer,parameter :: nElemPerDir = 4 ! elements per direction
  integer,parameter :: controlDegree = 7 ! polynomial degree (Gauss-Lobatto)
  integer,parameter :: targetDegree = 14 ! plotting/interpolation degree
  real(prec),parameter :: Lbox = 2.0_prec*pi ! domain length per direction [m]
  real(prec),parameter :: de = Lbox/real(nElemPerDir,prec) ! element width [m]

  ! Time integration. The acoustic CFL (set by c ~ 347 m/s) is restrictive while
  ! the advective turnover time (~ 2*pi/|u| ~ 4 s) is much longer; the explicit
  ! step is therefore acoustic-limited. The defaults below integrate only a
  ! short interval (10 steps) so the example runs quickly, as with the other
  ! examples in this directory. For a genuine preservation study, build with an
  ! optimized configuration and increase nReports (and/or iointerval) to cover
  ! several advective turnover times; the diagnostics below are what you watch.
  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  real(prec),parameter :: dt = 5.0_prec*10.0_prec**(-5) ! time step [s]
  real(prec),parameter :: iointerval = 1.0_prec*10.0_prec**(-4) ! diagnostic interval [s]
  integer,parameter :: nReports = 5 ! number of diagnostic intervals (2 steps each)

  type(ESAtmo3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  ! Velocity (u,v,w) on the collocation grid, host-resident. Diagnostics derive
  ! vorticity from this array using the interpolant derivative matrix and the
  ! (constant) geometry metric terms, so the whole diagnostics path is plain
  ! host Fortran that compiles and runs identically on CPU and GPU builds.
  real(prec),allocatable :: uvw(:,:,:,:,:) ! (i,j,k,iEl,1:3) -> (u,v,w)
  real(prec) :: e0,h0 ! baseline invariants E(0), H(0)
  real(prec) :: errL2,Eabs,Habs,Rinf
  integer :: iout

  ! Create the interpolant (Gauss-Lobatto control nodes for the EC split form)
  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS_LOBATTO, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  ! Create a fully triply-periodic structured mesh on [0, 2*pi]^3.
  ! Domain-boundary faces are wired to the element on the opposite side, so the
  ! ABC field sees a genuine periodic box (no walls).
  call mesh%PeriodicStructuredMesh(nxPerTile=nElemPerDir, &
                                   nyPerTile=nElemPerDir, &
                                   nzPerTile=nElemPerDir, &
                                   nTileX=1,nTileY=1,nTileZ=1, &
                                   dx=de,dy=de,dz=de)

  ! Generate geometry (metric terms) from the mesh elements
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  ! Initialize the model
  call modelobj%Init(mesh,geometry)
  modelobj%tecplot_enabled = .false.

  ! No gravity (constant geopotential is enforced in the initial condition; a
  ! zero gravitational acceleration makes the source term vanish identically).
  modelobj%g = 0.0_prec

  ! Velocity workspace for the vorticity diagnostics
  allocate(uvw(1:controlDegree+1,1:controlDegree+1,1:controlDegree+1, &
               1:mesh%nElem,1:3))

  ! Initial condition: ABC velocity over a constant rho/theta/Phi background
  call set_abc_initial_condition()

  ! Baseline invariants for the drift diagnostics
  call compute_diagnostics(errL2,Eabs,Habs,Rinf)
  e0 = Eabs
  h0 = Habs

  ! Write the initial condition
  call modelobj%WriteModel()
  call modelobj%IncrementIOCounter()

  ! Diagnostic table header + t = 0 row
  write(*,'(A)') ''
  write(*,'(A)') '# ABC flow preservation diagnostics (triply periodic box)'
  write(*,'(A12,4A16)') 't','E_err','dE/E0','dH/H0','R_inf'
  call print_diag_row(0.0_prec,errL2,Eabs,Habs,Rinf)

  call modelobj%SetTimeIntegrator(integrator)

  ! March one iointerval at a time, reporting diagnostics after each interval.
  ! ForwardStep advances the model from its internal time this%t to the target
  ! time in iointerval chunks; passing target = (iout + 1/2)*iointerval
  ! guarantees exactly one chunk per call regardless of floating-point rounding
  ! (int(1.5) = 1, and the leftover half-interval is never consumed because the
  ! chunk advances this%t by exactly iointerval).
  do iout = 1,nReports
    call modelobj%ForwardStep((real(iout,prec)+0.5_prec)*iointerval,dt,iointerval)
    call modelobj%solution%UpdateHost()
    call compute_diagnostics(errL2,Eabs,Habs,Rinf)
    call print_diag_row(real(iout,prec)*iointerval,errL2,Eabs,Habs,Rinf)
  enddo

  ! Guard against blow-up
  if(Eabs /= Eabs .or. Rinf /= Rinf) then
    print*,"Error: diagnostics are inf or nan"
    stop 1
  endif

  ! Write the final state
  call modelobj%WriteModel()

  ! Clean up
  deallocate(uvw)
  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

contains

  subroutine set_abc_initial_condition()
    !! Fill the solution interior with the ABC velocity field carried as momentum
    !! density (rho * u) over a constant rho/theta/Phi background.
    implicit none
    integer :: i,j,k,iEl
    real(prec) :: xq,yq,zq,uu,vv,ww

    do concurrent(i=1:interp%N+1,j=1:interp%N+1, &
                  k=1:interp%N+1,iEl=1:mesh%nElem)

      xq = geometry%x%interior(i,j,k,iEl,1,1)
      yq = geometry%x%interior(i,j,k,iEl,1,2)
      zq = geometry%x%interior(i,j,k,iEl,1,3)

      uu = A_abc*sin(zq)+C_abc*cos(yq)
      vv = B_abc*sin(xq)+A_abc*cos(zq)
      ww = C_abc*sin(yq)+B_abc*cos(xq)

      modelobj%solution%interior(i,j,k,iEl,1) = rho0 ! rho
      modelobj%solution%interior(i,j,k,iEl,2) = rho0*uu ! rho*u
      modelobj%solution%interior(i,j,k,iEl,3) = rho0*vv ! rho*v
      modelobj%solution%interior(i,j,k,iEl,4) = rho0*ww ! rho*w
      modelobj%solution%interior(i,j,k,iEl,5) = rho0*theta0 ! rho*theta
      modelobj%solution%interior(i,j,k,iEl,6) = phi0 ! geopotential

    enddo

    call modelobj%solution%UpdateDevice()

  endsubroutine set_abc_initial_condition

  subroutine compute_diagnostics(errL2rel,Eabs,Habs,Rinf)
    !! Compute the four ABC diagnostics over the (rank-local) collocation grid:
    !!   errL2rel = ||u - u_ABC||_2 / ||u_ABC||_2
    !!   Eabs     = (1/2) integral |u|^2 dV
    !!   Habs     = integral u.omega dV
    !!   Rinf     = max_x |u x omega|
    !! Vorticity omega = curl(u) is obtained from the strong-form mapped gradient
    !! of the velocity field, evaluated inline (see velgrad) from the interpolant
    !! derivative matrix and the constant geometry metric terms. All arrays read
    !! here are host-resident (solution synced via UpdateHost; geometry/interp
    !! are constant after setup), so this routine compiles and runs identically
    !! on CPU and GPU builds. Integrals use Gauss-Lobatto quadrature weights and
    !! the metric Jacobian (cf. CalculateEntropy). Single-rank: no MPI reduction.
    implicit none
    real(prec),intent(out) :: errL2rel,Eabs,Habs,Rinf
    ! Local
    integer :: i,j,k,iEl
    real(prec) :: rho,uu,vv,ww,jac,wq
    real(prec) :: wx,wy,wz,cx,cy,cz,rloc
    real(prec) :: xq,yq,zq,ua,va,wa
    real(prec) :: Eint,Hint,errnum,errden

    call modelobj%solution%UpdateHost()

    ! Primitive velocity (u,v,w) = (rho*u, rho*v, rho*w) / rho on the host grid
    do concurrent(i=1:interp%N+1,j=1:interp%N+1, &
                  k=1:interp%N+1,iEl=1:mesh%nElem)
      rho = modelobj%solution%interior(i,j,k,iEl,1)
      uvw(i,j,k,iEl,1) = modelobj%solution%interior(i,j,k,iEl,2)/rho
      uvw(i,j,k,iEl,2) = modelobj%solution%interior(i,j,k,iEl,3)/rho
      uvw(i,j,k,iEl,3) = modelobj%solution%interior(i,j,k,iEl,4)/rho
    enddo

    Eint = 0.0_prec
    Hint = 0.0_prec
    errnum = 0.0_prec
    errden = 0.0_prec
    Rinf = 0.0_prec

    do iEl = 1,mesh%nElem
      do k = 1,interp%N+1
        do j = 1,interp%N+1
          do i = 1,interp%N+1

            jac = abs(geometry%J%interior(i,j,k,iEl,1))
            wq = interp%qWeights(i)*interp%qWeights(j)*interp%qWeights(k)

            uu = uvw(i,j,k,iEl,1)
            vv = uvw(i,j,k,iEl,2)
            ww = uvw(i,j,k,iEl,3)

            ! Vorticity omega = curl(u); velgrad(ivar,idir,...) = d u_ivar / d x_idir
            wx = velgrad(3,2,i,j,k,iEl)-velgrad(2,3,i,j,k,iEl) ! dw/dy - dv/dz
            wy = velgrad(1,3,i,j,k,iEl)-velgrad(3,1,i,j,k,iEl) ! du/dz - dw/dx
            wz = velgrad(2,1,i,j,k,iEl)-velgrad(1,2,i,j,k,iEl) ! dv/dx - du/dy

            ! Energy and helicity integrands
            Eint = Eint+0.5_prec*(uu*uu+vv*vv+ww*ww)*jac*wq
            Hint = Hint+(uu*wx+vv*wy+ww*wz)*jac*wq

            ! Beltrami residual u x omega (== 0 for an exact Beltrami field)
            cx = vv*wz-ww*wy
            cy = ww*wx-uu*wz
            cz = uu*wy-vv*wx
            rloc = sqrt(cx*cx+cy*cy+cz*cz)
            Rinf = max(Rinf,rloc)

            ! Pointwise error against the analytic ABC field
            xq = geometry%x%interior(i,j,k,iEl,1,1)
            yq = geometry%x%interior(i,j,k,iEl,1,2)
            zq = geometry%x%interior(i,j,k,iEl,1,3)
            ua = A_abc*sin(zq)+C_abc*cos(yq)
            va = B_abc*sin(xq)+A_abc*cos(zq)
            wa = C_abc*sin(yq)+B_abc*cos(xq)
            errnum = errnum+((uu-ua)**2+(vv-va)**2+(ww-wa)**2)*jac*wq
            errden = errden+(ua*ua+va*va+wa*wa)*jac*wq

          enddo
        enddo
      enddo
    enddo

    Eabs = Eint
    Habs = Hint
    errL2rel = sqrt(errnum/errden)

  endsubroutine compute_diagnostics

  function velgrad(ivar,idir,i,j,k,iEl) result(dfdx)
    !! Strong-form physical-space derivative d(u_ivar)/d(x_idir) of the velocity
    !! field uvw at collocation node (i,j,k) of element iEl. This is the
    !! per-node form of the mapped strong gradient (cf. MappedGradient):
    !!
    !!   d f / d x_idir = (1/J) * sum_dim sum_ii D(ii,.) f dsdx(.,idir,dim)
    !!
    !! contracting the reference-space derivative against the contravariant
    !! metric vectors. Uses only host-resident, run-constant geometry/interp
    !! data, so it is backend-agnostic (CPU and GPU builds alike).
    implicit none
    integer,intent(in) :: ivar,idir,i,j,k,iEl
    real(prec) :: dfdx
    integer :: ii
    real(prec) :: s

    s = 0.0_prec
    do ii = 1,interp%N+1
      s = s+interp%dMatrix(ii,i)*uvw(ii,j,k,iEl,ivar)* &
          geometry%dsdx%interior(ii,j,k,iEl,1,idir,1)
    enddo
    do ii = 1,interp%N+1
      s = s+interp%dMatrix(ii,j)*uvw(i,ii,k,iEl,ivar)* &
          geometry%dsdx%interior(i,ii,k,iEl,1,idir,2)
    enddo
    do ii = 1,interp%N+1
      s = s+interp%dMatrix(ii,k)*uvw(i,j,ii,iEl,ivar)* &
          geometry%dsdx%interior(i,j,ii,iEl,1,idir,3)
    enddo

    dfdx = s/geometry%J%interior(i,j,k,iEl,1)

  endfunction velgrad

  subroutine print_diag_row(tt,errL2rel,Eabs,Habs,Rinf)
    !! Print one row of the diagnostic table, converting absolute energy/helicity
    !! to relative drifts against the baseline values e0, h0.
    implicit none
    real(prec),intent(in) :: tt,errL2rel,Eabs,Habs,Rinf
    real(prec) :: dE,dH

    dE = (Eabs-e0)/e0
    dH = 0.0_prec
    if(abs(h0) > 0.0_prec) dH = (Habs-h0)/h0

    write(*,'(ES12.4,4ES16.6)') tt,errL2rel,dE,dH,Rinf

  endsubroutine print_diag_row

endprogram esatmo3d_abc_flow
