# Linear Euler (2D)

## Definition
The [`SELF_LinearEuler2D_t` module](../ford/module/self_lineareuler2d_t.html) defines the [`LinearEuler2D_t` class](ford/type/lineareuler2d_t.html). In SELF, models are posed in the form of a generic conservation law

$$
  \vec{s}_t + \nabla \cdot \overleftrightarrow{f} = \vec{q}
$$

where $\vec{s}$ is a vector of solution variables, $\overleftrightarrow{f}$ is the conservative flux, and $\vec{q}$ are non-conservative source terms. 

For the Linear Euler equations in 2-D

$$
    \vec{s} = 
    \begin{pmatrix}
    u \\ 
    v \\ 
    p \\
    c \\
    \rho_0 \\
    \sigma
    \end{pmatrix}
$$

where $u$ and $v$ are the $x$ and $y$ components of the fluid velocity (respectively), $p$ is the pressure, $c$ is the speed of sound, $\rho_0$ is the background density, and $\sigma$ is the relaxation rate (the inverse of a local damping time scale) used to build sponge/damping layers. The density anomaly is *not* carried as a solution variable: for a motionless background state it is slaved to the pressure through the acoustic relation $\rho = p/c^2$ and never feeds back into the velocity or pressure dynamics, so only the velocity components and the pressure are forward-stepped. If the density anomaly is needed as a diagnostic, it can be recovered pointwise as $\rho = p/c^2$. The sound speed $c$, the background density $\rho_0$, and the relaxation rate $\sigma$ are carried as per-node solution variables so that they can vary in space (e.g. across material interfaces, or through the thickness of a sponge layer); their flux and source are identically zero, so $c$, $\rho_0$, and $\sigma$ are held fixed in time at each spatial location. This is entropy-stable for piecewise-constant material regions aligned with element boundaries: element interiors have $\nabla \rho_0 = \nabla c = 0$ (so the flux-divergence form is exact) and the impedance-matched Riemann flux handles the jumps at faces.

When we assume an ideal gas, and a motionless background state, the conservative fluxes are

$$
    \overleftrightarrow{f} = 
    \begin{pmatrix}
    \frac{p}{\rho_0} \hat{x} \\
    \frac{p}{\rho_0} \hat{y} \\
    \rho_0 c^2 (u \hat{x} + v \hat{y}) \\
    \vec{0} \\
    \vec{0} \\
    \vec{0}
    \end{pmatrix}
$$

The source term is the linear relaxation ("sponge"/damping) term

$$
    \vec{q} = -\sigma
    \begin{pmatrix}
    u \\
    v \\
    p \\
    0 \\
    0 \\
    0
    \end{pmatrix}.
$$

With $\sigma = 0$ everywhere - the default, since the solution array is zero-initialized - the source vanishes identically and the model reduces to the undamped linear Euler system. See [Sponge layers and damping](#sponge-layers-and-damping) below.

To track stability of the Euler equation, the total entropy function is

$$
    e = \frac{1}{2} \int_V \rho_0\,(u^2 + v^2) + \frac{p^2}{\rho_0 c^2} \hspace{1mm} dV.
$$

## Implementation
The Linear Euler 2D model is implemented as a type extension of the [`DGModel2D` class](../ford/type/dgmodel2d_t.html). The [`LinearEuler2D_t` class](../ford/type/lineareuler2d_t.html) keeps a scalar `rho0` attribute (used as the reference value that fills variable 5 in the built-in initial conditions) and overrides `SetNumberOfVariables` (to declare `nvar = 6` with `nstepped = 3`, so the last three variables are static), `SetMetadata`, `AdditionalInit`, `entropy_func`, `flux2d`, `riemannflux2d`, and `SourceMethod`. The sound speed lives in `solution(:,:,:,4)`, the background density in `solution(:,:,:,5)`, and the relaxation rate in `solution(:,:,:,6)`; all three can be set independently per node when initializing the simulation.

### Riemann Solver
The `LinearEuler2D` class is defined using the conservative form of the conservation law. The interface flux is the exact upwind (Godunov) solver for the linearized acoustic system obtained by characteristic decomposition. The normal-flux Jacobian has eigenstructure

* $+c$: right-going acoustic mode, $W_+ = \rho_0 c\, u_n + p$
* $-c$: left-going  acoustic mode, $W_- = -\rho_0 c\, u_n + p$
* $0$: tangential vorticity mode, $u_t$

Upwinding each mode by its characteristic direction at the face gives the impedance-matched interface state

$$
    u_n^* = \frac{Z_L u_{n,L} + Z_R u_{n,R} + (p_L - p_R)}{Z_L + Z_R}, \quad
    p^*   = \frac{Z_R p_L + Z_L p_R + Z_L Z_R (u_{n,L} - u_{n,R})}{Z_L + Z_R},
$$

with the per-side acoustic impedance $Z = \rho_0 c$, where each side uses its own background density $\rho_0$ and sound speed $c$. This reduces correctly to Fresnel reflection / transmission across an impedance jump (e.g. $Z_L \neq Z_R$ at a material interface, whether from a density jump, a sound-speed jump, or both). The interface flux is then

$$
    \overleftrightarrow{f}^* \cdot \hat{n} =
    \begin{pmatrix}
    p^*\, n_x / \overline{\rho_0} \\
    p^*\, n_y / \overline{\rho_0} \\
    \overline{\rho_0}\,\overline{c^2}\,u_n^* \\
    0 \\
    0
    \end{pmatrix}, \qquad
    \overline{\rho_0} = \tfrac{1}{2}(\rho_{0,L} + \rho_{0,R}), \quad
    \overline{c^2} = \tfrac{1}{2}(c_L^2 + c_R^2).
$$

The reconstructed momentum/pressure fluxes use the face-averaged $\overline{\rho_0}$ and $\overline{c^2}$ as a pragmatic treatment of the non-conservative products $p/\rho_0$ and $\rho_0 c^2 \nabla \cdot \vec{v}$ at a face where the coefficients jump; the physically important reflection/transmission is carried exactly by the per-side impedances above. The previously used local Lax-Friedrichs solver was found to over-dissipate the tangential mode and to fail to stably handle impedance mismatch at high polynomial order (aliasing instability at material interfaces), and has been replaced by the characteristic flux above. Details are in [`self_lineareuler2d_t.f90`](../ford/sourcefile/self_lineareuler2d_t.f90.html).

### Boundary conditions
Boundary conditions are managed through the [extensible boundary-condition system](../Learning/BoundaryConditions.md). Each model registers its hyperbolic boundary conditions inside `AdditionalInit` by calling `hyperbolicBCs%RegisterBoundaryCondition(id, name, fn)`. `LinearEuler2D_t` registers both a no-normal-flow and a radiation handler out of the box (CPU); the GPU build overwrites both registrations with equivalent device kernels. Prescribed boundary conditions are registered by user code (or by an example subclass) so that the external state can be set as a function of position and time.

The built-in boundary identifiers used with the mesh generators are

* `SELF_BC_RADIATION` — set the external state on the boundary to zero in the Riemann solver (open/non-reflecting).
* `SELF_BC_NONORMALFLOW` — reflect the velocity vector about the boundary normal and prolong $p$, $c$, $\rho_0$, and $\sigma$. This produces a reflecting (free-slip) wall and works for arbitrarily oriented normals.
* `SELF_BC_PRESCRIBED` — use a user-registered handler to fill the external state.

As an example, when using the built-in structured mesh generator,

```fortran

type(Mesh2D),target :: mesh
integer :: bcids(1:4)

  bcids(1:4) = (/&
                  SELF_BC_NONORMALFLOW,& ! South boundary condition
                  SELF_BC_RADIATION,&    ! East boundary condition
                  SELF_BC_PRESCRIBED,&   ! North boundary condition
                  SELF_BC_RADIATION &    ! West boundary condition
                /)
  call mesh%StructuredMesh(nxPerTile=5,nyPerTile=5,&
                            nTileX=2,nTileY=2,&
                            dx=0.1_prec,dy=0.1_prec,bcids)

```

!!! note
    See the [Structured Mesh documentation](../MeshGeneration/StructuredMesh.md) for details on using the `structuredmesh` procedure, and the [Boundary Condition System](../Learning/BoundaryConditions.md) for how to register new BC handlers.

!!! note
    To set a prescribed state as a function of position and time, create a type-extension of `LinearEuler2D` and register a custom BC method against `SELF_BC_PRESCRIBED` from `AdditionalInit`. Remember that your handler must also fill `solution%extBoundary(:,:,:,4)` (sound speed) and `solution%extBoundary(:,:,:,5)` (background density) with the appropriate values at the boundary — the planewave examples show this pattern. The relaxation rate carries no flux, so `solution%extBoundary(:,:,:,6)` never enters the tendency and does not need to be filled.

### Setting the sound speed and background density

Because $c$ and $\rho_0$ are solution variables, you initialize them the same way you initialize $u$, $v$, and $p$:

```fortran
this%solution%interior(i,j,iel,4) = c_value_at_this_node
this%solution%interior(i,j,iel,5) = rho0_value_at_this_node
```

For a uniform background, set every node to the same constants. For a piecewise-constant medium (e.g. bone embedded in marrow), assign the local material's $c$ and $\rho_0$. The `SphericalSoundWave` initializer takes the (uniform) sound speed as an explicit argument and fills the background density from the scalar `this%rho0`; its `rhoprime` argument sets the pressure-pulse amplitude through the acoustic relation $p = \rho' c^2$:

```fortran
call model%SphericalSoundWave(rhoprime=1.0e-2_prec, Lr=0.1_prec, &
                              x0=0.5_prec, y0=0.5_prec, c=1.0_prec)
```


### Sponge layers and damping

The relaxation rate $\sigma$ (solution variable 6, units s$^{-1}$) is the inverse of a local damping time scale. Wherever $\sigma > 0$, the source term

$$
    \vec{q} = -\sigma\,(u, v, p, 0, 0, 0)^T
$$

pulls the acoustic state back toward the motionless, unperturbed background. Because the *same* rate is applied to the velocity components and to the pressure, the local acoustic energy density decays as $e^{-2\sigma t}$ and the local characteristic structure is untouched: the damping introduces no impedance mismatch of its own between a damped and an undamped region. That is what makes it usable as a **sponge layer** — a band of elements next to an open boundary in which outgoing waves are absorbed before they can reflect back into the region of interest.

$\sigma$ is set exactly the way $c$ and $\rho_0$ are set — it is a static, per-node solution variable with zero flux and zero source of its own, and it is excluded from the time-stepped variables, so once written it is preserved bitwise for the whole run:

```fortran
this%solution%interior(i,j,iel,6) = sigma_value_at_this_node
```

$\sigma = 0$ (the default from the zero-initialized solution array) recovers the undamped model exactly, so existing setups need no changes.

#### Choosing a profile

Two rules of thumb:

* **Ramp $\sigma$ smoothly from zero.** A jump in $\sigma$ at the inner edge of the layer partially reflects the incoming wave. Start the layer at $\sigma = 0$ and grow it toward the boundary, e.g. quadratically or cubically in the normalized depth into the layer.
* **Make the layer a few elements thick and pick $\sigma_{max}$ from the transit time.** A wave crossing a layer of thickness $L$ at speed $c$ is attenuated in amplitude by $\exp\left(-\int \sigma\, \mathrm{d}t\right)$; for a quadratic ramp that integral is $\sigma_{max} L / (3c)$. Aim for a value of a few (an amplitude reduction of $10^{-2}$ to $10^{-3}$) rather than for the largest $\sigma$ that is stable — an over-strong layer reflects more than it absorbs. Keep $\sigma \Delta t \lesssim 1$ so that the explicit time integrator resolves the relaxation.

A layer of thickness `layer` next to every domain boundary, with a quadratic ramp up to `sigma_max`, is written as

```fortran
do iel = 1,mesh%nElem
  do j = 1,modelobj%solution%N+1
    do i = 1,modelobj%solution%N+1
      x = modelobj%geometry%x%interior(i,j,iel,1,1)
      y = modelobj%geometry%x%interior(i,j,iel,1,2)

      ! distance to the nearest domain boundary
      d = min(x-xmin,xmax-x,y-ymin,ymax-y)

      if(d < layer) then
        ! zero at the inner edge of the layer, sigma_max at the boundary
        modelobj%solution%interior(i,j,iel,6) = sigma_max*((layer-d)/layer)**2
      else
        modelobj%solution%interior(i,j,iel,6) = 0.0_prec
      endif
    enddo
  enddo
enddo
call modelobj%solution%UpdateDevice()
```

Pair the layer with `SELF_BC_RADIATION` on the outer boundaries: the sponge removes most of the outgoing energy and the radiation condition handles the remainder. `SELF_BC_NONORMALFLOW` also works — the sponge then absorbs the wave before it can reflect off the wall — which is useful when the mesh has no open boundary.

The relaxation is not restricted to boundary layers: a spatially uniform $\sigma$ acts as a bulk absorber, and a $\sigma$ supported on an arbitrary subregion damps only there.

!!! note
    For a boundary treatment that is *exactly* non-reflecting for the linearized Euler equations rather than merely strongly absorbing, see the [Linear Euler 2D PML model](linear-euler-2d-pml-model.md), which implements the Hu (2001) unsplit perfectly matched layer. The sponge layer described here is cheaper, but not free: $\sigma$ raises `nvar` from 5 to 6, which adds one variable slot to every model buffer sized by `nvar` (`solution`, `dSdt`, `workSol`, `solutionGradient`, `flux`, `source`, and their boundary arrays) whether or not a sponge is used. The PML instead carries three time-stepped auxiliary variables (`nvar = 7`) *plus* two extra per-node `MappedScalar2D` fields for $\sigma_x$ and $\sigma_y$, and its auxiliaries are integrated in time rather than held static. Prefer the sponge when strong absorption is enough; prefer the PML when even small reflections matter.

!!! warning
    $\sigma$ is expected to be non-negative. A negative $\sigma$ amplifies the solution and is not a supported configuration.

## GPU Acceleration
When building SELF with GPU acceleration enabled, the Linear Euler (2-D) model overrides the following `DGModel2D` type-bound procedures

* `BoundaryFlux`
* `FluxMethod` 
* `SourceMethod`
* `SetBoundaryCondition`
* `SetGradientBoundaryCondition`

These methods are one-level above the usual `pure function` type-bound procedures used to define the riemann solver, flux, source terms, and boundary conditions. These procedures need to be overridden with calls to GPU accelerated kernels to make the solver fully resident on the GPU. 

The relaxation source term is fully device-resident: `SourceMethod` launches the `sourcemethod_LinearEuler2D_gpu` kernel, which reads $\sigma$ from solution variable 6 on the device, so enabling a sponge layer adds no host-device traffic.

Out-of-the-box, the no-normal-flow and radiation boundary conditions are GPU accelerated. However, prescribed boundary conditions are CPU-only. We have opted to keep the prescribed boundary conditions CPU-only so that their implementation remains easy-to-use. This implies that some data is copied between host and device every iteration when prescribed boundary conditions are enabled. 

!!! note
    In simulations where no prescribed boundaries are used, or your prescribed boundaries are time independent, you can disable prescribed boundary conditions by explicitly setting `modelobj % prescribed_bcs_enabled = .false.`. This can improve the time-to-solution for your simulation by avoiding unnecessary host-device memory movement. An example of this feature is shown in [`examples/lineareuler2d_spherical_soundwave_closeddomain.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_spherical_soundwave_closeddomain.f90)


## Example usage

For examples, see any of the following

* [`examples/lineareuler2d_spherical_soundwave_closeddomain.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_spherical_soundwave_closeddomain.f90) - Simulation with a gaussian pressure anomaly as an initial condition in a domain with no-normal-flow boundary conditions on all sides. Demonstrates uniform sound speed via the `SphericalSoundWave` initializer.
* [`examples/linear_euler2d_planewave_propagation.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_planewave_propagation.f90) - Gaussian plane wave that propagates at a $45^\circ$ angle through a square domain. The initial and boundary conditions are an exact plane-wave solution to the linear Euler equations. The example subclass carries its own `c` attribute and writes it into `solution(...,4)` for both the initial condition and the prescribed boundary state.
* [`examples/linear_euler2d_planewave_reflection.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_planewave_reflection.f90) - Gaussian plane wave reflected off a wall at $x=1$ via the method of images. Combines prescribed boundary conditions with no-normal-flow on the reflecting side.
* [`test/lineareuler2d_sponge_damping.f90`](https://github.com/FluidNumerics/SELF/blob/main/test/lineareuler2d_sponge_damping.f90) - Sponge-layer regression test. A Gaussian pulse is released at the center of a square domain with radiation boundaries on all sides; the outermost two rings of elements carry a quadratically ramped $\sigma$. The same problem is run twice, with $\sigma = 0$ and with the sponge, and the two are compared. Shows the profile-setting pattern described in [Sponge layers and damping](#sponge-layers-and-damping).
* [`examples/linear_euler2d_boneandmarrow.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_boneandmarrow.f90) - Heterogeneous-medium test on a HOHQMesh ISM-MM mesh tagged with three materials (muscle/bone/marrow). Each material is mapped to a representative sound speed and background density, written into `solution(...,4)` and `solution(...,5)`, and a Gaussian acoustic pulse refracts and reflects at the material interfaces. Exercises the impedance-matched Riemann solver across $\rho_0 c$ (impedance) discontinuities.