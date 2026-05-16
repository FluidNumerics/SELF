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
    \rho \\ 
    u \\ 
    v \\ 
    p \\
    c
    \end{pmatrix}
$$

where $\rho$ is a density anomaly, referenced to the density $\rho_0$, $u$ and $v$ are the $x$ and $y$ components of the fluid velocity (respectively), $p$ is the pressure, and $c$ is the speed of sound. The sound speed is carried as a per-node solution variable so that it can vary in space (e.g. across material interfaces); its flux and source are identically zero, so $c$ is held fixed in time at each spatial location.

When we assume an ideal gas, and a motionless background state, the conservative fluxes are

$$
    \overleftrightarrow{f} = 
    \begin{pmatrix}
    \rho_0(u \hat{x} + v \hat{y}) \\
    \frac{p}{\rho_0} \hat{x} \\
    \frac{p}{\rho_0} \hat{y} \\
    \rho_0 c^2 (u \hat{x} + v \hat{y}) \\
    \vec{0}
    \end{pmatrix}
$$

The source term is zero,

$$
    \vec{q} = \vec{0}.
$$ 

To track stability of the Euler equation, the total entropy function is

$$
    e = \frac{1}{2} \int_V \rho_0\,(u^2 + v^2) + \frac{p^2}{\rho_0 c^2} \hspace{1mm} dV.
$$

## Implementation
The Linear Euler 2D model is implemented as a type extension of the [`DGModel2D` class](../ford/type/dgmodel2d_t.html). The [`LinearEuler2D_t` class](../ford/type/lineareuler2d_t.html) adds a parameter for the reference density $\rho_0$ and overrides `SetNumberOfVariables` (to declare `nvar = 5`), `SetMetadata`, `AdditionalInit`, `entropy_func`, `flux2d`, and `riemannflux2d`. The sound speed is no longer a scalar attribute of the model — it lives in `solution(:,:,:,5)` and can be set independently per node when initializing the simulation.

### Riemann Solver
The `LinearEuler2D` class is defined using the conservative form of the conservation law. The interface flux is the exact upwind (Godunov) solver for the linearized acoustic system obtained by characteristic decomposition. The normal-flux Jacobian has eigenstructure

* $+c$: right-going acoustic mode, $W_+ = \rho_0 c\, u_n + p$
* $-c$: left-going  acoustic mode, $W_- = -\rho_0 c\, u_n + p$
* $0$: entropy density mode, $W_0 = \rho - p/c^2$
* $0$: tangential vorticity mode, $u_t$

Upwinding each mode by its characteristic direction at the face gives the impedance-matched interface state

$$
    u_n^* = \frac{Z_L u_{n,L} + Z_R u_{n,R} + (p_L - p_R)}{Z_L + Z_R}, \quad
    p^*   = \frac{Z_R p_L + Z_L p_R + Z_L Z_R (u_{n,L} - u_{n,R})}{Z_L + Z_R},
$$

with acoustic impedance $Z = \rho_0 c$. This reduces correctly to Fresnel reflection / transmission across an impedance jump (e.g. $c_L \neq c_R$ at a material interface). The interface flux is then

$$
    \overleftrightarrow{f}^* \cdot \hat{n} =
    \begin{pmatrix}
    \rho_0\,u_n^* \\
    p^*\, n_x / \rho_0 \\
    p^*\, n_y / \rho_0 \\
    \rho_0\,\overline{c^2}\,u_n^* \\
    0
    \end{pmatrix}, \qquad
    \overline{c^2} = \tfrac{1}{2}(c_L^2 + c_R^2).
$$

The pressure-flux uses an averaged $c^2$ as a pragmatic treatment of the non-conservative product $\rho_0 c^2 \nabla \cdot \vec{v}$ at a face where $c$ jumps. The previously used local Lax-Friedrichs solver was found to over-dissipate the tangential and entropy modes and to fail to stably handle impedance mismatch at high polynomial order (aliasing instability at material interfaces), and has been replaced by the characteristic flux above. Details are in [`self_lineareuler2d_t.f90`](../ford/sourcefile/self_lineareuler2d_t.f90.html).

### Boundary conditions
Boundary conditions are managed through the [extensible boundary-condition system](../Learning/BoundaryConditions.md). Each model registers its hyperbolic boundary conditions inside `AdditionalInit` by calling `hyperbolicBCs%RegisterBoundaryCondition(id, name, fn)`. `LinearEuler2D_t` registers a no-normal-flow handler out of the box; the GPU build additionally registers a radiation handler that runs on the device. Prescribed boundary conditions are registered by user code (or by an example subclass) so that the external state can be set as a function of position and time.

The built-in boundary identifiers used with the mesh generators are

* `SELF_BC_RADIATION` — set the external state on the boundary to zero in the Riemann solver (open/non-reflecting).
* `SELF_BC_NONORMALFLOW` — reflect the velocity vector about the boundary normal and prolong $\rho$, $p$, and $c$. This produces a reflecting (free-slip) wall and works for arbitrarily oriented normals.
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
    To set a prescribed state as a function of position and time, create a type-extension of `LinearEuler2D` and register a custom BC method against `SELF_BC_PRESCRIBED` from `AdditionalInit`. Remember that your handler must also fill `solution%extBoundary(:,:,:,5)` with the appropriate sound-speed value at the boundary — the planewave examples show this pattern.

### Setting the sound speed

Because $c$ is a solution variable, you initialize it the same way you initialize $\rho$, $u$, $v$, and $p$:

```fortran
this%solution%interior(i,j,iel,5) = c_value_at_this_node
```

For a uniform background, set every node to the same constant. For a piecewise-constant medium (e.g. bone embedded in marrow), assign the local material's $c$. The `SphericalSoundWave` initializer takes the (uniform) sound speed as an explicit argument:

```fortran
call model%SphericalSoundWave(rhoprime=1.0e-2_prec, Lr=0.1_prec, &
                              x0=0.5_prec, y0=0.5_prec, c=1.0_prec)
```


## GPU Acceleration
When building SELF with GPU acceleration enabled, the Linear Euler (2-D) model overrides the following `DGModel2D` type-bound procedures

* `BoundaryFlux`
* `FluxMethod` 
* `SourceMethod`
* `SetBoundaryCondition`
* `SetGradientBoundaryCondition`

These methods are one-level above the usual `pure function` type-bound procedures used to define the riemann solver, flux, source terms, and boundary conditions. These procedures need to be overridden with calls to GPU accelerated kernels to make the solver fully resident on the GPU. 

Out-of-the-box, the no-normal-flow and radiation boundary conditions are GPU accelerated. However, prescribed boundary conditions are CPU-only. We have opted to keep the prescribed boundary conditions CPU-only so that their implementation remains easy-to-use. This implies that some data is copied between host and device every iteration when prescribed boundary conditions are enabled. 

!!! note
    In simulations where no prescribed boundaries are used, or your prescribed boundaries are time independent, you can disable prescribed boundary conditions by explicitly setting `modelobj % prescribed_bcs_enabled = .false.`. This can improve the time-to-solution for your simulation by avoiding unnecessary host-device memory movement. An example of this feature is shown in [`examples/lineareuler2d_spherical_soundwave_closeddomain.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_spherical_soundwave_closeddomain.f90)


## Example usage

For examples, see any of the following

* [`examples/lineareuler2d_spherical_soundwave_closeddomain.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_spherical_soundwave_closeddomain.f90) - Simulation with a gaussian pressure and density anomaly as an initial condition in a domain with no-normal-flow boundary conditions on all sides. Demonstrates uniform sound speed via the `SphericalSoundWave` initializer.
* [`examples/linear_euler2d_planewave_propagation.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_planewave_propagation.f90) - Gaussian plane wave that propagates at a $45^\circ$ angle through a square domain. The initial and boundary conditions are an exact plane-wave solution to the linear Euler equations. The example subclass carries its own `c` attribute and writes it into `solution(...,5)` for both the initial condition and the prescribed boundary state.
* [`examples/linear_euler2d_planewave_reflection.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_planewave_reflection.f90) - Gaussian plane wave reflected off a wall at $x=1$ via the method of images. Combines prescribed boundary conditions with no-normal-flow on the reflecting side.
* [`examples/linear_euler2d_boneandmarrow.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_boneandmarrow.f90) - Heterogeneous-medium test on a HOHQMesh ISM-MM mesh tagged with three materials (muscle/bone/marrow). Each material is mapped to a representative sound speed written into `solution(...,5)`, and a Gaussian acoustic pulse refracts and reflects at the material interfaces. Exercises the impedance-matched Riemann solver across $c$ discontinuities.