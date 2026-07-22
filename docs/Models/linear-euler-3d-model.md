# Linear Euler (3D)

## Definition
The [`SELF_LinearEuler3D_t` module](../ford/module/self_lineareuler3d_t.html) defines the [`LinearEuler3D_t` class](ford/type/lineareuler3d_t.html). In SELF, models are posed in the form of a generic conservation law

$$
  \vec{s}_t + \nabla \cdot \overleftrightarrow{f} = \vec{q}
$$

where $\vec{s}$ is a vector of solution variables, $\overleftrightarrow{f}$ is the conservative flux, and $\vec{q}$ are non-conservative source terms. 

For the Linear Euler equations in 3-D



$$
    \vec{s} = 
    \begin{pmatrix}
    u \\ 
    v \\ 
    w \\
    p \\
    c \\
    \rho_0
    \end{pmatrix}
$$

where $u$, $v$, and $w$ are the $x$, $y$, and $z$ components of the fluid velocity (respectively), and $p$ is the pressure. The density anomaly is *not* carried as a solution variable: for a motionless background state it is slaved to the pressure through the acoustic relation $\rho = p/c^2$ and never feeds back into the velocity or pressure dynamics, so only the velocity components and the pressure are forward-stepped. If the density anomaly is needed as a diagnostic, it can be recovered pointwise as $\rho = p/c^2$. Both the sound speed $c$ and the background density $\rho_0$ are carried as solution variables so that heterogeneous (spatially varying) media are supported; they have zero flux and zero source, so they are static in time and set by the initial condition (mirroring the 2-D model). This is entropy-stable for piecewise-constant material regions aligned with element boundaries. When we assume an ideal gas, and a motionless background state, the conservative fluxes are

$$
    \overleftrightarrow{f} = 
    \begin{pmatrix}
    \frac{p}{\rho_0} \hat{x} \\
    \frac{p}{\rho_0} \hat{y} \\
    \frac{p}{\rho_0} \hat{z} \\
    \rho_0c^2(u \hat{x} + v \hat{y} + w \hat{z}) \\
    0 \\
    0
    \end{pmatrix}
$$

The source term is set to zero.

$$
    \vec{q} = \vec{0}
$$ 

To track stability of the Euler equation, the total entropy function is

$$
    e = \frac{1}{2} \int_V \rho_0 (u^2 + v^2 + w^2) + \frac{p^2}{\rho_0 c^2} \hspace{1mm} dV
$$

## Implementation
The Linear Euler 3D model is implemented as a type extension of the [`DGModel3D` class](../ford/type/dgmodel3d_t.html). The [`LinearEuler3D_t` class](../ford/type/lineareuler3d_t.html) keeps scalar `rho0` and `c` attributes (used as the reference values that fill variables 6 and 5 in the built-in initial conditions) and overrides `SetNumberOfVariables` (to declare `nvar = 6` with `nstepped = 4`, so the last two variables are static), `SetMetadata`, `AdditionalInit`, `entropy_func`, `flux3d`, and `riemannflux3d`. The sound speed lives in `solution(:,:,:,:,5)` and the background density in `solution(:,:,:,:,6)`; both can be set independently per node when initializing the simulation.

### Riemann Solver
The `LinearEuler3D` class is defined using the conservative form of the conservation law. The Riemann solver for the hyperbolic part of the Euler equation is the impedance-matched (characteristic/Godunov) upwind flux, identical in form to the 2-D model's. With the per-side acoustic impedances $Z_L = \rho_{0,L} c_L$ and $Z_R = \rho_{0,R} c_R$ evaluated from the per-node background density and sound speed on either side of the face, the interface normal velocity and pressure are

$$
    u_n^* = \frac{Z_L u_{n,L} + Z_R u_{n,R} + (p_L - p_R)}{Z_L + Z_R}, \qquad
    p^* = \frac{Z_R p_L + Z_L p_R + Z_L Z_R (u_{n,L} - u_{n,R})}{Z_L + Z_R}
$$

and the normal flux is

$$
    \overleftrightarrow{f}_h^* \cdot \hat{n} = 
    \begin{pmatrix}
    p^* n_x / \overline{\rho_0} \\
    p^* n_y / \overline{\rho_0} \\
    p^* n_z / \overline{\rho_0} \\
    \overline{\rho_0} \overline{c^2} u_n^* \\
    0 \\
    0
    \end{pmatrix}, \qquad
    \overline{\rho_0} = \frac{1}{2}(\rho_{0,L} + \rho_{0,R}), \quad
    \overline{c^2} = \frac{1}{2}(c_L^2 + c_R^2)
$$

Because the interface states are resolved with the per-side impedances, material interfaces in a heterogeneous field (density and/or sound-speed jumps) produce the physically correct transmission and reflection. The face-averaged $\overline{\rho_0}$ and $\overline{c^2}$ are used to reconstruct the momentum/pressure fluxes. The details for this implementation can be found in [`self_lineareuler3d_t.f90`](../ford/sourcefile/self_lineareuler3d_t.f90.html)

### Boundary conditions
When initializing the mesh for your Euler 3D equation solver, you can change the boundary conditions to 

* `SELF_BC_Radiation` to set the external state on model boundaries to 0 in the Riemann solver
* `SELF_BC_NoNormalFlow` to set the external normal velocity to the negative of the interior normal velocity and prolong the pressure, tangential velocity, sound speed, and background density (free slip). This effectively creates a reflecting boundary condition.
* `SELF_BC_Prescribed` to set a prescribed external state.


As an example, when using the built-in structured mesh generator, you can do the following

```fortran

type(Mesh3D),target :: mesh
integer :: bcids(1:6)

  bcids(1:6) = (/&
                  SELF_NONORMALFLOW,& ! Bottom boundary condition
                  SELF_NONORMALFLOW,& ! South boundary condition
                  SELF_RADIATION,&    ! East boundary condition
                  SELF_PRESCRIBED,&   ! North boundary condition
                  SELF_RADIATION &    ! West boundary condition
                  SELF_NONORMALFLOW,& ! Top boundary condition
                /)   
  call mesh%StructuredMesh(nxPerTile=5,nyPerTile=5,nzPerTile=5,&
                            nTileX=2,nTileY=2,nTileZ=2,&
                            dx=0.1_prec,dy=0.1_prec,dz=0.1_prec,bcids)

```

!!! note
    See the [Structured Mesh documentation](../MeshGeneration/StructuredMesh.md) for details on using the `structuredmesh` procedure

!!! note
    To set a prescribed state as a function of position and time, you can create a type-extension of the `LinearEuler3D` class and override the [`hbc3d_Prescribed`](../ford/proc/hbc3d_prescribed_model.html) 

#### The no-normal-flow boundary condition
To set the no-normal-flow boundary condition in SELF, we set the external state that is used as input to a Riemann solver. To determine the three components of the velocity field, we use the following conditions

* $\vec{u}_{ext}\cdot \hat{n} = -\vec{u}_{in}\cdot \hat{n}$ 
* $\vec{u}_{ext}\cdot \hat{t}_1 = \vec{u}_{in}\cdot \hat{t}_1$ 
* $\vec{u}_{ext}\cdot \hat{t}_2 = \vec{u}_{in}\cdot \hat{t}_2$ 

where $\hat{n}$ is the outward pointing unit normal vector and $\hat{t}_1$ and $\hat{t}_2$ are mutually orthogonal vectors that are tangent to the boundary surface.

## GPU Acceleration
When building SELF with GPU acceleration enabled, the Linear Euler (3-D) model overrides the following `DGModel3D` type-bound procedures

* `BoundaryFlux`
* `FluxMethod` 
* `SourceMethod`
* `SetBoundaryCondition`
* `SetGradientBoundaryCondition`

These methods are one-level above the usual `pure function` type-bound procedures used to define the riemann solver, flux, source terms, and boundary conditions. These procedures need to be overridden with calls to GPU accelerated kernels to make the solver fully resident on the GPU. 

Out-of-the-box, the no-normal-flow and radiation boundary conditions are GPU accelerated. However, prescribed boundary conditions are CPU-only. We have opted to keep the prescribed boundary conditions CPU-only so that their implementation remains easy-to-use. This implies that some data is copied between host and device every iteration when prescribed boundary conditions are enabled. 

!!! note
    In simulations where no prescribed boundaries are used, or your prescribed boundaries are time independent, you can disable prescribed boundary conditions by explicitly setting `modelobj % prescribed_bcs_enabled = .false.`. This can improve the time-to-solution for your simulation by avoiding unnecessary host-device memory movement. An example of this feature is shown in [`examples/lineareuler3d_spherical_soundwave_radiation.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler3d_spherical_soundwave_radiation.f90)


## Example usage

For examples, see any of the following

* [`examples/linear_euler3d_spherical_soundwave_radiation.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler3d_spherical_soundwave_radiation.f90) - Implements a simulation with a gaussian pressure anomaly as an initial condition in a domain with radiation boundary conditions on all sides. Uses the `SphericalSoundWave` initializer, which fills the uniform sound speed and background density fields.
