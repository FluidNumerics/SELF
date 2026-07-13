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
    \rho \\ 
    u \\ 
    v \\ 
    w \\
    p \\
    c
    \end{pmatrix}
$$

where $\rho$ is a density anomaly, referenced to the density $\rho_0$, $u$, $v$, and $w$ are the $x$, $y$, and $z$ components of the fluid velocity (respectively), and $p$ is the pressure. The sound speed $c$ is carried as a solution variable so that heterogeneous (spatially varying) media are supported; it has zero flux and zero source, so it is static in time and is set by the initial condition (mirroring the 2-D model). When we assume an ideal gas, and a motionless background state, the conservative fluxes are

$$
    \overleftrightarrow{f} = 
    \begin{pmatrix}
    \rho_0(u \hat{x} + v \hat{y} + w \hat{z}) \\
    p \hat{x} \\
    p \hat{y} \\
    p \hat{z} \\
    \rho_0c^2(u \hat{x} + v \hat{y} + w \hat{z}) \\
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
The Linear Euler 3D model is implemented as a type extension of the [`DGModel3D` class](../ford/type/dgmodel3d_t.html). The [`LinearEuler3D_t` class](../ford/type/lineareuler3d_t.html) adds parameters for the reference density and the reference speed of sound and overrides the `SetMetadata`, `entropy_func`, `flux3d`, and `riemannflux3d` type-bound procedures.

### Riemann Solver
The `LinearEuler3D` class is defined using the conservative form of the conservation law. The Riemann solver for the hyperbolic part of the Euler equation is the impedance-matched (characteristic/Godunov) upwind flux, identical in form to the 2-D model's. With the acoustic impedances $Z_L = \rho_0 c_L$ and $Z_R = \rho_0 c_R$ evaluated from the per-node sound speed on either side of the face, the interface normal velocity and pressure are

$$
    u_n^* = \frac{Z_L u_{n,L} + Z_R u_{n,R} + (p_L - p_R)}{Z_L + Z_R}, \qquad
    p^* = \frac{Z_R p_L + Z_L p_R + Z_L Z_R (u_{n,L} - u_{n,R})}{Z_L + Z_R}
$$

and the normal flux is

$$
    \overleftrightarrow{f}_h^* \cdot \hat{n} = 
    \begin{pmatrix}
    \rho_0 u_n^* \\
    p^* n_x / \rho_0 \\
    p^* n_y / \rho_0 \\
    p^* n_z / \rho_0 \\
    \rho_0 \overline{c^2} u_n^* \\
    0
    \end{pmatrix}, \qquad \overline{c^2} = \frac{1}{2}(c_L^2 + c_R^2)
$$

Because the interface states are resolved with the one-sided impedances, material interfaces in a heterogeneous sound-speed field produce the physically correct transmission and reflection. The details for this implementation can be found in [`self_lineareuler3d_t.f90`](../ford/sourcefile/self_lineareuler3d_t.f90.html)

### Boundary conditions
When initializing the mesh for your Euler 3D equation solver, you can change the boundary conditions to 

* `SELF_BC_Radiation` to set the external state on model boundaries to 0 in the Riemann solver
* `SELF_BC_NoNormalFlow` to set the external normal velocity to the negative of the interior normal velocity and prolong the density, pressure, and tangential velocity (free slip). This effectively creates a reflecting boundary condition.
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

* [`examples/lineareuler3d_spherical_soundwave_closeddomain.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler3d_spherical_soundwave_closeddomain.f90) - Implements a simulation with a gaussian pressure and density anomaly as an initial condition in a domain with no normal flow boundary conditions on all sides.
* [`examples/lineareuler3d_spherical_soundwave_radiation.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler3d_spherical_soundwave_radiation.f90) - Implements a simulation with a gaussian pressure and density anomaly as an initial condition in a domain with radiation boundary conditions on all sides.
* [`examples/linear_euler3d_planewave_propagation.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler3d_planewave_propagation.f90) - Implements a simulation with a gaussian plane wave that propagates at a $45^\circ$ angle through a square domain. The initial and boundary conditions are all taken as an exact plane wave solution to the Linear Euler equations in 3D. This provides an example for using prescribed boundary conditions.
