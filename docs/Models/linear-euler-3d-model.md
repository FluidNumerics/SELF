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
    \rho_0 \\
    \sigma
    \end{pmatrix}
$$

where $u$, $v$, and $w$ are the $x$, $y$, and $z$ components of the fluid velocity (respectively), $p$ is the pressure, $c$ is the speed of sound, $\rho_0$ is the background density, and $\sigma$ is the relaxation rate (the inverse of a local damping time scale) used to build sponge/damping layers. The density anomaly is *not* carried as a solution variable: for a motionless background state it is slaved to the pressure through the acoustic relation $\rho = p/c^2$ and never feeds back into the velocity or pressure dynamics, so only the velocity components and the pressure are forward-stepped. If the density anomaly is needed as a diagnostic, it can be recovered pointwise as $\rho = p/c^2$. The sound speed $c$, the background density $\rho_0$, and the relaxation rate $\sigma$ are carried as solution variables so that heterogeneous (spatially varying) media and sponge/damping layers are supported; they have zero flux and zero source, so they are static in time and set by the initial condition (mirroring the 2-D model). This is entropy-stable for piecewise-constant material regions aligned with element boundaries. When we assume an ideal gas, and a motionless background state, the conservative fluxes are

$$
    \overleftrightarrow{f} = 
    \begin{pmatrix}
    \frac{p}{\rho_0} \hat{x} \\
    \frac{p}{\rho_0} \hat{y} \\
    \frac{p}{\rho_0} \hat{z} \\
    \rho_0c^2(u \hat{x} + v \hat{y} + w \hat{z}) \\
    0 \\
    0 \\
    0
    \end{pmatrix}
$$

The source term is the linear relaxation ("sponge"/damping) term

$$
    \vec{q} = -\sigma\,(u, v, w, p, 0, 0, 0)^T
$$

With $\sigma = 0$ everywhere - the default, since the solution array is zero-initialized - the source vanishes identically and the model reduces to the undamped linear Euler system. See [Sponge layers and damping](#sponge-layers-and-damping) below.

To track stability of the Euler equation, the total entropy function is

$$
    e = \frac{1}{2} \int_V \rho_0 (u^2 + v^2 + w^2) + \frac{p^2}{\rho_0 c^2} \hspace{1mm} dV
$$

## Implementation
The Linear Euler 3D model is implemented as a type extension of the [`DGModel3D` class](../ford/type/dgmodel3d_t.html). The [`LinearEuler3D_t` class](../ford/type/lineareuler3d_t.html) keeps scalar `rho0` and `c` attributes (used as the reference values that fill variables 6 and 5 in the built-in initial conditions) and overrides `SetNumberOfVariables` (to declare `nvar = 7` with `nstepped = 4`, so the last three variables are static), `SetMetadata`, `AdditionalInit`, `entropy_func`, `flux3d`, `riemannflux3d`, and `SourceMethod`. The sound speed lives in `solution(:,:,:,:,5)`, the background density in `solution(:,:,:,:,6)`, and the relaxation rate in `solution(:,:,:,:,7)`; all three can be set independently per node when initializing the simulation.

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
* `SELF_BC_NoNormalFlow` to set the external normal velocity to the negative of the interior normal velocity and prolong the pressure, tangential velocity, sound speed, background density, and relaxation rate (free slip). This effectively creates a reflecting boundary condition.
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

### Sponge layers and damping

The relaxation rate $\sigma$ (solution variable 7, units s$^{-1}$) is the inverse of a local damping time scale. Wherever $\sigma > 0$, the source term

$$
    \vec{q} = -\sigma\,(u, v, w, p, 0, 0, 0)^T
$$

pulls the acoustic state back toward the motionless, unperturbed background. Because the *same* rate is applied to all three velocity components and to the pressure, the local acoustic energy density decays as $e^{-2\sigma t}$ and the local characteristic structure is untouched: the damping introduces no impedance mismatch of its own between a damped and an undamped region. That is what makes it usable as a **sponge layer** — a shell of elements next to an open boundary in which outgoing waves are absorbed before they can reflect back into the region of interest.

$\sigma$ is set exactly the way $c$ and $\rho_0$ are set — it is a static, per-node solution variable with zero flux and zero source of its own, and it is excluded from the time-stepped variables, so once written it is preserved bitwise for the whole run:

```fortran
this%solution%interior(i,j,k,iel,7) = sigma_value_at_this_node
```

$\sigma = 0$ (the default from the zero-initialized solution array) recovers the undamped model exactly, so existing setups need no changes.

#### Choosing a profile

Two rules of thumb:

* **Ramp $\sigma$ smoothly from zero.** A jump in $\sigma$ at the inner edge of the layer partially reflects the incoming wave. Start the layer at $\sigma = 0$ and grow it toward the boundary, e.g. quadratically or cubically in the normalized depth into the layer.
* **Make the layer a few elements thick and pick $\sigma_{max}$ from the transit time.** A wave crossing a layer of thickness $L$ at speed $c$ is attenuated in amplitude by $\exp\left(-\int \sigma\, \mathrm{d}t\right)$; for a quadratic ramp that integral is $\sigma_{max} L / (3c)$. Aim for a value of a few (an amplitude reduction of $10^{-2}$ to $10^{-3}$) rather than for the largest $\sigma$ that is stable — an over-strong layer reflects more than it absorbs. Keep $\sigma \Delta t \lesssim 1$ so that the explicit time integrator resolves the relaxation.

A layer of thickness `layer` next to every domain boundary, with a quadratic ramp up to `sigma_max`, is written as

```fortran
do iel = 1,mesh%nElem
  do k = 1,modelobj%solution%N+1
    do j = 1,modelobj%solution%N+1
      do i = 1,modelobj%solution%N+1
        x = modelobj%geometry%x%interior(i,j,k,iel,1,1)
        y = modelobj%geometry%x%interior(i,j,k,iel,1,2)
        z = modelobj%geometry%x%interior(i,j,k,iel,1,3)

        ! distance to the nearest domain boundary
        d = min(x-xmin,xmax-x,y-ymin,ymax-y,z-zmin,zmax-z)

        if(d < layer) then
          ! zero at the inner edge of the layer, sigma_max at the boundary
          modelobj%solution%interior(i,j,k,iel,7) = sigma_max*((layer-d)/layer)**2
        else
          modelobj%solution%interior(i,j,k,iel,7) = 0.0_prec
        endif
      enddo
    enddo
  enddo
enddo
call modelobj%solution%UpdateDevice()
```

Pair the layer with `SELF_BC_RADIATION` on the outer boundaries: the sponge removes most of the outgoing energy and the radiation condition handles the remainder. `SELF_BC_NONORMALFLOW` also works — the sponge then absorbs the wave before it can reflect off the wall — which is useful when the mesh has no open boundary.

The relaxation is not restricted to boundary layers: a spatially uniform $\sigma$ acts as a bulk absorber, and a $\sigma$ supported on an arbitrary subregion damps only there.

!!! warning
    $\sigma$ is expected to be non-negative. A negative $\sigma$ amplifies the solution and is not a supported configuration.

## GPU Acceleration
When building SELF with GPU acceleration enabled, the Linear Euler (3-D) model overrides the following `DGModel3D` type-bound procedures

* `BoundaryFlux`
* `FluxMethod` 
* `SourceMethod`
* `SetBoundaryCondition`
* `SetGradientBoundaryCondition`

These methods are one-level above the usual `pure function` type-bound procedures used to define the riemann solver, flux, source terms, and boundary conditions. These procedures need to be overridden with calls to GPU accelerated kernels to make the solver fully resident on the GPU. 

The relaxation source term is fully device-resident: `SourceMethod` launches the `sourcemethod_LinearEuler3D_gpu` kernel, which reads $\sigma$ from solution variable 7 on the device, so enabling a sponge layer adds no host-device traffic.

Out-of-the-box, the no-normal-flow and radiation boundary conditions are GPU accelerated. However, prescribed boundary conditions are CPU-only. We have opted to keep the prescribed boundary conditions CPU-only so that their implementation remains easy-to-use. This implies that some data is copied between host and device every iteration when prescribed boundary conditions are enabled. 

!!! note
    In simulations where no prescribed boundaries are used, or your prescribed boundaries are time independent, you can disable prescribed boundary conditions by explicitly setting `modelobj % prescribed_bcs_enabled = .false.`. This can improve the time-to-solution for your simulation by avoiding unnecessary host-device memory movement. An example of this feature is shown in [`examples/lineareuler3d_spherical_soundwave_radiation.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler3d_spherical_soundwave_radiation.f90)


## Example usage

For examples, see any of the following

* [`examples/linear_euler3d_spherical_soundwave_radiation.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler3d_spherical_soundwave_radiation.f90) - Implements a simulation with a gaussian pressure anomaly as an initial condition in a domain with radiation boundary conditions on all sides. Uses the `SphericalSoundWave` initializer, which fills the uniform sound speed and background density fields.
* [`test/lineareuler3d_sponge_damping.f90`](https://github.com/FluidNumerics/SELF/blob/main/test/lineareuler3d_sponge_damping.f90) - Sponge-layer regression test. A Gaussian pulse is released at the center of a cubic domain with radiation boundaries on all sides; the outermost shell of elements carries a quadratically ramped $\sigma$. The same problem is run twice, with $\sigma = 0$ and with the sponge, and the two are compared. Shows the profile-setting pattern described in [Sponge layers and damping](#sponge-layers-and-damping).
