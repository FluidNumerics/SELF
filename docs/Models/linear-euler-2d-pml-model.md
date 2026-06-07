# Linear Euler (2D) with Perfectly Matched Layer

## Definition

The [`SELF_LinearEuler2D_PML_t` module](../ford/module/self_lineareuler2d_pml_t.html) defines the [`LinearEuler2D_PML_t` class](../ford/type/lineareuler2d_pml_t.html), a type extension of [`LinearEuler2D_t`](./linear-euler-2d-model.md) that adds a Perfectly Matched Layer (PML) absorbing region for open-domain acoustic simulations.

The interior dynamics are identical to the parent linear Euler 2D model. Inside the PML region the governing equations are augmented with damping coefficients $\sigma_x(x)$, $\sigma_y(y)$ and four auxiliary variables $\phi_\rho, \phi_u, \phi_v, \phi_P$ that together drive outgoing waves to zero with minimal reflection at the interior/PML interface.

### Formulation

We follow the unsplit auxiliary-differential-equation (ADE) form of the PML for the linearised Euler equations, after Hu (2001, *JCP* 173:455–480). The modified system is

$$
    \frac{\partial \vec{q}}{\partial t}
    + A\,\frac{\partial \vec{q}}{\partial x}
    + B\,\frac{\partial \vec{q}}{\partial y}
    + (\sigma_x + \sigma_y)\,\vec{q}
    + \sigma_x\,\sigma_y\,\vec{\phi}
    = 0,
    \qquad
    \frac{\partial \vec{\phi}}{\partial t} = \vec{q}
$$

where $\vec{q} = (\rho, u, v, p)$ is the acoustic state, $\vec{\phi} = (\phi_\rho, \phi_u, \phi_v, \phi_P)$ is the auxiliary "memory" vector that integrates $\vec{q}$ in time, and $A$, $B$ are the linear-Euler flux Jacobians (see the [Linear Euler 2D reference](./linear-euler-2d-model.md)).

In the interior we set $\sigma_x = \sigma_y = 0$ and the system reduces exactly to the standard linear Euler 2D model. The auxiliary $\vec{\phi}$ still integrates $\vec{q}$, but with the coupling coefficient $\sigma_x \sigma_y = 0$ it never feeds back into the dynamics and starts at $\vec{\phi}(t=0) = \vec{0}$.

### Solution vector

To accommodate the auxiliary $\vec{\phi}$, the PML model carries `nvar = 9`:

| index | name      | meaning                                      |
|-------|-----------|----------------------------------------------|
| 1     | `rho`     | density perturbation (acoustic)              |
| 2     | `u`       | $x$-velocity                                 |
| 3     | `v`       | $y$-velocity                                 |
| 4     | `P`       | pressure perturbation                        |
| 5     | `c`       | sound speed (per-node, static)               |
| 6     | `phi_rho` | $\int_0^t \rho\,dt$ (PML auxiliary)          |
| 7     | `phi_u`   | $\int_0^t u\,dt$ (PML auxiliary)             |
| 8     | `phi_v`   | $\int_0^t v\,dt$ (PML auxiliary)             |
| 9     | `phi_P`   | $\int_0^t p\,dt$ (PML auxiliary)             |

The auxiliary variables carry **identically zero flux** in both spatial directions. They are evolved exclusively by the source term above. Sound speed $c$ remains static (zero flux, zero source), as in the parent model.

### Damping fields $\sigma_x$, $\sigma_y$

`LinearEuler2D_PML_t` stores the per-node damping coefficients in two `MappedScalar2D` fields, `sigma_x` and `sigma_y` (one variable each). They are populated once by `SetPMLProfile` from a user-supplied interior bounding box, PML thickness, peak value, and polynomial ramp exponent, and are pushed to the GPU automatically when the GPU build is in use.

The default ramp is a cubic polynomial that rises from zero at the interior/PML interface up to the requested $\sigma_\mathrm{max}$ at the outer edge of the layer:

$$
    \sigma_x(x) = \sigma_\mathrm{max}\,\left(\frac{d_x}{L}\right)^p,
    \qquad
    d_x = \max\bigl(0,\ x_\mathrm{min} - x,\ x - x_\mathrm{max}\bigr)
$$

and analogously for $\sigma_y$, where $L$ is `pml_width`, $p$ is `ramp_exponent` (default 3), and $[x_\mathrm{min}, x_\mathrm{max}] \times [y_\mathrm{min}, y_\mathrm{max}]$ is the interior box you want to leave undamped. Nodes outside elements tagged as PML always receive $\sigma = 0$, regardless of position.

### Identifying PML elements

PML elements are picked out of the mesh by material name. `SetPMLProfile` walks `mesh%elemMaterial` and `mesh%materialNames` and applies the damping ramp only to elements whose material name **starts with the configured prefix** (default `"pml"`). All other elements remain at $\sigma_x = \sigma_y = 0$.

You can supply that tagging in two equivalent ways:

* **Programmatic** — write directly to `mesh%elemMaterial` after calling `mesh%StructuredMesh`. This is the path used by the bundled example and test.
* **HOHQMesh ISM-MM** — generate a `.mesh` file in which the outer regions are emitted with material name `pml` (or `pml_east`, `pml_corner_NE`, etc.). `mesh%Read_HOHQMesh` populates `elemMaterial`/`materialNames` from the file automatically.

Both paths are covered step-by-step in the [PML tutorial](../Tutorials/LinearEuler2D/PerfectlyMatchedLayer.md).

## Implementation

`LinearEuler2D_PML_t` is implemented in [`src/SELF_LinearEuler2D_PML_t.f90`](https://github.com/FluidNumerics/SELF/blob/main/src/SELF_LinearEuler2D_PML_t.f90) as a type extension of [`LinearEuler2D_t`](./linear-euler-2d-model.md). It overrides

* `SetNumberOfVariables` — declares `nvar = 9`.
* `SetMetadata` — registers names and units for the four PML auxiliaries (plus `sigma_x`, `sigma_y`).
* `AdditionalInit` — allocates `sigma_x` and `sigma_y` and registers PML-aware no-normal-flow and radiation boundary handlers.
* `AdditionalFree` — releases the damping fields.
* `flux2d`, `riemannflux2d` — reuse the parent linear-Euler flux for variables 1–5 and return zero for variables 6–9 (auxiliaries carry no flux).
* `sourcemethod` — implements the ADE source term node-by-node. We override `sourcemethod` rather than `source2d` because the pure `source2d(s, dsdx)` signature has no access to the per-node $\sigma_x(i,j,iel)$, $\sigma_y(i,j,iel)$ lookups.

A new procedure `SetPMLProfile(x_interior_min, x_interior_max, y_interior_min, y_interior_max, pml_width, sigma_max, ramp_exponent)` populates `sigma_x` and `sigma_y` from the per-element material tags and the geometric ramp.

### Boundary conditions

`LinearEuler2D_PML_t` re-registers both the **no-normal-flow** and **radiation** boundary handlers so that the auxiliary variables 6–9 are also handled at the boundary. Variables 1–5 are treated exactly as in the parent model:

* `SELF_BC_NONORMALFLOW` — reflect $\vec{v}$ about the boundary normal; copy $\rho$, $p$, $c$; zero $\vec{\phi}$ in `extBoundary`.
* `SELF_BC_RADIATION` — set $\rho = u = v = p = 0$ in the exterior state; copy $c$ from the interior side so the Riemann solver sees a consistent sound speed; zero $\vec{\phi}$.

Because the auxiliary variables carry zero Riemann flux, their `extBoundary` value is mathematically irrelevant — we zero it purely for diagnostic cleanliness.

`SELF_BC_PRESCRIBED` is *not* automatically registered. If you need it (e.g. driving a planewave in from one side while absorbing on the other), follow the same pattern as the prescribed-BC examples for the parent model and remember to fill `extBoundary(:,:,:,6:9)` (zero is the natural choice).

## GPU acceleration

When SELF is built with GPU support, the GPU-backed `LinearEuler2D_PML` type (in [`src/gpu/SELF_LinearEuler2D_PML.f90`](https://github.com/FluidNumerics/SELF/blob/main/src/gpu/SELF_LinearEuler2D_PML.f90) and [`src/gpu/SELF_LinearEuler2D_PML.cpp`](https://github.com/FluidNumerics/SELF/blob/main/src/gpu/SELF_LinearEuler2D_PML.cpp)) overrides

* `BoundaryFlux` — launches a 9-variable PML-aware Riemann kernel.
* `FluxMethod` — launches a 9-variable PML-aware interior-flux kernel.
* `SourceMethod` — launches the PML source kernel that reads `sigma_x_gpu`, `sigma_y_gpu`, and `solution_gpu`.

The GPU no-normal-flow and radiation BC handlers are re-registered to point to device-resident kernels (`hbc2d_*_lineareuler2d_pml_gpu`). `sigma_x` and `sigma_y` live on the device after `SetPMLProfile` calls `UpdateDevice`, so no host/device traffic is needed during time stepping unless prescribed BCs are enabled.

## Example usage

A complete worked example with on-the-fly programmatic PML tagging is provided at

* [`examples/linear_euler2d_pml_planewave.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_pml_planewave.f90)

and an absorption regression test (planar pulse, east-side PML, asserts >95% absorption) at

* [`test/lineareuler2d_pml_absorption.f90`](https://github.com/FluidNumerics/SELF/blob/main/test/lineareuler2d_pml_absorption.f90).

See the [PML tutorial](../Tutorials/LinearEuler2D/PerfectlyMatchedLayer.md) for an end-to-end walkthrough including both the programmatic and HOHQMesh ISM-MM workflows for tagging PML elements.
