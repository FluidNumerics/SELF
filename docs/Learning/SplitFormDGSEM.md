# Split-Form Entropy-Conserving DGSEM

This page documents the split-form discontinuous Galerkin spectral element method (DGSEM) implemented in SELF's `ECDGModel2D_t` and `ECDGModel3D_t` classes. The formulation follows Gassner, Winters, and Kopriva (2016) and is aligned with the [Trixi.jl](https://github.com/trixi-framework/Trixi.jl) implementation for curved meshes.

## Mathematical Background

### Conservation Law
We consider a system of conservation laws on a domain $\Omega$:

$$
\mathbf{s}_t + \vec{\nabla} \cdot \vec{\mathbf{f}}(\mathbf{s}) = \mathbf{q}
$$

where $\mathbf{s}$ is the state vector, $\vec{\mathbf{f}}$ is the flux, and $\mathbf{q}$ is a source term.

### Entropy Conservation
A convex entropy function $\eta(\mathbf{s})$ satisfies the entropy conservation law $\eta_t + \vec{\nabla} \cdot \vec{\Psi} = 0$ (see [Provable Stability](ProvableStability.md)). A numerical scheme is **entropy-conserving** if it satisfies a discrete analogue of this conservation law, and **entropy-stable** if entropy can only decrease (mimicking physical dissipation at discontinuities).

### Two-Point Entropy-Conserving Flux
An entropy-conserving two-point flux $\mathbf{f}^\#(\mathbf{s}_L, \mathbf{s}_R)$ is a symmetric function of two states satisfying the Tadmor condition:

$$
(\mathbf{w}_R - \mathbf{w}_L)^T \mathbf{f}^\#(\mathbf{s}_L, \mathbf{s}_R) = \Psi_R - \Psi_L
$$

where $\mathbf{w} = \partial \eta / \partial \mathbf{s}$ are the entropy variables and $\Psi$ is the entropy flux potential.

**Example:** For linear advection $f = a s$ with entropy $\eta = s^2/2$, the arithmetic mean $f^\#(s_L, s_R) = a(s_L + s_R)/2$ is entropy-conserving.

## The Split-Form Derivative Operator

### Standard Derivative Matrix $D$
The standard DGSEM derivative matrix $D$ satisfies $D \cdot \mathbf{1} = 0$ (derivative of a constant is zero) and the SBP property:

$$
M D + D^T M = B
$$

where $M = \text{diag}(w_0, \ldots, w_N)$ is the mass matrix (quadrature weights) and $B$ is the boundary operator.

### Split-Form Matrix $D_\text{split}$
The split-form derivative matrix is defined as:

$$
D_\text{split} = D - \tfrac{1}{2} M^{-1} B
$$

**Key property:** $D_\text{split}$ is skew-symmetric under the $M$ inner product:

$$
M D_\text{split} + D_\text{split}^T M = 0
$$

This follows directly from the SBP property: $M D_\text{split} + D_\text{split}^T M = (M D + D^T M) - B = B - B = 0$.

**$D_\text{split}$ is NOT an SBP operator.** Its weighted symmetric part is zero (not the boundary term $B$). This means the volume integral formed with $D_\text{split}$ contributes **zero** to the entropy rate. All entropy change passes through the surface term.

For Gauss-Lobatto-Legendre (GLL) nodes, $D_\text{split}$ reduces to $(D - D^T)/2$, and SELF requires GLL nodes for `TwoPointVector` types.

In SELF, the split-form matrix is stored as `interp%dSplitMatrix` in the Lagrange class:

```fortran
dSplitMatrix(ii,i) = dMatrix(ii,i) &
  - 0.5*(bMatrix(i,2)*bMatrix(ii,2) - bMatrix(i,1)*bMatrix(ii,1)) / qWeights(i)
```

### Relationship to Trixi.jl
Trixi.jl defines `derivative_split = 2D - M^{-1}B` and applies it **without** a factor of 2. SELF defines `dSplitMatrix = D - 0.5 M^{-1}B` and multiplies by 2 in the divergence formula. These are algebraically identical:

$$
2 \cdot D_\text{split}^\text{SELF} = D_\text{split}^\text{Trixi}
$$


## The EC-DGSEM Semidiscretization

### Full Tendency
The `ECDGModel2D_t` computes the time tendency as:

$$
\frac{d\mathbf{s}}{dt} = \mathbf{q} - \frac{1}{J} \left[ 2 \sum_n D_{\text{split},n,i} \, \tilde{F}^r(n,i,j) \;+\; M^{-1} B^T \mathbf{f}_\text{Riemann} \right]
$$

where:

- $\tilde{F}^r(n,i,j)$ is a **scalar contravariant** two-point flux for direction $r$
- $\mathbf{f}_\text{Riemann}$ is the surface Riemann flux (e.g., Local Lax-Friedrichs)
- $J$ is the element Jacobian

### Equivalence to the Strong-Form Penalty
Using $D_\text{split}$ for the volume combined with the plain Riemann flux on the surface is algebraically identical to using the standard $D$ for the volume with a **penalty** $(f_\text{Riemann} - f_\text{local})$ on the surface:

$$
2 D_\text{split} \tilde{F} + M^{-1} B^T f_\text{Riemann} = 2 D \tilde{F} + M^{-1} B^T (f_\text{Riemann} - \tilde{F}_\text{boundary})
$$

SELF uses the left-hand form (matching Trixi.jl's `SurfaceIntegralWeakForm`).

## Contravariant Two-Point Fluxes on Curved Meshes

### Convention
Following Trixi.jl for curved meshes, `interior(n,i,j,iEl,iVar,r)` stores a **scalar** contravariant two-point flux for the $r$-th computational direction. Each direction uses its own node pairing and its own averaged metric:

| Direction $r$ | Node pair | Averaged metric |
|:---:|---|---|
| $\xi^1$ | $(i,j) \leftrightarrow (n,j)$ | $\overline{J\mathbf{a}^1} = \tfrac{1}{2}\left(J\mathbf{a}^1_{i,j} + J\mathbf{a}^1_{n,j}\right)$ |
| $\xi^2$ | $(i,j) \leftrightarrow (i,n)$ | $\overline{J\mathbf{a}^2} = \tfrac{1}{2}\left(J\mathbf{a}^2_{i,j} + J\mathbf{a}^2_{i,n}\right)$ |
| $\xi^3$ | $(i,j,k) \leftrightarrow (i,j,n)$ | $\overline{J\mathbf{a}^3} = \tfrac{1}{2}\left(J\mathbf{a}^3_{i,j,k} + J\mathbf{a}^3_{i,j,n}\right)$ |

The contravariant flux is the dot product of the averaged metric with the physical EC flux:

$$
\tilde{F}^r(n,i,j) = \sum_d \overline{J a^r_d} \; f^\#_d(\mathbf{s}_L, \mathbf{s}_R)
$$

where $f^\#_d$ is the $d$-th physical component of the entropy-conserving two-point flux evaluated between the direction-specific node pair.

### Why Not Physical-Space Storage?
A single `interior(n,i,j,...,d)` array with $d$ indexing physical directions cannot correctly serve both the $\xi^1$ and $\xi^2$ sums, because the same index $n$ refers to different physical node partners in each direction. Storing pre-projected scalar contravariant fluxes (one per computational direction) avoids this cross-contamination entirely.

### Implementation
`TwoPointFluxMethod` in `ECDGModel2D_t` computes the contravariant projection:

```fortran
! xi^1: pair (i,j)-(nn,j), project onto avg(Ja^1)
sR = solution(nn,j,iel,:)
Fphys = twopointflux2d(sL, sR)         ! physical EC flux (nvar x 2)
Fc = sum_d avg(Ja^1_d) * Fphys(:,d)    ! scalar contravariant
twoPointFlux%interior(nn,i,j,...,1) = Fc

! xi^2: pair (i,j)-(i,nn), project onto avg(Ja^2)
sR = solution(i,nn,iel,:)
Fphys = twopointflux2d(sL, sR)
Fc = sum_d avg(Ja^2_d) * Fphys(:,d)
twoPointFlux%interior(nn,i,j,...,2) = Fc
```

`MappedDivergence` then simply applies `Divergence / J` -- no metric operations inside.

## Implementing a Concrete EC Model

To create an entropy-conserving model, extend `ECDGModel2D_t` (or `ECDGModel3D_t`) and override:

| Procedure | Purpose |
|-----------|---------|
| `SetNumberOfVariables` | Set `this%nvar` |
| `twopointflux2d(sL, sR)` | Return the physical EC two-point flux (nvar x 2 array) |
| `riemannflux2d(sL, sR, dsdx, nhat)` | Return the surface Riemann flux (nvar array). Use LLF for symmetric dissipation. |
| `entropy_func(s)` | Return the scalar entropy $\eta(\mathbf{s})$ for diagnostics |
| `SetMetadata` | (Optional) name and units for each variable |

The base class handles the split-form volume integral, metric averaging, surface correction, time integration (inherited from `DGModel2D_t`), and I/O.

**Example:** `ECAdvection2D_t` implements entropy-conserving linear advection with:
- `twopointflux2d`: arithmetic mean $f^\# = \mathbf{a}(s_L + s_R)/2$
- `riemannflux2d`: Local Lax-Friedrichs $f_\text{LLF} = \tfrac{1}{2}\left[u_n(s_L+s_R) - |\mathbf{a}|(s_R-s_L)\right]$
- `entropy_func`: $\eta = s^2/2$


## References

1. Gassner, G. J., Winters, A. R., & Kopriva, D. A. (2016). Split form nodal discontinuous Galerkin schemes with summation-by-parts property for the compressible Euler equations. *Journal of Computational Physics*, 327, 39-66.

2. Winters, A. R., Kopriva, D. A., Gassner, G. J., & Hindenlang, F. (2021). Construction of Modern Robust Nodal Discontinuous Galerkin Spectral Element Methods for the Compressible Navier-Stokes Equations. *Lecture Notes in Computational Science and Engineering*, Springer.

3. Ranocha, H., Schlottke-Lakemper, M., Winters, A. R., Faulhaber, E., Chan, J., & Gassner, G. J. (2022). Adaptive numerical simulations with Trixi.jl: A case study of Julia for scientific computing. *Proceedings of the JuliaCon Conferences*.
