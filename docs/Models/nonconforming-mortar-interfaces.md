# Nonconforming (Mortar) Interfaces in 2-D

SELF supports 2:1 nonconforming quadrilateral interfaces — "mortars" — in its 2-D
discontinuous Galerkin models. At a mortar interface, one **big** element edge coincides
with the two half-edges of two refined (**small**) neighbor elements. This is the
standard interface treatment used for quadtree-style adaptive meshes.

## Method

SELF implements the projection-based DGSEM mortar method (Maday, Mavriplis & Patera,
1989; Kopriva, 1996):

1. **Big → small (solution traces).** The big side's edge trace, a degree-\(N\)
   polynomial on \([-1,1]\), is restricted to the sub-edges \([-1,0]\) and \([0,1]\)
   with precomputed restriction matrices (`Lagrange % mortarR`). Restriction of a
   polynomial is exact.
2. **Riemann solve on the small sides.** Each small side computes its surface flux with
   its own outward normal and boundary metric, using the model's `riemannflux2d` —
   no model code changes are involved.
3. **Small → big (flux integrands).** The big side's surface-flux integrand
   \(f^* \cdot \hat{n}\, n_{scale}\) is replaced by the \(L^2\) projection of the two
   small sides' integrands (`Lagrange % mortarP`, applied by `MortarFluxCollect`), with
   the sign flipped for the opposing normal.

The projection matrices are built with the exact one-dimensional mass matrix, evaluated
with an internal Gauss rule that is exact for the polynomial integrands. Two discrete
identities follow, both enforced by the `mortarprojection_identity` test:

- **Forward–backward consistency:** restricting a trace to both sub-edges and
  projecting back recovers it exactly, for Legendre–Gauss *and*
  Legendre–Gauss–Lobatto control points.
- **Conservation:** the big side's discrete surface integral equals minus the sum of
  the small sides' discrete surface integrals, to roundoff. The mortar interface
  neither creates nor destroys the conserved quantities.

For models with parabolic terms (`gradient_enabled`), the solution-gradient traces are
exchanged through the same restriction/projection operators and the interface gradient
is the average of projected traces (project-then-average BR1).

## Mesh representation

Mortar connectivity lives on `Mesh2D` alongside the conforming `sideInfo`:

- `mesh % nMortars` — number of mortar interfaces (0 for conforming meshes; all
  mortar code paths are skipped and conforming meshes execute exactly as before).
- `mesh % mortarInfo(1:8, 1:nMortars)` — for each mortar: the big element and local
  side, the two small elements with their `10*side + flip` codes (sub-edge 1 covers
  big-edge coordinate \([-1,0]\), sub-edge 2 covers \([0,1]\)), and a global side id
  per sub-edge used for MPI message tags. The table is replicated on all ranks;
  element ids are global.

Sides that participate in a mortar carry `sideInfo(3) = 0` (no conforming neighbor)
and `sideInfo(5) = 0` (no boundary condition), so the conforming side exchange,
flip recovery, and boundary-condition mapping all skip them automatically;
`sideInfo(1)` records the mortar index.

### Building a nonconforming mesh

Two built-in constructors create nonconforming meshes. `SimpleMortarMesh` is the
smallest useful example — one 2dx-by-2dx element joined to two dx-by-dx elements
across a single mortar interface — and `DoubleMortarMesh` is a six-element,
two-mortar configuration whose second mortar includes a reversed-orientation
(`flip = 1`) small element, used by the validation suite to exercise every trace
reorientation and MPI message pattern:

```fortran
type(Mesh2D),target :: mesh
integer :: bcids(1:4)

bcids = [SELF_BC_RADIATION,SELF_BC_RADIATION, &
         SELF_BC_RADIATION,SELF_BC_RADIATION] ! south, east, north, west
call mesh%SimpleMortarMesh(0.1_prec,bcids)
```

Everything downstream is unchanged from a conforming workflow: generate the geometry
with `geometry % GenerateFromMesh(mesh)`, initialize a model with
`modelobj % Init(mesh,geometry)`, and call `ForwardStep`. Meshes with many mortars can
be constructed the same way this constructor does it: fill `nodeCoords`,
`globalNodeIDs`, and `sideInfo` for all elements, then populate `nMortars` and
`mortarInfo` (see `SimpleMortarMesh_Mesh2D_t` in `src/SELF_Mesh_2D_t.f90` as a
template).

## Backends and parallelism

- **CPU and GPU** are both supported. On GPU builds, traces are staged, reoriented,
  restricted, and projected entirely in device memory (`src/gpu/SELF_Mortar.cpp`);
  the operator matrices are resident on the device.
- **MPI:** mortar interfaces may straddle rank boundaries in any way. The big-side
  rank exchanges traces point-to-point with each remote small-side rank in a single
  message round per exchange, using the sub-edge global side ids as tags (the same
  convention as the conforming `SideExchange`). On GPU builds, messages are posted on
  device memory (GPU-aware MPI), matching the conforming exchange.

## Limitations

- 2-D only; 3-D (hex faces with four sub-faces and general orientation) is future
  work.
- 2:1 refinement ratio only, and h-nonconformity only (all elements share one
  polynomial degree; no p-mortars).
- The entropy-conserving split-form models (`ECDGModel2D` and its descendants) do
  **not** support mortar meshes: the plain \(L^2\) projection would break their
  provable entropy estimate, and entropy-stable mortar operators (Friedrich et al.,
  2018) are not implemented. Their initialization raises an error on a mesh with
  `nMortars > 0`.
- `mortarInfo` is replicated across ranks, which is appropriate for the modest mortar
  counts of static nonconforming meshes.

## Tests

| Test | What it checks |
|---|---|
| `mortarprojection_identity` | Operator identities: \(\sum_k P_k R_k = I\) and discrete conservation, Gauss and Gauss–Lobatto, to a near-roundoff tolerance |
| `mappedscalarmortarexchange_2d_linear` (+`_mpi`) | Exact trace recovery in `extBoundary` on all sides of both mortars of the `DoubleMortarMesh` (including the `flip = 1` sub-edge), with distinct degree-N fields per variable |
| `mappedvectordgdivergence_2d_mortar` (+`_mpi`) | DG divergence of \(\vec{f}=(x,y)\) equals 2 to roundoff across both mortars; per-interface conservation defect at roundoff |
| `lineareuler2d_mortar_soundwave` (+`_mpi`) | Acoustic pulse crossing the interface: stability (entropy non-increase) and NaN-free evolution |
| `advection_diffusion_2d_mortar` (+`_mpi`) | Parabolic (BR1) terms across mortars: the solution-gradient mortar exchange and project-then-average interface gradient, with entropy non-increase |
| `ec_advection_2d_mortar_guard` | EC split-form models abort on nonconforming meshes (registered `WILL_FAIL`) |

The `_mpi` variants run on two ranks. On the `DoubleMortarMesh`, mortar 1 splits
big/small elements across the rank boundary while mortar 2 places both small
elements remote from the big element's rank, so every mortar message pattern —
including independent delivery of the big-side trace to each sub-edge — is
exercised.

## References

- Y. Maday, C. Mavriplis, A. Patera, "Nonconforming mortar element methods:
  application to spectral discretizations," *Domain Decomposition Methods*, 1989.
- D. A. Kopriva, "A conservative staggered-grid Chebyshev multidomain method for
  compressible flows. II. A semi-structured method," *JCP* 128 (1996).
- L. Friedrich, A. R. Winters, D. C. Del Rey Fernández, G. J. Gassner, M. H.
  Carpenter, "An entropy stable h/p non-conforming discontinuous Galerkin method with
  the summation-by-parts property," *JSC* 77 (2018).
