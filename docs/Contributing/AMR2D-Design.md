# Design Assessment: Conforming Adaptive Mesh Refinement (AMR) for 2D Models

**Status:** Proposal / scoping document — no implementation yet.
**Scope:** 2D only. Conforming quadrilateral h-refinement (no mortars, no hanging
nodes), spectral (Legendre modal) refinement indicator, CPU and GPU backends,
MPI-parallel.

This document records a full audit of the current SELF architecture as it
relates to adaptivity, the key design decisions (with alternatives and
trade-offs), the proposed new components, and a phased work breakdown. Several
decisions need maintainer sign-off before implementation begins; these are
collected in [Open questions](#12-open-questions-for-maintainers).

---

## 1. Requirements

1. A new mesh type supporting adaptive refinement/coarsening of 2D quad meshes.
2. **Conforming elements only** — every interior side is shared by exactly two
   elements; no mortared/nonconforming interfaces, no hanging nodes.
3. New model type(s) that leverage the AMR mesh (starting from `DGModel2D`).
4. Refinement indicator built from **Legendre polynomial projection
   coefficients** of the element-local solution: refine when the modal spectrum
   is "too shallow" (slow decay) or when high-order modes carry too much
   energy.
5. CPU and GPU backends.

---

## 2. Current-state audit

What exists today (references are to files at the time of writing):

| Capability | Status | Where |
|---|---|---|
| Conforming static quad mesh, 1:1 side pairing, `flip ∈ {0,1}` | ✅ | `Mesh2D_t`, `src/SELF_Mesh_2D_t.f90:102-136` |
| Side-connectivity reconstruction from corner-node pairs (hash matching) | ✅ | HOHQMesh reader, `src/SELF_Mesh_2D_t.f90:886-941` |
| Side orientation recovery (incl. MPI) | ✅ | `RecalculateFlip`, `src/SELF_Mesh_2D_t.f90:1089-1235` |
| Geometry/metric generation from any conforming quad mesh | ✅ | `GenerateFromMesh_SEMQuad`, `src/SELF_Geometry_2D.f90:110-145` |
| Block domain decomposition + halo exchange | ✅ | `SELF_DomainDecomposition_t.f90`, `SideExchange` in `src/SELF_MappedScalar_2D_t.f90:221-282` |
| Nodal interpolation between grids (`GridInterp`, CPU + hipBLAS GPU) | ✅ | `src/SELF_Scalar_2D_t.f90:176-198`, `src/gpu/SELF_Scalar_2D.f90:177-190` |
| Legendre polynomial evaluation | ✅ | `LegendrePolynomial`, `src/SELF_Quadrature.f90:360` |
| Nodal→modal (Legendre) transform, spectral indicator | ❌ none | — |
| L2 projection / coarsening operators | ❌ none | — |
| Mortars, nonconforming sides, h/p-refinement, element trees | ❌ none | — |
| Dynamic `nElem` (reallocatable fields, resizable GPU mirrors) | ❌ none | all data objects sized at `Init` |

Load-bearing consequences:

- **Everything downstream of the mesh already works for *any* conforming quad
  mesh.** `sideInfo`, `RecalculateFlip`, `GenerateFromMesh_SEMQuad`,
  `SideExchange`, and the DG operators make no structured-grid assumptions.
  If AMR produces a valid conforming quad mesh, the entire existing solver
  stack runs on it unchanged. This is the single biggest reuse win and drives
  the "regenerate the mesh each adaptation" design below.
- **Nothing is sized dynamically.** `Mesh2D`, `SEMQuad`, and all
  `MappedScalar2D`/`MappedVector2D` fields (plus their `_gpu` device mirrors,
  hand-rolled `c_ptr` + `hipMalloc`) are allocated once for a fixed `nElem`.
  Adaptation therefore implies free/re-init of mesh, geometry, and model
  storage — acceptable if adaptation is infrequent relative to time steps.
- **BC bookkeeping is element-indexed.** `MapBoundaryConditions_DGModel2D_t`
  (`src/SELF_DGModel2D_t.f90:453-530`) caches `(element, side)` lists per BC,
  with device copies. These must be rebuilt after every adaptation.

---

## 3. Central design decision: conforming refinement without mortars

Isotropic quad splitting (1 → 2×2) creates hanging nodes wherever a refined
element neighbors an unrefined one. DGSEM codes usually reconcile this with
mortars; those are excluded here. To stay conforming we must close the
refinement boundary with **transition templates** — extra quads that connect
trisected/bisected edges to whole edges.

### 3.1 Why 3-refinement (1 → 3×3), not 2-refinement (1 → 2×2)

A counting argument (face–vertex incidence) shows a quad with a single
bisected edge **cannot** be partitioned into quads without either adding nodes
to its other edges (which propagates refinement to neighbors) or leaving a
non-quad face. This is the classical result behind Schneiders' work: 2-refinement
in all-quad meshes requires sheet/buffer propagation and has no local
single-element closure.

**3-refinement does have local closures.** With edges trisected, a quad with
one refined edge splits into 4 quads with 2 interior nodes:

```
D ______________ C          D ______________ C
 |              |            |              |
 |              |            |   P1____P2   |     Q1 = (A,T1,P1,D)
 |              |    --->    |   |      |   |     Q2 = (T1,T2,P2,P1)
 |              |            |   |      |   |     Q3 = (T2,B,C,P2)
 |______________|            |___|______|___|     Q4 = (P1,P2,C,D)
A                B          A    T1     T2   B

P1 = (1/3, 1/2), P2 = (2/3, 1/2) in the element reference square
```

All four cells are convex quads; the trisected edge (T1, T2) mates exactly with
the three edge-children of the refined neighbor. This is the standard
template-based 3-refinement of Schneiders (*"Refining quadrilateral and
hexahedral element meshes"*, 1996/2000), also used in Sandia meshing tools.

**Recommendation:** ternary (3-refinement) quadtree per element. Each
refinement level divides h by 3 (a refined element becomes 9 children). The
factor-9 growth per level is the price of the conforming constraint; it should
be stated clearly to users (mortars would give factor 4 — see §3.4).

### 3.2 Minimal template set via flag promotion

Rather than implementing the full Schneiders template catalogue (1 edge, 2
adjacent, 2 opposite, 3, 4 edges), we use **flag promotion** to need only two
cases:

1. **Full refinement** — 3×3 split (9 children), when the element itself is
   flagged.
2. **Single-edge transition** — the 4-quad template above (up to rotation),
   when exactly one edge abuts a refined neighbor.
3. Trivial extra case: **two opposite refined edges** → 3 vertical/horizontal
   strip quads (both trisected edges mate; no interior nodes needed).

Promotion rule, iterated to a fixed point: any element with **two or more
refined edges that are not an opposite pair** is promoted to full refinement.
Iteration converges (worst case: uniform refinement). Combined with standard
**level grading** (edge-neighbor tree levels differ by at most 1, enforced by
flag smoothing before template generation), transitions only ever bridge one
level.

In 2D, corner-only (diagonal) adjacency needs no special treatment: trisection
introduces new nodes only on edge interiors, so a diagonal neighbor at a
different level shares just the original corner node — conforming by
construction.

### 3.3 Tree + regenerated leaf mesh (the "two-layer" mesh design)

The authoritative adaptive structure is a **forest of ternary quadtrees**, one
tree per element of the user-supplied base mesh ("root elements"). Transition
quads are **ephemeral**: they are *not* tree nodes, but are regenerated from
the red (fully-refined) tree at every adaptation, in the closure pass.

Each adaptation cycle produces a complete, ordinary conforming `Mesh2D` (the
"leaf mesh") from scratch:

1. Evaluate the spectral indicator per element → refine/coarsen flags.
2. Smooth flags: enforce level grading ≤ 1; apply promotion (§3.2);
   coarsen only complete 9-sibling groups whose removal keeps grading valid.
3. Update the red trees (split/merge tree nodes).
4. Emit leaf quads: red leaves + closure transitions; assign global corner-node
   IDs (hash on parametric position within root, quantized) and global element
   ordering (roots in base-mesh order, depth-first within each tree —
   deterministic, locality-preserving).
5. Build `sideInfo` by corner-pair hash matching (reuse the algorithm at
   `src/SELF_Mesh_2D_t.f90:886-941`, lifted into a shared helper), then
   `RecalculateFlip`.
6. Regenerate `SEMQuad` geometry via the existing `GenerateFromMesh` path.

Why regeneration instead of incremental connectivity surgery: it reuses the
entire, already-MPI-correct construction path; it is far easier to verify; and
adaptation cost is amortized over many time steps. Incremental updates are a
later optimization if profiling demands it.

**Geometry of children.** Every leaf carries its bilinear map (4 corner
points) in the **root element's reference square**, plus `(rootElem, level)`.
For `nGeo = 1` base meshes, child `nodeCoords` are the bilinear image of the
template corners — exact. For curved elements (`nGeo > 1`, e.g. HOHQMesh ISM),
child geometry nodes are sampled from the parent's isoparametric map at the
template-mapped points, preserving boundary geometry to the parent's order.
Positivity of each child Jacobian is checked at generation time (transition
quads are non-affine; a validity check is mandatory).

**Boundary conditions.** A leaf edge lying on a root-element boundary edge
inherits that edge's `sideInfo(5)` BC id from the base mesh; interior leaf
edges get 0. `BCType`/`BCNames` are copied from the base mesh unchanged.

### 3.4 Trade-off vs. mortars (flagged for maintainer awareness)

The conforming constraint buys: unchanged DG kernels, unchanged `SideExchange`,
no new interface operators in the time loop, and provable conservation at faces
comes for free (watertight faces, single shared Riemann flux). It costs:
factor-9 refinement granularity, transition elements with reduced quality
(bounded, fixed template shapes — their reference-space aspect/skew is known a
priori), a more complex mesh-generation step, and coarsening constrained to
full sibling groups. Mortar-based 2:1 AMR would invert these trade-offs but
touches the numerics (mortar projections in the flux path) — explicitly out of
scope per the requirements.

---

## 4. Spectral refinement indicator (Legendre modal analysis)

### 4.1 Nodal → modal transform

New machinery (currently absent — only `LegendrePolynomial` at
`src/SELF_Quadrature.f90:360` exists as a building block). For control points
`ξ_i` with quadrature weights `w_i` (Gauss or Gauss–Lobatto, degree N), define
the discrete Legendre coefficients of a nodal field `u_i = u(ξ_i)`:

```
â_k = (1/γ_k) Σ_i w_i u_i L_k(ξ_i),   γ_k = Σ_i w_i L_k(ξ_i)²
```

Using the *discrete* norms `γ_k` (rather than the analytic `2/(2k+1)`) makes
the transform exactly invertible on the nodal space for both Gauss and
Gauss–Lobatto points (for Gauss, `γ_k = 2/(2k+1)` exactly for k ≤ N; for GLL,
`γ_N = 2/N`). In matrix form `â = P u` with `P = Γ⁻¹ Lᵀ W`; the 2D transform is
the tensor application `Â = P U Pᵀ` per element and variable — i.e. exactly the
two-pass batched-GEMM shape already used by `GridInterp`/derivative operators
(`self_blas_matrixop_dim{1,2}_2d`, `src/gpu/SELF_GPUBLAS.f90`).

Implementation home: add `pMatrix(0:N,0:N)` (nodal→modal) and its inverse to
`Lagrange_t` (`src/SELF_Lagrange_t.f90`), built in `Init_Lagrange_t`, with
`pMatrix_gpu` mirrors added in `src/gpu/SELF_Lagrange.f90` (and `src/apu/`)
following the existing pattern for `dgMatrix` etc. This also unlocks future
modal filtering/de-aliasing work.

### 4.2 Element indicator

Per element and indicator variable, form directional spectra from the 2D modal
matrix `Â_{mn}`:

```
qξ_m = sqrt( Σ_n Â_{mn}² ),    qη_n = sqrt( Σ_m Â_{mn}² )
```

Two criteria, matching the stated requirement:

1. **Spectral decay rate ("too shallow")** — Mavriplis-style least-squares fit
   `q_m ≈ C e^(−σ m)` over the upper half of the spectrum
   (`m = max(1, N/2) … N`, skipping near-zero modes below a floor to avoid
   log-of-noise). Small σ ⇒ under-resolved.
2. **Tail energy ("too much energy in high modes")** — Persson–Peraire-style
   ratio using the top two modes (two, to be robust to odd/even cancellation):
   `η = (q_N² + q_{N−1}²) / Σ_m q_m²`.

Decision logic with hysteresis (two threshold pairs, so elements don't flip
refine/coarsen every cycle):

```
refine  if  min(σξ, ση) < σ_refine   or  max(ηξ, ηη) > η_refine
coarsen if  min(σξ, ση) > σ_coarsen  and max(ηξ, ηη) < η_coarsen
          (and level > 0, and all 8 siblings also vote coarsen)
```

Defaults to be calibrated during Phase 1 against manufactured cases (smooth
Gaussian: no flags; near-discontinuity: flags); reasonable starting points are
`σ_refine ≈ 1`, `η_refine ≈ 10⁻⁴ · (N+1)⁻⁴`-scaled, with coarsen thresholds a
factor ~4 apart from refine thresholds. All thresholds, the adapt cadence, the
max refinement level, and the set of indicator variables are user-facing
parameters; the indicator variable set defaults to the `nstepped` prognostic
variables (max over variables), overridable per model (e.g. pressure only for
`LinearEuler2D`).

Note the indicator is computed in **reference space**, so transition-element
skew does not bias it.

### 4.3 Backends

- **CPU:** `do concurrent` over elements; two small matmuls + O(N²) reduction
  per element. Negligible cost relative to one RHS evaluation.
- **GPU:** modal transform via the existing hipBLAS batched pattern with
  `pMatrix_gpu`; one custom kernel (`SpectralIndicator_2D_gpu` in a new
  `src/gpu/SELF_AMR.cpp`) computes `q`, the decay fit, tail ratio, and writes
  per-element `int` flags; a single small D2H `hipMemcpy` brings flags to the
  host, where all adaptation logic runs. No new host↔device traffic patterns
  inside time stepping — the indicator runs only at adaptation checkpoints.

---

## 5. Solution transfer between meshes

Every leaf stores its bilinear map into the root reference square, so any two
generations of leaves relate through ancestor polynomial spaces. Rules:

- **Unchanged leaf:** direct copy (element index remap only).
- **Refinement (parent → children):** evaluate the parent interpolant at each
  child's control points mapped into the parent's reference square. This is
  exact (children's polynomial spaces contain the restricted parent
  polynomial), so it is conservative and cannot create new extrema beyond
  interpolation effects. Implemented as per-template precomputed interpolation
  matrices (`(N+1)² × (N+1)²`, one per child position per template — a small,
  static set), applied as batched GEMMs grouped by template class.
- **Coarsening (children → parent):** discrete **L2 projection**: for parent
  mode `(m,n)`, `â^P_{mn} = (1/γ_m γ_n) Σ_children Σ_ij w_i w_j |J_c(ξ_ij)|/|J_P|
  · u_c(ξ_ij) L_m(ξ^P) L_n(η^P)` with child quadrature points mapped into the
  parent reference square. L2 projection preserves the element mean (mode 0)
  exactly for affine subcells, i.e. the transfer is conservative up to
  quadrature/geometry error on curved elements.
- **Transition regions:** transitions are ephemeral (§3.3), so an old
  transition group may be replaced by a full 3×3 (or revert to the parent).
  Transfer routes through the common parent space: L2-project the old covering
  leaves onto the parent polynomial, then interpolate/project to the new
  leaves. The parent-space detour loses only what the parent space cannot
  represent — acceptable at resolution boundaries, and conservative.

A `MeshTransfer2D` object encapsulates old→new leaf correspondence and the
matrix set; it is built on the host during adaptation and (Phase 4) applied on
device for leaves resident on the GPU.

**Numerical-correctness note (per CLAUDE.md §3):** solution transfer changes
solution values by design; it is an approximation-theoretic operation
(interpolation/L2 projection), not a change to the discretization. Entropy may
increase slightly at coarsening (projection is L2-, not entropy-stable);
this is standard for AMR and will be measured in validation. If needed later,
entropy-bounded projection is a contained follow-up.

---

## 6. New components and types

New files (no existing modules renamed/moved; public APIs untouched):

```
src/SELF_ModalTransform.f90        ! pMatrix construction + spectra helpers (used by Lagrange_t)
src/SELF_QuadTemplates_2D.f90      ! template geometry tables (corner coords, edge-child maps)
src/SELF_AMRMesh_2D_t.f90          ! AMRMesh2D_t: forest, flags, leaf-mesh generation
src/SELF_MeshTransfer_2D_t.f90     ! old→new transfer operator
src/SELF_AMRController_2D_t.f90    ! indicator + adapt orchestration for class(DGModel2D_t)
src/cpu/SELF_AMRMesh_2D.f90        ! thin backend extensions (pattern-matching existing cpu/)
src/gpu/SELF_AMRMesh_2D.f90
src/gpu/SELF_AMR.cpp               ! SpectralIndicator_2D_gpu, transfer-apply kernels
                                   ! + interfaces added to src/gpu/SELF_GPUInterfaces.f90
```

### 6.1 `AMRMesh2D_t` (containment, not inheritance)

```fortran
type :: AMRMesh2D_t
  type(Mesh2D) :: base       ! user-supplied root mesh (never modified)
  type(Mesh2D) :: leaf       ! regenerated conforming mesh; models point here
  ! forest state
  integer :: nRoots, nLeaves, maxLevel
  integer, allocatable :: leafRoot(:)      ! root element of each leaf
  integer, allocatable :: leafLevel(:)     ! tree level (transitions: parent level)
  integer, allocatable :: leafKind(:)      ! red child / transition template id + orientation
  real(prec), allocatable :: leafMap(:,:,:) ! (2,4,nLeaves) corners in root ref square
  type(...), allocatable :: trees(:)        ! per-root ternary quadtree (packed arrays)
contains
  procedure :: Init      ! from an existing Mesh2D (structured or file-read)
  procedure :: Adapt     ! flags in -> new leaf mesh + transfer maps out
  procedure :: Free
endtype
```

Containment is chosen over `extends(Mesh2D)` because `Mesh2D%Init` has a
different lifecycle (regenerated repeatedly) and because models only need a
`type(Mesh2D), pointer` — pointing at `amrMesh%leaf` keeps `DGModel2D_t`,
`SideExchange`, geometry, and I/O completely unaware of AMR.

### 6.2 Model integration: controller + thin AMR model types

Adaptation is model-agnostic given `class(DGModel2D_t)`; physics lives in
overridden `flux2d`/`riemannflux2d`/BCs and is untouched. An
`AMRController2D` holds the `AMRMesh2D`, thresholds, cadence, and performs:

```
AdaptModel(model):
  1. indicator (device or host)             -> flags
  2. amrMesh%Adapt(flags)                   -> new leaf Mesh2D + MeshTransfer2D
  3. regenerate SEMQuad geometry            (existing GenerateFromMesh)
  4. capture old solution; Free + re-Init model fields (solution, flux, dSdt, ...)
     re-AssociateGeometry; rebuild BC lists (MapBoundaryConditions; device copies)
  5. apply MeshTransfer2D to fill solution  (+ UpdateDevice on GPU builds)
  6. refresh decomposition (MPI: repartition + redistribute; §7)
```

New model types are thin wrappers that own a controller and trigger it from the
existing hook system (`PostStepHook`, every `nStepsPerAdapt` steps — between RK
steps, never inside a stage), e.g.:

```fortran
type, extends(advection_diffusion_2d) :: amr_advection_diffusion_2d
  type(AMRController2D) :: amr
contains
  procedure :: PostStepHook => AdaptHook_amr_advection_diffusion_2d
endtype
```

Pilot models: `amr_advection_diffusion_2d` and `AMRLinearEuler2D` (CPU + GPU
variants following the existing backend file pattern). Metadata, equation
parsers, and integrator state carry over; `workSol` is rebuilt (zeroed) at
adaptation, which is correct because adaptation happens between complete steps.

**Time step:** explicit CFL limit shrinks with the smallest h (factor 3 per
level). v1 keeps the user-controlled global `dt` contract: the controller
exposes `hMin()` and the examples rescale `dt` accordingly on adaptation. Local
time stepping is explicitly out of scope.

### 6.3 I/O

`Write_DGModel2D_t` writes per-element arrays and already handles arbitrary
conforming meshes; output continues to work, but element counts change between
files. Add the leaf mesh (or `leafRoot`/`leafLevel`/`leafMap`) to the output so
post-processing (pyself) can reconstruct geometry per snapshot. Existing
outputs of static-mesh models are untouched (no reference-file changes).

---

## 7. MPI strategy

- **Ordering & partition:** global leaf order = base-mesh root order,
  depth-first within each tree (deterministic; preserves locality similarly to
  an SFC). Reuse `GenerateDecomposition`'s contiguous block partition
  (`src/SELF_DomainDecomposition_t.f90:129-190`) on the new `nLeaves`. Equal
  element count is a valid balance measure because all elements share the same
  polynomial degree N.
- **Consistent flags:** refine/coarsen flags for ghost-adjacent elements must
  agree across ranks before smoothing. v1: `MPI_Allgatherv` of flags (one int
  per element) so every rank runs identical smoothing/promotion — deterministic
  and simple; the tree metadata is tiny compared to field data. Rank-local
  smoothing with ghost-layer exchange is a scalability follow-up.
- **Redistribution:** old and new partitions are both contiguous over the same
  deterministic ordering, so each rank computes exact send/recv ranges
  (`MPI_Alltoallv` on packed solution payloads keyed by parent/root). Transfer
  operators that need children resident on one rank (coarsening groups,
  transition groups) are applied on the receiving side after redistribution of
  parent-space coefficients.
- **Constraint compliance:** all collectives happen inside the adaptation
  checkpoint (between time steps), never inside `CalculateTendency` or RK
  stages, per CLAUDE.md §5. Existing halo-exchange patterns (`MPIExchangeAsync`
  tags on global side IDs, `src/SELF_MappedScalar_2D_t.f90:121-172`) are reused
  untouched — regeneration produces fresh valid global side IDs.
- **MPI testing:** required on ≥ 2 ranks (see §9).

---

## 8. GPU strategy

Division of labor: **numerics on device, adaptation logic on host.**

Runs on device (new or reused):
- Spectral indicator: `pMatrix_gpu` batched GEMMs + `SpectralIndicator_2D_gpu`
  reduction kernel → per-element flags (single small D2H copy).
- Solution transfer apply (Phase 4): children/parent gather kernels batched by
  template class, using precomputed matrices resident on device.
- Everything already on device (RHS, exchange, RK updates) — unchanged.

Runs on host:
- Tree updates, flag smoothing/promotion, template emission, corner hashing,
  `RecalculateFlip`, partition computation.

**Memory churn:** each adaptation frees and re-`hipMalloc`s every field's
device mirror (existing `Init`/`Free`/`UpdateDevice` machinery). v1 accepts
this (adaptation is infrequent). Phase 5 optimization: capacity-based
allocation (`nElemCapacity = growthFactor * nElem`) inside the data-object
`Init`s so most adaptations skip reallocation — contained change to the gpu
backend `Init` routines, no interface change.

GPU-aware MPI paths are unaffected: redistribution stages through host in v1
(adaptation checkpoint, not hot path).

The `apu/` backend gets `pMatrix` plumbing for `Lagrange` (5-file backend is
already partial; full APU AMR support tracked separately).

---

## 9. Testing & validation plan

Unit (serial, both precisions, tolerances precision-aware per existing
convention, e.g. `test/mappedscalargradient_2d_linear.f90:52-56`):
- `modaltransform_roundtrip`: `P⁻¹ P u = u` to roundoff; Gauss and GL nodes.
- `spectralindicator_smooth` / `spectralindicator_front`: manufactured fields
  (low-degree polynomial ⇒ no refine flags; `tanh` front ⇒ flags exactly in
  front-crossing elements).
- `quadtemplates_validity`: every template/orientation yields positive
  Jacobians and correct edge-child adjacency.
- `amrmesh2d_conforming`: after randomized refine/coarsen cycles, every
  interior side matches exactly one neighbor (audit `sideInfo` both ways),
  flips verified via `RecalculateFlip` idempotence, mesh is watertight
  (shared-edge node coincidence).
- `meshtransfer_exactness`: refine transfer exact for polynomial data of
  degree ≤ N; coarsen preserves `∫ u dA` to roundoff (affine geometry).

Integration:
- `amr_advection_2d_gaussian`: advect a Gaussian; AMR tracks it; error vs. a
  uniform-fine reference within tolerance and conservation error at roundoff.
- `amr_lineareuler2d_wave`: acoustic pulse; entropy non-increase across the run
  (existing entropy machinery), refinement follows the wavefront.
- `amr_*_mpi` variants on 2 ranks (registered via `add_mpi_fortran_tests`).
- GPU CI: same integration tests in the Buildkite MI210/V100 pipelines;
  CPU/GPU agreement of indicator flags and post-adapt solutions to tolerance.

No existing reference outputs are modified; all current regression tests must
remain green (AMR is purely additive).

---

## 10. Work breakdown & phasing

| Phase | Content | New/changed code (est.) |
|---|---|---|
| 0. Modal infrastructure | `pMatrix` in `Lagrange_t` + GPU/APU mirrors, `SELF_ModalTransform.f90`, indicator math + unit tests | ~600 LOC |
| 1. Serial CPU AMR core | templates, forest, leaf-mesh generation (connectivity-hash reuse), `AMRMesh2D_t`, transfer operators, conforming/validity/transfer tests | ~2,500 LOC — the long pole |
| 2. Model integration (CPU) | `AMRController2D`, field re-init + BC rebuild path, `amr_advection_diffusion_2d`, `AMRLinearEuler2D`, integration tests, example programs | ~1,200 LOC |
| 3. MPI | flag allgather, repartition + `Alltoallv` redistribution, 2-rank tests | ~800 LOC |
| 4. GPU | indicator kernels + interfaces, device transfer apply, BC device-list rebuild, GPU CI tests | ~1,000 LOC |
| 5. Hardening/perf | capacity-based device allocation, adaptation-cost profiling, threshold calibration docs, user docs (this page → user guide) | ~500 LOC |

Rough total ~6–7k LOC plus tests. Phases 0–2 deliver a usable serial CPU
capability; 3 and 4 are independent of each other once 2 lands.

---

## 11. Compliance with repository constraints (CLAUDE.md)

- No change to formulation, discretization order, basis, quadrature, nodal
  ordering, flux routines, or time integrators — AMR acts *between* steps.
- No renamed/moved modules; all existing public interfaces preserved; new
  functionality is additive (`Lagrange_t` gains members, no signature changes).
- No collectives inside time-stepping loops; halo patterns unchanged.
- No external dependencies (forest is home-grown; no p4est).
- Element-local operations remain element-local; the leaf mesh is a plain
  `Mesh2D`.
- Fortran 2008 only; formatting via `fprettify`.

---

## 12. Open questions for maintainers

1. **Accept 3-refinement (factor 9/level) as the price of "conforming, no
   mortars"?** 2-refinement provably has no local conforming all-quad closure;
   the alternatives are mortars (excluded) or propagating sheet refinement.
2. **Adaptation cadence & placement:** `PostStepHook` every `nStepsPerAdapt`
   steps — acceptable, or should adaptation align with IO intervals in
   `ForwardStep`?
3. **`dt` policy on adaptation:** keep user-controlled global `dt` (examples
   rescale via `hMin()`), or auto-rescale inside the controller?
4. **Indicator variables:** default to all prognostic (`nstepped`) variables
   with max-combining — confirm, and confirm per-model override hook.
5. **Curved base meshes (`nGeo > 1`) in v1**, or restrict v1 to `nGeo = 1`
   (structured + straight-sided ISM) and defer curved-element child geometry?
6. **Coarsening policy:** strict 9-sibling groups only (proposed), or also
   allow transition-region relaxation when the refined neighbor coarsens?
7. **Output format for time-varying meshes:** embed leaf mesh per snapshot
   (proposed) vs. separate mesh series files — coordinate with pyself.

---

## 13. References

- R. Schneiders, "Refining quadrilateral and hexahedral element meshes,"
  *Numerical Grid Generation in Computational Field Simulations*, 1996; and
  "Octree-based hexahedral mesh generation," *IJCGA*, 2000 (3-refinement
  templates).
- P.-O. Persson, J. Peraire, "Sub-Cell Shock Capturing for Discontinuous
  Galerkin Methods," AIAA 2006-112 (modal tail-energy sensor).
- C. Mavriplis, "Adaptive mesh strategies for the spectral element method,"
  *CMAME* 116 (1994) (spectral decay-rate estimator).
- D. A. Kopriva, *Implementing Spectral Methods for Partial Differential
  Equations* (algorithms referenced throughout `SELF_Lagrange_t.f90` /
  `SELF_Quadrature.f90`).
