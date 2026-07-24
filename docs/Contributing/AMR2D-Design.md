# Design Assessment: Adaptive Mesh Refinement (AMR) for 2D Models

**Status:** Proposal / scoping document — no implementation yet.
**Scope:** 2D quadrilateral meshes, quadtree (2:1) h-refinement built on the
implemented mortar interface support, spectral (Legendre modal) refinement
indicator, CPU and GPU backends, MPI-parallel.

**Revision note:** the original version of this document scoped AMR under a
"conforming elements only, no mortars" requirement, which forced a
template-based 3-refinement design (9 children per level, transition-element
closures). That constraint has since been dropped: 2:1 nonconforming mortar
interfaces are now implemented in SELF (see
[Nonconforming (Mortar) Interfaces](../Models/nonconforming-mortar-interfaces.md)
and the companion [mortar design document](Mortar2D-Design.md)). This document
has been rewritten around the mortar-based quadtree design, which is the
industry-standard AMR structure (FLEXI, Trixi.jl, Nektar++). The retired
conforming-template analysis is preserved in the git history of this file for
reference.

---

## 1. Requirements

1. A new mesh type supporting adaptive refinement/coarsening of 2D quad meshes,
   built on the existing 2:1 mortar interface machinery.
2. New model type(s) that leverage the AMR mesh (starting from `DGModel2D`).
3. Refinement indicator built from **Legendre polynomial projection
   coefficients** of the element-local solution: refine when the modal spectrum
   is "too shallow" (slow decay) or when high-order modes carry too much
   energy.
4. CPU and GPU backends.

---

## 2. Current-state audit

What exists today:

| Capability | Status | Where |
|---|---|---|
| Conforming quad mesh, side pairing, geometry, halo exchange | ✅ | `Mesh2D_t`, `SEMQuad`, mapped data classes |
| **2:1 mortar interfaces** (restriction/L2-projection operators, `mortarInfo` connectivity, `MortarExchange`/`MortarFluxCollect`, CPU + GPU + MPI) | ✅ | `Lagrange_t%mortarR/mortarP`, `src/SELF_Mesh_2D_t.f90`, mapped 2-D data classes, `src/gpu/SELF_Mortar.cpp` |
| Nonconforming test mesh constructor | ✅ | `SimpleMortarMesh` |
| Mortar validation tests (operators, exchange, conservation, dynamics; serial + MPI) | ✅ | `test/mortar*`, `test/*mortar*` |
| Aggregated GPU halo exchange with per-field packed buffers and persistent requests | ✅ | `SideExchangeStart/Finish`, `decomp%halo_*` tables |
| Legendre polynomial evaluation | ✅ | `LegendrePolynomial`, `src/SELF_Quadrature.f90` |
| Nodal→modal (Legendre) transform, spectral indicator | ❌ none | — |
| Quadtree/forest structures, refine/coarsen operators | ❌ none | — |
| Dynamic `nElem` (reallocatable fields, resizable GPU mirrors) | ❌ none | all data objects sized at `Init` |

Load-bearing consequences:

- **The hard interface problem is solved.** With mortars in place, an adapted
  quadtree mesh is directly representable: coarse/fine element interfaces
  become `mortarInfo` entries, and the entire solver stack (side exchange,
  mortar exchange, flux collect, DG operators) already runs on such meshes on
  both backends with MPI. AMR reduces to *mesh management + indicator +
  solution transfer*.
- **Nothing is sized dynamically.** `Mesh2D`, `SEMQuad`, and all field objects
  (host arrays and `c_ptr` device mirrors) are allocated once for a fixed
  `nElem`. Adaptation implies free/re-init of mesh, geometry, and model
  storage — acceptable when adaptation is infrequent relative to time steps.
- **Per-field exchange caches must be invalidated on adaptation.** The
  aggregated halo exchange lazily builds packed device buffers, persistent MPI
  requests (`halo_reqs`), and shared side tables (`decomp%halo_*`); the mortar
  exchange lazily allocates `mortarBuff`. All of these are sized by the
  current mesh/partition and must be torn down and rebuilt when the mesh
  changes.
- **BC bookkeeping is element-indexed.** `MapBoundaryConditions` caches
  `(element, side)` lists per BC with device copies; these must be rebuilt
  after every adaptation.

---

## 3. Mesh design: forest of quadtrees over a base mesh

The authoritative adaptive structure is a **forest of binary-refinement
quadtrees**, one tree per element of the user-supplied base mesh ("root
elements"). A refined element splits into 4 congruent children (2×2 in
reference space); each level halves h.

Each adaptation cycle produces an ordinary `Mesh2D` (the "leaf mesh") plus its
mortar table:

1. Evaluate the spectral indicator per element → refine/coarsen flags.
2. Smooth flags: enforce **2:1 balance** (edge-neighbor tree levels differ by
   at most 1); coarsen only complete 4-sibling groups whose removal keeps the
   balance valid.
3. Update the trees (split/merge nodes).
4. Emit leaves: assign global corner-node IDs (hash on quantized parametric
   position within the root) and a deterministic global ordering (roots in
   base-mesh order, Z-order/depth-first within each tree — locality-preserving
   for the block partition).
5. Build connectivity: equal-level edge pairings become conforming `sideInfo`
   entries (reusing the corner-pair hash matching already in
   `src/SELF_Mesh_2D_t.f90`, lifted into a shared helper, plus
   `RecalculateFlip`); level-differing edge pairings become `mortarInfo`
   entries, with sub-edge order and flips fixed by the tree geometry.
6. Regenerate `SEMQuad` geometry via the existing `GenerateFromMesh` path.

Design points:

- **Regeneration over incremental surgery.** Rebuilding the leaf mesh wholesale
  reuses the existing, MPI-correct construction path and is far easier to
  verify; adaptation cost is amortized over many steps. Incremental updates
  are a later optimization if profiling demands it.
- **No transition elements, no promotion passes.** Unlike the retired
  conforming design, tree leaves *are* the mesh; closure is the 2:1 balance
  constraint alone. Coarsening is a local 4-sibling merge.
- **Geometry of children.** Every leaf carries its affine map (corner points)
  in the root element's reference square plus `(rootElem, level)`. For
  `nGeo = 1` base meshes, child corners are the bilinear image of the dyadic
  sub-squares — exact, and the mortar interface's watertightness requirement
  is satisfied by construction. For curved elements (`nGeo > 1`), child
  geometry nodes are sampled from the parent's isoparametric map, preserving
  boundary geometry to the parent's order (both sides of a mortar then derive
  their metrics from the same parent map, which keeps the discrete
  free-stream property).
- **Boundary conditions.** A leaf edge on a root-element boundary edge
  inherits that edge's `sideInfo(5)` BC id; interior leaf edges get 0.
  `BCType`/`BCNames` copy from the base mesh unchanged.
- **Containment, not inheritance:**

```fortran
type :: AMRMesh2D_t
  type(Mesh2D) :: base       ! user-supplied root mesh (never modified)
  type(Mesh2D) :: leaf       ! regenerated mesh (+ mortarInfo); models point here
  ! forest state
  integer :: nRoots, nLeaves, maxLevel
  integer, allocatable :: leafRoot(:)      ! root element of each leaf
  integer, allocatable :: leafLevel(:)     ! tree level of each leaf
  real(prec), allocatable :: leafMap(:,:)  ! (2,nLeaves) dyadic offset in root ref square
  type(...), allocatable :: trees(:)       ! per-root quadtree (packed arrays)
contains
  procedure :: Init      ! from an existing Mesh2D (structured or file-read)
  procedure :: Adapt     ! flags in -> new leaf mesh + transfer maps out
  procedure :: Free
endtype
```

Models only need a `type(Mesh2D), pointer` at `amrMesh%leaf`, keeping
`DGModel2D`, the exchanges, geometry, and I/O unaware of AMR.

---

## 4. Spectral refinement indicator (Legendre modal analysis)

### 4.1 Nodal → modal transform

New machinery (currently absent — only `LegendrePolynomial` exists as a
building block). For control points `ξ_i` with quadrature weights `w_i` (Gauss
or Gauss–Lobatto, degree N), define the discrete Legendre coefficients of a
nodal field `u_i = u(ξ_i)`:

```
â_k = (1/γ_k) Σ_i w_i u_i L_k(ξ_i),   γ_k = Σ_i w_i L_k(ξ_i)²
```

Using the *discrete* norms `γ_k` makes the transform exactly invertible on the
nodal space for both node families (for Gauss, `γ_k = 2/(2k+1)` exactly for
k ≤ N; for GLL, `γ_N = 2/N`). In matrix form `â = P u` with `P = Γ⁻¹ Lᵀ W`;
the 2D transform is the tensor application `Â = P U Pᵀ` per element and
variable — the same two-pass batched shape as the existing derivative and
interpolation operators, so the GPU implementation follows the established
kernel pattern.

Implementation home: add `pMatrix(0:N,0:N)` (nodal→modal) and its inverse to
`Lagrange_t`, built in `Init_Lagrange_t` next to the mortar matrices, with
device mirrors following the `mortarR_gpu`/`mortarP_gpu` pattern. This also
unlocks future modal filtering/de-aliasing work.

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
          (and level > 0, and all 3 siblings also vote coarsen)
```

Defaults to be calibrated against manufactured cases (smooth Gaussian: no
flags; near-discontinuity: flags). All thresholds, the adapt cadence, the max
refinement level, and the indicator-variable set are user-facing parameters;
the variable set defaults to the `nstepped` prognostic variables (max over
variables), overridable per model (e.g. pressure only for `LinearEuler2D`).

The indicator is computed in reference space, so element size and mapping do
not bias it.

### 4.3 Backends

- **CPU:** `do concurrent` over elements; two small matmuls + O(N²) reduction
  per element. Negligible cost relative to one RHS evaluation.
- **GPU:** modal transform with `pMatrix_gpu` kernels; one custom kernel
  (`SpectralIndicator_2D_gpu`, following the `SELF_Mortar.cpp` file pattern)
  computes `q`, the decay fit, tail ratio, and writes per-element `int` flags;
  a single small D2H copy brings flags to the host, where all adaptation
  logic runs. The indicator runs only at adaptation checkpoints — no new
  host↔device traffic inside time stepping.

---

## 5. Solution transfer between meshes

Dyadic refinement makes the transfer operators tensor products of **the 1-D
mortar matrices that already exist**:

- **Refinement (parent → 4 children):** each child's nodal values are the
  parent polynomial evaluated on a half-interval in each direction — exactly
  the sub-edge restriction `mortarR(:,:,k)` applied per direction:
  `U_child(kx,ky) = R_kx U_parent R_kyᵀ`. Exact for polynomial data, hence
  conservative.
- **Coarsening (4 children → parent):** the volume L2 projection is the tensor
  product of the 1-D trace projections:
  `U_parent = Σ_{kx,ky} P_kx U_child(kx,ky) P_kyᵀ`, where `mortarP` already
  carries the 1/2 sub-interval Jacobian per direction (1/4 per child in 2-D,
  as required). The same exactness argument as the mortar operators gives
  `Σ P R = I` (refine-then-coarsen is the identity) and preservation of the
  element mean to roundoff on affine children.
- **Unchanged leaf:** direct copy (element index remap only).

Because refinement level changes by one per adaptation cycle (2:1 balance is
enforced before splitting), transfer never needs to traverse more than one
tree level per cycle.

A `MeshTransfer2D` object encapsulates old→new leaf correspondence; it is
built on the host during adaptation and applied as batched tensor-product
operations (host first; device application is a Phase-5 optimization, and the
operator matrices are already resident on the GPU).

**Numerical-correctness note (CLAUDE.md §3):** transfer is an
approximation-theoretic operation (interpolation/L2 projection), not a change
to the discretization. Entropy may increase slightly at coarsening (L2
projection is not entropy-stable); this is standard for AMR and will be
measured in validation. The same caveat already applies to the mortar flux
projection, and the same v1 policy is adopted: standard `DGModel2D` models
only — the EC split-form models reject nonconforming meshes at init and are
likewise out of scope for AMR v1.

---

## 6. Model integration

Adaptation is model-agnostic given `class(DGModel2D_t)`; physics lives in
overridden `flux2d`/`riemannflux2d`/BCs and is untouched. An
`AMRController2D` holds the `AMRMesh2D`, thresholds, cadence, and performs:

```
AdaptModel(model):
  1. indicator (device or host)             -> flags
  2. amrMesh%Adapt(flags)                   -> new leaf Mesh2D (+ mortarInfo) + MeshTransfer2D
  3. regenerate SEMQuad geometry            (existing GenerateFromMesh)
  4. capture old solution; Free + re-Init model fields; re-AssociateGeometry;
     rebuild BC lists (MapBoundaryConditions + device copies);
     drop per-field exchange caches (packed halo buffers, persistent halo
     requests, mortarBuff) so they rebuild lazily against the new mesh
  5. apply MeshTransfer2D to fill solution  (+ UpdateDevice on GPU builds)
  6. refresh decomposition (MPI: repartition + redistribute; §7)
```

New model types are thin wrappers that own a controller and trigger it from
the existing hook system (`PostStepHook`, every `nStepsPerAdapt` steps —
between complete steps, never inside an RK stage):

```fortran
type, extends(LinearEuler2D) :: AMRLinearEuler2D
  type(AMRController2D) :: amr
contains
  procedure :: PostStepHook => AdaptHook_AMRLinearEuler2D
endtype
```

Pilot models: `amr_advection_diffusion_2d` and `AMRLinearEuler2D` (CPU + GPU
variants following the existing backend file pattern). Metadata, equation
parsers, and integrator state carry over; `workSol` is rebuilt (zeroed) at
adaptation, which is correct because adaptation happens between complete
steps.

**Time step:** explicit CFL shrinks with the smallest h (factor 2 per level).
v1 keeps the user-controlled global `dt` contract: the controller exposes
`hMin()` and the examples rescale `dt` on adaptation. Local time stepping is
out of scope.

**I/O:** `Write_DGModel2D_t` already handles arbitrary meshes, but element
counts change between files. Add the leaf metadata
(`leafRoot`/`leafLevel`/`leafMap`) to the output so post-processing (pyself)
can reconstruct geometry per snapshot. Static-mesh outputs are untouched.

---

## 7. MPI strategy

- **Ordering & partition:** global leaf order = base-mesh root order, Z-order
  within each tree (deterministic, locality-preserving). Reuse the contiguous
  block partition (`GenerateDecomposition`) on the new `nLeaves`; equal
  element count is a valid balance measure because all elements share one
  polynomial degree.
- **Consistent flags:** refine/coarsen flags for ghost-adjacent elements must
  agree across ranks before 2:1 balancing. v1: `MPI_Allgatherv` of flags (one
  int per element) so every rank runs identical balancing — deterministic and
  simple; tree metadata is tiny compared to field data. Rank-local balancing
  with ghost-layer exchange is a scalability follow-up.
- **Redistribution:** old and new partitions are both contiguous over the same
  deterministic ordering, so each rank computes exact send/recv ranges
  (`MPI_Alltoallv` on packed solution payloads). Coarsening groups are applied
  on the receiving side after redistribution so all four siblings are
  resident.
- **Exchange table rebuild:** the aggregated halo tables (`decomp%halo_*`) and
  every field's persistent halo requests are partition-specific and are
  rebuilt lazily after step 4/6 above. The mortar sub-edge global side ids in
  the regenerated `mortarInfo` keep the MPI tag convention valid without any
  new machinery.
- **Constraint compliance:** all collectives happen inside the adaptation
  checkpoint (between time steps), never inside `CalculateTendency` or RK
  stages, per CLAUDE.md §5.

---

## 8. GPU strategy

Division of labor: **numerics on device, adaptation logic on host** — the same
split the mortar implementation uses.

Runs on device:
- Spectral indicator (`pMatrix_gpu` transforms + `SpectralIndicator_2D_gpu`
  reduction → per-element flags, one small D2H copy).
- All existing solver machinery on the adapted mesh (RHS, conforming and
  mortar exchanges, RK updates) — unchanged; this is already exercised by the
  mortar GPU tests.
- Solution transfer apply (Phase 5): batched tensor-product kernels using the
  device-resident `mortarR_gpu`/`mortarP_gpu`.

Runs on host:
- Tree updates, flag balancing, leaf emission, corner hashing,
  `RecalculateFlip`, mortar-table emission, partition computation.

**Memory churn:** each adaptation frees and re-allocates every field's device
mirrors plus the lazily-built exchange buffers (existing `Init`/`Free`/
`UpdateDevice` machinery). v1 accepts this (adaptation is infrequent). Phase-5
optimization: capacity-based allocation (`nElemCapacity = growthFactor*nElem`)
inside the data-object `Init`s so most adaptations skip reallocation.

One environment constraint to respect (and a reason to land issue #151
first if convenient): device selection currently happens in the mesh's
domain-decomposition init, so the AMR controller must never allocate device
resources before the base mesh exists — the same ordering rule the mortar
tests document.

---

## 9. Testing & validation plan

Unit (serial, both precisions, precision-aware tolerances):
- `modaltransform_roundtrip`: `P⁻¹ P u = u` to roundoff; Gauss and GL nodes.
- `spectralindicator_smooth` / `spectralindicator_front`: manufactured fields
  (low-degree polynomial ⇒ no refine flags; `tanh` front ⇒ flags exactly in
  front-crossing elements).
- `amrmesh2d_adapt`: randomized refine/coarsen cycles preserve mesh validity —
  every interior side is either a matched conforming pair or a well-formed
  `mortarInfo` entry; 2:1 balance holds; corner hashing watertight;
  `RecalculateFlip` idempotent.
- `meshtransfer_exactness`: refine transfer exact for polynomial data of
  degree ≤ N; refine-then-coarsen is the identity; coarsen preserves
  `∫ u dA` to roundoff (affine geometry).

Integration:
- `amr_advection_2d_gaussian`: advected Gaussian tracked by the indicator;
  error vs. a uniform-fine reference within tolerance; conservation error at
  roundoff.
- `amr_lineareuler2d_wave`: acoustic pulse; entropy non-increase across
  multiple adaptation cycles; refinement follows the wavefront.
- `amr_*_mpi` variants on 2 ranks, including at least one adaptation that
  changes the partition (elements migrate between ranks).
- GPU CI: the same integration tests on the Buildkite MI210/V100 pipelines;
  CPU/GPU agreement of indicator flags and post-adapt solutions to tolerance.

The mortar interface itself is already covered by the existing mortar test
suite and needs no re-validation here — AMR tests target the tree logic,
indicator, transfer, and re-initialization plumbing.

---

## 10. Work breakdown & phasing

| Phase | Content | New/changed code (est.) |
|---|---|---|
| 0. Modal infrastructure | `pMatrix` in `Lagrange_t` + device mirrors, indicator math + unit tests | ~600 LOC |
| 1. Quadtree mesh core | forest, 2:1 balancing, leaf + `sideInfo`/`mortarInfo` emission, `AMRMesh2D_t`, transfer operators (tensor products of mortar matrices), validity/transfer tests | ~1,500 LOC |
| 2. Model integration (CPU) | `AMRController2D`, field re-init + BC/cache rebuild path, `amr_advection_diffusion_2d`, `AMRLinearEuler2D`, integration tests, examples | ~1,200 LOC |
| 3. MPI | flag allgather, repartition + `Alltoallv` redistribution, halo/persistent-request invalidation, 2-rank tests | ~800 LOC |
| 4. GPU | indicator kernels + interfaces, device transfer apply, GPU CI tests | ~800 LOC |
| 5. Hardening/perf | capacity-based device allocation, adaptation-cost profiling, threshold calibration docs, user docs | ~500 LOC |

Rough total ~5–5.5k LOC plus tests — meaningfully smaller than the retired
conforming-template design (~6–7k), because the interface treatment, its GPU
kernels, its MPI paths, and its tests already exist. Phases 0–2 deliver a
usable serial CPU capability; 3 and 4 are independent once 2 lands.

---

## 11. Compliance with repository constraints (CLAUDE.md)

- No change to formulation, discretization order, basis, quadrature, nodal
  ordering, flux routines, or time integrators — AMR acts *between* steps; the
  interface discretization it relies on (mortars) is already merged and
  validated.
- No renamed/moved modules; all existing public interfaces preserved; new
  functionality is additive.
- No collectives inside time-stepping loops; halo and mortar exchange patterns
  unchanged — adaptation rebuilds their tables between steps.
- No external dependencies (forest is home-grown; no p4est).
- EC split-form models excluded (consistent with the existing mortar guard).
- Fortran 2008 only; formatting via `fprettify`.

---

## 12. Open questions for maintainers

1. **Adaptation cadence & placement:** `PostStepHook` every `nStepsPerAdapt`
   steps — acceptable, or should adaptation align with IO intervals in
   `ForwardStep`?
2. **`dt` policy on adaptation:** keep user-controlled global `dt` (examples
   rescale via `hMin()`), or auto-rescale inside the controller?
3. **Indicator variables:** default to all prognostic (`nstepped`) variables
   with max-combining — confirm, and confirm the per-model override hook.
4. **Curved base meshes (`nGeo > 1`) in v1**, or restrict v1 to `nGeo = 1`
   and defer curved-element child geometry?
5. **Coarsening policy:** strict 4-sibling groups only (proposed), or a more
   aggressive region-based merge?
6. **Output format for time-varying meshes:** embed leaf metadata per snapshot
   (proposed) vs. separate mesh-series files — coordinate with pyself.
7. **Sequencing with issue #151** (library-level MPI/device init): not a hard
   dependency, but landing it first removes the initialization-order
   constraint the AMR controller would otherwise have to respect.

---

## 13. References

- Y. Maday, C. Mavriplis, A. Patera, "Nonconforming mortar element methods:
  application to spectral discretizations," *Domain Decomposition Methods*,
  1989.
- D. A. Kopriva, "A conservative staggered-grid Chebyshev multidomain method
  for compressible flows. II. A semi-structured method," *JCP* 128 (1996)
  (the DGSEM 2:1 mortar algorithm underlying the implemented interfaces).
- P.-O. Persson, J. Peraire, "Sub-Cell Shock Capturing for Discontinuous
  Galerkin Methods," AIAA 2006-112 (modal tail-energy sensor).
- C. Mavriplis, "Adaptive mesh strategies for the spectral element method,"
  *CMAME* 116 (1994) (spectral decay-rate estimator).
- D. A. Kopriva, *Implementing Spectral Methods for Partial Differential
  Equations* (algorithms referenced throughout `SELF_Lagrange_t.f90` /
  `SELF_Quadrature.f90`).
