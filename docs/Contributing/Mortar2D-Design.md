# Design Assessment: Mortar (Nonconforming) Interfaces for 2D Models

**Status:** Implemented (v1) — see the user guide page
[Nonconforming (Mortar) Interfaces](../Models/nonconforming-mortar-interfaces.md).
This document is retained as the design rationale; the section on HOPR mortar-mesh
reading and quadtree AMR generation remains future work.
**Scope:** 2D quadrilateral meshes, 2:1 h-nonconforming interfaces ("h-mortars"),
CPU and GPU backends, MPI-parallel. p-mortars (mixed polynomial degree) are
noted as a future extension only.

Companion to the [conforming AMR design assessment](AMR2D-Design.md), which
excluded mortars by requirement. This document answers: *what would it take to
implement mortar elements instead*, so the two routes can be compared on equal
footing. The one-line summary of the comparison: **mortars move the complexity
out of mesh generation and into the core numerics and communication layer.**
The conforming-AMR mesh machinery gets dramatically simpler (plain 2:1
quadtree, factor-4 refinement, no transition templates), but the flux path,
the side-exchange machinery, the gradient (BR1) path, the MPI pattern, and the
GPU kernels all acquire a second interface type.

---

## 1. What a 2D h-mortar is

At a 2:1 nonconforming interface, one "big" element edge coincides with two
"small" element edges. The mortar method (Maday–Mavriplis–Patera 1989; DGSEM
form by Kopriva 1996) computes the interface flux in the *small* trace spaces
and transfers data with 1D L2 projections:

1. **Big → small (solution):** restrict the big edge trace, a polynomial of
   degree N on `[-1,1]`, to the sub-edges `[-1,0]` and `[0,1]`, sampled at each
   small side's quadrature points. Restriction of a polynomial is exact, so
   this is plain interpolation with two precomputed `(N+1)×(N+1)` matrices.
2. **Riemann solve on each small sub-edge** using the small side's own
   `nHat`/`nScale` — the standard two-point flux, unchanged.
3. **Small → big (flux):** the big side's surface integrand is the L2
   projection of the two small-side integrands back onto the big trace space
   (with a sign flip for the opposing normal). Projecting the *integrand*
   (flux × `nScale`) rather than the pointwise flux is what makes the scheme
   discretely conservative: the big element's surface integral equals the sum
   of the small elements' integrals to roundoff.

Only step 3 is mathematically new to SELF; steps 1–2 reuse existing machinery.

---

## 2. Why SELF's flux path is unusually mortar-friendly

Two existing design facts do most of the work:

- **Redundant per-element Riemann evaluation.** `BoundaryFlux_DGModel2D_t`
  (`src/SELF_DGModel2D_t.f90:403-431`) is a flat `do concurrent` over
  `(point, side, element)`: every element computes its own flux from its
  `boundary` (interior trace) and `extBoundary` (neighbor trace). There is no
  shared face loop. Therefore **small sides need no code changes at all** —
  they just need `extBoundary` filled with the restricted big trace, and the
  existing loop produces the correct mortar fluxes.
- **The stored quantity is already the projectable one.** Line 427 stores
  `boundaryNormal = riemannflux2d(sL,sR,dsdx,nhat) * nScale` — the surface
  integrand. The only new step in the tendency pipeline is a small post-pass
  that overwrites the *big* side's `boundaryNormal` with the projected sum of
  the two small sides' `boundaryNormal` values (sign-flipped). Everything
  downstream (`MappedDGDivergence`,
  `src/SELF_MappedVector_2D_t.f90:338-389`) is untouched.

The modified `CalculateTendency` order (cf. `src/SELF_DGModel2D_t.f90:570-603`):

```
1. solution%BoundaryInterp                      (unchanged)
2. solution%SideExchange                        (unchanged for conforming sides)
2m. MortarSolutionExchange   <- NEW: fill small extBoundary from big traces,
                                     big extBoundary from projected small traces
3-6. hooks / BCs / gradient / source            (see §4 for gradient)
7. BoundaryFlux                                 (unchanged - computes all sides)
7m. MortarFluxCollect        <- NEW: big boundaryNormal := -P(sum small integrands)
8-10. FluxMethod / MappedDGDivergence / dSdt    (unchanged)
```

Filling the big side's `extBoundary` with the projected small traces (step 2m)
is not required for the flux (7m overwrites the big side's integrand anyway)
but keeps `AverageSides`, BC hooks, and diagnostics well-defined on mortar
sides.

---

## 3. Component-by-component work list

### 3.1 Connectivity: `mortarInfo` alongside `sideInfo`

`sideInfo(1,:,:)` — the HOPR "Side Type" field — is currently unused in SELF
(noted at `src/SELF_Mesh_2D_t.f90:77-88`) and is the natural flag: `0` =
conforming, `>0` = index into a new mortar table. Proposed addition to
`Mesh2D_t`:

```fortran
integer :: nMortars
integer, pointer :: mortarInfo(:,:)   ! (1:8, 1:nMortars)
  ! 1: big element (global), 2: big local side
  ! 3: small element 1,      4: 10*side + flip   (sub-edge on big's xi in [-1,0])
  ! 5: small element 2,      6: 10*side + flip   (sub-edge on big's xi in [ 0,1])
  ! 7: global side id (MPI tag base), 8: reserved (future p-mortar / 3D)
```

Sub-edge ordering is defined in the **big side's edge coordinate**, so the
existing binary `flip` convention still suffices in 2D. `RecalculateFlip`
(`src/SELF_Mesh_2D_t.f90:1089-1235`) gains a mortar branch that matches each
small edge's corner pair against the big edge's endpoint + midpoint node.

Sources of nonconforming meshes:
- **Quadtree AMR** (the main driver): with mortars allowed, the AMR mesh layer
  in the companion document simplifies to a plain 1→2×2 quadtree with 2:1
  balance — no promotion pass, no transition templates, no ephemeral elements,
  clean 4-sibling coarsening. Mortar interfaces fall directly out of the tree
  (big side = coarser leaf edge, smalls = the two finer leaf edges).
- **Static nonconforming meshes** (optional input path): HOPR encodes mortars
  in its `SideInfo` (MortarType / negative neighbor indices). The 2D reader
  (`Read_HOPr_Mesh2D_t`, `src/SELF_Mesh_2D_t.f90:480-648`) could translate
  these; useful for testing mortars independently of AMR, but not required.

### 3.2 Projection operators (on `Lagrange_t`)

Four precomputed `(N+1)×(N+1)` matrices, built in `Init_Lagrange_t` next to
`iMatrix`/`dgMatrix`/`bMatrix` and mirrored to device like the others
(`src/gpu/SELF_Lagrange.f90:150-185` pattern):

- `mortarR1, mortarR2` — restriction of the big trace to sub-edges
  `[-1,0]`, `[0,1]` (Lagrange interpolation at scaled points; exact).
- `mortarP1, mortarP2` — L2 projections from each sub-edge space back to the
  big space: `P_k = M⁻¹ R_kᵀ M_k · (1/2)`, with `M = diag(qWeights)` and the
  1/2 the sub-edge-to-big Jacobian. Discretely, `P_k(i,j) =
  (1/2) w_j R_k(j,i) / w_i`. The identity `P₁R₁ + P₂R₂ = I` (forward-backward
  consistency) is a build-time assertion and a unit test.

With Gauss points the discrete L2 projection is exact; with Gauss–Lobatto it
is the standard lumped-GLL projection (same choice as `dgMatrix`), consistent
with the rest of the discretization.

### 3.3 Data-plane changes (`MappedScalar2D` / `MappedVector2D`)

- `MortarSolutionExchange`: loop over mortars; small `extBoundary(:,s,e,:) =
  R_k · bigBoundary` (flip-aware); big `extBoundary = P₁·small₁ + P₂·small₂`.
  CPU: `do concurrent` over `(point, mortar, var)`; GPU: one kernel in the
  mold of `SideExchange_2D_gpu` (`src/gpu/SELF_MappedData.cpp`).
- `MortarFluxCollect` (operates on `flux%boundaryNormal`): big-side integrand
  := `-(P₁·f₁ + P₂·f₂)` with each small integrand reversed per its flip. The
  sign flip accounts for the opposing outward normals.
- `AverageSides` (`src/SELF_Scalar_2D_t.f90:159-174`) works unchanged *given*
  step 2m filled mortar `extBoundary` consistently.

No changes to array shapes: `boundary`/`extBoundary` stay
`(N+1, 4, nElem, nVar)`. Mortars need no extra trace storage in 2D because the
big side's projected data lives in its ordinary `extBoundary` slot; a small
scratch buffer `(N+1, 2, nMortars, nVar)` is only needed transiently inside the
two mortar routines (allocated once at init, device-mirrored).

### 3.4 Gradient / BR1 path (parabolic terms)

`CalculateSolutionGradient` (`src/SELF_DGModel2D_t.f90:331`) consumes
`avgBoundary`; with 2m in place, small sides average against restricted big
traces and the big side against projected small traces. This is the standard
"project-then-average" nonconforming BR1 and is the v1 choice. Caveat to state
in the docs: naive BR1 across nonconforming interfaces can lose the provable
stability of the conforming scheme; the corrected schemes (e.g. Friedrich et
al. 2018, Kopriva-style lifting on mortar spaces) can be added later behind
the same two hook points. Advection-dominated and linear models are
unaffected in practice.

### 3.5 Entropy-conservative models — explicitly deferred

`ECDGModel2D` / split-form models (`dSplitMatrix`, two-point volume fluxes)
require *entropy-stable* mortar operators (Friedrich, Winters, Del Rey
Fernández, Gassner et al., JSC 2018) — projections alone do not preserve the
entropy estimate. v1 restricts mortars to the standard `DGModel2D` family and
raises a runtime error for EC models on nonconforming meshes.

### 3.6 Geometry and watertightness

Small sides use their own `nHat`/`nScale` from `SEMQuad` — already computed
per element. Two consistency requirements:

- **Watertight subdivision:** small edges must lie exactly on the big edge.
  For quadtree children this holds by construction at `nGeo = 1`; for curved
  parents (`nGeo > 1`), child geometry must be sampled from the parent's
  isoparametric map (same rule as the conforming design) so the traces
  coincide and free-stream preservation survives. A discrete free-stream test
  guards this.
- **Metric identity on the interface:** with watertight geometry,
  `nScale_big(ξ) = 2·nScale_small(ξ_k)` pointwise for affine faces, and the
  projection identity guarantees `∮ big = Σ ∮ small` for the integrands. For
  curved faces the identity holds at the discrete level because both sides'
  metrics derive from the same parent map.

### 3.7 MPI

Mortar interfaces can straddle ranks (big on one rank, smalls on one or two
others). Design: **redundant computation, single exchange round** — the same
philosophy as the current conforming exchange, where both elements compute the
shared flux independently:

- `MPIExchangeAsync` (`src/SELF_MappedScalar_2D_t.f90:121-172`) gains a mortar
  branch: the big rank sends its big trace to each remote small rank; each
  small rank sends its small trace to the big rank (and nothing else — smalls
  never need each other). Tags reuse the global-side-id scheme with the
  sub-edge index folded in: `tag = globalSideId*2 + k + 2*nUniqueSides*(ivar-1)`.
- Each rank then runs `MortarSolutionExchange` / `MortarFluxCollect` locally
  on the traces it holds; the big-side integrand is computed identically
  (same operations, same order) wherever it is needed, preserving determinism
  and requiring **no second communication round** and no flux send-back inside
  the RK stage.
- One structural change to size bookkeeping: `maxMsg` accounting (currently
  `nUniqueSides`-based) must count mortar sub-messages.

This keeps the "halo exchange pattern" qualitatively intact (async
point-to-point, posted per side, one round per stage) but it *is* a
modification of a fixed pattern per CLAUDE.md §5/§11 — maintainer approval
required.

### 3.8 GPU

New device pieces, all following existing patterns:
- `mortarInfo_gpu` on the gpu `Mesh2D` (mirrors `sideInfo_gpu`,
  `src/gpu/SELF_Mesh_2D.f90`).
- `mortarR/P` matrices on the gpu `Lagrange`.
- Kernels in a new `src/gpu/SELF_Mortar.cpp` + interfaces in
  `SELF_GPUInterfaces.f90`: `MortarSolutionExchange_2D_gpu`,
  `MortarFluxCollect_2D_gpu` (one thread per point per mortar per var; the
  `(N+1)×(N+1)` matvecs are tiny — hand-rolled loops as in
  `BoundaryInterp_2D_gpu`, no hipBLAS needed at these sizes).
- The existing `SideExchange_2D_gpu` / `ApplyFlip_2D_gpu` kernels skip sides
  flagged as mortar (branch on `sideInfo(1)`).

No new host↔device traffic in the time loop; mortar traces stay resident.

### 3.9 Model layer & API surface

No changes to any concrete model (`flux2d`, `riemannflux2d`, BCs untouched).
`CalculateTendency_DGModel2D_t` (and its gpu override) gain the two mortar
calls, gated on `mesh%nMortars > 0` so conforming meshes execute byte-identical
code paths — existing regression tests and bitwise reproducibility on
conforming meshes are unaffected. BC mapping (`MapBoundaryConditions`,
`src/SELF_DGModel2D_t.f90:453-530`) is untouched: mortar sides are interior
sides (`sideInfo(5) = 0`).

---

## 4. Interaction with the AMR design

If mortars are accepted, the companion AMR design changes as follows:

| Aspect | Conforming (3-refinement) | With mortars (2:1 quadtree) |
|---|---|---|
| Refinement granularity | 9 children / level | 4 children / level |
| Mesh generation | Templates + promotion + closure passes | Trivial (tree leaves are the mesh) |
| Transition elements | Yes (quality-reduced, ephemeral) | None |
| Coarsening | 9-sibling groups, closure rebuild | 4-sibling groups, local |
| Solution transfer | Through parent space (template maps) | Standard dyadic restriction/prolongation |
| Core numerics | Untouched | Flux path + exchange + BR1 modified |
| MPI exchange | Untouched | Mortar branch in async exchange |
| GPU kernels | Indicator + transfer only | + 2 mortar kernels, exchange branch |
| EC/split-form models | Work as-is | Deferred (entropy-stable mortars) |
| Formulation risk | None (pure meshing) | Interface discretization added |

The spectral indicator (Legendre modal transform, decay/tail criteria) and the
adaptation controller are **identical in both designs** and should be built
first either way (Phase 0 of the AMR document).

---

## 5. Testing & validation

- `mortarops_identity`: `P₁R₁ + P₂R₂ = I` to roundoff; restriction exactness
  for degree-N polynomials (Gauss and GL).
- `mortar_conservation`: two-element big/small patch; surface integrals of
  `boundaryNormal` on big side equal minus the sum over small sides to
  roundoff, for a nonlinear flux.
- `mortar_freestream`: constant state on a nonconforming (and curved, when
  `nGeo>1` lands) mesh stays constant to roundoff over many steps.
- `mortar_convergence`: manufactured solution on a fixed 2:1 mesh; design
  order preserved (mortars are known not to degrade order at 2:1).
- `mortar_*_mpi`: same patches split so mortars straddle 2 ranks; results
  bitwise-match the serial run.
- GPU parity tests on the Buildkite pipelines, matching CPU to tolerance.
- All existing conforming-mesh regression tests must remain bitwise unchanged
  (the `nMortars = 0` gate makes this checkable).

---

## 6. Effort estimate

| Piece | Est. |
|---|---|
| Projection operators on `Lagrange_t` (+ gpu/apu mirrors, unit tests) | ~350 LOC |
| `Mesh2D_t` mortar connectivity (+ flip logic, quadtree emission, optional HOPR translation) | ~800 LOC |
| CPU mortar exchange + flux collect (scalar + vector, gradient path) | ~700 LOC |
| `CalculateTendency` integration + gating + EC guard | ~200 LOC |
| MPI mortar branch + message accounting | ~500 LOC |
| GPU kernels + interfaces + device mirrors | ~700 LOC |
| Tests (conservation, free-stream, convergence, MPI, GPU) | ~700 LOC |
| **Total (mortar capability)** | **~4k LOC** |

For the full AMR-with-mortars picture: this ~4k replaces roughly 2–2.5k LOC of
template/closure/transfer complexity in the conforming design's Phase 1, so
**total effort for AMR is similar on both routes (~6–7k LOC)**. The difference
is *where the risk sits*: conforming keeps all risk in mesh generation
(verifiable by connectivity audits, zero numerics risk); mortars put risk in
the interface discretization, the exchange layer, and stability theory (BR1,
EC models), but buy factor-4 refinement, simpler meshing, and the
industry-standard AMR structure (FLEXI, Trixi.jl, Nektar++ all use mortars).

---

## 7. Constraint compliance & required approvals (CLAUDE.md)

Mortars, unlike the conforming route, **do touch guarded areas** and need
explicit maintainer sign-off on three points:

1. **Mathematical formulation (§3):** a new interface discretization
   (projection-based mortar flux) is added. It provably preserves
   conservation and design order at 2:1 interfaces, but it is a formulation
   extension, not a pure meshing change.
2. **Halo exchange pattern (§5/§11):** the async point-to-point pattern is
   preserved (one round per stage, no new collectives), but mortar sub-face
   messages are added to it.
3. **Split-form/EC models (§3):** excluded from mortar meshes in v1 pending
   entropy-stable mortar operators.

Everything else stays within bounds: no module renames, no public-interface
changes, no external dependencies, conforming code paths bitwise-identical.

---

## 8. References

- Y. Maday, C. Mavriplis, A. Patera, "Nonconforming mortar element methods:
  application to spectral discretizations," *Domain Decomposition Methods*,
  1989.
- D. A. Kopriva, "A conservative staggered-grid Chebyshev multidomain method
  for compressible flows. II. A semi-structured method," *JCP* 128 (1996)
  (the DGSEM 2:1 mortar algorithm).
- L. Friedrich, A. R. Winters, D. C. Del Rey Fernández, G. J. Gassner,
  M. H. Carpenter, "An entropy stable h/p non-conforming discontinuous
  Galerkin method with the summation-by-parts property," *JSC* 77 (2018).
- F. Hindenlang, G. Gassner et al., FLEXI / FLUXO implementations of 2:1
  mortars in curvilinear DGSEM (reference production implementations).
