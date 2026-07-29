# Adaptive Mesh Refinement (AMR)

This page describes SELF's approach to adaptive mesh refinement for 2-D discontinuous
Galerkin spectral element models. It documents the **refinement trigger** that is
implemented today and lays out a **staged design** for the remaining mesh-adaptation
machinery so that the invasive parts (dynamic mesh mutation, solution transfer, MPI
re-partitioning, and GPU re-allocation) can be reviewed and landed incrementally.

AMR in SELF builds directly on the 2:1 nonconforming (mortar) interface support
(see [Nonconforming (Mortar) Interfaces](../Models/nonconforming-mortar-interfaces.md)).
Refining a quadrilateral element into four children produces exactly the 2:1 hanging-node
configurations that the mortar operators already couple correctly, so the flux-coupling
layer that adaptive refinement needs is already in place and tested.

---

## 1. Scope and status

| Component | Status |
| --- | --- |
| Refinement **trigger** (Legendre modal-decay indicator, CPU + GPU) | **Implemented** |
| `h`-refinement primitives (isoparametric subdivision + refined connectivity) | **Implemented** |
| Uniform `h`-refinement (conforming, serial) | **Implemented** |
| Adaptive quad-forest (flagged refine / coarsen, level tracking) | **Implemented** |
| Solution transfer (prolongation / restriction, conservative) | **Implemented** |
| Forest face-neighbour queries + 2:1 balancing | **Implemented** |
| Hanging-node / mortar-table + `Mesh2D_t` emission | **Implemented** |
| Adaptation-epoch transfer plan (old-leaf → new-leaf mapping) | **Implemented** |
| MPI dynamic re-partitioning / load balancing | Designed (Stage 5) |
| GPU device re-allocation for a changing element count | Designed (Stage 6) |

The full **serial** adaptive-refinement loop is wired into the library today: flag with the
Stage-1 indicator, mutate the forest (Stage 2b), transfer the solution (Stage 3), balance and
emit a runnable nonconforming `Mesh2D_t` (Stage 4). What remains is **scaling** that loop -
dynamic MPI re-partitioning across ranks (Stage 5) and GPU device re-allocation for a changing
element count (Stage 6) - each deferred behind this design so it lands as a self-contained,
reviewable piece.

---

## 2. The refinement trigger (implemented)

The trigger answers a single, element-local question: *is the solution on this element
well-resolved by the current polynomial degree, or does it need more resolution?* It never
mutates the mesh; it only produces a per-element flag. This makes it cheap, embarrassingly
parallel, free of MPI communication, and identical in structure on CPU and GPU.

### 2.1 Legendre modal-decay (spectral) indicator

For a nodal solution field \(u_{ij}\) on an element, we form its tensor-product Legendre
modal expansion

$$
u(\xi,\eta) = \sum_{p=0}^{N}\sum_{q=0}^{N} \hat u_{pq}\,
\tilde L_p(\xi)\,\tilde L_q(\eta),
$$

where \(\tilde L_p\) are the \(L^2\)-normalized Legendre polynomials on \([-1,1]\)
(\(\int_{-1}^1 \tilde L_p \tilde L_q\,d\xi = \delta_{pq}\)). The nodal-to-modal transform is the
exact inverse of the Legendre Vandermonde \(V_{ip} = \tilde L_p(\xi_i)\) built from the
interpolant control points, so the indicator is **independent of the control-node type**
(Legendre–Gauss or Legendre–Gauss–Lobatto) and is exact for polynomial data. With the
normalized basis the total modal energy equals the exact \(L^2\) energy on the reference
element,

$$
E_{\text{tot}} = \sum_{p,q}\hat u_{pq}^2 = \lVert u\rVert_{L^2([-1,1]^2)}^2 .
$$

The smoothness of the field is measured by how much of that energy sits in the highest modes.
Defining the clipped energies that drop the highest one and two modes in each direction,

$$
E_{\text{clip1}} = \sum_{p,q\le N-1}\hat u_{pq}^2, \qquad
E_{\text{clip2}} = \sum_{p,q\le N-2}\hat u_{pq}^2,
$$

the smoothness ratio is the larger of the top-shell and next-shell energy fractions,

$$
S_e = \max\!\left(\frac{E_{\text{tot}}-E_{\text{clip1}}}{E_{\text{tot}}},\;
                  \frac{E_{\text{clip1}}-E_{\text{clip2}}}{E_{\text{clip1}}}\right),
\qquad \sigma_e = \log_{10} S_e .
$$

A smooth, well-resolved field has energy concentrated in the low modes, so \(S_e\) is tiny and
\(\sigma_e\) is very negative; an under-resolved field (a steep front or a discontinuity)
leaves energy in the top modes, driving \(\sigma_e\) toward 0.

This is the modal-energy indicator of **Persson & Peraire (2006)**. The second (next-shell)
term is the robustification introduced by **Hennemann, Rueda-Ramírez, Hindenlang & Gassner
(2021)** — the shock indicator used in Trixi.jl — which guards against the odd/even parity
dropouts a single-mode measure suffers for symmetric data. It is only active for \(N\ge 2\).

### 2.2 Refine / keep / coarsen semantics

With user-supplied thresholds \(\sigma_{\text{refine}} > \sigma_{\text{coarsen}}\), each element
is flagged

| Condition | Flag |
| --- | --- |
| \(\sigma_e > \sigma_{\text{refine}}\) | `SELF_AMR_REFINE` (+1) |
| \(\sigma_e < \sigma_{\text{coarsen}}\) | `SELF_AMR_COARSEN` (−1) |
| otherwise | `SELF_AMR_KEEP` (0) |

Two thresholds (rather than one) give a hysteresis band that prevents elements from thrashing
between refine and coarsen on successive checks. Recommended starting values for double
precision are \(\sigma_{\text{refine}}\approx-3\) and \(\sigma_{\text{coarsen}}\approx-8\); the
best values are problem-dependent and are deliberately required arguments rather than hidden
defaults.

### 2.3 API

The indicator lives in `SELF_RefinementIndicator_2D` and follows the standard SELF backend
pattern: a portable base type with a `do concurrent` implementation
(`SELF_RefinementIndicator_2D_t`), a thin CPU extension, and a GPU extension that runs a
device kernel (`SELF_Refinement.cpp`) and copies the resulting flags back to the host, where
the mesh-adaptation logic runs.

```fortran
use SELF_RefinementIndicator_2D

type(RefinementIndicator2D) :: amr

! interp : the model interpolant (Lagrange); nElem : rank-local element count
call amr%Init(interp, nElem, refineThreshold=-3.0_prec, coarsenThreshold=-8.0_prec)

! solution : the model's MappedScalar2D (or any Scalar2D); ivar selects the driving
! variable, or SELF_AMR_ALLVARS (=0) reduces over all variables (most conservative).
call amr%Estimate(solution, ivar=1)

! amr%indicator(1:nElem)  -- per-element sigma_e (host)
! amr%flag(1:nElem)       -- per-element SELF_AMR_REFINE / KEEP / COARSEN (host)
n = amr%CountFlagged(SELF_AMR_REFINE)
```

The GPU device kernel uses per-thread scratch bounded to \(N \le 15\); higher degrees fail
loudly at `Init` with a pointer to the bound (`AMR2D_MAXNP` in `SELF_Refinement.cpp`).

### 2.4 Validation

`test/refinement_indicator_2d_spectraldecay.f90` checks the indicator against fields with known
spectra, on both Gauss and Gauss–Lobatto nodes:

- a **constant** field (energy only in mode 0) floors out → `COARSEN`;
- a **pure highest tensor mode** \(\tilde L_N(\xi)\tilde L_N(\eta)\) has \(S_e=1\) so
  \(\sigma_e = 0\) to roundoff → `REFINE` (this also verifies the transform recovers a single
  mode exactly);
- a **low-degree polynomial** (degree \(\le N-2\)) has no top-shell energy → `COARSEN`;
- a **steep tanh front** is under-resolved → `REFINE`.

---

## 2.5 `h`-refinement primitives and uniform refinement (implemented)

The mesh-mutation core of Stage 2 is implemented as two additive modules that leave every
existing mesh type and interface untouched:

- **`SELF_RefinementPrimitives_2D`** — the element-local, dependency-light pieces:
    - `SubdivideNodeCoords` performs **exact isoparametric subdivision** of one element's
      geometry into its four children. The parent geometry is the degree-`nGeo` Lagrange
      interpolant through the element's geometry nodes; each child node coordinate is that
      interpolant evaluated at the corresponding point of the parent reference square. This is
      exact for straight-sided and curved (isoparametric) elements alike and for any control
      node type, so refinement never perturbs the represented domain.
    - `RefineConnectivity` builds the refined mesh's connectivity by **pure deterministic
      integer bookkeeping**. Sibling faces interior to a parent are same-orientation (flip 0); a
      child face on a parent boundary inherits the parent face's neighbor, side, and flip, with
      the sub-position pairing across the face taken from the parent flip. Because refinement is
      orientation preserving, this reproduces exactly the flips a corner-node matching pass would
      compute, without requiring globally consistent node ids (which structured meshes do not
      guarantee) and without any coordinate hashing.

- **`SELF_MeshRefinement_2D`** — `UniformRefineMesh(meshIn, meshOut)` assembles a fully-formed
  `Mesh2D_t` with 4× the elements from those primitives: child geometry by isoparametric
  subdivision, connectivity and flips inherited from the base, boundary-condition and material
  metadata carried over, and a (serial) domain decomposition. The result is conforming (no
  hanging nodes), so it needs neither mortars nor 2:1 balancing and is immediately usable for
  geometry generation and time stepping. It is a genuine capability on its own (e.g. grid
  convergence studies) and it exercises all of the subdivision, connectivity, and array
  machinery that adaptive refinement will reuse.

Uniform refinement is currently **serial** (it replicates the input mesh's single-rank MPI
state rather than re-decomposing); multi-rank refinement is Stage 5.

**Validation.** The primitives are unit-tested in isolation (exact child geometry for affine and
curved elements; neighbor / side / flip / global-side-id reciprocity of the refined
connectivity). The end-to-end `UniformRefineMesh` is covered by
`test/mesh2d_uniform_refine.f90`: refining a structured mesh quadruples the element count,
conserves the domain area (integral of the Jacobian) to roundoff, keeps every Jacobian strictly
positive, produces reciprocal interior connectivity, and doubles the physical-boundary side
count.

## 2.6 Adaptive quad-forest (implemented)

`SELF_QuadTreeMesh_2D` provides the adaptive mesh-mutation data structure of Stage 2b. Each base
element is the root of a quadtree; `QuadTreeMesh2D` stores a growable node pool (level, parent,
per-node child pointers, originating root) plus the base geometry, and maintains the active leaf
set.

- `Init(mesh)` seeds one root per base element (all leaves at level 0).
- `AdaptFromFlags(flag)` consumes a per-leaf flag array **indexed exactly like the Stage-1
  indicator's `flag(:)`** (`+1` refine, `-1` coarsen, `0` keep): flagged leaves are subdivided
  into four children, and a family of four leaf siblings is merged back only when all four are
  flagged coarsen (the standard de-refinement rule). Refine and coarsen are resolved against the
  same pre-adaptation snapshot so they never interfere.
- `LeafCoords(i, geomInterp, coords)` regenerates any leaf's physical geometry by repeated exact
  isoparametric subdivision of its root along the quadtree path (reusing the §2.5 primitive), so
  the forest never stores redundant coordinates.
- Node storage grows by amortized doubling; the leaf set is always recovered by traversal from
  the roots, which makes nodes orphaned by coarsening invisible without an explicit free list.

This wires directly to the trigger: `indicator%Estimate(solution, ivar)` then
`forest%AdaptFromFlags(indicator%flag)` performs one adaptation step. What the forest deliberately
does **not** yet do is answer *face-neighbour* queries, enforce 2:1 balance, or emit a
solver-ready `Mesh2D_t` with hanging-node mortars - an adaptively refined forest is generally
nonconforming, and turning it into a runnable mesh is Stage 4.

The forest is unit-tested (refine/coarsen leaf and level bookkeeping, the four-sibling coarsening
guard, amortized capacity growth, and leaf geometry matching direct subdivision at multiple
levels).

## 2.7 Solution transfer (implemented)

`SELF_SolutionTransfer_2D` moves the prognostic solution with the mesh when an element is refined
or coarsened, reusing the mortar operators already built on `Lagrange`:

- `ProlongToChildren(interp, nVar, uParent, uChildren)` samples the parent's degree-N nodal
  polynomial onto its four children - a tensor product of `mortarR` in each direction. Exact
  interpolation, no loss.
- `RestrictFromChildren(interp, nVar, uChildren, uParent)` L2-projects the four children back onto
  the parent - a tensor product of `mortarP` (the adjoint of `mortarR`, already carrying the 1/2
  per-direction sub-edge Jacobian for solution traces).

Both are element-local and portable (host `do concurrent`); the driver maps forest parent/child
relations onto the element index ranges. From the 1-D mortar identities they inherit exactly:

- **Reversibility** - `RestrictFromChildren(ProlongToChildren(u)) = u` to roundoff (refine then
  immediately coarsen leaves the solution unchanged), from `sum_k P_k R_k = I` per direction.
- **Conservation** - the reference-cell integral of the restricted parent equals the sum of the
  children's, i.e. `sum_ij w_i w_j u_parent = (1/4) sum_c sum_ij w_i w_j u_child_c`; weighted by
  the geometry Jacobian this is conservation of the cell-integrated quantity.

Validated at two levels: a unit test drives the tensor operators with the exact mortar matrices
(prolong reproduces a degree-N polynomial at the child nodes; prolong-then-restrict is the
identity to roundoff; discrete conservation defect is zero), and `test/solution_transfer_2d.f90`
checks the end-to-end behaviour against Stage-2 geometry - prolonging a coarse field onto a
`UniformRefineMesh` and back is reversible, and `int u dA` (with the geometry Jacobian) is
identical on the coarse and refined meshes.

## 2.8 Forest face-neighbours and 2:1 balancing (implemented)

`SELF_QuadTreeMesh_2D` also answers *face-neighbour* queries on the forest and enforces the 2:1
balance condition - the first part of Stage 4 and the prerequisite for mortar generation. The
forest stores the base mesh's root face connectivity (`rootNbr` / `rootNbrSide` / `rootFlip`,
from `sideInfo`; a conforming base is assumed).

- `FaceNeighbor(node, s, nbr, ns, nf)` returns the equal-or-larger neighbour across local side s
  by the classic quadtree ascend/descend search: cross to a sibling when the face is interior to
  the parent, otherwise ascend to the parent's neighbour and descend one level, matching
  sub-positions across the face through the base flip. `nbr` is either a leaf at any level
  `<= level(node)` or an internal node at exactly `level(node)`; `ns` / `nf` are the neighbour's
  facing side and the edge flip. Consequently a 2:1 hanging face is precisely "`nbr` is a leaf
  with `level(nbr) = level(node)-1`" and finer neighbours are precisely "`nbr` is internal" -
  exactly the classification Stage 4b needs to emit `mortarInfo`.
- `Balance2to1()` iterates to a fixed point: any leaf whose equal-or-larger neighbour is a leaf
  two or more levels coarser refines that neighbour, and the (possibly rippling) refinement
  repeats until no face violates the condition.
- `MaxLevelJump()` reports the largest level difference across any leaf face (0 conforming, 1 for
  a balanced adaptive forest) - a cheap invariant for tests and drivers.

Unit-tested standalone: level-0 neighbour queries on a structured base (including boundaries and
directional reciprocity); a uniformly refined forest is conforming (`MaxLevelJump = 0`); an
adaptive refinement that creates a two-level jump is reduced to one level by `Balance2to1`
(rippling into the coarse neighbour); and equal-level leaf-neighbour reciprocity holds across the
balanced forest.

## 2.9 Mesh emission (implemented) — closing the loop

`SELF_AdaptiveMesh_2D`'s `EmitMesh(forest, baseMesh, outMesh)` turns a 2:1-balanced forest into a
solver-ready `Mesh2D_t`. Each leaf becomes an element (leaf-list order); for every leaf face the
Stage-4a `FaceNeighbor` classification drives the emitted connectivity:

- **domain boundary** → `sideInfo(3)=0`, `sideInfo(5)` = the base element's BC id on that side;
- **same-level leaf** → a conforming interior side (`sideInfo(3)=neighbour`, `(4)=10*side+flip`,
  a shared global side id);
- **one-level-finer neighbour** → this leaf is the **big** side of a 2:1 mortar; the two small
  elements are the finer neighbour node's children on the shared face, with the big edge
  coordinate `[-1,0]`/`[0,1]` mapped to the correct child through the face flip;
- **one-level-coarser neighbour** → a **small** side, filled when its big side is processed.

Mortar sides carry `sideInfo(1)=mortar index` and `sideInfo(3)=sideInfo(5)=0`, and the emitted
`mortarInfo(1:8, :)` follows the exact layout of the hand-built `SimpleMortarMesh` (big elem/side;
small elem + `10*side+flip` per sub-edge; two sub-edge global side ids), so the existing mortar
side-exchange, projection, and flux machinery consume it unchanged. Leaf geometry comes from
`LeafCoords`; the output uses a serial decomposition (Stage 5 will re-partition).

Validated at two levels. A standalone structural test confirms - on an adaptively refined,
balanced forest - the side classification is exclusive, conforming sides are reciprocal, every
mortar's big and small sides reference the same mortar index, and (geometrically) each big face is
exactly tiled by its two small faces. The full-pipeline CI test `test/adaptive_mortar_2d.f90`
refines a structured mesh, balances, emits the mesh, builds its geometry (strictly positive
Jacobians; total area equal to the base mesh), and runs the real `SideExchange` + `MortarExchange`
on a globally linear field: the external trace matches the interior trace to roundoff on every
conforming **and** 2:1 mortar side - the same criterion the hand-built mortar-mesh tests use,
now on a mesh produced entirely by the AMR pipeline.

With this, the serial loop closes: `indicator → forest.AdaptFromFlags → (transfer solution) →
forest.Balance2to1 → EmitMesh` yields a runnable adaptive mesh.

## 2.10 Adaptation-epoch transfer plan (implemented)

`SELF_TransferPlan_2D` is the driver layer of Stage 3: it connects the element-local transfer
operators (§2.7) to an actual forest mutation. One *adaptation epoch* is: snapshot the leaf
list (`nOld`, `oldLeaf`), mutate the forest (at most one `AdaptFromFlags`, then any number of
refinements — `Balance2to1`, `RefineNode`), then `BuildTransferPlan(forest, nOld, oldLeaf,
plan)`. The plan records, for every new leaf in leaf-list order (the element ordering `EmitMesh`
produces), where its solution comes from in the old element ordering:

- **copy** — the leaf survived unchanged;
- **prolong** — the leaf descends from an old leaf; the old polynomial is interpolated down the
  quadtree path, one step per level, so a fresh child re-refined by balance ripple in the same
  epoch is handled by depth > 1;
- **restrict** — the leaf is (an ancestor of) a coarsened family; the four old children are
  L2-projected onto their parent, then prolonged down any further steps (depth > 0 occurs when
  a just-coarsened parent is immediately re-refined by balancing).

Reconstruction after the fact is possible because forest node ids are stable: refinement
appends nodes and coarsening only detaches children, whose `level`/`parent`/`quadrant` entries
persist. Each new leaf ascends its parent chain until it meets an old leaf or a complete
old-leaf family; a snapshot that cannot explain a leaf fails loudly. `ApplyTransferPlan`
executes the plan on nodal data in the `MappedScalar2D%interior` layout and inherits the §2.7
identities (exact prolongation, conservative restriction, exact refine-coarsen round trips).

`test/transfer_plan_2d.f90` validates one epoch that simultaneously coarsens (then re-refines)
a family, refines a leaf whose child is refined again (depth-2 prolongation), and lets 2:1
balancing ripple into an untouched root: the classification multiplicities are checked exactly;
a bilinear field is reproduced at the emitted new mesh's nodes to roundoff through all three
transfer kinds; the Jacobian-weighted global integral of a non-polynomial field is conserved to
roundoff; and refine-everything/coarsen-everything across two epochs is the identity.

---

## 3. Comparison with Trixi.jl

Trixi.jl offers two families of indicators that drive both shock capturing and AMR:

- **`IndicatorHennemannGassner`** — the Persson–Peraire modal-energy indicator with the
  next-shell robustification, mapped through a logistic function to a blending coefficient
  \(\alpha\in[0,1]\). SELF's indicator uses the **same modal-energy quantity** (§2.1); it stops
  at \(\sigma_e = \log_{10} S_e\) and thresholds directly rather than forming \(\alpha\),
  because AMR needs a discrete refine/keep/coarsen decision, not a continuous blend.
- **`IndicatorLöhner`** — a normalized second-difference (curvature) estimate of a nodal
  quantity. It is cheaper (no modal transform) but noisier and less tightly coupled to spectral
  resolution. It is a natural future alternative behind the same `flag` interface.

Two design points worth noting for anyone extending this:

- Trixi additionally offers a **`IndicatorMax`**-style controller and clip/smoothing passes
  across element neighbors. SELF's indicator is strictly element-local; neighbor smoothing (to
  avoid isolated refined elements) belongs in the Stage-4 balancing step, not the trigger.
- Trixi's `ControllerThreeLevel` maps indicator values to target refinement **levels** with
  hysteresis. SELF's two-threshold refine/keep/coarsen flag is the minimal equivalent; a
  level-target controller can be layered on top once the tree (Stage 2) exists.

An alternative spectral trigger is the **Mavriplis (1994)** decay-rate estimator, which
least-squares-fits \(\lvert\hat u_p\rvert \sim c\,e^{-\sigma p}\) to the modal coefficients and
estimates the truncation error from the fitted decay rate and tail. It is more informative for
`p`-adaptivity (choosing *how much* to refine) but more fragile than the energy-fraction
measure; it is a candidate for a future `p`/`hp`-adaptive extension.

---

## 4. Staged plan for the remaining machinery

The current mesh (`Mesh2D_t`) is **statically allocated**: `elemInfo(1:6,1:nElem)`,
`sideInfo(1:5,1:4,1:nElem)`, `nodeCoords(...,1:nElem)`, and a hand-built `mortarInfo` table.
Dynamic AMR requires making the element set mutable while preserving every invariant the solver
relies on. The following stages are each independently reviewable and testable.

### Stage 2 — `h`-refinement mesh mutation

- **(done)** Generate child `nodeCoords` from the parent geometry map by exact isoparametric
  subdivision, and build refined connectivity deterministically (`SubdivideNodeCoords`,
  `RefineConnectivity`). Metric terms/Jacobians for children come from the existing geometry
  routines with the refined `nodeCoords` as input — **no change to the geometry algorithms**.
- **(done)** Uniform (conforming) refinement end to end (`UniformRefineMesh`), serial.
- **(2b, done)** An explicit **quadtree / forest-of-quadtrees** parent–child structure
  (`SELF_QuadTreeMesh_2D`): each leaf is an active element; adaptive refinement replaces a flagged
  leaf with four children (each spanning a reference sub-quadrant), coarsening merges four
  siblings back to their parent, with amortized-capacity node growth and traversal-based leaf
  enumeration. Driven directly by the Stage-1 indicator flags via `AdaptFromFlags`.
- **(next)** Storage compaction: reclaim nodes orphaned by coarsening (currently the node pool
  grows monotonically).

### Stage 3 — Solution transfer (prolongation / restriction) — **done**

Implemented in `SELF_SolutionTransfer_2D` (see §2.7):

- **Prolongation** (parent → 4 children): tensor product of the mortar restriction operator
  `Lagrange%mortarR` — exact interpolation of the parent polynomial onto the children.
- **Restriction** (4 children → parent): tensor product of the \(L^2\)-projection adjoint
  `Lagrange%mortarP`; conservative by construction (\(\sum_k P_k R_k = I\) and discrete
  conservation), so coarsening preserves cell-integrated quantities.
- The transfer is a separate operator applied between time steps and does not touch the solver's
  floating-point reductions.

### Stage 4 — Mortar regeneration and 2:1 balancing

- **(4a, done)** **Face-neighbour navigation** on the forest (`FaceNeighbor`, ascend/descend
  quadtree search across the base root connectivity, honouring base side pairings and flips) and
  **2:1 balance** (`Balance2to1`: neighbours more than one level coarser than a leaf are refined,
  rippling to a fixed point). See §2.8.
- **(4b, done)** `EmitMesh` (`SELF_AdaptiveMesh_2D`) rebuilds `sideInfo` + `mortarInfo` from the
  balanced tree and emits a solver-ready `Mesh2D_t`: every face where `FaceNeighbor` sees a
  finer/coarser neighbour becomes a 2:1 mortar (the configuration the solver already handles);
  same-level faces are conforming; leaf geometry comes from `LeafCoords`. See §2.9.
- **(next)** Optional **neighbour smoothing** of the trigger flags (avoid isolated refined
  elements), on top of the balance pass.

### Stage 5 — MPI dynamic re-partitioning

- Refinement changes per-rank element counts and thus load balance. Introduce a
  **space-filling-curve (Morton/Hilbert) ordering** of leaves and repartition by equal-weight
  arc segments after each adaptation, matching the existing domain-decomposition ownership
  model.
- Migrate element data for reassigned leaves, then rebuild the halo-exchange tables. Preserve
  the existing rank-local ownership and halo patterns — the exchange kernels are unchanged; only
  the tables they consume are rebuilt.

### Stage 6 — GPU device re-allocation

- On adaptation, device arrays sized by `nElem` must be reallocated and re-uploaded. Reuse the
  amortized-capacity scheme from Stage 2 so device reallocation is infrequent, and copy
  survivor data device-to-device where possible to avoid host round-trips.
- Keep memory access patterns identical to the static case so kernel performance is unaffected
  between adaptation steps.

### Driver integration

Adaptation runs **between** time steps (e.g. every _k_ steps): estimate → flag → (optional
neighbor smoothing) → refine/coarsen tree → transfer solution → regenerate mortars/balance →
repartition (MPI) → re-upload (GPU). The time-integration loop itself is unchanged; the solver
sees only a (possibly) different, still-valid mesh at the top of the next step.

---

## 5. Application roadmap: `LinearEuler2D` ultrasound point source under AMR

This section maps the serial AMR machinery above onto a first *running application*: a single
point-source wavelet in the ultrasound frequency range, propagating through a ~1 m × 1 m
domain with the 2-D linear Euler model, on a dynamically adapting mesh. Serial CPU and single
GPU are the initial targets (Stage 5 MPI repartitioning is not required). Each gap below is an
additive, independently testable piece; existing model/mesh interfaces stay untouched.

### 5.1 What already works (no changes needed)

- **Flux coupling on adaptive meshes** — `EmitMesh` produces the same `mortarInfo` layout the
  `LinearEuler2D` mortar tests (`test/lineareuler2d_mortar_soundwave.f90`) already exercise.
- **Per-epoch time step** — `ForwardStep(tn, dt, ioInterval)` takes `dt` on every call, so a
  driver that re-computes `dt` after each adaptation needs no time-integrator changes.
- **Output of a changing mesh** — every HDF5 snapshot written by `WriteModel` carries its own
  `/controlgrid/geometry` alongside the solution, so per-snapshot meshes are already
  representable; `pyself` reads them file-by-file.
- **Initial condition** — `SphericalSoundWave` (Gaussian pressure pulse) generates an outgoing
  wavelet whose spectral content is set by the pulse half-width `Lr`; choosing `Lr` of a few
  millimetres puts the dominant wavelength in the ultrasound band with **zero model changes**.

### 5.2 Gap 1 — Transfer plan: old-leaf → new-leaf solution mapping — **implemented**

*Status: implemented in `SELF_TransferPlan_2D` (see §2.10), validated by
`test/transfer_plan_2d.f90`.*

`AdaptFromFlags` + `Balance2to1` mutate the forest but record no correspondence between the
pre- and post-adaptation leaf lists, which the Stage-3 transfer operators need. Because node
ids are stable in the pool (coarsening only orphans nodes), the plan can be built after the
fact from a snapshot of the old leaf array:

- new leaf **is** an old leaf → copy;
- new leaf **descends from** an old leaf → prolong along the quadtree path
  (`Balance2to1` ripple can refine a fresh child again, so prolongation must handle **multiple
  levels**, applying `ProlongToChildren` per step of the path);
- new leaf **is the parent of four** old leaves → restrict (always exactly one level:
  `AdaptFromFlags` coarsens one level per call and balancing never coarsens).

Deliverable: a `BuildTransferPlan` (forest + saved old-leaf list → typed plan) and an
`ApplyTransfer` driver mapping `solution%interior(:,:,oldIdx,:)` to the new element ordering.
Tests: adapt→transfer conservation of `∫u dA` (Jacobian-weighted), refine-then-coarsen
reversibility through a full plan, and a balanced two-level ripple case.

### 5.3 Gap 2 — Model regrid: rebinding a live `DGModel2D` to an emitted mesh

`DGModel2D` storage (7 `MappedScalar/Vector` objects) is sized by `nElem` at `Init`, and
`Init` is `intent(out)` — it resets `t`, model parameters (`rho0`), BC registrations, and the
IO counter. Rather than making the model mutable, add an **external AMR controller** module
(e.g. `SELF_AMRController_2D`) that owns the forest, the indicator, and double-buffered
`Mesh2D`/`SEMQuad` instances, and performs one adaptation epoch:

1. `indicator%Estimate` on the current solution (driving variable: pressure, `ivar=3`);
2. **flag halo expansion** — grow refine flags to face-neighbours of flagged leaves (via
   `FaceNeighbor`) so the moving wavefront cannot outrun the refined band between epochs;
   this is the "neighbour smoothing" already anticipated in §4;
3. snapshot old leaves → `AdaptFromFlags` → `Balance2to1` → `BuildTransferPlan` → `EmitMesh`;
4. new `SEMQuad` geometry (`Init` + `GenerateFromMesh`), free the old buffer;
5. save model scalars (`t`, `rho0`, integrator choice, flags), `Free` + `Init` the model on
   the new mesh/geometry, restore scalars, apply the transferred solution, `UpdateDevice`.

Step 5 works on GPU today because a fresh `Init` allocates correctly sized device arrays;
host-side transfer with an `UpdateHost`/`UpdateDevice` round-trip per epoch is acceptable at
demo cadence (Stage-6 device-side transfer remains the later optimization). Background fields
`c` (var 4) and `rho0` (var 5) ride the same prolong/restrict — exact for uniform media.

### 5.4 Gap 3 — Time-step control

- **Now (required):** level-based global `dt`. For a quadtree, the fine-level element scale is
  exactly `h_root / 2^maxLevel`, so `dt_epoch = dt_base / 2^maxLevel` with `dt_base` chosen
  for the base mesh by the usual explicit-DG bound `dt ≈ C·h/(c·N²)`. Deterministic, free,
  and no geometry reduction is needed.
- **Later (optional):** local time stepping (LTS) — leaves at level ℓ subcycle with
  `dt/2^ℓ`. This touches the RK update and requires time-interpolated interface/mortar data
  between levels, i.e. exactly the time-integration and flux-exchange machinery that is
  frozen by policy; it needs its own design + review round (call it Stage 7). Cost analysis
  for this demo says it is not needed to start: with the refined band confined to the
  wavefront annulus, a global fine `dt` costs `nElem_total × fine-step-count`, and most
  elements are coarse *and cheap*; LTS buys roughly `2^maxLevel×` on the coarse bulk — worth
  having, not blocking.

### 5.5 Gap 4 — The example and its CI-scale test

`examples/linear_euler2d_amr_pointsource.f90` (plus a reduced `test/` variant):

- Domain `[0,1]²` m via `mesh%StructuredMesh`, radiation BCs on all four sides; source at the
  centre.
- Medium: water (`c = 1500 m/s`, `rho0 = 1000 kg/m³`) with `f₀ ≈ 100 kHz` → `λ = 15 mm`
  (an air / 40 kHz variant, `λ ≈ 8.6 mm`, also fits but needs one more refinement level).
- Resolution: base 16×16 (`h₀ = 62.5 mm`), `N = 7`, max level 3 (`h = 7.8 mm`), giving
  ≈ 15 points per wavelength on the fine level — comfortable for wave propagation; the coarse
  bulk intentionally under-resolves the front so the indicator *must* refine to keep σ below
  threshold.
- `dt ≈ 3×10⁻⁸ s` at level 3; an end time of ~0.3 ms (front travels 45 cm) is ~10⁴ steps —
  seconds-to-minutes serial CPU, trivial on one GPU.
- Adaptation cadence: regrid every k steps with `k·dt·c ≤` one fine element (`k ≈ 100` at the
  numbers above, with the §5.3 halo providing the safety margin).
- CI assertions: solution NaN-free; entropy finite and non-increasing (upwind flux +
  radiation BCs are dissipative); refinement actually occurs (`forest%MaxLevel() > 0`, leaf
  count grows) **and** coarsening occurs behind the front (leaf count later shrinks);
  Jacobian-weighted transfer conservation defect at machine precision per epoch.

### 5.6 Gap 5 — Visualization: pressure field + mesh skeleton

A companion `examples/linear_euler2d_amr_plot.py` (pyself + matplotlib/pyvista) that, per
snapshot: renders the pressure field from `/controlgrid/solution` and overlays the **element
wireframe** traced from each element's four edges in `/controlgrid/geometry`, then assembles
PNG frames into a movie. Because each file carries its own geometry, frames with different
element counts need no special handling. This is the artifact that *shows* the refinement
band tracking the expanding wavefront.

### 5.7 Deferred / follow-on

- **Time-dependent transducer source.** A true point *forcing* (e.g. Ricker wavelet at `f₀`)
  rather than an initial pulse: `source2d` currently has no access to position or time, so
  this needs a localized-forcing hook plus per-epoch source relocation (the containing
  element changes identity on regrid). Physically nicer (continuous-wave and pulse-train
  experiments); not required for the first demo.
- **LTS (Stage 7)** as scoped in §5.4, and **device-side transfer (Stage 6)**.
- **Storage compaction** of orphaned forest nodes on long runs (§4, Stage 2 "next").

### 5.8 Suggested PR sequence

| PR | Content | Depends on |
| --- | --- | --- |
| 1 | Transfer plan (old→new leaf map, multi-level prolong) + tests | — |
| 2 | AMR controller (halo flags, regrid orchestration, level-based dt) + adapt-epoch soundwave test | 1 |
| 3 | Ultrasound example, CI-scale test, plotting script, docs | 2 |
| 4 | GPU epoch test; device-side transfer if profiling justifies it | 2 |
| 5 | (design first) LTS; time-dependent point forcing | 3 |

---

## 6. References

- P.-O. Persson and J. Peraire, *Sub-cell shock capturing for discontinuous Galerkin methods*,
  AIAA 2006-112 (2006).
- S. Hennemann, A. M. Rueda-Ramírez, F. J. Hindenlang, and G. J. Gassner, *A provably entropy
  stable subcell shock capturing approach for high order split form DG for the compressible
  Euler equations*, J. Comput. Phys. 426, 109935 (2021).
- C. Mavriplis, *Adaptive mesh strategies for the spectral element method*, Comput. Methods
  Appl. Mech. Engrg. 116 (1994) 77–86.
- M. Schlottke-Lakemper et al., *Trixi.jl: Adaptive high-order numerical simulations of
  hyperbolic PDEs in Julia*, and the Trixi.jl `IndicatorHennemannGassner` / `IndicatorLöhner`
  documentation.
