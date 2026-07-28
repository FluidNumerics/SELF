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
| Face-neighbour queries + 2:1 balancing + hanging-node/mortar emission | Designed (Stage 4) |
| MPI dynamic re-partitioning / load balancing | Designed (Stage 5) |
| GPU device re-allocation for a changing element count | Designed (Stage 6) |

The trigger and the conforming (uniform) `h`-refinement path are wired into the library today
(see §2.5). Adaptive refinement with hanging nodes reuses the same subdivision primitive but
additionally needs the mortar regeneration and 2:1 balancing of Stage 4; dynamic MPI
re-partitioning is Stage 5. These remain deferred behind this design so each lands as a
self-contained, reviewable piece.

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

- Add **face-neighbour navigation** on the forest (ascend/descend quadtree search across the base
  mesh's root connectivity, honouring base side pairings and flips). This is the prerequisite for
  both balancing and mortar detection and is the natural first piece of Stage 4.
- After a refine/coarsen sweep, rebuild `mortarInfo` from the tree: every face between elements
  of different levels becomes a 2:1 mortar (the configuration the solver already handles).
- Enforce **2:1 balance** (no face may separate elements differing by more than one level) by
  propagating refinement: neighbors of a refined element that would violate the balance are
  themselves refined. This is where optional **neighbor smoothing** of the trigger flags lives.
- Reuse the existing conforming-side connectivity generation for same-level faces.

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

## 5. References

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
