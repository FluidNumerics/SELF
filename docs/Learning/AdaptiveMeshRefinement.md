# Adaptive Mesh Refinement (AMR)

This page describes SELF's approach to adaptive mesh refinement for discontinuous Galerkin
spectral element models. Sections 1-5 document the **2-D** implementation, which came first: the
refinement trigger, and the staged design along which the invasive parts (dynamic mesh mutation,
solution transfer, MPI re-partitioning, and GPU re-allocation) were reviewed and landed.
[Section 6](#6-three-dimensional-amr-octrees-face-mortars) covers the **3-D** stack - octrees and
2:1 face mortars - which is a deliberate transcription of the 2-D one and records only what
differs.

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
| Model regrid (`DGModel2D%Regrid`) + AMR controller (serial, CPU/GPU) | **Implemented** |
| Ultrasound point-source example + AMR visualization script | **Implemented** |
| MPI dynamic re-partitioning / load balancing (v2: replicated forest, point-to-point migration) | **Implemented** |
| GPU device re-allocation for a changing element count | **Implemented** (amortized high-water-mark storage, Stage 6b) |
| Device-side solution transfer (no host round trip, any rank count) | **Implemented** (Stage 6a; multi-rank in §6.7) |
| Geometry reuse for unchanged elements across an epoch | **Implemented** (Stage 6c) |

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

### 2.1.1 The amplitude gate (relative energy floor)

\(S_e\) is a **ratio**, so it is scale-free by construction: it reports how an element's energy is
distributed across modes and says nothing about how much energy the element carries. Grid-scale
residue at \(10^{-5}\) of the field's amplitude has exactly the same \(S_e\approx 1\) as a
full-amplitude discontinuity.

Guarding the ratio with an *absolute* floor at machine epsilon does not fix this, because
\(E_{\text{tot}}\) carries the **square** of the field amplitude while \(\varepsilon\) is a fixed
\(2.2\times10^{-16}\): a wake at \(10^{-5}\) of the peak sits some ten orders of magnitude above the
floor and is still judged on shape alone. Such elements never reach
\(\sigma_{\text{coarsen}}\), so the mesh behind a passing front is never released and the refined
region grows to the whole area the wave has swept (issue #162).

Each element therefore also carries a **gate energy**

$$
g_e = \sum_v w_v\,E_{\text{tot},e,v},
$$

a weighted sum of the per-variable modal energies. Since \(E_{\text{tot},e,v}\) is the exact \(L^2\)
energy of variable \(v\) on the reference element, \(g_e\) is a convex quadratic functional of the
state — a **discrete entropy integral** over the element when \(w\) is taken from the diagonal of an
entropy Hessian. It is compared against a field-wide energy scale,

$$
g_e \;\le\; \max\!\left(\varepsilon,\; \texttt{relativeEnergyFloor}\cdot \texttt{energyScale}\right)
\quad\Longrightarrow\quad S_e := 0 ,
$$

so a quiescent element is reported perfectly resolved (\(\sigma_e = \log_{10}\varepsilon\), hence
`COARSEN`) whatever its modal shape.

Because energy is amplitude squared, `relativeEnergyFloor` is \(10^{dB/10}\) in amplitude terms.
The default is \(10^{-12}\): amplitudes below \(10^{-6}\) (−120 dB) of the field scale count as
quiescent. It sits below \(\varepsilon\) in `real32`, so in a single-precision build the absolute
term dominates and the gate is inactive unless a caller raises it — the safe direction, since a
`real32` field cannot represent −120 dB structure anyway.

#### The energy hysteresis band

As described so far the gate is a **single hard cut**, which reintroduces on the energy axis exactly
the thrashing the two \(\sigma\) thresholds exist to prevent: an element whose energy drifts across
that one value flips `COARSEN` ↔ `REFINE` on successive epochs, churning the mesh without the
solution changing meaningfully. `significantEnergyFloor` opens a band instead:

| gate energy \(g_e\) | flag |
| --- | --- |
| \(\le\) quiescent floor | `COARSEN`, whatever the spectrum says |
| between the floors | `KEEP` — hold the mesh |
| \(>\) significant floor | the spectrum decides, as before |

The middle zone reads *too weak to be worth spending levels on, too strong to declare resolved*,
and is stable under a small change in amplitude. Inside the band \(\sigma_e\) still reports the
**true** spectrum — the band holds the flag, it does not falsify the diagnostic — so for those
elements `flag` is deliberately not what thresholding \(\sigma_e\) would give.

It defaults to the quiescent floor, collapsing the band to the single cut, so behaviour is
unchanged unless a caller opts in. `test/refinement_indicator_2d_energy_hysteresis.f90` covers the
three zones and asserts that an in-band energy drifting by a factor of ten either way does not move
the flag, while the same drift across a degenerate band does.

**Keep the band narrow.** Its upper edge is functionally *"do not spend levels below this energy"* —
the same knob as the floor, with the same cost in refinement depth. On the ultrasound benchmark's
initial adaptation, with a quiescent floor of \(10^{-12}\):

| significant floor | 1e-12 | 3e-12 | 1e-11 | 3e-11 | 1e-10 |
| --- | --- | --- | --- | --- | --- |
| elements / max level | 328 / 2 | 328 / 2 | **268 / 1** | 268 / 1 | 268 / 1 |

A band up to ~3× the floor leaves depth untouched; at 10× it collapses a level, exactly as raising
the floor to \(10^{-10}\) does. The band is a thrash damper, not a savings knob — widen it only as
far as is needed to stop flags oscillating.

**Raising the floor is not free.** The gate cannot distinguish residue from the low-amplitude
*flank* of a feature that is genuinely under-resolved, so an aggressive floor buys element count at
the cost of refinement **depth**. Measured on the ultrasound point-source benchmark
(`examples/linear_euler2d_amr_ultrasound_pointsource`, 30 epochs, `maxLevel = 2`):

| `relativeEnergyFloor` | initial adaptation | mean elements | final |
| --- | --- | --- | --- |
| 0 (no gate) | 328, level 2 | 1482 | 1732 |
| **1e-12 (default)** | **328, level 2** | **1057 (1.40x fewer)** | **1528** |
| 1e-8 | 268, level **1** | 2210 | 3304 |

At `1e-8` the gate also suppresses the source pulse's skirt, so the forest never reaches the level
cap; `RecommendedTimeStep` is tied to the max level, so `dt` doubles, and the resulting
time-integration error drives enough later refinement to roughly *double* the final mesh. Measure
before raising this. A related trap: a 2-D cylindrical pulse leaves a genuine algebraically-decaying
tail behind its front (Huygens' principle fails in even dimensions), which is physics rather than
residue and should not be gated away.

`energyScale` defaults to the largest \(g_e\) over the elements, reduced with `MPI_MAX` when the
indicator is given a communicator, so the gate follows a decaying front. Two consequences worth
knowing:

- the automatic scale makes the flags depend on a floating-point reduction, which can differ at
  round-off between rank counts — `SetEnergyScale` pins it and is the deterministic escape hatch;
- once a wave has left the domain entirely the scale collapses onto the residue's own peak and the
  gate stops discriminating; the absolute \(\varepsilon\) floor is what bounds that case.

Note that raising \(\sigma_{\text{coarsen}}\) is **not** a substitute. In a low-amplitude wake
\(\sigma_e\) is near \(0\), so any threshold high enough to coarsen the wake is also high enough to
stop the front from refining; the two decisions cannot be separated by thresholds alone, which is
what motivates an amplitude-based gate.

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

! Amplitude gate (all optional; these are the defaults unless set).
call amr%SetRelativeEnergyFloor(1.0e-12_prec) ! 0 disables the gate
call amr%SetRelativeEnergyFloor(1.0e-12_prec, &  ! ... or open a hysteresis band
                                significantEnergyFloor=1.0e-10_prec)
call amr%SetEnergyScale(pRef**2)              ! pin the scale; ClearEnergyScale undoes it
call amr%SetEnergyWeights(w)                  ! per-variable gate weights, e.g. from the entropy

! solution : the model's MappedScalar2D (or any Scalar2D); ivar selects the driving
! variable, or SELF_AMR_ALLVARS (=0) reduces over all variables (most conservative).
! comm : optional, makes the automatic energy scale global (MPI_MAX).
! gate : optional, a caller-computed per-element gate energy replacing the weighted sum.
call amr%Estimate(solution, ivar=1)

! amr%indicator(1:nElem)  -- per-element sigma_e (host)
! amr%flag(1:nElem)       -- per-element SELF_AMR_REFINE / KEEP / COARSEN (host)
! amr%gate(1:nElem)       -- per-element gate energy actually used (host, diagnostic)
n = amr%CountFlagged(SELF_AMR_REFINE)
```

The estimate runs in two phases: an element-local, parallel (device) pass that forms the spectra,
the raw \(S_e\) and \(g_e\); then a host pass that reduces for the energy scale, applies the gate,
takes the \(\log_{10}\) and sets the flags. The second phase is host-side because the scale needs a
reduction that no element-local pass can supply, and it is shared by both backends, so CPU and GPU
produce identical flags. `AMRController2D` forwards `relativeEnergyFloor` and `energyWeights` from
its own `Init` and re-applies them after each epoch's indicator resize; it also supplies the
communicator so the scale is global.

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

## 2.11 Model regrid and the AMR controller (implemented)

Two pieces close the loop around a *live, time-stepping model*:

- **`DGModel2D%Regrid(mesh, geometry)`** rebinds a model to a new mesh/geometry pair: the
  mesh-sized solution storage is reallocated and the boundary-condition registrations and maps
  are rebuilt (mirroring the mesh-sized portion of `Init`/`Free`, including the GPU backend's
  BC device arrays), while everything that is not mesh-sized is preserved — the time state
  (`t`, `dt`, entropy, IO counter), the time-integrator selection, configuration flags, and
  model-specific parameters, all of which a fresh `Init` (`intent(out)`) would reset. The
  solution interior is left for the caller to fill via the §2.10 transfer.

- **`SELF_AMRController_2D`** owns the forest, the indicator, and the meshes/geometries it
  emits (double-buffered), and performs one adaptation epoch per `Adapt(model, adapted)` call:
  estimate → cap refine flags at a configurable `maxLevel` → spread refine flags to face
  neighbours for `nHalo` passes (so a feature moving at speed \(c\) stays inside the refined
  band when the adaptation cadence satisfies \(k\,\Delta t\,c \le n_{halo} h_{fine}\)) →
  `AdaptFromFlags` + `Balance2to1` (a no-op epoch leaves the model untouched) →
  `BuildTransferPlan` + `EmitMesh` + new `SEMQuad` → `model%Regrid`, apply the transferred
  solution, upload to device. `RecommendedTimeStep(dtBase) = dtBase / 2^{MaxLevel}` gives the
  level-based explicit-stability bound to pass to `ForwardStep` after each epoch — exact for
  the quadtree, whose children are exact half-scale subdivisions.

`test/lineareuler2d_amr_soundwave.f90` runs the full loop on a deliberately under-resolved
acoustic pulse (LinearEuler2D, radiation boundaries, RK3): the initial adaptation refines
around the pulse up to the level cap; every mid-run adaptation conserves the Jacobian-weighted
global integral of each prognostic variable to roundoff; the model's time and parameters
survive regridding; the mesh evolves as the wave propagates; and the acoustic energy stays
finite and non-increasing across the whole adaptive run.

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

### Stage 5 — MPI dynamic re-partitioning — **implemented (v2)**

*Status: implemented. `QuadTreeMesh2D%InitGlobal` builds the rank-replicated forest from
allgathered global base tables; `EmitMesh` builds the global connectivity/mortar tables on
every rank and stores only its contiguous slice of a freshly generated decomposition
(`sideInfo(3)` global ids, global `nUniqueSides`, fully replicated `mortarInfo`, exactly the
invariants `SideExchange`/`MortarExchange` require); the controller allgathers the indicator
flags per epoch, and migration is point-to-point: each rank receives only the window of the old
field its own new element range reads (5d below). Validated by
`test/lineareuler2d_amr_soundwave_mpi.f90` at 2 and 4 ranks: the global element trajectory and
entropy history match the serial run, transfers conserve globally, and a leaf-list checksum
confirms forest replication — plus `SELF_AMR_MIGRATE_VERIFY=1` re-runs, which assert the
migrated window matches the v1 allgathered field value for value, and the serial
`test/transfer_plan_{2,3}d_window.f90` routing tables. The v1 gather-then-slice path is retained
behind `SELF_AMR_MIGRATE_GATHER=1`. What remains open is the distributed forest, whose memory is
still O(global elements) per rank (issue #167 Stage 2).*

Two observations make a correct first version tractable:

- The forest's leaf list is **already a space-filling curve**: root-major depth-first traversal
  is Morton order within each quadtree, so "SFC partitioning" is just contiguous ranges of the
  existing leaf list — the same contiguous-ownership model (`offsetElem`) the domain
  decomposition already uses.
- The forest is cheap (a few integers per node), so it can be **replicated on every rank**.
  If all ranks apply identical flags, they compute identical adapted/balanced forests, transfer
  plans, and emitted global connectivity — deterministically, with no communication beyond the
  flags themselves.

Sub-stages:

- **(5a) Global flags + replicated mutation.** The controller allgathers the rank-local
  indicator flags (by the decomposition's element offsets) into a global per-leaf flag array;
  every rank then runs the same cap/halo/`AdaptFromFlags`/`Balance2to1` sequence on its forest
  copy. One small collective per epoch, outside the time-stepping loop.
- **(5b) Multi-rank `EmitMesh`.** Every rank builds the same global `sideInfo`/`mortarInfo`
  (deterministic from the forest), then decomposes it exactly the way the existing mesh
  constructors do for `nRanks > 1`, so `SideExchange`/`MortarExchange` consume the result
  unchanged. Repartitioning is implicit: each epoch's emitted mesh is re-decomposed over the
  new leaf list, so equal-count SFC arcs move with the refinement.
- **(5c) Solution migration, v1 — gather-then-slice, retained as a fallback.** Allgatherv the
  old rank-local solutions into a global old field, then `ApplyTransferPlan` only for the new
  rank-local elements. Correct and simple; memory and per-rank traffic are one global solution
  copy per adapting epoch, which is fine at single-node scale and prohibitive beyond it. Still
  reachable with `SELF_AMR_MIGRATE_GATHER=1`, and exercised in CI so it cannot rot.
- **(5d) Solution migration, v2 — point-to-point windows (issue #167).** The default. Four
  things make it simple:

    1. *The window.* The old elements a rank's new range `[eFirst,eLast]` reads through the plan
       — `sourceElem`, or all four (2-D) / eight (3-D) `family` members of a coarsened family —
       lie in a contiguous window `[wFirst,wLast]` of the old element list, because both
       partitions are contiguous ranges of the same leaf order. `PlanWindows` takes the min/max
       hull of those indices in one pass over the plan. Correctness does not depend on the
       ordering being monotone: a non-monotone ordering only widens the hull, and the worst case
       is the whole old list — exactly what v1 moves. Locality, not correctness, is what the
       space-filling curve buys.
    2. *Routing needs no communication.* The plan and both `offsetElem` tables are replicated,
       so a rank computes its own window **and every peer's**, and derives its send runs
       (its old range ∩ a peer's window) and receive runs (its window ∩ a peer's old range)
       locally. Both ends of a pair evaluate the same `max`/`min` on the same integers, so the
       schedules match by construction — no count exchange, no handshake, no collective.
    3. *Nothing is packed.* A run of elements at a fixed variable is contiguous in both
       `solution%interior` and the window buffer, so each `MPI_Isend`/`MPI_Irecv` reads and
       writes the real storage: one message per (peer, variable), a handful of peers, no pack
       buffer and no unpack loop. (A rank talks to a couple of peers for a balanced repartition;
       nothing bounds it to that, and a pathological repartition simply produces more, smaller
       messages.) Receives are posted first and the exchange closes with one
       `MPI_Waitall`, the shape every other SELF exchange uses. It runs *before* `Regrid`, so
       the sends read the still-live pre-regrid solution and no traffic is in flight while
       `Regrid` rebuilds mesh state on the same communicator.
    4. *`ApplyTransferPlanWindow`* is what consumes it: `ApplyTransferPlanRange` (v1) is that
       routine over the whole old field, so the two paths share one body and are **bit-identical**
       — the same operands, the same operators, the same order, no reductions. That is a
       guarantee, not a tolerance, and `SELF_AMR_MIGRATE_VERIFY=1` asserts it in-process. On a GPU
       build the consumer is the kernel instead, reading a device-resident window with an
       `oldFirst` offset (§6.7); keep the two questions apart when reading a verification failure.
       **Migration** is byte movement — device-to-device copies, MPI transfers — so it stays exact
       and `MIGRATE_VERIFY` stays bitwise. The **apply** is arithmetic, and on the device agrees
       with the host only to round-off, because the compiler contracts its multiply-accumulates
       into FMAs; that has its own switch, `SELF_AMR_TRANSFER_VERIFY`, and its own tolerance.
       Do not merge them — a single switch would have to take the looser bar.

  A rank owning no new elements gets an empty window, posts no receives, and skips the apply,
  but still serves its peers. A coarsened family whose children were owned by different ranks is
  fully covered by the hull and projected on the receiving rank in child-octant order, so the
  result does not depend on the partition. The per-rank migration footprint becomes the window —
  local elements plus the overhang — instead of a full global field, in a controller-owned
  grow-only buffer, so a settled adapting run allocates nothing here.
- Tests on ≥ 2 ranks: forest determinism across ranks (identical leaf checksums after an
  epoch), global conservation of the transferred solution (mpi_allreduce), and the AMR
  soundwave regression run distributed.

### Stage 6 — GPU device re-allocation (implemented)

Both parts are implemented and measured on one MI300X (gfx942, ROCm 6.4.3, exclusive node),
ultrasound example, 20 epochs, mean of three runs. The starting point was the "correct but
unamortized" form: every adapting epoch freed and re-initialized all model storage, re-uploaded
mesh/geometry, and moved the solution through a host round trip.

| | adaptation total | AMR share of epoch loop |
| --- | --- | --- |
| before | 0.998 s | 51.7% |
| + 6a device-side transfer | 0.830 s | 47.0% |
| + 6b amortized capacity | **0.616 s** | **39.4%** |

Time integration is unchanged throughout (0.93-0.95 s), so the whole gain is in adaptation:
**38% off the cost of an epoch**. On the curved multi-material bone-and-marrow case the AMR
share falls from 19.7% to 14.2%.

- **(6a) Device-side transfer — implemented.** `StageSolutionForTransfer` and
  `ApplyTransferPlan` are type-bound on the model, using the same backend split as `Regrid`.
  The GPU override stages the pre-regrid solution device-to-device and applies the plan in
  `TransferSolution_2D_gpu` (`src/gpu/SELF_SolutionTransfer.cpp`), so a single-GPU adapting run
  moves no solution data across the host link. The kernel applies only the mortar operator pair
  of the child on the recorded path rather than prolonging to all four and discarding three, so
  device and host agree to round-off while conservation stays exact. Multi-rank runs took the
  host path until #172, which assembles the migrated window in device memory and gives the kernel a
  window offset, so the device transfer now serves any rank count (§6.7).
- **(6b) Amortized capacity — implemented.** The prediction recorded here previously was that
  "at typical cadences allocation cost is expected to be far below the re-upload and
  geometry-generation cost". **The measurement contradicted it:** the free/re-initialize cycle
  was the single largest component of an adaptation at 51.5%, ahead of both the transfer and
  geometry regeneration. The cost was not the device allocator (`hipMalloc` + `hipFree` were
  only ~8% of an adaptation) but the host-side work a fresh allocation drags along - allocating,
  zeroing, rebuilding metadata and equation parsers, and uploading the zeros. `Resize` on the
  data classes now reuses storage via rank-1 pools with pointer remapping (see
  `src/SELF_DataPool.f90` for why this avoids changing any kernel stride).

**Behaviour change to be aware of.** On a GPU build the transferred solution is left on the
device, so `solution%interior` (the host mirror) is stale after `Adapt`. This matches the rest
of the time loop, where the device is authoritative and a caller wanting host data calls
`UpdateHost()` first, as `Write_DGModel2D_t` does. It used to be incidentally fresh because the
transfer ran on the host.

**How this scales with refinement depth.** Depth is the knob that reaches ultrasound length
scales (the base mesh stays 16x16; `Lr` and the epoch length track `maxLevel`). Adaptation cost
grows faster than time integration, so the AMR share rises with depth - but it rises from a
lower base and more slowly than before:

| maxLevel | f0 | elements | baseline | + 6a/6b | + 6c | geometry reused |
| --- | --- | --- | --- | --- | --- | --- |
| 3 | 0.10 MHz | 2,092 | 51.7% | 40.1% | **38.1%** | 72.8% |
| 4 | 0.20 MHz | 4,156 | 56.8% | 43.1% | **42.1%** | 73.0% |
| 5 | 0.40 MHz | 7,528 | 58.9% | 45.9% | **45.1%** | 73.9% |
| 6 | 0.80 MHz | 14,896 | 61.6% | 48.4% | **47.6%** | 74.6% |
| 7 | 1.60 MHz | 30,760 | 62.9% | 50.6% | **49.9%** | 75.4% |
| 8 | 3.20 MHz | 66,616 | *out of memory* | 51.8% | **51.0%** | 76.5% |

maxLevel 8 could not run at all before this work: it exhausted the device after ~555
adaptations (see the leak below). It now completes, which is the first datapoint at a
production-relevant frequency. The reused fraction rises with depth, as expected - the refined
annulus is a thinner shell relative to the whole mesh - so geometry reuse pays off more, not
less, at production resolution. Per-epoch element counts are identical at every depth before and
after 6c, which is the check that the incremental geometry changed no physics.

**A leak this work exposed.** `boundarynormal_gpu` was allocated in `Init_MappedScalar2D` and
never freed, leaking 12.8 kB per element per adaptation at N=7, nvar=5 across the model's five
`MappedScalar2D` fields. It aborted a 640-epoch maxLevel-8 run with an out-of-memory error
after ~555 adaptations on a 192 GiB device. Fixed, with a device-memory regression test
(`test/data_2d_device_memory.f90`) since nothing in the suite had been watching device memory.

### Stage 6c - geometry (implemented)

Re-profiling after 6a/6b put geometry at the top of an adaptation: 46.8%, ahead of `Regrid`, split
between allocate/zero/free churn (`Init_SEMQuad` + `Free_SEMQuad`, 21.6%) and real computation
(`GenerateFromMesh_SEMQuad`, 25.2%). Both are addressed:

- **Persistent geometry.** `Tensor2D` joins the amortized-storage scheme (it was the last class on
  plain `allocate`/`hipMalloc`, and `dxds`/`dsdx` are the largest arrays in `SEMQuad`), `SEMQuad`
  gains `Resize`, and the controller owns two long-lived geometry buffers instead of allocating one
  per epoch. The `nGeo -> N` interpolant and node-coordinate staging that `GenerateFromMesh` used to
  rebuild on every call are cached, since only the element count varies between epochs.
- **Incremental generation.** Most elements do not change in an epoch, and the transfer plan already
  names them: `SELF_TRANSFER_COPY` is set only when the new and old element are the *same forest
  node*. Because `LeafCoords` is a pure function of root coordinates, level and quadrant path, and
  per-element geometry generation is strictly element-local, such an element's whole geometry block
  is reproducible **bit for bit** - so it is copied forward rather than recomputed. The changed
  elements are generated compacted into a scratch geometry and scattered into place, which keeps the
  existing generation loops untouched. The two buffers alternate so the previous epoch's geometry is
  readable while the new one is assembled.
- **Multi-rank.** `sourceElem` indexes the global old element list while a rank holds only its slice,
  so reuse additionally requires the source to be locally owned - a range test against the previous
  decomposition. On one rank that is always true; on several, migrated elements are regenerated. No
  communication is added. (Geometry is therefore still *regenerated* rather than migrated when an
  element changes rank, even though its solution now migrates point-to-point.)

Measured: 72.8% of elements reused per epoch at maxLevel 3, rising to 76.5% at maxLevel 8. Adaptation
cost 0.616 s -> 0.585 s, and the bone-and-marrow case 14.2% -> 10.7% AMR share.

Correctness is checked numerically rather than argued. `test/geometry_2d_reuse.f90` requires every
kept element to match its previous-epoch source bit for bit, and adds a positive-Jacobian check plus
the discrete metric identity `sum_i d(J dsdx(i,j))/dxi_i = 0`, which the suite did not previously
test. `test/geometry_2d_reuse_mpi.f90` covers the multi-rank branch, asserting both that reuse and
regeneration actually fire and that the resulting mixed assembly equals a from-scratch generation
exactly. Three env-gated switches (`SELF_AMR_GEOM_NO_REUSE`, `SELF_AMR_GEOM_FULL`,
`SELF_AMR_GEOM_VERIFY`) allow a suspected geometry problem to be localized in one run without a
rebuild. Two more do the same for multi-rank solution migration: `SELF_AMR_MIGRATE_GATHER=1`
selects the v1 allgather path, and `SELF_AMR_MIGRATE_VERIFY=1` runs both migrations and compares
them value for value (§ Stage 5, 5d).

### A correctness bug the amortized scheme introduced

Worth recording, because the reasoning error is easy to repeat. Stage 6b removed the `UpdateDevice`
from `Resize` on the grounds that uploading freshly zeroed arrays is waste when the caller
overwrites them. That is true of the solution but not of every field, and it silently dropped an
invariant `Init` had provided: after `Init` a field's *device* buffer was defined, because the zeros
had been uploaded. After `Resize` it held whatever was previously in that memory, and a fresh
`hipMalloc` is uninitialized.

The low-storage Runge-Kutta update reads its accumulator before writing it -
`grk = rk_a*grk + dSdt` with `rk3_a(1) = 0`, where `grk` is `workSol%interior_gpu`. Multiplying by
zero annihilates any finite leftover, which is why this survived initial testing, but `0*NaN` and
`0*Inf` are NaN. A resized `workSol` buffer that happened to contain a NaN bit pattern poisoned the
solution on the next step; the wave died, the indicator found nothing to refine, and the mesh
coarsened back to the base mesh. It presented as an intermittent, allocation-history-dependent
failure - the same source on the same node passing at one time and failing at another - and no test
caught it, because none runs enough GPU adaptation epochs for a stale buffer to land on such a
pattern. `EnsureDeviceBuffer` now zeroes on resize, restoring the invariant unconditionally; a
device-side fill costs HBM bandwidth rather than a PCIe transfer, so the point of skipping
`UpdateDevice` is preserved.

The general lesson: "the caller overwrites it before reading" is a per-field claim, not a property
of an amortized-storage scheme.

**Still open.** The geometry copies and the device upload are both still done the simple way: element
blocks are copied on the host (roughly eleven array-slice copies per element), and the whole geometry
is re-uploaded each epoch regardless of how little changed. That is why 6c yields ~5% rather than the
~35% a naive reading of the 46.8% attribution would predict. The fix is a device-side indexed-copy
kernel - structurally the 6a transfer kernel without the operator - plus uploading only the compacted
changed elements. There is now a measurement justifying it rather than an estimate. Mesh emission
(`EmitMesh`, ~5% of an adaptation) is also still a full rebuild per epoch.

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

### 5.3 Gap 2 — Model regrid: rebinding a live `DGModel2D` to an emitted mesh — **implemented**

*Status: implemented as `DGModel2D%Regrid` + `SELF_AMRController_2D` (see §2.11), validated by
`test/lineareuler2d_amr_soundwave.f90`. The level-based time step of §5.4 ("now") is
`RecommendedTimeStep`.*

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

### 5.5 Gap 4 — The example and its CI-scale test — **implemented**

*Status: implemented as `examples/linear_euler2d_amr_ultrasound_pointsource.f90` (water,
`c₀ = 1500 m/s`, `f₀ ≈ 100 kHz`, 16×16 base at `N = 7`, level-3 cap, radiation boundaries).
The example is registered as a CI test at a 6-epoch (30 µs) default;
`SELF_AMR_ULTRASOUND_EPOCHS=60` extends it to a full-domain movie run. The generic AMR-loop
mechanics are separately covered by `test/lineareuler2d_amr_soundwave.f90` (§2.11).*

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

### 5.6 Gap 5 — Visualization: pressure field + mesh skeleton — **implemented**

*Status: implemented as `examples/linear_euler2d_amr_plot.py` (h5py + numpy + matplotlib
only). Each snapshot's field and geometry are interpolated from the Gauss control points to a
uniform per-element grid including the element edges (barycentric Lagrange, exact for the
polynomial data), rendered as a filled pressure field with the element-outline wireframe
overlaid, one PNG per snapshot plus an MP4 when ffmpeg is available.*

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
| 4 | GPU epoch test; device-side transfer (Stage 6a) and amortized capacity (Stage 6b) - both done, profiling justified them | 2 |
| 5 | (design first) LTS; time-dependent point forcing | 3 |

---

## 6. Three-dimensional AMR (octrees + face mortars)

Everything above describes the 2-D implementation, which came first and remains the more
thoroughly exercised of the two. The 3-D stack is a deliberate transcription of it - octrees in
place of quadtrees, 2:1 **face** mortars in place of edge mortars, eight children in place of four
- and is implemented and validated end to end. This section records what is the same, what is
genuinely different, and what is not done yet.

The design record is [AMR (3D) Design](../Contributing/AMR3D-Design.md); it is a delta document,
and where a 3-D module is a mechanical transcription the 2-D design remains the authoritative
rationale.

### 6.1 Status

| Component | 2-D | 3-D |
| --- | --- | --- |
| Modal-decay indicator (CPU + GPU) | Implemented | Implemented (`SELF_RefinementIndicator_3D`) |
| `h`-refinement primitives, uniform refinement | Implemented | Implemented (`SELF_RefinementPrimitives_3D`) |
| Adaptive forest, level tracking, 2:1 balancing | Quadtree | Octree (`SELF_OctreeMesh_3D`) |
| Nonconforming mesh emission | Implemented | Implemented (face mortars) |
| Transfer plan + solution transfer | Implemented | Implemented (`SELF_TransferPlan_3D`) |
| Model regrid + AMR controller | Implemented | Implemented (`SELF_AMRController_3D`) |
| MPI re-partitioning (replicated forest, point-to-point migration) | Implemented (v2) | Implemented (v2) |
| Migration diagnostics (`SELF_AMR_MIGRATE_*`) | Implemented | Implemented |
| Device-side solution transfer | Implemented (Stage 6a) | **Implemented** (§6.4) |
| Device-side transfer on more than one rank | Implemented (§6.7) | Implemented (§6.7) |
| Amortized high-water-mark storage | Implemented (Stage 6b) | Inherited via the shared data classes |
| Geometry reuse across an epoch | Implemented (Stage 6c) | Implemented (`SELF_AMR_GEOM_*` switches) |
| Example | ultrasound point source | `linear_euler3d_amr_spherical_soundwave` |
| Coarsen-wake regression | Implemented | **Not yet** |

The pilot model is `LinearEuler3D`. Note it carries `nvar = 7`: `u, v, w, P` are advanced in
time, while `c`, `rho0` and `sigma` (the sponge-layer relaxation rate) are spatially varying but
time-constant background fields. They are still solution variables, so the transfer must carry all
seven - a transfer that moved only the prognostic variables would silently blank the medium on
every newly created element.

### 6.2 What is genuinely different from 2-D

- **Eight children, not four.** `transferAxc/Ayc/Azc(1:8)` map an octant to its
  `(x,y,z)` half, in CGNS corner order: children 1-4 walk the x/y quadrants of the lower-z layer,
  5-8 the same quadrants of the upper-z layer.
- **Triple tensor products.** Prolongation is
  `U_child = (R_kx ⊗ R_ky ⊗ R_kz) U_parent`; restriction is
  `U_parent = Σ_c (P_kx ⊗ P_ky ⊗ P_kz) U_child(c)`. Each 1-D `P_k` carries the half-interval
  Jacobian, so the triple product carries **1/8** rather than 1/4, which is what makes the 3-D
  restriction conservative. The identities are the same: `Restrict(Prolong(u)) = u`, and
  `Σ w³ u_parent = (1/8) Σ_c Σ w³ u_child`.
- **Cost scales far more steeply.** A level of refinement multiplies the element count by 8
  rather than 4, and each element carries `(N+1)³` nodes rather than `(N+1)²`. Adaptation is
  therefore a much larger share of a 3-D epoch than of a 2-D one at comparable settings: in the
  maxLevel-3 case measured in §6.4, adaptation is 89% of the epoch loop before this change and 88%
  after, against roughly 50% for the 2-D ultrasound case.
- **Balancing is face-based.** The forest balances across faces and tolerates 2-level edge and
  corner jumps, because DG face mortars carry all interface data and hanging edges/corners carry
  none. `EmitMesh` guards `MaxLevelJump() <= 1`.

### 6.3 The transfer protocol and where the solution lives

The plan (`SELF_TransferPlan_3D`) is built on the host after the last forest mutation, and
classifies every new leaf as `COPY`, `PROLONG` or `RESTRICT`, with a `depth` and an octant `path`
for any further descent. It is applied around the regrid as a three-step, type-bound protocol on
the model, so the backend split that already selects `Regrid` selects the transfer too:

```fortran
call model%StageSolutionForTransfer()   ! preserve the field; Regrid may then free it
call model%Regrid(newMesh,newGeom)
call model%ApplyTransferPlan(plan,interp,eFirst,eLast)
```

Staging exists because `Regrid` reallocates (or resizes) the storage the solution lives in. The
portable implementation stages into a host array whose lifetime is exactly stage-to-apply;
applying without a preceding stage is a hard error, pinned by
`test/dgmodel3d_guard_transfer_unstaged.f90`.

`ApplyTransferPlan` takes an optional `uGlobal` plus an optional `oldFirst`: the **multi-rank
path**. On more than one rank the controller hands in host-side old-field data and each rank fills
exactly its own new element range. By default `uGlobal` is the rank's *window* of the old field
and `oldFirst` is that window's first global old element index, so the plan's global
`sourceElem`/`family` indices still resolve (Stage-5 v2, § Stage 5); under
`SELF_AMR_MIGRATE_GATHER=1` it is the whole allgathered global old field, which is the same call
with `oldFirst = 1`.

Since #172 the default multi-rank path no longer uses `uGlobal` at all: the window is migrated into
model state by `MigrateOldWindow` and `ApplyTransferPlan` reads it from there, which is what lets the
GPU backend hold it in device memory and apply the plan with a kernel on any rank count (§6.7).
`uGlobal` remains for the retained `SELF_AMR_MIGRATE_GATHER=1` path and for external callers.

`[eFirst,eLast]` must be the rank's WHOLE new element range on both backends. The portable apply
writes an exactly-shaped `uNew`, so a sub-range of a larger field would place every variable after
the first at the wrong stride; the kernel is indifferent because it is told the field's element
stride separately. Do not rely on that difference.

**Where the solution lives after `Adapt`.** On a GPU build the transferred solution is left on the
**device** — at any rank count since #172 — so `solution%interior` (the host mirror) is **stale**. This matches the
rest of the time loop, where the device is authoritative and a caller wanting host data calls
`UpdateHost()` first, as `Write_DGModel3D_t` does. Before the device transfer existed the mirror
happened to be fresh; do not rely on that. On CPU builds the two are the same storage.

### 6.4 Device-side transfer (implemented)

Until this landed, a single-rank GPU adaptation epoch moved the entire solution across the host
link twice - `UpdateHost` before the regrid, `UpdateDevice` after - and ran the triple
tensor-product interpolation on the CPU in between. Both steps are now overridden on the GPU
backend: staging is a device-to-device copy into a model-owned buffer, and the plan is applied by
`TransferSolution_3D_gpu` (`src/gpu/SELF_SolutionTransfer.cpp`).

**Kernel.** One workgroup per new element with `(N+1)²` threads, where thread `(a,b)` owns one
`(a,b)` pencil and loops the free index through each of the three directional passes - the same
decomposition the 3-D modal indicator uses. An `(N+1)³` thread block is not an option: it is 4096
threads at `N = 15`, past the block limit, and it would push the working buffers into per-thread
scratch. The three `(N+1)³` working buffers are static `__shared__` arrays sized `AMR3D_MAXNP³`
with `AMR3D_MAXNP = 12`, i.e. 41 kB reserved per block at any degree; the Fortran caller guards
`N+1 <= 12` so a larger degree fails loudly rather than overrunning.

A dynamic-shared variant sized to the degree actually in use (12 kB at `N = 7`, the idiom
`SELF_MatrixMultiply.cpp` uses) was implemented and measured, and was **1-3% slower** on a B300
across three runs at two refinement depths - the SM shared budget there is large enough that 41 kB
per block does not gate occupancy, so the launch-time sizing was pure overhead. It was therefore
not adopted. The argument for it has not gone away on a 64 kB-LDS AMD device, where the static
footprint admits one workgroup per CU rather than five; that case has not been measured, and since
the alternative on every platform is the host round trip this kernel replaces, the static form
cannot regress anything as it stands.

**One deliberate divergence from the host.** The host descent calls `ProlongToChildren`, which
forms all eight children and discards seven; the kernel applies only the operator triple of the
child actually on the recorded path - an eighth of the work per descent step. What is dropped is
the seven discarded children, which are independent computations writing disjoint slices, and not
any term of the retained value: each contraction sums the same products against the same mortar
column in the same ascending index order as the host loop, so the reduction order is preserved.
Device and host nonetheless agree to round-off rather than bitwise, because the device compiler
contracts these multiply-accumulates into FMAs.

`test/solution_transfer_3d_device.f90` pins the result value by value against the portable
`ApplyTransferPlanRange`, on an epoch built to contain all three transfer kinds and prolongation
depths 0-2, with a non-polynomial field that differs in every variable. That is the check the AMR
regressions cannot make: conservation and entropy non-growth are invariants a kernel with a
transposed direction could still satisfy.

**Measurements.** One B300 (sm_103, CUDA 13, fp64), `linear_euler3d_amr_spherical_soundwave`,
20 epochs, mean of three runs, before and after differing **only** in the three files that
implement the transfer. `tAdapt` is reported per adapting epoch, because epochs whose leaf set is
unchanged cost only an indicator pass and mixing the two makes the metric depend on how many
epochs happened to adapt:

| maxLevel | elements | adapting epochs | before | after | |
| --- | --- | --- | --- | --- | --- |
| 1 | 512 | 1 | 119.5 ms | **101.6 ms** | -15.0% |
| 2 | 4,096 | 10 | 571.6 ms | **516.3 ms** | -9.7% |
| 3 | 20,784 | 20 | 2182.6 ms | **1959.6 ms** | -10.2% |

Time integration is unchanged, as it must be - nothing on this path runs inside `ForwardStep`.
`tForwardStep` over the same runs was 0.672 -> 0.615 s, 2.152 -> 2.116 s and 5.199 -> 5.176 s,
all inside the run-to-run spread (which is under 1% on `tAdapt`).

**Where the saving actually comes from, which is not where you would guess.** It is tempting to
credit the two eliminated full-field host/device copies, and the 2-D Stage 6a note above leads
with them. Do the arithmetic for the maxLevel-3 case: 20,784 elements x 125 nodes x 6 variables x
8 B is a 125 MB field, so the round trip is 250 MB, which at PCIe Gen5 rates is roughly 5 ms -
**0.2%** of a 2183 ms adaptation. The copies are not the story at this scale.

The 223 ms actually saved is the **host-side interpolation**: the portable
`ApplyTransferPlanRange` runs the triple tensor-product contractions on one CPU core for every
new element, and that work moves onto the GPU. The kernel also does an eighth of the descent work
the host does, because it prolongs onto the child on the recorded path instead of forming all
eight and discarding seven. Both effects scale with element count, which is why the saving holds
at roughly 10% of adaptation across a 40x range in mesh size while the copy time would have
faded to nothing.

The remaining ~90% of an adaptation is geometry generation, regrid and the indicator - untouched
by this change, and where the next 3-D optimization work belongs.

### 6.5 What is not done in 3-D

- **A coarsening regression.** 2-D has `lineareuler2d_amr_coarsen_wake`, which pins the amplitude
  gate's ability to release the mesh behind a passing front (§2.1.1). There is no 3-D analogue, so
  3-D coarsening is exercised only incidentally, by the soundwave regression and by the transfer
  tests' hand-built epochs.
- **A host-staged fallback for the migration.** The multi-rank device transfer (#172, §6.7) posts
  its sends and receives on device pointers, so it requires a GPU-aware MPI and there is no
  host-staged path behind it. That is not new — the per-step aggregated halo exchange has always
  done the same, with no fallback either, and `SELF_REQUIRE_GPU_AWARE_MPI` makes the configure-time
  check fatal by default — but it does mean a site without it has no degraded mode to fall back to,
  only `SELF_AMR_TRANSFER_HOST=1`, which gives up the device apply as well as the device receive.
- **A cheaper `SELF_AMR_MIGRATE_VERIFY` on the device path.** The diagnostic still downloads the
  whole window to run its bitwise comparison on the host, which is exactly the traffic §6.7 exists
  to remove. It is off by default and out of the time loop, so this is a cost only under the switch.
- **A distributed forest.** The forest and the transfer plan are still replicated on every rank,
  so forest memory is O(global elements) per rank and the per-epoch flag allgather remains. This is
  issue #167's Stage 2 (the p4est-style design), deliberately deferred until measurement justifies
  the much larger change: 2:1 balancing becomes a neighbourhood iteration, and the routing argument
  that makes migration free of handshakes depends on the plan being replicated.
- **Weighted partitioning.** The repartition is equal element count per rank, which is the right
  balance measure only because every element carries the same degree.
- **A 3-D visualization script.** 2-D ships `examples/linear_euler2d_amr_plot.py`. The 3-D
  example writes the same self-describing HDF5 snapshots (solution + geometry per file, so the
  changing mesh needs no special handling downstream) but no renderer is provided.

### 6.6 Multi-rank migration in 3-D

The mechanism is dimension-independent (§ Stage 5, 5d); what differs in 3-D is the price of getting
it wrong. An element's payload is `(N+1)^3 * nvar` reals rather than `(N+1)^2 * nvar`, so a full
field is large: at `maxLevel = 3` the example reaches 20,784 elements, which at 125 nodes and 6
variables in double precision is **125 MB**. Under the v1 gather every rank received every element
it did not own, once per adapting epoch, and allocated a second full field alongside the live one.
Under v2 a rank receives only the old elements its new range reads from a peer, which for a
balanced repartition is the load-imbalance slack at the ends of its range.

The examples report this directly. `BENCH_migrateBytesRecv`, `BENCH_migrateBytesSent` and
`BENCH_migrateElemRemote` are per-rank run totals, and both migration paths count them, so one
binary measures both sides:

```shell
# after (default) vs before (v1), same binary, same mesh trajectory
mpirun -n 4 ./examples/linear_euler3d_amr_spherical_soundwave
SELF_AMR_MIGRATE_GATHER=1 mpirun -n 4 ./examples/linear_euler3d_amr_spherical_soundwave
```

Compare `BENCH_tAdapt_s / BENCH_nAdaptEpochs` - the cost of an *adapting* epoch - rather than the
whole-epoch figure, which is diluted by indicator-only no-op epochs; and treat the comparison as
void unless `BENCH_nAdaptEpochs` and `BENCH_nElemFinal` agree between the two runs, since
otherwise the two took different mesh trajectories. `BENCH_tForwardStep_s` is the control: nothing
was added to the time-stepping loop.

#### What was measured

**Communication volume**, from the 3-D spherical soundwave and the 2-D ultrasound point source at
`maxLevel = 2` / `maxLevel = 4`, six epochs (five adapting), on a 12-thread CPU box. Per-rank
totals over the run; both runs of a pair have identical `nAdaptEpochs` and `nElemFinal`, and the
final pressure extrema agree to every printed digit. Note this is a *smaller* mesh than the 125 MB
example above, which is an illustration of the full-field size at a deeper level, not a measurement
of these runs.

| case | ranks | received per rank, v1 | received per rank, v2 | remote elements v1 -> v2 |
|---|---|---|---|---|
| 3-D | 2 | 35.9 MB | **0** | 5988 -> 0 |
| 3-D | 4 | 53.9 MB | **0.95 MB** | 8982 -> 159 |
| 3-D | 8 | 62.9 MB | **1.80 MB** | 10479 -> 300 |
| 2-D | 2 | 10.9 MB | **0** | 4242 -> 0 |
| 2-D | 4 | 16.3 MB | **0.55 MB** | 6363 -> 215 |
| 2-D | 8 | 19.0 MB | **0.28 MB** | 7420 -> 108 |

These count **bytes arriving in a rank's own memory** - the quantity that sets the migration
footprint - not total network traffic: what a collective moves internally depends on the algorithm
the MPI implementation selects, which is not something to model from the outside. Read two things
from the table. First, v1's received volume grows with rank count in both dimensions (each rank
approaches "everyone else's whole field"), while v2's is set by how far the partition boundary
actually moves. Second, at 2 ranks v2's traffic is exactly **zero**: a contiguous
space-filling-curve repartition of a symmetric refinement leaves each rank's new range sourced
entirely from its own old range - the windows printed under `SELF_AMR_MIGRATE_VERIFY=1` are
`[1,256]` and `[257,512]`, precisely the two old ranges. That is a property of the refinement
pattern, not of the mechanism: `test/amr_migrate_{2,3}d_window_mpi.f90` refines asymmetrically on
purpose and asserts that elements do cross ranks (17 of them at 2 ranks, and at 4 ranks one rank's
window is entirely remote while another's spans three peers), so a regression that silently made
the exchange a no-op would fail there.

Peak RSS is **not** a usable signal at this problem size: the deltas (16 MB at 2 ranks, 1 MB at 4,
7 MB at 8) are smaller than and not proportional to the buffer removed, because the footprint is
dominated by the replicated forest, the geometry double-buffers and the halo tables. The buffer
that went away matters at scale, where it is O(global field) per rank; here it is not what sets the
peak.

**Adaptation time**, on one B300 node with four GPUs (CUDA sm_103, one rank per GPU), `maxLevel = 3`,
ten adapting epochs, four repetitions per configuration with the two migration paths interleaved
per repetition so warm-up cannot land preferentially on one of them. Median seconds per adapting
epoch, with the observed range:

| ranks | v1 gather | v2 point-to-point | received per rank, v1 -> v2 |
|---|---|---|---|
| 1 | 0.982 [0.914, 0.985] | 0.965 [0.961, 0.988] | 0 -> 0 (nothing migrates) |
| 2 | 0.538 [0.525, 0.558] | 0.553 [0.528, 0.558] | 123.5 MB -> **0** |
| 4 | 0.332 [0.296, 0.340] | 0.287 [0.286, 0.302] | 185.2 MB -> **1.5 MB** |

At one rank the two are the same code path, and they measure the same - which is the check that
the windowed apply and its bounds guard cost nothing on the path that does not migrate. At two
ranks eliminating 123 MB per rank of host collective buys no measurable time: an allgather within
one node is fast, and adaptation is dominated by geometry and the indicator. At four ranks the
ranges no longer overlap and the point-to-point path is about 13% faster per adapting epoch. The
volume reduction is the durable result; the time saved is what that volume happened to cost on this
machine at this scale.

`BENCH_tForwardStep_s` is unchanged between the two paths at every rank count, as it must be.
(It does grow steeply with rank count on this node - 1.2 s, 4.7 s, 23.8 s - because the halo
exchange runs through host memory in this container configuration. That is a pre-existing property
of multi-rank GPU runs, identical in both paths, and unrelated to migration.)

### 6.7 The device-resident multi-rank transfer

§6.4 moved the transfer onto the device on ONE rank. On several it stayed on the host for a reason
that had nothing to do with the kernel: `ApplyTransferPlan_DGModel{2,3}D` fell back to the portable
path whenever a caller supplied host old-field data, and after #167 a multi-rank caller always did,
because the migrated window was received into host memory. So multi-rank GPU adaptation ran the
tensor-product interpolation on one CPU core and pushed the field across the host link twice per
adapting epoch. #172 closes that, in both dimensions.

**The kernel contract.** `TransferSolution_{2,3}D_gpu` takes an `oldFirst0` argument, and `nOld`
now means `uOld`'s per-variable element STRIDE rather than the global old-leaf count. A plan's
global source index is rebased as `src - oldFirst0` at the four places the kernel reads `uOld`.
The whole-field case is exactly `oldFirst0 = 0`, which is what the single-rank path passes, so its
indexing is unchanged arithmetic and the §6.4 measurements still describe it. This is stated as a
contract because it is what the tests pin, not an implementation detail.

**Where the window lives.** It became model state (`MigrateOldWindow` / `DownloadOldWindow` on
`DGModel{2,3}D_t`) rather than controller state, which is what lets each backend supply it its own
way while the controller stays portable. The portable implementation migrates into a flat host
buffer; the GPU override assembles the window in device memory — its own run by a device-to-device
copy out of the live pre-regrid solution, the peers' runs `MPI_Irecv`'d straight into device memory.
The schedule itself is shared: it moved into `SELF_SolutionMigration`, on flat
`(perElem, nElem, nvar)` buffers, which is also why one copy now serves both dimensions where there
were two.

**Device pointers into MPI is not a new dependency.** `SideExchangeStart` has always handed
`halo_sendbuf_gpu`/`halo_recvbuf_gpu` to `MPI_SEND_INIT`/`MPI_RECV_INIT`, with no host-staged
fallback anywhere, and `SELF_REQUIRE_GPU_AWARE_MPI` is a fatal configure-time check by default. A
multi-rank GPU run without a GPU-aware MPI could never have taken a time step, so the migration
requires nothing the time loop did not already require.

**The synchronization, which is easy to get wrong.** The local run is copied by a C++ launcher,
`MigrateWindowLocal_gpu`, not by a loop of Fortran `hipMemcpy` calls, because it also has to
synchronize the device before the caller posts MPI against `solution%interior_gpu` — and it does so
UNCONDITIONALLY, outside the "is there anything to copy" guard. Letting the copies carry the
synchronization implicitly fails twice: a rank whose local run is empty, every old element it owned
claimed by peers, would post MPI against unsynchronized device state; and it would rest on an
unstated assumption that every kernel writing the solution shares the copies' stream.
`HaloPack_{2,3}D_gpu` calls `hipDeviceSynchronize` before its own MPI for exactly this reason, and
that is the precedent followed here.

**What the kernel cannot check, the caller does.** On the host, `ApplyTransferPlanWindow` verifies
that every element's sources lie inside the window and stops if not. A kernel cannot: it has no way
to report, and an out-of-window index is a silent out-of-bounds device read rather than a visibly
wrong answer. So `ApplyTransferPlan` scans the plan over its own element range before launching.
Integer comparisons, once per adapting epoch.

**Two verification bars, kept apart deliberately.** `SELF_AMR_MIGRATE_VERIFY` compares the migrated
window against an allgathered reference BIT FOR BIT, and `VerifyMigration` is unchanged — on the
device path the window is downloaded and the same exact check runs. Everything that touched those
values is data movement (a device-to-device copy, MPI byte transfers, the download), so no
arithmetic has been applied and exactness is available, therefore it is demanded. The transfer
APPLY is a different question: it agrees between host and device only to round-off, because the
device compiler contracts the multiply-accumulates into FMAs and the kernel's descent prolongs only
onto the child on the recorded path. That gets its own switch, `SELF_AMR_TRANSFER_VERIFY`, and its
own tolerance. Merging the two would force the strict one down to the loose bar, and naming them
separately is the defence against someone later doing precisely that.

**`SELF_AMR_TRANSFER_HOST=1`** routes a GPU build back through the portable windowed path. It is
both an escape hatch and the A/B switch this path is timed with, so that before and
after share one binary, one build and one mesh trajectory — the condition an earlier #165 baseline
violated by changing the example between runs.

**Coverage, in two layers that fail for different reasons.** `test/amr_migrate_{2,3}d_device_apply_mpi`
drives the whole path through the model at 2 and 4 ranks against an independently built whole-field
reference. It uses ASYMMETRIC refinement, and must: under a symmetric refinement the contiguous
space-filling-curve repartition leaves every rank's new range sourced from its own old range, so not
a single byte moves and only the local copy runs — the AMR soundwave MPI regressions pass without
ever exercising a message. The test also asserts its own premise, that some rank's window begins
past 1 AND away from that rank's own first old element, because a window beginning exactly where the
rank's old range begins cannot distinguish a correct `oldFirst0` from a dropped one.

The offset is ALSO covered serially, at a lower level, by
`test/solution_transfer_{2,3}d_device_offset`, which calls `TransferSolution_{2,3}D_gpu` directly
with a synthetic window and a nonzero `oldFirst0`. Going under the model is what makes that
possible: through the model the offset is inherently multi-rank, because a rank's window base is
nonzero only when some other rank owned the elements it reads, and at one rank the hull of all new
elements is always `[1, nOld]`. The launcher takes `oldFirst0` as an ordinary argument, so it is not
bound by the model's whole-range contract. The value of having both is separation — the direct test
fails only for a kernel indexing error, the MPI test only for a migration or routing error, and each
asserts its own premise (`wFirst > 1`, `nWinElem < nOld`, `nWinElem >= 2`) so it cannot quietly stop
testing what it exists for. The direct tests are the one GPU-gated pair in this change, and they
print SKIP on a CPU build, because there the kernel is absent and there is nothing to compare. They
are not thereby untested: Buildkite runs on-hardware coverage pipelines for MI210 (HIP) and V100
(CUDA) on every pull-request branch, so both device backends build and run the full suite per PR.
GitHub Actions is the CPU-only half of the estate. Note that every job in both halves uses
gfortran - see `AMR3D-Design.md` §6 - so CI says nothing about ifx, nvfortran or amdflang.

On the host side `transfer_plan_{2,3}d_window` already pins the offset arithmetic bitwise over five
partition tables, including offset, all-remote and empty windows. And
`test/solution_transfer_2d_device` closes a separate pre-existing gap: the 2-D kernel had no
value-level test at all, only 3-D did.


**Measurements.** B300 (CUDA sm_103), one dedicated node, 4 GPUs, 4 repetitions with the two modes
interleaved WITHIN each repetition. `SELF_AMR_TRANSFER_HOST=1` selects the host side, so both sides
are one binary, one build and one mesh trajectory. Figures are seconds per ADAPTING epoch
(`BENCH_tAdapt_s / BENCH_nAdaptEpochs`), taken as the max over ranks — an adaptation is collective,
so its slowest rank gates it — then the median over repetitions. The percentage is the median of the
PAIRED per-repetition ratios, not the ratio of the medians: interleaving exists to make each
repetition a matched pair, and comparing medians throws that away. On 2-D at 4 ranks the two
differ by 13 points, which is how much the pairing is worth.

3-D, `linear_euler3d_amr_spherical_soundwave` at `maxLevel = 3`, 20 epochs, 20,784 elements — the
same configuration §6.4 measured, so the rows are directly comparable to it:

| ranks | elem/rank | host | device | paired median | per rep | recv/rank |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 20,784 | 2.006 s | 2.003 s | **-0.5%** | +0/-1/-1/-0% | 0 |
| 2 | 10,392 | 1.112 s | 1.003 s | **-10.3%** | -12/-11/-10/-9% | 0 |
| 4 | 5,196 | 0.601 s | 0.537 s | **-10.3%** | -10/-10/-10/-11% | 8.9 MB |

**One rank is a control, not a datapoint.** The branch is gated on `nRanks > 1`, so the two labels
execute identical code there and must agree; -0.5% is the harness reporting that it is measuring
what it claims to. A #167 sweep that looped the mode OUTER and repetitions inner put the whole
warm-up penalty on whichever mode ran first and manufactured a 1.7x difference at one rank, which is
why the control is kept and why the loop order is what it is.

**The 2-rank row is the interesting one.** It migrates ZERO bytes — under a symmetric refinement the
contiguous space-filling-curve repartition leaves every rank's new range sourced from its own old
range — and still gains 10.3%. So the saving cannot be the eliminated host copies; it is the
per-element tensor-product interpolation leaving the CPU, exactly as §6.4 concluded for one rank.
The 4-rank row does move 8.9 MB per rank, directly into device memory, and gains the same 10.3%.
The figure tracks rank-local element count rather than communication volume.

`BENCH_tForwardStep_s` is the control on everything else, and its paired drift is within ±0.3% at
every rank count. Watch the per-repetition values rather than a magnitude, though: at 2 ranks they
run -9/+36/+1/-1%, which CHANGES SIGN and is therefore run-to-run variance in the halo exchange
(which goes through host memory in this container configuration, per §6.6) rather than an effect of
the switch. Reported as magnitudes it would have looked like a consistent 12% regression.

2-D, `linear_euler2d_amr_ultrasound_pointsource` at `maxLevel = 6`, 20 epochs, 1,648 elements. The
gains are larger in relative terms and **should not be quoted precisely**: a 2-D adaptation here is
15-25 ms, small enough that the per-repetition spread swamps the effect. The direction is
consistent with 3-D; the magnitude is not measurable on this node at this scale.

| ranks | elem/rank | host | device | paired median | per rep |
| --- | --- | --- | --- | --- | --- |
| 1 | 1,648 | 0.0231 s | 0.0231 s | **+1.0%** (control) | -1/+1/+1/+3% |
| 2 | 824 | 0.0212 s | 0.0148 s | -28% | -34/-50/-23/**+9**% |
| 4 | 412 | 0.0146 s | 0.0117 s | -15% | -2/-16/-15/-37% |

The `+9%` repetition at 2 ranks is left in the table deliberately. A summary statistic that hid it
would imply a reliability the measurement does not have at this problem size.

Validity gates, all of which held: `BENCH_nAdaptEpochs`, `BENCH_nElemFinal` and `BENCH_nSteps`
identical between modes at every rank count, and `BENCH_migrateBytesRecv` identical too - this
change moves where bytes land, not how many, so a difference there would have meant the routing had
changed and voided the comparison.

---

## 7. References

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
