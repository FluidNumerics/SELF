# Design: Adaptive Mesh Refinement (AMR) for 3D Models — Octrees + Face Mortars

**Status:** Implemented (Stages 1–6a) — mirrors the implemented 2D design.
Stages 1–5 landed with #164; the device-resident solution transfer (Stage 6a's
3D analogue) landed with #165. What is still outstanding from §6 is noted in
that table. Sections written in the future tense below have been updated where
the delivered code differs from the plan; where it does not, the plan text is
the record of what was built.
**Scope:** 3D hexahedral meshes, octree (2:1) h-refinement built on new 3D
2:1 face-mortar interface support, spectral (Legendre modal) refinement
indicator, CPU and GPU backends, MPI-parallel. Pilot model: `LinearEuler3D`.

Companion documents: [Mortar2D-Design.md](Mortar2D-Design.md) and
[AMR2D-Design.md](AMR2D-Design.md). The 2D stack (mortars → quadtree forest →
transfer plan → controller → MPI → GPU) is implemented and validated; this
document records only the 3D-specific decisions. Where a 3D module is a
mechanical transcription of its 2D counterpart, the 2D document remains the
authoritative rationale.

---

## 1. Structural mapping from 2D

| 2D (implemented) | 3D (this plan) |
|---|---|
| `Mesh2D_t%mortarInfo(1:8,:)` (1 big + 2 small edges) | `Mesh3D_t%mortarInfo(1:14,:)` (1 big + 4 small faces) |
| `MortarExchange` / `MortarFluxCollect` on 2D mapped data | same on `MappedScalar3D_t` / `MappedVector3D_t` |
| 1D `mortarR/mortarP` applied per edge | same 1D matrices applied as 2D tensor products per face |
| flip ∈ {0,1}; sub-edge pairing `t` vs `3-t` | flip ∈ 0..7 (dihedral group); quadrant permutation table |
| flux collect factor −2 (sub-edge Jacobian ½) | factor −4 (sub-face Jacobian ¼) |
| `SELF_QuadTreeMesh_2D` (4 children) | `SELF_OctreeMesh_3D` (8 children) |
| `childOfSide(1:2,1:4)` | `childOfFace(1:4,1:6)` |
| `RefineConnectivity` 4 inner face pairs | 12 inner face pairs |
| transfer: 2-pass tensor product, 4 children, ¼ | 3-pass tensor product, 8 children, ⅛ |
| indicator: 2D modal transform `Â = P U Pᵀ` | 3D modal transform (three passes), same criteria |
| `mortarBuff(N+1,1:4,m,v)` | `mortarBuff(N+1,N+1,1:8,m,v)` |
| `dt = dtBase/2**MaxLevel` | identical |

New modules (each mirrors the 2D file of the same name, 2D→3D):
`SELF_OctreeMesh_3D.f90`, `SELF_AdaptiveMesh_3D.f90`,
`SELF_MeshRefinement_3D.f90`, `SELF_RefinementPrimitives_3D.f90`,
`SELF_RefinementIndicator_3D_t.f90` (+ `cpu`/`gpu` backends),
`SELF_SolutionTransfer_3D.f90`, `SELF_TransferPlan_3D.f90`,
`SELF_AMRController_3D.f90`.

---

## 2. Conventions (normative)

### 2.1 Face trace coordinates

`boundary(i,j,side,elem)` face coordinates, from `BoundaryInterp`:

| sides | face (i,j) |
|---|---|
| 1 (Bottom), 6 (Top) | (ξ₁, ξ₂) |
| 2 (South), 4 (North) | (ξ₁, ξ₃) |
| 3 (East), 5 (West) | (ξ₂, ξ₃) |

### 2.2 Children

Children are ordered by CGNS corner-node order:
`axc = [0,1,1,0,0,1,1,0]`, `ayc = [0,0,1,1,0,0,1,1]`,
`azc = [0,0,0,0,1,1,1,1]`; child c covers `ξ_d ∈ [a-1, a]` per direction.
Subdivision is orientation-preserving: a child face lying on parent face `s`
is that child's local face `s`, so neighbor/flip data inherits verbatim
(same load-bearing fact as in 2D).

### 2.3 Face quadrants and `childOfFace`

A big face is split into 4 sub-faces ("quadrants") indexed in the **big
face's** trace coordinates: `q = kx + 2*(ky-1)` with `kx, ky ∈ {1,2}` the
half-interval index along the face's i and j coordinate (1 = [-1,0],
2 = [0,1]). Restriction to quadrant q is the tensor product
`R_kx ⊗ R_ky` of the existing 1D `mortarR`; projection back is
`P_kx ⊗ P_ky` of `mortarP`.

`childOfFace(q,s)` — the child of a refined hex whose face covers quadrant
`q` of parent face `s` (derived from §2.1 + §2.2):

| s | face | childOfFace(1:4,s) |
|---|---|---|
| 1 | Bottom | 1, 2, 4, 3 |
| 2 | South | 1, 2, 5, 6 |
| 3 | East | 2, 3, 6, 7 |
| 4 | North | 4, 3, 8, 7 |
| 5 | West | 1, 4, 5, 8 |
| 6 | Top | 5, 6, 8, 7 |

Octree helpers: `oc_opposite = [6,4,5,2,3,1]`;
`oc_reflect(s,c)`: s∈{1,6} → `[5,6,7,8,1,2,3,4]`,
s∈{2,4} → `[4,3,2,1,8,7,6,5]`, s∈{3,5} → `[2,1,4,3,6,5,8,7]`;
`oc_subpos(c,s)` is the inverse of `childOfFace`.

### 2.4 Flips and quadrant permutation

3D flips follow the existing 8-state convention
(`SELF_Mesh_3D_t.f90`): `(i2,j2) = F_f(i1,j1)` maps receiver-face indices to
donor-face indices. The same maps applied to half-interval indices give the
quadrant permutation `mortarQuadPerm(q,f)` — the donor-face quadrant that
coincides with receiver-face quadrant q:

| f | perm(1:4) | | f | perm(1:4) |
|---|---|---|---|---|
| 0 | 1,2,3,4 | | 4 | 1,3,2,4 |
| 1 | 2,1,4,3 | | 5 | 2,4,1,3 |
| 2 | 4,3,2,1 | | 6 | 4,2,3,1 |
| 3 | 3,4,1,2 | | 7 | 3,1,4,2 |

On octree-emitted meshes, mortars interior to a root tree always carry
flip 0; mortars across a root face inherit the base mesh's face flip. The
mortar exchange implements all 8 flips so unstructured (HOHQMesh/HOPr) base
meshes are supported.

### 2.5 `mortarInfo(1:14, 1:nMortars)`

Replicated on all ranks; element ids are global. Rows:

```
 1    big element id            2    big local face id
 3,4  small elem, 10*face+flip  (sub-face covering big-face quadrant 1)
 5,6  same for quadrant 2
 7,8  same for quadrant 3
 9,10 same for quadrant 4
11:14 global side ids of sub-faces 1..4 (MPI message tags)
```

As in 2D: mortar faces carry `sideInfo(3) = 0` and `sideInfo(5) = 0`
(conforming machinery and BC mapping skip them) and `sideInfo(1)` = mortar
index; the flip stored with each small face maps **big-face quadrant
coordinates to the small face's own coordinates** (the same receiver→donor
convention as `sideInfo(4)`).

### 2.6 Mortar data plane

`mortarBuff(1:N+1, 1:N+1, 1:8, 1:nMortars, 1:nvar [,1:3])` — slots 1–4 hold
the big-face trace (replicated per quadrant; MPI receives land per-quadrant),
slots 5–8 hold the small-face traces re-oriented into the big face's
coordinates. Exchange algorithm, per mortar:

1. Stage rank-local big trace into slots 1–4 and (flip-reoriented) small
   traces into slots 4+q; MPI posts fill remote slots (small ranks receive
   the big trace into slot q; the big rank receives small traces into 4+q),
   with received small traces flip-reoriented after the wait.
2. Small `extBoundary(i,j)` at face coords `F_f(i,j)` ← `(R_kx ⊗ R_ky)`
   applied to slot q (exact restriction).
3. Big `extBoundary` ← `Σ_q (P_kx ⊗ P_ky)` applied to slots 4+q
   (L2 projection; `mortarP` carries ½ per direction ⇒ ¼ total, the correct
   solution-trace projection).

`MortarFluxCollect` (on `boundaryNormal`): big-face integrand
`:= −4 Σ_q (P_kx ⊗ P_ky) g_q`. The ¼ from the tensor-product `mortarP`
combines with the ×4 integrand conversion so that the big face's discrete
surface integral equals minus the sum of the four small faces' integrals to
roundoff — the same conservation argument as 2D with ½/×2.

MPI messages use tags `globalSideId + nUniqueSides*(ivar-1)` per sub-face,
exactly the 2D scheme; the message pattern remains async point-to-point,
one round per stage.

### 2.7 Balance and hanging edges

2:1 balance is enforced across **faces only**, the direct analogue of 2D
(which balances across edges and tolerates 2-level corner jumps). DG face
mortars carry all interface data; hanging edges/corners carry none, so
face-2:1 balance is sufficient for well-posed connectivity. `EmitMesh`
guards `MaxLevelJump() <= 1` as in 2D.

### 2.8 Solution transfer

**Operators.** Prolongation: `U_child = (R_kx ⊗ R_ky ⊗ R_kz) U_parent` (exact),
`(kx,ky,kz) = (axc(c)+1, ayc(c)+1, azc(c)+1)`. Restriction:
`U_parent = Σ_c (P_kx ⊗ P_ky ⊗ P_kz) U_child(c)` — conservative, and
`Restrict(Prolong(u)) = u` to roundoff via `Σ_k P_k R_k = I` per direction.
Reference-cell conservation: `Σ w³ u_parent = (1/8) Σ_c Σ w³ u_child`. Each 1D
`P_k` carries the half-interval Jacobian, so the triple product carries 1/8.
Implemented in `SELF_SolutionTransfer_3D.f90` as three sequential directional
contractions; the loop order there is normative, since CLAUDE.md forbids
reordering floating-point reductions.

**Epoch protocol.** The plan (`SELF_TransferPlan_3D`) is built on the host after
the last forest mutation and classifies every new leaf as `COPY`, `PROLONG` or
`RESTRICT`, plus a `depth` and an octant `path` for any further descent. It is
applied around the regrid as a three-step, type-bound protocol on the model, so
the backend split that already selects `Regrid` selects the transfer too:

```fortran
call model%StageSolutionForTransfer()   ! preserve the field; Regrid may then free it
call model%Regrid(newMesh,newGeom)
call model%ApplyTransferPlan(plan,interp,eFirst,eLast)
```

`StageSolutionForTransfer` exists because `Regrid` reallocates (or resizes) the
storage the solution lives in, so the pre-regrid field must be preserved first.
The portable implementation stages into `DGModel3D_t%transferStage`, a host
array whose lifetime is exactly stage-to-apply — it is allocated by the first
call, consumed and deallocated by the second, and released in `Free`. Calling
`ApplyTransferPlan` without a preceding stage is a hard error (`stop 1`), pinned
by `test/dgmodel3d_guard_transfer_unstaged.f90`, because the alternative is
reading unallocated storage and failing somewhere less obvious.

`ApplyTransferPlan` takes an optional `uGlobal` argument: the **multi-rank
escape hatch**. The forest is rank-replicated and migration is gather-then-slice
(2D Stage-5 v1), so on more than one rank the controller allgathers the old
field on the host and passes it through `uGlobal`; each rank then fills exactly
its own new contiguous element range and elements that changed ranks are
migrated by construction. Because the allgather is inherently a host operation,
that branch stays on the portable host transfer on every backend.

**Device path.** On a single-rank GPU build both steps are overridden
(`src/gpu/SELF_DGModel3D.f90`): staging is a device-to-device copy into a
model-owned buffer, and the plan is applied by `TransferSolution_3D_gpu`, so an
adapting run moves no solution data across the host link at all. Consequence
worth knowing: the transferred solution is then left on the device and
`solution%interior` (the host mirror) is **stale** after `Adapt`. That matches
the rest of the time loop, where the device is authoritative, but it is a
behaviour change from the host-only implementation, where the mirror happened to
be fresh. See §4 for the kernel, and the Learning page for the measurements.

### 2.9 Indicator

3D tensor Legendre modal transform (three passes with the same `Pmodal`);
`E_tot`, `E_clip1`, `E_clip2` computed over the 3D mode cube with the same
top-mode clipping per direction; identical `S_e`/`σ_e` definition, amplitude
gate, relative energy floor (default `1e-12`), hysteresis band, and the
two-phase (device/host) split that keeps CPU and GPU flags identical.

### 2.10 EC / split-form models

`ECDGModel3D` raises a runtime error on meshes with `nMortars > 0`
(entropy-stable mortar operators out of scope), mirroring the 2D guard.

---

## 3. What carries over unchanged

- `Lagrange_t%mortarR/mortarP` (+ existing device mirrors) — no new operator
  math anywhere in the plan.
- `DomainDecomposition_t` — contiguous SFC partition of the Morton-ordered
  leaf list per epoch, `elemToRank` locality tests, tag conventions.
- Controller architecture: rank-replicated forest, flag allgather,
  level cap + halo passes, `AdaptFromFlags` (coarsen-then-refine, snapshot
  based), `Balance2to1` fixed-point sweeps, no-op detection,
  `BuildTransferPlan` walk over stable node ids, `EmitMesh` regeneration,
  geometry double-buffering with `CopyElements` reuse
  (`sourceKind == COPY` ⇒ same forest node ⇒ bit-identical coordinates),
  `Regrid` via `Resize` (not Free+Init), Stage-5 v1 gather-then-slice
  migration, `RecommendedTimeStep = dtBase/2**MaxLevel`, and the
  `SELF_AMR_GEOM_{FULL,VERIFY,NO_REUSE}` diagnostics.

New 3D-side prerequisites (absent today, required by the controller):
`Resize` on `Scalar3D`/`Vector3D`/`Tensor3D`/mapped variants,
`SEMHex%Resize/CopyElements/GenerateFromNodeCoords/UploadGeometry` + cached
scratch + `nElem = 0` default, `DGModel3D_t%Regrid/
StageSolutionForTransfer/ApplyTransferPlan` + `transferStage`.

Geometry reuse matters *more* in 3D: `CalculateContravariantBasis_SEMHex`
uses the curl-invariant metric form and is substantially more expensive than
its 2D counterpart.

---

## 4. GPU strategy

Same division of labor as 2D: numerics on device, adaptation logic on host.
Kernels follow the existing files' pattern: 3D mortar gather/scatter and
flux-scatter kernels (`SELF_Mortar.cpp`), a 3D indicator kernel
(`SELF_Refinement.cpp`, `AMR3D_MAXNP` bound), and the 3D transfer kernel
(`SELF_SolutionTransfer.cpp`). All are implemented.

The transfer kernel is `TransferSolution_3D_gpu`. As planned, the 2D kernel's
one-block-per-element with `Np²` threads and one thread per node becomes
one-block-per-element with `Np²` threads where thread `(a,b)` owns one `(a,b)`
pencil and loops the free index through each of the three directional passes — the
same decomposition `RefinementIndicator_3D_gpukernel` uses. A `Np³` thread block
is not an option: it is 4096 threads at N=15, past the block limit, and it would
push the working buffers into per-thread scratch.

Three `Np³` `__shared__` buffers do fit, at a bound: `AMR3D_MAXNP` is 12, giving
3·12³·8 B = 41,472 B, inside both the 48 KB CUDA static shared-memory limit and
the 64 KB gfx90a/gfx942 LDS budget. The allocation is static, so it is 41,472 B
at every degree, not `(N+1)³`; a dynamic-shared variant sized to the degree in
use measured 1-3% slower on a B300 and was not adopted (see the Learning page's
§6.4). The Fortran caller guards `interp%N+1 <= 12` so a larger degree fails
loudly rather than overrunning.

One deliberate divergence from the host reference, inherited from 2D: the host
descent calls `ProlongToChildren`, which forms all eight children and discards
seven, whereas the kernel applies only the operator triple of the child actually
on the recorded path — an eighth of the work per descent step. What is dropped is
the seven discarded children, not any part of the retained value: each contraction
sums the same terms against the same mortar column in the same ascending index
order as the host loop, so the reduction order CLAUDE.md protects is preserved.
Device and host nonetheless agree only to round-off rather than bitwise, because
the device compiler contracts these multiply-accumulates into FMAs.
`test/solution_transfer_3d_device.f90` pins the kernel's output against the host
reference value by value; the AMR regressions assert the invariants (conservation,
entropy non-growth), which conservation arguments make exact either way.

---

## 5. Testing & validation plan

Serial/CPU (mirroring the 2D suite):
- `mappedscalarmortarexchange_3d_linear`: linear field reproduced exactly
  across a hand-built 3D mortar mesh (all interior faces conforming or
  mortar), including a rotated-neighbor variant covering nonzero flips.
- `mappedvectordgdivergence_3d_mortar`: divergence of a linear flux on a
  mortar mesh matches the conforming result; conservation of the big-face
  surface integral to roundoff.
- `lineareuler3d_mortar_soundwave`: acoustic pulse crossing a static mortar
  interface; conservation + entropy checks (a conforming-mesh
  `lineareuler3d_soundwave` baseline lands first — LinearEuler3D currently
  has no dedicated regression test).
- `octree_balance_3d`, `mesh3d_uniform_refine`, `adaptive_mortar_3d`,
  `solution_transfer_3d`, `transfer_plan_3d`, `refinement_indicator_3d_*`:
  transcriptions of the 2D unit tests (validity invariants, transfer
  exactness/identity/conservation, indicator smooth/front + guards).
- `lineareuler3d_amr_soundwave`: the 2D soundwave AMR regression in 3D —
  adapt-to-convergence, conservation per epoch, dt halving, entropy
  non-increase, forest replication checksum.
- Guard tests (`WILL_FAIL`) for every error stop added, following the 2D
  guard-test inventory (controller decomp/maxlevel/wrongmesh, unbalanced
  emit, transfer-plan misuse, indicator argument guards, EC mortar guard).

MPI (2 ranks): `mappedscalarmortarexchange_3d_linear_mpi`,
`lineareuler3d_mortar_soundwave_mpi`, `lineareuler3d_amr_soundwave_mpi`
(partition-changing adaptation), `geometry_3d_reuse_mpi` + env-gated verify.

GPU: the same integration tests through the Buildkite MI210/V100 pipelines;
CPU/GPU flag and post-adapt solution agreement enforced by construction
(shared host `FinalizeIndicator`) and by tolerance respectively.

---

## 6. Work breakdown & phasing

| Stage | Content | Mirrors |
|---|---|---|
| 1 | 3D face mortars: `mortarInfo`, `SimpleMortarMesh3D`(+rotated), scalar/vector exchange + flux collect, model integration + EC guard, serial+MPI tests | PR #145 |
| 2 | Octree forest, refinement primitives, uniform refine, `EmitMesh`, validity tests | AMR Stages 2/4a/4b |
| 3 | Indicator 3D, transfer plan/solution transfer, `Resize`/geometry infrastructure, `Regrid`, controller, unit tests | AMR Stages 1/3/4 |
| 4 | `lineareuler3d_amr_soundwave` (+ coarsen-wake) integration tests, example | #156/#162 |
| 5 | MPI: forest-from-decomposed-mesh, flag allgather, gather-then-slice migration, 2-rank tests | Stage 5 |
| 6 | GPU: mortar/indicator/transfer kernels, device staging, geometry upload paths | Stage 6 |

Delivery status: stages 1–3 and 5 landed with #164, together with the
`lineareuler3d_amr_soundwave` serial and 2-rank integration tests of stage 4.
Stage 6's mortar and indicator kernels landed with #164; the device-resident
solution transfer and its staging landed with #165 (`TransferSolution_3D_gpu`,
`StageSolutionForTransfer_DGModel3D`, `ApplyTransferPlan_DGModel3D`), which also
delivered stage 4's example, `linear_euler3d_amr_spherical_soundwave`.

Still outstanding: the 3D `coarsen-wake` regression (the analogue of
`lineareuler2d_amr_coarsen_wake`), and the 3D analogues of 2D Stages 6b/6c
(amortized high-water-mark storage and geometry reuse) to the extent they are
not already inherited from the shared data classes.

Stages are landed as separate reviewable commits on one PR (or stacked PRs
at maintainer preference), each leaving the suite green.

---

## 7. Compliance with repository constraints (CLAUDE.md)

Identical posture to the 2D work, which established the precedents:
- The mortar interface discretization is the same formulation extension
  already approved and merged in 2D (#145), extended dimensionally.
- Conforming meshes execute byte-identical code paths (`nMortars > 0`
  gates); all existing 3D regression tests must pass bitwise.
- No collectives inside time stepping; adaptation runs between steps.
- No module renames/moves; additive public API only; Fortran 2008;
  no external dependencies (octree is home-grown); `fprettify` formatting.
- EC split-form 3D models excluded from nonconforming meshes (guarded).
