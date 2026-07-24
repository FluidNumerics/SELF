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
| `h`-refinement mesh mutation (element subdivision / coarsening) | Designed (Stage 2) |
| Solution transfer (prolongation / restriction) | Designed (Stage 3) |
| Mortar regeneration + 2:1 balancing | Designed (Stage 4) |
| MPI dynamic re-partitioning / load balancing | Designed (Stage 5) |
| GPU device re-allocation for a changing element count | Designed (Stage 6) |

Only the trigger is wired into the library today. Everything from Stage 2 onward changes the
statically-allocated mesh model and is intentionally deferred behind this design so it can be
landed and reviewed in self-contained pieces.

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

- Add an explicit **quadtree / forest-of-quadtrees** parent–child structure alongside the flat
  element arrays: each leaf is an active element; refinement replaces one leaf with four
  children (each spanning a reference sub-quadrant), coarsening merges four siblings back to
  their parent.
- Grow the mesh arrays via an amortized-capacity scheme (allocate spare capacity, reallocate in
  chunks) so refinement does not reallocate every step. Keep the flat solver arrays contiguous
  by compacting active leaves.
- Generate child `nodeCoords` from the parent geometry map (isoparametric subdivision), and
  recompute metric terms/Jacobians for children via the existing geometry routines — **no
  change to the geometry algorithms**, only new inputs.

### Stage 3 — Solution transfer (prolongation / restriction)

- **Prolongation** (parent → 4 children): interpolate the parent nodal solution onto each
  child's nodes. Because a child occupies a reference sub-quadrant, this is exactly the
  one-sided mortar restriction operator `Lagrange%mortarR` applied in each direction — the
  operator is already built and tested.
- **Restriction** (4 children → parent): the \(L^2\)-projection adjoint `Lagrange%mortarP`,
  again already available, applied per direction. Restriction is conservative by construction
  (the mortar operators satisfy \(\sum_k P_k R_k = I\) and discrete conservation), so coarsening
  preserves cell-integrated quantities.
- This stage must preserve **conservation** and **not reorder floating-point reductions** in the
  solver; the transfer is a separate operator applied between time steps.

### Stage 4 — Mortar regeneration and 2:1 balancing

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
