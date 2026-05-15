# Off-Grid Point Sampling

SELF stores fields on element-local, per-quadrature-node arrays
(`MappedScalar2D%interior(i,j,iEl,iVar)` /
`MappedScalar3D%interior(i,j,k,iEl,iVar)`). For probes, observation
operators, particle trackers, or post-processing along arbitrary tracks
you often need to sample a field at a user-supplied physical
coordinate that does **not** sit on a quadrature node. The
`SELF_Points` module provides this primitive.

A `Points` instance holds a cloud of physical coordinates and offers
two operations:

1. **`LocatePoints`** — find which (rank-local) element contains each
   point and recover its reference coordinates
   $(\xi_1, \xi_2[, \xi_3]) \in [-1, 1]^d$ via a spatial-hash candidate
   pre-filter and a Newton inverse-map.
2. **`EvaluateScalar`** — sample a `MappedScalar2D` or `MappedScalar3D`
   at the located points via tensor-product Lagrange interpolation,
   returning an array shaped `(nPoints, nVar)`.

## Data layout

```fortran
type, public :: Points_t
  integer :: nPoints           ! number of points in the cloud
  integer :: nDim              ! 2 or 3
  real(prec), allocatable :: x(:,:)            ! (nPoints, nDim) physical coords (input)
  integer,    allocatable :: elements(:)       ! (nPoints) rank-local element id; 0 = not found
  real(prec), allocatable :: coordinates(:,:)  ! (nPoints, nDim) reference coords (s,t[,u])
end type
```

After a successful `LocatePoints` call, the module also caches per-point
Lagrange basis values (one row of `lS_cache(0:N, p)` per point per
reference direction) so that subsequent `EvaluateScalar` calls skip
basis evaluation entirely; see [Cached basis](#cached-basis) below.

## Algorithm

### Spatial hash

For each element $e$, an axis-aligned bounding box is computed over the
control-grid nodes `geometry%x%interior(:,:[,:],e,1,1:d)`. A uniform
grid of cells covers the global bounding box, sized so each cell
contains $\mathcal{O}(1)$ elements on average. Every element is
inserted into all cells its bounding box overlaps, producing a
CSR-style table

```text
cellStart(0:nCells)         ! prefix-summed offsets
cellElems(1:totalEntries)   ! element ids per cell
```

For each query point we read off the owning cell, walk its candidate
list, apply a cheap bounding-box reject, and run the Newton solver on
the survivors.

### Newton inverse-map

For element $e$ with the high-order isoparametric map
$X^e(\boldsymbol{\xi})$, we solve $X^e(\boldsymbol{\xi}) =
\mathbf{x}_{\text{target}}$ for $\boldsymbol{\xi}$.

At each iteration $k$:

1. Evaluate the 1D Lagrange basis at each component of $\boldsymbol{\xi}^{(k)}$
   via `Lagrange%CalculateLagrangePolynomials` (barycentric form,
   Kopriva 2009, Algorithm 34).
2. Tensor-product evaluate the current physical coordinate

    \begin{equation}
      X^e(\boldsymbol{\xi}^{(k)}) = \sum_{i,j[,k]}
        \ell_i(\xi_1)\,\ell_j(\xi_2)\,[\ell_k(\xi_3)]\;
        \mathbf{x}^e_{i,j[,k]} .
    \end{equation}

3. Tensor-product evaluate the Jacobian $J = \partial X^e / \partial \boldsymbol{\xi}$
   from the precomputed covariant basis tensor
   `geometry%dxds%interior(:,:[,:],e,1,r,c)` — index convention
   $J_{r,c} = \partial x_r / \partial \xi_c$.
4. Solve the $d \times d$ linear system $J \cdot \delta\boldsymbol{\xi}
   = \mathbf{x}_{\text{target}} - X^e(\boldsymbol{\xi}^{(k)})$ in
   closed form ($2{\times}2$ Cramer / $3{\times}3$ cofactor).
5. Update $\boldsymbol{\xi}^{(k+1)} = \boldsymbol{\xi}^{(k)} +
   \delta\boldsymbol{\xi}$.

The iteration is declared converged when $\|\delta\boldsymbol{\xi}\|_\infty
< $ `newtonTolerance` (1e-8). It is rejected when the converged
$\boldsymbol{\xi}$ lies outside $[-1, 1]^d$ (with a small slack) or
when `newtonMax` (500) iterations are exhausted. Rejected candidates
fall through to the next element in the spatial-hash cell.

### Cached basis

Because the converged $\boldsymbol{\xi}$ is a function of the
reference coordinates only, the per-point Lagrange basis values
$\ell_i(\xi_1), \ell_j(\xi_2), [\ell_k(\xi_3)]$ are computed **once**
at the end of `LocatePoints` and stored on the `Points` object:

```fortran
lS_cache(0:N, 1:nPoints)
lT_cache(0:N, 1:nPoints)
lU_cache(0:N, 1:nPoints)   ! 3D only
```

`EvaluateScalar` checks `nCached == scalar%interp%N` and, if so, uses
the cached basis directly — for probes sampled every timestep this
removes the dominant per-call cost. If the cache is invalid
(different polynomial degree, or no `LocatePoints` has run),
`EvaluateScalar` falls back transparently to on-the-fly basis
evaluation.

### Sampling

For each located point $p$, the sampled value of variable $v$ is

\begin{equation}
  s_v(p) = \sum_{i,j[,k]}
    \ell_i(\xi_1)\,\ell_j(\xi_2)\,[\ell_k(\xi_3)]\;
    \text{scalar.interior}(i,j[,k],e(p),v) .
\end{equation}

Points with `elements(p) == 0` (off-mesh or owned by another rank)
sample to zero.

## Rank-local search

Each MPI rank owns a subset of elements
(`mesh%decomp%elemToRank`). `LocatePoints` runs only over the
elements the rank's `geometry` object holds — it does **not** perform
any cross-rank communication. Points that lie inside a non-local
element produce `elements(p) = 0` and `EvaluateScalar` returns zero
for them, exactly as for points outside the mesh.

This rank-local design matches the existing "rank-local data
ownership" rule used throughout SELF. If a global point-to-element
map is needed, callers can post-process by collecting the local
`elements` arrays across ranks (e.g. with `MPI_Allreduce` choosing
the non-zero value).

## Usage

### Basic CPU example (2D)

```fortran
use SELF_Constants
use SELF_Lagrange
use SELF_Mesh_2D
use SELF_Geometry_2D
use SELF_MappedScalar_2D
use SELF_Points

type(Lagrange), target :: interp
type(Mesh2D),   target :: mesh
type(SEMQuad),  target :: geometry
type(MappedScalar2D)   :: f
type(Points)           :: probes

integer, parameter :: nProbes = 100
real(prec) :: xProbes(nProbes, 2)
real(prec) :: values(nProbes, 1)

call interp%Init(N=7, controlNodeType=GAUSS, M=16, targetNodeType=UNIFORM)
call mesh%Read_HOPr("mymesh.h5")
call geometry%Init(interp, mesh%nElem)
call geometry%GenerateFromMesh(mesh)

call f%Init(interp, nvar=1, nelem=mesh%nElem)
call f%AssociateGeometry(geometry)
call f%SetEquation(1, 'f = sin(x)*cos(y)')
call f%SetInteriorFromEquation(geometry, 0.0_prec)

! Fill xProbes(:,1) and xProbes(:,2) with the physical coordinates
! you want to sample at...

call probes%Init(nProbes, 2)
call probes%SetPoints(xProbes)
call probes%LocatePoints(geometry)   ! one-shot setup; fills cache

call probes%EvaluateScalar(f, values)  ! cheap; reuses cache
```

`LocatePoints` is the expensive step; call it once per cloud.
`EvaluateScalar` is cheap and is safe to call every timestep.

### 3D and multi-variable

3D usage is identical — call `probes%Init(nPoints, 3)`, pass a
`type(SEMHex)` geometry, and use a `MappedScalar3D`. For
multi-variable fields, `EvaluateScalar` writes a column per variable:

```fortran
real(prec) :: values(nProbes, nvar)
call probes%EvaluateScalar(f, values)
```

### GPU backend

When SELF is built with `SELF_ENABLE_GPU=ON`, `Points` is a
device-resident object: it carries `elements_gpu`, `coordinates_gpu`,
and `lS_cache_gpu` / `lT_cache_gpu` / `lU_cache_gpu` pointers in
addition to the host arrays. `LocatePoints` still runs on the host
(the point search is irregular and one-shot) and pushes the result
to the device via `UpdateDevice`. **`EvaluateScalar` on GPU consumes
only the cached basis** — there is no Lagrange-polynomial evaluation
in the kernel, only the tensor-product contraction.

The GPU `EvaluateScalar` is dispatched off the same generic name as
the CPU one. Pass a caller-allocated device pointer instead of a
Fortran array:

```fortran
use SELF_GPU

type(c_ptr)        :: values_dev
integer(c_size_t)  :: nBytes

nBytes = int(nProbes*nvar, c_size_t) * int(prec, c_size_t)
call gpuCheck(hipMalloc(values_dev, nBytes))

call probes%EvaluateScalar(f, values_dev)   ! kernel launch; result on device

! ... use values_dev directly with other device kernels, or copy back ...
call gpuCheck(hipMemcpy(c_loc(values), values_dev, nBytes, hipMemcpyDeviceToHost))
call gpuCheck(hipFree(values_dev))
```

The same call site `probes%EvaluateScalar(f, ...)` resolves to either
the host-output specific (when the second argument is a Fortran
array) or the device-output specific (when it is a `type(c_ptr)`).

## Off-mesh points

Points whose physical coordinate lies outside every element this rank
owns receive the sentinel `elements(p) = 0`. `coordinates(p,:)` is
undefined for such points; `EvaluateScalar` returns zero for them.
This is the same path that's taken when a point lives in a non-local
element under MPI.

## Limitations and future work

- The point search is rank-local; a global gather-scatter must be
  handled by the caller.
- The GPU path runs `EvaluateScalar` only — `LocatePoints` is still
  host-side. The host search is one-shot per cloud, so this is rarely
  the bottleneck for the use case (sampling many timesteps with a
  fixed cloud).
- Sampling of `MappedVector2D/3D` and `MappedTensor2D/3D` is not yet
  implemented; the same machinery applies.
