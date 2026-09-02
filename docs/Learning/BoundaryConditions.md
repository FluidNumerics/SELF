# Boundary Condition System

This page describes the internal design of SELF's extensible boundary condition (BC) system. It is intended for developers who want to understand the architecture, add new BC types to built-in models, or extend the system itself.

For a practical guide to using boundary conditions in your own applications, see [Boundary Conditions (User Guide)](../Models/boundary-conditions.md).

## Design Goals

The BC system is designed to:

- Allow each model to register only the BC types it supports
- Dispatch BC application through procedure pointers, avoiding monolithic `select case` blocks
- Pre-compute which mesh faces belong to each BC type, so BC routines loop only over relevant faces
- Support GPU-accelerated BC kernels with per-BC device arrays

## Architecture Overview

The system is built on three components:

1. **`BoundaryCondition`** -- a node in a linked list that holds a BC identifier, a procedure pointer to the BC implementation, and arrays of element/side indices
2. **`BoundaryConditionList`** -- a doubly-linked list that manages registration, lookup, and iteration over `BoundaryCondition` nodes
3. **`SELF_bcMethod`** -- an abstract interface that all BC implementations must satisfy

These are defined in `src/SELF_BoundaryConditions.f90`.

### The `BoundaryCondition` Type

```fortran
type BoundaryCondition
  procedure(SELF_bcMethod), pointer :: bcMethod => null()
  integer :: bcid
  character(SELF_BCNAME_LENGTH) :: bcname
  integer :: nBoundaries
  integer, allocatable :: elements(:)
  integer, allocatable :: sides(:)
  type(c_ptr) :: elements_gpu = c_null_ptr
  type(c_ptr) :: sides_gpu = c_null_ptr
  type(BoundaryCondition), pointer :: next => null()
  type(BoundaryCondition), pointer :: prev => null()
endtype
```

| Field | Purpose |
|-------|---------|
| `bcMethod` | Procedure pointer to the BC implementation |
| `bcid` | Integer constant identifying the BC type (e.g., `SELF_BC_PRESCRIBED`) |
| `bcname` | Human-readable name for diagnostics |
| `nBoundaries` | Number of boundary faces that carry this BC |
| `elements(:)` | Element indices for each boundary face |
| `sides(:)` | Local side indices for each boundary face |
| `elements_gpu`, `sides_gpu` | Device pointers used by GPU kernels |

### The `SELF_bcMethod` Interface

Every BC implementation must match this signature:

```fortran
subroutine SELF_bcMethod(this, mymodel)
  class(BoundaryCondition), intent(in) :: this
  class(Model), intent(inout) :: mymodel
endsubroutine
```

The BC receives itself (providing access to `elements`, `sides`, and `nBoundaries`) and the model (providing access to solution data). Implementations use `select type` to downcast `mymodel` to the concrete model type.

### The `BoundaryConditionList` Type

Each `DGModel` carries two lists:

- `hyperbolicBCs` -- for boundary conditions on the solution (used by `SetBoundaryCondition`)
- `parabolicBCs` -- for boundary conditions on the solution gradient (used by `SetGradientBoundaryCondition`)

Key methods:

| Method | Purpose |
|--------|---------|
| `Init()` | Initialize an empty list |
| `Free()` | Deallocate all nodes |
| `RegisterBoundaryCondition(bcid, bcname, bcfunc)` | Add a new BC or update an existing one |
| `GetBCForID(bcid)` | Return the node for a given `bcid`, or `null()` |
| `PopulateBoundaries(bcid, nBoundaries, elements, sides)` | Fill element/side arrays after mesh scanning |

If `RegisterBoundaryCondition` is called with a `bcid` that is already registered, it updates the procedure pointer without creating a new node. This is how GPU model variants override CPU implementations.

## BC Identifiers

A `bcid` is an arbitrary integer. The system never interprets it: `GetBCForID` compares it for
equality against the registered nodes, and `MapBoundaryConditions` compares it against
`sideInfo(5,...)`. Nothing else looks at the value, so any integer is legal as long as the mesh
tagging and the registration agree.

`src/SELF_Mesh.f90` defines five parameters:

```fortran
! Conditions on the solution
integer, parameter :: SELF_BC_PRESCRIBED = 100
integer, parameter :: SELF_BC_RADIATION = 101
integer, parameter :: SELF_BC_NONORMALFLOW = 102

! Conditions on the solution gradients
integer, parameter :: SELF_BC_PRESCRIBED_STRESS = 200
integer, parameter :: SELF_BC_NOSTRESS = 201
```

These are **deprecated as a fixed enumeration**. They are the ids the built-in models happen to
register, retained because the existing models, tests and examples reference them. They are not
a list of the conditions SELF supports, and no new ones will be added — a model defines the ids
it needs in its own module.

The same five values are duplicated as `#define`s in `src/gpu/SELF_GPU_Macros.h` for device
code. Nothing checks that the two lists agree, so they are kept in sync by hand. Device kernels
dispatch from a per-BC element/side list rather than by comparing against these macros, so a new
`bcid` does not need an entry there.

## Initialization Flow

The BC system is initialized as part of model creation. The sequence is:

```
Model%Init(mesh, geometry)
  |
  +-- hyperbolicBCs%Init()
  +-- parabolicBCs%Init()
  +-- AdditionalInit()          <-- subclass registers BCs here
  +-- MapBoundaryConditions()   <-- scans mesh, populates element/side arrays
```

### Step 1: Register BCs in `AdditionalInit`

Each model subclass overrides `AdditionalInit` to register its supported BC types:

```fortran
subroutine AdditionalInit_ECAdvection2D_t(this)
  class(ECAdvection2D_t), intent(inout) :: this
  procedure(SELF_bcMethod), pointer :: bcfunc

  bcfunc => hbc2d_NoNormalFlow_ECAdvection2D
  call this%hyperbolicBCs%RegisterBoundaryCondition( &
    SELF_BC_NONORMALFLOW, "no_normal_flow", bcfunc)
endsubroutine
```

### Step 2: Map Mesh Faces in `MapBoundaryConditions`

After registration, `MapBoundaryConditions` scans the mesh `sideInfo` array. For each registered BC, it performs two passes:

1. **Count** how many boundary faces carry that `bcid`
2. **Collect** the element and side indices into arrays

These arrays are stored in the `BoundaryCondition` node via `PopulateBoundaries`. A boundary face is identified by `sideInfo(3,j,iEl) == 0` (no neighbor element) and `sideInfo(5,j,iEl) == bcid`.

### Step 3: The Reverse Scan

Both mapping loops iterate over *registrations*, not over mesh faces, so a boundary face whose
`bcid` matches no registration is never enumerated by them. `MapBoundaryConditions` therefore
finishes with a scan in the opposite direction — over mesh faces, looking each `bcid` up in the
**hyperbolic** list — and records what it finds on the model:

| Field | Meaning |
|-------|---------|
| `nUnmappedBoundaries` | boundary faces with no hyperbolic registration, summed over all ranks |
| `unmappedBoundaryID` | the first such `bcid` found; meaningful only when the count is nonzero, since `-1` is itself a legal `bcid` |
| `unmappedBoundariesReported` | one-shot latch, cleared on every `MapBoundaryConditions` |

Only the hyperbolic list decides whether a face is handled. `SetBoundaryCondition` dispatches
that list alone, and it is what writes `solution%extBoundary` — the trace the Riemann solver
consumes. The parabolic list writes `solutionGradient%extBoundary` by way of
`SetGradientBoundaryCondition`, so a `bcid` registered only parabolically still leaves the
solution trace unwritten, and is reported.

Three details matter in that scan:

- **Mortar sides are excluded.** They carry `sideInfo(3) = 0` exactly like a physical boundary,
  with `sideInfo(1)` holding the mortar index. The exclusion keys on `sideInfo(1)`, and only on
  a mesh with `nMortars > 0` — the HOPr readers copy `sideInfo(1)` verbatim from the file, where
  it is the HOPr side type and may be nonzero on an ordinary face. `ResetBoundaryConditionType`
  in `SELF_Mesh_{2,3}D_t` applies the identical test, for the identical reason: a reset that
  tagged a mortar side would put an interior face into a boundary condition's element/side
  list, and that condition would then overwrite the exterior state the mortar exchange had
  just written. The `nMortars` gate is load-bearing in both places — without it, a reset on a
  HOPr mesh silently tags nothing.
- **The count is reduced across ranks.** Each rank owns a slice of the mesh, so a `bcid` absent
  locally may be present elsewhere. `MapBoundaryConditions` `mpi_allreduce`s the count
  (`MPI_SUM`). The representative id is *not* reduced: a bcid is any integer, so no value can
  mark "absent". Ranks instead agree by `MPI_MIN` on the lowest-numbered rank that actually
  holds an offender — a rank without one bids `nRanks` and can never win — and that rank
  broadcasts its first offending id. Every rank reaches the routine on both the
  `Init` and the `Regrid` path, so the collective cannot deadlock.
- **1-D differs.** `Mesh1D` is replicated on every rank, so no reduction is needed, and
  `SetBoundaryCondition` seeds a periodic default before dispatching, so an unmapped endpoint
  falls back to periodic rather than to a zero exterior state. A `bcid` of `0` in 1-D is that
  periodic default and never counts.

`ForwardStep` calls `ReportUnmappedBoundaries` before its first step. The base `Model`
implementation is a no-op — `Model` is abstract and carries neither a mesh nor the BC lists —
and each `DGModel{1,2,3}D_t` overrides it to print a warning, once, from rank 0. It is a warning
and not an error on purpose: a zero exterior state is a meaningful condition for some systems.

## Runtime Dispatch

During time integration, `SetBoundaryCondition` iterates through the linked list and calls each registered BC:

```fortran
subroutine SetBoundaryCondition(this)
  class(DGModel2D_t), intent(inout) :: this
  type(BoundaryCondition), pointer :: bc
  procedure(SELF_bcMethod), pointer :: apply_bc

  bc => this%hyperbolicBCs%head
  do while (associated(bc))
    apply_bc => bc%bcMethod
    call apply_bc(bc, this)
    bc => bc%next
  enddo
endsubroutine
```

The same pattern applies to `SetGradientBoundaryCondition` for parabolic BCs.

## GPU Acceleration

GPU-enabled models follow a layered pattern:

1. The CPU model class (e.g., `ECAdvection2D_t`) registers a Fortran BC implementation in its `AdditionalInit`
2. The GPU model class (e.g., `ECAdvection2D` in `src/gpu/`) extends the CPU class and:
    - Calls the parent `AdditionalInit` to register the CPU version
    - Re-registers the same `bcid` with a GPU wrapper function, which replaces the procedure pointer
3. During `Init`, the GPU class uploads `elements` and `sides` arrays to device memory (`elements_gpu`, `sides_gpu`)
4. During `Free`, the GPU class deallocates device arrays

### GPU Wrapper Pattern

A GPU wrapper is a Fortran subroutine matching `SELF_bcMethod` that calls a C/C++ kernel:

```fortran
subroutine hbc2d_Mirror_ECAdvection2D_GPU_wrapper(bc, mymodel)
  class(BoundaryCondition), intent(in) :: bc
  class(Model), intent(inout) :: mymodel

  select type (m => mymodel)
  class is (ECAdvection2D)
    if (bc%nBoundaries > 0) then
      call hbc2d_mirror_ecadvection2d_gpu( &
        m%solution%extBoundary_gpu, &
        m%solution%boundary_gpu, &
        bc%elements_gpu, bc%sides_gpu, &
        bc%nBoundaries, m%solution%interp%N, &
        m%solution%nElem, m%solution%nvar)
    endif
  endselect
endsubroutine
```

The C++ kernel receives device pointers and iterates over the pre-filtered boundary face list:

```cpp
__global__ void hbc2d_mirror_ecadvection2d_kernel(
    real *extBoundary, real *boundary,
    int *elements, int *sides,
    int nBoundaries, int N, int nel, int nvar)
{
  // Thread indexing over DOFs, boundary faces, and variables
  // elements[n] and sides[n] identify which face to process
}
```

### GPU Memory Lifecycle

```
Init_ECAdvection2D(mesh, geometry)
  |
  +-- Init_ECDGModel2D_t()    (parent: registers CPU BCs, maps mesh)
  +-- for each BC in hyperbolicBCs:
        hipMalloc(elements_gpu)
        hipMemcpy(elements -> elements_gpu)
        hipMalloc(sides_gpu)
        hipMemcpy(sides -> sides_gpu)

Free_ECAdvection2D()
  |
  +-- for each BC in hyperbolicBCs:
        hipFree(elements_gpu)
        hipFree(sides_gpu)
  +-- Free_ECDGModel2D_t()    (parent: frees BC list nodes)
```

## Adding a New BC Type to a Built-in Model

To add a new BC type (e.g., an inflow condition) to an existing model:

1. **Define a BC ID in the model's own module**, not in `src/SELF_Mesh.f90`:

    ```fortran
    integer, parameter :: MYMODEL_BC_INFLOW = 103
    ```

    Keeping the parameter next to the model that registers it is the point: the id is private
    to the agreement between that model and the meshes it is run on, and adding to the
    `SELF_BC_*` list would imply a global enumeration that does not exist.

2. **Write the BC implementation** in the model's `_t` source file, matching the `SELF_bcMethod` interface

3. **Register it** in the model's `AdditionalInit`:

    ```fortran
    bcfunc => hbc2d_Inflow_MyModel
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      MYMODEL_BC_INFLOW, "inflow", bcfunc)
    ```

    `AdditionalInit` is the required place: `Regrid` frees both BC lists and re-runs
    `AdditionalInit` on every adaptation epoch, so a registration made anywhere else is lost
    after the first regrid.

4. **For GPU models**, write a C++ kernel and Fortran wrapper, then re-register in the GPU class `AdditionalInit`

5. **Tag mesh faces** with the new BC ID in the mesh setup (e.g., via `ResetBoundaryConditionType` or by setting `sideInfo(5,:,:)` appropriately)

## Key Source Files

| File | Contents |
|------|----------|
| `src/SELF_BoundaryConditions.f90` | `BoundaryCondition`, `BoundaryConditionList`, `SELF_bcMethod` interface |
| `src/SELF_Mesh.f90` | The deprecated built-in BC id parameters (`SELF_BC_PRESCRIBED`, etc.) |
| `src/SELF_DGModel{1D,2D,3D}_t.f90` | `MapBoundaryConditions`, `SetBoundaryCondition`, `SetGradientBoundaryCondition` |
| `src/SELF_Model.f90` | Base `AdditionalInit` / `AdditionalFree` / `ReportUnmappedBoundaries` stubs, and `ForwardStep` |
| `src/gpu/SELF_ECAdvection2D.f90` | Example GPU BC wrapper pattern |
| `src/gpu/SELF_ECAdvection2D.cpp` | Example GPU BC kernel |
