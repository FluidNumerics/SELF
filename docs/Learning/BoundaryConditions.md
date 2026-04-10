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

## BC Identifier Constants

BC type identifiers are integer parameters defined in `src/SELF_Mesh.f90`:

```fortran
! Conditions on the solution
integer, parameter :: SELF_BC_PRESCRIBED = 100
integer, parameter :: SELF_BC_RADIATION = 101
integer, parameter :: SELF_BC_NONORMALFLOW = 102

! Conditions on the solution gradients
integer, parameter :: SELF_BC_PRESCRIBED_STRESS = 200
integer, parameter :: SELF_BC_NOSTRESS = 201
```

These constants are used both when tagging mesh faces (in `sideInfo`) and when registering BCs.

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

1. **Define a BC ID** in `src/SELF_Mesh.f90` if one does not already exist:

    ```fortran
    integer, parameter :: SELF_BC_INFLOW = 103
    ```

2. **Write the BC implementation** in the model's `_t` source file, matching the `SELF_bcMethod` interface

3. **Register it** in the model's `AdditionalInit`:

    ```fortran
    bcfunc => hbc2d_Inflow_MyModel
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_INFLOW, "inflow", bcfunc)
    ```

4. **For GPU models**, write a C++ kernel and Fortran wrapper, then re-register in the GPU class `AdditionalInit`

5. **Tag mesh faces** with the new BC ID in the mesh setup (e.g., via `ResetBoundaryConditionType` or by setting `sideInfo(5,:,:)` appropriately)

## Key Source Files

| File | Contents |
|------|----------|
| `src/SELF_BoundaryConditions.f90` | `BoundaryCondition`, `BoundaryConditionList`, `SELF_bcMethod` interface |
| `src/SELF_Mesh.f90` | BC ID constants (`SELF_BC_PRESCRIBED`, etc.) |
| `src/SELF_DGModel{1D,2D,3D}_t.f90` | `MapBoundaryConditions`, `SetBoundaryCondition`, `SetGradientBoundaryCondition` |
| `src/SELF_Model.f90` | Base `AdditionalInit` / `AdditionalFree` stubs |
| `src/gpu/SELF_ECAdvection2D.f90` | Example GPU BC wrapper pattern |
| `src/gpu/SELF_ECAdvection2D.cpp` | Example GPU BC kernel |
