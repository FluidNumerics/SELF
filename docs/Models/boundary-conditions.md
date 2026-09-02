# Boundary Conditions

SELF uses an extensible boundary condition system that lets you register custom boundary conditions for your models. This page explains how to use boundary conditions in practice.

For details on the internal architecture, see [Boundary Condition System (Developer Guide)](../Learning/BoundaryConditions.md).

## Overview

Every DG model in SELF maintains two boundary condition lists:

- **Hyperbolic BCs** -- conditions on the solution at domain boundaries
- **Parabolic BCs** -- conditions on the solution gradient at domain boundaries

You populate these lists by overriding the `AdditionalInit` method in your model subclass. SELF then automatically determines which mesh faces match each BC and calls your BC routines during time integration.

## BC identifiers are just integers

A boundary condition is identified by an integer, the `bcid`. **Any** integer works. A `bcid`
is meaningful only because two places agree on it:

1. the value stored in the mesh — `sideInfo(5,:,:)` in 2-D and 3-D, `mesh%bcid(1:2)` in 1-D
2. the value passed to `RegisterBoundaryCondition`

Nothing else interprets it. There is no registry of legal values, no `select case` on it, and
no requirement to pick from a predefined set.

### The ids the built-in models use

The following constants are defined in `SELF_Mesh`. They exist because the built-in models
use them, and they are what you tag your mesh with when you use a built-in model:

| Constant | Value | Category | Description |
|----------|-------|----------|-------------|
| `SELF_BC_PRESCRIBED` | 100 | Hyperbolic | Dirichlet-type: set the exterior state to a prescribed value |
| `SELF_BC_RADIATION` | 101 | Hyperbolic | Radiation / non-reflecting outflow |
| `SELF_BC_NONORMALFLOW` | 102 | Hyperbolic | Mirror / no-normal-flow (wall) |
| `SELF_BC_PRESCRIBED_STRESS` | 200 | Parabolic | Prescribed stress on the gradient |
| `SELF_BC_NOSTRESS` | 201 | Parabolic | Zero-stress (free-slip) |

!!! warning "These are not an enumeration of what SELF supports"
    They are being phased out as a fixed list, and no new ones will be added. Do not treat
    them as the set of boundary conditions available to you, and do not add to them when you
    need a new condition — define your own, as below.

### Defining your own ids

Declare the ids your model needs in your own module, next to the model that registers them.
That keeps the id and its implementation together, and does not require touching `SELF_Mesh`:

```fortran
module my_inflow_model
  use self_lineareuler2d
  use SELF_BoundaryConditions
  implicit none

  ! This model's own boundary condition id. Any integer will do; it only has to match what
  ! the mesh is tagged with. Every id you tag needs its own registration.
  integer, parameter :: MYMODEL_BC_INFLOW = 1

  type, extends(LinearEuler2D) :: my_inflow
  contains
    procedure :: AdditionalInit => AdditionalInit_my_inflow
  endtype

contains

  subroutine AdditionalInit_my_inflow(this)
    class(my_inflow), intent(inout) :: this
    procedure(SELF_bcMethod), pointer :: bcfunc

    ! Keep the parent's registrations (no_normal_flow, radiation). Dispatch through the
    ! parent COMPONENT rather than naming a routine: the CPU LinearEuler2D inherits
    ! AdditionalInit_LinearEuler2D_t, while the GPU one overrides it to register device
    ! kernels. Calling the _t routine by name installs the host handlers on a GPU build and
    ! the device kernels never replace them, leaving the boundary conditions inert there.
    call this%LinearEuler2D%AdditionalInit()

    bcfunc => hbc2d_Inflow_my_inflow
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      MYMODEL_BC_INFLOW, "inflow", bcfunc)
  endsubroutine
endmodule
```

Small ids like `1` and `2` are convenient here for another reason: a mesh read from a
HOHQMesh `.mesh` file carries a 1-based index into `mesh%BCNames` in `sideInfo(5,...)`, so
the first boundary name in the file is id `1`, the second is `2`, and so on. Choosing your
ids to match lets you tag boundaries by name in the mesh generator and register against those
same ids in Fortran, with no remapping step.

## Diagnosing unhandled boundaries

A boundary face whose `bcid` has no boundary condition registered on `hyperbolicBCs` is
**not** applied. A registration on `parabolicBCs` alone does not count — that list is
dispatched by `SetGradientBoundaryCondition` and writes `solutionGradient%extBoundary`, not
the solution trace — so such a face is reported too.
Nothing writes `extBoundary` on that face, so the Riemann solver reads whatever is there:
zero on the first step, and the previous step's values afterwards. For some systems a zero
exterior state is a meaningful condition (for linear Euler it is effectively radiation), so
SELF does not treat this as an error — but it does not hide it either.

`MapBoundaryConditions` counts those faces at `Init` and after every `Regrid`, and the first
`ForwardStep` prints a warning naming the count and one of the unregistered ids:

```
 src/SELF_DGModel2D_t.f90 : Warning :  40  mesh boundary edges carry a bcid with no registered boundary condition.
 src/SELF_DGModel2D_t.f90 : Warning : One of the unregistered bcids is  999
```

The count is summed over all MPI ranks and printed once, by rank 0. It is also readable from
your own code as `model%nUnmappedBoundaries`, with a representative id in
`model%unmappedBoundaryID`.

Common causes:

- a typo in a `bcids(:)` entry, or a `bcid` the model does not register (for instance
  `SELF_BC_NONORMALFLOW` on `LinearEuler3D`, which registers radiation only)
- a mesh read from a HOHQMesh `.mesh` file, where `sideInfo(5,...)` carries `BCNames`
  indices, without a `ResetBoundaryConditionType` call or matching registrations
- registering boundary conditions somewhere other than `AdditionalInit`. `Regrid` frees the
  BC lists and re-runs `AdditionalInit` on every adaptation epoch, so registrations made
  elsewhere are lost after the first one.

In 1-D the fallback differs: an unregistered endpoint keeps the periodic default that
`SetBoundaryCondition` seeds, rather than a zero exterior state. A `bcid` of `0` in 1-D *is*
the periodic default and is never reported.

## Workflow

Setting up boundary conditions involves three steps:

1. **Tag mesh faces** with BC identifiers
2. **Write BC subroutines** that set exterior state values
3. **Register BCs** in your model's `AdditionalInit`

### Step 1: Tag Mesh Faces

After creating your mesh, assign BC identifiers to boundary faces. The simplest way is `ResetBoundaryConditionType`:

```fortran
! Set left boundary to prescribed, right boundary to prescribed
call mesh%ResetBoundaryConditionType(SELF_BC_PRESCRIBED, SELF_BC_PRESCRIBED)
```

For 2D and 3D meshes, boundary face tagging is typically set via the mesh file or by modifying `sideInfo(5,:,:)` after mesh creation.

### Step 2: Write BC Subroutines

Each BC subroutine must match the `SELF_bcMethod` interface:

```fortran
subroutine my_bc(bc, mymodel)
  use SELF_BoundaryConditions
  use SELF_Model
  class(BoundaryCondition), intent(in) :: bc
  class(Model), intent(inout) :: mymodel
endsubroutine
```

Inside the subroutine, use `select type` to access your model's data, then loop over `bc%nBoundaries` to set the exterior state on each boundary face:

```fortran
subroutine hbc1d_Prescribed_mymodel(bc, mymodel)
  class(BoundaryCondition), intent(in) :: bc
  class(Model), intent(inout) :: mymodel
  integer :: n, iEl, s

  select type (m => mymodel)
  class is (my_model_type)
    do n = 1, bc%nBoundaries
      iEl = bc%elements(n)   ! element index
      s = bc%sides(n)        ! local side index
      ! Set the exterior boundary state
      m%solution%extBoundary(s, iEl, 1:m%nvar) = ...
    enddo
  endselect
endsubroutine
```

!!! note
    The `extBoundary` array holds the exterior (ghost) state used by the Riemann solver at domain boundaries. Setting `extBoundary = boundary` (the interior state) produces a mirror/no-normal-flow condition. Setting it to a specific value produces a Dirichlet condition.

### Step 3: Register BCs in `AdditionalInit`

Override `AdditionalInit` in your model subclass and register each BC:

```fortran
subroutine AdditionalInit_mymodel(this)
  class(my_model_type), intent(inout) :: this
  procedure(SELF_bcMethod), pointer :: bcfunc

  ! Register a hyperbolic BC
  bcfunc => hbc1d_Prescribed_mymodel
  call this%hyperbolicBCs%RegisterBoundaryCondition( &
    SELF_BC_PRESCRIBED, "prescribed", bcfunc)

  ! Register a parabolic BC (if needed)
  bcfunc => pbc1d_Prescribed_mymodel
  call this%parabolicBCs%RegisterBoundaryCondition( &
    SELF_BC_PRESCRIBED, "prescribed", bcfunc)
endsubroutine
```

SELF calls `AdditionalInit` during model initialization, before scanning the mesh. After your BCs are registered, SELF automatically determines which mesh faces belong to each BC type and stores the element/side arrays.

## Complete Example: Traveling Shock (Burgers 1D)

This example shows how to extend the built-in `burgers1D` model with prescribed boundary conditions for a traveling shock solution.

### Model Module

```fortran
module burgers1d_shock_model
  use self_Burgers1D
  use SELF_BoundaryConditions
  implicit none

  type, extends(burgers1D) :: burgers1d_shock
    real(prec) :: ul = 1.0_prec   ! Left state
    real(prec) :: ur = 0.0_prec   ! Right state
    real(prec) :: x0 = 0.1_prec   ! Initial shock position
  contains
    procedure :: AdditionalInit => AdditionalInit_burgers1d_shock
  endtype

contains

  subroutine AdditionalInit_burgers1d_shock(this)
    class(burgers1d_shock), intent(inout) :: this
    procedure(SELF_bcMethod), pointer :: bcfunc

    ! Register prescribed BC for the solution
    bcfunc => hbc1d_Prescribed_burgers1d_shock
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_PRESCRIBED, "prescribed", bcfunc)

    ! Register prescribed BC for the gradient
    bcfunc => pbc1d_Prescribed_burgers1d_shock
    call this%parabolicBCs%RegisterBoundaryCondition( &
      SELF_BC_PRESCRIBED, "prescribed", bcfunc)
  endsubroutine

  subroutine hbc1d_Prescribed_burgers1d_shock(bc, mymodel)
    class(BoundaryCondition), intent(in) :: bc
    class(Model), intent(inout) :: mymodel
    integer :: n, iEl, s
    real(prec) :: x, jump, spd

    select type (m => mymodel)
    class is (burgers1d_shock)
      do n = 1, bc%nBoundaries
        iEl = bc%elements(n)
        s = bc%sides(n)
        x = m%geometry%x%boundary(s, iEl, 1)
        jump = m%ul - m%ur
        spd = 0.5_prec*(m%ul + m%ur)
        m%solution%extBoundary(s, iEl, 1) = &
          spd - 0.5_prec*tanh((x - spd*m%t - m%x0)*jump &
                               / (4.0_prec*m%nu))
      enddo
    endselect
  endsubroutine

  subroutine pbc1d_Prescribed_burgers1d_shock(bc, mymodel)
    class(BoundaryCondition), intent(in) :: bc
    class(Model), intent(inout) :: mymodel
    integer :: n, iEl, s
    real(prec) :: x, jump, spd, r, drdx

    select type (m => mymodel)
    class is (burgers1d_shock)
      do n = 1, bc%nBoundaries
        iEl = bc%elements(n)
        s = bc%sides(n)
        x = m%geometry%x%boundary(s, iEl, 1)
        jump = m%ul - m%ur
        spd = 0.5_prec*(m%ul + m%ur)
        r = (x - spd*m%t - m%x0)*jump/(4.0_prec*m%nu)
        drdx = jump/(4.0_prec*m%nu)
        m%solutionGradient%extBoundary(s, iEl, 1) = &
          -0.5_prec*drdx*(2.0_prec/(exp(r) + exp(-r)))**2
      enddo
    endselect
  endsubroutine

endmodule
```

### Main Program

```fortran
program traveling_shock
  use self_data
  use burgers1d_shock_model
  implicit none

  type(burgers1d_shock) :: modelobj
  type(Lagrange), target :: interp
  type(Mesh1D), target :: mesh
  type(Geometry1D), target :: geometry

  ! Create mesh
  call mesh%StructuredMesh(nElem=10, x=(/0.0_prec, 1.0_prec/))

  ! Tag both boundaries as prescribed
  call mesh%ResetBoundaryConditionType(SELF_BC_PRESCRIBED, SELF_BC_PRESCRIBED)

  ! Create interpolant and geometry
  call interp%Init(N=7, controlNodeType=GAUSS, M=10, targetNodeType=UNIFORM)
  call geometry%Init(interp, mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  ! Initialize model (this calls AdditionalInit and MapBoundaryConditions)
  call modelobj%Init(mesh, geometry)
  modelobj%gradient_enabled = .true.
  modelobj%nu = 0.01_prec

  ! Set initial condition, run, clean up...
  call modelobj%SetTimeIntegrator('rk3')
  call modelobj%ForwardStep(2.0_prec, 1.0e-5_prec, 0.05_prec)
  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()
endprogram
```

## Extending a Built-in Model's BCs

When your model extends a built-in SELF model that already registers BCs, call the parent's `AdditionalInit` first, then register your additional BCs:

```fortran
subroutine AdditionalInit_my_extended_model(this)
  class(my_extended_model), intent(inout) :: this
  procedure(SELF_bcMethod), pointer :: bcfunc

  ! Register the parent's BCs first, through the parent COMPONENT. Naming
  ! AdditionalInit_LinearEuler2D_t directly installs the host handlers on a GPU build and
  ! stops the device kernels from ever replacing them.
  call this%LinearEuler2D%AdditionalInit()

  ! Then register additional BCs
  bcfunc => hbc2d_Prescribed_my_extended_model
  call this%hyperbolicBCs%RegisterBoundaryCondition( &
    SELF_BC_PRESCRIBED, "prescribed", bcfunc)
endsubroutine
```

## Replacing a BC Implementation

If a BC with a given ID is already registered and you call `RegisterBoundaryCondition` with the same ID again, the procedure pointer is updated to point to your new implementation. The element/side arrays are preserved. This mechanism is used internally by GPU model variants to replace CPU BC routines with GPU-accelerated kernels, but you can also use it to override a parent model's BC behavior.

## Tips

- **One BC per ID per list.** Each `bcid` can appear at most once in the hyperbolic list and once in the parabolic list.
- **Mesh tagging must match registration.** If you register a BC for `SELF_BC_PRESCRIBED` but no mesh faces carry that ID, the BC will have `nBoundaries = 0` and its subroutine will never be called. The reverse — mesh faces carrying an ID nothing is registered for — is reported by a warning on the first `ForwardStep`; see [Diagnosing unhandled boundaries](#diagnosing-unhandled-boundaries).
- **Use `select type` in BC routines.** The `SELF_bcMethod` interface receives `class(Model)`, so you must downcast to access your model's fields.
- **BC routines are called every time step.** Keep them efficient. For GPU models, use device kernels rather than host-side loops.
