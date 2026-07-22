# Linear Euler 2D — Perfectly Matched Layer Tutorial

This tutorial shows how to use the [`LinearEuler2D_PML`](../../Models/linear-euler-2d-pml-model.md) model to absorb outgoing acoustic waves with a Perfectly Matched Layer. It covers the two supported ways of declaring where the PML lives in your mesh:

1. **Programmatic tagging** on a `StructuredMesh` (no external mesh tool needed).
2. **Material tagging in a HOHQMesh ISM-MM `.mesh` file**, for unstructured geometries.

The numerics and PML formulation are summarised in the [PML model reference](../../Models/linear-euler-2d-pml-model.md); this page focuses on the workflow.

## What `LinearEuler2D_PML` adds

`LinearEuler2D_PML` is a type extension of `LinearEuler2D`. Compared to the parent model it

* extends `nvar` from 4 to 7 (the original acoustic state plus three PML auxiliary variables `phi_u`, `phi_v`, `phi_P`);
* carries two per-node fields `sigma_x` and `sigma_y` that hold the PML damping coefficients;
* registers PML-aware no-normal-flow and radiation boundary handlers automatically.

You declare the PML region by tagging mesh elements with a material name that starts with `"pml"` (the prefix is configurable via `pml_material_prefix`). Inside those elements, `SetPMLProfile` builds a polynomial ramp so $\sigma$ rises smoothly from zero at the interior/PML interface to its peak value at the outer boundary.

## Common setup

Both workflows share the same initialisation skeleton:

```fortran
use self_data
use self_lineareuler2d_pml
use self_mesh_2d

type(LinearEuler2D_PML) :: modelobj
type(Lagrange),target   :: interp
type(Mesh2D),target     :: mesh
type(SEMQuad),target    :: geometry

! ... build the mesh (the two workflows differ here) ...
! ... ensure mesh%materialNames contains an entry beginning with "pml"
!     and mesh%elemMaterial(iel) points to it for every PML element ...

call interp%Init(N=5, controlNodeType=GAUSS, &
                 M=10, targetNodeType=UNIFORM)

call geometry%Init(interp, mesh%nElem)
call geometry%GenerateFromMesh(mesh)

call modelobj%Init(mesh, geometry)
modelobj%rho0 = 1.0_prec

call modelobj%SetPMLProfile(x_interior_min = 0.0_prec, &
                            x_interior_max = 4.0_prec, &
                            y_interior_min = 0.0_prec, &
                            y_interior_max = 2.0_prec, &
                            pml_width      = 2.0_prec, &
                            sigma_max      = 20.0_prec, &
                            ramp_exponent  = 3)
```

`SetPMLProfile` only writes non-zero $\sigma_x, \sigma_y$ into elements whose material name starts with the configured prefix; every other element is left at $\sigma = 0$. The interior bounding box `[x_interior_min, x_interior_max] × [y_interior_min, y_interior_max]` defines where you want the wave equation untouched — the ramp grows with distance *outside* that box.

After this point the PML is fully active. Initial-condition setup, time integration, and IO follow the same patterns as the [parent `LinearEuler2D`](../../Models/linear-euler-2d-model.md).

!!! note
    Always set the three PML auxiliary variables to zero in your initial condition:

    ```fortran
    this%solution%interior(i,j,iel,5:7) = 0.0_prec
    ```

    Their initial value is the integral $\int_0^t \vec{q}\,dt$, so starting them at zero corresponds to "the simulation has just begun".

## Workflow 1: programmatic tagging on a structured mesh

This is the easiest path: build a tile mesh with `StructuredMesh`, then walk the element list and decide which elements are PML based on their centroid.

```fortran
integer  :: bcids(1:4), iel
real(prec) :: xc
real(prec),parameter :: x_interior_max = 4.0_prec

bcids = [SELF_BC_NONORMALFLOW, & ! south
         SELF_BC_NONORMALFLOW, & ! east (sits behind the PML)
         SELF_BC_NONORMALFLOW, & ! north
         SELF_BC_NONORMALFLOW]   ! west

! 12 elements in x by 4 in y, each 0.5 wide. Interior is x in [0,4],
! the east-side PML occupies x in [4,6] (4 elements, 2 units thick).
call mesh%StructuredMesh(12, 4, 1, 1, 0.5_prec, 0.5_prec, bcids)

! Replace the default single-material table with our own.
if (allocated(mesh%materialNames)) deallocate(mesh%materialNames)
mesh%nMaterials = 2
allocate(mesh%materialNames(1:2))
mesh%materialNames(1) = "interior"
mesh%materialNames(2) = "pml"

! Tag any element whose centroid sits east of x_interior_max as PML.
do iel = 1, mesh%nElem
  xc = sum(mesh%nodeCoords(1,:,:,iel)) / &
       real(size(mesh%nodeCoords,2) * size(mesh%nodeCoords,3), prec)
  if (xc > x_interior_max) then
    mesh%elemMaterial(iel) = 2  ! pml
  else
    mesh%elemMaterial(iel) = 1  ! interior
  end if
end do
```

Things to know:

* `StructuredMesh` initialises `materialNames` to a single entry `"default"` and points every element at it. Reallocating the table and rewriting `elemMaterial` is safe — nothing else in `LinearEuler2D_PML` depends on the original table.
* The centroid is computed as a corner-node average. For the default linear (`nGeo = 1`) structured mesh this is exact; if you ever increase `nGeo` you should still average across the geometry nodes shown above.
* For PML on multiple sides, just generalise the centroid test (e.g. `if (xc < x_min .or. xc > x_max .or. yc < y_min .or. yc > y_max) elemMaterial = 2`). All PML elements share one material name; `SetPMLProfile` figures out the per-direction damping from the node positions.

The complete file is at [`examples/linear_euler2d_pml_planewave.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_pml_planewave.f90).

### Outer boundaries behind the PML

You still need a boundary condition on the outer edge of the PML elements — the PML attenuates but does not eliminate the outgoing wave entirely, and whatever is left will hit whatever you put behind it. The two natural choices are

* `SELF_BC_NONORMALFLOW` (default in the example above) — anything that survives the PML reflects back into the PML, where it is damped on the return trip. Cheap and robust.
* `SELF_BC_RADIATION` — sets the exterior acoustic state to zero, so the residual signal is mostly transmitted out instead of reflected. Slightly better absorption for thinner PMLs.

A well-designed PML (a few elements thick with $\sigma_\mathrm{max}$ tuned to give one or two decades of decay) makes the choice essentially invisible.

## Workflow 2: HOHQMesh ISM-MM mesh

When you need a curved or unstructured geometry, generate the mesh with HOHQMesh and emit it in the **ISM-MM** format. ISM-MM associates a material-name string with every element, exactly as needed by `SetPMLProfile`.

### Authoring the control file

HOHQMesh control files assign one material per geometry block. Mark the outer absorbing region with a name that starts with `pml` (the prefix `LinearEuler2D_PML` looks for); any other name is treated as interior.

A minimal sketch of a HOHQMesh control file with an interior disk surrounded by a PML annulus:

```
\begin{MODEL}
   \begin{OUTER_BOUNDARY}
      \begin{PARAMETRIC_EQUATION_CURVE}
         name     = pml_outer
         xEqn     = f(t) = R_out*cos(2*pi*t)
         yEqn     = f(t) = R_out*sin(2*pi*t)
         zEqn     = f(t) = 0.0
      \end{PARAMETRIC_EQUATION_CURVE}
   \end{OUTER_BOUNDARY}
   \begin{INNER_BOUNDARIES}
      \begin{CHAIN}
         name     = pml_interface       % shared edge between interior and PML
         \begin{PARAMETRIC_EQUATION_CURVE}
            name = pml_interface
            xEqn = f(t) = R_in*cos(2*pi*t)
            yEqn = f(t) = R_in*sin(2*pi*t)
            zEqn = f(t) = 0.0
         \end{PARAMETRIC_EQUATION_CURVE}
      \end{CHAIN}
   \end{INNER_BOUNDARIES}
\end{MODEL}

\begin{MATERIALS}
   \begin{MATERIAL}
      material name = interior          % anything not starting with "pml"
      material id   = 1
   \end{MATERIAL}
   \begin{MATERIAL}
      material name = pml               % MUST start with the configured prefix
      material id   = 2
   \end{MATERIAL}
\end{MATERIALS}
```

When HOHQMesh writes the resulting `.mesh` file in ISM-MM mode, every element block ends with the material name of the region it falls in. The reader picks it up automatically:

```
ISM-MM
  <nNodes> <nEdges> <nElem> <polyOrder>
  ... node coordinates ...
  <c1> <c2> <c3> <c4>   pml          % corner IDs + material name
  <curveFlag1> ... <curveFlag4>
  ... (curve sample points for any curved sides) ...
  <bdyName1> <bdyName2> <bdyName3> <bdyName4>
  ...
```

You can also mark the outer-edge boundary segments (those that survive after the PML attenuates the wave) with a single name like `outer` and remap them in your Fortran setup with `mesh%ResetBoundaryConditionType(SELF_BC_NONORMALFLOW)` or `SELF_BC_RADIATION`, exactly like the [BoneAndMarrow example](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_boneandmarrow.f90) does for its single `outer` boundary.

### Loading and configuring

The Fortran setup is even shorter than the structured-mesh version, because `mesh%Read_HOHQMesh` already populates `materialNames` and `elemMaterial`:

```fortran
character(LEN=255) :: WORKSPACE

call get_environment_variable("WORKSPACE", WORKSPACE)
call mesh%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/MyDomain/MyDomain.mesh")

! All physical boundaries are behind the PML. Map them to a single
! BC (radiation here; no-normal-flow is equally valid).
call mesh%ResetBoundaryConditionType(SELF_BC_RADIATION)

call interp%Init(N=5, controlNodeType=GAUSS, M=10, targetNodeType=UNIFORM)
call geometry%Init(interp, mesh%nElem)
call geometry%GenerateFromMesh(mesh)

call modelobj%Init(mesh, geometry)
modelobj%rho0 = 1.0_prec

call modelobj%SetPMLProfile(x_interior_min = -R_in, x_interior_max = R_in, &
                            y_interior_min = -R_in, y_interior_max = R_in, &
                            pml_width      = R_out - R_in, &
                            sigma_max      = 20.0_prec)
```

Things to know:

* `Read_HOHQMesh` auto-detects ISM vs ISM-MM from the file header. Plain ISM files have no material strings, so every element ends up tagged as `"default"` and no PML is applied — make sure your mesh writer is set to ISM-MM.
* You can use multiple `pml*` materials in the same mesh (`pml_east`, `pml_corner_NE`, etc.). `LinearEuler2D_PML` only cares about the prefix; it does not distinguish between them. The per-direction strength is set by the geometric position of each node, not by the material name.
* The interior bounding box you pass to `SetPMLProfile` is purely a geometric construct. For an annular PML around a disk it is convenient to use a circumscribed square, as in the snippet above: $d_x = d_y = 0$ in the interior, and they grow as nodes move outward.

## Verifying that the PML works

Two quick sanity checks before running production simulations:

1. **Stability.** Run the simulation long enough for the pulse to fully enter the PML and check that `entropy` is finite and non-increasing. The bundled regression test at [`test/lineareuler2d_pml_absorption.f90`](https://github.com/FluidNumerics/SELF/blob/main/test/lineareuler2d_pml_absorption.f90) does exactly this and is a good template.

2. **Absorption.** Pull the solution back to the host and measure `max(abs(p))` over the interior elements (those with `xc <= x_interior_max`, etc.). After the wave has had time to traverse the PML, the residual should be a small fraction of the initial peak.

If absorption is weaker than expected:

* Make the PML thicker — three to four elements is a sensible starting point at $N = 5$.
* Increase `sigma_max`. RK3 is stable for `dt * sigma_max < ~2.5`, so you can usually push $\sigma_\mathrm{max}$ up to $\mathcal{O}(10)$ relative to the inverse PML traversal time.
* Use a smoother ramp (`ramp_exponent = 3` is the default; bump to 4 if you see interface reflections).

If the simulation goes unstable, check first that you did *not* leave `phi_*` set to anything other than zero in the initial condition, that `dt * sigma_max` is well below the RK stability limit, and that the PML thickness is at least a couple of elements.

## Running the bundled example

The bundled example uses the programmatic tagging workflow. After [installing SELF](../../GettingStarted/install.md):

```shell
${SELF_ROOT}/examples/linear_euler2d_pml_planewave
```

This launches a 2D Gaussian acoustic pulse near the western interior, propagates it eastward into a PML zone, and writes `solution.*.h5` snapshots every 0.5 time units. Open them in PySELF or ParaView to watch the pulse decay smoothly through the PML.
