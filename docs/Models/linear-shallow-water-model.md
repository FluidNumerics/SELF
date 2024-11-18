# Linear Shallow Water Equations

## Definition
The [`SELF_LinearShallowWater2D_t` module](../ford/sourcefile/self_LinearShallowWater2D_t.f90.html) defines the [`LinearShallowWater2D_t` class](ford/type/LinearShallowWater2D_t.html). In SELF, models are posed in the form of a generic conservation law

$$
  \vec{s}_t + \nabla \cdot \overleftrightarrow{f} = \vec{q}
$$

where $\vec{s}$ is a vector of solution variables, $\overleftrightarrow{f}$ is the conservative flux, and $\vec{q}$ are non-conservative source terms. 

For the linear shallow water equations in 2-D

$$
    \vec{s} = 
    \begin{pmatrix}
        u \\ 
        v \\ 
        \eta
    \end{pmatrix}
$$

where $u$ and $v$ are the x and y components of the barotropic velocity ($\vec{u} =  u \hat{x} + v \hat{y}$) and $\eta$ is the deviation of the fluid free surface relative to the resting fluid.

$$
    \overleftrightarrow{f} = 
    \begin{pmatrix}
        g \eta \hat{x} \\ 
        g \eta \hat{y} \\ 
        H \vec{u}
    \end{pmatrix}
$$

where $g$ is acceleration due to gravity and $H$ is uniform resting fluid depth. The source term is set to zero.

$$
    \vec{q} = \vec{0}
$$

To track stability of the Euler equation, the total entropy function is

$$
    e = \frac{1}{2} \int_V H u^2 + H v^2 + g \eta^2 \hspace{1mm} dV
$$

## Implementation
The 2D Linear Shallow Water model is implemented as a type extension of the `DGModel2d` class. The `LinearShallowWater2D_t` class adds parameters for acceleration due to gravity and the uniform resting fluid depth. It also overrides `SetMetaData`, `entropy_func`, `flux2d`, and `riemannflux2d` type-bound procedures.

### Riemann Solver
The `LinearShallowWater2D` class is defined using the advective form.
The Riemann solver for the hyperbolic part of the shallow water equations is the local Lax-Friedrichs upwind Riemann solver

$$
    \overleftrightarrow{f} \cdot \hat{n} =
    \frac{1}{2}
    (\overleftrightarrow{f}_L \cdot \hat{n} + 
    \overleftrightarrow{f}_R \cdot \hat{n} +
    c(\vec{s}_L - \vec{s}_R))
$$

where $c = \sqrt{gH}$, and

$$
    \overleftrightarrow{f}_L \cdot \hat{n} =
    \begin{pmatrix}
        g \eta_L n_x \\ 
        g \eta_L n_y \\ 
        H \vec{u}_L \cdot \hat{n}
    \end{pmatrix}
$$

$$
    \overleftrightarrow{f}_R \cdot \hat{n} =
    \begin{pmatrix}
        g \eta_R n_x \\ 
        g \eta_R n_y \\ 
        H \vec{u}_R \cdot \hat{n}
    \end{pmatrix}
$$

Together, this becomes

$$
    \overleftrightarrow{f} \cdot \hat{n} =
    \frac{1}{2}
    \begin{pmatrix}
        (g \eta_L + g \eta_R + c(\vec{u}_L \cdot \hat{n} - \vec{u}_R \cdot \hat{n}))n_x \\
        (g \eta_L + g \eta_R + c(\vec{u}_L \cdot \hat{n} - \vec{u}_R \cdot \hat{n}))n_y \\
        H \vec{u}_L \cdot \hat{n} + H \vec{u}_R \cdot \hat{n} + c(\eta_L - \eta_R)
    \end{pmatrix}
$$

The details for this implementation can be found in [self_LinearShallowWater2D_t.f90](../ford/sourcefile/self_LinearShallowWater2D_t.f90.html).


### Boundary Conditions
When initializing the mesh for your 2D Linear Shallow Water Equations solver, you can change the boundary conditions to 

* `SELF_BC_Radiation` to set the external state on model boundaries to 0 in the Riemann solver
* `SELF_BC_NoNormalFlow` to set the external normal velocity to the negative of the interior normal velocity and prolong the density, pressure, and tangential velocity (free slip). This effectively creates a reflecting boundary condition.
* `SELF_BC_Prescribed` to set a prescribed external state.

As an example, when using the built-in structured mesh generator, you can do the following

```fortran

type(Mesh2D),target :: mesh
integer :: bcids(1:4)

bcids(1:4) = (/&
                SELF_NONORMALFLOW,& ! South boundary condition
                SELF_RADIATION,&    ! East boundary condition
                SELF_PRESCRIBED,&   ! North boundary condition
                SELF_RADIATION &    ! West boundary condition
              /)   
                            
call mesh%StructuredMesh(nxPerTile=5,nyPerTile=5,&
                         nTileX=2,nTileY=2,&
                         dx=0.1_prec,dy=0.1_prec,bcids)

```

!!! note
    See the [Structured Mesh documentation](../MeshGeneration/StructuredMesh.md) for details on using the `structuredmesh` procedure

!!! note
    To set a prescribed state as a function of position and time, you can create a type-extension of the `LinearShallowWater2D` class and override the [`hbc2d_Prescribed`](../ford/proc/hbc2d_prescribed_model.html) 

## Example usage

For examples, see any of the following

* [`examples/LinearShallowWater2D.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/LinearShallowWater2D.f90) - Implements the 2D shallow water equations with no normal flow boundary conditions