# Linear Euler (2D)

## Definition
The [`SELF_LinearEuler2D_t` module](../ford/module/self_lineareuler2d_t.html) defines the [`LinearEuler2D_t` class](ford/type/lineareuler2d_t.html). In SELF, models are posed in the form of a generic conservation law

\begin{equation}
\vec{s}_t + \nabla \cdot \overleftrightarrow{f} = \vec{q}
\end{equation}

where $\vec{s}$ is a vector of solution variables, $\overleftrightarrow{f}$ is the conservative flux, and $\vec{q}$ are non-conservative source terms. 

For the Linear Euler equations in 2-D

\begin{equation}
\vec{s} = \begin{pmatrix}
\rho \\ 
u \\ 
v \\ 
p
\end{pmatrix}
\end{equation}
where $\rho$ is a density anomaly, referenced to the density $\rho_0$, $u$ and $v$ are the $x$ and $y$ components of the fluid velocity (respectively), and $p$ is the pressure. When we assume an ideal gas, and a motionless background state, the conservative fluxes are

\begin{equation}
\overleftrightarrow{f} = \begin{pmatrix}
\rho_0(u \hat{x} + v \hat{y}) \\
p \hat{x} \\
p \hat{y} \\
\rho_0c^2(u \hat{x} + v \hat{y})
\end{pmatrix}
\end{equation}
where $c$ is the (constant) speed of sound. The source term is set to zero.

\begin{equation}
\vec{q} = \vec{0}
\end{equation}

To track stability of the Euler equation, the total entropy function is

\begin{equation}
e = \frac{1}{2} \int_V u^2 + v^2 + \frac{p}{\rho_0 c^2} \hspace{1mm} dV
\end{equation}

## Implementation
The Linear Euler 2D model is implemented as a type extension of the [`DGModel2D` class](../ford/type/dgmodel2d_t.html). The [`LinearEuler2D_t` class](../ford/type/lineareuler2d_t.html) adds parameters for the reference density and the speed speed of sound and overrides the `SetMetadata`, `entropy_func`, `flux2d`, and `riemannflux2d` type-bound procedures.

### Riemann Solver
The `LinearEuler2D` class is defined using the conservative form of the conservation law. The Riemman solver for the hyperbolic part of Euler equation is the local Lax Friedrichs upwind riemann solver

\begin{equation}
\overleftrightarrow{f}_h^* \cdot \hat{n} = \frac{1}{2}( \overleftrightarrow{f}_L \cdot \hat{n}  + \overleftrightarrow{f}_R \cdot \hat{n}  + c (\vec{s}_L - \vec{s}_R))
\end{equation}

where 
\begin{equation}
f_L = 
\end{equation}

\begin{equation}
f_R = 
\end{equation}

The details for this implementation can be found in [`self_lineareuler2d_t.f90`](../ford/sourcefile/self_lineareuler2d_t.f90.html)

### Boundary conditions
When initializing the mesh for your Euler 2D equation solver, you can change the boundary conditions to 

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
    To set a prescribed state as a function of position and time, you can create a type-extension of the `LinearEuler2D` class and override the [`hbc2d_Prescribed`](../ford/proc/hbc2d_prescribed_model.html) 


## Example usage

For examples, see any of the following

* [`examples/lineareuler2d_spherical_soundwave_closeddomain.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/linear_euler2d_spherical_soundwave_closeddomain.f90) - Implements a simulation with a gaussian pressure and density anomaly as an initial condition in a domain with no normal flow boundary conditions on all sides.