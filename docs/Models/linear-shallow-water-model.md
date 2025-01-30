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

where $g$ is acceleration due to gravity and $H$ is uniform resting fluid depth. The source term includes a coriolis force

$$
    \vec{q} = 
        \begin{pmatrix}
        -fv - C_d u\\ 
        fu - C_d v\\ 
        0
    \end{pmatrix}
$$

where $f$ is the coriolis parameter and $C_d$ is the linear drag coefficient.

To track stability of the shallow water equations, the total entropy function is taken to be the total (kinetic plus potential) energy

$$
    e = \frac{1}{2} \int_V H u^2 + H v^2 + g \eta^2 \hspace{1mm} dV
$$

## Implementation
The 2D Linear Shallow Water model is implemented as a type extension of the `DGModel2d` class. The `LinearShallowWater2D_t` class adds parameters for acceleration due to gravity and the uniform resting fluid depth. It also overrides `SetMetaData`, `entropy_func`, `flux2d`, and `riemannflux2d` type-bound procedures.

### Defining the coriolis parameter and geostrophic velocities
The `LinearShallowWater2D` class has a generic method (`SetCoriolis`) that can be used for defining the coriolis parameter at each location in the model domain. The `SetCoriolis` method can be used for either setting an $f$ or $beta$ plane.

#### Setting up an f-plane
Assuming you've created interpolant ,mesh, geometry objects, and model objects you can define a constant value for the coriolis parameter using the following
```fortran
type(LinearShallowWater2D) :: modelobj
real(prec), parameter :: f0 = 10.0_prec*(-4)
...

  call modelobj%SetCoriolis(f0)

```

#### Setting up a beta-plane
Assuming you've created interpolant ,mesh, geometry objects, and model objects you can define the coriolis so that it varies with the `y` coordinate in the geometry using
```fortran
type(LinearShallowWater2D) :: modelobj
real(prec), parameter :: f0 = 10.0_prec*(-4)
real(prec), parameter :: beta = 10.0_prec*(-11) 
...

  call modelobj%SetCoriolis(f0,beta)

```

#### Setting arbitrary spatially varying coriolis parameter
Perhaps you find that f-plane and beta-plane scenarios are just too boring, or their not an appropriate model for what you're considering. In this case, you can easily set the `fCori%interior` attribute of the `LinearShallowWater2D` class directly


```fortran
type(LinearShallowWater2D) :: modelobj
integer :: iel
integer :: i
integer :: j

do concurrent(i=1:modelobj%solution%N+1,j=1:modelobj%solution%N+1,iel=1:modelobj%mesh%nElem)
    x = modelobj%geometry%x%interior(i,j,1,iel,1) ! Get the x-coordinate
    y = modelobj%geometry%x%interior(i,j,1,iel,2) ! Get the y-coordinate
    this%fCori%interior(i,j,1,iel) =  ! Define the coriolis parameter here as a function of x and y
enddo
call this%fCori%UpdateDevice()
```

#### Defining Geostrophic velocities
With the `fCori` attribute defined, you can define geostrophic velocities from an initial condition for the free-surface height.

!!! note
    Setting geostrophic velocities is only valid when $f \neq 0$ everywhere in the domain.

To define geostrophic velocities, you can simply use the `DiagnoseGeostrophicVelocity` procedure. This will 

* Reset the velocity field to zero,
* Calculate the free-surface height gradients using the `CalculateTendency` method
* Compute the `u` and `v` variables using geostrophic balance

As an example,

```fortran
type(LinearShallowWater2D) :: modelobj
real(prec), parameter :: f0 = 10.0_prec*(-4)
real(prec), parameter :: beta = 10.0_prec*(-11) 
...

  call modelobj%SetCoriolis(f0,beta)  

  ! Set the free-surface height using an equation string
  call modelobj%solution%SetEquation(3,'f = 0.01*exp( -( (x-500000.0)^2 + (y-500000.0)^2 )/(2.0*(10.0^10)) )')
  call modelobj%solution%SetInteriorFromEquation(geometry,0.0_prec)

  ! Calculate u and v from the free surface height using geostrophy
  call modelobj%DiagnoseGeostrophicVelocity()

```

### Setting the Drag coefficient
Assuming you've created interpolant ,mesh, geometry objects, and model objects you can define a constant value for the linear drag coefficient by setting the constant parameter `Cd`, e.g. 

```fortran
type(LinearShallowWater2D) :: modelobj
real(prec), parameter :: fCd = 0.25
...

  modelobj % Cd = Cd ! Set the drag coefficient

```
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

* [Gravity waves in closed square domain](../Tutorials/LinearShallowWater/LinearShallowWater.md)
* [Kelvin waves in a closed circular rotating domain (f-plane)](../Tutorials/LinearShallowWater/KelvinWaves.md)
* [Planetary Rossby waves in an open square domain (beta-plane)](../Tutorials/LinearShallowWater/PlanetaryRossbyWave.md)