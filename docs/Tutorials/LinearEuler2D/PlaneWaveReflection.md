# Linear Euler 2D - Plane Wave Reflection Tutorial
This tutorial will walk you through using an example program that uses the `LinearEuler2D` class to run a simulation with the linear Euler equations for an ideal gas in 2-D. This example is configured using the built in structured mesh generator with prescribed boundary conditions on north, west, and south boundaries and a no-normal-flow boundary condition on the east boundary.

## Problem statement

### Equations Solved
In this example, we are solving the linear Euler equations in 2-D for an ideal gas, given by

$$
  \vec{s}_t + \nabla \cdot \overleftrightarrow{f} = \vec{q}
$$

where


$$
    \vec{s} = 
    \begin{pmatrix}
    \rho \\ 
    u \\ 
    v \\ 
    p
    \end{pmatrix}
$$

and

$$
    \overleftrightarrow{f} = 
    \begin{pmatrix}
    \rho_0(u \hat{x} + v \hat{y}) \\
    p \hat{x} \\
    p \hat{y} \\
    \rho_0c^2(u \hat{x} + v \hat{y})
    \end{pmatrix}
$$

$$
    \vec{q} = \vec{0}
$$ 



The variables are defined as follows

* $\rho$ is a density anomaly referenced to the density $\rho_0$
* $u$ and $v$ are the $x$ and $y$ components of the fluid velocity (respectively)
* $p$ is the pressure
* $c$ is the (constant) speed of sound. 

### Model Domain
The physical domain is defined by $\vec{x} \in [0, 1]\times[0,1]$. We use the `StructuredMesh` routine to create a domain with 20 × 20 elements that are dimensioned 0.05 × 0.05 . Model boundaries on the south, north, and west edges of the domain are all tagged with the `SELF_BC_PRESCRIBED` flag to implement prescribed boundary conditions; the east boundary is tagged with `SELF_BC_NONORMALFLOW` to implement no-normal-flow boundary conditions.

Within each element, all variables are approximated by a Lagrange interpolating polynomial of degree 7. The interpolation knots are the Legendre-Gauss points.


### Initial and Boundary Conditions
The initial and prescribed boundary conditions are set using an exact solution. The exact solution is found using the method of images where a no-normal-flow wall is placed at $x=1$ . We define the solution as the sum of an incident wave and a reflecting wave

$$
\vec{s} = \vec{s}_i + \vec{s}_r
$$

where

$$
    \vec{s}_i = 
    \begin{pmatrix}
    \frac{1}{c^2} \\
    \frac{k_x}{c} \\ 
    \frac{k_y}{c} \\ 
    1
    \end{pmatrix} \bar{p} e^{-\left( (\frac{k_x(x-x_0) + k_y(y-y_0) - ct)^2}{L^2} \right)}
$$

is the incident wave, and

$$
    \vec{s}_r = 
    \begin{pmatrix}
    \frac{1}{c^2} \\
    -\frac{k_x}{c} \\ 
    \frac{k_y}{c} \\ 
    1
    \end{pmatrix} \bar{p} e^{-\left( \frac{(-k_x(x-(2-x_0)) + k_y(y-y_0) - ct)^2}{L^2} \right)}
$$

is the reflecting wave.

The parameters used in the exact solution are as follows : 

* $\bar{p} = 10^{-4}$ is the amplitude of the sound wave
* $x_0 = y_0 = 0.2$ defines the center of the initial sound wave
* $L = \frac{0.2}{2\sqrt{\ln{2}}}$ is the half-width of the sound wave
* $k_x = k_y = \frac{\sqrt{2}}{2}$ are the $x$ and $y$ components of the wave number
* $c = 1$ is the speed of sound.

The model domain
<figure markdown>
![Plane wave initial condition](./img/planewave_p_init.png){ align=left }
  <figcaption>Pressure field for the initial condition, showing a plane wave with a front oriented at 45 degrees</figcaption>
</figure>

<figure markdown>
![Plane wave during initial reflection](./img/planewave_r_t05.png){ align=left }
  <figcaption>Pressure field at t=0.5 computed with SELF</figcaption>
</figure>

<figure markdown>
![Plane wave reflection later in the simulation](./img/planewave_r_t075.png){ align=left }
  <figcaption>Pressure field at t=0.75 computed with SELF</figcaption>
</figure>

## How we implement this
You can find the example file for this demo in the `examples/linear_euler2d_planewave_propagation.f90` file. This file defines the `lineareuler2d_planewave_model` module in addition to a program that runs the propagating plane wave simulation.


The `lineareuler2d_planewave_model` module defines the `lineareuler2d_planewave` class, which is a type extension of the `lineareuler2d` class that is provided by SELF. We make this type extension so that we can 

* add attributes ( `kx` and `ky` ) for the x and y components of the plane-wave wave number
* add attributes ( `x0` and `y0` ) for the initial center position of the plane-wave
* add an attribute ( `p` ) for the pressure amplitude of the wave
* ad an attribute ( `L` ) for the half-width of the plane wave
* override the `hbc1d_Prescribed` type-bound procedure to set the boundary condition to the exact solution

## Running this example

!!! note
    To run this example, you must first [install SELF](../../GettingStarted/install.md) . We assume that SELF is installed in path referenced by the `SELF_ROOT` environment variable.


To run this example, simply execute

```shell
${SELF_ROOT}/examples/linear_euler2d_planewave_propagation
```

This will run the simulation from $t=0$ to $t=1.0$ and write model output at intervals of $Δ t_{io} = 0.05$.

During the simulation, tecplot (`solution.*.tec`) files are generated which can easily be visualized with paraview.