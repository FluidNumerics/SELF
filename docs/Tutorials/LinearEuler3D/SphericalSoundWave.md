# Linear Euler 2D - Spherical Sound Wave Example
This tutorial will walk you through using an example program that uses the `LinearEuler2D` class to run a simulation with the linear Euler equations for an ideal gas in 2-D. This example is configured using the built in structured mesh generator with no-normal-flow boundary conditions on all domain boundaries.

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
The physical domain is defined by $\vec{x} \in [0, 1]\times[0,1]$. We use the `StructuredMesh` routine to create a domain with 20 × 20 elements that are dimensioned 0.05 × 0.05 . All model boundaries are all tagged with the `SELF_BC_NONORMALFLOW` flag to implement no-normal-flow boundary conditions.

Within each element, all variables are approximated by a Lagrange interpolating polynomial of degree 7. The interpolation knots are the Legendre-Gauss points.


### Initial Conditions
The initial condition is set using


$$
    \begin{pmatrix}
    ρ \\ 
    u \\ 
    v \\ 
    P
    \end{pmatrix} = 
    \begin{pmatrix}
    \frac{1}{c^2} \\
    0 \\ 
    0 \\ 
    1
    \end{pmatrix} \bar{p} e^{-\left( \frac{ (x-x_0)^2 + (y-y_0)^2 }{L^2} \right)}
$$


The parameters used in the exact solution are as follows : 

* $\bar{p} = 10^{-4}$ is the amplitude of the sound wave
* $x_0 = y_0 = 0.5$ defines the center of the initial sound wave
* $L = \frac{0.06}{\sqrt{\ln{2}}}$ is the half-width of the sound wave
* $c = 1$ is the speed of sound.

This initial condition is similar to the spherical sound wave on pg. 218 of [Kopriva (2009), "Implementing Spectral Methods for Partial Differential Equations"](https://link.springer.com/book/10.1007/978-90-481-2261-5)


<figure markdown>
![Spherical sound-wave initial condition](./img/sphericalwave_r_init.png){ align=left }
  <figcaption>Pressure field for the initial condition, showing a spherical sound wave centered in the domain. </figcaption>
</figure>

<figure markdown>
![Spherical sound-wave halfway through the simulation](./img/sphericalwave_r_t05.png){ align=left }
  <figcaption>Pressure field at t=0.5 computed with SELF</figcaption>
</figure>

<figure markdown>
![Spherical sound-wave at the end of the simulation](./img/sphericalwave_r_t1.png){ align=left }
  <figcaption>Pressure field at t=1 computed with SELF</figcaption>
</figure>

## How we implement this
You can find the example file for this demo in the `examples/linear_euler2d_spherical_soundwave_closeddomain.f90` file. This examples using the `LinearEuler2D` class as-is to write a simple program that simulates the expansion of a spherical sound-wave in a closed domain.

## Running this example

!!! note
    To run this example, you must first [install SELF](../../GettingStarted/install.md) . We assume that SELF is installed in path referenced by the `SELF_ROOT` environment variable.


To run this example, simply execute

```shell
${SELF_ROOT}/examples/linear_euler2d_spherical_soundwave_closeddomain
```

This will run the simulation from $t=0$ to $t=1.0$ and write model output at intervals of $Δ t_{io} = 0.1$.

During the simulation, tecplot (`solution.*.tec`) files are generated which can easily be visualized with paraview.