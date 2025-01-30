# Planetary Rossby Wave
This experiment is designed to show how a geostrophic monopole evolves on a $\beta$-plane. 

## Configuration

### Equations

The equations solved are the linear shallow water equations, given by
$$
    u_t - fv = -g \eta_x
$$
$$
    v_t + fu = -g \eta_y
$$
$$
    \eta_t + (Hu)_x + (Hv)_y = 0
$$

where $\vec{u} =  u \hat{x} + v \hat{y}$ is the barotropic velocity, $g$ is the acceleration of gravity, $H$ is a uniform resting fluid depth, and $\eta$ is the deviation of the fluid free surface relative to the resting fluid. In this model, the $x$ direction is similar to longitude and $y$ is similar to latitude. 

A $\beta$-plane, in geophysical fluid dynamics, is an approximation that accounts for first order variability in the (vertical component of the) coriolis parameter with latitude, 

$$
    f = f_0 + \beta y
$$

The background variation in the planetary vorticity supports Rossby waves, which propagate "westward" with higher potential vorticity to the right of phase propagation.

### Domain Discretization
In this problem, the domain is a square with $(x,y) \in [-500km, 500km]^2$. The model domain is divided into $10\times 10$ elements of uniform size. Within each element, the solution is approximated as a Lagrange interpolating polynomial of degree 7, using the Legendre-Gauss quadrature points as interpolating knots. To exchange momentum and mass fluxes between neighboring elements, we use a local upwind (Lax-Friedrich's) Riemann solver.

### Initial Condition
The initial condition is defined by setting the free surface height to a Gaussian, centered at the origin, with a half width of 10 km and a height of 1 cm.
$$
    \eta(t=0) = 0.01e^{ -( (x^2 + y^2 )/(2.0*10.0^{10}) )}
$$

The initial velocity field is calculated by using the pressure gradient force and using geostrophic balance; in SELF, this is handled by the `LinearShallowWater % DiagnoseGeostrophicVelocity` type bound procedure after setting the initial free surface height.

### Boundary Conditions
Radiation boundary conditions are applied by setting the external state to a motionless fluid with no free surface height variation  ( $u=v=0, \eta = 0$). The model is integrated forward in time using Williamson's $3^{rd}$ order low storage Runge-Kutta, with a time step of $\Delta t = 0.5 s$. 

### Physical Parameters
The remaining parameters for the problem are as follows

* $g = 10 m s^{-2}$
* $f_0 = 10^{-4} s^{-1}$
* $\beta = 10^{-11} m^{-1} s^{-1}$
* $H = 1000 m$


