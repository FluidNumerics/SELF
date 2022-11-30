# Skew Symmetric Formulations for the Shallow Water Equations

The nonlinear shallow water equations can be written as


When approximated with spectral element methods, the nonlinear terms in the flux function introduce an aliasing error that stems from an inability to satisfy the chain rule. In the shallow water equations, this error unfortunately generates spurious flows in the presence of variable bottom topography, even in one dimension.

For example, consider the following initial conditions 

\begin{subequations}
  \begin{align}
     Hu &= 0 \\
     H &= h(x) = 1-0.2e^{-\frac{(x-0.5)^2}{0.01}}
  \end{align}
\end{subequations}

This initial condition sets the initial velocity to null, free surface height to null, and the bottom topography to be a Gaussian hill of $0.2$. For the continous form of the 1-D shallow water equations, the exact solution is identical to the initial condition (the "lake at rest" solution).


If we integrate the model forward using a standard nodal discontinuous galerkin formulation by a single time step, we can see the change in the velocity field immediately.


