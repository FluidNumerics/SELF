# Compressible Navier-Stokes


## Hydrostatic Adjustment

Modeling compressible fluids in the presence of potential forces (such as gravity) often requires defining initial conditions that have already adjusted to the potential force. For some simple configurations, such as those with constant gravitational acceleration, it is easy to write down the fluid state for a hydrostatic compressible fluid. More complicated gravitational potentials pose a challenge. 

SELF's Compressible Ideal Gas modules come equipped with the `HydrostaticAdjustment` method, which can be used to compute a hydrostatic fluid state, given the potential function and an initial density and energy field. This method currently works by forward-stepping the equations of motion with an artifical momentum drag term until the fluid momentum reaches a specified tolerance. Ideally, we would solve this system using an implicit time stepping scheme. Because SELF currently only provides explicit time stepping schemes, we brute force our way to equilibrium, The method chooses a time step so that the maximum CFL number based on the maximum initial sound wave speed is 0.75. 


### Choosing the artificial momentum drag
To explain how we choose the artificial momentum drag coefficient, consider the compressible euler equations, without the momentum advection terms and with the additional momentum drag

\begin{align}
 (\rho \vec{u})_t &= -\nabla p - C_d \rho \vec{u}  \\
 \rho_t + \nabla \cdot ( \rho \vec{u} ) &= 0
\end{align}

Taking the divergence of the momentum equation and substituting into the time derivative of the mass conservation equation gives

\begin{equation}
\rho_{tt} - \nabla^2 p - C_d \nabla (\rho \vec{u}) = 0 
\end{equation}

which we can also write as

\begin{equation}
\rho_{tt} - \nabla^2 p = - C_d \rho_t 
\end{equation}


Using the equation of state, we can write the density as a function of the pressure $\rho = \rho(p)$. By applying the chain rule, we can write

\begin{equation}
\rho_t = \rho_p p_{tt}
\end{equation}

where we have made the assumption that $(\rho_p)_t = 0$. Using these asssumptions, we can write a single equation for the pressure

\begin{equation}
p_{tt} - c^2 \nabla^2 p = - C_d p_t
\end{equation}

where $c = (\rho_p)^{-1/2}$ is the speed of sound. Using Fourier solutions for the pressure 

\begin{equation}
p = \hat{p} e^{i(\vec{k}\cdot\vec{x} - \sigma t)}
\end{equation}

we obtain the following dispersion relation

\begin{equation}
\sigma^2 + i C_d \sigma - c^2 | \vec{k} |^2 = 0
\end{equation}

which has roots

\begin{equation}
\sigma = \frac{1}{2}( -i C_d \pm \sqrt{ 4 c^2 |\vec{k}|^2 - C_d^2} )
\end{equation}

The frequency becomes purely complex, exhibiting no oscillatory motions, when

\begin{equation}
C_d \geq 2 c ||\vec{k}||
\end{equation}

For a numerical method, the largest that the right hand side can be is when $\vec{k}$ is associated with the shortest resolvable wave mode; usually, this is the Nyquist mode, which has a wavelength of $2\Delta x$, so that

\begin{equation}
C_d \geq \frac{c}{\Delta x}
\end{equation}

From this brief analysis, we  gain some insight into how the momentum drag can influence the evolution of sound waves.In adjusting the fluid to a hydrostatic state, imbalances in the potential forces and the pressure gradient force lead to an erruption of sound waves. The momentum drag acts to damp out these disturbances and we can choose the drag coefficient so that a wide range of wavelengths are damped, with no oscillation. We found a condition where all modes are damped; this corresponds to making the momentum drag as important as sound wave propagation at the grid scale. So long as the model is stepped forward in a stable manner, we can integrate the compressible equations until the fluid momentum is near zero.
