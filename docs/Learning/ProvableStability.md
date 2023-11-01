# Provable stability

## Motivation
When we develop numerical methods, we want to know if the method we are going to spend time on is supposed to be **stable**. Numerical methods for conservation laws will iteratively apply the same operations over and over. Successive iterations form a sequence of values that we ascribe physical meaning to through a mathematical model. It's important that the generated sequence does not "blow up" and the iterations don't result in arbitrarily large values. In some sense, this is what many of us mean when we say a method is unstable.

Hyperbolic systems, like one-way advection, or the two-way wave equation are models for information transmission between locations. They do not amplify or deamplify the original amount of energy stored in a system; instead, they move it around. Some hyperbolic systems, like variable advection and the nonlinear "Burgers Equation", can develop discontinuities, independent of the structure of the initial condition. For nonlinear hyperbolic systems, there can be a number of possible solutions to the PDE when discontinuities are present. Physically, this is often interpreted as a breakdown of the assumptions of the original model and we know that some other dynamics needs to take over to "smooth out" the discontinuities. In the cases where smoothing occurs, the model incorporates some mechanism to dissipate energy that accumulates in the vicinity of a discontinuity.

With this kind of thinking, we have the idea that the "energy" associated with a hyperbolic conservation law ought to be bounded by the initial amount of energy that we started with. When the solution does not develop discontinuites, the total energy remains constant. When discontinuities do arise, the total energy decreases. To formulate this mathematically, we start with the definition of our conservation law, which is generally,

\begin{equation}
  \mathbf{s}_t + \vec{\nabla} \cdot \vec{\mathbf{f}} = 0 
\end{equation}

where $\mathbf{s}$ is a state vector containing the prognostic variables we are modeling and $\vec{\mathbf{f}}$ is a "block vector" containing the physical flux vector for each prognostic variable. We've explicitly ignored non-conservative source terms, for the sake of the discussion here.


## Stability in terms of Mathematical Entropy
In spectral element literature, researchers will focus on stability criteria with respect to **convex (mathematical) entropy** functions. The term "entropy" is used because its use is motivated by thermodynamic entropy in the context of computational fluid dynamics. An entropy function is a nonlinear scalar function of the prognostic variables

\begin{equation}
  q = q(\mathbf{s})
\end{equation}

The entropy function is called convex if

\begin{equation}
  \mathbf{v}^T \frac{\partial^2 q}{\partial \mathbf{s}^2}\mathbf{v} > 0, \hspace{2mm} \forall \mathbf{v} \neq 0
\end{equation}


We can create another conservation law for the entropy function that is consistent with the original conservation law. First, we define the "entropy variables" as the derivatives of the entropy function with respect to each prognostic variable

\begin{equation}
  \mathbf{w} = \begin{pmatrix}
                \frac{\partial q}{\partial s_1} \\
                \frac{\partial q}{\partial s_2} \\
                \frac{\partial q}{\partial s_3} \\
                 ... \\
                \frac{\partial q}{\partial s_M} \\
                \end{pmatrix}
\end{equation}

With this defintion in mind, the time derivative of the entropy function can be written as follows :
\begin{equation}
q_t = \sum_{m} \frac{\partial q}{\partial s_m} \frac{\partial s_m}{\partial t} = \mathbf{w}^T \mathbf{s}_t
\end{equation} 

Performing a dot-product of the conservation law with the entropy variables gives

\begin{equation}
  \mathbf{w}^T \mathbf{s}_t + \mathbf{w}^T \vec{\nabla} \cdot \vec{\mathbf{f}} = q_t + \mathbf{w}^T \vec{\nabla} \cdot \vec{\mathbf{f}} = 0
\end{equation}

Approaching this from another perspective, we want the entropy to be conservative so that the entropy magnitude is bounded by its initial value and the boundary contributions. The entropy function obeys a conservation law of the form

\begin{equation}
q_t + \vec{\nabla} \cdot \vec{\Psi} = 0
\end{equation}

where $\vec{\Psi}$ is the "entropy flux". The entropy function is usually defined in such a way that

\begin{equation}
 \mathbf{w}^T \vec{\nabla} \cdot \vec{\mathbf{f}} = \vec{\nabla}\cdot\vec{\Psi}
\end{equation}

This contraction is important, as it allows us to treat the entropy as a conserved quantity. Integration over space and time for the entropy then results in

\begin{equation}
  \int_{V} q(T) \hspace{1mm} dV = \int_{V} q(0) \hspace{1mm} dV + \int_{0}^T \oint_{\partial V} \vec{\Psi} \cdot \hat{n} \hspace{1mm} dS 
\end{equation}

which indicates that the total entropy at a later time is equal to the initial total entropy plus any additional entropy brought in through boundary conditions. A closed system is one in which no entropy is introduced at the boundaries, so that the total entropy at a later time is equivalent to the total entropy at the initial time. As we alluded to previously, when discontinuities develop in the solution, we expect that smoothing of the solution is necessary in such a way that the entropy is reduced (our previous discussion used "energy" as a hopefully more familiar analogy). 

To summarize, we express stability as 

\begin{equation}
  \int_{V} q(T) \hspace{1mm} dV \leq \int_{V} q(0) \hspace{1mm} dV + \int_{0}^T \oint_{\partial V} \vec{\Psi} \cdot \hat{n} \hspace{1mm} dS 
\end{equation}

When developing a provably stable numerical methods, we need to do the following :

1. Define a nonlinear convex entropy function that satisfies the entropy flux contraction $\mathbf{w}^T \vec{\nabla} \cdot \vec{\mathbf{f}} = \vec{\nabla}\cdot\vec{\Psi} \rightarrow \vec{\mathbf{f}}^T \cdot \vec{\nabla} \mathbf{w} = 0$
2. Define a discretization for the original conservation law that discretely satisfies entropy conservation


### Example : 1-D Advection
To help conceptualize this, let's consider a simple 1-D advection problem,

\begin{equation}
  s_t + \frac{\partial}{\partial x}(u s)  = 0
\end{equation}

where $s$ is the solution and $u$ is a velocity field. A convex entropy function here could be 

\begin{equation}
  q(s) = s^2
\end{equation}


## Split form equations


## Two-point flux

