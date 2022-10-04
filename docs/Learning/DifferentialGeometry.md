# Differential Geometry

## Problem Statement
When working with Spectral Element Methods, we are able to handle complex geometry through the use of an unstructured mesh. Additionally, each element in the mesh is allowed to have curved sides. To make this possible, we define a mapping from physical space to a reference "computational" element. The mapping is then approximated by an interpolant whose knots coincide with those used to represent the solution; this is called an isoparameteric mapping.

In SELF, all of the elements are tensor-product element; in 1-D we work  with line segments, in 2-D we have quadrilaterals, and in 3-D we have hexahedra. The reference space is then

\begin{equation}
    \vec{\xi} \in [-1,1]^D
\end{equation}

where $\vec{\xi} = \sum_{i=1}^D \xi_i \hat{\xi_i}$ and $D$ is the number of spatial dimensions. The mapping from physical space to computational space is defined as

\begin{equation}
    \vec{x}=\vec{x}(\vec{\xi})
\end{equation}

Recall that we are interested in approximating solutions to the conservation law

\begin{equation}
    \mathbf{s}_t + \nabla \cdot \vec{\mathbf{f}} = \mathbf{q}
\end{equation}

The problem we have is that we want to compute the divergence of the conservative flux in physical space, but we only know the coordinate transformation and how to estimate derivatives in computational space. For example, consider calculating $\frac{\partial s}{\partial x}$. Using chain rule, we have

\begin{equation}
    \frac{\partial s}{\partial x} = \sum_{i=1}^D\frac{\partial s}{\partial \xi_i}\frac{\partial \xi_i}{\partial x}
\end{equation}

The problem here is that we know how to compute $\frac{\partial \vec{x}}{\partial \xi_i}$, but not the other way around.

These notes explain how we can calculate derivative operators in physical space given that we have the ability to calculate derivatives in computational space and that we know the coordinate transformation from computational to physical space.

## Div, Grad, and Curl in mapped coordinates
In this section, we derive formulas for computing derivatives, divergence, gradient, and curl under coordinate transformations. We then use the formulas to develop algorithms for implementing spectral element methods under coordinate transformations.

We begin with the premise that we have a logically hexahedral space. The physical positions in this space is represented by a vector with three components

\begin{equation}
    \vec{x} = x_j \hat{x}_j \hspace{3mm} j = 1:3
\end{equation}

We assume that the physical positions can be calculated in terms of “computational” positions, $\vec{\xi}_i \in [-1,1]^3$

\begin{equation}
    \vec{x} = \vec{x}(\vec{\xi})
\end{equation}

![Coordinate System Mapping Schematic](./coordinate-mapping.png){ align=center }

When working in mapped coordinates, we want to understand how to calculate divergence, gradient, and curl of scalar and vector quantities. We are usually given a coordinate transformation $\vec{x} =\vec{x}(\vec{\xi})$ and functions that depend of the computational coordinates $\vec{\xi}$.

As an example, to calculate the derivative of a function $f$, we write

\begin{equation}
    \begin{split}
        \frac{\partial f}{\partial x} &= \frac{\partial f}{\partial \xi_1}\frac{\partial \xi_1}{\partial x} + \frac{\partial f}{\partial \xi_2}\frac{\partial \xi_2}{\partial x} + \frac{\partial f}{\partial \xi_3}\frac{\partial \xi_3}{\partial x} \\
        &= \nabla_\xi f \cdot \frac{\partial \vec{\xi}}{\partial x}
    \end{split}
\end{equation}

From the above equation, you can show that the divergence of a vector is

\begin{equation}
    \nabla \cdot \vec{f} = \sum_{i=1}^{3} \frac{\partial \vec{f}}{\partial \xi_i} \cdot \nabla \xi_i
\end{equation}

where

\begin{equation}
    \nabla \xi_i = \frac{\partial \xi_i}{\partial x_1} \hat{x}_1 + \frac{\partial \xi_i}{\partial x_2} \hat{x}_2 + \frac{\partial \xi_i}{\partial x_3} \hat{x}_3 \label{eq:contravariant_definition}
\end{equation}

Similarly, the gradient of a scalar function is

\begin{equation}
    \nabla f= \sum_{i=1}^{3} \frac{\partial f}{\partial \xi_i} \nabla \xi_i
\end{equation}

and the curl of a vector is

\begin{equation}
    \nabla \times \vec{f}= \sum_{i=1}^{3} \frac{\partial \vec{f}}{\partial \xi_i} \times \nabla \xi_i
\end{equation}

The vectors are called the **contravariant basis vectors** and are usually denoted $\vec{a}^i = \nabla \xi_i$.


The problem with these formulations is that we usually know the transformation from computational space to physical space ( $\vec{x} = \vec{x}(\xi)$ ) and not the other way around. Because of this, it is not readily apparent how to evaluate the contravariant basis vectors. However, we'll show that we can calculate the contravariant basis vectors using quantities we can readily calculate.

Let's first see how we can construct measures of length area and volume. We know that we can relate a small change in $\vec{x}$ to changes in $\vec{\xi}$ through the coordinate transformation.

\begin{equation}
    \begin{split}
        d\vec{x} &= \frac{\partial \vec{x}}{\partial \xi_1} \Delta \xi_1 + \frac{\partial \vec{x}}{\partial \xi_2} \Delta \xi_2 + \frac{\partial \vec{x}}{\partial \xi_3} \Delta \xi_3 \\
        &= \sum_{i=1}^3 \frac{\partial \vec{x}}{\partial \xi_i} \Delta \xi_i
    \end{split}
\end{equation}

The vectors $\frac{\partial \vec{x}}{\partial \xi_i}$ are called the **covariant basis vectors** and are denoted $\vec{a}_i = \frac{\partial \vec{x}}{\partial \xi_i}$. This equation can be used to calculate a measure of arc length

\begin{equation}
    \begin{split}
        ds = | d\vec{x} | &= \sqrt{ \sum_{i=1}^3 \sum_{j=1}^3 \vec{a}_i \cdot \vec{a}_j \Delta \xi_i \Delta \xi_j }\\
        &= \sqrt{ \sum_{i=1}^3 \sum_{j=1}^3 g_{i,j} \Delta \xi_i \Delta \xi_j }
    \end{split}
\end{equation}

The quantities $g_{i,j} = \vec{a}_i \cdot \vec{a}_j$ are the elements of the **covariant metric tensor**.

A small unit area can be calculated by taking the cross product of two units of length.

\begin{equation}
    \begin{split}
        d\vec{A}_i &= \left( \frac{\partial \vec{x}}{\partial \xi_j} \Delta \xi_j \right) \times \left( \frac{\partial \vec{x}}{\partial \xi_k} \Delta \xi_k \right)\\
        &= ( \vec{a}_j \times \vec{a}_k ) \Delta \xi_j \Delta \xi_k \hspace{4mm}  i,j,k \hspace{1mm} cyclic 
    \end{split}
\end{equation}

A differential volume element is calculated by projecting a differential area in the direction of $\vec{a}_i$

\begin{equation}
    \begin{split}
        dV &= d\vec{A}_i \cdot \vec{a}_i \Delta \xi_i \\
        &= (\vec{a}_j \times \vec{a}_k) \cdot \vec{a}_i \Delta \xi_i \Delta \xi_j \Delta \xi_k\\
        &= J \Delta \xi_i \Delta \xi_j \Delta \xi_k \hspace{4mm}  i,j,k \hspace{1mm} cyclic 
    \end{split}
\end{equation}

The quantity $(\vec{a}_j \times \vec{a}_k) \cdot \vec{a}_i$ is called the **Jacobian** of the coordinate transformation.


Let’s now look at how to calculate the divergence. The divergence is defined as

\begin{equation}
    \nabla \cdot \vec{f} = \lim_{\Delta V \rightarrow 0} \frac{1}{\Delta V} \oint_{\partial \Delta V} \vec{f} \cdot d\vec{A}
\end{equation}

The integral on the right hand side is a boundary integral over the faces of the volume element $\Delta V$.

\begin{equation}
    \begin{split}
        \oint_{\partial \Delta V} \vec{f} \cdot d\vec{A} &= \left( \vec{f}\cdot (\vec{a}_2 \times \vec{a}_3)|_{(\xi_1 + \Delta \xi_1,\xi_2,\xi_3)} - \vec{f}\cdot (\vec{a}_2 \times \vec{a}_3)|_{(\xi_1,\xi_2,\xi_3)} \right) \Delta \xi_2 \Delta \xi_3\\
        &+ \left( \vec{f}\cdot (\vec{a}_3 \times \vec{a}_1)|_{(\xi_1,\xi_2 + \Delta \xi_2,\xi_3)} - \vec{f}\cdot (\vec{a}_3 \times \vec{a}_1)|_{(\xi_1,\xi_2,\xi_3)} \right)\Delta \xi_3 \Delta \xi_1 \\
        &+ \left( \vec{f}\cdot (\vec{a}_1 \times \vec{a}_2)|_{(\xi_1,\xi_2,\xi_3 + \Delta \xi_3)} - \vec{f}\cdot (\vec{a}_1 \times \vec{a}_2)|_{(\xi_1,\xi_2,\xi_3)} \right) \Delta \xi_1 \Delta \xi_2
    \end{split}
\end{equation}

With $\Delta V = J \Delta \xi_1 \Delta \xi_2 \Delta \xi_3$, the divergence becomes

\begin{equation}
    \begin{split}
        \nabla \cdot \vec{f} &= \frac{1}{J}\lim_{\Delta V \rightarrow 0} \frac{ \vec{f}\cdot (\vec{a}_2 \times \vec{a}_3)|_{(\xi_1 + \Delta \xi_1,\xi_2,\xi_3)} - \vec{f}\cdot (\vec{a}_2 \times \vec{a}_3)|_{(\xi_1,\xi_2,\xi_3)} }{\Delta \xi_1}\\
        &+ \frac{\vec{f}\cdot (\vec{a}_3 \times \vec{a}_1)|_{(\xi_1,\xi_2 + \Delta \xi_2,\xi_3)} - \vec{f}\cdot (\vec{a}_3 \times \vec{a}_1)|_{(\xi_1,\xi_2,\xi_3)} }{\Delta \xi_2} \\
        &+ \frac{\vec{f}\cdot (\vec{a}_1 \times \vec{a}_2)|_{(\xi_1,\xi_2,\xi_3 + \Delta \xi_3)} - \vec{f}\cdot (\vec{a}_1 \times \vec{a}_2)|_{(\xi_1,\xi_2,\xi_3)}}{\Delta \xi_3}
    \end{split}
\end{equation}

Upon taking the limit as the volume shrinks to zero, we have

\begin{equation}
    \nabla \cdot \vec{f} = \frac{1}{J}\sum_{i=1}^{3}\left[\frac{\partial}{\partial \xi_i}\left( \vec{f} \cdot (\vec{a}_j \times \vec{a}_k) \right)\right] \hspace{4mm}  i,j,k \hspace{1mm} cyclic \label{eq:mapped_divergence_covariant}
\end{equation}


Now, we know that the divergence of a constant vector is equal to zero. This implies that

\begin{equation}
    \nabla \cdot \vec{c} = 0 = \frac{1}{J} \vec{c} \cdot \left[\sum_{i=1}^{3}\frac{\partial}{\partial \xi_i}\left(\vec{a}_j \times \vec{a}_k\right)\right] \hspace{4mm}  i,j,k \hspace{1mm} cyclic \label{eq:constant_divergence}
\end{equation}

In order for the divergence of a constant vector to be valid for any arbitrary constant, it must be that 

\begin{equation}
    \sum_{i=1}^{3}\frac{\partial}{\partial \xi_i}\left( \vec{a}_j \times \vec{a}_k \right)= 0 \hspace{4mm}  i,j,k \hspace{1mm} cyclic \label{eq:metric_identities}
\end{equation}

This equation is known as the **metric identities**. Using the metric identities allows us to write the divergence as

\begin{equation}
    \nabla \cdot \vec{f} = \frac{1}{J}\sum_{i=1}^{3} \frac{\partial\vec{f}}{\partial \xi_i} \cdot (\vec{a}_j \times \vec{a}_k) \hspace{4mm}  i,j,k \hspace{1mm} cyclic \label{eq:mapped_divergence_non_conservative_covariant}
\end{equation}

In comparison with our original formulations for the divergence, we can now relate the contravariant basis vectors to the covariant basis vectors

\begin{equation}
    J\vec{a}^i = \vec{a}_j \times \vec{a}_k \hspace{4mm}  i,j,k \hspace{1mm} cyclic \label{eq:covariant_to_contravariant}
\end{equation}

Similar to the definition of the covariant metric tensor, the **contravariant metric tensor** is define as

\begin{equation}
    g^{i,j} = \vec{a}^i \cdot \vec{a}^j \label{eq:contravariant_metric_tensor}
\end{equation}

With this, we can now define the divergence, gradient, and curl under a coordinate transformation.

Operator | Conservative Form  | Non-Conservative Form |
-------- | ------------------ | --------------------- |
Divergence | \begin{equation} \frac{1}{J}\sum_{i=1}^{3} \frac{\partial}{\partial \xi_i} \cdot (J\vec{a}^i \cdot \vec{f}) \end{equation} | \begin{equation} \sum_{i=1}^{3} \frac{\partial\vec{f}}{\partial \xi_i} \cdot \vec{a}^i \end{equation} |
Gradient | \begin{equation} \frac{1}{J}\sum_{i=1}^{3} \frac{\partial}{\partial \xi_i} \cdot (J\vec{a}^i f) \end{equation} | \begin{equation} \sum_{i=1}^{3} \frac{\partial f}{\partial \xi_i} \vec{a}^i \end{equation} |
Curl | \begin{equation} \frac{1}{J}\sum_{i=1}^{3} \frac{\partial}{\partial \xi_i} \cdot (J\vec{a}^i \times \vec{f}) \end{equation} | \begin{equation} \sum_{i=1}^{3} \frac{\partial\vec{f}}{\partial \xi_i} \times \vec{a}^i \end{equation} |

## Free Stream Preservation
Free stream preservation is a property of a numerical method that refers to the ability of the method to satisfy the metric identities.

With collocation style spectral element methods, we approximate the Jacobian-weighted contravariant basis vectors as a polynomial of degree $N$. When the coordinate mapping is a polynomial of degree $N$, the covariant basis vectors are also a polynomial of degree $N$ and the Jacobian weighted contravariant basis vectors are of degree $2N$. This introduces an aliasing error that causes a lack of satisfaction of the metric identities.

Kopriva (2006) demonstrated the generation of spurious sound waves in an Euler equations solver when the metric identities are not satisfied discretely. Further, he showed that we can use alternate formulations for the contravariant basis vectors that eliminate the aliasing error through basic vector calculus identities.

First, recall that we can write

\begin{equation}
    \nabla u \times \nabla v = - \nabla \times ( v \nabla u ) = \nabla \times( u \nabla v)
\end{equation}

These relations follow from the product rule and the fact that the curl of a gradient is identically zero. From this identity, we can write

\begin{equation}
    \nabla_\xi x_m \times \nabla_\xi x_l = - \nabla \times ( x_l \nabla x_m )
\end{equation}

Note that 
\begin{equation}
\begin{split}
    \nabla_\xi x_m \times \nabla_\xi x_l &= \left( \frac{\partial x_m}{\partial \xi_2}\frac{\partial x_l}{\partial \xi_3} -  \frac{\partial x_m}{\partial \xi_3}\frac{\partial x_l}{\partial \xi_2}\right) \hat{x}_1 \\
    &- \left( \frac{\partial x_m}{\partial \xi_1}\frac{\partial x_l}{\partial \xi_3} - -  \frac{\partial x_m}{\partial \xi_3}\frac{\partial x_l}{\partial \xi_1} \right) \hat{x}_2 \\
    &+ \left( \frac{\partial x_m}{\partial \xi_1}\frac{\partial x_l}{\partial \xi_2} - -  \frac{\partial x_m}{\partial \xi_2}\frac{\partial x_l}{\partial \xi_1} \right) \hat{x}_3
\end{split}
\end{equation}

Further, the contravariant basis vectors can be written

\begin{equation}
\begin{split}
    J\vec{a}^i &= \frac{\partial \vec{x}}{\partial \xi_j} \times \frac{\partial \vec{x}}{\partial \xi_k} \\
    &= \left( \frac{\partial x_2}{\partial \xi_j}\frac{\partial x_3}{\partial \xi_k} -  \frac{\partial x_3}{\partial \xi_j}\frac{\partial x_2}{\partial \xi_k}\right) \hat{x}_1 \\
    &- \left( \frac{\partial x_1}{\partial \xi_j}\frac{\partial x_3}{\partial \xi_k} -  \frac{\partial x_3}{\partial \xi_j}\frac{\partial x_1}{\partial \xi_k}\right) \hat{x}_2 \\
    &+\left( \frac{\partial x_1}{\partial \xi_j}\frac{\partial x_2}{\partial \xi_k} -  \frac{\partial x_2}{\partial \xi_j}\frac{\partial x_1}{\partial \xi_k}\right) \hat{x}_3
\end{split}
\end{equation}

If we let $i=1,2,3$, and $i,j,k$ be cyclic, and let $n=1,2,3$ with $n,m,l$ cyclic we have that

\begin{equation}
    Ja^i_n = \hat{x}_i \cdot ( \nabla_\xi x_m \times \nabla_\xi x_l ) = -\hat{x}_i \cdot ( \nabla_\xi \times ( x_l \nabla_\xi x_m )) = \hat{x}_i \cdot ( \nabla_\xi \times ( x_m \nabla_\xi x_l ))
\end{equation}

For example, when $n=1$, we have that 

\begin{equation}
\begin{split}
    \nabla_\xi x_2 \times \nabla_\xi x_3 &= \left( \frac{\partial x_2}{\partial \xi_2}\frac{\partial x_3}{\partial \xi_3} -  \frac{\partial x_2}{\partial \xi_3}\frac{\partial x_3}{\partial \xi_2}\right) \hat{x}_1\\
    &- \left( \frac{\partial x_2}{\partial \xi_1}\frac{\partial x_3}{\partial \xi_3} - -  \frac{\partial x_2}{\partial \xi_3}\frac{\partial x_3}{\partial \xi_1} \right) \hat{x}_2\\
    &+ \left( \frac{\partial x_2}{\partial \xi_1}\frac{\partial x_3}{\partial \xi_2} - -  \frac{\partial x_2}{\partial \xi_2}\frac{\partial x_3}{\partial \xi_1} \right) \hat{x}_3
\end{split}
\end{equation}

and when $i=1$, we have

\begin{equation}
\begin{split}
    J\vec{a}^1 &= \frac{\partial \vec{x}}{\partial \xi_2} \times \frac{\partial \vec{x}}{\partial \xi_3}\\
    &= \left( \frac{\partial x_2}{\partial \xi_2}\frac{\partial x_3}{\partial \xi_3} -  \frac{\partial x_3}{\partial \xi_2}\frac{\partial x_3}{\partial \xi_3}\right) \hat{x}_1\\
    &- \left( \frac{\partial x_1}{\partial \xi_2}\frac{\partial x_3}{\partial \xi_3} -  \frac{\partial x_3}{\partial \xi_2}\frac{\partial x_1}{\partial \xi_3}\right) \hat{x}_2\\
    &+\left( \frac{\partial x_1}{\partial \xi_2}\frac{\partial x_2}{\partial \xi_3} -  \frac{\partial x_2}{\partial \xi_2}\frac{\partial x_1}{\partial \xi_3}\right) \hat{x}_3
\end{split}
\end{equation}

By comparison, we can see that $Ja^1_1 = \hat{x}_1 \cdot ( \nabla_\xi x_2 \times \nabla_\xi x_3 )$.

Through basic vector calculus, we've shown that we can write the Jacobian weighted contravariant basis vectors as either

\begin{equation}
    Ja^i_n = -\hat{x}_i \cdot ( \nabla_\xi \times ( x_l \nabla_\xi x_m )) 
\end{equation}

or 

\begin{equation}
    Ja^i_n = \hat{x}_i \cdot ( \nabla_\xi \times ( x_m \nabla_\xi x_l ))
\end{equation}

Both formulations are guaranteed to be Divergence free and are known as the "conservative forms". The "curl invariant" form is obtained by averaging the two, to get

\begin{equation}
    Ja^i_n = \frac{1}{2} \hat{x}_i \cdot \nabla_\xi \times [ x_m \nabla_\xi x_l  - x_l \nabla_\xi x_m ]
\end{equation}

When calculating the metric terms in a numerical method, we opt to use the curl form of the metric terms, since we can guarantee that the metric terms will be satisfied.

