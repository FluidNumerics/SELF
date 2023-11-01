# Software Architecture

The Spectral Element Library in Fortran (SELF) uses a pattern driven design that starts with a few core classes and bootstraps its way to discrete calculus operations in complex geometry and workflows for PDE solvers. In this section of the documentation, you will learn about 

* The conceptual design of SELF
* The project scaffolding

## Universal constants
Every universe has constants, so why shouldn't SELF ? 

The `SELF_Constants` module provides definitions for floating point precision across the entirety of the SELF. Additionally, there are some other useful constants, like $\pi$, machine precision, and various `CHARACTER` lengths. Floating point precision is set using a C-preprocessing conditional. If the CPP flag `DOUBLE_PRECISION` is defined at build time, then `prec=real64` (from the `ISO_FORTRAN_ENV` module); alternatively `prec=real32`. Throughout the rest of SELF, `REAL` types are declared with `prec` precision, e.g. `REAL(prec)`.

## Memory Management
SELF is designed to provide a simplified interface for managing memory on GPU accelerated systems. Because of this, the bottom end of SELF is a set of classes that manage arrays of 32-bit integers, 64-bit integers, and floating point values (floating point precision is determined at build-time). Classes in `SELF_Memory` define structure for 1-D through 7-D arrays for each type, giving 21 classes. 

Each class is named according to the underlying data-type and the rank of the underlying array. The convention we use is

```
hfTYPE_rN
```

where `TYPE` is `int32`, `int64`, or `real` and `N` is the rank and varies between one and 7. The examples below should help you better understand the naming convention

* `hfInt32_r3` - A structure for working with a rank 3 arrays of 32-bit integers on CPU and GPU
* `hfReal_r7` - A structure for working with rank 7 arrays of floating point values on CPU and GPU
* `hfInt64_r1` - A structure for working with rank 1 arrays of 64-bit floating point values on CPU and GPU

Each class has a `hostData` and `deviceData` attribute. Data within these classes that are stored on the CPU side (`hostData`) are Fortran `POINTER`s, while data stored on the GPU side (`deviceData`) side are `TYPE(c_ptr)`s from the `ISO_C_BINDING` module. Memory allocation and deallocation on the host are handled by Fortran's intrinsic `ALLOCATE` and `DEALLOCATE` methods. On the device, GPU memory is allocated and deallocated using `hipMalloc` and `hipFree`, respectively, from AMD's HIP/HIPFort. 

From the user's perspective, you can work with these classes using their type-bound procedures for

* Memory allocation
* Memory deallocation
* Memory copy from host to device
* Memory copy from device to host

Accessing data is fairly straightforward using Fortran syntax for access class attributes.


```
USE SELF_Memory

TYPE(hfReal_r4) :: mydata

CALL mydata % Alloc(loBound=(/1,1,1,1/),&
                    upBound=(/10,10,10,10/))

mydata % hostData = 1.0_prec

CALL myData % UpdateDevice()

CALL myData % Free()
```

## Quadrature
Polynomial-based spectral and spectral element methods use discrete quadrature to approximate integrals. In nodal spectral element methods, like those implemented in SELF, functions are represented as interpolants whose interpolation knots are chosen to match specific quadrature points. The quadrature points and weights are chosen to preserve discrete orthogonality for specific polynomial basis, like Chebyshev or Legendre polynomials.

The `SELF_Quadrature` module provides routines for generating Gauss, Gauss-Lobatto, and Gauss-Radau quadrature for both Chebyshev and Legendre basis. This module is primarily used by the `SELF_Lagrange` module for building Lagrange interpolating polynomials that exhibit spectral accuracy.

## Lagrange Polynomials



## Data Types

### Scalars

### Vectors

### Two-point vectors
Two-point vectors are used in the construction of Entropy-Conservative nodal DGSEM (ECDGSEM).

Using notation similar to [Winters et al. 2020], we write the flux vector in a conservation law as a block vector : 

\begin{equation}
  \overleftrightarrow{\mathbf{F}} = \mathbf{F}^1 \hat{\xi_1} + \mathbf{F}^2 \hat{\xi_2} + \mathbf{F}^3 \hat{\xi_3}
\end{equation}

where $\hat{\xi_i}$ is the $i^{th}$ computational coordinate direction and $\mathbf{F}^i$ is a state vector for the flux in the $i^{th}$ computational direction. When a two-point vector is used to represent the flux, the discrete divergence operation can be written

\begin{equation}
 \nabla \cdot \overleftrightarrow{\mathbf{F}} = 2 \sum_{n=0}^{N} \left( \mathbf{F}^1_{(i,n),j,k} D_{i,n} + \mathbf{F}^2_{i,(j,n),k} D_{j,n} + \mathbf{F}^3_{i,j,(k,n)} D_{k,n} \right)
\end{equation}

The details of how the two-point flux are calculated is best reserved for discussions in the development of specific conservation law solvers. However, to illustrate how they are used in practice, we consider a linear flux

\begin{equation}
  \overleftrightarrow{\mathbf{F}} = A^1 \mathbf{s} \hat{\xi_1} + A^2 \mathbf{s} \hat{\xi_2} + A^3 \mathbf{s} \hat{\xi_3}
\end{equation}

where the $A^i$ are scalar functions that do not depend on $\mathbf{s}$. Often the two-point fluxes are defined so that they provide equivalent discretizations that would be obtained if one started from the split form of the PDE. The two-point flux in this case is taken as the product of averages of the $A^i$ and $\vec{s}$, e.g.

\begin{equation}
\mathbf{F}^1_{(i,n),j,k} = \left( \frac{A^1_{i,j,k} + A^1_{n,j,k}}{2} \right) \left(\frac{\mathbf{s}_{i,j,k} + \mathbf{s}_{n,j,k}}{2} \right)
\end{equation}

\begin{equation}
\mathbf{F}^2_{i,(j,n),k} = \left( \frac{A^2_{i,j,k} + A^2_{i,n,k}}{2} \right) \left(\frac{\mathbf{s}_{i,j,k} + \mathbf{s}_{i,n,k}}{2} \right)
\end{equation}

\begin{equation}
\mathbf{F}^3_{i,j,(k,n)} = \left( \frac{A^3_{i,j,k} + A^3_{i,j,n}}{2} \right) \left(\frac{\mathbf{s}_{i,j,k} + \mathbf{s}_{i,j,n}}{2} \right)
\end{equation}

This has been shown to be an equivalent discretization for the split form of the advective operator.

As a practical aspect, we see that the two-point flux adds a dimension in the computational coordinates, relative to the standard "vectors" implemented in SELF. Because of this, we use the following memory layout for two-point vectors :

In 2-D 
```fortran
F(1:2,0:N,0:N,0:N,1:nVar,1:nEl)
```

and in 3-D
```fortran
F(1:3,0:N,0:N,0:N,0:N,1:nVar,1:nEl)
```

The first dimension refers to the directional component (physical or computational) of the vector. The second dimension is a computational coordinate direction and variations in this dimensions are equivalent to looping over $n$ in the two-point vector representation shown above. The remaining array dimensions of size `0:N` are computational coordinate directions and are the usual `(i,j,k)` dimensions, in this order. This is followed by the dimension over the solution variables and then the elements. 


### Tensors
