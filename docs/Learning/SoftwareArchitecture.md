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


