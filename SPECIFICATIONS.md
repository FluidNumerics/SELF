# Spectral Element Library in Fortran (SELF)
Copyright 2017-2020 Fluid Numerics LLC

## 1. Introduction

### 1.1 Purpose
The Spectral Element Library in Fortran (SELF) is intended to provide an object-oriented Fortran API for performing numerical differentiation, interpolation, and integration using high order spectral and spectral element methods. This interface will provide access to GPU accelerated routines and methods for handling multi-GPU/multi-core/multi-node acceleration with MPI+, allowing scientific software developers to focus more on their algorithms rather than on the middleware integration with Exascale Era compute hardware. 

SELF will also leverage modern cloud computing infrastructure to support a robust testing infrastructure for each public API call that tests SELF on each commit to develop and main branches.

### 1.2 Intended Audience
The SELF Software Requirements Specification (SSRS) is intended for current and prospective developers and maintainers of SELF and current and prospective developers of applications that leverage SELF (users).

SELF developers and maintainers have core competencies in Fortran, MPI+GPU application development, and Spectral Element theory.

SELF users are expected to come from all levels of experience, from students learning the basics of numerical analysis and discrete ODEs/PDEs to experienced researchers familiar with GPU and multi-GPU computing and applications of Spectral Element Methods.

### 1.3 Intended Use
The SSRS is intended to guide the development and continued maintenance of SELF.

SELF maintainers and developers should use this document to ensure that development activity is kept within scope, meets user needs, and delivers the system features and requirements documented in this SSRS.

SELF users should use this document to hold maintainers and developers accountable on meeting their needs and delivering the system features and requirements specified within this SSRS.

Maintainers, developers, and users should use this document to discuss changes to the scope, user needs, and system features and requirements.

### 1.4 Scope
SELF provides an object-oriented API in Fortran for creating ordinary and partial differential equation solvers using spectral element methods on unstructured grids. All routines provide options for portable GPU acceleration so that end-users can leverage Exascale era compute platforms. 

SELF includes a complete testing infrastructure for each public routine. This infrastructure records errors (numerical + roundoff) and runtime on an array of CPU, GPU, and multi-GPU platforms. SELF tests are incorporated into a continuous integration/continuous benchmarking (CI/CB) platform with a publicly available dashboard for vetting application results. SELF, with the inclusion of its complete testing infrastructure, provides a demonstrably reliable API for numerical model development leveraging spectral element methods.

This SSRS is included with SELF and is meant to evolve through interactions and discussions with users. The SSRS is used to guide developers on user needs and required features and specifications. Users and developers are encouraged to have discussions to help grow SELF and this SSRS for the benefit of the community.

### 1.5 Definitions and Acronyms


## 2. Description
SELF is the culmination of years of prototype and experimental development that has led to this specification of a pattern driven API whose routines are designed to work on single core, multi-core, multi-node, GPU, and multi-GPU platforms.

SELF is originally planned to be a core component of SELF-Fluids, a computational fluid dynamics application based on the Nodal Discontinuous Galerkin Spectral Element method for compressible, incompressible, and hydrostatic primitive equations of motion for fluids. However, the SELF API is written so that it can be leveraged by other developers to implement spectral collocation, continuous galerkin, and discontinuous galerkin spectral element methods based on Legendre-Gauss, Legendre-Gauss-Lobatto, and Legendre-Gauss-Radau quadratures.

### 2.1 User Needs
SELF developers and maintainers assume that users will leverage the SELF API to solve conservation laws. In general, this implies that the coupled systems of partial differential equations being solved can be expressed in terms of finding the time rate of change of some prognostic variables as the divergence of a flux plus source terms. 

SELF users will require methods for calculating divergence, gradient, and curl in one, two, or three spatial dimensions on scalar, vector, and tensor functions. We anticipate that users are interested in performing these operations in complex domains that require unstructured isoparametric elements. SELF users are expected to be able to benefit from built in structured mesh generation and mesh read/write routines. Supported mesh I/O formats will be integrated as developer capacity allows and as the user-base demands.

### 2.2 Assumptions and Dependencies

We assume that end users will have access to mesh generation and pre-processor tools (e.g. Gmsh, HOPR) to specify unstructured meshes when building their applications that incorporate SELF.

To leverage multi-core and multi-GPU platforms, SELF depends on a modern implementation of MPI.
