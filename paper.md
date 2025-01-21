---
title: 'Spectral Element Library in Fortran : A portable, spectrally accurate, multi-GPU accelerated conservation law solver written in Fortran'
tags:
  - Fortran
  - Computational Fluid Dynamics
  - Conservation Laws
  - High performance computing
authors:
  - name: Joseph A. Schoonover
    orcid: 0000-0001-5650-7095
    equal-contrib: true
    affiliation: 1
  - name: Garrett Byrd
    equal-contrib: true
    affiliation: 1
  - name: Siddhartha Bishnu
  
affiliations:
 - name: Fluid Numerics, United States
   index: 1
date: 13 August 2025
bibliography: paper.bib

---

# Summary
Hyperbolic and parabolic partial differential equations can be cast as a generic conservation law, defined by an array of conservative fluxes and non-conservative source terms. A variety of physical phenomena can be modeled using such conservation laws, which allows for the development of a software framework that can be used to solve a variety of conservation laws. This allows for the development of a flexible framework that can be used to define partial differential equation solvers to model a variety of physical phenomena. The Spectral Element Library in Fortran (SELF) is an object-oriented implementation of the spectral element method applied to generic conservation laws. SELF consists of core algorithms for computing weak and strong form differential operations on unstructured isoparametric grids and higher level abstract models. The abstract models provide a suite of explicit Runge-Kutta time integrators and a standard forward stepping pipeline. End users of SELF can create their own solvers by defining a minimal set of pure fortran functions, which can reduce the time-to-science. All methods included in SELF are implemented on CPU and GPU platforms. Domain decomposition is built in to SELF, which allows for easy scaling, simply by running SELF programs with multiple MPI ranks. GPU acceleration is supported by HIP or CUDA, which permits portability across venddor platforms.


# Statement of need

Simulating complex physical phenomena often requires solving partial differential equations (PDEs) cast as conservation laws. These equations underpin a variety of scientific and engineering disciplines, including fluid dynamics, physical oceanography, magnetohydrodynamics, and aerospace engineering. Despite the widespread use of conservation laws in modeling, researchers and practitioners face significant barriers in developing efficient and scalable solvers tailored to their specific applications. These barriers include time-intensive algorithm development, challenges in implementing scalable parallel solvers, and the complexity of targeting both CPU and GPU platforms.

The Spectral Element Library in Fortran (SELF) addresses these challenges by providing a flexible, modular, and performant framework for solving generic conservation laws using spectral element methods. SELF is designed to reduce the "time-to-science" by enabling users to quickly prototype and implement PDE solvers using a minimal set of pure Fortran functions. Its object-oriented architecture supports a wide range of applications by providing:

1. Core algorithms for differential operations on unstructured isoparametric grids, ensuring geometric flexibility and high-order accuracy.
2. Abstract models that include explicit Runge-Kutta time integrators and a standardized forward-stepping pipeline.
3. Scalable parallelism through built-in support for domain decomposition and MPI.
4. Cross-platform GPU acceleration via HIP or CUDA, offering portability across hardware from different vendors.

SELF is particularly valuable for researchers and developers who require highly efficient solvers that can be scaled to modern heterogeneous computing systems. It caters to a diverse user base, from academic researchers developing novel numerical methods to engineers and domain scientists modeling complex real-world phenomena. By abstracting low-level implementation details and providing robust, reusable components, SELF enables users to focus on scientific discovery rather than software development.

SELF’s explicit support for both CPUs and GPUs, its portability across vendor platforms, and its ability to seamlessly scale from personal computers to high-performance computing clusters make it an indispensable tool for modern computational science. It represents a significant advancement in the development of open-source, high-performance tools for solving PDEs.

# Examples
SELF includes example models to demonstrate its usage :

* Burgers Equation in 1-D
* Linear Shallow Water in 2-D
* Linear Euler Equations in 2-D and 3-D


Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge Jorge Galvez-Vallejo and Jonathan Moore for their support of the Spectral Element Library in Fortran.

# References