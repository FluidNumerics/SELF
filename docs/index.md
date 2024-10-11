# Spectral Element Library in Fortran

## SELF
The Spectral Element Library in Fortran (SELF) is an object oriented suite of Fortran modules that enable researchers to implement conservation law solvers using Spectral Element Methods. One goal of this project is to help dispel the notion that Fortran codes are ancient monoliths designed for a single (potentially unknown) purpose. The other goal is to make it easy to explore high order conservation law solvers on a broad range of compute platforms.

SELF is built with modern compute architecture in mind. Our intention is to enable scientific software developers to build SELF for personal workstations, GPU and multi-GPU accelerated platforms, and HPC clusters. SELF can be built with optional GPU acceleration provided by either CUDA & CUBLAS or HIP & HIPBLAS. With GPU-Aware MPI installed, you can also run SELF applications on multi-GPU platforms.


## Models
SELF provides a few pre-canned models for 1-D, 2-D, and 3-D simulations. The specific models based on a generic Discontinuous Galerkin framework that is used to solve a conservation law of the form

\begin{equation}
\vec{s}_t + \nabla \cdot \overleftrightarrow{f} = \vec{q}
\end{equation}

where $\vec{s}$ is a vector of solution variables, $\overleftrightarrow{f}$ is the conservative flux, and $\vec{q}$ are non-conservative source terms. The conservative fluxes are assumed to be functions of the solution variables and their gradients. Effectively, this means that SELF can be used to solve hyperbolic and weakly parabolic partial differential equations. 

Building your own model is done by making a type extension of one of SELF's template classes and overriding the necessary type-bound procedures to fit your needs. Often, it's as simple as defining a `pure function` for the flux, source terms, riemann solver, and boundary conditions. You can use the following template classes to easily build your own models:

* `DGModel1D`
* `DGModel2D`
* `DGModel3D`

Alternatively, you can build on top of one of the pre-canned models:

* `Burgers1D`


[**Learn more about building your own models**](./Tutorials/MakingYourOwnModel.md)


## Support

SELF is made available as a free and open-source software to give you the opportunity to work with the code at your own pace and expense. You can report bugs, issues, and feature requests using the repository's issue tracker. Fluid Numerics, the developers of SELF, will prioritize tasks based on paying customer demand (first) and community demand (second).

Fluid Numerics LLC provides a range of professional services for any of the following

* Implementing new features
* Adding new pre-canned models
* Resolving bugs
* Porting SELF to new hardware
* On call support
* Mathematical Modeling and Numerical Analysis
* SELF usage training

You can reach out for professional support at [support@fluidnumerics.com](mailto:support@fluidnumerics.com)
