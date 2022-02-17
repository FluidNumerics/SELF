# Spectral Element Libraries in Fortran (SELF)
Copyright 2020-2021 Fluid Numerics LLC

[![Documentation Status](https://readthedocs.org/projects/self/badge/?version=latest)](https://self.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/FluidNumerics/SELF/branch/master/graph/badge.svg)](https://codecov.io/gh/FluidNumerics/SELF)


SELF is licensed for use under the [Anti-Capitalist Software License](./LICENSE). For other licensure, reach out to support@fluidnumerics.com.

## Support this project
SELF is open-source made freely available to the public. Fluid Numerics is a nano-business (less than 5 people) that is funded solely through service engagements with customers; we are not funded by venture capital in any way. Developing SELF and the associated documentation and delivering live coding sessions and tutorials all require time & labor as well as compute resources. 

If you value this repository and the activities of the SELF developers, you can [support the continued development and maintenance SELF on Open Collective](https://opencollective.com/higher-order-methods/projects/fluid-self), where we transparently share our operational budget and expenses for this project. Through the [Higher Order Methods collective](https://opencollective.com/higher-order-methods) you can support a variety of activities all focused on providing publicly accessible tools and training for using spectral element methods for physical modeling.

## About
SELF is an object-oriented Fortran library that support the implementation of Spectral Element Methods for solving partial differential equations.

The SELF API is designed based on the assumption that SEM developers and researchers need to be able to implement derivatives in 1-D and divergence, gradient, and curl in 2-D and 3-D on scalar, vector, and tensor functions using spectral collocation, continuous galerkin, and discontinuous galerkin spectral element methods. Additionally, as we enter the exascale era, we are currently faced with a zoo of compute hardware that is available. Because of this, SELF routines provide support for GPU acceleration through AMD's HIP and support for multi-core, multi-node, and multi-GPU platforms with MPI.

See the [Specifications](./SPECIFICATIONS.md) for more details.

To learn more about the SELF API and the software layout, [check out the SELF documentation](https://higherordermethods.github.io/SELF/)
## Getting started
[Join the HigherOrderMethods Discord Server](https://discord.gg/57aNxcpYMW) to get community support.

### Dependencies
Before installing SELF, make sure your system has the following dependencies installed : 
* Fortran compiler (e.g. gfortran)
* [HIP](https://github.com/ROCm-Developer-Tools/HIP) (For GPU Accelerated Builds)
* [HIPFort](https://github.com/ROCmSoftwarePlatform/hipfort) (For GPU Accelerated Builds)
* [stdlib](https://github.com/fortran-lang/stdlib)
* [json-fortran](https://github.com/jacobwilliams/json-fortran)
* [FLAP](https://github.com/szaghi/FLAP)
* [feq-parse](https://github.com/FluidNumerics/feqparse)
* [HDF5]()

For GPU accelerated builds, we recommend installing `rocm-dev` on either CentOS or Ubuntu OS following the [ROCm installation guide](https://rocmdocs.amd.com/en/latest/Installation_Guide/Installation-Guide.html). You should install HIPfort from source, using the Fortran compiler you will build SELF with.

### Building SELF

**Serial CPU Builds**
```
BUILD=release \
FC=gfortran \
make
```
The SELF make system assumes that `feqparse`, `FLAP`, and `HDF5` are installed under `/opt/view`. If these dependencies are installed elseswhere you can use the following environment variables specify the linker and includes flags.
```
     SELF_FEQPARSE_LIBS     Set the linker flags for feq-parse (Default: -L/opt/view/lib -lfeqparse)
     SELF_FEQPARSE_INC      Set the includes flags for feq-parse (Default: -I/opt/view/include)
     SELF_FLAP_LIBS         Set the linker flags for FLAP (Default: -L/opt/view/lib/ -lFLAP -lFACE -lPENF) 
     SELF_FLAP_INC          Set the includes flags for FLAP (Default: -I/opt/view/include/FLAP -I/opt/view/include/PENF -I/opt/view/include/FACE)
     SELF_HDF5_LIBS         Set the linker flags for HDF5 (Default: -L/opt/view/lib/ -lhdf5_fortran -lhdf5 -lz -lm) 
     SELF_HDF5_INC          Set the includes flags for HDF5 (Default: -I/opt/view/include/shared)
```
By default, self will install under `/opt/self`. To change the install path, use the `SELF_PREFIX` environment variable.

Upon installation, you will have
* `${SELF_PREFIX}/lib/libself.a` : Static library for linking into other application
* `${SELF_PREFIX}/include/*.mod` : Fortran module files required for building other applications that depend on SELF
* `${SELF_PREFIX}/bin/self` : A test program that can be used to test each SELF routine
* `${SELF_PREFIX}/util` : A directory with useful utilities for pre and post-processing data produced by SELF

## Maintainers
* [Joseph Schoonover, Fluid Numerics LLC](https://fluidnumerics.com/people/joe-schoonover)

Want to become a maintainer ? Reach out to support@fluidnumerics.com

If you'd like to contribute, see [CONTRIBUTING.md](./CONTRIBUTING.md) to get started.

If you need help, [open an issue](https://github.com/FluidNumerics/SELF/issues/new)

