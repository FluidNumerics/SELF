# Spectral Element Libraries in Fortran (SELF)
Copyright 2020 Fluid Numerics LLC

SELF is licensed for use under the [Anti-Capitalist Software License](./LICENSE). For other licensure, reach out to support@fluidnumerics.com.

## About
SELF is an object-oriented Fortran library that support the implementation of Spectral Element Methods for solving partial differential equations.

The SELF API is designed based on the assumption that SEM developers and researchers need to be able to implement derivatives in 1-D and divergence, gradient, and curl in 2-D and 3-D on scalar, vector, and tensor functions using spectral collocation, continuous galerkin, and discontinuous galerkin spectral element methods. Additionally, as we enter the exascale era, we are currently faced with a zoo of compute hardware that is available. Because of this, SELF routines provide support for GPU acceleration through AMD's HIP and support for multi-core, multi-node, and multi-GPU platforms with MPI.

See the [Specifications](./SPECIFICATIONS.md) for more details.

Currently SELF is being developed from a refactoring of SELF-Fluids.

## Getting started
### Dependencies
Before installing SELF, make sure your system has the following dependencies installed : 
* Fortran compiler (e.g. gfortran)
* [HIP](https://github.com/ROCm-Developer-Tools/HIP) (For GPU Accelerated Builds)
* [HIPFort](https://github.com/ROCmSoftwarePlatform/hipfort) (For GPU Accelerated Builds)
* [FLAP](https://github.com/szaghi/FLAP)
* [feqparse](https://github.com/FluidNumerics/feqparse)


### Building SELF

**Serial CPU Builds**
```
BUILD=release \
FC=gfortran \
make
```
The SELF make system assumes that `feqparse` and `FLAP` are installed under `/opt/feqparse` and `/opt/FLAP` respectively. If these dependencies are installed elseswhere you can use the following environment variables specify the linker and includes flags.
```
     SELF_FEQPARSE_LIBS     Set the linker flags for feq-parse (Default: -L/opt/feqparse/lib -lfeqparse)
     SELF_FEQPARSE_INC      Set the includes flags for feq-parse (Default: -I/opt/feqparse/include)
     SELF_FLAP_LIBS         Set the linker flags for FLAP (Default: -L/opt/FLAP/lib/ -lFLAP -lFACE -lPENF) 
     SELF_FLAP_INC          Set the includes flags for FLAP (Default: -I/opt/FLAP/include/FLAP -I/opt/FLAP/include/PENF -I/opt/self/include/FACE)
```
By default, self will install under `/opt/self`. To change the install path, use the `SELF_PREFIX` environment variable.

Upon installation, you will have
* `${SELF_PREFIX}/lib/libself.a` : Static library for linking into other application
* `${SELF_PREFIX}/include/*.mod` : Fortran module files required for building other applications that depend on SELF
* `${SELF_PREFIX}/bin/self` : A test program that can be used to test each SELF routine
* `${SELF_PREFIX}/bin/ci.sh` : A shell script to run all of the CI tests for SELF

## Maintainers
* [Joseph Schoonover, Fluid Numerics LLC](https://fluidnumerics.com/people/joe-schoonover)

Want to become a maintainer ? Reach out to support@fluidnumerics.com

If you'd like to contribute, see [CONTRIBUTING.md](./CONTRIBUTING.md) to get started.

If you need help, [open an issue](https://github.com/FluidNumerics/SELF/issues/new)

