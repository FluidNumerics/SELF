# Spectral Element Libraries in Fortran (SELF)
Copyright 2020-2023 Fluid Numerics LLC

[![codecov](https://codecov.io/gh/FluidNumerics/SELF/branch/main/graph/badge.svg?token=AKKSL5CWK6)](https://codecov.io/gh/FluidNumerics/SELF)

[![linux-gnu-cmake](https://github.com/FluidNumerics/SELF/actions/workflows/linux-gnu-cmake.yml/badge.svg)](https://github.com/FluidNumerics/SELF/actions/workflows/linux-gnu-cmake.yml)

[![linux-gnu-multithreaded-cmake](https://github.com/FluidNumerics/SELF/actions/workflows/linux-gnu-multithreaded-cmake.yml/badge.svg)](https://github.com/FluidNumerics/SELF/actions/workflows/linux-gnu-multithreaded-cmake.yml)

[![linux-amdflang-cmake](https://github.com/FluidNumerics/SELF/actions/workflows/linux-amdflang-cmake.yaml/badge.svg)](https://github.com/FluidNumerics/SELF/actions/workflows/linux-amdflang-cmake.yaml)

## Licensing
SELF is licensed for use under a [3 Clause BSD with attribution license](./LICENSE). [Fluid Numerics](https://www.fluidnumerics.com) is a small family-owned business. We want to make SELF available to folks who want to build conservation laws that run on a wide range of platforms. Under the license, you can use, modify, and redistribute SELF so long as attribution is given to Fluid Numerics. 

## How to support this repository
Continued support of SELF relies on users and customers. Here's a few ways you can support this project:

* Give this repository a star on Github
* [Open an issue](https://github.com/FluidNumerics/SELF/issues/new/choose) if you have a question or want to report a problem. We want to help!
* [Subscribe to Fluid Numerics on Youtube](https://www.youtube.com/@FluidNumerics?sub_confirmation=1)
* [Sponsor this project on Open Collective](https://opencollective.com/opensource-fluidnumerics)
* [Work with us](https://www.fluidnumerics.com/services)

## About
SELF is an object-oriented Fortran library that support the implementation of Spectral Element Methods for solving partial differential equations.

The SELF API is designed based on the assumption that SEM developers and researchers need to be able to implement derivatives in 1-D and divergence, gradient, and curl in 2-D and 3-D on scalar, vector, and tensor functions using spectral collocation, continuous galerkin, and discontinuous galerkin spectral element methods. Additionally, as we enter the exascale era, we are currently faced with a zoo of compute hardware that is available. Because of this, SELF routines provide support for GPU acceleration through AMD's HIP and support for multi-core, multi-node, and multi-GPU platforms with MPI.


### Prerequisites
All you need is a Fortran compiler that is compliant with the Fortran 2008 standard and supports C interoperability. You can see which compilers are regularly tested on the [Github actions page](https://github.com/FluidNumerics/feq-parse/actions/workflows/ci.yml). Additionally, the table below lists the [supported compilers](#supported-compilers)

## Supported Compilers, Operating Systems, and software stacks

The following combinations are tested on the main branch of self :

Name | Version | Platform | Build System | Stack | Architecture
--- | --- | --- | --- | --- | --- |
GNU Fortran `gfortran` | 13.2.0 | Ubuntu 22.04.2 LTS | `cmake` | openmpi/5.0.1, feq-parse/2.2.2, hdf5/1.14.3 | x86_64 - gfx90a (MI210)
GNU Fortran `gfortran` | 13.2.0 | Ubuntu 22.04.2 LTS | `cmake` | openmpi/5.0.1, feq-parse/2.2.2, hdf5/1.14.3 | x86_64
GNU Fortran `gfortran` | 12.3.0 | Ubuntu 22.04.2 LTS | `cmake` | openmpi/5.0.1, feq-parse/2.2.2, hdf5/1.14.3 | x86_64
GNU Fortran `gfortran` | 12.3.0 | Ubuntu 22.04.2 LTS | `cmake` | openmpi/5.0.3 (ucx+rocm), feq-parse/2.2.2, hdf5/1.14.3 | x86_64 - gfx90a (MI210)
GNU Fortran `gfortran` | 11.4.0 | Ubuntu 22.04.2 LTS | `cmake` | openmpi/5.0.1, feq-parse/2.2.2, hdf5/1.14.3 | x86_64
GNU Fortran `gfortran` | 10.5.0 | Ubuntu 22.04.2 LTS | `cmake` | openmpi/5.0.1, feq-parse/2.2.2, hdf5/1.14.3 | x86_64
GNU Fortran `gfortran` | 9.5.0 | Ubuntu 22.04.2 LTS | `cmake` | openmpi/5.0.1, feq-parse/2.2.2, hdf5/1.14.3 | x86_64
AOMP `amdflang` | 6.1.2 | Ubuntu 22.04.2 LTS | `cmake` | openmpi/5.0.1, feq-parse/2.2.2, hdf5/1.14.3 | x86_64 - gfx90a (MI210)
AOMP `amdflang` | 6.1.2 | Ubuntu 22.04.2 LTS | `cmake` | openmpi/5.0.1, feq-parse/2.2.2, hdf5/1.14.3 | x86_64


"Supported" for us means that we test `self` regularly on the platforms listed. Of course, we want to have `self` working on as many platforms as possible; [open an issue](https://github.com/FluidNumerics/SELF/issues/new/choose) if you encounter any problems installing or running `self` on your own platform.

## Support

### Documentation

* [**User & Developer Documentation**](https://fluidnumerics.github.io/SELF)
* [**API Documentation**](https://fluidnumerics.github.io/SELF/ford/)


If you'd like to contribute, see [CONTRIBUTING.md](./CONTRIBUTING.md) to get started.

If you need help, [open an issue](https://github.com/FluidNumerics/SELF/issues/new)

