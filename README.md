# Spectral Element Libraries in Fortran (SELF)
Copyright 2020-2023 Fluid Numerics LLC

[![codecov](https://codecov.io/gh/FluidNumerics/SELF/branch/main/graph/badge.svg?token=AKKSL5CWK6)](https://codecov.io/gh/FluidNumerics/SELF)


## Licensing
SELF is licensed for use under a [non-commercial use visible-source license](./LICENSE). Fluid Numerics is a small family-owned business and wants to make SELF available to researchers for academic use. Under the license, you can use, modify, and redistribute SELF so long as attribution is given to Fluid Numerics. However, since we are interested in protecting our time-and-effort investment in SELF, sale and commercial-use of SELF is prohibited under the license.

If you are interested in commercial licensure and would like support from Fluid Numerics, reach out to support@fluidnumerics.com

## About
SELF is an object-oriented Fortran library that support the implementation of Spectral Element Methods for solving partial differential equations.

The SELF API is designed based on the assumption that SEM developers and researchers need to be able to implement derivatives in 1-D and divergence, gradient, and curl in 2-D and 3-D on scalar, vector, and tensor functions using spectral collocation, continuous galerkin, and discontinuous galerkin spectral element methods. Additionally, as we enter the exascale era, we are currently faced with a zoo of compute hardware that is available. Because of this, SELF routines provide support for GPU acceleration through AMD's HIP and support for multi-core, multi-node, and multi-GPU platforms with MPI.

## Installation
`self` can be installed using CMake on Linux platforms that have the following packages already installed

* CMake (3.21-3.27)
* 2008 standard compliant Fortran Compiler
* ROCm 5.7.0 or greater
* CUDA 11 or greater (if building for Nvidia GPU)
* HDF5
* [FEQParse](https://github.com/fluidnumerics/feq-parse)

### Prerequisites
All you need is a Fortran compiler that is compliant with the Fortran 2008 standard and supports C interoperability. You can see which compilers are regularly tested on the [Github actions page](https://github.com/FluidNumerics/feq-parse/actions/workflows/ci.yml). Additionally, the table below lists the [supported compilers](#supported-compilers)

### CMake
For a quick installation to `${HOME}/.local/self`,
```
mkdir build/
cd build/
cmake ../ -DCMAKE_INSTALL_PREFIX=${HOME}/.local/self
make
sudo make install
```
If you'd like to run the provided tests to verify your installation, use `ctest` to run the provided tests from within the `build/` directory
```
ctest .
```

The above steps install
```
${HOME}/.local/self/lib/libself-static.a
${HOME}/.local/self/lib/libself.so
${HOME}/.local/self/include/*.mod
${HOME}/.local/self/example/
${HOME}/.local/self/test
```

## Supported Compilers, Operating Systems, and software stacks

The following combinations are tested on the main branch of self :

Name | Version | Platform | Build System | Stack | Architecture
--- | --- | --- | --- | --- | --- |
GNU Fortran | 13.2.0 | Ubuntu 22.04.2 LTS | `cmake` | openmpi/5.0.0, feq-parse/2.0.3, hdf5/1.12.2 | x86_64 - gfx90a
GNU Fortran | 13.2.0 | Ubuntu 22.04.2 LTS | `cmake` | openmpi/5.0.0, feq-parse/2.0.3, hdf5/1.12.2 | x86_64 - gfx906

"Supported" for us means that we test `self` regularly on the platforms listed. Of course, we want to have `self` working on as many platforms as possible; [open an issue](https://github.com/FluidNumerics/SELF/issues/new/choose) if you encounter any problems installing or running `self` on your own platform.

## Support

### Documentation

* [**User & Developer Documentation**](https://fluidnumerics.github.io/SELF)
* [**API Documentation**](https://fluidnumerics.github.io/SELF/ford/)


If you'd like to contribute, see [CONTRIBUTING.md](./CONTRIBUTING.md) to get started.

If you need help, [open an issue](https://github.com/FluidNumerics/SELF/issues/new)

