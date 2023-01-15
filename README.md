# Spectral Element Libraries in Fortran (SELF)
Copyright 2020-2022 Fluid Numerics LLC

[![codecov](https://codecov.io/gh/FluidNumerics/SELF/branch/main/graph/badge.svg?token=AKKSL5CWK6)](https://codecov.io/gh/FluidNumerics/SELF)


## Licensing
SELF is licensed for use under a [non-commercial use visible-source license](./LICENSE). Fluid Numerics is a small business (with only two owner-employees!) and wants to make SELF available to researchers for academic use. Under the license, you can use, modify, and redistribute SELF so long as attribution is given to Fluid Numerics. However, since we are interested in protecting our time-and-effort investment in SELF, sale and commercial-use of SELF is prohibited under the license.

If you are interested in commercial licensure and would like support from Fluid Numerics, reach out to support@fluidnumerics.com

## About
SELF is an object-oriented Fortran library that support the implementation of Spectral Element Methods for solving partial differential equations.

The SELF API is designed based on the assumption that SEM developers and researchers need to be able to implement derivatives in 1-D and divergence, gradient, and curl in 2-D and 3-D on scalar, vector, and tensor functions using spectral collocation, continuous galerkin, and discontinuous galerkin spectral element methods. Additionally, as we enter the exascale era, we are currently faced with a zoo of compute hardware that is available. Because of this, SELF routines provide support for GPU acceleration through AMD's HIP and support for multi-core, multi-node, and multi-GPU platforms with MPI.

## Support

### Documentation
[**User & Developer Documentation**](https://fluidnumerics.github.io/SELF)
[**API Documentation**](https://fluidnumerics.github.io/SELF/ford/)

### Community


### Maintainers
* [Joseph Schoonover, Fluid Numerics LLC](https://fluidnumerics.com/people)

If you'd like to contribute, see [CONTRIBUTING.md](./CONTRIBUTING.md) to get started.
If you need help, [open an issue](https://github.com/FluidNumerics/SELF/issues/new)

