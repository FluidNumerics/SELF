# Spectral Element Libraries in Fortran (SELF)
Copyright 2020-2022 Fluid Numerics LLC

[![Documentation Status](https://readthedocs.org/projects/self/badge/?version=latest)](https://self.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/FluidNumerics/SELF/branch/master/graph/badge.svg)](https://codecov.io/gh/FluidNumerics/SELF)
[![Discord Chat](https://img.shields.io/discord/308323056592486420.svg)](https://discord.gg/57aNxcpYMW)  
[![Youtube](https://img.shields.io/youtube/channel/views/UCpd92vU2HjjTPup-AIN0pkg?style=social)](https://www.youtube.com/playlist?list=PLRO4xf5MdhAv9CNTETor75rANZtBqPVgQ)


SELF is licensed for use under the [Anti-Corporatist Software License](./LICENSE). For other licensure, reach out to support@fluidnumerics.com.

## About
SELF is an object-oriented Fortran library that support the implementation of Spectral Element Methods for solving partial differential equations.

The SELF API is designed based on the assumption that SEM developers and researchers need to be able to implement derivatives in 1-D and divergence, gradient, and curl in 2-D and 3-D on scalar, vector, and tensor functions using spectral collocation, continuous galerkin, and discontinuous galerkin spectral element methods. Additionally, as we enter the exascale era, we are currently faced with a zoo of compute hardware that is available. Because of this, SELF routines provide support for GPU acceleration through AMD's HIP and support for multi-core, multi-node, and multi-GPU platforms with MPI.

## Support

### Documentation
* [**API Documentation**](https://fluidnumerics.github.io/SELF/ford/)
* [**ReadTheDocs** *(Work in Progress)*](https://self.readthedocs.io/en/latest/)

### Community

#### Open Collective
SELF is part of the Higher Order Methods Collective, which is fiscally hosted by [WATERCHaNGE](https://www.waterchange.org).
You can keep track of updates and announcements for livestreams and training events at the [**Higher Order Methods Open Collective **](https://opencollective.com/higher-order-methods).

You can support SELF and related educational activities focused on numerical analysis and higher order methods for solving conservation laws by contributing to the Open Collective.
[![Open Collective](https://github.com/opencollective/opencollective-images/blob/main/src/static/images/contribute.svg)](https://opencollective.com/higher-order-methods/contribute)


### Maintainers
* [Joseph Schoonover, Fluid Numerics LLC](https://fluidnumerics.com/people/joe-schoonover)
* **You** Want to become a maintainer ? Reach out to support@fluidnumerics.com

If you'd like to contribute, see [CONTRIBUTING.md](./CONTRIBUTING.md) to get started.

If you need help, [open an issue](https://github.com/FluidNumerics/SELF/issues/new)

