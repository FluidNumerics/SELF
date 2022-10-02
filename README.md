# Spectral Element Libraries in Fortran (SELF)
Copyright 2020-2022 Fluid Numerics LLC

[![Documentation Status](https://readthedocs.org/projects/self/badge/?version=latest)](https://self.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/FluidNumerics/SELF/branch/main/graph/badge.svg?token=AKKSL5CWK6)](https://codecov.io/gh/FluidNumerics/SELF)
[![Youtube](https://img.shields.io/youtube/channel/subscribers/UCW5e-TavOnw1AABGH-VMbRg?style=social)](https://www.youtube.com/channel/UCW5e-TavOnw1AABGH-VMbRg?sub_confirmation=1)
[![Reddit](https://img.shields.io/reddit/subreddit-subscribers/fluidnumerics?style=social)](https://www.reddit.com/r/FluidNumerics/)

Join the [Higher Order Methods Slack group](https://join.slack.com/t/higherordermethods/shared_invite/zt-1da6fpyjo-c4yNNXD_o0F3Yrxe8isgJg)

SELF is licensed for use under the [Researcher Software License](./LICENSE). For other licensure, reach out to support@fluidnumerics.com.

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

