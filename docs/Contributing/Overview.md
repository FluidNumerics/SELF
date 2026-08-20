# Contributing

Fluid Numerics developers are the primary maintainers of the Spectral Element Library in Fortran. We welcome contributions from the broader commmunity. Before getting too deep into developing code, it's best to take a look at the [public issue tracker on Github for SELF](https://github.com/FluidNumerics/SELF/issues). This is the first spot where bugs and feature requests are raised and discussed. Discussion here can help scope software development activities and ensure that any proposed changes will be accepted into SELF down the road.

Although SELF is written to be supported on AMD and Nvidia GPU hardware, you do not need to have this hardware on your systems in order to contribute. Fluid Numerics maintains build actions that work with our [Galapagos Cluster](https://galapagos.fluidnumerics.com/), which enables us to test on a variety of platforms automatically.

## Setting up the developer environment
To help enable developers, we provide a registry of Docker container images that can be deployed on local workstations or HPC environmnents using `docker`, `singularity`/`apptainer`, `podman`, or `enroot`. If you prefer bare-metal configurations, we also provide [`spack`](https://spack.io) environment files so that you can easily set up a [spack environment](https://spack.readthedocs.io/en/latest/environments.html). 

This section of the documentation will help you get starred with setting up your software development environment.

### Bare-metal environments with Spack


