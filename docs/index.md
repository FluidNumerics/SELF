# Spectral Element Library in Fortran

## SELF
The Spectral Element Library in Fortran (SELF) is an object oriented suite of Fortran modules that enable researchers to implement conservation law solvers using Spectral Element Methods. One goal of this project is to help dispel the notion that Fortran codes are ancient monoliths designed for a single (potentially unknown) purpose. 

SELF is built with modern compute architecture in mind. Our intention is to enable scientific software developers to build SELF for personal workstations, GPU and multi-GPU accelerated platforms, and HPC clusters. Portable GPU acceleration in SELF is enabled by AMD's HIP and additional coarse-grained parallelism is implemented through domain decomposition with MPI.

 A main program (`self`) is able to digest JSON input files to select and drive ready-to-go conservation law solvers, including 
 
 * Viscous Burger's equation in 1-D
 * Linear Shallow Water Equations in 2-D
 * Non-linear Shallow Water Equations in 2-D
 * Compressible Navier Stokes in 2-D

## PySELF
If you prefer working in Python, our team is actively working on building out a [Python interface for SELF](https://fluidnumerics.github.io/pyself). The goal of [PySELF](https://fluidnumerics.github.io/pyself) is to help you with visualizing and analyzing SELF output. Ultimately, our goal is to provide a python API for you to manage complete simulation workflows from Python. 


## Support

SELF is made available as a free and open-source software to give you the opportunity to work with the code at your own pace and expense. You can report bugs, issues, and feature requests using the repository's issue tracker. Fluid Numerics, the developers of SELF, will prioritize tasks based on paying customer demand (first) and community demand (second).

Fluid Numerics LLC provides a range of professional services for any of the following

* Implementing new features
* Resolving bugs
* Porting SELF to new hardware
* On call support
* Mathematical Modeling and Numerical Analysis
* SELF usage training

You can reach out for professional support at [support@fluidnumerics.com](mailto:support@fluidnumerics.com)
