.. SELF documentation master file, created by
   sphinx-quickstart on Thu Jul 22 15:17:21 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=======
SELF
=======
The Spectral Element Library in Fortran (SELF) is a library of Fortran classes that assist in the construction of PDE solvers using Spectral Element Methods. SELF is developed using a pattern driven design based on interpolation and derivative operations with Lagrange interpolating polynomials. Two and three-dimensional algorithms are provided for working with logically quadrilateral (2-D) and hexahedral (3-D) elements in an unstructured mesh. All elements can be isoparametric so that high order representations of physical boundaries can be preserved. A memory management class is provided to support CPU and GPU memory on both AMD and Nvidia GPUs and all interpolation and differentiation methods have implementations for both CPU and GPU. Additionally, SELF can be built using an MPI+ model so that multi-GPU systems can be targeted.

SELF is meant to be used as a library to develop PDE solvers. To support this, you can obtain Docker containers, Singularity containers, or VM Images for Google Cloud that provide you with a complete development environment to get started.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started:
   QuickStart/shallow_water_equations
   QuickStart/maxwells_equations_2d
   QuickStart/maxwells_equations_3d

.. toctree::
   :maxdepth: 2
   :caption: For Developers:
   Developers/ci_build_system



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
