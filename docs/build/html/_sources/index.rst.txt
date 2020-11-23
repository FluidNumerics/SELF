Spectral Element Libraries in Fortran (SELF)
===============================================

About
*****
SELF is an object-oriented Fortran library that support the implementation of Spectral Element Methods for solving partial differential equations.

The SELF API is designed based on the assumption that SEM developers and researchers need to be able to implement derivatives in 1-D and divergence, gradient, and curl in 2-D and 3-D on scalar, vector, and tensor functions using spectral collocation, continuous galerkin, and discontinuous galerkin spectral element methods. Additionally, as we enter the exascale era, we are currently faced with a zoo of compute hardware that is available. Because of this, SELF routines provide support for GPU acceleration through AMD's HIP and support for multi-core, multi-node, and multi-GPU platforms with MPI.

Guide
******
.. toctree::
   :maxdepth: 2

   getting-started
   software-design
   license



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
