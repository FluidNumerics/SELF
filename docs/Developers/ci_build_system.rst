.. SELF documentation master file, created by
   sphinx-quickstart on Thu Jul 22 15:17:21 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. ci_build_system:

====================================
Continuous Integration Build System
====================================
A good continuous integration system for a scientific code answers the following questions

* Does the code build ?
* Does the code run without failure ?
* Does the code get the right answer ?
* Does the code perform well ?

Further, since SELF is meant to be run on a variety of platforms, these questions need to be answered unders a variety of build and runtime configurations to ensure portability.

The SELF repository contains Fortran code and a make system that is used to create the SELF library and a command line interface for exercising each of the classes and their type-bound procedures. The ``ci/`` directory contains infrastructure as code for an image bakery (hosted on Google Cloud) that produces Docker and Singularity container images and GCE VM images. Each of these artifacts are subjected to tests that call the ``self`` CLI to verify that all of the routines run without failure and obtain the correct answer. To understand each subroutine's performance, our CI tests run the CLI commands underneath the hpc-toolkit profiler and process the profiles so that they can be saved in a database (BigQuery) for analysis by the development team. This provides us the ability to ensure performance regressions are not introduced over time and to associate performance changes with specific commits.

To determine correctness, our CI tests run each of the interpolation and differentiation routines with a variety of mesh inputs, over a range of polynomial degrees, and using all of the supported quadrature types. Each CLI call produces HDF5 output that is analyzed to provide a measure of the numerical error. For each CLI call, a measure of the error is stored in our CI database, alongside exit codes and profiles. This provides the development team the ability to ensure that correctness errors are not introduced over time and to associate numerical error cahnges with specific commits.

