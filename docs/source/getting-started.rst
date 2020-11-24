###############
Getting Started
###############

**************************
Installation with CMake
**************************
This section of the documentation will walk you through installing SELF on your own system


Dependencies
============
SELF requires a 2008 compliant Fortran compiler in addition to the following package dependencies

* HIP
* hipfort
* JSON-fortran
* FEQParse

Optionally, for parallel builds leveraging domain decomposition, an MPI implementation is required


Build
==============
To build and install SELF with GPU support (recommended)::

 FC="/path/to/hipfc" \
 CC="/path/to/hipfc" \
 CXX="/path/to/hipfc" \
 FFLAGS="-v -DGPU -ffree-line-length-none" \
 CXXFLAGS="-v" \
 cmake -DCMAKE_INSTALL_PREFIX="/path/to/install"
 make
 make install


******
Docker
******

***********
Singularity
***********


