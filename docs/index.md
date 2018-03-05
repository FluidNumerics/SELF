## About
 The Spectral Element Libraries in Fortran (SELF) are a set of Fortran modules that define data structures and type-bound procedures that can be used to develop Spectral Element discretizations for partial differential equations. SELF-Fluids is a an extended set of modules built on top of SELF to solve Compressible Navier Stokes using the Spectral Element Method on CPUs and GPUs.
 
SELF-Fluids is accelerated on GPUs with CUDA Fortran. CUDA, MPI only, and MPI+CUDA flavors of SELF-Fluids executables are possible.
 
 
 ### Verification
 
 The heart of the SELF-Fluids solvers are the divergence and gradient operators in mapped coordinate systems.
 
