## About
 The Spectral Element Libraries in Fortran (SELF) are a set of Fortran modules that define data structures and type-bound procedures that can be used to develop Spectral Element discretizations for partial differential equations. SELF-Fluids is a an extended set of modules built on top of SELF to solve Compressible Navier Stokes using the Spectral Element Method on CPUs and GPUs.
 
SELF-Fluids is accelerated on GPUs with CUDA Fortran. CUDA, MPI only, and MPI+CUDA flavors of SELF-Fluids executables are possible.
 
 
 
 ## Thermal Bubble Demonstration
 This example can be found in the `examples/thermalbubble/` directory of the SELF-Fluids repository. The initial conditions consist of a warm motionless ball of fluid in a neutrally stable background environment. As warm fluid begins to rise, a ring vortex (like a smoke ring) forms and accelerates the fluid upwards.
 
 [![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/NVToKGeOy94/0.jpg)](https://www.youtube.com/watch?v=NVToKGeOy94)
 
