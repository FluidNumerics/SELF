## About
 The Spectral Element Libraries in Fortran (SELF) are a set of Fortran modules that define data structures and type-bound procedures that can be used to develop Spectral Element discretizations for partial differential equations. SELF-Fluids is a an extended set of modules built on top of SELF to solve Compressible Navier Stokes using the Spectral Element Method on CPUs and GPUs.
 
SELF-Fluids is accelerated on GPUs with CUDA Fortran. CUDA, MPI only, and MPI+CUDA flavors of SELF-Fluids executables are possible.
 
 
 
 ## Thermal Bubble Demonstration
 This example can be found in the `examples/thermalbubble/` directory of the SELF-Fluids repository. The initial conditions consist of a warm motionless ball of fluid in a neutrally stable background environment. As warm fluid begins to rise, a ring vortex (like a smoke ring) forms and accelerates the fluid upwards.
 
 [![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/NVToKGeOy94/0.jpg)](https://www.youtube.com/watch?v=NVToKGeOy94)
 
## Strong Scaling with and without GPUs

Table 1: The wall times, speedup, and scaling efficiency are shown for the single GPU and select MPI configurations. The system used for this study has sixteen cores per node and two Tesla K40 GPU's per node. For the cases with GPU's, two MPI ranks per node were used, with one rank assigned to each GPU. Identical affinity is used in the MPI-only configurations for 2, 4, and 8 rank configurations. These results are based on the ten averages of instrumented wall times for computing 1,000 simulation time steps.

| No. Ranks	| GPU	| Wall Time |	Speedup	| Scaling Efficiency |
| --- | --- | --- | --- | --- |
| 1 |	yes	| 554.419 |	34.350	| – |
| 2	| yes	| 301.585	| 63.148	| 91.92 % |
| 4	| yes	| 152.702	| 124.717	| 90.77 % |
| 8	| yes	| 77.8016	| 244.785	| 89.08 % |
| 2	| no	| 10127.85	| 1.880	| 94.02 % |
| 4	| no	| 5201.868	| 3.661	| 91.53 % |
| 8	| no	| 2502.976	| 7.609	| 95.11 % |
| 64	| no |	514.882	| 36.988 |	30.73 % |

 ## Boundary Layer Turbulence
 This example can be found in the `examples/boundarylayer/` directory of the SELF-Fluids repository. The initial conditions consist of a neutrally stable fluid in a doubly periodic domain moving uniformly at 10 m/s. At the bottom of the domain, a drag force slows down the fluid resulting in an unstable shear. This evenutally becomes turbulent as depicted in the video below, showing the spatially and temporally varying vertical velocity component.
 
 [![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/74vc87gVgWk/0.jpg)](https://www.youtube.com/watch?v=74vc87gVgWk)
