# SELF-Fluids


## Spectral Element Libraries in Fortran for Modeling Fluids
SELF-Fluids is a modular framework for working with spectral element methods that has 
been organized into a master class ( Fluid ) and a program ( `sfluid` ) that can be 
used to solve the Compressible Navier-Stokes equations.

The Fluid class is employs the Nodal Discontinuous Galerkin Spectral Element method
for the spatial discretization. The inter-element fluxes are numerically estimated
using the upwind Lax-Friedrich's approximate Riemman solver. Currently, the model
is forward-stepped using Williamson's 3rd Order Low-Storage Runge-Kutta time
integrator. There are plans to develop and support additional explicit and implicit
time-stepping schemes.

## Getting Started
SELF-Fluids ships with a number of verified examples to help you get started quickly.

In this documentation, we'll walk through the planetary boundary layer example, which
can be found under `SELF-Fluids/examples/boundarylayer/`. This demo solves the fluid
equations in a doubly periodic domain with no-normal-flow boundary conditions at the
bottom of the domain, and prescribed boundary conditions at the top of the domain. 
The flow is set to 10 m/s in the x-direction initially and at the prescribed boundary.
At the bottom boundary, a drag layer is used to slow the fluid and generate a vertical
shear. Eventually, this shear layer becomes unstable ( Kelvin-Helmholtz ) and turbulence
ensues.


## Compiling and Installing
SELF-Fluids is built and installed using autotools. 
For a serial build, with single precision, no diagnostics, no timers, and without GPU 
acceleration, installation is as simple as
```
./configure --prefix=/path/to/install/directory
make
make install
```

### MPI Support
To enable building of the sfluid executable with MPI support, you can add the flag
`--enable-mpi` at the configure stage.
```
./configure --enable-mpi --prefix=/path/to/install/directory
make
make install
```
This causes the configure script to seek out the mpif90 compiler. Because of this,
it is necessary to have an MPI library, like OpenMPI or MPICH, installed.


### OpenMP Multi-Threaded Acceleration
To enable building of the sfluid executable with multi-threading support via OpenMP,
 you can add the flag `--enable-openmp` at the configure stage.
```
./configure --enable-openmp --prefix=/path/to/install/directory
make
make install
```
This adds the appropriate flag to `FCFLAGS` so that OpenMP directives within the
SELF-Fluids source are interpreted. Note that OpenMP support can be used with MPI, 
but cannot be enabled when GPU accelerations are enabled.

At runtime, you will need to set the environment variable `OMP_NUM_THREADS` to
the desired number of threads.


### GPU Acceleration
To enable building of the sfluid executable with GPU acceleration via CUDA-Fortran,
 you can add the flag `--enable-cuda` at the configure stage. You must also specify
a GPU architecture using the `GPU_ARCH` flag ( e.g `GPU_ARCH=cc60` specifies
compilation for devices with compute capability 6.0 ).
```
./configure --enable-cuda --prefix=/path/to/install/directory GPU_ARCH=<compute-capability>
make
make install
```
Note that CUDA support can be used with MPI. To use the GPU acceleration, you must use
the PGI compilers ( https://www.pgroup.com/products/community.htm )



## Running
At run-time, the sfluid executable uses the `runtime.params` file for determining the length of 
the run, mesh size (if structured), and other numerical and physical parameters. Additionally,
the initial conditions and drag forces are set in `self.equations`. 

The first time the `sfluid`is executed for an example, it will generate a mesh consistent with 
the settings in `runtime.params` and initial conditions consistent with those specified in 
`self.equations`. To get started quickly, run
```
./sfluid
```

If you would just like to generate the mesh, or generate a new mesh from scratch, run
```
./sfluid meshgen
```

If you would like to proceed up to initial condition generation, and not forward step the model,
you can run
```
./sfluid init
```

Note that if you are running with MPI enabled, you will need to prepend `mpirun -np XX`, replacing
`XX` with the number of MPI processes.


### Output
The sfluid executable will generate some files specific to SELF-Fluids ( .mesh, .bcm, and .pickup )
and tecplot output ( .tec ) that can be visualized in VisIt, Paraview, or TecPlot. There are currently
plans to switch to VTK output to replace the tecplot output in later versions.


