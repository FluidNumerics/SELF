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


First, to get started, change directories into the boundary layer example
```
cd SELF-Fluids/examples/boundarylayer
```

There are a variety of flavors of the sfluid executable that can be built dependent on 
environment variables that control things like the fortran compiler used, gpu acceleration,
in-situ timers, integrated diagnostics, passive tracers, etc. Currently yhe environment
variables that control these features are set in a file called `SELF_environment_settings`.
To build with various options, you'll usually, change this file, source it, and recompile.

### Set up for serial CPU
Open `SELF_environment_settings` in a text editor. Make sure that `MPI=no`, `OpenMP=no`, and
`CUDA=no`. Also check to make sure `FC` is set to a fortran compiler in your path. 

### Set up with MPI support
Open `SELF_environment_settings` in a text editor. Make sure that `MPI=yes`, `OpenMP=no`, and
`CUDA=no`. Also check to make sure `MPIFC` and is set to an MPI-Fortran compiler in your path. 

### Set up with CUDA-Fortran support
Open `SELF_environment_settings` in a text editor. Make sure that `MPI=no`, `OpenMP=no`, and
`CUDA=yes`. Also check to make sure `FC=pgfortran`. 

Note that to use CUDA-Fortran support, you must have the PGI Compilers installed.

### Set up with OpenMP support
Open `SELF_environment_settings` in a text editor. Make sure that `MPI=no`, `OpenMP=yes`, and
`CUDA=no`. Also check to make sure `FC` is set to a fortran compiler in your path. 

### Set up with MPI+CUDA support
Open `SELF_environment_settings` in a text editor. Make sure that `MPI=yes`, `OpenMP=no`, and
`CUDA=yes`. Also check to make sure `FC=pgfortran` and `MPIFC=mpif90`. 

### Set up with MPI+OpenMP support
Open `SELF_environment_settings` in a text editor. Make sure that `MPI=yes`, `OpenMP=yes`, and
`CUDA=no`. Also check to make sure `FC` is set to a fortran compiler in your path. 

Note that OpenMP cannot be used when CUDA is enabled.

### Compiling
Once your environment settings file is set up, do
```
source SELF_environment_settings
make sfluid
```

### Running
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


