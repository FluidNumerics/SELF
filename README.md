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
./configure [OPTIONS] --prefix=/path/to/install/directory
make
make install
```
This installs the `sfluid` binary in the directory specified by the `prefix` argument
at the configure stage. To make this binary visible, add
```
export PATH=${PATH}:/path/to/install/directory/bin
```
to your `.bashrc` file.

Additional options for the configure step include
```
  --enable-mpi
  --enable-openmp
  --enable-cuda
  --enable-double-precision
  --enable-timing
```


### MPI Support
To enable building of the sfluid executable with MPI support, you can add the flag
`--enable-mpi` at the configure stage.
```
./configure --enable-mpi --prefix=/path/to/install/directory
make
make install
```
To use MPI, it is necessary to have an MPI library, like OpenMPI or MPICH, installed
and have binaries in your path.


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

The first time the `sfluid` is executed for an example, it will generate a mesh consistent with 
the settings in `runtime.params` and initial conditions consistent with those specified in 
`self.equations`. If `runtime.params` is not present, a sample will be generated for you.
If `self.equations` is not present, a sample equations file will be generated
for you as well.

To get started quickly, run
```
sfluid
```

If you would just like to generate the mesh, or generate a new mesh from scratch, run
```
sfluid meshgen
```

If you would like to proceed up to initial condition generation, and not forward step the model,
you can run
```
sfluid init
```

Note that if you are running with MPI enabled, you will need to prepend `mpirun -np XX`, replacing
`XX` with the number of MPI processes.


### Output
The sfluid executable will generate some files specific to SELF-Fluids ( .mesh, .bcm, and .pickup )
and tecplot output ( .tec ) that can be visualized in VisIt, Paraview, or TecPlot. There are currently
plans to switch to VTK output to replace the tecplot output in later versions.


## Input
The sfluid binary manages mesh generation, initial condition generation,
and forward execution of the compressible fluids model. Command line options
are given as follows : 

```
sfluid [tool] [options]**
```
**[tool]** can be :


 **help**

   Display a help message



 **meshgen**

   Run only the mesh generator to generate a structured mesh.
   The structured mesh is built using nXElem, nYElem, and nZElem
   specified in runtime.params.

   Topography can be set using an equation like


         h = exp( -(x-500.0)^2/200.0 )

   in the self.equations file. This will result in a terrain-following
   structured mesh.


 **init**

   Run up to the initial condition generation and do not forward
   step the model. The initial conditions are read in from the 
   self.equations file. 


 **[options]** can be :
  **--param-file /path/to/param/file**
     Specifies the full path to a file with namelist settings for
     the sfluid application. If not provided, runtime.params in  
     your current directory is assumed.                          

 **--equation-file /path/to/equation/file**
     Specifies the full path to an equation file for setting the 
     initial conditions, topography shape (for structured mesh), 
     and the drag field. If not provided, self.equations in your 
     current directory is assumed.

## Output

### ExtComm.*.bcm
"External Communications" file. This file is an ASCII file generated in
`src/self/BoundaryCommunicator_Class.F90`. Each line of the file corresponds
to a boundary face for a process's domain. A boundary face is a face that
is associated with only one element in the process's mesh. The boundary face
can be associated with a physical boundary condition or an MPI communication
boundary. 

Each line has the unique boundary ID, the external process ID, and 
a mapping back to the local face ID in the owning process's mesh. 
One file is generated for each MPI rank. Between "box" and "geom" is a 0-padded integer
corresponding to the MPI rank ID ( if MPI is not used the rank ID is assumed to be 0 ).

This file is created during mesh generation and is read-only for initial
condition generation and model forward integration.
	
### box.*.mesh
Mesh file. This file is a binary file generated in `src/self/HexMesh_Class.F90`.
This file lists the number of nodes, number of elements, number of faces, polynomial
degree, the element-node connectivity, element-node-face connectivity, face-boundary mapping,
and element neighbor orientation.

One file is generated for each MPI rank. Between "box" and "geom" is a 0-padded integer
corresponding to the MPI rank ID ( if MPI is not used the rank ID is assumed to be 0 ).

This file is created during mesh generation and is read-only for initial
condition generation and model forward integration.

### box.*.geom
Geometry file. This file is a binary file in `src/self/HexMesh_Class.F90`.
This file stores the mesh node positions in addition to the nodal locations
of quadrature points and metric terms at all quadrature points for each element.
One file is generated for each MPI rank. Between "box" and "geom" is a 0-padded integer
corresponding to the MPI rank ID ( if MPI is not used the rank ID is assumed to be 0 ).

This file is created during mesh generation and is read-only for initial
condition generation and model forward integration.


## State.*.h5
State file. This file is an HDF5 file containing the fluid state at each time level. It is
generated in `src/fluid/Fluid_Class.F90`. 
The fluid state consists of the momentum vector, density, density-weighted temperature, pressure, and a
passive density-weighted tracer. Each state variable is specified at the all of the quadrature points
in each element of the mesh. 

The index between "State" and "h5" is a 0-padded integer corresponding to the time stamp DDDDHHMMSSmmm,
where D = day, H = hour, M = minute, S = second, m = millisecond

These files are created during initial condition generation and forward integration. State files
are read in during the start of forward integration
