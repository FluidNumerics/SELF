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

### Building with Google Cloud Build
To get started, you'll need to build a dependency container for SELF-Fluids. For instructions
on building the dependency container with [Google Cloud Build](https://cloud.google.com/cloud-build/),
see the [cloudbuild README](cloudbuild/README.md).

Once you have a dependency container stored in [Google Container Registry](https://cloud.google.com/container-registry/),
you can build a SELF-Fluids singularity containter by executing
```
gcloud builds submit . --substitutions=_BUILD_BASE=<BUILD BASE>
```
where `<BUILD BASE>` is replaced with the prefix of the dependency container.

### Autotools builds
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
export PATH=$PATH}:/path/to/install/directory/bin
```
to your `.bashrc` file.

Additional options for the configure step include
```
  --enable-mpi
  --enable-openmp
  --enable-cuda
  --enable-double-precision
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



## Input
At run-time, the sfluid executable uses the `runtime.params` file for determining the length of 
the run, mesh size (if structured), and other numerical and physical parameters. Additionally,
the initial conditions, drag forces, and boundary conditions and topography (for structured mesh) 
are set in `self.equations`. Together, these two files can completely define a fluid simulation

## Running
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

If you have initial conditions and a mesh file, you can run
```
sfluid integrate
```

Note that if you are running with MPI enabled, you will need to prepend `mpirun -np XX`, replacing
`XX` with the number of MPI processes.

## Output

### Mesh Files
The mesh files are stored in HDF5 format. The mesh file has the following data layout
```
   GROUP "mesh" 
      GROUP "global" 
         GROUP "elements" 
            DATASET "node-ids" 
         GROUP "faces" 
            DATASET "boundary-ids" 
            DATASET "element-ids" 
            DATASET "element-sides" 
            DATASET "i-inc" 
            DATASET "i-start" 
            DATASET "j-inc" 
            DATASET "j-start" 
            DATASET "node-ids" 
            DATASET "swap-dimension" 
         GROUP "geometry" 
            DATASET "boundary-positions" 
            DATASET "positions" 
         GROUP "nodes" 
            DATASET "positions" 

      GROUP "decomp" 
         DATASET "element-to-blockid" 
         GROUP "proc000000" 
            DATASET "element-ids" 
            DATASET "face-ids" 
            DATASET "node-ids" 
         GROUP "proc000001" 
            DATASET "element-ids" 
            DATASET "face-ids" 
         GROUP "proc000002" 
            DATASET "element-ids" 
            DATASET "face-ids" 
            DATASET "node-ids" 
         GROUP "proc000003" 
            DATASET "element-ids" 
            DATASET "face-ids" 
            DATASET "node-ids" 
```
The `global` group contains the global mesh structure. SELF-Fluids works with a hexahedral unstructured nodal spectral element isoparametric mesh.

#### Nodes, Elements, and Faces
The hexahedral unstructured characteristic of the mesh means the we store the mesh as logical cube elements with 8 corner nodes. The nodes in the mesh
are defined by their physical position ( /mesh/global/nodes/positions ). The elements (connectivity) are defined by 8 corner nodes (/mesh/global/elements/node-ids).

Although its not required to completely define the unstructured mesh, we also store faces in the mesh in this file (/mesh/global/faces). Faces are defined by
four corner nodes. This definition leads to other defining features, such as the elements that share the face (/mesh/global/faces/element-ids). The relative orientation
of the two elements sharing the face is important for flux calculations in spectral element methods. This orientiation information is stored in
* /mesh/global/faces/element-sides
* /mesh/global/faces/i-inc
* /mesh/global/faces/i-start
* /mesh/global/faces/j-inc
* /mesh/global/faces/j-start
* /mesh/global/faces/swap-dimensions
For faces that are not shared by two elements, they are marked as boundary faces. Along these faces, boundary conditions are typically applied. SELF-Fluids contains solution
storage data structure attributes specifically for handling external boundary states. To index those attributes, the faces are also aligned with a boundary ID (/mesh/global/faces/boundary-ids).

#### Element Geometry
The isoparametric characteristic of the mesh means that the elements can have curved boundaries. This geometry is defined by the physical positions (/mesh/global/geometry/positions) of the
element at the quadrature points within the element. Further, we opt to store the boundary positions of the element faces (/mesh/global/geometry/boundary-positions). The physical positions 
at quadrature points define a mapping from physical space to a computational space. When the mesh is read in from file, we can calculate the covariant and contravariant basis vectors in 
addition to boundary normals and the Jacobian of the mapping.


#### Domain Decomposition
SELF-Fluids is parallelized using data-parallelism. Since the mesh determines the amount of work and (in part) the communication patterns with Spectral Element PDE solvers, the mesh is decomposed
into blocks and assigned to MPI ranks. The decomposition is defined by a mapping of the elements to blocks ( /mesh/decomp/element-to-blockid ). From this, we can derive a group for each MPI rank that
contains
* global element ID's (/mesh/decomp/procXXXXXX/element-ids)
* face ID's (/mesh/decomp/procXXXXXX/face-ids)
* node ID's (/mesh/decomp/procXXXXXX/node-ids)
owned by that rank. The XXXXXX, in practice, are replaced with a 0 padded integer for each rank.


### Fluid State Files
The fluid state files are HDF5 files containing the fluid state at each time level. The fluid state consists of the momentum vector, density, density-weighted temperature, pressure, and a
passive density-weighted tracer. Each state variable is specified at the all of the quadrature points in each element of the mesh. 

The index between "State" and "h5" is a 0-padded integer corresponding to the time stamp DDDDHHMMSSmmm,
where D = day, H = hour, M = minute, S = second, m = millisecond

These files are created during initial condition generation and forward integration. State files are read in during the start of forward integration



## Future Development Plans

* Variable viscosity (LES)
* Implicit time integration
* Topography input to mesh
