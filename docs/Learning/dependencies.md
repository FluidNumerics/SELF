# Dependencies

This documentation provides more detailed information about SELF's dependencies.


## Message Passing Interface (MPI)

### Background
In order to support scaling of conservation law solvers to large supercomputers, SELF exposes additional parallelism through domain decomposition. This strategy divides a physical domain into multiple subdomains and gives ownership of each subdomain to a single process. These processes are often called "MPI ranks" or "MPI processes". So long as the SELF program you write enables domain decomposition during mesh creation, launching with `mpirun -np 4 /path/to/your/program` will provide you with domain decomposition across four MPI ranks. 

With domain decomposition, each process has its own private memory. However, your conservation law likely demands that neighboring elements in your mesh participate in flux calculations together. When two neighboring elements reside on two separate processes, information must be shared. In SELF, this is achieved as part of the `SideExchange` methods for each data type when domain decomposition is enabled.

On CPU-only platforms, the `boundary` attributes of the `SELF_Mapped*` types are passed to `MPI_ISEND` (asynchronous send) to send data to neighboring processes. On the flip size, the `extboundary` attributes are passed to the `MPI_IRECV` (asynchronous recv) to receive data from the neighboring processes. With each `MPI_ISEND/IRECV` pair, a unique tag for the message is calculated using the global edge/face ID and variable ID. This ensures the messages get stored in the correct address locations when they are received.

On GPU and multi-GPU accelerated platforms, we assume that you are working with a GPU-Aware installation of MPI. In this case, the GPU pointer for the `boundary` attribute (`boundary_gpu`) is passed to `MPI_ISEND`. Similarly, `extboundary_gpu` is passed to `MPI_IRECV`. If your installation of MPI is not GPU aware and you are using a GPU accelerated build of SELF, your program will halt during the initialization of the `DomainDecomposition` class. At this stage, we check the results of `MPIX_Query_rocm_support` and `MPIX_Query_cuda_support` to determine if you have GPU aware MPI enabled at runtime.

### What is supported for MPI+GPU ?
At the moment, we have only tested SELF with OpenMPI with GPU aware support on AMD and Nvidia GPU platforms. As requested by users, we can work to test other MPI flavors ( [Open an issue](https://github.com/FluidNumerics/SELF/issues/new/choose) ). To find if your OpenMPI is built with GPU-Aware Support, you can use `ompi_info` and look for the `MPI extensions` attribute. For AMD platforms, you should see `rocm` , and for Nvidia platforms, you should see `cuda`; if you do not see these, your OpenMPI installation is not GPU aware
```shell
ompi_info | grep "MPI extensions"
       MPI extensions: affinity, cuda, ftmpi, rocm
```

### How to install OpenMPI with ROCm for AMD GPUs
Fluid Numerics recommends building OpenMPI with the Unified Communications Framework (UCX) with ROCm enabled. At the time we are writing this documentation, Spack does not have such a build option available for OpenMPI (though we're hoping it happens soon!). Instead, we can provide this build script below that we have found to work on Debian based linux platforms.

```shell
#!/bin/bash

export ROCM_PATH=/opt/rocm-6.0.2 ## Set this to the path of your ROCm installation
export INSTALL_ROOT=/apps/experimental/rocm-6.0.2 ## Set this to the path you want UCX and OpenMPI installed to
export BUILD_DIR=${INSTALL_ROOT}/tmp/
export UCX_DIR=${INSTALL_ROOT}/ucx
export OMPI_DIR=${INSTALL_ROOT}/ompi_ucx

module load rocm/6.0.2
module load gcc/12.3.0

rm -rf ${BUILD_DIR}
mkdir -p ${BUILD_DIR}


# Install ucx
cd $BUILD_DIR
git clone https://github.com/openucx/ucx.git -b v1.15.x
cd ucx
./autogen.sh
mkdir build
cd build
../configure -prefix=${UCX_DIR} \
    --with-rocm=${ROCM_PATH}
make -j $(nproc)
make -j $(nproc) install

# Install openmpi
cd $BUILD_DIR
git clone --recursive https://github.com/open-mpi/ompi.git \
    -b v5.0.x
cd ompi
./autogen.pl
mkdir build
cd build
../configure --prefix=$OMPI_DIR \
	     --with-ucx=$UCX_DIR \
             --with-rocm=${ROCM_PATH} \
	     --with-hwloc=internal \
             --with-slurm=/usr/local
make -j $(nproc)
make install
```


## HDF5
SELF is set up so that you can focus on building your conservation law solver easily and not worry about other more mundane coding tasks. When you run a simulation, you need to get data out somehow. SELF uses [HDF5](https://www.hdfgroup.org/solutions/hdf5/) to support serial and parallel file IO. All of SELF's parent model classes provide a simple `WriteModel` API call that will write the `solution` data (on the quadrature grid) to file, even if you are using domain decomposition. These files are particularly useful for doing "pickup" runs or for doing analysis offline in another toolkit. 

In order to provide serial and parallel IO functionality, SELF requires HDF5 built with Fortran and MPI bindings. When using spack, this can easily be installed using `spack install hdf5+fortran+mpi`.


## FEQ-Parse
SELF conveniently provides a simple way to set boundary conditions and initial conditions for your model using equations stored in Fortran `CHARACTER` variables. By using the `SetEquation` API calls, you can easily define initial and boundary conditions for your model. Under the hood, SELF uses Fluid Numerics' [`feq-parse`](https://github.com/fluidnumerics/feq-parse) to parse the `CHARACTER` equations and evaluate them on your model grid.

Fluid Numerics maintains the `feq-parse` spack package, which you can easily install using `spack install feq-parse`.

