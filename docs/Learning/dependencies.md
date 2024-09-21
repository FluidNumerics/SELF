# Dependencies

This documentation provides more detailed information about SELF's dependencies.


## GPU-Aware Message Passing Interface (MPI)
In order to support scaling of conservation law solvers to large supercomputers, SELF exposes additional parallelism through domain decomposition. This strategy divides a physical domain into multiple subdomains and gives ownership of each subdomain to a single process. These processes are often called "MPI ranks" or "MPI processes". So long as the SELF program you write enables domain decomposition during mesh creation, launching with `mpirun -np 4 /path/to/your/program` will provide you with domain decomposition across four MPI ranks. 

With domain decomposition, each process has its own private memory. However, your conservation law likely demands that neighboring elements in your mesh participate in flux calculations together. When two neighboring elements reside on two separate processes, information must be shared. In SELF, this is achieved as part of the `SideExchange` methods for each data type when domain decomposition is enabled.

On CPU-only platforms, the `boundary` attributes of the `SELF_Mapped*` types are passed to `MPI_ISEND` (asynchronous send) to send data to neighboring processes. On the flip size, the `extboundary` attributes are passed to the `MPI_IRECV` (asynchronous recv) to receive data from the neighboring processes. With each `MPI_ISEND/IRECV` pair, a unique tag for the message is calculated using the global edge/face ID and variable ID. This ensures the messages get stored in the correct address locations when they are received.

On GPU and multi-GPU accelerated platforms, we assume that you are working with a GPU-Aware installation of MPI. In this case, the GPU pointer for the `boundary` attribute (`boundary_gpu`) is passed to `MPI_ISEND`. Similarly, `extboundary_gpu` is passed to `MPI_IRECV`. If your installation of MPI is not GPU aware and you are using a GPU accelerated build of SELF, your program will halt during the initialization of the `DomainDecomposition` class. At this stage, we check the results of `MPIX_Query_rocm_support` and `MPIX_Query_cuda_support` to determine if you have GPU aware MPI enabled at runtime.

At the moment, we have only tested SELF with OpenMPI with GPU aware support on AMD and Nvidia GPU platforms. As requested by users, we can work to test other MPI flavors ( [Open an issue](https://github.com/FluidNumerics/SELF/issues/new/choose) ). To find if your OpenMPI is built with GPU-Aware Support
```shell
ompi_info --parsable --all | grep mpi_built_with_cuda_support:value 
mca:mpi:base:param:mpi_built_with_cuda_support:value:true
```