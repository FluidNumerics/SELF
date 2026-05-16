#include <stdio.h>
#include <mpi.h>
#include <mpi-ext.h> /* Needed for ROCm-aware check */

extern "C"
{
    int check_gpu_aware_support()
    {
        int gpuaware = 0;

#if defined(OMPI_HAVE_MPI_EXT_ROCM) && OMPI_HAVE_MPI_EXT_ROCM
        gpuaware = (int) MPIX_Query_rocm_support();
        printf("Query rocm support");
#endif

#if defined(OMPI_HAVE_MPI_EXT_CUDA) && OMPI_HAVE_MPI_EXT_CUDA
        gpuaware = (int) MPIX_Query_cuda_support();
#endif

        return gpuaware;

    }
}