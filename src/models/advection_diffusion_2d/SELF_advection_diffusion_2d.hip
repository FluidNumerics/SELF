#include <hip/hip_runtime.h>
#include "../../SELF_HIP_Macros.h"

__global__ void SetBoundaryCondition_advection_diffusion_2d_gpu_kernel(int *sideInfo, real *extBoundary, int N, int nEl){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*nEl*4;
  uint32_t i1 = idof % (N+1);
  uint32_t s1 = (idof/(N+1)) % 4;
  uint32_t e1 = idof/(N+1)/4;
  size_t ivar = blockIdx.y;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];

  if(e2 == 0){
    extBoundary[idof + ndof*ivar] = 0.0;
  }
  
}

extern "C"
{
  void SetBoundaryCondition_advection_diffusion_2d_gpu(int *sideInfo, real *extBoundary, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*4*nEl;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nVar,1);
    dim3 nthreads(threads_per_block,1,1);
    SetBoundaryCondition_advection_diffusion_2d_gpu_kernel<<<nblocks,nthreads>>>(sideInfo, extBoundary, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}