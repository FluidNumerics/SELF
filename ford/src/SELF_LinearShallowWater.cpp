#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"
#include <cstdio>

__global__ void Flux_LinearShallowWater_gpu(real *flux, real *solution, real g, real H, int N, int nVar){

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    if( iVar == 0 ){
      flux[VE_2D_INDEX(1,i,j,iVar,iEl,N,1)] = g*solution[SC_2D_INDEX(i,j,2,iEl,N,nVar)];
      flux[VE_2D_INDEX(2,i,j,iVar,iEl,N,1)] = 0.0;
    } else if ( iVar == 1) {
      flux[VE_2D_INDEX(1,i,j,iVar,iEl,N,1)] = 0.0;
      flux[VE_2D_INDEX(2,i,j,iVar,iEl,N,1)] = g*solution[SC_2D_INDEX(i,j,2,iEl,N,nVar)];
    } else if ( iVar == 2) {
      flux[VE_2D_INDEX(1,i,j,iVar,iEl,N,1)] = H*solution[SC_2D_INDEX(i,j,0,iEl,N,nVar)];
      flux[VE_2D_INDEX(2,i,j,iVar,iEl,N,1)] = H*solution[SC_2D_INDEX(i,j,1,iEl,N,nVar)];
    }
}

extern "C"
{
  void Flux_LinearShallowWater_gpu_wrapper(real **flux, real **solution, real g, real H, int N, int nVar, int nEl)
  {

    // Block size is set to match the size of the element exactly
    // Grid size is set to ( number of tracers X number of elements )
    // DGSEM is beautiful
    Flux_LinearShallowWater_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*flux, *solution, g, H, N, nVar);
  }
}
