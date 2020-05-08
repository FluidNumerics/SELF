#include <hip/hip_runtime.h>
#include "SELF_Macros.h"


__global__ void GridInterpolate_1D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;

  real fm = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    fm += f[SP_1D_INDEX(ii,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
  }
  fInterp[SP_1D_INDEX(i,iVar,iEl,M,nVar)] = fm;

}

extern "C"
{
  void GridInterpolate_1D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((GridInterpolate_1D_gpu), dim3(nVar,nEl,1), dim3(M+1,1,1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
  } 
}
