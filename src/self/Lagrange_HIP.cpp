#include <hip/hip_runtime.h>
#include "SELF_Macros.h"


// GridInterpolate_1D //
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


// GridDerivative_1D //
__global__ void GridDerivative_1D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;

  real fm = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    fm += f[SP_1D_INDEX(ii,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
  }
  df[SP_1D_INDEX(i,iVar,iEl,N,nVar)] = fm;

}

extern "C"
{
  void GridDerivative_1D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((GridDerivative_1D_gpu), dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0, *dMatrix, *f, *df, N, nVar);
  } 
}