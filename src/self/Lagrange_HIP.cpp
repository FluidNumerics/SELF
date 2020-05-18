#include <hip/hip_runtime.h>
#include "SELF_Macros.h"


// ScalarGridInterp_1D //
__global__ void ScalarGridInterp_1D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;

  real fm = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    fm += f[SC_1D_INDEX(ii,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
  }
  fInterp[SC_1D_INDEX(i,iVar,iEl,M,nVar)] = fm;

}

extern "C"
{
  void ScalarGridInterp_1D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((ScalarGridInterp_1D_gpu), dim3(nVar,nEl,1), dim3(M+1,1,1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
  } 
}

// ScalarGridInterp_2D //
__global__ void ScalarGridInterp_2D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;

  real fij = 0.0;
  for (int jj=0; jj<N+1; jj++) {
    real fi = 0.0;
    for (int ii=0; ii<N+1; ii++) {
      fi += f[SC_2D_INDEX(ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
    }
    fij += fi*iMatrix[jj+j*(N+1)];
  }
  fInterp[SC_2D_INDEX(i,j,iVar,iEl,M,nVar)] = fij;

}

extern "C"
{
  void ScalarGridInterp_2D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((ScalarGridInterp_2D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
  } 
}

// ScalarGridInterp_3D //
__global__ void ScalarGridInterp_3D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t k = hipThreadIdx_z;

  real fijk = 0.0;
  for (int kk=0; kk<N+1; kk++) {
    real fij = 0.0;
    for (int jj=0; jj<N+1; jj++) {
      real fi = 0.0;
      for (int ii=0; ii<N+1; ii++) {
        fi += f[SC_3D_INDEX(ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      }
      fij += fi*iMatrix[jj+j*(N+1)];
    }
    fijk += fij*iMatrix[kk+k*(N+1)];
  }
  fInterp[SC_3D_INDEX(i,j,k,iVar,iEl,M,nVar)] = fijk;

}

extern "C"
{
  void ScalarGridInterp_3D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((ScalarGridInterp_3D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,M+1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
  } 
}


// GridDerivative_1D //
__global__ void GridDerivative_1D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;

  real fm = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    fm += f[SC_1D_INDEX(ii,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
  }
  df[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = fm;

}

extern "C"
{
  void GridDerivative_1D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((GridDerivative_1D_gpu), dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0, *dMatrix, *f, *df, N, nVar);
  } 
}
