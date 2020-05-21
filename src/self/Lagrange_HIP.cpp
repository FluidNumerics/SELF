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

// VectorGridInterp_2D //
__global__ void VectorGridInterp_2D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;

  real fij1 = 0.0;
  real fij2 = 0.0;
  for (int jj=0; jj<N+1; jj++) {
    real fi1 = 0.0;
    real fi2 = 0.0;
    for (int ii=0; ii<N+1; ii++) {
      fi1 += f[VE_2D_INDEX(1,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      fi2 += f[VE_2D_INDEX(2,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
    }
    fij1 += fi1*iMatrix[jj+j*(N+1)];
    fij2 += fi2*iMatrix[jj+j*(N+1)];
  }
  fInterp[VE_2D_INDEX(1,i,j,iVar,iEl,M,nVar)] = fij1;
  fInterp[VE_2D_INDEX(2,i,j,iVar,iEl,M,nVar)] = fij2;

}

extern "C"
{
  void VectorGridInterp_2D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((VectorGridInterp_2D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
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

// VectorGridInterp_3D //
__global__ void VectorGridInterp_3D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t k = hipThreadIdx_z;

  real fijk1 = 0.0;
  real fijk2 = 0.0;
  real fijk3 = 0.0;
  for (int kk=0; kk<N+1; kk++) {
    real fij1 = 0.0;
    real fij2 = 0.0;
    real fij3 = 0.0;
    for (int jj=0; jj<N+1; jj++) {
      real fi1 = 0.0;
      real fi2 = 0.0;
      real fi3 = 0.0;
      for (int ii=0; ii<N+1; ii++) {
        fi1 += f[VE_3D_INDEX(1,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi2 += f[VE_3D_INDEX(2,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi3 += f[VE_3D_INDEX(3,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      }
      fij1 += fi1*iMatrix[jj+j*(N+1)];
      fij2 += fi2*iMatrix[jj+j*(N+1)];
      fij3 += fi3*iMatrix[jj+j*(N+1)];
    }
    fijk1 += fij1*iMatrix[kk+k*(N+1)];
    fijk2 += fij2*iMatrix[kk+k*(N+1)];
    fijk3 += fij3*iMatrix[kk+k*(N+1)];
  }
  fInterp[VE_3D_INDEX(1,i,j,k,iVar,iEl,M,nVar)] = fijk1;
  fInterp[VE_3D_INDEX(2,i,j,k,iVar,iEl,M,nVar)] = fijk2;
  fInterp[VE_3D_INDEX(3,i,j,k,iVar,iEl,M,nVar)] = fijk3;

}

extern "C"
{
  void VectorGridInterp_3D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((VectorGridInterp_3D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,M+1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
  } 
}


// Derivative_1D //
__global__ void Derivative_1D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

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
  void Derivative_1D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((Derivative_1D_gpu), dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0, *dMatrix, *f, *df, N, nVar);
  } 
}
