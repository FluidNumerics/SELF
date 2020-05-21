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
  real fi = 0.0;
  for (int jj=0; jj<N+1; jj++) {
    fi = 0.0;
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

  real fij[2] = {0.0};
  real fi[2] = {0.0};
  for (int jj=0; jj<N+1; jj++) {
    fi[1] = 0.0;
    fi[2] = 0.0;
    for (int ii=0; ii<N+1; ii++) {
      fi[1] += f[VE_2D_INDEX(1,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      fi[2] += f[VE_2D_INDEX(2,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
    }
    fij[1] += fi[1]*iMatrix[jj+j*(N+1)];
    fij[2] += fi[2]*iMatrix[jj+j*(N+1)];
  }
  fInterp[VE_2D_INDEX(1,i,j,iVar,iEl,M,nVar)] = fij[1];
  fInterp[VE_2D_INDEX(2,i,j,iVar,iEl,M,nVar)] = fij[2];

}

extern "C"
{
  void VectorGridInterp_2D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((VectorGridInterp_2D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
  } 
}

// TensorGridInterp_2D //
__global__ void TensorGridInterp_2D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;

  real fij[4] = {0.0};
  real fi[4] = {0.0};
  for (int jj=0; jj<N+1; jj++) {
    fi[1] = 0.0;
    fi[2] = 0.0;
    fi[3] = 0.0;
    fi[4] = 0.0;
    for (int ii=0; ii<N+1; ii++) {
      fi[1] += f[TE_2D_INDEX(1,1,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      fi[2] += f[TE_2D_INDEX(2,1,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      fi[3] += f[TE_2D_INDEX(1,2,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      fi[4] += f[TE_2D_INDEX(2,2,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
    }
    fij[1] += fi[1]*iMatrix[jj+j*(N+1)];
    fij[2] += fi[2]*iMatrix[jj+j*(N+1)];
    fij[3] += fi[3]*iMatrix[jj+j*(N+1)];
    fij[4] += fi[4]*iMatrix[jj+j*(N+1)];
  }
  fInterp[TE_2D_INDEX(1,1,i,j,iVar,iEl,M,nVar)] = fij[1];
  fInterp[TE_2D_INDEX(2,1,i,j,iVar,iEl,M,nVar)] = fij[2];
  fInterp[TE_2D_INDEX(1,2,i,j,iVar,iEl,M,nVar)] = fij[3];
  fInterp[TE_2D_INDEX(2,2,i,j,iVar,iEl,M,nVar)] = fij[4];

}

extern "C"
{
  void TensorGridInterp_2D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((TensorGridInterp_2D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
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
  real fij = 0.0;
  real fi = 0.0;
  for (int kk=0; kk<N+1; kk++) {
    fij = 0.0;
    for (int jj=0; jj<N+1; jj++) {
      fi = 0.0;
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

  real fijk[3] = {0.0};
  real fij[3] = {0.0};
  real fi[3] = {0.0};
  for (int kk=0; kk<N+1; kk++) {
    fij[1] = 0.0;
    fij[2] = 0.0;
    fij[3] = 0.0;
    for (int jj=0; jj<N+1; jj++) {
      fi[1] = 0.0;
      fi[2] = 0.0;
      fi[3] = 0.0;
      for (int ii=0; ii<N+1; ii++) {
        fi[1] += f[VE_3D_INDEX(1,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[2] += f[VE_3D_INDEX(2,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[3] += f[VE_3D_INDEX(3,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      }
      fij[1] += fi[1]*iMatrix[jj+j*(N+1)];
      fij[2] += fi[2]*iMatrix[jj+j*(N+1)];
      fij[3] += fi[3]*iMatrix[jj+j*(N+1)];
    }
    fijk[1] += fij[1]*iMatrix[kk+k*(N+1)];
    fijk[2] += fij[2]*iMatrix[kk+k*(N+1)];
    fijk[3] += fij[3]*iMatrix[kk+k*(N+1)];
  }
  fInterp[VE_3D_INDEX(1,i,j,k,iVar,iEl,M,nVar)] = fijk[1];
  fInterp[VE_3D_INDEX(2,i,j,k,iVar,iEl,M,nVar)] = fijk[2];
  fInterp[VE_3D_INDEX(3,i,j,k,iVar,iEl,M,nVar)] = fijk[3];

}

extern "C"
{
  void VectorGridInterp_3D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((VectorGridInterp_3D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,M+1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
  } 
}

// TensorGridInterp_3D //
__global__ void TensorGridInterp_3D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t k = hipThreadIdx_z;

  real fijk[9] = {0.0};
  real fij[9] = {0.0};
  real fi[9] = {0.0};
  for (int kk=0; kk<N+1; kk++) {
    fij[1] = 0.0;
    fij[2] = 0.0;
    fij[3] = 0.0;
    fij[4] = 0.0;
    fij[5] = 0.0;
    fij[6] = 0.0;
    fij[7] = 0.0;
    fij[8] = 0.0;
    fij[9] = 0.0;
    for (int jj=0; jj<N+1; jj++) {
      fi[1] = 0.0;
      fi[2] = 0.0;
      fi[3] = 0.0;
      fi[4] = 0.0;
      fi[5] = 0.0;
      fi[6] = 0.0;
      fi[7] = 0.0;
      fi[8] = 0.0;
      fi[9] = 0.0;
      for (int ii=0; ii<N+1; ii++) {
        fi[1] += f[TE_3D_INDEX(1,1,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[2] += f[TE_3D_INDEX(2,1,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[3] += f[TE_3D_INDEX(3,1,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[4] += f[TE_3D_INDEX(1,2,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[5] += f[TE_3D_INDEX(2,2,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[6] += f[TE_3D_INDEX(3,2,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[7] += f[TE_3D_INDEX(1,3,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[8] += f[TE_3D_INDEX(2,3,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[9] += f[TE_3D_INDEX(3,3,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      }
      fij[1] += fi[1]*iMatrix[jj+j*(N+1)];
      fij[2] += fi[2]*iMatrix[jj+j*(N+1)];
      fij[3] += fi[3]*iMatrix[jj+j*(N+1)];
      fij[4] += fi[4]*iMatrix[jj+j*(N+1)];
      fij[5] += fi[5]*iMatrix[jj+j*(N+1)];
      fij[6] += fi[6]*iMatrix[jj+j*(N+1)];
      fij[7] += fi[7]*iMatrix[jj+j*(N+1)];
      fij[8] += fi[8]*iMatrix[jj+j*(N+1)];
      fij[9] += fi[9]*iMatrix[jj+j*(N+1)];
    }
    fijk[1] += fij[1]*iMatrix[kk+k*(N+1)];
    fijk[2] += fij[2]*iMatrix[kk+k*(N+1)];
    fijk[3] += fij[3]*iMatrix[kk+k*(N+1)];
    fijk[4] += fij[4]*iMatrix[kk+k*(N+1)];
    fijk[5] += fij[5]*iMatrix[kk+k*(N+1)];
    fijk[6] += fij[6]*iMatrix[kk+k*(N+1)];
    fijk[7] += fij[7]*iMatrix[kk+k*(N+1)];
    fijk[8] += fij[8]*iMatrix[kk+k*(N+1)];
    fijk[9] += fij[9]*iMatrix[kk+k*(N+1)];
  }
  fInterp[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,M,nVar)] = fijk[1];
  fInterp[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,M,nVar)] = fijk[2];
  fInterp[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,M,nVar)] = fijk[3];
  fInterp[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,M,nVar)] = fijk[4];
  fInterp[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,M,nVar)] = fijk[5];
  fInterp[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,M,nVar)] = fijk[6];
  fInterp[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,M,nVar)] = fijk[7];
  fInterp[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,M,nVar)] = fijk[8];
  fInterp[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,M,nVar)] = fijk[9];

}

extern "C"
{
  void TensorGridInterp_3D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((TensorGridInterp_3D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,M+1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
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
