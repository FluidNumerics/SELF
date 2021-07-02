#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

/*
// Template
__global__ void Template_{1D|2D|3D}_gpu( , int N, int nVar, int nEl){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;


  // How to access scalars 
  f[SC_1D_INDEX(i,iVar,iEl,N,nVar)];
  f[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];
  f[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

  // How to access vectors (2d)
  f[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)];
  f[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)];


  // How to access vectors (3d)
  f[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)];
  f[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)];
  f[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)];

  // How to access tensors (2d)
  f[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)];
  f[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)];
  f[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)];
  f[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)];


  // How to access tensors (3d)
  f[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,N,nVar)];

}

extern "C"
{
  void Template_{1D|2D|3D}_gpu_wrapper(int N, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((Template_{1D|2D|3D}_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0,  N, M, nVar);
  } 
}


*/

// ScalarGridInterp_1D //
__global__ void ScalarGridInterp_1D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

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
	  //hipLaunchKernelGGL((ScalarGridInterp_1D_gpu), dim3(nVar,nEl,1), dim3(M+1,1,1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
	  ScalarGridInterp_1D_gpu<<<dim3(nVar,nEl,1), dim3(M+1,1,1), 0, 0>>>(*iMatrix, *f, *fInterp, N, M, nVar);
  } 
}

// ScalarGridInterp_2D //
__global__ void ScalarGridInterp_2D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

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
	  //hipLaunchKernelGGL((ScalarGridInterp_2D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
	  ScalarGridInterp_2D_gpu<<<dim3(nVar,nEl,1), dim3(M+1,M+1,1), 0, 0>>>(*iMatrix, *f, *fInterp, N, M, nVar);
  } 
}

// VectorGridInterp_2D //
__global__ void VectorGridInterp_2D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real fij[2] = {0.0};
  real fi[2] = {0.0};
  for (int jj=0; jj<N+1; jj++) {
    fi[0] = 0.0;
    fi[1] = 0.0;
    for (int ii=0; ii<N+1; ii++) {
      fi[0] += f[VE_2D_INDEX(1,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      fi[1] += f[VE_2D_INDEX(2,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
    }
    fij[0] += fi[0]*iMatrix[jj+j*(N+1)];
    fij[1] += fi[1]*iMatrix[jj+j*(N+1)];
  }
  fInterp[VE_2D_INDEX(1,i,j,iVar,iEl,M,nVar)] = fij[0];
  fInterp[VE_2D_INDEX(2,i,j,iVar,iEl,M,nVar)] = fij[1];

}

extern "C"
{
  void VectorGridInterp_2D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorGridInterp_2D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
	  VectorGridInterp_2D_gpu<<<dim3(nVar,nEl,1), dim3(M+1,M+1,1), 0, 0>>>(*iMatrix, *f, *fInterp, N, M, nVar);
  } 
}

// TensorGridInterp_2D //
__global__ void TensorGridInterp_2D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real fij[4] = {0.0};
  real fi[4] = {0.0};
  for (int jj=0; jj<N+1; jj++) {
    fi[0] = 0.0;
    fi[1] = 0.0;
    fi[2] = 0.0;
    fi[3] = 0.0;
    for (int ii=0; ii<N+1; ii++) {
      fi[0] += f[TE_2D_INDEX(1,1,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      fi[1] += f[TE_2D_INDEX(2,1,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      fi[2] += f[TE_2D_INDEX(1,2,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      fi[3] += f[TE_2D_INDEX(2,2,ii,jj,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
    }
    fij[0] += fi[0]*iMatrix[jj+j*(N+1)];
    fij[1] += fi[1]*iMatrix[jj+j*(N+1)];
    fij[2] += fi[2]*iMatrix[jj+j*(N+1)];
    fij[3] += fi[3]*iMatrix[jj+j*(N+1)];
  }
  fInterp[TE_2D_INDEX(1,1,i,j,iVar,iEl,M,nVar)] = fij[0];
  fInterp[TE_2D_INDEX(2,1,i,j,iVar,iEl,M,nVar)] = fij[1];
  fInterp[TE_2D_INDEX(1,2,i,j,iVar,iEl,M,nVar)] = fij[2];
  fInterp[TE_2D_INDEX(2,2,i,j,iVar,iEl,M,nVar)] = fij[3];

}

extern "C"
{
  void TensorGridInterp_2D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((TensorGridInterp_2D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
	  TensorGridInterp_2D_gpu<<<dim3(nVar,nEl,1), dim3(M+1,M+1,1), 0, 0>>>(*iMatrix, *f, *fInterp, N, M, nVar);
  } 
}

// ScalarGridInterp_3D //
__global__ void ScalarGridInterp_3D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

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
	  //hipLaunchKernelGGL((ScalarGridInterp_3D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,M+1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
	  ScalarGridInterp_3D_gpu<<<dim3(nVar,nEl,1), dim3(M+1,M+1,M+1), 0, 0>>>(*iMatrix, *f, *fInterp, N, M, nVar);
  } 
}

// VectorGridInterp_3D //
__global__ void VectorGridInterp_3D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  real fijk[3] = {0.0};
  real fij[3] = {0.0};
  real fi[3] = {0.0};
  for (int kk=0; kk<N+1; kk++) {
    fij[0] = 0.0;
    fij[1] = 0.0;
    fij[2] = 0.0;
    for (int jj=0; jj<N+1; jj++) {
      fi[0] = 0.0;
      fi[1] = 0.0;
      fi[2] = 0.0;
      for (int ii=0; ii<N+1; ii++) {
        fi[0] += f[VE_3D_INDEX(1,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[1] += f[VE_3D_INDEX(2,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[2] += f[VE_3D_INDEX(3,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      }
      fij[0] += fi[0]*iMatrix[jj+j*(N+1)];
      fij[1] += fi[1]*iMatrix[jj+j*(N+1)];
      fij[2] += fi[2]*iMatrix[jj+j*(N+1)];
    }
    fijk[0] += fij[0]*iMatrix[kk+k*(N+1)];
    fijk[1] += fij[1]*iMatrix[kk+k*(N+1)];
    fijk[2] += fij[2]*iMatrix[kk+k*(N+1)];
  }
  fInterp[VE_3D_INDEX(1,i,j,k,iVar,iEl,M,nVar)] = fijk[0];
  fInterp[VE_3D_INDEX(2,i,j,k,iVar,iEl,M,nVar)] = fijk[1];
  fInterp[VE_3D_INDEX(3,i,j,k,iVar,iEl,M,nVar)] = fijk[2];

}

extern "C"
{
  void VectorGridInterp_3D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorGridInterp_3D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,M+1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
	  VectorGridInterp_3D_gpu<<<dim3(nVar,nEl,1), dim3(M+1,M+1,M+1), 0, 0>>>(*iMatrix, *f, *fInterp, N, M, nVar);
  } 
}

// TensorGridInterp_3D //
__global__ void TensorGridInterp_3D_gpu(real *iMatrix, real *f, real *fInterp, int N, int M, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  real fijk[9] = {0.0};
  real fij[9] = {0.0};
  real fi[9] = {0.0};
  for (int kk=0; kk<N+1; kk++) {
    fij[0] = 0.0;
    fij[1] = 0.0;
    fij[2] = 0.0;
    fij[3] = 0.0;
    fij[4] = 0.0;
    fij[5] = 0.0;
    fij[6] = 0.0;
    fij[7] = 0.0;
    fij[8] = 0.0;
    for (int jj=0; jj<N+1; jj++) {
      fi[0] = 0.0;
      fi[1] = 0.0;
      fi[2] = 0.0;
      fi[3] = 0.0;
      fi[4] = 0.0;
      fi[5] = 0.0;
      fi[6] = 0.0;
      fi[7] = 0.0;
      fi[8] = 0.0;
      for (int ii=0; ii<N+1; ii++) {
        fi[0] += f[TE_3D_INDEX(1,1,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[1] += f[TE_3D_INDEX(2,1,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[2] += f[TE_3D_INDEX(3,1,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[3] += f[TE_3D_INDEX(1,2,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[4] += f[TE_3D_INDEX(2,2,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[5] += f[TE_3D_INDEX(3,2,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[6] += f[TE_3D_INDEX(1,3,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[7] += f[TE_3D_INDEX(2,3,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
        fi[8] += f[TE_3D_INDEX(3,3,ii,jj,kk,iVar,iEl,N,nVar)]*iMatrix[ii+i*(N+1)];
      }
      fij[0] += fi[0]*iMatrix[jj+j*(N+1)];
      fij[1] += fi[1]*iMatrix[jj+j*(N+1)];
      fij[2] += fi[2]*iMatrix[jj+j*(N+1)];
      fij[3] += fi[3]*iMatrix[jj+j*(N+1)];
      fij[4] += fi[4]*iMatrix[jj+j*(N+1)];
      fij[5] += fi[5]*iMatrix[jj+j*(N+1)];
      fij[6] += fi[6]*iMatrix[jj+j*(N+1)];
      fij[7] += fi[7]*iMatrix[jj+j*(N+1)];
      fij[8] += fi[8]*iMatrix[jj+j*(N+1)];
    }
    fijk[0] += fij[0]*iMatrix[kk+k*(N+1)];
    fijk[1] += fij[1]*iMatrix[kk+k*(N+1)];
    fijk[2] += fij[2]*iMatrix[kk+k*(N+1)];
    fijk[3] += fij[3]*iMatrix[kk+k*(N+1)];
    fijk[4] += fij[4]*iMatrix[kk+k*(N+1)];
    fijk[5] += fij[5]*iMatrix[kk+k*(N+1)];
    fijk[6] += fij[6]*iMatrix[kk+k*(N+1)];
    fijk[7] += fij[7]*iMatrix[kk+k*(N+1)];
    fijk[8] += fij[8]*iMatrix[kk+k*(N+1)];
  }
  fInterp[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,M,nVar)] = fijk[0];
  fInterp[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,M,nVar)] = fijk[1];
  fInterp[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,M,nVar)] = fijk[2];
  fInterp[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,M,nVar)] = fijk[3];
  fInterp[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,M,nVar)] = fijk[4];
  fInterp[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,M,nVar)] = fijk[5];
  fInterp[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,M,nVar)] = fijk[6];
  fInterp[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,M,nVar)] = fijk[7];
  fInterp[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,M,nVar)] = fijk[8];

}

extern "C"
{
  void TensorGridInterp_3D_gpu_wrapper(real **iMatrix, real **f, real **fInterp, int N, int M, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((TensorGridInterp_3D_gpu), dim3(nVar,nEl,1), dim3(M+1,M+1,M+1), 0, 0, *iMatrix, *f, *fInterp, N, M, nVar);
	  TensorGridInterp_3D_gpu<<<dim3(nVar,nEl,1), dim3(M+1,M+1,M+1), 0, 0>>>(*iMatrix, *f, *fInterp, N, M, nVar);
  } 
}

// ScalarBoundaryInterp_1D //
__global__ void ScalarBoundaryInterp_1D_gpu(real *bMatrix, real *f, real *fBound, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t bid = threadIdx.x;

  real fb = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    fb += f[SC_1D_INDEX(ii,iVar,iEl,N,nVar)]*bMatrix[ii+bid*(N+1)];
  }
  fBound[SCB_1D_INDEX(iVar,bid,iEl,N,nVar)] = fb;
}

extern "C"
{
  void ScalarBoundaryInterp_1D_gpu_wrapper(real **bMatrix, real **f, real **fBound, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((ScalarBoundaryInterp_1D_gpu), dim3(nVar,nEl,1), dim3(2,1,1), 0, 0, *bMatrix, *f, *fBound, N, nVar);
	  ScalarBoundaryInterp_1D_gpu<<<dim3(nVar,nEl,1), dim3(2,1,1), 0, 0>>>(*bMatrix, *f, *fBound, N, nVar);
  } 
}
// ScalarBoundaryInterp_2D //
__global__ void ScalarBoundaryInterp_2D_gpu(real *bMatrix, real *f, real *fBound, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

  real fb[4] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
    fb[0] += f[SC_2D_INDEX(i,ii,iVar,iEl,N,nVar)]*bMatrix[ii]; // South
    fb[1] += f[SC_2D_INDEX(ii,i,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // East
    fb[2] += f[SC_2D_INDEX(i,ii,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // North
    fb[3] += f[SC_2D_INDEX(ii,i,iVar,iEl,N,nVar)]*bMatrix[ii]; // West
  }
  fBound[SCB_2D_INDEX(i,iVar,1,iEl,N,nVar)] = fb[0];
  fBound[SCB_2D_INDEX(i,iVar,2,iEl,N,nVar)] = fb[1];
  fBound[SCB_2D_INDEX(i,iVar,3,iEl,N,nVar)] = fb[2];
  fBound[SCB_2D_INDEX(i,iVar,4,iEl,N,nVar)] = fb[3];
}

extern "C"
{
  void ScalarBoundaryInterp_2D_gpu_wrapper(real **bMatrix, real **f, real **fBound, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((ScalarBoundaryInterp_2D_gpu), dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0, *bMatrix, *f, *fBound, N, nVar);
	  ScalarBoundaryInterp_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*bMatrix, *f, *fBound, N, nVar);
  } 
}

// VectorBoundaryInterp_2D //
__global__ void VectorBoundaryInterp_2D_gpu(real *bMatrix, real *f, real *fBound, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t idir = threadIdx.x+1;
  size_t i = threadIdx.y;

  real fb[4] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
      fb[0] += f[VE_2D_INDEX(idir,i,ii,iVar,iEl,N,nVar)]*bMatrix[ii]; // South
      fb[1] += f[VE_2D_INDEX(idir,ii,i,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // East
      fb[2] += f[VE_2D_INDEX(idir,i,ii,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // North
      fb[3] += f[VE_2D_INDEX(idir,ii,i,iVar,iEl,N,nVar)]*bMatrix[ii]; // West
  }
  fBound[VEB_2D_INDEX(idir,i,iVar,1,iEl,N,nVar)] = fb[0];
  fBound[VEB_2D_INDEX(idir,i,iVar,2,iEl,N,nVar)] = fb[1];
  fBound[VEB_2D_INDEX(idir,i,iVar,3,iEl,N,nVar)] = fb[2];
  fBound[VEB_2D_INDEX(idir,i,iVar,4,iEl,N,nVar)] = fb[3];
}

extern "C"
{
  void VectorBoundaryInterp_2D_gpu_wrapper(real **bMatrix, real **f, real **fBound, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorBoundaryInterp_2D_gpu), dim3(nVar,nEl,1), dim3(2,N+1,1), 0, 0, *bMatrix, *f, *fBound, N, nVar);
	  VectorBoundaryInterp_2D_gpu<<<dim3(nVar,nEl,1), dim3(2,N+1,1), 0, 0>>>(*bMatrix, *f, *fBound, N, nVar);
  } 
}

// TensorBoundaryInterp_2D //
__global__ void TensorBoundaryInterp_2D_gpu(real *bMatrix, real *f, real *fBound, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

  real fb[16] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
    fb[0] += f[TE_2D_INDEX(1,1,i,ii,iVar,iEl,N,nVar)]*bMatrix[ii]; // South
    fb[1] += f[TE_2D_INDEX(1,1,ii,i,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // East
    fb[2] += f[TE_2D_INDEX(1,1,i,ii,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // North
    fb[3] += f[TE_2D_INDEX(1,1,ii,i,iVar,iEl,N,nVar)]*bMatrix[ii]; // West

    fb[4] += f[TE_2D_INDEX(2,1,i,ii,iVar,iEl,N,nVar)]*bMatrix[ii]; // South
    fb[5] += f[TE_2D_INDEX(2,1,ii,i,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // East
    fb[6] += f[TE_2D_INDEX(2,1,i,ii,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // North
    fb[7] += f[TE_2D_INDEX(2,1,ii,i,iVar,iEl,N,nVar)]*bMatrix[ii]; // West

    fb[8] += f[TE_2D_INDEX(1,2,i,ii,iVar,iEl,N,nVar)]*bMatrix[ii]; // South
    fb[9] += f[TE_2D_INDEX(1,2,ii,i,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // East
    fb[10] += f[TE_2D_INDEX(1,2,i,ii,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // North
    fb[11] += f[TE_2D_INDEX(1,2,ii,i,iVar,iEl,N,nVar)]*bMatrix[ii]; // West

    fb[12] += f[TE_2D_INDEX(2,2,i,ii,iVar,iEl,N,nVar)]*bMatrix[ii]; // South
    fb[13] += f[TE_2D_INDEX(2,2,ii,i,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // East
    fb[14] += f[TE_2D_INDEX(2,2,i,ii,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // North
    fb[15] += f[TE_2D_INDEX(2,2,ii,i,iVar,iEl,N,nVar)]*bMatrix[ii]; // West
  }
  fBound[TEB_2D_INDEX(1,1,i,iVar,1,iEl,N,nVar)] = fb[0];
  fBound[TEB_2D_INDEX(1,1,i,iVar,2,iEl,N,nVar)] = fb[1];
  fBound[TEB_2D_INDEX(1,1,i,iVar,3,iEl,N,nVar)] = fb[2];
  fBound[TEB_2D_INDEX(1,1,i,iVar,4,iEl,N,nVar)] = fb[3];

  fBound[TEB_2D_INDEX(2,1,i,iVar,1,iEl,N,nVar)] = fb[4];
  fBound[TEB_2D_INDEX(2,1,i,iVar,2,iEl,N,nVar)] = fb[5];
  fBound[TEB_2D_INDEX(2,1,i,iVar,3,iEl,N,nVar)] = fb[6];
  fBound[TEB_2D_INDEX(2,1,i,iVar,4,iEl,N,nVar)] = fb[7];

  fBound[TEB_2D_INDEX(1,2,i,iVar,1,iEl,N,nVar)] = fb[8];
  fBound[TEB_2D_INDEX(1,2,i,iVar,2,iEl,N,nVar)] = fb[9];
  fBound[TEB_2D_INDEX(1,2,i,iVar,3,iEl,N,nVar)] = fb[10];
  fBound[TEB_2D_INDEX(1,2,i,iVar,4,iEl,N,nVar)] = fb[11];

  fBound[TEB_2D_INDEX(2,2,i,iVar,1,iEl,N,nVar)] = fb[12];
  fBound[TEB_2D_INDEX(2,2,i,iVar,2,iEl,N,nVar)] = fb[13];
  fBound[TEB_2D_INDEX(2,2,i,iVar,3,iEl,N,nVar)] = fb[14];
  fBound[TEB_2D_INDEX(2,2,i,iVar,4,iEl,N,nVar)] = fb[15];

}

extern "C"
{
  void TensorBoundaryInterp_2D_gpu_wrapper(real **bMatrix, real **f, real **fBound, int N, int nVar, int nEl)
  {
    //hipLaunchKernelGGL((TensorBoundaryInterp_2D_gpu), dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0, *bMatrix, *f, *fBound, N, nVar);
    TensorBoundaryInterp_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*bMatrix, *f, *fBound, N, nVar);
  }
}
// ScalarBoundaryInterp_3D //
__global__ void ScalarBoundaryInterp_3D_gpu(real *bMatrix, real *f, real *fBound, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real fb[6] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
    fb[0] += f[SC_3D_INDEX(i,j,ii,iVar,iEl,N,nVar)]*bMatrix[ii]; // Bottom
    fb[1] += f[SC_3D_INDEX(i,ii,j,iVar,iEl,N,nVar)]*bMatrix[ii]; // South
    fb[2] += f[SC_3D_INDEX(ii,i,j,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // East
    fb[3] += f[SC_3D_INDEX(i,ii,j,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // North
    fb[4] += f[SC_3D_INDEX(ii,i,j,iVar,iEl,N,nVar)]*bMatrix[ii]; // West
    fb[5] += f[SC_3D_INDEX(i,j,ii,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // Top
  }
  fBound[SCB_3D_INDEX(i,j,iVar,1,iEl,N,nVar)] = fb[0];
  fBound[SCB_3D_INDEX(i,j,iVar,2,iEl,N,nVar)] = fb[1];
  fBound[SCB_3D_INDEX(i,j,iVar,3,iEl,N,nVar)] = fb[2];
  fBound[SCB_3D_INDEX(i,j,iVar,4,iEl,N,nVar)] = fb[3];
  fBound[SCB_3D_INDEX(i,j,iVar,5,iEl,N,nVar)] = fb[4];
  fBound[SCB_3D_INDEX(i,j,iVar,6,iEl,N,nVar)] = fb[5];
}

extern "C"
{
  void ScalarBoundaryInterp_3D_gpu_wrapper(real **bMatrix, real **f, real **fBound, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((ScalarBoundaryInterp_3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *bMatrix, *f, *fBound, N, nVar);
	  ScalarBoundaryInterp_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*bMatrix, *f, *fBound, N, nVar);
  } 
}

// VectorBoundaryInterp_3D //
__global__ void VectorBoundaryInterp_3D_gpu(real *bMatrix, real *f, real *fBound, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t idir = threadIdx.x+1;
  size_t i = threadIdx.y;
  size_t j = threadIdx.z;

  real fb[6] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
    fb[0] += f[VE_3D_INDEX(idir,i,j,ii,iVar,iEl,N,nVar)]*bMatrix[ii]; // Bottom
    fb[1] += f[VE_3D_INDEX(idir,i,ii,j,iVar,iEl,N,nVar)]*bMatrix[ii]; // South
    fb[2] += f[VE_3D_INDEX(idir,ii,i,j,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // East
    fb[3] += f[VE_3D_INDEX(idir,i,ii,j,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // North
    fb[4] += f[VE_3D_INDEX(idir,ii,i,j,iVar,iEl,N,nVar)]*bMatrix[ii]; // West
    fb[5] += f[VE_3D_INDEX(idir,i,j,ii,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // Top
  }
  fBound[VEB_3D_INDEX(idir,i,j,iVar,1,iEl,N,nVar)] = fb[0];
  fBound[VEB_3D_INDEX(idir,i,j,iVar,2,iEl,N,nVar)] = fb[1];
  fBound[VEB_3D_INDEX(idir,i,j,iVar,3,iEl,N,nVar)] = fb[2];
  fBound[VEB_3D_INDEX(idir,i,j,iVar,4,iEl,N,nVar)] = fb[3];
  fBound[VEB_3D_INDEX(idir,i,j,iVar,5,iEl,N,nVar)] = fb[4];
  fBound[VEB_3D_INDEX(idir,i,j,iVar,6,iEl,N,nVar)] = fb[5];
}

extern "C"
{
  void VectorBoundaryInterp_3D_gpu_wrapper(real **bMatrix, real **f, real **fBound, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorBoundaryInterp_3D_gpu), dim3(nVar,nEl,1), dim3(3,N+1,N+1), 0, 0, *bMatrix, *f, *fBound, N, nVar);
	  VectorBoundaryInterp_3D_gpu<<<dim3(nVar,nEl,1), dim3(3,N+1,N+1), 0, 0>>>(*bMatrix, *f, *fBound, N, nVar);
  } 
}

// TensorBoundaryInterp_3D //
__global__ void TensorBoundaryInterp_3D_gpu(real *bMatrix, real *f, real *fBound, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  for (int col=1; col<=3; col++){
    for (int row=1; row<=3; row++){
      fBound[TEB_3D_INDEX(row,col,i,j,iVar,1,iEl,N,nVar)] = 0.0;
      fBound[TEB_3D_INDEX(row,col,i,j,iVar,2,iEl,N,nVar)] = 0.0;
      fBound[TEB_3D_INDEX(row,col,i,j,iVar,3,iEl,N,nVar)] = 0.0;
      fBound[TEB_3D_INDEX(row,col,i,j,iVar,4,iEl,N,nVar)] = 0.0;
      fBound[TEB_3D_INDEX(row,col,i,j,iVar,5,iEl,N,nVar)] = 0.0;
      fBound[TEB_3D_INDEX(row,col,i,j,iVar,6,iEl,N,nVar)] = 0.0;
    }
  }
  for (int ii=0; ii<N+1; ii++) {
    for (int col=1; col<=3; col++){
      for (int row=1; row<=3; row++){
        fBound[TEB_3D_INDEX(row,col,i,j,iVar,1,iEl,N,nVar)] += f[TE_3D_INDEX(row,col,i,j,ii,iVar,iEl,N,nVar)]*bMatrix[ii]; // Bottom
        fBound[TEB_3D_INDEX(row,col,i,j,iVar,2,iEl,N,nVar)] += f[TE_3D_INDEX(row,col,i,ii,j,iVar,iEl,N,nVar)]*bMatrix[ii]; // South
        fBound[TEB_3D_INDEX(row,col,i,j,iVar,3,iEl,N,nVar)] += f[TE_3D_INDEX(row,col,ii,i,j,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // East
        fBound[TEB_3D_INDEX(row,col,i,j,iVar,4,iEl,N,nVar)] += f[TE_3D_INDEX(row,col,i,ii,j,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // North
        fBound[TEB_3D_INDEX(row,col,i,j,iVar,5,iEl,N,nVar)] += f[TE_3D_INDEX(row,col,ii,i,j,iVar,iEl,N,nVar)]*bMatrix[ii]; // West
        fBound[TEB_3D_INDEX(row,col,i,j,iVar,6,iEl,N,nVar)] += f[TE_3D_INDEX(row,col,i,j,ii,iVar,iEl,N,nVar)]*bMatrix[ii+(N+1)]; // Top
      }
    }
  }

}

extern "C"
{
  void TensorBoundaryInterp_3D_gpu_wrapper(real **bMatrix, real **f, real **fBound, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((TensorBoundaryInterp_3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *bMatrix, *f, *fBound, N, nVar);
	  TensorBoundaryInterp_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*bMatrix, *f, *fBound, N, nVar);
  } 
}

// Derivative_1D //
__global__ void Derivative_1D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

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
	  //hipLaunchKernelGGL((Derivative_1D_gpu), dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0, *dMatrix, *f, *df, N, nVar);
	  Derivative_1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*dMatrix, *f, *df, N, nVar);
  } 
}

// DGDerivative_1D //
__global__ void DGDerivative_1D_gpu(real *dMatrix, real *bMatrix, real *qWeight, real *f, real *bf, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

  real fm = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    fm += f[SC_1D_INDEX(ii,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
  }

  fm += (bMatrix[i+(N+1)]*bf[SCB_1D_INDEX(iVar,1,iEl,N,nVar)]-
	 bMatrix[i]*bf[SCB_1D_INDEX(iVar,0,iEl,N,nVar)])/qWeight[i];

  df[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = fm;

}

extern "C"
{
  void DGDerivative_1D_gpu_wrapper(real **dMatrix, real **bMatrix, real **qWeight, real **f, real **bf, real **df, int N, int nVar, int nEl)
  {
//	  hipLaunchKernelGGL((DGDerivative_1D_gpu), dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0,
//                        *dMatrix, *bMatrix, *qWeight, *f, *bf, *df, N, nVar);

	  DGDerivative_1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*dMatrix, *bMatrix, *qWeight, *f, *bf, *df, N, nVar);
  }
}

// ScalarGradient_2D //
__global__ void ScalarGradient_2D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real fs = 0.0;
  real fp = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    fs += f[SC_2D_INDEX(ii,j,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
    fp += f[SC_2D_INDEX(i,ii,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)];
  }
  df[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] = fs;
  df[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] = fp;

}

extern "C"
{
  void ScalarGradient_2D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((ScalarGradient_2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *dMatrix, *f, *df, N, nVar);
	  ScalarGradient_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*dMatrix, *f, *df, N, nVar);
  } 
}

__global__ void ScalarDGGradient_2D_gpu(real *dgMatrix, real *bMatrix, real *qWeights, real *f, real *bf, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real fs = 0.0;
  real fp = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    fs += f[SC_2D_INDEX(ii,j,iVar,iEl,N,nVar)]*dgMatrix[ii+i*(N+1)];
    fp += f[SC_2D_INDEX(i,ii,iVar,iEl,N,nVar)]*dgMatrix[ii+j*(N+1)];
  }

  fs += (bf[SCB_2D_INDEX(j,iVar,2,iEl,N,nVar)]*bMatrix[i+(N+1)]-
         bf[SCB_2D_INDEX(j,iVar,4,iEl,N,nVar)]*bMatrix[i])/qWeights[i];

  fp += (bf[SCB_2D_INDEX(i,iVar,3,iEl,N,nVar)]*bMatrix[j+(N+1)]-
         bf[SCB_2D_INDEX(i,iVar,1,iEl,N,nVar)]*bMatrix[j])/qWeights[j];

  df[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] = fs;
  df[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] = fp;

}

extern "C"
{
  void ScalarDGGradient_2D_gpu_wrapper(real **dgMatrix, real **bMatrix, real **qWeights, real **f, real **bf, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((ScalarDGGradient_2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *dgMatrix, *bMatrix, *qWeights, *f, *bf, *df, N, nVar);
	  ScalarDGGradient_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*dgMatrix, *bMatrix, *qWeights, *f, *bf, *df, N, nVar);
  } 
}

// VectorGradient_2D //
__global__ void VectorGradient_2D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real dfloc[4] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
    dfloc[0] += f[VE_2D_INDEX(1,ii,j,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
    dfloc[1] += f[VE_2D_INDEX(2,ii,j,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
    dfloc[2] += f[VE_2D_INDEX(1,i,ii,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)];
    dfloc[3] += f[VE_2D_INDEX(2,i,ii,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)];
  }
  df[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)] = dfloc[0];
  df[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)] = dfloc[1];
  df[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)] = dfloc[2];
  df[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)] = dfloc[3];

}

extern "C"
{
  void VectorGradient_2D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorGradient_2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *dMatrix, *f, *df, N, nVar);
	  VectorGradient_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*dMatrix, *f, *df, N, nVar);
  } 
}

// VectorDGGradient_2D //
__global__ void VectorDGGradient_2D_gpu(real *dgMatrix, real *bMatrix, real *qWeights, real *f, real *bf, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real dfloc[4] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
    dfloc[0] += f[VE_2D_INDEX(1,ii,j,iVar,iEl,N,nVar)]*dgMatrix[ii+i*(N+1)];
    dfloc[1] += f[VE_2D_INDEX(2,ii,j,iVar,iEl,N,nVar)]*dgMatrix[ii+i*(N+1)];
    dfloc[2] += f[VE_2D_INDEX(1,i,ii,iVar,iEl,N,nVar)]*dgMatrix[ii+j*(N+1)];
    dfloc[3] += f[VE_2D_INDEX(2,i,ii,iVar,iEl,N,nVar)]*dgMatrix[ii+j*(N+1)];
  }

  dfloc[0] += (bf[VEB_2D_INDEX(1,j,iVar,2,iEl,N,nVar)]*bMatrix[i+(N+1)]-
               bf[VEB_2D_INDEX(1,j,iVar,4,iEl,N,nVar)]*bMatrix[i])/qWeights[i];
  dfloc[1] += (bf[VEB_2D_INDEX(2,j,iVar,2,iEl,N,nVar)]*bMatrix[i+(N+1)]-
               bf[VEB_2D_INDEX(2,j,iVar,4,iEl,N,nVar)]*bMatrix[i])/qWeights[i];
  dfloc[2] += (bf[VEB_2D_INDEX(1,i,iVar,3,iEl,N,nVar)]*bMatrix[j+(N+1)]-
               bf[VEB_2D_INDEX(1,i,iVar,1,iEl,N,nVar)]*bMatrix[j])/qWeights[j];
  dfloc[3] += (bf[VEB_2D_INDEX(2,i,iVar,3,iEl,N,nVar)]*bMatrix[j+(N+1)]-
               bf[VEB_2D_INDEX(2,i,iVar,1,iEl,N,nVar)]*bMatrix[j])/qWeights[j];

  df[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)] = dfloc[0];
  df[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)] = dfloc[1];
  df[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)] = dfloc[2];
  df[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)] = dfloc[3];

}

extern "C"
{
  void VectorDGGradient_2D_gpu_wrapper(real **dMatrix, real **bMatrix, real **qWeights, real **f, real **bf, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorDGGradient_2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *dMatrix, *bMatrix, *qWeights, *f, *bf, *df, N, nVar);
	  VectorDGGradient_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*dMatrix, *bMatrix, *qWeights, *f, *bf, *df, N, nVar);
  } 
}

// VectorDivergence_2D //
__global__ void VectorDivergence_2D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real dfloc = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    dfloc += f[VE_2D_INDEX(1,ii,j,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)] 
            +f[VE_2D_INDEX(2,i,ii,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)];
  }
  df[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = dfloc; 

}

extern "C"
{
  void VectorDivergence_2D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorDivergence_2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *dMatrix, *f, *df, N, nVar);
	  VectorDivergence_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*dMatrix, *f, *df, N, nVar);
  } 
}

// VectorDGDivergence_2D //
__global__ void VectorDGDivergence_2D_gpu(real *dgMatrix, real *bMatrix, real *qWeights, real *f, real *bf, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  df[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    df[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] += f[VE_2D_INDEX(1,ii,j,iVar,iEl,N,nVar)]*dgMatrix[ii+i*(N+1)]+
             f[VE_2D_INDEX(2,i,ii,iVar,iEl,N,nVar)]*dgMatrix[ii+j*(N+1)];
  }

  df[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] += (bf[VEB_2D_INDEX(1,j,iVar,2,iEl,N,nVar)]*bMatrix[i+(N+1)] +
            bf[VEB_2D_INDEX(1,j,iVar,4,iEl,N,nVar)]*bMatrix[i])/
           qWeights[i] +
           (bf[VEB_2D_INDEX(2,i,iVar,3,iEl,N,nVar)]*bMatrix[j+(N+1)] +
            bf[VEB_2D_INDEX(2,i,iVar,1,iEl,N,nVar)]*bMatrix[j])/
           qWeights[j];


}

extern "C"
{
  void VectorDGDivergence_2D_gpu_wrapper(real **dgMatrix, real **bMatrix, real **qWeights, real **f, real **bf, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorDGDivergence_2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *dgMatrix, *bMatrix, *qWeights, *f, *bf, *df, N, nVar);
	  VectorDGDivergence_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*dgMatrix, *bMatrix, *qWeights, *f, *bf, *df, N, nVar);
  } 
}

// VectorCurl_2D //
__global__ void VectorCurl_2D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real dfloc = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    dfloc += f[VE_2D_INDEX(2,i,ii,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)] 
            -f[VE_2D_INDEX(1,ii,j,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
  }
  df[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = dfloc; 

}

extern "C"
{
  void VectorCurl_2D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorCurl_2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *dMatrix, *f, *df, N, nVar);
	  VectorCurl_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*dMatrix, *f, *df, N, nVar);
  } 
}

// TensorDivergence_2D //
__global__ void TensorDivergence_2D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real df1 = 0.0;
  real df2 = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    df1 += f[TE_2D_INDEX(1,1,ii,j,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)] 
          +f[TE_2D_INDEX(2,1,i,ii,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)];

    df2 += f[TE_2D_INDEX(1,2,ii,j,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)] 
          +f[TE_2D_INDEX(2,2,i,ii,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)];
  }
  df[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] = df1; 
  df[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] = df2; 

}

extern "C"
{
  void TensorDivergence_2D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((TensorDivergence_2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *dMatrix, *f, *df, N, nVar);
	  TensorDivergence_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*dMatrix, *f, *df, N, nVar);
  } 
}

// TensorDGDivergence_2D //
__global__ void TensorDGDivergence_2D_gpu(real *dgMatrix, real *bMatrix, real *qWeights, real *f, real *bf, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  df[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] = 0.0;
  df[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    df[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] += f[TE_2D_INDEX(1,1,ii,j,iVar,iEl,N,nVar)]*dgMatrix[ii+i*(N+1)]+
           f[TE_2D_INDEX(2,1,i,ii,iVar,iEl,N,nVar)]*dgMatrix[ii+j*(N+1)];

    df[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] += f[TE_2D_INDEX(1,2,ii,j,iVar,iEl,N,nVar)]*dgMatrix[ii+i*(N+1)]+
           f[TE_2D_INDEX(2,2,i,ii,iVar,iEl,N,nVar)]*dgMatrix[ii+j*(N+1)];
  }

  df[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] += (bf[TEB_2D_INDEX(1,1,j,iVar,2,iEl,N,nVar)]*bMatrix[i+(N+1)]+
          bf[TEB_2D_INDEX(1,1,j,iVar,4,iEl,N,nVar)]*bMatrix[i])/qWeights[i]+
         (bf[TEB_2D_INDEX(2,1,i,iVar,3,iEl,N,nVar)]*bMatrix[j+(N+1)]+
          bf[TEB_2D_INDEX(2,1,i,iVar,1,iEl,N,nVar)]*bMatrix[j])/qWeights[j];

  df[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] += (bf[TEB_2D_INDEX(1,2,j,iVar,2,iEl,N,nVar)]*bMatrix[i+(N+1)]+
          bf[TEB_2D_INDEX(1,2,j,iVar,4,iEl,N,nVar)]*bMatrix[i])/qWeights[i]+
         (bf[TEB_2D_INDEX(2,2,i,iVar,3,iEl,N,nVar)]*bMatrix[j+(N+1)]+
          bf[TEB_2D_INDEX(2,2,i,iVar,1,iEl,N,nVar)]*bMatrix[j])/qWeights[j];

}

extern "C"
{
  void TensorDGDivergence_2D_gpu_wrapper(real **dgMatrix, real **bMatrix, real **qWeights, real **f, real **bf, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((TensorDGDivergence_2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *dgMatrix, *bMatrix, *qWeights, *f, *bf, *df, N, nVar);
	  TensorDGDivergence_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*dgMatrix, *bMatrix, *qWeights, *f, *bf, *df, N, nVar);
  }
}

// ScalarGradient_3D //
__global__ void ScalarGradient_3D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  real dfloc[3] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
    dfloc[0] += f[SC_3D_INDEX(ii,j,k,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
    dfloc[1] += f[SC_3D_INDEX(i,ii,k,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)];
    dfloc[2] += f[SC_3D_INDEX(i,j,ii,iVar,iEl,N,nVar)]*dMatrix[ii+k*(N+1)];
  }
  df[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] = dfloc[0];
  df[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] = dfloc[1];
  df[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] = dfloc[2];
}

extern "C"
{
  void ScalarGradient_3D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((ScalarGradient_3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *dMatrix, *f, *df, N, nVar);
	  ScalarGradient_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*dMatrix, *f, *df, N, nVar);
  } 
}

// VectorGradient_3D //
__global__ void VectorGradient_3D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  real dfloc[9] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
    dfloc[0] += f[VE_3D_INDEX(1,ii,j,k,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
    dfloc[1] += f[VE_3D_INDEX(2,ii,j,k,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
    dfloc[2] += f[VE_3D_INDEX(3,ii,j,k,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
    dfloc[3] += f[VE_3D_INDEX(1,i,ii,k,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)];
    dfloc[4] += f[VE_3D_INDEX(2,i,ii,k,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)];
    dfloc[5] += f[VE_3D_INDEX(3,i,ii,k,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)];
    dfloc[6] += f[VE_3D_INDEX(1,i,j,ii,iVar,iEl,N,nVar)]*dMatrix[ii+k*(N+1)];
    dfloc[7] += f[VE_3D_INDEX(2,i,j,ii,iVar,iEl,N,nVar)]*dMatrix[ii+k*(N+1)];
    dfloc[8] += f[VE_3D_INDEX(3,i,j,ii,iVar,iEl,N,nVar)]*dMatrix[ii+k*(N+1)];
  }
  df[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,N,nVar)] = dfloc[0];
  df[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,N,nVar)] = dfloc[1];
  df[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,N,nVar)] = dfloc[2];
  df[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)] = dfloc[3];
  df[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,N,nVar)] = dfloc[4];
  df[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)] = dfloc[5];
  df[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)] = dfloc[6];
  df[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)] = dfloc[7];
  df[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,N,nVar)] = dfloc[8];
}

extern "C"
{
  void VectorGradient_3D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorGradient_3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *dMatrix, *f, *df, N, nVar);
	  VectorGradient_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*dMatrix, *f, *df, N, nVar);
  } 
}

// VectorDivergence_3D //
__global__ void VectorDivergence_3D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  df[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)] = 0.0; 
  for (int ii=0; ii<N+1; ii++) {
    df[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)] += f[VE_3D_INDEX(1,ii,j,k,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)] 
            +f[VE_3D_INDEX(2,i,ii,k,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)]
            +f[VE_3D_INDEX(3,i,j,ii,iVar,iEl,N,nVar)]*dMatrix[ii+k*(N+1)];
  }

}

extern "C"
{
  void VectorDivergence_3D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorDivergence_3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *dMatrix, *f, *df, N, nVar);
	  VectorDivergence_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*dMatrix, *f, *df, N, nVar);
  } 
}

// VectorDGDivergence_3D //
__global__ void VectorDGDivergence_3D_gpu(real *dgMatrix, real *bMatrix, real *qWeights, real *f, real *bf, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  df[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)] = 0.0; 
  for (int ii=0; ii<N+1; ii++) {
    df[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)] += f[VE_3D_INDEX(1,ii,j,k,iVar,iEl,N,nVar)]*dgMatrix[ii+i*(N+1)]+
             f[VE_3D_INDEX(2,i,ii,k,iVar,iEl,N,nVar)]*dgMatrix[ii+j*(N+1)]+
             f[VE_3D_INDEX(3,i,j,ii,iVar,iEl,N,nVar)]*dgMatrix[ii+k*(N+1)];
  }
  df[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)] += (bf[VEB_3D_INDEX(1,j,k,iVar,3,iEl,N,nVar)]*bMatrix[i+(N+1)]+
            bf[VEB_3D_INDEX(1,j,k,iVar,5,iEl,N,nVar)]*bMatrix[i])/
           qWeights[i]+
           (bf[VEB_3D_INDEX(2,i,k,iVar,4,iEl,N,nVar)]*bMatrix[j+(N+1)]+
            bf[VEB_3D_INDEX(2,i,k,iVar,2,iEl,N,nVar)]*bMatrix[j])/
           qWeights[j]+
           (bf[VEB_3D_INDEX(3,i,j,iVar,6,iEl,N,nVar)]*bMatrix[k+(N+1)]+
            bf[VEB_3D_INDEX(3,i,j,iVar,1,iEl,N,nVar)]*bMatrix[k])/
           qWeights[k];


}

extern "C"
{
  void VectorDGDivergence_3D_gpu_wrapper(real **dgMatrix, real **bMatrix, real **qWeights, real **f, real **bf, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorDGDivergence_3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *dgMatrix, *bMatrix, *qWeights, *f, *bf, *df, N, nVar);
	  VectorDGDivergence_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*dgMatrix, *bMatrix, *qWeights, *f, *bf, *df, N, nVar);
  } 
}

// VectorCurl_3D //
__global__ void VectorCurl_3D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  real dfloc[3] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
    dfloc[0] += f[VE_3D_INDEX(3,i,ii,k,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)] 
               -f[VE_3D_INDEX(2,i,j,ii,iVar,iEl,N,nVar)]*dMatrix[ii+k*(N+1)];
    dfloc[1] += f[VE_3D_INDEX(1,i,j,ii,iVar,iEl,N,nVar)]*dMatrix[ii+k*(N+1)] 
               -f[VE_3D_INDEX(3,ii,j,k,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
    dfloc[2] += f[VE_3D_INDEX(2,i,ii,k,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)] 
               -f[VE_3D_INDEX(1,ii,j,k,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)];
  }
  df[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] = dfloc[0]; 
  df[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] = dfloc[1]; 
  df[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] = dfloc[2]; 

}

extern "C"
{
  void VectorCurl_3D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((VectorCurl_3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *dMatrix, *f, *df, N, nVar);
	  VectorCurl_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*dMatrix, *f, *df, N, nVar);
  } 
}

// TensorDivergence_3D //
__global__ void TensorDivergence_3D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  real df1 = 0.0;
  real df2 = 0.0;
  real df3 = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    df1 += f[TE_3D_INDEX(1,1,ii,j,k,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)] 
          +f[TE_3D_INDEX(2,1,i,ii,k,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)]
          +f[TE_3D_INDEX(3,1,i,j,ii,iVar,iEl,N,nVar)]*dMatrix[ii+k*(N+1)];

    df2 += f[TE_3D_INDEX(1,2,ii,j,k,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)] 
          +f[TE_3D_INDEX(2,2,i,ii,k,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)]
          +f[TE_3D_INDEX(3,2,i,j,ii,iVar,iEl,N,nVar)]*dMatrix[ii+k*(N+1)];

    df3 += f[TE_3D_INDEX(1,3,ii,j,k,iVar,iEl,N,nVar)]*dMatrix[ii+i*(N+1)] 
          +f[TE_3D_INDEX(2,3,i,ii,k,iVar,iEl,N,nVar)]*dMatrix[ii+j*(N+1)]
          +f[TE_3D_INDEX(3,3,i,j,ii,iVar,iEl,N,nVar)]*dMatrix[ii+k*(N+1)];
  }
  df[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] = df1; 
  df[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] = df2; 
  df[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] = df3; 

}

extern "C"
{
  void TensorDivergence_3D_gpu_wrapper(real **dMatrix, real **f, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((TensorDivergence_3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *dMatrix, *f, *df, N, nVar);
	  TensorDivergence_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*dMatrix, *f, *df, N, nVar);
  } 
}

// TensorDGDivergence_3D //
__global__ void TensorDGDivergence_3D_gpu(real *dgMatrix, real *bMatrix, real *qWeights, real *f, real *bf, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  df[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] = 0.0;
  df[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] = 0.0;
  df[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    df[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] += f[TE_3D_INDEX(1,1,ii,j,k,iVar,iEl,N,nVar)]*dgMatrix[ii+i*(N+1)] 
       +f[TE_3D_INDEX(2,1,i,ii,k,iVar,iEl,N,nVar)]*dgMatrix[ii+j*(N+1)]
       +f[TE_3D_INDEX(3,1,i,j,ii,iVar,iEl,N,nVar)]*dgMatrix[ii+k*(N+1)];

    df[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] += f[TE_3D_INDEX(1,2,ii,j,k,iVar,iEl,N,nVar)]*dgMatrix[ii+i*(N+1)] 
       +f[TE_3D_INDEX(2,2,i,ii,k,iVar,iEl,N,nVar)]*dgMatrix[ii+j*(N+1)]
       +f[TE_3D_INDEX(3,2,i,j,ii,iVar,iEl,N,nVar)]*dgMatrix[ii+k*(N+1)];

    df[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] += f[TE_3D_INDEX(1,3,ii,j,k,iVar,iEl,N,nVar)]*dgMatrix[ii+i*(N+1)] 
          +f[TE_3D_INDEX(2,3,i,ii,k,iVar,iEl,N,nVar)]*dgMatrix[ii+j*(N+1)]
          +f[TE_3D_INDEX(3,3,i,j,ii,iVar,iEl,N,nVar)]*dgMatrix[ii+k*(N+1)];
  }

  df[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] += (bf[TEB_3D_INDEX(1,1,j,k,iVar,3,iEl,N,nVar)]*bMatrix[i+(N+1)]+
          bf[TEB_3D_INDEX(1,1,j,k,iVar,5,iEl,N,nVar)]*bMatrix[i])/
         qWeights[i]+
         (bf[TEB_3D_INDEX(2,1,i,k,iVar,4,iEl,N,nVar)]*bMatrix[j+(N+1)]+
          bf[TEB_3D_INDEX(2,1,i,k,iVar,2,iEl,N,nVar)]*bMatrix[j])/
         qWeights[j]+
         (bf[TEB_3D_INDEX(3,1,i,j,iVar,6,iEl,N,nVar)]*bMatrix[k+(N+1)]+
          bf[TEB_3D_INDEX(3,1,i,j,iVar,1,iEl,N,nVar)]*bMatrix[k])/
         qWeights[k];

  df[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] += (bf[TEB_3D_INDEX(1,2,j,k,iVar,3,iEl,N,nVar)]*bMatrix[i+(N+1)]+
          bf[TEB_3D_INDEX(1,2,j,k,iVar,5,iEl,N,nVar)]*bMatrix[i])/
         qWeights[i]+
         (bf[TEB_3D_INDEX(2,2,i,k,iVar,4,iEl,N,nVar)]*bMatrix[j+(N+1)]+
          bf[TEB_3D_INDEX(2,2,i,k,iVar,2,iEl,N,nVar)]*bMatrix[j])/
         qWeights[j]+
         (bf[TEB_3D_INDEX(3,2,i,j,iVar,6,iEl,N,nVar)]*bMatrix[k+(N+1)]+
          bf[TEB_3D_INDEX(3,2,i,j,iVar,1,iEl,N,nVar)]*bMatrix[k])/
         qWeights[k];

  df[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] += (bf[TEB_3D_INDEX(1,3,j,k,iVar,3,iEl,N,nVar)]*bMatrix[i+(N+1)]+
          bf[TEB_3D_INDEX(1,3,j,k,iVar,5,iEl,N,nVar)]*bMatrix[i])/
         qWeights[i]+
         (bf[TEB_3D_INDEX(2,3,i,k,iVar,4,iEl,N,nVar)]*bMatrix[j+(N+1)]+
          bf[TEB_3D_INDEX(2,3,i,k,iVar,2,iEl,N,nVar)]*bMatrix[j])/
         qWeights[j]+
         (bf[TEB_3D_INDEX(3,3,i,j,iVar,6,iEl,N,nVar)]*bMatrix[k+(N+1)]+
          bf[TEB_3D_INDEX(3,3,i,j,iVar,1,iEl,N,nVar)]*bMatrix[k])/
         qWeights[k];


}

extern "C"
{
  void TensorDGDivergence_3D_gpu_wrapper(real **dgMatrix, real **bMatrix, real **qWeights, real **f, real **bf, real **df, int N, int nVar, int nEl)
  {
	  //hipLaunchKernelGGL((TensorDGDivergence_3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *dgMatrix, *bMatrix, *qWeights, *f, *bf, *df, N, nVar);
	  TensorDGDivergence_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*dgMatrix, *bMatrix, *qWeights, *f, *bf, *df, N, nVar);
  } 
}
