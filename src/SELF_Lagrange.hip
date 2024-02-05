#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"



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
  void ScalarBoundaryInterp_2D_gpu_wrapper(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
	  ScalarBoundaryInterp_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(bMatrix, f, fBound, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  void VectorBoundaryInterp_2D_gpu_wrapper(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
	  VectorBoundaryInterp_2D_gpu<<<dim3(nVar,nEl,1), dim3(2,N+1,1), 0, 0>>>(bMatrix, f, fBound, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  void TensorBoundaryInterp_2D_gpu_wrapper(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
    TensorBoundaryInterp_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(bMatrix, f, fBound, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  void ScalarBoundaryInterp_3D_gpu_wrapper(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
	  ScalarBoundaryInterp_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(bMatrix, f, fBound, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  void VectorBoundaryInterp_3D_gpu_wrapper(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
	  VectorBoundaryInterp_3D_gpu<<<dim3(nVar,nEl,1), dim3(3,N+1,N+1), 0, 0>>>(bMatrix, f, fBound, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  void TensorBoundaryInterp_3D_gpu_wrapper(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
	  TensorBoundaryInterp_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(bMatrix, f, fBound, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}


// DGDerivative_1D //
__global__ void DGDerivative_1D_gpu(real *dMatrix, real *bMatrix, real *qWeight, real *f, real *bf, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

  df[SC_1D_INDEX(i,iVar,iEl,N,nVar)] += (bMatrix[i+(N+1)]*bf[SCB_1D_INDEX(iVar,1,iEl,N,nVar)]+
	 bMatrix[i]*bf[SCB_1D_INDEX(iVar,0,iEl,N,nVar)])/qWeight[i];

}

extern "C"
{
  void DGDerivative_1D_gpu_wrapper(real *dMatrix, real *bMatrix, real *qWeight, real *f, real *bf, real *df, int N, int nVar, int nEl)
  {
	  DGDerivative_1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(dMatrix, bMatrix, qWeight, f, bf, df, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  void ScalarDGGradient_2D_gpu_wrapper(real *dgMatrix, real *bMatrix, real *qWeights, real *f, real *bf, real *df, int N, int nVar, int nEl)
  {
	  ScalarDGGradient_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(dgMatrix, bMatrix, qWeights, f, bf, df, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  void VectorGradient_2D_gpu_wrapper(real *dMatrix, real *f, real *df, int N, int nVar, int nEl)
  {
	  VectorGradient_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(dMatrix, f, df, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  void VectorDGGradient_2D_gpu_wrapper(real *dMatrix, real *bMatrix, real *qWeights, real *f, real *bf, real *df, int N, int nVar, int nEl)
  {
	  VectorDGGradient_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(dMatrix, bMatrix, qWeights, f, bf, df, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  void VectorDivergence_2D_gpu_wrapper(real *dMatrix, real *f, real *df, int N, int nVar, int nEl)
  {
	  VectorDivergence_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(dMatrix, f, df, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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

  df[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] += (bf[SCB_2D_INDEX(j,iVar,2,iEl,N,nVar)]*bMatrix[i+(N+1)] +
            bf[SCB_2D_INDEX(j,iVar,4,iEl,N,nVar)]*bMatrix[i])/
           qWeights[i] +
           (bf[SCB_2D_INDEX(i,iVar,3,iEl,N,nVar)]*bMatrix[j+(N+1)] +
            bf[SCB_2D_INDEX(i,iVar,1,iEl,N,nVar)]*bMatrix[j])/
           qWeights[j];


}

extern "C"
{
  void VectorDGDivergence_2D_gpu_wrapper(real *dgMatrix, real *bMatrix, real *qWeights, real *f, real *bf, real *df, int N, int nVar, int nEl)
  {
	  VectorDGDivergence_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(dgMatrix, bMatrix, qWeights, f, bf, df, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

// P2VectorDivergence_2D //
__global__ void P2VectorDivergence_2D_gpu(real *dMatrix, real *f, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real dfloc = 0.0;
  for (int n=0; n<N+1; n++) {
    dfloc += f[P2VE_2D_INDEX(1,n,i,j,iVar,iEl,N,nVar)]*dMatrix[n+i*(N+1)]+ 
             f[P2VE_2D_INDEX(2,n,i,j,iVar,iEl,N,nVar)]*dMatrix[n+j*(N+1)];
  }
  df[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = dfloc; 

}

extern "C"
{
  void P2VectorDivergence_2D_gpu_wrapper(real *dMatrix, real *f, real *df, int N, int nVar, int nEl)
  {
	  P2VectorDivergence_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(dMatrix, f, df, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

// VectorDGDivergence_2D //
__global__ void P2VectorDGDivergence_2D_gpu(real *dgMatrix, real *bMatrix, real *qWeights, real *f, real *bf, real *df, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real dfloc = 0.0;
  for (int n =0; n<N+1; n++) {
    dfloc += f[P2VE_2D_INDEX(1,n,i,j,iVar,iEl,N,nVar)]*dgMatrix[n+i*(N+1)]+
             f[P2VE_2D_INDEX(2,n,i,j,iVar,iEl,N,nVar)]*dgMatrix[n+j*(N+1)];
  }

  dfloc += (bf[SCB_2D_INDEX(j,iVar,2,iEl,N,nVar)]*bMatrix[i+(N+1)] +
            bf[SCB_2D_INDEX(j,iVar,4,iEl,N,nVar)]*bMatrix[i])/
           qWeights[i] +
           (bf[SCB_2D_INDEX(i,iVar,3,iEl,N,nVar)]*bMatrix[j+(N+1)] +
            bf[SCB_2D_INDEX(i,iVar,1,iEl,N,nVar)]*bMatrix[j])/
           qWeights[j];

  df[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = dfloc; 

}

extern "C"
{
  void P2VectorDGDivergence_2D_gpu_wrapper(real *dgMatrix, real *bMatrix, real *qWeights, real *f, real *bf, real *df, int N, int nVar, int nEl)
  {
	  P2VectorDGDivergence_2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(dgMatrix, bMatrix, qWeights, f, bf, df, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  void ScalarGradient_3D_gpu_wrapper(real *dMatrix, real *f, real *df, int N, int nVar, int nEl)
  {
	  ScalarGradient_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(dMatrix, f, df, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  void VectorGradient_3D_gpu_wrapper(real *dMatrix, real *f, real *df, int N, int nVar, int nEl)
  {
	  VectorGradient_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(dMatrix, f, df, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  void VectorDivergence_3D_gpu_wrapper(real *dMatrix, real *f, real *df, int N, int nVar, int nEl)
  {
	  VectorDivergence_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(dMatrix, f, df, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
  df[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)] += (bf[SCB_3D_INDEX(j,k,iVar,3,iEl,N,nVar)]*bMatrix[i+(N+1)]+
            bf[SCB_3D_INDEX(j,k,iVar,5,iEl,N,nVar)]*bMatrix[i])/
           qWeights[i]+
           (bf[SCB_3D_INDEX(i,k,iVar,4,iEl,N,nVar)]*bMatrix[j+(N+1)]+
            bf[SCB_3D_INDEX(i,k,iVar,2,iEl,N,nVar)]*bMatrix[j])/
           qWeights[j]+
           (bf[SCB_3D_INDEX(i,j,iVar,6,iEl,N,nVar)]*bMatrix[k+(N+1)]+
            bf[SCB_3D_INDEX(i,j,iVar,1,iEl,N,nVar)]*bMatrix[k])/
           qWeights[k];


}

extern "C"
{
  void VectorDGDivergence_3D_gpu_wrapper(real *dgMatrix, real *bMatrix, real *qWeights, real *f, real *bf, real *df, int N, int nVar, int nEl)
  {
	  VectorDGDivergence_3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(dgMatrix, bMatrix, qWeights, f, bf, df, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}