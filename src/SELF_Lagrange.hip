#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"



// ScalarBoundaryInterp_2D //
__global__ void ScalarBoundaryInterp_2D_gpu(real *bMatrix, real *f, real *fBound, int N){

  uint32_t iVar = blockIdx.y;
  uint32_t iEl = blockIdx.x;
  uint32_t nEl = gridDim.x;
  uint32_t i = threadIdx.x;

  real fb[4] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
    fb[0] += f[SC_2D_INDEX(i,ii,iEl,iVar,N,nEl)]*bMatrix[ii]; // South
    fb[1] += f[SC_2D_INDEX(ii,i,iEl,iVar,N,nEl)]*bMatrix[ii+(N+1)]; // East
    fb[2] += f[SC_2D_INDEX(i,ii,iEl,iVar,N,nEl)]*bMatrix[ii+(N+1)]; // North
    fb[3] += f[SC_2D_INDEX(ii,i,iEl,iVar,N,nEl)]*bMatrix[ii]; // West
  }
  fBound[SCB_2D_INDEX(i,0,iEl,iVar,N,nEl)] = fb[0];
  fBound[SCB_2D_INDEX(i,1,iEl,iVar,N,nEl)] = fb[1];
  fBound[SCB_2D_INDEX(i,2,iEl,iVar,N,nEl)] = fb[2];
  fBound[SCB_2D_INDEX(i,3,iEl,iVar,N,nEl)] = fb[3];
}

extern "C"
{
  void ScalarBoundaryInterp_2D_gpu_wrapper(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
	  ScalarBoundaryInterp_2D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,1,1), 0, 0>>>(bMatrix, f, fBound, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

// ScalarBoundaryInterp_3D //
__global__ void ScalarBoundaryInterp_3D_gpu(real *bMatrix, real *f, real *fBound, int N){

  uint32_t iVar = blockIdx.y;
  uint32_t iEl = blockIdx.x;
  uint32_t nEl = gridDim.x;
  uint32_t i = threadIdx.x;
  uint32_t j = threadIdx.y;

  real fb[6] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
    fb[0] += f[SC_3D_INDEX(i,j,ii,iEl,iVar,N,nEl)]*bMatrix[ii]; // Bottom
    fb[1] += f[SC_3D_INDEX(i,ii,j,iEl,iVar,N,nEl)]*bMatrix[ii]; // South
    fb[2] += f[SC_3D_INDEX(ii,i,j,iEl,iVar,N,nEl)]*bMatrix[ii+(N+1)]; // East
    fb[3] += f[SC_3D_INDEX(i,ii,j,iEl,iVar,N,nEl)]*bMatrix[ii+(N+1)]; // North
    fb[4] += f[SC_3D_INDEX(ii,i,j,iEl,iVar,N,nEl)]*bMatrix[ii]; // West
    fb[5] += f[SC_3D_INDEX(i,j,ii,iEl,iVar,N,nEl)]*bMatrix[ii+(N+1)]; // Top
  }
  fBound[SCB_3D_INDEX(i,j,0,iEl,iVar,N,nEl)] = fb[0];
  fBound[SCB_3D_INDEX(i,j,1,iEl,iVar,N,nEl)] = fb[1];
  fBound[SCB_3D_INDEX(i,j,2,iEl,iVar,N,nEl)] = fb[2];
  fBound[SCB_3D_INDEX(i,j,3,iEl,iVar,N,nEl)] = fb[3];
  fBound[SCB_3D_INDEX(i,j,4,iEl,iVar,N,nEl)] = fb[4];
  fBound[SCB_3D_INDEX(i,j,5,iEl,iVar,N,nEl)] = fb[5];
}

extern "C"
{
  void ScalarBoundaryInterp_3D_gpu_wrapper(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
	  ScalarBoundaryInterp_3D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(bMatrix, f, fBound, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}


// VectorBoundaryInterp_2D //
__global__ void VectorBoundaryInterp_2D_gpu(real *bMatrix, real *f, real *fBound, int N){

  uint32_t iVar = blockIdx.y;
  uint32_t nVar = gridDim.y;
  uint32_t iEl = blockIdx.x;
  uint32_t nEl = gridDim.x;
  uint32_t idir = blockIdx.z;
  uint32_t i = threadIdx.x;

  real fb[4] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
      fb[0] += f[VE_2D_INDEX(i,ii,iEl,iVar,idir,N,nEl,nVar)]*bMatrix[ii]; // South
      fb[1] += f[VE_2D_INDEX(ii,i,iEl,iVar,idir,N,nEl,nVar)]*bMatrix[ii+(N+1)]; // East
      fb[2] += f[VE_2D_INDEX(i,ii,iEl,iVar,idir,N,nEl,nVar)]*bMatrix[ii+(N+1)]; // North
      fb[3] += f[VE_2D_INDEX(ii,i,iEl,iVar,idir,N,nEl,nVar)]*bMatrix[ii]; // West
  }
  fBound[VEB_2D_INDEX(i,0,iEl,iVar,idir,N,nEl,nVar)] = fb[0];
  fBound[VEB_2D_INDEX(i,1,iEl,iVar,idir,N,nEl,nVar)] = fb[1];
  fBound[VEB_2D_INDEX(i,2,iEl,iVar,idir,N,nEl,nVar)] = fb[2];
  fBound[VEB_2D_INDEX(i,3,iEl,iVar,idir,N,nEl,nVar)] = fb[3];
}

extern "C"
{
  void VectorBoundaryInterp_2D_gpu_wrapper(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
	  VectorBoundaryInterp_2D_gpu<<<dim3(nEl,nVar,2), dim3(N+1,1,1), 0, 0>>>(bMatrix, f, fBound, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}


// VectorBoundaryInterp_3D //
__global__ void VectorBoundaryInterp_3D_gpu(real *bMatrix, real *f, real *fBound, int N){

  uint32_t iVar = blockIdx.y;
  uint32_t nVar = gridDim.y;
  uint32_t iEl = blockIdx.x;
  uint32_t nEl = gridDim.x;
  uint32_t idir = blockIdx.z;
  uint32_t i = threadIdx.x;
  uint32_t j = threadIdx.y;

  real fb[6] = {0.0};
  for (int ii=0; ii<N+1; ii++) {
      fb[0] += f[VE_3D_INDEX(i,j,ii,iEl,iVar,idir,N,nEl,nVar)]*bMatrix[ii]; // Bottom
      fb[1] += f[VE_3D_INDEX(i,ii,j,iEl,iVar,idir,N,nEl,nVar)]*bMatrix[ii]; // South
      fb[2] += f[VE_3D_INDEX(ii,i,j,iEl,iVar,idir,N,nEl,nVar)]*bMatrix[ii+(N+1)]; // East
      fb[3] += f[VE_3D_INDEX(i,ii,j,iEl,iVar,idir,N,nEl,nVar)]*bMatrix[ii+(N+1)]; // North
      fb[4] += f[VE_3D_INDEX(ii,i,j,iEl,iVar,idir,N,nEl,nVar)]*bMatrix[ii]; // West
      fb[5] += f[VE_3D_INDEX(i,j,ii,iEl,iVar,idir,N,nEl,nVar)]*bMatrix[ii+(N+1)]; // Top
  }
  fBound[VEB_3D_INDEX(i,j,0,iEl,iVar,idir,N,nEl,nVar)] = fb[0];
  fBound[VEB_3D_INDEX(i,j,1,iEl,iVar,idir,N,nEl,nVar)] = fb[1];
  fBound[VEB_3D_INDEX(i,j,2,iEl,iVar,idir,N,nEl,nVar)] = fb[2];
  fBound[VEB_3D_INDEX(i,j,3,iEl,iVar,idir,N,nEl,nVar)] = fb[3];
  fBound[VEB_3D_INDEX(i,j,4,iEl,iVar,idir,N,nEl,nVar)] = fb[4];
  fBound[VEB_3D_INDEX(i,j,5,iEl,iVar,idir,N,nEl,nVar)] = fb[5];
}

extern "C"
{
  void VectorBoundaryInterp_3D_gpu_wrapper(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
	  VectorBoundaryInterp_3D_gpu<<<dim3(nEl,nVar,3), dim3(N+1,N+1,1), 0, 0>>>(bMatrix, f, fBound, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

extern "C"
{
  void TensorBoundaryInterp_2D_gpu_wrapper(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
    // Intentionally using the vector boundary interpolation method with the 
    // gridDim.z = 4 to account for the 2x2 structure at each point
    VectorBoundaryInterp_2D_gpu<<<dim3(nVar,nEl,4), dim3(N+1,1,1), 0, 0>>>(bMatrix, f, fBound, N);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

extern "C"
{
  void TensorBoundaryInterp_3D_gpu_wrapper(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
    // Intentionally using the vector boundary interpolation method with the 
    // gridDim.z = 9 to account for the 3x3 structure at each point
	  VectorBoundaryInterp_3D_gpu<<<dim3(nVar,nEl,9), dim3(N+1,N+1,1), 0, 0>>>(bMatrix, f, fBound, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

__global__ void DGDerivative_BoundaryContribution_1D(real *bMatrix, real *qWeight, real *bf, real *df, int N, int nEl){

  
  uint32_t ivar = blockIdx.y;
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;

  uint32_t ndof=(N+1)*nEl;
  if( idof < ndof ){
    uint32_t iel = idof/(N+1); // Calculate the element ID
    uint32_t i = idof - iel*(N+1); // Calculate the quadrature point node id
    df[idof + ndof*ivar] += (bMatrix[i+(N+1)]*bf[SCB_1D_INDEX(1,iel,ivar,N,nEl)]+
	                           bMatrix[i]*bf[SCB_1D_INDEX(0,iel,ivar,N,nEl)])/qWeight[i];
  }

}

extern "C"
{
  void DGDerivative_BoundaryContribution_1D_gpu(real *bMatrix, real *qWeight,real *bf, real *df, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*nEl;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
	  DGDerivative_BoundaryContribution_1D<<<dim3(nblocks_x,nVar,1), dim3(threads_per_block,1,1), 0, 0>>>(bMatrix, qWeight, bf, df, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// VectorDGDivergence_2D //
__global__ void VectorDGDivergence_BoundaryContribution_2D(real *bMatrix, real *qWeights, real *bf, real *df, int N){

  uint32_t nel = blockDim.x;
  uint32_t ivar = blockIdx.y;
  uint32_t iel = blockIdx.x;
  uint32_t i = threadIdx.x;
  uint32_t j = threadIdx.y;

  df[SC_2D_INDEX(i,j,iel,ivar,N,nel)] += (bf[SCB_2D_INDEX(j,1,iel,ivar,N,nel)]*bMatrix[i+(N+1)] + // east
                                          bf[SCB_2D_INDEX(j,3,iel,ivar,N,nel)]*bMatrix[i])/       // west
                                         qWeights[i] +
                                         (bf[SCB_2D_INDEX(i,2,iel,ivar,N,nel)]*bMatrix[j+(N+1)] + // north
                                          bf[SCB_2D_INDEX(i,0,iel,ivar,N,nel)]*bMatrix[j])/       // south
                                         qWeights[j];


}

extern "C"
{
  void VectorDGDivergence_BoundaryContribution_2D_gpu(real *bMatrix, real *qWeights, real *bf, real *df, int N, int nVar, int nEl)
  {
	  VectorDGDivergence_BoundaryContribution_2D<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(bMatrix, qWeights, bf, df, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

__global__ void VectorDGDivergence_BoundaryContribution_3D(real *bMatrix, real *qWeights, real *bf, real *df, int N){

  uint32_t nel = blockDim.y;
  uint32_t iq = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t iel = blockIdx.y;
  uint32_t ivar = blockIdx.z;

  uint32_t i = iq % (N+1);
  uint32_t j = (iq/(N+1)) % (N+1);
  uint32_t k = iq/(N+1)/(N+1);

  if( iq < (N+1)*(N+1)*(N+1) ){

    df[SC_3D_INDEX(i,j,k,iel,ivar,N,nel)] += (bf[SCB_3D_INDEX(j,k,2,iel,ivar,N,nel)]*bMatrix[i+(N+1)] +
                                              bf[SCB_3D_INDEX(j,k,4,iel,ivar,N,nel)]*bMatrix[i])/
                                            qWeights[i] +
                                            (bf[SCB_3D_INDEX(i,k,3,iel,ivar,N,nel)]*bMatrix[j+(N+1)] +
                                              bf[SCB_3D_INDEX(i,k,1,iel,ivar,N,nel)]*bMatrix[j])/
                                            qWeights[j] +
                                            (bf[SCB_3D_INDEX(i,j,5,iel,ivar,N,nel)]*bMatrix[k+(N+1)] +
                                              bf[SCB_3D_INDEX(i,j,0,iel,ivar,N,nel)]*bMatrix[k])/
                                            qWeights[k];
  }


}

extern "C"
{
  void VectorDGDivergence_BoundaryContribution_3D_gpu(real *bMatrix, real *qWeights, real *bf, real *df, int N, int nVar, int nEl)
  {
    uint32_t nquad_pts = (N+1)*(N+1)*(N+1);
    uint32_t nthreads = 256;
    uint32_t nblocks_x = nquad_pts/nthreads + 1;

	  VectorDGDivergence_BoundaryContribution_3D<<<dim3(nblocks_x,nEl,nVar), dim3(nthreads,1,1), 0, 0>>>(bMatrix, qWeights, bf, df, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

__global__ void TensorDGDivergence_BoundaryContribution_2D(real *bMatrix, real *qWeights, real *bf, real *df, int N){

  uint32_t nel = blockDim.x;
  uint32_t nvar = blockDim.y;

  uint32_t iel = blockIdx.x;
  uint32_t ivar = blockIdx.y;

  uint32_t i = threadIdx.x;
  uint32_t j = threadIdx.y;

  df[VE_2D_INDEX(i,j,iel,ivar,0,N,nel,nvar)] += (bf[TEB_2D_INDEX(j,1,iel,ivar,0,0,N,nel,nvar)]*bMatrix[i+(N+1)] + // east
                                                    bf[TEB_2D_INDEX(j,3,iel,ivar,0,0,N,nel,nvar)]*bMatrix[i])/    // west
                                                   qWeights[i] +
                                                   (bf[TEB_2D_INDEX(i,2,iel,ivar,0,1,N,nel,nvar)]*bMatrix[j+(N+1)] + // north
                                                    bf[TEB_2D_INDEX(i,0,iel,ivar,0,1,N,nel,nvar)]*bMatrix[j])/       // south
                                                   qWeights[j];

  df[VE_2D_INDEX(i,j,iel,ivar,1,N,nel,nvar)] += (bf[TEB_2D_INDEX(j,1,iel,ivar,1,0,N,nel,nvar)]*bMatrix[i+(N+1)] + // east
                                                    bf[TEB_2D_INDEX(j,3,iel,ivar,1,0,N,nel,nvar)]*bMatrix[i])/    // west
                                                   qWeights[i] +
                                                   (bf[TEB_2D_INDEX(i,2,iel,ivar,1,1,N,nel,nvar)]*bMatrix[j+(N+1)] + // north
                                                    bf[TEB_2D_INDEX(i,0,iel,ivar,1,1,N,nel,nvar)]*bMatrix[j])/       // south
                                                   qWeights[j];

}

extern "C"
{
  void TensorDGDivergence_BoundaryContribution_2D_gpu(real *bMatrix, real *qWeights, real *bf, real *df, int N, int nVar, int nEl)
  {
	  TensorDGDivergence_BoundaryContribution_2D<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(bMatrix, qWeights, bf, df, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

__global__ void TensorDGDivergence_BoundaryContribution_3D(real *bMatrix, real *qWeights, real *bf, real *df, int N){

  uint32_t nel = blockDim.y;
  uint32_t nvar = blockDim.z;

  uint32_t iq = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t iel = blockIdx.y;
  uint32_t ivar = blockIdx.z;

  uint32_t i = iq % (N+1);
  uint32_t j = (iq/(N+1)) % (N+1);
  uint32_t k = iq/(N+1)/(N+1);

  if( iq < (N+1)*(N+1)*(N+1) ){
    df[VE_3D_INDEX(i,j,k,iel,ivar,0,N,nel,nvar)] += (bf[TEB_3D_INDEX(j,k,2,iel,ivar,0,0,N,nel,nvar)]*bMatrix[i+(N+1)] + // east
                                                      bf[TEB_3D_INDEX(j,k,4,iel,ivar,0,0,N,nel,nvar)]*bMatrix[i])/      // west
                                                    qWeights[i] +
                                                    (bf[TEB_3D_INDEX(i,k,3,iel,ivar,0,1,N,nel,nvar)]*bMatrix[j+(N+1)] + // north
                                                      bf[TEB_3D_INDEX(i,k,1,iel,ivar,0,1,N,nel,nvar)]*bMatrix[j])/      // south
                                                    qWeights[j] +
                                                    (bf[TEB_3D_INDEX(i,j,5,iel,ivar,0,2,N,nel,nvar)]*bMatrix[k+(N+1)] + // top
                                                      bf[TEB_3D_INDEX(i,j,0,iel,ivar,0,2,N,nel,nvar)]*bMatrix[k])/      // bottom
                                                    qWeights[k];

    df[VE_3D_INDEX(i,j,k,iel,ivar,1,N,nel,nvar)] += (bf[TEB_3D_INDEX(j,k,2,iel,ivar,1,0,N,nel,nvar)]*bMatrix[i+(N+1)] + // east
                                                      bf[TEB_3D_INDEX(j,k,4,iel,ivar,1,0,N,nel,nvar)]*bMatrix[i])/      // west
                                                    qWeights[i] +
                                                    (bf[TEB_3D_INDEX(i,k,3,iel,ivar,1,1,N,nel,nvar)]*bMatrix[j+(N+1)] + // north
                                                      bf[TEB_3D_INDEX(i,k,1,iel,ivar,1,1,N,nel,nvar)]*bMatrix[j])/      // south
                                                    qWeights[j] +
                                                    (bf[TEB_3D_INDEX(i,j,5,iel,ivar,1,2,N,nel,nvar)]*bMatrix[k+(N+1)] + // top
                                                      bf[TEB_3D_INDEX(i,j,0,iel,ivar,1,2,N,nel,nvar)]*bMatrix[k])/      // bottom
                                                    qWeights[k];

    df[VE_3D_INDEX(i,j,k,iel,ivar,2,N,nel,nvar)] += (bf[TEB_3D_INDEX(j,k,2,iel,ivar,2,0,N,nel,nvar)]*bMatrix[i+(N+1)] + // east
                                                      bf[TEB_3D_INDEX(j,k,4,iel,ivar,2,0,N,nel,nvar)]*bMatrix[i])/      // west
                                                    qWeights[i] +
                                                    (bf[TEB_3D_INDEX(i,k,3,iel,ivar,2,1,N,nel,nvar)]*bMatrix[j+(N+1)] + // north
                                                      bf[TEB_3D_INDEX(i,k,1,iel,ivar,2,1,N,nel,nvar)]*bMatrix[j])/      // south
                                                    qWeights[j] +
                                                    (bf[TEB_3D_INDEX(i,j,5,iel,ivar,2,2,N,nel,nvar)]*bMatrix[k+(N+1)] + // top
                                                      bf[TEB_3D_INDEX(i,j,0,iel,ivar,2,2,N,nel,nvar)]*bMatrix[k])/      // bottom
                                                    qWeights[k];
  }

}

extern "C"
{
  void TensorDGDivergence_BoundaryContribution_3D_gpu(real *bMatrix, real *qWeights, real *bf, real *df, int N, int nVar, int nEl)
  {
    uint32_t nquad_pts = (N+1)*(N+1)*(N+1);
    uint32_t nthreads = 256;
    uint32_t nblocks_x = nquad_pts/nthreads + 1;
	  TensorDGDivergence_BoundaryContribution_3D<<<dim3(nblocks_x,nEl,nVar), dim3(nthreads,1,1), 0, 0>>>(bMatrix, qWeights, bf, df, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}