#include "SELF_GPU_Macros.h"

__global__ void Average(real *avgf, real *f1, real *f2, int ndof){

  uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;

  if( i < ndof ){
    avgf[i] =0.5*(f1[i]+f2[i]);
  }
}

extern "C"
{
  void Average_gpu(real *f, real *f1, real *f2, int ndof)
  {
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,1,1);
    dim3 nthreads(threads_per_block,1,1);

    Average<<<nblocks, nthreads, 0, 0>>>(f, f1, f2, ndof);
  }
}

__global__ void BoundaryInterp_2D_gpukernel(real *bMatrix, real *f, real *fBound, int N, int nEl){

  int iq = threadIdx.x + blockIdx.x*blockDim.x;
  int i = iq % (N+1);
  int iEl = (iq/(N+1)) % (nEl);
  int iVar = iq/(N+1)/(nEl);

  real fbl = 0.0;
  real fbr = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    fbl += f[SC_2D_INDEX(i,ii,iEl,iVar,N,nEl)]*bMatrix[ii]; // South
    fbr += f[SC_2D_INDEX(i,ii,iEl,iVar,N,nEl)]*bMatrix[ii+(N+1)]; // North
  }
  fBound[SCB_2D_INDEX(i,0,iEl,iVar,N,nEl)] = fbl; // South
  fBound[SCB_2D_INDEX(i,2,iEl,iVar,N,nEl)] = fbr; // North


  fbl = 0.0;
  fbr = 0.0;
  for (int ii=0; ii<N+1; ii++) {
    fbl += f[SC_2D_INDEX(ii,i,iEl,iVar,N,nEl)]*bMatrix[ii]; // West
    fbr += f[SC_2D_INDEX(ii,i,iEl,iVar,N,nEl)]*bMatrix[ii+(N+1)]; // East
  }
  fBound[SCB_2D_INDEX(i,3,iEl,iVar,N,nEl)] = fbl; // West
  fBound[SCB_2D_INDEX(i,1,iEl,iVar,N,nEl)] = fbr; // East

}

extern "C"
{
  void BoundaryInterp_2D_gpu(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*nEl*nVar;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
	  BoundaryInterp_2D_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(bMatrix, f, fBound, N, nEl);
  } 
}

__global__ void BoundaryInterp_3D_gpukernel(real *bMatrix, real *f, real *fBound, int N, int nEl){
  
  int iq = threadIdx.x + blockIdx.x*blockDim.x;
  int i = iq % (N+1);
  int j = (iq/(N+1))%(N+1);
  int iEl = (iq/(N+1)/(N+1)) % (nEl);
  int iVar = iq/(N+1)/(N+1)/(nEl);

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
  void BoundaryInterp_3D_gpu(real *bMatrix, real *f, real *fBound, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(N+1)*nEl*nVar;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
	  BoundaryInterp_3D_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(bMatrix, f, fBound, N, nEl);
  } 
}
