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

__global__ void BoundaryInterp_2D_gpukernel(real *bMatrix, real *f, real *fBound, int N, int nel, int nvar){
  int ndof = (N+1)*nel*nvar;
  int iq = threadIdx.x + blockIdx.x*blockDim.x;
  if(iq < ndof){
    int i = iq % (N+1);
    int iEl = (iq/(N+1)) % (nel);
    int iVar = iq/(N+1)/(nel);

    real fbl = 0.0;
    real fbr = 0.0;
    for (int ii=0; ii<N+1; ii++) {
      fbl += f[SC_2D_INDEX(i,ii,iEl,iVar,N,nel)]*bMatrix[ii]; // South
      fbr += f[SC_2D_INDEX(i,ii,iEl,iVar,N,nel)]*bMatrix[ii+(N+1)]; // North
    }
    fBound[SCB_2D_INDEX(i,0,iEl,iVar,N,nel)] = fbl; // South
    fBound[SCB_2D_INDEX(i,2,iEl,iVar,N,nel)] = fbr; // North


    fbl = 0.0;
    fbr = 0.0;
    for (int ii=0; ii<N+1; ii++) {
      fbl += f[SC_2D_INDEX(ii,i,iEl,iVar,N,nel)]*bMatrix[ii]; // West
      fbr += f[SC_2D_INDEX(ii,i,iEl,iVar,N,nel)]*bMatrix[ii+(N+1)]; // East
    }
    fBound[SCB_2D_INDEX(i,3,iEl,iVar,N,nel)] = fbl; // West
    fBound[SCB_2D_INDEX(i,1,iEl,iVar,N,nel)] = fbr; // East
  }

}

extern "C"
{
  void BoundaryInterp_2D_gpu(real *bMatrix, real *f, real *fBound, int N, int nvar, int nel)
  {
    int ndof = (N+1)*nel*nvar;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
	  BoundaryInterp_2D_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(bMatrix, f, fBound, N, nel, nvar);
  } 
}

__global__ void BoundaryInterp_3D_gpukernel(real *bMatrix, real *f, real *fBound, int N, int nel, int nvar){
  int ndof = (N+1)*(N+1)*nel*nvar;
  int iq = threadIdx.x + blockIdx.x*blockDim.x;
  if( iq < ndof ){
    int i = iq % (N+1);
    int j = (iq/(N+1))%(N+1);
    int iEl = (iq/(N+1)/(N+1)) % (nel);
    int iVar = iq/(N+1)/(N+1)/(nel);

    real fb[6] = {0.0};
    for (int ii=0; ii<N+1; ii++) {
      fb[0] += f[SC_3D_INDEX(i,j,ii,iEl,iVar,N,nel)]*bMatrix[ii]; // Bottom
      fb[1] += f[SC_3D_INDEX(i,ii,j,iEl,iVar,N,nel)]*bMatrix[ii]; // South
      fb[2] += f[SC_3D_INDEX(ii,i,j,iEl,iVar,N,nel)]*bMatrix[ii+(N+1)]; // East
      fb[3] += f[SC_3D_INDEX(i,ii,j,iEl,iVar,N,nel)]*bMatrix[ii+(N+1)]; // North
      fb[4] += f[SC_3D_INDEX(ii,i,j,iEl,iVar,N,nel)]*bMatrix[ii]; // West
      fb[5] += f[SC_3D_INDEX(i,j,ii,iEl,iVar,N,nel)]*bMatrix[ii+(N+1)]; // Top
    }
    fBound[SCB_3D_INDEX(i,j,0,iEl,iVar,N,nel)] = fb[0];
    fBound[SCB_3D_INDEX(i,j,1,iEl,iVar,N,nel)] = fb[1];
    fBound[SCB_3D_INDEX(i,j,2,iEl,iVar,N,nel)] = fb[2];
    fBound[SCB_3D_INDEX(i,j,3,iEl,iVar,N,nel)] = fb[3];
    fBound[SCB_3D_INDEX(i,j,4,iEl,iVar,N,nel)] = fb[4];
    fBound[SCB_3D_INDEX(i,j,5,iEl,iVar,N,nel)] = fb[5];
  }
}

extern "C"
{
  void BoundaryInterp_3D_gpu(real *bMatrix, real *f, real *fBound, int N, int nvar, int nel)
  {
    int ndof = (N+1)*(N+1)*nel*nvar;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
	  BoundaryInterp_3D_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(bMatrix, f, fBound, N, nel, nvar);
  } 
}

template <int blockSize>
__global__ void __launch_bounds__(256) Divergence_2D_gpukernel(real *f, real *df, real *dmatrix, int nq, int N){

    uint32_t iq = threadIdx.x;
    if( iq < nq ){
        
        uint32_t iel = blockIdx.x;
        uint32_t nel = gridDim.x;
        uint32_t ivar = blockIdx.y;
        uint32_t nvar = gridDim.y;
        uint32_t i = iq % (N+1);
        uint32_t j = (iq/(N+1));

        __shared__ real f1[blockSize];
        __shared__ real f2[blockSize];
        __shared__ real dmloc[blockSize];
        f1[iq] = f[iq + nq*(iel + nel*(ivar))]; // x-component
        f2[iq] = f[iq + nq*(iel + nel*(ivar + nvar))]; // y-component
        dmloc[iq] = dmatrix[iq];
        __syncthreads();

        real dfloc = 0.0;
        for(int ii = 0; ii<N+1; ii++){
            dfloc += dmloc[ii+(N+1)*i]*f1[ii+(N+1)*(j)]+
                     dmloc[ii+(N+1)*j]*f2[i+(N+1)*(ii)];
        }
        df[iq + nq*(iel + nel*ivar)] = dfloc;
    }

}

extern "C"
{
  void Divergence_2D_gpu(real *f, real *df, real *dmatrix, int N, int nvar, int nel){
    int nq = (N+1)*(N+1);

    if( N <= 7 ){
      Divergence_2D_gpukernel<64><<<dim3(nel,nvar,1), dim3(64,1,1), 0, 0>>>(f,df,dmatrix,nq,N);
    } else {
      Divergence_2D_gpukernel<256><<<dim3(nel,nvar,1), dim3(256,1,1), 0, 0>>>(f,df,dmatrix,nq,N);
    }

  }
}

__global__ void __launch_bounds__(256) DG_BoundaryContribution_2D_gpukernel(real *bMatrix, real *qWeights, real *bf, real *df, int N, int nq){

  uint32_t iq = threadIdx.x;

  if( iq < nq ){
    uint32_t i = iq % (N+1);
    uint32_t j = (iq/(N+1));
    uint32_t iel = blockIdx.x;
    uint32_t nel = gridDim.x;
    uint32_t ivar = blockIdx.y;

    df[iq + nq*(iel + nel*ivar)] += (bMatrix[i+(N+1)]*bf[SCB_2D_INDEX(j,1,iel,ivar,N,nel)] + // east
                                           bMatrix[i]*bf[SCB_2D_INDEX(j,3,iel,ivar,N,nel)])/ // west
                                         qWeights[i];

    df[iq + nq*(iel + nel*ivar)] += (bMatrix[j+(N+1)]*bf[SCB_2D_INDEX(i,2,iel,ivar,N,nel)] + // north
                                           bMatrix[j]*bf[SCB_2D_INDEX(i,0,iel,ivar,N,nel)])/  // south
                                         qWeights[j];
  }

}

extern "C"
{
  void DG_BoundaryContribution_2D_gpu(real *bMatrix, real *qWeights, real *bf, real *df, int N, int nvar, int nel)
  {
    int nq = (N+1)*(N+1);

    if( N <= 7 ){
      DG_BoundaryContribution_2D_gpukernel<<<dim3(nel,nvar,1), dim3(64,1,1), 0, 0>>>(bMatrix, qWeights, bf, df, N, nq);
    } else {
      DG_BoundaryContribution_2D_gpukernel<<<dim3(nel,nvar,1), dim3(256,1,1), 0, 0>>>(bMatrix, qWeights, bf, df, N, nq);
    }
  } 
}

template<int blockSize, int matSize>
__global__ void __launch_bounds__(512) Divergence_3D_gpukernel(double *f, double *df, double *dmatrix, int nq, int N, int nel, int nvar){

    uint32_t idof = threadIdx.x;
    if( idof < nq ){
        
        uint32_t iel = blockIdx.x;
        uint32_t ivar = blockIdx.y;
        uint32_t i = idof % (N+1);
        uint32_t j = (idof/(N+1)) % (N+1);
        uint32_t k = (idof/(N+1)/(N+1));

        __shared__ double f1[blockSize];
        __shared__ double f2[blockSize];
        __shared__ double f3[blockSize];
        __shared__ double dmloc[matSize];
        f1[i+(N+1)*(j+(N+1)*k)] = f[i+(N+1)*(j+(N+1)*(k + (N+1)*(iel + nel*(ivar))))];
        f2[i+(N+1)*(j+(N+1)*k)] = f[i+(N+1)*(j+(N+1)*(k + (N+1)*(iel + nel*(ivar + nvar))))];
        f3[i+(N+1)*(j+(N+1)*k)] = f[i+(N+1)*(j+(N+1)*(k + (N+1)*(iel + nel*(ivar + 2*nvar))))];
        if( k == 0 ){
            dmloc[i+(N+1)*j] = dmatrix[i+(N+1)*j];
        }
        __syncthreads();

        double dfloc = 0.0;

        for(int ii = 0; ii<N+1; ii++){
            dfloc += dmloc[ii+(N+1)*i]*f1[ii+(N+1)*(j+(N+1)*(k))]+
                     dmloc[ii+(N+1)*j]*f2[i+(N+1)*(ii+(N+1)*(k))]+
                     dmloc[ii+(N+1)*k]*f3[i+(N+1)*(j+(N+1)*(ii))];
        }
        df[idof + nq*(iel + nel*ivar)] = dfloc;
    }

}

extern "C"
{
  void Divergence_3D_gpu(double *f, double *df, double *dmatrix, int N, int nel, int nvar){
    int nq = (N+1)*(N+1)*(N+1);
    if( N < 4 ){
        Divergence_3D_gpukernel<64,16><<<dim3(nel,nvar,1), dim3(64,1,1), 0, 0>>>(f,df,dmatrix,nq,N,nel,nvar);

    } else if( N >= 4 && N < 8 ){
        Divergence_3D_gpukernel<512,64><<<dim3(nel,nvar,1), dim3(512,1,1), 0, 0>>>(f,df,dmatrix,nq,N,nel,nvar);
    }
  }
}

__global__ void __launch_bounds__(512) DG_BoundaryContribution_3D_gpukernel(real *bMatrix, real *qWeights, real *bf, real *df, int N, int nel){

  uint32_t iq = threadIdx.x;
  uint32_t nq = (N+1)*(N+1)*(N+1);

  if( iq < nq ){
    uint32_t i = iq % (N+1);
    uint32_t j = (iq/(N+1))%(N+1);
    uint32_t k = (iq/(N+1)/(N+1))%(N+1);
    uint32_t iel = blockIdx.x;
    uint32_t ivar = blockIdx.y;
    df[iq + nq*(iel + nel*ivar)] += (bf[SCB_3D_INDEX(i,j,5,iel,ivar,N,nel)]*bMatrix[k+(N+1)] + // top
                                              bf[SCB_3D_INDEX(i,j,0,iel,ivar,N,nel)]*bMatrix[k])/       // bottom
                                              qWeights[k];

    df[iq + nq*(iel + nel*ivar)] += (bf[SCB_3D_INDEX(j,k,2,iel,ivar,N,nel)]*bMatrix[i+(N+1)] + // east
                                              bf[SCB_3D_INDEX(j,k,4,iel,ivar,N,nel)]*bMatrix[i])/       // west
                                              qWeights[i];

    df[iq + nq*(iel + nel*ivar)] += (bf[SCB_3D_INDEX(i,k,3,iel,ivar,N,nel)]*bMatrix[j+(N+1)] + // north
                                              bf[SCB_3D_INDEX(i,k,1,iel,ivar,N,nel)]*bMatrix[j])/       // south
                                              qWeights[j];
  }

}

extern "C"
{
  void DG_BoundaryContribution_3D_gpu(real *bMatrix, real *qWeights, real *bf, real *df, int N, int nvar, int nel)
  {

    int nq = (N+1)*(N+1)*(N+1);
    if( N < 4 ){
        DG_BoundaryContribution_3D_gpukernel<<<dim3(nel,nvar,1), dim3(64,1,1), 0, 0>>>(bMatrix, qWeights, bf, df, N, nel);

    } else if( N >= 4 && N < 8 ){
        DG_BoundaryContribution_3D_gpukernel<<<dim3(nel,nvar,1), dim3(512,1,1), 0, 0>>>(bMatrix, qWeights, bf, df, N, nel);
    }
  } 
}