#include "SELF_GPU_Macros.h"

// JacobianWeight functions
// The functions take in an array of data and divide by the jacobian.
__global__ void JacobianWeight(real *f, real *jacobian, int ndof){

  size_t ivar = blockIdx.y;
  size_t idof = threadIdx.x + blockIdx.x*blockDim.x;

  if( idof < ndof ){
    f[idof + ndof*ivar] = f[idof + ndof*ivar]/jacobian[idof];
  }
}

extern "C"
{
  void JacobianWeight_1D_gpu(real *f, real *jacobian, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(nEl);
    int threads_per_block = 256;
    int nblocksx = ndof/threads_per_block+1; 
    JacobianWeight<<<dim3(nblocksx,nVar,1), dim3(threads_per_block,1,1), 0, 0>>>(f, jacobian, ndof);
  }
}
extern "C"
{
  void JacobianWeight_2D_gpu(real *f, real *jacobian, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(N+1)*(nEl);
    int threads_per_block = 256;
    int nblocksx = ndof/threads_per_block+1; 
    JacobianWeight<<<dim3(nblocksx,nVar,1), dim3(threads_per_block,1,1), 0, 0>>>(f, jacobian, ndof);
  }
}
extern "C"
{
  void JacobianWeight_3D_gpu(real *f, real *jacobian, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(N+1)*(N+1)*(nEl);
    int threads_per_block = 256;
    int nblocksx = ndof/threads_per_block+1; 
    JacobianWeight<<<dim3(nblocksx,nVar,1), dim3(threads_per_block,1,1), 0, 0>>>(f, jacobian, ndof);
  }
}

__global__ void DGDerivative_BoundaryContribution_1D(real *bMatrix, real *qWeight, real *bf, real *df, int N, int nEl){

  
  uint32_t ivar = blockIdx.y;
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;

  uint32_t ndof=(N+1)*nEl;
  if( idof < ndof ){
    uint32_t iel = idof/(N+1); // Calculate the element ID
    uint32_t i = idof - iel*(N+1); // Calculate the quadrature point node id
    df[idof + ndof*ivar] += (bMatrix[i+(N+1)]*bf[SCB_1D_INDEX(1,iel,ivar,nEl)]+
	                        bMatrix[i]*bf[SCB_1D_INDEX(0,iel,ivar,nEl)])/qWeight[i];
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
  }
}


__global__ void SideExchange_2D(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nEl){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*nEl*4;
  uint32_t i1 = idof % (N+1);
  uint32_t s1 = (idof/(N+1)) % 4;
  uint32_t e1 = idof/(N+1)/4;
  size_t ivar = blockIdx.y;
  
  if(idof < ndof){
    int e2Global = sideInfo[INDEX3(2,s1,e1,5,4)]-1;
    int e2 = e2Global - offset;
    int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
    int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
    int bcid = sideInfo[INDEX3(4,s1,e1,5,4)];
    
    if(s2 > 0 || bcid == 0){
      int neighborRank = elemToRank[e2Global];
      if( neighborRank == rankId ){
        if(flip == 0){
          extBoundary[idof + ndof*ivar] = boundary[SCB_2D_INDEX(i1,s2-1,e2,ivar,N,nEl)];
        }
        else if(flip == 1){
          int i2 = N-i1;
          extBoundary[idof + ndof*ivar] = boundary[SCB_2D_INDEX(i2,s2-1,e2,ivar,N,nEl)];
        }
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_2D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*4*nEl;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nVar,1);
    dim3 nthreads(threads_per_block,1,1);
    SideExchange_2D<<<nblocks,nthreads>>>(extBoundary, boundary, sideInfo, elemToRank, rankId, offset, N, nEl);
  }
}
