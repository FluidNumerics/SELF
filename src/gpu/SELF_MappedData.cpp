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
  uint32_t ivar = blockIdx.y;
  
  if(idof < ndof){
    int e2Global = sideInfo[INDEX3(2,s1,e1,5,4)];
    int e2 = e2Global - offset;
    int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
    int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
    
    if(e2Global != 0){
      int neighborRank = elemToRank[e2Global-1];
      if( neighborRank == rankId ){
        if(flip == 0){
          extBoundary[idof + ndof*ivar] = boundary[SCB_2D_INDEX(i1,s2-1,e2-1,ivar,N,nEl)];
        }
        else if(flip == 1){
          int i2 = N-i1;
          extBoundary[idof + ndof*ivar] = boundary[SCB_2D_INDEX(i2,s2-1,e2-1,ivar,N,nEl)];
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

__global__ void ContravariantWeight_2D_gpukernel(real *scalar, real *dsdx, real *tensor, int ndof){

  uint32_t ivar = blockIdx.y; // variable dimension
  uint32_t nvar = gridDim.y; // number of variables
  uint32_t tdim = blockIdx.z; // tensor dimension (flattened index for the rows and columns of the tensor)
  uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;

  if( i < ndof ){
    tensor[i+ndof*(ivar + nvar*tdim)] = dsdx[i+ndof*tdim]*scalar[i+ndof*ivar];
  }

}

extern "C"
{
  void ContravariantWeight_2D_gpu(real *scalar, real *dsdx, real *tensor, int N, int nvar, int nel)
  {
    int ndof = (N+1)*(N+1)*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nvar,4);
    dim3 nthreads(threads_per_block,1,1);

    ContravariantWeight_2D_gpukernel<<<nblocks, nthreads, 0, 0>>>(scalar, dsdx, tensor, ndof);

  }
}

__global__ void NormalWeight_2D_gpukernel(real *fb, real *nhat, real *nscale, real *fbn, int N, int nvar, int nel){

  uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*4*nel;

  if( i < ndof ){
    uint32_t ivar = blockIdx.y;
    real f = fb[i+ndof*ivar];
    real nmag = nscale[i];
    fbn[i+ndof*ivar] = f*nhat[i]*nmag; // x-direction
    fbn[i+ndof*(ivar+nvar)] = f*nhat[i+ndof]*nmag; // y-direction
  }

}

extern "C"
{
  void NormalWeight_2D_gpu(real *fb, real *nhat, real *nscale, real *fbn, int N, int nvar, int nel)
  {
    int ndof = (N+1)*4*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nvar,1);
    dim3 nthreads(threads_per_block,1,1);

    NormalWeight_2D_gpukernel<<<nblocks, nthreads, 0, 0>>>(fb,nhat,nscale,fbn,N,nvar,nel);

  }
}
__global__ void DG_BoundaryContribution_2D_gpukernel(real *bMatrix, real *qWeights, real *bf, real *df, int N, int nel){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ivar = blockIdx.y;
  uint32_t ndof = nel*(N+1)*(N+1);
  if( idof < ndof ){
    uint32_t i = idof % (N+1);
    uint32_t j = (idof/(N+1))%(N+1);
    uint32_t iel = (idof/(N+1)/(N+1));
    df[SC_2D_INDEX(i,j,iel,ivar,N,nel)] += (bf[SCB_2D_INDEX(j,1,iel,ivar,N,nel)]*bMatrix[i+(N+1)] + // east
                                          bf[SCB_2D_INDEX(j,3,iel,ivar,N,nel)]*bMatrix[i])/       // west
                                         qWeights[i] +
                                         (bf[SCB_2D_INDEX(i,2,iel,ivar,N,nel)]*bMatrix[j+(N+1)] + // north
                                          bf[SCB_2D_INDEX(i,0,iel,ivar,N,nel)]*bMatrix[j])/       // south
                                         qWeights[j];
  }

}

extern "C"
{
  void DG_BoundaryContribution_2D_gpu(real *bMatrix, real *qWeights, real *bf, real *df, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(N+1)*nEl;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nVar,1);
    dim3 nthreads(threads_per_block,1,1);

	  DG_BoundaryContribution_2D_gpukernel<<<nblocks,nthreads, 0, 0>>>(bMatrix, qWeights, bf, df, N, nEl);
  } 
}

// __global__ void VectorDG_BoundaryContribution_2D_gpukernel(real *bMatrix, real *qWeights, real *bf, real *df, int N, int nel){

//   uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
//   uint32_t ivar = blockIdx.y;
//   uint32_t ndof = nel*(N+1)*(N+1);
//   if( idof < ndof ){
//     uint32_t i = idof % (N+1);
//     uint32_t j = (idof/(N+1))%(N+1);
//     uint32_t iel = (idof/(N+1)/(N+1));
//     df[SC_2D_INDEX(i,j,iel,ivar,N,nel)] += (bf[SCB_2D_INDEX(j,1,iel,ivar,N,nel)]*bMatrix[i+(N+1)] + // east
//                                           bf[SCB_2D_INDEX(j,3,iel,ivar,N,nel)]*bMatrix[i])/       // west
//                                          qWeights[i] +
//                                          (bf[SCB_2D_INDEX(i,2,iel,ivar,N,nel)]*bMatrix[j+(N+1)] + // north
//                                           bf[SCB_2D_INDEX(i,0,iel,ivar,N,nel)]*bMatrix[j])/       // south
//                                          qWeights[j];
//   }

// }

// extern "C"
// {
//   void VectorDG_BoundaryContribution_2D_gpu(real *bMatrix, real *qWeights, real *bf, real *df, int N, int nVar, int nEl)
//   {
//     int ndof = (N+1)*(N+1)*nEl;
//     int threads_per_block = 256;
//     int nblocks_x = ndof/threads_per_block + 1;

//     dim3 nblocks(nblocks_x,nVar,1);
//     dim3 nthreads(threads_per_block,1,1);

// 	  VectorDG_BoundaryContribution_2D_gpukernel<<<nblocks,nthreads, 0, 0>>>(bMatrix, qWeights, bf, df, N, nEl);
//   } 
// }

__global__ void ContravariantProjection_2D_gpukernel(real *vector, real *dsdx, int N, int ndof){

    uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;
    uint32_t ivar = blockIdx.y;
    uint32_t nvar = blockDim.y;

    real Fx = vector[i+ndof*ivar];
    real Fy = vector[i+ndof*(ivar + nvar)];

    if( i < ndof ){
      vector[i+ndof*ivar] = dsdx[i]*Fx + dsdx[i + ndof]*Fy;
      vector[i+ndof*(ivar + nvar)] = dsdx[i + 2*ndof]*Fx + dsdx[i + 3*ndof]*Fy;
    }

}

extern "C"
{
  void ContravariantProjection_2D_gpu(real *vector, real *dsdx, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(N+1)*nEl;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nVar,1);
    dim3 nthreads(threads_per_block,1,1);
    ContravariantProjection_2D_gpukernel<<<nblocks,nthreads, 0, 0>>>(vector, dsdx, N, ndof);

  } 
}


__global__ void ContravariantProjection_3D_gpukernel(real *vector, real *dsdx, int N, int ndof){

    uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;
    uint32_t ivar = blockIdx.y;
    uint32_t nvar = blockDim.y;

    real Fx = vector[i+ndof*ivar];
    real Fy = vector[i+ndof*(ivar + nvar)];
    real Fz = vector[i+ndof*(ivar + 2*nvar)];

    if( i < ndof ){
      vector[i+ndof*ivar] = dsdx[i]*Fx + dsdx[i + ndof]*Fy + dsdx[i + 2*ndof]*Fz;
      vector[i+ndof*(ivar + nvar)] = dsdx[i + 3*ndof]*Fx + dsdx[i + 4*ndof]*Fy + + dsdx[i + 5*ndof]*Fz;
      vector[i+ndof*(ivar + 2*nvar)] = dsdx[i + 6*ndof]*Fx + dsdx[i + 7*ndof]*Fy + + dsdx[i + 8*ndof]*Fz;
    }

}

extern "C"
{
  void ContravariantProjection_3D_gpu(real *vector, real *dsdx, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(N+1)*(N+1)*nEl;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nVar,1);
    dim3 nthreads(threads_per_block,1,1);
    ContravariantProjection_3D_gpukernel<<<nblocks,nthreads, 0, 0>>>(vector, dsdx, N, ndof);
  } 
}