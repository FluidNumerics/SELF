#include "SELF_GPU_Macros.h"

// JacobianWeight functions
// The functions take in an array of data and divide by the jacobian.
__global__ void JacobianWeight(real *f, real *jacobian, int ndof){

 uint32_t ivar = blockIdx.y;
 uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;

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

__global__ void ApplyFlip_2D(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nVar, int nEl){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = nVar*nEl*4;
  uint32_t s1 = (idof) % 4;
  uint32_t e1 = (idof/4) % nEl;
  uint32_t ivar = idof/4/nEl;
  
  if(idof < ndof){
    int e2Global = sideInfo[INDEX3(2,s1,e1,5,4)];
    int e2 = e2Global - offset;
    int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
    int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
    real buff[24]; // warning : set fixed buffer size for applying flip. This limits the polynomial degree to 23 

    if(e2Global != 0){
      int neighborRank = elemToRank[e2Global-1];
      if( neighborRank != rankId ){
        if(flip == 1){
          for( int i1 = 0; i1<N+1; i1++){
            int i2 = N-i1;
            buff[i1] = extBoundary[SCB_2D_INDEX(i2,s1,e1,ivar,N,nEl)];
          }
          for( int i1 = 0; i1<N+1; i1++){
            extBoundary[SCB_2D_INDEX(i1,s1,e1,ivar,N,nEl)] = buff[i1];
          }
        }
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_2D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    int ndof = 4*nEl*nVar;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,1,1);
    dim3 nthreads(threads_per_block,1,1);
    ApplyFlip_2D<<<nblocks,nthreads>>>(extBoundary, sideInfo, elemToRank, rankId, offset, N, nVar, nEl);
  }
}

__global__ void SideExchange_3D(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nEl){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*(N+1)*nEl*6;

  if(idof < ndof){

    uint32_t s1 = (idof/(N+1)/(N+1)) % 6;
    uint32_t e1 = idof/(N+1)/(N+1)/6;

    int e2Global = sideInfo[INDEX3(2,s1,e1,5,6)];
    
    if(e2Global != 0){

      int neighborRank = elemToRank[e2Global-1];

      if( neighborRank == rankId ){
        uint32_t i1 = idof % (N+1);
        uint32_t j1 = (idof/(N+1)) % (N+1);
        int e2 = e2Global - offset;
        int s2 = sideInfo[INDEX3(3,s1,e1,5,6)]/10;
        int flip = sideInfo[INDEX3(3,s1,e1,5,6)]-s2*10; 
        uint32_t ivar = blockIdx.y;
        if(flip == 0){
          extBoundary[idof + ndof*ivar] = boundary[SCB_3D_INDEX(i1,j1,s2-1,e2-1,ivar,N,nEl)];
        }
        else if(flip == 1){
          int i2 = N-i1;
          int j2 = j1;
          extBoundary[idof + ndof*ivar] = boundary[SCB_3D_INDEX(i2,j2,s2-1,e2-1,ivar,N,nEl)];
        }
        else if(flip == 2){
          int i2 = N-i1;
          int j2 = N-j1;
          extBoundary[idof + ndof*ivar] = boundary[SCB_3D_INDEX(i2,j2,s2-1,e2-1,ivar,N,nEl)];
        }
        else if(flip == 3){
          int i2 = i1;
          int j2 = N-j1;
          extBoundary[idof + ndof*ivar] = boundary[SCB_3D_INDEX(i2,j2,s2-1,e2-1,ivar,N,nEl)];
        }
        else if(flip == 4){
          int i2 = j1;
          int j2 = i1;
          extBoundary[idof + ndof*ivar] = boundary[SCB_3D_INDEX(i2,j2,s2-1,e2-1,ivar,N,nEl)];
        }
        else if(flip == 5){
          int i2 = N-j1;
          int j2 = i1;
          extBoundary[idof + ndof*ivar] = boundary[SCB_3D_INDEX(i2,j2,s2-1,e2-1,ivar,N,nEl)];
        }
        else if(flip == 6){
          int i2 = N-j1;
          int j2 = N-i1;
          extBoundary[idof + ndof*ivar] = boundary[SCB_3D_INDEX(i2,j2,s2-1,e2-1,ivar,N,nEl)];
        }
        else if(flip == 7){
          int i2 = j1;
          int j2 = N-i1;
          extBoundary[idof + ndof*ivar] = boundary[SCB_3D_INDEX(i2,j2,s2-1,e2-1,ivar,N,nEl)];
        }
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_3D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(N+1)*6*nEl;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nVar,1);
    dim3 nthreads(threads_per_block,1,1);
    SideExchange_3D<<<nblocks,nthreads>>>(extBoundary, boundary, sideInfo, elemToRank, rankId, offset, N, nEl);
  }
}

__global__ void ApplyFlip_3D(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nVar, int nEl){

  uint32_t s1 = blockIdx.x;
  uint32_t e1 = blockIdx.y;
  uint32_t ivar = blockIdx.z;
  uint32_t i = threadIdx.x;
  uint32_t j = threadIdx.y;

  __shared__ real extBuff[256];

  extBuff[i+(N+1)*j] = extBoundary[SCB_3D_INDEX(i,j,s1,e1,ivar,N,nEl)];

  __syncthreads();
  
  int e2Global = sideInfo[INDEX3(2,s1,e1,5,4)];
  int e2 = e2Global - offset;
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;

  if(e2Global != 0){
    int neighborRank = elemToRank[e2Global-1];
    if( neighborRank != rankId ){

      if(flip == 1){
        int i2 = N-i;
        int j2 = j;
        extBoundary[SCB_3D_INDEX(i,j,s1,e1,ivar,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 2){
        int i2 = N-i;
        int j2 = N-j;
        extBoundary[SCB_3D_INDEX(i,j,s1,e1,ivar,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 3){
        int i2 = i;
        int j2 = N-j;
        extBoundary[SCB_3D_INDEX(i,j,s1,e1,ivar,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 4){
        int i2 = j;
        int j2 = i;
        extBoundary[SCB_3D_INDEX(i,j,s1,e1,ivar,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 5){
        int i2 = N-j;
        int j2 = i;
        extBoundary[SCB_3D_INDEX(i,j,s1,e1,ivar,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 6){
        int i2 = N-j;
        int j2 = N-i;
        extBoundary[SCB_3D_INDEX(i,j,s1,e1,ivar,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 7){
        int i2 = j;
        int j2 = N-i;
        extBoundary[SCB_3D_INDEX(i,j,s1,e1,ivar,N,nEl)] = extBuff[i2+(N+1)*j2];    
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_3D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    dim3 nblocks(6,nEl,nVar);
    dim3 nthreads(N+1,N+1,1);
    ApplyFlip_3D<<<nblocks,nthreads>>>(extBoundary, sideInfo, elemToRank, rankId, offset, N, nVar, nEl);
  }
}

__global__ void ContravariantWeight_gpukernel(real *scalar, real *dsdx, real *tensor, int ndof){

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

    ContravariantWeight_gpukernel<<<nblocks, nthreads, 0, 0>>>(scalar, dsdx, tensor, ndof);

  }
}

extern "C"
{
  void ContravariantWeight_3D_gpu(real *scalar, real *dsdx, real *tensor, int N, int nvar, int nel)
  {
    int ndof = (N+1)*(N+1)*(N+1)*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nvar,9);
    dim3 nthreads(threads_per_block,1,1);

    ContravariantWeight_gpukernel<<<nblocks, nthreads, 0, 0>>>(scalar, dsdx, tensor, ndof);

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

__global__ void NormalWeight_3D_gpukernel(real *fb, real *nhat, real *nscale, real *fbn, int N, int nvar, int nel){

  uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*(N+1)*6*nel;

  if( i < ndof ){
    uint32_t ivar = blockIdx.y;
    real f = fb[i+ndof*ivar];
    real nmag = nscale[i];
    fbn[i+ndof*ivar] = f*nhat[i]*nmag; // x-direction
    fbn[i+ndof*(ivar+nvar)] = f*nhat[i+ndof]*nmag; // y-direction
    fbn[i+ndof*(ivar+2*nvar)] = f*nhat[i+2*ndof]*nmag; // z-direction
  }

}

extern "C"
{
  void NormalWeight_3D_gpu(real *fb, real *nhat, real *nscale, real *fbn, int N, int nvar, int nel)
  {
    int ndof = (N+1)*(N+1)*6*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nvar,1);
    dim3 nthreads(threads_per_block,1,1);

    NormalWeight_3D_gpukernel<<<nblocks, nthreads, 0, 0>>>(fb,nhat,nscale,fbn,N,nvar,nel);
  }
}

__global__ void ContravariantProjection_2D_gpukernel(real *vector, real *dsdx, int nq){

    uint32_t idof = threadIdx.x;

    if( idof < nq ){
      uint32_t iel = blockIdx.x;
      uint32_t nel = gridDim.x;
      uint32_t ivar = blockIdx.y;
      uint32_t nvar = gridDim.y;
      real Fx = vector[idof + nq*(iel + nel*(ivar))];
      real Fy = vector[idof + nq*(iel + nel*(ivar + nvar))];
            
      vector[idof + nq*(iel + nel*(ivar))] = dsdx[idof + nq*(iel)]*Fx + // dsdx(...,0,0)*Fx
                                             dsdx[idof + nq*(iel+nel)]*Fy; // dsdx(...,1,0)*Fy;

      vector[idof + nq*(iel + nel*(ivar+nvar))] = dsdx[idof + nq*(iel+nel*2)]*Fx + //dsdx(...,0,1)*Fx
                                                  dsdx[idof + nq*(iel+nel*3)]*Fy;  //dsdx(...,1,1)*Fy
    }

}

extern "C"
{
  void ContravariantProjection_2D_gpu(real *vector, real *dsdx, int N, int nVar, int nEl)
  {
    int nq = (N+1)*(N+1);

    if( N <= 7 ){
      ContravariantProjection_2D_gpukernel<<<dim3(nEl,nVar,1), dim3(64,1,1), 0, 0>>>(vector, dsdx, nq);
    } else {
      ContravariantProjection_2D_gpukernel<<<dim3(nEl,nVar,1), dim3(256,1,1), 0, 0>>>(vector, dsdx, nq);
    }
  } 
}

__global__ void ContravariantProjection_3D_gpukernel(real *vector, real *dsdx, int N, int ndof){

    uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;

    if( i < ndof ){
      uint32_t ivar = blockIdx.y;
      uint32_t nvar = blockDim.y;
      real Fx = vector[i+ndof*ivar];
      real Fy = vector[i+ndof*(ivar + nvar)];
      real Fz = vector[i+ndof*(ivar + 2*nvar)];
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