#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

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
  void JacobianWeight_1D_gpu(real *f, real *dxds, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(nEl);
    int threads_per_block = 256;
    int nblocksx = ndof/threads_per_block+1; 
    JacobianWeight<<<dim3(nblocksx,nVar,1), dim3(threads_per_block,1,1), 0, 0>>>(f, dxds, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

extern "C"
{
  void JacobianWeight_2D_gpu(real *scalar, real *jacobian, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(N+1)*(nEl);
    int threads_per_block = 256;
    int nblocksx = ndof/threads_per_block+1; 

    dim3 nblocks(nblocksx,nVar,1);
    dim3 nthreads(threads_per_block,1,1);

    JacobianWeight<<<nblocks, nthreads, 0, 0>>>(scalar, jacobian, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

extern "C"
{
  void JacobianWeight_3D_gpu(real *scalar, real *jacobian, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(N+1)*(N+1)*(nEl);
    int threads_per_block = 256;
    int nblocksx = ndof/threads_per_block+1; 

    dim3 nblocks(nblocksx,nVar,1);
    dim3 nthreads(threads_per_block,1,1);

    JacobianWeight<<<nblocks, nthreads, 0, 0>>>(scalar, jacobian, ndof);
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void SideExchange_3D(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nEl){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*(N+1)*nEl*6;
  uint32_t i1 = idof % (N+1);
  uint32_t j1 = (idof/(N+1)) % (N+1);
  uint32_t s1 = (idof/(N+1)/(N+1)) % 6;
  uint32_t e1 = idof/(N+1)/(N+1)/6;
  size_t ivar = blockIdx.y;
  

  if(idof < ndof){
    int e2Global = sideInfo[INDEX3(2,s1,e1,5,6)]-1;
    int e2 = e2Global - offset;
    int s2 = sideInfo[INDEX3(3,s1,e1,5,6)]/10;
    int flip = sideInfo[INDEX3(3,s1,e1,5,6)]-s2*10;
    int bcid = sideInfo[INDEX3(4,s1,e1,5,6)];

    if(s2 > 0 || bcid == 0){
      int neighborRank = elemToRank[e2Global];
      if( neighborRank == rankId ){
        if(flip == 0){
          extBoundary[SCB_3D_INDEX(i1,j1,s1,e1,ivar,N,nEl)] = boundary[SCB_3D_INDEX(i1,j1,s2-1,e2,ivar,N,nEl)];
        }
        else if(flip == 1){
          int i2 = j1;
          int j2 = N-i1;
          extBoundary[SCB_3D_INDEX(i1,j1,s1,e1,ivar,N,nEl)] = boundary[SCB_3D_INDEX(i2,j2,s2-1,e2,ivar,N,nEl)];
        }
        else if(flip == 2){
          int i2 = N-i1;
          int j2 = N-j1;
          extBoundary[SCB_3D_INDEX(i1,j1,s1,e1,ivar,N,nEl)] = boundary[SCB_3D_INDEX(i2,j2,s2-1,e2,ivar,N,nEl)];
        }
        else if(flip == 3){
          int i2 = N-j1;
          int j2 = i1;
          extBoundary[SCB_3D_INDEX(i1,j1,s1,e1,ivar,N,nEl)] = boundary[SCB_3D_INDEX(i2,j2,s2-1,e2,ivar,N,nEl)];
        }
        else if(flip == 4){
          int i2 = j1;
          int j2 = i1;
          extBoundary[SCB_3D_INDEX(i1,j1,s1,e1,ivar,N,nEl)] = boundary[SCB_3D_INDEX(i2,j2,s2-1,e2,ivar,N,nEl)];
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
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void ApplyFlip_2D(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nEl){

  size_t s1 = blockIdx.x;
  size_t e1 = blockIdx.y;
  size_t ivar = blockIdx.z;
  size_t i1 = threadIdx.x;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
  int i2 = N-i1;
  int neighborRank = elemToRank[e2];
  int bcid = sideInfo[INDEX3(4,s1,e1,5,4)];


  if(s2 > 0 || bcid == 0){ // Interior Element
    if( neighborRank != rankId ){ // Side shared with another rank
      if(flip == 1){ // Neighboring elements have different orientation

        __shared__ real extBuff[16];
        extBuff[i1] = extBoundary[SCB_2D_INDEX(i2,s1,e1,ivar,N,nEl)];
        __syncthreads();

        extBoundary[SCB_2D_INDEX(i1,s1,e1,ivar,N,nEl)] = extBuff[i1];
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_2D_gpu(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_2D<<<dim3(4,nEl,nVar), dim3(N+1,1,1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void ApplyFlip_3D(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nEl){


  size_t s1 = blockIdx.x;
  size_t e1 = blockIdx.y;
  size_t ivar = blockIdx.z;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,6)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,6)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,6)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1,e1,5,6)];

  __shared__ real extBuff[256];

  extBuff[i1+(N+1)*j1] = extBoundary[SCB_3D_INDEX(i1,j1,s1,e1,ivar,N,nEl)];

  __syncthreads();

  if(s2 > 0 || bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank != rankId ){
      if(flip == 1){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[SCB_3D_INDEX(i1,j1,s1,e1,ivar,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 2){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[SCB_3D_INDEX(i1,j1,s1,e1,ivar,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 3){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[SCB_3D_INDEX(i1,j1,s1,e1,ivar,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = i1;
        extBoundary[SCB_3D_INDEX(i1,j1,s1,e1,ivar,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
    }
  }
}

extern "C"
{
  void ApplyFlip_3D_gpu(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_3D<<<dim3(6,nEl,nVar), dim3(N+1,N+1,1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void BassiRebaySides_gpu(real *avgBoundary, real *boundary, real *extBoundary, int ndof){

  uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;

  if( i < ndof ){
    avgBoundary[i] =0.5*(boundary[i]+extBoundary[i]);
  }
}

extern "C"
{
  void BassiRebaySides_2D_gpu(real *avgBoundary, real *boundary, real *extBoundary, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*4*nEl*nVar;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,1,1);
    dim3 nthreads(threads_per_block,1,1);

    BassiRebaySides_gpu<<<nblocks, nthreads, 0, 0>>>(avgBoundary, boundary, extBoundary, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

extern "C"
{
  void BassiRebaySides_3D_gpu(real *avgBoundary, real *boundary, real *extBoundary, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(N+1)*6*nEl*nVar;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,1,1);
    dim3 nthreads(threads_per_block,1,1);

    BassiRebaySides_gpu<<<nblocks, nthreads, 0, 0>>>(avgBoundary, boundary, extBoundary, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void ContravariantWeight_MappedScalar_gpu(real *scalar, real *dsdx, real *tensor, int ndof){

  uint32_t ivar = blockIdx.y; // variable dimension
  uint32_t nvar = gridDim.y; // number of variables
  uint32_t tdim = blockIdx.z; // tensor dimension (flattened index for the rows and columns of the tensor)
  uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;
  //uint32_t tid = threadIdx.x + blockIdx.x*blockDim.x;
  //uint32_t stride = blockDim.x*gridDim.x;

  //for (uint32_t i = tid; i < ndof; i += stride){
  if( i < ndof ){
    tensor[i+ndof*(ivar + nvar*tdim)] = dsdx[i+ndof*tdim]*scalar[i+ndof*ivar];
  }

}

extern "C"
{
  void ContravariantWeight_MappedScalar2D_gpu(real *scalar, real *dsdx, real *tensor, int isize, int jsize, int nvar, int nel)
  {
    int ndof = isize*jsize*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nvar,4);
    dim3 nthreads(threads_per_block,1,1);

    ContravariantWeight_MappedScalar_gpu<<<nblocks, nthreads, 0, 0>>>(scalar, dsdx, tensor, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

extern "C"
{
  void ContravariantWeight_MappedScalar3D_gpu(real *scalar, real *dsdx, real *tensor, int isize, int jsize, int ksize, int nvar, int nel)
  {
    int ndof = isize*jsize*ksize*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nvar,9);
    dim3 nthreads(threads_per_block,1,1);

    ContravariantWeight_MappedScalar_gpu<<<nblocks, nthreads, 0, 0>>>(scalar, dsdx, tensor, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void ContravariantProjection_MappedVector2D_gpu(real *vector, real *dsdx, int N, int ndof){

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
  void ContravariantProjection_MappedVector2D_gpu_wrapper(real *vector, real *dsdx, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(N+1)*nEl;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nVar,1);
    dim3 nthreads(threads_per_block,1,1);
    ContravariantProjection_MappedVector2D_gpu<<<nblocks,nthreads, 0, 0>>>(vector, dsdx, N, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}


__global__ void ContravariantProjection_MappedVector3D_gpu(real *vector, real *dsdx, int N, int ndof){

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
  void ContravariantProjection_MappedVector3D_gpu_wrapper(real *vector, real *dsdx, int N, int nVar, int nEl)
  {
    int ndof = (N+1)*(N+1)*(N+1)*nEl;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nVar,1);
    dim3 nthreads(threads_per_block,1,1);
    ContravariantProjection_MappedVector3D_gpu<<<nblocks,nthreads, 0, 0>>>(vector, dsdx, N, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}