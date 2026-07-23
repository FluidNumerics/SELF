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

// Aggregated MPI halo pack/unpack for the 2-D and 3-D side exchanges.
//
// haloSides holds (element, side, flip) integer triplets, with 0-based element
// and side indices, for every locally-owned side whose neighbor element lives
// on another rank. Entries are grouped by neighbor rank and sorted by global
// side id within each group so that the send and receive buffers on the two
// ranks of an interface enumerate sides in the same order; no per-message
// metadata is required. Buffer layout for entry n, variable ivar:
//   buf[i + (N+1)*(j + (N+1)*(ivar + nVar*n))]
// The sender packs its boundary trace in its native orientation; the
// receiver applies its side's flip permutation during the unpack.

__global__ void HaloPack_2D(real *boundary, real *sendBuf, int *haloSides, int N, int nVar, int nEl, int nHalo){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*nHalo;

  if(idof < ndof){
    uint32_t i = idof % (N+1);
    uint32_t n = idof/(N+1);
    uint32_t ivar = blockIdx.y;
    int e1 = haloSides[3*n];
    int s1 = haloSides[3*n+1];

    sendBuf[i + (N+1)*(ivar + nVar*n)] = boundary[SCB_2D_INDEX(i,s1,e1,ivar,N,nEl)];
  }

}

extern "C"
{
  void HaloPack_2D_gpu(real *boundary, real *sendBuf, int *haloSides, int N, int nVar, int nEl, int nHalo)
  {
    int ndof = (N+1)*nHalo;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nVar,1);
    dim3 nthreads(threads_per_block,1,1);
    HaloPack_2D<<<nblocks,nthreads>>>(boundary, sendBuf, haloSides, N, nVar, nEl, nHalo);
    // The packed send buffer must be complete before MPI_Isend is posted on
    // the host; synchronize the device before returning.
#ifdef __HIP_PLATFORM_AMD__
    CHECK(hipDeviceSynchronize());
#else
    CHECK(cudaDeviceSynchronize());
#endif
  }
}

__global__ void HaloUnpack_2D(real *recvBuf, real *extBoundary, int *haloSides, int N, int nVar, int nEl, int nHalo){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*nHalo;

  if(idof < ndof){
    uint32_t i = idof % (N+1);
    uint32_t n = idof/(N+1);
    uint32_t ivar = blockIdx.y;
    int e1 = haloSides[3*n];
    int s1 = haloSides[3*n+1];
    int flip = haloSides[3*n+2];
    int i2 = i;

    // extBoundary(i) receives the neighbor's trace at (i2), reversing the
    // trace when the neighbor's coordinate runs in the opposite direction.
    if(flip == 1){
      i2 = N-i;
    }

    extBoundary[SCB_2D_INDEX(i,s1,e1,ivar,N,nEl)] = recvBuf[i2 + (N+1)*(ivar + nVar*n)];
  }

}

extern "C"
{
  void HaloUnpack_2D_gpu(real *recvBuf, real *extBoundary, int *haloSides, int N, int nVar, int nEl, int nHalo)
  {
    int ndof = (N+1)*nHalo;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nVar,1);
    dim3 nthreads(threads_per_block,1,1);
    HaloUnpack_2D<<<nblocks,nthreads>>>(recvBuf, extBoundary, haloSides, N, nVar, nEl, nHalo);
  }
}

__global__ void HaloPack_3D(real *boundary, real *sendBuf, int *haloSides, int N, int nVar, int nEl, int nHalo){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*(N+1)*nHalo;

  if(idof < ndof){
    uint32_t i = idof % (N+1);
    uint32_t j = (idof/(N+1)) % (N+1);
    uint32_t n = idof/(N+1)/(N+1);
    uint32_t ivar = blockIdx.y;
    int e1 = haloSides[3*n];
    int s1 = haloSides[3*n+1];

    sendBuf[i + (N+1)*(j + (N+1)*(ivar + nVar*n))] = boundary[SCB_3D_INDEX(i,j,s1,e1,ivar,N,nEl)];
  }

}

extern "C"
{
  void HaloPack_3D_gpu(real *boundary, real *sendBuf, int *haloSides, int N, int nVar, int nEl, int nHalo)
  {
    int ndof = (N+1)*(N+1)*nHalo;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nVar,1);
    dim3 nthreads(threads_per_block,1,1);
    HaloPack_3D<<<nblocks,nthreads>>>(boundary, sendBuf, haloSides, N, nVar, nEl, nHalo);
    // The packed send buffer must be complete before MPI_Isend is posted on
    // the host; synchronize the device before returning.
#ifdef __HIP_PLATFORM_AMD__
    CHECK(hipDeviceSynchronize());
#else
    CHECK(cudaDeviceSynchronize());
#endif
  }
}

__global__ void HaloUnpack_3D(real *recvBuf, real *extBoundary, int *haloSides, int N, int nVar, int nEl, int nHalo){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*(N+1)*nHalo;

  if(idof < ndof){
    uint32_t i = idof % (N+1);
    uint32_t j = (idof/(N+1)) % (N+1);
    uint32_t n = idof/(N+1)/(N+1);
    uint32_t ivar = blockIdx.y;
    int e1 = haloSides[3*n];
    int s1 = haloSides[3*n+1];
    int flip = haloSides[3*n+2];
    int i2 = i;
    int j2 = j;

    // extBoundary(i,j) receives the neighbor's trace at (i2,j2), matching the
    // relative orientation (flip) of the two element faces on the interface.
    if(flip == 1){
      i2 = N-i;
      j2 = j;
    }
    else if(flip == 2){
      i2 = N-i;
      j2 = N-j;
    }
    else if(flip == 3){
      i2 = i;
      j2 = N-j;
    }
    else if(flip == 4){
      i2 = j;
      j2 = i;
    }
    else if(flip == 5){
      i2 = N-j;
      j2 = i;
    }
    else if(flip == 6){
      i2 = N-j;
      j2 = N-i;
    }
    else if(flip == 7){
      i2 = j;
      j2 = N-i;
    }

    extBoundary[SCB_3D_INDEX(i,j,s1,e1,ivar,N,nEl)] = recvBuf[i2 + (N+1)*(j2 + (N+1)*(ivar + nVar*n))];
  }

}

extern "C"
{
  void HaloUnpack_3D_gpu(real *recvBuf, real *extBoundary, int *haloSides, int N, int nVar, int nEl, int nHalo)
  {
    int ndof = (N+1)*(N+1)*nHalo;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nVar,1);
    dim3 nthreads(threads_per_block,1,1);
    HaloUnpack_3D<<<nblocks,nthreads>>>(recvBuf, extBoundary, haloSides, N, nVar, nEl, nHalo);
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

    uint32_t iq= threadIdx.x;

    if( iq< nq ){
      uint32_t iel = blockIdx.x;
      uint32_t nel = gridDim.x;
      uint32_t ivar = blockIdx.y;
      uint32_t nvar = gridDim.y;
      real Fx = vector[iq+ nq*(iel + nel*(ivar))];
      real Fy = vector[iq+ nq*(iel + nel*(ivar + nvar))];
            
      vector[iq+ nq*(iel + nel*(ivar))] = dsdx[iq+ nq*(iel)]*Fx + // dsdx(...,0,0)*Fx
                                             dsdx[iq+ nq*(iel+nel)]*Fy; // dsdx(...,1,0)*Fy;

      vector[iq+ nq*(iel + nel*(ivar+nvar))] = dsdx[iq+ nq*(iel+nel*2)]*Fx + //dsdx(...,0,1)*Fx
                                                  dsdx[iq+ nq*(iel+nel*3)]*Fy;  //dsdx(...,1,1)*Fy
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

__global__ void ContravariantProjection_3D_gpukernel(real *vector, real *dsdx, int nq){

    // Grid-strided over quadrature points so that any polynomial degree N is
    // supported (nq = (N+1)^3 may exceed the thread-block size).
    for(uint32_t iq = threadIdx.x; iq < nq; iq += blockDim.x){
      uint32_t iel = blockIdx.x;
      uint32_t nel = gridDim.x;
      uint32_t ivar = blockIdx.y;
      uint32_t nvar = gridDim.y;

      real Fx = vector[iq+ nq*(iel + nel*(ivar))];
      real Fy = vector[iq+ nq*(iel + nel*(ivar + nvar))];
      real Fz = vector[iq+ nq*(iel + nel*(ivar + 2*nvar))];

      vector[iq+ nq*(iel + nel*(ivar))] = dsdx[iq+ nq*iel]*Fx + 
                                          dsdx[iq+ nq*(iel+nel)]*Fy + 
                                          dsdx[iq+ nq*(iel+2*nel)]*Fz;

      vector[iq+ nq*(iel + nel*(ivar + nvar))] = dsdx[iq+ nq*(iel+3*nel)]*Fx + 
                                                 dsdx[iq+ nq*(iel+4*nel)]*Fy + 
                                                 dsdx[iq+ nq*(iel+5*nel)]*Fz;

      vector[iq+ nq*(iel + nel*(ivar + 2*nvar))] = dsdx[iq+ nq*(iel+6*nel)]*Fx + 
                                                   dsdx[iq+ nq*(iel+7*nel)]*Fy + 
                                                   dsdx[iq+ nq*(iel+8*nel)]*Fz;
    }

}

extern "C"
{
  void ContravariantProjection_3D_gpu(real *vector, real *dsdx, int N, int nvar, int nel)
  {

    int nq = (N+1)*(N+1)*(N+1);
    // Grid-strided single launch, valid for any N. The mapped-vector divergence
    // reaches this only as the high-N fallback of
    // MappedContravariantDivergence_3D_gpu (N>=13); it remains a correct
    // standalone contravariant projection for any degree.
    ContravariantProjection_3D_gpukernel<<<dim3(nel,nvar,1), dim3(256,1,1), 0, 0>>>(vector, dsdx, nq);
  }
}

// Fused contravariant-projection + interior divergence for 3-D mapped vectors.
//
// Profiling of the LinearEuler3D forward step on MI300X showed the separate
// ContravariantProjection_3D kernel is the #1 hotspot (memory-bound, ~28% L2
// hit) because it streams the physical flux + 9 metric terms through global
// memory, writes the projected field back, and VectorDivergence_3D then re-reads
// that field from global. This kernel fuses the two: each node's contravariant
// components are computed ONCE during a shared-memory staging load (physical
// flux * dsdx), then the tensor contraction reads them from LDS. This removes
// the global write-back of the projected field and converts the contraction's
// repeated global re-reads into LDS reads. It leaves the input `f` unmodified
// (unlike the in-place ContravariantProjection).
//
// The arithmetic is IDENTICAL to ContravariantProjection_3D followed by
// VectorDivergence_3D (same products, same summation order), so `df` is
// bitwise-identical to the unfused chain. LDS use is (3*(N+1)^3 + (N+1)^2)
// reals; the launcher only dispatches this kernel when that fits the LDS budget
// (the Fortran caller falls back to the two-kernel path otherwise).
__global__ void MappedContravariantDivergence_3D_gpukernel(real *dsdx, real *A, real *f, real *df,
                                                           real *jacobian, int N, int nel, int nvar){

  int iel = blockIdx.x;
  int ivar = blockIdx.y;
  int nq = (N+1)*(N+1)*(N+1);

  extern __shared__ real s[];
  real *c1 = s;              // contravariant component 0, size nq
  real *c2 = s + nq;         // component 1
  real *c3 = s + 2*nq;       // component 2
  real *sA = s + 3*nq;       // derivative matrix, size (N+1)*(N+1)

  for(int idx = threadIdx.x; idx < (N+1)*(N+1); idx += blockDim.x){
    sA[idx] = A[idx];
  }

  // Stage the contravariant-projected flux into LDS (projection computed once
  // per node, identical order to ContravariantProjection_3D_gpukernel).
  for(int iq = threadIdx.x; iq < nq; iq += blockDim.x){
    real Fx = f[iq + nq*(iel + nel*(ivar))];
    real Fy = f[iq + nq*(iel + nel*(ivar + nvar))];
    real Fz = f[iq + nq*(iel + nel*(ivar + 2*nvar))];
    c1[iq] = dsdx[iq + nq*iel]*Fx + dsdx[iq + nq*(iel+nel)]*Fy + dsdx[iq + nq*(iel+2*nel)]*Fz;
    c2[iq] = dsdx[iq + nq*(iel+3*nel)]*Fx + dsdx[iq + nq*(iel+4*nel)]*Fy + dsdx[iq + nq*(iel+5*nel)]*Fz;
    c3[iq] = dsdx[iq + nq*(iel+6*nel)]*Fx + dsdx[iq + nq*(iel+7*nel)]*Fy + dsdx[iq + nq*(iel+8*nel)]*Fz;
  }
  __syncthreads();

  // Tensor contraction from LDS (identical order to VectorDivergence_3D_gpukernel).
  for(int iq = threadIdx.x; iq < nq; iq += blockDim.x){
    int i = iq % (N+1);
    int j = (iq/(N+1)) % (N+1);
    int k = iq/(N+1)/(N+1);
    real acc = 0.0;
    for(int a = 0; a < N+1; a++){
      acc += sA[a + (N+1)*i]*c1[a + (N+1)*(j + (N+1)*k)]
           + sA[a + (N+1)*j]*c2[i + (N+1)*(a + (N+1)*k)]
           + sA[a + (N+1)*k]*c3[i + (N+1)*(j + (N+1)*a)];
    }
    // Optional Jacobian weight folded into the write (strong-form path passes a
    // non-null jacobian; the DG path passes null and applies /J in the fused
    // boundary-contribution kernel instead). Bitwise-identical to a separate
    // JacobianWeight pass.
    if(jacobian){
      df[SC_3D_INDEX(i,j,k,iel,ivar,N,nel)] = acc/jacobian[iq + nq*iel];
    }else{
      df[SC_3D_INDEX(i,j,k,iel,ivar,N,nel)] = acc;
    }
  }
}

// Defined in SELF_MatrixMultiply.cpp; used as the high-N fallback below.
extern "C" void VectorDivergence_3D_gpu(real *A, real *f, real *df, int N, int nvar, int nel);

// Maximum dynamic shared memory (bytes) a plain <<<...,smem,...>>> launch may
// request on the active device WITHOUT an opt-in (cudaFuncSetAttribute /
// hipFuncSetAttribute). This is 48 KiB on NVIDIA sm_70/sm_80 and 64 KiB on AMD
// gfx942; a launch requesting more fails (silently for sm_70). Queried once and
// cached, so the fused divergence kernel is used up to the device's real limit
// (falling back to the two-kernel path above it) rather than a fixed guess.
static size_t maxDynamicSharedBytes(){
  static size_t cached = 0;
  if(cached) return cached;
  int dev = 0, v = 0;
#ifdef __HIP_PLATFORM_AMD__
  CHECK(hipGetDevice(&dev));
  CHECK(hipDeviceGetAttribute(&v, hipDeviceAttributeMaxSharedMemoryPerBlock, dev));
#else
  CHECK(cudaGetDevice(&dev));
  CHECK(cudaDeviceGetAttribute(&v, cudaDevAttrMaxSharedMemoryPerBlock, dev));
#endif
  cached = (size_t)v;
  return cached;
}

extern "C"
{
  // If jacobian is non-null, the /J weight is folded into the epilogue write
  // (used by the strong-form path, which then skips the separate JacobianWeight
  // call). The DG path passes null and applies /J in
  // DG_BoundaryContribution_JacobianWeight_3D_gpu after the boundary terms.
  void MappedContravariantDivergence_3D_gpu(real *dsdx, real *A, real *f, real *df,
                                            real *jacobian, int N, int nvar, int nel)
  {
    int nq = (N+1)*(N+1)*(N+1);
    size_t smem = (size_t)(3*nq + (N+1)*(N+1))*sizeof(real);
    // Use the device's real per-block dynamic-shared limit (48 KiB on sm_70,
    // 64 KiB on gfx942); requests above it fall back to the two-kernel path.
    const size_t maxLDS = maxDynamicSharedBytes();
    if( smem <= maxLDS ){
      MappedContravariantDivergence_3D_gpukernel<<<dim3(nel,nvar,1), dim3(256,1,1), smem, 0>>>(dsdx,A,f,df,jacobian,N,nel,nvar);
    } else {
      // Fallback when the fused kernel's LDS request exceeds the ceiling:
      // in-place contravariant projection followed by the grid-strided
      // divergence, then the Jacobian weight if requested. Numerically identical
      // to the fused path.
      ContravariantProjection_3D_gpu(f, dsdx, N, nvar, nel);
      VectorDivergence_3D_gpu(A, f, df, N, nvar, nel);
      if(jacobian){
        JacobianWeight_3D_gpu(df, jacobian, N, nvar, nel);
      }
    }
  }
}