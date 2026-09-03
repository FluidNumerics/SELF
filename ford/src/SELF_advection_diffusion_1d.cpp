#include "SELF_GPU_Macros.h"


__global__ void setboundarycondition_advection_diffusion_1d_gpukernel(real *extBoundary, real *boundary, int nel, int nvar){

  uint32_t ivar = threadIdx.x + blockIdx.x*blockDim.x;
  if(ivar < nvar){
    extBoundary[SCB_1D_INDEX(0,0,ivar,nel)] = boundary[SCB_1D_INDEX(1,nel-1,ivar,nel)];
    extBoundary[SCB_1D_INDEX(1,nel-1,ivar,nel)] = boundary[SCB_1D_INDEX(0,0,ivar,nel)];
  }
  
}

extern "C" 
{
  void setboundarycondition_advection_diffusion_1d_gpu(real *extBoundary, real *boundary, int nel, int nvar){
    int threads_per_block = 64;
    int nblocks_x = nvar/threads_per_block +1;
	  setboundarycondition_advection_diffusion_1d_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(extBoundary,boundary,nel,nvar);
  }
}

__global__ void fluxmethod_advection_diffusion_1d_gpukernel(real *solution, real *solutiongradient, real *flux, real u, real nu, int ndof){
  uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;

  if( i < ndof ){
    flux[i] = u*solution[i] - nu*solutiongradient[i];
  }

}
extern "C"
{
  void fluxmethod_advection_diffusion_1d_gpu(real *solution, real *solutiongradient, real *flux, real u, real nu, int ndof){
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
    fluxmethod_advection_diffusion_1d_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(solution,solutiongradient,flux,u,nu,ndof);
  }

}

__global__ void boundaryflux_advection_diffusion_1d_gpukernel(real *fb, real *fextb, real *dfavg, real *flux, real u, real nu, int ndof){
  uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;
  // when i is even, we are looking at the left side of the element and the boundary normal is negative
  // when i is odd, we are looking at the right side of the element and boundary normal is positive
  //
  //   i%2 is 0 when i is even
  //          1 when i is odd
  //
  //   2*(i%2) is 0 when i is even
  //              2 when i is odd
  //
  //   2*(i%2)-1 is -1 when i is even
  //                 1 when i is odd
  real nhat = 2.0*(i%2)-1.0;

  if( i < ndof ){
    flux[i] = 0.5*(u*nhat*(fb[i]+fextb[i]) + fabsf(u*nhat)*(fb[i]-fextb[i])) - nu*dfavg[i]*nhat;
  }

}
extern "C"
{
  void boundaryflux_advection_diffusion_1d_gpu(real *fb, real *fextb, real *dfavg, real *flux, real u, real nu, int ndof){
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
    boundaryflux_advection_diffusion_1d_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(fb,fextb,dfavg,flux,u,nu,ndof);
  }

}

