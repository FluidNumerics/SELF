#include "SELF_GPU_Macros.h"


__global__ void setboundarycondition_shallow_water_2d_gpukernel(real *extBoundary, real *boundary, int *sideInfo, real *nhat, int N, int nEl, int nvar){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*4*nEl;

  if(idof < ndof){
    uint32_t i = idof % (N+1);
    uint32_t s1 = (idof/(N+1)) % 4;
    uint32_t e1 = idof/(N+1)/4;
    uint32_t e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
    uint32_t bcid = sideInfo[INDEX3(4,s1,e1,5,4)];
    if( e2 == 0){
      // if(bcid == SELF_BC_RADIATION){
        extBoundary[SCB_2D_INDEX(i,s1,e1,0,N,nEl)] = 0.0;
        extBoundary[SCB_2D_INDEX(i,s1,e1,1,N,nEl)] = 0.0;
        extBoundary[SCB_2D_INDEX(i,s1,e1,2,N,nEl)] = 0.0;
      // }
      // else if (bcid == SELF_BC_NONORMALFLOW){
        // real nx = nhat[VEB_2D_INDEX(i,s1,e1,0,0,N,nEl,1)];
        // real ny = nhat[VEB_2D_INDEX(i,s1,e1,0,1,N,nEl,1)];
        // real u = boundary[SCB_2D_INDEX(i,s1,e1,0,N,nEl)];
        // real v = boundary[SCB_2D_INDEX(i,s1,e1,1,N,nEl)];
        // real eta = boundary[SCB_2D_INDEX(i,s1,e1,2,N,nEl)];
        
        // extBoundary[SCB_2D_INDEX(i,s1,e1,0,N,nEl)] = (ny * ny - nx * nx) * u - 2 * nx * ny * v;
        // extBoundary[SCB_2D_INDEX(i,s1,e1,1,N,nEl)] = (nx * nx - ny * ny) * v - 2 * nx * ny * u;
        // extBoundary[SCB_2D_INDEX(i,s1,e1,2,N,nEl)] = eta;
      // }
    }
  }
}

extern "C" 
{
  void setboundarycondition_shallow_water_2d_gpu(real *extBoundary, real *boundary, int *sideInfo, real *nhat, int N, int nEl, int nvar){
    int threads_per_block = 256;
    int ndof = (N+1)*4*nEl;
    int nblocks_x = ndof/threads_per_block +1;

    dim3 nblocks(nblocks_x,1,1);
    dim3 nthreads(threads_per_block,1,1);

	setboundarycondition_shallow_water_2d_gpukernel<<<nblocks,nthreads, 0, 0>>>(extBoundary,boundary,sideInfo,nhat,N,nEl,nvar);
  }
}

__global__ void fluxmethod_shallow_water_2d_gpukernel(real *solution, real *flux, real g, real H, int ndof){
  uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;

  if( i < ndof ){
    real u = solution[i];
    real v = solution[i + ndof];
    real eta = solution[i + 2*ndof];
    flux[i] = g*eta; // x-component of u
    flux[i+ndof] = 0.0; // x-component of v
    flux[i+2*ndof] = H*u; // x-component of eta
    flux[i+3*ndof] = 0.0; // y-component of u
    flux[i+4*ndof] = g*eta; // y-component of v
    flux[i+5*ndof] = H*v; // y-component of eta

    // flux[i] = 0.0; // x-component of u
    // flux[i+ndof] = 0.0; // x-component of v
    // flux[i+2*ndof] = 0.0; // x-component of eta
    // flux[i+3*ndof] = 0.0; // y-component of u
    // flux[i+4*ndof] = 0.0; // y-component of v
    // flux[i+5*ndof] = 0.0; // y-component of eta
  }

}
extern "C"
{
  void fluxmethod_shallow_water_2d_gpu(real *solution, real *flux, real g, real H, int N, int nEl, int nvar){
    int ndof = (N+1)*(N+1)*nEl;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
    fluxmethod_shallow_water_2d_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(solution,flux,g,H,ndof);
  }
}

__global__ void riemannsolver_shallow_water_2d_gpukernel(real *fb, real *fextb, real *nhat, real *nscale, real *flux, real g, real H, int N, int nEl, int nvar){
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*4*nEl;

  if( idof < ndof ){

    // uint32_t i = idof % (N+1);
    // uint32_t j = (idof/(N+1)) % 4;
    // uint32_t iel = idof/(N+1)/4;

    real nx = nhat[idof];
    real ny = nhat[idof+ndof];
    real nmag = nscale[idof];

    real uL = fb[idof];
    real vL = fb[idof + ndof];
    real etaL = fb[idof + 2*ndof];
    real uR = fextb[idof];
    real vR = fextb[idof + ndof];
    real etaR = fextb[idof + 2*ndof];

    real unL = uL * nx + vL * ny;
    real unR = uR * nx + vR * ny;

    real c = sqrt(g * H);

    flux[idof]          = 0.5 * (g * (etaL + etaR) + c * (unL - unR)) * nx * nmag;
    flux[idof + ndof]   = 0.5 * (g * (etaL + etaR) + c * (unL - unR)) * ny * nmag;
    flux[idof + 2*ndof] = 0.5 * (H * (unL + unR) + c * (etaL - etaR)) * nmag;
    // flux[idof] = 0.0;
    // flux[idof + ndof] = 0.0;
    // flux[idof + 2*ndof] = 0.0;
  }

}
extern "C"
{
  void riemannsolver_shallow_water_2d_gpu(real *fb, real *fextb, real *nhat, real *nscale, real *flux, real g, real H, int N, int nEl, int nvar){
    int threads_per_block = 256;
    uint32_t ndof = (N+1)*4*nEl;
    int nblocks_x = ndof/threads_per_block +1;

    dim3 nblocks(nblocks_x,1,1);
    dim3 nthreads(threads_per_block,1,1);

    riemannsolver_shallow_water_2d_gpukernel<<<nblocks,nthreads>>>(fb,fextb,nhat,nscale,flux,g,H,N,nEl,nvar);
  }

}

