#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"
#include <cstdio>


__global__ void SetBoundaryCondition_LinearShallowWater_gpu( real *solution, real *extBoundary, real *nHat, int *sideInfo, int N, int nVar ){

  size_t iSide = blockIdx.x+1;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

  int bcid = sideInfo[4+5*(iSide-1 +4*(iEl))];
  int e2 = sideInfo[2+5*(iSide-1 + 4*(iEl))];

  if( e2 == 0 ){

    if( bcid == SELF_BC_RADIATION ){

      extBoundary[SCB_2D_INDEX(i,0,iSide,iEl,N,nVar)] = 0.0;
      extBoundary[SCB_2D_INDEX(i,1,iSide,iEl,N,nVar)] = 0.0;
      extBoundary[SCB_2D_INDEX(i,2,iSide,iEl,N,nVar)] = 0.0;

    } else if ( bcid == SELF_BC_NONORMALFLOW ){

      real nx = nHat[VEB_2D_INDEX(1,i,0,iSide,iEl,N,1)];
      real ny = nHat[VEB_2D_INDEX(2,i,0,iSide,iEl,N,1)];
      real u = solution[SCB_2D_INDEX(i,0,iSide,iEl,N,nVar)];
      real v = solution[SCB_2D_INDEX(i,1,iSide,iEl,N,nVar)];
      real eta = solution[SCB_2D_INDEX(i,2,iSide,iEl,N,nVar)];

      extBoundary[SCB_2D_INDEX(i,0,iSide,iEl,N,nVar)] = (ny*ny-nx*ny)*u-2.0*nx*ny*v;
      extBoundary[SCB_2D_INDEX(i,0,iSide,iEl,N,nVar)] = (nx*nx-ny*ny)*v-2.0*nx*ny*u;
      extBoundary[SCB_2D_INDEX(i,0,iSide,iEl,N,nVar)] = eta;

    } else {

      extBoundary[SCB_2D_INDEX(i,0,iSide,iEl,N,nVar)] = 0.0;
      extBoundary[SCB_2D_INDEX(i,1,iSide,iEl,N,nVar)] = 0.0;
      extBoundary[SCB_2D_INDEX(i,2,iSide,iEl,N,nVar)] = 0.0;

    }
  }

}

extern "C"
{
  void SetBoundaryCondition_LinearShallowWater_gpu_wrapper( real **solution, real **extBoundary, real **nHat, int **sideInfo, int N, int nVar, int nEl)
  {
    SetBoundaryCondition_LinearShallowWater_gpu<<<dim3(4,nEl,1), dim3(N+1,1,1), 0, 0>>>(*solution, *extBoundary, *nHat, *sideInfo, N, nVar);
  }
}

__global__ void Source_LinearShallowWater_gpu(real *source, real *solution, real *f, int N, int nVar){

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    if( iVar == 0 ){
      source[SC_2D_INDEX(i,j,iVar,iEl,N,1)] = f[SC_2D_INDEX(i,j,1,iEl,N,1)]*solution[SC_2D_INDEX(i,j,1,iEl,N,nVar)];
    } else if ( iVar == 1) {
      source[SC_2D_INDEX(i,j,iVar,iEl,N,1)] = -f[SC_2D_INDEX(i,j,1,iEl,N,1)]*solution[SC_2D_INDEX(i,j,0,iEl,N,nVar)];
    } else if ( iVar == 2) {
      source[SC_2D_INDEX(i,j,iVar,iEl,N,1)] = 0.0;
    }
}

extern "C"
{
  void Source_LinearShallowWater_gpu_wrapper(real **source, real **solution, real **f, int N, int nVar, int nEl)
  {
    Source_LinearShallowWater_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*source, *solution, *f, N, nVar);
  }
}


__global__ void Flux_LinearShallowWater_gpu(real *flux, real *solution, real g, real H, int N, int nVar){

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    if( iVar == 0 ){
      flux[VE_2D_INDEX(1,i,j,iVar,iEl,N,1)] = g*solution[SC_2D_INDEX(i,j,2,iEl,N,nVar)];
      flux[VE_2D_INDEX(2,i,j,iVar,iEl,N,1)] = 0.0;
    } else if ( iVar == 1) {
      flux[VE_2D_INDEX(1,i,j,iVar,iEl,N,1)] = 0.0;
      flux[VE_2D_INDEX(2,i,j,iVar,iEl,N,1)] = g*solution[SC_2D_INDEX(i,j,2,iEl,N,nVar)];
    } else if ( iVar == 2) {
      flux[VE_2D_INDEX(1,i,j,iVar,iEl,N,1)] = H*solution[SC_2D_INDEX(i,j,0,iEl,N,nVar)];
      flux[VE_2D_INDEX(2,i,j,iVar,iEl,N,1)] = H*solution[SC_2D_INDEX(i,j,1,iEl,N,nVar)];
    }
}

extern "C"
{
  void Flux_LinearShallowWater_gpu_wrapper(real **flux, real **solution, real g, real H, int N, int nVar, int nEl)
  {
    Flux_LinearShallowWater_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*flux, *solution, g, H, N, nVar);
  }
}

__global__ void RiemannSolver_LinearShallowWater_gpu( real *flux, real *solution, real *extBoundary, real *nHat, real *nScale, real g, real H, int N, int nVar ){

  size_t iSide = blockIdx.x+1;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

  real nx = nHat[VEB_2D_INDEX(1,i,0,iSide,iEl,N,1)];
  real ny = nHat[VEB_2D_INDEX(2,i,0,iSide,iEl,N,1)];
	
  real unL = solution[SCB_2D_INDEX(i,0,iSide,iEl,N,nVar)]*nx+
             solution[SCB_2D_INDEX(i,1,iSide,iEl,N,nVar)]*ny;

  real unR = extBoundary[SCB_2D_INDEX(i,0,iSide,iEl,N,nVar)]*nx+
            extBoundary[SCB_2D_INDEX(i,1,iSide,iEl,N,nVar)]*ny;

  real etaL = solution[SCB_2D_INDEX(i,2,iSide,iEl,N,nVar)];
  real etaR = extBoundary[SCB_2D_INDEX(i,2,iSide,iEl,N,nVar)];

  real c = sqrt(g*H);

  real wL = 0.5*(unL/g + etaL/c);
  real wR = 0.5*(unR/g - etaR/c);

  real nmag = nScale[SCB_2D_INDEX(i,0,iSide,iEl,N,1)];

  flux[SCB_2D_INDEX(i,0,iSide,iEl,N,nVar)] = g*c*(wL-wR)*nx*nmag;
  flux[SCB_2D_INDEX(i,1,iSide,iEl,N,nVar)] = g*c*(wL-wR)*ny*nmag;
  flux[SCB_2D_INDEX(i,2,iSide,iEl,N,nVar)] = c*c*(wL+wR)*nmag;

}

extern "C"
{
  void RiemannSolver_LinearShallowWater_gpu_wrapper( real **flux, real **solution, real **extBoundary, real **nHat, real **nScale, real g, real H, int N, int nVar, int nEl)
  {
    RiemannSolver_LinearShallowWater_gpu<<<dim3(4,nEl,1), dim3(N+1,1,1), 0, 0>>>(*flux, *solution, *extBoundary, *nHat, *nScale, g, H, N, nVar);
  }
}
