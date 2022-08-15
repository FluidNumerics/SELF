#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

__global__ void UpdateSolution_Model1D_gpu(real *solution, real *dSdt, real dt, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

    solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)] += dt*dSdt[SC_1D_INDEX(i,iVar,iEl,N,nVar)];

}

extern "C"
{
  void UpdateSolution_Model1D_gpu_wrapper(real **solution, real **dSdt, real dt, int N, int nVar, int nEl)
  {
    UpdateSolution_Model1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*solution, *dSdt, dt, N, nVar);
  }
}

__global__ void UpdateGRK_Model1D_gpu(real *grk, real *solution, real *dSdt, real rk_a, real rk_g, real dt, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

    grk[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = rk_a*grk[SC_1D_INDEX(i,iVar,iEl,N,nVar)] + dSdt[SC_1D_INDEX(i,iVar,iEl,N,nVar)];
    solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)] += rk_g*dt*grk[SC_1D_INDEX(i,iVar,iEl,N,nVar)];

}

extern "C"
{
  void UpdateGRK_Model1D_gpu_wrapper(real **grk, real **solution, real **dSdt, real rk_a, real rk_g, real dt, int N, int nVar, int nEl)
  {
    UpdateGRK_Model1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*grk, *solution, *dSdt, rk_a, rk_g, dt, N, nVar);
  }
}

__global__ void CalculateDSDt_Model1D_gpu(real *fluxDivergence, real *source, real *dSdt, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

    dSdt[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = source[SC_1D_INDEX(i,iVar,iEl,N,nVar)]-
	    fluxDivergence[SC_1D_INDEX(i,iVar,iEl,N,nVar)];

}

extern "C"
{
  void CalculateDSDt_Model1D_gpu_wrapper(real **fluxDivergence, real **source, real **dSdt, int N, int nVar, int nEl)
  {
    CalculateDSDt_Model1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*fluxDivergence, *source, *dSdt, N, nVar);
  }
}

__global__ void UpdateGRK_Model2D_gpu(real *grk, real *solution, real *dSdt, real rk_a, real rk_g, real dt, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    grk[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = rk_a*grk[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] + dSdt[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];
    solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] += rk_g*dt*grk[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

}

extern "C"
{
  void UpdateGRK_Model2D_gpu_wrapper(real **grk, real **solution, real **dSdt, real rk_a, real rk_g, real dt, int N, int nVar, int nEl)
  {
    UpdateGRK_Model2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*grk, *solution, *dSdt, rk_a, rk_g, dt, N, nVar);
  }
}

__global__ void UpdateSolution_Model2D_gpu(real *solution, real *dSdt, real dt, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] += dt*dSdt[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

}

extern "C"
{
  void UpdateSolution_Model2D_gpu_wrapper(real **solution, real **dSdt, real dt, int N, int nVar, int nEl)
  {
    UpdateSolution_Model2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*solution, *dSdt, dt, N, nVar);
  }
}

__global__ void CalculateDSDt_Model2D_gpu(real *fluxDivergence, real *source, real *dSdt, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    dSdt[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = source[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]-
	    fluxDivergence[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

}

extern "C"
{
  void CalculateDSDt_Model2D_gpu_wrapper(real **fluxDivergence, real **source, real **dSdt, int N, int nVar, int nEl)
  {
    CalculateDSDt_Model2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*fluxDivergence, *source, *dSdt, N, nVar);
  }
}

__global__ void CalculateDSDt_Model3D_gpu(real *fluxDivergence, real *source, real *dSdt, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    dSdt[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)] = source[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]-
	    fluxDivergence[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

}

extern "C"
{
  void CalculateDSDt_Model3D_gpu_wrapper(real **fluxDivergence, real **source, real **dSdt, int N, int nVar, int nEl)
  {
    CalculateDSDt_Model3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*fluxDivergence, *source, *dSdt, N, nVar);
  }
}
