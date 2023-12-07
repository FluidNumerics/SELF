#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

__global__ void UpdateSolution_Model1D_gpu(real *solution, real *dSdt, real dt, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;

    solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)] += dt*dSdt[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

}

extern "C"
{
  void UpdateSolution_Model1D_gpu_wrapper(real **solution, real **dSdt, real dt, int N, int nVar, int nEl)
  {
    UpdateSolution_Model1D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,1,1), 0, 0>>>(*solution, *dSdt, dt, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void UpdateGAB2_Model1D_gpu(real *prevsol, real *solution, int m, int nPrev, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;

  if( m == 0 ){

    prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

  }
  else if( m == 1 ) {

    solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

  }
  else if( m == 2 ) {

    prevsol[SC_1D_INDEX(i,iEl,nVar+iVar,N,nEl)] = prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

    prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

    solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = 1.5*prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)]-
	    0.5*prevsol[SC_1D_INDEX(i,iEl,nVar+iVar,N,nEl)];

  }

}

extern "C"
{
  void UpdateGAB2_Model1D_gpu_wrapper(real **prevsol, real **solution, int m, int nPrev, int N, int nVar, int nEl)
  {
    UpdateGAB2_Model1D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,1,1), 0, 0>>>(*prevsol, *solution, m, nPrev, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void UpdateGAB3_Model1D_gpu(real *prevsol, real *solution, int m, int nPrev, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;

  if( m == 0 ){

    prevsol[SC_1D_INDEX(i,iEl,nVar+iVar,N,nEl)] = solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

  }
  else if( m == 1 ){

    prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

  }
  else if( m == 2 ) {

    solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

  }
  else {

    prevsol[SC_1D_INDEX(i,iEl,2*nVar+iVar,N,nEl)] = prevsol[SC_1D_INDEX(i,iEl,nVar+iVar,N,nEl)];

    prevsol[SC_1D_INDEX(i,iEl,nVar+iVar,N,nEl)] = prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

    solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = (23.0*prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)]-
	    16.0*prevsol[SC_1D_INDEX(i,iEl,nVar+iVar,N,nEl)] +
	    5.0*prevsol[SC_1D_INDEX(i,iEl,2*nVar+iVar,N,nEl)])/12.0;

  }

}

extern "C"
{
  void UpdateGAB3_Model1D_gpu_wrapper(real **prevsol, real **solution, int m, int nPrev, int N, int nVar, int nEl)
  {
    UpdateGAB3_Model1D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,1,1), 0, 0>>>(*prevsol, *solution, m, nPrev, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void UpdateGAB4_Model1D_gpu(real *prevsol, real *solution, int m, int nPrev, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;

  if( m == 0 ){

    prevsol[SC_1D_INDEX(i,iEl,2*nVar+iVar,N,nEl)] = solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

  }
  else if( m == 1 ){

    prevsol[SC_1D_INDEX(i,iEl,nVar+iVar,N,nEl)] = solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

  }
  else if( m == 2 ){

    prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

  }
  else if( m == 3 ) {

    solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

  }
  else {


    prevsol[SC_1D_INDEX(i,iEl,3*nVar+iVar,N,nEl)] = prevsol[SC_1D_INDEX(i,iEl,2*nVar+iVar,N,nEl)];

    prevsol[SC_1D_INDEX(i,iEl,2*nVar+iVar,N,nEl)] = prevsol[SC_1D_INDEX(i,iEl,nVar+iVar,N,nEl)];

    prevsol[SC_1D_INDEX(i,iEl,nVar+iVar,N,nEl)] = prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

    solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = (55.0*prevsol[SC_1D_INDEX(i,iEl,iVar,N,nEl)]-
	    59.0*prevsol[SC_1D_INDEX(i,iEl,nVar+iVar,N,nEl)] +
	    37.0*prevsol[SC_1D_INDEX(i,iEl,2*nVar+iVar,N,nEl)]-
	    9.0*prevsol[SC_1D_INDEX(i,iEl,3*nVar+iVar,N,nEl)])/24.0;

  }

}

extern "C"
{
  void UpdateGAB4_Model1D_gpu_wrapper(real **prevsol, real **solution, int m, int nPrev, int N, int nVar, int nEl)
  {
    UpdateGAB4_Model1D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,1,1), 0, 0>>>(*prevsol, *solution, m, nPrev, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void UpdateGRK_Model1D_gpu(real *grk, real *solution, real *dSdt, real rk_a, real rk_g, real dt, int nWork, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;

    grk[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = rk_a*grk[SC_1D_INDEX(i,iEl,iVar,N,nEl)] + dSdt[SC_1D_INDEX(i,iEl,iVar,N,nEl)];
    solution[SC_1D_INDEX(i,iEl,iVar,N,nEl)] += rk_g*dt*grk[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

}

extern "C"
{
  void UpdateGRK_Model1D_gpu_wrapper(real **grk, real **solution, real **dSdt, real rk_a, real rk_g, real dt, int nWork, int N, int nVar, int nEl)
  {
    UpdateGRK_Model1D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,1,1), 0, 0>>>(*grk, *solution, *dSdt, rk_a, rk_g, dt, nWork, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void CalculateDSDt_Model1D_gpu(real *fluxDivergence, real *source, real *dSdt, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;

    dSdt[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = source[SC_1D_INDEX(i,iEl,iVar,N,nEl)]-
	    fluxDivergence[SC_1D_INDEX(i,iEl,iVar,N,nEl)];

}

extern "C"
{
  void CalculateDSDt_Model1D_gpu_wrapper(real **fluxDivergence, real **source, real **dSdt, int N, int nVar, int nEl)
  {
    CalculateDSDt_Model1D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,1,1), 0, 0>>>(*fluxDivergence, *source, *dSdt, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}
