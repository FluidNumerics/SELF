#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

__global__ void UpdateSolution_DGModel1D_gpu(real *solution, real *dSdt, real dt, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

    solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)] += dt*dSdt[SC_1D_INDEX(i,iVar,iEl,N,nVar)];

}

extern "C"
{
  void UpdateSolution_DGModel1D_gpu_wrapper(real **solution, real **dSdt, real dt, int N, int nVar, int nEl)
  {
    UpdateSolution_DGModel1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*solution, *dSdt, dt, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void UpdateGAB2_DGModel1D_gpu(real *prevsol, real *solution, int m, int nPrev, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

  if( m == 0 ){

    prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)] = solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)];

  }
  else if( m == 1 ) {

    solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)];

  }
  else if( m == 2 ) {

    prevsol[SC_1D_INDEX(i,nVar+iVar,iEl,N,nPrev)] = prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)];

    prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)] = solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)];

    solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = 1.5*prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)]-
	    0.5*prevsol[SC_1D_INDEX(i,nVar+iVar,iEl,N,nPrev)];

  }

}

extern "C"
{
  void UpdateGAB2_DGModel1D_gpu_wrapper(real **prevsol, real **solution, int m, int nPrev, int N, int nVar, int nEl)
  {
    UpdateGAB2_DGModel1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*prevsol, *solution, m, nPrev, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void UpdateGAB3_DGModel1D_gpu(real *prevsol, real *solution, int m, int nPrev, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

  if( m == 0 ){

    prevsol[SC_1D_INDEX(i,nVar+iVar,iEl,N,nPrev)] = solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)];

  }
  else if( m == 1 ){

    prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)] = solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)];

  }
  else if( m == 2 ) {

    solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)];

  }
  else {

    prevsol[SC_1D_INDEX(i,2*nVar+iVar,iEl,N,nPrev)] = prevsol[SC_1D_INDEX(i,nVar+iVar,iEl,N,nPrev)];

    prevsol[SC_1D_INDEX(i,nVar+iVar,iEl,N,nPrev)] = prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)];

    solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = (23.0*prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)]-
	    16.0*prevsol[SC_1D_INDEX(i,nVar+iVar,iEl,N,nPrev)] +
	    5.0*prevsol[SC_1D_INDEX(i,2*nVar+iVar,iEl,N,nPrev)])/12.0;

  }

}

extern "C"
{
  void UpdateGAB3_DGModel1D_gpu_wrapper(real **prevsol, real **solution, int m, int nPrev, int N, int nVar, int nEl)
  {
    UpdateGAB3_DGModel1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*prevsol, *solution, m, nPrev, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void UpdateGAB4_DGModel1D_gpu(real *prevsol, real *solution, int m, int nPrev, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

  if( m == 0 ){

    prevsol[SC_1D_INDEX(i,2*nVar+iVar,iEl,N,nPrev)] = solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)];

  }
  else if( m == 1 ){

    prevsol[SC_1D_INDEX(i,nVar+iVar,iEl,N,nPrev)] = solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)];

  }
  else if( m == 2 ){

    prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)] = solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)];

  }
  else if( m == 3 ) {

    solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)];

  }
  else {


    prevsol[SC_1D_INDEX(i,3*nVar+iVar,iEl,N,nPrev)] = prevsol[SC_1D_INDEX(i,2*nVar+iVar,iEl,N,nPrev)];

    prevsol[SC_1D_INDEX(i,2*nVar+iVar,iEl,N,nPrev)] = prevsol[SC_1D_INDEX(i,nVar+iVar,iEl,N,nPrev)];

    prevsol[SC_1D_INDEX(i,nVar+iVar,iEl,N,nPrev)] = prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)];

    solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = (55.0*prevsol[SC_1D_INDEX(i,iVar,iEl,N,nPrev)]-
	    59.0*prevsol[SC_1D_INDEX(i,nVar+iVar,iEl,N,nPrev)] +
	    37.0*prevsol[SC_1D_INDEX(i,2*nVar+iVar,iEl,N,nPrev)]-
	    9.0*prevsol[SC_1D_INDEX(i,3*nVar+iVar,iEl,N,nPrev)])/24.0;

  }

}

extern "C"
{
  void UpdateGAB4_DGModel1D_gpu_wrapper(real **prevsol, real **solution, int m, int nPrev, int N, int nVar, int nEl)
  {
    UpdateGAB4_DGModel1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*prevsol, *solution, m, nPrev, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void UpdateGRK_DGModel1D_gpu(real *grk, real *solution, real *dSdt, real rk_a, real rk_g, real dt, int nWork, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

    grk[SC_1D_INDEX(i,iVar,iEl,N,nWork)] = rk_a*grk[SC_1D_INDEX(i,iVar,iEl,N,nWork)] + dSdt[SC_1D_INDEX(i,iVar,iEl,N,nVar)];
    solution[SC_1D_INDEX(i,iVar,iEl,N,nVar)] += rk_g*dt*grk[SC_1D_INDEX(i,iVar,iEl,N,nWork)];

}

extern "C"
{
  void UpdateGRK_DGModel1D_gpu_wrapper(real **grk, real **solution, real **dSdt, real rk_a, real rk_g, real dt, int nWork, int N, int nVar, int nEl)
  {
    UpdateGRK_DGModel1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*grk, *solution, *dSdt, rk_a, rk_g, dt, nWork, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void CalculateDSDt_DGModel1D_gpu(real *fluxDivergence, real *source, real *dSdt, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

    dSdt[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = source[SC_1D_INDEX(i,iVar,iEl,N,nVar)]-
	    fluxDivergence[SC_1D_INDEX(i,iVar,iEl,N,nVar)];

}

extern "C"
{
  void CalculateDSDt_DGModel1D_gpu_wrapper(real **fluxDivergence, real **source, real **dSdt, int N, int nVar, int nEl)
  {
    CalculateDSDt_DGModel1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*fluxDivergence, *source, *dSdt, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }
}
