#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

__global__ void UpdateGAB2_Model2D_gpu(real *prevsol, real *solution, int m, int nPrev, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  if( m == 0 ){

    prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)] = solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

  }
  else if( m == 1 ) {

    solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)];

  }
  else {

    prevsol[SC_2D_INDEX(i,j,nVar+iVar,iEl,N,nPrev)] = prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)];

    prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)] = solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

    solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = 1.5*prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)]-
	    0.5*prevsol[SC_2D_INDEX(i,j,nVar+iVar,iEl,N,nPrev)];

  }

}

extern "C"
{
  void UpdateGAB2_Model2D_gpu_wrapper(real **prevsol, real **solution, int m, int nPrev, int N, int nVar, int nEl)
  {
    UpdateGAB2_Model2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*prevsol, *solution, m, nPrev, N, nVar);
  }
}

__global__ void UpdateGAB3_Model2D_gpu(real *prevsol, real *solution, int m, int nPrev, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  if( m == 0 ){

    prevsol[SC_2D_INDEX(i,j,nVar+iVar,iEl,N,nPrev)] = solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

  }
  else if( m == 1 ){

    prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)] = solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

  }
  else if( m == 2 ) {

    solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)];

  }
  else {

    prevsol[SC_2D_INDEX(i,j,2*nVar+iVar,iEl,N,nPrev)] = prevsol[SC_2D_INDEX(i,j,nVar+iVar,iEl,N,nPrev)];

    prevsol[SC_2D_INDEX(i,j,nVar+iVar,iEl,N,nPrev)] = prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)];

    solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = (23.0*prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)]-
	    16.0*prevsol[SC_2D_INDEX(i,j,nVar+iVar,iEl,N,nPrev)] +
	    5.0*prevsol[SC_2D_INDEX(i,j,2*nVar+iVar,iEl,N,nPrev)])/12.0;

  }

}

extern "C"
{
  void UpdateGAB3_Model2D_gpu_wrapper(real **prevsol, real **solution, int m, int nPrev, int N, int nVar, int nEl)
  {
    UpdateGAB3_Model2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*prevsol, *solution, m, nPrev, N, nVar);
  }
}

__global__ void UpdateGAB4_Model2D_gpu(real *prevsol, real *solution, int m, int nPrev, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  if( m == 0 ){

    prevsol[SC_2D_INDEX(i,j,2*nVar+iVar,iEl,N,nPrev)] = solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

  }
  else if( m == 1 ){

    prevsol[SC_2D_INDEX(i,j,nVar+iVar,iEl,N,nPrev)] = solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

  }
  else if( m == 2 ){

    prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)] = solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

  }
  else if( m == 3 ) {

    solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)];

  }
  else {


    prevsol[SC_2D_INDEX(i,j,3*nVar+iVar,iEl,N,nPrev)] = prevsol[SC_2D_INDEX(i,j,2*nVar+iVar,iEl,N,nPrev)];

    prevsol[SC_2D_INDEX(i,j,2*nVar+iVar,iEl,N,nPrev)] = prevsol[SC_2D_INDEX(i,j,nVar+iVar,iEl,N,nPrev)];

    prevsol[SC_2D_INDEX(i,j,nVar+iVar,iEl,N,nPrev)] = prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)];

    solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = (55.0*prevsol[SC_2D_INDEX(i,j,iVar,iEl,N,nPrev)]-
	    59.0*prevsol[SC_2D_INDEX(i,j,nVar+iVar,iEl,N,nPrev)] +
	    37.0*prevsol[SC_2D_INDEX(i,j,2*nVar+iVar,iEl,N,nPrev)]-
	    9.0*prevsol[SC_2D_INDEX(i,j,3*nVar+iVar,iEl,N,nPrev)])/24.0;

  }

}

extern "C"
{
  void UpdateGAB4_Model2D_gpu_wrapper(real **prevsol, real **solution, int m, int nPrev, int N, int nVar, int nEl)
  {
    UpdateGAB4_Model2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*prevsol, *solution, m, nPrev, N, nVar);
  }
}

__global__ void UpdateGRK_Model2D_gpu(real *grk, real *solution, real *dSdt, real rk_a, real rk_g, real dt, int nWork, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    grk[SC_2D_INDEX(i,j,iVar,iEl,N,nWork)] = rk_a*grk[SC_2D_INDEX(i,j,iVar,iEl,N,nWork)] + dSdt[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];
    solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] += rk_g*dt*grk[SC_2D_INDEX(i,j,iVar,iEl,N,nWork)];

}

extern "C"
{
  void UpdateGRK_Model2D_gpu_wrapper(real **grk, real **solution, real **dSdt, real rk_a, real rk_g, real dt, int nWork, int N, int nVar, int nEl)
  {
    UpdateGRK_Model2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*grk, *solution, *dSdt, rk_a, rk_g, dt, nWork, N, nVar);
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
