#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

__global__ void UpdateSolution_Model(real *solution, real *dSdt, real dt, uint32_t ndof){

  size_t i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i < ndof ){
    solution[i] += dt*dSdt[i];
  }

}

__global__ void UpdateGAB2_Model(real *prevsol, real *solution, int m, uint32_t ndof){

  size_t i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i < ndof ){

    if( m == 0 ){

      prevsol[i] = solution[i];

    }
    else if( m == 1 ) {

      solution[i] = prevsol[i];

    }
    else if( m == 2 ) {

      prevsol[i + ndof] = prevsol[i];
      prevsol[i] = solution[i];
      solution[i] = 1.5*prevsol[i]-0.5*prevsol[i+ndof];

    }
  }

}

__global__ void UpdateGAB3_Model(real *prevsol, real *solution, int m, uint32_t ndof){

  size_t i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i < ndof ){

    if( m == 0 ){

      prevsol[i+ndof] = solution[i];

    }
    else if( m == 1 ){

      prevsol[i] = solution[i];

    }
    else if( m == 2 ) {

      solution[i] = prevsol[i];

    }
    else {

      prevsol[i+2*ndof] = prevsol[i+ndof];
      prevsol[i+ndof] = prevsol[i];
      solution[i] = (23.0*prevsol[i]-16.0*prevsol[i+ndof] + 5.0*prevsol[i+2*ndof])/12.0;

    }
  }

}

__global__ void UpdateGAB4_Model(real *prevsol, real *solution, int m, uint32_t ndof){

  size_t i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i < ndof ){

    if( m == 0 ){

      prevsol[i+2*ndof] = solution[i];

    }
    else if( m == 1 ){

      prevsol[i+ndof] = solution[i];

    }
    else if( m == 2 ){

      prevsol[i] = solution[i];

    }
    else if( m == 3 ) {

      solution[i] = prevsol[i];

    }
    else {


      prevsol[i+3*ndof] = prevsol[i+2*ndof];
      prevsol[i+2*ndof] = prevsol[i+ndof];
      prevsol[i+ndof] = prevsol[i];
      solution[i] = (55.0*prevsol[i]-59.0*prevsol[i+ndof]+37.0*prevsol[i+2*ndof]-9.0*prevsol[i+3*ndof])/24.0;

    }
  }

}

__global__ void UpdateGRK_Model(real *grk, real *solution, real *dSdt, real rk_a, real rk_g, real dt, uint32_t ndof){

  size_t i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i < ndof ){
    grk[i] = rk_a*grk[i] + dSdt[i];
    solution[i] += rk_g*dt*grk[i];
  }

}

__global__ void CalculateDSDt_Model(real *fluxDivergence, real *source, real *dSdt, uint32_t ndof){

  size_t i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i < ndof ){
    dSdt[i] = source[i]-fluxDivergence[i];
  }

}