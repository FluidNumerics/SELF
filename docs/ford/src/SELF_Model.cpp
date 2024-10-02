#include "SELF_GPU_Macros.h"

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

extern "C"
{
  void UpdateSolution_gpu(real *solution, real *dSdt, real dt, int ndof)
  {
    uint32_t nthreads = 256;
    uint32_t nblocks_x = ndof/nthreads + 1;
    UpdateSolution_Model<<<dim3(nblocks_x,1), dim3(nthreads,1,1), 0, 0>>>(solution, dSdt, dt, ndof);
  }
}

extern "C"
{
  void UpdateGRK_gpu(real *grk, real *solution, real *dSdt, real rk_a, real rk_g, real dt, int ndof)
  {
    uint32_t nthreads = 256;
    uint32_t nblocks_x = ndof/nthreads + 1;
    UpdateGRK_Model<<<dim3(nblocks_x,1), dim3(nthreads,1,1), 0, 0>>>(grk, solution, dSdt, rk_a, rk_g, dt, ndof);
  }
}

extern "C"
{
  void CalculateDSDt_gpu(real *fluxDivergence, real *source, real *dSdt, int ndof)
  {
    uint32_t nthreads = 256;
    uint32_t nblocks_x = ndof/nthreads + 1;
    CalculateDSDt_Model<<<dim3(nblocks_x,1), dim3(nthreads,1,1), 0, 0>>>(fluxDivergence, source, dSdt, ndof);
  }
}

__global__ void GradientNormal_1d_gpukernel(real *fbn, real *fbavg, int ndof){
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
    fbn[i] = nhat*fbavg[i];
  }

}
extern "C"
{
  void GradientNormal_1d_gpu(real *fbn, real *fbavg, int ndof){
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
    GradientNormal_1d_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(fbn,fbavg,ndof);
  }

}