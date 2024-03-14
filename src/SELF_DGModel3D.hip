#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"
#include "SELF_Model.h"



extern "C"
{
  void UpdateSolution_DGModel3D_gpu(real *solution, real *dSdt, real dt, int N, int nVar, int nEl)
  {
    uint32_t ndof = (N+1)*(N+1)*(N+1)*(nEl)*(nVar);
    uint32_t nthreads = 256;
    uint32_t nblocks_x = ndof/nthreads + 1;
    UpdateSolution_Model<<<dim3(nblocks_x,1), dim3(nthreads,1,1), 0, 0>>>(solution, dSdt, dt, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

extern "C"
{
  void UpdateGAB2_DGModel3D_gpu(real *prevsol, real *solution, int m, int nPrev, int N, int nVar, int nEl)
  {
    uint32_t ndof = (N+1)*(N+1)*(N+1)*(nEl)*(nVar);
    uint32_t nthreads = 256;
    uint32_t nblocks_x = ndof/nthreads + 1;
    UpdateGAB2_Model<<<dim3(nblocks_x,1), dim3(nthreads,1,1), 0, 0>>>(prevsol, solution, m, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

extern "C"
{
  void UpdateGAB3_DGModel3D_gpu(real *prevsol, real *solution, int m, int nPrev, int N, int nVar, int nEl)
  {
    uint32_t ndof = (N+1)*(N+1)*(N+1)*(nEl)*(nVar);
    uint32_t nthreads = 256;
    uint32_t nblocks_x = ndof/nthreads + 1;
    UpdateGAB3_Model<<<dim3(nblocks_x,1), dim3(nthreads,1,1), 0, 0>>>(prevsol, solution, m, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

extern "C"
{
  void UpdateGAB4_DGModel3D_gpu(real *prevsol, real *solution, int m, int nPrev, int N, int nVar, int nEl)
  {
    uint32_t ndof = (N+1)*(N+1)*(N+1)*(nEl)*(nVar);
    uint32_t nthreads = 256;
    uint32_t nblocks_x = ndof/nthreads + 1;
    UpdateGAB4_Model<<<dim3(nblocks_x,1), dim3(nthreads,1,1), 0, 0>>>(prevsol, solution, m, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

extern "C"
{
  void UpdateGRK_DGModel3D_gpu(real *grk, real *solution, real *dSdt, real rk_a, real rk_g, real dt, int nWork, int N, int nVar, int nEl)
  {
    uint32_t ndof = (N+1)*(N+1)*(N+1)*(nEl)*(nVar);
    uint32_t nthreads = 256;
    uint32_t nblocks_x = ndof/nthreads + 1;
    UpdateGRK_Model<<<dim3(nblocks_x,1), dim3(nthreads,1,1), 0, 0>>>(grk, solution, dSdt, rk_a, rk_g, dt, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

extern "C"
{
  void CalculateDSDt_DGModel3D_gpu(real *fluxDivergence, real *source, real *dSdt, int N, int nVar, int nEl)
  {
    uint32_t ndof = (N+1)*(N+1)*(N+1)*(nEl)*(nVar);
    uint32_t nthreads = 256;
    uint32_t nblocks_x = ndof/nthreads + 1;
    CalculateDSDt_Model<<<dim3(nblocks_x,1), dim3(nthreads,1,1), 0, 0>>>(fluxDivergence, source, dSdt, ndof);
    HIP_SAFE_CALL(hipGetLastError());
  }
}
