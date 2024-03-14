#ifndef SELF_MODEL_H
#define SELF_MODEL_H

#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

__global__ void UpdateSolution_Model(real *solution, real *dSdt, real dt, uint32_t ndof);
__global__ void UpdateGAB2_Model(real *prevsol, real *solution, int m, uint32_t ndof);
__global__ void UpdateGAB3_Model(real *prevsol, real *solution, int m, uint32_t ndof);
__global__ void UpdateGAB4_Model(real *prevsol, real *solution, int m, uint32_t ndof);
__global__ void UpdateGRK_Model(real *grk, real *solution, real *dSdt, real rk_a, real rk_g, real dt, uint32_t ndof);
__global__ void CalculateDSDt_Model(real *fluxDivergence, real *source, real *dSdt, uint32_t ndof);

#endif
