#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"
#include <cstdio>

__global__ void Source_CompressibleIdealGas2D_gpu(real *source, real *solution, real *environmentalsGradient, int N, int nVar)
{

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real rhou = solution[SC_2D_INDEX(i, j, 0, iEl, N, nVar)];
  real rhov = solution[SC_2D_INDEX(i, j, 1, iEl, N, nVar)];
  real rho = solution[SC_2D_INDEX(i, j, 2, iEl, N, nVar)];
  real gx = environmentalsGradient[VE_2D_INDEX(1, i, j, 0, iEl, N, 1)];
  real gy = environmentalsGradient[VE_2D_INDEX(2, i, j, 0, iEl, N, 1)];

  source[SC_2D_INDEX(i, j, 0, iEl, N, nVar)] = -rho * gx;
  source[SC_2D_INDEX(i, j, 1, iEl, N, nVar)] = -rho * gy;
  source[SC_2D_INDEX(i, j, 2, iEl, N, nVar)] = 0.0;
  source[SC_2D_INDEX(i, j, 3, iEl, N, nVar)] = -rhou * gx - rhou * gy;
}

extern "C"
{
  void Source_CompressibleIdealGas2D_gpu_wrapper(real **source, real **solution, real **environmentalsGradient, int N, int nVar, int nEl)
  {
    Source_CompressibleIdealGas2D_gpu<<<dim3(nVar, nEl, 1), dim3(N + 1, N + 1, 1), 0, 0>>>(*source, *solution, *environmentalsGradient, N, nVar);
  }
}

__global__ void ConservativeFlux_CompressibleIdealGas2D_gpu(real *flux, real *solution, real *primitive, int N, int nVar)
{

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  /*
  // Use shared memory to hold
  // u, v, rhou, rhov, rhoE, p
  __shared__ real u[64];
  __shared__ real v[64];
  __shared__ real rhou[64];
  __shared__ real rhov[64];
  __shared__ real rhoE[64];
  __shared__ real p[64];


  u[i+j*(N+1)] = primitive[SC_2D_INDEX(i,j,0,iEl,N,nVar)];
  v[i+j*(N+1)] = primitive[SC_2D_INDEX(i,j,1,iEl,N,nVar)];
  p[i+j*(N+1)] = primitive[SC_2D_INDEX(i,j,3,iEl,N,nVar)];
  rhou[i+j*(N+1)] = solution[SC_2D_INDEX(i,j,0,iEl,N,nVar)];
  rhov[i+j*(N+1)] = solution[SC_2D_INDEX(i,j,1,iEl,N,nVar)];
  rhoE[i+j*(N+1)] = solution[SC_2D_INDEX(i,j,3,iEl,N,nVar)];

  __syncthreads();

  */
  if (iVar == 0)
  { // rho*u

    for (int n = 0; n <= N; n++)
    {
      // rho*u*u + p
      flux[P2PVE_2D_INDEX(1, 1, n, i, j, iVar, iEl, N, nVar)] =
          primitive[SC_2D_INDEX(n, j, 0, iEl, N, nVar)] *       // u
              solution[SC_2D_INDEX(n, j, iVar, iEl, N, nVar)] + // rho*u
          primitive[SC_2D_INDEX(n, j, 3, iEl, N, nVar)];

      // rho*u*v
      flux[P2PVE_2D_INDEX(2, 1, n, i, j, iVar, iEl, N, nVar)] =
          primitive[SC_2D_INDEX(n, j, 1, iEl, N, nVar)] *  // v
          solution[SC_2D_INDEX(n, j, iVar, iEl, N, nVar)]; // rho*u
                                                           //
      // rho*u*u + p
      flux[P2PVE_2D_INDEX(1, 2, n, i, j, iVar, iEl, N, nVar)] =
          primitive[SC_2D_INDEX(i, n, 0, iEl, N, nVar)] *       // u
              solution[SC_2D_INDEX(i, n, iVar, iEl, N, nVar)] + // rho*u
          primitive[SC_2D_INDEX(n, j, 3, iEl, N, nVar)];

      // rho*u*v
      flux[P2PVE_2D_INDEX(2, 2, n, i, j, iVar, iEl, N, nVar)] =
          primitive[SC_2D_INDEX(i, n, 1, iEl, N, nVar)] *  // v
          solution[SC_2D_INDEX(i, n, iVar, iEl, N, nVar)]; // rho*u
    }
  }
  else if (iVar == 1) // rho*v
  {
    for (int n = 0; n <= N; n++)
    {
      // rho*v*u
      flux[P2PVE_2D_INDEX(1, 1, n, i, j, iVar, iEl, N, nVar)] =
          primitive[SC_2D_INDEX(n, j, 0, iEl, N, nVar)] *  // u
          solution[SC_2D_INDEX(n, j, iVar, iEl, N, nVar)]; // rho*v

      // rho*v*v + p
      flux[P2PVE_2D_INDEX(2, 1, n, i, j, iVar, iEl, N, nVar)] =
          primitive[SC_2D_INDEX(n, j, 1, iEl, N, nVar)] *       // v
              solution[SC_2D_INDEX(n, j, iVar, iEl, N, nVar)] + // rho*v
          primitive[SC_2D_INDEX(n, j, 3, iEl, N, nVar)];

      // rho*v*u
      flux[P2PVE_2D_INDEX(1, 2, n, i, j, iVar, iEl, N, nVar)] =
          primitive[SC_2D_INDEX(i, n, 0, iEl, N, nVar)] *  // u
          solution[SC_2D_INDEX(i, n, iVar, iEl, N, nVar)]; // rho*v

      // rho*v*v + p
      flux[P2PVE_2D_INDEX(2, 2, n, i, j, iVar, iEl, N, nVar)] =
          primitive[SC_2D_INDEX(i, n, 1, iEl, N, nVar)] *       // v
              solution[SC_2D_INDEX(i, n, iVar, iEl, N, nVar)] + // rho*v
          primitive[SC_2D_INDEX(i, n, 3, iEl, N, nVar)];
    }
  }
  else if (iVar == 2)
  { // density
    for (int n = 0; n <= N; n++)
    {
      flux[P2PVE_2D_INDEX(1, 1, n, i, j, iVar, iEl, N, nVar)] =
          solution[SC_2D_INDEX(n, j, 0, iEl, N, nVar)]; // rho*u

      flux[P2PVE_2D_INDEX(2, 1, n, i, j, iVar, iEl, N, nVar)] =
          solution[SC_2D_INDEX(n, j, 1, iEl, N, nVar)]; // rho*v
                                                        //
      flux[P2PVE_2D_INDEX(1, 2, n, i, j, iVar, iEl, N, nVar)] =
          solution[SC_2D_INDEX(i, n, 0, iEl, N, nVar)]; // rho*u

      flux[P2PVE_2D_INDEX(2, 2, n, i, j, iVar, iEl, N, nVar)] =
          solution[SC_2D_INDEX(i, n, 1, iEl, N, nVar)]; // rho*v
    }
  }
  else if (iVar == 3)
  { // total energy (rho*u*H)

    for (int n = 0; n <= N; n++)
    {
      // Calculate Enthalpy
      real H = (solution[SC_2D_INDEX(n, j, iVar, iEl, N, nVar)] +
                primitive[SC_2D_INDEX(n, j, 3, iEl, N, nVar)]);

      flux[P2PVE_2D_INDEX(1, 1, n, i, j, iVar, iEl, N, nVar)] =
          primitive[SC_2D_INDEX(n, j, 0, iEl, N, nVar)] * H;

      flux[P2PVE_2D_INDEX(2, 1, n, i, j, iVar, iEl, N, nVar)] =
          primitive[SC_2D_INDEX(n, j, 1, iEl, N, nVar)] * H;

      // Calculate Enthalpy
      H = (solution[SC_2D_INDEX(i, n, iVar, iEl, N, nVar)] +
           primitive[SC_2D_INDEX(i, n, 3, iEl, N, nVar)]);

      flux[P2PVE_2D_INDEX(1, 1, n, i, j, iVar, iEl, N, nVar)] =
          primitive[SC_2D_INDEX(i, n, 0, iEl, N, nVar)] * H;

      flux[P2PVE_2D_INDEX(2, 1, n, i, j, iVar, iEl, N, nVar)] =
          primitive[SC_2D_INDEX(i, n, 1, iEl, N, nVar)] * H;
    }
  }
}

extern "C"
{
  void ConservativeFlux_CompressibleIdealGas2D_gpu_wrapper(real **flux, real **solution, real **primitive, int N, int nVar, int nEl)
  {
    ConservativeFlux_CompressibleIdealGas2D_gpu<<<dim3(nVar, nEl, 1), dim3(N + 1, N + 1, 1), 0, 0>>>(*flux, *solution, *primitive, N, nVar);
  }
}

__global__ void LocalLaxFriedrichs_CompressibleIdealGas2D_gpu(real *flux,
                                                    real *solution,
                                                    real *extSolution,
                                                    real *primitive,
                                                    real *extPrimitive,
                                                    real *diagnostics,
                                                    real *extDiagnostics,
                                                    real *nHat,
                                                    real *nScale,
                                                    int N,
                                                    int nVar,
                                                    int nDiag)
{

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iSide = blockIdx.y + 1;
  size_t iEl = blockIdx.z;
  size_t i = threadIdx.x;

  // Get the boundary normals on cell edges from the mesh geometry
  real nx = nHat[VEB_2D_INDEX(1, i, 0, iSide, iEl, N, 1)];
  real ny = nHat[VEB_2D_INDEX(2, i, 0, iSide, iEl, N, 1)];

  // Calculate the normal velocity at the cell edges
  real unL = primitive[SCB_2D_INDEX(i, 0, iSide, iEl, N, nVar)] * nx +
             primitive[SCB_2D_INDEX(i, 1, iSide, iEl, N, nVar)] * ny;

  real unR = extPrimitive[SCB_2D_INDEX(i, 0, iSide, iEl, N, nVar)] * nx +
             extPrimitive[SCB_2D_INDEX(i, 1, iSide, iEl, N, nVar)] * ny;

  real fluxL = 0.0;
  real fluxR = 0.0;
  if (iVar == 0)
  {

    fluxL = unL * solution[SCB_2D_INDEX(i, iVar, iSide, iEl, N, nVar)] +
            primitive[SCB_2D_INDEX(i, 3, iSide, iEl, N, nVar)] * nx;

    fluxR = unR * extSolution[SCB_2D_INDEX(i, iVar, iSide, iEl, N, nVar)] +
            extPrimitive[SCB_2D_INDEX(i, 3, iSide, iEl, N, nDiag)] * nx;
  }
  else if (iVar == 1)
  {
    fluxL = unL * solution[SCB_2D_INDEX(i, iVar, iSide, iEl, N, nVar)] +
            primitive[SCB_2D_INDEX(i, 3, iSide, iEl, N, nVar)] * ny;

    fluxR = unR * extSolution[SCB_2D_INDEX(i, iVar, iSide, iEl, N, nVar)] +
            extPrimitive[SCB_2D_INDEX(i, 3, iSide, iEl, N, nVar)] * ny;
  }
  else if (iVar == 2)
  {
    fluxL = solution[SCB_2D_INDEX(i, 0, iSide, iEl, N, nVar)] * nx +
            solution[SCB_2D_INDEX(i, 1, iSide, iEl, N, nVar)] * ny;

    fluxR = extSolution[SCB_2D_INDEX(i, 0, iSide, iEl, N, nVar)] * nx +
            extSolution[SCB_2D_INDEX(i, 1, iSide, iEl, N, nVar)] * ny;
  }
  else if (iVar == 3)
  {
    fluxL = unL * (solution[SCB_2D_INDEX(i, 3, iSide, iEl, N, nVar)] +
                   primitive[SCB_2D_INDEX(i, 3, iSide, iEl, N, nVar)]);
    fluxR = unR * (extSolution[SCB_2D_INDEX(i, 3, iSide, iEl, N, nVar)] +
                   extPrimitive[SCB_2D_INDEX(i, 3, iSide, iEl, N, nVar)]);
  }

  real cL = diagnostics[SCB_2D_INDEX(i, 2, iSide, iEl, N, nDiag)];
  real cR = extDiagnostics[SCB_2D_INDEX(i, 2, iSide, iEl, N, nDiag)];

#ifdef DOUBLE_PRECISION
  real alpha = fmax(fabs(unL), fabs(unR)) +
               fmax(fabs(cL), fabs(cR));
#else
  real alpha = fmaxf(fabs(unL), fabs(unR)) +
               fmaxf(fabs(cL), fabs(cR));
#endif

  real jump = solution[SCB_2D_INDEX(i, iVar, iSide, iEl, N, nVar)] -
              extSolution[SCB_2D_INDEX(i, iVar, iSide, iEl, N, nVar)];
  real nmag = nScale[SCB_2D_INDEX(i, 0, iSide, iEl, N, 1)];

  // Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
  flux[SCB_2D_INDEX(i, iVar, iSide, iEl, N, nVar)] = 0.5 * (fluxL + fluxR + alpha * jump) * nmag;
}

extern "C"
{
  void LocalLaxFriedrichs_CompressibleIdealGas2D_gpu_wrapper(real **flux,
                                                   real **solution,
                                                   real **extSolution,
                                                   real **primitive,
                                                   real **extPrimitive,
                                                   real **diagnostics,
                                                   real **extDiagnostics,
                                                   real **nHat,
                                                   real **nScale,
                                                   int N,
                                                   int nVar,
                                                   int nDiag,
                                                   int nEl)
  {
    LocalLaxFriedrichs_CompressibleIdealGas2D_gpu<<<dim3(nVar, 4, nEl), dim3(N + 1, 1, 1), 0, 0>>>(*flux, *solution, *extSolution, *primitive, *extPrimitive, *diagnostics, *extDiagnostics, *nHat, *nScale, N, nVar, nDiag);
  }
}

__global__ void SetSolutionBoundaryCondition_CompressibleIdealGas2D_gpu(real *solution,
                                                                real *extSolution,
                                                                real *prescribedSolution,
                                                                real *nHat,
                                                                int *sideInfo,
                                                                int N, int nVar)
{

  // Get the array indices from the GPU thread IDs
  size_t iSide = blockIdx.x + 1;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

  int bcid = sideInfo[4 + 5 * (iSide - 1 + 4 * (iEl))];
  int e2 = sideInfo[2 + 5 * (iSide - 1 + 4 * (iEl))];

  if (e2 == 0)
  {
    if (bcid == SELF_BC_NONORMALFLOW)
    {
      real nx = nHat[VEB_2D_INDEX(1, i, 0, iSide, iEl, N, 1)];
      real ny = nHat[VEB_2D_INDEX(2, i, 0, iSide, iEl, N, 1)];
      real rhoU = solution[SCB_2D_INDEX(i, 0, iSide, iEl, N, nVar)];
      real rhoV = solution[SCB_2D_INDEX(i, 1, iSide, iEl, N, nVar)];

      extSolution[SCB_2D_INDEX(i, 0, iSide, iEl, N, nVar)] = (ny * ny - nx * nx) * rhoU - 2.0 * nx * ny * rhoV;
      extSolution[SCB_2D_INDEX(i, 1, iSide, iEl, N, nVar)] = (nx * nx - ny * ny) * rhoV - 2.0 * nx * ny * rhoU;
      extSolution[SCB_2D_INDEX(i, 2, iSide, iEl, N, nVar)] = solution[SCB_2D_INDEX(i, 2, iSide, iEl, N, nVar)];
      extSolution[SCB_2D_INDEX(i, 3, iSide, iEl, N, nVar)] = solution[SCB_2D_INDEX(i, 3, iSide, iEl, N, nVar)];
    }
    else if (bcid == SELF_BC_PRESCRIBED || bcid == SELF_BC_RADIATION)
    {
      extSolution[SCB_2D_INDEX(i, 0, iSide, iEl, N, nVar)] =
          prescribedSolution[SCB_2D_INDEX(i, 0, iSide, iEl, N, nVar)];

      extSolution[SCB_2D_INDEX(i, 1, iSide, iEl, N, nVar)] =
          prescribedSolution[SCB_2D_INDEX(i, 1, iSide, iEl, N, nVar)];

      extSolution[SCB_2D_INDEX(i, 2, iSide, iEl, N, nVar)] =
          prescribedSolution[SCB_2D_INDEX(i, 2, iSide, iEl, N, nVar)];

      extSolution[SCB_2D_INDEX(i, 3, iSide, iEl, N, nVar)] =
          prescribedSolution[SCB_2D_INDEX(i, 3, iSide, iEl, N, nVar)];
    }
  }
}

extern "C"
{
  void SetSolutionBoundaryCondition_CompressibleIdealGas2D_gpu_wrapper(real **solution,
                                                               real **extSolution,
                                                               real **prescribedSolution,
                                                               real **nHat,
                                                               int **sideInfo,
                                                               int N,
                                                               int nVar,
                                                               int nDiag,
                                                               int nEl)
  {
    SetSolutionBoundaryCondition_CompressibleIdealGas2D_gpu<<<dim3(4, nEl, 1), dim3(N + 1, 1, 1), 0, 0>>>(*solution,
                                                                                                  *extSolution,
                                                                                                  *prescribedSolution,
                                                                                                  *nHat,
                                                                                                  *sideInfo,
                                                                                                  N, nVar);
  }
}


__global__ void SetDiagBoundaryCondition_CompressibleIdealGas2D_gpu(real *diagnostics, real *extDiagnostics, real *prescribedDiagnostics, int *sideInfo, int N, int nDiag)
{

  // Get the array indices from the GPU thread IDs
  size_t iSide = blockIdx.x + 1;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t iVar = threadIdx.y;

  int bcid = sideInfo[4 + 5 * (iSide - 1 + 4 * (iEl))];
  int e2 = sideInfo[2 + 5 * (iSide - 1 + 4 * (iEl))];
  if (e2 == 0)
  {
    if (bcid == SELF_BC_NONORMALFLOW)
    {
      // Prolong the diagnostic values to the external state
      extDiagnostics[SCB_2D_INDEX(i, iVar, iSide, iEl, N, nDiag)] =
          diagnostics[SCB_2D_INDEX(i, iVar, iSide, iEl, N, nDiag)];
    }
    else if (bcid == SELF_BC_PRESCRIBED || bcid == SELF_BC_RADIATION)
    {
      // Prolong the diagnostic values to the external state
      extDiagnostics[SCB_2D_INDEX(i, iVar, iSide, iEl, N, nDiag)] =
          prescribedDiagnostics[SCB_2D_INDEX(i, iVar, iSide, iEl, N, nDiag)];
    }
  }
}

extern "C"
{
  void SetDiagBoundaryCondition_CompressibleIdealGas2D_gpu_wrapper(real **diagnostics,
                                                               real **extDiagnostics,
                                                               real **prescribedDiagnostics,
                                                               real **nHat,
                                                               int **sideInfo,
                                                               int N,
                                                               int nVar,
                                                               int nDiag,
                                                               int nEl)
  {

    SetDiagBoundaryCondition_CompressibleIdealGas2D_gpu<<<dim3(4, nEl, 1), dim3(N + 1, nDiag, 1), 0, 0>>>(*diagnostics,
                                                                                                                 *extDiagnostics,
                                                                                                                 *prescribedDiagnostics,
                                                                                                                 *sideInfo,
                                                                                                                 N, nDiag);
  }
}

__global__ void CalculateDiagnostics_CompressibleIdealGas2D_gpu(real *solution,
                                                                real *diagnostics,
                                                                real expansionFactor,
                                                                real R,
                                                                int N,
                                                                int nVar,
                                                                int nDiag)
{

  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real rhoU = solution[SC_2D_INDEX(i, j, 0, iEl, N, nVar)];
  real rhoV = solution[SC_2D_INDEX(i, j, 1, iEl, N, nVar)];
  real rho = solution[SC_2D_INDEX(i, j, 2, iEl, N, nVar)];
  real rhoE = solution[SC_2D_INDEX(i, j, 3, iEl, N, nVar)];
  real rhoKE = 0.5 * (rhoU * rhoU + rhoV * rhoV) / rho;
  real p = (expansionFactor - 1.0) * (rhoE - rhoKE);

  // Kinetic Energy
  diagnostics[SC_2D_INDEX(i, j, 0, iEl, N, nDiag)] = rhoKE;

  // Enthalpy
  diagnostics[SC_2D_INDEX(i, j, 1, iEl, N, nDiag)] = rhoE + p;

  // Speed of sound
#ifdef DOUBLE_PRECSISION
  diagnostics[SC_2D_INDEX(i, j, 2, iEl, N, nDiag)] = sqrt(expansionFactor * p / rho);
#else
  diagnostics[SC_2D_INDEX(i, j, 2, iEl, N, nDiag)] = sqrtf(expansionFactor * p / rho);
#endif

  // In-Situ Temperature
  diagnostics[SC_2D_INDEX(i, j, 3, iEl, N, nDiag)] = (2.0 / 3.0) * ((rhoE - rhoKE) / rho) / R;
}

extern "C"
{
  void CalculateDiagnostics_CompressibleIdealGas2D_gpu_wrapper(real **solution,
                                                               real **diagnostics,
                                                               real expansionFactor,
                                                               real R,
                                                               int N,
                                                               int nVar,
                                                               int nDiag,
                                                               int nEl)
  {
    CalculateDiagnostics_CompressibleIdealGas2D_gpu<<<dim3(nEl, 1, 1), dim3(N + 1, N + 1, 1), 0, 0>>>(*solution, *diagnostics, expansionFactor, R, N, nVar, nDiag);
  }
}

__global__ void ConservativeToPrimitive_CompressibleIdealGas2D_gpu(real *solution,
                                                                   real *primitive,
                                                                   real expansionFactor,
                                                                   int N,
                                                                   int nVar)
{

  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real rhoU = solution[SC_2D_INDEX(i, j, 0, iEl, N, nVar)];
  real rhoV = solution[SC_2D_INDEX(i, j, 1, iEl, N, nVar)];
  real rho = solution[SC_2D_INDEX(i, j, 2, iEl, N, nVar)];
  real rhoE = solution[SC_2D_INDEX(i, j, 3, iEl, N, nVar)];
  real rhoKE = 0.5 * (rhoU * rhoU + rhoV * rhoV) / rho;
  real p = (expansionFactor - 1.0) * (rhoE - rhoKE);

  primitive[SC_2D_INDEX(i, j, 0, iEl, N, nVar)] = rhoU / rho;
  primitive[SC_2D_INDEX(i, j, 1, iEl, N, nVar)] = rhoV / rho;
  primitive[SC_2D_INDEX(i, j, 2, iEl, N, nVar)] = rho;
  primitive[SC_2D_INDEX(i, j, 3, iEl, N, nVar)] = p;
}

extern "C"
{
  void ConservativeToPrimitive_CompressibleIdealGas2D_gpu_wrapper(real **solution,
                                                                  real **primitive,
                                                                  real expansionFactor,
                                                                  int N,
                                                                  int nVar,
                                                                  int nEl)
  {
    ConservativeToPrimitive_CompressibleIdealGas2D_gpu<<<dim3(nEl, 1, 1), dim3(N + 1, N + 1, 1), 0, 0>>>(*solution, *primitive, expansionFactor, N, nVar);
  }
}

__global__ void ConservativeToEntropy_CompressibleIdealGas2D_gpu(real *solution,
                                                                 real *entropy,
                                                                 real expansionFactor,
                                                                 int N,
                                                                 int nVar)
{

  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real rhoU = solution[SC_2D_INDEX(i, j, 0, iEl, N, nVar)];
  real rhoV = solution[SC_2D_INDEX(i, j, 1, iEl, N, nVar)];
  real rho = solution[SC_2D_INDEX(i, j, 2, iEl, N, nVar)];
  real E = solution[SC_2D_INDEX(i, j, 3, iEl, N, nVar)];
  real u = rhoU / rho;
  real v = rhoV / rho;
  real KE = 0.5 * (u * rhoU + v * rhoV);
  real p = (expansionFactor - 1.0) * (E - KE);
  real s = log(p) - expansionFactor * log(rho);

  entropy[SC_2D_INDEX(i, j, 0, iEl, N, nVar)] = u * rho / p;
  entropy[SC_2D_INDEX(i, j, 1, iEl, N, nVar)] = v * rho / p;
  entropy[SC_2D_INDEX(i, j, 2, iEl, N, nVar)] = (expansionFactor - s) / (expansionFactor - 1.0) - KE / p;
  entropy[SC_2D_INDEX(i, j, 3, iEl, N, nVar)] = -rho / p;
}

extern "C"
{
  void ConservativeToEntropy_CompressibleIdealGas2D_gpu_wrapper(real **solution,
                                                                real **entropy,
                                                                real expansionFactor,
                                                                int N,
                                                                int nVar,
                                                                int nEl)
  {
    ConservativeToEntropy_CompressibleIdealGas2D_gpu<<<dim3(nEl, 1, 1), dim3(N + 1, N + 1, 1), 0, 0>>>(*solution, *entropy, expansionFactor, N, nVar);
  }
}
