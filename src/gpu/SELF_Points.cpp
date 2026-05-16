#include "SELF_GPU_Macros.h"

// EvalScalarPoints_2D_kernel
// One thread per (point, variable). Each thread:
//   - reads the rank-local element id from elements[p] (1-based, 0 = unlocated)
//   - reads the cached 1D Lagrange basis values for this point at xi(1), xi(2)
//   - performs the tensor-product contraction against scalar[i,j,iEl,iVar]
//   - writes values[p + iVar*nPoints]
//
// Memory layouts (column-major, mirrors Fortran):
//   elements:   int [nPoints]                            (1-based; 0 means not found)
//   lS, lT:     real [nPoints][N+1]                      (indexed as p*(N+1)+i)
//   scalar:     real [nVar][nElem][N+1][N+1]             (use SC_2D_INDEX macro)
//   values:     real [nVar][nPoints]                     (indexed as p + iVar*nPoints)

__global__ void EvalScalarPoints_2D_kernel(real *values, int *elements,
                                           real *lS, real *lT, real *scalar,
                                           int N, int nPoints, int nElem, int nVar) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int totalWork = nPoints * nVar;
  if (idx >= totalWork) return;

  int p = idx % nPoints;
  int iVar = idx / nPoints;

  int iEl = elements[p];
  if (iEl <= 0) {
    values[p + iVar * nPoints] = (real)0;
    return;
  }
  iEl -= 1; // Fortran 1-based -> C 0-based

  // Per-point basis slice (length N+1) starts at p*(N+1).
  real *lSp = &lS[p * (N + 1)];
  real *lTp = &lT[p * (N + 1)];

  real fij = (real)0;
  for (int j = 0; j <= N; j++) {
    real fi = (real)0;
    for (int i = 0; i <= N; i++) {
      fi += lSp[i] * scalar[SC_2D_INDEX(i, j, iEl, iVar, N, nElem)];
    }
    fij += lTp[j] * fi;
  }
  values[p + iVar * nPoints] = fij;
}

extern "C" {
void EvalScalarPoints_2D_gpu(real *values, int *elements,
                             real *lS, real *lT, real *scalar,
                             int N, int nPoints, int nElem, int nVar) {
  int totalWork = nPoints * nVar;
  if (totalWork <= 0) return;
  int threadsPerBlock = 256;
  int nBlocks = (totalWork + threadsPerBlock - 1) / threadsPerBlock;
  EvalScalarPoints_2D_kernel<<<dim3(nBlocks, 1, 1), dim3(threadsPerBlock, 1, 1), 0, 0>>>(
      values, elements, lS, lT, scalar, N, nPoints, nElem, nVar);
}
}

__global__ void EvalScalarPoints_3D_kernel(real *values, int *elements,
                                           real *lS, real *lT, real *lU,
                                           real *scalar,
                                           int N, int nPoints, int nElem, int nVar) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int totalWork = nPoints * nVar;
  if (idx >= totalWork) return;

  int p = idx % nPoints;
  int iVar = idx / nPoints;

  int iEl = elements[p];
  if (iEl <= 0) {
    values[p + iVar * nPoints] = (real)0;
    return;
  }
  iEl -= 1;

  real *lSp = &lS[p * (N + 1)];
  real *lTp = &lT[p * (N + 1)];
  real *lUp = &lU[p * (N + 1)];

  real fijk = (real)0;
  for (int k = 0; k <= N; k++) {
    real fij = (real)0;
    for (int j = 0; j <= N; j++) {
      real fi = (real)0;
      for (int i = 0; i <= N; i++) {
        fi += lSp[i] * scalar[SC_3D_INDEX(i, j, k, iEl, iVar, N, nElem)];
      }
      fij += lTp[j] * fi;
    }
    fijk += lUp[k] * fij;
  }
  values[p + iVar * nPoints] = fijk;
}

extern "C" {
void EvalScalarPoints_3D_gpu(real *values, int *elements,
                             real *lS, real *lT, real *lU, real *scalar,
                             int N, int nPoints, int nElem, int nVar) {
  int totalWork = nPoints * nVar;
  if (totalWork <= 0) return;
  int threadsPerBlock = 256;
  int nBlocks = (totalWork + threadsPerBlock - 1) / threadsPerBlock;
  EvalScalarPoints_3D_kernel<<<dim3(nBlocks, 1, 1), dim3(threadsPerBlock, 1, 1), 0, 0>>>(
      values, elements, lS, lT, lU, scalar, N, nPoints, nElem, nVar);
}
}
