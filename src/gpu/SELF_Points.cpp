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

// DiracDelta_2D_kernel
// One thread per (i, j, p). Each thread that owns a located point computes
//
//   J_0(p) = sum_{m,n} lS[p][m] * lT[p][n] * J[m,n,iEl-1,0]
//   scalar[i,j,iEl-1,p] = lS[p][i] * lT[p][j] / (qW[i] * qW[j] * J_0(p))
//
// for iEl = elements[p]. Threads for points with elements[p] <= 0 do nothing;
// the array is zeroed by the host launcher before the kernel runs so all other
// (i,j,iEl',p) entries (iEl' != elements[p]) end up zero.
//
// Memory layouts (column-major, mirrors Fortran):
//   scalar:    real [nVar=nPoints][nElem][N+1][N+1]    (use SC_2D_INDEX)
//   J:         real [1][nElem][N+1][N+1]               (J%interior, iVar=0)
//   lS, lT:    real [nPoints][N+1]                     (p*(N+1)+i)
//   qW:        real [N+1]                              (LGL weights)
//   elements:  int  [nPoints]                          (1-based; 0 means not found)

__global__ void DiracDelta_2D_kernel(real *scalar, real *J,
                                     real *lS, real *lT, real *qW,
                                     int *elements,
                                     int N, int nPoints, int nElem) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int Np1 = N + 1;
  int totalWork = Np1 * Np1 * nPoints;
  if (idx >= totalWork) return;

  int i = idx % Np1;
  int j = (idx / Np1) % Np1;
  int p = idx / (Np1 * Np1);

  int iEl = elements[p];
  if (iEl <= 0) return; // array was pre-zeroed by the host launcher
  iEl -= 1;

  real *lSp = &lS[p * Np1];
  real *lTp = &lT[p * Np1];

  real J0 = (real)0;
  for (int n = 0; n < Np1; n++) {
    real ln = lTp[n];
    for (int m = 0; m < Np1; m++) {
      J0 += lSp[m] * ln * J[SC_2D_INDEX(m, n, iEl, 0, N, nElem)];
    }
  }

  scalar[SC_2D_INDEX(i, j, iEl, p, N, nElem)] =
      lSp[i] * lTp[j] / (qW[i] * qW[j] * J0);
}

extern "C" {
void DiracDelta_2D_gpu(real *scalar, real *J,
                       real *lS, real *lT, real *qW,
                       int *elements,
                       int N, int nPoints, int nElem, int nVar) {
  // Zero the whole field first so points with elements==0 and all other-element
  // entries (iEl' != elements[p]) are zero on exit.
  size_t nBytes = (size_t)(N + 1) * (N + 1) * (size_t)nElem * (size_t)nVar * sizeof(real);
#ifdef __HIP_PLATFORM_AMD__
  CHECK(hipMemset(scalar, 0, nBytes));
#else
  CHECK(cudaMemset(scalar, 0, nBytes));
#endif

  int totalWork = (N + 1) * (N + 1) * nPoints;
  if (totalWork <= 0) return;
  int threadsPerBlock = 256;
  int nBlocks = (totalWork + threadsPerBlock - 1) / threadsPerBlock;
  DiracDelta_2D_kernel<<<dim3(nBlocks, 1, 1), dim3(threadsPerBlock, 1, 1), 0, 0>>>(
      scalar, J, lS, lT, qW, elements, N, nPoints, nElem);
}
}

// DiracDelta_3D_kernel: see DiracDelta_2D_kernel.

__global__ void DiracDelta_3D_kernel(real *scalar, real *J,
                                     real *lS, real *lT, real *lU, real *qW,
                                     int *elements,
                                     int N, int nPoints, int nElem) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int Np1 = N + 1;
  int totalWork = Np1 * Np1 * Np1 * nPoints;
  if (idx >= totalWork) return;

  int i = idx % Np1;
  int j = (idx / Np1) % Np1;
  int k = (idx / (Np1 * Np1)) % Np1;
  int p = idx / (Np1 * Np1 * Np1);

  int iEl = elements[p];
  if (iEl <= 0) return;
  iEl -= 1;

  real *lSp = &lS[p * Np1];
  real *lTp = &lT[p * Np1];
  real *lUp = &lU[p * Np1];

  real J0 = (real)0;
  for (int l = 0; l < Np1; l++) {
    real ll = lUp[l];
    for (int n = 0; n < Np1; n++) {
      real ln = lTp[n];
      for (int m = 0; m < Np1; m++) {
        J0 += lSp[m] * ln * ll * J[SC_3D_INDEX(m, n, l, iEl, 0, N, nElem)];
      }
    }
  }

  scalar[SC_3D_INDEX(i, j, k, iEl, p, N, nElem)] =
      lSp[i] * lTp[j] * lUp[k] / (qW[i] * qW[j] * qW[k] * J0);
}

extern "C" {
void DiracDelta_3D_gpu(real *scalar, real *J,
                       real *lS, real *lT, real *lU, real *qW,
                       int *elements,
                       int N, int nPoints, int nElem, int nVar) {
  size_t Np1 = (size_t)(N + 1);
  size_t nBytes = Np1 * Np1 * Np1 * (size_t)nElem * (size_t)nVar * sizeof(real);
#ifdef __HIP_PLATFORM_AMD__
  CHECK(hipMemset(scalar, 0, nBytes));
#else
  CHECK(cudaMemset(scalar, 0, nBytes));
#endif

  int totalWork = (N + 1) * (N + 1) * (N + 1) * nPoints;
  if (totalWork <= 0) return;
  int threadsPerBlock = 256;
  int nBlocks = (totalWork + threadsPerBlock - 1) / threadsPerBlock;
  DiracDelta_3D_kernel<<<dim3(nBlocks, 1, 1), dim3(threadsPerBlock, 1, 1), 0, 0>>>(
      scalar, J, lS, lT, lU, qW, elements, N, nPoints, nElem);
}
}
