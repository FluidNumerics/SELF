#include "SELF_GPU_Macros.h"

// ============================================================
// TwoPointVectorDivergence_2D
//
// Computes the split-form (two-point) divergence in 2-D on the
// reference element:
//
//   df(i,j,iel,ivar) = 2 * sum_{n=0}^{N} [ D_{n,i} * F^0(n,i,j,iel,ivar)
//                                         + D_{n,j} * F^1(n,i,j,iel,ivar) ]
//
// Input  f  : two-point vector, layout TPV_2D_INDEX(n,i,j,iel,ivar,idir,...)
// Output df : scalar, layout SC_2D_INDEX / iq + nq*(iel + nel*ivar)
// ============================================================

template <int blockSize>
__global__ void __launch_bounds__(256) TwoPointVectorDivergence_2D_gpukernel(
    real *f, real *df, real *dmatrix, int nq, int N, int nvar)
{
  uint32_t iq = threadIdx.x;
  if (iq < nq) {
    uint32_t iel  = blockIdx.x;
    uint32_t nel  = gridDim.x;
    uint32_t ivar = blockIdx.y;
    uint32_t i    = iq % (N+1);
    uint32_t j    = iq / (N+1);

    // Load D-matrix into shared memory (nq entries, same as block size)
    __shared__ real dmloc[blockSize];
    dmloc[iq] = dmatrix[iq];
    __syncthreads();

    // nq3 = (N+1)^3 : stride between idir slices in the two-point array
    int nq3 = nq*(N+1);

    real dfLoc = 0.0;
    for (int nn = 0; nn < N+1; nn++) {
      // f^0(n,i,j,...) at idir=0  :  TPV_2D_INDEX reduces to
      //   nn + (N+1)*iq + nq3*(iel + nel*ivar)
      real f0 = f[nn + (N+1)*iq + nq3*(iel + nel*ivar)];
      // f^1(n,i,j,...) at idir=1  :  nq3*(iel + nel*(ivar + nvar))
      real f1 = f[nn + (N+1)*iq + nq3*(iel + nel*(ivar + nvar))];

      dfLoc += dmloc[nn + (N+1)*i]*f0 + dmloc[nn + (N+1)*j]*f1;
    }
    df[iq + nq*(iel + nel*ivar)] = 2.0*dfLoc;
  }
}

extern "C"
{
  void TwoPointVectorDivergence_2D_gpu(real *f, real *df, real *dmatrix, int N, int nvar, int nel)
  {
    int nq = (N+1)*(N+1);
    if (N <= 7) {
      TwoPointVectorDivergence_2D_gpukernel<64><<<dim3(nel,nvar,1), dim3(64,1,1), 0, 0>>>(
        f, df, dmatrix, nq, N, nvar);
    } else {
      TwoPointVectorDivergence_2D_gpukernel<256><<<dim3(nel,nvar,1), dim3(256,1,1), 0, 0>>>(
        f, df, dmatrix, nq, N, nvar);
    }
  }
}

// ============================================================
// TwoPointVectorDivergence_3D
//
// Computes the split-form divergence in 3-D on the reference element:
//
//   df(i,j,k,iel,ivar) = 2 * sum_n [ D_{n,i} F^0(n,i,j,k)
//                                   + D_{n,j} F^1(n,i,j,k)
//                                   + D_{n,k} F^2(n,i,j,k) ]
// ============================================================

template <int blockSize, int matSize>
__global__ void __launch_bounds__(512) TwoPointVectorDivergence_3D_gpukernel(
    real *f, real *df, real *dmatrix, int nq, int N, int nvar)
{
  uint32_t iq = threadIdx.x;
  if (iq < nq) {
    uint32_t iel  = blockIdx.x;
    uint32_t nel  = gridDim.x;
    uint32_t ivar = blockIdx.y;
    uint32_t i    = iq % (N+1);
    uint32_t j    = (iq/(N+1)) % (N+1);
    uint32_t k    = iq/(N+1)/(N+1);

    // Load D-matrix (matSize = (N+1)^2 entries) into shared memory.
    // Only threads with k==0 contribute; all others are idle for this load.
    __shared__ real dmloc[matSize];
    if (k == 0) {
      dmloc[i + (N+1)*j] = dmatrix[i + (N+1)*j];
    }
    __syncthreads();

    // nq4 = (N+1)^4 : stride between idir slices
    int nq4 = nq*(N+1);

    real dfLoc = 0.0;
    for (int nn = 0; nn < N+1; nn++) {
      real f0 = f[nn + (N+1)*iq + nq4*(iel + nel*ivar)];
      real f1 = f[nn + (N+1)*iq + nq4*(iel + nel*(ivar + nvar))];
      real f2 = f[nn + (N+1)*iq + nq4*(iel + nel*(ivar + 2*nvar))];

      dfLoc += dmloc[nn + (N+1)*i]*f0
             + dmloc[nn + (N+1)*j]*f1
             + dmloc[nn + (N+1)*k]*f2;
    }
    df[iq + nq*(iel + nel*ivar)] = 2.0*dfLoc;
  }
}

extern "C"
{
  void TwoPointVectorDivergence_3D_gpu(real *f, real *df, real *dmatrix, int N, int nvar, int nel)
  {
    int nq = (N+1)*(N+1)*(N+1);
    if (N < 4) {
      TwoPointVectorDivergence_3D_gpukernel<64,16><<<dim3(nel,nvar,1), dim3(64,1,1), 0, 0>>>(
        f, df, dmatrix, nq, N, nvar);
    } else if (N >= 4 && N < 8) {
      TwoPointVectorDivergence_3D_gpukernel<512,64><<<dim3(nel,nvar,1), dim3(512,1,1), 0, 0>>>(
        f, df, dmatrix, nq, N, nvar);
    }
  }
}

// ============================================================
// MappedTwoPointVectorDivergence_2D
//
// Fused kernel for the physical-space split-form divergence on a
// curvilinear 2-D mesh following Winters, Kopriva, Gassner & Hindenlang.
//
// For each node (i,j) the contravariant two-point fluxes are formed
// by projecting the physical-space two-point fluxes onto *averaged*
// metric terms:
//
//   F~^r_{(i,n),j} = sum_d (Ja^r_d(i,j) + Ja^r_d(n,j))/2 * f^d(n,i,j)
//
// where dsdx[iq + nq*(iel + nel*(d + 2*r))] = J*a^r_d  (0-based d,r).
//
// The physical divergence is:
//   df = (2/J) * sum_n [ D_{n,i} F~^1 + D_{n,j} F~^2 ]
// ============================================================

template <int blockSize>
__global__ void __launch_bounds__(256) MappedTwoPointVectorDivergence_2D_gpukernel(
    real *f, real *df, real *dmatrix, real *dsdx, real *jacobian, int nq, int N, int nvar)
{
  uint32_t iq = threadIdx.x;
  if (iq < nq) {
    uint32_t iel  = blockIdx.x;
    uint32_t nel  = gridDim.x;
    uint32_t ivar = blockIdx.y;
    uint32_t i    = iq % (N+1);
    uint32_t j    = iq / (N+1);

    __shared__ real dmloc[blockSize];
    dmloc[iq] = dmatrix[iq];
    __syncthreads();

    int nq3 = nq*(N+1);

    // Pre-load the four metric components at the current node (i,j) to
    // avoid repeated global-memory reads inside the nn loop.
    real m00 = dsdx[iq + nq*(iel + nel*0)]; // J*a^1_x  at (i,j)
    real m10 = dsdx[iq + nq*(iel + nel*1)]; // J*a^1_y  at (i,j)
    real m01 = dsdx[iq + nq*(iel + nel*2)]; // J*a^2_x  at (i,j)
    real m11 = dsdx[iq + nq*(iel + nel*3)]; // J*a^2_y  at (i,j)

    real dfLoc = 0.0;
    for (int nn = 0; nn < N+1; nn++) {
      // Node (nn,j) for the xi^1 averaged metric
      uint32_t iq_nn_j = nn + (N+1)*j;
      // Node (i,nn) for the xi^2 averaged metric
      uint32_t iq_i_nn = i + (N+1)*nn;

      // Physical-space two-point fluxes (same value used for both sums)
      real fx = f[nn + (N+1)*iq + nq3*(iel + nel*ivar)];
      real fy = f[nn + (N+1)*iq + nq3*(iel + nel*(ivar + nvar))];

      // Contravariant two-point flux in xi^1: metric averaged between (i,j)-(nn,j)
      real Fc1 = 0.5*(m00 + dsdx[iq_nn_j + nq*(iel + nel*0)])*fx
               + 0.5*(m10 + dsdx[iq_nn_j + nq*(iel + nel*1)])*fy;

      // Contravariant two-point flux in xi^2: metric averaged between (i,j)-(i,nn)
      real Fc2 = 0.5*(m01 + dsdx[iq_i_nn + nq*(iel + nel*2)])*fx
               + 0.5*(m11 + dsdx[iq_i_nn + nq*(iel + nel*3)])*fy;

      dfLoc += dmloc[nn + (N+1)*i]*Fc1 + dmloc[nn + (N+1)*j]*Fc2;
    }
    df[iq + nq*(iel + nel*ivar)] = 2.0*dfLoc/jacobian[iq + nq*iel];
  }
}

extern "C"
{
  void MappedTwoPointVectorDivergence_2D_gpu(
      real *f, real *df, real *dmatrix, real *dsdx, real *jacobian, int N, int nvar, int nel)
  {
    int nq = (N+1)*(N+1);
    if (N <= 7) {
      MappedTwoPointVectorDivergence_2D_gpukernel<64><<<dim3(nel,nvar,1), dim3(64,1,1), 0, 0>>>(
        f, df, dmatrix, dsdx, jacobian, nq, N, nvar);
    } else {
      MappedTwoPointVectorDivergence_2D_gpukernel<256><<<dim3(nel,nvar,1), dim3(256,1,1), 0, 0>>>(
        f, df, dmatrix, dsdx, jacobian, nq, N, nvar);
    }
  }
}

// ============================================================
// MappedTwoPointVectorDivergence_3D
//
// Fused kernel for the physical-space split-form divergence on a
// curvilinear 3-D mesh.
//
// dsdx[iq + nq*(iel + nel*(d + 3*r))] = J*a^r_d  (0-based d in {0,1,2},
//                                                   r in {0,1,2}).
//
// Three averaging pairs:
//   xi^1: metric averaged between (i,j,k)-(nn,j,k)
//   xi^2: metric averaged between (i,j,k)-(i,nn,k)
//   xi^3: metric averaged between (i,j,k)-(i,j,nn)
// ============================================================

template <int blockSize, int matSize>
__global__ void __launch_bounds__(512) MappedTwoPointVectorDivergence_3D_gpukernel(
    real *f, real *df, real *dmatrix, real *dsdx, real *jacobian, int nq, int N, int nvar)
{
  uint32_t iq = threadIdx.x;
  if (iq < nq) {
    uint32_t iel  = blockIdx.x;
    uint32_t nel  = gridDim.x;
    uint32_t ivar = blockIdx.y;
    uint32_t i    = iq % (N+1);
    uint32_t j    = (iq/(N+1)) % (N+1);
    uint32_t k    = iq/(N+1)/(N+1);

    __shared__ real dmloc[matSize];
    if (k == 0) {
      dmloc[i + (N+1)*j] = dmatrix[i + (N+1)*j];
    }
    __syncthreads();

    int nq4 = nq*(N+1);

    // Pre-load the nine metric components at the current node (i,j,k)
    real m00 = dsdx[iq + nq*(iel + nel*0)]; // J*a^1_x
    real m10 = dsdx[iq + nq*(iel + nel*1)]; // J*a^1_y
    real m20 = dsdx[iq + nq*(iel + nel*2)]; // J*a^1_z
    real m01 = dsdx[iq + nq*(iel + nel*3)]; // J*a^2_x
    real m11 = dsdx[iq + nq*(iel + nel*4)]; // J*a^2_y
    real m21 = dsdx[iq + nq*(iel + nel*5)]; // J*a^2_z
    real m02 = dsdx[iq + nq*(iel + nel*6)]; // J*a^3_x
    real m12 = dsdx[iq + nq*(iel + nel*7)]; // J*a^3_y
    real m22 = dsdx[iq + nq*(iel + nel*8)]; // J*a^3_z

    real dfLoc = 0.0;
    for (int nn = 0; nn < N+1; nn++) {
      // Neighbour nodes for each coordinate direction's averaging pair
      uint32_t iq_nn_jk = nn + (N+1)*(j + (N+1)*k); // (nn,j,k) for xi^1
      uint32_t iq_i_nn_k = i + (N+1)*(nn + (N+1)*k); // (i,nn,k) for xi^2
      uint32_t iq_ij_nn  = i + (N+1)*(j + (N+1)*nn); // (i,j,nn) for xi^3

      // Physical two-point fluxes
      real fx = f[nn + (N+1)*iq + nq4*(iel + nel*ivar)];
      real fy = f[nn + (N+1)*iq + nq4*(iel + nel*(ivar + nvar))];
      real fz = f[nn + (N+1)*iq + nq4*(iel + nel*(ivar + 2*nvar))];

      // Contravariant two-point flux in xi^1: (i,j,k)-(nn,j,k)
      real Fc1 = 0.5*(m00 + dsdx[iq_nn_jk + nq*(iel + nel*0)])*fx
               + 0.5*(m10 + dsdx[iq_nn_jk + nq*(iel + nel*1)])*fy
               + 0.5*(m20 + dsdx[iq_nn_jk + nq*(iel + nel*2)])*fz;

      // Contravariant two-point flux in xi^2: (i,j,k)-(i,nn,k)
      real Fc2 = 0.5*(m01 + dsdx[iq_i_nn_k + nq*(iel + nel*3)])*fx
               + 0.5*(m11 + dsdx[iq_i_nn_k + nq*(iel + nel*4)])*fy
               + 0.5*(m21 + dsdx[iq_i_nn_k + nq*(iel + nel*5)])*fz;

      // Contravariant two-point flux in xi^3: (i,j,k)-(i,j,nn)
      real Fc3 = 0.5*(m02 + dsdx[iq_ij_nn + nq*(iel + nel*6)])*fx
               + 0.5*(m12 + dsdx[iq_ij_nn + nq*(iel + nel*7)])*fy
               + 0.5*(m22 + dsdx[iq_ij_nn + nq*(iel + nel*8)])*fz;

      dfLoc += dmloc[nn + (N+1)*i]*Fc1
             + dmloc[nn + (N+1)*j]*Fc2
             + dmloc[nn + (N+1)*k]*Fc3;
    }
    df[iq + nq*(iel + nel*ivar)] = 2.0*dfLoc/jacobian[iq + nq*iel];
  }
}

extern "C"
{
  void MappedTwoPointVectorDivergence_3D_gpu(
      real *f, real *df, real *dmatrix, real *dsdx, real *jacobian, int N, int nvar, int nel)
  {
    int nq = (N+1)*(N+1)*(N+1);
    if (N < 4) {
      MappedTwoPointVectorDivergence_3D_gpukernel<64,16><<<dim3(nel,nvar,1), dim3(64,1,1), 0, 0>>>(
        f, df, dmatrix, dsdx, jacobian, nq, N, nvar);
    } else if (N >= 4 && N < 8) {
      MappedTwoPointVectorDivergence_3D_gpukernel<512,64><<<dim3(nel,nvar,1), dim3(512,1,1), 0, 0>>>(
        f, df, dmatrix, dsdx, jacobian, nq, N, nvar);
    }
  }
}
