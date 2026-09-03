#include "SELF_GPU_Macros.h"

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// //
//
// SELF_MatrixMultiply.cpp
//
// Hand-written CUDA/HIP tensor-product contraction kernels that replace the
// cuBLAS/hipBLAS (Xgemm / XgemvStridedBatched) matrix operators previously
// used for spectral interpolation and differentiation.
//
// All operators act on nodal, column-major (Fortran) data. The operator
// matrix A is stored column-major with leading dimension (Nc+1), so that
//     A[k + (Nc+1)*a]
// maps control node k -> target node a. This matches the storage used for
// interp (iMatrix), strong-derivative (dMatrix), DG-derivative (dgMatrix) and
// boundary (bMatrix) matrices constructed in SELF_Lagrange_t.f90.
//
// Every kernel sums its contractions in an explicit, fixed loop order that
// matches the CPU reference algorithms in the *_t.f90 templates. Results may
// therefore differ from the previous BLAS path at the round-off level, but are
// consistent with the CPU backend.
//
// Kernels are grid-strided over the quadrature/target points so that they
// support arbitrary polynomial degree N (control) and M (target) without the
// fixed-thread-count limitations present in some earlier kernels.
//
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// //

#define SELF_MATOP_THREADS 256

// ---------------------------------------------------------------------------
// 1-D operator application:  Af = A^T f
//
//   Af[a,c] = sum_{k=0}^{opAcols-1} A[k + opAcols*a] * f[k + opAcols*c]
//
// with a in [0,opArows), c in [0,ncol). One thread-block handles one column c
// (a single (element,variable) pair); the operator matrix is staged in shared
// memory and reused across all output rows of that column. Covers 1-D grid
// interpolation (A=iMatrix, opArows=M+1), strong/DG derivative (A=dMatrix or
// dgMatrix, opArows=N+1) and boundary interpolation (A=bMatrix, opArows=2).
// ---------------------------------------------------------------------------
__global__ void MatrixOp_1D_gpukernel(real *A, real *f, real *Af, int opArows, int opAcols, int ncol){

  int c = blockIdx.x; // column index (element*variable)
  if( c >= ncol ) return;

  extern __shared__ real sA[]; // opAcols*opArows entries
  for(int idx = threadIdx.x; idx < opAcols*opArows; idx += blockDim.x){
    sA[idx] = A[idx];
  }
  __syncthreads();

  for(int a = threadIdx.x; a < opArows; a += blockDim.x){
    real acc = 0.0;
    for(int k = 0; k < opAcols; k++){
      acc += sA[k + opAcols*a]*f[k + opAcols*c];
    }
    Af[a + opArows*c] = acc;
  }
}

extern "C"
{
  void MatrixOp_1D_gpu(real *A, real *f, real *Af, int opArows, int opAcols, int ncol){
    int nthreads = SELF_MATOP_THREADS;
    size_t smem = (size_t)opArows*opAcols*sizeof(real);
    MatrixOp_1D_gpukernel<<<dim3(ncol,1,1), dim3(nthreads,1,1), smem, 0>>>(A,f,Af,opArows,opAcols,ncol);
  }
}

// ---------------------------------------------------------------------------
// 2-D grid interpolation (tensor product), control degree N -> target degree M
//
//   fInterp[a,b] = sum_{i,j} A[i+(N+1)*a] * A[j+(N+1)*b] * f[i,j]
//
// One block per (element,variable); the interpolation matrix is staged in
// shared memory. Threads grid-stride over the (M+1)^2 target points.
// ---------------------------------------------------------------------------
__global__ void GridInterp_2D_gpukernel(real *A, real *f, real *fInterp, int N, int M, int nel, int nvar){

  int iel = blockIdx.x;
  int ivar = blockIdx.y;
  int nq_in = (N+1)*(N+1);
  int nq_out = (M+1)*(M+1);
  int fbase = nq_in*(iel + nel*ivar);
  int gbase = nq_out*(iel + nel*ivar);

  extern __shared__ real sA[]; // (N+1)*(M+1)
  for(int idx = threadIdx.x; idx < (N+1)*(M+1); idx += blockDim.x){
    sA[idx] = A[idx];
  }
  __syncthreads();

  for(int iq = threadIdx.x; iq < nq_out; iq += blockDim.x){
    int a = iq % (M+1);
    int b = iq/(M+1);
    real acc = 0.0;
    for(int i = 0; i < N+1; i++){
      real tmp = 0.0;
      for(int j = 0; j < N+1; j++){
        tmp += sA[j + (N+1)*b]*f[i + (N+1)*j + fbase];
      }
      acc += sA[i + (N+1)*a]*tmp;
    }
    fInterp[a + (M+1)*b + gbase] = acc;
  }
}

extern "C"
{
  void GridInterp_2D_gpu(real *A, real *f, real *fInterp, int N, int M, int nvar, int nel){
    size_t smem = (size_t)(N+1)*(M+1)*sizeof(real);
    GridInterp_2D_gpukernel<<<dim3(nel,nvar,1), dim3(SELF_MATOP_THREADS,1,1), smem, 0>>>(A,f,fInterp,N,M,nel,nvar);
  }
}

// ---------------------------------------------------------------------------
// 3-D grid interpolation (tensor product), control degree N -> target degree M
//
//   fInterp[a,b,c] = sum_{i,j,k} A[i+(N+1)*a]*A[j+(N+1)*b]*A[k+(N+1)*c]*f[i,j,k]
// ---------------------------------------------------------------------------
__global__ void GridInterp_3D_gpukernel(real *A, real *f, real *fInterp, int N, int M, int nel, int nvar){

  int iel = blockIdx.x;
  int ivar = blockIdx.y;
  int nq_in = (N+1)*(N+1)*(N+1);
  int nq_out = (M+1)*(M+1)*(M+1);
  int fbase = nq_in*(iel + nel*ivar);
  int gbase = nq_out*(iel + nel*ivar);

  extern __shared__ real sA[]; // (N+1)*(M+1)
  for(int idx = threadIdx.x; idx < (N+1)*(M+1); idx += blockDim.x){
    sA[idx] = A[idx];
  }
  __syncthreads();

  for(int iq = threadIdx.x; iq < nq_out; iq += blockDim.x){
    int a = iq % (M+1);
    int b = (iq/(M+1)) % (M+1);
    int c = iq/(M+1)/(M+1);
    real acc = 0.0;
    for(int k = 0; k < N+1; k++){
      real tk = 0.0;
      for(int j = 0; j < N+1; j++){
        real tj = 0.0;
        for(int i = 0; i < N+1; i++){
          tj += sA[i + (N+1)*a]*f[i + (N+1)*(j + (N+1)*k) + fbase];
        }
        tk += sA[j + (N+1)*b]*tj;
      }
      acc += sA[k + (N+1)*c]*tk;
    }
    fInterp[a + (M+1)*(b + (M+1)*c) + gbase] = acc;
  }
}

extern "C"
{
  void GridInterp_3D_gpu(real *A, real *f, real *fInterp, int N, int M, int nvar, int nel){
    size_t smem = (size_t)(N+1)*(M+1)*sizeof(real);
    GridInterp_3D_gpukernel<<<dim3(nel,nvar,1), dim3(SELF_MATOP_THREADS,1,1), smem, 0>>>(A,f,fInterp,N,M,nel,nvar);
  }
}

// ---------------------------------------------------------------------------
// 2-D scalar gradient (strong form), control degree N -> N.
//
//   df[i,j,0] = sum_a A[a+(N+1)*i] * f[a,j]   (d/dxi^1)
//   df[i,j,1] = sum_a A[a+(N+1)*j] * f[i,a]   (d/dxi^2)
//
// The two directional derivatives are written to separate direction slots
// (VE_2D layout). Used for scalar gradients (A=dMatrix) and, with an effective
// variable count of 2*nvar, vector gradients.
// ---------------------------------------------------------------------------
__global__ void ScalarGradient_2D_gpukernel(real *A, real *f, real *df, int N, int nel, int nvar){

  int iel = blockIdx.x;
  int ivar = blockIdx.y;
  int nq = (N+1)*(N+1);
  int base = nq*(iel + nel*ivar);

  extern __shared__ real sA[]; // (N+1)*(N+1)
  for(int idx = threadIdx.x; idx < (N+1)*(N+1); idx += blockDim.x){
    sA[idx] = A[idx];
  }
  __syncthreads();

  for(int iq = threadIdx.x; iq < nq; iq += blockDim.x){
    int i = iq % (N+1);
    int j = iq/(N+1);
    real d1 = 0.0;
    real d2 = 0.0;
    for(int a = 0; a < N+1; a++){
      d1 += sA[a + (N+1)*i]*f[a + (N+1)*j + base];
      d2 += sA[a + (N+1)*j]*f[i + (N+1)*a + base];
    }
    // VE_2D_INDEX(i,j,iel,ivar,idir,N,nel,nvar)
    df[VE_2D_INDEX(i,j,iel,ivar,0,N,nel,nvar)] = d1;
    df[VE_2D_INDEX(i,j,iel,ivar,1,N,nel,nvar)] = d2;
  }
}

extern "C"
{
  void ScalarGradient_2D_gpu(real *A, real *f, real *df, int N, int nvar, int nel){
    size_t smem = (size_t)(N+1)*(N+1)*sizeof(real);
    ScalarGradient_2D_gpukernel<<<dim3(nel,nvar,1), dim3(SELF_MATOP_THREADS,1,1), smem, 0>>>(A,f,df,N,nel,nvar);
  }
}

// ---------------------------------------------------------------------------
// 3-D scalar gradient (strong form), control degree N -> N.
//
//   df[i,j,k,0] = sum_a A[a+(N+1)*i]*f[a,j,k]
//   df[i,j,k,1] = sum_a A[a+(N+1)*j]*f[i,a,k]
//   df[i,j,k,2] = sum_a A[a+(N+1)*k]*f[i,j,a]
// ---------------------------------------------------------------------------
__global__ void ScalarGradient_3D_gpukernel(real *A, real *f, real *df, int N, int nel, int nvar){

  int iel = blockIdx.x;
  int ivar = blockIdx.y;
  int nq = (N+1)*(N+1)*(N+1);
  int base = nq*(iel + nel*ivar);

  extern __shared__ real sA[]; // (N+1)*(N+1)
  for(int idx = threadIdx.x; idx < (N+1)*(N+1); idx += blockDim.x){
    sA[idx] = A[idx];
  }
  __syncthreads();

  for(int iq = threadIdx.x; iq < nq; iq += blockDim.x){
    int i = iq % (N+1);
    int j = (iq/(N+1)) % (N+1);
    int k = iq/(N+1)/(N+1);
    real d1 = 0.0;
    real d2 = 0.0;
    real d3 = 0.0;
    for(int a = 0; a < N+1; a++){
      d1 += sA[a + (N+1)*i]*f[a + (N+1)*(j + (N+1)*k) + base];
      d2 += sA[a + (N+1)*j]*f[i + (N+1)*(a + (N+1)*k) + base];
      d3 += sA[a + (N+1)*k]*f[i + (N+1)*(j + (N+1)*a) + base];
    }
    df[VE_3D_INDEX(i,j,k,iel,ivar,0,N,nel,nvar)] = d1;
    df[VE_3D_INDEX(i,j,k,iel,ivar,1,N,nel,nvar)] = d2;
    df[VE_3D_INDEX(i,j,k,iel,ivar,2,N,nel,nvar)] = d3;
  }
}

extern "C"
{
  void ScalarGradient_3D_gpu(real *A, real *f, real *df, int N, int nvar, int nel){
    size_t smem = (size_t)(N+1)*(N+1)*sizeof(real);
    ScalarGradient_3D_gpukernel<<<dim3(nel,nvar,1), dim3(SELF_MATOP_THREADS,1,1), smem, 0>>>(A,f,df,N,nel,nvar);
  }
}

// ---------------------------------------------------------------------------
// 2-D vector divergence (strong or DG form), control degree N -> N.
//
//   df[i,j] = sum_a A[a+(N+1)*i]*f[a,j,dir=0]
//           + sum_a A[a+(N+1)*j]*f[i,a,dir=1]
//
// f is a vector field laid out with a trailing direction index (VE_2D); the
// scalar result df uses SC_2D layout. Used for the mapped (DG-)gradient, where
// f is the contravariant-weighted field jas. A=dMatrix (strong) or dgMatrix.
// ---------------------------------------------------------------------------
__global__ void VectorDivergence_2D_gpukernel(real *A, real *f, real *df, int N, int nel, int nvar){

  int iel = blockIdx.x;
  int ivar = blockIdx.y;
  int nq = (N+1)*(N+1);

  extern __shared__ real sA[]; // (N+1)*(N+1)
  for(int idx = threadIdx.x; idx < (N+1)*(N+1); idx += blockDim.x){
    sA[idx] = A[idx];
  }
  __syncthreads();

  for(int iq = threadIdx.x; iq < nq; iq += blockDim.x){
    int i = iq % (N+1);
    int j = iq/(N+1);
    real acc = 0.0;
    for(int a = 0; a < N+1; a++){
      acc += sA[a + (N+1)*i]*f[VE_2D_INDEX(a,j,iel,ivar,0,N,nel,nvar)]
           + sA[a + (N+1)*j]*f[VE_2D_INDEX(i,a,iel,ivar,1,N,nel,nvar)];
    }
    df[SC_2D_INDEX(i,j,iel,ivar,N,nel)] = acc;
  }
}

extern "C"
{
  void VectorDivergence_2D_gpu(real *A, real *f, real *df, int N, int nvar, int nel){
    size_t smem = (size_t)(N+1)*(N+1)*sizeof(real);
    VectorDivergence_2D_gpukernel<<<dim3(nel,nvar,1), dim3(SELF_MATOP_THREADS,1,1), smem, 0>>>(A,f,df,N,nel,nvar);
  }
}

// ---------------------------------------------------------------------------
// 3-D vector divergence (strong or DG form), control degree N -> N.
//
//   df[i,j,k] = sum_a A[a+(N+1)*i]*f[a,j,k,dir=0]
//             + sum_a A[a+(N+1)*j]*f[i,a,k,dir=1]
//             + sum_a A[a+(N+1)*k]*f[i,j,a,dir=2]
// ---------------------------------------------------------------------------
__global__ void VectorDivergence_3D_gpukernel(real *A, real *f, real *df, int N, int nel, int nvar){

  int iel = blockIdx.x;
  int ivar = blockIdx.y;
  int nq = (N+1)*(N+1)*(N+1);

  extern __shared__ real sA[]; // (N+1)*(N+1)
  for(int idx = threadIdx.x; idx < (N+1)*(N+1); idx += blockDim.x){
    sA[idx] = A[idx];
  }
  __syncthreads();

  for(int iq = threadIdx.x; iq < nq; iq += blockDim.x){
    int i = iq % (N+1);
    int j = (iq/(N+1)) % (N+1);
    int k = iq/(N+1)/(N+1);
    real acc = 0.0;
    for(int a = 0; a < N+1; a++){
      acc += sA[a + (N+1)*i]*f[VE_3D_INDEX(a,j,k,iel,ivar,0,N,nel,nvar)]
           + sA[a + (N+1)*j]*f[VE_3D_INDEX(i,a,k,iel,ivar,1,N,nel,nvar)]
           + sA[a + (N+1)*k]*f[VE_3D_INDEX(i,j,a,iel,ivar,2,N,nel,nvar)];
    }
    df[SC_3D_INDEX(i,j,k,iel,ivar,N,nel)] = acc;
  }
}

extern "C"
{
  void VectorDivergence_3D_gpu(real *A, real *f, real *df, int N, int nvar, int nel){
    size_t smem = (size_t)(N+1)*(N+1)*sizeof(real);
    VectorDivergence_3D_gpukernel<<<dim3(nel,nvar,1), dim3(SELF_MATOP_THREADS,1,1), smem, 0>>>(A,f,df,N,nel,nvar);
  }
}
