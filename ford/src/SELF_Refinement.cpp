#include "SELF_GPU_Macros.h"

// Maximum supported (N+1) for the device modal-indicator kernel. The per-thread scratch
// buffers below are sized to this bound; the Fortran GPU Init guards N so a larger degree
// fails loudly rather than overrunning. Raise this (and rebuild) to support higher degrees.
#define AMR2D_MAXNP 16

// Absolute (machine-epsilon) guard on the modal energy RATIOS below - it only keeps 0/0 out of
// r1 and r2. The amplitude gate is RELATIVE to a field-wide energy scale and is applied on the
// host in FinalizeIndicator, because the scale is a reduction (an MPI collective when the mesh
// is decomposed) that an element-local device thread cannot supply.
#ifdef DOUBLE_PRECISION
  #define AMR_ENERGY_FLOOR DBL_EPSILON
#else
  #define AMR_ENERGY_FLOOR FLT_EPSILON
#endif

// RefinementIndicator_2D_gpukernel
//
//   One device thread per element. Computes the tensor-product Legendre modal spectrum of the
//   driving field(s) on the element via two matrix passes with the precomputed nodal->modal
//   transform Pmodal (Pmodal[i + (N+1)*p] = contribution of node i to mode p), forms the total
//   and clipped modal energies, and writes the raw smoothness ratio S_e together with the gate
//   energy g = sum_v w[v]*etot(v). The mathematics mirror the portable implementation in
//   SELF_RefinementIndicator_2D_t; the log10, the amplitude gate and the refine/keep/coarsen
//   thresholds are applied host-side by the shared FinalizeIndicator, so both backends produce
//   identical flags.
__global__ void RefinementIndicator_2D_gpukernel(real *Pmodal, real *f, real *w,
                                                 real *ratio, real *gate,
                                                 int N, int nvar, int ivar, int nel){

  int iel = threadIdx.x + blockIdx.x*blockDim.x;
  if( iel < nel ){
    int Np = N+1;
    int v0, v1;
    if( ivar == 0 ){ v0 = 0; v1 = nvar-1; } // SELF_AMR_ALLVARS: reduce over all variables
    else { v0 = ivar-1; v1 = ivar-1; }      // Fortran 1-based -> C 0-based

    real tmp[AMR2D_MAXNP*AMR2D_MAXNP];
    real uhat[AMR2D_MAXNP*AMR2D_MAXNP];

    real semax = 0.0;
    real g = 0.0;
    for(int v=0; v<nvar; v++){

      // A variable is transformed if it drives the indicator, or if it carries a non-zero gate
      // weight (its energy is then needed for the gate). With the default weights these
      // coincide, so the work is exactly what it was before the gate existed.
      real wv = w[v];
      bool needSe = ( v >= v0 && v <= v1 );
      bool needG  = ( wv != 0.0 );
      if( !(needSe || needG) ) continue;

      // Pass 1 (xi): tmp[p + Np*j] = sum_i Pmodal[i + Np*p] * u(i,j)
      for(int j=0; j<Np; j++){
        for(int p=0; p<Np; p++){
          real acc = 0.0;
          for(int i=0; i<Np; i++){
            acc += Pmodal[i + Np*p]*f[SC_2D_INDEX(i,j,iel,v,N,nel)];
          }
          tmp[p + Np*j] = acc;
        }
      }
      // Pass 2 (eta): uhat[p + Np*q] = sum_j Pmodal[j + Np*q] * tmp[p + Np*j]
      for(int q=0; q<Np; q++){
        for(int p=0; p<Np; p++){
          real acc = 0.0;
          for(int jj=0; jj<Np; jj++){
            acc += Pmodal[jj + Np*q]*tmp[p + Np*jj];
          }
          uhat[p + Np*q] = acc;
        }
      }

      // Total and clipped modal energies (clip drops the highest one / two modes per direction).
      real etot = 0.0, eclip1 = 0.0, eclip2 = 0.0;
      for(int q=0; q<Np; q++){
        for(int p=0; p<Np; p++){
          real e = uhat[p + Np*q]*uhat[p + Np*q];
          etot += e;
          if( p <= Np-2 && q <= Np-2 ) eclip1 += e;
          if( p <= Np-3 && q <= Np-3 ) eclip2 += e;
        }
      }

      // Gate energy: the discrete (quadratic) entropy integral over the element when the weights
      // come from an entropy Hessian.
      if( needG ) g += wv*etot;

      if( needSe ){
        real se;
        if( etot <= (real)AMR_ENERGY_FLOOR ){
          se = 0.0;
        } else {
          real r1 = (etot - eclip1)/etot;
          real r2 = 0.0;
          if( N >= 2 && eclip1 > (real)AMR_ENERGY_FLOOR ){
            r2 = (eclip1 - eclip2)/eclip1;
          }
          se = (r1 > r2) ? r1 : r2;
        }
        if( se > semax ) semax = se;
      }
    }

    ratio[iel] = semax; // raw ratio; the host gates it, takes the log10 and sets the flag
    gate[iel]  = g;
  }
}

extern "C"
{
  void RefinementIndicator_2D_gpu(real *Pmodal, real *f, real *w, real *ratio, real *gate,
                                  int N, int nvar, int ivar, int nel)
  {
    int threads_per_block = 256;
    int nblocks_x = nel/threads_per_block + 1;
    RefinementIndicator_2D_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(
      Pmodal, f, w, ratio, gate, N, nvar, ivar, nel);
  }
}

// Maximum supported (N+1) for the 3-D device modal-indicator kernel. Unlike the 2-D kernel's
// per-thread scratch, the 3-D kernel holds two (N+1)^3 buffers in per-block __shared__ memory
// (a 16^3 per-thread buffer would be 32 KB of local memory per thread, far past any budget).
// Two 12^3 double-precision shared buffers are 2*13824 B = 27648 B, inside the 48 KB static
// shared-memory limit; the Fortran GPU Init guards N so a larger degree fails loudly rather
// than overrunning. Raise this only as far as shared memory allows (and rebuild).
#define AMR3D_MAXNP 12

// RefinementIndicator_3D_gpukernel
//
//   One device block per element, with (N+1)*(N+1) threads. Computes the tensor-product Legendre
//   modal spectrum of the driving field(s) on the element via three matrix passes with the
//   precomputed nodal->modal transform Pmodal (Pmodal[i + (N+1)*p] = contribution of node i to
//   mode p), forms the total and clipped modal energies, and writes the raw smoothness ratio S_e
//   together with the gate energy g = sum_v w[v]*etot(v). The mathematics mirror the portable
//   implementation in SELF_RefinementIndicator_3D_t; the log10, the amplitude gate and the
//   refine/keep/coarsen thresholds are applied host-side by the shared FinalizeIndicator, so
//   both backends produce identical flags.
//
//   Decomposition: thread (a,b) = (tid % (N+1), tid / (N+1)) owns one (a,b) pencil and loops the
//   third index in each pass; the two shared buffers alternate (tmp -> uhat -> tmp), exactly as
//   in the portable implementation, so the full modal coefficients land in tmp after pass 3.
//   __syncthreads() separates the passes. The O((N+1)^3) energy sums are then done by thread 0
//   alone - simple and correct, and cheap because N is small - which also keeps the summation
//   order fixed.
__global__ void RefinementIndicator_3D_gpukernel(real *Pmodal, real *f, real *w,
                                                 real *ratio, real *gate,
                                                 int N, int nvar, int ivar, int nel){

  int iel = blockIdx.x;
  int tid = threadIdx.x;
  int Np = N+1;
  int a = tid % Np; // this thread's first pencil index
  int b = tid / Np; // this thread's second pencil index
  int v0, v1;
  if( ivar == 0 ){ v0 = 0; v1 = nvar-1; } // SELF_AMR_ALLVARS: reduce over all variables
  else { v0 = ivar-1; v1 = ivar-1; }      // Fortran 1-based -> C 0-based

  __shared__ real tmp[AMR3D_MAXNP*AMR3D_MAXNP*AMR3D_MAXNP];
  __shared__ real uhat[AMR3D_MAXNP*AMR3D_MAXNP*AMR3D_MAXNP];

  // Only thread 0's accumulators are ever read; every thread carries them so the loop below
  // stays uniform across the block (the __syncthreads() calls require that).
  real semax = 0.0;
  real g = 0.0;
  for(int v=0; v<nvar; v++){

    // A variable is transformed if it drives the indicator, or if it carries a non-zero gate
    // weight (its energy is then needed for the gate). With the default weights these
    // coincide, so the work is exactly what it was before the gate existed. The condition is
    // uniform across the block, so skipping a variable skips its barriers uniformly too.
    real wv = w[v];
    bool needSe = ( v >= v0 && v <= v1 );
    bool needG  = ( wv != 0.0 );
    if( !(needSe || needG) ) continue;

    // Pass 1 (xi): tmp[p + Np*(j + Np*k)] = sum_i Pmodal[i + Np*p] * u(i,j,k), (p,j) = (a,b)
    for(int k=0; k<Np; k++){
      real acc = 0.0;
      for(int i=0; i<Np; i++){
        acc += Pmodal[i + Np*a]*f[SC_3D_INDEX(i,b,k,iel,v,N,nel)];
      }
      tmp[a + Np*(b + Np*k)] = acc;
    }
    __syncthreads();
    // Pass 2 (eta): uhat[p + Np*(q + Np*k)] = sum_j Pmodal[j + Np*q] * tmp[p + Np*(j + Np*k)],
    // (p,q) = (a,b)
    for(int k=0; k<Np; k++){
      real acc = 0.0;
      for(int j=0; j<Np; j++){
        acc += Pmodal[j + Np*b]*tmp[a + Np*(j + Np*k)];
      }
      uhat[a + Np*(b + Np*k)] = acc;
    }
    __syncthreads();
    // Pass 3 (zeta): tmp[p + Np*(q + Np*r)] = sum_k Pmodal[k + Np*r] * uhat[p + Np*(q + Np*k)],
    // (p,q) = (a,b): tmp now holds the full 3-D modal coefficients uhat(p,q,r).
    for(int r=0; r<Np; r++){
      real acc = 0.0;
      for(int k=0; k<Np; k++){
        acc += Pmodal[k + Np*r]*uhat[a + Np*(b + Np*k)];
      }
      tmp[a + Np*(b + Np*r)] = acc;
    }
    __syncthreads();

    // Total and clipped modal energies (clip drops the highest one / two modes per direction),
    // summed by thread 0 alone in a fixed order.
    if( tid == 0 ){
      real etot = 0.0, eclip1 = 0.0, eclip2 = 0.0;
      for(int r=0; r<Np; r++){
        for(int q=0; q<Np; q++){
          for(int p=0; p<Np; p++){
            real e = tmp[p + Np*(q + Np*r)]*tmp[p + Np*(q + Np*r)];
            etot += e;
            if( p <= Np-2 && q <= Np-2 && r <= Np-2 ) eclip1 += e;
            if( p <= Np-3 && q <= Np-3 && r <= Np-3 ) eclip2 += e;
          }
        }
      }

      // Gate energy: the discrete (quadratic) entropy integral over the element when the weights
      // come from an entropy Hessian.
      if( needG ) g += wv*etot;

      if( needSe ){
        real se;
        if( etot <= (real)AMR_ENERGY_FLOOR ){
          se = 0.0;
        } else {
          real r1 = (etot - eclip1)/etot;
          real r2 = 0.0;
          if( N >= 2 && eclip1 > (real)AMR_ENERGY_FLOOR ){
            r2 = (eclip1 - eclip2)/eclip1;
          }
          se = (r1 > r2) ? r1 : r2;
        }
        if( se > semax ) semax = se;
      }
    }
    // Thread 0 reads tmp above while the other threads wait here; without this barrier the next
    // variable's pass 1 would overwrite tmp under thread 0's energy sums.
    __syncthreads();
  }

  if( tid == 0 ){
    ratio[iel] = semax; // raw ratio; the host gates it, takes the log10 and sets the flag
    gate[iel]  = g;
  }
}

extern "C"
{
  void RefinementIndicator_3D_gpu(real *Pmodal, real *f, real *w, real *ratio, real *gate,
                                  int N, int nvar, int ivar, int nel)
  {
    int Np = N+1;
    RefinementIndicator_3D_gpukernel<<<dim3(nel,1,1), dim3(Np*Np,1,1), 0, 0>>>(
      Pmodal, f, w, ratio, gate, N, nvar, ivar, nel);
  }
}
