#include "SELF_GPU_Macros.h"

// Maximum supported (N+1) for the device modal-indicator kernel. The per-thread scratch
// buffers below are sized to this bound; the Fortran GPU Init guards N so a larger degree
// fails loudly rather than overrunning. Raise this (and rebuild) to support higher degrees.
#define AMR2D_MAXNP 16

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
//   and clipped modal energies, and writes the log10 smoothness ratio and the refine/keep/coarsen
//   flag. The mathematics mirror the portable implementation in SELF_RefinementIndicator_2D_t.
__global__ void RefinementIndicator_2D_gpukernel(real *Pmodal, real *f, real *indicator, int *flag,
                                                 real refineThreshold, real coarsenThreshold,
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
    for(int v=v0; v<=v1; v++){

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

    real sigma = log10( (semax > (real)AMR_ENERGY_FLOOR) ? semax : (real)AMR_ENERGY_FLOOR );
    indicator[iel] = sigma;
    if( sigma > refineThreshold )      flag[iel] =  1; // SELF_AMR_REFINE
    else if( sigma < coarsenThreshold) flag[iel] = -1; // SELF_AMR_COARSEN
    else                               flag[iel] =  0; // SELF_AMR_KEEP
  }
}

extern "C"
{
  void RefinementIndicator_2D_gpu(real *Pmodal, real *f, real *indicator, int *flag,
                                  real refineThreshold, real coarsenThreshold,
                                  int N, int nvar, int ivar, int nel)
  {
    int threads_per_block = 256;
    int nblocks_x = nel/threads_per_block + 1;
    RefinementIndicator_2D_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(
      Pmodal, f, indicator, flag, refineThreshold, coarsenThreshold, N, nvar, ivar, nel);
  }
}
