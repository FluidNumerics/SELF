#include "SELF_GPU_Macros.h"
#include <math.h>

// ============================================================
// Mirror boundary condition: extBoundary = boundary at domain faces
// ============================================================

__global__ void setboundarycondition_ecadvection2d_gpukernel(
    real *extBoundary, real *boundary, int *sideInfo, int N, int nel, int nvar)
{
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*4*nel;

  if (idof < ndof) {
    uint32_t i  = idof % (N+1);
    uint32_t s1 = (idof/(N+1)) % 4;
    uint32_t e1 = idof/(N+1)/4;
    uint32_t e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
    if (e2 == 0) {
      uint32_t ivar = blockIdx.y;
      extBoundary[SCB_2D_INDEX(i,s1,e1,ivar,N,nel)] =
        boundary[SCB_2D_INDEX(i,s1,e1,ivar,N,nel)];
    }
  }
}

extern "C"
{
  void setboundarycondition_ecadvection2d_gpu(
      real *extBoundary, real *boundary, int *sideInfo, int N, int nel, int nvar)
  {
    int ndof = (N+1)*4*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    setboundarycondition_ecadvection2d_gpukernel<<<dim3(nblocks_x,nvar,1),
      dim3(threads_per_block,1,1), 0, 0>>>(extBoundary, boundary, sideInfo, N, nel, nvar);
  }
}

// ============================================================
// LLF boundary flux for linear advection
//   flux = 0.5*(un*(sL+sR) - lambda*(sR-sL)) * nScale
//   un = u*nx + v*ny,  lambda = sqrt(u^2+v^2)
// ============================================================

__global__ void boundaryflux_ecadvection2d_gpukernel(
    real *fb, real *fextb, real *nhat, real *nscale, real *flux,
    real u, real v, real lam, int N, int nel, int nvar)
{
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*4*nel;

  if (idof < ndof) {
    uint32_t i   = idof % (N+1);
    uint32_t j   = (idof/(N+1)) % 4;
    uint32_t iel = idof/(N+1)/4;
    uint32_t ivar = blockIdx.y;

    real nx   = nhat[VEB_2D_INDEX(i,j,iel,0,0,N,nel,1)];
    real ny   = nhat[VEB_2D_INDEX(i,j,iel,0,1,N,nel,1)];
    real un   = u*nx + v*ny;
    real nmag = nscale[SCB_2D_INDEX(i,j,iel,0,N,nel)];

    real sL = fb[idof + ivar*ndof];
    real sR = fextb[idof + ivar*ndof];

    flux[idof + ivar*ndof] = 0.5*(un*(sL+sR) - lam*(sR-sL)) * nmag;
  }
}

extern "C"
{
  void boundaryflux_ecadvection2d_gpu(
      real *fb, real *fextb, real *nhat, real *nscale, real *flux,
      real u, real v, int N, int nel, int nvar)
  {
    real lam = sqrt(u*u + v*v);
    int ndof = (N+1)*4*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    boundaryflux_ecadvection2d_gpukernel<<<dim3(nblocks_x,nvar,1),
      dim3(threads_per_block,1,1), 0, 0>>>(fb, fextb, nhat, nscale, flux, u, v, lam, N, nel, nvar);
  }
}

// ============================================================
// TwoPointFluxMethod for EC linear advection on a curvilinear mesh
//
// Computes contravariant two-point fluxes:
//   Fc^r(nn,i,j) = avg(Ja^r . a) * (s(node_L) + s(node_R)) / 2
//
// where a = (u,v) is the constant advection velocity and the metric
// averaging uses the correct partner per direction.
//
// f   : output, TPV_2D_INDEX layout
// s   : solution interior, SC_2D_INDEX layout
// dsdx: metric terms, TE_2D_INDEX layout (nVar=1 for geometry)
// ============================================================

template <int blockSize>
__global__ void __launch_bounds__(256) twopointfluxmethod_ecadvection2d_gpukernel(
    real *f, real *s, real *dsdx, real u, real v, int nq, int N, int nvar)
{
  uint32_t iq = threadIdx.x;
  if (iq < nq) {
    uint32_t iel  = blockIdx.x;
    uint32_t nel  = gridDim.x;
    uint32_t ivar = blockIdx.y;
    uint32_t i    = iq % (N+1);
    uint32_t j    = iq / (N+1);

    real s_ij = s[iq + nq*(iel + nel*ivar)];
    int nq3 = nq*(N+1);  // stride between idir slices in TPV layout

    for (int nn = 0; nn < N+1; nn++) {
      uint32_t iq_nn_j = nn + (N+1)*j;    // node (nn,j) for xi^1
      uint32_t iq_i_nn = i  + (N+1)*nn;   // node (i,nn) for xi^2

      // xi^1: pair (i,j)-(nn,j), project onto avg(Ja^1)
      real s_nn_j = s[iq_nn_j + nq*(iel + nel*ivar)];
      real savg1 = 0.5*(s_ij + s_nn_j);
      real un_contra1 = 0.5*(dsdx[iq + nq*(iel + nel*0)] + dsdx[iq_nn_j + nq*(iel + nel*0)])*u
                      + 0.5*(dsdx[iq + nq*(iel + nel*1)] + dsdx[iq_nn_j + nq*(iel + nel*1)])*v;
      f[nn + (N+1)*iq + nq3*(iel + nel*ivar)] = un_contra1 * savg1;

      // xi^2: pair (i,j)-(i,nn), project onto avg(Ja^2)
      real s_i_nn = s[iq_i_nn + nq*(iel + nel*ivar)];
      real savg2 = 0.5*(s_ij + s_i_nn);
      real un_contra2 = 0.5*(dsdx[iq + nq*(iel + nel*2)] + dsdx[iq_i_nn + nq*(iel + nel*2)])*u
                      + 0.5*(dsdx[iq + nq*(iel + nel*3)] + dsdx[iq_i_nn + nq*(iel + nel*3)])*v;
      f[nn + (N+1)*iq + nq3*(iel + nel*(ivar + nvar))] = un_contra2 * savg2;
    }
  }
}

extern "C"
{
  void twopointfluxmethod_ecadvection2d_gpu(
      real *f, real *s, real *dsdx, real u, real v, int N, int nvar, int nel)
  {
    int nq = (N+1)*(N+1);
    if (N <= 7) {
      twopointfluxmethod_ecadvection2d_gpukernel<64><<<dim3(nel,nvar,1),
        dim3(64,1,1), 0, 0>>>(f, s, dsdx, u, v, nq, N, nvar);
    } else {
      twopointfluxmethod_ecadvection2d_gpukernel<256><<<dim3(nel,nvar,1),
        dim3(256,1,1), 0, 0>>>(f, s, dsdx, u, v, nq, N, nvar);
    }
  }
}
