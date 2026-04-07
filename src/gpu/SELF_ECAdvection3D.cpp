#include "SELF_GPU_Macros.h"
#include <math.h>

// ============================================================
// Mirror boundary condition: extBoundary = boundary at domain faces
// 3D: 6 sides, boundary arrays indexed (i,j,side,iel,ivar)
// ============================================================

__global__ void setboundarycondition_ecadvection3d_gpukernel(
    real *extBoundary, real *boundary, int *sideInfo, int N, int nel, int nvar)
{
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*(N+1)*6*nel;

  if (idof < ndof) {
    uint32_t s1 = (idof/(N+1)/(N+1)) % 6;
    uint32_t e1 = idof/(N+1)/(N+1)/6;
    uint32_t e2 = sideInfo[INDEX3(2,s1,e1,5,6)];
    if (e2 == 0) {
      uint32_t ivar = blockIdx.y;
      extBoundary[idof + ndof*ivar] = boundary[idof + ndof*ivar];
    }
  }
}

extern "C"
{
  void setboundarycondition_ecadvection3d_gpu(
      real *extBoundary, real *boundary, int *sideInfo, int N, int nel, int nvar)
  {
    int ndof = (N+1)*(N+1)*6*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    setboundarycondition_ecadvection3d_gpukernel<<<dim3(nblocks_x,nvar,1),
      dim3(threads_per_block,1,1), 0, 0>>>(extBoundary, boundary, sideInfo, N, nel, nvar);
  }
}

// ============================================================
// LLF boundary flux for 3-D linear advection
//   flux = 0.5*(un*(sL+sR) - lambda*(sR-sL)) * nScale
//   un = u*nx + v*ny + w*nz,  lambda = sqrt(u^2+v^2+w^2)
// ============================================================

__global__ void boundaryflux_ecadvection3d_gpukernel(
    real *fb, real *fextb, real *nhat, real *nscale, real *flux,
    real u, real v, real w, real lam, int N, int nel, int nvar)
{
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*(N+1)*6*nel;

  if (idof < ndof) {
    uint32_t i   = idof % (N+1);
    uint32_t j   = (idof/(N+1)) % (N+1);
    uint32_t s1  = (idof/(N+1)/(N+1)) % 6;
    uint32_t iel = idof/(N+1)/(N+1)/6;
    uint32_t ivar = blockIdx.y;

    real nx   = nhat[VEB_3D_INDEX(i,j,s1,iel,0,0,N,nel,1)];
    real ny   = nhat[VEB_3D_INDEX(i,j,s1,iel,0,1,N,nel,1)];
    real nz   = nhat[VEB_3D_INDEX(i,j,s1,iel,0,2,N,nel,1)];
    real un   = u*nx + v*ny + w*nz;
    real nmag = nscale[SCB_3D_INDEX(i,j,s1,iel,0,N,nel)];

    real sL = fb[idof + ivar*ndof];
    real sR = fextb[idof + ivar*ndof];

    flux[idof + ivar*ndof] = 0.5*(un*(sL+sR) - lam*(sR-sL)) * nmag;
  }
}

extern "C"
{
  void boundaryflux_ecadvection3d_gpu(
      real *fb, real *fextb, real *nhat, real *nscale, real *flux,
      real u, real v, real w, int N, int nel, int nvar)
  {
    real lam = sqrt(u*u + v*v + w*w);
    int ndof = (N+1)*(N+1)*6*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    boundaryflux_ecadvection3d_gpukernel<<<dim3(nblocks_x,nvar,1),
      dim3(threads_per_block,1,1), 0, 0>>>(fb, fextb, nhat, nscale, flux, u, v, w, lam, N, nel, nvar);
  }
}

// ============================================================
// TwoPointFluxMethod for EC linear advection on a curvilinear 3-D mesh
//
// Computes contravariant two-point fluxes:
//   Fc^r(nn,i,j,k) = avg(Ja^r . a) * (s(node_L) + s(node_R)) / 2
//
// f    : output, TPV_3D_INDEX layout
// s    : solution interior, SC_3D_INDEX layout
// dsdx : metric terms, TE_3D_INDEX layout (nVar=1)
//        dsdx[iq + nq*(iel + nel*(d + 3*r))]  for physical dir d, comp dir r
// ============================================================

template <int blockSize, int matSize>
__global__ void __launch_bounds__(512) twopointfluxmethod_ecadvection3d_gpukernel(
    real *f, real *s, real *dsdx, real u, real v, real w, int nq, int N, int nvar)
{
  uint32_t iq = threadIdx.x;
  if (iq < nq) {
    uint32_t iel  = blockIdx.x;
    uint32_t nel  = gridDim.x;
    uint32_t ivar = blockIdx.y;
    uint32_t ii   = iq % (N+1);
    uint32_t jj   = (iq/(N+1)) % (N+1);
    uint32_t kk   = iq/(N+1)/(N+1);

    real s_ijk = s[iq + nq*(iel + nel*ivar)];
    int nq4 = nq*(N+1);  // stride between idir slices in TPV_3D layout

    for (int nn = 0; nn < N+1; nn++) {
      // xi^1: pair (i,j,k)-(nn,j,k)
      uint32_t iq_nn = nn + (N+1)*(jj + (N+1)*kk);
      real savg1 = 0.5*(s_ijk + s[iq_nn + nq*(iel + nel*ivar)]);
      real uc1 = 0.5*(dsdx[iq + nq*(iel + nel*0)] + dsdx[iq_nn + nq*(iel + nel*0)])*u
               + 0.5*(dsdx[iq + nq*(iel + nel*1)] + dsdx[iq_nn + nq*(iel + nel*1)])*v
               + 0.5*(dsdx[iq + nq*(iel + nel*2)] + dsdx[iq_nn + nq*(iel + nel*2)])*w;
      f[nn + (N+1)*iq + nq4*(iel + nel*ivar)] = uc1 * savg1;

      // xi^2: pair (i,j,k)-(i,nn,k)
      uint32_t iq_nn2 = ii + (N+1)*(nn + (N+1)*kk);
      real savg2 = 0.5*(s_ijk + s[iq_nn2 + nq*(iel + nel*ivar)]);
      real uc2 = 0.5*(dsdx[iq + nq*(iel + nel*3)] + dsdx[iq_nn2 + nq*(iel + nel*3)])*u
               + 0.5*(dsdx[iq + nq*(iel + nel*4)] + dsdx[iq_nn2 + nq*(iel + nel*4)])*v
               + 0.5*(dsdx[iq + nq*(iel + nel*5)] + dsdx[iq_nn2 + nq*(iel + nel*5)])*w;
      f[nn + (N+1)*iq + nq4*(iel + nel*(ivar + nvar))] = uc2 * savg2;

      // xi^3: pair (i,j,k)-(i,j,nn)
      uint32_t iq_nn3 = ii + (N+1)*(jj + (N+1)*nn);
      real savg3 = 0.5*(s_ijk + s[iq_nn3 + nq*(iel + nel*ivar)]);
      real uc3 = 0.5*(dsdx[iq + nq*(iel + nel*6)] + dsdx[iq_nn3 + nq*(iel + nel*6)])*u
               + 0.5*(dsdx[iq + nq*(iel + nel*7)] + dsdx[iq_nn3 + nq*(iel + nel*7)])*v
               + 0.5*(dsdx[iq + nq*(iel + nel*8)] + dsdx[iq_nn3 + nq*(iel + nel*8)])*w;
      f[nn + (N+1)*iq + nq4*(iel + nel*(ivar + 2*nvar))] = uc3 * savg3;
    }
  }
}

extern "C"
{
  void twopointfluxmethod_ecadvection3d_gpu(
      real *f, real *s, real *dsdx, real u, real v, real w, int N, int nvar, int nel)
  {
    int nq = (N+1)*(N+1)*(N+1);
    if (N < 4) {
      twopointfluxmethod_ecadvection3d_gpukernel<64,16><<<dim3(nel,nvar,1),
        dim3(64,1,1), 0, 0>>>(f, s, dsdx, u, v, w, nq, N, nvar);
    } else if (N >= 4 && N < 8) {
      twopointfluxmethod_ecadvection3d_gpukernel<512,64><<<dim3(nel,nvar,1),
        dim3(512,1,1), 0, 0>>>(f, s, dsdx, u, v, w, nq, N, nvar);
    }
  }
}
