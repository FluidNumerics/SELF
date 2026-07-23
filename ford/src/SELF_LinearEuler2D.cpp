/*
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2024 Fluid Numerics LLC
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
!    the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
*/

#include "SELF_GPU_Macros.h"


__global__ void boundaryflux_LinearEuler2D_kernel(real *fb, real *extfb, real *nhat, real *nmag, real *flux, int ndof){
  // Characteristic-decomposition (impedance-matched) interface flux for
  // linear acoustics with possibly discontinuous sound speed and background
  // density. See the CPU-side Fortran subroutine riemannflux2d_LinearEuler2D_t
  // for the derivation. This replaces LLF, which over-dissipates the tangential
  // mode and fails to stably handle the impedance mismatch
  // at high polynomial order (aliasing instability). The per-side background
  // density rho0 is read from variable 5 (index 4); the impedances use each
  // side's own rho0, and the reconstructed fluxes use the face average.
  // Variable layout (0-based): 0=u, 1=v, 2=p, 3=c, 4=rho0.
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;

  if( idof < ndof ){
    real nx  = nhat[idof];
    real ny  = nhat[idof+ndof];
    real nm  = nmag[idof];

    real cL     = fb[idof + 3*ndof];
    real cR     = extfb[idof + 3*ndof];
    real rho0L  = fb[idof + 4*ndof];
    real rho0R  = extfb[idof + 4*ndof];
    real rho0_avg = 0.5*(rho0L + rho0R);
    real ZL  = rho0L*cL;
    real ZR  = rho0R*cR;

    real unL = fb[idof]*nx + fb[idof + ndof]*ny;
    real unR = extfb[idof]*nx + extfb[idof + ndof]*ny;
    real pL  = fb[idof + 2*ndof];
    real pR  = extfb[idof + 2*ndof];

    real un_star = (ZL*unL + ZR*unR + (pL - pR)) / (ZL + ZR);
    real p_star  = (ZR*pL  + ZL*pR  + ZL*ZR*(unL - unR)) / (ZL + ZR);
    real c2_avg  = 0.5*(cL*cL + cR*cR);

    flux[idof]          = (p_star*nx/rho0_avg) * nm;        // u
    flux[idof + ndof]   = (p_star*ny/rho0_avg) * nm;        // v
    flux[idof + 2*ndof] = (rho0_avg*c2_avg*un_star) * nm;   // pressure
    flux[idof + 3*ndof] = 0.0;                              // sound speed
    flux[idof + 4*ndof] = 0.0;                              // background density
  }
}

extern "C"
{
  void boundaryflux_LinearEuler2D_gpu(real *fb, real *extfb,real *nhat, real *nmag, real *flux, int N, int nel, int nvar){
    int threads_per_block = 256;
    uint32_t ndof = (N+1)*4*nel;
    int nblocks_x = ndof/threads_per_block +1;

    dim3 nblocks(nblocks_x,nvar,1);
    dim3 nthreads(threads_per_block,1,1);

    boundaryflux_LinearEuler2D_kernel<<<nblocks,nthreads>>>(fb,extfb,nhat,nmag,flux,ndof);
  }
}

  __global__ void fluxmethod_LinearEuler2D_gpukernel(real *solution, real *flux, int ndof, int nvar){
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;

  if( idof < ndof ){
    real u = solution[idof];
    real v = solution[idof + ndof];
    real p = solution[idof + 2*ndof];
    real c = solution[idof + 3*ndof];
    real rho0 = solution[idof + 4*ndof];

    flux[idof + ndof*(0 + nvar*0)] = p/rho0; // x-velocity, x flux; p/rho0
    flux[idof + ndof*(0 + nvar*1)] = 0.0; // x-velocity, y flux; 0

    flux[idof + ndof*(1 + nvar*0)] = 0.0; // y-velocity, x flux; 0
    flux[idof + ndof*(1 + nvar*1)] = p/rho0; // y-velocity, y flux; p/rho0

    flux[idof + ndof*(2 + nvar*0)] = c*c*rho0*u; // pressure, x flux : rho0*c^2*u
    flux[idof + ndof*(2 + nvar*1)] = c*c*rho0*v; // pressure, y flux : rho0*c^2*v

    flux[idof + ndof*(3 + nvar*0)] = 0.0; // sound speed, x flux; 0 (c held fixed in time)
    flux[idof + ndof*(3 + nvar*1)] = 0.0; // sound speed, y flux; 0 (c held fixed in time)

    flux[idof + ndof*(4 + nvar*0)] = 0.0; // background density, x flux; 0 (rho0 held fixed in time)
    flux[idof + ndof*(4 + nvar*1)] = 0.0; // background density, y flux; 0 (rho0 held fixed in time)
  }

}
extern "C"
{
  void fluxmethod_LinearEuler2D_gpu(real *solution, real *flux, int N, int nel, int nvar){
    int ndof = (N+1)*(N+1)*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
    fluxmethod_LinearEuler2D_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(solution,flux,ndof,nvar);
  }

}
// ============================================================
// No-normal-flow BC kernel for 2D Linear Euler
// Operates on pre-filtered boundary faces via elements/sides arrays
// ============================================================
__global__ void hbc2d_nonormalflow_lineareuler2d_kernel(
    real *extBoundary, real *boundary, real *nhat,
    int *elements, int *sides,
    int nBoundaries, int N, int nel)
{
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t total_dofs = nBoundaries * (N+1);

  if(idof < total_dofs){
    uint32_t i  = idof % (N+1);
    uint32_t n  = idof / (N+1);
    uint32_t e1 = elements[n] - 1; // Fortran 1-based to C 0-based
    uint32_t s1 = sides[n] - 1;

    real u  = boundary[SCB_2D_INDEX(i,s1,e1,0,N,nel)];
    real v  = boundary[SCB_2D_INDEX(i,s1,e1,1,N,nel)];
    real nx = nhat[VEB_2D_INDEX(i,s1,e1,0,0,N,nel,1)];
    real ny = nhat[VEB_2D_INDEX(i,s1,e1,0,1,N,nel,1)];

    extBoundary[SCB_2D_INDEX(i,s1,e1,0,N,nel)] = (ny*ny-nx*nx)*u - 2.0*nx*ny*v; // u
    extBoundary[SCB_2D_INDEX(i,s1,e1,1,N,nel)] = (nx*nx-ny*ny)*v - 2.0*nx*ny*u; // v
    extBoundary[SCB_2D_INDEX(i,s1,e1,2,N,nel)] = boundary[SCB_2D_INDEX(i,s1,e1,2,N,nel)]; // pressure
    extBoundary[SCB_2D_INDEX(i,s1,e1,3,N,nel)] = boundary[SCB_2D_INDEX(i,s1,e1,3,N,nel)]; // c
    extBoundary[SCB_2D_INDEX(i,s1,e1,4,N,nel)] = boundary[SCB_2D_INDEX(i,s1,e1,4,N,nel)]; // rho0
  }
}

extern "C"
{
  void hbc2d_nonormalflow_lineareuler2d_gpu(
      real *extBoundary, real *boundary, real *nhat,
      int *elements, int *sides,
      int nBoundaries, int N, int nel)
  {
    int threads_per_block = 256;
    int total_dofs = nBoundaries * (N+1);
    int nblocks_x = total_dofs/threads_per_block + 1;
    hbc2d_nonormalflow_lineareuler2d_kernel<<<dim3(nblocks_x,1,1),
      dim3(threads_per_block,1,1), 0, 0>>>(extBoundary, boundary, nhat,
        elements, sides, nBoundaries, N, nel);
  }
}

// ============================================================
// Radiation BC kernel for 2D Linear Euler
// Sets u/v/p extBoundary = 0 on pre-filtered boundary
// faces. The sound speed (index 3) and background density (index 4)
// are copied from the interior side so that face Riemann fluxes see
// a consistent c and rho0.
// ============================================================
__global__ void hbc2d_radiation_lineareuler2d_kernel(
    real *extBoundary, real *boundary,
    int *elements, int *sides,
    int nBoundaries, int N, int nel)
{
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t total_dofs = nBoundaries * (N+1);

  if(idof < total_dofs){
    uint32_t i  = idof % (N+1);
    uint32_t n  = idof / (N+1);
    uint32_t e1 = elements[n] - 1;
    uint32_t s1 = sides[n] - 1;

    extBoundary[SCB_2D_INDEX(i,s1,e1,0,N,nel)] = 0.0;
    extBoundary[SCB_2D_INDEX(i,s1,e1,1,N,nel)] = 0.0;
    extBoundary[SCB_2D_INDEX(i,s1,e1,2,N,nel)] = 0.0;
    extBoundary[SCB_2D_INDEX(i,s1,e1,3,N,nel)] = boundary[SCB_2D_INDEX(i,s1,e1,3,N,nel)]; // c preserved
    extBoundary[SCB_2D_INDEX(i,s1,e1,4,N,nel)] = boundary[SCB_2D_INDEX(i,s1,e1,4,N,nel)]; // rho0 preserved
  }
}

extern "C"
{
  void hbc2d_radiation_lineareuler2d_gpu(
      real *extBoundary, real *boundary,
      int *elements, int *sides,
      int nBoundaries, int N, int nel)
  {
    int threads_per_block = 256;
    int total_dofs = nBoundaries * (N+1);
    int nblocks_x = total_dofs/threads_per_block + 1;
    hbc2d_radiation_lineareuler2d_kernel<<<dim3(nblocks_x,1,1),
      dim3(threads_per_block,1,1), 0, 0>>>(extBoundary, boundary,
        elements, sides, nBoundaries, N, nel);
  }
}
