/*
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2026 Fluid Numerics LLC
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
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
*/

#include "SELF_GPU_Macros.h"


// ============================================================
// Riemann boundary-flux kernel (PML variant).
// Acoustic variables 0..4 use the impedance-matched solver from
// the parent LinearEuler2D kernel; auxiliary variables 5..8
// (the PML phi fields) carry zero flux and are evolved purely by
// the source term.
// ============================================================
__global__ void boundaryflux_LinearEuler2D_PML_kernel(real *fb, real *extfb, real *nhat, real *nmag, real *flux, real rho0, int ndof){
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;

  if( idof < ndof ){
    real nx  = nhat[idof];
    real ny  = nhat[idof+ndof];
    real nm  = nmag[idof];

    real cL  = fb[idof + 4*ndof];
    real cR  = extfb[idof + 4*ndof];
    real ZL  = rho0*cL;
    real ZR  = rho0*cR;

    real unL = fb[idof +     ndof]*nx + fb[idof + 2*ndof]*ny;
    real unR = extfb[idof +  ndof]*nx + extfb[idof + 2*ndof]*ny;
    real pL  = fb[idof + 3*ndof];
    real pR  = extfb[idof + 3*ndof];

    real un_star = (ZL*unL + ZR*unR + (pL - pR)) / (ZL + ZR);
    real p_star  = (ZR*pL  + ZL*pR  + ZL*ZR*(unL - unR)) / (ZL + ZR);
    real c2_avg  = 0.5*(cL*cL + cR*cR);

    flux[idof]          = (rho0*un_star) * nm;          // density
    flux[idof + ndof]   = (p_star*nx/rho0) * nm;        // u
    flux[idof + 2*ndof] = (p_star*ny/rho0) * nm;        // v
    flux[idof + 3*ndof] = (rho0*c2_avg*un_star) * nm;   // pressure
    flux[idof + 4*ndof] = 0.0;                          // sound speed
    flux[idof + 5*ndof] = 0.0;                          // phi_rho
    flux[idof + 6*ndof] = 0.0;                          // phi_u
    flux[idof + 7*ndof] = 0.0;                          // phi_v
    flux[idof + 8*ndof] = 0.0;                          // phi_p
  }
}

extern "C"
{
  void boundaryflux_LinearEuler2D_PML_gpu(real *fb, real *extfb, real *nhat, real *nmag, real *flux, real rho0, int N, int nel, int nvar){
    int threads_per_block = 256;
    uint32_t ndof = (N+1)*4*nel;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,1,1);
    dim3 nthreads(threads_per_block,1,1);

    boundaryflux_LinearEuler2D_PML_kernel<<<nblocks,nthreads>>>(fb,extfb,nhat,nmag,flux,rho0,ndof);
  }
}

// ============================================================
// Interior-flux kernel (PML variant). Acoustic variables 0..4
// match the parent LinearEuler2D kernel; auxiliary variables 5..8
// have zero flux in both directions.
// ============================================================
__global__ void fluxmethod_LinearEuler2D_PML_kernel(real *solution, real *flux, real rho0, int ndof, int nvar){
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;

  if( idof < ndof ){
    real u = solution[idof +     ndof];
    real v = solution[idof + 2*ndof];
    real p = solution[idof + 3*ndof];
    real c = solution[idof + 4*ndof];

    flux[idof + ndof*(0 + nvar*0)] = rho0*u;   // rho, x-flux
    flux[idof + ndof*(0 + nvar*1)] = rho0*v;   // rho, y-flux

    flux[idof + ndof*(1 + nvar*0)] = p/rho0;   // u, x-flux
    flux[idof + ndof*(1 + nvar*1)] = 0.0;      // u, y-flux

    flux[idof + ndof*(2 + nvar*0)] = 0.0;      // v, x-flux
    flux[idof + ndof*(2 + nvar*1)] = p/rho0;   // v, y-flux

    flux[idof + ndof*(3 + nvar*0)] = c*c*rho0*u; // p, x-flux
    flux[idof + ndof*(3 + nvar*1)] = c*c*rho0*v; // p, y-flux

    flux[idof + ndof*(4 + nvar*0)] = 0.0;      // c, x-flux
    flux[idof + ndof*(4 + nvar*1)] = 0.0;      // c, y-flux

    // PML auxiliary variables 5..8 (phi_rho, phi_u, phi_v, phi_p)
    // carry zero flux; they are evolved exclusively by the source.
    for(int ivar = 5; ivar < 9; ivar++){
      flux[idof + ndof*(ivar + nvar*0)] = 0.0;
      flux[idof + ndof*(ivar + nvar*1)] = 0.0;
    }
  }
}

extern "C"
{
  void fluxmethod_LinearEuler2D_PML_gpu(real *solution, real *flux, real rho0, int N, int nel, int nvar){
    int ndof = (N+1)*(N+1)*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    fluxmethod_LinearEuler2D_PML_kernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(solution,flux,rho0,ndof,nvar);
  }
}

// ============================================================
// Source-term kernel (Hu 2001 unsplit PML, ADE form).
//
// For acoustic variables (0..3) the PML adds:
//   source = -(sigma_x + sigma_y) q - sigma_x sigma_y phi
// For sound speed (4): source = 0 (held fixed in time).
// For auxiliaries (5..8): source = q (so dphi/dt = q).
//
// sigma_x and sigma_y are single-variable MappedScalar2D fields,
// so their device storage is exactly ndof reals each (one value
// per spatial DOF).
// ============================================================
__global__ void sourcemethod_LinearEuler2D_PML_kernel(real *source, real *solution, real *sigma_x, real *sigma_y, int ndof){
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;

  if( idof < ndof ){
    real sx = sigma_x[idof];
    real sy = sigma_y[idof];
    real rho     = solution[idof +     0*ndof];
    real u       = solution[idof +     1*ndof];
    real v       = solution[idof +     2*ndof];
    real p       = solution[idof +     3*ndof];
    // solution[idof + 4*ndof] is c (sound speed); not used here.
    real phi_rho = solution[idof +     5*ndof];
    real phi_u   = solution[idof +     6*ndof];
    real phi_v   = solution[idof +     7*ndof];
    real phi_p   = solution[idof +     8*ndof];

    real sxy = sx + sy;
    real sxsy = sx * sy;

    source[idof + 0*ndof] = -sxy*rho - sxsy*phi_rho;
    source[idof + 1*ndof] = -sxy*u   - sxsy*phi_u;
    source[idof + 2*ndof] = -sxy*v   - sxsy*phi_v;
    source[idof + 3*ndof] = -sxy*p   - sxsy*phi_p;
    source[idof + 4*ndof] = 0.0;
    source[idof + 5*ndof] = rho;
    source[idof + 6*ndof] = u;
    source[idof + 7*ndof] = v;
    source[idof + 8*ndof] = p;
  }
}

extern "C"
{
  void sourcemethod_LinearEuler2D_PML_gpu(real *source, real *solution, real *sigma_x, real *sigma_y, int N, int nel, int nvar){
    int ndof = (N+1)*(N+1)*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    sourcemethod_LinearEuler2D_PML_kernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(source,solution,sigma_x,sigma_y,ndof);
  }
}

// ============================================================
// No-normal-flow BC (PML variant).
// Variables 0..4 follow the parent LinearEuler2D no-normal-flow
// behaviour; variables 5..8 (PML auxiliaries) are zeroed in
// extBoundary for cleanliness (they have zero Riemann flux, so
// their exterior state is mathematically irrelevant).
// ============================================================
__global__ void hbc2d_nonormalflow_lineareuler2d_pml_kernel(
    real *extBoundary, real *boundary, real *nhat,
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

    real u  = boundary[SCB_2D_INDEX(i,s1,e1,1,N,nel)];
    real v  = boundary[SCB_2D_INDEX(i,s1,e1,2,N,nel)];
    real nx = nhat[VEB_2D_INDEX(i,s1,e1,0,0,N,nel,1)];
    real ny = nhat[VEB_2D_INDEX(i,s1,e1,0,1,N,nel,1)];

    extBoundary[SCB_2D_INDEX(i,s1,e1,0,N,nel)] = boundary[SCB_2D_INDEX(i,s1,e1,0,N,nel)]; // density
    extBoundary[SCB_2D_INDEX(i,s1,e1,1,N,nel)] = (ny*ny-nx*nx)*u - 2.0*nx*ny*v;            // u
    extBoundary[SCB_2D_INDEX(i,s1,e1,2,N,nel)] = (nx*nx-ny*ny)*v - 2.0*nx*ny*u;            // v
    extBoundary[SCB_2D_INDEX(i,s1,e1,3,N,nel)] = boundary[SCB_2D_INDEX(i,s1,e1,3,N,nel)]; // pressure
    extBoundary[SCB_2D_INDEX(i,s1,e1,4,N,nel)] = boundary[SCB_2D_INDEX(i,s1,e1,4,N,nel)]; // c
    extBoundary[SCB_2D_INDEX(i,s1,e1,5,N,nel)] = 0.0;                                     // phi_rho
    extBoundary[SCB_2D_INDEX(i,s1,e1,6,N,nel)] = 0.0;                                     // phi_u
    extBoundary[SCB_2D_INDEX(i,s1,e1,7,N,nel)] = 0.0;                                     // phi_v
    extBoundary[SCB_2D_INDEX(i,s1,e1,8,N,nel)] = 0.0;                                     // phi_p
  }
}

extern "C"
{
  void hbc2d_nonormalflow_lineareuler2d_pml_gpu(
      real *extBoundary, real *boundary, real *nhat,
      int *elements, int *sides,
      int nBoundaries, int N, int nel)
  {
    int threads_per_block = 256;
    int total_dofs = nBoundaries * (N+1);
    int nblocks_x = total_dofs/threads_per_block + 1;
    hbc2d_nonormalflow_lineareuler2d_pml_kernel<<<dim3(nblocks_x,1,1),
      dim3(threads_per_block,1,1), 0, 0>>>(extBoundary, boundary, nhat,
        elements, sides, nBoundaries, N, nel);
  }
}

// ============================================================
// Radiation BC (PML variant).
// Zeros acoustic perturbations; preserves c from the interior
// side so the Riemann solver sees a consistent sound speed;
// zeros PML auxiliaries.
// ============================================================
__global__ void hbc2d_radiation_lineareuler2d_pml_kernel(
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
    extBoundary[SCB_2D_INDEX(i,s1,e1,3,N,nel)] = 0.0;
    extBoundary[SCB_2D_INDEX(i,s1,e1,4,N,nel)] = boundary[SCB_2D_INDEX(i,s1,e1,4,N,nel)]; // c preserved
    extBoundary[SCB_2D_INDEX(i,s1,e1,5,N,nel)] = 0.0;
    extBoundary[SCB_2D_INDEX(i,s1,e1,6,N,nel)] = 0.0;
    extBoundary[SCB_2D_INDEX(i,s1,e1,7,N,nel)] = 0.0;
    extBoundary[SCB_2D_INDEX(i,s1,e1,8,N,nel)] = 0.0;
  }
}

extern "C"
{
  void hbc2d_radiation_lineareuler2d_pml_gpu(
      real *extBoundary, real *boundary,
      int *elements, int *sides,
      int nBoundaries, int N, int nel)
  {
    int threads_per_block = 256;
    int total_dofs = nBoundaries * (N+1);
    int nblocks_x = total_dofs/threads_per_block + 1;
    hbc2d_radiation_lineareuler2d_pml_kernel<<<dim3(nblocks_x,1,1),
      dim3(threads_per_block,1,1), 0, 0>>>(extBoundary, boundary,
        elements, sides, nBoundaries, N, nel);
  }
}
