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


// Impedance-matched (characteristic/Godunov) Riemann flux with a per-node
// sound speed carried as solution variable 6 (0-based index 5), mirroring the
// 2-D model. The interface states are resolved with the acoustic impedances
// Z = rho0*c on either side:
//   un* = (ZL*unL + ZR*unR + (pL - pR)) / (ZL + ZR)
//   p*  = (ZR*pL + ZL*pR + ZL*ZR*(unL - unR)) / (ZL + ZR)
// Variable layout (0-based): 0=rho, 1=u, 2=v, 3=w, 4=p, 5=c. The sound-speed
// variable carries zero flux.
__global__ void boundaryflux_LinearEuler3D_kernel(real *fb, real *extfb, real *nhat, real *nmag, real *flux, real rho0, int ndof){
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;

  if( idof < ndof ){

    real nx = nhat[idof];
    real ny = nhat[idof+ndof];
    real nz = nhat[idof+2*ndof];
    real nm = nmag[idof];

    real cL = fb[idof + 5*ndof];
    real cR = extfb[idof + 5*ndof];
    real ZL = rho0*cL;
    real ZR = rho0*cR;

    real unL = fb[idof + ndof]*nx + fb[idof + 2*ndof]*ny + fb[idof + 3*ndof]*nz;
    real unR = extfb[idof + ndof]*nx + extfb[idof + 2*ndof]*ny + extfb[idof + 3*ndof]*nz;
    real pL = fb[idof + 4*ndof];
    real pR = extfb[idof + 4*ndof];

    real un_star = (ZL*unL + ZR*unR + (pL - pR))/(ZL + ZR);
    real p_star = (ZR*pL + ZL*pR + ZL*ZR*(unL - unR))/(ZL + ZR);
    real c2_avg = 0.5*(cL*cL + cR*cR);

    flux[idof] = (rho0*un_star)*nm;              // density
    flux[idof+ndof] = (p_star*nx/rho0)*nm;       // u
    flux[idof+2*ndof] = (p_star*ny/rho0)*nm;     // v
    flux[idof+3*ndof] = (p_star*nz/rho0)*nm;     // w
    flux[idof+4*ndof] = (rho0*c2_avg*un_star)*nm; // pressure
    flux[idof+5*ndof] = 0.0;                     // sound speed (static)
  }
}

extern "C"
{
  void boundaryflux_LinearEuler3D_gpu(real *fb, real *extfb,real *nhat, real *nmag, real *flux, real rho0, int N, int nel, int nvar){
    int threads_per_block = 256;
    uint32_t ndof = (N+1)*(N+1)*6*nel;
    int nblocks_x = ndof/threads_per_block +1;

    dim3 nblocks(nblocks_x,1,1);
    dim3 nthreads(threads_per_block,1,1);

    boundaryflux_LinearEuler3D_kernel<<<nblocks,nthreads>>>(fb,extfb,nhat,nmag,flux,rho0,ndof);
  }
}

  __global__ void fluxmethod_LinearEuler3D_gpukernel(real *solution, real *flux, real rho0, int ndof, int nvar){
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;

  if( idof < ndof ){
    real u = solution[idof + ndof];
    real v = solution[idof + 2*ndof];
    real w = solution[idof + 3*ndof];
    real p = solution[idof + 4*ndof];
    real c = solution[idof + 5*ndof];

    flux[idof + ndof*(0 + nvar*0)] = rho0*u; // density, x flux ; rho0*u
    flux[idof + ndof*(0 + nvar*1)] = rho0*v; // density, y flux ; rho0*v
    flux[idof + ndof*(0 + nvar*2)] = rho0*w; // density, z flux ; rho0*w

    flux[idof + ndof*(1 + nvar*0)] = p/rho0; // x-velocity, x flux; p/rho0
    flux[idof + ndof*(1 + nvar*1)] = 0.0; // x-velocity, y flux; 0
    flux[idof + ndof*(1 + nvar*2)] = 0.0; // x-velocity, z flux; 0

    flux[idof + ndof*(2 + nvar*0)] = 0.0; // y-velocity, x flux; 0
    flux[idof + ndof*(2 + nvar*1)] = p/rho0; // y-velocity, y flux; p/rho0
    flux[idof + ndof*(2 + nvar*2)] = 0.0; // y-velocity, z flux; 0

    flux[idof + ndof*(3 + nvar*0)] = 0.0; // z-velocity, x flux; 0
    flux[idof + ndof*(3 + nvar*1)] = 0.0; // z-velocity, y flux; 0
    flux[idof + ndof*(3 + nvar*2)] = p/rho0; // z-velocity, z flux; p/rho0

    flux[idof + ndof*(4 + nvar*0)] = c*c*rho0*u; // pressure, x flux : rho0*c^2*u
    flux[idof + ndof*(4 + nvar*1)] = c*c*rho0*v; // pressure, y flux : rho0*c^2*v
    flux[idof + ndof*(4 + nvar*2)] = c*c*rho0*w; // pressure, z flux : rho0*c^2*w

    flux[idof + ndof*(5 + nvar*0)] = 0.0; // sound speed, x flux; 0 (c held fixed in time)
    flux[idof + ndof*(5 + nvar*1)] = 0.0; // sound speed, y flux; 0
    flux[idof + ndof*(5 + nvar*2)] = 0.0; // sound speed, z flux; 0

  }

}
extern "C"
{
  void fluxmethod_LinearEuler3D_gpu(real *solution, real *flux, real rho0, int N, int nel, int nvar){
    int ndof = (N+1)*(N+1)*(N+1)*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
    fluxmethod_LinearEuler3D_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(solution,flux,rho0,ndof,nvar);
  }

}
// Radiation BC kernel for 3D Linear Euler
// Zeroes the acoustic perturbation (rho,u,v,w,p) in extBoundary on
// pre-filtered boundary faces; the sound speed (variable 5, 0-based) is
// copied from the interior side so face Riemann fluxes see a consistent c.
__global__ void hbc3d_radiation_lineareuler3d_kernel(
    real *extBoundary, real *boundary,
    int *elements, int *sides,
    int nBoundaries, int N, int nel)
{
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t dofs_per_face = (N+1)*(N+1);
  uint32_t total_dofs = nBoundaries * dofs_per_face;

  if(idof < total_dofs){
    uint32_t i  = idof % (N+1);
    uint32_t j  = (idof / (N+1)) % (N+1);
    uint32_t n  = idof / dofs_per_face;
    uint32_t e1 = elements[n] - 1;
    uint32_t s1 = sides[n] - 1;

    extBoundary[SCB_3D_INDEX(i,j,s1,e1,0,N,nel)] = 0.0;
    extBoundary[SCB_3D_INDEX(i,j,s1,e1,1,N,nel)] = 0.0;
    extBoundary[SCB_3D_INDEX(i,j,s1,e1,2,N,nel)] = 0.0;
    extBoundary[SCB_3D_INDEX(i,j,s1,e1,3,N,nel)] = 0.0;
    extBoundary[SCB_3D_INDEX(i,j,s1,e1,4,N,nel)] = 0.0;
    extBoundary[SCB_3D_INDEX(i,j,s1,e1,5,N,nel)] = boundary[SCB_3D_INDEX(i,j,s1,e1,5,N,nel)]; // c preserved
  }
}

extern "C"
{
  void hbc3d_radiation_lineareuler3d_gpu(
      real *extBoundary, real *boundary,
      int *elements, int *sides,
      int nBoundaries, int N, int nel)
  {
    int threads_per_block = 256;
    int dofs_per_face = (N+1)*(N+1);
    int total_dofs = nBoundaries * dofs_per_face;
    int nblocks_x = total_dofs/threads_per_block + 1;
    hbc3d_radiation_lineareuler3d_kernel<<<dim3(nblocks_x,1,1),
      dim3(threads_per_block,1,1), 0, 0>>>(extBoundary, boundary,
        elements, sides, nBoundaries, N, nel);
  }
}
