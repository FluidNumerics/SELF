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

__global__ void boundaryflux_LinearShallowWater2D_kernel(real *fb, real *extfb, real *nhat, real *nmag, real *flux, real g, real H, int ndof){
    uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;

    if( idof < ndof){

        real nx = nhat[idof];
        real ny = nhat[idof+ndof];
        real nm = nmag[idof];

        real fl[3];
        fl[0] = fb[idof];           // uL
        fl[1] = fb[idof + ndof];    // vL
        fl[2] = fb[idof + 2*ndof];  // etaL

        real fr[3];
        fr[0] = extfb[idof];           // uR
        fr[1] = extfb[idof + ndof];    // vR
        fr[2] = extfb[idof + 2*ndof];  // etaR
        
        real unL = fl[0] * nx + fl[1] * ny;
        real unR = fr[0] * nx + fr[1] * ny;

        real c = sqrt(g * H);

        flux[idof]          = 0.5 * (g * (fl[2] + fr[2]) + c * (unL - unR)) * nx * nm;
        flux[idof + ndof]   = 0.5 * (g * (fl[2] + fr[2]) + c * (unL - unR)) * ny * nm;
        flux[idof + 2*ndof] = 0.5 * (H * (unL + unR) + c * (fl[2] - fr[2])) * nm;
    }
}

extern "C"
{
    void boundaryflux_LinearShallowWater2D_gpu(real *fb, real *extfb, real *nhat, real *nmag, real *flux, real g, real H, int N, int nel, int nvar){
        int threads_per_block = 256;
        uint32_t ndof = (N+1)*4*nel;
        int nblocks_x = ndof/threads_per_block + 1;

        dim3 nblocks(nblocks_x, nvar, 1);
        dim3 nthreads(threads_per_block, 1, 1);

        boundaryflux_LinearShallowWater2D_kernel<<<nblocks,nthreads>>>(fb,extfb,nhat,nmag,flux,g,H,ndof);
    }
}

__global__ void fluxmethod_LinearShallowWater2D_gpukernel(real *solution, real *flux, real g, real H, int ndof, int nvar){
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;

  if( idof < ndof ){
    real u = solution[idof];
    real v = solution[idof + ndof];
    real eta = solution[idof + 2*ndof];

    flux[idof + ndof*(0 + nvar*0)] = g*eta; // x-component of u
    flux[idof + ndof*(0 + nvar*1)] = 0.0;   // y-component of u
    flux[idof + ndof*(1 + nvar*0)] = 0.0;   // x-component of v
    flux[idof + ndof*(1 + nvar*1)] = g*eta; // y-component of v
    flux[idof + ndof*(2 + nvar*0)] = H*u;   // x-component of eta
    flux[idof + ndof*(2 + nvar*1)] = H*v;   // y-component of eta

  }

}
extern "C"
{
  void fluxmethod_LinearShallowWater2D_gpu(real *solution, real *flux, real g, real H, int N, int nel, int nvar){
    int ndof = (N+1)*(N+1)*nel;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
    fluxmethod_LinearShallowWater2D_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(solution,flux,g,H,ndof,nvar);
  }
}

__global__ void setboundarycondition_LinearShallowWater2D_gpukernel(real *extBoundary, real *boundary, int *sideInfo, real *nhat, int N, int nEl, int nvar){
    uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
    uint32_t ndof = (N+1)*4*nEl;

    if(idof < ndof){
        uint32_t i = idof % (N+1);
        uint32_t s1 = (idof/(N+1)) % 4;
        uint32_t e1 = idof/(N+1)/4;
        uint32_t e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
        uint32_t bcid = sideInfo[INDEX3(4,s1,e1,5,4)];
        if( e2 == 0){
            if( bcid == SELF_BC_NONORMALFLOW){
                real u   = boundary[SCB_2D_INDEX(i,s1,e1,0,N,nEl)];
                real v   = boundary[SCB_2D_INDEX(i,s1,e1,1,N,nEl)];
                real eta = boundary[SCB_2D_INDEX(i,s1,e1,2,N,nEl)];
                real nx      = nhat[VEB_2D_INDEX(i,s1,e1,0,0,N,nEl,1)];
                real ny      = nhat[VEB_2D_INDEX(i,s1,e1,0,1,N,nEl,1)];

                extBoundary[SCB_2D_INDEX(i,s1,e1,0,N,nEl)] = (ny * ny - nx * nx) * u - 2 * nx * ny * v;
                extBoundary[SCB_2D_INDEX(i,s1,e1,1,N,nEl)] = (nx * nx - ny * ny) * v - 2 * nx * ny * u;
                extBoundary[SCB_2D_INDEX(i,s1,e1,2,N,nEl)] = eta; 
            } else if ( bcid == SELF_BC_RADIATION){
                extBoundary[SCB_2D_INDEX(i,s1,e1,0,N,nEl)] = 0.0;
                extBoundary[SCB_2D_INDEX(i,s1,e1,1,N,nEl)] = 0.0;
                extBoundary[SCB_2D_INDEX(i,s1,e1,2,N,nEl)] = 0.0; 
            }
        }
    }
}

extern "C" 
{
  void setboundarycondition_LinearShallowWater2D_gpu(real *extBoundary, real *boundary, int *sideInfo, real *nhat, int N, int nel, int nvar){
    int threads_per_block = 256;
    int ndof = (N+1)*4*nel;
    int nblocks_x = ndof/threads_per_block +1;

    dim3 nblocks(nblocks_x,1,1);
    dim3 nthreads(threads_per_block,1,1);

	setboundarycondition_LinearShallowWater2D_gpukernel<<<nblocks,nthreads, 0, 0>>>(extBoundary,boundary,sideInfo,nhat,N,nel,nvar);
  }
}