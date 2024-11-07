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


__global__ void boundaryflux_LinearEuler2D_kernel(real *fb, real *extfb, real *nhat, real *nmag, real *flux, real rho0, real c, int ndof){
  uint32_t idof = threadIdx.x;
   
  if( idof < ndof ){

    real fl[4];
    real nx = nhat[idof];
    real ny = nhat[idof+ndof];
    real un = fb[idof + ndof]*nx + fb[idof + 2*ndof]*ny;
    real p = fb[idof + 3*ndof]; 

    fl[0] = rho0*un; // density flux
    fl[1] = p*nx/rho0; // x-momentum flux
    fl[2] = p*ny/rho0; // y-momentum flux
    fl[3] = rho0*c*c*un; // pressure flux

    real fr[4];
    un = extfb[idof + ndof]*nx + extfb[idof + 2*ndof]*ny;
    p = extfb[idof + 3*ndof]; 
    
    fr[0] = rho0*un; // density flux
    fr[1] = p*nx/rho0; // x-momentum flux
    fr[2] = p*ny/rho0; // y-momentum flux
    fr[3] = rho0*c*c*un; // pressure flux

    real nm = nmag[idof];
    flux[idof] = (0.5*(fl[0]+fr[0])+c*(fb[idof]-extfb[idof]))*nm; // density
    flux[idof+ndof] = (0.5*(fl[0]+fr[0])+c*(fb[idof+ndof]-extfb[idof+ndof]))*nm; // u
    flux[idof+2*ndof] = (0.5*(fl[0]+fr[0])+c*(fb[idof+2*ndof]-extfb[idof+2*ndof]))*nm; // v
    flux[idof+3*ndof] = (0.5*(fl[0]+fr[0])+c*(fb[idof+3*ndof]-extfb[idof+3*ndof]))*nm; // p
  }
}

extern "C"
{
  void boundaryflux_LinearEuler2D_gpu(real *fb, real *extfb,real *nhat, real *nmag, real *flux, real rho0, real c, int N, int nel, int nvar){
    int threads_per_block = 256;
    uint32_t ndof = (N+1)*4*nel;
    int nblocks_x = ndof/threads_per_block +1;

    dim3 nblocks(nblocks_x,nvar,1);
    dim3 nthreads(threads_per_block,1,1);

    boundaryflux_LinearEuler2D_kernel<<<nblocks,nthreads>>>(fb,extfb,nhat,nmag,flux,rho0,c,ndof);
  }

}
