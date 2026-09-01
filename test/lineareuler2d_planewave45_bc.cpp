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
// Time-dependent prescribed BC for the 45-degree plane wave used by
// lineareuler2d_planewave45_dtconvergence.
//
// The exterior state on every prescribed face is set to the exact plane wave
//
//   phase = kx*x + ky*y - omega*t
//   P = amp*sin(phase)
//   u = amp*(kx/|k|)*sin(phase)/(rho0*c)
//   v = amp*(ky/|k|)*sin(phase)/(rho0*c)
//
// evaluated at the current model time. Variable layout (0-based):
// 0=u, 1=v, 2=p, 3=c, 4=rho0, 5=sigma. This mirrors
// hbc2d_Prescribed_planewave45 on the host; the two must agree to round-off.
// ============================================================
__global__ void hbc2d_planewave45_lineareuler2d_kernel(
    real *extBoundary, real *xBoundary,
    int *elements, int *sides,
    int nBoundaries, int N, int nel,
    real kx, real ky, real omega, real amp, real rho0, real c, real t)
{
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t total_dofs = nBoundaries * (N+1);

  if(idof < total_dofs){
    uint32_t i  = idof % (N+1);
    uint32_t n  = idof / (N+1);
    uint32_t e1 = elements[n] - 1; // Fortran 1-based to C 0-based
    uint32_t s1 = sides[n] - 1;

    real x = xBoundary[VEB_2D_INDEX(i,s1,e1,0,0,N,nel,1)];
    real y = xBoundary[VEB_2D_INDEX(i,s1,e1,0,1,N,nel,1)];

    real kmag   = sqrt(kx*kx + ky*ky);
    real sphase = sin(kx*x + ky*y - omega*t);

    extBoundary[SCB_2D_INDEX(i,s1,e1,0,N,nel)] = amp*(kx/kmag)*sphase/(rho0*c); // u
    extBoundary[SCB_2D_INDEX(i,s1,e1,1,N,nel)] = amp*(ky/kmag)*sphase/(rho0*c); // v
    extBoundary[SCB_2D_INDEX(i,s1,e1,2,N,nel)] = amp*sphase;                    // p
    extBoundary[SCB_2D_INDEX(i,s1,e1,3,N,nel)] = c;                             // sound speed
    extBoundary[SCB_2D_INDEX(i,s1,e1,4,N,nel)] = rho0;                          // background density
    extBoundary[SCB_2D_INDEX(i,s1,e1,5,N,nel)] = 0.0;                           // relaxation rate
  }
}

extern "C"
{
  void hbc2d_planewave45_lineareuler2d_gpu(
      real *extBoundary, real *xBoundary,
      int *elements, int *sides,
      int nBoundaries, int N, int nel,
      real kx, real ky, real omega, real amp, real rho0, real c, real t)
  {
    int threads_per_block = 256;
    int total_dofs = nBoundaries * (N+1);
    int nblocks_x = total_dofs/threads_per_block + 1;
    hbc2d_planewave45_lineareuler2d_kernel<<<dim3(nblocks_x,1,1),
      dim3(threads_per_block,1,1), 0, 0>>>(extBoundary, xBoundary,
        elements, sides, nBoundaries, N, nel,
        kx, ky, omega, amp, rho0, c, t);
  }
}
