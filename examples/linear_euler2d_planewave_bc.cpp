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
// Device-resident prescribed boundary conditions for the Gaussian plane wave
// examples (linear_euler2d_planewave_propagation and
// linear_euler2d_planewave_reflection).
//
// Incident wave:
//   phi   = wx*(x-x0) + wy*(y-y0) - c*t
//   shi   = exp(-phi^2/L^2)
//
// Reflected wave (method of images about the eastern wall at x = 2):
//   phr   = -wx*(x+x0-2) + wy*(y-y0) - c*t
//   shr   = exp(-phr^2/L^2)
//
// With reflect = 0 only the incident wave is applied. Variable layout (0-based):
// 0=u, 1=v, 2=p, 3=c, 4=rho0, 5=sigma.
//
// These mirror the host routines in the two examples. The exterior state that the
// Riemann solver reads lives in solution%extBoundary_gpu, so on a GPU build the
// boundary condition has to be written here: filling the host extBoundary array
// leaves the device state untouched and the boundary silently does nothing.
// ============================================================
__global__ void hbc2d_planewave_gaussian_lineareuler2d_kernel(
    real *extBoundary, real *xBoundary,
    int *elements, int *sides,
    int nBoundaries, int N, int nel,
    real wx, real wy, real amp, real x0, real y0, real L,
    real rho0, real c, real t, int reflect)
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

    real u = amp*wx/c;
    real v = amp*wy/c;

    real phi = wx*(x-x0) + wy*(y-y0) - c*t;
    real shi = exp(-phi*phi/(L*L));

    real shr = 0.0;
    if(reflect != 0){
      real phr = -wx*(x+x0-2.0) + wy*(y-y0) - c*t;
      shr = exp(-phr*phr/(L*L));
    }

    extBoundary[SCB_2D_INDEX(i,s1,e1,0,N,nel)] = u*(shi-shr); // u
    extBoundary[SCB_2D_INDEX(i,s1,e1,1,N,nel)] = v*(shi+shr); // v
    extBoundary[SCB_2D_INDEX(i,s1,e1,2,N,nel)] = amp*(shi+shr); // pressure
    extBoundary[SCB_2D_INDEX(i,s1,e1,3,N,nel)] = c; // sound speed
    extBoundary[SCB_2D_INDEX(i,s1,e1,4,N,nel)] = rho0; // background density
    extBoundary[SCB_2D_INDEX(i,s1,e1,5,N,nel)] = 0.0; // relaxation rate
  }
}

extern "C"
{
  void hbc2d_planewave_gaussian_lineareuler2d_gpu(
      real *extBoundary, real *xBoundary,
      int *elements, int *sides,
      int nBoundaries, int N, int nel,
      real wx, real wy, real amp, real x0, real y0, real L,
      real rho0, real c, real t, int reflect)
  {
    int threads_per_block = 256;
    int total_dofs = nBoundaries * (N+1);
    int nblocks_x = total_dofs/threads_per_block + 1;
    hbc2d_planewave_gaussian_lineareuler2d_kernel<<<dim3(nblocks_x,1,1),
      dim3(threads_per_block,1,1), 0, 0>>>(extBoundary, xBoundary,
        elements, sides, nBoundaries, N, nel,
        wx, wy, amp, x0, y0, L, rho0, c, t, reflect);
  }
}
