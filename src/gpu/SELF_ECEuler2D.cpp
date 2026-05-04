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
#include <math.h>

/*
 * Equation of state: p = p0 * (rho * Rd * theta / p0)^gamma
 */
__device__ inline real eos_pressure_2d(real rho, real theta, real p0, real Rd, real gamma)
{
  return p0 * pow(rho * Rd * theta / p0, gamma);
}

/*
 * Logarithmic mean used by the Souza et al. (2023, JAMES)
 * entropy-conservative two-point flux for compressible Euler with
 * potential temperature:
 *
 *   <a>_log = (a - b) / (ln(a) - ln(b))   for a != b
 *   <a>_log = (a + b) / 2                 for a == b
 *
 * Implemented with a Taylor series in u = ((a-b)/(a+b))^2 near a == b
 * to avoid 0/0.
 */
__device__ inline real log_mean_dev_2d(real a, real b)
{
  real zeta = a / b;
  real f = (zeta - (real)1.0) / (zeta + (real)1.0);
  real u = f * f;
  real F_F;
  if (u < (real)1.0e-2) {
    F_F = (real)1.0 + u/(real)3.0 + u*u/(real)5.0 + u*u*u/(real)7.0;
  } else {
    F_F = log(zeta) / ((real)2.0 * f);
  }
  return (a + b) / ((real)2.0 * F_F);
}

// ============================================================
// No-normal-flow boundary condition for EC Euler 2D
//
// Mirrors density (var 0), rho*theta (var 3), and Phi (var 4).
// Reflects momentum: (rho*v)_ext = (rho*v)_int - 2*((rho*v).n)*n
// ============================================================

__global__ void hbc2d_nonormalflow_eceuler2d_kernel(
    real *extBoundary, real *boundary, real *nhat,
    int *elements, int *sides,
    int nBoundaries, int N, int nel, int nvar)
{
  uint32_t idof = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t dofs_per_face = (N + 1);
  uint32_t total_dofs = nBoundaries * dofs_per_face;

  if (idof < total_dofs) {
    uint32_t i  = idof % (N + 1);
    uint32_t n  = idof / dofs_per_face;
    uint32_t e1 = elements[n] - 1;
    uint32_t s1 = sides[n] - 1;

    // Load unit normal (2 components)
    real nx = nhat[VEB_2D_INDEX(i, s1, e1, 0, 0, N, nel, 1)];
    real ny = nhat[VEB_2D_INDEX(i, s1, e1, 0, 1, N, nel, 1)];

    // Mirror density
    extBoundary[SCB_2D_INDEX(i, s1, e1, 0, N, nel)] =
        boundary[SCB_2D_INDEX(i, s1, e1, 0, N, nel)];

    // Load interior momentum (2 components)
    real rhou = boundary[SCB_2D_INDEX(i, s1, e1, 1, N, nel)];
    real rhov = boundary[SCB_2D_INDEX(i, s1, e1, 2, N, nel)];

    // Normal momentum component
    real rhovn = rhou * nx + rhov * ny;

    // Reflect: ext = int - 2*(int.n)*n
    extBoundary[SCB_2D_INDEX(i, s1, e1, 1, N, nel)] = rhou - (real)2.0 * rhovn * nx;
    extBoundary[SCB_2D_INDEX(i, s1, e1, 2, N, nel)] = rhov - (real)2.0 * rhovn * ny;

    // Mirror rho*theta
    extBoundary[SCB_2D_INDEX(i, s1, e1, 3, N, nel)] =
        boundary[SCB_2D_INDEX(i, s1, e1, 3, N, nel)];

    // Mirror geopotential (Phi depends on y only)
    extBoundary[SCB_2D_INDEX(i, s1, e1, 4, N, nel)] =
        boundary[SCB_2D_INDEX(i, s1, e1, 4, N, nel)];
    (void)nvar;
  }
}

extern "C"
{
  void hbc2d_nonormalflow_eceuler2d_gpu(
      real *extBoundary, real *boundary, real *nhat,
      int *elements, int *sides,
      int nBoundaries, int N, int nel, int nvar)
  {
    int dofs_per_face = (N + 1);
    int total_dofs = nBoundaries * dofs_per_face;
    int threads_per_block = 256;
    int nblocks_x = total_dofs / threads_per_block + 1;
    hbc2d_nonormalflow_eceuler2d_kernel<<<dim3(nblocks_x, 1, 1),
        dim3(threads_per_block, 1, 1), 0, 0>>>(
        extBoundary, boundary, nhat,
        elements, sides, nBoundaries, N, nel, nvar);
  }
}

// ============================================================
// Parabolic no-stress / no-heat-flux BC.
// Reflects the normal component of the solution gradient at every
// wall node so BR1 averaging gives avgGrad . n = 0.
//   grad_ext = grad_int - 2*(grad_int . n)*n
// ============================================================

__global__ void pbc2d_nostress_eceuler2d_kernel(
    real *extGrad, real *grad, real *nhat,
    int *elements, int *sides,
    int nBoundaries, int N, int nel, int nvar)
{
  uint32_t idof = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t dofs_per_face = (N + 1);
  uint32_t total_dofs = nBoundaries * dofs_per_face;

  if (idof < total_dofs) {
    uint32_t i  = idof % (N + 1);
    uint32_t n  = idof / dofs_per_face;
    uint32_t e1 = elements[n] - 1;
    uint32_t s1 = sides[n] - 1;

    real nx = nhat[VEB_2D_INDEX(i, s1, e1, 0, 0, N, nel, 1)];
    real ny = nhat[VEB_2D_INDEX(i, s1, e1, 0, 1, N, nel, 1)];

    for (int iVar = 0; iVar < nvar; iVar++) {
      real gx = grad[VEB_2D_INDEX(i, s1, e1, iVar, 0, N, nel, nvar)];
      real gy = grad[VEB_2D_INDEX(i, s1, e1, iVar, 1, N, nel, nvar)];
      real gn = gx * nx + gy * ny;
      extGrad[VEB_2D_INDEX(i, s1, e1, iVar, 0, N, nel, nvar)] = gx - (real)2.0 * gn * nx;
      extGrad[VEB_2D_INDEX(i, s1, e1, iVar, 1, N, nel, nvar)] = gy - (real)2.0 * gn * ny;
    }
  }
}

extern "C"
{
  void pbc2d_nostress_eceuler2d_gpu(
      real *extGrad, real *grad, real *nhat,
      int *elements, int *sides,
      int nBoundaries, int N, int nel, int nvar)
  {
    int dofs_per_face = (N + 1);
    int total_dofs = nBoundaries * dofs_per_face;
    int threads_per_block = 256;
    int nblocks_x = total_dofs / threads_per_block + 1;
    pbc2d_nostress_eceuler2d_kernel<<<dim3(nblocks_x, 1, 1),
        dim3(threads_per_block, 1, 1), 0, 0>>>(
        extGrad, grad, nhat,
        elements, sides, nBoundaries, N, nel, nvar);
  }
}

// ============================================================
// LMARS boundary flux for EC Euler 2D
//
// Variables: [rho, rhou, rhov, rhotheta, Phi]
// Geopotential (var 4) carries no flux. Gravity is folded into
// SourceMethod via the Souza non-conservative term so no
// hydrostatic pressure split is applied here.
// ============================================================

__global__ void boundaryflux_eceuler2d_kernel(
    real *fb, real *fextb, real *nhat, real *nscale, real *flux,
    real p0, real Rd, real gamma, int N, int nel)
{
  uint32_t idof = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t ndof = (N + 1) * 4 * nel;

  if (idof < ndof) {
    real nx   = nhat[idof];
    real ny   = nhat[idof + ndof];
    real nmag = nscale[idof];

    // Left state (interior)
    real rhoL    = fb[idof];
    real rhouL   = fb[idof + ndof];
    real rhovL   = fb[idof + 2 * ndof];
    real rthetaL = fb[idof + 3 * ndof];

    real uL = rhouL / rhoL;
    real vL = rhovL / rhoL;
    real thetaL = rthetaL / rhoL;
    real pL = eos_pressure_2d(rhoL, thetaL, p0, Rd, gamma);
    real unL = uL * nx + vL * ny;
    real cL = sqrt(gamma * pL / rhoL);

    // Right state (exterior)
    real rhoR    = fextb[idof];
    real rhouR   = fextb[idof + ndof];
    real rhovR   = fextb[idof + 2 * ndof];
    real rthetaR = fextb[idof + 3 * ndof];

    real uR = rhouR / rhoR;
    real vR = rhovR / rhoR;
    real thetaR = rthetaR / rhoR;
    real pR = eos_pressure_2d(rhoR, thetaR, p0, Rd, gamma);
    real unR = uR * nx + vR * ny;
    real cR = sqrt(gamma * pR / rhoR);

    // LMARS interface velocity and pressure.
    real rho_bar = (real)0.5 * (rhoL + rhoR);
    real c_bar   = (real)0.5 * (cL + cR);
    real rc      = rho_bar * c_bar;

    real un_star = (real)0.5 * (unL + unR) - (pR - pL) / ((real)2.0 * rc);
    real p_star  = (real)0.5 * (pL + pR) - (real)0.5 * rc * (unR - unL);

    // Upwind on un_star
    real s0, s1, s2, s3;
    if (un_star >= (real)0) {
      s0 = rhoL;   s1 = rhouL;  s2 = rhovL;  s3 = rthetaL;
    } else {
      s0 = rhoR;   s1 = rhouR;  s2 = rhovR;  s3 = rthetaR;
    }

    flux[idof]            = (s0 * un_star) * nmag;
    flux[idof + ndof]     = (s1 * un_star + p_star * nx) * nmag;
    flux[idof + 2 * ndof] = (s2 * un_star + p_star * ny) * nmag;
    flux[idof + 3 * ndof] = (s3 * un_star) * nmag;
    flux[idof + 4 * ndof] = (real)0;  // geopotential — no flux
  }
}

extern "C"
{
  void boundaryflux_eceuler2d_gpu(
      real *fb, real *fextb, real *nhat, real *nscale, real *flux,
      real p0, real Rd, real gamma, int N, int nel)
  {
    int ndof = (N + 1) * 4 * nel;
    int threads_per_block = 256;
    int nblocks_x = ndof / threads_per_block + 1;
    boundaryflux_eceuler2d_kernel<<<dim3(nblocks_x, 1, 1),
        dim3(threads_per_block, 1, 1), 0, 0>>>(
        fb, fextb, nhat, nscale, flux, p0, Rd, gamma, N, nel);
  }
}

// ============================================================
// Souza et al. (2023, JAMES) entropy-conservative two-point flux for
// EC Euler 2D with potential temperature on curvilinear mesh.
//
// Physical flux:
//   f_d(rho)     = <rho>_log * <v_d>
//   f_d(rho*v_i) = <rho>_log * <v_i> * <v_d> + <p> * delta_{id}
//   f_d(rho*th)  = <rho*theta>_log * <v_d>
//   f_d(Phi)     = 0
//
// nvar = 5: (rho, rho*u, rho*v, rho*theta, Phi). Phi has zero flux;
// gravity is handled by the Souza non-conservative source term.
//
// dsdx layout: dsdx[iq + nq*(iel + nel*(d + 2*r))] for Ja^r_d
//   slot 0 = (d=0, r=0) -> Ja^1_x
//   slot 1 = (d=1, r=0) -> Ja^1_y
//   slot 2 = (d=0, r=1) -> Ja^2_x
//   slot 3 = (d=1, r=1) -> Ja^2_y
// ============================================================

template <int blockSize, int matSize>
__global__ void __launch_bounds__(matSize) twopointfluxmethod_eceuler2d_kernel(
    real *f, real *s, real *dsdx,
    real p0, real Rd, real gamma,
    int nq, int N, int nvar)
{
  uint32_t iq = threadIdx.x;
  if (iq < (uint32_t)nq) {
    uint32_t iel = blockIdx.x;
    uint32_t nel = gridDim.x;
    uint32_t ii  = iq % (N + 1);
    uint32_t jj  = iq / (N + 1);

    // Load state at (i,j)
    real rho_ijk    = s[iq + nq * (iel + nel * 0)];
    real rhou_ijk   = s[iq + nq * (iel + nel * 1)];
    real rhov_ijk   = s[iq + nq * (iel + nel * 2)];
    real rtheta_ijk = s[iq + nq * (iel + nel * 3)];

    real u_ijk     = rhou_ijk / rho_ijk;
    real v_ijk     = rhov_ijk / rho_ijk;
    real theta_ijk = rtheta_ijk / rho_ijk;
    real p_ijk     = eos_pressure_2d(rho_ijk, theta_ijk, p0, Rd, gamma);

    int nq3 = nq * (N + 1); // stride between idir slices in TPV layout

    for (int nn = 0; nn < N + 1; nn++) {

      // ============== xi^1: pair (i,j)-(nn,j) ==============
      uint32_t iq1 = nn + (N + 1) * jj;
      real rho1    = s[iq1 + nq * (iel + nel * 0)];
      real rhou1   = s[iq1 + nq * (iel + nel * 1)];
      real rhov1   = s[iq1 + nq * (iel + nel * 2)];
      real rtheta1 = s[iq1 + nq * (iel + nel * 3)];
      real u1      = rhou1 / rho1;
      real v1      = rhov1 / rho1;
      real theta1  = rtheta1 / rho1;
      real p1      = eos_pressure_2d(rho1, theta1, p0, Rd, gamma);

      real rho_log = log_mean_dev_2d(rho_ijk, rho1);
      real rth_log = log_mean_dev_2d(rtheta_ijk, rtheta1);
      real u_a     = (real)0.5 * (u_ijk + u1);
      real v_a     = (real)0.5 * (v_ijk + v1);
      real p_a     = (real)0.5 * (p_ijk + p1);

      real Ja1_x = (real)0.5 * (dsdx[iq  + nq * (iel + nel * (0 + 2 * 0))] +
                                 dsdx[iq1 + nq * (iel + nel * (0 + 2 * 0))]);
      real Ja1_y = (real)0.5 * (dsdx[iq  + nq * (iel + nel * (1 + 2 * 0))] +
                                 dsdx[iq1 + nq * (iel + nel * (1 + 2 * 0))]);

      for (int ivar = 0; ivar < 4; ivar++) {
        real Fphys_x, Fphys_y;
        if (ivar == 0) {
          Fphys_x = rho_log * u_a;
          Fphys_y = rho_log * v_a;
        } else if (ivar == 1) {
          Fphys_x = rho_log * u_a * u_a + p_a;
          Fphys_y = rho_log * u_a * v_a;
        } else if (ivar == 2) {
          Fphys_x = rho_log * v_a * u_a;
          Fphys_y = rho_log * v_a * v_a + p_a;
        } else {
          Fphys_x = rth_log * u_a;
          Fphys_y = rth_log * v_a;
        }
        f[nn + (N + 1) * iq + nq3 * (iel + nel * (ivar + nvar * 0))] =
            Ja1_x * Fphys_x + Ja1_y * Fphys_y;
      }
      // Geopotential (var 4) carries no flux.
      f[nn + (N + 1) * iq + nq3 * (iel + nel * (4 + nvar * 0))] = (real)0;

      // ============== xi^2: pair (i,j)-(i,nn) ==============
      uint32_t iq2 = ii + (N + 1) * nn;
      real rho2    = s[iq2 + nq * (iel + nel * 0)];
      real rhou2   = s[iq2 + nq * (iel + nel * 1)];
      real rhov2   = s[iq2 + nq * (iel + nel * 2)];
      real rtheta2 = s[iq2 + nq * (iel + nel * 3)];
      real u2      = rhou2 / rho2;
      real v2      = rhov2 / rho2;
      real theta2  = rtheta2 / rho2;
      real p2      = eos_pressure_2d(rho2, theta2, p0, Rd, gamma);

      rho_log = log_mean_dev_2d(rho_ijk, rho2);
      rth_log = log_mean_dev_2d(rtheta_ijk, rtheta2);
      u_a     = (real)0.5 * (u_ijk + u2);
      v_a     = (real)0.5 * (v_ijk + v2);
      p_a     = (real)0.5 * (p_ijk + p2);

      real Ja2_x = (real)0.5 * (dsdx[iq  + nq * (iel + nel * (0 + 2 * 1))] +
                                 dsdx[iq2 + nq * (iel + nel * (0 + 2 * 1))]);
      real Ja2_y = (real)0.5 * (dsdx[iq  + nq * (iel + nel * (1 + 2 * 1))] +
                                 dsdx[iq2 + nq * (iel + nel * (1 + 2 * 1))]);

      for (int ivar = 0; ivar < 4; ivar++) {
        real Fphys_x, Fphys_y;
        if (ivar == 0) {
          Fphys_x = rho_log * u_a;
          Fphys_y = rho_log * v_a;
        } else if (ivar == 1) {
          Fphys_x = rho_log * u_a * u_a + p_a;
          Fphys_y = rho_log * u_a * v_a;
        } else if (ivar == 2) {
          Fphys_x = rho_log * v_a * u_a;
          Fphys_y = rho_log * v_a * v_a + p_a;
        } else {
          Fphys_x = rth_log * u_a;
          Fphys_y = rth_log * v_a;
        }
        f[nn + (N + 1) * iq + nq3 * (iel + nel * (ivar + nvar * 1))] =
            Ja2_x * Fphys_x + Ja2_y * Fphys_y;
      }
      f[nn + (N + 1) * iq + nq3 * (iel + nel * (4 + nvar * 1))] = (real)0;
    }
  }
  (void)blockSize;
}

extern "C"
{
  void twopointfluxmethod_eceuler2d_gpu(
      real *f, real *s, real *dsdx,
      real p0, real Rd, real gamma,
      int N, int nvar, int nel)
  {
    int nq = (N + 1) * (N + 1);
    if (N < 7) {
      twopointfluxmethod_eceuler2d_kernel<64, 64><<<dim3(nel, 1, 1),
          dim3(64, 1, 1), 0, 0>>>(f, s, dsdx, p0, Rd, gamma, nq, N, nvar);
    } else {
      twopointfluxmethod_eceuler2d_kernel<256, 256><<<dim3(nel, 1, 1),
          dim3(256, 1, 1), 0, 0>>>(f, s, dsdx, p0, Rd, gamma, nq, N, nvar);
    }
  }
}

// ============================================================
// Souza et al. (2023) non-conservative gravity flux differencing.
// At each receiving node (i,j,iel) compute
//   src(rhov)|_i = -(1/J_i) sum_r sum_j D_split^r[i_r, j]
//                  * 0.5*(Ja^r_y(i)+Ja^r_y(j))
//                  * <rho>_log(s_i, s_j) * (Phi_j - Phi_i)
// using the geopotential Phi carried as state variable index 4
// (0-indexed). All other source components are set to zero.
//
// In 2D, gravity acts on the rho*v equation (variable index 2,
// 0-indexed) since y is the vertical coordinate.
// ============================================================
template <int blockSize, int matSize>
__global__ void __launch_bounds__(matSize) sourcemethod_eceuler2d_kernel(
    real *source, real *solution, real *dsdx, real *J, real *dSplit,
    int nq, int N, int nvar)
{
  uint32_t iq = threadIdx.x;
  if (iq < (uint32_t)nq) {
    uint32_t iel = blockIdx.x;
    uint32_t nel = gridDim.x;
    uint32_t ii  = iq % (N + 1);
    uint32_t jj  = iq / (N + 1);

    real rho_ijk = solution[iq + nq * (iel + nel * 0)];
    real phi_ijk = solution[iq + nq * (iel + nel * 4)];

    real acc = (real)0;
    for (int nn = 0; nn < N + 1; nn++) {
      // xi^1 partner (nn, jj) — uses Ja^1_y (slot d=1, r=0 -> idx 1)
      uint32_t iq1 = nn + (N + 1) * jj;
      real rho_p = solution[iq1 + nq * (iel + nel * 0)];
      real phi_p = solution[iq1 + nq * (iel + nel * 4)];
      real rho_log = log_mean_dev_2d(rho_ijk, rho_p);
      real Ja_y_avg = (real)0.5 *
                      (dsdx[iq  + nq * (iel + nel * (1 + 2 * 0))]
                     + dsdx[iq1 + nq * (iel + nel * (1 + 2 * 0))]);
      real D_r = dSplit[nn + (N + 1) * ii];
      acc += D_r * Ja_y_avg * rho_log * (phi_p - phi_ijk);

      // xi^2 partner (ii, nn) — Ja^2_y slot d=1, r=1 -> idx 3
      uint32_t iq2 = ii + (N + 1) * nn;
      rho_p = solution[iq2 + nq * (iel + nel * 0)];
      phi_p = solution[iq2 + nq * (iel + nel * 4)];
      rho_log = log_mean_dev_2d(rho_ijk, rho_p);
      Ja_y_avg = (real)0.5 *
                 (dsdx[iq  + nq * (iel + nel * (1 + 2 * 1))]
                + dsdx[iq2 + nq * (iel + nel * (1 + 2 * 1))]);
      D_r = dSplit[nn + (N + 1) * jj];
      acc += D_r * Ja_y_avg * rho_log * (phi_p - phi_ijk);
    }

    real jac = J[iq + nq * iel];
    source[iq + nq * (iel + nel * 0)] = (real)0;          // rho
    source[iq + nq * (iel + nel * 1)] = (real)0;          // rhou
    source[iq + nq * (iel + nel * 2)] = -acc / jac;       // rhov
    source[iq + nq * (iel + nel * 3)] = (real)0;          // rhotheta
    source[iq + nq * (iel + nel * 4)] = (real)0;          // Phi
    (void)nvar;
  }
  (void)blockSize;
}

extern "C"
{
  void sourcemethod_eceuler2d_gpu(
      real *source, real *solution, real *dsdx, real *J, real *dSplit,
      int N, int nvar, int nel)
  {
    int nq = (N + 1) * (N + 1);
    if (N < 7) {
      sourcemethod_eceuler2d_kernel<64, 64><<<dim3(nel, 1, 1), dim3(64, 1, 1), 0, 0>>>(
          source, solution, dsdx, J, dSplit, nq, N, nvar);
    } else {
      sourcemethod_eceuler2d_kernel<256, 256><<<dim3(nel, 1, 1), dim3(256, 1, 1), 0, 0>>>(
          source, solution, dsdx, J, dSplit, nq, N, nvar);
    }
  }
}

// ============================================================
// Diffusive flux for ECEuler2D — constant-coefficient Laplacian
// (Bassi-Rebay 1) with separate momentum (nu) and thermal (kappa)
// diffusivities.
//
// Interior fill:
//   F_d(rho)       = 0
//   F_d(rho*v_i)   = -nu    * d(rho*v_i)/dx_d   for i = 1,2
//   F_d(rho*theta) = -kappa * d(rho*theta)/dx_d
//   F_d(Phi)       = 0
//
// Variable indexing (0-indexed): 0=rho, 1=rhou, 2=rhov, 3=rhotheta, 4=Phi
// ============================================================

__global__ void diffusiveflux_eceuler2d_kernel(
    real *diffFlux, real *grad, real nu, real kappa,
    uint32_t ndof_node, uint32_t nvar)
{
  uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof_total = ndof_node * nvar * 2;
  if (i < ndof_total) {
    uint32_t var = (i / ndof_node) % nvar;
    real coeff;
    if (var == 0 || var == 4) coeff = (real)0;
    else if (var == 3) coeff = kappa;
    else coeff = nu;
    diffFlux[i] = -coeff * grad[i];
  }
}

extern "C"
{
  void diffusiveflux_eceuler2d_gpu(
      real *diffFlux, real *grad, real nu, real kappa,
      int N, int nvar, int nel)
  {
    uint32_t ndof_node = (N+1) * (N+1) * (uint32_t)nel;
    uint32_t ndof_total = ndof_node * (uint32_t)nvar * 2;
    uint32_t nthreads = 256;
    uint32_t nblocks  = (ndof_total + nthreads - 1) / nthreads;
    diffusiveflux_eceuler2d_kernel<<<dim3(nblocks,1,1), dim3(nthreads,1,1), 0, 0>>>(
        diffFlux, grad, nu, kappa, ndof_node, (uint32_t)nvar);
  }
}

// ============================================================
// Boundary BR1 + SIPG flux for the parabolic terms:
//   fluxN = -coeff*(avgGrad . n)*nmag + tau*(uL - uR)*nmag
// tau_nu, tau_kappa precomputed on host.
// ============================================================

__global__ void diffusiveboundaryflux_eceuler2d_kernel(
    real *fluxN, real *avgGrad, real *uBnd, real *uExt,
    real *nhat, real *nscale,
    real nu, real kappa, real tau_nu, real tau_kappa,
    uint32_t ndof_face, uint32_t nvar)
{
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  if (idof < ndof_face) {
    real nx   = nhat[idof];
    real ny   = nhat[idof + ndof_face];
    real nmag = nscale[idof];

    for (uint32_t var = 0; var < nvar; var++) {
      real coeff, tau;
      if (var == 0 || var == 4) {
        coeff = (real)0;
        tau   = (real)0;
      } else if (var == 3) {
        coeff = kappa;
        tau   = tau_kappa;
      } else {
        coeff = nu;
        tau   = tau_nu;
      }

      // avgGrad components for this variable: stride between d-slabs
      // is ndof_face*nvar.
      real gx = avgGrad[idof + ndof_face*(var + nvar*0)];
      real gy = avgGrad[idof + ndof_face*(var + nvar*1)];
      real gn = gx*nx + gy*ny;
      real uL = uBnd[idof + ndof_face*var];
      real uR = uExt[idof + ndof_face*var];
      fluxN[idof + ndof_face*var] = (-coeff*gn + tau*(uL - uR)) * nmag;
    }
  }
}

extern "C"
{
  void diffusiveboundaryflux_eceuler2d_gpu(
      real *fluxN, real *avgGrad, real *uBnd, real *uExt,
      real *nhat, real *nscale,
      real nu, real kappa, real tau_nu, real tau_kappa,
      int N, int nvar, int nel)
  {
    uint32_t ndof_face = (N+1) * 4 * (uint32_t)nel;
    uint32_t nthreads = 256;
    uint32_t nblocks  = (ndof_face + nthreads - 1) / nthreads;
    diffusiveboundaryflux_eceuler2d_kernel<<<dim3(nblocks,1,1), dim3(nthreads,1,1), 0, 0>>>(
        fluxN, avgGrad, uBnd, uExt, nhat, nscale,
        nu, kappa, tau_nu, tau_kappa, ndof_face, (uint32_t)nvar);
  }
}
