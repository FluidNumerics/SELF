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
__device__ real eos_pressure(real rho, real theta, real p0, real Rd, real gamma)
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
__device__ inline real log_mean_dev(real a, real b)
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
// No-normal-flow boundary condition for EC Euler 3D
//
// Mirrors density (var 0) and rho*theta (var 4).
// Reflects momentum: (rho*v)_ext = (rho*v)_int - 2*((rho*v).n)*n
// ============================================================

__global__ void hbc3d_nonormalflow_eceuler3d_kernel(
    real *extBoundary, real *boundary, real *nhat,
    int *elements, int *sides,
    int nBoundaries, int N, int nel, int nvar)
{
  uint32_t idof = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t dofs_per_face = (N + 1) * (N + 1);
  uint32_t total_dofs = nBoundaries * dofs_per_face;

  if (idof < total_dofs) {
    uint32_t i  = idof % (N + 1);
    uint32_t j  = (idof / (N + 1)) % (N + 1);
    uint32_t n  = idof / dofs_per_face;
    uint32_t e1 = elements[n] - 1;
    uint32_t s1 = sides[n] - 1;

    // Load unit normal
    real nx = nhat[VEB_3D_INDEX(i, j, s1, e1, 0, 0, N, nel, 1)];
    real ny = nhat[VEB_3D_INDEX(i, j, s1, e1, 0, 1, N, nel, 1)];
    real nz = nhat[VEB_3D_INDEX(i, j, s1, e1, 0, 2, N, nel, 1)];

    // Mirror density
    extBoundary[SCB_3D_INDEX(i, j, s1, e1, 0, N, nel)] =
        boundary[SCB_3D_INDEX(i, j, s1, e1, 0, N, nel)];

    // Load interior momentum
    real rhou = boundary[SCB_3D_INDEX(i, j, s1, e1, 1, N, nel)];
    real rhov = boundary[SCB_3D_INDEX(i, j, s1, e1, 2, N, nel)];
    real rhow = boundary[SCB_3D_INDEX(i, j, s1, e1, 3, N, nel)];

    // Normal momentum component
    real rhovn = rhou * nx + rhov * ny + rhow * nz;

    // Reflect: ext = int - 2*(int.n)*n
    extBoundary[SCB_3D_INDEX(i, j, s1, e1, 1, N, nel)] = rhou - 2.0 * rhovn * nx;
    extBoundary[SCB_3D_INDEX(i, j, s1, e1, 2, N, nel)] = rhov - 2.0 * rhovn * ny;
    extBoundary[SCB_3D_INDEX(i, j, s1, e1, 3, N, nel)] = rhow - 2.0 * rhovn * nz;

    // Mirror rho*theta
    extBoundary[SCB_3D_INDEX(i, j, s1, e1, 4, N, nel)] =
        boundary[SCB_3D_INDEX(i, j, s1, e1, 4, N, nel)];
  }
}

extern "C"
{
  void hbc3d_nonormalflow_eceuler3d_gpu(
      real *extBoundary, real *boundary, real *nhat,
      int *elements, int *sides,
      int nBoundaries, int N, int nel, int nvar)
  {
    int dofs_per_face = (N + 1) * (N + 1);
    int total_dofs = nBoundaries * dofs_per_face;
    int threads_per_block = 256;
    int nblocks_x = total_dofs / threads_per_block + 1;
    hbc3d_nonormalflow_eceuler3d_kernel<<<dim3(nblocks_x, 1, 1),
        dim3(threads_per_block, 1, 1), 0, 0>>>(
        extBoundary, boundary, nhat,
        elements, sides, nBoundaries, N, nel, nvar);
  }
}

// ============================================================
// LLF (Rusanov) boundary flux for EC Euler 3D
//
// F* = 0.5*(fL + fR) - 0.5*lambda_max*(sR - sL)
// lambda_max = max(|unL|+cL, |unR|+cR)
// c = sqrt(gamma * p / rho)
//
// Variables: [rho, rhou, rhov, rhow, rhotheta]
// ============================================================

__global__ void boundaryflux_eceuler3d_kernel(
    real *fb, real *fextb, real *nhat, real *nscale, real *flux,
    real *p_hyd_b, real *p_hyd_extb,
    real p0, real Rd, real gamma, int N, int nel)
{
  uint32_t idof = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t ndof = (N + 1) * (N + 1) * 6 * nel;

  if (idof < ndof) {
    // Load unit normal and scale
    real nx   = nhat[idof];
    real ny   = nhat[idof + ndof];
    real nz   = nhat[idof + 2 * ndof];
    real nmag = nscale[idof];
    // Hydrostatic pressure on each side of the face (well-balanced split).
    // Reading both interior and exterior buffers lets us support jump
    // discontinuities in the hydrostatic profile across element edges
    // (e.g., for stretched or non-conforming meshes); for a smooth
    // profile both values agree by construction.
    real p_hyd_L = p_hyd_b[idof];
    real p_hyd_R = p_hyd_extb[idof];

    // Left state (interior)
    real rhoL   = fb[idof];
    real rhouL  = fb[idof + ndof];
    real rhovL  = fb[idof + 2 * ndof];
    real rhowL  = fb[idof + 3 * ndof];
    real rthetaL = fb[idof + 4 * ndof];

    real uL = rhouL / rhoL;
    real vL = rhovL / rhoL;
    real wL = rhowL / rhoL;
    real thetaL = rthetaL / rhoL;
    real pL = eos_pressure(rhoL, thetaL, p0, Rd, gamma);
    real unL = uL * nx + vL * ny + wL * nz;
    real cL = sqrt(gamma * pL / rhoL);

    // Right state (exterior)
    real rhoR   = fextb[idof];
    real rhouR  = fextb[idof + ndof];
    real rhovR  = fextb[idof + 2 * ndof];
    real rhowR  = fextb[idof + 3 * ndof];
    real rthetaR = fextb[idof + 4 * ndof];

    real uR = rhouR / rhoR;
    real vR = rhovR / rhoR;
    real wR = rhowR / rhoR;
    real thetaR = rthetaR / rhoR;
    real pR = eos_pressure(rhoR, thetaR, p0, Rd, gamma);
    real unR = uR * nx + vR * ny + wR * nz;
    real cR = sqrt(gamma * pR / rhoR);

    // Maximum wave speed
    real lam = max(fabs(unL) + cL, fabs(unR) + cR);

    // Normal fluxes with WB pressure split: each side subtracts its
    // own p_hyd so that conservation holds even with discontinuities.
    real fL0 = rhoL * unL;
    real fR0 = rhoR * unR;

    real fL1 = rhouL * unL + (pL - p_hyd_L) * nx;
    real fR1 = rhouR * unR + (pR - p_hyd_R) * nx;

    real fL2 = rhovL * unL + (pL - p_hyd_L) * ny;
    real fR2 = rhovR * unR + (pR - p_hyd_R) * ny;

    real fL3 = rhowL * unL + (pL - p_hyd_L) * nz;
    real fR3 = rhowR * unR + (pR - p_hyd_R) * nz;

    real fL4 = rthetaL * unL;
    real fR4 = rthetaR * unR;

    // LLF flux * nScale
    flux[idof]              = (0.5 * (fL0 + fR0) - 0.5 * lam * (rhoR - rhoL)) * nmag;
    flux[idof + ndof]       = (0.5 * (fL1 + fR1) - 0.5 * lam * (rhouR - rhouL)) * nmag;
    flux[idof + 2 * ndof]   = (0.5 * (fL2 + fR2) - 0.5 * lam * (rhovR - rhovL)) * nmag;
    flux[idof + 3 * ndof]   = (0.5 * (fL3 + fR3) - 0.5 * lam * (rhowR - rhowL)) * nmag;
    flux[idof + 4 * ndof]   = (0.5 * (fL4 + fR4) - 0.5 * lam * (rthetaR - rthetaL)) * nmag;
  }
}

extern "C"
{
  void boundaryflux_eceuler3d_gpu(
      real *fb, real *fextb, real *nhat, real *nscale, real *flux,
      real *p_hyd_b, real *p_hyd_extb,
      real p0, real Rd, real gamma, int N, int nel)
  {
    int ndof = (N + 1) * (N + 1) * 6 * nel;
    int threads_per_block = 256;
    int nblocks_x = ndof / threads_per_block + 1;
    boundaryflux_eceuler3d_kernel<<<dim3(nblocks_x, 1, 1),
        dim3(threads_per_block, 1, 1), 0, 0>>>(
        fb, fextb, nhat, nscale, flux, p_hyd_b, p_hyd_extb, p0, Rd, gamma, N, nel);
  }
}

// ============================================================
// Souza et al. (2023, JAMES) entropy-conservative two-point flux for
// EC Euler 3D with potential temperature on curvilinear mesh, with
// well-balanced hydrostatic pressure split.
//
// Physical flux:
//   f_d(rho)     = <rho>_log * <v_d>
//   f_d(rho*v_i) = <rho>_log * <v_i> * <v_d> + (<p> - <p_hyd>) * delta_{id}
//   f_d(rho*th)  = <rho*theta>_log * <v_d>
//
// where <a>_log is the logarithmic mean, <a> the arithmetic mean.
// The hydrostatic pressure subtraction makes the volume tendency for
// (rho*u, rho*v, rho*w) vanish exactly in the hydrostatic state.
//
// For each computational direction r and node pair, the contravariant
// flux is the parent's split-form projection:
//
//   Fc^r = avg(Ja^r_d) * Fphys_d  (summed over physical dirs d=1,2,3)
//
// Memory layout:
//   s:     SC_3D_INDEX(i,j,k,iel,ivar,N,nel)
//   dsdx:  dsdx[iq + nq*(iel + nel*(d + 3*r))] for Ja^r_d
//   p_hyd: p_hyd[iq + nq*iel]   (Scalar3D nVar=1, interior buffer)
//   f:     TPV_3D_INDEX(nn,i,j,k,iel,ivar,idir,N,nel,nvar)
// ============================================================

template <int blockSize>
__global__ void __launch_bounds__(512) twopointfluxmethod_eceuler3d_kernel(
    real *f, real *s, real *dsdx, real *p_hyd,
    real p0, real Rd, real gamma,
    int nq, int N, int nvar)
{
  uint32_t iq = threadIdx.x;
  if (iq < (uint32_t)nq) {
    uint32_t iel = blockIdx.x;
    uint32_t nel = gridDim.x;
    uint32_t ii  = iq % (N + 1);
    uint32_t jj  = (iq / (N + 1)) % (N + 1);
    uint32_t kk  = iq / (N + 1) / (N + 1);

    // Load state at (i,j,k)
    real rho_ijk   = s[iq + nq * (iel + nel * 0)];
    real rhou_ijk  = s[iq + nq * (iel + nel * 1)];
    real rhov_ijk  = s[iq + nq * (iel + nel * 2)];
    real rhow_ijk  = s[iq + nq * (iel + nel * 3)];
    real rtheta_ijk = s[iq + nq * (iel + nel * 4)];

    real u_ijk     = rhou_ijk / rho_ijk;
    real v_ijk     = rhov_ijk / rho_ijk;
    real w_ijk     = rhow_ijk / rho_ijk;
    real theta_ijk = rtheta_ijk / rho_ijk;
    real p_ijk     = eos_pressure(rho_ijk, theta_ijk, p0, Rd, gamma);
    real ph_ijk    = p_hyd[iq + nq * iel];

    int nq4 = nq * (N + 1); // stride between idir slices in TPV layout

    for (int nn = 0; nn < N + 1; nn++) {

      // ============== xi^1: pair (i,j,k)-(nn,j,k) ==============
      uint32_t iq1 = nn + (N + 1) * (jj + (N + 1) * kk);
      real rho1   = s[iq1 + nq * (iel + nel * 0)];
      real rhou1  = s[iq1 + nq * (iel + nel * 1)];
      real rhov1  = s[iq1 + nq * (iel + nel * 2)];
      real rhow1  = s[iq1 + nq * (iel + nel * 3)];
      real rtheta1 = s[iq1 + nq * (iel + nel * 4)];
      real u1     = rhou1 / rho1;
      real v1     = rhov1 / rho1;
      real w1     = rhow1 / rho1;
      real theta1 = rtheta1 / rho1;
      real p1     = eos_pressure(rho1, theta1, p0, Rd, gamma);
      real ph1    = p_hyd[iq1 + nq * iel];

      real rho_log = log_mean_dev(rho_ijk, rho1);
      real rth_log = log_mean_dev(rtheta_ijk, rtheta1);
      real u_a     = (real)0.5 * (u_ijk + u1);
      real v_a     = (real)0.5 * (v_ijk + v1);
      real w_a     = (real)0.5 * (w_ijk + w1);
      real p_a     = (real)0.5 * (p_ijk + p1) - (real)0.5 * (ph_ijk + ph1);

      for (int ivar = 0; ivar < 5; ivar++) {
        real Fphys_x, Fphys_y, Fphys_z;
        if (ivar == 0) {
          Fphys_x = rho_log * u_a;
          Fphys_y = rho_log * v_a;
          Fphys_z = rho_log * w_a;
        } else if (ivar == 1) {
          Fphys_x = rho_log * u_a * u_a + p_a;
          Fphys_y = rho_log * u_a * v_a;
          Fphys_z = rho_log * u_a * w_a;
        } else if (ivar == 2) {
          Fphys_x = rho_log * v_a * u_a;
          Fphys_y = rho_log * v_a * v_a + p_a;
          Fphys_z = rho_log * v_a * w_a;
        } else if (ivar == 3) {
          Fphys_x = rho_log * w_a * u_a;
          Fphys_y = rho_log * w_a * v_a;
          Fphys_z = rho_log * w_a * w_a + p_a;
        } else {
          Fphys_x = rth_log * u_a;
          Fphys_y = rth_log * v_a;
          Fphys_z = rth_log * w_a;
        }
        real Ja1_x = (real)0.5 * (dsdx[iq + nq * (iel + nel * (0 + 3 * 0))] +
                                   dsdx[iq1 + nq * (iel + nel * (0 + 3 * 0))]);
        real Ja1_y = (real)0.5 * (dsdx[iq + nq * (iel + nel * (1 + 3 * 0))] +
                                   dsdx[iq1 + nq * (iel + nel * (1 + 3 * 0))]);
        real Ja1_z = (real)0.5 * (dsdx[iq + nq * (iel + nel * (2 + 3 * 0))] +
                                   dsdx[iq1 + nq * (iel + nel * (2 + 3 * 0))]);
        f[nn + (N + 1) * iq + nq4 * (iel + nel * (ivar + nvar * 0))] =
            Ja1_x * Fphys_x + Ja1_y * Fphys_y + Ja1_z * Fphys_z;
      }

      // ============== xi^2: pair (i,j,k)-(i,nn,k) ==============
      uint32_t iq2 = ii + (N + 1) * (nn + (N + 1) * kk);
      real rho2   = s[iq2 + nq * (iel + nel * 0)];
      real rhou2  = s[iq2 + nq * (iel + nel * 1)];
      real rhov2  = s[iq2 + nq * (iel + nel * 2)];
      real rhow2  = s[iq2 + nq * (iel + nel * 3)];
      real rtheta2 = s[iq2 + nq * (iel + nel * 4)];
      real u2     = rhou2 / rho2;
      real v2     = rhov2 / rho2;
      real w2     = rhow2 / rho2;
      real theta2 = rtheta2 / rho2;
      real p2     = eos_pressure(rho2, theta2, p0, Rd, gamma);
      real ph2    = p_hyd[iq2 + nq * iel];

      rho_log = log_mean_dev(rho_ijk, rho2);
      rth_log = log_mean_dev(rtheta_ijk, rtheta2);
      u_a     = (real)0.5 * (u_ijk + u2);
      v_a     = (real)0.5 * (v_ijk + v2);
      w_a     = (real)0.5 * (w_ijk + w2);
      p_a     = (real)0.5 * (p_ijk + p2) - (real)0.5 * (ph_ijk + ph2);

      for (int ivar = 0; ivar < 5; ivar++) {
        real Fphys_x, Fphys_y, Fphys_z;
        if (ivar == 0) {
          Fphys_x = rho_log * u_a;
          Fphys_y = rho_log * v_a;
          Fphys_z = rho_log * w_a;
        } else if (ivar == 1) {
          Fphys_x = rho_log * u_a * u_a + p_a;
          Fphys_y = rho_log * u_a * v_a;
          Fphys_z = rho_log * u_a * w_a;
        } else if (ivar == 2) {
          Fphys_x = rho_log * v_a * u_a;
          Fphys_y = rho_log * v_a * v_a + p_a;
          Fphys_z = rho_log * v_a * w_a;
        } else if (ivar == 3) {
          Fphys_x = rho_log * w_a * u_a;
          Fphys_y = rho_log * w_a * v_a;
          Fphys_z = rho_log * w_a * w_a + p_a;
        } else {
          Fphys_x = rth_log * u_a;
          Fphys_y = rth_log * v_a;
          Fphys_z = rth_log * w_a;
        }
        real Ja2_x = (real)0.5 * (dsdx[iq + nq * (iel + nel * (0 + 3 * 1))] +
                                   dsdx[iq2 + nq * (iel + nel * (0 + 3 * 1))]);
        real Ja2_y = (real)0.5 * (dsdx[iq + nq * (iel + nel * (1 + 3 * 1))] +
                                   dsdx[iq2 + nq * (iel + nel * (1 + 3 * 1))]);
        real Ja2_z = (real)0.5 * (dsdx[iq + nq * (iel + nel * (2 + 3 * 1))] +
                                   dsdx[iq2 + nq * (iel + nel * (2 + 3 * 1))]);
        f[nn + (N + 1) * iq + nq4 * (iel + nel * (ivar + nvar * 1))] =
            Ja2_x * Fphys_x + Ja2_y * Fphys_y + Ja2_z * Fphys_z;
      }

      // ============== xi^3: pair (i,j,k)-(i,j,nn) ==============
      uint32_t iq3 = ii + (N + 1) * (jj + (N + 1) * nn);
      real rho3   = s[iq3 + nq * (iel + nel * 0)];
      real rhou3  = s[iq3 + nq * (iel + nel * 1)];
      real rhov3  = s[iq3 + nq * (iel + nel * 2)];
      real rhow3  = s[iq3 + nq * (iel + nel * 3)];
      real rtheta3 = s[iq3 + nq * (iel + nel * 4)];
      real u3     = rhou3 / rho3;
      real v3     = rhov3 / rho3;
      real w3     = rhow3 / rho3;
      real theta3 = rtheta3 / rho3;
      real p3     = eos_pressure(rho3, theta3, p0, Rd, gamma);
      real ph3    = p_hyd[iq3 + nq * iel];

      rho_log = log_mean_dev(rho_ijk, rho3);
      rth_log = log_mean_dev(rtheta_ijk, rtheta3);
      u_a     = (real)0.5 * (u_ijk + u3);
      v_a     = (real)0.5 * (v_ijk + v3);
      w_a     = (real)0.5 * (w_ijk + w3);
      p_a     = (real)0.5 * (p_ijk + p3) - (real)0.5 * (ph_ijk + ph3);

      for (int ivar = 0; ivar < 5; ivar++) {
        real Fphys_x, Fphys_y, Fphys_z;
        if (ivar == 0) {
          Fphys_x = rho_log * u_a;
          Fphys_y = rho_log * v_a;
          Fphys_z = rho_log * w_a;
        } else if (ivar == 1) {
          Fphys_x = rho_log * u_a * u_a + p_a;
          Fphys_y = rho_log * u_a * v_a;
          Fphys_z = rho_log * u_a * w_a;
        } else if (ivar == 2) {
          Fphys_x = rho_log * v_a * u_a;
          Fphys_y = rho_log * v_a * v_a + p_a;
          Fphys_z = rho_log * v_a * w_a;
        } else if (ivar == 3) {
          Fphys_x = rho_log * w_a * u_a;
          Fphys_y = rho_log * w_a * v_a;
          Fphys_z = rho_log * w_a * w_a + p_a;
        } else {
          Fphys_x = rth_log * u_a;
          Fphys_y = rth_log * v_a;
          Fphys_z = rth_log * w_a;
        }
        real Ja3_x = (real)0.5 * (dsdx[iq + nq * (iel + nel * (0 + 3 * 2))] +
                                   dsdx[iq3 + nq * (iel + nel * (0 + 3 * 2))]);
        real Ja3_y = (real)0.5 * (dsdx[iq + nq * (iel + nel * (1 + 3 * 2))] +
                                   dsdx[iq3 + nq * (iel + nel * (1 + 3 * 2))]);
        real Ja3_z = (real)0.5 * (dsdx[iq + nq * (iel + nel * (2 + 3 * 2))] +
                                   dsdx[iq3 + nq * (iel + nel * (2 + 3 * 2))]);
        f[nn + (N + 1) * iq + nq4 * (iel + nel * (ivar + nvar * 2))] =
            Ja3_x * Fphys_x + Ja3_y * Fphys_y + Ja3_z * Fphys_z;
      }
    }
  }
}

extern "C"
{
  void twopointfluxmethod_eceuler3d_gpu(
      real *f, real *s, real *dsdx, real *p_hyd,
      real p0, real Rd, real gamma,
      int N, int nvar, int nel)
  {
    int nq = (N + 1) * (N + 1) * (N + 1);
    if (N < 4) {
      twopointfluxmethod_eceuler3d_kernel<64><<<dim3(nel, 1, 1),
          dim3(64, 1, 1), 0, 0>>>(f, s, dsdx, p_hyd, p0, Rd, gamma, nq, N, nvar);
    } else if (N >= 4 && N < 8) {
      twopointfluxmethod_eceuler3d_kernel<512><<<dim3(nel, 1, 1),
          dim3(512, 1, 1), 0, 0>>>(f, s, dsdx, p_hyd, p0, Rd, gamma, nq, N, nvar);
    }
  }
}

// ============================================================
// Gravitational source term with the well-balanced split:
//   S = [0, 0, 0, -(rho - rho_hyd)*g, 0]
// In the hydrostatic state rho == rho_hyd, so the source vanishes
// exactly (matching the WB volume + boundary fluxes).
// Fully device-resident — reads solution and rho_hyd, writes source.
// ============================================================

__global__ void sourcemethod_eceuler3d_kernel(
    real *source, real *solution, real *rho_hyd, real g, int ndof)
{
  uint32_t idof = threadIdx.x + blockIdx.x * blockDim.x;

  if (idof < (uint32_t)ndof) {
    source[idof]              = 0.0; // rho
    source[idof + ndof]       = 0.0; // rhou
    source[idof + 2 * ndof]   = 0.0; // rhov
    source[idof + 3 * ndof]   = -(solution[idof] - rho_hyd[idof]) * g; // rhow
    source[idof + 4 * ndof]   = 0.0; // rhotheta
  }
}

extern "C"
{
  void sourcemethod_eceuler3d_gpu(
      real *source, real *solution, real *rho_hyd, real g, int N, int nel)
  {
    int ndof = (N + 1) * (N + 1) * (N + 1) * nel;
    int threads_per_block = 256;
    int nblocks_x = ndof / threads_per_block + 1;
    sourcemethod_eceuler3d_kernel<<<dim3(nblocks_x, 1, 1),
        dim3(threads_per_block, 1, 1), 0, 0>>>(source, solution, rho_hyd, g, ndof);
  }
}

// ============================================================
// Diffusive flux for ECEuler3D — constant-coefficient Laplacian
// (Bassi-Rebay 1) with separate momentum (nu) and thermal (kappa)
// diffusivities.
//
// Interior fill:
//   F_d(rho)       = 0
//   F_d(rho*v_i)   = -nu    * d(rho*v_i)/dx_d   for i = 1,2,3
//   F_d(rho*theta) = -kappa * d(rho*theta)/dx_d
//
// solutionGradient layout (Vector3D):
//   dsdx[(idof + ndof_node*ivar) + ndof_per_dir*d]
// where ndof_node = (N+1)^3 * nel  and  ndof_per_dir = ndof_node*nvar.
//
// diffFlux is also a Vector3D with the same layout. Using the
// MappedDGDivergence pipeline downstream gives the BR1 weak-form
// diffusive divergence in a single MappedDGDivergence call.
// ============================================================

__global__ void diffusiveflux_eceuler3d_kernel(
    real *diffFlux, real *grad, real nu, real kappa,
    uint32_t ndof_node, uint32_t nvar)
{
  uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof_total = ndof_node * nvar * 3;
  if (i < ndof_total) {
    // Decode (node, ivar, d). Layout: outer d, middle ivar, inner node.
    uint32_t node = i % ndof_node;
    uint32_t var  = (i / ndof_node) % nvar;
    // Choose coefficient by variable index: 0 -> rho (no diffusion);
    // 1,2,3 -> momentum (nu); 4 -> rho*theta (kappa).
    real coeff;
    if (var == 0) coeff = (real)0;
    else if (var == 4) coeff = kappa;
    else coeff = nu;
    diffFlux[i] = -coeff * grad[i];
    (void)node;
  }
}

extern "C"
{
  void diffusiveflux_eceuler3d_gpu(
      real *diffFlux, real *grad, real nu, real kappa,
      int N, int nvar, int nel)
  {
    uint32_t ndof_node = (N+1) * (N+1) * (N+1) * (uint32_t)nel;
    uint32_t ndof_total = ndof_node * (uint32_t)nvar * 3;
    uint32_t nthreads = 256;
    uint32_t nblocks  = (ndof_total + nthreads - 1) / nthreads;
    diffusiveflux_eceuler3d_kernel<<<dim3(nblocks,1,1), dim3(nthreads,1,1), 0, 0>>>(
        diffFlux, grad, nu, kappa, ndof_node, (uint32_t)nvar);
  }
}

// ============================================================
// Boundary BR1 central flux for the parabolic terms:
//   f_R^diff(iVar) = -coeff(iVar) * (avg_grad(iVar,d) . n_d) * nmag
//
// avgGrad layout (Vector3D, boundary side, with the 3 components
// stored as separate "vars" of the 3*nvar ordering used internally):
//   avgGrad[idof_face + ndof_face*(var + nvar*d)]
// where idof_face = (i + (N+1)*j + (N+1)^2*side + 6*(N+1)^2*iel).
// nhat layout (Vector3D, nVar=1, boundary):
//   nhat[idof_face]                  : x-comp
//   nhat[idof_face + ndof_face]      : y-comp
//   nhat[idof_face + 2*ndof_face]    : z-comp
// nscale layout (Scalar3D, nVar=1):  nscale[idof_face]
// fluxN layout (Scalar3D, boundary): fluxN[idof_face + ndof_face*var]
// ============================================================

__global__ void diffusiveboundaryflux_eceuler3d_kernel(
    real *fluxN, real *avgGrad, real *uBnd, real *uExt,
    real *nhat, real *nscale,
    real nu, real kappa, real tau_nu, real tau_kappa,
    uint32_t ndof_face, uint32_t nvar)
{
  // SIPG-stabilised BR1 boundary flux:
  //   fluxN = -coeff*(avgGrad . n)*nmag + tau*(uL - uR)*nmag
  // tau_nu, tau_kappa precomputed on host:
  //   tau = eta_penalty * coeff * (N+1)^2 / length_scale.
  // Solution buffer layout (Scalar3D, boundary):
  //   uBnd[idof + ndof_face*var]   = solution%boundary(... ,var)
  //   uExt[idof + ndof_face*var]   = solution%extBoundary(... ,var)
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  if (idof < ndof_face) {
    real nx   = nhat[idof];
    real ny   = nhat[idof + ndof_face];
    real nz   = nhat[idof + 2*ndof_face];
    real nmag = nscale[idof];

    for (uint32_t var = 0; var < nvar; var++) {
      real coeff, tau;
      if (var == 0) {
        coeff = (real)0;
        tau   = (real)0;
      } else if (var == 4) {
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
      real gz = avgGrad[idof + ndof_face*(var + nvar*2)];
      real gn = gx*nx + gy*ny + gz*nz;
      real uL = uBnd[idof + ndof_face*var];
      real uR = uExt[idof + ndof_face*var];
      fluxN[idof + ndof_face*var] = (-coeff*gn + tau*(uL - uR)) * nmag;
    }
  }
}

extern "C"
{
  void diffusiveboundaryflux_eceuler3d_gpu(
      real *fluxN, real *avgGrad, real *uBnd, real *uExt,
      real *nhat, real *nscale,
      real nu, real kappa, real tau_nu, real tau_kappa,
      int N, int nvar, int nel)
  {
    uint32_t ndof_face = (N+1) * (N+1) * 6 * (uint32_t)nel;
    uint32_t nthreads = 256;
    uint32_t nblocks  = (ndof_face + nthreads - 1) / nthreads;
    diffusiveboundaryflux_eceuler3d_kernel<<<dim3(nblocks,1,1), dim3(nthreads,1,1), 0, 0>>>(
        fluxN, avgGrad, uBnd, uExt, nhat, nscale,
        nu, kappa, tau_nu, tau_kappa, ndof_face, (uint32_t)nvar);
  }
}
