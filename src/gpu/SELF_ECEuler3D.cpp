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

    // Normal fluxes
    real fL0 = rhoL * unL;
    real fR0 = rhoR * unR;

    real fL1 = rhouL * unL + pL * nx;
    real fR1 = rhouR * unR + pR * nx;

    real fL2 = rhovL * unL + pL * ny;
    real fR2 = rhovR * unR + pR * ny;

    real fL3 = rhowL * unL + pL * nz;
    real fR3 = rhowR * unR + pR * nz;

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
      real p0, real Rd, real gamma, int N, int nel)
  {
    int ndof = (N + 1) * (N + 1) * 6 * nel;
    int threads_per_block = 256;
    int nblocks_x = ndof / threads_per_block + 1;
    boundaryflux_eceuler3d_kernel<<<dim3(nblocks_x, 1, 1),
        dim3(threads_per_block, 1, 1), 0, 0>>>(
        fb, fextb, nhat, nscale, flux, p0, Rd, gamma, N, nel);
  }
}

// ============================================================
// Kennedy-Gruber two-point flux for EC Euler 3D on curvilinear mesh
//
// For each computational direction r and node pair, computes:
//   Fc^r = avg(Ja^r_d) * Fphys_d  (summed over physical dirs d=1,2,3)
//
// Physical flux (Kennedy-Gruber split form):
//   f_d(rho)    = {{rho}} * {{v_d}}
//   f_d(rho*vi) = {{rho}} * {{vi}} * {{v_d}} + {{p}} * delta_{id}
//   f_d(rho*th) = {{rho}} * {{theta}} * {{v_d}}
//
// Memory layout:
//   s:    SC_3D_INDEX(i,j,k,iel,ivar,N,nel)
//   dsdx: stored as dsdx[iq + nq*(iel + nel*(d + 3*r))] for Ja^r_d
//   f:    TPV_3D_INDEX(nn,i,j,k,iel,ivar,idir,N,nel,nvar)
// ============================================================

template <int blockSize>
__global__ void __launch_bounds__(512) twopointfluxmethod_eceuler3d_kernel(
    real *f, real *s, real *dsdx,
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

    int nq4 = nq * (N + 1); // stride between idir slices in TPV layout

    for (int nn = 0; nn < N + 1; nn++) {

      // === xi^1: pair (i,j,k)-(nn,j,k) ===
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

      real rho_a = 0.5 * (rho_ijk + rho1);
      real u_a   = 0.5 * (u_ijk + u1);
      real v_a   = 0.5 * (v_ijk + v1);
      real w_a   = 0.5 * (w_ijk + w1);
      real th_a  = 0.5 * (theta_ijk + theta1);
      real p_a   = 0.5 * (p_ijk + p1);

      // Physical flux dotted with averaged Ja^1
      // dsdx layout: dsdx[iq + nq*(iel + nel*(d + 3*r))]
      // For r=0 (xi^1): Ja^1_d at d=0,1,2
      for (int ivar = 0; ivar < 5; ivar++) {
        real Fphys_x, Fphys_y, Fphys_z;
        if (ivar == 0) {
          Fphys_x = rho_a * u_a;
          Fphys_y = rho_a * v_a;
          Fphys_z = rho_a * w_a;
        } else if (ivar == 1) {
          Fphys_x = rho_a * u_a * u_a + p_a;
          Fphys_y = rho_a * u_a * v_a;
          Fphys_z = rho_a * u_a * w_a;
        } else if (ivar == 2) {
          Fphys_x = rho_a * v_a * u_a;
          Fphys_y = rho_a * v_a * v_a + p_a;
          Fphys_z = rho_a * v_a * w_a;
        } else if (ivar == 3) {
          Fphys_x = rho_a * w_a * u_a;
          Fphys_y = rho_a * w_a * v_a;
          Fphys_z = rho_a * w_a * w_a + p_a;
        } else {
          Fphys_x = rho_a * th_a * u_a;
          Fphys_y = rho_a * th_a * v_a;
          Fphys_z = rho_a * th_a * w_a;
        }
        real Ja1_x = 0.5 * (dsdx[iq + nq * (iel + nel * (0 + 3 * 0))] +
                             dsdx[iq1 + nq * (iel + nel * (0 + 3 * 0))]);
        real Ja1_y = 0.5 * (dsdx[iq + nq * (iel + nel * (1 + 3 * 0))] +
                             dsdx[iq1 + nq * (iel + nel * (1 + 3 * 0))]);
        real Ja1_z = 0.5 * (dsdx[iq + nq * (iel + nel * (2 + 3 * 0))] +
                             dsdx[iq1 + nq * (iel + nel * (2 + 3 * 0))]);
        f[nn + (N + 1) * iq + nq4 * (iel + nel * (ivar + nvar * 0))] =
            Ja1_x * Fphys_x + Ja1_y * Fphys_y + Ja1_z * Fphys_z;
      }

      // === xi^2: pair (i,j,k)-(i,nn,k) ===
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

      rho_a = 0.5 * (rho_ijk + rho2);
      u_a   = 0.5 * (u_ijk + u2);
      v_a   = 0.5 * (v_ijk + v2);
      w_a   = 0.5 * (w_ijk + w2);
      th_a  = 0.5 * (theta_ijk + theta2);
      p_a   = 0.5 * (p_ijk + p2);

      for (int ivar = 0; ivar < 5; ivar++) {
        real Fphys_x, Fphys_y, Fphys_z;
        if (ivar == 0) {
          Fphys_x = rho_a * u_a;
          Fphys_y = rho_a * v_a;
          Fphys_z = rho_a * w_a;
        } else if (ivar == 1) {
          Fphys_x = rho_a * u_a * u_a + p_a;
          Fphys_y = rho_a * u_a * v_a;
          Fphys_z = rho_a * u_a * w_a;
        } else if (ivar == 2) {
          Fphys_x = rho_a * v_a * u_a;
          Fphys_y = rho_a * v_a * v_a + p_a;
          Fphys_z = rho_a * v_a * w_a;
        } else if (ivar == 3) {
          Fphys_x = rho_a * w_a * u_a;
          Fphys_y = rho_a * w_a * v_a;
          Fphys_z = rho_a * w_a * w_a + p_a;
        } else {
          Fphys_x = rho_a * th_a * u_a;
          Fphys_y = rho_a * th_a * v_a;
          Fphys_z = rho_a * th_a * w_a;
        }
        real Ja2_x = 0.5 * (dsdx[iq + nq * (iel + nel * (0 + 3 * 1))] +
                             dsdx[iq2 + nq * (iel + nel * (0 + 3 * 1))]);
        real Ja2_y = 0.5 * (dsdx[iq + nq * (iel + nel * (1 + 3 * 1))] +
                             dsdx[iq2 + nq * (iel + nel * (1 + 3 * 1))]);
        real Ja2_z = 0.5 * (dsdx[iq + nq * (iel + nel * (2 + 3 * 1))] +
                             dsdx[iq2 + nq * (iel + nel * (2 + 3 * 1))]);
        f[nn + (N + 1) * iq + nq4 * (iel + nel * (ivar + nvar * 1))] =
            Ja2_x * Fphys_x + Ja2_y * Fphys_y + Ja2_z * Fphys_z;
      }

      // === xi^3: pair (i,j,k)-(i,j,nn) ===
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

      rho_a = 0.5 * (rho_ijk + rho3);
      u_a   = 0.5 * (u_ijk + u3);
      v_a   = 0.5 * (v_ijk + v3);
      w_a   = 0.5 * (w_ijk + w3);
      th_a  = 0.5 * (theta_ijk + theta3);
      p_a   = 0.5 * (p_ijk + p3);

      for (int ivar = 0; ivar < 5; ivar++) {
        real Fphys_x, Fphys_y, Fphys_z;
        if (ivar == 0) {
          Fphys_x = rho_a * u_a;
          Fphys_y = rho_a * v_a;
          Fphys_z = rho_a * w_a;
        } else if (ivar == 1) {
          Fphys_x = rho_a * u_a * u_a + p_a;
          Fphys_y = rho_a * u_a * v_a;
          Fphys_z = rho_a * u_a * w_a;
        } else if (ivar == 2) {
          Fphys_x = rho_a * v_a * u_a;
          Fphys_y = rho_a * v_a * v_a + p_a;
          Fphys_z = rho_a * v_a * w_a;
        } else if (ivar == 3) {
          Fphys_x = rho_a * w_a * u_a;
          Fphys_y = rho_a * w_a * v_a;
          Fphys_z = rho_a * w_a * w_a + p_a;
        } else {
          Fphys_x = rho_a * th_a * u_a;
          Fphys_y = rho_a * th_a * v_a;
          Fphys_z = rho_a * th_a * w_a;
        }
        real Ja3_x = 0.5 * (dsdx[iq + nq * (iel + nel * (0 + 3 * 2))] +
                             dsdx[iq3 + nq * (iel + nel * (0 + 3 * 2))]);
        real Ja3_y = 0.5 * (dsdx[iq + nq * (iel + nel * (1 + 3 * 2))] +
                             dsdx[iq3 + nq * (iel + nel * (1 + 3 * 2))]);
        real Ja3_z = 0.5 * (dsdx[iq + nq * (iel + nel * (2 + 3 * 2))] +
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
      real *f, real *s, real *dsdx,
      real p0, real Rd, real gamma,
      int N, int nvar, int nel)
  {
    int nq = (N + 1) * (N + 1) * (N + 1);
    if (N < 4) {
      twopointfluxmethod_eceuler3d_kernel<64><<<dim3(nel, 1, 1),
          dim3(64, 1, 1), 0, 0>>>(f, s, dsdx, p0, Rd, gamma, nq, N, nvar);
    } else if (N >= 4 && N < 8) {
      twopointfluxmethod_eceuler3d_kernel<512><<<dim3(nel, 1, 1),
          dim3(512, 1, 1), 0, 0>>>(f, s, dsdx, p0, Rd, gamma, nq, N, nvar);
    }
  }
}

// ============================================================
// Gravitational source term: S = [0, 0, 0, -rho*g, 0]
// Fully device-resident — reads solution, writes source on GPU.
// ============================================================

__global__ void sourcemethod_eceuler3d_kernel(
    real *source, real *solution, real g, int ndof)
{
  uint32_t idof = threadIdx.x + blockIdx.x * blockDim.x;

  if (idof < (uint32_t)ndof) {
    source[idof]              = 0.0; // rho
    source[idof + ndof]       = 0.0; // rhou
    source[idof + 2 * ndof]   = 0.0; // rhov
    source[idof + 3 * ndof]   = -solution[idof] * g; // rhow = -rho*g
    source[idof + 4 * ndof]   = 0.0; // rhotheta
  }
}

extern "C"
{
  void sourcemethod_eceuler3d_gpu(
      real *source, real *solution, real g, int N, int nel)
  {
    int ndof = (N + 1) * (N + 1) * (N + 1) * nel;
    int threads_per_block = 256;
    int nblocks_x = ndof / threads_per_block + 1;
    sourcemethod_eceuler3d_kernel<<<dim3(nblocks_x, 1, 1),
        dim3(threads_per_block, 1, 1), 0, 0>>>(source, solution, g, ndof);
  }
}
