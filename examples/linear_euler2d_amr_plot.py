# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
#
# Maintainers : support@fluidnumerics.com
# Official Repository : https://github.com/FluidNumerics/self/
#
# Copyright © 2024 Fluid Numerics LLC
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
"""Render 2-D AMR model output: pressure field + mesh skeleton.

Reads the HDF5 snapshots written by SELF 2-D models (e.g. the
linear_euler2d_amr_ultrasound_pointsource example), renders the pressure field with the
element skeleton overlaid, and writes one PNG per snapshot (plus an MP4 movie when ffmpeg is
available). Because every snapshot carries its own geometry, snapshots with different element
counts - an adapting mesh - need no special handling: the wireframe IS the refinement pattern.

Only h5py, numpy, and matplotlib are required (no pyvista / GL). Each element's field and
geometry are interpolated from the Gauss control points to a uniform per-element grid that
includes the element edges (barycentric Lagrange interpolation, exact for the polynomial
data), so both the filled field and the element outlines are faithful to the high-order
representation.

Usage:
    python linear_euler2d_amr_plot.py <directory-with-solution.*.h5> [options]

Options:
    --var NAME        solution variable to render (default: P)
    --nplot M         per-element plot resolution, including edges (default: 8)
    --vmax V          symmetric color scale; default: 0.5 * max|field| of the first snapshot
    --frame-time DT   seconds of simulation time per snapshot, used for frame titles
    --fps FPS         movie frame rate (default: 10)
    --outdir DIR      output directory for frames/movie (default: <input directory>/plots)
    --no-movie        skip the ffmpeg movie assembly
"""

import argparse
import glob
import os
import subprocess
import sys

import h5py
import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection


def lagrange_matrix(control_points, target_points):
    """Barycentric Lagrange interpolation matrix L: L @ f(control) = f(target).

    Exact for polynomials of degree <= len(control_points)-1. Rows for target points that
    coincide with a control point reduce to a Kronecker delta (handled explicitly to avoid
    division by zero).
    """
    xc = np.asarray(control_points, dtype=np.float64)
    xt = np.asarray(target_points, dtype=np.float64)
    n = xc.size
    # Barycentric weights
    w = np.ones(n)
    for j in range(n):
        w[j] = 1.0 / np.prod(xc[j] - np.delete(xc, j))
    L = np.zeros((xt.size, n))
    for i, x in enumerate(xt):
        d = x - xc
        hit = np.where(np.abs(d) < 1e-14)[0]
        if hit.size > 0:
            L[i, hit[0]] = 1.0
        else:
            t = w / d
            L[i, :] = t / np.sum(t)
    return L


def load_snapshot(fname, varname):
    """Load control points, per-element geometry, and one solution variable.

    Returns (xi, x, y, f) with x, y, f of shape (nelem, Np, Np). Fortran writes (i,j,iel), so
    h5py sees (iel, j, i); the tensor interpolation below is symmetric in the two node
    indices, so no transposition is needed as long as x, y, and f are treated identically.
    """
    with h5py.File(fname, "r") as f:
        xi = f["interp/controlpoints"][()]
        x = f["controlgrid/geometry/x_dim1"][()]
        y = f["controlgrid/geometry/x_dim2"][()]
        u = f[f"controlgrid/solution/{varname}"][()]
    return xi, x, y, u


def upsample(L, arr):
    """Tensor-product interpolation of (nelem, Np, Np) data to (nelem, m, m)."""
    # (m,Np) @ (nelem,Np,Np) @ (Np,m) -> (nelem,m,m), applied to both node indices
    return L @ arr @ L.T


def element_triangulation(xu, yu):
    """One global Triangulation over all elements' upsampled nodes (2 triangles per cell)."""
    nelem, m, _ = xu.shape
    pts_per_el = m * m
    tris = []
    for e in range(nelem):
        base = e * pts_per_el
        for j in range(m - 1):
            for i in range(m - 1):
                n0 = base + j * m + i
                n1 = base + j * m + i + 1
                n2 = base + (j + 1) * m + i + 1
                n3 = base + (j + 1) * m + i
                tris.append((n0, n1, n2))
                tris.append((n0, n2, n3))
    return mtri.Triangulation(xu.ravel(), yu.ravel(), np.asarray(tris))


def element_outlines(xu, yu):
    """Element-boundary polylines (the mesh skeleton) from the upsampled geometry."""
    segs = []
    for e in range(xu.shape[0]):
        xe, ye = xu[e], yu[e]
        for xs, ys in (
            (xe[0, :], ye[0, :]),
            (xe[-1, :], ye[-1, :]),
            (xe[:, 0], ye[:, 0]),
            (xe[:, -1], ye[:, -1]),
        ):
            segs.append(np.column_stack([xs, ys]))
    return segs


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("directory", help="directory containing solution.*.h5 snapshots")
    ap.add_argument("--var", default="P", help="solution variable to render (default: P)")
    ap.add_argument("--nplot", type=int, default=8,
                    help="per-element plot resolution incl. edges (default: 8)")
    ap.add_argument("--vmax", type=float, default=None,
                    help="symmetric color scale (default: 0.8*max|field| of the last "
                         "snapshot, so the decayed propagating wave stays visible; earlier, "
                         "higher-amplitude frames saturate)")
    ap.add_argument("--frame-time", type=float, default=None,
                    help="simulation seconds per snapshot, for frame titles")
    ap.add_argument("--fps", type=int, default=10, help="movie frame rate (default: 10)")
    ap.add_argument("--outdir", default=None,
                    help="output directory (default: <directory>/plots)")
    ap.add_argument("--no-movie", action="store_true", help="skip ffmpeg movie assembly")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.directory, "solution.*.h5")))
    if not files:
        print(f"No solution.*.h5 files found in {args.directory}")
        return 1

    outdir = args.outdir or os.path.join(args.directory, "plots")
    os.makedirs(outdir, exist_ok=True)

    m = args.nplot
    target = np.linspace(-1.0, 1.0, m)

    vmax = args.vmax
    if vmax is None:
        # Fix the scale across frames to the (decayed) final-snapshot amplitude so the
        # propagating wave stays visible; the initial pulse saturates, which reads naturally.
        _, _, _, ulast = load_snapshot(files[-1], args.var)
        vmax = 0.8 * float(np.max(np.abs(ulast)))
        if vmax == 0.0:
            vmax = 1.0

    for k, fname in enumerate(files):
        xi, x, y, u = load_snapshot(fname, args.var)
        L = lagrange_matrix(xi, target)
        xu, yu, uu = upsample(L, x), upsample(L, y), upsample(L, u)

        fig, ax = plt.subplots(figsize=(8, 7), dpi=130)
        tri = element_triangulation(xu, yu)
        tpc = ax.tripcolor(tri, uu.ravel(), shading="gouraud",
                           cmap="RdBu_r", vmin=-vmax, vmax=vmax, rasterized=True)
        ax.add_collection(LineCollection(element_outlines(xu, yu),
                                         colors="black", linewidths=0.4, alpha=0.5))
        ax.set_aspect("equal")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        nelem = x.shape[0]
        title = f"{args.var}  |  frame {k}  |  {nelem} elements"
        if args.frame_time is not None:
            title += f"  |  t = {k*args.frame_time*1e6:.1f} µs"
        ax.set_title(title)
        fig.colorbar(tpc, ax=ax, label=args.var)
        frame = os.path.join(outdir, f"amr_{args.var}_{k:05d}.png")
        fig.savefig(frame, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {frame}  ({nelem} elements)")

    if not args.no_movie and len(files) > 1:
        movie = os.path.join(outdir, f"amr_{args.var}.mp4")
        cmd = ["ffmpeg", "-y", "-framerate", str(args.fps),
               "-i", os.path.join(outdir, f"amr_{args.var}_%05d.png"),
               "-pix_fmt", "yuv420p",
               "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
               movie]
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            print(f"wrote {movie}")
        except (FileNotFoundError, subprocess.CalledProcessError) as err:
            print(f"movie assembly skipped ({err}); frames are in {outdir}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
