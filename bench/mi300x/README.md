# MI300X (gfx942) benchmark + profiling harness

Single-GPU LinearEuler3D forward-step benchmarking and GPU profiling for the SELF
3-D DG divergence kernel chain, on the `jschoonover@iad-mj-login` cluster
(8×MI300X, slurm + enroot/pyxis). These are the ROCm counterparts of the
NVIDIA `ncu` scripts in the repo root.

## What gets measured

The benchmark binary is `examples/bench_lineareuler3d_scaling`
(`[nex ney nez ntx nty ntz nsteps cdegree]`; single-GPU = `1 1 1` tiles). Its
RK3 forward step exercises the divergence chain
`ContravariantProjection_3D_gpukernel → VectorDivergence_3D_gpukernel →
DG_BoundaryContribution_3D_gpukernel → JacobianWeight` plus boundary interp,
side exchange, flux, and DSDt kernels. (Gradient kernels are not exercised by
this model — `gradient_enabled=.false.`.)

Profile at **N=7** for the baseline. After the N≥8 kernel fix is validated on
the image, extend to N=8–12 (pass `CDEGREE=…`).

## Tool mapping (vs the NVIDIA `ncu` scripts)

| Purpose | NVIDIA | MI300X (here) |
|---|---|---|
| per-kernel time (hotspots) | `ncu --csv --metrics gpu__time_duration.sum` | `rocprofv3 --kernel-trace` |
| HW events (occupancy/VALU/LDS/L2) | `ncu --set full` | `rocprofv3 --pmc <counters>` |
| kernel filter | `-k regex:…` | `--kernel-include-regex` |
| device pin | `NVIDIA_VISIBLE_DEVICES` | `ROCR_VISIBLE_DEVICES=0` |

**Profiler install.** The base image ships ROCm 6.4.3 but NO profiler CLI.
`install_profiler.slurm` adds `rocprofv3` from the **`rocprofiler-sdk`** package
via dnf from repo.radeon.com (Rocky 9), producing `self-gfx942-dp-prof.sqsh`.
(Do not confuse `rocprofiler-sdk` — which provides `rocprofv3` — with the
deprecated `rocprofiler` v2 package; and `rocprof-compute`/omniperf installs but
currently python-errors on launch, so we use `rocprofv3` directly.)

**Critical gotcha:** `rocprofv3` writes a `.rocprofv3` scratch dir in the current
working directory, so it MUST run from a writable path — the harness uses
`--container-workdir=/results` (the bind-mounted output dir). Running from the
read-only `/opt/self/build` makes the app abort (SIGABRT). rocprofv3 collects a
few counters per pass, so `--pmc` lists are kept short (split across passes).
rocprofv3 CSV kernel names contain quoted commas — aggregate with a quote-aware
parser (gawk `FPAT`), not plain `-F,`.

## Scripts

| Script | Role |
|---|---|
| `sync_source.sh` | rsync the working tree to `~/self-src` on the cluster (run locally) |
| `build_gfx942_sqsh.slurm` | in-container rebuild → `self-gfx942-dp.sqsh` (`DP=OFF` → fp32) |
| `install_profiler.slurm` | add `rocprofv3` (rocprofiler-sdk) → `self-gfx942-dp-prof.sqsh` |
| `validate_ctest.slurm` | full serial ctest (126/126 dp) + the new N≥8 regressions |
| `hotspot_sweep.slurm` | per-kernel time over 4³/8³/16³/24³ (rocprofv3 --kernel-trace) |
| `hwevents_divergence.slurm` | occupancy/VALU/LDS/L2 per chain kernel (rocprofv3 --pmc) |
| `rocprof_profile.slurm` | legacy fallback: rocprof v1 `--hip-trace --stats` hotspots |
| `ab_compare.slurm` | baseline vs modified image: per-kernel deltas + DOF/s |

## Iterate loop

1. `bench/mi300x/sync_source.sh` — ship source to the cluster.
2. `sbatch bench/mi300x/build_gfx942_sqsh.slurm` — build `self-gfx942-dp.sqsh`.
3. `sbatch bench/mi300x/install_profiler.slurm` — add rocprofv3 → `self-gfx942-dp-prof.sqsh` (once per image).
4. `sbatch bench/mi300x/validate_ctest.slurm` — confirm correctness.
5. `sbatch bench/mi300x/hotspot_sweep.slurm` — rank kernels; find the top one.
6. `sbatch bench/mi300x/hwevents_divergence.slurm` — bottleneck class of that kernel.
6. Freeze the baseline: `cp ~/self-gfx942-dp.sqsh ~/self-gfx942-dp-baseline.sqsh`.
7. Make a code change locally → repeat 1–2 (rebuild into `self-gfx942-dp.sqsh`).
8. `sbatch bench/mi300x/ab_compare.slurm` — per-kernel Δns/% and DOF/s vs baseline.
9. Record in a running ledger (kernel, size, precision, mean ns, SOL%, occupancy,
   L2 hit, HBM GB/s, DOF/s); loop.

Then build `self-gfx942-fp32.sqsh` (`DP=OFF`) and A/B double-vs-fp32 through the
same harness (fp32 changes results — report it as a lever alongside, not
replacing, the validated double baseline).

All jobs mount `~/le3d-results/<tagged-dir>` into the container as `/results`, so
outputs persist on the login node.

## Cluster notes

- slurm binaries at `/opt/slurm-24.05.7/bin` (may not be on PATH).
- Base image `~/self-base-gfx942-dp.sqsh` has the ROCm toolchain + `/opt/self`.
- `node006` has broken container GPU access → all scripts `--exclude=node006`.
- Container OpenMPI is `--without-slurm`; single-rank runs strip
  `PMIX_*/PMI_*/SLURM_*/OMPI_*/PRTE_*` and run as a singleton (done in-script).
- ctest needs `--container-writable` and `WORKSPACE=/opt/self/src`.
