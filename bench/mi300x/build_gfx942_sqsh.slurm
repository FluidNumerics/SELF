#!/bin/bash
# Build/refresh a gfx942 (MI300X) SELF enroot image from rsync'd source.
#
# Mirrors the documented ~70s in-container rebuild: start from the base image
# ~/self-base-gfx942-dp.sqsh (full ROCm toolchain + OpenMPI 5 --with-rocm +
# /opt/self populated), copy the freshly-synced source over /opt/self/src, run
# cmake+make+install, and save the result with --container-save.
#
# Submit AFTER bench/mi300x/sync_source.sh has populated ~/self-src.
#
# Environment overrides (export before sbatch, or edit defaults):
#   DP        ON  -> double precision (self-gfx942-dp.sqsh, validated)
#             OFF -> single precision (self-gfx942-fp32.sqsh, tuning lever)
#   BASE      base image                 (default ~/self-base-gfx942-dp.sqsh)
#   SRC       synced source dir          (default ~/self-src)
#   OUT_SQSH  output image               (default derived from DP)
#SBATCH -J self-gfx942-build
#SBATCH -p main
#SBATCH -n1
#SBATCH --cpus-per-task=48
#SBATCH --mem=220G
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --exclude=node006
#SBATCH -o self-gfx942-build-%j.out
#SBATCH -e self-gfx942-build-%j.err
set -e

export ENROOT_TEMP_PATH="$HOME/.enroot/tmp" ENROOT_CACHE_PATH="$HOME/.enroot/cache" ENROOT_DATA_PATH="$HOME/.enroot/data" ENROOT_RUNTIME_PATH="$HOME/.enroot/run"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_CACHE_PATH" "$ENROOT_DATA_PATH" "$ENROOT_RUNTIME_PATH"
export MELLANOX_VISIBLE_DEVICES=none

DP="${DP:-ON}"
BASE="${BASE:-$HOME/self-base-gfx942-dp.sqsh}"
SRC="${SRC:-$HOME/self-src}"
if [ "$DP" = "ON" ]; then OUT_SQSH="${OUT_SQSH:-$HOME/self-gfx942-dp.sqsh}"; else OUT_SQSH="${OUT_SQSH:-$HOME/self-gfx942-fp32.sqsh}"; fi

echo "BASE=$BASE  SRC=$SRC  DP=$DP  OUT_SQSH=$OUT_SQSH"

srun --gres=gpu:1 --mpi=none --exclude=node006 \
     --container-writable --container-image="$BASE" \
     --container-mounts="${SRC}:/mnt/src" --container-save="$OUT_SQSH" \
  bash -lc '
    set -e
    source /opt/spack-environment/activate.sh 2>/dev/null || true
    rocminfo | grep -m1 gfx || true
    # rsync is not present inside the container image; use cp to fully refresh
    # /opt/self/src from the mounted host source (already pruned of build/.git).
    rm -rf /opt/self/src
    mkdir -p /opt/self/src
    cp -a /mnt/src/. /opt/self/src/
    cmake -S /opt/self/src -B /opt/self/build \
      -DSELF_ENABLE_GPU=ON -DSELF_GPU_BACKEND=HIP \
      -DCMAKE_HIP_ARCHITECTURES=gfx942 -DGPU_TARGETS=gfx942 \
      -DSELF_ENABLE_EXAMPLES=ON -DSELF_ENABLE_TESTING=ON \
      -DSELF_ENABLE_DOUBLE_PRECISION='"$DP"'
    cmake --build /opt/self/build -j48
    cmake --install /opt/self/build
    # Ensure the ROCm compute profiler is present (this is the only writable step).
    command -v rocprofv3 || echo "WARNING: rocprofv3 not found"
    command -v rocprof-compute >/dev/null 2>&1 || command -v omniperf >/dev/null 2>&1 || \
      pip install --quiet rocprofiler-compute 2>/dev/null || echo "NOTE: rocprof-compute unavailable; hwevents script will use rocprofv3 --pmc fallback"
    echo "build complete"
  '
echo "Saved image: $OUT_SQSH"
