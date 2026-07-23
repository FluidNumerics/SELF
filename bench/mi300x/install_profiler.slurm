#!/bin/bash
# Add the ROCm profiler (rocprofv3 from rocprofiler-sdk, and rocprof-compute if
# available) to the gfx942 SELF image. The base image ships ROCm 6.4.3 under the
# official /opt/rocm-6.4.3 layout but WITHOUT the profiler CLIs, so we install
# them from the AMD RHEL9 repo (the container is Rocky 9 with internet access)
# into a writable copy and save a new image.
#
# In:  IN_SQSH  (default ~/self-gfx942-dp.sqsh)
# Out: OUT_SQSH (default ~/self-gfx942-dp-prof.sqsh)
#SBATCH -J self-gfx942-prof
#SBATCH -p main
#SBATCH -n1
#SBATCH --cpus-per-task=16
#SBATCH --mem=220G
#SBATCH --gres=gpu:1
#SBATCH --time=00:40:00
#SBATCH --exclude=node006
#SBATCH -o self-gfx942-prof-%j.out
#SBATCH -e self-gfx942-prof-%j.err
set -e

export ENROOT_TEMP_PATH="$HOME/.enroot/tmp" ENROOT_CACHE_PATH="$HOME/.enroot/cache" ENROOT_DATA_PATH="$HOME/.enroot/data" ENROOT_RUNTIME_PATH="$HOME/.enroot/run"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_CACHE_PATH" "$ENROOT_DATA_PATH" "$ENROOT_RUNTIME_PATH"
export MELLANOX_VISIBLE_DEVICES=none

IN_SQSH="${IN_SQSH:-$HOME/self-gfx942-dp.sqsh}"
OUT_SQSH="${OUT_SQSH:-$HOME/self-gfx942-dp-prof.sqsh}"
ROCM_VER="${ROCM_VER:-6.4.3}"
echo "IN=$IN_SQSH  OUT=$OUT_SQSH  ROCM_VER=$ROCM_VER"

# --container-remap-root: dnf must run as root; enroot otherwise maps the calling
# (unprivileged) user into the container.
srun --gres=gpu:1 --mpi=none --exclude=node006 \
     --container-writable --container-remap-root \
     --container-image="$IN_SQSH" --container-save="$OUT_SQSH" \
  bash -lc '
    set -e
    cat > /etc/yum.repos.d/rocm.repo <<REPO
[ROCm-'"$ROCM_VER"']
name=ROCm'"$ROCM_VER"'
baseurl=https://repo.radeon.com/rocm/rhel9/'"$ROCM_VER"'/main
enabled=1
priority=50
gpgcheck=1
gpgkey=https://repo.radeon.com/rocm/rocm.gpg.key
REPO
    dnf clean all
    # rocprofv3 CLI is provided by the rocprofiler-sdk package (verified via
    # dnf provides */bin/rocprofv3). NOTE: do NOT pair it with a non-existent
    # package name in one transaction or the whole install aborts.
    dnf install -y --nogpgcheck rocprofiler-sdk rocprofiler-sdk-roctx 2>&1 | tail -25 || true
    # rocprof-compute (omniperf) is optional and may fail to run (python deps).
    dnf install -y --nogpgcheck rocprofiler-compute 2>&1 | tail -8 || true
    export PATH=/opt/rocm/bin:/opt/rocm-'"$ROCM_VER"'/bin:$PATH
    echo "=== verify ==="
    command -v rocprofv3 && rocprofv3 --version 2>&1 | head -3 || echo "rocprofv3 MISSING"
    echo "--- counter availability (head) ---"
    rocprofv3 --list-avail 2>&1 | head -15 || rocprofv3 -L 2>&1 | head -15 || true
    command -v rocprof-compute && rocprof-compute --version 2>&1 | head -2 || echo "rocprof-compute MISSING"
  '
echo "Saved image: $OUT_SQSH"
