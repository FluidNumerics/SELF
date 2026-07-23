#!/bin/bash
# A/B comparison of two SELF gfx942 images on the single-GPU LinearEuler3D
# forward step. Because each code change is rebuilt into a new .sqsh, the two
# "variants" are two images (baseline vs modified). For each image we capture
# the headline DOF/s and a per-kernel hotspot trace; then join offline for
# per-kernel delta ns / %.
#
# This also serves double-vs-fp32 A/B: point BASE_SQSH/MOD_SQSH at
# self-gfx942-dp.sqsh and self-gfx942-fp32.sqsh.
#
# Env overrides:
#   BASE_SQSH  frozen baseline image (default ~/self-gfx942-dp-baseline.sqsh)
#   MOD_SQSH   modified image        (default ~/self-gfx942-dp.sqsh)
#   OUT        host results dir       (default ~/le3d-results/ab_<timestamp>)
#   NSTEPS, CDEGREE
#SBATCH -J self-le3d-ab
#SBATCH -p main
#SBATCH -n1
#SBATCH --cpus-per-task=16
#SBATCH --mem=220G
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --exclude=node006
#SBATCH -o self-le3d-ab-%j.out
#SBATCH -e self-le3d-ab-%j.err
set -e

export ENROOT_TEMP_PATH="$HOME/.enroot/tmp" ENROOT_CACHE_PATH="$HOME/.enroot/cache" ENROOT_DATA_PATH="$HOME/.enroot/data" ENROOT_RUNTIME_PATH="$HOME/.enroot/run"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_CACHE_PATH" "$ENROOT_DATA_PATH" "$ENROOT_RUNTIME_PATH"
export MELLANOX_VISIBLE_DEVICES=none

BASE_SQSH="${BASE_SQSH:-$HOME/self-gfx942-dp-baseline.sqsh}"
MOD_SQSH="${MOD_SQSH:-$HOME/self-gfx942-dp.sqsh}"
OUT="${OUT:-$HOME/le3d-results/ab_$(date +%Y%m%d_%H%M%S)}"; mkdir -p "$OUT"
NSTEPS="${NSTEPS:-20}"
CDEGREE="${CDEGREE:-7}"

echo "BASE=$BASE_SQSH  MOD=$MOD_SQSH  OUT=$OUT  N=$CDEGREE"

run_variant () {
  local V="$1" SQSH="$2"
  echo "########## variant=$V  image=$SQSH ##########"
  srun --gres=gpu:1 --mpi=none --exclude=node006 --container-image="$SQSH" \
       --container-mounts="${OUT}:/results" --container-workdir=/opt/self/build \
    bash -lc '
      source /opt/spack-environment/activate.sh 2>/dev/null || true
      export LD_LIBRARY_PATH=/opt/self/build:/opt/self/install/lib:${LD_LIBRARY_PATH}
      export WORKSPACE=/opt/self/src
      export ROCR_VISIBLE_DEVICES=0
      for v in $(env | grep -E "^(PMIX_|PMI_|SLURM_|OMPI_|PRTE_)" | cut -d= -f1); do unset "$v"; done
      export OMPI_MCA_btl=self,vader OMPI_MCA_pml=ob1
      BIN=/opt/self/build/examples/bench_lineareuler3d_scaling
      for S in "8 8 8" "16 16 16"; do
        TAG=$(echo $S | tr " " x)
        "$BIN" $S 1 1 1 '"$NSTEPS"' '"$CDEGREE"' | tee /results/'"$V"'_${TAG}_bench.txt
        rocprofv3 --kernel-trace --output-format csv \
          -d /results/'"$V"'_${TAG} --output-file trace \
          -- "$BIN" $S 1 1 1 '"$NSTEPS"' '"$CDEGREE"'
      done
      # focused divergence SOL if rocprof-compute is available
      RPC=$(command -v rocprof-compute || command -v omniperf || echo "")
      if [ -n "$RPC" ]; then
        "$RPC" profile -n '"$V"'_div16 --path /results/omniperf \
          -k Divergence_3D_gpukernel -- "$BIN" 16 16 16 1 1 1 12 '"$CDEGREE"' || true
      fi
    '
}

run_variant base "$BASE_SQSH"
run_variant mod  "$MOD_SQSH"

echo "== Join base vs mod per-kernel totals offline: =="
echo "   for each */*kernel_trace*.csv aggregate by kernel name (see hotspot_sweep.slurm awk),"
echo "   then diff base_<size> vs mod_<size> tables -> delta_ns, delta_pct."
echo "Results under: $OUT"
