#!/bin/bash
# Reconfirm correctness on the built gfx942 image before trusting profiles.
# Runs the full serial ctest suite (expected 124/124 in double precision) plus,
# separately, the new high-order N>=8 3-D divergence/gradient regressions that
# exercise the grid-strided kernel fix.
#
# Requires --container-writable (ctest writes Testing/) and WORKSPACE set so the
# mesh files under share/ resolve.
#
# Env overrides: SQSH (image), OUT (host results dir).
#SBATCH -J self-le3d-ctest
#SBATCH -p main
#SBATCH -n1
#SBATCH --cpus-per-task=16
#SBATCH --mem=220G
#SBATCH --gres=gpu:1
#SBATCH --time=00:40:00
#SBATCH --exclude=node006
#SBATCH -o self-le3d-ctest-%j.out
#SBATCH -e self-le3d-ctest-%j.err
set -e

export ENROOT_TEMP_PATH="$HOME/.enroot/tmp" ENROOT_CACHE_PATH="$HOME/.enroot/cache" ENROOT_DATA_PATH="$HOME/.enroot/data" ENROOT_RUNTIME_PATH="$HOME/.enroot/run"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_CACHE_PATH" "$ENROOT_DATA_PATH" "$ENROOT_RUNTIME_PATH"
export MELLANOX_VISIBLE_DEVICES=none

SQSH="${SQSH:-$HOME/self-gfx942-dp.sqsh}"
OUT="${OUT:-$HOME/le3d-results/ctest_$(date +%Y%m%d_%H%M%S)}"; mkdir -p "$OUT"

echo "SQSH=$SQSH  OUT=$OUT"

srun --gres=gpu:1 --mpi=none --exclude=node006 --container-writable --container-image="$SQSH" \
     --container-mounts="${OUT}:/results" --container-workdir=/opt/self/build \
  bash -lc '
    source /opt/spack-environment/activate.sh 2>/dev/null || true
    export LD_LIBRARY_PATH=/opt/self/build:/opt/self/install/lib:${LD_LIBRARY_PATH}
    export WORKSPACE=/opt/self/src
    export ROCR_VISIBLE_DEVICES=0
    for v in $(env | grep -E "^(PMIX_|PMI_|SLURM_|OMPI_|PRTE_)" | cut -d= -f1); do unset "$v"; done
    export OMPI_MCA_btl=self,vader OMPI_MCA_pml=ob1

    echo "== full serial suite =="
    ctest --output-on-failure -E "_mpi" 2>&1 | tee /results/ctest_full.log || true

    echo "== new high-order N>=8 regressions =="
    ctest --output-on-failure -R "3d_linear_n(8|12)" 2>&1 | tee /results/ctest_highN.log || true
  '
echo "Logs under: $OUT"
