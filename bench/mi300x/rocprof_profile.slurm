#!/bin/bash
# Hotspot + hardware-event profile of the LinearEuler3D forward step on MI300X
# using rocprof (v1, from the rocprofiler-2.0 package installed into the *-prof
# image). rocprofv3 is not available on this image; rocprof v1 covers both:
#   * --stats            -> per-kernel call count / total / average duration
#   * -i <counters.txt>  -> hardware counters (occupancy, VALU/mem busy, cache,
#                           LDS bank conflicts) per kernel dispatch
#
# Profile at 16^3 / N=7 (fills the CUs). Env: SQSH (prof image), OUT, NELEM,
# NSTEPS, CDEGREE.
#SBATCH -J self-le3d-rocprof
#SBATCH -p main
#SBATCH -n1
#SBATCH --cpus-per-task=16
#SBATCH --mem=220G
#SBATCH --gres=gpu:1
#SBATCH --time=00:40:00
#SBATCH --exclude=node006
#SBATCH -o self-le3d-rocprof-%j.out
#SBATCH -e self-le3d-rocprof-%j.err
set -e

export ENROOT_TEMP_PATH="$HOME/.enroot/tmp" ENROOT_CACHE_PATH="$HOME/.enroot/cache" ENROOT_DATA_PATH="$HOME/.enroot/data" ENROOT_RUNTIME_PATH="$HOME/.enroot/run"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_CACHE_PATH" "$ENROOT_DATA_PATH" "$ENROOT_RUNTIME_PATH"
export MELLANOX_VISIBLE_DEVICES=none

SQSH="${SQSH:-$HOME/self-gfx942-dp-prof.sqsh}"
OUT="${OUT:-$HOME/le3d-results/rocprof_$(date +%Y%m%d_%H%M%S)}"; mkdir -p "$OUT"
NELEM="${NELEM:-16}"
NSTEPS="${NSTEPS:-12}"
CDEGREE="${CDEGREE:-7}"
echo "SQSH=$SQSH OUT=$OUT NELEM=$NELEM NSTEPS=$NSTEPS CDEGREE=$CDEGREE"

srun --gres=gpu:1 --mpi=none --exclude=node006 --container-image="$SQSH" \
     --container-mounts="${OUT}:/results" --container-workdir=/opt/self/build \
  bash -lc '
    export PATH=/opt/rocm/bin:$PATH
    export LD_LIBRARY_PATH=/opt/self/build:/opt/self/install/lib:${LD_LIBRARY_PATH}
    export WORKSPACE=/opt/self/src ROCR_VISIBLE_DEVICES=0
    for v in $(env | grep -E "^(PMIX_|PMI_|SLURM_|OMPI_|PRTE_)" | cut -d= -f1); do unset "$v"; done
    export OMPI_MCA_btl=self,vader OMPI_MCA_pml=ob1
    BIN=/opt/self/build/examples/bench_lineareuler3d_scaling
    ARGS="'"$NELEM"' '"$NELEM"' '"$NELEM"' 1 1 1 '"$NSTEPS"' '"$CDEGREE"'"

    echo "===== hotspots: rocprof --stats ====="
    rocprof --stats --timestamp on -o /results/hotspots.csv "$BIN" $ARGS 2>/results/hotspots.log || true
    echo "-- per-kernel stats (sorted by total duration) --"
    if [ -f /results/hotspots.stats.csv ]; then
      sort -t, -k4 -nr /results/hotspots.stats.csv | head -20
    fi

    echo "===== hardware counters: rocprof -i counters.txt ====="
    cat > /results/counters.txt <<CNT
pmc: Wavefronts VALUInsts SALUInsts VALUUtilization VALUBusy SALUBusy
pmc: MemUnitBusy MemUnitStalled WriteUnitStalled L2CacheHit
pmc: LDSBankConflict GPUBusy FetchSize WriteSize
CNT
    rocprof -i /results/counters.txt --timestamp on -o /results/hwcounters.csv "$BIN" $ARGS 2>/results/hwcounters.log || true
    echo "-- divergence-chain counter rows --"
    if [ -f /results/hwcounters.csv ]; then
      head -1 /results/hwcounters.csv
      grep -E "Divergence_3D|ContravariantProjection_3D|DG_BoundaryContribution_3D|JacobianWeight|VectorDivergence_3D" /results/hwcounters.csv | head -40
    fi
    echo "== outputs =="; ls -l /results/*.csv 2>/dev/null
  '
echo "Results under: $OUT"
