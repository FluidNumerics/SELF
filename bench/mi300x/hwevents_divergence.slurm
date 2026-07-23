#!/bin/bash
# Hardware-event profile of the 3-D divergence GPU kernel chain on MI300X using
# rocprofv3 (rocprofiler-sdk). Collects per-dispatch hardware counters for the
# divergence chain kernels: occupancy (SQ_WAVES / SQ_BUSY_CYCLES), VALU work
# (SQ_INSTS_VALU and the FP64 breakdown), LDS pressure (SQ_INSTS_LDS), and L2
# cache behaviour (TCC_HIT / TCC_MISS -> hit rate). rocprofv3 also reports
# VGPR/SGPR/LDS block size and grid/workgroup size per dispatch.
#
# Requires the *-prof image (rocprofv3 from the rocprofiler-sdk package, added by
# install_profiler.slurm). IMPORTANT: run from a WRITABLE workdir (/results,
# bind-mounted) — rocprofv3 writes a .rocprofv3 scratch dir in CWD and aborts on
# the read-only /opt/self/build.
#
# rocprofv3 collects at most a few counters per pass; keep the --pmc list short
# (or split into multiple --pmc lines / runs). Env: SQSH, OUT, NELEM, CDEGREE.
#SBATCH -J self-le3d-hwevents
#SBATCH -p main
#SBATCH -n1
#SBATCH --cpus-per-task=16
#SBATCH --mem=220G
#SBATCH --gres=gpu:1
#SBATCH --time=00:40:00
#SBATCH --exclude=node006
#SBATCH -o self-le3d-hwevents-%j.out
#SBATCH -e self-le3d-hwevents-%j.err
set -e

export ENROOT_TEMP_PATH="$HOME/.enroot/tmp" ENROOT_CACHE_PATH="$HOME/.enroot/cache" ENROOT_DATA_PATH="$HOME/.enroot/data" ENROOT_RUNTIME_PATH="$HOME/.enroot/run"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_CACHE_PATH" "$ENROOT_DATA_PATH" "$ENROOT_RUNTIME_PATH"
export MELLANOX_VISIBLE_DEVICES=none

SQSH="${SQSH:-$HOME/self-gfx942-dp-prof.sqsh}"
OUT="${OUT:-$HOME/le3d-results/hwevents_$(date +%Y%m%d_%H%M%S)}"; mkdir -p "$OUT"
NELEM="${NELEM:-16}"
CDEGREE="${CDEGREE:-7}"
echo "SQSH=$SQSH OUT=$OUT NELEM=$NELEM CDEGREE=$CDEGREE"

srun --gres=gpu:1 --mpi=none --exclude=node006 --container-image="$SQSH" \
     --container-mounts="${OUT}:/results" --container-workdir=/results bash -lc '
    export PATH=/opt/rocm/bin:$PATH
    export LD_LIBRARY_PATH=/opt/self/build:/opt/self/install/lib:${LD_LIBRARY_PATH}
    export WORKSPACE=/opt/self/src ROCR_VISIBLE_DEVICES=0
    for v in $(env | grep -E "^(PMIX_|PMI_|SLURM_|OMPI_|PRTE_)" | cut -d= -f1); do unset "$v"; done
    export OMPI_MCA_btl=self,vader OMPI_MCA_pml=ob1
    BIN=/opt/self/build/examples/bench_lineareuler3d_scaling
    ARGS="'"$NELEM"' '"$NELEM"' '"$NELEM"' 1 1 1 6 '"$CDEGREE"'"
    KRE="Divergence_3D_gpukernel|ContravariantProjection_3D_gpukernel|DG_BoundaryContribution_3D_gpukernel|JacobianWeight|VectorDivergence_3D_gpukernel"

    # Occupancy + VALU (incl FP64) + LDS + L2 cache. Two passes keep each within
    # the hardware counter-slot budget.
    rocprofv3 --pmc SQ_WAVES SQ_BUSY_CYCLES GRBM_GUI_ACTIVE SQ_INSTS_VALU SQ_INSTS_LDS \
      --kernel-include-regex "$KRE" --output-format csv -d /results/pmc_occ --output-file occ \
      -- "$BIN" $ARGS >/results/occ.log 2>&1 || echo "occ rc=$?"
    rocprofv3 --pmc SQ_INSTS_VALU_FMA_F64 SQ_INSTS_VALU_ADD_F64 SQ_INSTS_VALU_MUL_F64 TCC_HIT TCC_MISS \
      --kernel-include-regex "$KRE" --output-format csv -d /results/pmc_mem --output-file mem \
      -- "$BIN" $ARGS >/results/mem.log 2>&1 || echo "mem rc=$?"

    echo "== counter CSVs =="; find /results -name "*counter_collection*.csv"
    # NOTE: kernel names contain commas inside quotes; aggregate with a
    # quote-aware parser (gawk FPAT or python), e.g.:
    #   gawk "BEGIN{FPAT=\"([^,]*)|(\\\"[^\\\"]*\\\")\"} ..." <csv>
    # Per kernel: L2 hit = TCC_HIT/(TCC_HIT+TCC_MISS); occupancy proxy =
    # SQ_WAVES vs peak; FP64 share from the *_F64 instruction counts.
  '
echo "Results under: $OUT"
