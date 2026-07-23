#!/bin/bash
# Hotspot sweep: per-kernel wall time for a single-GPU LinearEuler3D forward step
# across a few problem sizes on MI300X, using rocprofv3 kernel tracing.
#
# rocprofv3 --kernel-trace is single-pass and records one row per kernel dispatch
# (name + start/end timestamps). We aggregate per kernel name to get total time,
# launch count, and mean dispatch time -> the per-step hotspot breakdown.
#
# Benchmark at N=7 (cdegree default) until the N>=8 3-D path is confirmed on the
# built image; degrees are the last CLI arg if you want to override.
#
# Env overrides: SQSH (image), OUT (host results dir), NSTEPS, CDEGREE.
#SBATCH -J self-le3d-hotspot
#SBATCH -p main
#SBATCH -n1
#SBATCH --cpus-per-task=16
#SBATCH --mem=220G
#SBATCH --gres=gpu:1
#SBATCH --time=00:40:00
#SBATCH --exclude=node006
#SBATCH -o self-le3d-hotspot-%j.out
#SBATCH -e self-le3d-hotspot-%j.err
set -e

export ENROOT_TEMP_PATH="$HOME/.enroot/tmp" ENROOT_CACHE_PATH="$HOME/.enroot/cache" ENROOT_DATA_PATH="$HOME/.enroot/data" ENROOT_RUNTIME_PATH="$HOME/.enroot/run"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_CACHE_PATH" "$ENROOT_DATA_PATH" "$ENROOT_RUNTIME_PATH"
export MELLANOX_VISIBLE_DEVICES=none

# Use the *-prof image (has rocprofv3 from rocprofiler-sdk).
SQSH="${SQSH:-$HOME/self-gfx942-dp-prof.sqsh}"
OUT="${OUT:-$HOME/le3d-results/hotspot_$(date +%Y%m%d_%H%M%S)}"; mkdir -p "$OUT"
NSTEPS="${NSTEPS:-20}"
CDEGREE="${CDEGREE:-7}"

echo "SQSH=$SQSH  OUT=$OUT  NSTEPS=$NSTEPS  CDEGREE=$CDEGREE"

srun --gres=gpu:1 --mpi=none --exclude=node006 --container-image="$SQSH" \
     --container-mounts="${OUT}:/results" --container-workdir=/results \
  bash -lc '
    source /opt/spack-environment/activate.sh 2>/dev/null || true
    export LD_LIBRARY_PATH=/opt/self/build:/opt/self/install/lib:${LD_LIBRARY_PATH}
    export WORKSPACE=/opt/self/src
    export ROCR_VISIBLE_DEVICES=0
    for v in $(env | grep -E "^(PMIX_|PMI_|SLURM_|OMPI_|PRTE_)" | cut -d= -f1); do unset "$v"; done
    export OMPI_MCA_btl=self,vader OMPI_MCA_pml=ob1
    rocminfo | grep -m1 gfx || true
    echo "rocprofv3: $(command -v rocprofv3 || echo MISSING)"

    BIN=/opt/self/build/examples/bench_lineareuler3d_scaling
    for S in "4 4 4" "8 8 8" "16 16 16" "24 24 24"; do
      TAG=$(echo $S | tr " " x)
      echo "===== elements=$TAG  N='"$CDEGREE"' ====="
      # Headline timing/DOF-per-s from the benchmark itself
      "$BIN" $S 1 1 1 '"$NSTEPS"' '"$CDEGREE"' | tee /results/bench_${TAG}.txt
      # Per-dispatch kernel trace
      rocprofv3 --kernel-trace --output-format csv \
        -d /results/hotspot_${TAG} --output-file trace \
        -- "$BIN" $S 1 1 1 '"$NSTEPS"' '"$CDEGREE"'
    done

    echo "== aggregating per-kernel totals =="
    # gawk FPAT is quote-aware: rocprofv3 Kernel_Name contains commas inside
    # quotes (the C++ signature), which a plain -F, would split incorrectly.
    for CSV in $(find /results -name "*kernel_trace*.csv"); do
      OUTAGG="${CSV%.csv}.agg.csv"
      gawk "
        BEGIN{ FPAT=\"([^,]*)|(\\\"[^\\\"]*\\\")\" }
        NR==1{
          for(i=1;i<=NF;i++){ h=\$i; gsub(/\"/,\"\",h);
            if(h==\"Kernel_Name\") kn=i;
            if(h==\"Start_Timestamp\") st=i;
            if(h==\"End_Timestamp\") en=i; }
          next }
        { name=\$kn; gsub(/\"/,\"\",name); sub(/\\(.*/,\"\",name);
          d=\$en-\$st; tot[name]+=d; cnt[name]++ }
        END{ print \"kernel,count,total_ns,mean_ns\";
             for(k in tot) printf \"%s,%d,%.0f,%.1f\n\", k, cnt[k], tot[k], tot[k]/cnt[k] }
      " "$CSV" | sort -t, -k3 -nr > "$OUTAGG"
      echo "-> $OUTAGG"; head -12 "$OUTAGG"
    done
  '
echo "Results under: $OUT"
