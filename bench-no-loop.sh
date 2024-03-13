#!/bin/bash
#SBATCH -n1 
#SBATCH --gpus-per-task=1
#SBATCH --sockets=1 
#SBATCH --cpus-per-task=16 
#SBATCH --partition=mi210
#SBATCH -o stdout
#SBATCH -e stdout

export OPERATION=OP_N
export PRECISION=64
export ROWS=11
export COLUMNS=10
export BATCHCOUNT=100
export PROFILE_DIR=blas_results/${OPERATION}_${PRECISION}/${ROWS}_${COLUMNS}_${BATCHCOUNT}
export FILENAME=build/blas/gemvstridedbatched_op_n_64
mkdir -p $PROFILE_DIR

source ~/.bashrc
source /etc/profile.d/z11_lmod.sh
module avail
module load gcc/13.2.0
module load hip/5.7.3

# Flat profile
rocprof --timestamp on $FILENAME $ROWS $COLUMNS
mv results.csv $PROFILE_DIR/

# Hotspot profile
rocprof --stats $FILENAME $ROWS $COLUMNS
mv results.stats.csv $PROFILE_DIR/

# # Trace Profile
rocprof --sys-trace $FILENAME $ROWS $COLUMNS
mv results.json $PROFILE_DIR/

# Hardware events profile (for bandwidth estimates and L2 Cache hit)
rocprof -i events.txt $FILENAME $ROWS $COLUMNS
mv events.csv $PROFILE_DIR/

mv results.* $PROFILE_DIR/
mv stdout $PROFILE_DIR/