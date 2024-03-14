#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --gpus-per-task=1
#SBATCH --sockets=1 
#SBATCH --cpus-per-task=16 
#SBATCH --partition=mi210
#SBATCH -o stdout
#SBATCH -e stderr

for ROWS in $(seq 8 8 128) 
do
    for COLUMNS in $(seq 8 8 128)
    do
        export OPERATION=op_t # use lowercase (op_n, op_t)
        export PRECISION=64
        export BATCHCOUNT=1000
        export PROFILE_DIR=../blas_results/${OPERATION}_${PRECISION}/${ROWS}_${COLUMNS}_${BATCHCOUNT}
        export FILENAME=../build/blas/gemvstridedbatched_${OPERATION}_${PRECISION}
        # export FILENAME=build/blas/gemvstridedbatched_op_n_64
        if [ ! -d "$PROFILE_DIR" ]; then
            mkdir -p $PROFILE_DIR

            source ~/.bashrc
            source /etc/profile.d/z11_lmod.sh
            module avail
            module load gcc/13.2.0
            module load hip/5.7.3

            # # Flat profile
            # rocprof --timestamp on $FILENAME $ROWS $COLUMNS $BATCHCOUNT
            # mv results.csv $PROFILE_DIR/

            # # Hotspot profile
            # rocprof --stats $FILENAME $ROWS $COLUMNS $BATCHCOUNT
            # mv results.stats.csv $PROFILE_DIR/

            # # # Trace Profile
            # rocprof --sys-trace $FILENAME $ROWS $COLUMNS $BATCHCOUNT
            # mv results.json $PROFILE_DIR/

            # Hardware events profile (for bandwidth estimates and L2 Cache hit)
            rocprof -i events.txt $FILENAME $ROWS $COLUMNS $BATCHCOUNT
            mv events.csv $PROFILE_DIR/

            mv results.* $PROFILE_DIR/
            mv stdout $PROFILE_DIR/
        fi
    done
done