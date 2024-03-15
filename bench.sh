#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --gpus-per-task=1
#SBATCH --sockets=1 
#SBATCH --cpus-per-task=16 
#SBATCH --partition=mi210
#SBATCH -o stdout
#SBATCH -e stderr
#SBATCH --job-name=op_n_32

export BATCHCOUNT=1000

for OP in "op_n" "op_t"
do
    export OPERATION=$OP # use lowercase (op_n, op_t)
    for PREC in 32 64
    do
        export PRECISION=$PREC
        for ROWS in $(seq 8 8 256) 
        do
            for COLUMNS in $(seq 8 8 256)
            do
                export PROFILE_DIR=blas_results/${OPERATION}_${PRECISION}/${ROWS}_${COLUMNS}_${BATCHCOUNT}
                export FILENAME=build/blas/gemvstridedbatched_${OPERATION}_${PRECISION}
                #if [ ! -d "$PROFILE_DIR" ]; then
                if [ ! -d "$PROFILE_DIR" ] || [ "$(ls -l "$PROFILE_DIR" | grep "^-" | wc -l)" -ne 9 ]; then
                    mkdir -p $PROFILE_DIR

                    source ~/.bashrc
                    source /etc/profile.d/z11_lmod.sh
                    module avail
                    module load gcc/13.2.0
                    module load hip/5.7.3

                    # Flat profile
                    rocprof --timestamp on $FILENAME $ROWS $COLUMNS $BATCHCOUNT
                    mv results.csv $PROFILE_DIR/

                    # Hotspot profile
                    rocprof --stats $FILENAME $ROWS $COLUMNS $BATCHCOUNT
                    mv results.stats.csv $PROFILE_DIR/

                    # Trace Profile
                    rocprof --sys-trace $FILENAME $ROWS $COLUMNS $BATCHCOUNT
                    mv results.json $PROFILE_DIR/

                    # Hardware events profile (for bandwidth estimates and L2 Cache hit)
                    rocprof -i events.txt $FILENAME $ROWS $COLUMNS $BATCHCOUNT
                    mv events.csv $PROFILE_DIR/

                    mv results.* $PROFILE_DIR/
                    mv stdout $PROFILE_DIR/
                fi
            done
        done
    done
done