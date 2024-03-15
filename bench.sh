#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --gpus-per-task=1
#SBATCH --sockets=1 
#SBATCH --cpus-per-task=16 
#SBATCH --partition=mi210
#SBATCH -o stdout
#SBATCH -e stderr

export BATCHCOUNT=1000

# Use lowercase "op_n" for normal, "op_t" for transposed.
for OP in "op_n" "op_t"; do
    export OPERATION=$OP
    # Use 32 for single/32-bit, 64 for double/64-bit.
    for PREC in 32 64; do 
        export PRECISION=$PREC
        # Use the following to set the bounds and stepsize for number of rows/columns.
        # It is usually expected to make START equal to STEPSIZE.
        # If rows/columns will not have the same START/STEPSIZE/STOP, this will need
        # to be changed in their respective for loops.
        export STEPSIZE=8
        export START=$STEPSIZE
        export STOP=256
        for ROWS in $(seq $START $STEPSIZE $STOP); do
            for COLUMNS in $(seq $START $STEPSIZE $STOP); do
                export PROFILE_DIR=blas_results/${OPERATION}_${PRECISION}/${ROWS}_${COLUMNS}_${BATCHCOUNT}
                export FILENAME=build/blas/gemvstridedbatched_${OPERATION}_${PRECISION}
                # If file does not already exist or if the number of files in the folder is not the expected value (9).
                # The second condition is implemented to handle easy scancel/sbatch from the user.
                # I.e., if a user stops the job early, you end up with a folder that exists, but has incomplete/missing files.
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