#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --gpus-per-task=1
#SBATCH --sockets=1 
#SBATCH --cpus-per-task=16 
#SBATCH --partition=gpu
#SBATCH --constraint=mi210
#SBATCH -o stdout
#SBATCH -e stderr

# Conversions
gb_to_mb=1000
gb_to_kb=1000000
gb_to_b=1000000000
# mb_to_gb=0.001
mb_to_kb=1000
mb_to_b=1000000
# kb_to_gb=0.000001
# kb_to_mb=0.001
kb_to_b=1000

# Get available disk space and assign it to `available_disk_space`
available_disk_space=$(df -BM . | awk 'NR==2 {print $4}')
available_disk_space=$(echo "$available_disk_space" | sed 's/M//')

echo $available_disk_space

# `wiggle_room`: The amount of disk space you want to remain available
# during the execution of this script.
wiggle_room=$((10 * $gb_to_mb))

# Max directory size of 80MB has only been qualitatively determined based
# on 100 runs per kernel. I.e., number of rows in results.csv.
max_directory_size=80

echo $max_directory_size

if [[ $available_disk_space -lt $wiggle_room ]]; then
    echo "uh oh"
fi

########
# axpy #
########

export SUBROUTINE="axpy"

# Use 32 for single/32-bit, 64 for double/64-bit.
for PREC in 32 64; do 
    export PRECISION=$PREC
    export STEPSIZE=8
    export START=$STEPSIZE
    export STOP=16384
    for ROWS in $(seq $START $STEPSIZE $STOP); do
        export PROFILE_DIR=blas_results/${SUBROUTINE}_${PRECISION}/${ROWS}
        export FILENAME=build/blas/${SUBROUTINE}_${PRECISION}
        export EXPECTEDFILES=4
        if [ ! -d "$PROFILE_DIR" ] || [ "$(ls -l "$PROFILE_DIR" | grep "^-" | wc -l)" -ne $EXPECTEDFILES ]; then
            mkdir -p $PROFILE_DIR

            source ~/.bashrc
            source /etc/profile.d/z11_lmod.sh
            module --default avail
            module load gcc/13.2.0
            module load hip/5.7.3

            # Flat profile
            rocprof --timestamp on $FILENAME $ROWS
            mv results.csv $PROFILE_DIR/

            # Hotspot profile
            rocprof --stats $FILENAME $ROWS
            mv results.stats.csv $PROFILE_DIR/

            # Trace Profile
            rocprof --sys-trace $FILENAME $ROWS
            mv results.json $PROFILE_DIR/

            # Hardware events profile (for bandwidth estimates and L2 Cache hit)
            rocprof -i events_axpy.txt $FILENAME $ROWS
            mv events_axpy.csv $PROFILE_DIR/
            mv $PROFILE_DIR/events_axpy.csv $PROFILE_DIR/events.csv
        fi
    done
done

######################
# gemvstridedbatched #
######################

export SUBROUTINE="gemvstridedbatched"

export BATCHCOUNT=1000

for OP in "op_n" "op_t"; do
    export OPERATION=$OP
    for PREC in 32 64; do 
        export PRECISION=$PREC
        export STEPSIZE=8
        export START=$STEPSIZE
        export STOP=256
        for ROWS in $(seq $START $STEPSIZE $STOP); do
            for COLUMNS in $(seq $START $STEPSIZE $STOP); do
                export PROFILE_DIR=blas_results/${SUBROUTINE}_${OPERATION}_${PRECISION}/${ROWS}_${COLUMNS}_${BATCHCOUNT}
                export FILENAME=build/blas/${SUBROUTINE}_${OPERATION}_${PRECISION}
                export EXPECTEDFILES=9
                if [ ! -d "$PROFILE_DIR" ] || [ "$(ls -l "$PROFILE_DIR" | grep "^-" | wc -l)" -ne $EXPECTEDFILES ]; then
                    mkdir -p $PROFILE_DIR

                    source ~/.bashrc
                    source /etc/profile.d/z11_lmod.sh
                    module --default avail
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
                    rocprof -i events_gemvsb.txt $FILENAME $ROWS $COLUMNS $BATCHCOUNT
                    mv events_gemvsb.csv $PROFILE_DIR/
                    mv $PROFILE_DIR/events_gemvsb.csv $PROFILE_DIR/events.csv

                    mv results.* $PROFILE_DIR/
                    mv stdout $PROFILE_DIR/
                fi
            done
        done
    done
done





########
# gemm #
########

export SUBROUTINE="gemm"

for OP in "op_n" "op_t"; do
    export OPERATION=$OP
    for PREC in 32 64; do 
        export PRECISION=$PREC
        export STEPSIZE=4
        export START=$STEPSIZE
        export STOP=16
        for ROWS in $(seq $START $STEPSIZE $STOP); do
            for COLUMNS_FACTOR in 1000 10000 100000; do
                let COLUMNS=ROWS*COLUMNS_FACTOR
                export COLUMNS
                export PROFILE_DIR=blas_results/${SUBROUTINE}_${OPERATION}_${PRECISION}/${ROWS}_${COLUMNS}
                export FILENAME=build/blas/${SUBROUTINE}_${OPERATION}_${PRECISION}
                export EXPECTEDFILES=4
                if [ ! -d "$PROFILE_DIR" ] || [ "$(ls -l "$PROFILE_DIR" | grep "^-" | wc -l)" -ne $EXPECTEDFILES ]; then
                    mkdir -p $PROFILE_DIR

                    source ~/.bashrc
                    source /etc/profile.d/z11_lmod.sh
                    module --default avail
                    module load gcc/13.2.0
                    module load hip/5.7.3

                    # Flat profile
                    rocprof --timestamp on $FILENAME $ROWS $COLUMNS
                    mv results.csv $PROFILE_DIR/

                    # Hotspot profile
                    rocprof --stats $FILENAME $ROWS $COLUMNS
                    mv results.stats.csv $PROFILE_DIR/

                    # Trace Profile
                    rocprof --sys-trace $FILENAME $ROWS $COLUMNS
                    mv results.json $PROFILE_DIR/

                    # Hardware events profile (for bandwidth estimates and L2 Cache hit)
                    rocprof -i events_gemm.txt $FILENAME $ROWS $COLUMNS
                    mv events_gemm.csv $PROFILE_DIR/
                    mv $PROFILE_DIR/events_gemm.csv $PROFILE_DIR/events.csv
                fi
            done
        done
    done
done

