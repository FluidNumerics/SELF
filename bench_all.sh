#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --gpus-per-task=1
#SBATCH --sockets=1 
#SBATCH --cpus-per-task=16 
#SBATCH --partition=gpu
#SBATCH --constraint=mi210
#SBATCH -o stdout
#SBATCH -e stderr

# # Conversions
# gb_to_mb=1000
# gb_to_kb=1000000
# gb_to_b=1000000000
# # mb_to_gb=0.001
# mb_to_kb=1000
# mb_to_b=1000000
# # kb_to_gb=0.000001
# # kb_to_mb=0.001
# kb_to_b=1000

# unused block; intention is to check available disk space
###############################################################################
# # Get available disk space and assign it to `available_disk_space`
# available_disk_space=$(df -BM . | awk 'NR==2 {print $4}')
# available_disk_space=$(echo "$available_disk_space" | sed 's/M//')

# echo $available_disk_space

# # `wiggle_room`: The amount of disk space you want to remain available
# # during the execution of this script.
# wiggle_room=$((10 * $gb_to_mb))

# # Max directory size of 80MB has only been qualitatively determined based
# # on 100 runs per kernel. I.e., number of rows in results.csv.
# max_directory_size=80

# echo $max_directory_size

# if ([ $available_disk_space -lt $wiggle_room ]); then
#     echo "uh oh"
# fi
###############################################################################

# flags
# default behavior: perform all rocprof, omniperf and valgrind must be set
# manually. currently no way to check if omniperf is done; i.e., omniperf will 
# be performed every time if the flag is given
rocprof_all_flag='true'
rocprof_flat_flag=''
rocprof_hotspot_flag=''
rocprof_trace_flag=''
rocprof_events_flag=''
omniperf_flag=''
valgrind_flag=''
all_flag=''
EXPECTEDFILES=0
FILENUMBER=0

while getopts 'rfhteov' flag; do
  case "${flag}" in
    r) rocprof_all_flag='true' 
        let "EXPECTEDFILES=EXPECTEDFILES+1111";;
    f) rocprof_flat_flag='true' 
        let "EXPECTEDFILES=EXPECTEDFILES+1";;
    h) rocprof_hotspot_flag='true' 
        let "EXPECTEDFILES=EXPECTEDFILES+10";;
    t) rocprof_trace_flag='true' 
        let "EXPECTEDFILES=EXPECTEDFILES+100";;
    e) rocprof_events_flag='true' 
        let "EXPECTEDFILES=EXPECTEDFILES+1000";;
    o) omniperf_flag='true' ;;
        # let "EXPECTEDFILES=EXPECTEDFILES+10000";;
    v) valgrind_flag='true' 
        let "EXPECTEDFILES=EXPECTEDFILES+10000";;
    *) error "Unexpected option ${flag}" ;;
  esac
done

if (([ "$rocprof_flat_flag" = true ] || [ "$rocprof_hotspot_flag" = true ] || [ "$rocprof_trace_flag" = true ] || [ "$rocprof_events_flag" = true ] || [ "$omniperf_flag" = true ] || [ "$valgrind_flag" = true ])) ; then
    rocprof_all_flag='false'
fi

# if rocprof_all is used in combination with other rocprof flags
EXPECTEDFILES_STRING="$EXPECTEDFILES"
EXPECTEDFILES_STRING="${EXPECTEDFILES_STRING//2/1}"
EXPECTEDFILES=$((EXPECTEDFILES_STRING))

get_file_number () {
    FILENUMBER=0
    if [ -f $1/results.csv ] ; then
        FILENUMBER=FILENUMBER+1
    fi
    if [ -f $1/resultsstats.csv ] ; then
        FILENUMBER=FILENUMBER+10
    fi
    if [ -f $1/results.json ] ; then
        FILENUMBER=FILENUMBER+100
    fi
    if [ -f $1/results.csv ] ; then
        FILENUMBER=FILENUMBER+1000
    fi
    if [ -f $1/callgrind.out.* ] ; then
        FILENUMBER=FILENUMBER+10000
    fi
}

########
# axpy #
########

export SUBROUTINE="axpy"

# Use 32 for single/32-bit, 64 for double/64-bit.
for PREC in 32 64; do 
    export PRECISION=$PREC
    export STEPSIZE=8
    export START=$STEPSIZE
    export STOP=4096
    for ROWS in $(seq $START $STEPSIZE $STOP); do
        export PROFILE_DIR=blas_results/${SUBROUTINE}_${PRECISION}/${ROWS}
        export FILENAME=build/blas/${SUBROUTINE}_${PRECISION}

        get_file_number $PROFILE_DIR 
        if [ $FILENUMBER = $EXPECTEDFILES ]; then
            mkdir -p $PROFILE_DIR

            if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_flat_flag" = true ] || [ "$rocprof_hotspot_flag" = true ] || [ "$rocprof_trace_flag" = true ] || [ "$rocprof_events_flag" = true ]) ; then
                source ~/.bashrc
                source /etc/profile.d/z11_lmod.sh
                module --default avail
                module load gcc/13.2.0
                module load hip/6.1.2
            fi
            if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_flat_flag" = true ]) ; then
                # Flat profile
                rocprof --timestamp on $FILENAME $ROWS
                mv results.csv $PROFILE_DIR/
            fi
            if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_hotspot_flag" = true ]) ; then
                # Hotspot profile
                rocprof --stats $FILENAME $ROWS
                mv results.stats.csv $PROFILE_DIR/
            fi
            if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_trace_flag" = true ]) ; then
                # Trace Profile
                rocprof --sys-trace $FILENAME $ROWS
                mv results.json $PROFILE_DIR/
            fi
            if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_events_flag" = true ]) ; then
                # Hardware events profile (for bandwidth estimates and L2 Cache hit)
                rocprof -i events_axpy.txt $FILENAME $ROWS
                mv events_axpy.csv $PROFILE_DIR/events.csv
                # mv $PROFILE_DIR/events_axpy.csv $PROFILE_DIR/events.csv
            fi
            if ([ "$all_flag" = true ] || [ "$omniperf_flag" = true]) ; then
                module load omniperf
                omniperf profile --name ${SUBROUTINE}_${PRECISION} -- .$FILENAME $ROWS
            fi
            if ([ "$all_flag" = true ] || [ "$valgrind_flag" = true]) ; then
                valgrind --tool=callgrind $FILENAME $ROWS
                mv callgrind.out.* $PROFILE_DIR/callgrind.out.*
            fi
        fi
    done
done

########
# gemv #
########

export SUBROUTINE="gemv"

for OP in "op_n" "op_t"; do
    export OPERATION=$OP
    for PREC in 32 64; do 
        export PRECISION=$PREC
        export STEPSIZE=8
        export START=$STEPSIZE
        export STOP=256
        for ROWS in $(seq $START $STEPSIZE $STOP); do
            for COLUMNS in $(seq $START $STEPSIZE $STOP); do
                export PROFILE_DIR=blas_results/${SUBROUTINE}_${OPERATION}_${PRECISION}/${ROWS}_${COLUMNS}
                export FILENAME=build/blas/${SUBROUTINE}_${OPERATION}_${PRECISION}

                get_file_number $PROFILE_DIR
                if [ $FILENUMBER = $EXPECTEDFILES ]; then
                    mkdir -p $PROFILE_DIR

                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_flat_flag" = true ] || [ "$rocprof_hotspot_flag" = true ] || [ "$rocprof_trace_flag" = true ] || [ "$rocprof_events_flag" = true ]) ; then
                        source ~/.bashrc
                        source /etc/profile.d/z11_lmod.sh
                        module --default avail
                        module load gcc/13.2.0
                        module load hip/6.1.2
                    fi
                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_flat_flag" = true ]) ; then
                        # Flat profile
                        rocprof --timestamp on $FILENAME $ROWS $COLUMNS
                        mv results.csv $PROFILE_DIR/
                    fi
                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_hotspot_flag" = true ]) ; then
                        # Hotspot profile
                        rocprof --stats $FILENAME $ROWS $COLUMNS
                        mv results.stats.csv $PROFILE_DIR/
                    fi
                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_trace_flag" = true ]) ; then
                        # Trace Profile
                        rocprof --sys-trace $FILENAME $ROWS $COLUMNS
                        mv results.json $PROFILE_DIR/
                    fi
                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_events_flag" = true ]) ; then
                        # Hardware events profile (for bandwidth estimates and L2 Cache hit)
                        rocprof -i events_gemvsb.txt $FILENAME $ROWS $COLUMNS
                        mv events_gemvsb.csv $PROFILE_DIR/events.csv
                        # mv $PROFILE_DIR/events_gemvsb.csv $PROFILE_DIR/events.csv
                    fi
                    if ([ "$all_flag" = true ] || [ "$omniperf_flag" = true]) ; then
                        module load omniperf
                        omniperf profile --name ${SUBROUTINE}_${OPERATION}_${PRECISION} -- .${FILENAME} $ROWS $COLUMNS
                    fi
                    if ([ "$all_flag" = true ] || [ "$valgrind_flag" = true]) ; then
                        valgrind --tool=callgrind $FILENAME $ROWS $COLUMNS
                        mv callgrind.out.* $PROFILE_DIR/callgrind.out.*
                    fi
                fi
            done
        done
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

                get_file_number $PROFILE_DIR
                if [ $FILENUMBER = $EXPECTEDFILES ]; then
                    mkdir -p $PROFILE_DIR

                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_flat_flag" = true ] || [ "$rocprof_hotspot_flag" = true ] || [ "$rocprof_trace_flag" = true ] || [ "$rocprof_events_flag" = true ]) ; then
                        source ~/.bashrc
                        source /etc/profile.d/z11_lmod.sh
                        module --default avail
                        module load gcc/13.2.0
                        module load hip/6.1.2
                    fi
                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_flat_flag" = true ]) ; then
                        # Flat profile
                        rocprof --timestamp on $FILENAME $ROWS $COLUMNS $BATCHCOUNT
                        mv results.csv $PROFILE_DIR/
                    fi
                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_hotspot_flag" = true ]) ; then
                        # Hotspot profile
                        rocprof --stats $FILENAME $ROWS $COLUMNS $BATCHCOUNT
                        mv results.stats.csv $PROFILE_DIR/
                    fi
                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_trace_flag" = true ]) ; then
                        # Trace Profile
                        rocprof --sys-trace $FILENAME $ROWS $COLUMNS $BATCHCOUNT
                        mv results.json $PROFILE_DIR/
                    fi
                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_events_flag" = true ]) ; then
                        # Hardware events profile (for bandwidth estimates and L2 Cache hit)
                        rocprof -i events_gemvsb.txt $FILENAME $ROWS $COLUMNS $BATCHCOUNT
                        mv events_gemvsb.csv $PROFILE_DIR/events.csv
                        # mv $PROFILE_DIR/events_gemvsb.csv $PROFILE_DIR/events.csv
                    fi
                    if ([ "$all_flag" = true ] || [ "$omniperf_flag" = true]) ; then
                        module load omniperf
                        omniperf profile --name ${SUBROUTINE}_${OPERATION}_${PRECISION} -- .${FILENAME} $ROWS $COLUMNS $BATCHCOUNT
                    fi
                    if ([ "$all_flag" = true ] || [ "$valgrind_flag" = true]) ; then
                        valgrind --tool=callgrind $FILENAME $ROWS $COLUMNS $BATCHCOUNT
                        mv callgrind.out.* $PROFILE_DIR/callgrind.out.*
                    fi
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

                get_file_number $PROFILE_DIR
                if [ $FILENUMBER = $EXPECTEDFILES ]; then
                    mkdir -p $PROFILE_DIR

                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_flat_flag" = true ] || [ "$rocprof_hotspot_flag" = true ] || [ "$rocprof_trace_flag" = true ] || [ "$rocprof_events_flag" = true ]) ; then
                        source ~/.bashrc
                        source /etc/profile.d/z11_lmod.sh
                        module --default avail
                        module load gcc/13.2.0
                        module load hip/6.1.2
                    fi
                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_flat_flag" = true ]) ; then
                        # Flat profile
                        rocprof --timestamp on $FILENAME $ROWS $COLUMNS
                        mv results.csv $PROFILE_DIR/
                    fi
                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_hotspot_flag" = true ]) ; then
                        # Hotspot profile
                        rocprof --stats $FILENAME $ROWS $COLUMNS
                        mv results.stats.csv $PROFILE_DIR/
                    fi
                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_trace_flag" = true ]) ; then
                        # Trace Profile
                        rocprof --sys-trace $FILENAME $ROWS $COLUMNS
                        mv results.json $PROFILE_DIR/
                    fi
                    if ([ "$all_flag" = true ] || [ "$rocprof_all_flag" = true ] || [ "$rocprof_events_flag" = true ]) ; then
                        # Hardware events profile (for bandwidth estimates and L2 Cache hit)
                        rocprof -i events_gemm.txt $FILENAME $ROWS $COLUMNS
                        mv events_gemm.csv $PROFILE_DIR/events.csv
                        # mv $PROFILE_DIR/events_gemm.csv $PROFILE_DIR/events.csv
                    fi
                    if ([ "$all_flag" = true ] || [ "$omniperf_flag" = true]) ; then
                        module load omniperf
                        omniperf profile --name ${SUBROUTINE}_${OPERATION}_${PRECISION} -- .${FILENAME} $ROWS $COLUMNS
                    fi
                    if ([ "$all_flag" = true ] || [ "$valgrind_flag" = true]) ; then
                        valgrind --tool=callgrind $FILENAME $ROWS $COLUMNS
                        mv callgrind.out.* $PROFILE_DIR/callgrind.out.*
                    fi
                fi
            done
        done
    done
done

