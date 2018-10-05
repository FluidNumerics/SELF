#!/bin/bash
# 
# How many nodes are needed
#SBATCH --nodes=4
#
#SBATCH --ntasks=4
#
# How many MPI tasks per node
#SBATCH --ntasks-per-node=1  
#
# How many physical CPU’s per task
#SBATCH --cpus-per-task=1
#  
# How long the job is anticipated to run
#SBATCH --time=10:00:00
# 
# The name of the job
#SBATCH --job-name=self_testing 
#

pwd; hostname

module purge
module use ${HOME}/apps/modulefiles
export EX_DIR="${HOME}/SELF-Fluids/examples/boundarylayer"
export RUN_DIR="${HOME}/SELF-Fluids/testing"
export OMP_NUM_THREADS=8


run_model(){

  sfluid --param-file ${EX_DIR}/runtime.params \
         --equation-file ${EX_DIR}/self.equations

  if [ ! -d "$1" ]; then

    mkdir -p $1

  fi

  mv Pipeline.* State.* box.* ExtComm.* Diagnostics.* mesh.* $1
  mv Timing.stats $1

}

run_mpi_model(){

  mpirun -np ${SLURM_NTASKS} sfluid --param-file ${EX_DIR}/runtime.params \
                                    --equation-file ${EX_DIR}/self.equations

  if [ ! -d "$1" ]; then

    mkdir -p $1

  fi

  mv  State.* box.* ExtComm.* Diagnostics.* mesh.* $1
  mv Timing.stats $1

}

## PGI (18.4) Single Process
#module load pgi/18.4 hdf5/1.10.2/serial/pgi/18.4

## Serial
#echo "Running PGI 18.4 Serial"
#module load self-fluids/dev/pgi/18.4/serial
#run_model ${RUN_DIR}/pgi/18.4/serial
#module unload self-fluids
#
#echo "Running PGI 18.4 Serial (Double)"
#module load self-fluids/dev/pgi/18.4/serial_double
#run_model ${RUN_DIR}/pgi/18.4/serial_double
#module unload self-fluids
#
## CUDA
#echo "Running PGI 18.4 CUDA"
#module load self-fluids/dev/pgi/18.4/cuda
#run_model ${RUN_DIR}/pgi/18.4/cuda
#module unload self-fluids
#
#echo "Running PGI 18.4 CUDA (Double)"
#module load self-fluids/dev/pgi/18.4/cuda_double
#run_model ${RUN_DIR}/pgi/18.4/cuda_double
#module unload self-fluids


## OpenMP
#echo "Running PGI 18.4 OpenMP"
#module load self-fluids/dev/pgi/18.4/openmp
#run_model ${RUN_DIR}/pgi/18.4/openmp
#module unload self-fluids

#echo "Running PGI 18.4 OpenMP (Double)"
#module load self-fluids/dev/pgi/18.4/openmp_double
#run_model ${RUN_DIR}/pgi/18.4/openmp_double
#module unload self-fluids


## PGI (18.4) MPI
#module unload hdf5

#module load pgi/18.4
#module load openmpi/3.1.2/pgi/18.4 hdf5/1.10.2/parallel/pgi/18.4/openmpi/3.1.2

#echo "Running PGI 18.4 MPI"
#module load self-fluids/dev/pgi/18.4/openmpi_3.1.2
#run_mpi_model ${RUN_DIR}/pgi/18.4/openmp_3.1.2
#module unload self-fluids

#echo "Running PGI 18.4 MPI (Double)"
#module load self-fluids/dev/pgi/18.4/openmpi_3.1.2_double
#run_mpi_model ${RUN_DIR}/pgi/18.4/openmp_3.1.2_double
#module unload self-fluids
#
#echo "Running PGI 18.4 MPI+CUDA"
#module load self-fluids/dev/pgi/18.4/cuda+openmpi_3.1.2
#run_mpi_model ${RUN_DIR}/pgi/18.4/cuda+openmp_3.1.2
#module unload self-fluids

#echo "Running PGI 18.4 MPI+CUDA (Double)"
#module load self-fluids/dev/pgi/18.4/cuda+openmpi_3.1.2_double
#run_mpi_model ${RUN_DIR}/pgi/18.4/cuda+openmp_3.1.2_double
#module unload self-fluids

#echo "Running PGI 18.4 MPI+OpenMP"
#module load self-fluids/dev/pgi/18.4/openmp+openmpi_3.1.2
#run_mpi_model ${RUN_DIR}/pgi/18.4/openmp+openmp_3.1.2
#module unload self-fluids

#echo "Running PGI 18.4 MPI+OpenMP (Double)"
#module load self-fluids/dev/pgi/18.4/openmp+openmpi_3.1.2_double
#run_mpi_model ${RUN_DIR}/pgi/18.4/openmp+openmp_3.1.2_double
#module unload self-fluids

#module unload openmpi hdf5 pgi

# GCC (8.1.0) Single Process
module load gcc/8.1.0 hdf5/1.10.2/serial/gcc/8.1.0

# Serial
echo "Running GCC 8.1.0 Serial"
module load self-fluids/dev/gcc/8.1.0/serial
run_model ${RUN_DIR}/gcc/8.1.0/serial
module unload self-fluids

echo "Running GCC 8.1.0 Serial (Double)"
module load self-fluids/dev/gcc/8.1.0/serial_double
run_model ${RUN_DIR}/gcc/8.1.0/serial_double
module unload self-fluids

# OpenMP
echo "Running GCC 8.1.0 OpenMP"
module load self-fluids/dev/gcc/8.1.0/openmp
run_model ${RUN_DIR}/gcc/8.1.0/openmp
module unload self-fluids

echo "Running GCC 8.1.0 OpenMP (Double)"
module load self-fluids/dev/gcc/8.1.0/openmp_double
run_model ${RUN_DIR}/gcc/8.1.0/openmp_double
module unload self-fluids


# GCC (8.1.0) MPI
module unload hdf5

module load openmpi/3.1.2/gcc/8.1.0 hdf5/1.10.2/parallel/gcc/8.1.0/openmpi/3.1.2

echo "Running GCC 8.1.0 MPI"
module load self-fluids/dev/gcc/8.1.0/openmpi_3.1.2
run_mpi_model ${RUN_DIR}/gcc/8.1.0/openmp_3.1.2
module unload self-fluids

echo "Running GCC 8.1.0 MPI (Double)"
module load self-fluids/dev/gcc/8.1.0/openmpi_3.1.2_double
run_mpi_model ${RUN_DIR}/gcc/8.1.0/openmp_3.1.2_double
module unload self-fluids

echo "Running GCC 8.1.0 MPI+OpenMP"
module load self-fluids/dev/gcc/8.1.0/openmp+openmpi_3.1.2
run_mpi_model ${RUN_DIR}/gcc/8.1.0/openmp+openmp_3.1.2
module unload self-fluids

echo "Running GCC 8.1.0 MPI+OpenMP (Double)"
module load self-fluids/dev/gcc/8.1.0/openmp+openmpi_3.1.2_double
run_mpi_model ${RUN_DIR}/gcc/8.1.0/openmp+openmp_3.1.2_double
module unload self-fluids


module unload openmpi hdf5
module unload gcc
