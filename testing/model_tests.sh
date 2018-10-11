#!/bin/bash
# 
# How many nodes are needed
#SBATCH --nodes=2
#
#SBATCH --ntasks=2
#
# How many MPI tasks per node
#SBATCH --ntasks-per-node=1  
#
# How many physical CPU’s per task
#SBATCH --cpus-per-task=8
#  
# How long the job is anticipated to run
#SBATCH --time=10:00:00
# 
# The name of the job
#SBATCH --job-name=self_testing 
#
# ///////////////////////////////////////////////////////////////////////////////// #
#
# In order to use this script the following environment variables need to be set
#
# SELFMODDIR - the directory where the SELF-Fluids modulefiles can be found
#
# SELFDIR - the path to the directory where SELF-Fluids is installed
#
# EX_DIR - the directory containing the example/test case to be run (this directory
#          contains the self.equations and runtime.params files)
#
# TEST_CASE - the name of the testcase - also the folder name under ${SELFDIR}/examples/
#
# RUN_DIR - the head directory where model output should be contained. Directory
#           $RUN_DIR/$TEST_CASE/$SELFMOD contains
#
# COMPILER - the environment module for the compiler used to build SELF-Fluids
#
# MPIMOD - the environment module for the MPI flavor used to build SELF-Fluids
#
# HDF5MOD - the environment module for the HDF5 flavor used to build SELF-Fluids
#
# SELFMOD - the environment module for loading a version of SELF-Fluids for testing
#
# OMP_NUM_THREADS - the number of OpenMP threads to use (if OpenMP is enabled)
#
# ///////////////////////////////////////////////////////////////////////////////// #

export SELFMODDIR="${HOME}/apps/modulefiles"
export SELFDIR="${HOME}/testing/SELF-Fluids"
export OMP_NUM_THREADS=8
export RUN_DIR="${SELFDIR}/testing"

job_wall_time="5:00:00"

test_cases=("boundarylayer" \
            "thermalbubble")

self_mods=("pgi/18.4/serial" \
           "pgi/18.4/cuda" \
           "pgi/18.4/openmp" \
           "pgi/18.4/openmpi_3.1.2" \
           "pgi/18.4/cuda+openmpi_3.1.2" \
           "pgi/18.4/openmp+openmpi_3.1.2" \
           "gcc/8.1.0/serial" \
           "gcc/8.1.0/openmp" \
           "gcc/8.1.0/openmpi_3.1.2" \
           "gcc/8.1.0/openmp+openmpi_3.1.2")

n_nodes=(1 1 1 2 2 2 1 1 2 2)
n_tasks=(1 1 1 16 2 2 1 1 16 2)
n_tasks_per_node=(1 1 1 8 1 1 1 1 8 1)
cpus_per_task=(1 1 1 1 1 1 1 1 1 1)
                   
compiler_mods=("pgi/18.4" \
               "pgi/18.4" \
               "pgi/18.4" \
               "pgi/18.4" \
               "pgi/18.4" \
               "pgi/18.4" \
               "gcc/8.1.0" \
               "gcc/8.1.0" \
               "gcc/8.1.0" \
               "gcc/8.1.0")

mpi_mods=("" \
          "" \
          "" \
          "openmpi/3.1.2/pgi/18.4" \
          "openmpi/3.1.2/pgi/18.4" \
          "openmpi/3.1.2/pgi/18.4" \
          "" \
          "" \
          "openmpi/3.1.2/gcc/8.1.0" \
          "openmpi/3.1.2/gcc/8.1.0")

hdf5_mods=("hdf5/1.10.2/serial/pgi/18.4"\
           "hdf5/1.10.2/serial/pgi/18.4"\
           "hdf5/1.10.2/serial/pgi/18.4"\
           "hdf5/1.10.2/parallel/pgi/18.4/openmpi/3.1.2"\
           "hdf5/1.10.2/parallel/pgi/18.4/openmpi/3.1.2"\
           "hdf5/1.10.2/parallel/pgi/18.4/openmpi/3.1.2"\
           "hdf5/1.10.2/serial/gcc/8.1.0"\
           "hdf5/1.10.2/serial/gcc/8.1.0"\
           "hdf5/1.10.2/parallel/gcc/8.1.0/openmpi/3.1.2"\
           "hdf5/1.10.2/parallel/gcc/8.1.0/openmpi/3.1.2")

n_demos=${#test_cases[@]}
n_builds=${#self_mods[@]}

for (( j=0; j<${n_demos}; j++ ));
do
  for (( i=0; i<${n_builds}; i++ ));
  do
    export TEST_CASE=${test_cases[$j]}
    export COMPILER=${compiler_mods[$i]}
    export MPIMOD=${mpi_mods[$i]}
    export HDF5MOD=${hdf5_mods[$i]}
    export SELFMOD=${self_mods[$i]}

    sbatch --export=ALL,\
           --nodes=${n_nodes[$i]} \
           --ntasks=${n_tasks[$i]} \
           --ntasks-per-node=${n_tasks_per_node[$i]} \
           --cpus-per-task=${cpus_per_task[$i]} \
           --time=${job_wall_time} \
           --job-name="self_testing.${i}.${j}" \
           self.modelrun.slurm

  done

done

unset SELFMODDIR
unset SELFDIR
unset RUN_DIR
unset OMP_NUM_THREADS
unset TEST_CASE
unset COMPILER
unset MPIMOD
unset HDF5MOD
unset SELFMOD
