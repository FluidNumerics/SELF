#!/bin/bash
# 
#
# ///////////////////////////////////////////////////////////////////////////////// #
#
# In order to use this script the following environment variables need to be set
#
# SELFMODDIR - the directory where the SELF-Fluids modulefiles can be found
#
# SELFDIR - the path to the directory where SELF-Fluids is installed
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
# Additionally, the SLURM variables must be set at command line : 
#
# How many nodes are needed
#  --nodes
#
#  --ntasks
#
# How many MPI tasks per node
#  --ntasks-per-node  
#
# How many physical CPU’s per task
#  --cpus-per-task
#  
# How long the job is anticipated to run
#  --time
# 
# The name of the job
#  --job-name
#
# ///////////////////////////////////////////////////////////////////////////////// #



pwd; hostname
compute_host=$(hostname)
module purge
module use ${SELFMODDIR}

run_model(){

  if [ ! -d "$1" ]; then

    mkdir -p $1

  fi

  cd $1

  echo " Compute node host(s)   : ${SLURM_JOB_NODELIST}" > logfile.txt
  echo " Slurm Job ID           : ${SLURM_JOB_ID}" >> logfile.txt
  echo " Slurm submission node  : ${SLURM_SUBMIT_HOST}" >> logfile.txt
  echo " Run directory          : $1" >> logfile.txt
  echo " SELF modules directory : ${SELFMODDIR} " >> logfile.txt
  echo " SELF install directory : ${SELFDIR} " >> logfile.txt
  echo " SELF test case         : ${TEST_CASE} " >> logfile.txt
  echo " SELF run directory     : ${RUN_DIR} " >> logfile.txt
  echo " Compiler               : ${COMPILER} " >> logfile.txt
  echo " MPI module             : ${MPIMOD} " >> logfile.txt
  echo " HDF5 module            : ${HDF5MOD} " >> logfile.txt
  echo " SELF module            : ${SELFMOD} " >> logfile.txt
  echo " Number of threads      : ${OMP_NUM_THREADS} " >> logfile.txt
  echo " Execution command      : sfluid --param-file ${SELFDIR}/examples/${TEST_CASE}/runtime.params --equation-file ${SELFDIR}/examples/${TEST_CASE}/self.equations " >> logfile.txt

  sfluid --param-file ${SELFDIR}/examples/${TEST_CASE}/runtime.params \
         --equation-file ${SELFDIR}/examples/${TEST_CASE}/self.equations

}

run_mpi_model(){

  if [ ! -d "$1" ]; then

    mkdir -p $1

  fi

  cd $1
  
  echo " Compute node host(s)   : ${SLURM_JOB_NODELIST}" > logfile.txt
  echo " Slurm Job ID           : ${SLURM_JOB_ID}" >> logfile.txt
  echo " Slurm submission node  : ${SLURM_SUBMIT_HOST}" >> logfile.txt
  echo " Run directory          : $1" >> logfile.txt
  echo " SELF modules directory : ${SELFMODDIR} " >> logfile.txt
  echo " SELF install directory : ${SELFDIR} " >> logfile.txt
  echo " SELF test case         : ${TEST_CASE} " >> logfile.txt
  echo " SELF run directory     : ${RUN_DIR} " >> logfile.txt
  echo " Compiler               : ${COMPILER} " >> logfile.txt
  echo " MPI module             : ${MPIMOD} " >> logfile.txt
  echo " HDF5 module            : ${HDF5MOD} " >> logfile.txt
  echo " SELF module            : ${SELFMOD} " >> logfile.txt
  echo " Number of threads      : ${OMP_NUM_THREADS} " >> logfile.txt
  echo " Execution command      : mpirun -np ${SLURM_NTASKS} sfluid --param-file ${SELFDIR}/examples/${TEST_CASE}/runtime.params --equation-file ${SELFDIR}/examples/${TEST_CASE}/self.equations " >> logfile.txt

  mpirun -np ${SLURM_NTASKS} sfluid --param-file ${SELFDIR}/examples/${TEST_CASE}/runtime.params \
                                    --equation-file ${SELFDIR}/examples/${TEST_CASE}/self.equations

}

if [ -z "$MPIMOD"];
then

  module load ${COMPILER} ${HDF5MOD} self-fluids/dev/${SELFMOD}
  run_model ${RUN_DIR}/${TEST_CASE}/${SELFMOD}
  module unload ${COMPILER} ${HDF5MOD} self-fluids/dev/${SELFMOD}

else

  module load ${COMPILER} ${MPIMOD} ${HDF5MOD} self-fluids/dev/${SELFMOD}
  run_mpi_model ${RUN_DIR}/${TEST_CASE}/${SELFMOD}
  module unload ${COMPILER} ${MPIMOD} ${HDF5MOD} self-fluids/dev/${SELFMOD}

fi

