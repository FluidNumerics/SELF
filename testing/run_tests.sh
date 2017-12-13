#!/bin/bash

export SELFDIR=$(pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')

# Compilers
export FC=pgfortran
export MPIFC=mpif90

# Support
export DEBUG=no
export TIMING=yes

# Precision
export DOUBLE_PRECISION=yes

# Parallelization schemes
export CUDA=no
export GPU_ARCH=cc60

export MPI=no
export MPI_LIB=
export MPI_INC=

export OMP=no

# ---------------------------------------------------------------------- #

export CUDA=no

make clean
make spectral_tests

./spectral_tests > cpu_spectral_tests.txt 

# ---------------------------------------------------------------------- #

export CUDA=yes

make clean
make spectral_tests

./spectral_tests > gpu_spectral_tests.txt 

# ---------------------------------------------------------------------- #

diff cpu_spectral_tests.txt cpu_spectral_tests.reference
diff gpu_spectral_tests.txt gpu_spectral_tests.reference
