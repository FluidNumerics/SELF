#!/bin/bash

export SELFDIR=$(pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')

source ${SELFDIR}/build/SELF_environment_settings

# ---------------------------------------------------------------------- #

export CUDA=no
export DOUBLE_PRECISION=yes

make -s clean
make -s spectral_tests

echo "spectral_tests : serial cpu double precision"

./spectral_tests > cpu_spectral_tests.double.txt 

# ---------------------------------------------------------------------- #

export DOUBLE_PRECISION=no

make -s clean
make -s spectral_tests

echo "spectral_tests : serial cpu single precision"

./spectral_tests > cpu_spectral_tests.single.txt 

# ---------------------------------------------------------------------- #

export CUDA=yes
export DOUBLE_PRECISION=yes

make -s clean
make -s spectral_tests

echo "spectral_tests : gpu double precision"

./spectral_tests > gpu_spectral_tests.double.txt 

# ---------------------------------------------------------------------- #

export CUDA=yes
export DOUBLE_PRECISION=no

make -s clean
make -s spectral_tests

echo "spectral_tests : gpu single precision"

./spectral_tests > gpu_spectral_tests.single.txt 

# ---------------------------------------------------------------------- #

diff -q cpu_spectral_tests.single.txt references/cpu_spectral_tests.single.reference
diff -q gpu_spectral_tests.single.txt references/gpu_spectral_tests.single.reference

diff -q cpu_spectral_tests.double.txt references/cpu_spectral_tests.double.reference
diff -q gpu_spectral_tests.double.txt references/gpu_spectral_tests.double.reference
