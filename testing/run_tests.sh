#!/bin/bash

export SELFDIR=$(pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')

source ${SELFDIR}/build/SELF_environment_settings

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
