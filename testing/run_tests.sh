#!/bin/bash

export SELFDIR=$(pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')

source ${SELFDIR}/build/SELF_environment_settings

# ---------------------------------------------------------------------- #

export CUDA=no
export DOUBLE_PRECISION=yes

make clean
make spectral_tests

./spectral_tests > cpu_spectral_tests.double.txt 

# ---------------------------------------------------------------------- #

export DOUBLE_PRECISION=no

make clean
make spectral_tests

./spectral_tests > cpu_spectral_tests.single.txt 

# ---------------------------------------------------------------------- #

export CUDA=yes
export DOUBLE_PRECISION=yes

make clean
make spectral_tests

./spectral_tests > gpu_spectral_tests.double.txt 

# ---------------------------------------------------------------------- #

export CUDA=yes
export DOUBLE_PRECISION=no

make clean
make spectral_tests

./spectral_tests > gpu_spectral_tests.single.txt 

# ---------------------------------------------------------------------- #

diff cpu_spectral_tests.single.txt references/cpu_spectral_tests.single.reference
diff gpu_spectral_tests.single.txt references/gpu_spectral_tests.single.reference

diff cpu_spectral_tests.double.txt references/cpu_spectral_tests.double.reference
diff gpu_spectral_tests.double.txt references/gpu_spectral_tests.double.reference
