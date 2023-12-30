#!/bin/bash

GPU_TARGET=gfx90a
WORKSPACE=$HOME/.local/workspace/
BUILD_TYPE=coverage
SRC_DIR=$(pwd)
BUILD_DIR=${SRC_DIR}/build

set -x
module load gcc/13.2.0
module load openmpi hdf5 feq-parse

# Clean out any old builds
rm -rf ${BUILD_DIR}
rm -rf ${WORKSPACE}/*

mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

FC=gfortran \
FFLAGS="-DDOUBLE_PRECISION" \
cmake -DCMAKE_PREFIX_PATH=/opt/rocm \
      -DCMAKE_HIP_ARCHITECTURES=${GPU_TARGET} \
      -DCMAKE_INSTALL_PREFIX=${WORKSPACE}/opt/self \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      ${SRC_DIR}
make VERBOSE=1 || exit 1
make install


# Set SELF_PREFIX for tests that require input mesh
export SELF_PREFIX=${SRC_DIR}

# Initialize coverage
mkdir -p ${WORKSPACE}/tmp/
lcov --no-external \
      --capture \
      --initial \
      --directory ${SRC_DIR} \
      --output-file ${WORKSPACE}/tmp/lcov_base.info

# Run ctests
ctest --test-dir ${BUILD_DIR}/test

# Compile coverage information
lcov --no-external \
    --capture \
    --directory ${SRC_DIR} \
    --output-file ${WORKSPACE}/tmp/lcov_test.info

lcov --add-tracefile ${WORKSPACE}/tmp/lcov_base.info \
     --add-tracefile ${WORKSPACE}/tmp/lcov_test.info \
     --output-file ${SRC_DIR}/lcov.info
