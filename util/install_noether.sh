#!/usr/bin/env -S bash -e

OMP_TARGET=none
#WORKSPACE_ROOT=$HOME/.local/workspace/
WORKSPACE_ROOT=/scratch/joe/workspace
BUILD_TYPE=coverage
SRC_DIR=$(pwd)
#BUILD_DIR=${SRC_DIR}/build
BUILD_DIR=/scratch/joe/build

module load gcc/13.2.0
module load openmpi/5.0.1 hdf5/1.12.3 feq-parse/2.2.2

# Clean out any old builds
rm -rf ${BUILD_DIR}
rm -rf ${WORKSPACE_ROOT}/*

mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

FC=gfortran \
cmake -DCMAKE_PREFIX_PATH=${ROCM_PATH} \
      -DOMP_TARGET=${OMP_TARGET} \
      -DCMAKE_INSTALL_PREFIX=${WORKSPACE_ROOT}/opt/self \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      ${SRC_DIR}

make VERBOSE=1 || exit 1

make install


# Set WORKSPACE for tests that require input mesh
# We use WORKSPACE so that we are consistent with 
# what we do for the superci tests
export WORKSPACE=${SRC_DIR}

# Initialize coverage
mkdir -p ${WORKSPACE_ROOT}/tmp/
lcov --no-external \
      --capture \
      --initial \
      --directory ${SRC_DIR} \
      --exclude '*/test/*' \
      --exclude '*/example/*' \
      --output-file ${WORKSPACE_ROOT}/tmp/lcov_base.info

# Run ctests
ctest --test-dir ${BUILD_DIR}

# Compile coverage information
lcov --no-external \
    --capture \
    --directory ${SRC_DIR} \
    --exclude '*/test/*' \
    --exclude '*/example/*' \
    --output-file ${WORKSPACE_ROOT}/tmp/lcov_test.info

lcov --add-tracefile ${WORKSPACE_ROOT}/tmp/lcov_base.info \
     --add-tracefile ${WORKSPACE_ROOT}/tmp/lcov_test.info \
     --exclude '*/test/*' \
     --exclude '*/example/*' \
     --output-file ${SRC_DIR}/lcov.info
