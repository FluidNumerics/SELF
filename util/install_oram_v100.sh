#!/usr/bin/env -S bash -e

OMP_TARGET=multicore
WORKSPACE_ROOT=/tmp/$(whoami)/workspace
BUILD_TYPE=coverage
SRC_DIR=$(pwd)
BUILD_DIR=/tmp/$(whoami)/build

module load cuda/12.5.0 cmake/3.29.6
module load openmpi/5.0.3 hdf5/1.12.3 feq-parse/2.2.2

# Clean out any old builds
rm -rf ${BUILD_DIR}
rm -rf ${WORKSPACE_ROOT}/*

mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

cmake -DOMP_TARGET=${OMP_TARGET} \
      -DCMAKE_INSTALL_PREFIX=${WORKSPACE_ROOT}/opt/self \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      ${SRC_DIR}

make VERBOSE=1 || exit 1

make install

export OMP_TARGET_OFFLOAD=DISABLED # Disable GPU offloading
export OMP_NUM_THREADS=2
# Set WORKSPACE for tests that require input mesh
# We use WORKSPACE so that we are consistent with 
# what we do for the superci tests
export WORKSPACE=${SRC_DIR}

# Initialize coverage
mkdir -p ${WORKSPACE_ROOT}/tmp/
# lcov --no-external \
#       --capture \
#       --initial \
#       --directory ${SRC_DIR} \
#       --exclude '*/test/*' \
#       --exclude '*/example/*' \
#       --output-file ${WORKSPACE_ROOT}/tmp/lcov_base.info

# Run ctests
ctest --test-dir ${BUILD_DIR}

# Compile coverage information
# lcov --no-external \
#     --capture \
#     --directory ${SRC_DIR} \
#     --exclude '*/test/*' \
#     --exclude '*/example/*' \
#     --output-file ${WORKSPACE_ROOT}/tmp/lcov_test.info

# lcov --add-tracefile ${WORKSPACE_ROOT}/tmp/lcov_base.info \
#      --add-tracefile ${WORKSPACE_ROOT}/tmp/lcov_test.info \
#      --exclude '*/test/*' \
#      --exclude '*/example/*' \
#      --output-file ${SRC_DIR}/lcov.info
