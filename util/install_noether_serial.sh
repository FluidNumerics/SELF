#!/usr/bin/env -S bash -e

WORKSPACE_ROOT=/scratch/$(whoami)/workspace
BUILD_TYPE=coverage
SRC_DIR=$(pwd)
BUILD_DIR=/scratch/$(whoami)/build
ENABLE_GPU=OFF
GPU_ARCH=gfx90a

module purge
module load gcc/12.3.0
module load openmpi/5.0.3 
module load hdf5/1.14.3 feq-parse/2.2.2
module load rocm/6.0.2
module list

# Clean out any old builds
rm -rf ${BUILD_DIR}
rm -rf ${WORKSPACE_ROOT}/*

mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

FC=gfortran \
CXX=hipcc \
cmake -DCMAKE_PREFIX_PATH=${ROCM_PATH} \
      -DCMAKE_INSTALL_PREFIX=${WORKSPACE_ROOT}/opt/self \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      -DSELF_ENABLE_GPU=${ENABLE_GPU} \
      -DCMAKE_HIP_ARCHITECTURE=${GPU_ARCH} \
      ${SRC_DIR}

make || exit 1
make install


# # Set WORKSPACE for tests that require input mesh
# # We use WORKSPACE so that we are consistent with 
# # what we do for the superci tests
export WORKSPACE=${SRC_DIR}

# # Initialize coverage
# mkdir -p ${WORKSPACE_ROOT}/tmp/
# lcov --no-external \
#       --capture \
#       --initial \
#       --directory ${SRC_DIR} \
#       --exclude '*/test/*' \
#       --exclude '*/example/*' \
#       --output-file ${WORKSPACE_ROOT}/tmp/lcov_base.info

# # Run ctests
ctest --test-dir ${BUILD_DIR}

# # Compile coverage information
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
