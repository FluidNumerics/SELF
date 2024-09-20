#!/usr/bin/env -S bash -e

WORKSPACE_ROOT=/scratch/$(whoami)/workspace
#BUILD_TYPE=coverage
BUILD_TYPE=release
SRC_DIR=$(pwd)
BUILD_DIR=/scratch/$(whoami)/build
ENABLE_GPU=ON
ENABLE_DOUBLE_PRECISION=ON
ENABLE_MULTITHREADING=OFF
NTHREADS=4
GPU_ARCH=gfx90a
GCOV=gcov-12

module purge
module load cmake/3.29.6
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
      -DSELF_ENABLE_MULTITHREADING=${ENABLE_MULTITHREADING} \
      -DSELF_MULTITHREADING_NTHREADS=${NTHREADS} \
      -DSELF_ENABLE_DOUBLE_PRECISION=${ENABLE_DOUBLE_PRECISION} \
      -DCMAKE_HIP_ARCHITECTURE=${GPU_ARCH} \
      ${SRC_DIR}

make VERBOSE=1 || exit 1
make install

# # Set WORKSPACE for tests that require input mesh
# # We use WORKSPACE so that we are consistent with 
# # what we do for the superci tests
export WORKSPACE=${SRC_DIR}

if [ "$BUILD_TYPE" = "coverage" ]; then
     lcov --capture \
          --initial \
          --directory ${BUILD_DIR}/src \
          --gcov=${GCOV} \
          --output-file ${SRC_DIR}/initial.info
fi

# # Run ctests
# ctest --test-dir ${BUILD_DIR}

if [ "$BUILD_TYPE" = "coverage" ]; then
     # Compile coverage information
     lcov --capture \
          --directory ${BUILD_DIR}/src \
          --gcov=${GCOV} \
          --output-file ${SRC_DIR}/ctest-capture.info

     lcov --add-tracefile ${SRC_DIR}/initial.info \
          --add-tracefile ${SRC_DIR}/ctest-capture.info \
          --gcov=${GCOV} \
          --output-file ${SRC_DIR}/coverage.info
          
     # Generate summary
     lcov --summary ${SRC_DIR}/coverage.info
fi
