#!/usr/bin/env bash
#
# Coverage build and test suite for the AMD MI210 (HIP, gfx90a) CI job.
#
# Runs inside the higherordermethods/selfish container, with the repository
# bind-mounted at /workspace. Driven by .github/workflows/amd-mi210-gpu-tests.yml.
#
# This lives in a file rather than being piped to `bash -s` from the workflow
# on purpose. When the script is read from stdin, the first MPI test to run
# under ctest consumes the remainder of it - Open MPI forwards stdin to rank 0 -
# so everything after the ctest line silently never executes and the step still
# exits 0.

set -euo pipefail

source /opt/spack-environment/activate.sh

# Unlike pyxis, which ran the container as the submitting user, the container
# here is uid 0. Open MPI refuses to launch as root without an explicit
# acknowledgement; 5.x reads the PRTE_ names, 4.x the OMPI_ ones.
export PRTE_ALLOW_RUN_AS_ROOT=1
export PRTE_ALLOW_RUN_AS_ROOT_CONFIRM=1
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

mkdir -p /workspace/build
cd /workspace/build

FC=gfortran cmake \
      -DCMAKE_INSTALL_PREFIX=/workspace/opt/self \
      -DCMAKE_BUILD_TYPE="coverage" \
      -DSELF_ENABLE_GPU=ON \
      -DSELF_GPU_BACKEND=HIP \
      -DSELF_ENABLE_TESTING=ON \
      -DCMAKE_HIP_ARCHITECTURES="gfx90a" \
      -DGPU_TARGETS="gfx90a" \
      -DSELF_ENABLE_EXAMPLES=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
      -DCMAKE_VERBOSE_MAKEFILE=ON \
      /workspace/

make -j"$(nproc)" VERBOSE=1

lcov --capture --initial \
     --directory /workspace/build/src/ \
     --output-file /workspace/initial.info

export WORKSPACE=/workspace
ctest --test-dir /workspace/build --output-on-failure

lcov --capture \
     --directory /workspace/build/src/ \
     --output-file /workspace/ctest-capture.info

lcov --add-tracefile /workspace/initial.info \
     --add-tracefile /workspace/ctest-capture.info \
     --output-file /workspace/coverage.info

# lcov records source paths as seen in here (/workspace/src/...). Strip the
# mount point so Codecov matches them against the repository tree.
sed -i 's|^SF:/workspace/|SF:|' /workspace/coverage.info

test -s /workspace/coverage.info
echo "coverage.info: $(wc -l < /workspace/coverage.info) lines"
