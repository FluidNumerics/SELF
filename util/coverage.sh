#!/bin/bash
 

export SELF_BUILD_DIR=/build/
export SELF_PREFIX=/opt/self

lcov --no-external \
     --capture \
     --initial \
     --directory ${SELF_BUILD_DIR} \
     --output-file /tmp/lcov_base.info

# Load dependency packages
. /etc/profile.d/z10_spack.sh

# Run the tests
ctest --test-dir ${SELF_BUILD_DIR}/build/test

rc=$?
if [ $rc -ne 0 ]; then
  echo "SELF tests failed!"
  exit $rc
fi

lcov --no-external \
     --capture \
     --directory ${SELF_BUILD_DIR} \
     --output-file /tmp/lcov_test.info

lcov --add-tracefile /tmp/lcov_base.info \
     --add-tracefile /tmp/lcov_test.info \ 
     --output-file /build/lcov.info
