#!/bin/bash


export SELF_BUILD_DIR=/build/self
export SELF_PREFIX=/opt/self

lcov --no-external \
     --capture \
     --initial \
     --directory ${SELF_BUILD_DIR} \
     --output-file /tmp/lcov_base.info

for file in $(ls $SELF_PREFIX/test/*); do

  echo "Running Test : $file"
  $file
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "$file failed!"
    exit $rc
  fi

  # GPU accelerated
  $file --gpu
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "$file failed!"
    exit $rc
  fi

done

lcov --no-external \
     --capture \
     --directory ${SELF_BUILD_DIR} \
     --output-file /tmp/lcov_test.info

lcov --add-tracefile /tmp/lcov_base.info \
     --add-tracefile /tmp/lcov_test.info \ 
     --output-file /build/lcov.info
