#!/bin/bash


export SELF_BUILD_DIR=/build
export SELF_PREFIX=/opt/self


moveGCNOFiles(){

  for file in $(ls $SELF_BUILD_DIR/src/*.gcno); do
    dest=$(echo $file | sed 's/\.f//g')
    cp $file $dest
  done

  for file in $(ls $SELF_BUILD_DIR/test/*.gcno); do
    dest=$(echo $file | sed 's/\.f90//g')
    cp $file $dest
  done

}

moveGCDAFiles(){
# This method is used to apply the correct file extension
# for gcovr to be able to match up the gcda file to the 
# source code.
#
# This is only needed because the SELF build system applies
# the .f.o extension to object files created from our .f90
# files.

  for file in $(ls $SELF_BUILD_DIR/src/*.gcda); do
    dest=$(echo $file | sed 's/\.f//g')
    cp $file $dest
  done

  for file in $(ls $SELF_BUILD_DIR/test/*.gcda); do
    dest=$(echo $file | sed 's/\.f90//g')
    cp $file $dest
  done

}

cleanupGCDAFiles(){
  rm $SELF_BUILD_DIR/test/*.gcda
  rm $SELF_BUILD_DIR/src/*.gcda
}

moveGCNOFiles

ls $SELF_PREFIX/test/*

# Reset counters
lcov --directory $SELF_BUILD_DIR --zerocounters
lcov --no-external --capture --initial --directory $SELF_BUILD_DIR --output-file lcov_base.info

for file in $(ls $SELF_PREFIX/test/*); do

  #covFile=$(echo $file | awk -F "/" '{print $NF}')

  echo "Running Test : $file"
  $file
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "$file failed!"
    exit $rc
  fi
  lcov --directory $SELF_BUILD_DIR --capture --output-file lcov.info

  # GPU accelerated
  $file --gpu
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "$file failed!"
    exit $rc
  fi
  lcov --directory $SELF_BUILD_DIR --capture --output-file lcov.info
  

  #moveGCDAFiles

  # Create a coverage file for this test.
  #moveGCDAFiles
  #gcovr -r ${SELF_BUILD_DIR} --json -o "${SELF_BUILD_DIR}/$covFile.json"
  #tomerge+="-a ${SELF_BUILD_DIR}/$covFile.json "
  #cleanupGCDAFiles

  

done

#gcovr $tomerge
#gcovr $tomerge -x -o ${SELF_BUILD_DIR}/coverage.xml
# Process coverage

[ -d "/workspace" ] && cp ${SELF_BUILD_DIR}/lcov.info /workspace
