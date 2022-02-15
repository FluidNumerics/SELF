#!/bin/bash


SELF_BUILD_DIR=/build
SELF_INSTALL_DIR=/opt/self

moveCovFiles(){

  for file in $(ls $SELF_BUILD_DIR/src/*.gcno); do
    dest=$(echo $file | sed 's/\.f//g')
    echo "$file -> $dest"
    cp $file $dest
  done

  for file in $(ls $SELF_BUILD_DIR/src/*.gcda); do
    dest=$(echo $file | sed 's/\.f//g')
    echo "$file -> $dest"
    cp $file $dest
  done

  for file in $(ls $SELF_BUILD_DIR/examples/*.gcno); do
    dest=$(echo $file | sed 's/\.f90//g')
    echo "$file -> $dest"
    cp $file $dest
  done

  for file in $(ls $SELF_BUILD_DIR/examples/*.gcda); do
    dest=$(echo $file | sed 's/\.f90//g')
    echo "$file -> $dest"
    cp $file $dest
  done

}

cleanupCovFiles(){
  rm $SELF_BUILD_DIR/examples/*.gcda
  rm $SELF_BUILD_DIR/src/*.gcda
}

tomerge=""
for file in $(ls $SELF_INSTALL_DIR/examples/*); do

  echo "Running Test : $file"
  $file

  moveCovFiles

#      for file in $(ls /build/src/*.f90); do
#        gcov --function-summaries \
#             --human-readable \
#             --demangled-names \
#             $file
#      done

  covFile=$(echo $file | awk -F "/" '{print $NF}')
  gcovr -r /build --json -o "/build/$covFile.json"
  tomerge+="-a /build/$covFile.json "

  cleanupCovFiles

done

  gcovr $tomerge
  gcovr $tomerge -x -o /build/coverage.xml
