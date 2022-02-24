#!/bin/bash
#
# Environment Variables
#
# VIEW - The path to the spack environment view.
# SELF_PREFIX - The path to install SELF. Defaults to $VIEW
# GPU_TARGET - GPU microarchitecture code to build for. Defaults to gfx900 (AMD MI100)
# PREC - Floating point precision to build with
# SELF_FFLAGS - compiler flags to build SELF


: "${VIEW:=${HOME}/view/self}"
: "${SELF_PREFIX:=${VIEW}}"

: "${GPU_TARGET:=gfx900}"
: "${PREC:=double}"
: "${SELF_FFLAGS:="-cpp -pg -g -O3"}"

if [ "${GPU_TARGET}" == "sm"* ];then
    HIP_PLATFORM=nvidia
else
    HIP_PLATFORM=amd
fi

export SELF_JSONF_LIBS="-L${VIEW}/lib -ljsonfortran"
export SELF_JSONF_INC="-I${VIEW}/include"
export SELF_FEQPARSE_LIBS="-L${VIEW}/lib -lfeqparse"
export SELF_FEQPARSE_INC="-I${VIEW}/include"
export SELF_FLAP_LIBS="-L${VIEW}/lib -lFLAP -lFACE -lPENF"
export SELF_FLAP_INC="-I${VIEW}/include/FLAP -I${VIEW}/include/PENF -I${VIEW}/include/FACE"
export SELF_HDF5_LIBS="-L${VIEW}/lib -lhdf5_fortran -lhdf5 -lz -lm"
export SELF_HDF5_INC="-I${VIEW}/include/shared"



rm -rf ${SELF_PREFIX}

spack env activate -d ./env/
make
