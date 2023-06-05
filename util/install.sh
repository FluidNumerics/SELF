#!/bin/bash
#
# Environment Variables
#
# VIEW - The path to the spack environment view.
# SELF_PREFIX - The path to install SELF. Defaults to $VIEW
# GPU_TARGET - GPU microarchitecture code to build for. Defaults to gfx900 (AMD MI25)
# PREC - Floating point precision to build with
# SELF_FFLAGS - compiler flags to build SELF


: "${VIEW:=/opt/view/self}"
: "${SELF_PREFIX:=/opt/self}"

: "${GPU_TARGET:=gfx906}"
: "${PREC:=single}"
: "${SELF_FFLAGS:="-cpp -pg -g -O0 -C -Wall -fbounds-check -fbacktrace --coverage -ffpe-trap=invalid,zero,overflow"}"

if [ "${GPU_TARGET}" == "sm"* ];then
    HIP_PLATFORM=nvidia
else
    HIP_PLATFORM=amd
fi

export SELF_JSONF_LIBS="-L${VIEW}/lib64 -ljsonfortran"
export SELF_JSONF_INC="-I${VIEW}/include"
export SELF_FEQPARSE_LIBS="-L${VIEW}/lib -lfeqparse"
export SELF_FEQPARSE_INC="-I${VIEW}/include"
export SELF_FLAP_LIBS="-L${VIEW}/lib64 -lFLAP -lFACE -lPENF"
export SELF_FLAP_INC="-I${VIEW}/include/FLAP -I${VIEW}/include/PENF -I${VIEW}/include/FACE"
export SELF_HDF5_LIBS="-L${VIEW}/lib -lhdf5_fortran -lhdf5 -lz -lm"
export SELF_HDF5_INC="-I${VIEW}/include/shared"

export SELF_PREFIX=${SELF_PREFIX}
export PREC=${PREC}
export SELF_FFLAGS=${SELF_FFLAGS}

make clean
rm -rf ${SELF_PREFIX}
make
