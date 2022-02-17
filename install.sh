#!/bin/bash


export HIP_PLATFORM=amd
export SELF_PREFIX=/home/joe/install/self
export PREC=double
export GPU_TARGET=gfx900
export VIEW=/home/joe/view/self

export SELF_JSONF_LIBS="-L${VIEW}/lib -ljsonfortran"
export SELF_JSONF_INC="-I${VIEW}/include"

export SELF_FEQPARSE_LIBS="-L${VIEW}/lib -lfeqparse"
export SELF_FEQPARSE_INC="-I${VIEW}/include"

export SELF_FLAP_LIBS="-L${VIEW}/lib -lFLAP -lFACE -lPENF"
export SELF_FLAP_INC="-I${VIEW}/include/FLAP -I${VIEW}/include/PENF -I${VIEW}/include/FACE"

export SELF_HDF5_LIBS="-L${VIEW}/lib -lhdf5_fortran -lhdf5 -lz -lm"
export SELF_HDF5_INC="-I${VIEW}/include/shared"

export SELF_FFLAGS="-cpp -pg -g -O3 -C -Wall -fbounds-check -fbacktrace --coverage -ffpe-trap=invalid,zero,overflow"


rm -rf ${SELF_PREFIX}

spack env activate -d ./env/
make
