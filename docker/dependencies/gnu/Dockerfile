# 
# HPC Base image
# 
# Contents:
#   CUDA version 10.0
#   HDF5 version 1.10.5
#   GNU compilers version 8.3
#   Python 2 and 3 (upstream)
# 

FROM nvidia/cuda:10.1-devel-centos7 AS devel

# Python
RUN yum install -y \
        python2 \
        python3 && \
    rm -rf /var/cache/yum/*

# GNU compiler
RUN yum install -y centos-release-scl && \
    yum install -y \
        devtoolset-8-gcc \
        devtoolset-8-gcc-c++ \
        devtoolset-8-gcc-gfortran && \
    rm -rf /var/cache/yum/*
ENV PATH=/opt/rh/devtoolset-8/root/usr/bin:$PATH

# HDF5 version 1.10.5
RUN yum install -y \
        bzip2 \
        file \
        make \
        wget \
        zlib-devel && \
    rm -rf /var/cache/yum/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.bz2 && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/hdf5-1.10.5.tar.bz2 -C /var/tmp -j && \
    cd /var/tmp/hdf5-1.10.5 &&  CC=gcc CXX=g++ F77=gfortran F90=gfortran FC=gfortran ./configure --prefix=/usr/local/hdf5 --enable-cxx --enable-fortran && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/hdf5-1.10.5.tar.bz2 /var/tmp/hdf5-1.10.5
ENV HDF5_DIR=/usr/local/hdf5 \
    LD_LIBRARY_PATH=/usr/local/hdf5/lib:$LD_LIBRARY_PATH \
    PATH=/usr/local/hdf5/bin:$PATH

RUN yum install -y cmake

RUN mkdir -p /var/tmp && \
    wget -q -nc --no-check-certificate -P /var/tmp http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz && \
    tar -xzf /var/tmp/metis-5.1.0.tar.gz -C /var/tmp && \
    cd /var/tmp/metis-5.1.0 && \
    make config prefix=/usr/local/metis && \
    make install

ENV LIB_METIS=/usr/local/metis/lib/libmetis.a
