# 
# HPC Base image
# 
# Contents:
#   CUDA version 10.0
#   HDF5 version 1.10.5
#   GNU compilers version 8.3
#   MVAPICH2 version 2.3.1
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

# Mellanox OFED version 4.6-1.0.1.1
RUN yum install -y \
        findutils \
        libnl \
        libnl3 \
        numactl-libs \
        wget && \
    rm -rf /var/cache/yum/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp http://content.mellanox.com/ofed/MLNX_OFED-4.6-1.0.1.1/MLNX_OFED_LINUX-4.6-1.0.1.1-rhel7.2-x86_64.tgz && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/MLNX_OFED_LINUX-4.6-1.0.1.1-rhel7.2-x86_64.tgz -C /var/tmp -z && \
    find /var/tmp/MLNX_OFED_LINUX-4.6-1.0.1.1-rhel7.2-x86_64 -regextype posix-extended -type f -regex ".*(libibmad|libibmad-devel|libibumad|libibumad-devel|libibverbs|libibverbs-devel|libibverbs-utils|libmlx4|libmlx4-devel|libmlx5|libmlx5-devel|librdmacm|librdmacm-devel)-[0-9].*x86_64.rpm" -not -path "*UPSTREAM*" -exec rpm --install {} + && \
    rm -rf /var/tmp/MLNX_OFED_LINUX-4.6-1.0.1.1-rhel7.2-x86_64.tgz /var/tmp/MLNX_OFED_LINUX-4.6-1.0.1.1-rhel7.2-x86_64

# MVAPICH2 version 2.3.1
RUN yum install -y \
        byacc \
        file \
        make \
        openssh-clients \
        wget && \
    rm -rf /var/cache/yum/*
RUN ln -s /usr/local/cuda/lib64/stubs/libnvidia-ml.so /usr/local/cuda/lib64/stubs/libnvidia-ml.so.1 && \
    ln -s /usr/local/cuda/lib64/stubs/libcuda.so /usr/local/cuda/lib64/stubs/libcuda.so.1 && \
    mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp http://mvapich.cse.ohio-state.edu/download/mvapich/mv2/mvapich2-2.3.1.tar.gz && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/mvapich2-2.3.1.tar.gz -C /var/tmp -z && \
    cd /var/tmp/mvapich2-2.3.1 &&  CC=gcc CXX=g++ F77=gfortran FC=gfortran LD_LIBRARY_PATH='/usr/local/cuda/lib64/stubs:$LD_LIBRARY_PATH' ./configure --prefix=/usr/local/mvapich2 --disable-mcast --enable-cuda --with-cuda=/usr/local/cuda && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/mvapich2-2.3.1.tar.gz /var/tmp/mvapich2-2.3.1
ENV LD_LIBRARY_PATH=/usr/local/mvapich2/lib:$LD_LIBRARY_PATH \
    PATH=/usr/local/mvapich2/bin:$PATH \
    PROFILE_POSTLIB="-L/usr/local/cuda/lib64/stubs -lnvidia-ml -lcuda"

ENV LD_LIBRARY_PATH=/usr/local/cuda/lib64/stubs:$LD_LIBRARY_PATH

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
    cd /var/tmp/hdf5-1.10.5 &&  CC=mpicc CXX=mpic++ F77=mpifort F90=mpifort FC=mpifort ./configure --prefix=/usr/local/hdf5 --enable-parallel --enable-threadsafe --enable-unsupported --enable-cxx --enable-fortran && \
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
