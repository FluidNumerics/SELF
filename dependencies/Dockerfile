FROM nvidia/cuda:10.1-devel-centos7 AS devel

## Install GCC-9
RUN yum install -y centos-release-scl && \
    yum install -y \
        wget \
        openssl-devel \
        devtoolset-9-gcc \
        devtoolset-9-gcc-c++ \
        devtoolset-9-gcc-gfortran && \
    rm -rf /var/cache/yum/*
ENV PATH=/opt/rh/devtoolset-9/root/usr/bin:$PATH

RUN source /opt/rh/devtoolset-9/enable

# Install rocm-3.8.0
RUN yum install -y epel-release && \
    yum install -y dkms kernel-headers-`uname -r` kernel-devel-`uname -r`

RUN echo -e "[ROCm]\nname=ROCm\nbaseurl=http://repo.radeon.com/rocm/yum/3.8\nenabled=1\ngpgcheck=1\ngpgkey=http://repo.radeon.com/rocm/rocm.gpg.key" > /etc/yum.repos.d/rocm.repo && yum -y update

RUN yum install -y rocm-dkms rocm-dev
RUN echo 'ADD_EXTRA_GROUPS=1' | tee -a /etc/adduser.conf && echo 'EXTRA_GROUPS=video' | tee -a /etc/adduser.conf

## Install CMake 3.18.4
RUN wget https://github.com/Kitware/CMake/releases/download/v3.18.4/cmake-3.18.4.tar.gz --directory-prefix=/tmp && \
    cd /tmp &&\
    tar -xvzf cmake-3.18.4.tar.gz &&\
    cd /tmp/cmake-3.18.4 &&\
    ./bootstrap --prefix=/usr/local &&\
    make &&\
    make install

## OpenMPI
RUN yum group install -y 'Development Tools' && \
    yum install -y perl-core \
        libtemplate-perl \
        zlib-devel \
        bzip2 \
        file \
        hwloc \
        make \
        numactl-devel \
        openssh-clients \
        perl \
        perl-Data-Dumper \
        autoconf \
        automake \
        libtool \
        libevent \
        tar && \
    rm -rf /var/cache/yum/*

RUN wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.5.tar.gz --directory-prefix=/tmp && \
    cd /tmp && tar -xvzf openmpi-4.0.5.tar.gz && \
    cd /tmp/openmpi-4.0.5 && CC=gcc CXX=g++ F77=gfortran F90=gfortran FC=gfortran ./configure --prefix=/usr/local/self --disable-getpwuid --enable-orterun-prefix-by-default --with-cuda --without-verbs && \
    make -j$(nproc) && \
    make -j$(nproc) install 
    
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
    cd /var/tmp/hdf5-1.10.5 && \
    source /opt/rh/devtoolset-9/enable && \
    CC=/usr/local/self/bin/mpicc CXX=/usr/local/self/bin/mpic++ F77=/usr/local/self/bin/mpifort F90=/usr/local/self/bin/mpifort FC=/usr/local/self/bin/mpifort ./configure --prefix=/usr/local/self --enable-parallel --enable-threadsafe --enable-unsupported --enable-cxx --enable-fortran && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/hdf5-1.10.5.tar.bz2 /var/tmp/hdf5-1.10.5

# BLAS and Lapack
RUN wget https://github.com/Reference-LAPACK/lapack/archive/v3.9.0.tar.gz --directory-prefix=/tmp/ &&\
    cd /tmp && tar -xvzf v3.9.0.tar.gz && mkdir -p /tmp/lapack-3.9.0/build && \
    cd /tmp/lapack-3.9.0/build && \
    source /opt/rh/devtoolset-9/enable && \
    /usr/local/bin/cmake -DCMAKE_INSTALL_PREFIX=/usr/local/self /tmp/lapack-3.9.0 &&  \
    make && make install &&\
    rm -rf /tmp/lapack-3.9.0

ENV HDF5_ROOT=/usr/local/self \
    LD_LIBRARY_PATH=/usr/local/self/lib:/opt/rocm/lib:$LD_LIBRARY_PATH \
    PATH=/usr/local/self/bin:/opt/rocm/bin:$PATH

#HOPR
RUN git clone https://github.com/FluidNumerics/hopr.git /tmp/hopr && \
    mkdir /tmp/hopr/build && \
    cd /tmp/hopr/build && \
    source /opt/rh/devtoolset-9/enable && \
    HDF5_ROOT=/usr/local/self /usr/local/bin/cmake -DHOPR_BUILD_HDF5=OFF -DCMAKE_INSTALL_PREFIX=/usr/local/self /tmp/hopr &&\
    make && make install &&\
    rm -rf /tmp/hopr

# FEQParse
RUN git clone https://github.com/FluidNumerics/feq-parse.git /tmp/feq-parse &&\
    mkdir -p /tmp/feq-parse/build && \
    cd /tmp/feq-parse/build && \
    source /opt/rh/devtoolset-9/enable && \
    FC=gfortran /usr/local/bin/cmake -DCMAKE_INSTALL_PREFIX=/usr/local/self /tmp/feq-parse && \
    make && make install && \
    rm -rf /tmp/feq-parse

# JSON-Fortran
RUN git clone https://github.com/jacobwilliams/json-fortran.git /tmp/json-fortran &&\
    mkdir -p /tmp/json-fortran/build && \
    cd /tmp/json-fortran/build && \
    source /opt/rh/devtoolset-9/enable && \
    FC=gfortran /usr/local/bin/cmake -DSKIP_DOC_GEN=True -DCMAKE_INSTALL_PREFIX=/usr/local/self /tmp/json-fortran && \
    make && make install &&\
    rm -rf /tmp/json-fortran

# HIPFort
RUN git clone https://github.com/ROCmSoftwarePlatform/hipfort.git /tmp/hipfort &&\
    mkdir -p /tmp/hipfort/build && \
    cd /tmp/hipfort/build && \
    source /opt/rh/devtoolset-9/enable && \
    FC=gfortran /usr/local/bin/cmake -DCMAKE_INSTALL_PREFIX=/opt/rocm/hipfort /tmp/hipfort &&\
    make && make install &&\
    rm -rf /tmp/hipfort


ENV HDF5_ROOT=/usr/local/self \
    LD_LIBRARY_PATH=/opt/rocm/hipfort/lib:$LD_LIBRARY_PATH \
    PATH=/opt/rocm/hipfort/bin:$PATH
