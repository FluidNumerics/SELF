FROM nvidia/cuda:10.1-devel-centos7 AS devel

#COPY . /build
#
## Install GCC-9
#RUN yum install -y centos-release-scl && \
#    yum install -y \
#        devtoolset-9-gcc \
#        devtoolset-9-gcc-c++ \
#        devtoolset-9-gcc-gfortran && \
#    rm -rf /var/cache/yum/*
#ENV PATH=/opt/rh/devtoolset-9/root/usr/bin:$PATH
#
#RUN  source /opt/rh/devtoolset-9/enable
#
## Install rocm-3.8.0
#RUN yum install -y epel-release && \
#    yum install -y dkms kernel-headers-`uname -r` kernel-devel-`uname -r`
#
#RUN echo -e "[ROCm]\nname=ROCm\nbaseurl=http://repo.radeon.com/rocm/yum/3.8\nenabled=1\ngpgcheck=1\ngpgkey=http://repo.radeon.com/rocm/rocm.gpg.key" > /etc/yum.repos.d/rocm.repo && yum -y update
#
#RUN yum install -y rocm-dkms rocm-dev
#RUN useradd -u 4000 self && usermod -a -G video self &&  echo 'ADD_EXTRA_GROUPS=1' | tee -a /etc/adduser.conf && echo 'EXTRA_GROUPS=video' | tee -a /etc/adduser.conf
#
## OpenMPI
#RUN yum group install -y 'Development Tools' && \
#    yum install -y perl-core \
#        libtemplate-perl \
#        zlib-devel \
#        bzip2 \
#        file \
#        hwloc \
#        make \
#        numactl-devel \
#        openssh-clients \
#        perl \
#        perl-Data-Dumper \
#        autoconf \
#        automake \
#        libtool \
#        libevent \
#        tar \
#        wget && \
#    rm -rf /var/cache/yum/*

RUN wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.5.tar.gz --directory-prefix=/tmp && \
    cd /tmp && tar -xvzf openmpi-4.0.5.tar.gz && \
    cd /tmp/openmpi-4.0.5 && CC=gcc CXX=g++ F77=gfortran F90=gfortran FC=gfortran ./configure --prefix=/usr/local/self --disable-getpwuid --enable-orterun-prefix-by-default --with-cuda --without-verbs && \
    make -j$(nproc) && \
    make -j$(nproc) install 
    
## HDF5 (from SELF submodules)
#RUN yum install -y \
#        bzip2 \
#        file \
#        make \
#        wget \
#        zlib-devel && \
#    rm -rf /var/cache/yum/*
#RUN cd /build/dependencies/hdf5 &&  CC=mpicc CXX=mpic++ F77=mpifort F90=mpifort FC=mpifort ./configure --prefix=/usr/local/self --enable-parallel --enable-threadsafe --enable-unsupported --enable-cxx --enable-fortran && \
#    make -j$(nproc) && \
#    make -j$(nproc) install
#
#ENV HDF5_DIR=/usr/local/self \
#    LD_LIBRARY_PATH=/usr/local/self/lib:/opt/rocm/lib:$LD_LIBRARY_PATH \
#    PATH=/usr/local/self/bin:/opt/rocm/bin:$PATH
#
##HOPR
#RUN mkdir /build/dependencies/hopr/build && cd /build/dependencies/hopr/build && cmake -DCMAKE_INSTALL_PREFIX=/usr/local/self ../ && make && make install
#
## FEQParse
#RUN mkdir /build/dependencies/feq-parse/build && cd /build/dependencies/feq-parse/build && cmake -DCMAKE_INSTALL_PREFIX=/usr/local/self ../ && make && make install
#
## JSON-Fortran
#RUN mkdir /build/dependencies/json-fortran/build && cd /build/dependencies/json-fortran/build && cmake -DSKIP_DOC_GEN=True -DCMAKE_INSTALL_PREFIX=/usr/local/self ../ && make && make install
#
## HIPFort
#RUN mkdir /build/dependencies/hipfort/build && cd /build/dependencies/hipfort/build && cmake -DCMAKE_INSTALL_PREFIX=/opt/rocm/hipfort ../ && make && make install
#
