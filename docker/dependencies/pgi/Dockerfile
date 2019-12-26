# 
# HPC Base image
# 
# Contents:
#   CUDA version 10.0
#   HDF5 version 1.10.5
#   PGI compilers version 19.10
# 

FROM nvidia/cuda:10.1-devel-centos7 AS devel


# PGI compiler version 19.10
RUN yum install -y \
        automake \
        autoconf \
        gcc \
        gcc-c++ \
        numactl-libs \
        perl \
        wget && \
    rm -rf /var/cache/yum/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -O /var/tmp/pgi-community-linux-x64-latest.tar.gz --referer https://www.pgroup.com/products/community.htm?utm_source=hpccm\&utm_medium=wgt\&utm_campaign=CE\&nvid=nv-int-14-39155 -P /var/tmp https://www.pgroup.com/support/downloader.php?file=pgi-community-linux-x64 && \
    mkdir -p /var/tmp/pgi && tar -x -f /var/tmp/pgi-community-linux-x64-latest.tar.gz -C /var/tmp/pgi -z && \
    cd /var/tmp/pgi && PGI_ACCEPT_EULA=accept PGI_INSTALL_DIR=/opt/pgi PGI_INSTALL_MPI=false PGI_INSTALL_NVIDIA=true PGI_MPI_GPU_SUPPORT=false PGI_SILENT=true ./install && \
    echo "variable LIBRARY_PATH is environment(LIBRARY_PATH);" >> /opt/pgi/linux86-64/19.10/bin/siterc && \
    echo "variable library_path is default(\$if(\$LIBRARY_PATH,\$foreach(ll,\$replace(\$LIBRARY_PATH,":",), -L\$ll)));" >> /opt/pgi/linux86-64/19.10/bin/siterc && \
    echo "append LDLIBARGS=\$library_path;" >> /opt/pgi/linux86-64/19.10/bin/siterc && \
    ln -sf /usr/lib64/libnuma.so.1 /opt/pgi/linux86-64/19.10/lib/libnuma.so && \
    ln -sf /usr/lib64/libnuma.so.1 /opt/pgi/linux86-64/19.10/lib/libnuma.so.1 && \
    rm -rf /var/tmp/pgi-community-linux-x64-latest.tar.gz /var/tmp/pgi
ENV LD_LIBRARY_PATH=/opt/pgi/linux86-64/19.10/lib:$LD_LIBRARY_PATH \
    PATH=/opt/pgi/linux86-64/19.10/bin:$PATH

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
    cd /var/tmp/hdf5-1.10.5 && \
    CPP=cpp CFLAGS="-fPIC -m64 -tp=px" CXXFLAGS="-fPIC -m64 -tp=px" FCFLAGS="-fPIC -m64 -tp=px" CC=pgcc CXX=pgc++ F77=pgfortran F90=pgfortran FC=pgfortran ./configure --prefix=/usr/local/hdf5 --enable-cxx --enable-fortran || cat config.log && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/hdf5-1.10.5.tar.bz2 /var/tmp/hdf5-1.10.5
ENV HDF5_DIR=/usr/local/hdf5 \
    LD_LIBRARY_PATH=/usr/local/hdf5/lib:$LD_LIBRARY_PATH \
    PATH=/usr/local/hdf5/bin:$PATH

RUN yum install -y \
        cmake
RUN mkdir -p /var/tmp && \
    wget -q -nc --no-check-certificate -P /var/tmp http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz && \
    tar -xzf /var/tmp/metis-5.1.0.tar.gz -C /var/tmp && \
    cd /var/tmp/metis-5.1.0 && \
    make config prefix=/usr/local/metis && \
    make install

ENV LIB_METIS=/usr/local/metis/lib/libmetis.a
