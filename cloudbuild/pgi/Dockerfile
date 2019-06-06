FROM launcher.gcr.io/google/centos7


# Install requirements (gcc & g++)
RUN yum -y update && \
    yum install -y gcc \
                   gcc-c++ \
                   gcc-gfortran \
                   autotools \
                   autoconf \
                   automake \
                   make \
                   perl \
                   wget

# Install PGI
ENV PGI_VERSION 18.10
ENV PGI_INSTALL_DIR /usr/local/pgi
ENV PGI_HOME    ${PGI_INSTALL_DIR}/linux86-64/${PGI_VERSION}
ENV PGI_BIN_DIR ${PGI_HOME}/bin
ENV PGI_LIB_DIR ${PGI_HOME}/lib
ENV PGI_MAN_DIR ${PGI_HOME}/man

ADD pgilinux-2018-1810-x86-64.tar.gz /tmp

RUN export PGI_SILENT=true && \
    export PGI_ACCEPT_EULA=accept && \
    export PGI_INSTALL_NVIDIA=true && \
    export PGI_INSTALL_MANAGED=true && \
    export PGI_INSTALL_AMD=false && \
    export PGI_INSTALL_JAVA=false && \
    export PGI_INSTALL_MPI=true && \
    export PGI_MPI_GPU_SUPPORT=true && \
    /tmp/install && \
    rm -rf /tmp/*

ADD hdf5-1.10.5.tar.gz /tmp

# HDF5
RUN ls /tmp/

RUN CPP=cpp CFLAGS="-fPIC -m64 -tp=px" CXXFLAGS="-fPIC -m64 -tp=px" FCFLAGS="-fPIC -m64 -tp=px" \
    CC=${PGI_INSTALL_DIR}/linux86-64/2018/bin/pgcc \
    CXX=${PGI_INSTALL_DIR}/linux86-64/2018/bin/pgc++ \
    FC=${PGI_INSTALL_DIR}/linux86-64/2018/bin/pgfortran \
    /tmp/hdf5-1.10.5/configure --enable-threadsafe --enable-cxx --enable-fortran --enable-unsupported --prefix=/usr/local/hdf5
RUN make
RUN make install

# METIS
ADD metis-5.1.0.tar.gz /tmp
RUN make config prefix=/usr/local/metis -f /tmp/metis-5.1.0/Makefile 
RUN make install -f /tmp/metis-5.1.0/Makefile 


RUN echo "${PGI_LIB_DIR}" >> /etc/ld.so.conf.d/pgi.conf

ENV PATH            ${PGI_BIN_DIR}:${PATH}
ENV LD_LIBRARY_PATH ${PGI_LIB_DIR}:${LD_LIBRARY_PATH}
ENV MANPATH         ${PGI_MAN_DIR}:${MANPATH}

