FROM launcher.gcr.io/google/centos7 as intermediate

ADD pgilinux-2018-1810-x86-64.tar.gz /tmp
ADD hdf5-1.10.5.tar.gz /tmp
ADD metis-5.1.0.tar.gz /tmp
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
                   cmake \
                   wget

# Install PGI
ENV PGI_VERSION 18.10
ENV PGI_INSTALL_DIR /usr/local/pgi
RUN export PGI_SILENT=true && \
    export PGI_ACCEPT_EULA=accept && \
    export PGI_INSTALL_NVIDIA=true && \
    export PGI_INSTALL_MANAGED=true && \
    export PGI_INSTALL_AMD=false && \
    export PGI_INSTALL_JAVA=false && \
    export PGI_INSTALL_MPI=true && \
    export PGI_MPI_GPU_SUPPORT=true && \
    /tmp/install && \
    echo "${PGI_INSTALL_DIR}/linux86-64/${PGI_VERSION}/bin" >> /etc/ld.so.conf.d/pgi.conf


# HDF5
RUN CPP=cpp CFLAGS="-fPIC -m64 -tp=px" CXXFLAGS="-fPIC -m64 -tp=px" FCFLAGS="-fPIC -m64 -tp=px" \
    CC=${PGI_INSTALL_DIR}/linux86-64/2018/bin/pgcc \
    CXX=${PGI_INSTALL_DIR}/linux86-64/2018/bin/pgc++ \
    FC=${PGI_INSTALL_DIR}/linux86-64/2018/bin/pgfortran \
    /tmp/hdf5-1.10.5/configure --enable-threadsafe --enable-cxx --enable-fortran --enable-unsupported --prefix=/usr/local/hdf5 &&\
    make && make install

# METIS
RUN cd /tmp/metis-5.1.0/ && make config prefix=/usr/local/metis && make install -f /tmp/metis-5.1.0/Makefile 


FROM launcher.gcr.io/google/centos7
COPY --from=intermediate /usr/local /usr/local

RUN yum -y update && \
    yum install -y gcc \
                   gcc-c++ \
                   gcc-gfortran \
                   autotools \
                   autoconf \
                   automake \
                   make \
                   perl \
                   cmake \
                   wget

# Set Environment variables
ENV PGI_VERSION 18.10
ENV PATH            ${PGI_INSTALL_DIR}/linux86-64/${PGI_VERSION}/bin:/usr/local/hdf5/bin:${PATH}
ENV LD_LIBRARY_PATH ${PGI_INSTALL_DIR}/linux86-64/${PGI_VERSION}/lib:/usr/local/hdf5/lib${LD_LIBRARY_PATH}
ENV MANPATH         ${PGI_INSTALL_DIR}/linux86-64/${PGI_VERSION}/man:${MANPATH}


