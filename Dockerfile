FROM debian:bullseye AS devel
ARG BUILD_TYPE=debug

COPY . /tmp

RUN mkdir -p /tmp/extern

# Install gcc, gfortran, and cmake
RUN apt update -y && \
    apt install -y gcc gfortran cmake git

# FEQParse
RUN git clone https://github.com/FluidNumerics/feq-parse.git /tmp/extern/feq-parse && \
    mkdir -p /tmp/extern/feq-parse/build && \
    cd /tmp/extern/feq-parse/build && \
    cmake -DCMAKE_INSTALL_PREFIX="/opt/self" /tmp/extern/feq-parse && \
    make && make install

# JSON-Fortran
RUN git clone https://github.com/jacobwilliams/json-fortran.git /tmp/extern/json-fortran && \
    mkdir -p /tmp/extern/json-fortran/build && \
    cd /tmp/extern/json-fortran/build && \
    cmake -DSKIP_DOC_GEN=True -DCMAKE_INSTALL_PREFIX="/opt/self" /tmp/extern/json-fortran && \
    make && make install

# FLAP
RUN git clone --recurse-submodules https://github.com/szaghi/FLAP.git /tmp/extern/FLAP && \
    mkdir -p /tmp/extern/FLAP/build && \
    cd /tmp/extern/FLAP/build && \
    FFLAGS=-cpp cmake -DCMAKE_INSTALL_PREFIX="/opt/self" /tmp/extern/FLAP && \
    make && make install

# Focal
RUN git clone https://github.com/LKedward/focal.git --depth 1 --branch v1.0.1 /tmp/extern/focal && \
    cd /tmp/extern/focal/ && \
    make && \
    mv obj/*.o /opt/self/lib/ && mv mod/*.mod /opt/self/include/

ENV BUILD_TYPE=${BUILD_TYPE}


RUN cd /tmp && \
    make install -f /tmp/build/self.make


FROM debian:bullseye

COPY --from=devel /opt/self /opt/self

# Install gcc, gfortran, cmake, git
RUN apt update -y && \
    apt install -y gcc gfortran cmake git

ENV LD_LIBRARY_PATH=/opt/self/lib:$LD_LIBRARY_PATH \
    PATH=/opt/self/bin:$PATH
