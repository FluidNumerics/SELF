FROM debian:bullseye AS devel
ARG BUILD_TYPE=release

COPY . /tmp

RUN mkdir -p /tmp/extern

# Install gcc, gfortran, and cmake
RUN apt update -y && \
    apt install -y gcc gfortran cmake git

# FEQParse
RUN git clone https://github.com/FluidNumerics/feq-parse.git /tmp/extern/feq-parse && \
    mkdir -p /tmp/extern/feq-parse/build && \
    cd /tmp/extern/feq-parse/build && \
    cmake -DCMAKE_INSTALL_PREFIX="/opt/feqparse" /tmp/extern/feq-parse && \
    make && make install

# FLAP
RUN git clone --recurse-submodules https://github.com/szaghi/FLAP.git /tmp/extern/FLAP && \
    mkdir -p /tmp/extern/FLAP/build && \
    cd /tmp/extern/FLAP/build && \
    FFLAGS=-cpp cmake -DCMAKE_INSTALL_PREFIX="/opt/FLAP" /tmp/extern/FLAP && \
    make && make install

RUN cd /tmp && \
    BUILD=${BUILD_TYPE} \
    SELF_PREFIX=/opt/self \
    FC=gfortran \
    make

FROM debian:bullseye

COPY --from=devel /opt /opt

# Install gcc, gfortran
RUN apt update -y && \
    apt install -y gcc gfortran

ENV LD_LIBRARY_PATH=/opt/self/lib:/opt/feqparse/lib:/opt/FLAP/lib:$LD_LIBRARY_PATH \
    PATH=/opt/self/bin:$PATH
