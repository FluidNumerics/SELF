FROM nvidia/cuda:11.8.0-devel-ubuntu22.04 as cuda

FROM us-docker.pkg.dev/fluidnumerics-self/self-base/ubuntu2204-rocm-dev-5.7

COPY --from=cuda /usr/local/cuda /usr/local/cuda
ENV PATH=$PATH:/usr/local/cuda/bin

ARG GPU_TARGET=sm_72
ARG BUILD_TYPE=Debug

# Install SELF
COPY . /build

RUN . /etc/profile.d/z10_spack_environment.sh && \
    mkdir /build/build && \
    cd /build/build && \
    FC=gfortran \
    FFLAGS="-DDOUBLE_PRECISION" \
    cmake -DCMAKE_PREFIX_PATH=/opt/rocm \
          -DCMAKE_HIP_ARCHITECTURES=${GPU_TARGET} \
          -DCMAKE_INSTALL_PREFIX=/opt/self \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          ../ && \
    make VERBOSE=1 && \
    make install
    
# Make input and output directories
RUN mkdir /input && \
    mkdir /output

ENV SELF_PREFIX=/opt/self
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/self/lib:/opt/self/lib64
ENV PATH=$PATH:/opt/self/bin:/opt/self/test

# Set the entrypoint to the main self program with the input file set to /self/input.json
ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l", "-c", "/opt/self/bin/self -i /self/input.json"]

WORKDIR /self

