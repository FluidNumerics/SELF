FROM gcr.io/self-fluids/self-base:nvidia-cuda11.2.1-rocm4.3 as builder
ARG GPU_TARGET=sm_72
ARG HIP_PLATFORM=nvidia
ARG PREC=double
ARG FFLAGS="-cpp -pg -g -O0 -C -Wall -fbounds-check -fbacktrace --coverage -ffpe-trap=invalid,zero,overflow"

# Install SELF
COPY . /build

ENV DEBIAN_FRONTEND=noninteractive   \
    LANGUAGE=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

RUN apt-get -yqq update \
 && apt-get -yqq install --no-install-recommends \
        build-essential \
 && rm -rf /var/lib/apt/lists/*

RUN . /etc/profile.d/z10_spack_environment.sh && \
    cd /build && \
    HIP_PLATFORM=${HIP_PLATFORM} \
    SELF_DIR=/build/ \
    SELF_PREFIX=/opt/self \
    PREC=${PREC} \
    GPU_TARGET=${GPU_TARGET} \
    SELF_FFLAGS=${FFLAGS} \
    make


# Bare OS image to run the installed executables
FROM nvidia/cuda:11.2.1-base

COPY --from=builder /build /build
COPY --from=builder /opt/rocm /opt/rocm
COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software
COPY --from=builder /opt/view /opt/view
COPY --from=builder /etc/profile.d/z10_spack_environment.sh /etc/profile.d/z10_spack_environment.sh

ENV DEBIAN_FRONTEND=noninteractive   \
    LANGUAGE=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

RUN apt-get -yqq update \
 && apt-get -yqq install --no-install-recommends \
        g++ \
        gcc \
        gfortran \
        libelf-dev \
        gcovr \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/self /opt/self

ENV SELF_PREFIX=/opt/self

ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l"]

