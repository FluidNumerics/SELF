FROM gcr.io/self-fluids/self-base:cuda-11.2.1 as builder
ARG GPU_TARGET=sm_72
ARG HIP_PLATFORM=nvidia
ARG PREC=double

# Install SELF
COPY . /workspace

RUN . /etc/profile.d/z10_spack_environment.sh && \
    cd /workspace && \
    HIP_PLATFORM=${HIP_PLATFORM} \
    SELF_DIR=/workspace/ \
    SELF_PREFIX=/opt/self \
    PREC=${PREC} \
    GPU_TARGET=${GPU_TARGET} \
    make

# Bare OS image to run the installed executables
# Bare OS image to run the installed executables
FROM nvidia/cuda:11.2.1-base

COPY --from=builder /workspace /workspace
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
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/self /opt/self

ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l"]

