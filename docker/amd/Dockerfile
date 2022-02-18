FROM gcr.io/self-fluids/self-base:amd-rocm4.3 as builder
ARG GPU_TARGET=gfx900
ARG HIP_PLATFORM=amd
ARG PREC=double

# Install SELF
COPY . /build

RUN . /etc/profile.d/z10_spack_environment.sh && \
    cd /build && \
    HIP_PLATFORM=${HIP_PLATFORM} \
    SELF_DIR=/build/ \
    SELF_PREFIX=/opt/self \
    PREC=${PREC} \
    GPU_TARGET=${GPU_TARGET} \
    make

# Bare OS image to run the installed executables
#FROM nvidia/cuda:11.2.1-base
#
#COPY --from=builder /build /build
#COPY --from=builder /opt/rocm /opt/rocm
#COPY --from=builder /opt/spack-environment /opt/spack-environment
#COPY --from=builder /opt/software /opt/software
#COPY --from=builder /opt/view /opt/view
#COPY --from=builder /etc/profile.d/z10_spack_environment.sh /etc/profile.d/z10_spack_environment.sh
#
#ENV DEBIAN_FRONTEND=noninteractive   \
#    LANGUAGE=en_US.UTF-8 \
#    LANG=en_US.UTF-8 \
#    LC_ALL=en_US.UTF-8
#
#RUN apt-get -yqq update \
# && apt-get -yqq install --no-install-recommends \
#        g++ \
#        gcc \
#        gfortran \
#        libelf-dev \
#        gcovr \
#
# && rm -rf /var/lib/apt/lists/*
#
#COPY --from=builder /opt/self /opt/self

ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l"]

