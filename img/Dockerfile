FROM gcr.io/self-fluids/self-base:latest as builder
ARG GPU_TARGET=sm_72

COPY . /tmp

RUN . /etc/profile.d/z10_spack_environment.sh && \
    cd /tmp && \
    BUILD=${BUILD_TYPE} \
    SELF_PREFIX=/opt/self \
    HIPFC=/opt/rocm/bin/hipfc \
    HIPFORT_COMPILER=h5pfc \
    PREC=double \
    HIPFORT_GPU=${GPU_TARGET} \
    make

FROM ubuntu:focal

COPY --from=builder /opt /opt

LABEL maintainer="joe@fluidnumerics.com"
