FROM gcr.io/self-fluids/self-base as builder
ARG GPU_TARGET=sm_70

RUN . /etc/profile.d/z10_spack_environment.sh && \
    cd /tmp && \
    BUILD=${BUILD_TYPE} \
    SELF_PREFIX=/opt/self \
    FC=gfortran \
    PREC=double \
    HIPFORT_GPU=${GPU_TARGET} \
    make

FROM ubuntu:focal

COPY --from=builder /opt /opt

LABEL maintainer="joe@fluidnumerics.com"
