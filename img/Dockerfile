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

FROM gcr.io/self-fluids/self-base:latest

COPY --from=builder /opt /opt
COPY --from=builder /etc/profile.d /etc/profile.d

ENV PATH=$PATH:/opt/view/bin
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/view/lib

LABEL maintainer="joe@fluidnumerics.com"
