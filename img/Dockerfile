FROM gcr.io/self-fluids/self-base:latest as builder
ARG GPU_TARGET=sm_72
ARG HIP_PLATFORM=nvidia

COPY . /tmp

RUN . /etc/profile.d/z10_spack_environment.sh && \
    cd /tmp && \
    HIP_PLATFORM=${HIP_PLATFORM} \
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

LABEL Maintainer="joe@fluidnumerics.com"
LABEL OS Ubuntu
LABEL Compiler gcc@9.2.0
LABEL MPI openmpi@4.0.5
LABEL Dependencies "hdf5@1.12.0+cxx+fortran+mpi, feq-parse@v1.1.0, flap@master"
LABEL Platform ${HIP_PLATFORM}
LABEL GPU ${GPU_TARGET}
