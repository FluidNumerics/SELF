FROM nvidia/cuda:11.8.0-devel-ubuntu20.04 as bootstrap

ARG ROCM_VERSION=5.2.3
ARG SPACK_VERSION=v0.19.2

ENV SPACK_ROOT=/opt/spack \
    CURRENTLY_BUILDING_DOCKER_IMAGE=1 \
    container=docker

ENV DEBIAN_FRONTEND=noninteractive   \
    LANGUAGE=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

RUN apt-get clean \
 && (apt-get update -y || apt-get update -y) \
 && apt-get install -y --no-install-recommends wget \
 && apt-get -yqq install --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        file \
        g++ \
        gcc \
        gfortran \
        git \
        gnupg2 \
        iproute2 \
        locales \
        lua-posix \
        make \
        python3 \
        python3-pip \
        python3-setuptools \
        unzip \
        wget \
        lcov \
 && locale-gen en_US.UTF-8 \
 && rm -rf /var/lib/apt/lists/*

RUN wget -q -O - https://repo.radeon.com/rocm/rocm.gpg.key | apt-key add - && \
    echo "deb [arch=amd64] https://repo.radeon.com/rocm/apt/${ROCM_VERSION}/ ubuntu main" | tee /etc/apt/sources.list.d/rocm.list &&\
    apt-get -yqq update && \
    apt-get -yqq install rocm-dev

RUN mkdir $SPACK_ROOT && cd $SPACK_ROOT && \
    git clone https://github.com/spack/spack.git . && git checkout ${SPACK_VERSION}  && \
    mkdir -p $SPACK_ROOT/opt/spack

RUN ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
          /usr/local/bin/docker-shell \
 && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
          /usr/local/bin/interactive-shell \
 && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
          /usr/local/bin/spack-env

RUN mkdir -p /root/.spack \
 && cp $SPACK_ROOT/share/spack/docker/modules.yaml \
        /root/.spack/modules.yaml \
 && rm -rf /root/*.* /run/nologin $SPACK_ROOT/.git

# [WORKAROUND]
# https://superuser.com/questions/1241548/
#     xubuntu-16-04-ttyname-failed-inappropriate-ioctl-for-device#1253889
RUN [ -f ~/.profile ]                                               \
 && sed -i 's/mesg n/( tty -s \\&\\& mesg n || true )/g' ~/.profile \
 || true


WORKDIR /root
SHELL ["docker-shell"]

# Creates the package cache
RUN spack spec hdf5+mpi

ENTRYPOINT ["/bin/bash", "/opt/spack/share/spack/docker/entrypoint.bash"]
CMD ["interactive-shell"]

# Build stage with Spack pre-installed and ready to be used
FROM bootstrap as builder
ARG ROCM_VERSION=5.2.3

# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir /opt/spack-environment \
&&  export ROCM_VERSION=${ROCM_VERSION} \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - hipfort@${ROCM_VERSION}" \
&&   echo "  - hdf5@1.12.0+cxx+fortran+mpi" \
&&   echo "  - json-fortran@7.1.0" \
&&   echo "  - feq-parse@1.1.0" \
&&   echo "  - flap@master" \
&&   echo "  packages:" \
&&   echo "    cuda:" \
&&   echo "      buildable: false" \
&&   echo "      externals:" \
&&   echo "      - spec: cuda@11.8.0" \
&&   echo "        prefix: /usr/local/cuda" \
&&   echo "    hip:" \
&&   echo "      buildable: false" \
&&   echo "      externals:" \
&&   echo "      - spec: hip@${ROCM_VERSION}" \
&&   echo "        prefix: /opt/rocm/" \
&&   echo "  config:" \
&&   echo "    install_tree: /opt/software" \
&&   echo "  view: /opt/view") > /opt/spack-environment/spack.yaml

# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && \
    spack env activate . && \
    spack install --fail-fast && \
    spack gc -y

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh && \
    echo "LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/opt/view/lib:/opt/view/lib64" >> /etc/profile.d/z10_spack_environment.sh

FROM builder as sbuild
ARG GPU_TARGET=sm_72
ARG HIP_PLATFORM=nvidia
ARG PREC=double
ARG FFLAGS="-cpp -pg -g -O0 -C -Wall -fbounds-check -fbacktrace --coverage -ffpe-trap=invalid,zero,overflow"

# Install SELF
COPY . /build

RUN . /etc/profile.d/z10_spack_environment.sh && \
    cd /build && \
    HIP_PLATFORM=${HIP_PLATFORM} \
    SELF_DIR=/build/ \
    SELF_PREFIX=/opt/self \
    PREC=${PREC} \
    GPU_TARGET=${GPU_TARGET} \
    SELF_FFLAGS=${FFLAGS} \
    make && \
    cp -r /build/util /opt/self/

ENV SELF_PREFIX=/opt/self

## Bare OS image to run the installed executables
#FROM nvidia/cuda:11.8.0-devel
#
#COPY --from=builder /opt/rocm /opt/rocm
#COPY --from=builder /opt/spack-environment /opt/spack-environment
#COPY --from=builder /opt/software /opt/software
#COPY --from=builder /opt/view /opt/view
#COPY --from=builder /etc/profile.d/z10_spack_environment.sh /etc/profile.d/z10_spack_environment.sh
#COPY --from=sbuild /opt/self /opt/self
#COPY --from=sbuild /workspace /workspace

#ENV DEBIAN_FRONTEND=noninteractive   \
#    LANGUAGE=en_US.UTF-8 \
#    LANG=en_US.UTF-8 \
#    LC_ALL=en_US.UTF-8
#
#RUN rm /etc/apt/sources.list.d/cuda.list \
# && rm /etc/apt/sources.list.d/nvidia-ml.list \
# && apt-key del 7fa2af80 \
# && apt-get -y update \
# && apt-get install -y --no-install-recommends wget \
# && wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-keyring_1.0-1_all.deb \
# && dpkg -i cuda-keyring_1.0-1_all.deb \
# && apt-get -yqq update \
# && apt-get -yqq install --no-install-recommends \
#        g++ \
#        gcc \
#        gfortran \
#        gcovr \
#        build-essential \
#        libelf-dev \
# && rm -rf /var/lib/apt/lists/*


ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l", "-c"]

