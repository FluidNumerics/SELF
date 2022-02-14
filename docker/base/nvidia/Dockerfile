FROM nvidia/cuda:11.2.1-devel as bootstrap

ENV SPACK_ROOT=/opt/spack \
    CURRENTLY_BUILDING_DOCKER_IMAGE=1 \
    container=docker

ENV DEBIAN_FRONTEND=noninteractive   \
    LANGUAGE=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

RUN apt-get -yqq update \
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
 && locale-gen en_US.UTF-8 \
 && pip3 install boto3 \
 && rm -rf /var/lib/apt/lists/*

RUN wget -q -O - https://repo.radeon.com/rocm/rocm.gpg.key | apt-key add - && \
    echo 'deb [arch=amd64] https://repo.radeon.com/rocm/apt/4.3/ xenial main' | tee /etc/apt/sources.list.d/rocm.list &&\
    apt-get -yqq update && \
    apt-get -yqq install rocm-dev rocm-libs

RUN mkdir $SPACK_ROOT && cd $SPACK_ROOT && \
    git clone https://github.com/spack/spack.git . && git checkout develop  && \
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
ARG GPU_TARGET=sm_72
ARG HIP_PLATFORM=nvidia
ARG PREC=double


# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir /opt/spack-environment \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - hipfort@4.3.0" \
&&   echo "  - hdf5@1.12.0+cxx+fortran+mpi" \
&&   echo "  - json-fortran@7.1.0" \
&&   echo "  - feq-parse@1.1.0" \
&&   echo "  - flap@master" \
&&   echo "  packages:" \
&&   echo "    cuda:" \
&&   echo "      buildable: false" \
&&   echo "      externals:" \
&&   echo "      - spec: cuda@11.2.1" \
&&   echo "        prefix: /usr/local/cuda" \
&&   echo "    hip:" \
&&   echo "      buildable: false" \
&&   echo "      externals:" \
&&   echo "      - spec: hip@4.3.0" \
&&   echo "        prefix: /opt/rocm/" \
&&   echo "  concretization: together" \
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
    spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh

# Bare OS image to run the installed executables
#FROM nvidia/cuda:11.2.1-base
#
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
# && rm -rf /var/lib/apt/lists/*


ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l"]

