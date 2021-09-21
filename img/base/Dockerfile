FROM ubuntu:focal as builder

RUN export DEBIAN_FRONTEND=noninteractive &&\
    apt-get update -y && \
    apt-get install -y git software-properties-common build-essential gfortran wget cmake

# Install rocm + hipfort
RUN export DEBIAN_FRONTEND=noninteractive \
&&  wget -q -O - https://repo.radeon.com/rocm/rocm.gpg.key | apt-key add - \
&&  echo 'deb [arch=amd64] https://repo.radeon.com/rocm/apt/4.3/ xenial main' | tee /etc/apt/sources.list.d/rocm.list \
&& apt-get update -y \
&& apt-get install -y rocm-dev rocm-libs \
&& git clone https://github.com/ROCmSoftwarePlatform/hipfort.git /tmp/hipfort \
&& mkdir /tmp/hipfort/build \
&& cd /tmp/hipfort/build \
&& FC=$(which gfortran) cmake -DCMAKE_INSTALL_PREFIX=/opt/rocm /tmp/hipfort \
&& make -j install \
&& echo "PATH=\${PATH}:/opt/rocm/bin" >> /etc/profile.d/z10_spack_environment.sh
&& echo "PATH=\${PATH}:/opt/rocm/hipfort/bin" >> /etc/profile.d/z10_spack_environment.sh

# Install cuda-compiler
RUN export DEBIAN_FRONTEND=noninteractive \
&& wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin \
&& mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600 \
&& wget https://developer.download.nvidia.com/compute/cuda/11.4.2/local_installers/cuda-repo-ubuntu2004-11-4-local_11.4.2-470.57.02-1_amd64.deb \
&& dpkg -i cuda-repo-ubuntu2004-11-4-local_11.4.2-470.57.02-1_amd64.deb \
&& apt-key add /var/cuda-repo-ubuntu2004-11-4-local/7fa2af80.pub \
&& apt-get update \
&& apt-get -y install cuda

RUN git clone https://github.com/FluidNumerics/spack.git /opt/spack && \
    echo "#!/bin/bash" > /etc/profile.d/spack.sh && \
    echo "SPACK_ROOT=/opt/spack" >> /etc/profile.d/spack.sh && \
    echo ". \${SPACK_ROOT}/share/spack/setup-env.sh" >> /etc/profile.d/spack.sh && \
    . /opt/spack/share/spack/setup-env.sh && \
    spack compiler find --scope=defaults

# Install HDF5
RUN mkdir /opt/spack-environment \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - cmake@3.18.4" \
&&   echo "  - openmpi@4.0.2" \
&&   echo "  - hdf5@1.12.0+cxx+fortran+mpi" \
&&   echo "  - feq-parse@v1.1.0" \
&&   echo "  - flap@master" \
&&   echo "  concretization: together" \
&&   echo "  config:" \
&&   echo "    install_tree: /opt/software" \
&&   echo "  view: /opt/view") > /opt/spack-environment/spack.yaml

# Install the software, remove unnecessary deps
RUN . /opt/spack/share/spack/setup-env.sh && \ 
    cd /opt/spack-environment && spack env activate . && spack install --fail-fast && spack gc -y

# Modifications to the environment that are necessary to run
RUN . /opt/spack/share/spack/setup-env.sh && \ 
    cd /opt/spack-environment && \
    spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh

LABEL maintainer="joe@fluidnumerics.com"

RUN mkdir /workspace
