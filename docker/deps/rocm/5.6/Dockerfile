FROM ubuntu:22.04 as rocm

ARG ROCM_VERSION=5.7
ARG PACKAGE_SET="rocm-dev"

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
        locales \
        wget \
        gnupg2 \ 
        gcc \
        gfortran \
 && locale-gen en_US.UTF-8 \
 && rm -rf /var/lib/apt/lists/*


# Install rocm
RUN mkdir --parents --mode=0755 /etc/apt/keyrings &&\
    wget https://repo.radeon.com/rocm/rocm.gpg.key -O - | \
    gpg --dearmor | tee /etc/apt/keyrings/rocm.gpg > /dev/null &&\
    echo "deb [arch=amd64 signed-by=/etc/apt/keyrings/rocm.gpg] https://repo.radeon.com/rocm/apt/${ROCM_VERSION}/ jammy main" | \
    tee /etc/apt/sources.list.d/rocm.list &&\
    echo "Package: *\nPin: release o=repo.radeon.com\nPin-Priority: 600" | \
    tee /etc/apt/preferences.d/rocm-pin-600 && \
    apt-get -y update && \
    apt-get -y install ${PACKAGE_SET}

# Build stage with Spack pre-installed and ready to be used
FROM spack/ubuntu-jammy:latest as builder
COPY --from=rocm /opt/rocm /opt/rocm
COPY --from=rocm /opt/rocm /opt/rocm-5.6.0

# 10/8/2023 - Check out fluidnumerics-joe/spack for updated feq-parse version
# When https://github.com/spack/spack/pull/40230 is merged, we won't need this
# anymore
RUN cd /opt/spack &&\
    git pull origin develop

# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir /opt/spack-environment \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - cmake@3.27.4" \
&&   echo "  - hdf5@1.12.0~cxx+fortran+mpi" \
&&   echo "  - json-fortran@7.1.0" \
&&   echo "  - feq-parse@2.0.3" \
&&   echo "  config:" \
&&   echo "    install_tree: /opt/software" \
&&   echo "    view: /opt/view" \
&&   echo "  concretizer:" \
&&   echo "    unify: true" \
&&   echo "  view: /opt/view") > /opt/spack-environment/spack.yaml

# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && spack env activate . && spack install --fail-fast && spack gc -y

# Strip all the binaries
RUN find -L /opt/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-archive\|x-sharedlib' | \
    awk -F: '{print $1}' | xargs strip -s

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh

# Bare OS image to run the installed executables
FROM rocm

COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software
COPY --from=builder /opt/._view /opt/._view
COPY --from=builder /opt/view /opt/view
COPY --from=builder /etc/profile.d/z10_spack_environment.sh /etc/profile.d/z10_spack_environment.sh

ENV PATH=$PATH:/opt/rocm/bin

RUN apt-get clean \
 && (apt-get update -y || apt-get update -y) \
 && apt-get install -y --no-install-recommends wget \
 && apt-get -yqq install --no-install-recommends \
        gcc \
        gfortran \ 
 && locale-gen en_US.UTF-8 \
 && rm -rf /var/lib/apt/lists/*

ENV CMAKE_MODULE_PATH=$CMAKE_MODULE_PATH:/opt/rocm/hip

ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l", "-c", "$*", "--" ]
CMD [ "/bin/bash" ]

