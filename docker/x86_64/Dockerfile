# Build stage with Spack pre-installed and ready to be used
FROM spack/rockylinux8:develop AS builder


# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)

RUN mkdir -p /opt/spack-environment && \
set -o noclobber \
&&  (echo spack: \
&&   echo '  specs:' \
&&   echo '  - cmake@3.29.6' \
&&   echo '  - feq-parse@2.2.2' \
&&   echo '  - openmpi@5.0.2' \
&&   echo '  - hdf5@1.14.3 +fortran +mpi' \
&&   echo '  concretizer:' \
&&   echo '    unify: when_possible' \
&&   echo '  config:' \
&&   echo '    install_tree: /opt/software' \
&&   echo '  view: /opt/views/view') > /opt/spack-environment/spack.yaml

# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && spack info openmpi && spack env activate . && spack install --fail-fast && spack gc -y

# Strip all the binaries
RUN find -L /opt/views/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-archive\|x-sharedlib' | \
    awk -F: '{print $1}' | xargs strip

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . > activate.sh

RUN spack env activate -d /opt/spack-environment && \
    mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=release  \
          -DSELF_ENABLE_MULTITHREADING=ON \
          -DCMAKE_INSTALL_PREFIX=/opt/views/view \
          ../ && \
    make VERBOSE=1 && \
    make install 

# Bare OS image to run the installed executables
FROM docker.io/rockylinux:8

COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software

# paths.view is a symlink, so copy the parent to avoid dereferencing and duplicating it
COPY --from=builder /opt/views /opt/views

RUN { \
      echo '#!/bin/sh' \
      && echo '.' /opt/spack-environment/activate.sh \
      && echo 'exec "$@"'; \
    } > /entrypoint.sh \
&& chmod a+x /entrypoint.sh \
&& ln -s /opt/views/view /opt/view


LABEL "app"="self"
LABEL "mpi"="openmpi@5.0.2"
LABEL "platform"="x86_64"
LABEL "gpu"="none"
LABEL "hdf5"="1.14.3"
LABEL "feqparse"="2.2.2"
ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash" ]

