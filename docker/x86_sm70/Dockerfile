ARG SELFISH_IMAGE=higherordermethods/selfish:latest-x86-cuda124-sm70
FROM ${SELFISH_IMAGE}

ARG SELFISH_SHA=""
LABEL org.opencontainers.image.base.name="${SELFISH_IMAGE}"
LABEL org.opencontainers.image.base.digest="${SELFISH_SHA}"

# Floating-point precision of the build. ON -> double precision (default),
# OFF -> single precision. Maps to the SELF_ENABLE_DOUBLE_PRECISION CMake option.
ARG SELF_DOUBLE_PRECISION=ON

COPY . /opt/self/src

RUN source /opt/spack-environment/activate.sh && \
    mkdir -p /opt/self/build && \
    cd /opt/self/build && \
    FC=gfortran cmake \
      -DCMAKE_INSTALL_PREFIX=/opt/self/install \
      -DCMAKE_BUILD_TYPE=Release \
      -DSELF_ENABLE_GPU=ON \
      -DSELF_GPU_BACKEND=CUDA \
      -DSELF_ENABLE_TESTING=ON \
      -DSELF_ENABLE_DOUBLE_PRECISION=${SELF_DOUBLE_PRECISION} \
      -DCMAKE_CUDA_ARCHITECTURES="70" \
      -DSELF_ENABLE_EXAMPLES=ON \
      /opt/self/src && \
    make -j$(nproc) && \
    make install

WORKDIR /opt/self/build
