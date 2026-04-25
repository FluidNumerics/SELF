ARG SELFISH_IMAGE=higherordermethods/selfish:latest-x86-cuda130-sm100
FROM ${SELFISH_IMAGE}

ARG SELFISH_SHA=""
LABEL org.opencontainers.image.base.name="${SELFISH_IMAGE}"
LABEL org.opencontainers.image.base.digest="${SELFISH_SHA}"

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
      -DCMAKE_CUDA_ARCHITECTURES="100" \
      -DSELF_ENABLE_EXAMPLES=ON \
      /opt/self/src && \
    make -j$(nproc) && \
    make install

WORKDIR /opt/self/build
