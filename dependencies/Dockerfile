FROM rocm/dev-ubuntu-20.04:latest AS devel

RUN DEBIAN_FRONTEND=noninteractive apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y git cmake gfortran wget

# Install CUDA Toolkit
RUN apt-get update &&\
    wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin &&\
    mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600 &&\
    wget https://developer.download.nvidia.com/compute/cuda/11.1.1/local_installers/cuda-repo-ubuntu1804-11-1-local_11.1.1-455.32.00-1_amd64.deb &&\
    dpkg -i cuda-repo-ubuntu1804-11-1-local_11.1.1-455.32.00-1_amd64.deb &&\
    apt-key add /var/cuda-repo-ubuntu1804-11-1-local/7fa2af80.pub &&\
    DEBIAN_FRONTEND=noninteractive apt-get update -y &&\
    DEBIAN_FRONTEND=noninteractive apt-get install -y cuda

# Install hipfort
RUN git clone  https://github.com/ROCmSoftwarePlatform/hipfort.git /tmp/hipfort &&\
    mkdir -p /tmp/hipfort/build &&\
    cd /tmp/hipfort/build &&\
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local /tmp/hipfort &&\
    make &&\
    make install



ENV PATH=${PATH}:/usr/local/hipfort/bin \
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/hipfort/lib
