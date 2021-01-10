#!/bin/bash

INSTALL_ROOT="/opt/self"

## Install the GCC-9 devtoolset
yum -y update
yum install -y centos-release-scl-rh
yum install -y devtoolset-9-toolchain
source /opt/rh/devtoolset-9/enable
cat > /etc/profile.d/gcc-9.sh <<EOL
#!/bin/bash
source /opt/rh/devtoolset-9/enable
EOL

## Install CMake 3.18.4
wget https://github.com/Kitware/CMake/releases/download/v3.18.4/cmake-3.18.4.tar.gz --directory-prefix=/tmp
cd /tmp
tar -xvzf cmake-3.18.4.tar.gz
cd /tmp/cmake-3.18.4
./bootstrap --prefix=/opt/cmake
make -j 8
make install

## Install hipfort
git clone https://github.com/ROCmSoftwarePlatform/hipfort.git /tmp/hipfort_src
cd /tmp/hipfort_src
git checkout origin/rocm-3.9.x -b rocm-3.9.x
mkdir /tmp/hipfort_src/build
cd /tmp/hipfort_src/build
/opt/cmake/bin/cmake -DHIPFORT_COMPILER=$(which gfortran) \
                     -DHIPFORT_INSTALL_DIR=/opt/hipfort \
                     ../
make
make install
cd /tmp/
rm -rf /tmp/hipfort_src

mkdir -p /tmp/extern

git clone https://github.com/FluidNumerics/feq-parse.git /tmp/extern/feq-parse && \
mkdir -p /tmp/extern/feq-parse/build && \
cd /tmp/extern/feq-parse/build && \
FC=/opt/rh/devtoolset-9/root/usr/bin/gfortran \
/opt/cmake/bin/cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_ROOT} /tmp/extern/feq-parse && \
make && make install

git clone https://github.com/jacobwilliams/json-fortran.git /tmp/extern/json-fortran && \
mkdir -p /tmp/extern/json-fortran/build && \
cd /tmp/extern/json-fortran/build && \
FC=/opt/rh/devtoolset-9/root/usr/bin/gfortran \
/opt/cmake/bin/cmake -DSKIP_DOC_GEN=True -DCMAKE_INSTALL_PREFIX=${INSTALL_ROOT} /tmp/extern/json-fortran && \
make && make install

git clone --recurse-submodules https://github.com/szaghi/FLAP.git /tmp/extern/FLAP && \
mkdir -p /tmp/extern/FLAP/build && \
cd /tmp/extern/FLAP/build && \
FC=/opt/rh/devtoolset-9/root/usr/bin/gfortran \
FFLAGS=-cpp /opt/cmake/bin/cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_ROOT} /tmp/extern/FLAP && \
make && make install

export HIP_PLATFORM=nvcc
export HIP_COMPILER=/usr/local/cuda/bin/nvcc

cd /tmp && make install -f /tmp/build/self.make

rm -rf /tmp/extern
