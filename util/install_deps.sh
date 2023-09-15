#!/bin/bash
#
#
# Maintainers : joe@fluidnumerics.com (@fluidnumerics_joe)
#
# //////////////////////////////////////////////////////////////// #

# Set any package variables here
SPACK_VERSION="${SPACK_VERSION:-v0.19.2}"
INSTALL_ROOT="${INSTALL_ROOT:-/opt}"
ROCM_VERSION="${ROCM_VERSION:-5.2.3}"
###

OS=$(cat /etc/os-release | grep ^ID= | awk -F "=" '{print $2}')
echo "Operating System : $OS"

if [[ $OS == "centos" ]];then

   PKGMGR=yum

elif [[ $OS == "debian" ]];then

   echo "ERROR : ROCm dependency not supported on Debian OS"
   exit 1

elif [[ $OS == "ubuntu" ]];then

   PKGMGR=apt-get
   export DEBIAN_FRONTEND=noninteractive
   apt-get remove -y unattended-upgrades

elif [[ $OS == "pop" ]];then

   echo "ERROR : ROCm dependency not supported on POP! OS"
   exit 1

fi

$PKGMGR update -y

# Install dependencies
if [[ $OS == "centos" ]];then
   $PKGMGR install -y gcc gcc-c++ gcc-gfortran git make wget patch valgrind valgrind-devel

   # CUDA
#   wget https://developer.download.nvidia.com/compute/cuda/11.8.0/local_installers/cuda_11.8.0_520.61.05_linux.run
#   sh cuda_11.8.0_520.61.05_linux.run --driver --toolkit --silent
 
   # ROCM
   yum install -y epel-release devtoolset-7 cmake

   cat > /etc/yum.repos.d/rocm.repo <<EOL
[rocm]
name=rocm
baseurl=https://repo.radeon.com/rocm/yum/rpm
enabled=1
gpgcheck=1
gpgkey=https://repo.radeon.com/rocm/rocm.gpg.key
EOL

   cat > /etc/yum.repos.d/amdgpu.repo <<EOL
[amdgpu]
name=amdgpu
baseurl=https://repo.radeon.com/amdgpu/latest/rhel/7.9/main/x86_64
enabled=1
gpgcheck=1
gpgkey=https://repo.radeon.com/rocm/rocm.gpg.key
EOL

   yum clean all -y
   yum update -y
   yum install -y rocm-dev

   cat > /etc/profile.d/z11_rocm.sh <<EOL
#!/bin/bash
export PATH=\${PATH}:/opt/rocm/bin
export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:/opt/rocm/lib:/opt/rocm/lib64
EOL


elif [[ $OS == "ubuntu" ]];then
   $PKGMGR install -y build-essential git wget libnuma-dev python3-dev python3-pip zip unzip valgrind gcovr

   # CUDA
#   wget https://developer.download.nvidia.com/compute/cuda/11.8.0/local_installers/cuda_11.8.0_520.61.05_linux.run
#   sh cuda_11.8.0_520.61.05_linux.run --driver --toolkit --silent

   # ROCM
   wget -q -O - https://repo.radeon.com/rocm/rocm.gpg.key | apt-key add - && \
   echo "deb [arch=amd64] https://repo.radeon.com/rocm/apt/${ROCM_VERSION}/ ubuntu main" | tee /etc/apt/sources.list.d/rocm.list &&\
   apt-get -yqq update && \
   apt-get -yqq install rocm-dev

elif [[ $OS == "pop" ]];then
   $PKGMGR install -y build-essential git wget libnuma-dev python3-dev python3-pip zip unzip valgrind gcovr

   # CUDA
#   wget https://developer.download.nvidia.com/compute/cuda/11.8.0/local_installers/cuda_11.8.0_520.61.05_linux.run
#   sh cuda_11.8.0_520.61.05_linux.run --driver --toolkit --silent

   # ROCM
   wget -q -O - https://repo.radeon.com/rocm/rocm.gpg.key | apt-key add - && \
   echo "deb [arch=amd64] https://repo.radeon.com/rocm/apt/${ROCM_VERSION}/ ubuntu main" | tee /etc/apt/sources.list.d/rocm.list &&\
   apt-get -yqq update && \
   apt-get -yqq install rocm-dev
fi


# Install spack
if [[ ! -d ${INSTALL_ROOT}/spack ]]; then
   git clone https://github.com/spack/spack.git ${INSTALL_ROOT}/spack
fi

source ${INSTALL_ROOT}/spack/share/spack/setup-env.sh
spack compiler find --scope site

# Install packages specified in the spack environment
cat ./env/spack.yaml
spack env activate -d ./env/
spack install --fail-fast --source
spack gc -y
spack env deactivate
spack env activate --sh -d ./env/ >> /etc/profile.d/z10_spack_environment.sh 
