# Install SELF
The Spectral Element Library in Fortran can be built provided the following dependencies are met

* Fortran 2008 compliant compiler
* MPI, e.g. [OpenMPI](https://www.open-mpi.org/)
* [GNU Make](https://www.gnu.org/software/make/)
* Fortran compiler ( `gfortran` recommended )
* [ROCm](https://rocm.docs.amd.com/en/latest/deploy/linux/quick_start.html), specifically the `rocm-hip-libraries`.
* (Optional) [CUDA Toolkit](), if you are building for Nvidia GPU hardware.
* [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
* [FluidNumerics/FEQParse](https://github.com/FluidNumerics/feq-parse)
* [jacobwilliams/JSON-Fortran](https://github.com/jacobwilliams/json-fortran)

!!! note
    Since ROCm is officially supported only on CentOS, RHEL, Ubuntu, and SLES operating systems, SELF currently can only be built on these operating systems. For deployment on Windows or MacOS systems, consider using the containerized builds.

You can install SELF in two possible ways

1. Bare Metal Installation
2. Docker image

A bare metal installation will require that you have a CentOS/RHEL 7 or 8, SLES, or Ubuntu 20.04 (focal) or 22.04 (jammy). You will also need to ensure that all of the dependencies are installed on your system before installing SELF.

A Docker image installation uses the Ubuntu 22.04 Docker image as a base and takes care of installing all of the dependencies for you. The resulting Docker image can be run using Docker or Singularity/Apptainer. 

This documentation will walk through steps to install SELF using bare metal installation and the Docker image approaches.


## Bare Metal Install

### Dependency installation

!!! note
    A bare metal installation will require that you have a CentOS/RHEL 7 or 8, SLES, or Ubuntu 20.04 (focal) or 22.04 (jammy).


On your system, make sure that you have a 2008 standard compliant Fortran compiler and GNU Make installed. On Ubuntu,

```
sudo apt-get install build-essential gcc gfortran
```

To help install the remainder of SELF's dependencies, we recommend that you use [Spack](https://spack.io). SELF comes with a Spack environment file that can be used to create a spack environment on your system.

First, install Spack

```
git clone https://github.com/spack/spack ~/spack
source ~/spack/share/spack/setup-env.sh
```

You can instruct Spack to use packages that are included with your operating system. This can help speed up the installation process.

```
spack external find --not-buildable
```

Download the SELF source code and activate the spack environment

```
git clone https://github.com/fluidnumerics/SELF ~/SELF
spack env activate -d ~/SELF/env
```

**Describe the SELF environment file and some modifications folks may want to make**

Once the environment is activated, show which packages will be installed and verify the output appears like what is shown below

```
$ spack find
==> In environment /home/joe/SELF/env
==> Root specs
feq-parse@1.1.0  flap@master  hdf5@1.12.0 +cxx+fortran+mpi  hipfort@4.5.2  json-fortran@7.1.0

==> 0 installed packages
```


Next, you can install the dependencies

```
spack install
```


### Install SELF

**Switch to make commands here and document the variables that control the build**
```
cd ~/SELF
./install.sh
```

This will install SELF under `${HOME}/view/self`. By default, this script is configured to build with serial, MPI, and GPU support with the target GPU set to AMD MI100 (`gfx900`). To change the behavior of the installation script, you can set the following environment variables before calling the script

* `VIEW` - The path to the spack environment view.
* `SELF_PREFIX` - The path to install SELF. Defaults to `$VIEW`
* `GPU_TARGET` - GPU microarchitecture code to build for. Defaults to `gfx900` (AMD MI100)
* `PREC` - Floating point precision to build with. Defaults to `double`. Change to `single` to build using 32-bit floating point arithmetic.
* `SELF_FFLAGS` - compiler flags to build SELF.


## Build a Docker Container
SELF comes with a Dockerfile to create builds that target specific GPU platforms. To build Docker containers, you will need to install [Docker](https://www.docker.com/). 

### Build with Cloud-Build-Local (Recommended)
To build a Docker container with SELF pre-installed, the SELF repository comes with a Cloud Build pipeline for use on your local system. This pipeline will execute `docker run` with the appropriate Dockerfile, depending on the target GPU architecture specified in the build substitutions. To use cloud-build-local, you will need to install [Docker](https://www.docker.com/), the [gcloud CLI, and google-cloud-sdk-cloud-build-local](https://cloud.google.com/sdk/docs/install).

Once installed, you can simply build SELF using the following command from the root of the SELF repository.

```
cloud-build-local --config=.cloudbuild/local.yaml --dryrun=false .
```

By default, this will build SELF with double precision floating point arithmetic, no optimizations (debug build), and with GPU kernels offloaded to Nvidia V100 GPUs. You can customize the behavior of the build process by using build substitutions. The following build substitution variables are currently available

* `_PREC` : The floating point precision to use in SELF; either `single` or `double`
* `_GPU_TARGET`: GPU microarchitecture code to build for. Defaults to `sm_72` (Nvidia V100)
* `_HIP_PLATFORM`: The value to set for the `HIP_PLATFORM` environment variable. Either `nvidia` or `amd`
* `_FFLAGS` : The compiler flags to send to the fortran compiler.

As an example, you can specify these substitution variables using something like the following

```
cloud-build-local --config=ci/cloudbuild.local.yaml --dryrun=false . --substitutions=_PREC=single,_GPU_TARGET=gfx906,_HIP_PLATFORM=amd
```
