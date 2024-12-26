# Install SELF


## Install with Spack
The easiest way to get started is to use the spack package manager. On a Linux platform, set up spack

```
git clone https://github.com/spack/spack ~/spack
source ~/spack/share/spack/setup-env.sh
```

Allow spack to locate your compilers (make sure you have C, C++, and Fortran compilers installed!)

```
spack compiler find
```

SELF comes with a spack environment file that defines the dependencies that are required for SELF. The versions listed in this environment file are the specific versions we regularly test against. To get this environment file, clone the SELF repository

```
git clone https://github.com/fluidnumerics/SELF/ ~/SELF/
```

**If you have a preferred compiler** you would like for spack to use, you can use `spack config add`, e.g.

```
spack -e ~/SELF/share/spack-env config add packages:all:require:["'%gcc@12.2.0'"]
```

The example above will force packages to be built with version 12.2.0 of gfortran from the `gcc` compiler set. 

!!! note
    If you do not set a preferred compiler, spack will pick one based on the available compilers found using `spack compiler find`

To reduce build time, import existing packages on your system
```
spack external find --not-buildable
```

Next, install SELF's dependencies (OpenMPI, HDF5, and feq-parse)
```
spack -e ~/SELF/share/spack-env install --no-check-signature
```

Then, install SELF
```
cd ~/SELF
spack env activate ~/SELF/share/spack-env
mkdir ~/SELF/build
cd ~/SELF/build
cmake -DCMAKE_INSTALL_PREFIX=${HOME}/opt/self ../
make
make install
```

If you'd like to run the tests included with SELF, to verify your installation, you can use `ctest`.

```
cd ${HOME}/opt/self/test
ctest
```

### Once v0.0.1 is released 
The easiest way to get started is to use the [spack package manager](https://spack.io). The spack package manager provides you with an easy command line interface to install research software from source code with all of its dependencies. 

!!! note
    Before proceeding, you will need to ensure that you have a 2008 compliant Fortran compiler.

On a Linux platform, set up spack : 

```shell
git clone https://github.com/spack/spack ~/spack
source ~/spack/share/spack/setup-env.sh
```

Allow spack to locate your compilers (make sure you have C, C++, and Fortran compilers installed!)

```shell
spack compiler find
```

The example above will force packages to be built with version 12.2.0 of gfortran from the `gcc` compiler set.

To reduce build time, import existing packages on your system

```shell
spack external find --not-buildable
```

Next, install SELF and it's dependencies

```shell
spack install self
```

By default, this will install SELF with the following features
* Double precision floating point arithmetic
* No unit tests and no examples
* No multi-threading, CPU-only

You can view documentation on all possible variants using

```shell
spack info self
```

### Enable Multithreading
Many of the computationally intensive methods in SELF are written using the `do concurrent` structure. We have provided the variant `+multithreading` which will enable multithreading for all `do concurrent` blocks. You can install SELF with multithreading using

```shell
spack install self+multithreading
```

If you are using the GNU compiler suite, the number of threads used for `do concurrent` blocks is determined during build time. Because of this, we have provided the `nthreads` option, which defaults to 4. You can change this option to a value more sensible for your platform, e.g.

```shell
spack install self+multithreading nthreads=16 % gcc
```

The `%gcc` here indicates that you intend to build SELF with the GNU compilers.

### Enable Nvidia GPU Acceleration
SELF provides GPU accelerated implementations of all methods that are used in forward stepping conservation law solvers. On Nvidia GPU platforms, you can take advantage of this using the `+cuda` variant : 

```shell
spack install self+cuda
```
This will also ensure that the MPI flavor that is used is GPU aware. You can specify the GPU architecture using the `gpu_arch` build option, e.g. for A100 GPUs

```shell 
spack install self+cuda gpu_arch=sm_80
```

!!! note
    AMD GPU-Aware MPI is currently not available in Spack. This means that these steps will not allow you to build SELF for multi-GPU platforms with AMD GPUs. See [Advanced Installation](#advanced-installation) for details on how to install for AMD GPU platforms.


## Advanced Installation

### Dependencies
The Spectral Element Library in Fortran can be built provided the following dependencies are met

* [Cmake (v3.21 or greater)](https://cmake.org/resources/)
* Fortran 2008 compliant compiler ( `gfortran` recommended )
* MPI, e.g. [OpenMPI](https://www.open-mpi.org/) with GPU-Aware Support
* [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
* [FluidNumerics/feq-parse](https://github.com/FluidNumerics/feq-parse)
* (Optional, AMD GPU Support) [ROCm v6.0.2 or greater](https://rocm.docs.amd.com/projects/install-on-linux/en/latest/)
* (Optional, Nvidia GPU Support) [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit)

[**Learn more about installing SELF's dependencies**](./dependencies.md)

### Installing SELF from source with CMake

!!! warning
    It is assumed that you have all of the necessary dependencies installed on your system and that they are discoverable by CMake.

SELF comes with a CMake build system that defines the build and installation process. 

1. Clone the SELF repository

```shell
git clone https://github.com/fluidnumerics/self ~/self/
```

2. Create a build directory

```shell
mkdir ~/self/build
cd ~/self/build
```

3. Run CMake to build the make system. You can set `CMAKE_INSTALL_PREFIX` to a path where you'd prefer to install SELF; here, we set it to `${HOME}/opt/self`

```shell
cmake -DCMAKE_INSTALL_PREFIX=${HOME}/opt/self ../
```

4. Build SELF and run the test suite to ensure everything is built properly

```shell
make
ctest --output-on-failure
```

5. Install SELF 

```shell
make install
```

When you install SELF, you will install the following artifacts

* `${CMAKE_INSTALL_PREFIX}/lib/libself-static.a` - A static library for SELF
* `${CMAKE_INSTALL_PREFIX}/lib/libself.so` - A shared object library for SELF
* `${CMAKE_INSTALL_PREFIX}/include/*.mod` - Module files generated by the Fortran compiler during the build process.
* `${CMAKE_INSTALL_PREFIX}/example/*` - A set of example programs
* `${CMAKE_INSTALL_PREFIX}/test/*` - A set of unit tests 

By default, this will install SELF with the following features
* Double precision floating point arithmetic
* No unit tests and no examples
* No multi-threading, CPU-only

There are a few CMake options that you can set to control the build features : 
* `SELF_ENABLE_MULTITHREADING`:  Option to enable CPU multithreading for `do concurrent` loop blocks. (Default: OFF)
* `SELF_ENABLE_TESTING`:  Option to enable build of tests. (Default: ON)
* `SELF_ENABLE_EXAMPLES`: Option to enable build of examples. (Default: ON)
* `SELF_ENABLE_GPU`: Option to enable GPU backend. Requires either CUDA or HIP. (Default: OFF)
* `SELF_ENABLE_DOUBLE_PRECISION` Option to enable double precision for floating point arithmetic. (Default: ON)

### Enabling Multithreading CPU support
Computationally heavy methods in SELF are expressed using Fortran's `do concurrent` loop blocks, which gives compilers the freedom to parallelize operations. Every Fortran compiler has their own set of compiler flags to enable parallelization of `do concurrent` blocks (see [this post on the Fortran-Lang discourse](https://fortran-lang.discourse.group/t/do-concurrent-compiler-flags-to-enable-parallelization/4300/6)). We have provided a single option in the CMake build system that allow you to enable parallelization. At the `cmake` stage of the build process, you can set `SELF_ENABLE_MULTITHREADING=ON`, e.g.

```shell
cmake -DSELF_ENABLE_MULTITHREADING=ON \
      -DCMAKE_INSTALL_PREFIX=${HOME}/opt/self \
       ../
```

If you are building with the GNU compilers (`gfortran`), the number of threads used for parallelization is determined at build time. By default, the SELF build system will set the number of threads to 4. You can override this setting at the `cmake` stage of the build process using the `SELF_MULTITHREADING_NTHREADS` build variable, e.g. to set the number of threads to 8 with `gfortran`,

```shell
cmake -DSELF_ENABLE_MULTITHREADING=ON \
      -DSELF_MULTITHREADING_NTHREADS=8 \
      -DCMAKE_INSTALL_PREFIX=${HOME}/opt/self \
       ../
```

The CMake build system will set the appropriate flags for multithreading for GNU, Intel (`ifort` and `ifx`), LLVM, and Nvidia HPC Compilers. If you are not using `gfortran`, you can set the number of threads for parallelism at runtime using the `OMP_NUM_THREADS` environment variable

### Enabling GPU Support 
SELF offers the option to use HIP or CUDA. Some of our "heavy-lifting" kernels, such as divergence, gradient, and grid interpolation operations are expressed using the BLAS API. For these, we use HIPBLAS or CUBLAS. GPU support is enabled in the CMake stage of the build by setting `SELF_ENABLE_GPU=ON`

The CMake build system will automatically search for HIP. If HIP is not found, then it will search for CUDA. If neither is found, the build process will fail.

#### HIP
HIP can be used to build for either AMD or Nvidia GPU's. If you have HIP installed and it is found, you can also set the `CMAKE_HIP_ARCHITECTURES` build variable to specify which GPU architecture you want to build for. Alternatively, if you are building SELF on the system that has a GPU installed, you can let HIP auto-detect the available GPU. At this time, we also advise setting the `CXX` enviornment variable to `hipcc`,e.g.

```shell
CXX=hipcc \
cmake -DSELF_ENABLE_GPU=ON \
      -DCMAKE_INSTALL_PREFIX=${HOME}/opt/self \
       ../
```

#### CUDA
SELF provides you the ability to use CUDA directly, in case you are on a system that does not have AMD's ROCm and HIP installed. If you have CUDA installed and it is found,  you can also set the `CMAKE_CUDA_ARCHITECTURES` build variable to specify which Nvidia GPU architecture you want to build for. At this time, we also advise setting the `CXX` enviornment variable to `nvcc`,e.g.

```shell
CXX=nvcc \
cmake -DSELF_ENABLE_GPU=ON \
      -DCMAKE_INSTALL_PREFIX=${HOME}/opt/self \
       ../
```


#### Reference table for GPU architecture codes

Vendor | Model | Architecture code(s) |
------ | ----- | -------------------- |
AMD    | Instinct MI100 | gfx908  |
AMD    | Instinct MI210 | gfx90a  |
AMD    | Instinct MI250 | gfx90a  |
AMD    | Instinct MI250x | gfx90a  |
AMD    | Radeon Pro W7900 | gfx1100  |
AMD    | Radeon Pro W7800 | gfx1100  |
Nvidia | Volta (V100) | sm_70, sm_72 |
Nvidia | Ampere (A100) | sm_80, sm_86, sm_87 |
Nvidia | Hopper (H100) | sm_90, sm_90a |


[If you encounter any problems, feel free to open an new issue](https://github.com/FluidNumerics/SELF/issues/new/choose)