#############
Build SELF
#############
The Spectral Element Library in Fortran can be built provided the following dependencies are met

* Fortran 2008 compliant compiler
* MPI, e.g. `OpenMPI <https://www.open-mpi.org/>`_
* `GNU Make <https://www.gnu.org/software/make/>`_
* `ROCm <https://rocmdocs.amd.com/en/latest/Installation_Guide/Installation-Guide.html>`_ ( >= 4.2 )
* `HIPFort <https://github.com/ROCmSoftwarePlatform/hipfort>`_ ( >= 4.2 ; Must be built with the Fortran compiler you will build SELF with. )
* `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_
* `FEQParse <https://github.com/FluidNumerics/feq-parse>`_
* `FLAP <https://github.com/szaghi/FLAP>`_
* `JSON-Fortran <https://github.com/jacobwilliams/json-fortran>`_


======================
Bare Metal Install
======================

Install SELF Dependencies
==========================
We recommend that you `install ROCm using using the instructions found in AMD's rocmdocs <https://rocmdocs.amd.com/en/latest/Installation_Guide/Installation-Guide.html>`_.

To help install the remainder of SELF's dependencies, we recommend that you use `Spack <https://spack.io>`_. SELF comes with a Spack environment file that can be used to create a spack environment on your system.

First, install Spack

.. code-block:: shell

   git clone https://github.com/spack/spack ~/spack
   source ~/spack/share/spack/setup-env.sh


Download the SELF source code and activate the spack environment

.. code-block:: shell

   git clone https://github.com/fluidnumerics/SELF ~/SELF
   spack env activate -d ~/SELF/env


Once the environment is activated, show which packages will be installed and verify the output appears like what is shown below


.. code-block:: shell

   $ spack find
   ==> In environment /home/joe/SELF/env
   ==> Root specs
   feq-parse@1.1.0  flap@master  hdf5@1.12.0 +cxx+fortran+mpi  hipfort@4.5.2  json-fortran@7.1.0
   
   ==> 0 installed packages

Next, you can install the dependencies

.. code-block:: shell

   $ spack install

Keep in mind, this installation process can take up to two hours.

Install SELF
=============
SELF comes with a `simple bash script <https://github.com/FluidNumerics/SELF/blob/main/install.sh>`_ that can be used to install SELF. Provided you followed the steps above to install SELF dependencies, you can simply run this script.


.. code-block:: shell

   $ cd ~/SELF
   $ ./install.sh

This will install SELF under :code:`${HOME}/view/self`. By default, this script is configured to build with serial, MPI, and GPU support with the target GPU set to AMD MI100 (:code:`gfx900`). To change the behavior of the installation script, you can set the following environment variables before calling the script

* :code:`VIEW` - The path to the spack environment view.
* :code:`SELF_PREFIX` - The path to install SELF. Defaults to :code:`$VIEW`
* :code:`GPU_TARGET` - GPU microarchitecture code to build for. Defaults to :code:`gfx900` (AMD MI100)
* :code:`PREC` - Floating point precision to build with. Defaults to :code:`double`. Change to :code:`single` to build using 32-bit floating point arithmetic.
* :code:`SELF_FFLAGS` - compiler flags to build SELF.


===========================
Build a Docker Container
===========================
SELF comes with Dockerfiles to create builds that target specific GPU platforms. To build Docker containers, you will need to install `Docker <https://www.docker.com/>`_. The SELF build process for Docker containers involves two steps

1. Build a base container with all of the necessary dependencies
2. Build SELF for the target architecture using the appropriate base container.


We've opted to split the build process into these two distinct stages, since installing SELF's dependencies takes 1-2 hours. Further, the Docker container builds are used as part of the continuous integration process and this approach reduces the costs associated with Google Cloud Build.

Installing Base Containers
============================
Base containers are hosted by Fluid Numerics at :code:`gcr.io/self-fluids/self-base`. Tags are used to mark the target platform in addition to the ROCm version and CUDA version (if applicable). For example, :code:`gcr.io/self-fluids/self-base:amd-rocm4.3` is the base container used for build SELF for AMD GPU's using ROCm 4.3. Below is a list of currently available base container

* :code:`gcr.io/self-fluids/self-base:amd-rocm4.3`
* :code:`gcr.io/self-fluids/self-base:nvidia-cuda11.2.1-rocm4.3`

To access these base containers, you need to be a `"SELF Member" under HigherOrderMethods <https://opencollective.com/higher-order-methods/contribute/self-member-35220/checkout>`_

Alternatively, you will need to build base container image. Container recipes for base containers can be found under :code:`SELF/docker/base/{HIP_PLATFORM}`, where :code:`{HIP_PLATFORM}` is one of :code:`nvidia` or :code:`amd`. To install the base containers locally, you can run

.. code-block:: shell

    docker build . -f docker/base/nvidia/Dockerfile -t gcr.io/self-fluids/self-base:nvidia-cuda11.2.1-rocm4.3
    docker build . -f docker/base/nvidia/Dockerfile -t gcr.io/self-fluids/self-base:amd-rocm4.3



Build with Cloud-Build-Local (Recommended)
============================================
To build a Docker container with SELF pre-installed, the SELF repository comes with a Cloud Build pipeline for use on your local system. This pipeline will execute :code:`docker run` with the appropriate Dockerfile, depending on the target GPU architecture specified in the build substitutions. To use cloud-build-local, you will need to install `Docker <https://www.docker.com/>`_, the `gcloud CLI, and google-cloud-sdk-cloud-build-local <https://cloud.google.com/sdk/docs/install>`_.

Once installed, you can simply build SELF using the following command from the root of the SELF repository.

.. code-block:: shell

   $ cloud-build-local --config=ci/cloudbuild.local.yaml --dryrun=false .

By default, this will build SELF with double precision floating point arithmetic and with GPU kernels offloaded to Nvidia V100 GPUs. You can customize the behavior of the build process by using build substitutions. The following build substitution variables are currently available

* :code:`_PREC` : The floating point precision to use in SELF; either :code:`single` or :code:`double`
* :code:`_GPU_TARGET`: GPU microarchitecture code to build for. Defaults to :code:`sm_72` (Nvidia V100)
* :code:`_HIP_PLATFORM`: The value to set for the :code:`HIP_PLATFORM` environment variable. Either `nvidia` or `amd`

As an example, you can specify these substitution variables using something like the following


.. code-block:: shell

   $ cloud-build-local --config=ci/cloudbuild.local.yaml --dryrun=false . --substitutions=_PREC=single,_GPU_TARGET=gfx900,_HIP_PLATFORM=amd


The Dockerfile that :code:`cloud-build-local` uses for the build recipe is under :code:`SELF/docker/${_HIP_PLATFORM}/Dockerfile`. The resulting Docker image when using :code:`cloud-build-local` is :code:`self:dev-${_PREC}-${_GPU_TARGET}`


Build with Docker
==================
If you prefer to build SELF using Docker directly, you can easily do so using the Dockerfile recipes under the :code:`SELF/docker/{HIP_PLATFORM}` directory. The :code:`{HIP_PLATFORM}` is either :code:`nvidia` or :code:`amd`. 

  args: ['build',
         '--build-arg','GPU_TARGET=${_GPU_TARGET}',
         '--build-arg','HIP_PLATFORM=${_HIP_PLATFORM}',
         '--build-arg','PREC=${_PREC}',
         '-f',
         'docker/${_HIP_PLATFORM}/Dockerfile',
         '.',
         '-t',
         'self:dev-${_PREC}-${_GPU_TARGET}']
