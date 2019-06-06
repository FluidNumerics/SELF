# cloudbuild for SELF-Fluids

This directory was set up for working with [Google Cloud Build](cloud.google.com/cloud-build/). To use Cloud Build, you must have the [gcloud SDK](https://cloud.google.com/sdk/) installed.
Cloud Build allows developers to build containers on private and secure resources hosted by Google Cloud Platform

## Directory Structure
SELF-Fluids depends on
* A fortran compiler
* HDF5 ( https://www.hdfgroup.org/solutions/hdf5/ )
* METIS ( http://glaros.dtc.umn.edu/gkhome/metis/ )
* An MPI flavor [Optional] (mvapich or openmpi preferred)
For GPU accelerated builds of SELF-Fluids, PGI compilers should be used for CUDA-Fortran

### Subdirectories
This directory contains subdirectories with Dockerfiles and cloudbuild.yaml instructions for building containers with all of SELF-Fluids dependencies.

#### singularity/
Contains Dockerfile and cloudbuild.yaml that installs Singularity on top of Google's Go container. The singularity container is used for building
of SELF-Fluids as Singularity images.

See https://cloud.google.com/community/tutorials/singularity-containers-with-cloud-build for more details.


#### pgi/
Contains Dockerfile and cloudbuild.yaml that installs
* PGI compilers
* HDF5 built for serial IO with PGI compilers
* METIS built with gcc
* Singularity definition file (self-fluids.def) that contains instructions for building SELF-Fluids with GPU acceleration from the associated Docker container
 
## Workflow

**1. Set your GCloud Project**

It is recommended that you run these commands from Google Cloud Shell within the project you plan to 
develop SELF-Fluids. Alternatively, you can use the gcloud SDK on your local system. In either case
set your project with
```
gcloud config set project <PROJECT ID>
```
where `<PROJECT ID>` is replaced with the desired project-id.


**2. Create Singularity Build Step**
Navigate to the singularity subdirectory
```
$ cd singularity
```
Then, build the build step
```
gcloud builds submit --config=cloudbuild.yaml
```
**3. Build a dependency container**
The dependency containers, after built, contain all of the dependencies necessary to building SELF-Fluids.
For example, to build the Docker image with PGI compilers, HDF5 (serial), and METIS
```
$ cd pgi/
$ gcloud builds submit . --timeout=2h
```
The 2 hour timeout is needed, particularly for building the HDF5 library from source code.
This step creates a Docker image with SELF-Fluids dependencies installed. This Docker image is stored in your project's private [Google Container Registry](https://cloud.google.com/container-registry/)


**4. Build SELF-Fluids**
Once a dependency container is built, you can then proceed with building SELF-Fluids into a singularity image from the root directory of this repository. 
From the root directory, 
```
$ gcloud builds submit . --substitutions=_BUILD_BASE=pgi
```
where `_BUILD_BASE` is an environment variable that controls which dependency container you want to build from.


