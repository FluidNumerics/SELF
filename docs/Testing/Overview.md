# Testing

The SELF repository provides mechanisms for you to be able to build and test SELF on your own system. The same mechanisms are also replicated in Google Cloud so that tests can be run on a variety of platforms before changes are made to the main branch.


## Pull Requests
SELF is tested when pull requests are submitted to update the main branch and after a repository owner approves a build. Tests are run using Google Cloud Build; the workflow is defined in `ci/cloudbuild.yaml`. The process builds a docker container image by installing SELF and its dependencies, including ROCm, CUDA, HDF5, FLAP, JSON-Fortran, and OpenMPI. 

Once built, the docker image is pushed to Google Container registry and the docker image is tested using Google Cloud Batch. Batch provides access to a variety of compute resources, including GPUs. This enables us to test serial, MPI, and GPU accelerated components of SELF.
