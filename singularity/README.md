# Base Builders
The base builders are singularity containers that have all of the prerequisities for
building SELF-Fluids. The base builders represent the supported compilers, MPI flavors and version,
and HDF5 versions that we support.

## Getting started
The base builders can be build on any platform that has singularity. There is support for building
on Google Cloud Platform with `gcloud builds` through cloudbuild yaml files. To build with cloud builds,
you will need to create a bucket as follows
```
gsutil mb gs://${PROJECT_ID}-singularity
```
All of the base builder Singularity Image Files (.sif) will be dumped to this location.
