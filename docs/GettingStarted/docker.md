# SELF Docker Images

## Fluid Numerics Artifact Registry
Fluid Numerics maintains a registry of Docker images on Google Cloud. We provide access to this registry through a subscription to the [Fluid Numerics Higher Order Methods Lab (homlab)](https://www.fluidnumerics.com/shop/p/higher-order-methods-lab). Once you have a [homlab subscription](https://www.fluidnumerics.com/shop/p/higher-order-methods-lab), you will have the ability to pull Docker images that we maintain and test regularly.

When you sign up for the subscription you will need a Gmail, Google Workspace, or Cloud Identity account that you use for Google Cloud.

## Build your own images
The SELF repository comes with docker recipes for CPU-only, ROCm, and CUDA platforms. Each recipe is organized under the `docker/` subdirectory:

* `docker/x86_64` - CPU only, MPI enabled, multithreaded
* `docker/x86_64_rocm` - MPI+HIP enabled, multithreaded
* `docker/x86_64_cuda` - MPI+CUDA enabled, multithreaded

Each subdirectory contains a `spack.yaml` and `Dockerfile`. The `spack.yaml` file is used to generate the Dockerfile using `spack containerize`. The `spack.yaml` files are retained for development purposes and are often not used by end-users. 

### CPU-only platforms
To build a Docker image for CPU only platforms :

1. Clone the SELF repository

```shell
git clone https://github.com/fluidnumerics/self ~/self/
cd ~/self/
```

2. Build the image with `docker build`. This process installs all of SELF's dependencies and SELF and usually takes about 30 minutes to complete.
```
docker build -f docker/x86_64/Dockerfile -t self-x86_64:latest .
```