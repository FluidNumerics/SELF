# SELF Container Images

SELF is distributed as Docker images from a registry hosted on Google Cloud. You can download and deploy SELF using Docker or Apptainer.

!!! note
    Currently access is granted to the Docker image registry is granted on a case-by-case basis. Email [support@fluidnumerics.com](mailto:support@fluidnumerics.com) to request access.


## Prerequisites

Fluid Numerics' Docker image registry for SELF is hosted on Google Cloud. Access is granted using a GMail or Google Workspace email address. To download Docker images from the registry, you will need to have the following installed on your workstation

* [gcloud CLI](https://cloud.google.com/sdk/docs/install) : After installing you will need to run `gcloud init` to authenticate using your GMail or Google Workspace email address.
* [Docker](https://docs.docker.com/engine/install/) or [Apptainer](https://apptainer.org/docs/admin/main/installation.html) (formerly Singularity)


Once installed and initialized, you will need to configure your Docker credentials (even if you are using Apptainer) using the gcloud CLI. You can do this by running the following in your workstation's terminal :

```
gcloud auth configure-docker \
    us-docker.pkg.dev
```

## Selecting an image

SELF is able to run on both Nvidia and AMD GPUs. This capability is enabled by AMD's open source HIPFort and HIP projects. Due to the nature of HIP and HIPFort, the target GPU architecture is determined at compile time. Additionally, container images need to have matching CUDA or ROCm versions as the host system they will be deployed on. 

If you are deploying the image on an AMD GPU platform, you will need to select a container image that has a matching ROCm version and target GPU architecture. Alternatively, if you are deploying the image on an Nvidia GPU platform, you will need to select a container image that has a matching CUDA version and target GPU architecture. If you plan to run SELF without GPU acceleration, you can use any of the container images available. 


The SELF images are hosted in the `us-docker.pkg.dev/self-fluids/self` registry; each image follows a naming convention ` base_rocm-{rocm_version}_cuda-{cuda_version}_{gpu_arch}_{prec}` that identifies the following

* `{rocm_version}` : The ROCm version
* `{cuda_version}` : The CUDA version
* `{gpu_arch}` : The target GPU architecture
* `{prec}` : The floating point precision (either `single` or `double`)

As an example,

```
apptainer pull self.sif docker://us-docker.pkg.dev/self-fluids/self/base_rocm-5.4.6_cuda-11.8.0_gfx906_double
```

will pull the latest version of SELF that was built with ROCm 5.4.6, targeting the gfx906 (MI50) GPU platform using, with double precision floating point arithmetic enabled. The resulting image will be stored in `self.sif`.

<center>
[Browse available images](https://console.cloud.google.com/artifacts/docker/self-fluids/us/self?project=self-fluids){ .md-button .md-button--primary }
</center>

### Selecting a specific version
Each of the SELF images posted to the registry are tagged with the corresponding Git SHA. This allows you to pinpoint a specific version of SELF, if you so desire. You can reference a specific version by using the `:{git_sha}` suffix on the image. For example

```
apptainer pull self.sif docker://us-docker.pkg.dev/self-fluids/self/base_rocm-5.4.6_cuda-11.8.0_gfx906_double:
```

## GPU Architectures Reference

  Vendor   |  GPU Model   |  Architecture Code  |
---------- | ------------ | ------------------- |
  AMD      |    MI25      |     gfx900          |
  AMD      |    MI50      |     gfx906          |
  AMD      |    MI100     |     gfx908          |
  AMD      |    MI210     |     gfx90a          |
  AMD      |    MI250     |     gfx90a          |
  AMD      |    MI250X    |     gfx90a          |
  Nvidia   |    P100      |     sm_62           |
  Nvidia   |    V100      |     sm_72           |
  Nvidia   |    A100      |     sm_86           |

## Running with Docker

```
docker run --mount type=bind,source="$(pwd)/input",target="/self" \
           self:test \
           "/opt/self/bin/self -i /opt/self/etc/schema/defaults/CompressibleIdealGas2D.json"
```