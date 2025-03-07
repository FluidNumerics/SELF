# SELF Docker Images

## Fluid Numerics Artifact Registry
Fluid Numerics maintains a registry of Docker images on Google Cloud. We provide access to this registry through a subscription to the [Fluid Numerics Higher Order Methods Lab (homlab)](https://www.fluidnumerics.com/shop/p/higher-order-methods-lab). Once you have a [homlab subscription](https://www.fluidnumerics.com/shop/p/higher-order-methods-lab), or if you are a [collaborator](https://opencollective.com/opensource-fluidnumerics/projects/spectral-element-library-in-fo) you will have the ability to pull Docker images that we maintain and test regularly.

### Available Images
| Image           | compiler | MPI flavor | GPU Software | Target GPU | Target CPU |
| --------------- | -------- | ------------- | --- | ----| ------ |
|**`self-x86_64`**| gfortran | openmpi@5.0.2 | N/A | N/A | x86_64 |
|**`self-x86_64-gfx90a_rocm6.2.1`**| gfortran | openmpi@5.0.2 | rocm@6.0.2 | gfx90a (MI200) | x86_64 |
|**`self-x86_64-gfx942_rocm6.2.1`**| gfortran | openmpi@5.0.2 | rocm@6.0.2 | gfx942 (MI300) | x86_64 |
|**`self-x86_64-sm72_cuda12.1`**| gfortran | openmpi@5.0.2 | cuda@12.1 | sm72 (V100) | x86_64 |
|**`self-x86_64-sm90_cuda12.1`**| gfortran | openmpi@5.0.2 | cuda@12.1 | sm90 (H100) | x86_64 |


This guide provides instructions for pulling our Docker images hosted on Google Cloud Artifact Registry using **Docker** and **Singularity/Apptainer**.


### **Prerequisites**

#### General Requirements
- Ensure you have access credentials for the Google Cloud project hosting the Docker image.

#### Docker-Specific Requirements
- [Docker](https://www.docker.com/) installed on your system (version 20.10 or later recommended).
- A Google Cloud service account or personal account with permissions to access the Artifact Registry.

#### Singularity/Apptainer-Specific Requirements
- [Singularity](https://sylabs.io/singularity/) or [Apptainer](https://apptainer.org/) installed on your system (version 3.8 or later recommended).
- Docker CLI available to pull images before converting them for use with Singularity/Apptainer.


### Pulling the Docker Image with Docker

1. **Authenticate with Google Cloud Artifact Registry**  
   Run the following command to authenticate Docker with Google Cloud:
   ```bash
   gcloud auth configure-docker us-docker.pkg.dev
   ```
   This configures Docker to use your Google Cloud credentials for accessing private images.

2. **Pull the Docker Image**  
   Use the `docker pull` command to fetch the image:
   ```bash
   docker pull us-docker.pkg.dev/fluidnumerics-research/self/self-x86_64:latest
   ```
   Replace `latest` with the desired tag if applicable.

3. **Verify the Image**  
   Confirm the image has been pulled successfully by listing it:
   ```bash
   docker images
   ```


### Pulling and Using the Image with Singularity/Apptainer

#### Option 1: Direct Pull Using `singularity pull`

Singularity/Apptainer can directly pull the image and convert it into a `.sif` file.

1. **Pull the Image**  
   Use the `singularity pull` or `apptainer pull` command:
   ```bash
   singularity pull docker://us-docker.pkg.dev/fluidnumerics-research/self/self-x86_64:latest
   ```
   This will fetch the image and save it as `self-x86_64_latest.sif` in your current directory.

2. **Verify the SIF File**  
   Check that the file was created:
   ```bash
   ls -lh self-x86_64_latest.sif
   ```

3. **Run the Image**  
   Execute commands using the Singularity container:
   ```bash
   singularity exec self-x86_64_latest.sif <your-command>
   ```

#### Option 2: Pull with Docker and Convert Locally

If you prefer to use Docker first and then convert the image for Singularity/Apptainer:

1. **Pull the Docker Image**  
   Follow the Docker instructions above to pull the image.

2. **Save the Docker Image as a `.tar` File**  
   Export the Docker image to a tarball:
   ```bash
   docker save us-docker.pkg.dev/fluidnumerics-research/self/self-x86_64:latest -o self-x86_64.tar
   ```

3. **Convert to a Singularity SIF File**  
   Use `singularity build` or `apptainer build` to create a SIF file:
   ```bash
   singularity build self-x86_64.sif docker-archive://self-x86_64.tar
   ```

4. **Run the Image**  
   Use the `singularity exec` command:
   ```bash
   singularity exec self-x86_64.sif <your-command>
   ```


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