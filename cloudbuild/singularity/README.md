# Singularity with gcloud builds

This directory contains files for creating a singularity build step in for Google Cloud Builds ( gcloud builds )
For developers looking to build new SELF-Fluids Singularity containers on Google Cloud Platform, creating
the singularity build step in the first step.

## Prerequisites
* [gcloud SDK](https://cloud.google.com/sdk/)

## Setup
It is recommended that you run these commands from Google Cloud Shell within the project you plan to 
develop SELF-Fluids. Alternatively, you can use the gcloud SDK on your local system. In either case
set your project with
```
gcloud config set project <PROJECT ID>
```
where `<PROJECT ID>` is replaced with the desired project-id.

Then, build the build step
```
gcloud builds submit --config=cloudbuild.yaml
```

After successfully building the build step, you can find build meta-data with
```
gcloud builds list
```
## Usage
This build step invokes singularity`commands in [Google Cloud Build](cloud.google.com/cloud-build/).

Arguments passed to this builder will be passed to `singularity` directly,
allowing callers to run [any singularity command](https://www.sylabs.io/guides/3.0/user-guide/cli.html).

