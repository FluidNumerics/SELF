terraform {
  backend "gcs" {
    bucket  = "self-fluids-terraform"
    prefix  = "self-build-triggers"
  }
}

// Configure the Google Cloud provider
provider "google-beta" {
}

resource "google_cloudbuild_trigger" "main" {
  provider = google-beta
  name = "self-main"
  project = var.primary_project
  description = "Builds the latest version of SELF from the main branch."
  github {
    owner = "FluidNumerics"
    name = "SELF"
    push {
      branch = "main"
    }
  }
  substitutions = {
    _IMAGE_TAG = "latest"
    _BUILD_TYPE = "release"
  }
  filename = "cloudbuild.yaml"
}

resource "google_cloudbuild_trigger" "dev" {
  provider = google-beta
  name = "self-dev"
  project = var.primary_project
  description = "Builds the latest version of SELF from the main branch."
  github {
    owner = "FluidNumerics"
    name = "SELF"
    push {
      branch = "develop"
    }
  }
  substitutions = {
    _IMAGE_TAG = "dev"
    _BUILD_TYPE = "release"
  }
  filename = "cloudbuild.yaml"
}
