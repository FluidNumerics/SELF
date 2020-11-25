terraform {
  backend "gcs" {
    bucket  = "self-fluids-terraform"
    prefix  = "self-build-triggers"
  }
}

// Configure the Google Cloud provider
provider "google" {
 version = "3.9"
}

resource "google_cloudbuild_trigger" "dev-builds" {
  name = "self-dev"
  github {
    owner = "FluidNumerics"
    name = "SELF"
    push {
      branch = "develop"
    }
  }
  substitutions {
    _IMAGE_TAG = "dev"
  }
  filename = "cloudbuild.yaml"
}

resource "google_cloudbuild_trigger" "main-builds" {
  name = "self-main"
  github {
    owner = "FluidNumerics"
    name = "SELF"
    push {
      branch = "main"
    }
  }
  substitutions {
    _IMAGE_TAG = "latest"
  }
  filename = "cloudbuild.yaml"
}
