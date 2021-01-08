terraform {
  backend "gcs" {
    bucket  = "self-fluids-terraform"
    prefix  = "self-build-triggers"
  }
}

// Configure the Google Cloud provider
provider "google-beta" {
}

resource "google_cloudbuild_trigger" "main-builds" {
  provider = google-beta
  name = "self-main"
  description = "Builds the latest version of SELF from the main branch"
  github {
    owner = "FluidNumerics"
    name = "SELF"
    push {
      branch = "main"
    }
  }
  substitutions = {
    _IMAGE_TAG = "latest"
  }
  filename = "cloudbuild.yaml"
}
