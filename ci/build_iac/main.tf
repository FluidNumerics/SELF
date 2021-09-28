terraform {
  backend "gcs" {
    bucket  = "self-fluids-terraform"
    prefix  = "fluid-cicb"
  }
}

provider "google" {
}

locals {
}

resource "google_cloudbuild_trigger" "builds" {
  count = length(var.builds)
  name = var.builds[count.index].name
  project = var.project
  description = var.builds[count.index].description 
  github {
    owner = var.github_owner
    name = var.github_repo
    push {
      branch = var.builds[count.index].branch
    }
  }
  substitutions = {
  _ZONE = var.builds[count.index].zone
  _PARTITIONS = var.builds[count.index].partitions
  _GPU_TARGET = var.builds[count.index].gpu_target
  }
  filename = "ci/cloudbuild.yaml"
}

module "fluid_cicb" {
  source = "github.com/FluidNumerics/fluid-run//tf/"
  bq_location = var.bq_location
  project = var.project
  subnet_cidr = var.subnet_cidr
  whitelist_ssh_ips = var.whitelist_ssh_ips
}
