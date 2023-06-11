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
    pull_request {
      branch = var.builds[count.index].branch
      comment_control = "COMMENTS_ENABLED"
    }
  }
  substitutions = {
  _ZONE = var.builds[count.index].zone
  _GPU_TARGET = var.builds[count.index].gpu_target
  _HIP_PLATFORM = var.builds[count.index].hip_platform
  _PREC = var.builds[count.index].prec
  }
  filename = "ci/cloudbuild.pr.yaml"
}

resource "google_cloudbuild_trigger" "branch_builds" {
  count = length(var.branch_builds)
  name = var.branch_builds[count.index].name
  project = var.project
  description = var.branch_builds[count.index].description 
  github {
    owner = var.github_owner
    name = var.github_repo
    push {
      branch = var.branch_builds[count.index].branch
    }
  }
  substitutions = {
  _ZONE = var.branch_builds[count.index].zone
  _GPU_TARGET = var.branch_builds[count.index].gpu_target
  _HIP_PLATFORM = var.branch_builds[count.index].hip_platform
  _PREC = var.branch_builds[count.index].prec
  }
  filename = "ci/cloudbuild.main.yaml"
}

module "fluid_cicb" {
  source = "github.com/FluidNumerics/fluid-run//tf/"
  bq_location = var.bq_location
  project = var.project
  subnet_cidr = var.subnet_cidr
  whitelist_ssh_ips = var.whitelist_ssh_ips
}
