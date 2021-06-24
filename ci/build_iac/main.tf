terraform {
  backend "gcs" {
    bucket  = "self-fluids-terraform"
    prefix  = "self-ci"
  }
}

provider "google" {
}

locals {
  region = trimsuffix(var.zone,substr(var.zone,-2,-2))
}

// Service account for CI tests
resource "google_service_account" "self_ci" {
  account_id = "self-cibot"
  display_name = "SELF Continuous Integration Service account"
  project = var.project
}

// **** Create the Shared VPC Network **** //
resource "google_compute_network" "vpc_network" {
  name = "self-ci"
  project = var.project
  auto_create_subnetworks = false
}

resource "google_compute_subnetwork" "self-ci" {
  name = "self-ci"
  ip_cidr_range = var.subnet_cidr
  region = local.region
  network = google_compute_network.vpc_network.self_link
  project = var.project
}

resource "google_compute_firewall" "default_ssh_firewall_rules" {
  name = "self-ci-ssh"
  network = google_compute_network.vpc_network.self_link
  target_tags = ["self-ci"]
  source_ranges = var.whitelist_ssh_ips
  project = var.project

  allow {
    protocol = "tcp"
    ports = ["22"]
  }
}

resource "google_compute_firewall" "default_internal_firewall_rules" {
  name = "self-ci-all-internal"
  network = google_compute_network.vpc_network.self_link
  source_tags = ["self-ci"]
  target_tags = ["self-ci"]
  project = var.project

  allow {
    protocol = "tcp"
    ports = ["0-65535"]
  }
  allow {
    protocol = "udp"
    ports = ["0-65535"]
  }
  allow {
    protocol = "icmp"
    ports = []
  }
}

resource "google_cloudbuild_trigger" "main-serial" {
  name = "self-main-serial-x86"
  project = var.project
  description = "Builds the latest version of SELF from the main branch."
  github {
    owner = "HigherOrderMethods"
    name = "SELF"
    push {
      branch = "main"
    }
  }
  substitutions = {
    _BUILD_TYPE = "release"
    _PLATFORM = "serial-x86"
    _SELF_NAME = "self"
    _SELF_ZONE = "us-west1-b"
    _SELF_MACHINE_TYPE = "n1-standard-2"
    _SELF_NODE_COUNT = "1"
    _SELF_IMAGE = "projects/hpc-apps/global/images/singularity-gcp-dev"
    _SELF_GPU_TYPE = "nvidia-tesla-v100"
    _SELF_GPU_COUNT = "0"
    _SELF_VPC_SUBNET = "${google_compute_subnetwork.self-ci.self_link}"
    _SELF_SERVICE_ACCOUNT = "${google_service_account.self_ci.email}"
    _SELF_TAGS = "self-ci"
  }
  filename = "ci/cloud-build/gce/cloudbuild.yaml"
}

resource "google_cloudbuild_trigger" "main-serial-x86-nvcc" {
  name = "self-main-serial-x86-nvcc"
  project = var.project
  description = "Builds the latest version of SELF from the main branch."
  github {
    owner = "HigherOrderMethods"
    name = "SELF"
    push {
      branch = "main"
    }
  }
  substitutions = {
    _BUILD_TYPE = "release"
    _PLATFORM = "serial-x86-nvcc"
    _SELF_NAME = "self"
    _SELF_ZONE = "us-west1-b"
    _SELF_MACHINE_TYPE = "n1-standard-8"
    _SELF_NODE_COUNT = "1"
    _SELF_IMAGE = "projects/hpc-apps/global/images/singularity-gcp-dev"
    _SELF_GPU_TYPE = "nvidia-tesla-v100"
    _SELF_GPU_COUNT = "1"
    _SELF_VPC_SUBNET = "${google_compute_subnetwork.self-ci.self_link}"
    _SELF_SERVICE_ACCOUNT = "${google_service_account.self_ci.email}"
    _SELF_TAGS = "self-ci"
  }
  filename = "ci/cloud-build/gce/cloudbuild.yaml"
}

resource "google_cloudbuild_trigger" "develop-serial" {
  name = "self-develop-serial-x86"
  project = var.project
  description = "Builds the latest version of SELF from the develop branch."
  github {
    owner = "HigherOrderMethods"
    name = "SELF"
    push {
      branch = "develop"
    }
  }
  substitutions = {
    _BUILD_TYPE = "release"
    _PLATFORM = "serial-x86"
    _SELF_NAME = "self"
    _SELF_ZONE = "us-west1-b"
    _SELF_MACHINE_TYPE = "n1-standard-2"
    _SELF_NODE_COUNT = "1"
    _SELF_IMAGE = "projects/hpc-apps/global/images/singularity-gcp-dev"
    _SELF_GPU_TYPE = "nvidia-tesla-v100"
    _SELF_GPU_COUNT = "0"
    _SELF_VPC_SUBNET = "${google_compute_subnetwork.self-ci.self_link}"
    _SELF_SERVICE_ACCOUNT = "${google_service_account.self_ci.email}"
    _SELF_TAGS = "self-ci"
  }
  filename = "ci/cloud-build/gce/cloudbuild.yaml"
}

resource "google_cloudbuild_trigger" "develop-serial-x86-nvcc" {
  name = "self-develop-serial-x86-nvcc"
  project = var.project
  description = "Builds the latest version of SELF from the develop branch."
  github {
    owner = "HigherOrderMethods"
    name = "SELF"
    push {
      branch = "develop"
    }
  }
  substitutions = {
    _BUILD_TYPE = "release"
    _PLATFORM = "serial-x86-nvcc"
    _SELF_NAME = "self"
    _SELF_ZONE = "us-west1-b"
    _SELF_MACHINE_TYPE = "n1-standard-8"
    _SELF_NODE_COUNT = "1"
    _SELF_IMAGE = "projects/hpc-apps/global/images/singularity-gcp-dev"
    _SELF_GPU_TYPE = "nvidia-tesla-v100"
    _SELF_GPU_COUNT = "1"
    _SELF_VPC_SUBNET = "${google_compute_subnetwork.self-ci.self_link}"
    _SELF_SERVICE_ACCOUNT = "${google_service_account.self_ci.email}"
    _SELF_TAGS = "self-ci"
  }
  filename = "ci/cloud-build/gce/cloudbuild.yaml"
}
