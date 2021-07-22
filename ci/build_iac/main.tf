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
    _SELF_IMAGE = var.gce_image
    _SELF_GPU_TYPE = "nvidia-tesla-v100"
    _SELF_GPU_COUNT = "0"
    _SELF_GPU_TARGET = "none"
    _SELF_VPC_SUBNET = google_compute_subnetwork.self-ci.self_link
    _SELF_SERVICE_ACCOUNT = google_service_account.self_ci.email
    _SELF_SINGULARITY_GPUFLAG = ""
    _SELF_TAGS = "self-ci"
  }
  filename = "ci/cloudbuild.yaml"
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
    _SELF_IMAGE = var.gce_image
    _SELF_GPU_TYPE = "nvidia-tesla-v100"
    _SELF_GPU_TARGET = "sm_72"
    _SELF_GPU_COUNT = "1"
    _SELF_VPC_SUBNET = google_compute_subnetwork.self-ci.self_link
    _SELF_SERVICE_ACCOUNT = google_service_account.self_ci.email
    _SELF_SINGULARITY_GPUFLAG = "--nv"
    _SELF_TAGS = "self-ci"
  }
  filename = "ci/cloudbuild.yaml"
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
    _BUILD_TYPE = "dev"
    _PLATFORM = "serial-x86"
    _SELF_NAME = "self"
    _SELF_ZONE = "us-west1-b"
    _SELF_MACHINE_TYPE = "n1-standard-2"
    _SELF_NODE_COUNT = "1"
    _SELF_IMAGE = var.gce_image
    _SELF_GPU_TYPE = "nvidia-tesla-v100"
    _SELF_GPU_COUNT = "0"
    _SELF_GPU_TARGET = "none"
    _SELF_VPC_SUBNET = google_compute_subnetwork.self-ci.self_link
    _SELF_SERVICE_ACCOUNT = google_service_account.self_ci.email
    _SELF_SINGULARITY_GPUFLAG = ""
    _SELF_TAGS = "self-ci"
  }
  filename = "ci/cloudbuild.yaml"
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
    _BUILD_TYPE = "dev"
    _PLATFORM = "serial-x86-nvcc"
    _SELF_NAME = "self"
    _SELF_ZONE = "us-west1-b"
    _SELF_MACHINE_TYPE = "n1-standard-8"
    _SELF_NODE_COUNT = "1"
    _SELF_IMAGE = var.gce_image
    _SELF_GPU_TYPE = "nvidia-tesla-v100"
    _SELF_GPU_COUNT = "1"
    _SELF_GPU_TARGET = "sm_72"
    _SELF_VPC_SUBNET = google_compute_subnetwork.self-ci.self_link
    _SELF_SERVICE_ACCOUNT = google_service_account.self_ci.email
    _SELF_SINGULARITY_GPUFLAG = "--nv"
    _SELF_TAGS = "self-ci"
  }
  filename = "ci/cloudbuild.yaml"
}

// Big Query Dataset for SELF CI Data
resource "google_bigquery_dataset" "self_ci" {
  dataset_id = "self_ci"
  friendly_name = "SELF CI/CB data"
  description = "A dataset containing build information for SELF from our CI/CB pipeline."
  location = var.bq_location
  project = var.project
}

resource "google_bigquery_table" "benchmarks" {
  dataset_id = google_bigquery_dataset.self_ci.dataset_id
  table_id = "cloud_build_data"
  project = var.project
  deletion_protection=false
  schema = <<EOF
[
  {
    "name": "cli_command",
    "type": "STRING",
    "mode": "REQUIRED",
    "description": "The CLI command being executed. Usually maps to a subroutine being exercised during testing."
  },
  {
    "name": "execution_command",
    "type": "STRING",
    "mode": "REQUIRED",
    "description": "The full command used to execute this benchmark"
  },
  {
    "name": "build_id",
    "type": "STRING",
    "mode": "REQUIRED",
    "description": "The Cloud Build build ID associated with this build."
  },
  {
    "name": "build_type",
    "type": "STRING",
    "mode": "REQUIRED",
    "description": "The build type passed to SELF make system during build. (_BUILD_TYPE Cloud Build Variable)"
  },
  {
    "name": "compiler",
    "type": "STRING",
    "mode": "NULLABLE",
    "description": "The compiler name and version (e.g. gcc@10.2.0) used to build SELF"
  },
  {
    "name": "container_platform",
    "type": "STRING",
    "mode": "REQUIRED",
    "description": "Name of the container platform used to build the application (e.g. docker, singularity). If not containerized, set to `none`"
  },
  {
    "name": "container_platform_runtime",
    "type": "STRING",
    "mode": "NULLABLE",
    "description": "Name of the container platform used to run the application (e.g. docker, singularity). If not containerized, set to `none`"
  },
  {
    "name": "mpi_provider",
    "type": "STRING",
    "mode": "NULLABLE",
    "description": "The MPI provider and version used to build the application (e.g. openmpi-4.0.2). If the application is not built with MPI, set to `none`."
  },
  {
    "name": "target_platform",
    "type": "STRING",
    "mode": "REQUIRED",
    "description": "Name of the target platform used in the build (_PLATFORM variable in Cloud Build)."
  },
  {
    "name": "machine_type",
    "type": "STRING",
    "mode": "NULLABLE",
    "description": "Node types as classified by the system provider."
  },
  {
    "name": "gpu_type",
    "type": "STRING",
    "mode": "NULLABLE",
    "description": "The vendor and model name of the GPU (e.g. nvidia-tesla-v100)"
  },
  {
    "name": "gpu_count",
    "type": "INT64",
    "mode": "NULLABLE",
    "description": "The number of GPUs, per compute node, on this compute system."
  },
  {
    "name": "node_count",
    "type": "INT64",
    "mode": "NULLABLE",
    "description": "The number of nodes used in testing SELF."
  },
  {
    "name": "datetime",
    "type": "DATETIME",
    "mode": "REQUIRED",
    "description" : "The UTC date and time of the build."
  },
  {
    "name": "exit_code",
    "type": "INT64",
    "mode": "REQUIRED",
    "description": "The system exit code thrown when executing benchmark_info.cli_command"
  },
  {
    "name": "git_sha",
    "type": "STRING",
    "mode": "REQUIRED",
    "description": "The git SHA associated with the version / commit being tested." 
  },
  {
    "name": "stderr",
    "type": "STRING",
    "mode": "NULLABLE",
    "description": "Standard error." 
  },
  {
    "name": "stdout",
    "type": "STRING",
    "mode": "NULLABLE",
    "description": "Standard output." 
  }
]
EOF
}
