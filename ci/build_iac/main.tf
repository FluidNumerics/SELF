terraform {
  backend "gcs" {
    bucket  = "self-fluids-terraform"
    prefix  = "fluid-cicb"
  }
}

provider "google" {
}

locals {
  region = trimsuffix(var.zone,substr(var.zone,-2,-2))
}

// Service account for CI tests
resource "google_service_account" "fluid_cicb" {
  account_id = "fluid-cicb"
  display_name = "Continuous Integration Service account"
  project = var.project
}

// Service Accounts //
resource "google_service_account" "slurm_controller" {
  account_id = "slurm-gcp-controller"
  display_name = "Slurm-GCP Controller Service Account"
  project = var.project
}

resource "google_service_account" "slurm_compute" {
  account_id = "slurm-gcp-compute"
  display_name = "Slurm-GCP Compute Service Account"
  project = var.project
}

// **** IAM Permissions **** ///
resource "google_project_iam_member" "project_compute_image_users" {
  project = var.project
  role = "roles/compute.imageUser"
  member = "serviceAccount:${google_service_account.slurm_controller.email}"
}

resource "google_project_iam_member" "project_compute_admins" {
  project = var.project
  role = "roles/compute.admin"
  member = "serviceAccount:${google_service_account.slurm_controller.email}"
}

resource "google_project_iam_member" "service_account_user" {
  project = var.project
  role = "roles/iam.serviceAccountUser"
  member = "serviceAccount:${google_service_account.slurm_controller.email}"
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
  _ZONE = var.zone
  _REGION = local.region
  _SLURM_CONTROLLER = module.slurm_cluster_controller.controller_node_name
  _BUILD_TYPE = var.builds[count.index].build_type
  _PLATFORM = var.builds[count.index].platform
  _PARTITIONS = var.builds[count.index].partitions
  _GPU_ACCEL = var.builds[count.index].gpu_accel
  }
  filename = var.cloudbuild_path
}

// Big Query Dataset for CICB Data
resource "google_bigquery_dataset" "fluid_cicb" {
  dataset_id = "fluid_cicb"
  friendly_name = "Fluid CI/CB data"
  description = "A dataset containing build information for the Fluid CI/CB pipeline."
  location = var.bq_location
  project = var.project
}

resource "google_bigquery_table" "benchmarks" {
  dataset_id = google_bigquery_dataset.fluid_cicb.dataset_id
  table_id = "app_runs"
  project = var.project
  deletion_protection=false
  schema = <<EOF
[
  {
    "name": "allocated_cpus",
    "type": "INT64",
    "mode": "NULLABLE",
    "description": "The number of CPUs that are allocated to run the execution_command."
  },
  {
    "name": "allocated_gpus",
    "type": "INT64",
    "mode": "NULLABLE",
    "description": "The number of GPUs that are allocated to run the execution_command."
  },
  {
    "name": "command_group",
    "type": "STRING",
    "mode": "REQUIRED",
    "description": "An identifier to allow grouping of execution_commands in reporting. This is particularly useful if you are exercising multiple options for the same CLI command and want to be able to group results and profile metrics for multiple execution commands."
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
    "description": "The number of nodes used in testing."
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
    "description": "The system exit code thrown when executing the execution_command"
  },
  {
    "name": "git_sha",
    "type": "STRING",
    "mode": "REQUIRED",
    "description": "The git SHA associated with the version / commit being tested." 
  },
  {
    "name": "max_memory_gb",
    "type": "FLOAT64",
    "mode": "NULLABLE",
    "description": "The maximum amount of memory used for the execution_command in GB."
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
  },
  {
    "name": "partition",
    "type": "STRING",
    "mode": "NULLABLE",
    "description": "(Optional) The name of the scheduler partition to run the job under. If provided, the execution_command is interpreted as the path to a batch script." 
  },
  {
    "name": "runtime",
    "type": "FLOAT64",
    "mode": "NULLABLE",
    "description": "The runtime for the execution_command in seconds."
  }
]
EOF
}

module "slurm_cluster_network" {
  source = "github.com/FluidNumerics/slurm-gcp//tf/modules/network"
  cluster_name                  = "fluid-cicb"
  disable_login_public_ips      = false
  disable_controller_public_ips = var.disable_controller_public_ips
  disable_compute_public_ips    = true
  network_name                  = null
  partitions                    = var.partitions
  shared_vpc_host_project       = null
  subnetwork_name               = null
  project = var.project
  region  = local.region
}

module "slurm_cluster_controller" {
  source = "github.com/FluidNumerics/slurm-gcp//tf/modules/controller"
  boot_disk_size                = var.controller_disk_size_gb
  boot_disk_type                = var.controller_disk_type
  image                         = var.controller_image
  instance_template             = var.controller_instance_template
  cluster_name                  = "fluid-cicb"
  compute_node_scopes           = ["https://www.googleapis.com/auth/cloud-platform"]
  compute_node_service_account  = google_service_account.slurm_compute.email
  disable_compute_public_ips    = true
  disable_controller_public_ips = var.disable_controller_public_ips
  labels                        = var.controller_labels
  login_network_storage         = []
  login_node_count              = 0
  machine_type                  = var.controller_machine_type
  munge_key                     = null
  jwt_key                       = null
  network_storage               = []
  partitions                    = var.partitions
  project                       = var.project
  region                        = local.region
  secondary_disk                = var.controller_secondary_disk
  secondary_disk_size           = var.controller_secondary_disk_size
  secondary_disk_type           = var.controller_secondary_disk_type
  shared_vpc_host_project       = var.project
  scopes                        = ["https://www.googleapis.com/auth/cloud-platform"]
  service_account               = google_service_account.slurm_controller.email
  subnet_depend                 = module.slurm_cluster_network.subnet_depend
  subnetwork_name               = null
  suspend_time                  = var.suspend_time
  zone                          = var.zone
}
