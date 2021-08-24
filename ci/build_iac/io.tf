
variable "builds" {
  type = list(object({
    name = string
    description = string
    branch = string
    build_type = string
    platform = string
    }))
  default = []
  description = "List of build triggers and their settings to configure"
}

variable "bq_location" {
  type = string
  description = "Valid location for Big Query Dataset. https://cloud.google.com/bigquery/docs/locations"
  default = "US"
}

variable "cloudbuild_path" {
  type = string
  description = "Path to the cloudbuild.yaml file in the Github repository"
  default = "ci/cloudbuild.yaml"
}

variable "controller_machine_type" {
  description = "Machine type to use for the controller instance"
  type        = string
  default     = "n1-standard-4"
}

variable "controller_disk_type" {
  description = "Disk type (pd-ssd or pd-standard) for controller."
  type        = string
  default     = "pd-standard"
}

variable "controller_image" {
  description = "Slurm image to use for the controller instance"
  type        = string
  default     = "projects/hpc-apps/global/images/family/fluid-cicb"
}

variable "controller_instance_template" {
  description = "Instance template to use to create controller instance"
  type        = string
  default     = null
}

variable "controller_disk_size_gb" {
  description = "Size of disk for the controller."
  type        = number
  default     = 50
}

variable "controller_labels" {
  description = "Labels to add to controller instance. List of key key, value pairs."
  type        = any
  default     = null
}

variable "controller_secondary_disk" {
  description = "Create secondary disk mounted to controller node"
  type        = bool
  default     = false
}

variable "controller_secondary_disk_size" {
  description = "Size of disk for the secondary disk"
  default     = 100
}

variable "controller_secondary_disk_type" {
  description = "Disk type (pd-ssd or pd-standard) for secondary disk"
  default     = "pd-ssd"
}

variable "controller_scopes" {
  description = "Scopes to apply to the controller"
  type        = list(string)
  default     = ["https://www.googleapis.com/auth/cloud-platform"]
}

variable "controller_service_account" {
  description = "Service Account for the controller"
  type        = string
  default     = null
}

variable "disable_controller_public_ips" {
  type    = bool
  default = true
}

variable "gce_image" {
  type = string
  description = "Google Compute Engine instance VM image for the test GCE cluster"
  default = "projects/hpc-apps/global/images/family/fluid-cicb"
}

variable "github_owner" {
  type = string
  description = "The owner of the Github repository"
}

variable "github_repo" {
  type = string
  description = "The name of the Github repository"
}

variable "partitions" {
  description = "An array of configurations for specifying multiple machine types residing in their own Slurm partitions."
  type = list(object({
    name                 = string,
    machine_type         = string,
    max_node_count       = number,
    zone                 = string,
    image                = string,
    image_hyperthreads   = bool,
    compute_disk_type    = string,
    compute_disk_size_gb = number,
    compute_labels       = any,
    cpu_platform         = string,
    gpu_type             = string,
    gpu_count            = number,
    gvnic                = bool,
    network_storage = list(object({
      server_ip    = string,
      remote_mount = string,
      local_mount  = string,
      fs_type      = string,
    mount_options = string })),
    preemptible_bursting = bool,
    vpc_subnet           = string,
    exclusive            = bool,
    enable_placement     = bool,
    regional_capacity    = bool,
    regional_policy      = any,
    instance_template    = string,
  static_node_count = number }))
}

variable "project" {
  type = string
  description = "GCP Project ID"
}

variable "subnet_cidr" {
  type = string
  description = "CIDR Range for the subnet (if vpc_subnet is not provided)."
  default = "10.10.0.0/16"
}

variable "suspend_time" {
  description = "Idle time (in sec) to wait before nodes go away"
  default     = 300
}

variable "whitelist_ssh_ips" {
  type = list(string)
  description = "IP addresses that should be added to a whitelist for ssh access"
  default = ["0.0.0.0/0"]
}

variable "zone" {
  type = string
  description = "GCP Zone to deploy your cluster cluster. Learn more at https://cloud.google.com/compute/docs/regions-zones"
}
