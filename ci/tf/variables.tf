variable "cluster_name" {
  type = string
  description = "Customer organization ID from the managed-fluid-slurm-gcp customers database"
}

variable "subnet_cidr" {
  type = string
  description = "CIDR Range for cluster VPC Subnet."
  default = "10.10.0.0/16"
}

variable "slurm_gcp_admins" {
  type = list(string)
  description = "A list of users that will serve as Linux System Administrators on your cluster. Set each element to 'user:someone@example.com' for users or 'group:somegroup@example.com' for groups"
}

variable "slurm_gcp_users" {
  type = list(string)
  description = "A list of users that will serve as Linux System Administrators on your cluster. Set each element to 'user:someone@example.com' for users or 'group:somegroup@example.com' for groups"
}

variable "controller_image" {
  type = string
  description = "Image to use for the fluid-slurm-gcp controller"
  default = "projects/fluid-cluster-ops/global/images/fluid-slurm-gcp-controller-centos"
}

variable "compute_image" {
  type = string
  description = "Image to use for the fluid-slurm-gcp compute instances (all partitions[].machines[])."
  default = "projects/fluid-cluster-ops/global/images/fluid-slurm-gcp-compute-centos"
}

variable "primary_project" {
  type = string
  description = "Main GCP project ID for the customer's managed solution"
}

variable "primary_zone" {
  type = string
  description = "Main GCP zone for the customer's managed solution"
}

variable "whitelist_ssh_ips" {
  type = list(string)
  description = "IP addresses that should be added to a whitelist for ssh access"
  default = ["0.0.0.0/0"]
}

variable "controller_machine_type" { 
  type = string
  description = "GCP Machine type to use for the controller node."
}

variable "default_partition" {
  type = string
  description = "Name of the default compute partition."
  default = ""
}

variable "builds" {
  type = list(object({
      name = string
      repo_owner = string
      repo_name = string
      repo_branch = string
      image_tag = string
      compute_partition = string
      cloudbuild_yaml = string
  }))
  description = "Settings for builds and pairing with compute partitions."
  default = []
}

variable "slurm_accounts" {
  type = list(object({
      name = string
      users = list(string)
      allowed_partitions = list(string)
  }))
  default = []
}

variable "munge_key" {
  type = string
  default = ""
}

variable "suspend_time" {
  type = number
  default = 300
}

variable "partitions" {
  type = list(object({
      name = string
      project = string
      max_time= string
      labels = map(string)
      machines = list(object({
        name = string
        gpu_count = number
        gpu_type = string
        machine_type=string
        max_node_count= number
        zone= string
      }))
  }))
  description = "Settings for partitions and compute instances available to the cluster."
  
  default = []
}
