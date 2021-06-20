
variable "project" {
  type = string
  description = "GCP Project ID"
}

variable "machine_type" {
  type = string
  description = "GCE instance types for your scheduler-less cluster. Learn more at https://cloud.google.com/compute/docs/machine-types."
  default = "n1-standard-8"
}

variable "node_count" {
  type = number
  description = "Number of VM instances in your cluster. "
  default = 1
}

variable "zone" {
  type = string
  description = "GCP Zone to deploy your cluster cluster. Learn more at https://cloud.google.com/compute/docs/regions-zones"
}

variable "tags" {
  type = list(string)
  description = "GCE Network tags to apply to the instance. Learn more at https://cloud.google.com/vpc/docs/add-remove-network-tags."
  default = []
}

variable "image" {
  type = string
  description = "VM image used to launch your HPC application"
  default = "projects/hpc-apps/global/images/gromacs-gcp-latest"
}

variable "disk_size_gb" {
  type = number
  description = "Size of each GCE instance boot disk in your cluster."
  default = 20
}

variable "disk_type" {
  type = string
  description = "Disk type for each GCE instance boot disk in your cluster. Learn more at https://cloud.google.com/compute/docs/disks"
  default = "pd-standard"
}

variable "labels" {
  type = map
  description = "Resource labels to apply to each GCE instance in your cluster."
  default = {"hpc-app"="gromacs"}
}

variable "name_prefix" {
  type = string
  description = "The name to prefix all GCE instances in your cluster"
  default = "gromacs"
}

variable "gpu_type" {
  type = string
  description = "Type of GPU Accelerator to attach to your GCE instance. Learn more at https://cloud.google.com/compute/docs/gpus"
  default = "nvidia-tesla-v100"
}

variable "gpu_count" {
  type = number
  description = "Number of GPU Accelerators to attach to your GCE instance."
  default = 0
}

variable "vpc_subnet" {
  type = string
  description = "VPC subnet (self-link) to use for your cluster."
  default = ""
}

variable "subnet_cidr" {
  type = string
  description = "CIDR Range for the subnet (if vpc_subnet is not provided)."
  default = "10.10.0.0/16"
}

variable "whitelist_ssh_ips" {
  type = list(string)
  description = "IP addresses that should be added to a whitelist for ssh access"
  default = ["0.0.0.0/0"]
}

output "head_node_ip_addr" {
  value = google_compute_instance.gce_nodes[0].network_interface[0].access_config[0].nat_ip
}
