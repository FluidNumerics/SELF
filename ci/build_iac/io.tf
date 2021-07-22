
variable "bq_location" {
  type = string
  description = "Valid location for Big Query Dataset. https://cloud.google.com/bigquery/docs/locations"
  default = "US"
}

variable "gce_image" {
  type = string
  description = "Google Compute Engine instance VM image for the test GCE cluster"
  default = "projects/hpc-apps/global/images/family/fluid-hpc-singularity-gcc-10-ompi-4"
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

variable "whitelist_ssh_ips" {
  type = list(string)
  description = "IP addresses that should be added to a whitelist for ssh access"
  default = ["0.0.0.0/0"]
}

variable "zone" {
  type = string
  description = "GCP Zone to deploy your cluster cluster. Learn more at https://cloud.google.com/compute/docs/regions-zones"
}
