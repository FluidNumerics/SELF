// Configure the Google Cloud provider
provider "google" {
}

resource "google_service_account" "gce_service_account" {
  account_id = "self-service"
  display_name = "SELF Service account"
  project = var.project
}

locals {
  region = trimsuffix(var.zone,substr(var.zone,-2,-2))
}


// **** Create the Shared VPC Network **** //
resource "google_compute_network" "vpc_network" {
  count = var.vpc_subnet == "" ? 1:0
  name = "${var.name_prefix}-network"
  project = var.project
  auto_create_subnetworks = false
}

resource "google_compute_subnetwork" "subnet" {
  count = var.vpc_subnet == "" ? 1:0
  name = "${var.name_prefix}-subnet"
  ip_cidr_range = var.subnet_cidr
  region = local.region
  network = google_compute_network.vpc_network[0].self_link
  project = var.project
}

resource "google_compute_firewall" "default_ssh_firewall_rules" {
  count = var.vpc_subnet == "" ? 1:0
  name = "${var.name_prefix}-ssh"
  network = google_compute_network.vpc_network[0].self_link
  target_tags = [var.name_prefix]
  source_ranges = var.whitelist_ssh_ips
  project = var.project

  allow {
    protocol = "tcp"
    ports = ["22"]
  }
}

resource "google_compute_firewall" "default_internal_firewall_rules" {
  count = var.vpc_subnet == "" ? 1:0
  name = "${var.name_prefix}-all-internal"
  network = google_compute_network.vpc_network[0].self_link
  source_tags = [var.name_prefix]
  target_tags = [var.name_prefix]
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


locals {
  subnet = var.vpc_subnet == "" ? google_compute_subnetwork.subnet[0].self_link:var.vpc_subnet
  tags = var.vpc_subnet == "" ? [var.name_prefix]: var.tags
}

// ***************************************** //
// Create the controller
resource "google_compute_instance" "gce_nodes" {
  count = var.node_count
  name = "${var.name_prefix}-${count.index}"
  project = var.project
  machine_type = var.machine_type
  zone = var.zone
  tags = local.tags
  boot_disk {
    auto_delete = true
    initialize_params {
      image = var.image
      size = var.disk_size_gb
      type = var.disk_type
    }
  }
  labels = var.labels 
  metadata_startup_script = file("${path.module}/startup-script.sh")
  metadata = {
    node_count = var.node_count
    name_prefix = var.name_prefix
    enable-oslogin = "TRUE"
  }
  network_interface {

    dynamic "access_config" {
      for_each = count.index == 0 ? [1] : []
        content {}
    }

    subnetwork  = local.subnet
  }

  dynamic "guest_accelerator" {
    for_each = var.gpu_count == 0 ? [] : [1]
    content{
      type = var.gpu_type
      count = var.gpu_count
    }
  }

  scheduling {
    on_host_maintenance = "TERMINATE"
  }


  service_account {
    email  = google_service_account.gce_service_account.email
    scopes = ["cloud-platform"]
  }
  lifecycle{
    ignore_changes = [metadata_startup_script]
  }
  allow_stopping_for_update = true
}
