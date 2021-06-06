// Configure the Google Cloud provider
provider "google" {
}

resource "google_service_account" "gce_service_account" {
  account_id = "gromacs"
  display_name = "GROMACS Service account"
  project = var.project
}


// ***************************************** //
// Create the controller
resource "google_compute_instance" "gce_nodes" {
  count = var.node_count
  name = "${var.name_prefix}-${count.index}"
  project = var.project
  machine_type = var.machine_type
  zone = var.zone
  tags = var.tags
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

    subnetwork  = var.vpc_subnet
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
