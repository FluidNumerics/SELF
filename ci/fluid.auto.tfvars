cluster_name = "<name>"
project = "<project>"
zone = "<zone>"
shared_vpc_host_project = "<project>"
subnetwork_name = "fluid-run"

suspend_time = 2

controller_image = "<image>"
disable_controller_public_ips = false
controller_machine_type = "n1-standard-8"
controller_disk_size_gb = 1024
controller_disk_type = "pd-ssd"

login_image = "<image>"
disable_login_public_ips = true
login_machine_type = "n1-standard-4"
login_node_count = 0


compute_node_scopes          = [
  "https://www.googleapis.com/auth/cloud-platform"
]
partitions = [
  { name                 = "c2-standard-4"
    machine_type         = "c2-standard-4"
    image                = "<image>"
    image_hyperthreads   = true
    static_node_count    = 0
    max_node_count       = 25
    zone                 = "<zone>"
    compute_disk_type    = "pd-standard"
    compute_disk_size_gb = 100
    compute_labels       = {}
    cpu_platform         = null
    gpu_count            = 0
    gpu_type             = null
    gvnic                = false
    network_storage      = []
    preemptible_bursting = false
    vpc_subnet           = "fluid-run"
    exclusive            = false
    enable_placement     = false
    regional_capacity    = false
    regional_policy      = null
    instance_template    = null
  }
]

# ** Uncomment to use CloudSQL as Slurm database ** #
#cloudsql_slurmdb = true
#cloudsql_enable_ipv4 = false
#cloudsql_name = slurmdb
#cloudsql_tier = "db-n1-standard-8"

create_filestore = false
filestore = { name = "filestore"
              zone = null
              tier = "PREMIUM"
              capacity_gb = 2048
              fs_name = "nfs"
              network = null
            }

create_lustre = false
lustre = { image = "projects/research-computing-cloud/global/images/family/rcc-lustre-centos-7"
           project = null
           zone = null
           vpc_subnet = null
           service_account = null
           network_tags = []
           name = "lustre-gcp"
           fs_name = "lustre"
           mds_node_count = 1
           mds_machine_type = "n2-standard-16"
           mds_boot_disk_type = "pd-standard"
           mds_boot_disk_size_gb = 100
           mdt_disk_type = "pd-ssd"
           mdt_disk_size_gb = 1024
           mdt_per_mds = 1
           oss_node_count = 2
           oss_machine_type = "n2-standard-16" 
           oss_boot_disk_type = "pd-standard"
           oss_boot_disk_size_gb = 100
           ost_disk_type = "local-ssd"
           ost_disk_size_gb = 1500 
           ost_per_oss = 1
           hsm_node_count = 0
           hsm_machine_type = "n2-standard-16"
           hsm_gcs_bucket = null
           hsm_gcs_prefix = null
         }


