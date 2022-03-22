cluster_name = "<name>"
project = "<project>"
zone = "<zone>"

controller_image = "<image>"
disable_controller_public_ips = false
controller_machine_type = "n1-standard-8"
controller_disk_size_gb = 100
controller_disk_type = "pd-standard"

login_image = "<image>"
disable_login_public_ips = true
login_machine_type = "n1-standard-4"
login_node_count = 0

suspend_time = 2

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
    vpc_subnet           = null
    exclusive            = false
    enable_placement     = false
    regional_capacity    = false
    regional_policy      = null
    instance_template    = null
  },
  { name                 = "v100"
    machine_type         = "n1-standard-8"
    image                = "<image>"
    image_hyperthreads   = true
    static_node_count    = 0
    max_node_count       = 25
    zone                 = "<zone>"
    compute_disk_type    = "pd-standard"
    compute_disk_size_gb = 100
    compute_labels       = {}
    cpu_platform         = null
    gpu_count            = 1
    gpu_type             = "nvidia-tesla-v100"
    gvnic                = false
    network_storage      = []
    preemptible_bursting = false
    vpc_subnet           = null
    exclusive            = false
    enable_placement     = false
    regional_capacity    = false
    regional_policy      = null
    instance_template    = null
  }
]
