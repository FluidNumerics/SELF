project = "self-fluids"
zone = "us-west1-b"
github_owner = "FluidNumerics"
github_repo = "SELF"
cloudbuild_path = "ci/cloudbuild.yaml"
disable_controller_public_ips = false

builds = [{name="main-serial-x86",
           description="Serial CPU only release branch builds for x86 platforms",
           branch="main",
           build_type="release",
           platform="serial-x86"},
          {name="main-serial-x86-nvcc",
           description="Serial CPU with single GPU accelerator release branch builds for x86+Nvidia platforms",
           branch="main",
           build_type="release",
           platform="serial-x86-nvcc"},
          {name="develop-serial-x86",
           description="Serial CPU only dev branch builds for x86 platforms",
           branch="develop",
           build_type="dev",
           platform="serial-x86"},
          {name="develop-serial-x86-nvcc",
           description="Serial CPU with single GPU accelerator dev branch builds for x86+Nvidia platforms",
           branch="develop",
           build_type="dev",
           platform="serial-x86-nvcc"},

]


partitions = [
  { name                 = "c2-standard-8"
    machine_type         = "c2-standard-8"
    static_node_count    = 0
    max_node_count       = 10
    zone                 = "us-west1-b"
    image                = "projects/research-computing-cloud/global/images/family/rcc-centos-7-v3"
    image_hyperthreads   = true
    compute_disk_type    = "pd-standard"
    compute_disk_size_gb = 50
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
    regional_policy      = {}
    instance_template    = null
  },
  { name                 = "n1-standard-4-v100"
    machine_type         = "n1-standard-4-v100"
    static_node_count    = 0
    max_node_count       = 10
    zone                 = "us-west1-b"
    image                = "projects/research-computing-cloud/global/images/family/rcc-centos-7-v3"
    image_hyperthreads   = true
    compute_disk_type    = "pd-standard"
    compute_disk_size_gb = 50
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
    regional_policy      = {}
    instance_template    = null
  }
]
