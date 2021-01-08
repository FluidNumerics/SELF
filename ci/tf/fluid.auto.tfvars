cluster_name = "self"
controller_machine_type = "n1-standard-1"
slurm_gcp_admins = ["group:support@fluidnumerics.com"]
slurm_gcp_users = ["user:support@fluidnumerics.com"]
primary_project = "self-fluids"
primary_zone = "us-west1-b"

partitions = [{name = "k80x1"
               project = "self-fluids"
               max_time = "1:00:00"
               labels = {"slurm-gcp"="compute"}
               machines = [{ name = "k80x1-usw1b"
                             gpu_count = 1
                             gpu_type = "nvidia-tesla-k80"
                             machine_type = "n1-standard-8"
                             max_node_count = 5
                             zone = "us-west1-b"
                          }]
               }]

builds = [{ name = "self-dev"
            repo_owner = "FluidNumerics"
            repo_name = "SELF"
            repo_branch = "develop"
            image_tag = "dev"
            compute_partition = "k80x1"
            cloudbuild_yaml = "cloudbuild.yaml"
         }]
