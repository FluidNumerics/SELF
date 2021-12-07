project = "self-fluids"
zone = "us-west1-b"
github_owner = "FluidNumerics"
github_repo = "SELF"

builds = [{name="main-v100",
           description="Latest main build targeting V100 GPUs",
           branch="main",
           gpu_target="sm_72",
           partitions="c2-standard-4",
           zone="us-west1-b"},
          {name="develop-v100",
           description="Latest develop buil targeting V100 GPUs",
           branch="develop",
           gpu_target="sm_72",
           zone="us-west1-b",
           partitions="c2-standard-4"}
]
