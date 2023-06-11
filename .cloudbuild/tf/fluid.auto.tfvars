project = "self-fluids"
zone = "us-west1-b"
github_owner = "FluidNumerics"
github_repo = "SELF"

branch_builds = [{name="main-branch-v100-double",
           description="Double precision build targeting Nvidia V100 GPU for main branch push",
           branch="main",
           gpu_target="sm_72",
           hip_platform="nvidia",
           prec="double"
           zone="us-west1-b"},
          {name="main-mi50-double",
           description="Double precision build targeting AMD MI50 GPU for main branch push",
           branch="main",
           gpu_target="gfx906",
           hip_platform="amd",
           prec="double"
           zone="us-west1-b"}
          ]

builds = [{name="main-v100-double",
           description="Double precision build targeting Nvidia V100 GPU for PR's to main branch",
           branch="main",
           gpu_target="sm_72",
           hip_platform="nvidia",
           prec="double"
           zone="us-west1-b"},
          {name="main-mi50-double",
           description="Double precision build targeting AMD MI50 GPU for PR's to main branch",
           branch="main",
           gpu_target="gfx906",
           hip_platform="amd",
           prec="double"
           zone="us-west1-b"}
          ]
