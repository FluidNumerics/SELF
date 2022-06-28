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
          {name="main-branch-v100-single",
           description="Single precision build targeting Nvidia V100 GPU for main branch push",
           branch="main",
           gpu_target="sm_72",
           hip_platform="nvidia",
           prec="single"
           zone="us-west1-b"},
          {name="main-branch-mi25-double",
           description="Double precision build targeting AMD MI25 GPU for main branch push",
           branch="main",
           gpu_target="gfx900",
           hip_platform="amd",
           prec="double"
           zone="us-west1-b"},
          {name="dev-branch-mi25-double",
           description="Double precision build targeting AMD MI25 GPU for dev branch push",
           branch="develop",
           gpu_target="gfx900",
           hip_platform="amd",
           prec="double"
           zone="us-west1-b"},
          ]

builds = [{name="main-v100-double",
           description="Double precision build targeting Nvidia V100 GPU for PR's to main branch",
           branch="main",
           gpu_target="sm_72",
           hip_platform="nvidia",
           prec="double"
           zone="us-west1-b"},
          {name="main-v100-single",
           description="Single precision build targeting Nvidia V100 GPU for PR's to main branch",
           branch="main",
           gpu_target="sm_72",
           hip_platform="nvidia",
           prec="single"
           zone="us-west1-b"},
          {name="main-mi25-double",
           description="Double precision build targeting AMD MI25 GPU for PR's to main branch",
           branch="main",
           gpu_target="gfx900",
           hip_platform="amd",
           prec="double"
           zone="us-west1-b"},
          {name="main-mi25-single",
           description="Single precision build targeting AMD MI25 GPU for PR's to main branch",
           branch="main",
           gpu_target="gfx900",
           hip_platform="amd",
           prec="single"
           zone="us-west1-b"}
          ]
