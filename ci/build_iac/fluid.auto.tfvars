project = "self-fluids"
zone = "us-west1-b"
github_owner = "FluidNumerics"
github_repo = "SELF"

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
          {name="dev-v100-double",
           description="Double precision build targeting Nvidia V100 GPU for PR's to develop branch",
           branch="develop",
           gpu_target="sm_72",
           hip_platform="nvidia",
           prec="double"
           zone="us-west1-b"},
          {name="develop-v100-single",
           description="Single precision build targeting Nvidia V100 GPU for PR's to develop branch",
           branch="main",
           gpu_target="sm_72",
           hip_platform="nvidia",
           prec="single"
           zone="us-west1-b"}

]
