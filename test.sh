#!/bin/sh
#
#  NAME
#    test.sh
#
#  DESCRIPTION
#    Execute the SELF CI tests for a Docker image, Singularity image, or GCE VM image
#    Tests are expected to be run on a slurm cluster under srun, e.g
#
#      srun -n1 --partition=PARTITION --gres=gpu:1 test.sh --artifact singularity --gpu yes
#
#  USAGE
#
#    test.sh [--artifact [{singularity|docker|gce|none}]] [--gpu [{yes|no}]]
#    [--image IMAGE]
#
#    By default, if invoked without any flags, test.sh will attempt to run 
#    the serial cpu-only tests without any artifacts (assuming SELF is installed
#    under /opt/self.
#


# START CLI Interface
while [ "$#" -ge "1" ]; do

  key="$1"
  case $key in
    --artifact)
      case "$2" in
        docker|Docker|DOCKER)
	  ARTIFACT='docker'
	  shift
	  ;;
	singularity|Singularity|SINGULARITY)
	  ARTIFACT='singularity'
	  shift
	  ;;
	gce|GCE)
	  ARTIFACT='gce'
	  shift
	  ;;
	none|None|NONE)
	  ARTIFACT='none'
	  shift
	  ;;
      esac
      ;;
    --gpu)
      case $2 in 
        yes|Yes|YES)
          CI_SH="ci.gpu.sh"
          GPU="yes"
          shift
          ;;
        no|No|NO)
          CI_SH="ci.sh"
          GPU="no"
          shift
          ;;
      esac
      ;;
    --image)
      IMAGE="$2"
      shift
      ;;
  esac
  shift
done

# END CLI Interface
# /////////////////////////////// #

set -x

if [ "$ARTIFACT" = 'singularity' ]; then

  source /apps/spack/share/spack/setup-env.sh
  spack load singularity

  if [ "$GPU" = 'yes' ];then
    singularity exec --nv $IMAGE ci.gpu.sh
  else
    singularity exec $IMAGE ci.sh
  fi

fi
