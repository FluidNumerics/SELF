#!/bin/bash
#
#
# 


help_msg () {
      cat << EOF
singularity_exec.sh
---------------
Description:
    A helper script for running SELF examples as benchmarks under fluid-run

Usage:
  singularity_exec.sh --example /path/to/self/binary --self_opts "[SELF OPTS]"

Options:
  --example          (Required) Path to the example binary that you want to run a benchmark for.
  --singularity_opts (Optional) Options to pass to singularity.
  --self_opts        (Optional) Command line options to use when running the benchmark.
  --help             Print this message.

Examples:

  singularity_exec.sh --example /path/to/self/binary --opts "[SELF OPTS]"

EOF
      exit 0

}


# /////////////////////////////////// #
# CLI Interface

while [ "$#" -ge "1" ]; do

  key="$1"
  case $key in

    --example)
      EXAMPLE="$2"
      shift
      ;;

    --singularity_opts)
      SINGULARITY_OPTS="$2"
      shift
      ;;

    --self_opts)
      SELF_OPTS="$2"
      shift
      ;;

    --help)
      help_msg
      ;;

  esac
  shift
done
# /////////////////////////////////// #
source /etc/profile.d/z10_spack_environment.sh

cd ${WORKSPACE}
mkdir -p ${WORKSPACE}/build

if [ -z ${EXAMPLE} ]; then
  echo "No example file set. Exiting"
  exit 1
fi

CMD="singularity exec $SINGULARITY_OPTS --env-file $ENV_FILE -B ${WORKSPACE}/build:/build $SINGULARITY_IMAGE $EXAMPLE $SELF_OPTS"
echo $CMD
$CMD
