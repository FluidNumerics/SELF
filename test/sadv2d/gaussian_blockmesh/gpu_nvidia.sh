#!/bin/bash

source /etc/profile.d/z10_spack_environment.sh
spack load singularity

mkdir -p ${WORKSPACE}/ci/test/codecov

singularity run --nv --bind ${WORKSPACE}/ci/test/codecov:/build \
                            ${SINGULARITY_IMAGE} /opt/self/bin/sadv2d \
			    -gpu \
			    -dt 0.001 \
			    -tn 1.0 \
			    -t0 0.0 \
			    -oi 0.5
