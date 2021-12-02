#!/bin/bash

source /etc/profile.d/z10_spack_environment.sh
spack load singularity

mkdir -p ${WORKSPACE}/codecov

singularity run --nv --bind ${WORKSPACE}/codecov:/build \
                            ${SINGULARITY_IMAGE} /opt/self/bin/sadv2d \
			    -gpu \
			    -dt 0.001 \
			    -tn 1.0 \
			    -t0 0.0 \
			    -oi 0.5
                            --mesh "${WORKSPACE}/test/sadv2d/gaussian_square_hohqmesh/Square.mesh"
