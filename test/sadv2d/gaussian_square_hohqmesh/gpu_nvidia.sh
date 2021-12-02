#!/bin/bash

source /etc/profile.d/z10_spack_environment.sh
spack load singularity

mkdir -p ${WORKSPACE}/codecov

singularity run --nv --bind ${WORKSPACE}/codecov:/build \
                            ${SINGULARITY_IMAGE} /opt/self/bin/sadv2d \
			    --gpu \
                            --mesh "${WORKSPACE}/test/sadv2d/gaussian_square_hohqmesh/Square.mesh"
