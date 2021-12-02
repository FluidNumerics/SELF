#!/bin/bash

source /etc/profile.d/z10_spack_environment.sh
spack load singularity

mkdir -p ${WORKSPACE}/codecov

singularity run --nv --bind ${WORKSPACE}/codecov:/build \
                            ${SINGULARITY_IMAGE} /opt/self/bin/sadv2d \
			    --gpu \
			    -dt 0.001 \
			    -tn 10.0 \
			    -t0 0.0 \
			    -oi 1.0 \
                            -vx "vx=1.0" \
			    -vy "vy=0.0" \
			    -ic "f=exp( -(x^2 + (y-2.5)^2)/0.5 )" \
                            -bc "f=0.0" \
                            --mesh "${WORKSPACE}/test/sadv2d/gaussian_spline_hohqmesh/SplineMesh.mesh"
