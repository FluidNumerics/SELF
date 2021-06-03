#!/bin/bash
COUNTER=0
trap '(( $? && ++errcount))' DEBUG

set -x

mkdir -p /workspace/v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
--output "/workspace/v3d_divergence/self.h5" \
v3d_divergence
echo $?

set +x


FAILCOUNT=$(( $errcount-1 ))

echo "============================="
echo "        Test Summary         "
echo "============================="
echo " "
echo " Total Tests : $COUNTER "
echo " Fail : $FAILCOUNT "

echo $FAILCOUNT > /workspace/failcount.txt