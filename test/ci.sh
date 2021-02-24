#!/bin/bash
COUNTER=0
trap '(( $? && ++errcount))' DEBUG

set -x

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "2.5E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "2.5E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "2.5E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "2.5E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "2.5E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "2.5E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=x" \
--vector-y "vy=y" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=x" \
--vector-y "vy=y" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=x" \
--vector-y "vy=y" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=x" \
--vector-y "vy=y" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=x" \
--vector-y "vy=y" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=x" \
--vector-y "vy=y" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=y" \
--vector-y "vy=-x" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=y" \
--vector-y "vy=-x" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=y" \
--vector-y "vy=-x" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=y" \
--vector-y "vy=-x" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=y" \
--vector-y "vy=-x" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=y" \
--vector-y "vy=-x" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
v2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.5E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.5E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.5E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.5E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.5E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.5E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
--vector-z "vz=exp(-z^2)" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
--vector-z "vz=exp(-z^2)" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
--vector-z "vz=exp(-z^2)" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
--vector-z "vz=exp(-z^2)" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
--vector-z "vz=exp(-z^2)" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
--vector-z "vz=exp(-z^2)" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-13 "t13=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t23=1.0" \
--tensor-31 "t31=1.0" \
--tensor-32 "t32=1.0" \
--tensor-33 "t33=1.0" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-13 "t13=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t23=1.0" \
--tensor-31 "t31=1.0" \
--tensor-32 "t32=1.0" \
--tensor-33 "t33=1.0" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-13 "t13=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t23=1.0" \
--tensor-31 "t31=1.0" \
--tensor-32 "t32=1.0" \
--tensor-33 "t33=1.0" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-13 "t13=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t23=1.0" \
--tensor-31 "t31=1.0" \
--tensor-32 "t32=1.0" \
--tensor-33 "t33=1.0" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-13 "t13=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t23=1.0" \
--tensor-31 "t31=1.0" \
--tensor-32 "t32=1.0" \
--tensor-33 "t33=1.0" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-13 "t13=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t23=1.0" \
--tensor-31 "t31=1.0" \
--tensor-32 "t32=1.0" \
--tensor-33 "t33=1.0" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.6E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-13 "t13=x*z" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
--tensor-23 "t23=y*z" \
--tensor-31 "t31=z*x" \
--tensor-32 "t32=z*y" \
--tensor-33 "t33=z" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.6E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-13 "t13=x*z" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
--tensor-23 "t23=y*z" \
--tensor-31 "t31=z*x" \
--tensor-32 "t32=z*y" \
--tensor-33 "t33=z" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.6E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-13 "t13=x*z" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
--tensor-23 "t23=y*z" \
--tensor-31 "t31=z*x" \
--tensor-32 "t32=z*y" \
--tensor-33 "t33=z" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.6E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-13 "t13=x*z" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
--tensor-23 "t23=y*z" \
--tensor-31 "t31=z*x" \
--tensor-32 "t32=z*y" \
--tensor-33 "t33=z" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.6E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-13 "t13=x*z" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
--tensor-23 "t23=y*z" \
--tensor-31 "t31=z*x" \
--tensor-32 "t32=z*y" \
--tensor-33 "t33=z" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.6E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-13 "t13=x*z" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
--tensor-23 "t23=y*z" \
--tensor-31 "t31=z*x" \
--tensor-32 "t32=z*y" \
--tensor-33 "t33=z" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-13 "t13=sin(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-23 "t23=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-31 "t31=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-32 "t32=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-33 "t33=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-13 "t13=sin(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-23 "t23=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-31 "t31=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-32 "t32=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-33 "t33=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-13 "t13=sin(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-23 "t23=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-31 "t31=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-32 "t32=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-33 "t33=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-13 "t13=sin(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-23 "t23=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-31 "t31=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-32 "t32=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-33 "t33=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-13 "t13=sin(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-23 "t23=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-31 "t31=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-32 "t32=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-33 "t33=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-13 "t13=sin(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-23 "t23=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-31 "t31=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-32 "t32=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-33 "t33=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
t3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "9.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "9.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "9.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "9.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "9.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "9.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.16" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.016" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0016" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.00016" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.6000000000000003e-05" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.6000000000000001e-06" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.25" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.025" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0025" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.00025" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "2.5e-05" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "2.5e-06" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=x" \
--vector-y "vy=y" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=x" \
--vector-y "vy=y" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=x" \
--vector-y "vy=y" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=x" \
--vector-y "vy=y" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=x" \
--vector-y "vy=y" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=x" \
--vector-y "vy=y" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=y" \
--vector-y "vy=-x" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=y" \
--vector-y "vy=-x" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=y" \
--vector-y "vy=-x" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=y" \
--vector-y "vy=-x" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=y" \
--vector-y "vy=-x" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=y" \
--vector-y "vy=-x" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
v2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
--vector-z "vz=exp(-z^2)" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
--vector-z "vz=exp(-z^2)" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
--vector-z "vz=exp(-z^2)" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
--vector-z "vz=exp(-z^2)" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
--vector-z "vz=exp(-z^2)" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "vx=-y*exp(-(x^2+y^2))" \
--vector-y "vy=x*exp(-(x^2+y^2))" \
--vector-z "vz=exp(-z^2)" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.19" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.019000000000000003" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0019" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.00019" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.9e-05" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.9000000000000002e-06" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-13 "t13=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t23=1.0" \
--tensor-31 "t31=1.0" \
--tensor-32 "t32=1.0" \
--tensor-33 "t33=1.0" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-13 "t13=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t23=1.0" \
--tensor-31 "t31=1.0" \
--tensor-32 "t32=1.0" \
--tensor-33 "t33=1.0" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-13 "t13=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t23=1.0" \
--tensor-31 "t31=1.0" \
--tensor-32 "t32=1.0" \
--tensor-33 "t33=1.0" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-13 "t13=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t23=1.0" \
--tensor-31 "t31=1.0" \
--tensor-32 "t32=1.0" \
--tensor-33 "t33=1.0" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-13 "t13=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t23=1.0" \
--tensor-31 "t31=1.0" \
--tensor-32 "t32=1.0" \
--tensor-33 "t33=1.0" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t12=1.0" \
--tensor-13 "t13=1.0" \
--tensor-21 "t21=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t23=1.0" \
--tensor-31 "t31=1.0" \
--tensor-32 "t32=1.0" \
--tensor-33 "t33=1.0" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-13 "t13=x*z" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
--tensor-23 "t23=y*z" \
--tensor-31 "t31=z*x" \
--tensor-32 "t32=z*y" \
--tensor-33 "t33=z" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-13 "t13=x*z" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
--tensor-23 "t23=y*z" \
--tensor-31 "t31=z*x" \
--tensor-32 "t32=z*y" \
--tensor-33 "t33=z" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-13 "t13=x*z" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
--tensor-23 "t23=y*z" \
--tensor-31 "t31=z*x" \
--tensor-32 "t32=z*y" \
--tensor-33 "t33=z" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-13 "t13=x*z" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
--tensor-23 "t23=y*z" \
--tensor-31 "t31=z*x" \
--tensor-32 "t32=z*y" \
--tensor-33 "t33=z" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-13 "t13=x*z" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
--tensor-23 "t23=y*z" \
--tensor-31 "t31=z*x" \
--tensor-32 "t32=z*y" \
--tensor-33 "t33=z" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t12=x*y" \
--tensor-13 "t13=x*z" \
--tensor-21 "t21=y*x" \
--tensor-22 "t22=y" \
--tensor-23 "t23=y*z" \
--tensor-31 "t31=z*x" \
--tensor-32 "t32=z*y" \
--tensor-33 "t33=z" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.3" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-13 "t13=sin(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-23 "t23=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-31 "t31=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-32 "t32=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-33 "t33=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.03" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-13 "t13=sin(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-23 "t23=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-31 "t31=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-32 "t32=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-33 "t33=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.003" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-13 "t13=sin(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-23 "t23=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-31 "t31=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-32 "t32=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-33 "t33=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0003" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-13 "t13=sin(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-23 "t23=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-31 "t31=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-32 "t32=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-33 "t33=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "3e-05" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-13 "t13=sin(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-23 "t23=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-31 "t31=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-32 "t32=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-33 "t33=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "3e-06" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--tensor-11 "t11=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-12 "t12=sin(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-13 "t13=sin(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-21 "t21=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-22 "t22=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-23 "t23=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
--tensor-31 "t31=cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--tensor-32 "t32=cos(2.0*pi*x)*cos(2.0*pi*y)" \
--tensor-33 "t33=cos(2.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z)" \
t3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "7.0E-5" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "7.0E-5" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "7.0E-5" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "7.0E-5" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "7.0E-5" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "7.0E-5" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "7.0E-5" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "7.0E-5" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "7.0E-5" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "7.0E-5" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "7.0E-5" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "7.0E-5" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.5E-5" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.5E-5" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.5E-5" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.5E-5" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.5E-5" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.5E-5" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.5E-5" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.5E-5" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.5E-5" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.5E-5" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.5E-5" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.5E-5" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.8" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
--derivative "df=2.0*pi*cos(2.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.08000000000000002" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
--derivative "df=2.0*pi*cos(2.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.008" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
--derivative "df=2.0*pi*cos(2.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0008" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
--derivative "df=2.0*pi*cos(2.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8e-05" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
--derivative "df=2.0*pi*cos(2.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.000000000000001e-06" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
--derivative "df=2.0*pi*cos(2.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.8" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
--derivative "df=2.0*pi*cos(2.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.08000000000000002" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
--derivative "df=2.0*pi*cos(2.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.008" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
--derivative "df=2.0*pi*cos(2.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0008" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
--derivative "df=2.0*pi*cos(2.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8e-05" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
--derivative "df=2.0*pi*cos(2.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.000000000000001e-06" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)" \
--derivative "df=2.0*pi*cos(2.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.2E-1" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.5" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.5" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.5" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-z "fy=2.0*pi*sin(2.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-z "fy=2.0*pi*sin(2.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-z "fy=2.0*pi*sin(2.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-z "fy=2.0*pi*sin(2.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-z "fy=2.0*pi*sin(2.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-z "fy=2.0*pi*sin(2.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.5" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-z "fy=2.0*pi*sin(2.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-z "fy=2.0*pi*sin(2.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-z "fy=2.0*pi*sin(2.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-z "fy=2.0*pi*sin(2.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-z "fy=2.0*pi*sin(2.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--function "f=sin(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-x "fx=2.0*pi*cos(2.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-y "fy=2.0*pi*sin(2.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z)" \
--vector-z "fy=2.0*pi*sin(2.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2))" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2))" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2))" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2))" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2))" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2))" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2))" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2))" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2))" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2))" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2))" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2))" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
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
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
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
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
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
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
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
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
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
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
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
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.08" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.008" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0008" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.000000000000001e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.08" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.008" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0008" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.000000000000001e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "8.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-21 "dfyx=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=1.0" \
--tensor-21 "dfyx=-1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=1.0" \
--tensor-21 "dfyx=-1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=1.0" \
--tensor-21 "dfyx=-1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=1.0" \
--tensor-21 "dfyx=-1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=1.0" \
--tensor-21 "dfyx=-1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=1.0" \
--tensor-21 "dfyx=-1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=1.0" \
--tensor-21 "dfyx=-1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=1.0" \
--tensor-21 "dfyx=-1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=1.0" \
--tensor-21 "dfyx=-1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=1.0" \
--tensor-21 "dfyx=-1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=1.0" \
--tensor-21 "dfyx=-1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "6.0E-2" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=1.0" \
--tensor-21 "dfyx=-1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.5" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.5" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
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
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfyz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfyz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfyz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfyz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfyz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfyz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfyz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfyz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfyz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfyz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfyz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfyz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfxz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfxz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfxz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfxz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfxz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfxz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfxz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfxz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfxz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfxz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfxz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-1" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfxy=0.0" \
--tensor-13 "dfxz=0.0" \
--tensor-21 "dfyx=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfxz=0.0" \
--tensor-31 "dfzx=0.0" \
--tensor-32 "dfzy=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.5" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.0" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.5" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 1 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient



set +x

echo "============================="
echo "        Test Summary         "
echo "============================="
echo " "
echo " Total Tests : $COUNTER "
echo " Fail : $errcount "
#test $errcount = 0
