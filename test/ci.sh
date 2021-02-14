#!/bin/bash
COUNTER=0
trap '(( $? && ++errcount))' DEBUG

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
s1d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
s2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
s3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--vector-x "vx=1.0" \
--vector-y "vy=1.0" \
--vector-z "vz=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "vx=-y" \
--vector-y "vy=x" \
--vector-z "vz=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
 t2d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-13 "t31=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t32=1.0" \
--tensor-31 "t13=1.0" \
--tensor-32 "t23=1.0" \
--tensor-33 "t33=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-13 "t31=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t32=1.0" \
--tensor-31 "t13=1.0" \
--tensor-32 "t23=1.0" \
--tensor-33 "t33=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-13 "t31=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t32=1.0" \
--tensor-31 "t13=1.0" \
--tensor-32 "t23=1.0" \
--tensor-33 "t33=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-13 "t31=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t32=1.0" \
--tensor-31 "t13=1.0" \
--tensor-32 "t23=1.0" \
--tensor-33 "t33=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-13 "t31=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t32=1.0" \
--tensor-31 "t13=1.0" \
--tensor-32 "t23=1.0" \
--tensor-33 "t33=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-13 "t31=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t32=1.0" \
--tensor-31 "t13=1.0" \
--tensor-32 "t23=1.0" \
--tensor-33 "t33=1.0" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-13 "t31=z*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
--tensor-23 "t32=z*y" \
--tensor-31 "t13=x*z" \
--tensor-32 "t23=y*z" \
--tensor-33 "t33=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-13 "t31=z*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
--tensor-23 "t32=z*y" \
--tensor-31 "t13=x*z" \
--tensor-32 "t23=y*z" \
--tensor-33 "t33=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-13 "t31=z*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
--tensor-23 "t32=z*y" \
--tensor-31 "t13=x*z" \
--tensor-32 "t23=y*z" \
--tensor-33 "t33=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-13 "t31=z*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
--tensor-23 "t32=z*y" \
--tensor-31 "t13=x*z" \
--tensor-32 "t23=y*z" \
--tensor-33 "t33=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-13 "t31=z*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
--tensor-23 "t32=z*y" \
--tensor-31 "t13=x*z" \
--tensor-32 "t23=y*z" \
--tensor-33 "t33=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-13 "t31=z*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
--tensor-23 "t32=z*y" \
--tensor-31 "t13=x*z" \
--tensor-32 "t23=y*z" \
--tensor-33 "t33=z" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-13 "t31=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-23 "t32=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-31 "t13=sin(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-32 "t23=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-33 "t33=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-13 "t31=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-23 "t32=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-31 "t13=sin(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-32 "t23=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-33 "t33=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-13 "t31=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-23 "t32=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-31 "t13=sin(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-32 "t23=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-33 "t33=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-13 "t31=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-23 "t32=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-31 "t13=sin(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-32 "t23=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-33 "t33=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-13 "t31=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-23 "t32=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-31 "t13=sin(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-32 "t23=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-33 "t33=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-13 "t31=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-23 "t32=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-31 "t13=sin(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-32 "t23=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-33 "t33=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
v3d_interp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
s1d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
s2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
s3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
 t2d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-13 "t31=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t32=1.0" \
--tensor-31 "t13=1.0" \
--tensor-32 "t23=1.0" \
--tensor-33 "t33=1.0" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-13 "t31=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t32=1.0" \
--tensor-31 "t13=1.0" \
--tensor-32 "t23=1.0" \
--tensor-33 "t33=1.0" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-13 "t31=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t32=1.0" \
--tensor-31 "t13=1.0" \
--tensor-32 "t23=1.0" \
--tensor-33 "t33=1.0" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-13 "t31=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t32=1.0" \
--tensor-31 "t13=1.0" \
--tensor-32 "t23=1.0" \
--tensor-33 "t33=1.0" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-13 "t31=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t32=1.0" \
--tensor-31 "t13=1.0" \
--tensor-32 "t23=1.0" \
--tensor-33 "t33=1.0" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=1.0" \
--tensor-12 "t21=1.0" \
--tensor-13 "t31=1.0" \
--tensor-21 "t12=1.0" \
--tensor-22 "t22=1.0" \
--tensor-23 "t32=1.0" \
--tensor-31 "t13=1.0" \
--tensor-32 "t23=1.0" \
--tensor-33 "t33=1.0" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-13 "t31=z*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
--tensor-23 "t32=z*y" \
--tensor-31 "t13=x*z" \
--tensor-32 "t23=y*z" \
--tensor-33 "t33=z" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-13 "t31=z*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
--tensor-23 "t32=z*y" \
--tensor-31 "t13=x*z" \
--tensor-32 "t23=y*z" \
--tensor-33 "t33=z" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-13 "t31=z*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
--tensor-23 "t32=z*y" \
--tensor-31 "t13=x*z" \
--tensor-32 "t23=y*z" \
--tensor-33 "t33=z" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-13 "t31=z*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
--tensor-23 "t32=z*y" \
--tensor-31 "t13=x*z" \
--tensor-32 "t23=y*z" \
--tensor-33 "t33=z" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-13 "t31=z*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
--tensor-23 "t32=z*y" \
--tensor-31 "t13=x*z" \
--tensor-32 "t23=y*z" \
--tensor-33 "t33=z" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=x" \
--tensor-12 "t21=y*x" \
--tensor-13 "t31=z*x" \
--tensor-21 "t12=x*y" \
--tensor-22 "t22=y" \
--tensor-23 "t32=z*y" \
--tensor-31 "t13=x*z" \
--tensor-32 "t23=y*z" \
--tensor-33 "t33=z" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-13 "t31=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-23 "t32=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-31 "t13=sin(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-32 "t23=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-33 "t33=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-13 "t31=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-23 "t32=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-31 "t13=sin(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-32 "t23=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-33 "t33=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-13 "t31=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-23 "t32=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-31 "t13=sin(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-32 "t23=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-33 "t33=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-13 "t31=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-23 "t32=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-31 "t13=sin(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-32 "t23=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-33 "t33=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-13 "t31=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-23 "t32=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-31 "t13=sin(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-32 "t23=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-33 "t33=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--tensor-11 "t11=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-12 "t21=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-13 "t31=cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--tensor-21 "t12=sin(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-22 "t22=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-23 "t32=cos(6.0*pi*x)*cos(6.0*pi*y)" \
--tensor-31 "t13=sin(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-32 "t23=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
--tensor-33 "t33=cos(6.0*pi*x)*cos(6.0*pi*y)*cos(6.0*pi*z)" \
v3d_binterp

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--derivative "df=0.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x" \
--derivative "df=1.0" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
--derivative "df=6.0*pi*cos(6.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
--derivative "df=6.0*pi*cos(6.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
--derivative "df=6.0*pi*cos(6.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
--derivative "df=6.0*pi*cos(6.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
--derivative "df=6.0*pi*cos(6.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
--derivative "df=6.0*pi*cos(6.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
--derivative "df=6.0*pi*cos(6.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
--derivative "df=6.0*pi*cos(6.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
--derivative "df=6.0*pi*cos(6.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
--derivative "df=6.0*pi*cos(6.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
--derivative "df=6.0*pi*cos(6.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)" \
--derivative "df=6.0*pi*cos(6.0*pi*x)" \
s1d_derivative

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y" \
--vector-x "fx=y" \
--vector-y "fy=x" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)" \
s2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=1.0" \
--vector-x "fx=0.0" \
--vector-y "fy=0.0" \
--vector-z "fz=0.0" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=x*y*z" \
--vector-x "fx=y*z" \
--vector-y "fy=x*z" \
--vector-z "fz=x*y" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-z "fy=6.0*pi*sin(6.0*pi*x)*sin(6.0*pi*y)*cos(6.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-z "fy=6.0*pi*sin(6.0*pi*x)*sin(6.0*pi*y)*cos(6.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-z "fy=6.0*pi*sin(6.0*pi*x)*sin(6.0*pi*y)*cos(6.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-z "fy=6.0*pi*sin(6.0*pi*x)*sin(6.0*pi*y)*cos(6.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-z "fy=6.0*pi*sin(6.0*pi*x)*sin(6.0*pi*y)*cos(6.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-z "fy=6.0*pi*sin(6.0*pi*x)*sin(6.0*pi*y)*cos(6.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-z "fy=6.0*pi*sin(6.0*pi*x)*sin(6.0*pi*y)*cos(6.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-z "fy=6.0*pi*sin(6.0*pi*x)*sin(6.0*pi*y)*cos(6.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-z "fy=6.0*pi*sin(6.0*pi*x)*sin(6.0*pi*y)*cos(6.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-z "fy=6.0*pi*sin(6.0*pi*x)*sin(6.0*pi*y)*cos(6.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-z "fy=6.0*pi*sin(6.0*pi*x)*sin(6.0*pi*y)*cos(6.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--function "f=sin(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-x "fx=6.0*pi*cos(6.0*pi*x)*sin(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-y "fy=6.0*pi*sin(6.0*pi*x)*cos(6.0*pi*y)*sin(6.0*pi*z)" \
--vector-z "fy=6.0*pi*sin(6.0*pi*x)*sin(6.0*pi*y)*cos(6.0*pi*z)" \
s3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--function "df=2.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--function "df=0.0" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--function "df=-2.0*(x+y)*exp(-(x^2+y^2))" \
v2d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--function "df=3.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--vector-z "fz=x*y" \
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
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
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
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y*exp(-(x^2+y^2+z^2))" \
--vector-y "fy=-x*exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2))" \
--function "df=0.0" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
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
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--vector-z "fz=exp(-(x^2+y^2+z^2))" \
--function "df=-2.0*(x+y+z)*exp(-(x^2+y^2+z^2))" \
v3d_divergence

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-21 "dfxy=0.0" \
 --tensor-22 "dfyy=1.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2))" \
--vector-y "fy=exp(-(x^2+y^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2))" \
 --tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2))" \
v2d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfyz=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfyz=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfyz=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfyz=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfyz=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfyz=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfyz=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfyz=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfyz=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfyz=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfyz=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=1.0" \
--vector-y "fy=1.0" \
--vector-z "fz=1.0" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=0.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfyz=0.0" \
--tensor-33 "dfzz=0.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfxz=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfxz=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfxz=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfxz=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfxz=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfxz=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfxz=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfxz=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfxz=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfxz=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfxz=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=x" \
--vector-y "fy=y" \
--vector-z "fz=z" \
--tensor-11 "dfxx=1.0" \
--tensor-12 "dfyx=0.0" \
--tensor-13 "dfzx=0.0" \
--tensor-21 "dfxy=0.0" \
--tensor-22 "dfyy=1.0" \
--tensor-23 "dfzy=0.0" \
--tensor-31 "dfxz=0.0" \
--tensor-32 "dfxz=0.0" \
--tensor-33 "dfzz=1.0" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=y" \
--vector-y "fy=-x" \
--tensor-11 "dfxx=0.0" \
--tensor-12 "dfyx=-1.0" \
--tensor-21 "dfxy=1.0" \
 --tensor-22 "dfyy=0.0" \
v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "strong" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.05" \
--control-quadrature "gauss" \
--control-degree 2 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.005000000000000001" \
--control-quadrature "gauss" \
--control-degree 3 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "0.0005" \
--control-quadrature "gauss" \
--control-degree 4 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-05" \
--control-quadrature "gauss" \
--control-degree 5 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5e-06" \
--control-quadrature "gauss" \
--control-degree 6 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient

let COUNTER++ 
${INSTALL_ROOT}/bin/self --tolerance "5.000000000000001e-07" \
--control-quadrature "gauss" \
--control-degree 7 \
--target-quadrature "gauss" \
--target-degree 7 \
--derivative-type "dg" \
--nelements 10 \
--nvar 5 \
--gpu-accel "false" \
--vector-x "fx=exp(-(x^2+y^2+z^2))" \
--vector-y "fy=exp(-(x^2+y^2+z^2))" \
--tensor-11 "dfxx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-12 "dfyx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-13 "dfzx=-2.0*x*exp(-(x^2+y^2+z^2))" \
--tensor-21 "dfxy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-22 "dfyy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-23 "dfzy=-2.0*y*exp(-(x^2+y^2+z^2))" \
--tensor-31 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-32 "dfxz=-2.0*z*exp(-(x^2+y^2+z^2))" \
--tensor-33 "dfzz=-2.0*z*exp(-(x^2+y^2+z^2))" \
 v3d_gradient



echo "============================="
echo "        Test Summary         "
echo "============================="
echo " "
echo " Total Tests : $COUNTER "
echo " Fail : $errcount "
test $errcount = 0
