#!/bin/bash

INSTALL_ROOT=/apps/self


${INSTALL_ROOT}/bin/self --tolerance "1.0E-5" blockmesh_1d
${INSTALL_ROOT}/bin/self --tolerance "1.0E-5" blockmesh_2d
${INSTALL_ROOT}/bin/self --tolerance "1.0E-5" blockmesh_3d

## 1-D (Scalar) Interpolation
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s1d_interp

# Exactness (Linear function)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --function "f=x" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s1d_interp

# Exponential error decay
${INSTALL_ROOT}/bin/self --tolerance "4.2E-1" \
                    --function "f=sin(6.0*pi*x)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s1d_interp

${INSTALL_ROOT}/bin/self --tolerance "4.2E-2" \
                    --function "f=sin(6.0*pi*x)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 3 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s1d_interp

${INSTALL_ROOT}/bin/self --tolerance "4.2E-3" \
                    --function "f=sin(6.0*pi*x)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 4 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s1d_interp

${INSTALL_ROOT}/bin/self --tolerance "4.2E-4" \
                    --function "f=sin(6.0*pi*x)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 5 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s1d_interp

${INSTALL_ROOT}/bin/self --tolerance "4.2E-5" \
                    --function "f=sin(6.0*pi*x)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 6 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s1d_interp

${INSTALL_ROOT}/bin/self --tolerance "4.2E-6" \
                    --function "f=sin(6.0*pi*x)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 7 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s1d_interp

## 2-D (Scalar) Interpolation
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s2d_interp

# Exactness linear
${INSTALL_ROOT}/bin/self --tolerance "6.0E-6" \
                    --function "f=x*y" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s2d_interp

# Exponential Error Decay
${INSTALL_ROOT}/bin/self --tolerance "6.8E-2" \
                    --function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s2d_interp

${INSTALL_ROOT}/bin/self --tolerance "8.0E-3" \
                    --function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 3 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s2d_interp

${INSTALL_ROOT}/bin/self --tolerance "6.8E-4" \
                    --function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 4 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s2d_interp

${INSTALL_ROOT}/bin/self --tolerance "6.8E-5" \
                    --function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 5 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s2d_interp

${INSTALL_ROOT}/bin/self --tolerance "6.8E-6" \
                    --function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 6 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s2d_interp

${INSTALL_ROOT}/bin/self --tolerance "0.0E0" \
                    --function "f=sin(6.0*pi*x)*sin(6.0*pi*y)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 7 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s2d_interp

## 3-D (Scalar) Interpolation

# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s3d_interp

# Exactness (Linear)
${INSTALL_ROOT}/bin/self --tolerance "1.1E-6" \
                    --function "f=x*y*z" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s3d_interp

# Exponential Error Decay
${INSTALL_ROOT}/bin/self --tolerance "8.6E-2" \
                    --function "f=sin(6.0*pi*x)sin(6.0*pi*y)sin(6.0*pi*z)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s3d_interp

${INSTALL_ROOT}/bin/self --tolerance "2.1E-2" \
                    --function "f=sin(6.0*pi*x)sin(6.0*pi*y)sin(6.0*pi*z)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 3 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s3d_interp

${INSTALL_ROOT}/bin/self --tolerance "5.7E-3" \
                    --function "f=sin(6.0*pi*x)sin(6.0*pi*y)sin(6.0*pi*z)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 4 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s3d_interp

${INSTALL_ROOT}/bin/self --tolerance "1.6E-3" \
                    --function "f=sin(6.0*pi*x)sin(6.0*pi*y)sin(6.0*pi*z)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 5 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s3d_interp

${INSTALL_ROOT}/bin/self --tolerance "2.0E-4" \
                    --function "f=sin(6.0*pi*x)sin(6.0*pi*y)sin(6.0*pi*z)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 6 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s3d_interp

${INSTALL_ROOT}/bin/self --tolerance "0.0E0" \
                    --function "f=sin(6.0*pi*x)sin(6.0*pi*y)sin(6.0*pi*z)" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 7 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s3d_interp

## 2-D (Vector) Interpolation

# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --vector-x "vx=1.0" \
                    --vector-y "vy=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    v2d_interp

# Exactness (Linear)
# Exponential Error Decay

## 3-D (Vector) Interpolation
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --vector-x "vx=1.0" \
                    --vector-y "vy=1.0" \
                    --vector-z "vz=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    v3d_interp

# Exactness (Linear)
# Exponential Error Decay

## 2-D (Tensor) Interpolation
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --tensor-11 "t11=1.0" \
                    --tensor-12 "t12=1.0" \
                    --tensor-21 "t21=1.0" \
                    --tensor-22 "t22=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    t2d_interp

# Exactness (Linear)
# Exponential Error Decay

## 3-D (Tensor) Interpolation
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --tensor-11 "t11=1.0" \
                    --tensor-12 "t12=1.0" \
                    --tensor-13 "t13=1.0" \
                    --tensor-21 "t21=1.0" \
                    --tensor-22 "t22=1.0" \
                    --tensor-23 "t23=1.0" \
                    --tensor-31 "t31=1.0" \
                    --tensor-32 "t32=1.0" \
                    --tensor-33 "t33=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    t3d_interp
# Exactness (Linear)
# Exponential Error Decay

## 1-D (Scalar) Boundary Interpolation
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s1d_binterp

# Exactness (Linear)
# Exponential Error Decay

## 2-D (Scalar) Boundary Interpolation
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s2d_binterp

# Exactness (Linear)
# Exponential Error Decay

## 3-D (Scalar) Boundary Interpolation
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s3d_binterp

# Exactness (Linear)
# Exponential Error Decay

## 2-D (Vector) Boundary Interpolation
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --vector-x "vx=1.0" \
                    --vector-y "vy=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    v2d_binterp

# Exactness (Linear)
# Exponential Error Decay

## 3-D (Vector) Boundary Interpolation
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --vector-x "vx=1.0" \
                    --vector-y "vy=1.0" \
                    --vector-z "vz=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    v3d_binterp

# Exactness (Linear)
# Exponential Error Decay

## 2-D (Tensor) Boundary Interpolation
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --tensor-11 "t11=1.0" \
                    --tensor-12 "t12=1.0" \
                    --tensor-21 "t21=1.0" \
                    --tensor-22 "t22=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    t2d_binterp

# Exactness (Linear)
# Exponential Error Decay

## 3-D (Tensor) Boundary Interpolation
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --tensor-11 "t11=1.0" \
                    --tensor-12 "t12=1.0" \
                    --tensor-13 "t13=1.0" \
                    --tensor-21 "t21=1.0" \
                    --tensor-22 "t22=1.0" \
                    --tensor-23 "t23=1.0" \
                    --tensor-31 "t31=1.0" \
                    --tensor-32 "t32=1.0" \
                    --tensor-33 "t33=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    t3d_binterp

# Exactness (Linear)
# Exponential Error Decay

## 1-D (Scalar) Derivative (Strong Form)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --derivative "df=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "strong" \
                    s1d_derivative

# Exactness (Linear)
# Exponential Error Decay

## 1-D (Scalar) Derivative (DG)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "9.6E-6" \
                    --function "f=1.0" \
                    --derivative "df=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "dg" \
                    s1d_derivative

# Exactness (Linear)
# Exponential Error Decay

## 2-D (Scalar) Gradient (Strong Form)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "2.0E-4" \
                    --function "f=1.0" \
                    --vector-x "gx=0.0" \
                    --vector-y "gy=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "strong" \
                    s2d_gradient

# Exactness (Linear)
# Exponential Error Decay

## 2-D (Scalar) Gradient (DG)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "5.8E-4" \
                    --function "f=1.0" \
                    --vector-x "gx=0.0" \
                    --vector-y "gy=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "dg" \
                    s2d_gradient

# Exactness (Linear)
# Exponential Error Decay

## 3-D (Scalar) Gradient (Strong)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "4.4E-4" \
                    --function "f=1.0" \
                    --vector-x "gx=0.0" \
                    --vector-y "gy=0.0" \
                    --vector-z "gz=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "strong" \
                    s3d_gradient

# Exactness (Linear)
# Exponential Error Decay

## 3-D (Scalar) Gradient (DG)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "7.9E-4" \
                    --function "f=1.0" \
                    --vector-x "gx=0.0" \
                    --vector-y "gy=0.0" \
                    --vector-z "gz=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "dg" \
                    s3d_gradient

# Exactness (Linear)
# Exponential Error Decay

## 2-D (Vector) Gradient (Strong Form)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.9E-4" \
                    --vector-x "vx=1.0" \
                    --vector-y "vy=1.0" \
                    --tensor-11 "vxx=0.0" \
                    --tensor-12 "vxy=0.0" \
                    --tensor-21 "vyx=0.0" \
                    --tensor-22 "vyy=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "strong" \
                    v2d_gradient

# Exactness (Linear)
# Exponential Error Decay

## 2-D (Vector) Gradient (DG)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "5.8E-4" \
                    --vector-x "vx=1.0" \
                    --vector-y "vy=1.0" \
                    --tensor-11 "vxx=0.0" \
                    --tensor-12 "vxy=0.0" \
                    --tensor-21 "vyx=0.0" \
                    --tensor-22 "vyy=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "dg" \
                    v2d_gradient

# Exactness (Linear)
# Exponential Error Decay

## 2-D (Vector) Divergence (Strong Form)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "2.78E-4" \
                    --vector-x "vx=1.0" \
                    --vector-y "vy=1.0" \
                    --function "divV=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "strong" \
                    v2d_divergence

# Exactness (Linear)
# Exponential Error Decay

## 2-D (Vector) Divergence (DG)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "2.75E-4" \
                    --vector-x "vx=1.0" \
                    --vector-y "vy=1.0" \
                    --function "divV=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "dg" \
                    v2d_divergence

# Exactness (Linear)
# Exponential Error Decay

## 3-D (Vector) Gradient (Strong Form)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "4.5E-4" \
                    --vector-x "vx=1.0" \
                    --vector-y "vy=1.0" \
                    --vector-z "vz=1.0" \
                    --tensor-11 "vxx=0.0" \
                    --tensor-12 "vxy=0.0" \
                    --tensor-13 "vxz=0.0" \
                    --tensor-21 "vyx=0.0" \
                    --tensor-22 "vyy=0.0" \
                    --tensor-23 "vyz=0.0" \
                    --tensor-31 "vzx=0.0" \
                    --tensor-32 "vzy=0.0" \
                    --tensor-33 "vzz=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "strong" \
                    v3d_gradient

# Exactness (Linear)
# Exponential Error Decay

## 3-D (Vector) Gradient (DG)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "7.9E-4" \
                    --vector-x "vx=1.0" \
                    --vector-y "vy=1.0" \
                    --vector-z "vz=1.0" \
                    --tensor-11 "vxx=0.0" \
                    --tensor-12 "vxy=0.0" \
                    --tensor-13 "vxz=0.0" \
                    --tensor-21 "vyx=0.0" \
                    --tensor-22 "vyy=0.0" \
                    --tensor-23 "vyz=0.0" \
                    --tensor-31 "vzx=0.0" \
                    --tensor-32 "vzy=0.0" \
                    --tensor-33 "vzz=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "dg" \
                    v3d_gradient

# Exactness (Linear)
# Exponential Error Decay

## 3-D (Vector) Divergence (Strong)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.0E-3" \
                    --vector-x "vx=1.0" \
                    --vector-y "vy=1.0" \
                    --vector-z "vz=1.0" \
                    --function "divV=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "strong" \
                    v3d_divergence

# Exactness (Linear)
# Exponential Error Decay

## 3-D (Vector) Divergence (DG)
# Exactness (Constant)
${INSTALL_ROOT}/bin/self --tolerance "1.2E-3" \
                    --vector-x "vx=1.0" \
                    --vector-y "vy=1.0" \
                    --vector-z "vz=1.0" \
                    --function "divV=0.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    --derivative-type "dg" \
                    v3d_divergence

# Exactness (Linear)
# Exponential Error Decay

