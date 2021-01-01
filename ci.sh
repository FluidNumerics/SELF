
/apps/self/bin/self --tolerance "1.0E-5" blockmesh_1d
/apps/self/bin/self --tolerance "1.0E-5" blockmesh_2d
/apps/self/bin/self --tolerance "1.0E-5" blockmesh_3d

/apps/self/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s1d_interp

/apps/self/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s2d_interp

/apps/self/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s3d_interp

/apps/self/bin/self --tolerance "1.0E-6" \
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

/apps/self/bin/self --tolerance "1.0E-6" \
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

/apps/self/bin/self --tolerance "1.0E-6" \
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

/apps/self/bin/self --tolerance "1.0E-6" \
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

/apps/self/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s1d_binterp

/apps/self/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s2d_binterp

/apps/self/bin/self --tolerance "1.0E-6" \
                    --function "f=1.0" \
                    --nvar 5 \
                    --nelements 10 \
                    --control-degree 2 \
                    --control-quadrature "gauss" \
                    --target-degree 7 \
                    --target-quadrature "gauss" \
                    --gpu-accel "false" \
                    s3d_binterp

/apps/self/bin/self --tolerance "1.0E-6" \
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

/apps/self/bin/self --tolerance "1.0E-6" \
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

/apps/self/bin/self --tolerance "1.0E-6" \
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

/apps/self/bin/self --tolerance "1.0E-6" \
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

/apps/self/bin/self --tolerance "1.0E-6" \
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

/apps/self/bin/self --tolerance "9.6E-6" \
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

/apps/self/bin/self --tolerance "2.0E-4" \
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

/apps/self/bin/self --tolerance "5.8E-4" \
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

/apps/self/bin/self --tolerance "4.4E-4" \
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
                    s3d_gradient


/apps/self/bin/self --tolerance "1.9E-4" \
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

/apps/self/bin/self --tolerance "5.8E-4" \
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

/apps/self/bin/self --tolerance "2.78E-4" \
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

/apps/self/bin/self --tolerance "2.75E-4" \
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

/apps/self/bin/self --tolerance "4.5E-4" \
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
                    v3d_gradient

/apps/self/bin/self --tolerance "1.0E-3" \
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

/apps/self/bin/self --tolerance "1.2E-3" \
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
