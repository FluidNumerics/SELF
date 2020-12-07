
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
                    s1d_derivative
