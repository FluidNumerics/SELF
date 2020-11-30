
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
