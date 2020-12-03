PROGRAM SELF

  USE SELF_Constants
  USE SELF_SupportRoutines
  USE SELF_Tests
  USE FLAP

  IMPLICIT NONE

  TYPE(COMMAND_LINE_INTERFACE) :: self_cli
  CHARACTER(20) :: cqTypeChar
  CHARACTER(20) :: tqTypeChar
  CHARACTER(240) :: functionChar
  CHARACTER(240) :: vx,vy,vz
  CHARACTER(240) :: tensorChar(1:3,1:3)
  REAL(prec) :: errorTolerance
  INTEGER :: cqType
  INTEGER :: tqType
  INTEGER :: cqDegree
  INTEGER :: tqDegree
  INTEGER :: nElem
  INTEGER :: nVar
  INTEGER :: error
  INTEGER :: errorCount

  CALL Parse_CLI(self_cli)

  CALL self_cli % get(val=cqTypeChar,switch='--control-quadrature')
  CALL self_cli % get(val=tqTypeChar,switch='--target-quadrature')
  CALL self_cli % get(val=cqDegree,switch='--control-degree')
  CALL self_cli % get(val=tqDegree,switch='--target-degree')
  CALL self_cli % get(val=nElem,switch='--nelements')
  CALL self_cli % get(val=nVar,switch='--nvar')
  CALL self_cli % get(val=functionChar,switch='--function')
  CALL self_cli % get(val=vx,switch='--vector-x')
  CALL self_cli % get(val=vy,switch='--vector-y')
  CALL self_cli % get(val=vz,switch='--vector-z')
  CALL self_cli % get(val=tensorChar(1,1),switch='--tensor-11')
  CALL self_cli % get(val=tensorChar(1,2),switch='--tensor-12')
  CALL self_cli % get(val=tensorChar(1,3),switch='--tensor-13')
  CALL self_cli % get(val=tensorChar(2,1),switch='--tensor-21')
  CALL self_cli % get(val=tensorChar(2,2),switch='--tensor-22')
  CALL self_cli % get(val=tensorChar(2,3),switch='--tensor-23')
  CALL self_cli % get(val=tensorChar(3,1),switch='--tensor-31')
  CALL self_cli % get(val=tensorChar(3,2),switch='--tensor-32')
  CALL self_cli % get(val=tensorChar(3,3),switch='--tensor-33')
  CALL self_cli % get(val=errorTolerance,switch='--tolerance')

  IF (TRIM(UpperCase(cqTypeChar)) == 'GAUSS') THEN
    cqType = GAUSS
  ELSEIF (TRIM(UpperCase(cqTypeChar)) == 'GAUSS-LOBATTO') THEN
    cqType = GAUSS_LOBATTO
  ELSE
    PRINT *, 'Invalid Control Quadrature'
    STOP 1
  END IF

  IF (TRIM(UpperCase(tqTypeChar)) == 'GAUSS') THEN
    tqType = GAUSS
  ELSEIF (TRIM(UpperCase(tqTypeChar)) == 'GAUSS-LOBATTO') THEN
    tqType = GAUSS_LOBATTO
  ELSE
    PRINT *, 'Invalid Target Quadrature'
    STOP 1
  END IF

  errorCount = 0
  IF (self_cli % run_command(group="ci-test")) THEN

    CALL BlockMesh1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,errorTolerance,error)
    errorCount = errorCount + error

    CALL BlockMesh2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,errorTolerance,error)
    errorCount = errorCount + error

    CALL BlockMesh3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,errorTolerance,error)
    errorCount = errorCount + error

    CALL ScalarInterp1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,functionChar,errorTolerance,error)
    errorCount = errorCount + error

    CALL ScalarBoundaryInterp1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,functionChar,errorTolerance,error)
    errorCount = errorCount + error

    CALL ScalarInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,functionChar,errorTolerance,error)
    errorCount = errorCount + error

    CALL ScalarInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,functionChar,errorTolerance,error)
    errorCount = errorCount + error

    CALL VectorInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,vx,vy,errorTolerance,error)
    errorCount = errorCount + error

    CALL VectorInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,vx,vy,vz,errorTolerance,error)
    errorCount = errorCount + error

    CALL TensorInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,tensorChar(1:2,1:2),errorTolerance,error)
    errorCount = errorCount + error

    CALL TensorInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,tensorChar,errorTolerance,error)
    errorCount = errorCount + error

  ELSEIF (self_cli % run_command(group="blockmesh_1d")) THEN

    CALL BlockMesh1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,errorTolerance,error)
    errorCount = errorCount + error

  ELSEIF (self_cli % run_command(group="blockmesh_2d")) THEN

    CALL BlockMesh2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,errorTolerance,error)
    errorCount = errorCount + error

  ELSEIF (self_cli % run_command(group="blockmesh_3d")) THEN

    CALL BlockMesh3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,errorTolerance,error)
    errorCount = errorCount + error

  ELSEIF (self_cli % run_command(group="s1d_interp")) THEN

    CALL ScalarInterp1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,functionChar,errorTolerance,error)
    errorCount = errorCount + error

  ELSEIF (self_cli % run_command(group="s1d_binterp")) THEN

    CALL ScalarBoundaryInterp1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,functionChar,errorTolerance,error)
    errorCount = errorCount + error

  ELSEIF (self_cli % run_command(group="s2d_interp")) THEN

    CALL ScalarInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,functionChar,errorTolerance,error)
    errorCount = errorCount + error

  ELSEIF (self_cli % run_command(group="s3d_interp")) THEN

    CALL ScalarInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,functionChar,errorTolerance,error)
    errorCount = errorCount + error

  ELSEIF (self_cli % run_command(group="v2d_interp")) THEN

    CALL VectorInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,vx,vy,errorTolerance,error)
    errorCount = errorCount + error

  ELSEIF (self_cli % run_command(group="t2d_interp")) THEN

    CALL TensorInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,tensorChar(1:2,1:2),errorTolerance,error)
    errorCount = errorCount + error

  ELSEIF (self_cli % run_command(group="v3d_interp")) THEN

    CALL VectorInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,vx,vy,vz,errorTolerance,error)
    errorCount = errorCount + error

  ELSEIF (self_cli % run_command(group="t3d_interp")) THEN

    CALL TensorInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,tensorChar,errorTolerance,error)
    errorCount = errorCount + error

  END IF

  CALL self_cli % free()

  STOP errorCount

CONTAINS

  SUBROUTINE Parse_CLI(cli)
    IMPLICIT NONE
    TYPE(COMMAND_LINE_INTERFACE),INTENT(out) :: cli

    CALL cli % init(progname="self", &
                    version="v0.0.0", &
                    description="Spectral Element Libraries in Fortran (SELF)", &
                    license="ANTI-CAPITALIST SOFTWARE LICENSE (v 1.4)", &
                    authors="Joseph Schoonover (Fluid Numerics LLC)")

    CALL cli % add(switch="--control-degree", &
                   switch_ab="-c", &
                   help="The polynomial degree of the control points."//NEW_LINE("A"), &
                   def="2", &
                   required=.FALSE.)

    CALL cli % add(switch="--target-degree", &
                   switch_ab="-t", &
                   help="The polynomial degree for the target points for interpolation."//NEW_LINE("A"), &
                   def="7", &
                   required=.FALSE.)

    CALL cli % add(switch="--control-quadrature", &
                   switch_ab="-cq", &
                   def="gauss", &
                   help="The quadrature type for control points."//NEW_LINE("A"), &
                   choices="gauss, gauss-lobatto, uniform", &
                   required=.FALSE.)

    CALL cli % add(switch="--target-quadrature", &
                   switch_ab="-tq", &
                   def="gauss", &
                   help="The quadrature type for target points."//NEW_LINE("A"), &
                   choices="gauss, gauss-lobatto, uniform", &
                   required=.FALSE.)

    CALL cli % add(switch="--nvar", &
                   switch_ab="-nv", &
                   help="The number of functions used to simulate increased workloads."//NEW_LINE("A"), &
                   def="5", &
                   required=.FALSE.)

    CALL cli % add(switch="--nelements", &
                   switch_ab="-ne", &
                   help="The number of elements to use on the test mesh. "// &
                   "In multi-dimensions, the number of elements in each direction."//NEW_LINE("A"), &
                   def="10", &
                   required=.FALSE.)

    CALL cli % add(switch="--gpu-accel", &
                   switch_ab="-gpu", &
                   help="Boolean flag for enabling or disabling GPU acceleration in tests."//NEW_LINE("A"), &
                   def="false", &
                   choices="true, false", &
                   required=.FALSE.)

    CALL cli % add(switch="--tolerance", &
                   switch_ab="-tol", &
                   help="Tolerance to use for determining if a test passes."//NEW_LINE("A"), &
                   def="1.0E-5", &
                   required=.FALSE.)

    CALL cli % add(switch="--function", &
                   switch_ab="-f", &
                   help="Function to interpolate from control points to target points"//NEW_LINE("A"), &
                   def="f=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--vector-x", &
                   switch_ab="-vx", &
                   help="x-component of the vector function to interpolate from control points "// &
                   "to target points"//NEW_LINE("A"), &
                   def="vx=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--vector-y", &
                   switch_ab="-vy", &
                   help="y-component of the vector function to interpolate from control points "// &
                   "to target points"//NEW_LINE("A"), &
                   def="vy=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--vector-z", &
                   switch_ab="-vz", &
                   help="z-component of the vector function to interpolate from control points "// &
                   "to target points"//NEW_LINE("A"), &
                   def="vz=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--tensor-11", &
                   switch_ab="-t11", &
                   help="Row 1 column 1 of the tensor function to interpolate from control points "// &
                   "to target points"//NEW_LINE("A"), &
                   def="t11=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--tensor-12", &
                   switch_ab="-t12", &
                   help="Row 1 column 2 of the tensor function to interpolate from control points "// &
                   "to target points"//NEW_LINE("A"), &
                   def="t12=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--tensor-13", &
                   switch_ab="-t13", &
                   help="Row 1 column 3 of the tensor function to interpolate from control points "// &
                   "to target points"//NEW_LINE("A"), &
                   def="t13=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--tensor-21", &
                   switch_ab="-t21", &
                   help="Row 2 column 1 of the tensor function to interpolate from control points "// &
                   "to target points"//NEW_LINE("A"), &
                   def="t21=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--tensor-22", &
                   switch_ab="-t22", &
                   help="Row 2 column 2 of the tensor function to interpolate from control points "// &
                   "to target points"//NEW_LINE("A"), &
                   def="t22=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--tensor-23", &
                   switch_ab="-t23", &
                   help="Row 2 column 3 of the tensor function to interpolate from control points "// &
                   "to target points"//NEW_LINE("A"), &
                   def="t23=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--tensor-31", &
                   switch_ab="-t31", &
                   help="Row 3 column 1 of the tensor function to interpolate from control points "// &
                   "to target points"//NEW_LINE("A"), &
                   def="t31=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--tensor-32", &
                   switch_ab="-t32", &
                   help="Row 3 column 2 of the tensor function to interpolate from control points "// &
                   "to target points"//NEW_LINE("A"), &
                   def="t32=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--tensor-33", &
                   switch_ab="-t33", &
                   help="Row 2 column 3 of the tensor function to interpolate from control points "// &
                   "to target points"//NEW_LINE("A"), &
                   def="t33=1.0", &
                   required=.FALSE.)

    CALL cli % add_group(group="ci-test", &
                         description="Run all CI-Tests")

    CALL cli % add_group(group="blockmesh_1d", &
                         description="Block Mesh generation in 1D")

    CALL cli % add_group(group="blockmesh_2d", &
                         description="Block Mesh generation in 2D")

    CALL cli % add_group(group="blockmesh_3d", &
                         description="Block Mesh generation in 3D")

    CALL cli % add_group(group="s1d_interp", &
                         description="Scalar 1D Grid Interpolation")

    CALL cli % add_group(group="s2d_interp", &
                         description="Scalar 2D Grid Interpolation")

    CALL cli % add_group(group="s3d_interp", &
                         description="Scalar 3D Grid Interpolation")

    CALL cli % add_group(group="v2d_interp", &
                         description="Vector 2D Grid Interpolation")

    CALL cli % add_group(group="v3d_interp", &
                         description="Vector 3D Grid Interpolation")

    CALL cli % add_group(group="t2d_interp", &
                         description="Tensor 2D Grid Interpolation")

    CALL cli % add_group(group="t3d_interp", &
                         description="Tensor 3D Grid Interpolation")

    CALL cli % add_group(group="s1d_binterp", &
                         description="Scalar 1D Boundary interpolation")

    CALL cli % add_group(group="s2d_binterp", &
                         description="Scalar 2D Boundary interpolation")

    CALL cli % add_group(group="s3d_binterp", &
                         description="Scalar 3D Boundary interpolation")

    CALL cli % add_group(group="v2d_binterp", &
                         description="Vector 2D Boundary interpolation")

    CALL cli % add_group(group="v3d_binterp", &
                         description="Vector 3D Boundary interpolation")

    CALL cli % add_group(group="t2d_binterp", &
                         description="Tensor 2D Boundary interpolation")

    CALL cli % add_group(group="t3d_binterp", &
                         description="Tensor 3D Boundary interpolation")

    CALL cli % add_group(group="s1d_derivative", &
                         description="Scalar 1D Derivative")

    CALL cli % add_group(group="s2d_gradient", &
                         description="Scalar 2D Gradient")

    CALL cli % add_group(group="s3d_gradient", &
                         description="Scalar 3D Gradient")

    CALL cli % add_group(group="v2d_gradient", &
                         description="Vector 2D Gradient")

    CALL cli % add_group(group="v2d_divergence", &
                         description="Vector 2D Divergence")

    CALL cli % add_group(group="v2d_curl", &
                         description="Vector 2D Curl")

    CALL cli % add_group(group="v3d_gradient", &
                         description="Vector 3D Gradient")

    CALL cli % add_group(group="v3d_divergence", &
                         description="Vector 3D Divergence")

    CALL cli % add_group(group="v3d_curl", &
                         description="Vector 3D Curl")


    CALL cli % parse()

  END SUBROUTINE Parse_CLI

END PROGRAM SELF
