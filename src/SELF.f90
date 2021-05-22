!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
PROGRAM SELF

  USE SELF_Constants
  USE SELF_SupportRoutines
  USE SELF_CLI
  USE FLAP

  IMPLICIT NONE

  TYPE(COMMAND_LINE_INTERFACE) :: selfCLI
  CHARACTER(240) :: meshFile
  CHARACTER(20) :: cqTypeChar
  CHARACTER(20) :: tqTypeChar
  CHARACTER(240) :: functionChar
  CHARACTER(240) :: derivativeChar
  CHARACTER(240) :: vectorChar(1:3)
  CHARACTER(240) :: tensorChar(1:3,1:3)
  CHARACTER(10) :: dFormChar
  CHARACTER(5) :: gpuAccelChar
  INTEGER :: dForm
  INTEGER :: cqType
  INTEGER :: tqType
  INTEGER :: cqDegree
  INTEGER :: tqDegree
  INTEGER :: nElem
  INTEGER :: nVar
  LOGICAL :: gpuAccel

  CALL Parse_CLI(selfCLI)

  CALL selfCLI % get(val=meshFile,switch='--mesh')
  CALL selfCLI % get(val=cqTypeChar,switch='--control-quadrature')
  CALL selfCLI % get(val=tqTypeChar,switch='--target-quadrature')
  CALL selfCLI % get(val=cqDegree,switch='--control-degree')
  CALL selfCLI % get(val=tqDegree,switch='--target-degree')
  CALL selfCLI % get(val=nElem,switch='--nelements')
  CALL selfCLI % get(val=nVar,switch='--nvar')
  CALL selfCLI % get(val=functionChar,switch='--function')
  CALL selfCLI % get(val=dFormChar,switch='--derivative-type')
  CALL selfCLI % get(val=derivativeChar,switch='--derivative')
  CALL selfCLI % get(val=gpuAccelChar,switch='--gpu-accel')
  CALL selfCLI % get(val=vectorChar(1),switch='--vector-x')
  CALL selfCLI % get(val=vectorChar(2),switch='--vector-y')
  CALL selfCLI % get(val=vectorChar(3),switch='--vector-z')
  CALL selfCLI % get(val=tensorChar(1,1),switch='--tensor-11')
  CALL selfCLI % get(val=tensorChar(1,2),switch='--tensor-12')
  CALL selfCLI % get(val=tensorChar(1,3),switch='--tensor-13')
  CALL selfCLI % get(val=tensorChar(2,1),switch='--tensor-21')
  CALL selfCLI % get(val=tensorChar(2,2),switch='--tensor-22')
  CALL selfCLI % get(val=tensorChar(2,3),switch='--tensor-23')
  CALL selfCLI % get(val=tensorChar(3,1),switch='--tensor-31')
  CALL selfCLI % get(val=tensorChar(3,2),switch='--tensor-32')
  CALL selfCLI % get(val=tensorChar(3,3),switch='--tensor-33')

  IF (TRIM(UpperCase(cqTypeChar)) == 'GAUSS') THEN
    cqType = GAUSS
  ELSEIF (TRIM(UpperCase(cqTypeChar)) == 'GAUSS-LOBATTO') THEN
    cqType = GAUSS_LOBATTO
  ELSE
    PRINT *, 'Invalid Control Quadrature'
    STOP -1
  END IF

  IF (TRIM(UpperCase(tqTypeChar)) == 'GAUSS') THEN
    tqType = GAUSS
  ELSEIF (TRIM(UpperCase(tqTypeChar)) == 'GAUSS-LOBATTO') THEN
    tqType = GAUSS_LOBATTO
  ELSE
    PRINT *, 'Invalid Target Quadrature'
    STOP -1
  END IF

  IF (TRIM(UpperCase(dFormChar)) == 'STRONG') THEN
    dForm = selfStrongForm
  ELSEIF (TRIM(UpperCase(dFormChar)) == 'DG') THEN
    dForm = selfWeakDGForm
  ELSEIF (TRIM(UpperCase(dFormChar)) == 'CG') THEN
    dForm = selfWeakCGForm
  ENDIF

  IF (TRIM(UpperCase(gpuAccelChar)) == 'TRUE') THEN
    gpuAccel = .TRUE.
  ELSE
    gpuAccel = .FALSE.
  END IF

  IF (selfCLI % run_command(group="blockmesh_1d")) THEN

    CALL BlockMesh1D_Test(cqType,tqType,cqDegree,tqDegree,nElem)
    

  ELSEIF (selfCLI % run_command(group="blockmesh_2d")) THEN

    CALL BlockMesh2D_Test(cqType,tqType,cqDegree,tqDegree,nElem)
    

  ELSEIF (selfCLI % run_command(group="blockmesh_3d")) THEN

    CALL BlockMesh3D_Test(cqType,tqType,cqDegree,tqDegree,nElem)
    

  ELSEIF (selfCLI % run_command(group="s1d_interp")) THEN

    CALL ScalarInterp1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                              functionChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="s2d_interp")) THEN

    CALL ScalarInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                              functionChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="s3d_interp")) THEN

    CALL ScalarInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                              functionChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="v2d_interp")) THEN

    CALL VectorInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                              vectorChar(1:2),gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="t2d_interp")) THEN

    CALL TensorInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                              tensorChar(1:2,1:2),gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="v3d_interp")) THEN

    CALL VectorInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                              vectorChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="t3d_interp")) THEN

    CALL TensorInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                              tensorChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="s1d_binterp")) THEN

    CALL ScalarBoundaryInterp1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                                      functionChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="s2d_binterp")) THEN

    CALL ScalarBoundaryInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                                      functionChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="s3d_binterp")) THEN

    CALL ScalarBoundaryInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                                      functionChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="v2d_binterp")) THEN

    CALL VectorBoundaryInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                                      vectorChar(1:2),gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="t2d_binterp")) THEN

    CALL TensorBoundaryInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                                      tensorChar(1:2,1:2),gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="v3d_binterp")) THEN

    CALL VectorBoundaryInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                                      vectorChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="t3d_binterp")) THEN

    CALL TensorBoundaryInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nVar,&
                                      tensorChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="s1d_derivative")) THEN

    CALL ScalarDerivative1D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nVar,&
                                  functionChar,derivativeChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="s2d_gradient")) THEN

    CALL ScalarGradient2D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nVar,&
                                functionChar,vectorChar(1:2),gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="s3d_gradient")) THEN

    CALL ScalarGradient3D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nVar,&
                                functionChar,vectorChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="v2d_gradient")) THEN

    CALL VectorGradient2D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nVar,&
                                vectorChar(1:2),tensorChar(1:2,1:2),gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="v2d_divergence")) THEN

    CALL VectorDivergence2D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nVar,&
                                vectorChar(1:2),functionChar,gpuAccel)
    


  ELSEIF (selfCLI % run_command(group="v3d_gradient")) THEN

    CALL VectorGradient3D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nVar,vectorChar,tensorChar,gpuAccel)
    

  ELSEIF (selfCLI % run_command(group="v3d_divergence")) THEN

    CALL VectorDivergence3D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nVar,&
                                vectorChar,functionChar,gpuAccel)
    

  END IF

  CALL selfCLI % free()

CONTAINS

  SUBROUTINE Parse_CLI(cli)
    IMPLICIT NONE
    TYPE(COMMAND_LINE_INTERFACE),INTENT(out) :: cli

    CALL cli % init(progname="self", &
                    version="v0.0.0", &
                    description="Spectral Element Libraries in Fortran (SELF)", &
                    license="ANTI-CAPITALIST SOFTWARE LICENSE (v 1.4)", &
                    authors="Joseph Schoonover (Fluid Numerics LLC)")

    CALL cli % add(switch="--mesh", &
                   switch_ab="-m", &
                   help="Path to a mesh file for control mesh."//NEW_LINE("A"), &
                   def="", &
                   required=.FALSE.)
    
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
                   choices="true,false", &
                   required=.FALSE.)

    CALL cli % add(switch="--derivative-type", &
                   switch_ab="-dtype", &
                   help="Flag to choose the type of derivative operator (weak/strong)."//NEW_LINE("A"), &
                   def="strong", &
                   choices="strong,dg", &
                   required=.FALSE.)

    CALL cli % add(switch="--function", &
                   switch_ab="-f", &
                   help="Function to interpolate from control points to target points"//NEW_LINE("A"), &
                   def="f=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--derivative", &
                   switch_ab="-df", &
                   help="Derivative of the test function; used for estimating errors."//NEW_LINE("A"), &
                   def="df=0.0", &
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
