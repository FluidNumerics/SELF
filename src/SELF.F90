PROGRAM SELF

USE SELF_SupportRoutines
USE SELF_Tests
USE FLAP

IMPLICIT NONE

  TYPE(COMMAND_LINE_INTERFACE) :: self_cli
  CHARACTER(20) :: cqTypeChar
  CHARACTER(20) :: tqTypeChar
  INTEGER :: cqType
  INTEGER :: tqType
  INTEGER :: nControlPoints
  INTEGER :: nTargetPoints
  INTEGER :: nElem
  INTEGER :: error


    CALL Parse_CLI(self_cli)

    CALL self_cli % get( val = cqTypeChar, switch = '--control-quadrature' )
    CALL self_cli % get( val = tqTypeChar, switch = '--target-quadrature' )
    CALL self_cli % get( val = nControlPoints, switch = '--control-points' )
    CALL self_cli % get( val = nTargetPoints, switch = '--target-points' )
    CALL self_cli % get( val = nElem, switch = '--nelements' )

    IF( TRIM(UpperCase(cqTypeChar)) == 'GAUSS' )THEN
      cqType = GAUSS
    ELSEIF( TRIM(UpperCase(cqTypeChar)) == 'GAUSS-LOBATTO' )THEN
      cqType = GAUSS_LOBATTO
    ELSE
      PRINT*, 'Invalid Control Quadrature'
      STOP 1
    ENDIF

    IF( TRIM(UpperCase(tqTypeChar)) == 'GAUSS' )THEN
      tqType = GAUSS
    ELSEIF( TRIM(UpperCase(tqTypeChar)) == 'GAUSS-LOBATTO' )THEN
      tqType = GAUSS_LOBATTO
    ELSE
      PRINT*, 'Invalid Target Quadrature'
      STOP 1
    ENDIF

    IF( self_cli % run_command( group = "blockmesh_1d" ) )THEN

      CALL BlockMesh1D_Test( cqType, tqType, nControlPoints, nTargetPoints, nElem, error )

    ENDIF

    CALL self_cli % free()

CONTAINS

  SUBROUTINE Parse_CLI(cli)
    IMPLICIT NONE
    TYPE(COMMAND_LINE_INTERFACE), INTENT(out) :: cli

    CALL cli % init( progname = "self", &
                     version = "v0.0.0", &
                     description = "Spectral Element Libraries in Fortran (SELF)", &
                     license = "ANTI-CAPITALIST SOFTWARE LICENSE (v 1.4)", &
                     authors = "Joseph Schoonover (Fluid Numerics LLC)")

    CALL cli % add( switch = "--control-points", &
                    switch_ab = "-c", &
                    help = "The number of control points on computational grid."//NEW_LINE("A"),&
                    def = "2", &
                    required = .FALSE. )

    CALL cli % add( switch = "--target-points", &
                    switch_ab = "-t", &
                    help = "The number of target points for grid interpolation methods."//NEW_LINE("A"),&
                    def = "7", &
                    required = .FALSE. )

    CALL cli % add( switch = "--control-quadrature", &
                    switch_ab = "-cq", &
                    def = "gauss", &
                    help = "The quadrature type for control points."//NEW_LINE("A"),&
                    choices = "gauss, gauss-lobatto, uniform", &
                    required = .FALSE. )

    CALL cli % add( switch = "--target-quadrature", &
                    switch_ab = "-tq", &
                    def = "gauss", &
                    help = "The quadrature type for target points."//NEW_LINE("A"),&
                    choices = "gauss, gauss-lobatto, uniform", &
                    required = .FALSE. )

    CALL cli % add( switch = "--nvar", &
                    switch_ab = "-nv", &
                    help = "The number of functions used to simulate increased workloads."//NEW_LINE("A"),&
                    def = "5", &
                    required = .FALSE. )

    CALL cli % add( switch = "--nelements", &
                    switch_ab = "-ne", &
                    help = "The number of elements to use on the test mesh. "//&
                           "In multi-dimensions, the number of elements in each direction."//NEW_LINE("A"),&
                    def = "10", &
                    required = .FALSE. )

    CALL cli % add( switch = "--gpu-accel", &
                    switch_ab = "-gpu", &
                    help = "Boolean flag for enabling or disabling GPU acceleration in tests."//NEW_LINE("A"),&
                    def = "false", &
                    choices = "true, false", &
                    required = .FALSE. )

    CALL cli % add_group( group = "blockmesh_1d", &
                          description = "Block Mesh generation in 1D" )

    CALL cli % add_group( group = "blockmesh_2d", &
                          description = "Block Mesh generation in 2D" )

    CALL cli % add_group( group = "blockmesh_3d", &
                          description = "Block Mesh generation in 3D" )

    CALL cli % add_group( group = "s1d_interp", &
                          description = "Scalar 1D Interpolation" )

    CALL cli % add( switch = "--function", &
                    switch_ab = "-f", &
                    group = "s1d_interp", &
                    help = "Function to interpolate from control points to target points", &
                    def = "f=1.0", &
                    required = .FALSE. )

    CALL cli % add_group( group = "s2d_interp", &
                          description = "Scalar 2D Interpolation" )

    CALL cli % add( switch = "--function", &
                    switch_ab = "-f", &
                    group = "s2d_interp", &
                    help = "Function to interpolate from control points to target points", &
                    def = "f=1.0", &
                    required = .FALSE. )

    CALL cli % add_group( group = "s3d_interp", &
                          description = "Scalar 3D Interpolation" )

    CALL cli % add( switch = "--function", &
                    switch_ab = "-f", &
                    group = "s3d_interp", &
                    help = "Function to interpolate from control points to target points", &
                    def = "f=1.0", &
                    required = .FALSE. )

    CALL cli % add_group( group = "v2d_interp", &
                          description = "Vector 2D Interpolation" )

    CALL cli % add( switch = "--vector-x", &
                    switch_ab = "-vx", &
                    group = "v2d_interp", &
                    help = "x-component of the vector function to interpolate from control points to target points", &
                    def = "vx=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--vector-y", &
                    switch_ab = "-vy", &
                    group = "v2d_interp", &
                    help = "y-component of the vector function to interpolate from control points to target points", &
                    def = "vy=1.0", &
                    required = .FALSE. )

    CALL cli % add_group( group = "v3d_interp", &
                          description = "Vector 3D Interpolation" )

    CALL cli % add( switch = "--vector-x", &
                    switch_ab = "-vx", &
                    group = "v3d_interp", &
                    help = "x-component of the vector function to interpolate from control points to target points", &
                    def = "vx=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--vector-y", &
                    switch_ab = "-vy", &
                    group = "v3d_interp", &
                    help = "y-component of the vector function to interpolate from control points to target points", &
                    def = "vy=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--vector-z", &
                    switch_ab = "-vz", &
                    group = "v3d_interp", &
                    help = "z-component of the vector function to interpolate from control points to target points", &
                    def = "vz=1.0", &
                    required = .FALSE. )

    CALL cli % add_group( group = "t2d_interp", &
                          description = "Tensor 2D Interpolation" )

    CALL cli % add( switch = "--tensor-11", &
                    switch_ab = "-t11", &
                    group = "t2d_interp", &
                    help = "Row 1 column 1 of the tensor function to interpolate from control points to target points", &
                    def = "t11=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--tensor-12", &
                    switch_ab = "-t12", &
                    group = "t2d_interp", &
                    help = "Row 1 column 2 of the tensor function to interpolate from control points to target points", &
                    def = "t12=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--tensor-21", &
                    switch_ab = "-t21", &
                    group = "t2d_interp", &
                    help = "Row 2 column 1 of the tensor function to interpolate from control points to target points", &
                    def = "t21=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--tensor-22", &
                    switch_ab = "-t22", &
                    group = "t2d_interp", &
                    help = "Row 2 column 2 of the tensor function to interpolate from control points to target points", &
                    def = "t22=1.0", &
                    required = .FALSE. )

    CALL cli % add_group( group = "t3d_interp", &
                          description = "Tensor 3D Interpolation" )

    CALL cli % add( switch = "--tensor-11", &
                    switch_ab = "-t11", &
                    group = "t3d_interp", &
                    help = "Row 1 column 1 of the tensor function to interpolate from control points to target points", &
                    def = "t11=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--tensor-12", &
                    switch_ab = "-t12", &
                    group = "t3d_interp", &
                    help = "Row 1 column 2 of the tensor function to interpolate from control points to target points", &
                    def = "t12=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--tensor-13", &
                    switch_ab = "-t13", &
                    group = "t3d_interp", &
                    help = "Row 1 column 3 of the tensor function to interpolate from control points to target points", &
                    def = "t13=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--tensor-21", &
                    switch_ab = "-t21", &
                    group = "t3d_interp", &
                    help = "Row 2 column 1 of the tensor function to interpolate from control points to target points", &
                    def = "t21=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--tensor-22", &
                    switch_ab = "-t22", &
                    group = "t3d_interp", &
                    help = "Row 2 column 2 of the tensor function to interpolate from control points to target points", &
                    def = "t22=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--tensor-23", &
                    switch_ab = "-t23", &
                    group = "t3d_interp", &
                    help = "Row 2 column 3 of the tensor function to interpolate from control points to target points", &
                    def = "t23=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--tensor-31", &
                    switch_ab = "-t31", &
                    group = "t3d_interp", &
                    help = "Row 3 column 1 of the tensor function to interpolate from control points to target points", &
                    def = "t31=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--tensor-32", &
                    switch_ab = "-t32", &
                    group = "t3d_interp", &
                    help = "Row 3 column 2 of the tensor function to interpolate from control points to target points", &
                    def = "t32=1.0", &
                    required = .FALSE. )

    CALL cli % add( switch = "--tensor-33", &
                    switch_ab = "-t33", &
                    group = "t3d_interp", &
                    help = "Row 2 column 3 of the tensor function to interpolate from control points to target points", &
                    def = "t33=1.0", &
                    required = .FALSE. )

    CALL cli % parse()

  END SUBROUTINE Parse_CLI

END PROGRAM SELF
