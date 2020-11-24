PROGRAM SELF


USE FLAP

IMPLICIT NONE

  TYPE(COMMAND_LINE_INTERFACE) :: cli

    CALL cli % init( progname = "SELF", &
                     version = "v0.0.0", &
                     description = "Spectral Element Libraries in Fortran (SELF)", &
                     license = "ANTI-CAPITALIST SOFTWARE LICENSE (v 1.4)", &
                     authors = "Joseph Schoonover (Fluid Numerics LLC)")

    CALL cli % add( switch = "--control-points", &
                    switch_ab = "-c", &
                    def = "2", &
                    required = .FALSE. )

    CALL cli % add( switch = "--target-points", &
                    switch_ab = "-t", &
                    def = "7", &
                    required = .FALSE. )

    CALL cli % add( switch = "--quadrature-type", &
                    switch_ab = "-q", &
                    def = "gauss", &
                    choices = "gauss, gauss-lobatto", &
                    required = .FALSE. )

    CALL cli % add( switch = "--nvar", &
                    switch_ab = "-nv", &
                    def = "5", &
                    required = .FALSE. )

    CALL cli % add( switch = "--nelements", &
                    switch_ab = "-ne", &
                    def = "10", &
                    required = .FALSE. )

    CALL cli % add_group( group = "s1d_interp", &
                          description = "Scalar 1D Interpolation" )

    CALL cli % add_group( group = "s2d_interp", &
                          description = "Scalar 2D Interpolation" )

    CALL cli % add_group( group = "s3d_interp", &
                          description = "Scalar 3D Interpolation" )

    CALL cli % parse()
    CALL cli % free()
END PROGRAM SELF
