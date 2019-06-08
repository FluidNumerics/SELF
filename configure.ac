AC_INIT([self-fluids],[1.2],[joe@fluidnumerics.com])
AM_INIT_AUTOMAKE(-Wall -Werror foreign)

# Set up build system to generate libSELF.a
AC_PROG_RANLIB
AM_PROG_AR

AC_CONFIG_MACRO_DIRS([m4])

# Find Fortran compilers
AC_PROG_FC([pgfortran gfortran])

AC_LANG_PUSH([Fortran])

# Check for CUDA-Fortran acceleration request
# Set --enable-cuda-fortran yes turns on the -DHAVE_CUDA CPP-flag if the PGI compilers are used.
AC_ARG_ENABLE( cuda, AS_HELP_STRING([--enable-cuda, default:no]))
AS_IF([test "x$enable_cuda" = "xyes"],[

  FCFLAGS+=" -DHAVE_CUDA -Mcuda=ptxinfo"

])

# Check for OpenMP
AC_ARG_ENABLE([openmp], AS_HELP_STRING([--enable-openmp], [default:no]))
AS_IF([test "x$enable_openmp" = "xyes"],[

  FCFLAGS+=" -DHAVE_OPENMP"
  AX_CHECK_COMPILE_FLAG([-mp], [FCFLAGS+=" -mp"])
  AX_CHECK_COMPILE_FLAG([-openmp], [FCFLAGS+=" -openmp"])
  AX_CHECK_COMPILE_FLAG([-fopenmp], [FCFLAGS+=" -fopenmp"])

])

AX_LIB_METIS()

# Check for HDF5
AX_LIB_HDF5()
AS_IF([test "x$with_hdf5" = "xyes"],[
  AC_MSG_NOTICE([found hdf5 $H5FC])
  FC="$H5FC"
],[
  AC_MSG_ERROR([Unable to find HDF5, we need HDF5.])
])

# Check for MPI
AC_CHECK_PROGS(MPIFC, mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c, $FC)
AC_ARG_ENABLE([mpi], AS_HELP_STRING([--enable-mpi], [default:no]))
AS_IF([test "x$enable_mpi" = "xyes"],[

  FCFLAGS+=" -DHAVE_MPI"
  AS_IF([test "x$HDF5_TYPE" = "xserial",],[
    AC_MSG_ERROR([HDF5 found is serial. We need parallel HDF5.])
  ])

])


# Check for potential temperature formulation request.
AC_ARG_ENABLE( pottemp, AS_HELP_STRING([--enable-pottemp, default:no]))
AS_IF([test "x$enable_pottemp" = "xyes"],[
  FCFLAGS+=" -DPOTENTIAL_TEMPERATURE"
])

# Check for Double Precision floating point arithmetic request.
AC_ARG_ENABLE( double, AS_HELP_STRING([--enable-double, default:no]))
AS_IF([test "x$enable_double" = "xyes"],[
  FCFLAGS+=" -DDOUBLE_PRECISION"
])


# Check for request to enable timers in Fluids_Class.f90
AC_ARG_ENABLE( timing, AS_HELP_STRING([--enable-timing, default:no]))

AS_IF( [test "x$enable_timing" = "xyes"],[
  FCFLAGS+=" -DTIMING"
])

# Check for request to enable diagnostics in Fluids_Class.f90
AC_ARG_ENABLE( diagnostics, AS_HELP_STRING([--enable-diagnostics, default:no]))

AS_IF( [test "x$enable_diagnostics" = "xyes"],[
  FCFLAGS+=" -DDIAGNOSTICS"
])

# Check for request to enable tecplot output in Fluids_Class.f90
AC_ARG_ENABLE( tecplot, AS_HELP_STRING([--enable-tecplot, default:no]))
AS_IF( [test "x$enable_tecplot" = "xyes"],[
  FCFLAGS+=" -DTECPLOT"
])

AC_ARG_ENABLE( unittest, AS_HELP_STRING([--enable-unittest, default:no]))
AS_IF( [test "x$enable_unittest" = "xyes"],[
  FCFLAGS+=" -DUNIT_TEST"
])


AC_CONFIG_FILES([
Makefile
src/self/Makefile
src/fluids/Makefile
src/cubed-sphere/Makefile
])

AC_OUTPUT
