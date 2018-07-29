AC_INIT([self-fluids],[1.2],[joe@fluidnumerics.com])
AM_INIT_AUTOMAKE(-Wall -Werror foreign)

# Set up build system to generate libSELF.a
AC_PROG_RANLIB
AM_PROG_AR

AC_CONFIG_MACRO_DIRS([m4])

# Find Fortran compilers
AC_PROG_FC([pgfortran gfortran])
AS_IF([test "x$FC" = "xpgfortran"],[
  FCFLAGS+=" -Mpreprocess"
])

AS_IF([test "x$FC" = "xgfortran"],[
  FCFLAGS+=" -cpp -ffree-line-length-none"
])

# Check for CUDA-Fortran acceleration request
# Set --enable-cuda-fortran yes turns on the -DHAVE_CUDA CPP-flag if the PGI compilers are used.
AC_ARG_ENABLE( cuda, AS_HELP_STRING([--enable-cuda-fortran, default:no]))
AM_CONDITIONAL(CUDA, test x"$cuda"=x"true")
AS_IF([test "x$cuda" = "xyes"],[

  if test "x$FC" = "xpgfortran"; then
    # GPU_ARCH needs to be set as an environment variable or passed, e.g.
    # GPU_ARCH=cc60 ./configure --enable-cuda
    FCFLAGS+=" -DHAVE_CUDA -Mcuda=$GPU_ARCH"
  else
    AC_MSG_ERROR("CUDA-Fortran can only be enabled with PGI compilers")
  fi

])

# Check for MPI
AC_ARG_ENABLE([mpi], AS_HELP_STRING([--enable-mpi], [default:no]))
AM_CONDITIONAL(MPI, test x"$mpi"=x"true")
AS_IF([test "x$enable_mpi" = "xyes"],[

  AX_MPI

  if test "$MPIFC"="mpif90"; then
    AC_MSG_NOTICE("Found $MPIFC")
    FC="$MPIFC"
    FCFLAGS+=" -DHAVE_MPI"
  else
    AC_MSG_ERROR('--enable-mpi set to yes but MPI-Fortran compiler not found.')
  fi

])

# Check for Double Precision floating point arithmetic request.
AC_ARG_ENABLE( double_precision, AS_HELP_STRING([--enable-double-precision, default:no]))
AM_CONDITIONAL(DOUBLE_PRECISION, test x"$double_precision"=x"true")
AS_IF([test "x$double_precision" = "xyes"],[
  FCFLAGS+=" -DDOUBLE_PRECISION"
])


# Check for request to enable timers in Fluids_Class.f90
AC_ARG_ENABLE( timing, AS_HELP_STRING([--enable-timing, default:no]))
AM_CONDITIONAL(TIMING, test x"$timing"=x"true")

AS_IF( [test "x$timing" = "xtrue"],[
  FCFLAGS+=" -DTIMING"
])

AC_CONFIG_FILES([
Makefile
src/self/Makefile
src/fluids/Makefile
])

AC_OUTPUT
