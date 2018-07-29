AC_INIT([self-fluids],[1.2],[joe@fluidnumerics.com],[https://gitlab.com/fluidnumerics/SELF-Fluids])
AM_INIT_AUTOMAKE(-Wall -Werror foreign)

# Set up build system to generate libSELF.a
AC_PROG_RANLIB
AM_PROG_AR

AC_CONFIG_MACRO_DIRS([m4])

# Find Fortran compilers
AC_PROG_FC([pgfortran gfortran])
if test "$FC"="pgfortran"; then
  FCFLAGS="-Mpreprocess"
fi

if test "$FC"="gfortran"; then
  FCFLAGS="-cpp -ffree-line-length-none"
fi


# Check for MPI
AC_ARG_ENABLE( mpi,
	      AS_HELP_STRING([--enable-mpi, default:no]),
[case "${enableval}" in
      yes) mpi=true ;;
      no) mpi=false ;;
      *) AC_MSG_ERROR([bad value ${enableval} for --enable-mpi]) ;;
esac],
[mpi=false] )
AM_CONDITIONAL(MPI, test x"$mpi"=x"true")

if test "$mpi"="true"; then
  AX_MPI(mpi=true, mpi=false)
  if test "$mpi"="true"; then
    FC=$MPIFC
    FCFLAGS+="-DHAVE_MPI"
  else
    AC_MSG_ERROR('--enable-mpi set to yes but MPI-Fortran compiler not found.')
  fi
fi

# Check for Double Precision floating point arithmetic request.
AC_ARG_ENABLE( double_precision,
	      AS_HELP_STRING([--enable-double-precision, default:no]),
[case "${enableval}" in
      yes) double_precision=true ;;
      no) double_precision=false ;;
      *) AC_MSG_ERROR([bad value ${enableval} for --enable-double-precision]) ;;
esac],
[double_precision=false] )
AM_CONDITIONAL(DOUBLE_PRECISION, test x"$double_precision"=x"true")

if test "$double_precision"="true"; then
  FCFLAGS+="-DDOUBLE_PRECISION"
fi


# Check for CUDA-Fortran acceleration request
# Set --enable-cuda-fortran yes turns on the -DHAVE_CUDA CPP-flag if the PGI compilers are used.
AC_ARG_ENABLE( cuda,
AS_HELP_STRING([--enable-cuda-fortran, default:no]),
[case "${enableval}" in
      yes) cuda=true ;;
      no) cuda=false ;;
      *) AC_MSG_ERROR([bad value ${enableval} for --enable-cuda-fortran]) ;;
esac],
[cuda=false] )
AM_CONDITIONAL(CUDA, test x"$cuda"=x"true")

if test "$cuda"="true"; then
  if test "$FC"="pgfortran"; then
    # GPU_ARCH needs to be set as an environment variable or passed, e.g.
    # GPU_ARCH=cc60 ./configure --enable-cuda
    FCFLAGS+="-DHAVE_CUDA -Mcuda=$GPU_ARCH"
  else
    AC_MSG_ERROR("CUDA-Fortran can only be enabled with PGI compilers")
  fi
fi

# Check for request to enable timers in Fluids_Class.f90
AC_ARG_ENABLE( timing,
AS_HELP_STRING([--enable-timing, default:no]),
[case "${enableval}" in
      yes) timing=true ;;
      no) timing=false ;;
      *) AC_MSG_ERROR([bad value ${enableval} for --enable-timing]) ;;
esac],
[timing=false] )
AM_CONDITIONAL(TIMING, test x"$timing"=x"true")

if test "$timing"="true"; then
  FCFLAGS+="-DTIMING"
fi

AC_CONFIG_FILES([
Makefile
src/self/Makefile
src/fluids/Makefile
])

AC_OUTPUT
