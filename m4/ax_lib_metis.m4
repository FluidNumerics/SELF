#
# SYNOPSIS
#
#   AX_LIB_METIS()
#
# DESCRIPTION
#
#   This macro provides tests of the availability of METIS
#
#   The macro adds a --with-hdf5 option accepting one of three values:
#
#     no   - do not check for the HDF5 library.
#     yes  - do check for HDF5 library in standard locations.
#     path - complete path to the HDF5 helper script h5cc or h5pcc.
#
#   If METIS is successfully found, this macro calls
#
#     AC_SUBST(METIS_LIBS)
#
#   and sets with_metis="yes".  Additionally, the macro sets
#
#   If METIS is disabled or not found, this macros sets with_metis="no".
#
#   To use the macro, one would code one of the following in "configure.ac"
#   before AC_OUTPUT:
#
#     1) dnl Check for HDF5 support
#        AX_LIB_METIS()
#
# LICENSE
#

AC_DEFUN([AX_LIB_METIS], [


dnl Add a default --with-hdf5 configuration option.
AC_ARG_WITH([metis],
  AS_HELP_STRING(
    [--with-metis=[PATH]],
    m4_case(m4_normalize([$1]),
            [location of libmetis.a])
  ),
  [if test "$withval" = ""; then
     with_metis="no"
   else
     with_metis="yes"
     METIS="$withval"
   fi],
   [with_metis="yes"]
)

dnl Set defaults to blank
METIS_LIB=""

dnl Try and find hdf5 compiler tools and options.
if test "$with_metis" = "yes"; then
  METIS_LIB=$METIS
  AC_SUBST([METIS_LIB])
  AC_DEFINE([HAVE_METIS], [1], [Defined if you have METIS])
fi
])
