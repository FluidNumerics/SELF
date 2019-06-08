bin_PROGRAMS=sfluid

AM_FCFLAGS = -I ../self/ $(HDF5_FFLAGS)
AM_LDFLAGS = $(HDF5_LDFLAGS)
AM_CPPFLAGS = $(HDF5_CPPLAGS)
sfluid_LDADD = ../self/libSELF.a ${METIS_LIB}

DISTCLEANFILES = fluid_equationparser_class.mod bodyforces_class.mod fluid_class.mod
include_HEADERS = fluid_equationparser_class.mod bodyforces_class.mod fluid_class.mod

sfluid_SOURCES= Fluid_EquationParser_Class.F90 BodyForces_Class.F90 Fluid_Class.F90 SELF_Fluids_Driver.F90

