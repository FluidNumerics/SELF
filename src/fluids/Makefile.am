bin_PROGRAMS=sfluid


if MPI
  sfluid_SOURCES=ModelParameters_Class.f90 MPILayer_Class.f90 Fluid_EquationParser_Class.f90 BodyForces_Class.f90 Fluid_Class.f90 SELF_Fluids_Driver.f90
else
  sfluid_SOURCES=ModelParameters_Class.f90 Fluid_EquationParser_Class.f90 BodyForces_Class.f90 Fluid_Class.f90 SELF_Fluids_Driver.f90
endif

sfluid_LDADD=../self/libself.a
