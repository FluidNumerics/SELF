FC=hipfc
GFLAGS=-v -g
FFLAGS=-v -g -DGPU -ffree-line-length-none -O0 -I/opt/json-fortran/lib -I/opt/feqparse/include
FLIBS=-L/opt/json-fortran/lib -ljsonfortran -L/opt/feqparse/lib -lfeqparse

test: SELF_Constants.o SELF_SupportRoutines.o SELF_Quadrature.o SELF_Lagrange.o SELF_Lagrange_Tests.o Lagrange_Test.o SysConf.o
	${FC} ${FLIBS} ${FFLAGS} *.o -o $@

#SELFPrecision.o: SELFPrecision.F90
#	${FC} ${FFLAGS} -c SELFPrecision.F90 -o $@

SELF_Constants.o: SELF_Constants.F90
	${FC} ${FFLAGS} -c SELF_Constants.F90 -o $@

SELF_Memory.o: SELF_Memory.F90 SELF_Constants.o
	${FC} ${FFLAGS} -c SELF_Memory.F90 -o $@

SysConf.o: SysConf.F90
	${FC} ${FFLAGS} -c SysConf.F90 -o $@

SELF_SupportRoutines.o: SELF_SupportRoutines.F90 SELF_Constants.o
	${FC} ${FFLAGS} -c SELF_SupportRoutines.F90 -o $@

SELF_Quadrature.o: SELF_Quadrature.F90 SELF_Constants.o
	${FC} ${FFLAGS} -c SELF_Quadrature.F90 -o $@

SELF_Lagrange.o : SELF_Constants.o SELF_Memory.o SELF_SupportRoutines.o SELF_Quadrature.o SELF_Lagrange_HIP.o SELF_Lagrange.F90
	${FC} ${FFLAGS} -c SELF_Lagrange.F90 -o $@

SELF_Lagrange_HIP.o : SELF_Lagrange_HIP.cpp SELF_Macros.h
	${FC} ${GFLAGS} -c SELF_Lagrange_HIP.cpp -o $@

Lagrange_Test.o : Lagrange_Test.F90 SELF_Constants.o SELF_SupportRoutines.o SELF_Lagrange.o SELF_Lagrange_Tests.o SELF_Lagrange_HIP.o
	${FC} ${FFLAGS} -c Lagrange_Test.F90 -o $@

SELF_Lagrange_Tests.o : SELF_Lagrange_Tests.F90 SELF_Constants.o SELF_SupportRoutines.o SELF_Lagrange.o SELF_Lagrange_HIP.o SysConf.o
	${FC} ${FFLAGS} -c SELF_Lagrange_Tests.F90 -o $@

SELF_Data.o: SELF_Data.F90 SELF_Constants.o SELF_Lagrange.o
	${FC} ${FFLAGS} -c SELF_Data.F90 -o $@

SELF_Mesh.o: SELF_Mesh.F90 SELF_Constants.o SELF_Lagrange.o SELF_Data.o
	${FC} ${FFLAGS} -c SELF_Mesh.F90 -o $@

SELF_MappedData_HIP.o: SELF_MappedData_HIP.cpp SELF_Mesh.o SELF_Macros.h
	${FC} ${GFLAGS} -c SELF_MappedData_HIP.cpp -o $@

SELF_MappedData.o: SELF_MappedData.F90 SELF_Mesh.o SELF_Constants.o SELF_Lagrange.o SELF_Data.o
	${FC} ${FFLAGS} -c SELF_MappedData.F90 -o $@





