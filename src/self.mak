FC=hipfc
GFLAGS=-v -g
FFLAGS=-v -g -DGPU -ffree-line-length-none -O0 -I/opt/json-fortran/jsonfortran-gnu-8.1.0/lib -I/opt/feqparse/include
FLIBS=-L/opt/json-fortran/jsonfortran-gnu-8.1.0/lib -ljsonfortran -L/opt/feqparse/lib -lfeqparse

test: SELFConstants.o CommonRoutines.o Quadrature.o Lagrange_Class.o Lagrange_Class_Tests.o Lagrange_Test.o SysConf.o
	${FC} ${FLIBS} ${FFLAGS} *.o -o $@

#SELFPrecision.o: SELFPrecision.F90
#	${FC} ${FFLAGS} -c SELFPrecision.F90 -o $@

SELFConstants.o: SELFConstants.F90
	${FC} ${FFLAGS} -c SELFConstants.F90 -o $@

SELFMemory.o: SELFMemory.F90 SELFConstants.o
	${FC} ${FFLAGS} -c SELFMemory.F90 -o $@

SysConf.o: SysConf.F90
	${FC} ${FFLAGS} -c SysConf.F90 -o $@

CommonRoutines.o: CommonRoutines.F90 SELFConstants.o
	${FC} ${FFLAGS} -c CommonRoutines.F90 -o $@

Quadrature.o: Quadrature.F90 SELFConstants.o
	${FC} ${FFLAGS} -c Quadrature.F90 -o $@

Lagrange_Class.o : SELFConstants.o SELFMemory.o CommonRoutines.o Quadrature.o Lagrange_HIP.o Lagrange_Class.F90
	${FC} ${FFLAGS} -c Lagrange_Class.F90 -o $@

Lagrange_HIP.o : Lagrange_HIP.cpp SELF_Macros.h
	${FC} ${GFLAGS} -c Lagrange_HIP.cpp -o $@

Lagrange_Test.o : Lagrange_Test.F90 SELFConstants.o CommonRoutines.o Lagrange_Class.o Lagrange_Class_Tests.o Lagrange_HIP.o
	${FC} ${FFLAGS} -c Lagrange_Test.F90 -o $@

Lagrange_Class_Tests.o : Lagrange_Class_Tests.F90 SELFConstants.o CommonRoutines.o Lagrange_Class.o Lagrange_HIP.o SysConf.o
	${FC} ${FFLAGS} -c Lagrange_Class_Tests.F90 -o $@

NodalSEMData.o: NodalSEMData.F90 SELFConstants.o Lagrange_Class.o
	${FC} ${FFLAGS} -c NodalSEMData.F90 -o $@

SEMMesh.o: SEMMesh.F90 SELFConstants.o Lagrange_Class.o NodalSEMData.o
	${FC} ${FFLAGS} -c SEMMesh.F90 -o $@

MappedNodalSEMData.o: MappedNodalSEMData.F90 SEMMesh.o SELFConstants.o Lagrange_Class.o NodalSEMData.o
	${FC} ${FFLAGS} -c MappedNodalSEMData.F90 -o $@





