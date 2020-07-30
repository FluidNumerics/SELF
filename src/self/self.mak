

FC=hipfc
GFLAGS=-v -g --fmad=false
FFLAGS=-v -g -DGPU -ffree-line-length-none -O0 -I/usr/local/jsonfortran-gnu-7.1.0/lib -I/opt/feqparse/include
FLIBS=-L/usr/local/jsonfortran-gnu-7.1.0/lib -ljsonfortran -L/opt/feqparse/lib -lfeqparse

test: ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Quadrature.o Lagrange_Class.o Lagrange_Class_Tests.o Lagrange_Test.o SysConf.o
	${FC} ${FLIBS} ${FFLAGS} *.o -o $@

ModelPrecision.o: ModelPrecision.F90
	${FC} ${FFLAGS} -c ModelPrecision.F90 -o $@

SysConf.o: SysConf.F90
	${FC} ${FFLAGS} -c SysConf.F90 -o $@

ConstantsDictionary.o: ConstantsDictionary.F90 ModelPrecision.o
	${FC} ${FFLAGS} -c ConstantsDictionary.F90 -o $@

CommonRoutines.o: CommonRoutines.F90 ModelPrecision.o ConstantsDictionary.o
	${FC} ${FFLAGS} -c CommonRoutines.F90 -o $@

Quadrature.o: Quadrature.F90 ModelPrecision.o ConstantsDictionary.o
	${FC} ${FFLAGS} -c Quadrature.F90 -o $@

Lagrange_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Quadrature.o Lagrange_HIP.o Lagrange_Class.F90
	${FC} ${FFLAGS} -c Lagrange_Class.F90 -o $@

Lagrange_HIP.o : Lagrange_HIP.cpp SELF_Macros.h
	${FC} ${GFLAGS} -c Lagrange_HIP.cpp -o $@

Lagrange_Test.o : Lagrange_Test.F90 ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Lagrange_Class.o Lagrange_Class_Tests.o Lagrange_HIP.o
	${FC} ${FFLAGS} -c Lagrange_Test.F90 -o $@

Lagrange_Class_Tests.o : Lagrange_Class_Tests.F90 ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Lagrange_Class.o Lagrange_HIP.o SysConf.o
	${FC} ${FFLAGS} -c Lagrange_Class_Tests.F90 -o $@



