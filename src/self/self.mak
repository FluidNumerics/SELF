

CXX=hipcc
#CXXFLAGS=
OPT=-O0 -g
FC=gfortran-7
FFLAGS=-ffree-line-length-none
LIB=-L/usr/local/jsonfortran-gnu-8.0.0/lib -ljsonfortran
INCLUDE=-I/usr/loca/jsonfortran-gnu-8.0.0/lib


LINK=${INCLUDE} ${LIB}


ModelPrecision.o: ModelPrecision.F03
	${FC} ${OPT} ${FFLAGS} -c ModelPrecision.F03 -o $@

ConstantsDictionary.o: ConstantsDictionary.F03 ModelPrecision.o
	${FC} ${OPT} ${FFLAGS} -c ConstantsDictionary.F03 -o $@

CommonRoutines.o: CommonRoutines.F03 ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} ${FFLAGS} -c CommonRoutines.F03 -o $@

Quadrature.o: Quadrature.F03 ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} ${FFLAGS} -c Quadrature.F03 -o $@

Lagrange_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Lagrange_Class.F03
	${FC} ${OPT} ${FFLAGS} -c Lagrange_Class.F03 -o $@




