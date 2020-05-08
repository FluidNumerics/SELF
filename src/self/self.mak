

CXX=hipcc
CXXFLAGS=-amdgpu-target=gfx900
OPT=-O0 -g
FC=gfortran-7
FFLAGS=-DGPU -ffree-line-length-none
LIB=${HIPFORTRAN_LIB} -L/usr/local/jsonfortran-gnu-8.0.0/lib -ljsonfortran -L/usr/lib/gcc/x86_64-linux-gnu/7 -lgfortran
INCLUDE=${HIPFORTRAN_INCLUDE} -I/usr/local/jsonfortran-gnu-8.0.0/lib


LINK=${INCLUDE} ${LIB}

test: ModelPrecision.o ConstantsDictionary.o CommonRoutines.o EquationParser_Class.o Quadrature.o Lagrange_Class.o Lagrange_Test.o
	${CXX} ${CXXFLAGS} *.o ${LINK} -o $@

ModelPrecision.o: ModelPrecision.F03
	${FC} ${OPT} ${FFLAGS} -c ModelPrecision.F03 -o $@

ConstantsDictionary.o: ConstantsDictionary.F03 ModelPrecision.o
	${FC} ${OPT} ${FFLAGS} -c ConstantsDictionary.F03 -o $@

CommonRoutines.o: CommonRoutines.F03 ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} ${FFLAGS} -c CommonRoutines.F03 -o $@

EquationParser_Class.o: EquationParser_Class.F03 ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} ${FFLAGS} -c EquationParser_Class.F03 -o $@

Quadrature.o: Quadrature.F03 ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} ${FFLAGS} -c Quadrature.F03 -o $@

Lagrange_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Lagrange_Class.F03
	${FC} ${OPT} ${FFLAGS} ${HIPFORTRAN_INCLUDE} -c Lagrange_Class.F03 -o $@

Lagrange_Test.o : Lagrange_Test.F03 ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Lagrange_Class.o EquationParser_Class.o
	${FC} ${OPT} ${FFLAGS} ${LINK} -c Lagrange_Test.F03 -o $@



