

CXX=hipcc
#CXXFLAGS=-amdgpu-target=gfx900
CXXFLAGS=
OPT=-O0 -g
#FC=gfortran-7
FC=gfortran-8
#FFLAGS=-ffree-line-length-none
FFLAGS=-DGPU -ffree-line-length-none
#LIB=${HIPFORTRAN_LIB} -L/usr/local/jsonfortran-gnu-8.0.0/lib -ljsonfortran -L/usr/lib/gcc/x86_64-linux-gnu/7 -lgfortran -lstdc++
LIB=${HIPFORTRAN_LIB} -L/apps/json-fortran/jsonfortran-gnu-7.1.0/lib -ljsonfortran -L/usr/lib/gcc/x86_64-linux-gnu/8 -lgfortran -lstdc++
#INCLUDE=${HIPFORTRAN_INCLUDE} -I/usr/local/jsonfortran-gnu-8.0.0/lib
INCLUDE=${HIPFORTRAN_INCLUDE} -I/apps/json-fortran/jsonfortran-gnu-7.1.0/lib


LINK=${INCLUDE} ${LIB}

test: ModelPrecision.o ConstantsDictionary.o CommonRoutines.o EquationParser_Class.o Quadrature.o Lagrange_Class.o Lagrange_Class_Tests.o Lagrange_Test.o
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

Lagrange_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Quadrature.o Lagrange_HIP.o Lagrange_Class.F03
	${FC} ${OPT} ${FFLAGS} ${HIPFORTRAN_INCLUDE} -c Lagrange_Class.F03 -o $@

Lagrange_HIP.o : Lagrange_HIP.cpp SELF_Macros.h
	${CXX} ${OPT} ${CXXFLAGS} -c ${HIPFORTRAN_INCLUDE} Lagrange_HIP.cpp -o $@

Lagrange_Test.o : Lagrange_Test.F03 ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Lagrange_Class.o Lagrange_Class_Tests.o EquationParser_Class.o Lagrange_HIP.o
	${FC} ${OPT} ${FFLAGS} ${LINK} -c Lagrange_Test.F03 -o $@

Lagrange_Class_Tests.o : Lagrange_Class_Tests.F03 ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Lagrange_Class.o EquationParser_Class.o Lagrange_HIP.o
	${FC} ${OPT} ${FFLAGS} ${LINK} -c Lagrange_Class_Tests.F03 -o $@



