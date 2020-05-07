
#FC=gfortran-8
#OPT=-O0 -g
LIB=-L/usr/local/jsonfortran/lib -ljsonfortran -L/usr/local/self-fluids/lib -lSELF
INC=-I/usr/local/jsonfortran/lib -I/usr/local/self-fluids/include


lagrange_test : Lagrange_Test.F90
	${FC} ${OPT} Lagrange_Test.F90 ${INC} ${LIB} -o $@
