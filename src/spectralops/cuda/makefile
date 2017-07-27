# makefile
#
#
#  Directives :
#
#    

SHELL=/bin/bash
FC=pgf90
LIB=
INC=
OPT=-fast -Mcuda -acc -ta=tesla:cc60 -Minfo=accel
WFLAGS= 
SELFDIR=~/Software/SELF-v3.0/



# determine the main directory for this experiment
DIR=$(shell pwd | rev | cut -d "/" -f2- | rev )
INIDIR=${DIR}/
SRCDIR=${SELFDIR}src/



TLOBJ = ModelPrecision.o \
        ConstantsDictionary.o \
        CommonRoutines.o \
        Timing.o \
        Quadrature.o \
        InterpolationSupportRoutines.o \
        Lagrange_Class.o \
        Lagrange_Cuda_Class.o
        
TNOBJ = ModelPrecision.o \
        ConstantsDictionary.o \
        CommonRoutines.o \
        Timing.o \
        Quadrature.o \
        InterpolationSupportRoutines.o \
        Lagrange_Class.o \
        Lagrange_Cuda_Class.o \
        NodalStorage_Cuda_Class.o \

TestLagrange_CUDA : ${TLOBJ} TestLagrangeClass_CUDA.o
	${FC} ${OPT} ${TLOBJ} TestLagrangeClass_CUDA.o -o $@

TestNodalStorage_CUDA : ${TNOBJ} TestNodalStorage_CUDA.o
	${FC} ${OPT} ${TNOBJ} TestNodalStorage_CUDA.o -o $@


# ---------------------------------------------------------------------- #
TestLagrangeClass_CUDA.o : Lagrange_Class.o Lagrange_Cuda_Class.o
	${FC} ${OPT} -c TestLagrangeClass_CUDA.cuf -o $@

TestNodalStorage_CUDA.o : Timing.o NodalStorage_Cuda_Class.o
	${FC} ${OPT} -c TestNodalStorage_CUDA.cuf -o $@

# ---------------------------------------------------------------------- #
#                           ~/src/common/                                #
# ---------------------------------------------------------------------- #

# Modules
ModelPrecision.o :
	${FC} ${OPT} -c ${SRCDIR}common/ModelPrecision.f90 -o $@

ConstantsDictionary.o : ModelPrecision.o
	${FC} ${OPT} -c ${SRCDIR}common/ConstantsDictionary.f90 -o $@

LinkedList_Class.o :
	${FC} ${OPT} -c ${SRCDIR}common/LinkedList_Class.f90 -o $@

KeyRing_Class.o : CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}common/KeyRing_Class.f90 -o $@

HashTable_Class.o : ConstantsDictionary.o LinkedList_Class.o
	${FC} ${OPT} -c ${SRCDIR}common/HashTable_Class.f90 -o $@

CommonRoutines.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}common/CommonRoutines.f90 -o $@

Timing.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}common/Timing.f90 -o $@

# Test programs
   
# ---------------------------------------------------------------------- #
#                           ~/src/spectralops/                           #
# ---------------------------------------------------------------------- #

# Modules
Quadrature.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}spectralops/Quadrature.f90 -o $@

InterpolationSupportRoutines.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o 
	${FC} ${OPT} -c ${SRCDIR}spectralops/InterpolationSupportRoutines.f90 -o $@

Lagrange_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o InterpolationSupportRoutines.o
	${FC} ${OPT} -c ${SRCDIR}spectralops/Lagrange_Class.f90 -o $@

Lagrange_Cuda_Class.o : ModelPrecision.o ConstantsDictionary.o Lagrange_Class.o
	${FC} ${OPT} -c ${SRCDIR}spectralops/cuda/Lagrange_Cuda_Class.cuf -o $@

NodalStorage_Class.o : ModelPrecision.o ConstantsDictionary.o  CommonRoutines.o \
                       Quadrature.o Lagrange_Class.o
	${FC} ${OPT} -c ${SRCDIR}spectralops/NodalStorage_Class.f90 -o $@

NodalStorage_Cuda_Class.o : ModelPrecision.o ConstantsDictionary.o  CommonRoutines.o \
                       Quadrature.o Lagrange_Cuda_Class.o
	${FC} ${OPT} -c ${SRCDIR}spectralops/cuda/NodalStorage_Cuda_Class.cuf -o $@


# ---------------------------------------------------------------------- #
#                          Other directives                              #
# ---------------------------------------------------------------------- #
.PHONY : clean

clean :
	rm *.o *.mod



















