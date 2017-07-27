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
OPT=-fast -Mcuda -ta=tesla:cc60
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
#                           ~/src/geom/                                  #
# ---------------------------------------------------------------------- #

#
BoundaryCommunicator_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}geom/BoundaryCommunicator_Class.f90 -o $@

Node_Class.o : ModelPrecision.o  ConstantsDictionary.o  LinkedList_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/Node_Class.f90 -o $@

Edge_Class.o : ModelPrecision.o  ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}geom/Edge_Class.f90 -o $@
	
Face_Class.o : ModelPrecision.o  ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}geom/Face_Class.f90 -o $@

Curve_Class.o : ModelPrecision.o ConstantsDictionary.o Lagrange_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/Curve_Class.f90 -o $@
	
Surface_Class.o : ModelPrecision.o ConstantsDictionary.o Lagrange_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/Surface_Class.f90 -o $@

QuadElement_Class.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}geom/QuadElement_Class.f90 -o $@
	
HexElement_Class.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}geom/HexElement_Class.f90 -o $@

MappedGeometry_2D_Class.o : ModelPrecision.o ConstantsDictionary.o  CommonRoutines.o \
                            Lagrange_Class.o Curve_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/MappedGeometry_2D_Class.f90 -o $@
	
MappedGeometry_3D_Class.o : ModelPrecision.o ConstantsDictionary.o  CommonRoutines.o \
                            Lagrange_Class.o Surface_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/MappedGeometry_3D_Class.f90 -o $@

QuadMesh_Class.o : ModelPrecision.o ConstantsDictionary.o  LinkedList_Class.o HashTable_Class.o \
                   Quadrature.o Curve_Class.o MappedGeometry_2D_Class.o QuadElement_Class.o \
                   Edge_Class.o Node_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/QuadMesh_Class.f90 -o $@

HexMesh_Class.o : ModelPrecision.o ConstantsDictionary.o  LinkedList_Class.o KeyRing_Class.o \
                  Quadrature.o Surface_Class.o MappedGeometry_3D_Class.o HexElement_Class.o \
                  Face_Class.o Node_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/HexMesh_Class.f90 -o $@

# /////////// CUDA Extensions /////////// #

BoundaryCommunicator_Cuda_Class.o : BoundaryCommunicator_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/BoundaryCommunicator_Cuda_Class.cuf -o $@
	
Node_Cuda_Class.o : ModelPrecision.o Node_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/Node_Cuda_Class.cuf -o $@
	
Edge_Cuda_Class.o : ModelPrecision.o Edge_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/Edge_Cuda_Class.cuf -o $@
	
Face_Cuda_Class.o : ModelPrecision.o Face_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/Face_Cuda_Class.cuf -o $@
	
Element_Cuda_Class.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/Element_Cuda_Class.cuf -o $@

MappedGeometry_2D_Cuda_Class.o : ModelPrecision.o MappedGeometry_2D_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/MappedGeometry_2D_Cuda_Class.cuf -o $@
	
MappedGeometry_3D_Cuda_Class.o : ModelPrecision.o MappedGeometry_3D_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/MappedGeometry_3D_Cuda_Class.cuf -o $@

QuadMesh_Cuda_Class.o : ModelPrecision.o ConstantsDictionary.o QuadMesh_Class.o \
                        Lagrange_Class.o Element_Cuda_Class.o Node_Cuda_Class.o Edge_Cuda_Class.o \
                        Curve_Class.o MappedGeometry_2D_Cuda_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/QuadMesh_Cuda_Class.cuf -o $@

HexMesh_Cuda_Class.o : ModelPrecision.o ConstantsDictionary.o HexMesh_Class.o \
                        Lagrange_Class.o Element_Cuda_Class.o Node_Cuda_Class.o Face_Cuda_Class.o \
                        Surface_Class.o MappedGeometry_3D_Cuda_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/HexMesh_Cuda_Class.cuf -o $@


# ---------------------------------------------------------------------- #
#                          Other directives                              #
# ---------------------------------------------------------------------- #
.PHONY : clean

clean :
	rm *.o *.mod



















