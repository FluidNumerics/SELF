# makefile
#
#
#  Directives :
#
#    

SHELL=/bin/bash

METISLIB=-L/usr/local/lib/ -lmetis

# determine the main directory for this experiment
SELFDIR=$(shell pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')
SRCDIR=${SELFDIR}/src/


ifeq (${FC},gfortran)

  $(info 'Compiling with GNU Compilers')

  CPPFLAGS=-cpp 

  ifeq (${DEBUG},yes)
    OPT+=-O0 -g -ffree-line-length-none -ffpe-trap=invalid -fbacktrace
  else
    OPT+=-Ofast -ffree-line-length-none
  endif

  ifeq (${OPENMP},yes)
    $(info '+OpenMP')
    OPT+=-fopenmp
    CPPFLAGS+=-DHAVE_OPENMP
  endif

  ifeq (${MPI},yes)
    $(info '+MPI')
    CPPFLAGS+=-DHAVE_MPI
    LIB+=${MPI_LIB}
    INC+=${MPI_INC}
  endif
 
    

else ifeq (${FC},pgfortran)

  $(info 'Compiling with PGI Compilers')
  CPPFLAGS=-Mpreprocess 

  ifeq (${DEBUG},yes)
    OPT+=-O0 -g -Kieee -traceback -Minfo=all -Mbounds
    CPPFLAGS+=-DDEBUG
  else
    OPT+=-fast
  endif

  ifeq (${CUDA},yes)
  
    $(info '+CUDA-Fortran')
    OPT+=-Mcuda=${GPU_ARCH}
    ifeq (${PTX_INFO},yes)
      OPT+=,ptxinfo
    endif

    CPPFLAGS+=-DHAVE_CUDA

    ifeq (${CUDA_DIRECT},yes)
      $(info '+CUDA-DIRECT')
      CPPFLAGS+=-DCUDA_DIRECT
    endif

  endif
  
  ifeq (${MPI},yes)
    $(info '+MPI')
    CPPFLAGS+=-DHAVE_MPI
    FC=${MPIFC}
    LIB+=${MPI_LIB}
    INC+=${MPI_INC}
  endif

  
endif
  

ifeq (${DOUBLE_PRECISION},yes)
  CPPFLAGS+=-DDOUBLE_PRECISION
endif

ifeq (${TIMING},yes)
  CPPFLAGS+=-DTIMING
endif

ifeq (${INSITU_VIZ},yes)
  ifeq (${MPI},yes)
    $(info '>>> MPI with in-situ graphics is currently not supported.')
  else
    CPPFLAGS+=-DINSITU_VIZ
    INC+=${LIBSIM_INC}
    LIB+=${LIBSIM_LIB}
  endif
endif

ifeq (${DIAGNOSTICS},yes)
  CPPFLAGS+=-DDIAGNOSTICS
endif

LIB+=${LAPACK_LIB}
INC+=${LAPACK_INC}

# ---------------------------------------------------------------------- #
#                           ~/src/common/                                #
# ---------------------------------------------------------------------- #

# Modules
ModelPrecision.o : ${SRCDIR}common/ModelPrecision.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}common/ModelPrecision.f90 -o $@

ConstantsDictionary.o : ModelPrecision.o ${SRCDIR}common/ConstantsDictionary.f90
	${FC} ${OPT} -c ${SRCDIR}common/ConstantsDictionary.f90 -o $@

LinkedList_Class.o : ${SRCDIR}common/LinkedList_Class.f90
	${FC} ${OPT} -c ${SRCDIR}common/LinkedList_Class.f90 -o $@

KeyRing_Class.o : CommonRoutines.o ${SRCDIR}common/KeyRing_Class.f90
	${FC} ${OPT} -c ${SRCDIR}common/KeyRing_Class.f90 -o $@

HashTable_Class.o : ConstantsDictionary.o LinkedList_Class.o ${SRCDIR}common/HashTable_Class.f90
	${FC} ${OPT} -c ${SRCDIR}common/HashTable_Class.f90 -o $@

CommonRoutines.o : ModelPrecision.o ConstantsDictionary.o ${SRCDIR}common/CommonRoutines.f90
	${FC} ${OPT} -c ${SRCDIR}common/CommonRoutines.f90 -o $@

Timing.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o ${SRCDIR}common/Timing.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}common/Timing.f90 -o $@

   
# ---------------------------------------------------------------------- #
#                           ~/src/spectral/                           #
# ---------------------------------------------------------------------- #

# Modules
Quadrature.o : ModelPrecision.o ConstantsDictionary.o ${SRCDIR}spectral/Quadrature.f90
	${FC} ${OPT} -c ${SRCDIR}spectral/Quadrature.f90 -o $@

Lagrange_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o ${SRCDIR}spectral/Lagrange_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}spectral/Lagrange_Class.f90 -o $@

NodalDG_Class.o : ModelPrecision.o ConstantsDictionary.o  CommonRoutines.o \
                  Quadrature.o Lagrange_Class.o ${SRCDIR}spectral/NodalDG_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}spectral/NodalDG_Class.f90 -o $@

NodalDGSolution_3D_Class.o : ModelPrecision.o ${SRCDIR}spectral/NodalDGSolution_3D_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}spectral/NodalDGSolution_3D_Class.f90 -o $@
	

# Test program

spectral_tests : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Quadrature.o Lagrange_Class.o NodalDG_Class.o NodalDGSolution_3D_Class.o Spectral_Tests.o
	${FC} ${OPT} ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Quadrature.o Lagrange_Class.o NodalDG_Class.o NodalDGSolution_3D_Class.o Spectral_Tests.o -o $@


Spectral_Tests.o : ModelPrecision.o CommonRoutines.o Lagrange_Class.o NodalDG_Class.o NodalDGSolution_3D_Class.o \
                   ${SRCDIR}spectral/testing/Spectral_Tests.F90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}spectral/testing/Spectral_Tests.F90 -o $@


# ---------------------------------------------------------------------- #
#                           ~/src/filters/                               #
# ---------------------------------------------------------------------- #

ModalCutoffFilter_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o ${SRCDIR}filters/ModalCutoffFilter_Class.f90
	${FC} ${OPT} -c ${SRCDIR}filters/ModalCutoffFilter_Class.f90 -o $@

RollOffFilter_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o ${SRCDIR}filters/RollOffFilter_Class.f90
	${FC} ${OPT} -c ${SRCDIR}filters/RollOffFilter_Class.f90 -o $@


# CUDA Extensions

RollOffFilter_Cuda_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o RollOffFilter_Class.o ${SRCDIR}filters/cuda/RollOffFilter_Cuda_Class.cuf
	${FC} ${OPT} -c ${SRCDIR}filters/cuda/RollOffFilter_Cuda_Class.cuf -o $@

   
# ---------------------------------------------------------------------- #
#                        ~/src/integration_testing/                      #
# ---------------------------------------------------------------------- #

cmdi : ModelPrecision.o CommonRoutines.o ModelDataInstances_Class.o CompareModelDataInstances.o
	${FC} ${OPT} ModelPrecision.o CommonRoutines.o ModelDataInstances_Class.o CompareModelDataInstances.o -o $@

ModelDataInstances_Class.o : ${SRCDIR}integration_testing/ModelDataInstances_Class.f90 ModelPrecision.o CommonRoutines.o
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}integration_testing/ModelDataInstances_Class.f90 -o $@

CompareModelDataInstances.o : ${SRCDIR}integration_testing/CompareModelDataInstances.f90 ModelDataInstances_Class.o ModelPrecision.o
	${FC} ${OPT} -c ${SRCDIR}integration_testing/CompareModelDataInstances.f90 -o $@

# ---------------------------------------------------------------------- #
#                            ~/src/geom/                                 #
# ---------------------------------------------------------------------- #

GOBJS = ModelPrecision.o \
        ConstantsDictionary.o \
        CommonRoutines.o \
        LinkedList_Class.o \
        HashTable_Class.o \
        Timing.o \
        Quadrature.o \
        Lagrange_Class.o \
        NodalStorage_Class.o \
        DGSEMSolutionStorage_2D_Class.o \
        Node_Class.o \
        Edge_Class.o \
        Curve_Class.o \
        MappedGeometry_2D_Class.o \
        QuadElement_Class.o \
        QuadMesh_Class.o \
        Params_Class.o \
        BoundaryCommunicator_Class.o

         
DecomposeQuadMesh : ${GOBJS} DecomposeQuadMesh.o
	${FC} ${OPT} ${GOBJS} DecomposeQuadMesh.o ${METISLIB} -o $@

GeneratePCMesh_2D : ${GOBJS} GenerateMeshFiles_2D.o 
	${FC} ${OPT} ${GOBJS} GenerateMeshFiles_2D.o  -o $@

GOBJ3 = ModelPrecision.o \
        ConstantsDictionary.o \
        CommonRoutines.o \
        LinkedList_Class.o \
        HashTable_Class.o \
        KeyRing_Class.o \
        Timing.o \
        Quadrature.o \
        Lagrange_Class.o \
        NodalStorage_Class.o \
        DGSEMSolutionStorage_2D_Class.o \
        Node_Class.o \
        Face_Class.o \
        Curve_Class.o \
        Surface_Class.o \
        MappedGeometry_3D_Class.o \
        HexElement_Class.o \
        HexMesh_Class.o \
        TopographicShapes.o \
        FluidParams_Class.o \
        BoundaryCommunicator_Class.o

DecomposeHexMesh : ${GOBJ3} DecomposeHexMesh.o
	${FC} ${OPT} ${GOBJ3} DecomposeHexMesh.o ${METISLIB} -o $@

meshgen : ${GOBJ3} DecomposeStructuredHexMesh.o
	${FC} ${OPT} ${GOBJ3} DecomposeStructuredHexMesh.o -o $@

GenerateHexMesh : ${GOBJ3} GenerateMeshFiles_3D.o 
	${FC} ${OPT} ${GOBJ3} GenerateMeshFiles_3D.o  -o $@


EQMOBJ = ModelPrecision.o \
         ConstantsDictionary.o \
         CommonRoutines.o \
         LinkedList_Class.o \
         HashTable_Class.o \
         KeyRing_Class.o \
         Timing.o \
         Quadrature.o \
         Lagrange_Class.o \
         NodalStorage_Class.o \
         Node_Class.o \
         Edge_Class.o \
         Face_Class.o \
         Curve_Class.o \
         Surface_Class.o \
         MappedGeometry_2D_Class.o \
         MappedGeometry_3D_Class.o \
         QuadElement_Class.o \
         HexElement_Class.o \
         QuadMesh_Class.o \
         HexMesh_Class.o \
         FluidParams_Class.o \
         BoundaryCommunicator_Class.o

ExtrudeQuadMesh : ${EQMOBJ} ExtrudeQuadMesh.o 
	${FC} ${OPT} ${EQMOBJ} ExtrudeQuadMesh.o  -o $@

## Programs
GenerateMeshFiles_2D.o : ModelPrecision.o NodalStorage_Class.o QuadMesh_Class.o \
                         BoundaryCommunicator_Class.o FluidParams_Class.o ${SRCDIR}geom/GenerateMeshFiles_2D.f90
	${FC} ${OPT} -c ${SRCDIR}geom/GenerateMeshFiles_2D.f90 -o $@

GenerateMeshFiles_3D.o : ModelPrecision.o NodalStorage_Class.o HexMesh_Class.o \
                         BoundaryCommunicator_Class.o TopographicShapes.o FluidParams_Class.o ${SRCDIR}geom/GenerateMeshFiles_3D.f90
	${FC} ${OPT} -c ${SRCDIR}geom/GenerateMeshFiles_3D.f90 -o $@

ExtrudeQuadMesh.o : ModelPrecision.o NodalStorage_Class.o QuadMesh_Class.o HexMesh_Class.o \
                    BoundaryCommunicator_Class.o FluidParams_Class.o ${SRCDIR}geom/ExtrudeQuadMesh.f90
	${FC} ${OPT} -c ${SRCDIR}geom/ExtrudeQuadMesh.f90 -o $@

DecomposeQuadMesh.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Timing.o \
                      NodalStorage_Class.o QuadMesh_Class.o FluidParams_Class.o ${SRCDIR}geom/DecomposeQuadMesh.f90
	${FC} ${OPT} -c ${SRCDIR}geom/DecomposeQuadMesh.f90 -o $@


DecomposeHexMesh.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Timing.o \
                      NodalStorage_Class.o HexMesh_Class.o TopographicShapes.o FluidParams_Class.o ${SRCDIR}geom/DecomposeHexMesh.f90
	${FC} ${OPT} -c ${SRCDIR}geom/DecomposeHexMesh.f90 -o $@

DecomposeStructuredHexMesh.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Timing.o \
                      NodalStorage_Class.o HexMesh_Class.o TopographicShapes.o FluidParams_Class.o ${SRCDIR}geom/DecomposeStructuredHexMesh.f90
	${FC} ${OPT} -c ${SRCDIR}geom/DecomposeStructuredHexMesh.f90 -o $@
	
## Modules 
#
BoundaryCommunicator_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o ${SRCDIR}geom/BoundaryCommunicator_Class.f90
	${FC} ${OPT} -c ${SRCDIR}geom/BoundaryCommunicator_Class.f90 -o $@

Node_Class.o : ModelPrecision.o  ConstantsDictionary.o  LinkedList_Class.o ${SRCDIR}geom/Node_Class.f90
	${FC} ${OPT} -c ${SRCDIR}geom/Node_Class.f90 -o $@

Edge_Class.o : ModelPrecision.o  ConstantsDictionary.o ${SRCDIR}geom/Edge_Class.f90
	${FC} ${OPT} -c ${SRCDIR}geom/Edge_Class.f90 -o $@
	
Face_Class.o : ModelPrecision.o  ConstantsDictionary.o ${SRCDIR}geom/Face_Class.f90
	${FC} ${OPT} -c ${SRCDIR}geom/Face_Class.f90 -o $@

Curve_Class.o : ModelPrecision.o ConstantsDictionary.o Lagrange_Class.o ${SRCDIR}geom/Curve_Class.f90
	${FC} ${OPT} -c ${SRCDIR}geom/Curve_Class.f90 -o $@
	
Surface_Class.o : ModelPrecision.o ConstantsDictionary.o Lagrange_Class.o ${SRCDIR}geom/Surface_Class.f90
	${FC} ${OPT} -c ${SRCDIR}geom/Surface_Class.f90 -o $@

QuadElement_Class.o : ModelPrecision.o ConstantsDictionary.o ${SRCDIR}geom/QuadElement_Class.f90
	${FC} ${OPT} -c ${SRCDIR}geom/QuadElement_Class.f90 -o $@
	
HexElement_Class.o : ModelPrecision.o ConstantsDictionary.o ${SRCDIR}geom/HexElement_Class.f90
	${FC} ${OPT} -c ${SRCDIR}geom/HexElement_Class.f90 -o $@

MappedGeometry_2D_Class.o : ModelPrecision.o ConstantsDictionary.o  CommonRoutines.o \
                            Lagrange_Class.o Curve_Class.o ${SRCDIR}geom/MappedGeometry_2D_Class.f90
	${FC} ${OPT} -c ${SRCDIR}geom/MappedGeometry_2D_Class.f90 -o $@
	
MappedGeometry_3D_Class.o : ModelPrecision.o ConstantsDictionary.o  CommonRoutines.o \
                            Lagrange_Class.o Surface_Class.o ${SRCDIR}geom/MappedGeometry_3D_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}geom/MappedGeometry_3D_Class.f90 -o $@

QuadMesh_Class.o : ModelPrecision.o ConstantsDictionary.o  LinkedList_Class.o HashTable_Class.o \
                   Quadrature.o Curve_Class.o MappedGeometry_2D_Class.o QuadElement_Class.o \
                   Edge_Class.o Node_Class.o ${SRCDIR}geom/QuadMesh_Class.f90
	${FC} ${OPT} -c ${SRCDIR}geom/QuadMesh_Class.f90 -o $@

HexMesh_Class.o : ModelPrecision.o ConstantsDictionary.o  LinkedList_Class.o KeyRing_Class.o \
                  Quadrature.o Surface_Class.o MappedGeometry_3D_Class.o HexElement_Class.o \
                  Face_Class.o Node_Class.o ${SRCDIR}geom/HexMesh_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}geom/HexMesh_Class.f90 -o $@

TopographicShapes.o : ModelPrecision.o ${SRCDIR}geom/TopographicShapes.f90 ${SRCDIR}geom/TopographicShapes.f90
	${FC} ${OPT} -c ${SRCDIR}geom/TopographicShapes.f90 -o $@ 

# /////////// CUDA Extensions /////////// #

BoundaryCommunicator_Cuda_Class.o : BoundaryCommunicator_Class.o ${SRCDIR}geom/cuda/BoundaryCommunicator_Cuda_Class.cuf
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/BoundaryCommunicator_Cuda_Class.cuf -o $@
	
Node_Cuda_Class.o : ModelPrecision.o Node_Class.o ${SRCDIR}geom/cuda/Node_Cuda_Class.cuf
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/Node_Cuda_Class.cuf -o $@
	
Edge_Cuda_Class.o : ModelPrecision.o Edge_Class.o ${SRCDIR}geom/cuda/Edge_Cuda_Class.cuf
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/Edge_Cuda_Class.cuf -o $@
	
Face_Cuda_Class.o : ModelPrecision.o Face_Class.o ${SRCDIR}geom/cuda/Face_Cuda_Class.cuf
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/Face_Cuda_Class.cuf -o $@
	
Element_Cuda_Class.o : ModelPrecision.o ConstantsDictionary.o ${SRCDIR}geom/cuda/Element_Cuda_Class.cuf
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/Element_Cuda_Class.cuf -o $@

MappedGeometry_2D_Cuda_Class.o : ModelPrecision.o MappedGeometry_2D_Class.o ${SRCDIR}geom/cuda/MappedGeometry_2D_Cuda_Class.cuf
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/MappedGeometry_2D_Cuda_Class.cuf -o $@
	
MappedGeometry_3D_Cuda_Class.o : ModelPrecision.o MappedGeometry_3D_Class.o ${SRCDIR}geom/cuda/MappedGeometry_3D_Cuda_Class.cuf
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/MappedGeometry_3D_Cuda_Class.cuf -o $@

QuadMesh_Cuda_Class.o : ModelPrecision.o ConstantsDictionary.o QuadMesh_Class.o \
                        Lagrange_Class.o Element_Cuda_Class.o Node_Cuda_Class.o Edge_Cuda_Class.o \
                        Curve_Class.o MappedGeometry_2D_Cuda_Class.o ${SRCDIR}geom/cuda/QuadMesh_Cuda_Class.cuf
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/QuadMesh_Cuda_Class.cuf -o $@

HexMesh_Cuda_Class.o : ModelPrecision.o ConstantsDictionary.o HexMesh_Class.o \
                        Lagrange_Class.o Lagrange_Cuda_Class.o Element_Cuda_Class.o Node_Cuda_Class.o Face_Cuda_Class.o \
                        Surface_Class.o MappedGeometry_3D_Cuda_Class.o ${SRCDIR}geom/cuda/HexMesh_Cuda_Class.cuf
	${FC} ${OPT} -c ${SRCDIR}geom/cuda/HexMesh_Cuda_Class.cuf -o $@

# Cuda Test Programs

# --------- QuadMesh ----------- #
TestQuadMesh_Cuda : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o LinkedList_Class.o HashTable_Class.o \
                   Quadrature.o Lagrange_Class.o NodalStorage_Class.o \
                   Node_Class.o Edge_Class.o QuadElement_Class.o Curve_Class.o MappedGeometry_2D_Class.o  QuadMesh_Class.o\
                   Node_Cuda_Class.o Edge_Cuda_Class.o Element_Cuda_Class.o MappedGeometry_2D_Cuda_Class.o \
                   QuadMesh_Cuda_Class.o FluidParams_Class.o TestQuadMesh_Cuda.o
	${FC} ${OPT} ModelPrecision.o ConstantsDictionary.o CommonRoutines.o LinkedList_Class.o HashTable_Class.o \
                   Quadrature.o Lagrange_Class.o NodalStorage_Class.o \
                   Node_Class.o Edge_Class.o QuadElement_Class.o Curve_Class.o MappedGeometry_2D_Class.o  QuadMesh_Class.o\
                   Node_Cuda_Class.o Edge_Cuda_Class.o Element_Cuda_Class.o MappedGeometry_2D_Cuda_Class.o \
                   QuadMesh_Cuda_Class.o FluidParams_Class.o TestQuadMesh_Cuda.o -o $@

TestQuadMesh_Cuda.o : ConstantsDictionary.o NodalStorage_Class.o FluidParams_Class.o QuadMesh_Cuda_Class.o ${SRCDIR}geom/tests/QuadMesh_Cuda/IOandDataTransfers/TestQuadMesh_Cuda.f90
	${FC} ${OPT} -c ${SELFDIR}/tests/QuadMesh_Cuda/IOandDataTransfers/TestQuadMesh_Cuda.f90 -o $@

# --------- HexMesh ----------- #
TestHexMesh_Cuda : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o LinkedList_Class.o KeyRing_Class.o \
                   Quadrature.o Lagrange_Class.o NodalStorage_Class.o \
                   Node_Class.o Face_Class.o HexElement_Class.o Surface_Class.o MappedGeometry_3D_Class.o  HexMesh_Class.o\
                   Node_Cuda_Class.o Face_Cuda_Class.o Element_Cuda_Class.o MappedGeometry_3D_Cuda_Class.o \
                   HexMesh_Cuda_Class.o FluidParams_Class.o TestHexMesh_Cuda.o
	${FC} ${OPT} ModelPrecision.o ConstantsDictionary.o CommonRoutines.o LinkedList_Class.o KeyRing_Class.o \
                   Quadrature.o Lagrange_Class.o NodalStorage_Class.o \
                   Node_Class.o Face_Class.o HexElement_Class.o Surface_Class.o MappedGeometry_3D_Class.o  HexMesh_Class.o\
                   Node_Cuda_Class.o Face_Cuda_Class.o Element_Cuda_Class.o MappedGeometry_3D_Cuda_Class.o \
                   HexMesh_Cuda_Class.o FluidParams_Class.o TestHexMesh_Cuda.o -o $@

TestHexMesh_Cuda.o : ConstantsDictionary.o NodalStorage_Class.o FluidParams_Class.o HexMesh_Cuda_Class.o ${SRCDIR}geom/tests/HexMesh_Cuda/IOandDataTransfers/TestHexMesh_Cuda.f90
	${FC} ${OPT} -c ${SELFDIR}/tests/HexMesh_Cuda/IOandDataTransfers/TestHexMesh_Cuda.f90 -o $@


# ---------------------------------------------------------------------- #
#                             ~/src/viz                                  #
# ---------------------------------------------------------------------- #

Ray_Class.o : ModelPrecision.o CommonRoutines.o Lagrange_Class.o BoundaryCommunicator_Class.o \
              MappedGeometry_3D_Class.o HexMesh_Class.o ${SRCDIR}viz/Ray_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}viz/Ray_Class.f90 -o $@

# ---------------------------------------------------------------------- #
#                             ~/src/fluid                                #
# ---------------------------------------------------------------------- #
	

ifeq (${CUDA},yes)

CNOBJ = ModelPrecision.o \
        ConstantsDictionary.o \
        CommonRoutines.o

ifeq (${TESTING},yes)
CNOBJ += ModelDataInstances_Class.o
endif

CNOBJ += Timing.o \
         LinkedList_Class.o \
         HashTable_Class.o \
         KeyRing_Class.o \
         Quadrature.o \
         Lagrange_Class.o \
         Lagrange_Cuda_Class.o \
         NodalStorage_Class.o \
         NodalStorage_Cuda_Class.o \
         RollOffFilter_Class.o \
         RollOffFilter_Cuda_Class.o \
         DGSEMSolutionStorage_3D_Class.o \
         DGSEMSolutionStorage_3D_Cuda_Class.o \
         Node_Class.o \
         Node_Cuda_Class.o \
         Face_Class.o \
         Face_Cuda_Class.o \
         Curve_Class.o \
         Surface_Class.o \
         MappedGeometry_3D_Class.o \
         MappedGeometry_3D_Cuda_Class.o \
         HexElement_Class.o \
         Element_Cuda_Class.o \
         HexMesh_Class.o \
         HexMesh_Cuda_Class.o \
         FluidParams_Class.o \
         BoundaryCommunicator_Class.o \
         BoundaryCommunicator_Cuda_Class.o \
         Fluid_Class.o
 
EUDEP = ModelPrecision.o ConstantsDictionary.o CommonRoutines.o \
                Lagrange_Class.o Lagrange_Cuda_Class.o NodalStorage_Class.o NodalStorage_Cuda_Class.o\
                RollOffFilter_Cuda_Class.o DGSEMSolutionStorage_3D_Cuda_Class.o \
                HexMesh_Class.o BoundaryCommunicator_Class.o \
                Face_Cuda_Class.o Element_Cuda_Class.o HexMesh_Cuda_Class.o BoundaryCommunicator_Cuda_Class.o \
                FluidParams_Class.o
                
else

ifeq ($(TESTING),yes)
CNOBJ = ModelPrecision.o \
        ConstantsDictionary.o \
        CommonRoutines.o \
        Timing.o \
        LinkedList_Class.o \
        HashTable_Class.o \
         KeyRing_Class.o \
         ModelDataInstances_Class.o \
         Quadrature.o \
         Lagrange_Class.o \
         NodalStorage_Class.o \
         RollOffFilter_Class.o \
         DGSEMSolutionStorage_3D_Class.o \
         Node_Class.o \
         Face_Class.o \
         Curve_Class.o \
         Surface_Class.o \
         MappedGeometry_3D_Class.o \
         HexElement_Class.o \
         HexMesh_Class.o \
         FluidParams_Class.o \
         BoundaryCommunicator_Class.o \
         Fluid_Class.o
else
CNOBJ = ModelPrecision.o \
        ConstantsDictionary.o \
        CommonRoutines.o \
        Timing.o \
        LinkedList_Class.o \
        HashTable_Class.o \
         KeyRing_Class.o \
         Quadrature.o \
         Lagrange_Class.o \
         NodalStorage_Class.o \
         RollOffFilter_Class.o \
         DGSEMSolutionStorage_3D_Class.o \
         Node_Class.o \
         Face_Class.o \
         Curve_Class.o \
         Surface_Class.o \
         MappedGeometry_3D_Class.o \
         HexElement_Class.o \
         HexMesh_Class.o \
         FluidParams_Class.o \
         BoundaryCommunicator_Class.o \
         Fluid_Class.o

endif
         
EUDEP = ModelPrecision.o ConstantsDictionary.o CommonRoutines.o \
                Lagrange_Class.o NodalStorage_Class.o \
                RollOffFilter_Class.o DGSEMSolutionStorage_3D_Class.o \
                HexMesh_Class.o BoundaryCommunicator_Class.o \
                Face_Class.o HexElement_Class.o HexMesh_Class.o BoundaryCommunicator_Class.o \
                FluidParams_Class.o
                
endif

FluidParams_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o ${SRCDIR}fluid/FluidParams_Class.f90
	${FC} ${OPT} -c ${SRCDIR}fluid/FluidParams_Class.f90 -o $@
	
Fluid_Class.o : ${EUDEP} ${SRCDIR}fluid/Fluid_Class.f90 ModelDataInstances_Class.o
	${FC} ${OPT} -c ${CPPFLAGS} ${SRCDIR}fluid/Fluid_Class.f90 ${LIB} -o $@

ifdef MODDIR
Fluid_InitialConditions.o : ModelPrecision.o FluidParams_Class.o Fluid_Class.o TopographicShapes.o ${MODDIR}/Fluid_InitialConditions.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${MODDIR}/Fluid_InitialConditions.f90 ${LIB} ${INC} -o $@
else
Fluid_InitialConditions.o : ModelPrecision.o FluidParams_Class.o Fluid_Class.o TopographicShapes.o ${SRCDIR}/fluid/Fluid_InitialConditions.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}fluid/Fluid_InitialConditions.f90 ${LIB} ${INC} -o $@
endif

Fluid_Driver.o : ModelPrecision.o FluidParams_Class.o Fluid_Class.o ${SRCDIR}fluid/Fluid_Driver.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}fluid/Fluid_Driver.f90 ${LIB} ${INC} -o $@

TestMPIStateExchange_Driver.o : ModelPrecision.o FluidParams_Class.o Fluid_Class.o ${SRCDIR}Fluid/TestMPIStateExchange_Driver.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}fluid/TestMPIStateExchange_Driver.f90 ${LIB} ${INC} -o $@
	
Fluid_PickupToTecplot.o : ModelPrecision.o FluidParams_Class.o Fluid_Class.o ${SRCDIR}Fluid/Fluid_PickupToTecplot.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}fluid/Fluid_PickupToTecplot.f90 ${LIB} ${INC} -o $@
	
SampleAlongRays.o : ModelPrecision.o Fluid_Class.o Ray_Class.o ${SRCDIR}fluid/SampleAlongRays.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}fluid/SampleAlongRays.f90 -o $@ 	

initialize: ${CNOBJ} Fluid_InitialConditions.o
	${FC} ${OPT} ${CNOBJ} Fluid_InitialConditions.o ${LIB} ${INC} -o $@

integrate: ${CNOBJ} Fluid_Driver.o
	${FC} ${OPT} ${CNOBJ} Fluid_Driver.o ${LIB} ${INC} -o $@
	
pickuptoptecplot : ${CNOBJ} Fluid_PickupToTecplot.o
	${FC} ${OPT} ${CNOBJ} Fluid_PickupToTecplot.o ${LIB} ${INC} -o $@
	
test_stateExchange: ${CNOBJ} TestMPIStateExchange_Driver.o
	${FC} ${OPT} ${CNOBJ} TestMPIStateExchange_Driver.o ${LIB} ${INC} -o $@
	
	
SOBJ = ModelPrecision.o \
       ConstantsDictionary.o \
       CommonRoutines.o \
       Lagrange_Class.o \
       Quadrature.o \
       NodalStorage_Class.o \
       RollOffFilter_Class.o \
       DGSEMSolutionStorage_3D_Class.o \
       LinkedList_Class.o \
       KeyRing_Class.o \
       Surface_Class.o \
       MappedGeometry_3D_Class.o \
       HexElement_Class.o \
       Face_Class.o \
       Node_Class.o \
       HexMesh_Class.o \
       BoundaryCommunicator_Class.o \
       FluidParams_Class.o \
       Fluid_Class.o \
       Ray_Class.o \
       SampleAlongRays.o

sars: ${SOBJ}
	${FC} ${OPT} ${CPPFLAGS} ${SOBJ} -o $@ 	

.PHONY : clean

clean :
	rm *.o *.mod



