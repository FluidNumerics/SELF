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

Quadrature.o : ModelPrecision.o ConstantsDictionary.o ${SRCDIR}spectral/Quadrature.f90
	${FC} ${OPT} -c ${SRCDIR}spectral/Quadrature.f90 -o $@

Lagrange_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o ${SRCDIR}spectral/Lagrange_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}spectral/Lagrange_Class.f90 -o $@

NodalDG_Class.o : ModelPrecision.o ConstantsDictionary.o  CommonRoutines.o \
                  Quadrature.o Lagrange_Class.o ${SRCDIR}spectral/NodalDG_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}spectral/NodalDG_Class.f90 -o $@

NodalDGSolution_3D_Class.o : ModelPrecision.o NodalDG_Class.o ${SRCDIR}spectral/NodalDGSolution_3D_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}spectral/NodalDGSolution_3D_Class.f90 -o $@
	

SpectralFilter_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Quadrature.o ${SRCDIR}spectral/SpectralFilter_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}spectral/SpectralFilter_Class.f90 -o $@



# Test program

spectral_tests : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Quadrature.o Lagrange_Class.o NodalDG_Class.o NodalDGSolution_3D_Class.o Spectral_Tests.o
	${FC} ${OPT} ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Quadrature.o Lagrange_Class.o NodalDG_Class.o NodalDGSolution_3D_Class.o Spectral_Tests.o -o $@


Spectral_Tests.o : ModelPrecision.o CommonRoutines.o Lagrange_Class.o NodalDG_Class.o NodalDGSolution_3D_Class.o \
                   ${SRCDIR}spectral/testing/Spectral_Tests.F90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}spectral/testing/Spectral_Tests.F90 -o $@


spectralfilter_tests : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Quadrature.o Lagrange_Class.o NodalDG_Class.o NodalDGSolution_3D_Class.o SpectralFilter_Class.o SpectralFilter_Tests.o
	${FC} ${OPT} ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Quadrature.o Lagrange_Class.o NodalDG_Class.o NodalDGSolution_3D_Class.o SpectralFilter_Class.o SpectralFilter_Tests.o -o $@


SpectralFilter_Tests.o : ModelPrecision.o CommonRoutines.o Lagrange_Class.o NodalDG_Class.o NodalDGSolution_3D_Class.o SpectralFilter_Class.o \
                   ${SRCDIR}spectral/testing/SpectralFilter_Tests.F90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}spectral/testing/SpectralFilter_Tests.F90 -o $@

# ---------------------------------------------------------------------- #
#                            ~/src/geom/                                 #
# ---------------------------------------------------------------------- #

Edges_Class.o : ModelPrecision.o ${SRCDIR}geom/Edges_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}geom/Edges_Class.f90 -o $@

Faces_Class.o : ModelPrecision.o  ConstantsDictionary.o ${SRCDIR}geom/Faces_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}geom/Faces_Class.f90 -o $@

Nodes_Class.o : ModelPrecision.o ${SRCDIR}geom/Nodes_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}geom/Nodes_Class.f90 -o $@
	
Surfaces_Class.o : ModelPrecision.o ConstantsDictionary.o Lagrange_Class.o ${SRCDIR}geom/Surfaces_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}geom/Surfaces_Class.f90 -o $@

HexElements_Class.o : ModelPrecision.o ConstantsDictionary.o  CommonRoutines.o \
                      Lagrange_Class.o Surfaces_Class.o ${SRCDIR}geom/HexElements_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}geom/HexElements_Class.f90 -o $@

BoundaryCommunicator_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o ${SRCDIR}geom/BoundaryCommunicator_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}geom/BoundaryCommunicator_Class.f90 -o $@

HexMesh_Class.o : ModelPrecision.o ConstantsDictionary.o  LinkedList_Class.o KeyRing_Class.o \
                  Quadrature.o Surfaces_Class.o HexElements_Class.o \
                  Faces_Class.o Nodes_Class.o ${SRCDIR}geom/HexMesh_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}geom/HexMesh_Class.f90 -o $@

TopographicShapes.o : ModelPrecision.o ${SRCDIR}geom/TopographicShapes.f90 ${SRCDIR}geom/TopographicShapes.f90
	${FC} ${OPT} -c ${SRCDIR}geom/TopographicShapes.f90 -o $@ 

MeshGenerator_3D.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Timing.o NodalDG_Class.o HexMesh_Class.o \
                     TopographicShapes.o FluidParams_Class.o BoundaryCommunicator_Class.o ${SRCDIR}geom/MeshGenerator_3D.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}geom/MeshGenerator_3D.f90

MeshGen_3D : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Timing.o Quadrature.o Lagrange_Class.o NodalDG_Class.o \
              LinkedList_Class.o KeyRing_Class.o Surfaces_Class.o HexElements_Class.o Faces_Class.o Nodes_Class.o HexMesh_Class.o \
              TopographicShapes.o FluidParams_Class.o BoundaryCommunicator_Class.o MeshGenerator_3D.o
	${FC} ${OPT}  ModelPrecision.o ConstantsDictionary.o CommonRoutines.o Timing.o Quadrature.o Lagrange_Class.o NodalDG_Class.o \
                      LinkedList_Class.o KeyRing_Class.o Surfaces_Class.o HexElements_Class.o Faces_Class.o Nodes_Class.o HexMesh_Class.o \
                      TopographicShapes.o FluidParams_Class.o BoundaryCommunicator_Class.o MeshGenerator_3D.o -o $@

	
# ---------------------------------------------------------------------- #
#                             ~/src/fluid                                #
# ---------------------------------------------------------------------- #

BodyForces_Class.o : ModelPrecision.o ${SRCDIR}fluid/BodyForces_Class.f90	
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}fluid/BodyForces_Class.f90

BoundaryConditions_Class.o : ModelPrecision.o ${SRCDIR}fluid/BoundaryConditions_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}fluid/BoundaryConditions_Class.f90

MPILayer_Class.o : ModelPrecision.o NodalDGSolution_3D_Class.o Faces_Class.o BoundaryConditions_Class.o BoundaryCommunicator_Class.o \
                   ${SRCDIR}fluid/MPILayer_Class.f90
	${FC} ${OPT} ${CPPFLAGS} -c ${SRCDIR}fluid/MPILayer_Class.f90

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



