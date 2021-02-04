
FC=gfortran
INC=-I/opt/self/include/FLAP -I/opt/self/include/PENF -I/opt/self/include/FACE -I/opt/self/include
FLIBS=-L/opt/self/lib/ -lFLAP -lFACE -lPENF -lfeqparse
PREFIX=/opt/self

ifeq ($(BUILD_TYPE),debug)
  FFLAGS=-O0 -g -pg
# -ftest-coverage -fprofile-arcs
else ifeq ($(BUILD_TYPE),benchmark)
  FFLAGS=-O3 -pg
else ifeq ($(BUILD_TYPE),release)
  FFLAGS=-O3
else
  FFLAGS=-O2 -g
endif

ifeq ($(GPU),yes)
    FFLAGS+= -DGPU
endif

install: libSELF.a self
	mkdir -p ${PREFIX}/bin
	mkdir -p ${PREFIX}/lib
	mkdir -p ${PREFIX}/include
	mv libSELF.a ${PREFIX}/lib/
	mv *.mod ${PREFIX}/include/
	sed -i 's/INSTALL_ROOT=.*/INSTALL_ROOT=\/opt\/self/g' test/ci.sh 
	sed -i 's/INSTALL_ROOT=.*/INSTALL_ROOT=\/opt\/self/g' test/ci.gpu.sh 
	cp test/ci.sh ${PREFIX}/bin/
	chmod 755 ${PREFIX}/bin/ci.sh
	cp test/ci.gpu.sh ${PREFIX}/bin/
	chmod 755 ${PREFIX}/bin/ci.gpu.sh
	cp src/*.h ${PREFIX}/include/
	mv self ${PREFIX}/bin/
	rm *.o

self: SELF.o
	${FC} ${INC} -I./src/ *.o  ${FLIBS} -o $@

SELF.o: libSELF.a
	${FC} -c ${FFLAGS} ${INC} -I./src/ ${FLIBS} src/SELF.F90 -o $@

libSELF.a: SELF_Constants.o SELF_Data.o SELF_Lagrange.o SELF_Geometry.o SELF_MappedData.o SELF_Memory.o SELF_Mesh.o SELF_Quadrature.o SELF_SupportRoutines.o SELF_Tests.o
	ar rcs $@ SELF_Constants.o SELF_Data.o SELF_Lagrange.o SELF_Geometry.o SELF_MappedData.o SELF_Memory.o SELF_Mesh.o SELF_Quadrature.o SELF_SupportRoutines.o SELF_Tests.o

SELF_Tests.o: SELF_Constants.o SELF_Data.o SELF_Lagrange.o SELF_MappedData.o SELF_Memory.o SELF_Mesh.o SELF_Quadrature.o SELF_SupportRoutines.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Tests.F90 -o $@

#SELF_MPILayer.o: SELF_MPILayer.F90

SELF_MappedData.o: SELF_Mesh.o SELF_Data.o SELF_Geometry.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_MappedData.F90 -o $@

#SELF_MappedData_HIP.o:
#	${CXX} -c ${CXXFLAGS} src/SELF_MappedData_HIP.cpp -o $@

SELF_Geometry.o: SELF_Data.o SELF_Lagrange.o SELF_Mesh.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Geometry.F90 -o $@

SELF_Mesh.o: SELF_Data.o SELF_Lagrange.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Mesh.F90 -o $@

SELF_Data.o: SELF_Lagrange.o SELF_Constants.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Data.F90 -o $@

SELF_Lagrange.o: SELF_Memory.o SELF_Quadrature.o SELF_Constants.o SELF_SupportRoutines.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Lagrange.F90 -o $@

#SELF_Lagrange_HIP.o:
#	${CXX} -c ${CXXFLAGS} src/SELF_Lagrange_HIP.cpp -o $@

SELF_Quadrature.o: SELF_Constants.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Quadrature.F90 -o $@

SELF_Memory.o:
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Memory.F90 -o $@

SELF_SupportRoutines.o: SELF_Constants.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_SupportRoutines.F90 -o $@

SELF_Constants.o:
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Constants.F90 -o $@
