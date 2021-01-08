

FC=gfortran
CXX=/opt/rocm/bin/hipcc
FFLAGS=-O0 -g -pg -DGPU
FLIBS=-L/apps/self/lib/ -lFLAP -lFACE -lPENF -lfeqparse
INC=-I/opt/hipfort/include/nvptx -I/apps/self/include/FLAP -I/apps/self/include/PENF -I/apps/self/include/FACE -I/apps/self/include
CXXFLAGS=
CXXLIBS=-O0 -g -pg -L/usr/local/cuda/lib64 -lcudart -L/opt/hipfort/lib -lhipfort-nvptx -lgfortran
PREFIX=/apps/self

install: libSELF.a self
	mkdir -p ${PREFIX}/bin
	mkdir -p ${PREFIX}/lib
	mkdir -p ${PREFIX}/include
	mv libSELF.a ${PREFIX}/lib/
	mv *.mod ${PREFIX}/include/
	cp ci.sh ${PREFIX}/bin
	cp src/*.h ${PREFIX}/include/
	mv self ${PREFIX}/bin/
	rm *.o

self: SELF.o
	${CXX} ${INC} -I./src/ *.o  ${FLIBS} ${CXXLIBS} -o $@

SELF.o: libSELF.a
	${FC} -c ${FFLAGS} ${INC} -I./src/ ${FLIBS} src/SELF.F90 -o $@

libSELF.a: SELF_Constants.o SELF_Data.o SELF_Lagrange_HIP.o SELF_Lagrange.o SELF_Geometry.o SELF_MappedData_HIP.o SELF_MappedData.o SELF_Memory.o SELF_Mesh.o SELF_Quadrature.o SELF_SupportRoutines.o SELF_Tests.o
	ar rcs $@ SELF_Constants.o SELF_Data.o SELF_Lagrange_HIP.o SELF_Lagrange.o SELF_Geometry.o SELF_MappedData_HIP.o SELF_MappedData.o SELF_Memory.o SELF_Mesh.o SELF_Quadrature.o SELF_SupportRoutines.o SELF_Tests.o

SELF_Tests.o: SELF_Constants.o SELF_Data.o SELF_Lagrange_HIP.o SELF_Lagrange.o SELF_MappedData_HIP.o SELF_MappedData.o SELF_Memory.o SELF_Mesh.o SELF_Quadrature.o SELF_SupportRoutines.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Tests.F90 -o $@

#SELF_MPILayer.o: SELF_MPILayer.F90

SELF_MappedData.o: SELF_Mesh.o SELF_Data.o SELF_Geometry.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_MappedData.F90 -o $@

SELF_MappedData_HIP.o:
	${CXX} -c ${CXXFLAGS} src/SELF_MappedData_HIP.cpp -o $@

SELF_Geometry.o: SELF_Data.o SELF_Lagrange.o SELF_Mesh.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Geometry.F90 -o $@

SELF_Mesh.o: SELF_Data.o SELF_Lagrange.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Mesh.F90 -o $@

SELF_Data.o: SELF_Lagrange.o SELF_Constants.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Data.F90 -o $@

SELF_Lagrange.o: SELF_Memory.o SELF_Quadrature.o SELF_Constants.o SELF_SupportRoutines.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Lagrange.F90 -o $@

SELF_Lagrange_HIP.o:
	${CXX} -c ${CXXFLAGS} src/SELF_Lagrange_HIP.cpp -o $@

SELF_Quadrature.o: SELF_Constants.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Quadrature.F90 -o $@

SELF_Memory.o:
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Memory.F90 -o $@

SELF_SupportRoutines.o: SELF_Constants.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_SupportRoutines.F90 -o $@

SELF_Constants.o:
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_Constants.F90 -o $@
