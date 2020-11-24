

FC=/opt/hipfort/bin/hipfc
CXX=/opt/rocm/bin/hipcc
FFLAGS=-DGPU -ffree-line-length-none
FLIBS=-L/apps/self/lib/ -lFLAP -lFACE -lPENF
INC=-I/apps/self/include/FLAP -I/apps/self/include/PENF -I/apps/self/include/FACE
CXXFLAGS=
PREFIX=/apps/self

install: libSELF.a self
	mkdir -p ${PREFIX}/bin
	mkdir -p ${PREFIX}/lib
	mkdir -p ${PREFIX}/include
	mv self ${PREFIX}/bin/
	mv libSELF.a ${PREFIX}/lib/
	mv *.mod ${PREFIX}/include/
	cp src/*.h ${PREFIX}/include/

self: SELF.o
	${FC} ${FFLAGS} ${INC} SELF.o ${FLIBS} -o $@

SELF.o:
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF.F90 -o $@

libSELF.a: SELF_Constants.o SELF_Data.o SELF_Lagrange_HIP.o SELF_Lagrange.o SELF_MappedData_HIP.o SELF_MappedData.o SELF_Memory.o SELF_Mesh.o SELF_Quadrature.o SELF_SupportRoutines.o
	ar rcs $@ SELF_Constants.o SELF_Data.o SELF_Lagrange_HIP.o SELF_Lagrange.o SELF_MappedData_HIP.o SELF_MappedData.o SELF_Memory.o SELF_Mesh.o SELF_Quadrature.o SELF_SupportRoutines.o 	

#SELF_MPILayer.o: SELF_MPILayer.F90

SELF_MappedData.o: SELF_Mesh.o SELF_Data.o
	${FC} -c ${FFLAGS} ${INC} ${FLIBS} src/SELF_MappedData.F90 -o $@

SELF_MappedData_HIP.o:
	${CXX} -c ${CXXFLAGS} src/SELF_MappedData_HIP.cpp -o $@

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
