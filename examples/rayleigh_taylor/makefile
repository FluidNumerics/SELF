
SHELL=/bin/bash
SELFDIR=$(shell pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')
INITDIR=$(shell pwd)

.PHONY: DecomposeStructuredHexMesh \
        GenerateHexMesh \
        initialize \
        integrate \
        test_stateExchange

meshgen :
	make --directory=${SELFDIR}/build/ CUDA=no MPI=no meshgen
	mv ${SELFDIR}/build/meshgen ./
	make --directory=${SELFDIR}/build/ clean


initialize :
	@echo 'Found SELF installation directory : '${SELFDIR}
	@echo 'Using modifications in directory : '${INITDIR}
	make -C ${SELFDIR}/build/ MODDIR=${INITDIR} initialize
	mv ${SELFDIR}/build/initialize ./

integrate :
	@echo 'Found SELF installation directory : '${SELFDIR}
	make -C ${SELFDIR}/build/ integrate
	mv ${SELFDIR}/build/integrate ./

clean :
	make -C ${SELFDIR}/build/ clean

