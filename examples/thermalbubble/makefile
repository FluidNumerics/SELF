
SHELL=/bin/bash
SELFDIR=$(shell pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')
INITDIR=$(shell pwd)

.PHONY: DecomposeStructuredHexMesh \
        GenerateHexMesh \
        initialize \
        integrate \
        test_stateExchange

DecomposeStructuredHexMesh :
	make --directory=${SELFDIR}/build/ DecomposeStructuredHexMesh
	mv ${SELFDIR}/build/DecomposeStructuredHexMesh ./
	make --directory=${SELFDIR}/build/ clean

GenerateHexMesh :
	make --directory=${SELFDIR}/build/ GenerateHexMesh
	mv ${SELFDIR}/build/GenerateHexMesh ./
	make --directory=${SELFDIR}/build/ clean

initialize :
	@echo 'Found SELF installation directory : '${SELFDIR}
	@echo 'Using modifications in directory : '${INITDIR}
	make -C ${SELFDIR}/build/ initialize
	mv ${SELFDIR}/build/initialize ./

integrate :
	@echo 'Found SELF installation directory : '${SELFDIR}
	make -C ${SELFDIR}/build/ integrate
	mv ${SELFDIR}/build/integrate ./

test_stateExchange :
	@echo 'Found SELF installation directory : '${SELFDIR}
	make -C ${SELFDIR}/build/ test_stateExchange
	mv ${SELFDIR}/build/test_stateExchange ./

clean :
	make -C ${SELFDIR}/build/ clean

