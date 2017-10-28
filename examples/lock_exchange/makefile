
SHELL=/bin/bash
SELFDIR=$(shell pwd | sed 's/\SELF.*/SELF/')
INITDIR=$(shell pwd)

.PHONY: DecomposeStructuredHexMesh \
        GenerateHexMesh \
        initialize \
        integrate 

DecomposeStructuredHexMesh :
	make --directory=${SELFDIR}/build/ DecomposeStructuredHexMesh
	mv ${SELFDIR}/build/DecomposeStructuredHexMesh ../run/
	make --directory=${SELFDIR}/build/ clean

GenerateHexMesh :
	make --directory=${SELFDIR}/build/ GenerateHexMesh
	mv ${SELFDIR}/build/GenerateHexMesh ../run/
	make --directory=${SELFDIR}/build/ clean

initialize :
	@echo 'Found SELF installation directory : '${SELFDIR}
	@echo 'Using modifications in directory : '${INITDIR}
	make MODDIR=${INITDIR} -C ${SELFDIR}/build/ initialize
	mv ${SELFDIR}/build/initialize ../run/

integrate :
	@echo 'Found SELF installation directory : '${SELFDIR}
	make -C ${SELFDIR}/build/ integrate
	mv ${SELFDIR}/build/integrate ../run/

clean :
	make -C ${SELFDIR}/build/ clean

