
SHELL=/bin/bash
SELFDIR=$(shell pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')
INITDIR=$(shell pwd)

.PHONY: meshgen \
        initialize \
        integrate 

meshgen :
	make --directory=${SELFDIR}/build/ -e CUDA=no MPI=no meshgen
	make --directory=${SELFDIR}/build/ clean
	mv ${SELFDIR}/build/meshgen ./

initialize :
	@echo 'Found SELF installation directory : '${SELFDIR}
	@echo 'Using modifications in directory : '${INITDIR}
	make MODDIR=${INITDIR} -C ${SELFDIR}/build/ initialize
	mv ${SELFDIR}/build/initialize ./

integrate :
	@echo 'Found SELF installation directory : '${SELFDIR}
	make -C ${SELFDIR}/build/ integrate
	mv ${SELFDIR}/build/integrate ./

clean :
	make -C ${SELFDIR}/build/ clean

