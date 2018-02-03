
SHELL=/bin/bash
SELFDIR=$(shell pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')
INITDIR=$(shell pwd)

.PHONY: SMeshGen_3D \
        initialize \
        integrate \
        test_stateExchange

SMeshGen_3D :
	make --directory=${SELFDIR}/build/ SMeshGen_3D
	mv ${SELFDIR}/build/SMeshGen_3D ./

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

