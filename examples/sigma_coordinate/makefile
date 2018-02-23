
SHELL=/bin/bash
SELFDIR=$(shell pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')
INITDIR=$(shell pwd)

.PHONY: mappedgradient 

mappedgradient :
	@echo 'Found SELF installation directory : '${SELFDIR}
	make -C ${SELFDIR}/build/ mappedgradient
	mv ${SELFDIR}/build/mappedgradient ./

clean :
	make -C ${SELFDIR}/build/ clean

