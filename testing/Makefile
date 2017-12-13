
SHELL=/bin/bash
SELFDIR=$(shell pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')
INITDIR=$(shell pwd)

.PHONY: spectral_tests \
        clean

spectral_tests :
	make --directory=${SELFDIR}/build/ spectral_tests
	mv ${SELFDIR}/build/spectral_tests ./

clean :
	make --directory=${SELFDIR}/build/ clean
