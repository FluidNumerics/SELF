
SHELL=/bin/bash
SELFDIR=$(shell pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')
INITDIR=$(shell pwd)

.PHONY: mappedgradient \
        MeshGen_3D 

MeshGen_3D :
	make --directory=${SELFDIR}/build/ MeshGen_3D
	mv ${SELFDIR}/build/MeshGen_3D ./

mappedgradient :
	@echo 'Found SELF installation directory : '${SELFDIR}
	make -C ${SELFDIR}/build/ mappedgradient
	mv ${SELFDIR}/build/mappedgradient ./

clean :
	make -C ${SELFDIR}/build/ clean

