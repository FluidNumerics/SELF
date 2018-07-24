
SHELL=/bin/bash
SELFDIR=$(shell pwd | sed 's/\SELF-Fluids.*/SELF-Fluids/')
INITDIR=$(shell pwd)

.PHONY: MeshGen_3D \
        sfluid

MeshGen_3D :
	make --directory=${SELFDIR}/build/ MeshGen_3D CUDA=no
	make clean
	mv ${SELFDIR}/build/MeshGen_3D ./

sfluid :
	@echo 'Found SELF installation directory : '${SELFDIR}
	source SELF_environment_settings
	make -C ${SELFDIR}/build/ sfluid
	mv ${SELFDIR}/build/sfluid ./

clean :
	make -C ${SELFDIR}/build/ clean

